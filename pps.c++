#include "Compression.h"
#include "G4Box.hh"
#include "G4EmParameters.hh"
#include "G4Event.hh"
#include "G4HadronicParameters.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4NuclearLevelData.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4PhysListFactory.hh"
#include "G4PVPlacement.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4RunManagerFactory.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4StepStatus.hh"
#include "G4SystemOfUnits.hh"
#include "G4Threading.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4UnitsTable.hh"
#include "G4UserEventAction.hh"
#include "G4UserRunAction.hh"
#include "G4UserSteppingAction.hh"
#include "G4UserTrackingAction.hh"
#include "G4VisExecutive.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VProcess.hh"
#include "G4VUserActionInitialization.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "Randomize.hh"
#include "ROOT/REntry.hxx"
#include "ROOT/RNTupleFillContext.hxx"
#include "ROOT/RNTupleFillStatus.hxx"
#include "ROOT/RNTupleModel.hxx"
#include "ROOT/RNTupleParallelWriter.hxx"
#include "ROOT/RNTupleWriteOptions.hxx"
#include "TFile.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

namespace PPS {

auto DefaultThreadCount() -> int {
    const auto coreCount = std::thread::hardware_concurrency();
    return coreCount > 0 ? static_cast<int>(coreCount) : 1;
}

struct Layer {
    std::string mMaterialName;
    double mThickness{0.0};
};

struct RadiationSource {
    std::string mParticleName;
    double mEnergy{0.0};
    double mIntensity{1.0};
    double mWeight{1.0};
};

struct Config {
    std::vector<Layer> mLayers;
    std::vector<RadiationSource> mSources;
    std::string mOutputFileName{"pps_output.root"};
    std::string mPhysicsListName{"QGSP_BIC_AllHP_EMZ"};
    int mThreads{DefaultThreadCount()};
    int mVerbose{0};
    bool mIncludeNeutrinos{false};
    int mEventCount{0};
    int mPrintProgress{1000};
    bool mUI{false};
    bool mOverwrite{false};
    bool mHelp{false};
    std::string mMacroFileName;

    auto TotalThickness() const -> double {
        auto total{0.0};
        for (const auto& layer : mLayers) {
            total += layer.mThickness;
        }
        return total;
    }
};

auto ParseQuantity(const std::string& text, const std::map<std::string, double>& units) -> double {
    auto stream = std::istringstream{text};
    auto value{0.0};
    auto suffix = std::string{};
    stream >> value;
    if (stream.fail()) {
        throw std::invalid_argument("cannot parse numeric value from '" + text + "'");
    }
    stream >> suffix;
    if (not suffix.empty()) {
        const auto unit = units.find(suffix);
        if (unit == units.end()) {
            throw std::invalid_argument("unknown unit '" + suffix + "' in '" + text + "'");
        }
        value *= unit->second;
    }
    // Reject NaN/inf and unit overflows (e.g. "1e308 MeV").
    if (not std::isfinite(value)) {
        throw std::invalid_argument("non-finite numeric value in '" + text + "'");
    }
    return value;
}

auto ParseLength(const std::string& text) -> double {
    static const auto units = std::map<std::string, double>{
        {"mm",       mm       },
        {"cm",       cm       },
        {"m",        meter    },
        {"km",       kilometer},
        {"um",       um       },
        {"nm",       nm       },
        {"angstrom", angstrom },
        {"fm",       fermi    }
    };
    return ParseQuantity(text, units);
}

auto ParseEnergy(const std::string& text) -> double {
    static const auto units = std::map<std::string, double>{
        {"eV",  eV },
        {"keV", keV},
        {"MeV", MeV},
        {"GeV", GeV},
        {"TeV", TeV},
        {"PeV", PeV}
    };
    return ParseQuantity(text, units);
}

auto Trim(const std::string& text) -> std::string {
    const auto first = text.find_first_not_of(" \t");
    if (first == std::string::npos) {
        return "";
    }
    const auto last = text.find_last_not_of(" \t");
    return text.substr(first, last - first + 1);
}

auto ParseLayers(const std::string& text) -> std::vector<Layer> {
    auto layers = std::vector<Layer>{};
    auto stream = std::istringstream{text};
    auto layerSpec = std::string{};
    while (std::getline(stream, layerSpec, ';')) {
        if (layerSpec.empty()) {
            throw std::invalid_argument("empty layer specification in '" + text + "'");
        }
        const auto colon = layerSpec.find(':');
        if (colon == std::string::npos) {
            throw std::invalid_argument(
                "layer '" + layerSpec + "' must have the form 'material:thickness', e.g. G4_Cu:10 cm");
        }
        auto layer = Layer{};
        layer.mMaterialName = Trim(layerSpec.substr(0, colon));
        layer.mThickness = ParseLength(Trim(layerSpec.substr(colon + 1)));
        if (layer.mMaterialName.empty()) {
            throw std::invalid_argument("layer '" + layerSpec + "' has an empty material name");
        }
        if (layer.mThickness <= 0.0) {
            throw std::invalid_argument("layer '" + layerSpec + "' has a non-positive thickness");
        }
        layers.push_back(std::move(layer));
    }
    if (layers.empty()) {
        throw std::invalid_argument("at least one material layer is required");
    }
    return layers;
}

auto ParseRadiationSources(const std::string& text) -> std::vector<RadiationSource> {
    auto sources = std::vector<RadiationSource>{};
    auto intensitySpecified = std::vector<bool>{};
    auto stream = std::istringstream{text};
    auto sourceSpec = std::string{};
    while (std::getline(stream, sourceSpec, ';')) {
        if (sourceSpec.empty()) {
            throw std::invalid_argument("empty source specification in '" + text + "'");
        }
        const auto firstColon = sourceSpec.find(':');
        if (firstColon == std::string::npos) {
            throw std::invalid_argument(
                "source '" + sourceSpec + "' must have the form 'particle:energy[:intensity]', e.g. e+:10 MeV:15");
        }
        const auto secondColon = sourceSpec.find(':', firstColon + 1);
        auto source = RadiationSource{};
        source.mParticleName = Trim(sourceSpec.substr(0, firstColon));
        source.mEnergy = ParseEnergy(Trim(sourceSpec.substr(firstColon + 1, secondColon - firstColon - 1)));
        if (source.mParticleName.empty()) {
            throw std::invalid_argument("source '" + sourceSpec + "' has an empty particle name");
        }
        if (source.mEnergy <= 0.0) {
            throw std::invalid_argument("source '" + sourceSpec + "' has a non-positive energy");
        }
        if (secondColon != std::string::npos) {
            const auto intensityText = Trim(sourceSpec.substr(secondColon + 1));
            auto parseIndex = std::size_t{0};
            try {
                source.mIntensity = std::stod(intensityText, &parseIndex);
            } catch (const std::exception&) {
                throw std::invalid_argument("cannot parse intensity '" + intensityText + "' in source '" + sourceSpec + "'");
            }
            if (parseIndex != intensityText.size()) {
                throw std::invalid_argument("cannot parse intensity '" + intensityText + "' in source '" + sourceSpec + "'");
            }
            if (source.mIntensity <= 0.0) {
                throw std::invalid_argument("source '" + sourceSpec + "' has a non-positive intensity");
            }
            if (not std::isfinite(source.mIntensity)) {
                throw std::invalid_argument("source '" + sourceSpec + "' has a non-finite intensity");
            }
            intensitySpecified.push_back(true);
        } else {
            intensitySpecified.push_back(false);
        }
        sources.push_back(std::move(source));
    }
    if (sources.empty()) {
        throw std::invalid_argument("at least one radiation source is required");
    }
    if (sources.size() > 1) {
        for (auto index{0}; index < static_cast<int>(sources.size()); ++index) {
            if (not intensitySpecified[index]) {
                throw std::invalid_argument(
                    "source '" + sources[index].mParticleName + "' is missing its relative intensity; "
                                                                "the intensity may be omitted only when the source list contains a single entry");
            }
        }
    }
    auto totalIntensity{0.0};
    for (const auto& source : sources) {
        totalIntensity += source.mIntensity;
    }
    for (auto& source : sources) {
        source.mWeight = source.mIntensity / totalIntensity;
    }
    return sources;
}

auto FormatEnergy(double energy) -> std::string {
    auto stream = std::ostringstream{};
    stream << G4BestUnit(energy, "Energy");
    return stream.str();
}

auto PrintUsage(const char* programName) -> void {
    G4cout
        << "usage: " << programName << " [options] [macroFile]\n"
        << "simulate one primary particle per event from a compound radiation source passing through a material\n"
        << "slab, recording per-event energy depositions inside the slab and particles penetrating or backscattering\n"
        << "into the world boundary; each source category is stored into its own RNTuple 'run{runId}_src{typeId}'\n"
        << "inside a ROOT file\n"
        << "\n"
        << "required options:\n"
        << "  -m, --material <spec>    material layers as semicolon-separated 'name:thickness' entries, built in\n"
        << "                           order from the source, e.g. \"G4_WATER:10 cm;G4_Cu:5 cm\"\n"
        << "  -s, --radiation-src <spec>\n"
        << "                           compound radiation source as semicolon-separated 'particle:energy[:intensity]'\n"
        << "                           entries; the intensity is the relative intensity of the source and may be\n"
        << "                           omitted only when the list contains a single source, e.g.\n"
        << "                           \"e+:10 MeV:15;gamma:3 MeV:12\" or \"neutron:1 MeV\"; intensities are\n"
        << "                           normalized to probabilities and every event emits one primary whose source\n"
        << "                           type is drawn from the resulting multinomial distribution\n"
        << "\n"
        << "optional options:\n"
        << "  -n, --n-event <count>    simulate <count> events in batch mode; may be combined with --ui to\n"
        << "                           pre-run events before the interactive session opens\n"
        << "  -i, --ui                 start an interactive UI session with visualization; if set, visualization\n"
        << "                           is enabled, otherwise the program runs without any UI\n"
        << "  -o, --output <file>      output ROOT file name (default: pps_output.root; never overwrites unless --force)\n"
        << "  -l, --physics <name>     reference physics list name (default: QGSP_BIC_AllHP_EMZ)\n"
        << "  -j, --threads <count>    worker thread count; 1 runs sequential, > 1 runs multithreaded (default: all CPU cores)\n"
        << "  -v, --verbose <level>    Geant4 verbosity level; 0 prints only the banner, progress, and summary (default: 0)\n"
        << "  -N, --neutrinos          include neutrino kinetic energy in the world boundary energy statistics;\n"
        << "                           neutrinos are ignored by default\n"
        << "  -f, --force              overwrite the output file if it already exists\n"
        << "  -h, --help               print this message\n"
        << "\n"
        << "a positional macroFile executes in batch mode and cannot be combined with --n-event or --ui;\n"
        << "without --ui, --n-event, or a macroFile the program refuses to start"
        << G4endl;
}

auto ParseCommandLine(int argc, char** argv) -> Config {
    auto config = Config{};
    auto positional = std::vector<std::string>{};
    auto nextValue = [&argv, argc](int& index, const std::string& option) -> std::string {
        if (index + 1 >= argc) {
            throw std::invalid_argument("missing value for option '" + option + "'");
        }
        ++index;
        return argv[index];
    };
    for (auto i{1}; i < argc; ++i) {
        const auto argument{std::string{argv[i]}};
        if (argument == "-h" or argument == "--help") {
            config.mHelp = true;
        } else if (argument == "-m" or argument == "--material") {
            config.mLayers = ParseLayers(nextValue(i, argument));
        } else if (argument == "-s" or argument == "--radiation-src") {
            config.mSources = ParseRadiationSources(nextValue(i, argument));
        } else if (argument == "-o" or argument == "--output") {
            config.mOutputFileName = nextValue(i, argument);
        } else if (argument == "-l" or argument == "--physics") {
            config.mPhysicsListName = nextValue(i, argument);
        } else if (argument == "-j" or argument == "--threads") {
            config.mThreads = std::stoi(nextValue(i, argument));
        } else if (argument == "-v" or argument == "--verbose") {
            config.mVerbose = std::stoi(nextValue(i, argument));
            if (config.mVerbose < 0) {
                throw std::invalid_argument("--verbose must be at least 0");
            }
        } else if (argument == "-n" or argument == "--n-event") {
            config.mEventCount = std::stoi(nextValue(i, argument));
            if (config.mEventCount < 1) {
                throw std::invalid_argument("--n-event must be at least 1");
            }
        } else if (argument == "-i" or argument == "--ui") {
            config.mUI = true;
        } else if (argument == "-N" or argument == "--neutrinos") {
            config.mIncludeNeutrinos = true;
        } else if (argument == "-f" or argument == "--force") {
            config.mOverwrite = true;
        } else if (not argument.empty() and argument.front() == '-') {
            throw std::invalid_argument("unknown option '" + argument + "'");
        } else {
            positional.push_back(argument);
        }
    }
    if (positional.size() > 1) {
        throw std::invalid_argument("at most one macro file may be given");
    }
    if (not positional.empty()) {
        config.mMacroFileName = positional.front();
    }
    if (not config.mHelp) {
        if (config.mLayers.empty()) {
            throw std::invalid_argument("required option --material is missing");
        }
        if (config.mSources.empty()) {
            throw std::invalid_argument("required option --radiation-src is missing");
        }
        if (config.mThreads < 1) {
            throw std::invalid_argument("--threads must be at least 1");
        }
        if (not config.mMacroFileName.empty() and (config.mUI or config.mEventCount > 0)) {
            throw std::invalid_argument("a macro file cannot be combined with --ui or --n-event");
        }
        if (config.mMacroFileName.empty() and not config.mUI and config.mEventCount == 0) {
            throw std::invalid_argument("specify --n-event <count>, --ui, or a macro file");
        }
    }
    return config;
}

auto ValidateMaterials(const Config& config) -> void {
    for (const auto& layer : config.mLayers) {
        if (G4NistManager::Instance()->FindOrBuildMaterial(layer.mMaterialName) == nullptr) {
            throw std::invalid_argument("material '" + layer.mMaterialName + "' not found in the NIST material database");
        }
    }
}

auto ValidateSources(const Config& config) -> void {
    for (const auto& source : config.mSources) {
        if (G4ParticleTable::GetParticleTable()->FindParticle(source.mParticleName) == nullptr) {
            throw std::invalid_argument("particle '" + source.mParticleName + "' not found in the particle table");
        }
    }
}

auto IsNeutrino(const G4ParticleDefinition* particle) -> bool {
    const auto pdg = std::abs(particle->GetPDGEncoding());
    return pdg == 12 or pdg == 14 or pdg == 16;
}

class OutputWriter {
public:
    static auto Instance() -> OutputWriter& {
        static auto instance = OutputWriter{};
        return instance;
    }

    OutputWriter(const OutputWriter&) = delete;
    auto operator=(const OutputWriter&) -> OutputWriter& = delete;

    auto Initialize(const std::string& fileName, bool overwrite) -> void {
        if (mInitialized) {
            return;
        }
        const auto exists = std::filesystem::exists(fileName);
        if (exists and not overwrite) {
            G4cerr << "error: output file '" << fileName
                   << "' already exists; pass --force to overwrite it" << G4endl;
            std::quick_exit(EXIT_FAILURE);
        }
        mFile = std::unique_ptr<TFile>{TFile::Open(fileName.c_str(), exists ? "RECREATE" : "UPDATE")};
        if (mFile == nullptr or mFile->IsZombie()) {
            G4cerr << "error: failed to open output file '" << fileName << "'" << G4endl;
            std::quick_exit(EXIT_FAILURE);
        }
        mInitialized = true;
    }

    auto BeginRun(int runId, int layerCount, int sourceCount) -> void {
        try {
            mWriters.clear();
            mWriters.reserve(static_cast<std::size_t>(sourceCount));
            auto options = ROOT::RNTupleWriteOptions{};
            options.SetCompression(ROOT::RCompressionSetting::EDefaults::kUseGeneralPurpose);
            options.SetUseImplicitMT(ROOT::RNTupleWriteOptions::EImplicitMT::kOff);
            for (auto type{0}; type < sourceCount; ++type) {
                auto model = CreateModel(layerCount);
                const auto name = "run" + std::to_string(runId) + "_src" + std::to_string(type);
                mWriters.push_back(ROOT::RNTupleParallelWriter::Append(std::move(model), name, *mFile, options));
            }
        } catch (const std::exception& error) {
            G4cerr << "error: failed to create RNTuples for run " << runId << ": " << error.what() << G4endl;
            std::quick_exit(EXIT_FAILURE);
        }
    }

    auto CreateFillContext(int sourceIndex) -> std::shared_ptr<ROOT::RNTupleFillContext> {
        return mWriters[sourceIndex]->CreateFillContext();
    }

    auto EndRun() -> void {
        for (auto& writer : mWriters) {
            if (writer != nullptr) {
                writer->CommitDataset();
            }
        }
        mWriters.clear();
    }

    auto Finalize() -> void {
        EndRun();
        if (mFile != nullptr) {
            mFile->Close();
            mFile.reset();
        }
        mInitialized = false;
    }

    // RNTupleParallelWriter instances that share one TFile are not synchronized with each other,
    // so every actual file write (cluster flush) must be serialized through this mutex.
    auto FileAccessMutex() -> std::mutex& {
        return mFileAccessMutex;
    }

private:
    auto CreateModel(int layerCount) -> std::unique_ptr<ROOT::RNTupleModel> {
        auto model = ROOT::RNTupleModel::CreateBare();
        model->MakeField<int>("event_id");
        model->MakeField<float>("total_e_pen");
        model->MakeField<std::vector<std::string>>("particle_pen");
        model->MakeField<std::vector<float>>("theta_pen");
        model->MakeField<std::vector<float>>("phi_pen");
        model->MakeField<std::vector<float>>("e_pen");
        model->MakeField<float>("total_e_dep");
        for (auto layerIndex{0}; layerIndex < layerCount; ++layerIndex) {
            const auto suffix = std::to_string(layerIndex);
            model->MakeField<float>("e_dep_" + suffix);
            model->MakeField<std::vector<std::string>>("particle_dep_" + suffix);
            model->MakeField<std::vector<float>>("x_dep_" + suffix);
            model->MakeField<std::vector<float>>("y_dep_" + suffix);
            model->MakeField<std::vector<float>>("z_dep_" + suffix);
            model->MakeField<std::vector<float>>("w_dep_" + suffix);
            model->MakeField<std::vector<std::string>>("proc_dep_" + suffix);
        }
        model->MakeField<float>("total_e_bsc");
        model->MakeField<std::vector<std::string>>("particle_bsc");
        model->MakeField<std::vector<float>>("theta_bsc");
        model->MakeField<std::vector<float>>("phi_bsc");
        model->MakeField<std::vector<float>>("e_bsc");
        return model;
    }

    OutputWriter() = default;

    std::unique_ptr<TFile> mFile;
    std::vector<std::unique_ptr<ROOT::RNTupleParallelWriter>> mWriters;
    std::mutex mFileAccessMutex;
    bool mInitialized{false};
};

class DetectorConstruction : public G4VUserDetectorConstruction {
public:
    explicit DetectorConstruction(const Config& config) :
        mConfig{config} {}
    ~DetectorConstruction() override = default;

    auto Construct() -> G4VPhysicalVolume* override {
        const auto nistManager = G4NistManager::Instance();
        const auto hydrogen = nistManager->FindOrBuildElement("H");
        const auto vacuumMaterial =
            new G4Material{"vacuum", 1e-16 * g / cm3, 1, kStateGas, 293.15 * kelvin, atmosphere};
        vacuumMaterial->AddElement(hydrogen, 1);

        const auto totalThickness = mConfig.TotalThickness();
        const auto xyHalfLength{500.0 * totalThickness};
        const auto worldZHalfLength{1.1 * totalThickness};

        const auto solidWorld = new G4Box{"World", xyHalfLength, xyHalfLength, worldZHalfLength};
        const auto logicalWorld = new G4LogicalVolume{solidWorld, vacuumMaterial, "World"};
        const auto physicalWorld = new G4PVPlacement{
            nullptr, G4ThreeVector{}, logicalWorld, "World", nullptr, false, 0,
            mConfig.mVerbose > 0};

        auto layerIndex{0};
        auto zPosition{0.0};
        for (const auto& layer : mConfig.mLayers) {
            const auto layerHalfLength = 0.5 * layer.mThickness;
            const auto solidLayer = new G4Box{"LayerSolid", xyHalfLength, xyHalfLength, layerHalfLength};
            const auto layerMaterial = nistManager->FindOrBuildMaterial(layer.mMaterialName);
            const auto logicalLayer = new G4LogicalVolume{solidLayer, layerMaterial, "Layer"};
            new G4PVPlacement{
                nullptr, G4ThreeVector{0.0, 0.0, zPosition + layerHalfLength},
                logicalLayer, "Layer",
                logicalWorld, false, layerIndex, mConfig.mVerbose > 0
            };
            mMaterialVolumes.emplace(logicalLayer, layerIndex);
            zPosition += layer.mThickness;
            ++layerIndex;
        }

        return physicalWorld;
    }

    auto GetMaterialVolumes() const -> const std::unordered_map<const G4LogicalVolume*, int>& {
        return mMaterialVolumes;
    }

private:
    const Config& mConfig;
    std::unordered_map<const G4LogicalVolume*, int> mMaterialVolumes;
};

class SimulationRun : public G4Run {
public:
    explicit SimulationRun(int layerCount, int sourceCount) :
        mSourceStatistics(sourceCount, SourceStatistics{layerCount}),
        mLayerCount{layerCount} {}

    auto Merge(const G4Run* other) -> void override {
        G4Run::Merge(other);
        const SimulationRun* otherSimulationRun = static_cast<const SimulationRun*>(other);
        for (auto sourceIndex{0}; sourceIndex < static_cast<int>(mSourceStatistics.size()); ++sourceIndex) {
            auto& target = mSourceStatistics[sourceIndex];
            const auto& source = otherSimulationRun->mSourceStatistics[sourceIndex];
            target.mEventCount += source.mEventCount;
            target.mPenetrationEnergySum += source.mPenetrationEnergySum;
            target.mPenetrationEnergySumSq += source.mPenetrationEnergySumSq;
            target.mDepositionEnergySum += source.mDepositionEnergySum;
            target.mDepositionEnergySumSq += source.mDepositionEnergySumSq;
            for (auto layerIndex{0}; layerIndex < static_cast<int>(target.mLayerDepositionEnergySums.size()); ++layerIndex) {
                target.mLayerDepositionEnergySums[layerIndex] += source.mLayerDepositionEnergySums[layerIndex];
                target.mLayerDepositionEnergySumsSq[layerIndex] += source.mLayerDepositionEnergySumsSq[layerIndex];
            }
            target.mBackscatteringEnergySum += source.mBackscatteringEnergySum;
            target.mBackscatteringEnergySumSq += source.mBackscatteringEnergySumSq;
            AddToMap(target.mPenetrationParticleCount, source.mPenetrationParticleCount);
            AddToMap(target.mPenetrationParticleEnergySum, source.mPenetrationParticleEnergySum);
            AddToMap(target.mPenetrationParticleEnergySumSq, source.mPenetrationParticleEnergySumSq);
            AddToMap(target.mBackscatteringParticleCount, source.mBackscatteringParticleCount);
            AddToMap(target.mBackscatteringParticleEnergySum, source.mBackscatteringParticleEnergySum);
            AddToMap(target.mBackscatteringParticleEnergySumSq, source.mBackscatteringParticleEnergySumSq);
        }
    }

    auto AddEventResult(int sourceIndex, double penetrationEnergy, const std::vector<float>& layerDepositionEnergies,
                        double backscatteringEnergy) -> void {
        auto& statistics = mSourceStatistics[sourceIndex];
        ++statistics.mEventCount;
        statistics.mPenetrationEnergySum += penetrationEnergy;
        statistics.mPenetrationEnergySumSq += penetrationEnergy * penetrationEnergy;
        auto totalDepositionEnergy{0.0};
        for (auto layerIndex{0}; layerIndex < static_cast<int>(layerDepositionEnergies.size()); ++layerIndex) {
            const auto layerDepositionEnergy = layerDepositionEnergies[layerIndex];
            statistics.mLayerDepositionEnergySums[layerIndex] += layerDepositionEnergy;
            statistics.mLayerDepositionEnergySumsSq[layerIndex] += layerDepositionEnergy * layerDepositionEnergy;
            totalDepositionEnergy += layerDepositionEnergy;
        }
        statistics.mDepositionEnergySum += totalDepositionEnergy;
        statistics.mDepositionEnergySumSq += totalDepositionEnergy * totalDepositionEnergy;
        statistics.mBackscatteringEnergySum += backscatteringEnergy;
        statistics.mBackscatteringEnergySumSq += backscatteringEnergy * backscatteringEnergy;
    }

    auto AddPenetrationEnergy(int sourceIndex, const std::string& particleName, double energy) -> void {
        auto& statistics = mSourceStatistics[sourceIndex];
        statistics.mPenetrationParticleCount[particleName] += 1.0;
        statistics.mPenetrationParticleEnergySum[particleName] += energy;
        statistics.mPenetrationParticleEnergySumSq[particleName] += energy * energy;
    }

    auto AddBackscatteringEnergy(int sourceIndex, const std::string& particleName, double energy) -> void {
        auto& statistics = mSourceStatistics[sourceIndex];
        statistics.mBackscatteringParticleCount[particleName] += 1.0;
        statistics.mBackscatteringParticleEnergySum[particleName] += energy;
        statistics.mBackscatteringParticleEnergySumSq[particleName] += energy * energy;
    }

    auto GetEventCount(int sourceIndex) const -> long long {
        return mSourceStatistics[sourceIndex].mEventCount;
    }
    auto GetPenetrationEnergySum(int sourceIndex) const -> double {
        return mSourceStatistics[sourceIndex].mPenetrationEnergySum;
    }
    auto GetPenetrationEnergySumSq(int sourceIndex) const -> double {
        return mSourceStatistics[sourceIndex].mPenetrationEnergySumSq;
    }
    auto GetDepositionEnergySum(int sourceIndex) const -> double {
        return mSourceStatistics[sourceIndex].mDepositionEnergySum;
    }
    auto GetDepositionEnergySumSq(int sourceIndex) const -> double {
        return mSourceStatistics[sourceIndex].mDepositionEnergySumSq;
    }
    auto GetLayerDepositionEnergySum(int sourceIndex, int layerIndex) const -> double {
        return mSourceStatistics[sourceIndex].mLayerDepositionEnergySums[layerIndex];
    }
    auto GetLayerDepositionEnergySumSq(int sourceIndex, int layerIndex) const -> double {
        return mSourceStatistics[sourceIndex].mLayerDepositionEnergySumsSq[layerIndex];
    }
    auto GetLayerCount() const -> int {
        return mLayerCount;
    }
    auto GetSourceCount() const -> int {
        return static_cast<int>(mSourceStatistics.size());
    }
    auto GetBackscatteringEnergySum(int sourceIndex) const -> double {
        return mSourceStatistics[sourceIndex].mBackscatteringEnergySum;
    }
    auto GetBackscatteringEnergySumSq(int sourceIndex) const -> double {
        return mSourceStatistics[sourceIndex].mBackscatteringEnergySumSq;
    }
    auto GetPenetrationParticleCounts(int sourceIndex) const -> const std::map<std::string, double>& {
        return mSourceStatistics[sourceIndex].mPenetrationParticleCount;
    }
    auto GetPenetrationParticleEnergySums(int sourceIndex) const -> const std::map<std::string, double>& {
        return mSourceStatistics[sourceIndex].mPenetrationParticleEnergySum;
    }
    auto GetPenetrationParticleEnergySumsSq(int sourceIndex) const -> const std::map<std::string, double>& {
        return mSourceStatistics[sourceIndex].mPenetrationParticleEnergySumSq;
    }
    auto GetBackscatteringParticleCounts(int sourceIndex) const -> const std::map<std::string, double>& {
        return mSourceStatistics[sourceIndex].mBackscatteringParticleCount;
    }
    auto GetBackscatteringParticleEnergySums(int sourceIndex) const -> const std::map<std::string, double>& {
        return mSourceStatistics[sourceIndex].mBackscatteringParticleEnergySum;
    }
    auto GetBackscatteringParticleEnergySumsSq(int sourceIndex) const -> const std::map<std::string, double>& {
        return mSourceStatistics[sourceIndex].mBackscatteringParticleEnergySumSq;
    }

private:
    struct SourceStatistics {
        explicit SourceStatistics(int layerCount) :
            mLayerDepositionEnergySums(layerCount, 0.0),
            mLayerDepositionEnergySumsSq(layerCount, 0.0) {}
        long long mEventCount{0};
        double mPenetrationEnergySum{0.0};
        double mPenetrationEnergySumSq{0.0};
        double mDepositionEnergySum{0.0};
        double mDepositionEnergySumSq{0.0};
        std::vector<double> mLayerDepositionEnergySums;
        std::vector<double> mLayerDepositionEnergySumsSq;
        double mBackscatteringEnergySum{0.0};
        double mBackscatteringEnergySumSq{0.0};
        std::map<std::string, double> mPenetrationParticleCount;
        std::map<std::string, double> mPenetrationParticleEnergySum;
        std::map<std::string, double> mPenetrationParticleEnergySumSq;
        std::map<std::string, double> mBackscatteringParticleCount;
        std::map<std::string, double> mBackscatteringParticleEnergySum;
        std::map<std::string, double> mBackscatteringParticleEnergySumSq;
    };

    auto AddToMap(std::map<std::string, double>& target, const std::map<std::string, double>& source) -> void {
        for (const auto& [key, value] : source) {
            target[key] += value;
        }
    }

    std::vector<SourceStatistics> mSourceStatistics;
    int mLayerCount{0};
};

class RunAction : public G4UserRunAction {
public:
    explicit RunAction(const Config& config) :
        mConfig{config} {}
    ~RunAction() override = default;

    auto GenerateRun() -> G4Run* override {
        SimulationRun* run = new SimulationRun{static_cast<int>(mConfig.mLayers.size()),
                                               static_cast<int>(mConfig.mSources.size())};
        mCurrentRun = run;
        return run;
    }

    auto BeginOfRunAction(const G4Run* run) -> void override {
        if (IsMaster()) {
            OutputWriter::Instance().BeginRun(run->GetRunID(), static_cast<int>(mConfig.mLayers.size()),
                                              static_cast<int>(mConfig.mSources.size()));
        }
        if (not IsMaster() or not G4Threading::IsMultithreadedApplication()) {
            mFillContexts.clear();
            mFillContexts.reserve(mConfig.mSources.size());
            for (auto sourceIndex{0}; sourceIndex < static_cast<int>(mConfig.mSources.size()); ++sourceIndex) {
                mFillContexts.push_back(OutputWriter::Instance().CreateFillContext(sourceIndex));
            }
        }
    }

    auto EndOfRunAction(const G4Run* run) -> void override {
        FlushAndReleaseFillState();
        if (IsMaster()) {
            OutputWriter::Instance().EndRun();
            PrintStatistics(static_cast<const SimulationRun&>(*run));
            G4cout << "run " << run->GetRunID() << " finished with " << run->GetNumberOfEvent() << " events" << G4endl;
        }
    }

    auto AddEventResult(int sourceIndex, double penetrationEnergy, const std::vector<float>& layerDepositionEnergies,
                        double backscatteringEnergy) -> void {
        mCurrentRun->AddEventResult(sourceIndex, penetrationEnergy, layerDepositionEnergies, backscatteringEnergy);
    }

    auto AddPenetrationEnergy(int sourceIndex, const std::string& particleName, double energy) -> void {
        mCurrentRun->AddPenetrationEnergy(sourceIndex, particleName, energy);
    }

    auto AddBackscatteringEnergy(int sourceIndex, const std::string& particleName, double energy) -> void {
        mCurrentRun->AddBackscatteringEnergy(sourceIndex, particleName, energy);
    }

    auto GetFillContext(int sourceIndex) const -> std::shared_ptr<ROOT::RNTupleFillContext> {
        return mFillContexts[sourceIndex];
    }

private:
    auto FlushAndReleaseFillState() -> void {
        if (mFillContexts.empty()) {
            return;
        }
        // Parallel writers that share one TFile are not synchronized with each other, so all actual
        // file writes must be serialized. Flush every remaining cluster under the file-access mutex
        // first; afterwards the fill-context destructors become no-ops (FlushCluster of an empty
        // context returns immediately) and no longer touch the file.
        {
            std::lock_guard<std::mutex> guard{OutputWriter::Instance().FileAccessMutex()};
            for (const auto& fillContext : mFillContexts) {
                fillContext->FlushCluster();
            }
        }
        mFillContexts.clear();
    }

    auto PrintStatistics(const SimulationRun& run) -> void {
        const auto eventCount = run.GetNumberOfEvent();
        if (eventCount < 1) {
            return;
        }
        G4cout << '\n';
        G4cout << "===============================================================================\n";
        auto printedAnySource{false};
        for (auto sourceIndex{0}; sourceIndex < run.GetSourceCount(); ++sourceIndex) {
            const auto& source = mConfig.mSources[sourceIndex];
            const auto sourceEventCount = run.GetEventCount(sourceIndex);
            if (sourceEventCount < 1) {
                continue;
            }
            if (printedAnySource) {
                // Separator between consecutive sources.
                G4cout << "-------------------------------------------------------------------------------\n";
            }
            printedAnySource = true;
            G4cout << " source " << sourceIndex << ": " << source.mParticleName << ' '
                   << G4BestUnit(source.mEnergy, "Energy") << ", intensity "
                   << std::fixed << std::setprecision(4) << source.mWeight << ", "
                   << sourceEventCount << " events:\n";
            G4cout << std::defaultfloat << std::setprecision(6);
            PrintEnergyRatio("penetration ratio (pen)", run.GetPenetrationEnergySum(sourceIndex),
                             run.GetPenetrationEnergySumSq(sourceIndex), sourceEventCount, source.mEnergy);
            PrintEnergyRatio("deposition ratio (dep)", run.GetDepositionEnergySum(sourceIndex),
                             run.GetDepositionEnergySumSq(sourceIndex), sourceEventCount, source.mEnergy);
            for (auto layerIndex{0}; layerIndex < run.GetLayerCount(); ++layerIndex) {
                const auto label = "  in layer " + std::to_string(layerIndex) + " (" +
                                   mConfig.mLayers[layerIndex].mMaterialName + ')';
                PrintEnergyRatio(label, run.GetLayerDepositionEnergySum(sourceIndex, layerIndex),
                                 run.GetLayerDepositionEnergySumSq(sourceIndex, layerIndex), sourceEventCount,
                                 source.mEnergy);
            }
            PrintEnergyRatio("back-scattering ratio (bsc)", run.GetBackscatteringEnergySum(sourceIndex),
                             run.GetBackscatteringEnergySumSq(sourceIndex), sourceEventCount, source.mEnergy);
            PrintParticleStatistics("penetration particles", run.GetPenetrationParticleCounts(sourceIndex),
                                    run.GetPenetrationParticleEnergySums(sourceIndex),
                                    run.GetPenetrationParticleEnergySumsSq(sourceIndex));
            PrintParticleStatistics("back-scattering particles", run.GetBackscatteringParticleCounts(sourceIndex),
                                    run.GetBackscatteringParticleEnergySums(sourceIndex),
                                    run.GetBackscatteringParticleEnergySumsSq(sourceIndex));
        }
        G4cout << "===============================================================================\n"
               << G4endl;
    }

    auto PrintParticleStatistics(const std::string& title, const std::map<std::string, double>& particleCounts,
                                 const std::map<std::string, double>& energySums,
                                 const std::map<std::string, double>& energySumsSq) -> void {
        if (particleCounts.empty()) {
            return;
        }
        G4cout << "   " << title << ":\n";
        G4cout << "     "
               << std::setw(14) << "particle"
               << std::setw(20) << "<E>"
               << std::setw(20) << "rms"
               << std::setw(8) << "n" << '\n';
        for (const auto& [particleName, count] : particleCounts) {
            const auto energySum = energySums.find(particleName)->second;
            const auto energySumSq = energySumsSq.find(particleName)->second;
            const auto meanEnergy = energySum / count;
            const auto variance = energySumSq / count - meanEnergy * meanEnergy;
            const auto rmsEnergy = variance > 0.0 ? std::sqrt(variance) : 0.0;
            G4cout << "     "
                   << std::setw(14) << particleName
                   << std::setw(20) << FormatEnergy(meanEnergy)
                   << std::setw(20) << FormatEnergy(rmsEnergy)
                   << std::setw(8) << static_cast<long long>(count) << '\n';
        }
    }

    auto PrintEnergyRatio(const std::string& label, double sumEnergy, double sumEnergySq, G4int eventCount,
                          double incidentEnergy) -> void {
        const auto meanEnergy{sumEnergy / eventCount};
        const auto ratio{meanEnergy / incidentEnergy};
        // Clamp to non-negative: with eventCount < 2 or exact cancellation the sample variance is
        // zero (or slightly negative from float rounding) and sqrt of a negative value would be NaN.
        const auto variance{std::max(0.0, (sumEnergySq - sumEnergy * sumEnergy / eventCount) / (eventCount - 1))};
        const auto ratioError{std::sqrt(variance / eventCount) / incidentEnergy};
        G4cout << "   " << std::left << std::setw(45) << label << '(' << 100.0 * ratio << " +/- "
               << 100.0 * ratioError << ") %" << G4endl;
    }

    const Config& mConfig;
    SimulationRun* mCurrentRun{nullptr};
    std::vector<std::shared_ptr<ROOT::RNTupleFillContext>> mFillContexts;
};

class EventAction : public G4UserEventAction {
public:
    explicit EventAction(RunAction& runAction, int layerCount, int sourceCount) :
        mRunAction{runAction},
        mLayerCount{layerCount},
        mSourceCount{sourceCount},
        mFillContexts(sourceCount),
        mEntries(sourceCount),
        mFillStatuses(sourceCount),
        mFields(sourceCount),
        mLayerEnergyDeposit(layerCount, 0.0F),
        mDepositionParticles(layerCount),
        mDepositionX(layerCount),
        mDepositionY(layerCount),
        mDepositionZ(layerCount),
        mDepositionWeight(layerCount),
        mDepositionProcess(layerCount) {}
    ~EventAction() override = default;

    auto SetSourceIndex(int sourceIndex) -> void {
        mCurrentSourceIndex = sourceIndex;
    }

    auto BeginOfEventAction(const G4Event*) -> void override {
        const auto runId = G4RunManager::GetRunManager()->GetCurrentRun()->GetRunID();
        if (runId != mRunId) {
            mRunId = runId;
            for (auto sourceIndex{0}; sourceIndex < mSourceCount; ++sourceIndex) {
                // Weak references on purpose: the fill contexts are owned by RunAction and must be
                // destroyed by the worker's EndOfRunAction so that RNTupleParallelWriter::CommitDataset()
                // (which requires all contexts to be gone) can run on the master. Holding shared_ptrs
                // here would keep the contexts alive and make CommitDataset throw.
                const auto fillContext{mRunAction.GetFillContext(sourceIndex)};
                mFillContexts[sourceIndex] = fillContext;
                mEntries[sourceIndex] = fillContext->CreateEntry();
                BindFields(*mEntries[sourceIndex], mFields[sourceIndex]);
            }
        }
        mTotalPenetrationEnergy = 0.0F;
        mPenetrationParticles.clear();
        mPenetrationTheta.clear();
        mPenetrationPhi.clear();
        mPenetrationEnergy.clear();
        mTotalDepositionEnergy = 0.0F;
        for (auto& layerEnergyDeposit : mLayerEnergyDeposit) {
            layerEnergyDeposit = 0.0F;
        }
        for (auto& depositionParticles : mDepositionParticles) {
            depositionParticles.clear();
        }
        for (auto& depositionX : mDepositionX) {
            depositionX.clear();
        }
        for (auto& depositionY : mDepositionY) {
            depositionY.clear();
        }
        for (auto& depositionZ : mDepositionZ) {
            depositionZ.clear();
        }
        for (auto& depositionWeight : mDepositionWeight) {
            depositionWeight.clear();
        }
        for (auto& depositionProcess : mDepositionProcess) {
            depositionProcess.clear();
        }
        mTotalBackscatteringEnergy = 0.0F;
        mBackscatteringParticles.clear();
        mBackscatteringTheta.clear();
        mBackscatteringPhi.clear();
        mBackscatteringEnergy.clear();
    }

    auto EndOfEventAction(const G4Event* event) -> void override {
        const auto sourceIndex = mCurrentSourceIndex;
        auto& fields = mFields[sourceIndex];
        *fields.mEventId = event->GetEventID();
        *fields.mTotalPenetrationEnergy = mTotalPenetrationEnergy;
        *fields.mPenetrationParticles = mPenetrationParticles;
        *fields.mPenetrationTheta = mPenetrationTheta;
        *fields.mPenetrationPhi = mPenetrationPhi;
        *fields.mPenetrationEnergy = mPenetrationEnergy;
        *fields.mTotalDepositionEnergy = mTotalDepositionEnergy;
        for (auto layerIndex{0}; layerIndex < mLayerCount; ++layerIndex) {
            *fields.mLayerEnergyDeposit[layerIndex] = mLayerEnergyDeposit[layerIndex];
            *fields.mDepositionParticles[layerIndex] = mDepositionParticles[layerIndex];
            *fields.mDepositionX[layerIndex] = mDepositionX[layerIndex];
            *fields.mDepositionY[layerIndex] = mDepositionY[layerIndex];
            *fields.mDepositionZ[layerIndex] = mDepositionZ[layerIndex];
            *fields.mDepositionWeight[layerIndex] = mDepositionWeight[layerIndex];
            *fields.mDepositionProcess[layerIndex] = mDepositionProcess[layerIndex];
        }
        *fields.mTotalBackscatteringEnergy = mTotalBackscatteringEnergy;
        *fields.mBackscatteringParticles = mBackscatteringParticles;
        *fields.mBackscatteringTheta = mBackscatteringTheta;
        *fields.mBackscatteringPhi = mBackscatteringPhi;
        *fields.mBackscatteringEnergy = mBackscatteringEnergy;
        // Fill the entry into the RNTuple of the source category drawn for this event. Filling only
        // buffers data in memory; the actual file write happens in the explicit FlushCluster call,
        // which is serialized across all parallel writers through the file-access mutex.
        const auto fillContext{mFillContexts[sourceIndex].lock()};
        if (fillContext == nullptr) {
            G4cerr << "error: RNTuple fill context for source " << sourceIndex << " is no longer available" << G4endl;
        } else {
            fillContext->FillNoFlush(*mEntries[sourceIndex], mFillStatuses[sourceIndex]);
            if (mFillStatuses[sourceIndex].ShouldFlushCluster()) {
                std::lock_guard<std::mutex> guard{OutputWriter::Instance().FileAccessMutex()};
                fillContext->FlushCluster();
            }
        }
        mRunAction.AddEventResult(sourceIndex, mTotalPenetrationEnergy, mLayerEnergyDeposit, mTotalBackscatteringEnergy);
    }

    auto AddPenetration(const std::string& particleName, const G4ThreeVector& direction, float energy) -> void {
        AddExitPoint(mPenetrationParticles, mPenetrationTheta, mPenetrationPhi, mPenetrationEnergy,
                     mTotalPenetrationEnergy, particleName, direction, energy);
        mRunAction.AddPenetrationEnergy(mCurrentSourceIndex, particleName, energy);
    }

    auto AddDeposition(int layerIndex, const std::string& particleName, const G4ThreeVector& position,
                       float energyDeposit, const std::string& processName) -> void {
        mTotalDepositionEnergy += energyDeposit;
        mLayerEnergyDeposit[layerIndex] += energyDeposit;
        mDepositionParticles[layerIndex].push_back(particleName);
        mDepositionX[layerIndex].push_back(static_cast<float>(position.x()));
        mDepositionY[layerIndex].push_back(static_cast<float>(position.y()));
        mDepositionZ[layerIndex].push_back(static_cast<float>(position.z()));
        mDepositionWeight[layerIndex].push_back(energyDeposit);
        mDepositionProcess[layerIndex].push_back(processName);
    }

    auto AddBackscattering(const std::string& particleName, const G4ThreeVector& direction, float energy) -> void {
        AddExitPoint(mBackscatteringParticles, mBackscatteringTheta, mBackscatteringPhi, mBackscatteringEnergy,
                     mTotalBackscatteringEnergy, particleName, direction, energy);
        mRunAction.AddBackscatteringEnergy(mCurrentSourceIndex, particleName, energy);
    }

private:
    struct SourceFields {
        std::shared_ptr<int> mEventId;
        std::shared_ptr<float> mTotalPenetrationEnergy;
        std::shared_ptr<std::vector<std::string>> mPenetrationParticles;
        std::shared_ptr<std::vector<float>> mPenetrationTheta;
        std::shared_ptr<std::vector<float>> mPenetrationPhi;
        std::shared_ptr<std::vector<float>> mPenetrationEnergy;
        std::shared_ptr<float> mTotalDepositionEnergy;
        std::vector<std::shared_ptr<float>> mLayerEnergyDeposit;
        std::vector<std::shared_ptr<std::vector<std::string>>> mDepositionParticles;
        std::vector<std::shared_ptr<std::vector<float>>> mDepositionX;
        std::vector<std::shared_ptr<std::vector<float>>> mDepositionY;
        std::vector<std::shared_ptr<std::vector<float>>> mDepositionZ;
        std::vector<std::shared_ptr<std::vector<float>>> mDepositionWeight;
        std::vector<std::shared_ptr<std::vector<std::string>>> mDepositionProcess;
        std::shared_ptr<float> mTotalBackscatteringEnergy;
        std::shared_ptr<std::vector<std::string>> mBackscatteringParticles;
        std::shared_ptr<std::vector<float>> mBackscatteringTheta;
        std::shared_ptr<std::vector<float>> mBackscatteringPhi;
        std::shared_ptr<std::vector<float>> mBackscatteringEnergy;
    };

    auto BindFields(ROOT::REntry& entry, SourceFields& fields) -> void {
        fields.mEventId = entry.GetPtr<int>("event_id");
        fields.mTotalPenetrationEnergy = entry.GetPtr<float>("total_e_pen");
        fields.mPenetrationParticles = entry.GetPtr<std::vector<std::string>>("particle_pen");
        fields.mPenetrationTheta = entry.GetPtr<std::vector<float>>("theta_pen");
        fields.mPenetrationPhi = entry.GetPtr<std::vector<float>>("phi_pen");
        fields.mPenetrationEnergy = entry.GetPtr<std::vector<float>>("e_pen");
        fields.mTotalDepositionEnergy = entry.GetPtr<float>("total_e_dep");
        fields.mLayerEnergyDeposit.clear();
        fields.mDepositionParticles.clear();
        fields.mDepositionX.clear();
        fields.mDepositionY.clear();
        fields.mDepositionZ.clear();
        fields.mDepositionWeight.clear();
        fields.mDepositionProcess.clear();
        for (auto layerIndex{0}; layerIndex < mLayerCount; ++layerIndex) {
            const auto suffix = std::to_string(layerIndex);
            fields.mLayerEnergyDeposit.push_back(entry.GetPtr<float>("e_dep_" + suffix));
            fields.mDepositionParticles.push_back(entry.GetPtr<std::vector<std::string>>("particle_dep_" + suffix));
            fields.mDepositionX.push_back(entry.GetPtr<std::vector<float>>("x_dep_" + suffix));
            fields.mDepositionY.push_back(entry.GetPtr<std::vector<float>>("y_dep_" + suffix));
            fields.mDepositionZ.push_back(entry.GetPtr<std::vector<float>>("z_dep_" + suffix));
            fields.mDepositionWeight.push_back(entry.GetPtr<std::vector<float>>("w_dep_" + suffix));
            fields.mDepositionProcess.push_back(entry.GetPtr<std::vector<std::string>>("proc_dep_" + suffix));
        }
        fields.mTotalBackscatteringEnergy = entry.GetPtr<float>("total_e_bsc");
        fields.mBackscatteringParticles = entry.GetPtr<std::vector<std::string>>("particle_bsc");
        fields.mBackscatteringTheta = entry.GetPtr<std::vector<float>>("theta_bsc");
        fields.mBackscatteringPhi = entry.GetPtr<std::vector<float>>("phi_bsc");
        fields.mBackscatteringEnergy = entry.GetPtr<std::vector<float>>("e_bsc");
    }

    auto AddExitPoint(std::vector<std::string>& particleNames, std::vector<float>& thetas, std::vector<float>& phis,
                      std::vector<float>& energies, float& totalEnergy, const std::string& particleName,
                      const G4ThreeVector& direction, float energy) -> void {
        totalEnergy += energy;
        particleNames.push_back(particleName);
        const auto z = std::clamp(direction.z(), -1.0, 1.0);
        thetas.push_back(static_cast<float>(std::acos(z)));
        phis.push_back(static_cast<float>(std::atan2(direction.y(), direction.x())));
        energies.push_back(energy);
    }

    RunAction& mRunAction;
    int mLayerCount{0};
    int mSourceCount{0};
    int mCurrentSourceIndex{0};
    int mRunId{-1};
    std::vector<std::weak_ptr<ROOT::RNTupleFillContext>> mFillContexts;
    std::vector<std::unique_ptr<ROOT::REntry>> mEntries;
    std::vector<ROOT::RNTupleFillStatus> mFillStatuses;
    std::vector<SourceFields> mFields;
    float mTotalPenetrationEnergy{0.0F};
    std::vector<std::string> mPenetrationParticles;
    std::vector<float> mPenetrationTheta;
    std::vector<float> mPenetrationPhi;
    std::vector<float> mPenetrationEnergy;
    float mTotalDepositionEnergy{0.0F};
    std::vector<float> mLayerEnergyDeposit;
    std::vector<std::vector<std::string>> mDepositionParticles;
    std::vector<std::vector<float>> mDepositionX;
    std::vector<std::vector<float>> mDepositionY;
    std::vector<std::vector<float>> mDepositionZ;
    std::vector<std::vector<float>> mDepositionWeight;
    std::vector<std::vector<std::string>> mDepositionProcess;
    float mTotalBackscatteringEnergy{0.0F};
    std::vector<std::string> mBackscatteringParticles;
    std::vector<float> mBackscatteringTheta;
    std::vector<float> mBackscatteringPhi;
    std::vector<float> mBackscatteringEnergy;
};

class SteppingAction : public G4UserSteppingAction {
public:
    explicit SteppingAction(EventAction& eventAction) :
        mEventAction{eventAction} {}
    ~SteppingAction() override = default;

    auto UserSteppingAction(const G4Step* step) -> void override {
        if (mMaterialVolumes.empty()) {
            const DetectorConstruction* detectorConstruction = static_cast<const DetectorConstruction*>(
                G4RunManager::GetRunManager()->GetUserDetectorConstruction());
            mMaterialVolumes = detectorConstruction->GetMaterialVolumes();
        }
        const auto logicalVolume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
        const auto layerIt = mMaterialVolumes.find(logicalVolume);
        if (layerIt == mMaterialVolumes.end()) {
            return;
        }
        const auto layerIndex = layerIt->second;
        const auto energyDeposit = step->GetTotalEnergyDeposit();
        if (energyDeposit <= 0.0) {
            return;
        }
        const auto process = step->GetPostStepPoint()->GetProcessDefinedStep();
        const auto processName = process != nullptr ? process->GetProcessName() : "<null>";
        mEventAction.AddDeposition(layerIndex, step->GetTrack()->GetDefinition()->GetParticleName(),
                                   step->GetPostStepPoint()->GetPosition(), static_cast<float>(energyDeposit), processName);
    }

private:
    EventAction& mEventAction;
    std::unordered_map<const G4LogicalVolume*, int> mMaterialVolumes;
};

class TrackingAction : public G4UserTrackingAction {
public:
    explicit TrackingAction(EventAction& eventAction, bool includeNeutrinos) :
        mEventAction{eventAction},
        mIncludeNeutrinos{includeNeutrinos} {}
    ~TrackingAction() override = default;

    auto PostUserTrackingAction(const G4Track* track) -> void override {
        const auto step = track->GetStep();
        if (step->GetPostStepPoint()->GetStepStatus() != fWorldBoundary) {
            return;
        }
        const auto definition = track->GetDefinition();
        if (not mIncludeNeutrinos and IsNeutrino(definition)) {
            return;
        }
        const auto direction = track->GetMomentumDirection();
        const auto energy = static_cast<float>(track->GetKineticEnergy());
        const auto particleName = definition->GetParticleName();
        if (direction.z() >= 0.0) {
            mEventAction.AddPenetration(particleName, direction, energy);
        } else {
            mEventAction.AddBackscattering(particleName, direction, energy);
        }
    }

private:
    EventAction& mEventAction;
    bool mIncludeNeutrinos{false};
};

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
    explicit PrimaryGeneratorAction(const Config& config, EventAction& eventAction) :
        mConfig{config},
        mEventAction{eventAction} {
        auto cumulative{0.0};
        mParticleDefinitions.reserve(config.mSources.size());
        for (const auto& source : config.mSources) {
            G4ParticleDefinition* particle = G4ParticleTable::GetParticleTable()->FindParticle(source.mParticleName);
            if (particle == nullptr) {
                throw std::invalid_argument("particle '" + source.mParticleName + "' not found in the particle table");
            }
            mParticleDefinitions.push_back(particle);
            // The intensities were already normalized into probabilities (mWeight); the cumulative
            // probabilities are used to draw the source type of each primary from the multinomial
            // distribution (one trial per event because every event carries a single primary particle)
            // with Geant4's per-thread random engine.
            cumulative += source.mWeight;
            mCumulativeProbabilities.push_back(cumulative);
        }
        // Guard against floating-point rounding in the last interval.
        mCumulativeProbabilities.back() = 1.0;
        mParticleGun = std::make_unique<G4ParticleGun>(1);
        mParticleGun->SetParticleMomentumDirection(G4ThreeVector{0.0, 0.0, 1.0});
    }
    ~PrimaryGeneratorAction() override = default;

    auto GeneratePrimaries(G4Event* event) -> void override {
        const auto sourceIndex = DrawSourceIndex();
        const auto& source = mConfig.mSources[sourceIndex];
        mEventAction.SetSourceIndex(sourceIndex);
        mParticleGun->SetParticleDefinition(mParticleDefinitions[sourceIndex]);
        mParticleGun->SetParticleEnergy(source.mEnergy);
        mParticleGun->SetParticlePosition(mSourcePosition);
        mParticleGun->GeneratePrimaryVertex(event);
    }

private:
    // Draw the source type for the current primary particle: a single trial of the multinomial
    // distribution over the normalized source intensities, sampled with Geant4's random engine
    // (each worker thread owns its own engine instance, so MT runs stay independent).
    auto DrawSourceIndex() const -> int {
        const auto uniformSample{G4UniformRand()};
        const auto it = std::upper_bound(mCumulativeProbabilities.begin(), mCumulativeProbabilities.end(),
                                         uniformSample);
        return it == mCumulativeProbabilities.end() ? static_cast<int>(mCumulativeProbabilities.size()) - 1 : static_cast<int>(it - mCumulativeProbabilities.begin());
    }

    const Config& mConfig;
    EventAction& mEventAction;
    std::unique_ptr<G4ParticleGun> mParticleGun;
    std::vector<G4ParticleDefinition*> mParticleDefinitions;
    std::vector<double> mCumulativeProbabilities;
    G4ThreeVector mSourcePosition{0.0, 0.0, 0.0};
};

class ActionInitialization : public G4VUserActionInitialization {
public:
    explicit ActionInitialization(const Config& config) :
        mConfig{config} {}
    ~ActionInitialization() override = default;

    auto BuildForMaster() const -> void override {
        SetUserAction(new RunAction{mConfig});
    }

    auto Build() const -> void override {
        RunAction* runAction = new RunAction{mConfig};
        SetUserAction(runAction);
        EventAction* eventAction = new EventAction{*runAction, static_cast<int>(mConfig.mLayers.size()),
                                                   static_cast<int>(mConfig.mSources.size())};
        SetUserAction(eventAction);
        SetUserAction(new PrimaryGeneratorAction{mConfig, *eventAction});
        SetUserAction(new SteppingAction{*eventAction});
        SetUserAction(new TrackingAction{*eventAction, mConfig.mIncludeNeutrinos});
    }

private:
    const Config& mConfig;
};

constexpr std::array defaultVisCommands = {
    "/vis/open OGL",
    "/vis/drawVolume",
    "/vis/filtering/trajectories/create/particleFilter neutrinoFilter",
    "/vis/filtering/trajectories/particleFilter/neutrinoFilter/add nu_e nu_mu nu_tau anti_nu_e anti_nu_mu anti_nu_tau",
    "/vis/filtering/trajectories/particleFilter/neutrinoFilter/invert true",
    "/vis/scene/add/trajectories smooth",
    "/vis/scene/endOfEventAction accumulate",
    "/vis/viewer/set/autoRefresh true",
    "/vis/viewer/set/viewpointThetaPhi 90 0",
};

} // namespace PPS

auto main(int argc, char** argv) -> int try {
    const auto config{PPS::ParseCommandLine(argc, argv)};
    if (config.mHelp) {
        PPS::PrintUsage(argv[0]);
        return EXIT_SUCCESS;
    }

    PPS::ValidateMaterials(config);

    const auto runManager =
        std::unique_ptr<G4RunManager>{G4RunManagerFactory::CreateRunManager(
            config.mThreads > 1 ? G4RunManagerType::MT : G4RunManagerType::Serial, config.mThreads)};
    runManager->SetVerboseLevel(config.mVerbose);

    runManager->SetUserInitialization(new PPS::DetectorConstruction{config});

    G4PhysListFactory physicsListFactory(config.mVerbose);
    if (not physicsListFactory.IsReferencePhysList(config.mPhysicsListName)) {
        G4cerr << "error: physics list '" << config.mPhysicsListName << "' is not available; available lists:";
        for (const auto& name : physicsListFactory.AvailablePhysLists()) {
            G4cerr << ' ' << name;
        }
        G4cerr << G4endl;
        return EXIT_FAILURE;
    }
    const auto physicsList = physicsListFactory.GetReferencePhysList(config.mPhysicsListName);
    physicsList->SetVerboseLevel(config.mVerbose);
    const auto emParameters = G4EmParameters::Instance();
    emParameters->SetMscMuHadStepLimitType(emParameters->MscStepLimitType());
    emParameters->SetVerbose(config.mVerbose);
    G4HadronicParameters::Instance()->SetVerboseLevel(config.mVerbose);
    G4NuclearLevelData::GetInstance()->GetParameters()->SetVerbose(config.mVerbose);
    runManager->SetUserInitialization(physicsList);

    PPS::ValidateSources(config);

    PPS::OutputWriter::Instance().Initialize(config.mOutputFileName, config.mOverwrite);
    runManager->SetUserInitialization(new PPS::ActionInitialization{config});
    const auto visManager = std::unique_ptr<G4VisExecutive>{new G4VisExecutive{"quiet"}};
    visManager->Initialize();
    runManager->Initialize();
    const auto uiManager = G4UImanager::GetUIpointer();
    uiManager->SetVerboseLevel(config.mVerbose);
    const auto printProgress =
        config.mEventCount > 0 ? std::max(1, config.mEventCount / 10) : config.mPrintProgress;
    uiManager->ApplyCommand("/run/printProgress " + std::to_string(printProgress));
    if (config.mUI) {
        std::unique_ptr<G4UIExecutive> uiExecutive{
            new G4UIExecutive{argc, argv}
        };
        for (const auto& command : PPS::defaultVisCommands) {
            uiManager->ApplyCommand(command);
        }
        if (config.mEventCount > 0) {
            uiManager->ApplyCommand("/run/beamOn " + std::to_string(config.mEventCount));
        }
        uiExecutive->SessionStart();
    } else if (not config.mMacroFileName.empty()) {
        uiManager->ApplyCommand("/control/execute " + config.mMacroFileName);
    } else {
        uiManager->ApplyCommand("/run/beamOn " + std::to_string(config.mEventCount));
    }
    PPS::OutputWriter::Instance().Finalize();

    return EXIT_SUCCESS;
} catch (const std::exception& error) {
    G4cerr << "error: " << error.what() << G4endl;
    PPS::OutputWriter::Instance().Finalize();
    std::quick_exit(EXIT_FAILURE);
}
