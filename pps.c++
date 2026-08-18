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
#include "ROOT/REntry.hxx"
#include "ROOT/RNTupleFillContext.hxx"
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

constexpr auto poissonUpperLimit{2.3};
constexpr auto ratioTolerance{1e-6};

struct Layer {
    std::string mMaterialName;
    double mThickness{0.0};
};

struct Config {
    std::vector<Layer> mLayers;
    double mEnergy{0.0};
    std::string mParticleName;
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
    if (suffix.empty()) {
        return value;
    }
    const auto unit = units.find(suffix);
    if (unit == units.end()) {
        throw std::invalid_argument("unknown unit '" + suffix + "' in '" + text + "'");
    }
    return value * unit->second;
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

auto FormatEnergy(double energy) -> std::string {
    auto stream = std::ostringstream{};
    stream << G4BestUnit(energy, "Energy");
    return stream.str();
}

auto PrintUsage(const char* programName) -> void {
    G4cout
        << "usage: " << programName << " [options] [macroFile]\n"
        << "simulate one particle per event passing through a material slab and record per-event\n"
        << "energy depositions inside the slab and particles penetrating or backscattering into an\n"
        << "RNTuple 'run{runId}' inside a ROOT file\n"
        << "\n"
        << "required options:\n"
        << "  -m, --material <spec>    material layers as semicolon-separated 'name:thickness' entries, built in\n"
        << "                           order from the source, e.g. \"G4_WATER:10 cm;G4_Cu:5 cm\"\n"
        << "  -e, --energy <value>     primary kinetic energy; bare value in MeV or with a unit, e.g. \"100 MeV\"\n"
        << "  -p, --particle <name>    primary particle name, e.g. proton, neutron, gamma (Geant4 particle table)\n"
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
        } else if (argument == "-e" or argument == "--energy") {
            config.mEnergy = ParseEnergy(nextValue(i, argument));
        } else if (argument == "-p" or argument == "--particle") {
            config.mParticleName = nextValue(i, argument);
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
        if (config.mEnergy <= 0.0) {
            throw std::invalid_argument("required option --energy is missing or not positive");
        }
        if (config.mParticleName.empty()) {
            throw std::invalid_argument("required option --particle is missing");
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

auto ValidateParticle(const Config& config) -> void {
    if (G4ParticleTable::GetParticleTable()->FindParticle(config.mParticleName) == nullptr) {
        throw std::invalid_argument("particle '" + config.mParticleName + "' not found in the particle table");
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

    auto BeginRun(int runId, int layerCount) -> void {
        try {
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
            auto options = ROOT::RNTupleWriteOptions{};
            options.SetCompression(ROOT::RCompressionSetting::EDefaults::kUseGeneralPurpose);
            options.SetUseImplicitMT(ROOT::RNTupleWriteOptions::EImplicitMT::kOff);
            mWriter = ROOT::RNTupleParallelWriter::Append(std::move(model), "run" + std::to_string(runId), *mFile, options);
        } catch (const std::exception& error) {
            G4cerr << "error: failed to create RNTuple 'run" << runId << "': " << error.what() << G4endl;
            std::quick_exit(EXIT_FAILURE);
        }
    }

    auto CreateFillContext() -> std::shared_ptr<ROOT::RNTupleFillContext> {
        return mWriter->CreateFillContext();
    }

    auto EndRun() -> void {
        if (mWriter != nullptr) {
            mWriter->CommitDataset();
            mWriter.reset();
        }
    }

    auto Finalize() -> void {
        EndRun();
        if (mFile != nullptr) {
            mFile->Close();
            mFile.reset();
        }
        mInitialized = false;
    }

private:
    OutputWriter() = default;

    std::unique_ptr<TFile> mFile;
    std::unique_ptr<ROOT::RNTupleParallelWriter> mWriter;
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
    explicit SimulationRun(int layerCount) :
        mLayerDepEnergySums(layerCount, 0.0),
        mLayerDepEnergySumsSq(layerCount, 0.0) {}

    auto Merge(const G4Run* other) -> void override {
        G4Run::Merge(other);
        const auto* otherSimulationRun = static_cast<const SimulationRun*>(other);
        mPenEnergySum += otherSimulationRun->mPenEnergySum;
        mPenEnergySumSq += otherSimulationRun->mPenEnergySumSq;
        mDepEnergySum += otherSimulationRun->mDepEnergySum;
        mDepEnergySumSq += otherSimulationRun->mDepEnergySumSq;
        for (auto layerIndex{0}; layerIndex < static_cast<int>(mLayerDepEnergySums.size()); ++layerIndex) {
            mLayerDepEnergySums[layerIndex] += otherSimulationRun->mLayerDepEnergySums[layerIndex];
            mLayerDepEnergySumsSq[layerIndex] += otherSimulationRun->mLayerDepEnergySumsSq[layerIndex];
        }
        mBscEnergySum += otherSimulationRun->mBscEnergySum;
        mBscEnergySumSq += otherSimulationRun->mBscEnergySumSq;
        AddToMap(mPenParticleCount, otherSimulationRun->mPenParticleCount);
        AddToMap(mPenParticleEnergySum, otherSimulationRun->mPenParticleEnergySum);
        AddToMap(mPenParticleEnergySumSq, otherSimulationRun->mPenParticleEnergySumSq);
        AddToMap(mBscParticleCount, otherSimulationRun->mBscParticleCount);
        AddToMap(mBscParticleEnergySum, otherSimulationRun->mBscParticleEnergySum);
        AddToMap(mBscParticleEnergySumSq, otherSimulationRun->mBscParticleEnergySumSq);
    }

    auto AddEventResult(double penetrationEnergy, const std::vector<float>& layerDepositionEnergies,
                        double backscatteringEnergy) -> void {
        mPenEnergySum += penetrationEnergy;
        mPenEnergySumSq += penetrationEnergy * penetrationEnergy;
        auto totalDepositionEnergy{0.0};
        for (auto layerIndex{0}; layerIndex < static_cast<int>(layerDepositionEnergies.size()); ++layerIndex) {
            const auto layerDepositionEnergy = layerDepositionEnergies[layerIndex];
            mLayerDepEnergySums[layerIndex] += layerDepositionEnergy;
            mLayerDepEnergySumsSq[layerIndex] += layerDepositionEnergy * layerDepositionEnergy;
            totalDepositionEnergy += layerDepositionEnergy;
        }
        mDepEnergySum += totalDepositionEnergy;
        mDepEnergySumSq += totalDepositionEnergy * totalDepositionEnergy;
        mBscEnergySum += backscatteringEnergy;
        mBscEnergySumSq += backscatteringEnergy * backscatteringEnergy;
    }

    auto AddPenetrationEnergy(const std::string& particleName, double energy) -> void {
        mPenParticleCount[particleName] += 1.0;
        mPenParticleEnergySum[particleName] += energy;
        mPenParticleEnergySumSq[particleName] += energy * energy;
    }

    auto AddBackscatteringEnergy(const std::string& particleName, double energy) -> void {
        mBscParticleCount[particleName] += 1.0;
        mBscParticleEnergySum[particleName] += energy;
        mBscParticleEnergySumSq[particleName] += energy * energy;
    }

    auto GetPenEnergySum() const -> double { return mPenEnergySum; }
    auto GetPenEnergySumSq() const -> double { return mPenEnergySumSq; }
    auto GetDepEnergySum() const -> double { return mDepEnergySum; }
    auto GetDepEnergySumSq() const -> double { return mDepEnergySumSq; }
    auto GetLayerDepEnergySum(int layerIndex) const -> double { return mLayerDepEnergySums[layerIndex]; }
    auto GetLayerDepEnergySumSq(int layerIndex) const -> double { return mLayerDepEnergySumsSq[layerIndex]; }
    auto GetLayerCount() const -> int { return static_cast<int>(mLayerDepEnergySums.size()); }
    auto GetBscEnergySum() const -> double { return mBscEnergySum; }
    auto GetBscEnergySumSq() const -> double { return mBscEnergySumSq; }
    auto GetPenParticleCounts() const -> const std::map<std::string, double>& { return mPenParticleCount; }
    auto GetPenParticleEnergySums() const -> const std::map<std::string, double>& { return mPenParticleEnergySum; }
    auto GetPenParticleEnergySumsSq() const -> const std::map<std::string, double>& { return mPenParticleEnergySumSq; }
    auto GetBscParticleCounts() const -> const std::map<std::string, double>& { return mBscParticleCount; }
    auto GetBscParticleEnergySums() const -> const std::map<std::string, double>& { return mBscParticleEnergySum; }
    auto GetBscParticleEnergySumsSq() const -> const std::map<std::string, double>& { return mBscParticleEnergySumSq; }

private:
    auto AddToMap(std::map<std::string, double>& target, const std::map<std::string, double>& source) -> void {
        for (const auto& [key, value] : source) {
            target[key] += value;
        }
    }

    double mPenEnergySum{0.0};
    double mPenEnergySumSq{0.0};
    double mDepEnergySum{0.0};
    double mDepEnergySumSq{0.0};
    std::vector<double> mLayerDepEnergySums;
    std::vector<double> mLayerDepEnergySumsSq;
    double mBscEnergySum{0.0};
    double mBscEnergySumSq{0.0};
    std::map<std::string, double> mPenParticleCount;
    std::map<std::string, double> mPenParticleEnergySum;
    std::map<std::string, double> mPenParticleEnergySumSq;
    std::map<std::string, double> mBscParticleCount;
    std::map<std::string, double> mBscParticleEnergySum;
    std::map<std::string, double> mBscParticleEnergySumSq;
};

class RunAction : public G4UserRunAction {
public:
    explicit RunAction(const Config& config) :
        mConfig{config} {}
    ~RunAction() override = default;

    auto GenerateRun() -> G4Run* override {
        return new SimulationRun{static_cast<int>(mConfig.mLayers.size())};
    }

    auto BeginOfRunAction(const G4Run* run) -> void override {
        mCurrentRun = const_cast<SimulationRun*>(static_cast<const SimulationRun*>(run));
        if (IsMaster()) {
            OutputWriter::Instance().BeginRun(run->GetRunID(), static_cast<int>(mConfig.mLayers.size()));
        }
        if (not IsMaster() or not G4Threading::IsMultithreadedApplication()) {
            mFillContext = OutputWriter::Instance().CreateFillContext();
            mEntry = mFillContext->CreateEntry();
        }
    }

    auto EndOfRunAction(const G4Run* run) -> void override {
        ReleaseFillState();
        if (IsMaster()) {
            OutputWriter::Instance().EndRun();
            PrintStatistics(static_cast<const SimulationRun&>(*run));
            G4cout << "run " << run->GetRunID() << " finished with " << run->GetNumberOfEvent() << " events" << G4endl;
        }
    }

    auto AddEventResult(double penetrationEnergy, const std::vector<float>& layerDepositionEnergies,
                        double backscatteringEnergy) -> void {
        mCurrentRun->AddEventResult(penetrationEnergy, layerDepositionEnergies, backscatteringEnergy);
    }

    auto AddPenetrationEnergy(const std::string& particleName, double energy) -> void {
        mCurrentRun->AddPenetrationEnergy(particleName, energy);
    }

    auto AddBackscatteringEnergy(const std::string& particleName, double energy) -> void {
        mCurrentRun->AddBackscatteringEnergy(particleName, energy);
    }

    auto GetFillContext() const -> std::shared_ptr<ROOT::RNTupleFillContext> {
        return mFillContext;
    }

private:
    auto ReleaseFillState() -> void {
        mEntry.reset();
        mFillContext.reset();
    }

    auto PrintStatistics(const SimulationRun& run) -> void {
        const auto eventCount = run.GetNumberOfEvent();
        if (eventCount < 1) {
            return;
        }
        G4cout << '\n';
        G4cout << "===============================================================================\n";
        G4cout << " energy ratio (incident " << G4BestUnit(mConfig.mEnergy, "Energy") << ' ' << mConfig.mParticleName
               << ", " << eventCount << " events):\n";
        PrintEnergyRatio("penetration ratio (pen)", run.GetPenEnergySum(), run.GetPenEnergySumSq(), eventCount);
        PrintEnergyRatio("deposition ratio (dep)", run.GetDepEnergySum(), run.GetDepEnergySumSq(), eventCount);
        for (auto layerIndex{0}; layerIndex < run.GetLayerCount(); ++layerIndex) {
            const auto label = "  in layer " + std::to_string(layerIndex) + " (" + mConfig.mLayers[layerIndex].mMaterialName + ')';
            PrintEnergyRatio(label, run.GetLayerDepEnergySum(layerIndex), run.GetLayerDepEnergySumSq(layerIndex), eventCount);
        }
        PrintEnergyRatio("back-scattering ratio (bsc)", run.GetBscEnergySum(), run.GetBscEnergySumSq(), eventCount);
        G4cout << "-------------------------------------------------------------------------------\n";
        PrintParticleStatistics("penetration particles", run.GetPenParticleCounts(), run.GetPenParticleEnergySums(), run.GetPenParticleEnergySumsSq());
        PrintParticleStatistics("back-scattering particles", run.GetBscParticleCounts(), run.GetBscParticleEnergySums(), run.GetBscParticleEnergySumsSq());
        G4cout << "===============================================================================\n"
               << G4endl;
    }

    auto PrintParticleStatistics(const std::string& title, const std::map<std::string, double>& particleCounts,
                                 const std::map<std::string, double>& energySums,
                                 const std::map<std::string, double>& energySumsSq) -> void {
        if (particleCounts.empty()) {
            return;
        }
        G4cout << ' ' << title << ":\n";
        G4cout << "   "
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
            G4cout << "   "
                   << std::setw(14) << particleName
                   << std::setw(20) << FormatEnergy(meanEnergy)
                   << std::setw(20) << FormatEnergy(rmsEnergy)
                   << std::setw(8) << static_cast<long long>(count) << '\n';
        }
    }

    auto PrintEnergyRatio(const std::string& label, double sumEnergy, double sumEnergySq, G4int eventCount) -> void {
        const auto meanEnergy = sumEnergy / eventCount;
        const auto ratio = meanEnergy / mConfig.mEnergy;
        const auto variance = (sumEnergySq - sumEnergy * sumEnergy / eventCount) / (eventCount - 1);
        const auto ratioError = std::sqrt(variance / eventCount) / mConfig.mEnergy;
        G4cout << "   " << std::left << std::setw(45) << label << '(' << 100.0 * ratio << " +/- "
               << 100.0 * ratioError << ") %" << G4endl;
    }

    const Config& mConfig;
    SimulationRun* mCurrentRun{nullptr};
    std::shared_ptr<ROOT::RNTupleFillContext> mFillContext;
    std::unique_ptr<ROOT::REntry> mEntry;
};

class EventAction : public G4UserEventAction {
public:
    explicit EventAction(RunAction* runAction, int layerCount) :
        mRunAction{runAction},
        mLayerCount{layerCount},
        mLayerEdep(layerCount, 0.0F),
        mParticleDep(layerCount),
        mXDep(layerCount),
        mYDep(layerCount),
        mZDep(layerCount),
        mWDep(layerCount),
        mProcDep(layerCount) {}
    ~EventAction() override = default;

    auto BeginOfEventAction(const G4Event*) -> void override {
        const auto runId = G4RunManager::GetRunManager()->GetCurrentRun()->GetRunID();
        if (runId != mRunId) {
            mRunId = runId;
            mFillContext = mRunAction->GetFillContext().get();
            mEntry = mFillContext->CreateEntry();
            mEventIdField = mEntry->GetPtr<int>("event_id");
            mTotalEPenField = mEntry->GetPtr<float>("total_e_pen");
            mParticlePenField = mEntry->GetPtr<std::vector<std::string>>("particle_pen");
            mThetaPenField = mEntry->GetPtr<std::vector<float>>("theta_pen");
            mPhiPenField = mEntry->GetPtr<std::vector<float>>("phi_pen");
            mEPenField = mEntry->GetPtr<std::vector<float>>("e_pen");
            mTotalEdepField = mEntry->GetPtr<float>("total_e_dep");
            mLayerEdepFields.clear();
            mParticleDepFields.clear();
            mXDepFields.clear();
            mYDepFields.clear();
            mZDepFields.clear();
            mWDepFields.clear();
            mProcDepFields.clear();
            for (auto layerIndex{0}; layerIndex < mLayerCount; ++layerIndex) {
                const auto suffix = std::to_string(layerIndex);
                mLayerEdepFields.push_back(mEntry->GetPtr<float>("e_dep_" + suffix));
                mParticleDepFields.push_back(mEntry->GetPtr<std::vector<std::string>>("particle_dep_" + suffix));
                mXDepFields.push_back(mEntry->GetPtr<std::vector<float>>("x_dep_" + suffix));
                mYDepFields.push_back(mEntry->GetPtr<std::vector<float>>("y_dep_" + suffix));
                mZDepFields.push_back(mEntry->GetPtr<std::vector<float>>("z_dep_" + suffix));
                mWDepFields.push_back(mEntry->GetPtr<std::vector<float>>("w_dep_" + suffix));
                mProcDepFields.push_back(mEntry->GetPtr<std::vector<std::string>>("proc_dep_" + suffix));
            }
            mTotalEBscField = mEntry->GetPtr<float>("total_e_bsc");
            mParticleBscField = mEntry->GetPtr<std::vector<std::string>>("particle_bsc");
            mThetaBscField = mEntry->GetPtr<std::vector<float>>("theta_bsc");
            mPhiBscField = mEntry->GetPtr<std::vector<float>>("phi_bsc");
            mEBscField = mEntry->GetPtr<std::vector<float>>("e_bsc");
        }
        mTotalEPen = 0.0F;
        mParticlePen.clear();
        mThetaPen.clear();
        mPhiPen.clear();
        mEPen.clear();
        mTotalEdep = 0.0F;
        for (auto& layerEdep : mLayerEdep) {
            layerEdep = 0.0F;
        }
        for (auto& particleDep : mParticleDep) {
            particleDep.clear();
        }
        for (auto& xDep : mXDep) {
            xDep.clear();
        }
        for (auto& yDep : mYDep) {
            yDep.clear();
        }
        for (auto& zDep : mZDep) {
            zDep.clear();
        }
        for (auto& wDep : mWDep) {
            wDep.clear();
        }
        for (auto& procDep : mProcDep) {
            procDep.clear();
        }
        mTotalEBsc = 0.0F;
        mParticleBsc.clear();
        mThetaBsc.clear();
        mPhiBsc.clear();
        mEBsc.clear();
    }

    auto EndOfEventAction(const G4Event* event) -> void override {
        *mEventIdField = event->GetEventID();
        *mTotalEPenField = mTotalEPen;
        *mParticlePenField = mParticlePen;
        *mThetaPenField = mThetaPen;
        *mPhiPenField = mPhiPen;
        *mEPenField = mEPen;
        *mTotalEdepField = mTotalEdep;
        for (auto layerIndex{0}; layerIndex < mLayerCount; ++layerIndex) {
            *mLayerEdepFields[layerIndex] = mLayerEdep[layerIndex];
            *mParticleDepFields[layerIndex] = mParticleDep[layerIndex];
            *mXDepFields[layerIndex] = mXDep[layerIndex];
            *mYDepFields[layerIndex] = mYDep[layerIndex];
            *mZDepFields[layerIndex] = mZDep[layerIndex];
            *mWDepFields[layerIndex] = mWDep[layerIndex];
            *mProcDepFields[layerIndex] = mProcDep[layerIndex];
        }
        *mTotalEBscField = mTotalEBsc;
        *mParticleBscField = mParticleBsc;
        *mThetaBscField = mThetaBsc;
        *mPhiBscField = mPhiBsc;
        *mEBscField = mEBsc;
        mFillContext->Fill(*mEntry);
        mRunAction->AddEventResult(mTotalEPen, mLayerEdep, mTotalEBsc);
    }

    auto AddPenetration(const std::string& particleName, const G4ThreeVector& direction, float energy) -> void {
        AddExitPoint(mParticlePen, mThetaPen, mPhiPen, mEPen, mTotalEPen, particleName, direction, energy);
        mRunAction->AddPenetrationEnergy(particleName, energy);
    }

    auto AddDeposition(int layerIndex, const std::string& particleName, const G4ThreeVector& position, float edep,
                       const std::string& processName) -> void {
        mTotalEdep += edep;
        mLayerEdep[layerIndex] += edep;
        mParticleDep[layerIndex].push_back(particleName);
        mXDep[layerIndex].push_back(static_cast<float>(position.x()));
        mYDep[layerIndex].push_back(static_cast<float>(position.y()));
        mZDep[layerIndex].push_back(static_cast<float>(position.z()));
        mWDep[layerIndex].push_back(edep);
        mProcDep[layerIndex].push_back(processName);
    }

    auto AddBackscattering(const std::string& particleName, const G4ThreeVector& direction, float energy) -> void {
        AddExitPoint(mParticleBsc, mThetaBsc, mPhiBsc, mEBsc, mTotalEBsc, particleName, direction, energy);
        mRunAction->AddBackscatteringEnergy(particleName, energy);
    }

private:
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

    RunAction* mRunAction;
    int mLayerCount{0};
    ROOT::RNTupleFillContext* mFillContext{nullptr};
    std::unique_ptr<ROOT::REntry> mEntry;
    int mRunId{-1};
    std::shared_ptr<int> mEventIdField;
    std::shared_ptr<float> mTotalEPenField;
    std::shared_ptr<std::vector<std::string>> mParticlePenField;
    std::shared_ptr<std::vector<float>> mThetaPenField;
    std::shared_ptr<std::vector<float>> mPhiPenField;
    std::shared_ptr<std::vector<float>> mEPenField;
    std::shared_ptr<float> mTotalEdepField;
    std::vector<std::shared_ptr<float>> mLayerEdepFields;
    std::vector<std::shared_ptr<std::vector<std::string>>> mParticleDepFields;
    std::vector<std::shared_ptr<std::vector<float>>> mXDepFields;
    std::vector<std::shared_ptr<std::vector<float>>> mYDepFields;
    std::vector<std::shared_ptr<std::vector<float>>> mZDepFields;
    std::vector<std::shared_ptr<std::vector<float>>> mWDepFields;
    std::vector<std::shared_ptr<std::vector<std::string>>> mProcDepFields;
    std::shared_ptr<float> mTotalEBscField;
    std::shared_ptr<std::vector<std::string>> mParticleBscField;
    std::shared_ptr<std::vector<float>> mThetaBscField;
    std::shared_ptr<std::vector<float>> mPhiBscField;
    std::shared_ptr<std::vector<float>> mEBscField;
    float mTotalEPen{0.0F};
    std::vector<std::string> mParticlePen;
    std::vector<float> mThetaPen;
    std::vector<float> mPhiPen;
    std::vector<float> mEPen;
    float mTotalEdep{0.0F};
    std::vector<float> mLayerEdep;
    std::vector<std::vector<std::string>> mParticleDep;
    std::vector<std::vector<float>> mXDep;
    std::vector<std::vector<float>> mYDep;
    std::vector<std::vector<float>> mZDep;
    std::vector<std::vector<float>> mWDep;
    std::vector<std::vector<std::string>> mProcDep;
    float mTotalEBsc{0.0F};
    std::vector<std::string> mParticleBsc;
    std::vector<float> mThetaBsc;
    std::vector<float> mPhiBsc;
    std::vector<float> mEBsc;
};

class SteppingAction : public G4UserSteppingAction {
public:
    explicit SteppingAction(EventAction* eventAction) :
        mEventAction{eventAction} {}
    ~SteppingAction() override = default;

    auto UserSteppingAction(const G4Step* step) -> void override {
        if (mMaterialVolumes.empty()) {
            const auto detectorConstruction = static_cast<const DetectorConstruction*>(
                G4RunManager::GetRunManager()->GetUserDetectorConstruction());
            mMaterialVolumes = detectorConstruction->GetMaterialVolumes();
        }
        const auto logicalVolume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
        const auto layerIt = mMaterialVolumes.find(logicalVolume);
        if (layerIt == mMaterialVolumes.end()) {
            return;
        }
        const auto layerIndex = layerIt->second;
        const auto edep = step->GetTotalEnergyDeposit();
        if (edep <= 0.0) {
            return;
        }
        const auto process = step->GetPostStepPoint()->GetProcessDefinedStep();
        const auto processName = process != nullptr ? process->GetProcessName() : "<null>";
        mEventAction->AddDeposition(layerIndex, step->GetTrack()->GetDefinition()->GetParticleName(),
                                    step->GetPostStepPoint()->GetPosition(), static_cast<float>(edep), processName);
    }

private:
    EventAction* mEventAction;
    std::unordered_map<const G4LogicalVolume*, int> mMaterialVolumes;
};

class TrackingAction : public G4UserTrackingAction {
public:
    explicit TrackingAction(EventAction* eventAction, bool includeNeutrinos) :
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
            mEventAction->AddPenetration(particleName, direction, energy);
        } else {
            mEventAction->AddBackscattering(particleName, direction, energy);
        }
    }

private:
    EventAction* mEventAction;
    bool mIncludeNeutrinos{false};
};

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
    explicit PrimaryGeneratorAction(const Config& config) :
        mConfig{config} {
        const auto particle = G4ParticleTable::GetParticleTable()->FindParticle(config.mParticleName);
        if (particle == nullptr) {
            throw std::invalid_argument("particle '" + config.mParticleName + "' not found in the particle table");
        }
        mParticleGun = std::make_unique<G4ParticleGun>(1);
        mParticleGun->SetParticleDefinition(particle);
        mParticleGun->SetParticleMomentumDirection(G4ThreeVector{0.0, 0.0, 1.0});
        mParticleGun->SetParticleEnergy(config.mEnergy);
    }
    ~PrimaryGeneratorAction() override = default;

    auto GeneratePrimaries(G4Event* event) -> void override {
        mParticleGun->SetParticlePosition(mSourcePosition);
        mParticleGun->GeneratePrimaryVertex(event);
    }

private:
    const Config& mConfig;
    std::unique_ptr<G4ParticleGun> mParticleGun;
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
        const auto runAction = new RunAction{mConfig};
        SetUserAction(runAction);
        const auto eventAction = new EventAction{runAction, static_cast<int>(mConfig.mLayers.size())};
        SetUserAction(eventAction);
        SetUserAction(new PrimaryGeneratorAction{mConfig});
        SetUserAction(new SteppingAction{eventAction});
        SetUserAction(new TrackingAction{eventAction, mConfig.mIncludeNeutrinos});
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

    PPS::ValidateParticle(config);

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
        std::unique_ptr<G4UIExecutive> ui(new G4UIExecutive{argc, argv});
        for (const auto& command : PPS::defaultVisCommands) {
            uiManager->ApplyCommand(command);
        }
        if (config.mEventCount > 0) {
            uiManager->ApplyCommand("/run/beamOn " + std::to_string(config.mEventCount));
        }
        ui->SessionStart();
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
