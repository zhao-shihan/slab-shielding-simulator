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

enum struct MaterialCompositionType {
    nist,
    atoms,
    massFractions,
    materials
};

struct MaterialComponent {
    std::string mName;
    double mQuantity{1.0};
    bool mQuantitySpecified{false};
};

struct LayerMaterial {
    std::string mName;
    MaterialCompositionType mKind{MaterialCompositionType::nist};
    std::vector<MaterialComponent> mComponents;
    double mDensity{0.0};
};

struct Layer {
    LayerMaterial mMaterial;
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
    bool mSave{false};
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

auto ParseDensity(const std::string& text) -> double {
    auto stream = std::istringstream{text};
    auto value{0.0};
    auto suffix = std::string{};
    stream >> value;
    if (stream.fail()) {
        throw std::invalid_argument("cannot parse numeric value from '" + text + "'");
    }
    stream >> suffix;
    // A bare number is interpreted in g/cm3, the conventional density unit.
    auto unit = g / cm3;
    if (suffix == "g/cm3" or suffix == "g/mL" or suffix == "g/ml") {
        unit = g / cm3;
    } else if (suffix == "kg/m3" or suffix == "kg/m^3") {
        unit = kg / m3;
    } else if (not suffix.empty()) {
        throw std::invalid_argument("unknown unit '" + suffix + "' in density '" + text + "'");
    }
    const auto density = value * unit;
    if (not std::isfinite(density) or density <= 0.0) {
        throw std::invalid_argument("invalid density '" + text + "'");
    }
    return density;
}

auto ParseMaterialComponents(const std::string& text) -> std::vector<MaterialComponent> {
    auto components = std::vector<MaterialComponent>{};
    auto stream = std::istringstream{text};
    auto componentSpec = std::string{};
    while (std::getline(stream, componentSpec, ',')) {
        if (componentSpec.empty()) {
            throw std::invalid_argument("empty component in material composition '" + text + "'");
        }
        auto component = MaterialComponent{};
        const auto colon = componentSpec.find(':');
        if (colon == std::string::npos) {
            component.mName = Trim(componentSpec);
        } else {
            component.mName = Trim(componentSpec.substr(0, colon));
            const auto quantityText = Trim(componentSpec.substr(colon + 1));
            auto parseIndex = std::size_t{0};
            try {
                component.mQuantity = std::stod(quantityText, &parseIndex);
            } catch (const std::exception&) {
                throw std::invalid_argument("cannot parse quantity '" + quantityText + "' in component '" + componentSpec + "'");
            }
            if (parseIndex != quantityText.size()) {
                throw std::invalid_argument("cannot parse quantity '" + quantityText + "' in component '" + componentSpec + "'");
            }
            component.mQuantitySpecified = true;
        }
        if (component.mName.empty()) {
            throw std::invalid_argument("empty component name in material composition '" + text + "'");
        }
        if (not std::isfinite(component.mQuantity) or component.mQuantity <= 0.0) {
            throw std::invalid_argument("component '" + componentSpec + "' has a non-positive quantity");
        }
        components.push_back(std::move(component));
    }
    if (components.empty()) {
        throw std::invalid_argument("material composition '" + text + "' is empty");
    }
    if (components.size() > 1) {
        for (const auto& component : components) {
            if (not component.mQuantitySpecified) {
                throw std::invalid_argument(
                    "component '" + component.mName + "' is missing its quantity; the quantity may be "
                                                      "omitted only when the composition has a single component");
            }
        }
    }
    return components;
}

auto ParseLayerSpec(const std::string& layerSpec) -> Layer {
    const auto firstColon = layerSpec.find(':');
    if (firstColon == std::string::npos) {
        throw std::invalid_argument(
            "layer '" + layerSpec + "' must be a NIST material 'name:thickness' or a custom material "
                                    "'name:[...|...|<...>]:density:thickness'");
    }
    auto layer = Layer{};
    layer.mMaterial.mName = Trim(layerSpec.substr(0, firstColon));
    if (layer.mMaterial.mName.empty()) {
        throw std::invalid_argument("layer '" + layerSpec + "' has an empty material name");
    }
    auto remainder = layerSpec.substr(firstColon + 1);
    const auto opening = remainder.empty() ? '\0' : remainder.front();
    if (opening == '[' or opening == '{' or opening == '<') {
        const auto closing = opening == '[' ? ']' : opening == '{' ? '}' :
                                                                     '>';
        const auto closingPos = remainder.find(closing);
        if (closingPos == std::string::npos) {
            throw std::invalid_argument("layer '" + layerSpec + "' has an unterminated composition list");
        }
        const auto compositionSpec = Trim(remainder.substr(1, closingPos - 1));
        remainder = remainder.substr(closingPos + 1);
        // The remainder has the form ':density:thickness'.
        if (remainder.empty() or remainder.front() != ':') {
            throw std::invalid_argument("layer '" + layerSpec + "' is missing its density or thickness");
        }
        const auto densityColon = remainder.find(':', 1);
        if (densityColon == std::string::npos) {
            throw std::invalid_argument("layer '" + layerSpec + "' is missing its density or thickness");
        }
        if (opening == '[') {
            layer.mMaterial.mKind = MaterialCompositionType::atoms;
        } else if (opening == '{') {
            layer.mMaterial.mKind = MaterialCompositionType::massFractions;
        } else {
            layer.mMaterial.mKind = MaterialCompositionType::materials;
        }
        layer.mMaterial.mComponents = ParseMaterialComponents(compositionSpec);
        if (layer.mMaterial.mKind == MaterialCompositionType::atoms) {
            for (const auto& component : layer.mMaterial.mComponents) {
                if (component.mQuantity != std::floor(component.mQuantity)) {
                    throw std::invalid_argument(
                        "component '" + component.mName + "' in layer '" + layerSpec + "' has a non-integer atom count");
                }
            }
        }
        layer.mMaterial.mDensity = ParseDensity(Trim(remainder.substr(1, densityColon - 1)));
        layer.mThickness = ParseLength(Trim(remainder.substr(densityColon + 1)));
    } else {
        layer.mMaterial.mKind = MaterialCompositionType::nist;
        layer.mThickness = ParseLength(Trim(remainder));
    }
    if (layer.mThickness <= 0.0) {
        throw std::invalid_argument("layer '" + layerSpec + "' has a non-positive thickness");
    }
    return layer;
}

auto ParseLayers(const std::string& text) -> std::vector<Layer> {
    auto layers = std::vector<Layer>{};
    auto stream = std::istringstream{text};
    auto layerSpec = std::string{};
    while (std::getline(stream, layerSpec, ';')) {
        if (layerSpec.empty()) {
            throw std::invalid_argument("empty layer specification in '" + text + "'");
        }
        layers.push_back(ParseLayerSpec(layerSpec));
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
        << "into the world boundary; with --output each source category is stored into its own RNTuple\n"
        << "'run{runId}_src{typeId}' inside a ROOT file, otherwise no ROOT file is produced\n"
        << "\n"
        << "required options:\n"
        << "  -m, --material <spec>    material layers as semicolon-separated entries, built in order from the\n"
        << "                           source; a layer is a NIST material 'name:thickness' or a custom material\n"
        << "                           'name:<composition>:density:thickness' with one of:\n"
        << "                             1. <composition> = [elem1:n1,elem2:n2,...]\n"
        << "                              elements by atom counts, e.g. water:[H:2,O:1]:1g/cm3:1mm\n"
        << "                             2. <composition> = {elem1:f1,elem2:f2,...}\n"
        << "                              elements by mass fractions, e.g. air:{N:0.7,O:0.3}:1g/cm3:1\n"
        << "                             3. <composition> = <mat1:f1,mat2:f2,...>\n"
        << "                              predefined materials by mass fractions, wc:<G4_CONCRETE:0.95,G4_WATER:0.05>:2.3:10 mm\n"
        << "                           a bare density number is in g/cm3; a single component may omit its atom count\n"
        << "                           or mass fraction\n"
        << "  -s, --rad-src <spec>     compound radiation source as semicolon-separated 'particle:energy[:intensity]'\n"
        << "                           entries; the intensity is the relative intensity of the source and may be\n"
        << "                           omitted only when the list contains a single source, e.g.\n"
        << "                           'e+:10MeV:15;gamma:3MeV:12' or 'neutron:1MeV'; intensities are\n"
        << "                           normalized to probabilities and every event emits one primary whose source\n"
        << "                           type is drawn from the resulting multinomial distribution\n"
        << "\n"
        << "optional options:\n"
        << "  -n, --n-event <count>    simulate <count> events in batch mode; may be combined with --ui to\n"
        << "                           pre-run events before the interactive session opens\n"
        << "  -i, --ui                 start an interactive UI session with visualization; if set, visualization\n"
        << "                           is enabled, otherwise the program runs without any UI\n"
        << "  -o, --output [<file>]    save per-event simulation results into a ROOT file, one RNTuple per\n"
        << "                           source category; the file name may be omitted (default: pps_output.root)\n"
        << "                           or given as a value, e.g. --output out.root; without this option no\n"
        << "                           ROOT file is produced (never overwrites unless --force)\n"
        << "  -l, --physics <name>     reference physics list name (default: QGSP_BIC_AllHP_EMZ)\n"
        << "  -j, --threads <count>    worker thread count; 1 runs sequential, > 1 runs multithreaded (default: all CPU cores)\n"
        << "  -v, --verbose <level>    Geant4 verbosity level; 0 prints only the banner, progress, and summary (default: 0)\n"
        << "  -N, --neutrinos          include neutrino kinetic energy in the world boundary energy statistics;\n"
        << "                           neutrinos are ignored by default\n"
        << "  -f, --force              overwrite the output file if it already exists\n"
        << "  -h, --help               print this message\n"
        << "\n"
        << "a positional macroFile executes in batch mode and cannot be combined with --n-event or --ui;\n"
        << "without --ui, --n-event, or a macroFile the program refuses to start" << G4endl;
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
        } else if (argument == "-s" or argument == "--rad-src") {
            config.mSources = ParseRadiationSources(nextValue(i, argument));
        } else if (argument == "-o" or argument == "--output") {
            config.mSave = true;
            if (i + 1 < argc and argv[i + 1][0] != '-') {
                config.mOutputFileName = nextValue(i, argument);
            }
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
            throw std::invalid_argument("required option --rad-src is missing");
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
    const auto nistManager = G4NistManager::Instance();
    for (const auto& layer : config.mLayers) {
        switch (layer.mMaterial.mKind) {
        case MaterialCompositionType::nist:
            if (nistManager->FindOrBuildMaterial(layer.mMaterial.mName) == nullptr) {
                throw std::invalid_argument(
                    "material '" + layer.mMaterial.mName + "' not found in the NIST material database");
            }
            break;
        case MaterialCompositionType::atoms:
        case MaterialCompositionType::massFractions:
            for (const auto& component : layer.mMaterial.mComponents) {
                if (nistManager->FindOrBuildElement(component.mName) == nullptr) {
                    throw std::invalid_argument(
                        "element '" + component.mName + "' of material '" + layer.mMaterial.mName + "' not found");
                }
            }
            break;
        case MaterialCompositionType::materials:
            for (const auto& component : layer.mMaterial.mComponents) {
                if (nistManager->FindOrBuildMaterial(component.mName) == nullptr) {
                    throw std::invalid_argument(
                        "material '" + component.mName + "' of material '" + layer.mMaterial.mName + "' not found in the NIST material database");
                }
            }
            break;
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
            const auto layerMaterial = BuildMaterial(layer.mMaterial);
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
    // Resolve a NIST material or build a custom material. Custom materials are cached by their
    // complete specification so that repeated layers reuse the same G4Material instead of creating
    // duplicate table entries. Custom materials are constructed as solids at NTP; the material
    // state is irrelevant for the transport processes used here.
    auto BuildMaterial(const LayerMaterial& layerMaterial) -> G4Material* {
        const auto nistManager = G4NistManager::Instance();
        if (layerMaterial.mKind == MaterialCompositionType::nist) {
            return nistManager->FindOrBuildMaterial(layerMaterial.mName);
        }
        auto key = std::ostringstream{};
        key << layerMaterial.mName << '|' << static_cast<int>(layerMaterial.mKind) << '|' << layerMaterial.mDensity;
        for (const auto& component : layerMaterial.mComponents) {
            key << '|' << component.mName << ':' << component.mQuantity;
        }
        const auto cached = mCustomMaterials.find(key.str());
        if (cached != mCustomMaterials.end()) {
            return cached->second;
        }
        auto material = new G4Material{layerMaterial.mName, layerMaterial.mDensity,
                                       static_cast<G4int>(layerMaterial.mComponents.size()), kStateSolid};
        if (layerMaterial.mKind == MaterialCompositionType::atoms) {
            for (const auto& component : layerMaterial.mComponents) {
                auto element = nistManager->FindOrBuildElement(component.mName);
                if (element == nullptr) {
                    throw std::invalid_argument("element '" + component.mName + "' not found");
                }
                material->AddElementByNumberOfAtoms(element, static_cast<G4int>(component.mQuantity));
            }
        } else if (layerMaterial.mKind == MaterialCompositionType::massFractions) {
            for (const auto& component : layerMaterial.mComponents) {
                auto element = nistManager->FindOrBuildElement(component.mName);
                if (element == nullptr) {
                    throw std::invalid_argument("element '" + component.mName + "' not found");
                }
                material->AddElementByMassFraction(element, component.mQuantity);
            }
        } else {
            for (const auto& component : layerMaterial.mComponents) {
                auto baseMaterial = nistManager->FindOrBuildMaterial(component.mName);
                if (baseMaterial == nullptr) {
                    throw std::invalid_argument("material '" + component.mName + "' not found");
                }
                material->AddMaterial(baseMaterial, component.mQuantity);
            }
        }
        mCustomMaterials.emplace(key.str(), material);
        return material;
    }

    const Config& mConfig;
    std::unordered_map<const G4LogicalVolume*, int> mMaterialVolumes;
    std::unordered_map<std::string, G4Material*> mCustomMaterials;
};

class SimulationRun : public G4Run {
public:
    struct SourceStatistics {
        explicit SourceStatistics(int layerCount) :
            mLayerDepositedEnergySums(layerCount, 0.0),
            mLayerDepositedEnergySumsSq(layerCount, 0.0),
            mLayerDepositedPrimaryCounts(layerCount, 0) {}
        long long mEventCount{0};
        double mPenetratingEnergySum{0.0};
        double mPenetratingEnergySumSq{0.0};
        double mDepositedEnergySum{0.0};
        double mDepositedEnergySumSq{0.0};
        std::vector<double> mLayerDepositedEnergySums;
        std::vector<double> mLayerDepositedEnergySumsSq;
        double mBackscatteredEnergySum{0.0};
        double mBackscatteredEnergySumSq{0.0};
        long long mPenetratingPrimaryCount{0};
        std::vector<long long> mLayerDepositedPrimaryCounts;
        long long mBackscatteredPrimaryCount{0};
        std::map<std::string, double> mPenetratingParticleCount;
        std::map<std::string, double> mPenetratingParticleEnergySum;
        std::map<std::string, double> mPenetratingParticleEnergySumSq;
        std::map<std::string, double> mBackscatteredParticleCount;
        std::map<std::string, double> mBackscatteredParticleEnergySum;
        std::map<std::string, double> mBackscatteredParticleEnergySumSq;
    };

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
            target.mPenetratingEnergySum += source.mPenetratingEnergySum;
            target.mPenetratingEnergySumSq += source.mPenetratingEnergySumSq;
            target.mDepositedEnergySum += source.mDepositedEnergySum;
            target.mDepositedEnergySumSq += source.mDepositedEnergySumSq;
            for (auto layerIndex{0}; layerIndex < static_cast<int>(target.mLayerDepositedEnergySums.size()); ++layerIndex) {
                target.mLayerDepositedEnergySums[layerIndex] += source.mLayerDepositedEnergySums[layerIndex];
                target.mLayerDepositedEnergySumsSq[layerIndex] += source.mLayerDepositedEnergySumsSq[layerIndex];
            }
            target.mBackscatteredEnergySum += source.mBackscatteredEnergySum;
            target.mBackscatteredEnergySumSq += source.mBackscatteredEnergySumSq;
            target.mPenetratingPrimaryCount += source.mPenetratingPrimaryCount;
            for (auto layerIndex{0}; layerIndex < static_cast<int>(target.mLayerDepositedPrimaryCounts.size()); ++layerIndex) {
                target.mLayerDepositedPrimaryCounts[layerIndex] += source.mLayerDepositedPrimaryCounts[layerIndex];
            }
            target.mBackscatteredPrimaryCount += source.mBackscatteredPrimaryCount;
            AddToMap(target.mPenetratingParticleCount, source.mPenetratingParticleCount);
            AddToMap(target.mPenetratingParticleEnergySum, source.mPenetratingParticleEnergySum);
            AddToMap(target.mPenetratingParticleEnergySumSq, source.mPenetratingParticleEnergySumSq);
            AddToMap(target.mBackscatteredParticleCount, source.mBackscatteredParticleCount);
            AddToMap(target.mBackscatteredParticleEnergySum, source.mBackscatteredParticleEnergySum);
            AddToMap(target.mBackscatteredParticleEnergySumSq, source.mBackscatteredParticleEnergySumSq);
        }
    }

    auto AddEventResult(int sourceIndex, double penetratingEnergy, const std::vector<double>& layerDepositedEnergies,
                        double backscatteredEnergy) -> void {
        auto& statistics = mSourceStatistics[sourceIndex];
        ++statistics.mEventCount;
        statistics.mPenetratingEnergySum += penetratingEnergy;
        statistics.mPenetratingEnergySumSq += penetratingEnergy * penetratingEnergy;
        auto totalDepositedEnergy{0.0};
        for (auto layerIndex{0}; layerIndex < static_cast<int>(layerDepositedEnergies.size()); ++layerIndex) {
            const auto layerDepositedEnergy = layerDepositedEnergies[layerIndex];
            statistics.mLayerDepositedEnergySums[layerIndex] += layerDepositedEnergy;
            statistics.mLayerDepositedEnergySumsSq[layerIndex] += layerDepositedEnergy * layerDepositedEnergy;
            totalDepositedEnergy += layerDepositedEnergy;
        }
        statistics.mDepositedEnergySum += totalDepositedEnergy;
        statistics.mDepositedEnergySumSq += totalDepositedEnergy * totalDepositedEnergy;
        statistics.mBackscatteredEnergySum += backscatteredEnergy;
        statistics.mBackscatteredEnergySumSq += backscatteredEnergy * backscatteredEnergy;
    }

    auto AddPenetratingEnergy(int sourceIndex, const std::string& particleName, double energy) -> void {
        auto& statistics = mSourceStatistics[sourceIndex];
        statistics.mPenetratingParticleCount[particleName] += 1.0;
        statistics.mPenetratingParticleEnergySum[particleName] += energy;
        statistics.mPenetratingParticleEnergySumSq[particleName] += energy * energy;
    }

    auto AddBackscatteredEnergy(int sourceIndex, const std::string& particleName, double energy) -> void {
        auto& statistics = mSourceStatistics[sourceIndex];
        statistics.mBackscatteredParticleCount[particleName] += 1.0;
        statistics.mBackscatteredParticleEnergySum[particleName] += energy;
        statistics.mBackscatteredParticleEnergySumSq[particleName] += energy * energy;
    }

    auto AddPrimaryPenetrating(int sourceIndex) -> void {
        ++mSourceStatistics[sourceIndex].mPenetratingPrimaryCount;
    }

    auto AddPrimaryDeposited(int sourceIndex, int layerIndex) -> void {
        ++mSourceStatistics[sourceIndex].mLayerDepositedPrimaryCounts[layerIndex];
    }

    auto AddPrimaryBackscattered(int sourceIndex) -> void {
        ++mSourceStatistics[sourceIndex].mBackscatteredPrimaryCount;
    }

    auto GetStatistics(int sourceIndex) const -> const SourceStatistics& {
        return mSourceStatistics[sourceIndex];
    }

    auto GetTotalStatistics() const -> SourceStatistics {
        auto total = SourceStatistics{mLayerCount};
        for (const auto& statistics : mSourceStatistics) {
            total.mEventCount += statistics.mEventCount;
            total.mPenetratingEnergySum += statistics.mPenetratingEnergySum;
            total.mPenetratingEnergySumSq += statistics.mPenetratingEnergySumSq;
            total.mDepositedEnergySum += statistics.mDepositedEnergySum;
            total.mDepositedEnergySumSq += statistics.mDepositedEnergySumSq;
            for (auto layerIndex{0}; layerIndex < mLayerCount; ++layerIndex) {
                total.mLayerDepositedEnergySums[layerIndex] += statistics.mLayerDepositedEnergySums[layerIndex];
                total.mLayerDepositedEnergySumsSq[layerIndex] += statistics.mLayerDepositedEnergySumsSq[layerIndex];
            }
            total.mBackscatteredEnergySum += statistics.mBackscatteredEnergySum;
            total.mBackscatteredEnergySumSq += statistics.mBackscatteredEnergySumSq;
            total.mPenetratingPrimaryCount += statistics.mPenetratingPrimaryCount;
            for (auto layerIndex{0}; layerIndex < mLayerCount; ++layerIndex) {
                total.mLayerDepositedPrimaryCounts[layerIndex] += statistics.mLayerDepositedPrimaryCounts[layerIndex];
            }
            total.mBackscatteredPrimaryCount += statistics.mBackscatteredPrimaryCount;
            AddToMap(total.mPenetratingParticleCount, statistics.mPenetratingParticleCount);
            AddToMap(total.mPenetratingParticleEnergySum, statistics.mPenetratingParticleEnergySum);
            AddToMap(total.mPenetratingParticleEnergySumSq, statistics.mPenetratingParticleEnergySumSq);
            AddToMap(total.mBackscatteredParticleCount, statistics.mBackscatteredParticleCount);
            AddToMap(total.mBackscatteredParticleEnergySum, statistics.mBackscatteredParticleEnergySum);
            AddToMap(total.mBackscatteredParticleEnergySumSq, statistics.mBackscatteredParticleEnergySumSq);
        }
        return total;
    }

    auto GetLayerCount() const -> int {
        return mLayerCount;
    }
    auto GetSourceCount() const -> int {
        return static_cast<int>(mSourceStatistics.size());
    }

private:
    static auto AddToMap(std::map<std::string, double>& target, const std::map<std::string, double>& source) -> void {
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
        if (IsMaster() and mConfig.mSave) {
            OutputWriter::Instance().BeginRun(run->GetRunID(), static_cast<int>(mConfig.mLayers.size()),
                                              static_cast<int>(mConfig.mSources.size()));
        }
        if (mConfig.mSave and (not IsMaster() or not G4Threading::IsMultithreadedApplication())) {
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
            if (mConfig.mSave) {
                OutputWriter::Instance().EndRun();
            }
            PrintStatistics(static_cast<const SimulationRun&>(*run));
            G4cout << "run " << run->GetRunID() << " finished with " << run->GetNumberOfEvent() << " events" << G4endl;
        }
    }

    auto AddEventResult(int sourceIndex, double penetratingEnergy, const std::vector<double>& layerDepositedEnergies,
                        double backscatteredEnergy) -> void {
        mCurrentRun->AddEventResult(sourceIndex, penetratingEnergy, layerDepositedEnergies, backscatteredEnergy);
    }

    auto AddPenetratingEnergy(int sourceIndex, const std::string& particleName, double energy) -> void {
        mCurrentRun->AddPenetratingEnergy(sourceIndex, particleName, energy);
    }

    auto AddBackscatteredEnergy(int sourceIndex, const std::string& particleName, double energy) -> void {
        mCurrentRun->AddBackscatteredEnergy(sourceIndex, particleName, energy);
    }

    auto AddPrimaryPenetrating(int sourceIndex) -> void {
        mCurrentRun->AddPrimaryPenetrating(sourceIndex);
    }

    auto AddPrimaryDeposited(int sourceIndex, int layerIndex) -> void {
        mCurrentRun->AddPrimaryDeposited(sourceIndex, layerIndex);
    }

    auto AddPrimaryBackscattered(int sourceIndex) -> void {
        mCurrentRun->AddPrimaryBackscattered(sourceIndex);
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
        // Effective incident energy of the compound source (intensity-weighted mean of the
        // per-source energies), used to normalize the aggregated total ratios.
        auto effectiveIncidentEnergy{0.0};
        for (const auto& source : mConfig.mSources) {
            effectiveIncidentEnergy += source.mWeight * source.mEnergy;
        }
        auto printedAnyBlock{false};
        for (auto sourceIndex{0}; sourceIndex < run.GetSourceCount(); ++sourceIndex) {
            const auto& source = mConfig.mSources[sourceIndex];
            const auto& statistics = run.GetStatistics(sourceIndex);
            if (statistics.mEventCount < 1) {
                continue;
            }
            if (printedAnyBlock) {
                // Separator between consecutive source blocks.
                G4cout << "-------------------------------------------------------------------------------\n";
            }
            printedAnyBlock = true;
            auto header = std::ostringstream{};
            header << "source " << sourceIndex << ": " << source.mParticleName << ' '
                   << G4BestUnit(source.mEnergy, "Energy") << ", intensity "
                   << std::fixed << std::setprecision(4) << source.mWeight << ", "
                   << statistics.mEventCount << " events";
            PrintSourceBlock(header.str(), statistics, source.mEnergy);
        }
        // Aggregated total block is printed last.
        G4cout << "-------------------------------------------------------------------------------\n";
        PrintSourceBlock("total: " + std::to_string(mConfig.mSources.size()) + " components, " + std::to_string(eventCount) + " events",
                         run.GetTotalStatistics(), effectiveIncidentEnergy);
        G4cout << "===============================================================================\n"
               << G4endl;
    }

    auto PrintSourceBlock(const std::string& header, const SimulationRun::SourceStatistics& statistics,
                          double incidentEnergy) -> void {
        G4cout << ' ' << header << ":\n";
        G4cout << std::defaultfloat << std::setprecision(6);
        PrintEnergyRatio("energy penetration ratio", statistics.mPenetratingEnergySum,
                         statistics.mPenetratingEnergySumSq, statistics.mEventCount, incidentEnergy);
        PrintEnergyRatio("energy deposition ratio", statistics.mDepositedEnergySum,
                         statistics.mDepositedEnergySumSq, statistics.mEventCount, incidentEnergy);
        for (auto layerIndex{0}; layerIndex < static_cast<int>(mConfig.mLayers.size()); ++layerIndex) {
            const auto label = "  in layer " + std::to_string(layerIndex) + " (" +
                               mConfig.mLayers[layerIndex].mMaterial.mName + ')';
            PrintEnergyRatio(label, statistics.mLayerDepositedEnergySums[layerIndex],
                             statistics.mLayerDepositedEnergySumsSq[layerIndex], statistics.mEventCount,
                             incidentEnergy);
        }
        PrintEnergyRatio("energy backscattering ratio", statistics.mBackscatteredEnergySum,
                         statistics.mBackscatteredEnergySumSq, statistics.mEventCount, incidentEnergy);
        // Primary-particle termination statistics: the fraction of primaries that end up before,
        // inside, or after the material. For unstable primaries the termination position includes
        // decay. The counts are mutually exclusive and sum to the event count of the block.
        auto totalDepositedPrimaryCount{0LL};
        for (auto layerIndex{0}; layerIndex < static_cast<int>(mConfig.mLayers.size()); ++layerIndex) {
            totalDepositedPrimaryCount += statistics.mLayerDepositedPrimaryCounts[layerIndex];
        }
        PrintPrimaryRatio("primary particle penetration ratio", statistics.mPenetratingPrimaryCount,
                          statistics.mEventCount);
        PrintPrimaryRatio("primary particle deposition ratio", totalDepositedPrimaryCount, statistics.mEventCount);
        for (auto layerIndex{0}; layerIndex < static_cast<int>(mConfig.mLayers.size()); ++layerIndex) {
            const auto label = "  in layer " + std::to_string(layerIndex) + " (" +
                               mConfig.mLayers[layerIndex].mMaterial.mName + ')';
            PrintPrimaryRatio(label, statistics.mLayerDepositedPrimaryCounts[layerIndex], statistics.mEventCount);
        }
        PrintPrimaryRatio("primary particle backscattering ratio", statistics.mBackscatteredPrimaryCount,
                          statistics.mEventCount);
        PrintParticleStatistics("penetrating particles", statistics.mPenetratingParticleCount,
                                statistics.mPenetratingParticleEnergySum,
                                statistics.mPenetratingParticleEnergySumSq);
        PrintParticleStatistics("backscattered particles", statistics.mBackscatteredParticleCount,
                                statistics.mBackscatteredParticleEnergySum,
                                statistics.mBackscatteredParticleEnergySumSq);
    }

    auto PrintPrimaryRatio(const std::string& label, long long count, long long eventCount) -> void {
        if (eventCount < 1) {
            return;
        }
        const auto proportion{static_cast<double>(count) / eventCount};
        // Standard error of the binomial proportion.
        const auto proportionError{std::sqrt(proportion * (1.0 - proportion) / eventCount)};
        G4cout << "   " << std::left << std::setw(45) << label + ':' << '(' << 100.0 * proportion << " +/- "
               << 100.0 * proportionError << ") %" << G4endl;
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
            const auto rmsEnergy = std::sqrt(variance);
            G4cout << "     "
                   << std::setw(14) << particleName
                   << std::setw(20) << FormatEnergy(meanEnergy)
                   << std::setw(20) << FormatEnergy(rmsEnergy)
                   << std::setw(8) << static_cast<long long>(count) << '\n';
        }
    }

    auto PrintEnergyRatio(const std::string& label, double sumEnergy, double sumEnergySq, long long eventCount,
                          double incidentEnergy) -> void {
        const auto meanEnergy{sumEnergy / eventCount};
        const auto ratio{meanEnergy / incidentEnergy};
        const auto variance{(sumEnergySq - sumEnergy * sumEnergy / eventCount) / (eventCount - 1)};
        const auto ratioError{std::sqrt(variance / eventCount) / incidentEnergy};
        G4cout << "   " << std::left << std::setw(45) << label + ':' << '(' << 100.0 * ratio << " +/- "
               << 100.0 * ratioError << ") %" << G4endl;
    }

    const Config& mConfig;
    SimulationRun* mCurrentRun{nullptr};
    std::vector<std::shared_ptr<ROOT::RNTupleFillContext>> mFillContexts;
};

class EventAction : public G4UserEventAction {
public:
    explicit EventAction(RunAction& runAction, int layerCount, int sourceCount, bool saveResults) :
        mRunAction{runAction},
        mLayerCount{layerCount},
        mSourceCount{sourceCount},
        mSaveResults{saveResults},
        mFillContexts(sourceCount),
        mEntries(sourceCount),
        mFillStatuses(sourceCount),
        mFields(sourceCount),
        mLayerEnergyDeposit(layerCount, 0.0),
        mDepositedParticles(layerCount),
        mDepositedX(layerCount),
        mDepositedY(layerCount),
        mDepositedZ(layerCount),
        mDepositedWeight(layerCount),
        mDepositedProcess(layerCount) {}
    ~EventAction() override = default;

    auto SetSourceIndex(int sourceIndex) -> void {
        mCurrentSourceIndex = sourceIndex;
    }

    auto BeginOfEventAction(const G4Event*) -> void override {
        const auto runId = G4RunManager::GetRunManager()->GetCurrentRun()->GetRunID();
        if (mSaveResults and runId != mRunId) {
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
        mTotalPenetratingEnergy = 0.0;
        mPenetratingParticles.clear();
        mPenetratingTheta.clear();
        mPenetratingPhi.clear();
        mPenetratingEnergy.clear();
        mTotalDepositedEnergy = 0.0;
        for (auto& layerEnergyDeposit : mLayerEnergyDeposit) {
            layerEnergyDeposit = 0.0;
        }
        for (auto& depositedParticles : mDepositedParticles) {
            depositedParticles.clear();
        }
        for (auto& depositedX : mDepositedX) {
            depositedX.clear();
        }
        for (auto& depositedY : mDepositedY) {
            depositedY.clear();
        }
        for (auto& depositedZ : mDepositedZ) {
            depositedZ.clear();
        }
        for (auto& depositedWeight : mDepositedWeight) {
            depositedWeight.clear();
        }
        for (auto& depositedProcess : mDepositedProcess) {
            depositedProcess.clear();
        }
        mTotalBackscatteredEnergy = 0.0;
        mBackscatteredParticles.clear();
        mBackscatteredTheta.clear();
        mBackscatteredPhi.clear();
        mBackscatteredEnergy.clear();
    }

    auto EndOfEventAction(const G4Event* event) -> void override {
        if (mSaveResults) {
            const auto sourceIndex = mCurrentSourceIndex;
            auto& fields = mFields[sourceIndex];
            *fields.mEventId = event->GetEventID();
            // All per-event computations are performed in double; the values are narrowed to the
            // float RNTuple fields exactly here, at the storage boundary.
            *fields.mTotalPenetratingEnergy = static_cast<float>(mTotalPenetratingEnergy);
            *fields.mPenetratingParticles = mPenetratingParticles;
            fields.mPenetratingTheta->assign(mPenetratingTheta.begin(), mPenetratingTheta.end());
            fields.mPenetratingPhi->assign(mPenetratingPhi.begin(), mPenetratingPhi.end());
            fields.mPenetratingEnergy->assign(mPenetratingEnergy.begin(), mPenetratingEnergy.end());
            *fields.mTotalDepositedEnergy = static_cast<float>(mTotalDepositedEnergy);
            for (auto layerIndex{0}; layerIndex < mLayerCount; ++layerIndex) {
                *fields.mLayerEnergyDeposit[layerIndex] = static_cast<float>(mLayerEnergyDeposit[layerIndex]);
                *fields.mDepositedParticles[layerIndex] = mDepositedParticles[layerIndex];
                fields.mDepositedX[layerIndex]->assign(mDepositedX[layerIndex].begin(), mDepositedX[layerIndex].end());
                fields.mDepositedY[layerIndex]->assign(mDepositedY[layerIndex].begin(), mDepositedY[layerIndex].end());
                fields.mDepositedZ[layerIndex]->assign(mDepositedZ[layerIndex].begin(), mDepositedZ[layerIndex].end());
                fields.mDepositedWeight[layerIndex]->assign(mDepositedWeight[layerIndex].begin(),
                                                            mDepositedWeight[layerIndex].end());
                *fields.mDepositedProcess[layerIndex] = mDepositedProcess[layerIndex];
            }
            *fields.mTotalBackscatteredEnergy = static_cast<float>(mTotalBackscatteredEnergy);
            *fields.mBackscatteredParticles = mBackscatteredParticles;
            fields.mBackscatteredTheta->assign(mBackscatteredTheta.begin(), mBackscatteredTheta.end());
            fields.mBackscatteredPhi->assign(mBackscatteredPhi.begin(), mBackscatteredPhi.end());
            fields.mBackscatteredEnergy->assign(mBackscatteredEnergy.begin(), mBackscatteredEnergy.end());
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
        }
        mRunAction.AddEventResult(mCurrentSourceIndex, mTotalPenetratingEnergy, mLayerEnergyDeposit,
                                  mTotalBackscatteredEnergy);
    }

    auto AddPenetratingParticle(const std::string& particleName, const G4ThreeVector& direction, double energy) -> void {
        AddExitPoint(mPenetratingParticles, mPenetratingTheta, mPenetratingPhi, mPenetratingEnergy,
                     mTotalPenetratingEnergy, particleName, direction, energy);
        mRunAction.AddPenetratingEnergy(mCurrentSourceIndex, particleName, energy);
    }

    auto AddDepositedEnergy(int layerIndex, const std::string& particleName, const G4ThreeVector& position,
                            double energyDeposit, const std::string& processName) -> void {
        mTotalDepositedEnergy += energyDeposit;
        mLayerEnergyDeposit[layerIndex] += energyDeposit;
        mDepositedParticles[layerIndex].push_back(particleName);
        mDepositedX[layerIndex].push_back(position.x());
        mDepositedY[layerIndex].push_back(position.y());
        mDepositedZ[layerIndex].push_back(position.z());
        mDepositedWeight[layerIndex].push_back(energyDeposit);
        mDepositedProcess[layerIndex].push_back(processName);
    }

    auto AddBackscatteredParticle(const std::string& particleName, const G4ThreeVector& direction, double energy) -> void {
        AddExitPoint(mBackscatteredParticles, mBackscatteredTheta, mBackscatteredPhi, mBackscatteredEnergy,
                     mTotalBackscatteredEnergy, particleName, direction, energy);
        mRunAction.AddBackscatteredEnergy(mCurrentSourceIndex, particleName, energy);
    }

    auto AddPrimaryPenetrating() -> void {
        mRunAction.AddPrimaryPenetrating(mCurrentSourceIndex);
    }

    auto AddPrimaryDeposited(int layerIndex) -> void {
        mRunAction.AddPrimaryDeposited(mCurrentSourceIndex, layerIndex);
    }

    auto AddPrimaryBackscattered() -> void {
        mRunAction.AddPrimaryBackscattered(mCurrentSourceIndex);
    }

private:
    struct SourceFields {
        std::shared_ptr<int> mEventId;
        std::shared_ptr<float> mTotalPenetratingEnergy;
        std::shared_ptr<std::vector<std::string>> mPenetratingParticles;
        std::shared_ptr<std::vector<float>> mPenetratingTheta;
        std::shared_ptr<std::vector<float>> mPenetratingPhi;
        std::shared_ptr<std::vector<float>> mPenetratingEnergy;
        std::shared_ptr<float> mTotalDepositedEnergy;
        std::vector<std::shared_ptr<float>> mLayerEnergyDeposit;
        std::vector<std::shared_ptr<std::vector<std::string>>> mDepositedParticles;
        std::vector<std::shared_ptr<std::vector<float>>> mDepositedX;
        std::vector<std::shared_ptr<std::vector<float>>> mDepositedY;
        std::vector<std::shared_ptr<std::vector<float>>> mDepositedZ;
        std::vector<std::shared_ptr<std::vector<float>>> mDepositedWeight;
        std::vector<std::shared_ptr<std::vector<std::string>>> mDepositedProcess;
        std::shared_ptr<float> mTotalBackscatteredEnergy;
        std::shared_ptr<std::vector<std::string>> mBackscatteredParticles;
        std::shared_ptr<std::vector<float>> mBackscatteredTheta;
        std::shared_ptr<std::vector<float>> mBackscatteredPhi;
        std::shared_ptr<std::vector<float>> mBackscatteredEnergy;
    };

    auto BindFields(ROOT::REntry& entry, SourceFields& fields) -> void {
        fields.mEventId = entry.GetPtr<int>("event_id");
        fields.mTotalPenetratingEnergy = entry.GetPtr<float>("total_e_pen");
        fields.mPenetratingParticles = entry.GetPtr<std::vector<std::string>>("particle_pen");
        fields.mPenetratingTheta = entry.GetPtr<std::vector<float>>("theta_pen");
        fields.mPenetratingPhi = entry.GetPtr<std::vector<float>>("phi_pen");
        fields.mPenetratingEnergy = entry.GetPtr<std::vector<float>>("e_pen");
        fields.mTotalDepositedEnergy = entry.GetPtr<float>("total_e_dep");
        fields.mLayerEnergyDeposit.clear();
        fields.mDepositedParticles.clear();
        fields.mDepositedX.clear();
        fields.mDepositedY.clear();
        fields.mDepositedZ.clear();
        fields.mDepositedWeight.clear();
        fields.mDepositedProcess.clear();
        for (auto layerIndex{0}; layerIndex < mLayerCount; ++layerIndex) {
            const auto suffix = std::to_string(layerIndex);
            fields.mLayerEnergyDeposit.push_back(entry.GetPtr<float>("e_dep_" + suffix));
            fields.mDepositedParticles.push_back(entry.GetPtr<std::vector<std::string>>("particle_dep_" + suffix));
            fields.mDepositedX.push_back(entry.GetPtr<std::vector<float>>("x_dep_" + suffix));
            fields.mDepositedY.push_back(entry.GetPtr<std::vector<float>>("y_dep_" + suffix));
            fields.mDepositedZ.push_back(entry.GetPtr<std::vector<float>>("z_dep_" + suffix));
            fields.mDepositedWeight.push_back(entry.GetPtr<std::vector<float>>("w_dep_" + suffix));
            fields.mDepositedProcess.push_back(entry.GetPtr<std::vector<std::string>>("proc_dep_" + suffix));
        }
        fields.mTotalBackscatteredEnergy = entry.GetPtr<float>("total_e_bsc");
        fields.mBackscatteredParticles = entry.GetPtr<std::vector<std::string>>("particle_bsc");
        fields.mBackscatteredTheta = entry.GetPtr<std::vector<float>>("theta_bsc");
        fields.mBackscatteredPhi = entry.GetPtr<std::vector<float>>("phi_bsc");
        fields.mBackscatteredEnergy = entry.GetPtr<std::vector<float>>("e_bsc");
    }

    auto AddExitPoint(std::vector<std::string>& particleNames, std::vector<double>& thetas, std::vector<double>& phis,
                      std::vector<double>& energies, double& totalEnergy, const std::string& particleName,
                      const G4ThreeVector& direction, double energy) -> void {
        totalEnergy += energy;
        particleNames.push_back(particleName);
        const auto z = std::clamp(direction.z(), -1.0, 1.0);
        thetas.push_back(std::acos(z));
        phis.push_back(std::atan2(direction.y(), direction.x()));
        energies.push_back(energy);
    }

    RunAction& mRunAction;
    int mLayerCount{0};
    int mSourceCount{0};
    bool mSaveResults{false};
    int mCurrentSourceIndex{0};
    int mRunId{-1};
    std::vector<std::weak_ptr<ROOT::RNTupleFillContext>> mFillContexts;
    std::vector<std::unique_ptr<ROOT::REntry>> mEntries;
    std::vector<ROOT::RNTupleFillStatus> mFillStatuses;
    std::vector<SourceFields> mFields;
    double mTotalPenetratingEnergy{0.0};
    std::vector<std::string> mPenetratingParticles;
    std::vector<double> mPenetratingTheta;
    std::vector<double> mPenetratingPhi;
    std::vector<double> mPenetratingEnergy;
    double mTotalDepositedEnergy{0.0};
    std::vector<double> mLayerEnergyDeposit;
    std::vector<std::vector<std::string>> mDepositedParticles;
    std::vector<std::vector<double>> mDepositedX;
    std::vector<std::vector<double>> mDepositedY;
    std::vector<std::vector<double>> mDepositedZ;
    std::vector<std::vector<double>> mDepositedWeight;
    std::vector<std::vector<std::string>> mDepositedProcess;
    double mTotalBackscatteredEnergy{0.0};
    std::vector<std::string> mBackscatteredParticles;
    std::vector<double> mBackscatteredTheta;
    std::vector<double> mBackscatteredPhi;
    std::vector<double> mBackscatteredEnergy;
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
        mEventAction.AddDepositedEnergy(layerIndex, step->GetTrack()->GetDefinition()->GetParticleName(),
                                        step->GetPostStepPoint()->GetPosition(), energyDeposit, processName);
    }

private:
    EventAction& mEventAction;
    std::unordered_map<const G4LogicalVolume*, int> mMaterialVolumes;
};

class TrackingAction : public G4UserTrackingAction {
public:
    explicit TrackingAction(EventAction& eventAction, const Config& config) :
        mEventAction{eventAction},
        mConfig{config} {}
    ~TrackingAction() override = default;

    auto PostUserTrackingAction(const G4Track* track) -> void override {
        const auto step = track->GetStep();
        const auto postStepPoint = step->GetPostStepPoint();
        // Primary-particle termination statistics (track id 1): classify where the primary ends up.
        // For unstable primaries this includes decay, because the final step position is the decay
        // vertex. The categories are mutually exclusive and exhaustive by the final z coordinate:
        //   z < 0                    -> terminated before the material (backscattered)
        //   0 <= z < totalThickness  -> terminated inside a material layer (deposited)
        //   z >= totalThickness      -> terminated after the material (penetrating)
        if (track->GetTrackID() == 1) {
            const auto z = postStepPoint->GetPosition().z();
            const auto totalThickness = mConfig.TotalThickness();
            if (z < 0.0) {
                mEventAction.AddPrimaryBackscattered();
            } else if (z >= totalThickness) {
                mEventAction.AddPrimaryPenetrating();
            } else {
                auto cumulative{0.0};
                auto layerIndex{0};
                for (auto i{0}; i < static_cast<int>(mConfig.mLayers.size()); ++i) {
                    cumulative += mConfig.mLayers[i].mThickness;
                    if (z < cumulative) {
                        layerIndex = i;
                        break;
                    }
                }
                mEventAction.AddPrimaryDeposited(layerIndex);
            }
        }
        if (postStepPoint->GetStepStatus() != fWorldBoundary) {
            return;
        }
        const auto definition = track->GetDefinition();
        if (not mConfig.mIncludeNeutrinos and IsNeutrino(definition)) {
            return;
        }
        const auto direction = track->GetMomentumDirection();
        const auto energy = track->GetKineticEnergy();
        const auto particleName = definition->GetParticleName();
        if (direction.z() >= 0.0) {
            mEventAction.AddPenetratingParticle(particleName, direction, energy);
        } else {
            mEventAction.AddBackscatteredParticle(particleName, direction, energy);
        }
    }

private:
    EventAction& mEventAction;
    const Config& mConfig;
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
                                                   static_cast<int>(mConfig.mSources.size()), mConfig.mSave};
        SetUserAction(eventAction);
        SetUserAction(new PrimaryGeneratorAction{mConfig, *eventAction});
        SetUserAction(new SteppingAction{*eventAction});
        SetUserAction(new TrackingAction{*eventAction, mConfig});
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

    if (config.mSave) {
        PPS::OutputWriter::Instance().Initialize(config.mOutputFileName, config.mOverwrite);
    }
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
    if (config.mSave) {
        PPS::OutputWriter::Instance().Finalize();
    }

    return EXIT_SUCCESS;
} catch (const std::exception& error) {
    G4cerr << "error: " << error.what() << G4endl;
    // Safe even when saving was disabled (no-op unless a ROOT file is open).
    PPS::OutputWriter::Instance().Finalize();
    std::quick_exit(EXIT_FAILURE);
}
