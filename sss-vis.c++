// sss-vis — read a ROOT file written by sss and build energy-deposition histograms.
//
// For every RNTuple in the input ROOT file (one per source category, named
// "run{runId}_src{typeId}"), this program builds a 3D histogram of the energy
// deposition positions weighted by the deposited energy, the 1D histogram of the
// energy deposition along z, and the three 2D projections onto the xy, xz, and
// yz planes. The histograms are stored into an output ROOT file.
//
// The data is read through ROOT::RDataFrame, which reads the RNTuple directly.
// The deposition positions of all material layers are combined into single
// columns with chained Define() calls, and the statistics, the fine xy
// histograms used to derive the automatic xy ranges, and the final histograms
// are booked as lazy RDataFrame actions: the statistics trigger one event loop,
// the fine xy histograms trigger a second one, and the final histograms (whose
// ranges depend on those intermediate results) trigger a third one when they are
// written to the output file. ROOT implicit multi-threading (IMT) is enabled so
// that the event loops run in parallel on all worker threads.

#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TFile.h"
#include "TKey.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TROOT.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <memory>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace SSSVis {

using PositionVector = ROOT::RVec<float>;

// Names of the columns created by the chained Define() calls that combine the
// per-layer deposition data of an RNTuple into a single view over all layers.
constexpr auto xPositionColumnName = "dep_x";
constexpr auto yPositionColumnName = "dep_y";
constexpr auto zPositionColumnName = "dep_z";
constexpr auto weightColumnName = "dep_w";
constexpr auto depositionCountColumnName = "dep_count";

// Bin count of the temporary fine 1D histograms used to derive the automatic
// energy-fraction xy ranges; 1000 bins resolve the range to 0.1% of the data
// extent.
constexpr auto fineHistogramBins = 1000;

// Upper bound on the number of material layers accepted from an input RNTuple;
// protects against crafted column names that would otherwise drive unbounded
// column-list allocation and Define chaining.
constexpr auto maximumLayerCount = 256;

struct Range {
    double mMinimum{0.0};
    double mMaximum{0.0};
};

struct ThreeDimensionalBinning {
    int mX{30};
    int mY{30};
    int mZ{30};
};

struct TwoDimensionalBinning {
    int mX{100};
    int mY{100};
};

struct Config {
    std::string mInputFileName{"sss_output.root"};
    std::string mOutputFileName;
    ThreeDimensionalBinning mThreeDimensionalBinning;
    TwoDimensionalBinning mTwoDimensionalBinning;
    int mOneDimensionalBins{300};
    std::optional<Range> mXRange;
    std::optional<Range> mYRange;
    std::optional<Range> mZRange;
    double mXyEnergyFraction{0.9};
    double mZExpandFactor{1.2};
    int mThreads{0};
    bool mForce{false};
    bool mHelp{false};
};

struct DepositionStatistics {
    double mWeightSum{0.0};
    long long mDepositionCount{0};
    double mXMinimum{0.0};
    double mXMaximum{0.0};
    double mYMinimum{0.0};
    double mYMaximum{0.0};
    double mZMinimum{0.0};
    double mZMaximum{0.0};
};

struct HistogramRanges {
    Range mX;
    Range mY;
    Range mZ;
};

struct DepositionHistograms {
    ROOT::RDF::RResultPtr<TH1D> mZ;
    ROOT::RDF::RResultPtr<TH3D> mThreeDimensional;
    ROOT::RDF::RResultPtr<TH2D> mXy;
    ROOT::RDF::RResultPtr<TH2D> mXz;
    ROOT::RDF::RResultPtr<TH2D> mYz;
};

auto DefaultOutputFileName(const std::string& inputFileName) -> std::string {
    const auto dotPosition = inputFileName.find_last_of('.');
    if (dotPosition == std::string::npos) {
        return inputFileName + "_vis";
    }
    return inputFileName.substr(0, dotPosition) + "_vis" + inputFileName.substr(dotPosition);
}

auto ParseDouble(const std::string& text) -> double {
    auto parsedIndex = std::size_t{0};
    auto value = 0.0;
    try {
        value = std::stod(text, &parsedIndex);
    } catch (const std::exception&) {
        throw std::invalid_argument("cannot parse number '" + text + "'");
    }
    if (parsedIndex != text.size()) {
        throw std::invalid_argument("cannot parse number '" + text + "'");
    }
    return value;
}

auto ParseInteger(const std::string& text) -> int {
    auto parsedIndex = std::size_t{0};
    auto value = 0;
    try {
        value = std::stoi(text, &parsedIndex);
    } catch (const std::exception&) {
        throw std::invalid_argument("cannot parse integer '" + text + "'");
    }
    if (parsedIndex != text.size()) {
        throw std::invalid_argument("cannot parse integer '" + text + "'");
    }
    return value;
}

auto ParseRange(const std::string& text) -> Range {
    const auto colon = text.find(':');
    if (colon == std::string::npos) {
        throw std::invalid_argument("range '" + text + "' must have the form 'min:max'");
    }
    const auto minimum = ParseDouble(text.substr(0, colon));
    const auto maximum = ParseDouble(text.substr(colon + 1));
    if (not std::isfinite(minimum) or not std::isfinite(maximum)) {
        throw std::invalid_argument("range '" + text + "' contains a non-finite value");
    }
    if (minimum >= maximum) {
        throw std::invalid_argument("range '" + text + "' has a minimum not below its maximum");
    }
    return Range{minimum, maximum};
}

auto ParseIntegerList(const std::string& text, std::size_t expectedCount) -> std::vector<int> {
    auto values = std::vector<int>{};
    auto stream = std::istringstream{text};
    auto token = std::string{};
    while (std::getline(stream, token, ':')) {
        auto parsedIndex = std::size_t{0};
        auto value = 0;
        try {
            value = std::stoi(token, &parsedIndex);
        } catch (const std::exception&) {
            throw std::invalid_argument("cannot parse integer '" + token + "' in '" + text + "'");
        }
        if (parsedIndex != token.size()) {
            throw std::invalid_argument("cannot parse integer '" + token + "' in '" + text + "'");
        }
        if (value < 1) {
            throw std::invalid_argument("bin count '" + token + "' in '" + text + "' must be positive");
        }
        values.push_back(value);
    }
    if (values.size() != expectedCount) {
        throw std::invalid_argument(
            "'" + text + "' must contain exactly " + std::to_string(expectedCount) + " positive integers");
    }
    return values;
}

auto ParseThreeDimensionalBinning(const std::string& text) -> ThreeDimensionalBinning {
    const auto values = ParseIntegerList(text, 3);
    return ThreeDimensionalBinning{values[0], values[1], values[2]};
}

auto ParseTwoDimensionalBinning(const std::string& text) -> TwoDimensionalBinning {
    const auto values = ParseIntegerList(text, 2);
    return TwoDimensionalBinning{values[0], values[1]};
}

auto ParseOneDimensionalBinning(const std::string& text) -> int {
    const auto values = ParseIntegerList(text, 1);
    return values[0];
}

auto PrintUsage(const char* programName) -> void {
    std::cout
        << "usage: " << programName << " [inputFile] [options]\n"
        << "read a ROOT file written by sss and, for every RNTuple in it, build a 3D histogram of the\n"
        << "energy deposition positions weighted by the deposited energy, the 1D histogram of the energy\n"
        << "deposition along z, and the three 2D projections onto the xy, xz, and yz planes; all histograms\n"
        << "are stored into an output ROOT file\n"
        << "\n"
        << "positional arguments:\n"
        << "  <inputFile>              input ROOT file written by sss (one RNTuple per source category);\n"
        << "                           default: sss_output.root (the default output file of sss)\n"
        << "\n"
        << "optional options:\n"
        << "  -o, --output <file>      output ROOT file (default: <inputFile> with '_vis' inserted before\n"
        << "                           the file extension)\n"
        << "  -3, --bins-3d <nx>:<ny>:<nz>\n"
        << "                           bin counts of the 3D histogram (default: 30:30:30)\n"
        << "  -2, --bins-2d <nx>:<ny>  bin counts of the 2D projection histograms (default: 100:100)\n"
        << "  -1, --bins-1d <count>    bin count of the 1D z histogram (default: 300)\n"
        << "  -x, --x-range <min>:<max>\n"
        << "                           manual x-axis range used by every histogram (default: narrowest\n"
        << "                           interval centred at 0 containing <xy-fraction> of the deposited\n"
        << "                           energy)\n"
        << "  -y, --y-range <min>:<max>\n"
        << "                           manual y-axis range (default: as for the x axis)\n"
        << "  -z, --z-range <min>:<max>\n"
        << "                           manual z-axis range (default: <z-expand> times the interval that\n"
        << "                           contains energy deposition, centered at z = 0)\n"
        << "  -s, --xy-fraction <value>\n"
        << "                           fraction of the deposited energy kept by the automatic xy range\n"
        << "                           (default: 0.9)\n"
        << "  -e, --z-expand <factor>  expansion factor of the automatic z range (default: 1.2)\n"
        << "  -j, --threads <count>    implicit-multithreading worker count (default: all CPU cores)\n"
        << "  -f, --force              overwrite the output file if it already exists\n"
        << "  -h, --help               print this message" << '\n';
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
        } else if (argument == "-o" or argument == "--output") {
            config.mOutputFileName = nextValue(i, argument);
        } else if (argument == "-3" or argument == "--bins-3d") {
            config.mThreeDimensionalBinning = ParseThreeDimensionalBinning(nextValue(i, argument));
        } else if (argument == "-2" or argument == "--bins-2d") {
            config.mTwoDimensionalBinning = ParseTwoDimensionalBinning(nextValue(i, argument));
        } else if (argument == "-1" or argument == "--bins-1d") {
            config.mOneDimensionalBins = ParseOneDimensionalBinning(nextValue(i, argument));
        } else if (argument == "-x" or argument == "--x-range") {
            config.mXRange = ParseRange(nextValue(i, argument));
        } else if (argument == "-y" or argument == "--y-range") {
            config.mYRange = ParseRange(nextValue(i, argument));
        } else if (argument == "-z" or argument == "--z-range") {
            config.mZRange = ParseRange(nextValue(i, argument));
        } else if (argument == "-s" or argument == "--xy-fraction") {
            config.mXyEnergyFraction = ParseDouble(nextValue(i, argument));
            if (not std::isfinite(config.mXyEnergyFraction) or config.mXyEnergyFraction <= 0.0
                or config.mXyEnergyFraction > 1.0) {
                throw std::invalid_argument("--xy-fraction must be in (0, 1]");
            }
        } else if (argument == "-e" or argument == "--z-expand") {
            config.mZExpandFactor = ParseDouble(nextValue(i, argument));
            if (not std::isfinite(config.mZExpandFactor) or config.mZExpandFactor <= 0.0) {
                throw std::invalid_argument("--z-expand must be a positive number");
            }
        } else if (argument == "-j" or argument == "--threads") {
            config.mThreads = ParseInteger(nextValue(i, argument));
            if (config.mThreads < 1) {
                throw std::invalid_argument("--threads must be at least 1");
            }
        } else if (argument == "-f" or argument == "--force") {
            config.mForce = true;
        } else if (not argument.empty() and argument.front() == '-') {
            throw std::invalid_argument("unknown option '" + argument + "'");
        } else {
            positional.push_back(argument);
        }
    }
    if (positional.size() > 1) {
        throw std::invalid_argument("at most one input file may be given");
    }
    if (not positional.empty()) {
        config.mInputFileName = positional.front();
    }
    return config;
}

auto CopyPositions(const PositionVector& values) -> PositionVector {
    return values;
}

auto ConcatenatePositions(const PositionVector& first, const PositionVector& second) -> PositionVector {
    auto result = PositionVector{};
    result.reserve(first.size() + second.size());
    for (const auto& value : first) {
        result.push_back(value);
    }
    for (const auto& value : second) {
        result.push_back(value);
    }
    return result;
}

// Chain Define() calls that concatenate the per-layer deposition columns into a
// single column named resultName. RDataFrame forbids defining a column with a
// name that already exists, so every intermediate concatenation gets a unique
// intermediate name and the final result name is attached as an alias.
auto CombineColumns(ROOT::RDF::RNode node, const std::vector<std::string>& columns,
                    const std::string& resultName) -> ROOT::RDF::RNode {
    if (columns.empty()) {
        throw std::invalid_argument("cannot combine an empty column list for '" + resultName + "'");
    }
    auto current = node.Define(resultName + "_c0", CopyPositions, {columns.front()});
    for (auto index{1}; index < static_cast<int>(columns.size()); ++index) {
        const auto previousName = resultName + "_c" + std::to_string(index - 1);
        const auto currentName = resultName + "_c" + std::to_string(index);
        current = current.Define(currentName, ConcatenatePositions,
                                 {previousName, columns[static_cast<std::size_t>(index)]});
    }
    const auto finalIntermediateName = resultName + "_c" + std::to_string(columns.size() - 1);
    return current.Alias(resultName, finalIntermediateName);
}

// Scalar statistics helper. The deposition count is accumulated per entry into a
// scalar column instead of materializing a full-size intermediate vector.
auto CountDepositions(const PositionVector& values) -> double {
    return static_cast<double>(values.size());
}

// Add the scalar column holding the number of energy depositions per entry.
auto AddStatisticsColumns(ROOT::RDF::RNode node) -> ROOT::RDF::RNode {
    return node.Define(depositionCountColumnName, CountDepositions, {weightColumnName});
}

// Book the lazy statistics actions and trigger the first event loop by reading
// their results. The energy deposition weight is always positive, so the plain
// minima and maxima of the combined positions bound the data extent, which is
// used as the binning range of the fine xy histograms.
auto ComputeDepositionStatistics(ROOT::RDF::RNode& node) -> DepositionStatistics {
    auto weightSumResult = node.Sum<PositionVector>(weightColumnName);
    auto xMinimumResult = node.Min<PositionVector>(xPositionColumnName);
    auto xMaximumResult = node.Max<PositionVector>(xPositionColumnName);
    auto yMinimumResult = node.Min<PositionVector>(yPositionColumnName);
    auto yMaximumResult = node.Max<PositionVector>(yPositionColumnName);
    auto zMinimumResult = node.Min<PositionVector>(zPositionColumnName);
    auto zMaximumResult = node.Max<PositionVector>(zPositionColumnName);
    auto depositionCountResult = node.Sum<double>(depositionCountColumnName);

    auto statistics = DepositionStatistics{};
    statistics.mWeightSum = static_cast<double>(weightSumResult.GetValue());
    statistics.mDepositionCount = static_cast<long long>(depositionCountResult.GetValue());
    statistics.mXMinimum = static_cast<double>(xMinimumResult.GetValue());
    statistics.mXMaximum = static_cast<double>(xMaximumResult.GetValue());
    statistics.mYMinimum = static_cast<double>(yMinimumResult.GetValue());
    statistics.mYMaximum = static_cast<double>(yMaximumResult.GetValue());
    statistics.mZMinimum = static_cast<double>(zMinimumResult.GetValue());
    statistics.mZMaximum = static_cast<double>(zMaximumResult.GetValue());
    return statistics;
}

// Widen a degenerate range (identical bounds) by a tiny symmetric padding so that
// the histograms always have a finite positive span.
auto EnsureFiniteSpan(const Range& range) -> Range {
    auto result = range;
    if (result.mMaximum <= result.mMinimum) {
        const auto magnitude = std::max({1.0, std::abs(result.mMinimum), std::abs(result.mMaximum)});
        const auto padding = 1.0e-6 * magnitude;
        result.mMinimum -= padding;
        result.mMaximum += padding;
    }
    return result;
}

// Return the narrowest bin interval of the given fine histogram that is centred
// at zero and contains the requested fraction of the histogram's total deposited
// energy. Long distribution tails therefore only extend the range by the energy
// they actually carry.
auto SymmetricEnergyInterval(TH1D& histogram, double fraction) -> Range {
    const auto binCount = histogram.GetNbinsX();
    auto prefix = std::vector<double>(static_cast<std::size_t>(binCount) + 1, 0.0);
    for (auto bin{1}; bin <= binCount; ++bin) {
        prefix[static_cast<std::size_t>(bin)] =
            prefix[static_cast<std::size_t>(bin) - 1] + histogram.GetBinContent(bin);
    }
    const auto total = prefix[static_cast<std::size_t>(binCount)];
    if (total <= 0.0) {
        return Range{histogram.GetXaxis()->GetXmin(), histogram.GetXaxis()->GetXmax()};
    }
    const auto target = total * fraction;
    auto bestUpperBin = binCount;
    for (auto upper{1}; upper <= binCount; ++upper) {
        // The symmetric interval [-L, +L] with L at the up edge of bin `upper`
        // covers bins binCount - upper + 1 .. upper.
        const auto lower = binCount - upper;
        if (lower >= upper) {
            continue;
        }
        const auto contained =
            prefix[static_cast<std::size_t>(upper)] - prefix[static_cast<std::size_t>(lower)];
        if (contained >= target) {
            bestUpperBin = upper;
            break;
        }
    }
    const auto halfWidth = histogram.GetXaxis()->GetBinUpEdge(bestUpperBin);
    return Range{-halfWidth, halfWidth};
}

// Determine the histogram ranges. For the xy axes the automatic range is the
// narrowest interval centred at zero that contains the configured fraction of the
// deposited energy; it is derived from temporary fine 1D histograms, which are
// booked first and produced in one event loop when their results are read below.
// The z range keeps the deposition-interval rule centred at z = 0.
auto DetermineHistogramRanges(ROOT::RDF::RNode& node, const DepositionStatistics& statistics,
                              const Config& config) -> HistogramRanges {
    auto ranges = HistogramRanges{};
    auto xFine = std::optional<ROOT::RDF::RResultPtr<TH1D>>{};
    auto yFine = std::optional<ROOT::RDF::RResultPtr<TH1D>>{};
    if (not config.mXRange.has_value()) {
        // Bin the fine histogram symmetrically around zero over the largest
        // absolute data coordinate, and pad the range slightly beyond it so that
        // depositions exactly at the extent bounds land inside the range.
        const auto xHalfExtent =
            std::max({0.0, std::abs(statistics.mXMinimum), std::abs(statistics.mXMaximum)});
        const auto xFineRange = EnsureFiniteSpan(Range{-xHalfExtent, xHalfExtent});
        const auto xFinePadding = 0.001 * (xFineRange.mMaximum - xFineRange.mMinimum);
        const auto xFineModel = ROOT::RDF::TH1DModel{
            "dep_x_fine", "dep_x_fine", fineHistogramBins,
            xFineRange.mMinimum - xFinePadding, xFineRange.mMaximum + xFinePadding};
        xFine = node.Histo1D<PositionVector, PositionVector>(xFineModel, xPositionColumnName, weightColumnName);
    }
    if (not config.mYRange.has_value()) {
        const auto yHalfExtent =
            std::max({0.0, std::abs(statistics.mYMinimum), std::abs(statistics.mYMaximum)});
        const auto yFineRange = EnsureFiniteSpan(Range{-yHalfExtent, yHalfExtent});
        const auto yFinePadding = 0.001 * (yFineRange.mMaximum - yFineRange.mMinimum);
        const auto yFineModel = ROOT::RDF::TH1DModel{
            "dep_y_fine", "dep_y_fine", fineHistogramBins,
            yFineRange.mMinimum - yFinePadding, yFineRange.mMaximum + yFinePadding};
        yFine = node.Histo1D<PositionVector, PositionVector>(yFineModel, yPositionColumnName, weightColumnName);
    }
    if (xFine.has_value()) {
        ranges.mX = EnsureFiniteSpan(SymmetricEnergyInterval(*xFine.value(), config.mXyEnergyFraction));
    } else {
        ranges.mX = config.mXRange.value();
    }
    if (yFine.has_value()) {
        ranges.mY = EnsureFiniteSpan(SymmetricEnergyInterval(*yFine.value(), config.mXyEnergyFraction));
    } else {
        ranges.mY = config.mYRange.value();
    }
    if (config.mZRange.has_value()) {
        ranges.mZ = config.mZRange.value();
    } else {
        const auto span = statistics.mZMaximum - statistics.mZMinimum;
        const auto halfSpan = 0.5 * config.mZExpandFactor * span;
        ranges.mZ = EnsureFiniteSpan(Range{-halfSpan, halfSpan});
    }
    return ranges;
}

// Book the lazy histogram actions on the combined node; the event loop is only
// triggered later when the results are accessed, i.e. when they are written.
auto BuildDepositionHistograms(ROOT::RDF::RNode& node, const std::string& ntupleName,
                               const HistogramRanges& ranges, const Config& config) -> DepositionHistograms {
    const auto zName = ntupleName + "_dep_z";
    const auto zTitle = ntupleName + ": energy deposition along z";
    const auto zModel = ROOT::RDF::TH1DModel{
        zName.c_str(), zTitle.c_str(), config.mOneDimensionalBins, ranges.mZ.mMinimum, ranges.mZ.mMaximum};
    auto z = node.Histo1D<PositionVector, PositionVector>(zModel, zPositionColumnName, weightColumnName);

    const auto threeDimensionalName = ntupleName + "_dep3d";
    const auto threeDimensionalTitle = ntupleName + ": energy deposition in xyz";
    const auto threeDimensionalModel = ROOT::RDF::TH3DModel{
        threeDimensionalName.c_str(), threeDimensionalTitle.c_str(),
        config.mThreeDimensionalBinning.mX, ranges.mX.mMinimum, ranges.mX.mMaximum,
        config.mThreeDimensionalBinning.mY, ranges.mY.mMinimum, ranges.mY.mMaximum,
        config.mThreeDimensionalBinning.mZ, ranges.mZ.mMinimum, ranges.mZ.mMaximum};
    auto threeDimensional = node.Histo3D<PositionVector, PositionVector, PositionVector, PositionVector>(
        threeDimensionalModel, xPositionColumnName, yPositionColumnName, zPositionColumnName, weightColumnName);

    const auto xyName = ntupleName + "_dep_xy";
    const auto xyTitle = ntupleName + ": energy deposition in xy";
    const auto xyModel = ROOT::RDF::TH2DModel{
        xyName.c_str(), xyTitle.c_str(),
        config.mTwoDimensionalBinning.mX, ranges.mX.mMinimum, ranges.mX.mMaximum,
        config.mTwoDimensionalBinning.mY, ranges.mY.mMinimum, ranges.mY.mMaximum};
    auto xy = node.Histo2D<PositionVector, PositionVector, PositionVector>(
        xyModel, xPositionColumnName, yPositionColumnName, weightColumnName);

    const auto xzName = ntupleName + "_dep_xz";
    const auto xzTitle = ntupleName + ": energy deposition in xz";
    const auto xzModel = ROOT::RDF::TH2DModel{
        xzName.c_str(), xzTitle.c_str(),
        config.mTwoDimensionalBinning.mX, ranges.mX.mMinimum, ranges.mX.mMaximum,
        config.mTwoDimensionalBinning.mY, ranges.mZ.mMinimum, ranges.mZ.mMaximum};
    auto xz = node.Histo2D<PositionVector, PositionVector, PositionVector>(
        xzModel, xPositionColumnName, zPositionColumnName, weightColumnName);

    const auto yzName = ntupleName + "_dep_yz";
    const auto yzTitle = ntupleName + ": energy deposition in yz";
    const auto yzModel = ROOT::RDF::TH2DModel{
        yzName.c_str(), yzTitle.c_str(),
        config.mTwoDimensionalBinning.mX, ranges.mY.mMinimum, ranges.mY.mMaximum,
        config.mTwoDimensionalBinning.mY, ranges.mZ.mMinimum, ranges.mZ.mMaximum};
    auto yz = node.Histo2D<PositionVector, PositionVector, PositionVector>(
        yzModel, yPositionColumnName, zPositionColumnName, weightColumnName);

    return DepositionHistograms{z, threeDimensional, xy, xz, yz};
}

// Accessing the histogram results triggers the histogram event loop; all five
// lazy actions booked on the combined node are produced in that single pass. The
// keys are written under each histogram's own name so the two can never drift
// apart.
auto WriteHistograms(TFile& outputFile, DepositionHistograms& histograms) -> void {
    outputFile.WriteTObject(histograms.mZ.GetPtr(), histograms.mZ->GetName(), "Overwrite");
    outputFile.WriteTObject(histograms.mThreeDimensional.GetPtr(), histograms.mThreeDimensional->GetName(), "Overwrite");
    outputFile.WriteTObject(histograms.mXy.GetPtr(), histograms.mXy->GetName(), "Overwrite");
    outputFile.WriteTObject(histograms.mXz.GetPtr(), histograms.mXz->GetName(), "Overwrite");
    outputFile.WriteTObject(histograms.mYz.GetPtr(), histograms.mYz->GetName(), "Overwrite");
}

// Extract the material layer index from a deposition column name of the form
// "x_dep_{i}"; non-deposition columns yield no value.
auto ParseLayerIndex(const std::string& columnName) -> std::optional<int> {
    constexpr auto layerColumnPrefix = std::string_view{"x_dep_"};
    if (not columnName.starts_with(layerColumnPrefix)) {
        return std::nullopt;
    }
    const auto suffix = std::string_view{columnName}.substr(layerColumnPrefix.size());
    if (suffix.empty() or not std::all_of(suffix.begin(), suffix.end(), [](unsigned char character) {
            return std::isdigit(character);
        })) {
        return std::nullopt;
    }
    return std::stoi(std::string{suffix});
}

auto DiscoverLayerCount(ROOT::RDataFrame& dataframe) -> int {
    auto maximumLayerIndex = -1;
    for (const auto& columnName : dataframe.GetColumnNames()) {
        const auto layerIndex = ParseLayerIndex(columnName);
        if (not layerIndex.has_value()) {
            continue;
        }
        if (layerIndex.value() >= maximumLayerCount) {
            throw std::runtime_error(
                "layer index " + std::to_string(layerIndex.value()) + " exceeds the supported maximum of "
                + std::to_string(maximumLayerCount));
        }
        maximumLayerIndex = std::max(maximumLayerIndex, layerIndex.value());
    }
    if (maximumLayerIndex < 0) {
        throw std::runtime_error("no energy deposition columns found in the RNTuple");
    }
    return maximumLayerIndex + 1;
}

auto BuildLayerColumnNames(const std::string& stem, int layerCount) -> std::vector<std::string> {
    auto columnNames = std::vector<std::string>{};
    columnNames.reserve(static_cast<std::size_t>(layerCount));
    for (auto layerIndex{0}; layerIndex < layerCount; ++layerIndex) {
        columnNames.push_back(stem + std::to_string(layerIndex));
    }
    return columnNames;
}

auto DiscoverRNTupleNames(const std::string& fileName) -> std::vector<std::string> {
    auto file = std::unique_ptr<TFile>{TFile::Open(fileName.c_str(), "READ")};
    if (file == nullptr or file->IsZombie()) {
        throw std::runtime_error("cannot open input file '" + fileName + "'");
    }
    auto ntupleNames = std::vector<std::string>{};
    for (const TObject* object : *file->GetListOfKeys()) {
        const TKey* key = static_cast<const TKey*>(object);
        if (std::string{key->GetClassName()} == "ROOT::RNTuple") {
            ntupleNames.emplace_back(key->GetName());
        }
    }
    std::sort(ntupleNames.begin(), ntupleNames.end());
    return ntupleNames;
}

auto ProcessNtuple(const std::string& ntupleName, const std::string& inputFileName, const Config& config,
                   TFile& outputFile) -> void {
    auto dataframe = ROOT::RDataFrame{ntupleName, inputFileName};
    const auto layerCount = DiscoverLayerCount(dataframe);
    const auto xColumns = BuildLayerColumnNames("x_dep_", layerCount);
    const auto yColumns = BuildLayerColumnNames("y_dep_", layerCount);
    const auto zColumns = BuildLayerColumnNames("z_dep_", layerCount);
    const auto weightColumns = BuildLayerColumnNames("w_dep_", layerCount);

    auto node = CombineColumns(dataframe, xColumns, xPositionColumnName);
    node = CombineColumns(std::move(node), yColumns, yPositionColumnName);
    node = CombineColumns(std::move(node), zColumns, zPositionColumnName);
    node = CombineColumns(std::move(node), weightColumns, weightColumnName);
    node = AddStatisticsColumns(std::move(node));

    const auto statistics = ComputeDepositionStatistics(node);
    if (statistics.mWeightSum <= 0.0) {
        std::cerr << "warning: RNTuple '" << ntupleName << "' contains no energy deposition; skipping" << '\n';
        return;
    }
    const auto ranges = DetermineHistogramRanges(node, statistics, config);
    auto histograms = BuildDepositionHistograms(node, ntupleName, ranges, config);
    WriteHistograms(outputFile, histograms);

    std::cout << "RNTuple '" << ntupleName << "': " << statistics.mDepositionCount
              << " depositions, total deposited energy " << statistics.mWeightSum << '\n';
    std::cout << "  x range = [" << ranges.mX.mMinimum << ", " << ranges.mX.mMaximum << "]\n";
    std::cout << "  y range = [" << ranges.mY.mMinimum << ", " << ranges.mY.mMaximum << "]\n";
    std::cout << "  z range = [" << ranges.mZ.mMinimum << ", " << ranges.mZ.mMaximum << "]\n";
    std::cout << "  wrote '" << ntupleName << "_dep_z', '" << ntupleName << "_dep3d', '" << ntupleName
              << "_dep_xy', '" << ntupleName << "_dep_xz', '" << ntupleName << "_dep_yz'" << '\n';
}

} // namespace SSSVis

auto main(int argc, char** argv) -> int try {
    const auto config = SSSVis::ParseCommandLine(argc, argv);
    if (config.mHelp) {
        SSSVis::PrintUsage(argv[0]);
        return EXIT_SUCCESS;
    }

    if (not std::filesystem::exists(config.mInputFileName)) {
        throw std::runtime_error("input file '" + config.mInputFileName + "' does not exist");
    }

    const auto outputFileName = config.mOutputFileName.empty()
                                    ? SSSVis::DefaultOutputFileName(config.mInputFileName)
                                    : config.mOutputFileName;
    if (std::filesystem::exists(outputFileName)) {
        if (not config.mForce) {
            throw std::runtime_error(
                "output file '" + outputFileName + "' already exists; pass --force to overwrite it");
        }
        if (std::filesystem::equivalent(config.mInputFileName, outputFileName)) {
            throw std::runtime_error("output file '" + outputFileName + "' is the same file as the input file");
        }
    }

    const auto ntupleNames = SSSVis::DiscoverRNTupleNames(config.mInputFileName);
    if (ntupleNames.empty()) {
        throw std::runtime_error("no RNTuple found in input file '" + config.mInputFileName + "'");
    }

    ROOT::EnableImplicitMT(static_cast<UInt_t>(config.mThreads));

    auto outputFile = std::unique_ptr<TFile>{TFile::Open(outputFileName.c_str(), "RECREATE")};
    if (outputFile == nullptr or outputFile->IsZombie()) {
        throw std::runtime_error("cannot open output file '" + outputFileName + "'");
    }

    for (const auto& ntupleName : ntupleNames) {
        SSSVis::ProcessNtuple(ntupleName, config.mInputFileName, config, *outputFile);
    }
    outputFile->Close();

    std::cout << "wrote " << ntupleNames.size() << " RNTuple visualization(s) to '" << outputFileName << "'\n";

    return EXIT_SUCCESS;
} catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    return EXIT_FAILURE;
}
