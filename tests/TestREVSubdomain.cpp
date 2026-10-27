#include "analysis/rev_analysis.h"
#include "analysis/rev_subdomain_flow.h"
#include "common/UnitTest.h"

#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <dirent.h>
#include <fstream>
#include <iterator>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>
#include <unistd.h>

namespace {
using Row = std::map<std::string, std::string>;

void require(bool condition, const std::string &message) {
    if (!condition)
        throw std::runtime_error(message);
}

std::vector<std::string> selectedRawFiles() {
    std::vector<std::string> files;
    DIR *directory = opendir(".");
    if (!directory)
        return files;
    while (const auto entry = readdir(directory)) {
        const std::string name(entry->d_name);
        if (name.compare(0, 10, "subdomain_") == 0 &&
            ((name.size() >= 4 && name.compare(name.size() - 4, 4, ".raw") == 0) ||
             (name.size() >= 8 && name.compare(name.size() - 8, 8, ".raw.tmp") == 0)))
            files.push_back(name);
    }
    closedir(directory);
    return files;
}

// Isolate generated fixtures and outputs from other tests and user files.
class TestDirectory {
public:
    TestDirectory() {
        char cwd[4096];
        require(getcwd(cwd, sizeof(cwd)) != nullptr, "Cannot read current directory");
        original = cwd;
        char path[] = "/tmp/lbpm-rev-selection-test-XXXXXX";
        char *created = mkdtemp(path);
        require(created != nullptr, "Cannot create test directory");
        temporary = created;
        require(chdir(temporary.c_str()) == 0, "Cannot enter test directory");
    }
    ~TestDirectory() {
        for (const auto &name : selectedRawFiles())
            std::remove(name.c_str());
        for (const char *name : {"input.db", "source.raw", "subdomain_selection.csv",
                                 "geometrical-rev-domain.db", "hydraulic-rev-domain.db",
                                 "geometrical-rev-domain.db.tmp", "hydraulic-rev-domain.db.tmp",
                                 "ID.00000", "stat_rev_geometric.csv", "stat_rev_hydraulic.csv"})
            std::remove(name);
        chdir(original.c_str());
        rmdir(temporary.c_str());
    }
private:
    std::string original, temporary;
};

std::vector<std::string> splitCSV(const std::string &line) {
    std::vector<std::string> fields;
    std::string field;
    bool quoted = false;
    for (size_t i = 0; i < line.size(); ++i) {
        const char c = line[i];
        if (c == '"') {
            if (quoted && i + 1 < line.size() && line[i + 1] == '"') {
                field += '"';
                ++i;
            } else {
                quoted = !quoted;
            }
        } else if (c == ',' && !quoted) {
            fields.push_back(field);
            field.clear();
        } else if (c != '\r') {
            field += c;
        }
    }
    fields.push_back(field);
    return fields;
}

std::map<std::string, Row> readSummary() {
    std::ifstream input("subdomain_selection.csv");
    require(input.good(), "Selection summary was not written");
    std::string line;
    std::getline(input, line);
    const auto header = splitCSV(line);
    require(header.size() == 45, "Summary must include both geometric and hydraulic score diagnostics");
    std::map<std::string, Row> result;
    while (std::getline(input, line)) {
        if (line.empty())
            continue;
        const auto values = splitCSV(line);
        require(values.size() == header.size(), "Summary CSV column count mismatch");
        Row row;
        for (size_t i = 0; i < header.size(); ++i)
            row[header[i]] = values[i];
        result[row.at("property")] = row;
    }
    return result;
}

std::vector<char> readBytes(const std::string &filename) {
    std::ifstream input(filename, std::ios::binary);
    require(input.good(), "Selected raw image was not written");
    return {std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>()};
}

std::string requireRawFile(const Row &row) {
    const std::string expected = "subdomain_" + row.at("property") + "_" + row.at("nx") +
                                 "x" + row.at("ny") + "x" + row.at("nz") + ".raw";
    require(row.at("raw_file") == expected,
            "Raw filename must contain the selected voxel dimensions in x/y/z order");
    return expected;
}

bool close(double actual, double expected) {
    return std::abs(actual - expected) < 1.0e-12;
}

std::string absolutePath(const std::string &filename) {
    char *resolved = realpath(filename.c_str(), nullptr);
    require(resolved != nullptr, "Cannot resolve source filename: " + filename);
    const std::string result(resolved);
    std::free(resolved);
    return result;
}

std::shared_ptr<Database> requireDomainFile(const Row &row) {
    const std::string expected_file = row.at("property") == "porosity"
        ? "geometrical-rev-domain.db" : "hydraulic-rev-domain.db";
    require(row.at("domain_file") == expected_file,
            "Selection summary did not identify its Domain database");
    Database selected(expected_file);
    require(selected.getAllKeys() == std::vector<std::string>{"Domain"},
            "Selected input database must contain only a Domain block");
    auto domain = selected.getDatabase("Domain");
    const auto verify_vector = [&](const char *key, const std::array<const char *, 3> &columns) {
        const auto values = domain->getVector<int>(key);
        require(values.size() == 3, std::string("Selected Domain.") + key + " must have three components");
        for (int d = 0; d < 3; ++d)
            require(values[d] == std::stoi(row.at(columns[d])),
                    std::string("Selected Domain.") + key + " differs from the selected source region");
    };
    verify_vector("n", {{"nx", "ny", "nz"}});
    verify_vector("N", {{"source_nx", "source_ny", "source_nz"}});
    verify_vector("offset", {{"source_x0", "source_y0", "source_z0"}});
    require(domain->getVector<int>("nproc") == std::vector<int>({1, 1, 1}),
            "Selected Domain must use the supported single-process layout");
    require(absolutePath(domain->getScalar<std::string>("Filename")) ==
                absolutePath(row.at("source_file")),
            "Selected Domain must reference the original full image");
    require(close(domain->getScalar<double>("voxel_length"), 2.5),
            "Selected Domain must preserve voxel length");
    if (domain->keyExists("L")) {
        const auto physical = domain->getVector<double>("L");
        const auto size = domain->getVector<int>("n");
        require(physical.size() == 3, "Selected Domain.L must have three components");
        for (int d = 0; d < 3; ++d)
            require(close(physical[d], size[d] * 2.5),
                    "Selected Domain.L must describe the selected physical extent");
    }
    require(domain->getScalar<int>("BC") == 0,
            "Selected Domain must reproduce periodic candidate geometry");
    for (const char *key : {"InletLayers", "OutletLayers"})
        if (domain->keyExists(key))
            require(domain->getVector<int>(key) == std::vector<int>({0, 0, 0}),
                    "Selected Domain must not reapply synthetic layers");
    if (row.at("read_type") == "signed8bit") {
        require(domain->getScalar<std::string>("ReadType") == "8bit",
                "ID input must be described using the supported raw reader");
        require(domain->getVector<int>("ReadValues") == domain->getVector<int>("WriteValues"),
                "Already-mapped ID labels must use an identity mapping");
    } else {
        require(domain->getScalar<std::string>("ReadType") == row.at("read_type"),
                "Selected Domain must preserve the original raw type");
        Database original("input.db");
        auto input_domain = original.getDatabase("Domain");
        for (const char *key : {"ReadValues", "WriteValues"})
            require(domain->getVector<int>(key) == input_domain->getVector<int>(key),
                    "Selected Domain must preserve raw label mapping");
    }
    return domain;
}

void requireDomainRoundTrip(const std::shared_ptr<Database> &selected,
                            const std::vector<signed char> &expected,
                            const Utilities::MPI &comm) {
    // Decomp writes ID.00000, so isolate the reload from the original image,
    // particularly when the original source itself was an ID.00000 file.
    auto reloaded = selected->cloneDatabase();
    const auto filename = absolutePath(selected->getScalar<std::string>("Filename"));
    reloaded->putScalar<std::string>("Filename", filename);
    TestDirectory reload_directory;
    Domain domain(reloaded, comm);
    domain.Decomp(filename);
    std::vector<signed char> interior;
    for (int z = 1; z < domain.Nz - 1; ++z)
        for (int y = 1; y < domain.Ny - 1; ++y)
            for (int x = 1; x < domain.Nx - 1; ++x)
                interior.push_back(domain.id[(z * domain.Ny + y) * domain.Nx + x]);
    require(interior == expected,
            "Reloading the emitted Domain did not reproduce the selected original-image labels");
}

void writeInput(int n, std::array<int, 3> source_size,
                std::array<int, 3> offset, const std::string &read_type,
                bool hydraulic = false, bool selection = true,
                bool synthetic_layers = false, const char *save_raw = "true",
                const char *subdomain_size = nullptr) {
    std::ofstream out("input.db");
    out << "Domain {\n"
        << "  Filename = \"source.raw\"\n"
        << "  ReadType = \"" << read_type << "\"\n"
        << "  N = " << source_size[0] << ", " << source_size[1] << ", " << source_size[2] << "\n"
        << "  n = " << n << ", " << n << ", " << n << "\n"
        << "  nproc = 1, 1, 1\n"
        << "  offset = " << offset[0] << ", " << offset[1] << ", " << offset[2] << "\n";
    if (synthetic_layers)
        out << "  InletLayers = 1, 1, 1\n  OutletLayers = 1, 1, 1\n";
    out << "  BC = 0\n  voxel_length = 2.5\n"
        << "  L = " << source_size[0] * 2.5 << ", " << source_size[1] * 2.5 << ", "
        << source_size[2] * 2.5 << "\n"
        << "  ReadValues = 17, 29\n  WriteValues = 0, 1\n}\n"
        << "REV {\n  Geometrical = true\n  Hydraulic = " << (hydraulic ? "true" : "false") << "\n"
        << "  SubdomainSelect = " << (selection ? "true" : "false") << "\n";
    if (save_raw)
        out << "  save_subdomain_raw = " << save_raw << "\n";
    if (subdomain_size)
        out << "  subdomain_size = " << subdomain_size << "\n";
    out << "}\n"
        << "MRT {\n  tau = 1.0\n  F = 0, 0, 0.00001\n"
        << "  timestepMax = 4000\n  analysis_interval = 20\n  tolerance = 0.00000001\n}\n";
    require(out.good(), "Cannot write test input database");
}

template<class T> void writeRaw(const std::vector<T> &data) {
    std::ofstream out("source.raw", std::ios::binary);
    out.write(reinterpret_cast<const char *>(data.data()), data.size() * sizeof(T));
    require(out.good(), "Cannot write test raw image");
}

struct Model {
    int Nx, Ny, Nz;
    std::shared_ptr<Domain> Dm;
    DoubleArray Distance, Velocity_x, Velocity_y, Velocity_z;
    double tau = 1.0, Fx = 0.0, Fy = 0.0, Fz = 1.0e-5;

    Model(int n, const Utilities::MPI &comm) : Model(std::array<int, 3>{{n, n, n}}, comm) {}
    Model(std::array<int, 3> size, const Utilities::MPI &comm)
        : Nx(size[0] + 2), Ny(size[1] + 2), Nz(size[2] + 2) {
        auto input = std::make_shared<Database>("input.db");
        Dm = std::make_shared<Domain>(input->getDatabase("Domain"), comm);
        Distance.resize(Nx, Ny, Nz);
        Distance.fill(-1.0);
        Velocity_x.resize(Nx, Ny, Nz);
        Velocity_y.resize(Nx, Ny, Nz);
        Velocity_z.resize(Nx, Ny, Nz);
        Velocity_x.fill(0.0);
        Velocity_y.fill(0.0);
        Velocity_z.fill(0.0);
    }
};

void setREV(ScaLBL_REV::StatisticalREVResult &result, int size, double mean,
            double surface_mean = 0.0, double tortuosity_mean = 1.0) {
    result.size = size;
    result.dimensions = {{size, size, size}};
    result.mean = mean;
    result.cv = 0.0;
    result.surface_mean = surface_mean;
    result.tortuosity_mean = tortuosity_mean;
}

void requireWeightedDiagnostics(const Row &row, double mean_surface,
                                double surface, double score) {
    require(close(std::stod(row.at("statistical_mean_surface")), mean_surface),
            "Summary lost the statistical surface mean at the porosity REV size");
    require(close(std::stod(row.at("selected_surface")), surface),
            "Summary reports incorrect selected specific surface area");
    const double surface_error = mean_surface != 0.0
        ? std::abs(surface - mean_surface) / std::abs(mean_surface)
        : (surface == 0.0 ? 0.0 : std::numeric_limits<double>::infinity());
    const double reported_error = std::stod(row.at("relative_deviation_surface"));
    require(reported_error == surface_error || close(reported_error, surface_error),
            "Summary reports incorrect relative surface deviation");
    const double reported_score = std::stod(row.at("selection_score"));
    require(reported_score == score || close(reported_score, score),
            "Selection score must be 60% porosity error plus 40% surface error");
    require(close(std::stod(row.at("porosity_weight")), 0.6) &&
                close(std::stod(row.at("surface_weight")), 0.4),
            "Summary does not record the requested 60/40 geometric weights");
    for (const char *key : {"statistical_mean_tortuosity", "selected_tortuosity",
                            "relative_deviation_tortuosity", "permeability_weight", "tortuosity_weight"})
        require(row.at(key).empty(), "Geometric scoring must not report hydraulic diagnostics");
}

void requireHydraulicDiagnostics(const Row &row, double mean_tortuosity,
                                 double selected_tortuosity) {
    require(close(std::stod(row.at("statistical_mean_tortuosity")), mean_tortuosity),
            "Summary lost the statistical tortuosity mean at the permeability REV size");
    require(close(std::stod(row.at("selected_tortuosity")), selected_tortuosity),
            "Summary reports incorrect hydraulic tortuosity");
    const double expected_error = std::abs(selected_tortuosity - mean_tortuosity) /
                                  std::abs(mean_tortuosity);
    require(close(std::stod(row.at("relative_deviation_tortuosity")), expected_error),
            "Summary reports incorrect relative tortuosity deviation");
    require(close(std::stod(row.at("selection_score")),
                  0.6 * std::stod(row.at("relative_deviation")) + 0.4 * expected_error),
            "Hydraulic selection score must be 60% permeability and 40% tortuosity error");
    require(close(std::stod(row.at("permeability_weight")), 0.6) &&
                close(std::stod(row.at("tortuosity_weight")), 0.4),
            "Summary must record the requested 60/40 hydraulic weights");
    for (const char *key : {"statistical_mean_surface", "selected_surface", "relative_deviation_surface",
                            "porosity_weight", "surface_weight"})
        require(row.at(key).empty(), "Hydraulic scoring must not report geometric diagnostics");
}

struct GeometryMeasurement {
    int pores = 0;
    int faces = 0;
    std::array<int, 3> face_axes{{0, 0, 0}};
    double porosity = 0.0;
    double surface = 0.0;
};

GeometryMeasurement measureWindow(const std::vector<uint8_t> &mask,
                                   const std::array<int, 3> &shape,
                                   const std::array<int, 3> &origin,
                                   const std::array<int, 3> &size, double h) {
    GeometryMeasurement measured;
    const auto pore = [&](int x, int y, int z) {
        return mask[x + shape[0] * (y + shape[1] * z)] != 0;
    };
    for (int z = origin[2]; z < origin[2] + size[2]; ++z)
        for (int y = origin[1]; y < origin[1] + size[1]; ++y)
            for (int x = origin[0]; x < origin[0] + size[0]; ++x) {
                const bool current = pore(x, y, z);
                measured.pores += current;
                if (x + 1 < origin[0] + size[0] && current != pore(x + 1, y, z))
                    ++measured.face_axes[0];
                if (y + 1 < origin[1] + size[1] && current != pore(x, y + 1, z))
                    ++measured.face_axes[1];
                if (z + 1 < origin[2] + size[2] && current != pore(x, y, z + 1))
                    ++measured.face_axes[2];
            }
    for (int count : measured.face_axes)
        measured.faces += count;
    const double volume = double(size[0]) * size[1] * size[2];
    measured.porosity = measured.pores / volume;
    measured.surface = measured.faces / (volume * h);
    return measured;
}

void testWeightedGeometrySelection(const Utilities::MPI &comm) {
    TestDirectory directory;
    writeInput(4, {{4, 4, 4}}, {{0, 0, 0}}, "8bit");
    const std::string fixture = "1111001001000110001111011011010100000000001111010011000111101011";
    std::vector<uint8_t> mask, source;
    Model model(4, comm);
    for (size_t i = 0; i < fixture.size(); ++i) {
        const bool pore = fixture[i] == '1';
        mask.push_back(pore ? 1 : 0);
        source.push_back(pore ? 29 : 17);
        const int x = i % 4, y = (i / 4) % 4, z = i / 16;
        model.Distance(x + 1, y + 1, z + 1) = pore ? 1.0 : -1.0;
    }
    writeRaw(source);
    const double mean_porosity = 0.5;
    const double mean_surface = 5.0 / (8.0 * 2.5);
    double best_score = std::numeric_limits<double>::infinity();
    std::array<int, 3> best_origin{{-1, -1, -1}};
    GeometryMeasurement best_measurement;
    int candidates = 0;
    for (int z = 0; z <= 2; ++z)
        for (int y = 0; y <= 2; ++y)
            for (int x = 0; x <= 2; ++x) {
                const auto measured = measureWindow(mask, {{4, 4, 4}}, {{x, y, z}}, {{2, 2, 2}}, 2.5);
                const double porosity_error = std::abs(measured.porosity - mean_porosity) / mean_porosity;
                const double surface_error = std::abs(measured.surface - mean_surface) / mean_surface;
                const double score = 0.6 * porosity_error + 0.4 * surface_error;
                ++candidates;
                if (score < best_score) {
                    best_score = score;
                    best_origin = {{x, y, z}};
                    best_measurement = measured;
                }
            }
    require(best_origin == std::array<int, 3>{{2, 0, 1}} &&
                best_measurement.pores == 3 && best_measurement.faces == 5 && close(best_score, 0.15),
            "Weighted regression fixture no longer has its expected compromise optimum");
    const auto porosity_only = measureWindow(mask, {{4, 4, 4}}, {{0, 0, 0}}, {{2, 2, 2}}, 2.5);
    require(close(porosity_only.porosity, mean_porosity) &&
                porosity_only.faces == 8 && best_measurement.porosity != mean_porosity,
            "Fixture must distinguish weighted selection from a porosity-only exact match");
    ScaLBL_REV rev(0, 1, comm);
    rev.ReadParams("input.db", false);
    setREV(rev.stat_poro_rev, 2, mean_porosity, mean_surface);
    rev.sub_domain_select(model);
    const Row row = readSummary().at("porosity");
    const std::array<int, 3> actual_origin{{std::stoi(row.at("x0")), std::stoi(row.at("y0")),
                                          std::stoi(row.at("z0"))}};
    require(actual_origin == best_origin,
            "Selection did not minimize the independently enumerated 60/40 geometry score");
    require(std::stoull(row.at("candidate_count")) == unsigned(candidates),
            "Weighted selection must still evaluate all valid windows");
    require(close(std::stod(row.at("selected_value")), best_measurement.porosity) &&
                close(std::stod(row.at("relative_deviation")), 0.25),
            "Geometric CSV must retain selected porosity and its separate relative deviation");
    requireWeightedDiagnostics(row, mean_surface, best_measurement.surface, best_score);
    std::vector<char> expected;
    for (int z = best_origin[2]; z < best_origin[2] + 2; ++z)
        for (int y = best_origin[1]; y < best_origin[1] + 2; ++y)
            for (int x = best_origin[0]; x < best_origin[0] + 2; ++x)
                expected.push_back(source[x + 4 * (y + 4 * z)]);
    require(readBytes(requireRawFile(row)) == expected,
            "Weighted-selected raw image does not contain the winning compromise geometry");
    requireDomainFile(row);
}

void testLastWindowAndOffset(const Utilities::MPI &comm) {
    TestDirectory directory;
    writeInput(4, {{6, 7, 8}}, {{1, 2, 3}}, "8bit");
    std::vector<uint8_t> source(6 * 7 * 8, 17);
    source[(6 * 7 * 6) + (6 * 5) + 4] = 29;
    writeRaw(source);
    Model model(4, comm);
    model.Distance(4, 4, 4) = 1.0;
    ScaLBL_REV rev(0, 1, comm);
    rev.ReadParams("input.db", false);
    setREV(rev.stat_poro_rev, 2, 0.125, 3.0 / (8.0 * 2.5));
    rev.sub_domain_select(model);
    const Row row = readSummary().at("porosity");
    requireWeightedDiagnostics(row, 3.0 / (8.0 * 2.5), 3.0 / (8.0 * 2.5), 0.0);
    for (const char *key : {"x0", "y0", "z0"})
        require(std::stoi(row.at(key)) == 2, "Final valid window was not considered");
    for (const char *key : {"x1", "y1", "z1"})
        require(std::stoi(row.at(key)) == 4, "Ending coordinate is not exclusive");
    require(std::stoull(row.at("candidate_count")) == 27, "Expected every valid window to be scanned");
    require(close(std::stod(row.at("selected_value")), 0.125), "Incorrect selected porosity");
    require(close(std::stod(row.at("relative_deviation")), 0.0), "Exact match should have zero error");
    require(close(std::stod(row.at("rev_size_um")), 5.0), "Incorrect physical REV size");
    require(std::stoi(row.at("source_x0")) == 3 && std::stoi(row.at("source_y0")) == 4 &&
                std::stoi(row.at("source_z0")) == 5, "Original-image offset was lost");
    require(std::stoi(row.at("source_x1")) == 5 && std::stoi(row.at("source_y1")) == 6 &&
                std::stoi(row.at("source_z1")) == 7, "Incorrect original-image exclusive bounds");
    const auto actual = readBytes(requireRawFile(row));
    std::vector<char> expected(8, 17);
    expected.back() = 29;
    require(actual == expected, "8-bit raw export changed labels, order, or offset");
    auto selected_domain = requireDomainFile(row);
    std::vector<signed char> mapped(8, 0);
    mapped.back() = 1;
    requireDomainRoundTrip(selected_domain, mapped, comm);
}

void testSixteenBitLabelsAndTies(const Utilities::MPI &comm) {
    TestDirectory directory;
    writeInput(3, {{3, 3, 3}}, {{0, 0, 0}}, "16bit");
    std::vector<uint16_t> source(27);
    Model model(3, comm);
    for (int z = 0; z < 3; ++z)
        for (int y = 0; y < 3; ++y)
            for (int x = 0; x < 3; ++x) {
                const bool pore = (x + y + z) % 2 != 0;
                source[x + 3 * (y + 3 * z)] = pore ? 60000 : 300;
                model.Distance(x + 1, y + 1, z + 1) = pore ? 1.0 : -1.0;
            }
    writeRaw(source);
    ScaLBL_REV rev(0, 1, comm);
    rev.ReadParams("input.db", false);
    setREV(rev.stat_poro_rev, 2, 0.5, 12.0 / (8.0 * 2.5));
    rev.sub_domain_select(model);
    const Row row = readSummary().at("porosity");
    requireWeightedDiagnostics(row, 12.0 / (8.0 * 2.5), 12.0 / (8.0 * 2.5), 0.0);
    require(std::stoi(row.at("x0")) == 0 && std::stoi(row.at("y0")) == 0 &&
                std::stoi(row.at("z0")) == 0, "Ties should retain the first window");
    require(std::stoull(row.at("candidate_count")) == 8, "An exact match must not stop the exhaustive scan");
    require(std::stoi(row.at("bytes_per_voxel")) == 2, "16-bit output type was lost");
    std::vector<uint16_t> expected;
    for (int z = 0; z < 2; ++z)
        for (int y = 0; y < 2; ++y)
            for (int x = 0; x < 2; ++x)
                expected.push_back(source[x + 3 * (y + 3 * z)]);
    const char *bytes = reinterpret_cast<const char *>(expected.data());
    require(readBytes(requireRawFile(row)) == std::vector<char>(bytes, bytes + 16),
            "16-bit raw export did not preserve original bytes and labels");
    requireDomainFile(row);
}

void testIndependentSizesAndZeroMeans(const Utilities::MPI &comm) {
    TestDirectory directory;
    writeInput(3, {{3, 3, 3}}, {{0, 0, 0}}, "8bit", true);
    writeRaw(std::vector<uint8_t>(27, 17));
    Model model(3, comm);
    ScaLBL_REV rev(0, 1, comm);
    rev.ReadParams("input.db", false);
    setREV(rev.stat_poro_rev, 2, 0.0);
    setREV(rev.stat_perm_rev, 1, 0.0);
    rev.sub_domain_select(model);
    const auto rows = readSummary();
    const Row porosity = rows.at("porosity"), permeability = rows.at("permeability");
    require(std::stoi(porosity.at("rev_size_vx")) == 2 &&
                std::stoi(permeability.at("rev_size_vx")) == 1, "Property REV sizes were coupled");
    require(std::stoull(porosity.at("candidate_count")) == 8 &&
                std::stoull(permeability.at("candidate_count")) == 27,
            "Each property must scan its own complete set of windows");
    requireWeightedDiagnostics(porosity, 0.0, 0.0, 0.0);
    for (const char *key : {"statistical_mean_surface", "selected_surface", "relative_deviation_surface",
                            "selection_score", "porosity_weight", "surface_weight"})
        require(permeability.at(key).empty(), "Hydraulic rows must not imply geometric scoring");
    require(close(std::stod(porosity.at("selected_value")), 0.0) &&
                close(std::stod(porosity.at("relative_deviation")), 0.0),
            "Solid geometry should exactly match zero geometric targets");
    require(permeability.at("status") == "no_valid_tortuosity",
            "Solid candidates must not fabricate a finite hydraulic tortuosity");
    require(readBytes(requireRawFile(porosity)).size() == 8,
            "Selected geometric dimensions were not independent of the hydraulic REV");
    require(!exists("subdomain_permeability_1x1x1.raw") && !exists("hydraulic-rev-domain.db"),
            "An undefined hydraulic score must not produce a selected raw image or Domain database");
    requireDomainFile(porosity);
}

void testNoConvergence(const Utilities::MPI &comm) {
    TestDirectory directory;
    writeInput(3, {{3, 3, 3}}, {{0, 0, 0}}, "8bit");
    writeRaw(std::vector<uint8_t>(27, 17));
    Model model(3, comm);
    ScaLBL_REV rev(0, 1, comm);
    rev.ReadParams("input.db", false);
    rev.sub_domain_select(model);
    require(readSummary().at("porosity").at("status") == "not_converged",
            "Missing statistical convergence should be explicit");
    require(selectedRawFiles().empty(), "An unconverged REV must not produce selected raw images");
    require(!exists("geometrical-rev-domain.db") && !exists("hydraulic-rev-domain.db"),
            "An unconverged property must not produce a Domain database");
    require(readSummary().at("porosity").at("domain_file").empty(),
            "An unconverged property must not advertise a Domain database");
}

void testFullDomainWindow(const Utilities::MPI &comm) {
    TestDirectory directory;
    writeInput(3, {{3, 3, 3}}, {{0, 0, 0}}, "8bit");
    std::vector<uint8_t> source(27, 17);
    source.back() = 29;
    writeRaw(source);
    Model model(3, comm);
    model.Distance(3, 3, 3) = 1.0;
    ScaLBL_REV rev(0, 1, comm);
    rev.ReadParams("input.db", false);
    setREV(rev.stat_poro_rev, 3, 1.0 / 27.0, 3.0 / (27.0 * 2.5));
    rev.sub_domain_select(model);
    const Row row = readSummary().at("porosity");
    require(std::stoull(row.at("candidate_count")) == 1, "A full-domain REV has exactly one candidate");
    require(std::stoi(row.at("x0")) == 0 && std::stoi(row.at("x1")) == 3,
            "Full-domain REV bounds are incorrect");
    const char *bytes = reinterpret_cast<const char *>(source.data());
    require(readBytes(requireRawFile(row)) == std::vector<char>(bytes, bytes + source.size()),
            "A full-domain selection should reproduce the original image");
}

void testUnmatchedZeroMean(const Utilities::MPI &comm) {
    TestDirectory directory;
    writeInput(3, {{3, 3, 3}}, {{0, 0, 0}}, "8bit");
    std::vector<uint8_t> source(27, 29);
    source.back() = 17;
    writeRaw(source);
    Model model(3, comm);
    model.Distance.fill(1.0);
    model.Distance(3, 3, 3) = -1.0;
    ScaLBL_REV rev(0, 1, comm);
    rev.ReadParams("input.db", false);
    setREV(rev.stat_poro_rev, 2, 0.0);
    rev.sub_domain_select(model);
    const Row row = readSummary().at("porosity");
    require(std::isinf(std::stod(row.at("relative_deviation"))),
            "Nonzero measurements relative to a zero mean should have infinite error");
    require(close(std::stod(row.at("selected_value")), 1.0), "A zero mean must not discard every candidate");
    require(std::stoi(row.at("x0")) == 0 && std::stoi(row.at("y0")) == 0 &&
                std::stoi(row.at("z0")) == 0,
            "Infinite-error ties should select the first candidate even when a later absolute error is smaller");
}

void testDisabled(const Utilities::MPI &comm) {
    TestDirectory directory;
    writeInput(3, {{3, 3, 3}}, {{0, 0, 0}}, "8bit", false, false);
    Model model(3, comm);
    ScaLBL_REV rev(0, 1, comm);
    rev.ReadParams("input.db", false);
    setREV(rev.stat_poro_rev, 2, 0.5, 12.0 / (8.0 * 2.5));
    // Deliberately omit the source image: disabled selection must not access it.
    rev.sub_domain_select(model);
    require(!exists("subdomain_selection.csv") && selectedRawFiles().empty(),
            "Disabled selection must not create outputs");
    require(!exists("geometrical-rev-domain.db") && !exists("hydraulic-rev-domain.db"),
            "Disabled selection must not create Domain databases");
}

void testOptionalRawExports(const Utilities::MPI &comm) {
    for (const char *setting : std::array<const char *, 3>{{nullptr, "false", "true"}}) {
        const bool save_raw = setting && std::string(setting) == "true";
        TestDirectory directory;
        writeInput(4, {{4, 4, 4}}, {{0, 0, 0}}, "8bit", true, true, false, setting);
        Model model(4, comm);
        std::vector<uint8_t> source(64, 17);
        Array<signed char> reference(4, 4, 4);
        reference.fill(0);
        for (int z = 0; z < 4; ++z)
            for (int y = 1; y <= 2; ++y)
                for (int x = 1; x <= 2; ++x) {
                    source[x + 4 * (y + 4 * z)] = 29;
                    reference(x, y, z) = 1;
                    model.Distance(x + 1, y + 1, z + 1) = 1.0;
                }
        writeRaw(source);
        const auto original = readBytes("source.raw");
        const auto hydraulic = rev_selection::simulateHydraulicProperties(reference, 2.5, 1.0,
                        {{0.0, 0.0, 1.0e-5}}, 4000, 20, 1.0e-8, comm);
        ScaLBL_REV rev(0, 1, comm);
        rev.ReadParams("input.db", false);
        setREV(rev.stat_poro_rev, 4, 0.25, 32.0 / (64.0 * 2.5));
        setREV(rev.stat_perm_rev, 4, hydraulic.permeability, 0.0, hydraulic.tortuosity);
        rev.sub_domain_select(model);
        const auto rows = readSummary();
        for (const char *property : {"porosity", "permeability"}) {
            const Row &row = rows.at(property);
            const std::string raw_file = std::string("subdomain_") + property + "_4x4x4.raw";
            require(row.at("status") == "selected" && std::stoull(row.at("candidate_count")) == 1,
                    "Optional raw export must preserve selection for both properties");
            requireDomainFile(row);
            if (save_raw) {
                require(row.at("raw_file") == raw_file && readBytes(raw_file) == original,
                        "Explicitly enabled raw export must preserve selected source bytes");
                require(!exists("subdomain_porosity.raw") && !exists("subdomain_permeability.raw"),
                        "Raw export must use dimensioned filenames instead of legacy names");
            } else {
                require(row.at("raw_file").empty(),
                        "Disabled raw export must not advertise a raw output in the summary");
                require(!exists(raw_file) && !exists(raw_file + ".tmp"),
                        "Disabled raw export must not create raw images or staging files");
            }
        }
        if (!save_raw) {
            std::vector<signed char> mapped;
            for (auto label : source)
                mapped.push_back(label == 29 ? 1 : 0);
            requireDomainRoundTrip(requireDomainFile(rows.at("porosity")), mapped, comm);
            // Files that raw export will not touch may alias the source. This
            // also verifies that existing raw images and staging files survive.
            for (const char *name : {"subdomain_porosity_4x4x4.raw", "subdomain_permeability_4x4x4.raw",
                                     "subdomain_porosity_4x4x4.raw.tmp", "subdomain_permeability_4x4x4.raw.tmp"})
                require(link("source.raw", name) == 0, "Cannot create unused raw output aliases");
            rev.sub_domain_select(model);
            for (const char *name : {"subdomain_porosity_4x4x4.raw", "subdomain_permeability_4x4x4.raw",
                                     "subdomain_porosity_4x4x4.raw.tmp", "subdomain_permeability_4x4x4.raw.tmp"})
                require(readBytes(name) == original,
                        "Disabled raw export must preserve existing raw images and staging files");
            require(readBytes("source.raw") == original,
                    "Selection without raw export must preserve the aliased original image");
            for (const auto &entry : readSummary()) {
                require(entry.second.at("raw_file").empty(),
                        "Existing unused raw files must not appear as newly exported outputs");
                requireDomainFile(entry.second);
            }
        }
    }
}

void testStatisticalRetention(const Utilities::MPI &comm) {
    TestDirectory directory;
    writeInput(4, {{4, 4, 4}}, {{0, 0, 0}}, "8bit", true, false);
    Model model(4, comm);
    model.Distance.fill(1.0);
    // Every two-voxel x interval has the same mean numerical velocity. A
    // one-voxel interval has dispersion, so hydraulic convergence occurs later
    // than porosity convergence and has a known mean at its retained size.
    for (int z = 1; z <= 4; ++z)
        for (int y = 1; y <= 4; ++y)
            for (int x = 1; x <= 4; ++x) {
                model.Velocity_z(x, y, z) = (x % 2 == 0 ? 3.0 : 1.0) * model.Fz;
                model.Velocity_x(x, y, z) = (x % 2 == 0 ? 4.0 : 0.0) * model.Fz;
            }
    ScaLBL_REV rev(0, 1, comm);
    rev.ReadParams("input.db", false);
    rev.average_pore_size = 0.5;
    rev.StatRevAnalysis(model);
    require(rev.stat_poro_rev.size == 1 && rev.stat_perm_rev.size == 2,
            "Statistical properties must retain their own first converged size");
    require(rev.stat_poro_rev.dimensions == std::array<int, 3>{{1, 1, 1}} &&
                rev.stat_perm_rev.dimensions == std::array<int, 3>{{2, 2, 2}},
            "Statistical sampling dimensions must match cubic REV sizes");
    require(close(rev.stat_poro_rev.mean, 1.0), "Wrong porosity mean retained at the statistical REV");
    require(close(rev.stat_poro_rev.surface_mean, 0.0),
            "An all-pore statistical REV must retain zero internal surface area");
    const double expected = 1013.0 * 2.5 * 2.5 * ((model.tau - 0.5) / 3.0) * 2.0;
    require(std::abs(rev.stat_perm_rev.mean - expected) < 1.0e-9,
            "Wrong permeability mean retained at its own statistical REV size");
    // Each two-voxel interval has speeds F and 5F and z velocities F and
    // 3F, so T=6/4. Three-voxel intervals have different tortuosities.
    require(close(rev.stat_perm_rev.tortuosity_mean, 1.5),
            "Tortuosity target must be retained at the permeability REV, not at a later size");
}

void testStatisticalSurfaceRetention(const Utilities::MPI &comm) {
    TestDirectory directory;
    writeInput(4, {{4, 4, 4}}, {{0, 0, 0}}, "8bit", false, false);
    Model model(4, comm);
    for (int z = 1; z <= 4; ++z)
        for (int y = 1; y <= 4; ++y)
            for (int x = 1; x <= 4; ++x)
                model.Distance(x, y, z) = (x + y + z) % 2 ? 1.0 : -1.0;
    ScaLBL_REV rev(0, 1, comm);
    rev.ReadParams("input.db", false);
    rev.average_pore_size = 0.5;
    rev.StatRevAnalysis(model);
    require(rev.stat_poro_rev.size == 2 && close(rev.stat_poro_rev.mean, 0.5),
            "Checkerboard porosity should first converge in a two-voxel cube");
    require(close(rev.stat_poro_rev.surface_mean, 12.0 / (8.0 * 2.5)),
            "Surface target must be retained at the porosity REV, not at a later statistical size");
}

void testSurfaceConventionsAndUnits(const Utilities::MPI &comm) {
    TestDirectory directory;
    writeInput(3, {{3, 3, 3}}, {{0, 0, 0}}, "8bit");
    std::vector<uint8_t> source(27, 17);
    source[13] = 29;
    writeRaw(source);
    Model model(3, comm);
    model.Distance(2, 2, 2) = 1.0;
    double previous_surface = 0.0;
    for (double h : {2.5, 1.25}) {
        model.Dm->voxel_length = h;
        ScaLBL_REV rev(0, 1, comm);
        rev.ReadParams("input.db", false);
        const double surface = 6.0 / (27.0 * h);
        setREV(rev.stat_poro_rev, 3, 1.0 / 27.0, surface);
        rev.sub_domain_select(model);
        const Row row = readSummary().at("porosity");
        requireWeightedDiagnostics(row, surface, surface, 0.0);
        if (previous_surface > 0.0)
            require(close(std::stod(row.at("selected_surface")), 2.0 * previous_surface),
                    "Specific surface must scale inversely with voxel length");
        previous_surface = std::stod(row.at("selected_surface"));
    }
    model.Dm->voxel_length = 2.5;
    // All interior voxels are pores but the surrounding halo is solid. The
    // artificial crop boundary must contribute no pore-solid surface faces.
    source.assign(27, 29);
    writeRaw(source);
    for (int z = 1; z <= 3; ++z)
        for (int y = 1; y <= 3; ++y)
            for (int x = 1; x <= 3; ++x)
                model.Distance(x, y, z) = 1.0;
    ScaLBL_REV full(0, 1, comm);
    full.ReadParams("input.db", false);
    setREV(full.stat_poro_rev, 3, 1.0, 0.0);
    full.sub_domain_select(model);
    requireWeightedDiagnostics(readSummary().at("porosity"), 0.0, 0.0, 0.0);

    // The existing statistical implementation requires a strict sign change;
    // an exact zero distance neither creates nor connects an interface face.
    source.assign(27, 17);
    source.front() = 29;
    writeRaw(source);
    model.Distance.fill(-1.0);
    model.Distance(1, 1, 1) = 1.0;
    model.Distance(2, 1, 1) = 0.0;
    ScaLBL_REV zero_distance(0, 1, comm);
    zero_distance.ReadParams("input.db", false);
    const double surface = 2.0 / (27.0 * 2.5);
    setREV(zero_distance.stat_poro_rev, 3, 1.0 / 27.0, surface);
    zero_distance.sub_domain_select(model);
    requireWeightedDiagnostics(readSummary().at("porosity"), surface, surface, 0.0);
}

void testRectangularUnitAxisSurface(const Utilities::MPI &comm) {
    TestDirectory directory;
    const std::array<int, 3> shape{{4, 3, 2}};
    writeInput(4, shape, {{0, 0, 0}}, "8bit");
    Database parameters("input.db");
    parameters.getDatabase("Domain")->putVector<int>("n", {4, 3, 2});
    {
        std::ofstream out("input.db");
        parameters.print(out);
    }
    Model model(shape, comm);
    std::vector<uint8_t> mask, source;
    for (int z = 0; z < shape[2]; ++z)
        for (int y = 0; y < shape[1]; ++y)
            for (int x = 0; x < shape[0]; ++x) {
                const bool pore = (x + y + z) % 2 != 0;
                mask.push_back(pore ? 1 : 0);
                source.push_back(pore ? 29 : 17);
                model.Distance(x + 1, y + 1, z + 1) = pore ? 1.0 : -1.0;
            }
    writeRaw(source);
    for (const auto size : {std::array<int, 3>{{3, 2, 1}},
                            std::array<int, 3>{{1, 2, 2}},
                            std::array<int, 3>{{2, 1, 2}}}) {
        const auto measured = measureWindow(mask, shape, {{0, 0, 0}}, size, 2.5);
        for (int d = 0; d < 3; ++d)
            if (size[d] == 1)
                require(measured.face_axes[d] == 0,
                        "A unit-width axis has no internal interface faces");
        ScaLBL_REV rev(0, 1, comm);
        rev.ReadParams("input.db", false);
        setREV(rev.stat_poro_rev, size[0], measured.porosity, measured.surface);
        rev.stat_poro_rev.dimensions = size;
        rev.sub_domain_select(model);
        const Row row = readSummary().at("porosity");
        requireWeightedDiagnostics(row, measured.surface, measured.surface, 0.0);
        require(std::stoi(row.at("x0")) == 0 && std::stoi(row.at("y0")) == 0 &&
                    std::stoi(row.at("z0")) == 0,
                "Equal rectangular checkerboard scores should retain the first window");
        const int expected = (shape[0] - size[0] + 1) * (shape[1] - size[1] + 1) *
                             (shape[2] - size[2] + 1);
        require(std::stoull(row.at("candidate_count")) == unsigned(expected),
                "Rectangular weighted selection must scan each unit-axis window");
        const std::string raw_file = "subdomain_porosity_" + std::to_string(size[0]) + "x" +
                                     std::to_string(size[1]) + "x" + std::to_string(size[2]) + ".raw";
        require(requireRawFile(row) == raw_file &&
                    readBytes(raw_file).size() == size_t(size[0] * size[1] * size[2]),
                "Rectangular raw filenames must use each selected voxel extent in x/y/z order");
        requireDomainFile(row);
    }
}

void testMissingStatisticalSurface(const Utilities::MPI &comm) {
    TestDirectory directory;
    writeInput(3, {{3, 3, 3}}, {{0, 0, 0}}, "8bit");
    writeRaw(std::vector<uint8_t>(27, 17));
    Model model(3, comm);
    for (double invalid : {std::numeric_limits<double>::quiet_NaN(),
                           std::numeric_limits<double>::infinity(), -1.0}) {
        ScaLBL_REV rev(0, 1, comm);
        rev.ReadParams("input.db", false);
        setREV(rev.stat_poro_rev, 2, 0.0, invalid);
        bool rejected = false;
        try {
            rev.sub_domain_select(model);
        } catch (const std::runtime_error &error) {
            rejected = std::string(error.what()).find("surface") != std::string::npos;
        }
        require(rejected, "Missing or invalid surface means must not create fabricated score targets");
        require(!exists("subdomain_porosity_2x2x2.raw") && !exists("geometrical-rev-domain.db"),
                "A score without a valid surface target must not export a selected domain");
    }
}

void testNumericalFlow(const Utilities::MPI &comm) {
    Array<signed char> channel(4, 4, 4);
    channel.fill(0);
    for (int z = 0; z < 4; ++z)
        for (int y = 1; y <= 2; ++y)
            for (int x = 1; x <= 2; ++x)
                channel(x, y, z) = 1;
    const auto solve = [&](double h, double force) {
        return rev_selection::simulatePermeability(channel, h, 1.0,
                    {{0.0, 0.0, force}}, 4000, 20, 1.0e-8, comm);
    };
    const double permeability = solve(1.0, 1.0e-5);
    require(std::isfinite(permeability) && permeability > 0.0,
            "An open numerical channel must have finite positive permeability");
    require(std::abs(solve(2.0, 1.0e-5) / permeability - 4.0) < 1.0e-8,
            "Numerical permeability must scale with voxel length squared");
    require(std::abs(solve(1.0, -1.0e-5) / permeability - 1.0) < 1.0e-6,
            "Reversing numerical flow direction must preserve channel permeability");
    const auto straight = rev_selection::simulateHydraulicProperties(channel, 1.0, 1.0,
                    {{0.0, 0.0, 1.0e-5}}, 4000, 20, 1.0e-8, comm);
    require(std::abs(straight.permeability / permeability - 1.0) < 1.0e-8 &&
                std::abs(straight.tortuosity - 1.0) < 1.0e-8,
            "A straight numerical channel must give unit hydraulic tortuosity from the same solve");
    Array<signed char> diagonal(4, 4, 4);
    diagonal.fill(0);
    for (int z = 0; z < 4; ++z)
        for (int y = 1; y <= 2; ++y)
            diagonal(z, y, z) = 1;
    const auto bent = rev_selection::simulateHydraulicProperties(diagonal, 1.0, 1.0,
                    {{0.0, 0.0, 1.0e-5}}, 4000, 20, 1.0e-8, comm);
    require(std::isfinite(bent.permeability) && bent.permeability > 0.0 &&
                std::abs(bent.tortuosity - std::sqrt(2.0)) < 1.0e-6,
            "A diagonal numerical channel must include its transverse velocity in tortuosity");
    Array<signed char> sealed(3, 3, 3);
    sealed.fill(0);
    sealed(1, 1, 1) = 1;
    require(std::abs(rev_selection::simulatePermeability(sealed, 1.0, 1.0,
                    {{0.0, 0.0, 1.0e-5}}, 4000, 20, 1.0e-8, comm)) < 1.0e-8,
            "A sealed numerical pore must converge to zero permeability");
    const auto no_flow = rev_selection::simulateHydraulicProperties(sealed, 1.0, 1.0,
                    {{0.0, 0.0, 1.0e-5}}, 4000, 20, 1.0e-8, comm);
    require(std::abs(no_flow.permeability) < 1.0e-8 && std::isnan(no_flow.tortuosity),
            "A sealed pore must report undefined tortuosity instead of a geometric substitute");
    bool rejected = false;
    try {
        rev_selection::simulatePermeability(channel, 1.0, 1.0,
                    {{0.0, 0.0, 1.0e-5}}, 2, 1, 1.0e-8, comm);
    } catch (const std::runtime_error &) {
        rejected = true;
    }
    require(rejected, "An unconverged numerical flow must not become a candidate permeability");
}

void testNumericalSelection(const Utilities::MPI &comm) {
    TestDirectory directory;
    writeInput(4, {{4, 4, 4}}, {{0, 0, 0}}, "8bit", true);
    Model model(4, comm);
    std::vector<uint8_t> source(64, 17);
    for (int z = 0; z < 4; ++z)
        for (int y = 1; y <= 2; ++y)
            for (int x = 1; x <= 2; ++x) {
                source[x + 4 * (y + 4 * z)] = 29;
                model.Distance(x + 1, y + 1, z + 1) = 1.0;
            }
    writeRaw(source);
    Array<signed char> reference(3, 3, 3);
    reference.fill(0);
    for (int z = 0; z < 3; ++z)
        for (int y = 1; y <= 2; ++y)
            for (int x = 1; x <= 2; ++x)
                reference(x, y, z) = 1;
    const double target = rev_selection::simulatePermeability(reference, 2.5, 1.0,
                    {{0.0, 0.0, 1.0e-5}}, 4000, 20, 1.0e-8, comm);
    ScaLBL_REV rev(0, 1, comm);
    rev.ReadParams("input.db", false);
    rev.geometrical = false;
    setREV(rev.stat_perm_rev, 3, target);
    // All saved model velocities are zero; a positive candidate measurement
    // must therefore come from a new numerical solve of the cropped geometry.
    rev.sub_domain_select(model);
    const Row row = readSummary().at("permeability");
    require(std::stoull(row.at("candidate_count")) == 8,
            "Numerical selection must evaluate all one-voxel window positions");
    require(std::stod(row.at("selected_value")) > 0.0 &&
                std::abs(std::stod(row.at("selected_value")) / target - 1.0) < 1.0e-8,
            "Hydraulic selection must use a fresh flow solve rather than saved local velocities");
    require(std::abs(std::stod(row.at("relative_deviation"))) < 1.0e-8,
            "Numerical selected permeability did not match its statistical target");
    require(requireRawFile(row) == "subdomain_permeability_3x3x3.raw" &&
                readBytes(requireRawFile(row)).size() == 27,
            "Numerical selection exported the wrong geometry size");
    requireHydraulicDiagnostics(row, 1.0, std::stod(row.at("selected_tortuosity")));
    require(std::abs(std::stod(row.at("selected_tortuosity")) - 1.0) < 1.0e-8,
            "A selected straight channel must use numerical unit tortuosity");
}

void testWeightedHydraulicSelection(const Utilities::MPI &comm) {
    TestDirectory directory;
    writeInput(5, {{5, 5, 5}}, {{0, 0, 0}}, "8bit", true);
    Model model(5, comm);
    const int shift[] = {1, 1, 2, 1, 2};
    std::vector<uint8_t> source(125, 17);
    for (int z = 0; z < 5; ++z)
        for (int y = 1; y <= 3; ++y)
            for (int x = shift[z]; x <= shift[z] + 1; ++x) {
                source[x + 5 * (y + 5 * z)] = 29;
                model.Distance(x + 1, y + 1, z + 1) = 1.0;
            }
    writeRaw(source);

    struct Candidate {
        std::array<int, 3> origin;
        rev_selection::HydraulicProperties measured;
    };
    std::vector<Candidate> candidates;
    Array<signed char> window(4, 4, 4);
    for (int z = 0; z <= 1; ++z)
        for (int y = 0; y <= 1; ++y)
            for (int x = 0; x <= 1; ++x) {
                for (int k = 0; k < 4; ++k)
                    for (int j = 0; j < 4; ++j)
                        for (int i = 0; i < 4; ++i)
                            window(i, j, k) = source[x + i + 5 * (y + j + 5 * (z + k))] == 29 ? 1 : 0;
                const auto measured = rev_selection::simulateHydraulicProperties(window, 2.5, 1.0,
                                {{0.0, 0.0, 1.0e-5}}, 4000, 20, 1.0e-8, comm);
                require(std::isfinite(measured.permeability) && measured.permeability > 0.0 &&
                            std::isfinite(measured.tortuosity) && measured.tortuosity >= 1.0 - 1.0e-10,
                        "Hydraulic compromise fixture must have resolved positive through-flow");
                candidates.push_back({{{x, y, z}}, measured});
            }

    // Select targets near a permeability midpoint, biased toward one flow
    // pattern, while matching the other pattern's tortuosity. Independently
    // enumerate all numerical candidates so the expected score cannot come
    // from the selector's saved velocity fields or a porosity approximation.
    int expected_index = -1;
    double mean_permeability = 0.0, mean_tortuosity = 0.0;
    double expected_score = std::numeric_limits<double>::infinity();
    for (size_t a = 0; a < candidates.size() && expected_index < 0; ++a) {
        for (size_t b = 0; b < candidates.size() && expected_index < 0; ++b) {
            if (std::abs(candidates[a].measured.tortuosity - candidates[b].measured.tortuosity) < 1.0e-5)
                continue;
            const double mean_k = 0.51 * candidates[a].measured.permeability +
                                  0.49 * candidates[b].measured.permeability;
            const double mean_t = candidates[b].measured.tortuosity;
            int permeability_only = -1, weighted = -1;
            double best_k_error = std::numeric_limits<double>::infinity();
            double best_score = std::numeric_limits<double>::infinity();
            for (size_t i = 0; i < candidates.size(); ++i) {
                const auto &measured = candidates[i].measured;
                const double k_error = std::abs(measured.permeability - mean_k) / mean_k;
                const double t_error = std::abs(measured.tortuosity - mean_t) / mean_t;
                const double score = 0.6 * k_error + 0.4 * t_error;
                if (k_error < best_k_error) {
                    best_k_error = k_error;
                    permeability_only = static_cast<int>(i);
                }
                if (score < best_score) {
                    best_score = score;
                    weighted = static_cast<int>(i);
                }
            }
            if (weighted != permeability_only) {
                expected_index = weighted;
                mean_permeability = mean_k;
                mean_tortuosity = mean_t;
                expected_score = best_score;
            }
        }
    }
    require(expected_index >= 0,
            "Hydraulic fixture must distinguish the 60/40 compromise from permeability-only selection");
    ScaLBL_REV rev(0, 1, comm);
    rev.ReadParams("input.db", false);
    rev.geometrical = false;
    setREV(rev.stat_perm_rev, 4, mean_permeability, 0.0, mean_tortuosity);
    const auto retained = rev.stat_perm_rev;
    rev.sub_domain_select(model);
    const Row row = readSummary().at("permeability");
    const auto &expected = candidates[expected_index];
    const std::array<int, 3> actual{{std::stoi(row.at("x0")), std::stoi(row.at("y0")),
                                    std::stoi(row.at("z0"))}};
    require(actual == expected.origin && std::stoull(row.at("candidate_count")) == candidates.size(),
            "Hydraulic selection did not minimize the independently enumerated weighted numerical score");
    require(std::abs(std::stod(row.at("selected_value")) / expected.measured.permeability - 1.0) < 1.0e-10 &&
                close(std::stod(row.at("selection_score")), expected_score),
            "Selected hydraulic property and score do not match independent numerical measurements");
    requireHydraulicDiagnostics(row, mean_tortuosity, expected.measured.tortuosity);
    std::vector<char> expected_raw;
    for (int z = actual[2]; z < actual[2] + 4; ++z)
        for (int y = actual[1]; y < actual[1] + 4; ++y)
            for (int x = actual[0]; x < actual[0] + 4; ++x)
                expected_raw.push_back(static_cast<char>(source[x + 5 * (y + 5 * z)]));
    require(readBytes(requireRawFile(row)) == expected_raw,
            "Hydraulic compromise output does not correspond to its winning numerical geometry");
    require(rev.stat_perm_rev.size == retained.size && rev.stat_perm_rev.mean == retained.mean &&
                rev.stat_perm_rev.tortuosity_mean == retained.tortuosity_mean,
            "Hydraulic scoring must reuse its retained statistical size and targets");
    requireDomainFile(row);
}

void testUnavailableStatisticalTortuosity(const Utilities::MPI &comm) {
    for (double target : {std::numeric_limits<double>::quiet_NaN(),
                           std::numeric_limits<double>::infinity(),
                           -std::numeric_limits<double>::infinity()}) {
        TestDirectory directory;
        writeInput(3, {{3, 3, 3}}, {{0, 0, 0}}, "8bit", true);
        writeRaw(std::vector<uint8_t>(27, 17));
        Model model(3, comm);
        ScaLBL_REV rev(0, 1, comm);
        rev.ReadParams("input.db", false);
        rev.geometrical = false;
        setREV(rev.stat_perm_rev, 2, 1.0, 0.0, target);
        rev.sub_domain_select(model);
        const Row row = readSummary().at("permeability");
        require(row.at("status") == "unavailable_tortuosity" &&
                    std::stoull(row.at("candidate_count")) == 0,
                "Missing or nonfinite tortuosity targets must be explicit instead of fabricated");
        require(!exists("subdomain_permeability_2x2x2.raw") && !exists("hydraulic-rev-domain.db"),
                "An unavailable statistical target must not export a scored hydraulic subdomain");
    }
    TestDirectory directory;
    writeInput(4, {{4, 4, 4}}, {{0, 0, 0}}, "8bit", true);
    writeRaw(std::vector<uint8_t>(64, 17));
    Model model(4, comm);
    ScaLBL_REV rev(0, 1, comm);
    rev.ReadParams("input.db", false);
    rev.geometrical = false;
    rev.average_pore_size = 0.5;
    rev.StatRevAnalysis(model);
    require(rev.stat_perm_rev.size == 1 && rev.stat_perm_rev.mean == 0.0 &&
                std::isnan(rev.stat_perm_rev.tortuosity_mean),
            "All-invalid statistical tortuosity samples must not turn into a displayed mean of one");
    rev.sub_domain_select(model);
    require(readSummary().at("permeability").at("status") == "unavailable_tortuosity",
            "Statistical zero-flow input must retain an unavailable hydraulic score target");
}

void testUndefinedCandidateTortuosity(const Utilities::MPI &comm) {
    TestDirectory directory;
    writeInput(5, {{5, 5, 5}}, {{0, 0, 0}}, "8bit", true);
    Model model(5, comm);
    std::vector<uint8_t> source(125, 17);
    for (int z = 0; z < 5; ++z) {
        source[25 * z] = 29;
        model.Distance(1, 1, z + 1) = 1.0;
    }
    writeRaw(source);
    Array<signed char> reference(3, 3, 3);
    reference.fill(0);
    for (int z = 0; z < 3; ++z)
        reference(0, 0, z) = 1;
    const auto target = rev_selection::simulateHydraulicProperties(reference, 2.5, 1.0,
                    {{0.0, 0.0, 1.0e-5}}, 4000, 20, 1.0e-8, comm);
    require(target.permeability > 0.0 && std::isfinite(target.tortuosity),
            "A resolved reference flow is needed to test undefined hydraulic candidates");
    ScaLBL_REV rev(0, 1, comm);
    rev.ReadParams("input.db", false);
    rev.geometrical = false;
    setREV(rev.stat_perm_rev, 3, target.permeability, 0.0, target.tortuosity);
    rev.sub_domain_select(model);
    const Row row = readSummary().at("permeability");
    require(row.at("status") == "selected" && std::stoull(row.at("candidate_count")) == 27,
            "Undefined zero-flow candidates must be counted while all valid windows are evaluated");
    require(std::stoi(row.at("x0")) == 0 && std::stoi(row.at("y0")) == 0 &&
                std::stoi(row.at("z0")) == 0,
            "Solid windows with undefined tortuosity must not displace a valid hydraulic optimum");
    requireHydraulicDiagnostics(row, target.tortuosity, target.tortuosity);
    require(close(std::stod(row.at("selection_score")), 0.0),
            "The first valid exact hydraulic match must retain its zero score");
}

void requireUnchangedREV(const ScaLBL_REV::StatisticalREVResult &actual,
                         const ScaLBL_REV::StatisticalREVResult &original) {
    require(actual.size == original.size && actual.dimensions == original.dimensions &&
                actual.mean == original.mean && actual.cv == original.cv &&
                (actual.surface_mean == original.surface_mean ||
                 (std::isnan(actual.surface_mean) && std::isnan(original.surface_mean))) &&
                (actual.tortuosity_mean == original.tortuosity_mean ||
                 (std::isnan(actual.tortuosity_mean) && std::isnan(original.tortuosity_mean))),
            "Selection must reuse the retained statistical REV unchanged");
}

void requireSelectionDimensions(const Row &row, std::array<int, 3> size, int reference_size) {
    require(std::stoi(row.at("rev_size_vx")) == reference_size &&
                close(std::stod(row.at("rev_size_um")), reference_size * 2.5),
            "Selection CSV must retain the statistical reference size");
    const std::array<const char *, 3> dimensions{{"nx", "ny", "nz"}};
    const std::array<const char *, 3> lengths{{"lx_um", "ly_um", "lz_um"}};
    for (int d = 0; d < 3; ++d)
        require(std::stoi(row.at(dimensions[d])) == size[d] &&
                    close(std::stod(row.at(lengths[d])), size[d] * 2.5),
                "Selection CSV must describe the overridden voxel and physical dimensions");
}

void testGeometricSizeOverride(const Utilities::MPI &comm) {
    TestDirectory directory;
    writeInput(4, {{6, 7, 8}}, {{1, 2, 3}}, "8bit", false, true, false, "true", "3");
    std::vector<uint8_t> source(6 * 7 * 8, 17);
    source[4 + 6 * (5 + 7 * 6)] = 29;
    writeRaw(source);
    const auto original = readBytes("source.raw");
    Model model(4, comm);
    model.Distance(4, 4, 4) = 1.0;
    ScaLBL_REV rev(0, 1, comm);
    rev.ReadParams("input.db", false);
    setREV(rev.stat_poro_rev, 2, 1.0 / 27.0, 3.0 / (27.0 * 2.5));
    const auto retained = rev.stat_poro_rev;
    for (const char *alias : {"subdomain_porosity_3x3x3.raw", "subdomain_porosity_3x3x3.raw.tmp"}) {
        require(link("source.raw", alias) == 0, "Cannot create override output alias");
        bool rejected = false;
        try {
            rev.sub_domain_select(model);
        } catch (const std::runtime_error &error) {
            rejected = std::string(error.what()).find("aliases") != std::string::npos;
        }
        require(rejected && !exists("subdomain_selection.csv") && readBytes("source.raw") == original,
                "Override output aliases must be rejected before writing any outputs");
        require(std::remove(alias) == 0, "Cannot remove override output alias");
    }
    rev.sub_domain_select(model);
    const Row row = readSummary().at("porosity");
    requireSelectionDimensions(row, {{3, 3, 3}}, 2);
    require(std::stoull(row.at("candidate_count")) == 8 &&
                std::stoi(row.at("x0")) == 1 && std::stoi(row.at("y0")) == 1 &&
                std::stoi(row.at("z0")) == 1 && std::stoi(row.at("x1")) == 4 &&
                std::stoi(row.at("y1")) == 4 && std::stoi(row.at("z1")) == 4,
            "Scalar override must scan cubic windows and select the exact retained-target match");
    require(close(std::stod(row.at("statistical_mean")), retained.mean) &&
                close(std::stod(row.at("selected_value")), retained.mean),
            "Geometry override must use the retained porosity target");
    requireWeightedDiagnostics(row, retained.surface_mean, retained.surface_mean, 0.0);
    std::vector<char> expected(27, 17);
    expected.back() = 29;
    require(requireRawFile(row) == "subdomain_porosity_3x3x3.raw" &&
                readBytes(requireRawFile(row)) == expected,
            "Override raw export must contain the selected three-voxel cube");
    require(std::stoi(row.at("source_x0")) == 2 && std::stoi(row.at("source_y0")) == 3 &&
                std::stoi(row.at("source_z0")) == 4,
            "Override selection must preserve the original-image offset");
    auto selected = requireDomainFile(row);
    std::vector<signed char> mapped(27, 0);
    mapped.back() = 1;
    requireDomainRoundTrip(selected, mapped, comm);
    requireUnchangedREV(rev.stat_poro_rev, retained);
}

void testRectangularSizeOverride(const Utilities::MPI &comm) {
    TestDirectory directory;
    writeInput(4, {{4, 4, 4}}, {{0, 0, 0}}, "8bit", false, true, false, "true", "3, 2, 1");
    Model model(4, comm);
    std::vector<uint8_t> source(64, 17);
    for (int z = 0; z < 4; ++z)
        for (int y = 0; y < 4; ++y)
            for (int x = 0; x < 4; ++x) {
                const bool pore = (x + y + z) % 2 != 0;
                source[x + 4 * (y + 4 * z)] = pore ? 29 : 17;
                model.Distance(x + 1, y + 1, z + 1) = pore ? 1.0 : -1.0;
            }
    writeRaw(source);
    ScaLBL_REV rev(0, 1, comm);
    rev.ReadParams("input.db", false);
    setREV(rev.stat_poro_rev, 2, 0.5, 12.0 / (8.0 * 2.5));
    const auto retained = rev.stat_poro_rev;
    rev.sub_domain_select(model);
    const Row row = readSummary().at("porosity");
    requireSelectionDimensions(row, {{3, 2, 1}}, 2);
    require(std::stoull(row.at("candidate_count")) == 24 &&
                std::stoi(row.at("x0")) == 0 && std::stoi(row.at("y0")) == 0 &&
                std::stoi(row.at("z0")) == 0,
            "Rectangular override must scan every valid window and preserve deterministic ties");
    const double surface = 7.0 / (6.0 * 2.5);
    const double score = 0.4 * std::abs(surface - retained.surface_mean) / retained.surface_mean;
    requireWeightedDiagnostics(row, retained.surface_mean, surface, score);
    std::vector<char> expected;
    for (int y = 0; y < 2; ++y)
        for (int x = 0; x < 3; ++x)
            expected.push_back(source[x + 4 * y]);
    require(requireRawFile(row) == "subdomain_porosity_3x2x1.raw" &&
                readBytes(requireRawFile(row)) == expected,
            "Rectangular override must control raw shape and x/y/z filename dimensions");
    requireDomainFile(row);
    requireUnchangedREV(rev.stat_poro_rev, retained);
}

void testHydraulicSizeOverride(const Utilities::MPI &comm) {
    TestDirectory directory;
    writeInput(4, {{4, 4, 4}}, {{0, 0, 0}}, "8bit", true, true, false, "true", "3");
    Model model(4, comm);
    std::vector<uint8_t> source(64, 17);
    Array<signed char> full(4, 4, 4), window(3, 3, 3);
    full.fill(0);
    window.fill(0);
    for (int z = 0; z < 4; ++z)
        for (int y = 1; y <= 2; ++y)
            for (int x = 1; x <= 2; ++x) {
                source[x + 4 * (y + 4 * z)] = 29;
                full(x, y, z) = 1;
                model.Distance(x + 1, y + 1, z + 1) = 1.0;
                if (z < 3)
                    window(x, y, z) = 1;
            }
    writeRaw(source);
    const auto simulate = [&](const Array<signed char> &geometry) {
        return rev_selection::simulateHydraulicProperties(geometry, 2.5, 1.0,
                        {{0.0, 0.0, 1.0e-5}}, 4000, 20, 1.0e-8, comm);
    };
    const auto target = simulate(full), expected = simulate(window);
    require(std::abs(expected.permeability / target.permeability - 1.0) > 0.1,
            "Hydraulic override fixture must distinguish the new measurement from its reference mean");
    ScaLBL_REV rev(0, 1, comm);
    rev.ReadParams("input.db", false);
    setREV(rev.stat_poro_rev, 2, 0.25, 32.0 / (64.0 * 2.5));
    setREV(rev.stat_perm_rev, 4, target.permeability, 0.0, target.tortuosity);
    const auto retained_poro = rev.stat_poro_rev, retained_perm = rev.stat_perm_rev;
    rev.sub_domain_select(model);
    const auto rows = readSummary();
    for (const char *property : {"porosity", "permeability"}) {
        const Row &row = rows.at(property);
        requireSelectionDimensions(row, {{3, 3, 3}}, std::string(property) == "porosity" ? 2 : 4);
        require(row.at("status") == "selected" && std::stoull(row.at("candidate_count")) == 8,
                "Both properties must use the shared override despite independent statistical sizes");
        require(readBytes(requireRawFile(row)).size() == 27,
                "Both property exports must use the overridden geometry size");
        requireDomainFile(row);
    }
    const Row &row = rows.at("permeability");
    require(close(std::stod(row.at("statistical_mean")), target.permeability) &&
                std::abs(std::stod(row.at("selected_value")) / expected.permeability - 1.0) < 1.0e-8,
            "Hydraulic override must simulate new candidate windows using the retained reference target");
    require(std::abs(std::stod(row.at("relative_deviation")) -
                     std::abs(expected.permeability - target.permeability) / target.permeability) < 1.0e-8,
            "Hydraulic override must report error relative to the original statistical mean");
    requireHydraulicDiagnostics(row, retained_perm.tortuosity_mean, expected.tortuosity);
    requireUnchangedREV(rev.stat_poro_rev, retained_poro);
    requireUnchangedREV(rev.stat_perm_rev, retained_perm);
}

void testSizeOverrideValidation(const Utilities::MPI &comm) {
    TestDirectory directory;
    for (const char *invalid : {"0", "-1", "2, 0, 3", "2, -1, 3", "1, 2", "1, 2, 3, 4", "2.5"}) {
        writeInput(4, {{4, 4, 4}}, {{0, 0, 0}}, "8bit", false, true, false, "true", invalid);
        ScaLBL_REV rev(0, 1, comm);
        bool rejected = false;
        try {
            rev.ReadParams("input.db", false);
        } catch (const std::runtime_error &error) {
            rejected = std::string(error.what()).find("REV.subdomain_size") != std::string::npos;
        }
        require(rejected, std::string("Invalid REV.subdomain_size was not diagnosed: ") + invalid);
    }
    for (const char *save_raw : {"false", "true"})
        for (const char *oversized : {"5", "2, 5, 3"}) {
            writeInput(4, {{4, 4, 4}}, {{0, 0, 0}}, "8bit", false, true, false, save_raw, oversized);
            Model model(4, comm);
            ScaLBL_REV rev(0, 1, comm);
            rev.ReadParams("input.db", false);
            setREV(rev.stat_poro_rev, 2, 0.5, 0.0);
            bool rejected = false;
            try {
                // No source file: bounds must be checked before reading it.
                rev.sub_domain_select(model);
            } catch (const std::runtime_error &error) {
                rejected = std::string(error.what()).find("REV.subdomain_size") != std::string::npos;
            }
            require(rejected && !exists("subdomain_selection.csv") && selectedRawFiles().empty() &&
                        !exists("geometrical-rev-domain.db") && !exists("hydraulic-rev-domain.db"),
                    "Oversized override must be rejected before all source reads and output writes");
        }
}

void testSizeOverrideRequiresStatisticalTargets(const Utilities::MPI &comm) {
    TestDirectory directory;
    writeInput(4, {{4, 4, 4}}, {{0, 0, 0}}, "8bit", true, true, false, "false", "3, 2, 1");
    Model model(4, comm);
    ScaLBL_REV rev(0, 1, comm);
    rev.ReadParams("input.db", false);
    // A size override cannot create converged references or require a source
    // image when neither statistical property has converged.
    rev.sub_domain_select(model);
    for (const auto &entry : readSummary())
        require(entry.second.at("status") == "not_converged" && entry.second.at("raw_file").empty() &&
                    entry.second.at("domain_file").empty(),
                "A size override must not bypass statistical convergence");
    require(rev.stat_poro_rev.size == 0 && rev.stat_perm_rev.size == 0 && selectedRawFiles().empty(),
            "A size override must not fabricate retained statistical results");
    writeRaw(std::vector<uint8_t>(64, 17));
    rev.geometrical = false;
    setREV(rev.stat_perm_rev, 2, 0.0, 0.0, std::numeric_limits<double>::quiet_NaN());
    rev.sub_domain_select(model);
    Row row = readSummary().at("permeability");
    require(row.at("status") == "unavailable_tortuosity", "Override must retain unavailable-target status");
    requireSelectionDimensions(row, {{3, 2, 1}}, 2);
    setREV(rev.stat_perm_rev, 2, 0.0);
    rev.sub_domain_select(model);
    row = readSummary().at("permeability");
    require(row.at("status") == "no_valid_tortuosity" && std::stoull(row.at("candidate_count")) == 24,
            "Invalid candidate status must retain the overridden scan count");
    requireSelectionDimensions(row, {{3, 2, 1}}, 2);
    require(row.at("raw_file").empty() && row.at("domain_file").empty() && selectedRawFiles().empty() &&
                !exists("geometrical-rev-domain.db") && !exists("hydraulic-rev-domain.db"),
            "Unavailable references and candidates must not export an override domain");
}

void testLayeredPorositySource(const Utilities::MPI &comm) {
    TestDirectory directory;
    writeInput(4, {{4, 4, 4}}, {{0, 0, 0}}, "8bit", false, true, true);
    std::vector<uint8_t> source(64, 17);
    source.back() = 29;
    writeRaw(source);
    Model model(4, comm);
    // Reservoir processing may turn original solid voxels into pore space.
    // The measured candidate and exported crop must both describe the source.
    model.Distance.fill(1.0);
    ScaLBL_REV rev(0, 1, comm);
    rev.ReadParams("input.db", false);
    setREV(rev.stat_poro_rev, 2, 0.125, 3.0 / (8.0 * 2.5));
    const auto retained = rev.stat_poro_rev;
    rev.sub_domain_select(model);
    const Row row = readSummary().at("porosity");
    requireWeightedDiagnostics(row, 3.0 / (8.0 * 2.5), 3.0 / (8.0 * 2.5), 0.0);
    for (const char *key : {"x0", "y0", "z0"})
        require(std::stoi(row.at(key)) == 2, "Layered selection did not scan the original final window");
    require(std::stoull(row.at("candidate_count")) == 27,
            "Layer settings must not exclude original-image window positions");
    require(close(std::stod(row.at("selected_value")), 0.125) &&
                close(std::stod(row.at("relative_deviation")), 0.0),
            "Layered selection measured synthetic geometry instead of the original source");
    std::vector<char> expected(8, 17);
    expected.back() = 29;
    require(readBytes(requireRawFile(row)) == expected,
            "Layered selection did not preserve the original 8-bit crop");
    const char *bytes = reinterpret_cast<const char *>(source.data());
    require(readBytes("source.raw") == std::vector<char>(bytes, bytes + source.size()),
            "Layer-compatible selection changed the source image");
    requireUnchangedREV(rev.stat_poro_rev, retained);
    auto selected_domain = requireDomainFile(row);
    std::vector<signed char> mapped(8, 0);
    mapped.back() = 1;
    requireDomainRoundTrip(selected_domain, mapped, comm);
}

void testLayeredSixteenBitSource(const Utilities::MPI &comm) {
    TestDirectory directory;
    writeInput(4, {{6, 7, 8}}, {{1, 2, 3}}, "16bit", false, true, true);
    // Domain::Decomp narrows 16-bit values before mapping: 273 -> 17 (solid),
    // 285 -> 29 (pore). Selection must preserve this analysis convention while
    // exporting the original two-byte labels without conversion.
    std::vector<uint16_t> source(6 * 7 * 8, 273);
    source[4 + 6 * (5 + 7 * 6)] = 285;
    writeRaw(source);
    Model model(4, comm);
    model.Distance.fill(1.0);
    ScaLBL_REV rev(0, 1, comm);
    rev.ReadParams("input.db", false);
    setREV(rev.stat_poro_rev, 2, 0.125, 3.0 / (8.0 * 2.5));
    const auto retained = rev.stat_poro_rev;
    rev.sub_domain_select(model);
    const Row row = readSummary().at("porosity");
    requireWeightedDiagnostics(row, 3.0 / (8.0 * 2.5), 3.0 / (8.0 * 2.5), 0.0);
    for (const char *key : {"x0", "y0", "z0"})
        require(std::stoi(row.at(key)) == 2, "Layered 16-bit source was not mapped like Domain::Decomp");
    require(std::stoi(row.at("source_x0")) == 3 && std::stoi(row.at("source_y0")) == 4 &&
                std::stoi(row.at("source_z0")) == 5,
            "Layered source mask lost Domain.offset");
    require(close(std::stod(row.at("selected_value")), 0.125) &&
                std::stoull(row.at("candidate_count")) == 27,
            "Layered 16-bit selection did not evaluate the complete original image");
    std::vector<uint16_t> expected(8, 273);
    expected.back() = 285;
    const char *bytes = reinterpret_cast<const char *>(expected.data());
    require(readBytes(requireRawFile(row)) == std::vector<char>(bytes, bytes + 16),
            "Layered 16-bit export changed source labels or byte order");
    requireUnchangedREV(rev.stat_poro_rev, retained);
    auto selected_domain = requireDomainFile(row);
    std::vector<signed char> mapped(8, 0);
    mapped.back() = 1;
    requireDomainRoundTrip(selected_domain, mapped, comm);
}

void testLayeredNumericalSelection(const Utilities::MPI &comm) {
    TestDirectory directory;
    writeInput(4, {{4, 4, 4}}, {{0, 0, 0}}, "8bit", true, true, true);
    std::vector<uint8_t> source(64, 17);
    for (int z = 0; z < 4; ++z)
        for (int y = 1; y <= 2; ++y)
            for (int x = 1; x <= 2; ++x)
                source[x + 4 * (y + 4 * z)] = 29;
    writeRaw(source);
    Model model(4, comm);
    // The model deliberately supplies no pore space. A positive permeability
    // can only come from a numerical solve using the original source mask.
    model.Distance.fill(-1.0);
    Array<signed char> reference(3, 3, 3);
    reference.fill(0);
    for (int z = 0; z < 3; ++z)
        for (int y = 1; y <= 2; ++y)
            for (int x = 1; x <= 2; ++x)
                reference(x, y, z) = 1;
    const double target = rev_selection::simulatePermeability(reference, 2.5, 1.0,
                    {{0.0, 0.0, 1.0e-5}}, 4000, 20, 1.0e-8, comm);
    ScaLBL_REV rev(0, 1, comm);
    rev.ReadParams("input.db", false);
    rev.geometrical = false;
    setREV(rev.stat_perm_rev, 3, target);
    const auto retained = rev.stat_perm_rev;
    rev.sub_domain_select(model);
    const Row row = readSummary().at("permeability");
    require(std::stoull(row.at("candidate_count")) == 8,
            "Layered permeability selection must evaluate every original-image window");
    require(std::stod(row.at("selected_value")) > 0.0 &&
                std::abs(std::stod(row.at("selected_value")) / target - 1.0) < 1.0e-8 &&
                std::abs(std::stod(row.at("relative_deviation"))) < 1.0e-8,
            "Layered permeability selection did not simulate the original rock geometry");
    requireHydraulicDiagnostics(row, 1.0, std::stod(row.at("selected_tortuosity")));
    require(std::abs(std::stod(row.at("selected_tortuosity")) - 1.0) < 1.0e-8,
            "Layered hydraulic scoring must measure tortuosity on the restored original channel");
    requireUnchangedREV(rev.stat_perm_rev, retained);
    const std::array<int, 3> origin{{std::stoi(row.at("x0")), std::stoi(row.at("y0")),
                                      std::stoi(row.at("z0"))}};
    std::vector<char> expected;
    for (int z = origin[2]; z < origin[2] + 3; ++z)
        for (int y = origin[1]; y < origin[1] + 3; ++y)
            for (int x = origin[0]; x < origin[0] + 3; ++x)
                expected.push_back(static_cast<char>(source[x + 4 * (y + 4 * z)]));
    require(readBytes(requireRawFile(row)) == expected,
            "Layered permeability export did not match the simulated source crop");
    requireDomainFile(row);
}

void testOriginalImageProtection(const Utilities::MPI &comm) {
    TestDirectory directory;
    writeInput(3, {{3, 3, 3}}, {{0, 0, 0}}, "8bit", true);
    writeRaw(std::vector<uint8_t>(27, 17));
    Model model(3, comm);
    const auto original = readBytes("source.raw");
    for (const char *reserved : {"subdomain_selection.csv", "subdomain_porosity_2x2x2.raw",
                                 "subdomain_permeability_1x1x1.raw", "subdomain_porosity_2x2x2.raw.tmp",
                                 "subdomain_permeability_1x1x1.raw.tmp", "geometrical-rev-domain.db",
                                 "hydraulic-rev-domain.db", "geometrical-rev-domain.db.tmp",
                                 "hydraulic-rev-domain.db.tmp"}) {
        require(link("source.raw", reserved) == 0,
                "Cannot create hard-linked fixture for original-image protection");
        ScaLBL_REV rev(0, 1, comm);
        rev.ReadParams("input.db", false);
        // The summary is also protected when no statistical REV converged:
        // writing only the status rows must never truncate an aliased source.
        if (std::string(reserved) != "subdomain_selection.csv") {
            setREV(rev.stat_poro_rev, 2, 0.0);
            setREV(rev.stat_perm_rev, 1, 0.0);
        }
        bool rejected = false;
        try {
            rev.sub_domain_select(model);
        } catch (const std::runtime_error &error) {
            rejected = std::string(error.what()).find("aliases") != std::string::npos;
        }
        require(rejected, std::string("Original-image alias was not rejected: ") + reserved);
        require(readBytes("source.raw") == original && readBytes(reserved) == original,
                std::string("Selection modified original bytes through an alias: ") + reserved);
        require(std::remove(reserved) == 0, "Cannot remove hard-linked test output");
    }
    for (const char *legacy : {"subdomain_porosity.raw", "subdomain_permeability.raw",
                               "subdomain_porosity.raw.tmp", "subdomain_permeability.raw.tmp"})
        require(link("source.raw", legacy) == 0, "Cannot create unused legacy output aliases");
    ScaLBL_REV rev(0, 1, comm);
    rev.ReadParams("input.db", false);
    setREV(rev.stat_poro_rev, 2, 0.0);
    rev.sub_domain_select(model);
    for (const char *legacy : {"subdomain_porosity.raw", "subdomain_permeability.raw",
                               "subdomain_porosity.raw.tmp", "subdomain_permeability.raw.tmp"})
        require(readBytes(legacy) == original,
                "Legacy raw filenames must remain untouched when dimensioned outputs are written");
}

void writeIDInput(int boundary_condition) {
    std::ofstream out("input.db");
    out << "Domain {\n  n = 3, 3, 3\n  nproc = 1, 1, 1\n"
        << "  BC = " << boundary_condition << "\n  voxel_length = 2.5\n}\n"
        << "REV {\n  Geometrical = true\n  Hydraulic = false\n  SubdomainSelect = true\n"
        << "  save_subdomain_raw = true\n}\n";
    require(out.good(), "Cannot write test ID input database");
}

void testIDSource(const Utilities::MPI &comm) {
    TestDirectory directory;
    writeIDInput(0);
    // ID.00000 stores signed labels and one halo cell per face. Distinct halo
    // labels expose accidental halo export or an incorrect source origin.
    std::vector<int8_t> padded(125, -100);
    for (int z = 1; z <= 3; ++z)
        for (int y = 1; y <= 3; ++y)
            for (int x = 1; x <= 3; ++x)
                padded[x + 5 * (y + 5 * z)] = -17;
    padded[3 + 5 * (3 + 5 * 3)] = 3;
    writeRaw(padded);
    require(std::rename("source.raw", "ID.00000") == 0, "Cannot create ID.00000 fixture");
    const auto original = readBytes("ID.00000");
    Model model(3, comm);
    model.Distance(3, 3, 3) = 1.0;
    ScaLBL_REV rev(0, 1, comm);
    rev.ReadParams("input.db", false);
    setREV(rev.stat_poro_rev, 2, 0.125, 3.0 / (8.0 * 2.5));
    rev.sub_domain_select(model);
    const Row row = readSummary().at("porosity");
    requireWeightedDiagnostics(row, 3.0 / (8.0 * 2.5), 3.0 / (8.0 * 2.5), 0.0);
    require(row.at("source_file") == "ID.00000" && row.at("read_type") == "signed8bit" &&
                std::stoi(row.at("bytes_per_voxel")) == 1,
            "ID source type and file metadata were not preserved");
    for (const char *key : {"x0", "y0", "z0"})
        require(std::stoi(row.at(key)) == 1, "ID selection coordinates must exclude halos");
    for (const char *key : {"source_x0", "source_y0", "source_z0"})
        require(std::stoi(row.at(key)) == 2, "ID source coordinates must include the stored halo offset");
    for (const char *key : {"source_x1", "source_y1", "source_z1"})
        require(std::stoi(row.at(key)) == 4, "ID source ending coordinates must be exclusive");
    for (const char *key : {"source_nx", "source_ny", "source_nz"})
        require(std::stoi(row.at(key)) == 5, "ID source dimensions must include stored halos");
    std::vector<char> expected(8, static_cast<char>(-17));
    expected.back() = 3;
    require(readBytes(requireRawFile(row)) == expected,
            "ID output must preserve signed labels while excluding stored halos");
    require(readBytes("ID.00000") == original, "ID selection must preserve the original image");
    auto selected_domain = requireDomainFile(row);
    std::vector<signed char> mapped(8, -17);
    mapped.back() = 3;
    requireDomainRoundTrip(selected_domain, mapped, comm);
    require(readBytes("ID.00000") == original,
            "Reloading the selected ID configuration changed the original image");

    writeIDInput(3);
    ScaLBL_REV pressure_bc(0, 1, comm);
    pressure_bc.ReadParams("input.db", false);
    setREV(pressure_bc.stat_poro_rev, 2, 0.125, 3.0 / (8.0 * 2.5));
    bool rejected = false;
    try {
        pressure_bc.sub_domain_select(model);
    } catch (const std::runtime_error &error) {
        rejected = std::string(error.what()).find("BC = 0") != std::string::npos;
    }
    require(rejected, "ID boundary conditions that introduce synthetic layers must be rejected");
    require(readBytes("ID.00000") == original,
            "Rejecting incompatible ID boundary conditions must preserve original bytes");
}
}

int main(int argc, char **argv) {
    Utilities::startup(argc, argv, false);
    int failures = 0;
    {
        Utilities::MPI comm(MPI_COMM_WORLD);
        ScaLBL_SetDevice(comm.getRank());
        UnitTest tests;
        const std::pair<const char *, void (*)(const Utilities::MPI &)> cases[] = {
            {"Weighted geometry compromise versus porosity-only selection", testWeightedGeometrySelection},
            {"Exhaustive boundary selection, coordinates, and source offset", testLastWindowAndOffset},
            {"Original 16-bit labels and deterministic ties", testSixteenBitLabelsAndTies},
            {"Independent REV sizes and zero statistical means", testIndependentSizesAndZeroMeans},
            {"Unconverged statistical REV", testNoConvergence},
            {"REV equal to the full domain", testFullDomainWindow},
            {"Nonzero local porosity relative to zero mean", testUnmatchedZeroMean},
            {"Selection boolean opt-out", testDisabled},
            {"Optional raw exports for geometric and hydraulic selections", testOptionalRawExports},
            {"Statistical first convergence and independent property means", testStatisticalRetention},
            {"Surface target retained at the porosity statistical REV", testStatisticalSurfaceRetention},
            {"Surface units, internal boundaries, and strict distance signs", testSurfaceConventionsAndUnits},
            {"Rectangular windows with each possible unit-width axis", testRectangularUnitAxisSurface},
            {"Missing and nonfinite statistical surface targets", testMissingStatisticalSurface},
            {"Numerical MRT flow, units, direction, and convergence", testNumericalFlow},
            {"Fresh numerical permeability across all sliding windows", testNumericalSelection},
            {"Weighted hydraulic compromise versus permeability-only selection", testWeightedHydraulicSelection},
            {"Unavailable statistical tortuosity targets", testUnavailableStatisticalTortuosity},
            {"Undefined hydraulic candidates counted without fabricated tortuosity", testUndefinedCandidateTortuosity},
            {"Scalar selection size override with retained targets and source protection", testGeometricSizeOverride},
            {"Rectangular selection size override and retained surface target", testRectangularSizeOverride},
            {"Shared selection size override with fresh hydraulic measurements", testHydraulicSizeOverride},
            {"Selection size override validation before output writes", testSizeOverrideValidation},
            {"Selection size override requires valid statistical references", testSizeOverrideRequiresStatisticalTargets},
            {"Original porosity geometry despite synthetic layer settings", testLayeredPorositySource},
            {"Original layered 16-bit labels, mapping, and offset", testLayeredSixteenBitSource},
            {"Original numerical flow geometry despite synthetic layer settings", testLayeredNumericalSelection},
            {"Original image protection against aliased result files", testOriginalImageProtection},
            {"Signed ID source, halo-free export, and boundary compatibility", testIDSource}};
        for (const auto &test : cases) {
            try {
                test.second(comm);
                tests.passes(test.first);
            } catch (const std::exception &error) {
                tests.failure(std::string(test.first) + ": " + error.what());
            }
        }
        tests.report();
        failures = static_cast<int>(tests.NumFailGlobal());
    }
    Utilities::shutdown();
    return failures;
}
