# sss — Slab Shielding Simulator

**sss** (Slab Shielding Simulator) is a single-file, command-line Monte Carlo simulation built on [Geant4](https://geant4.web.cern.ch/) that transports radiation through a multi-layer material slab. It emits one primary particle per event from a compound radiation source, tracks every energy deposition inside the slab, and records particles that penetrate or backscatter across the world boundary. Per-event results can optionally be saved to a [ROOT](https://root.cern/) file using the modern RNTuple columnar format.

The companion tool **sss-vis** reads such a ROOT file and builds, for every RNTuple in it, energy-deposition histograms — a 3D histogram of the deposition positions, a 1D histogram along z, and the three 2D projections — so the deposition pattern inside the slab can be inspected directly. See [sss-vis — visualization](#sss-vis--visualization).

The program is intended for shielding and dosimetry studies: estimate how much energy is deposited in each material layer, how much escapes the back face (penetration), how much returns toward the source (backscattering), and which particle species are involved. It is equally well suited for studying electromagnetic and hadronic shower development and the resulting energy deposition in multi-layer materials.

## Features

- **Multi-layer slab** — any number of material layers, built in order from the source, each with its own material and thickness.
- **Flexible material specification** — use a NIST material by name, or build custom materials from element atom counts, element mass fractions, or mixtures of predefined NIST materials.
- **Compound radiation source** — any mix of particles and energies (e.g. electrons, photons, neutrons), each with a relative intensity; the source type of every primary is drawn from the resulting multinomial distribution.
- **Per-event output** — optional ROOT RNTuple file, one ntuple per source category (`run{runId}_src{typeId}`), storing energy depositions, deposition positions, the responsible processes, and penetrating/backscattered particle spectra.
- **End-of-run statistics** — console report with penetration, deposition, and backscattering energy and primary-particle ratios, plus per-particle mean energy and RMS tables.
- **Multithreaded** — Geant4 MT support; by default one worker per CPU core.
- **Interactive mode** — optional GUI with trajectory visualization.
- **Macro support** — run Geant4 macro files in batch mode.

## Examples

The following examples each demonstrate one of the four material definition methods.

Simulate 10,000 10 MeV electrons through 5 cm of lead (NIST material):

```sh
sss -m Pb:5cm -s e-:10MeV -n 10000
```

Two-layer slab with a custom water layer defined by atom counts (H₂O) in front of lead, saving per-event results to a ROOT file:

```sh
sss -m "water:[H:2,O:1]:1g/cm3:1mm;Pb:5cm" -s e+:10MeV -n 10000 -o out.root
```

Custom borated polyethylene defined by mass fractions, exposed to 662 keV gammas:

```sh
sss -m "bpe:{H:0.1,C:0.8,B:0.1}:0.95g/cm3:5cm" -s gamma:0.662MeV -n 10000
```

Water-loaded concrete defined as a mixture of predefined NIST materials, shielding 1 MeV neutrons on 8 threads:

```sh
sss -m "wc:<G4_CONCRETE:0.95,G4_WATER:0.05>:2.3:10cm" -s neutron:1MeV -n 10000 -j 8
```

## Requirements

- CMake ≥ 3.16
- A C++20 compiler (GCC, Clang, or MSVC)
- [Geant4](https://geant4.web.cern.ch/) (installed and discoverable via `find_package`)
- [ROOT](https://root.cern/) with the `ROOTNTuple`, `ROOTNTupleUtil`, and `ROOTDataFrame` components

## Building

```sh
cmake -B build
cmake --build build -j
```

The executables are produced at `build/sss` and `build/sss-vis`.

## Usage

```
sss [options] [macroFile]
```

Run modes (at least one is required):

- `--n-event <count>` — simulate `<count>` events in batch mode.
- `--ui` — open an interactive UI session with visualization.
- `<macroFile>` — execute a Geant4 macro file in batch mode (cannot be combined with `--ui` or `--n-event`).

### Required options

| Option | Description |
| --- | --- |
| `-m, --material <spec>` | Material layers as semicolon-separated entries, built in order from the source. See [Material specification](#material-specification). |
| `-s, --rad-src <spec>` | Compound radiation source as semicolon-separated `particle:energy[:intensity]` entries. See [Source specification](#source-specification). |

### Optional options

| Option | Description |
| --- | --- |
| `-n, --n-event <count>` | Simulate `<count>` events in batch mode. May be combined with `--ui` to pre-run events before the interactive session opens. |
| `-l, --phys-list <name>` | Reference physics list name (default: `QBBC`). |
| `-N, --neutrinos` | Include neutrino kinetic energy in the world-boundary energy statistics (neutrinos are ignored by default). |
| `-j, --threads <count>` | Worker thread count; `1` runs sequentially, `> 1` runs multithreaded (default: all CPU cores). |
| `-v, --verbose <level>` | Geant4 verbosity level; `0` prints only the banner, progress, and summary (default: `0`). |
| `-i, --ui` | Start an interactive UI session with visualization. |
| `-o, --output [<file>]` | Save per-event results into a ROOT file, one RNTuple per source category. The file name may be omitted (default `sss_output.root`) or given as a value, e.g. `--output out.root`. Never overwrites an existing file unless `--force`. |
| `-f, --force` | Overwrite the output file if it already exists. |
| `-h, --help` | Print the usage message. |

### Material specification

A layer is either a NIST material or a custom material:

- **NIST material** — `name:thickness`, e.g. `Pb:5cm`.
- **Custom material** — `name:<composition>:density:thickness`, where `<composition>` is one of:
  1. **Atom counts** — `[elem1:n1,elem2:n2,...]`, e.g. `water:[H:2,O:1]:1g/cm3:1mm`.
  2. **Mass fractions** — `{elem1:f1,elem2:f2,...}`, e.g. `air:{N:0.7,O:0.3}:1g/cm3:1`.
  3. **Mixture of predefined materials** — `<mat1:f1,mat2:f2,...>`, e.g. `wc:<G4_CONCRETE:0.95,G4_WATER:0.05>:2.3:10mm`.

Layers are separated by semicolons. A bare density value is in `g/cm3`; `kg/m3` (or `kg/m^3`) is also accepted. Supported length units: `mm`, `cm`, `m`, `km`, `um`, `nm`, `angstrom`, `fm`. A single-component composition may omit its atom count or mass fraction.

### Source specification

A semicolon-separated list of `particle:energy[:intensity]` entries, e.g. `e+:10MeV:15;gamma:3MeV:12` or `neutron:1MeV`. The intensity is the *relative* intensity of the source and may be omitted only when the list contains a single source. Intensities are normalized to probabilities, and every event emits one primary whose source type is drawn from the resulting multinomial distribution. Supported energy units: `eV`, `keV`, `MeV`, `GeV`, `TeV`, `PeV`.

## Output

### ROOT file

With `--output`, one RNTuple per source category is written to the ROOT file, named `run{runId}_src{typeId}` (source type indices follow the order of the `--rad-src` list). Each event row contains:

| Field | Type | Description |
| --- | --- | --- |
| `event_id` | `int` | Event number. |
| `total_e_pen` | `float` | Total energy of particles penetrating the world boundary. |
| `particle_pen`, `theta_pen`, `phi_pen`, `e_pen` | vectors | For each penetrating particle: name, exit direction (polar/azimuthal angle), and kinetic energy. |
| `total_e_dep` | `float` | Total energy deposited in the slab. |
| `e_dep_{i}` | `float` | Energy deposited in layer `i` (per layer: `0 .. N-1`). |
| `particle_dep_{i}`, `x_dep_{i}`, `y_dep_{i}`, `z_dep_{i}`, `w_dep_{i}`, `proc_dep_{i}` | vectors | For each deposition in layer `i`: depositing particle name, deposition position, deposited energy, and the responsible process name. |
| `total_e_bsc` | `float` | Total energy of particles backscattered across the world boundary. |
| `particle_bsc`, `theta_bsc`, `phi_bsc`, `e_bsc` | vectors | For each backscattered particle: name, exit direction, and kinetic energy. |

The output file is never overwritten unless `--force` is given.

### Console statistics

At the end of each run, the program prints, for every source category and for the total:

- **Energy ratios** — penetration, total deposition, per-layer deposition, and backscattering, each as a percentage of the incident energy with a statistical error.
- **Primary-particle ratios** — penetration, deposition, and backscattering, each as a percentage of the event count with a statistical error.
- **Particle tables** — for penetrating and backscattered particles, sorted by count: particle name, mean energy `<E>`, energy RMS spreading, and count.

## sss-vis — visualization

`sss-vis` is the companion visualization tool of `sss`. It reads a ROOT file written with `--output` and, for every RNTuple in it (one per source category), builds histograms of the energy deposition positions weighted by the deposited energy:

- a **3D histogram** of the deposition positions in xyz,
- a **1D histogram** of the deposition along **z**,
- the three **2D projections** onto the **xy**, **xz**, and **yz** planes.

All histograms are stored in an output ROOT file (by default the input file name with `_vis` inserted before the extension, e.g. `sss_output.root` → `sss_output_vis.root`). The coordinates and the deposited-energy weights are in Geant4 native units (millimetres for positions, MeV for energy).

### Usage

```
sss-vis [inputFile] [options]
```

`inputFile` is the ROOT file written by `sss`; when omitted, the default is `sss_output.root` (the default output file of `sss`). The program enables ROOT implicit multi-threading and reads the RNTuples through `ROOT::RDataFrame`; all transformations and actions are lazy, so the event loops run only when a result is actually needed. RNTuples are processed in alphabetical order.

| Option | Description |
| --- | --- |
| `-o, --output <file>` | Output ROOT file (default: `inputFile` with `_vis` before the extension). |
| `-3, --bins-3d <nx>:<ny>:<nz>` | Bin counts of the 3D histogram (default: `30:30:30`). |
| `-2, --bins-2d <nx>:<ny>` | Bin counts of the 2D projection histograms (default: `100:100`). |
| `-1, --bins-1d <count>` | Bin count of the 1D z histogram (default: `300`). |
| `-x, --x-range <min>:<max>` | Manual x-axis range used by every histogram. |
| `-y, --y-range <min>:<max>` | Manual y-axis range used by every histogram. |
| `-z, --z-range <min>:<max>` | Manual z-axis range used by every histogram. |
| `-s, --xy-fraction <value>` | Fraction of the deposited energy kept by the automatic xy range (default: `0.9`). |
| `-e, --z-expand <factor>` | Expansion factor of the automatic z range (default: `1.2`). |
| `-j, --threads <count>` | Implicit-multithreading worker count (default: all CPU cores). |
| `-f, --force` | Overwrite the output file if it already exists. |
| `-h, --help` | Print the usage message. |

### Automatic histogram ranges

When a manual range is not given:

- the **x** and **y** ranges are the *narrowest intervals centred at 0* that contain `--xy-fraction` (default 0.9) of the total deposited energy, trimmed symmetrically from both tails; long distribution tails therefore extend the range only by the energy they actually carry;
- the **z** range is `--z-expand` (default 1.2) times the interval that contains energy deposition, centred at `z = 0`.

A degenerate range (all depositions at the same position) is widened by a tiny symmetric padding so that the histograms always have a finite span.

### Output histograms

For an RNTuple named `run0_src0`, the following histograms are written:

| Histogram name | Description |
| --- | --- |
| `run0_src0_dep_z` | 1D histogram of the energy deposition along z. |
| `run0_src0_dep3d` | 3D histogram of the energy deposition positions. |
| `run0_src0_dep_xy`, `run0_src0_dep_xz`, `run0_src0_dep_yz` | 2D projections of the energy deposition onto the xy, xz, and yz planes. |

The output file is never overwritten unless `--force` is given.

### Example

```sh
sss -m Pb:5cm -s e-:10MeV -n 10000 -o out.root
sss-vis out.root
# writes out_vis.root containing run0_src0_dep_z, run0_src0_dep3d,
# run0_src0_dep_xy, run0_src0_dep_xz, and run0_src0_dep_yz
```

## Physics

The default reference physics list is `QBBC`. Any list provided by the Geant4 physics-list factory may be selected with `--phys-list`.

### Geometry and scoring

- The slab is a stack of layers along **+z**; its front face sits at `z = -d/2`, where `d` is the total thickness.
- The world boundary is placed snugly around the slab, separated by a tiny vacuum gap.
- The primary starts just in front of the slab front face and moves along **+z**.
- A particle leaving the world boundary with `z ≥ 0` counts as **penetrating**; with `z < 0` it counts as **backscattered**.
- For unstable primaries, the final step position is the decay vertex, so decay is included automatically in the primary-termination classification.
- Neutrinos are excluded from the world-boundary energy statistics by default; pass `--neutrinos` to include them.
- When a step's energy deposition is caused by multiple Coulomb scattering (`msc`), the recorded deposition position is sampled uniformly on the segment connecting the pre-step and post-step positions, because `msc` deposits energy throughout the step; otherwise the post-step position is recorded.
