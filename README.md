# sss — Slab Shielding Simulator

**sss** (Slab Shielding Simulator) is a single-file, command-line Monte Carlo simulation built on [Geant4](https://geant4.web.cern.ch/) that transports radiation through a multi-layer material slab. It emits one primary particle per event from a compound radiation source, tracks every energy deposition inside the slab, and records particles that penetrate or backscatter across the world boundary. Per-event results can optionally be saved to a [ROOT](https://root.cern/) file using the modern RNTuple columnar format.

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
- [ROOT](https://root.cern/) with the `ROOTNTuple` and `ROOTNTupleUtil` components

## Building

```sh
cmake -B build
cmake --build build -j
```

The executable is produced at `build/sss`.

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

## Physics

The default reference physics list is `QBBC`. Any list provided by the Geant4 physics-list factory may be selected with `--phys-list`.

### Geometry and scoring

- The slab is a stack of layers along **+z**; its front face sits at `z = -d/2`, where `d` is the total thickness.
- The world boundary is placed snugly around the slab, separated by a tiny vacuum gap.
- The primary starts just in front of the slab front face and moves along **+z**.
- A particle leaving the world boundary with `z ≥ 0` counts as **penetrating**; with `z < 0` it counts as **backscattered**.
- For unstable primaries, the final step position is the decay vertex, so decay is included automatically in the primary-termination classification.
- Neutrinos are excluded from the world-boundary energy statistics by default; pass `--neutrinos` to include them.
