# REstretto
REstretto (REuse of sub-STRuctures as an Effective Technique for protein-ligand docking TOol): An virtual screening-oriented protein-ligand docking tool with reuse of fragments

## Installation

Dockerfile is at `.devcontainer/Dockerfile`. It explains the requirement and what packages are needed.

### Requirement
- Boost (latest)
- Open Babel 2.4.1
  - **OpenBabel 3.x (latest) is not supported.**

### Build

```sh
make all
```

## Usage

### Execution Sample
```sh
./atomgrid-gen testdata/testgrid.in       # make atom grid
./conformer-docking testdata/testgrid.in  # execute docking
```

Each program output progress to log file (not stdout).

## Configs

Sample Input: [testdata/testgrid.in](testdata/testgrid.in)
```
INNERBOX 10, 10, 10
OUTERBOX 20, 20, 20
BOX_CENTER 0.3826, 81.705, 109.195
SEARCH_PITCH 1, 1, 1
SCORING_PITCH 0.25, 0.25, 0.25
MEMORY_SIZE 8000
RECEPTOR testdata/2HU4_A_r.pdbqt
LIGAND testdata/G39.mol2
LIGAND testdata/conformers.mol2
OUTPUT testdata/G39_docked.sdf
GRID_FOLDER testdata/grid
```


### INNERBOX
`INNERBOX` is a region where the center of a compound will be put in.
The larger `INNERBOX` makes wider search but larger calculation cost.

`INNERBOX [X], [Y], [Z]`
- X: x-width of innerbox (search grid) [Å]
- Y: y-width of innerbox (search grid) [Å]
- Z: z-width of innerbox (search grid) [Å]

### OUTERBOX
`OUTERBOX` is a region where all atoms of a compound will be put in.

`OUTERBOX [X], [Y], [Z]`
- X: x-width of outerbox (atom grid) [Å]
- Y: y-width of outerbox (atom grid) [Å]
- Z: z-width of outerbox (atom grid) [Å]

### BOX_CENTER
`BOX_CENTER` is the center of the `INNERBOX` and `OUTERBOX`.

`BOX_CENTER [X], [Y], [Z]`
- X: x-coordinate of center of the pocket [Å]
- Y: y-coordinate of center of the pocket [Å]
- Z: z-coordinate of center of the pocket [Å]

### SEARCH_PITCH
`SEARCH_PITCH` is a search interval. Smaller `SEARCH_PITCH` results in finer search, but takes longer calculation time.

`SEARCH_PITCH [X], [Y], [Z]`
- X: search pitch of x-direction [Å]
- Y: search pitch of y-direction [Å]
- Z: search pitch of z-direction [Å]

### SCORING_PITCH
`SCORING_PITCH` is an scoring interval of atom grids and fragment grids. The finer pitch makes the finer scoring, but takes longer calculation time.

`SCORING_PITCH [X], [Y], [Z]`
- X: atom grid pitch of x-direction [Å]
- Y: atom grid pitch of y-direction [Å]
- Z: atom grid pitch of z-direction [Å]

### MEMORY_SIZE
It defines how much memory will be dedicated for storing fragment grids.

`MEMORY_SIZE [MEM]`
- MEM: memory size used for fragment grid [MB]

### RECEPTOR
A receptor (protein) file path.

`RECEPTOR [RECEPTOR]`
- RECEPTOR: file path of receptor

### LIGAND
Ligands (compounds) file pathes.

`LIGAND [LIGAND]`
- LIGAND: file path of ligand
  - You can specify more than one by multiple lines

### OUTPUT
File path for docked structures (the 3D-coordinates of compounds are only output).

`OUTPUT [OUTPUT]`
- OUTPUT: file path of docked structure

### GRID_FOLDER
A path to directory storing atom grids.

`GRID_FOLDER [GRID]`
- GRID: directory path of atomgrid


## Contributions & Support
Contributions and suggestions are welcome as well as bug reports. We greatly appreciate if you create a new issue and/or a PR (pull request) for these, especially including examples.

## Reference
- Yanagisawa K, Kubota R, Yoshikawa Y, Ohue M, Akiyama Y. (in prep)
