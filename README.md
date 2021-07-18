# REstretto
REstretto (REuse of sub-STRuctures as an Effective Technique for protein-ligand docking TOol): An virtual screening-oriented protein-ligand docking tool with reuse of fragments

## Installation

### Requirement
- Boost
- Open Babel 2.x

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


## Configs

### INNERBOX
`INNERBOX [X], [Y], [Z]`
- X: x-width of innerbox (search grid) [Å]
- Y: y-width of innerbox (search grid) [Å]
- Z: z-width of innerbox (search grid) [Å]

### OUTERBOX
`OUTERBOX [X], [Y], [Z]`
- X: x-width of outerbox (atom grid) [Å]
- Y: y-width of outerbox (atom grid) [Å]
- Z: z-width of outerbox (atom grid) [Å]

### BOX_CENTER
`BOX_CENTER [X], [Y], [Z]`
- X: x-coordinate of center of the pocket [Å]
- Y: y-coordinate of center of the pocket [Å]
- Z: z-coordinate of center of the pocket [Å]

### SEARCH_PITCH
`SEARCH_PITCH [X], [Y], [Z]`
- X: search pitch of x-direction [Å]
- Y: search pitch of y-direction [Å]
- Z: search pitch of z-direction [Å]

### SCORING_PITCH
`SCORING_PITCH [X], [Y], [Z]`
- X: atomgrid pitch of x-direction [Å]
- Y: atomgrid pitch of y-direction [Å]
- Z: atomgrid pitch of z-direction [Å]

### MEMORY_SIZE
`MEMORY_SIZE [MEM]`
- MEM: memory size used for fragment grid [MB]

### RECEPTOR
`RECEPTOR [RECEPTOR]`
- RECEPTOR: file path of receptor

### LIGAND
`LIGAND [LIGAND]`
- LIGAND: file path of ligand
  - You can specify more than one by multiple lines

### OUTPUT
`OUTPUT [OUTPUT]`
- OUTPUT: file path of docked structure

### GRID_FOLDER
`GRID_FOLDER [GRID]`
- GRID: folder path of atomgrid


## Contributions & Support
Contributions and suggestions are welcome as well as bug reports. We greatly appreciate if you create a new issue and/or a PR (pull request) for these, especially including examples.

## Reference
- Yanagisawa K, Kubota R, Yoshikawa Y, Ohue M, Akiyama Y. (in prep)
