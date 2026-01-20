# VASPsum for VS Code

Comprehensive syntax highlighting, tooltips, auto-completion, summarization and visualization for VASP (Vienna Ab initio Simulation Package) input and output files.

![img](./screenshot.png)

## Features

### Syntax Highlighting

Full syntax highlighting for all major VASP file types:

- **Input Files**: INCAR, POSCAR, KPOINTS
- **Structure Files**: POSCAR, CONTCAR, XDATCAR
- **Charge/Potential Files**: CHGCAR, CHG, LOCPOT, ELFCAR, PARCHG, AECCAR0/1/2
- **Output Files**: OUTCAR, OSZICAR, DOSCAR, EIGENVAL, PROCAR, IBZKPT, PCDAT

### Hover Tooltips

Hover over any INCAR/OUTCAR parameter to see:
- Parameter description
- Type information
- Default values
- Allowed values
- Direct link to official VASP Wiki

### Auto-completion

IntelliSense auto-completion for INCAR parameters with:
- Smart snippets for boolean parameters (`.TRUE.`/`.FALSE.`)
- Value suggestions for enumerated parameters
- Full documentation in completion items

### OUTCAR Navigation

Use the document outline (Ctrl+Shift+O) to navigate OUTCAR sections:
- Iteration steps
- Energy sections
- Force and stress tensors
- Convergence information

### System Info Command

Debug and support tool accessible via Command Palette:
- **Command**: `VASP: Show System Info`
- Displays extension, VS Code, and system information
- Detects VASP installation and environment variables
- Shows Python/ASE/pymatgen versions
- One-click copy to clipboard for sharing

### VASPsum - Unified Summary & Visualization

The main feature of this extension is the **VASPsum** command that provides a comprehensive view of your VASP calculation:

- **Command**: `VASPsum: Show Structure & Summary`
- **Shortcut**: `Ctrl+Shift+V` (Windows/Linux) / `Cmd+Shift+V` (Mac)
- **CodeLens**: Click the "VASPsum" link at the top of any VASP file

#### Structure Visualization (3Dmol.js)

When POSCAR or CONTCAR is present:
- Interactive 3D rotation, zoom, and pan
- Element-colored atoms using CPK color scheme
- Unit cell box display (toggle on/off)
- Switch between sphere and stick visualization styles
- Atom labels toggle
- **Unit cell repetition**: Expand the view in a, b, c directions (1-10x)
- **Constrained atom markers**: Atoms with selective dynamics constraints (F flags) are shown with:
  - Slight transparency (0.9 opacity)
  - Red 3D marker (X in XY plane + vertical line) visible from all angles
- **Theme-adaptive colors**: Background and unit cell lines automatically match your VS Code theme (light/dark)

#### Calculation Summary

When vasprun.xml is present:
- Total energy and energy per atom
- Fermi energy
- Convergence status (checks if relaxation completed before NSW steps)
- Maximum force and pressure
- Stress tensor (collapsible)
- Atomic forces table (collapsible, highlights forces > 0.05 eV/Å)
- Calculation parameters: XC functional, ENCUT, k-points, ISPIN, IBRION, NSW, EDIFF, EDIFFG

#### Relaxation Progress Viewer

For geometry optimizations (multiple ionic steps):
- Energy vs. ionic step graph (Chart.js)
- Interactive slider to step through the relaxation
- Red dot on graph shows current step
- Energy and max force display for selected step
- **Live structure update**: 3D viewer updates atom positions as you move the slider
- Constraint markers are preserved during playback

## Supported File Types

### Input Files

| File    | Names/Extensions                          | Description             |
|---------|-------------------------------------------|-------------------------|
| INCAR   | `INCAR`, `INCAR.orig`, `*.incar`          | Calculation parameters  |
| POSCAR  | `POSCAR`, `CONTCAR`, `*.poscar`, `*.vasp` | Atomic structure        |
| KPOINTS | `KPOINTS`, `*.kpoints`                    | K-point mesh definition |

### Output Files

| File     | Names/Extensions         | Description                               |
|----------|--------------------------|-------------------------------------------|
| OUTCAR   | `OUTCAR`, `*.outcar`     | Main output with all results              |
| OSZICAR  | `OSZICAR`, `*.oszicar`   | Iteration history (energies, convergence) |
| DOSCAR   | `DOSCAR`, `*.doscar`     | Density of states                         |
| EIGENVAL | `EIGENVAL`, `*.eigenval` | Band eigenvalues                          |
| PROCAR   | `PROCAR`, `*.procar`     | Projected DOS onto orbitals               |
| IBZKPT   | `IBZKPT`, `*.ibzkpt`     | Irreducible k-points                      |
| PCDAT    | `PCDAT`, `*.pcdat`       | Pair correlation function                 |

### Volumetric Data Files

| File   | Names/Extensions                | Description                    |
|--------|---------------------------------|--------------------------------|
| CHGCAR | `CHGCAR`, `CHG`, `*.chgcar`     | Charge density                 |
| LOCPOT | `LOCPOT`, `*.locpot`            | Local potential                |
| ELFCAR | `ELFCAR`, `*.elfcar`            | Electron localization function |
| PARCHG | `PARCHG`, `*.parchg`            | Partial charge density         |
| AECCAR | `AECCAR0`, `AECCAR1`, `AECCAR2` | All-electron charge density    |

### Trajectory Files

| File    | Names/Extensions       | Description             |
|---------|------------------------|-------------------------|
| XDATCAR | `XDATCAR`, `*.xdatcar` | MD trajectory positions |

## Installation

### From VS Code Marketplace

Search for "VASP Syntax" in the Extensions view (Ctrl+Shift+X).

### From Package

```bash
cd vscode-vasp-syntax
make install
```

### Development Mode

```bash
cd vscode-vasp-syntax
make install-dev
```

### Manual Installation

1. Build the extension:
   ```bash
   npm install
   npm run compile
   ```

2. Package as VSIX:
   ```bash
   npm install -g @vscode/vsce
   vsce package
   ```

3. Install in VS Code:
   ```bash
   code --install-extension vasp-syntax-1.0.0.vsix
   ```

## Usage

### Basic Usage

1. Open any VASP file (INCAR, POSCAR, CHGCAR, etc.)
2. The file will automatically be recognized and highlighted
3. Hover over parameters to see documentation
4. Use Ctrl+Space for auto-completion in INCAR files

### VASPsum - View Calculation Summary

Three ways to open the VASPsum panel:

1. **CodeLens**: Click the "VASPsum: Structure + Summary" link at the top of any VASP file
2. **Keyboard**: Press `Ctrl+Shift+V` (Windows/Linux) or `Cmd+Shift+V` (Mac)
3. **Command Palette**: Run "VASPsum: Show Structure & Summary"

The panel automatically detects available files in the directory:
- Shows 3D structure if POSCAR or CONTCAR exists
- Shows calculation summary if vasprun.xml exists
- Shows relaxation viewer with slider if multiple ionic steps

### Open VASP Wiki

1. Place cursor on a parameter name
2. Open Command Palette (Ctrl+Shift+P)
3. Run "VASP: Open Wiki Documentation"

### System Info (Debugging)

1. Open Command Palette (Ctrl+Shift+P)
2. Run "VASP: Show System Info"
3. Click "Copy All to Clipboard" to share for support

The System Info panel shows:
- Extension and VS Code versions
- Operating system details
- VASP environment variables (`VASP_PP_PATH`, `VASP_COMMAND`, etc.)
- VASP executable locations
- POTCAR directory detection
- Python/ASE/pymatgen versions

### OUTCAR Navigation

1. Open an OUTCAR file
2. Press Ctrl+Shift+O to open the outline
3. Click on sections to jump directly:
   - Iteration 1, 2, 3...
   - Free Energy
   - Total Force
   - Stress Tensor
   - Converged!

## Supported Parameters

Over 200 VASP parameters are documented with full descriptions:

| Category          | Parameters                                          |
|-------------------|-----------------------------------------------------|
| **Basic**         | ENCUT, PREC, ALGO, EDIFF, NELM, LREAL, GGA, METAGGA |
| **K-points**      | ISMEAR, SIGMA, KSPACING, KGAMMA                     |
| **Ionic**         | NSW, IBRION, ISIF, EDIFFG, POTIM, NFREE             |
| **Spin**          | ISPIN, MAGMOM, NUPDOWN                              |
| **Spin-Orbit**    | LSORBIT, LNONCOLLINEAR, SAXIS, LORBMOM              |
| **DFT+U**         | LDAU, LDAUTYPE, LDAUL, LDAUU, LDAUJ                 |
| **Hybrid**        | LHFCALC, HFSCREEN, AEXX, PRECFOCK                   |
| **Van der Waals** | IVDW, LUSE_VDW, VDW_S6, VDW_S8                      |
| **Output**        | LWAVE, LCHARG, LORBIT, NEDOS, LVTOT, LELF           |
| **Parallel**      | NCORE, NPAR, KPAR, LPLANE                           |
| **GW/BSE**        | NOMEGA, ENCUTGW, NBANDSO, NBANDSV                   |
| **MLFF**          | ML_LMLFF, ML_MODE, ML_RCUT1, ML_RCUT2               |
| **NEB**           | IMAGES, SPRING, LCLIMB, ICHAIN                      |
| **MD**            | MDALGO, SMASS, TEBEG, TEEND                         |
| **Mixing**        | AMIX, BMIX, AMIX_MAG, IMIX                          |

## Screenshots

### Syntax Highlighting

INCAR files are colorized by parameter category:

![INCAR highlighting](images/incar-highlight.png)

### Hover Tooltips

Hover over any parameter to see documentation:

![Hover tooltip](images/hover-tooltip.png)

### Auto-completion

Type to get suggestions with documentation:

![Completion](images/completion.png)

## Commands

| Command                             | Shortcut                       | Description                                                                |
|-------------------------------------|--------------------------------|----------------------------------------------------------------------------|
| `VASPsum: Show Structure & Summary` | `Ctrl+Shift+V` / `Cmd+Shift+V` | Unified view with 3D structure, calculation summary, and relaxation viewer |
| `VASP: Open Wiki Documentation`     | -                              | Open VASP Wiki for parameter under cursor                                  |
| `VASP: Show System Info`            | -                              | Display system/VASP environment information                                |


## Contributing

Contributions are welcome!

- Report issues and feature requests at the GitHub repository
- We also consider pull requests

## License

MIT

## Links

- [VASP Official Site](https://www.vasp.at/)
- [VASP Wiki](https://www.vasp.at/wiki/)
- [GitHub Repository](https://github.com/jkitchin/vscode-vaspsum)
