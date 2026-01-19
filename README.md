# VASP Syntax Highlighting for VS Code

Comprehensive syntax highlighting, tooltips, and auto-completion for VASP (Vienna Ab initio Simulation Package) input and output files.

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

## Supported File Types

### Input Files

| File | Names/Extensions | Description |
|------|------------------|-------------|
| INCAR | `INCAR`, `INCAR.orig`, `*.incar` | Calculation parameters |
| POSCAR | `POSCAR`, `CONTCAR`, `*.poscar`, `*.vasp` | Atomic structure |
| KPOINTS | `KPOINTS`, `*.kpoints` | K-point mesh definition |

### Output Files

| File | Names/Extensions | Description |
|------|------------------|-------------|
| OUTCAR | `OUTCAR`, `*.outcar` | Main output with all results |
| OSZICAR | `OSZICAR`, `*.oszicar` | Iteration history (energies, convergence) |
| DOSCAR | `DOSCAR`, `*.doscar` | Density of states |
| EIGENVAL | `EIGENVAL`, `*.eigenval` | Band eigenvalues |
| PROCAR | `PROCAR`, `*.procar` | Projected DOS onto orbitals |
| IBZKPT | `IBZKPT`, `*.ibzkpt` | Irreducible k-points |
| PCDAT | `PCDAT`, `*.pcdat` | Pair correlation function |

### Volumetric Data Files

| File | Names/Extensions | Description |
|------|------------------|-------------|
| CHGCAR | `CHGCAR`, `CHG`, `*.chgcar` | Charge density |
| LOCPOT | `LOCPOT`, `*.locpot` | Local potential |
| ELFCAR | `ELFCAR`, `*.elfcar` | Electron localization function |
| PARCHG | `PARCHG`, `*.parchg` | Partial charge density |
| AECCAR | `AECCAR0`, `AECCAR1`, `AECCAR2` | All-electron charge density |

### Trajectory Files

| File | Names/Extensions | Description |
|------|------------------|-------------|
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

| Category | Parameters |
|----------|------------|
| **Basic** | ENCUT, PREC, ALGO, EDIFF, NELM, LREAL, GGA, METAGGA |
| **K-points** | ISMEAR, SIGMA, KSPACING, KGAMMA |
| **Ionic** | NSW, IBRION, ISIF, EDIFFG, POTIM, NFREE |
| **Spin** | ISPIN, MAGMOM, NUPDOWN |
| **Spin-Orbit** | LSORBIT, LNONCOLLINEAR, SAXIS, LORBMOM |
| **DFT+U** | LDAU, LDAUTYPE, LDAUL, LDAUU, LDAUJ |
| **Hybrid** | LHFCALC, HFSCREEN, AEXX, PRECFOCK |
| **Van der Waals** | IVDW, LUSE_VDW, VDW_S6, VDW_S8 |
| **Output** | LWAVE, LCHARG, LORBIT, NEDOS, LVTOT, LELF |
| **Parallel** | NCORE, NPAR, KPAR, LPLANE |
| **GW/BSE** | NOMEGA, ENCUTGW, NBANDSO, NBANDSV |
| **MLFF** | ML_LMLFF, ML_MODE, ML_RCUT1, ML_RCUT2 |
| **NEB** | IMAGES, SPRING, LCLIMB, ICHAIN |
| **MD** | MDALGO, SMASS, TEBEG, TEEND |
| **Mixing** | AMIX, BMIX, AMIX_MAG, IMIX |

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

| Command | Description |
|---------|-------------|
| `VASP: Open Wiki Documentation` | Open VASP Wiki for parameter under cursor |
| `VASP: Show System Info` | Display system/VASP environment information |

## Development

```bash
# Install dependencies
npm install

# Build
npm run compile

# Watch for changes
npm run watch

# Lint
npm run lint

# Package
make package
```

## Contributing

Contributions are welcome!

- Add new parameters to `src/extension.ts` in the `VASP_PARAMETERS` object
- Add new file type support in `package.json` and create grammar in `syntaxes/`
- Report issues at the GitHub repository

## License

MIT

## Links

- [VASP Official Site](https://www.vasp.at/)
- [VASP Wiki](https://www.vasp.at/wiki/)
- [GitHub Repository](https://github.com/jkitchin/vscode-vasp-syntax)
