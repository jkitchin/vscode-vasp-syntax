import * as vscode from 'vscode';
import * as os from 'os';
import * as path from 'path';
import * as fs from 'fs';
import { execSync } from 'child_process';

// VASP parameter documentation with links to official VASP wiki
interface VaspParameter {
    description: string;
    type: string;
    default?: string;
    values?: string[];
    category: string;
    wikiPage: string;
}

const VASP_WIKI_BASE = 'https://www.vasp.at/wiki/index.php/';

const VASP_PARAMETERS: Record<string, VaspParameter> = {
    // ============================================================
    // Basic Parameters
    // ============================================================
    'SYSTEM': {
        description: 'Name of the system. Used as title in output files.',
        type: 'string',
        default: 'unknown system',
        category: 'Basic',
        wikiPage: 'SYSTEM'
    },
    'ENCUT': {
        description: 'Plane-wave energy cutoff in eV. Higher values give more accurate but slower calculations.',
        type: 'float',
        default: 'largest ENMAX in POTCAR',
        category: 'Basic',
        wikiPage: 'ENCUT'
    },
    'ENAUG': {
        description: 'Kinetic energy cutoff for the augmentation charges in eV.',
        type: 'float',
        default: 'EAUG from POTCAR',
        category: 'Basic',
        wikiPage: 'ENAUG'
    },
    'PREC': {
        description: 'Precision mode controlling FFT grids and other settings.',
        type: 'string',
        default: 'Normal',
        values: ['Low', 'Medium', 'High', 'Normal', 'Accurate', 'Single'],
        category: 'Basic',
        wikiPage: 'PREC'
    },
    'ALGO': {
        description: 'Electronic minimization algorithm.',
        type: 'string',
        default: 'Normal',
        values: ['Normal', 'VeryFast', 'Fast', 'Conjugate', 'All', 'Damped', 'Subrot', 'Eigenval', 'Exact', 'None', 'Nothing', 'CHI', 'GW0', 'GW', 'scGW0', 'scGW', 'EVGW0', 'EVGW', 'BSE'],
        category: 'Basic',
        wikiPage: 'ALGO'
    },
    'EDIFF': {
        description: 'Energy convergence criterion for electronic self-consistency in eV.',
        type: 'float',
        default: '1E-4',
        category: 'Basic',
        wikiPage: 'EDIFF'
    },
    'NELM': {
        description: 'Maximum number of electronic self-consistency steps.',
        type: 'integer',
        default: '60',
        category: 'Basic',
        wikiPage: 'NELM'
    },
    'NELMIN': {
        description: 'Minimum number of electronic self-consistency steps.',
        type: 'integer',
        default: '2',
        category: 'Basic',
        wikiPage: 'NELMIN'
    },
    'NELMDL': {
        description: 'Number of non-self-consistent steps at the beginning.',
        type: 'integer',
        default: '-5 (ISTART=0) or 0 (ISTART>0)',
        category: 'Basic',
        wikiPage: 'NELMDL'
    },
    'LREAL': {
        description: 'Determines whether projection operators are evaluated in real or reciprocal space.',
        type: 'string',
        default: '.FALSE.',
        values: ['.FALSE.', '.TRUE.', 'Auto', 'On', 'Off'],
        category: 'Basic',
        wikiPage: 'LREAL'
    },
    'ADDGRID': {
        description: 'Adds an additional support grid for augmentation charges.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Basic',
        wikiPage: 'ADDGRID'
    },
    'LASPH': {
        description: 'Include non-spherical contributions from gradient corrections inside PAW spheres.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Basic',
        wikiPage: 'LASPH'
    },
    'LPAW': {
        description: 'Use PAW (Projector Augmented Wave) method. Always .TRUE. for PAW potentials.',
        type: 'boolean',
        default: '.TRUE. (for PAW POTCAR)',
        category: 'Basic',
        wikiPage: 'LPAW'
    },
    'EAUG': {
        description: 'Augmentation charge cutoff from POTCAR (eV).',
        type: 'float',
        category: 'POTCAR',
        wikiPage: 'EAUG'
    },
    'TITEL': {
        description: 'Title/name of the pseudopotential from POTCAR.',
        type: 'string',
        category: 'POTCAR',
        wikiPage: 'POTCAR'
    },
    'LULTRA': {
        description: 'Use ultrasoft pseudopotential.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'POTCAR',
        wikiPage: 'POTCAR'
    },
    'IUNSCR': {
        description: 'Unscreening type: 0=linear, 1=nonlinear, 2=none.',
        type: 'integer',
        values: ['0', '1', '2'],
        category: 'POTCAR',
        wikiPage: 'POTCAR'
    },
    'RPACOR': {
        description: 'Partial core radius from POTCAR (a.u.).',
        type: 'float',
        category: 'POTCAR',
        wikiPage: 'POTCAR'
    },
    'RCORE': {
        description: 'Outmost cutoff radius from POTCAR (a.u.).',
        type: 'float',
        category: 'POTCAR',
        wikiPage: 'POTCAR'
    },
    'ENMAX': {
        description: 'Maximum recommended energy cutoff from POTCAR (eV).',
        type: 'float',
        category: 'POTCAR',
        wikiPage: 'POTCAR'
    },
    'ENMIN': {
        description: 'Minimum recommended energy cutoff from POTCAR (eV).',
        type: 'float',
        category: 'POTCAR',
        wikiPage: 'POTCAR'
    },
    'ICORE': {
        description: 'Core electron treatment from POTCAR.',
        type: 'integer',
        category: 'POTCAR',
        wikiPage: 'POTCAR'
    },
    'LCOR': {
        description: 'Correct augmentation charges.',
        type: 'boolean',
        category: 'POTCAR',
        wikiPage: 'POTCAR'
    },
    'VRHFIN': {
        description: 'Valence electron configuration from POTCAR.',
        type: 'string',
        category: 'POTCAR',
        wikiPage: 'POTCAR'
    },
    'GGA': {
        description: 'Type of GGA functional.',
        type: 'string',
        default: 'PE (PBE)',
        values: ['PE', 'RP', 'RE', 'PS', 'AM', 'MK', 'BO', 'OR', 'ML', 'B3', 'BF'],
        category: 'Basic',
        wikiPage: 'GGA'
    },
    'METAGGA': {
        description: 'Type of meta-GGA functional.',
        type: 'string',
        values: ['SCAN', 'TPSS', 'RTPSS', 'M06L', 'MBJ'],
        category: 'Basic',
        wikiPage: 'METAGGA'
    },
    'ISTART': {
        description: 'Determines how to start or restart a job.',
        type: 'integer',
        default: '0 (if no WAVECAR), 1 (if WAVECAR exists)',
        values: ['0', '1', '2'],
        category: 'Basic',
        wikiPage: 'ISTART'
    },
    'INIWAV': {
        description: 'Initial electronic wavefunction: 0=lowest eigenvalue of H, 1=random.',
        type: 'integer',
        default: '1',
        values: ['0', '1'],
        category: 'Basic',
        wikiPage: 'INIWAV'
    },
    'ICHARG': {
        description: 'Determines how the initial charge density is constructed.',
        type: 'integer',
        default: '2 (ISTART=0), 0 (ISTART>0)',
        values: ['0', '1', '2', '4', '10', '11', '12'],
        category: 'Basic',
        wikiPage: 'ICHARG'
    },

    // ============================================================
    // K-points and Smearing
    // ============================================================
    'ISMEAR': {
        description: 'Smearing method for partial occupancies. -5=tetrahedron, 0=Gaussian, 1-N=Methfessel-Paxton.',
        type: 'integer',
        default: '1',
        values: ['-5', '-4', '-3', '-2', '-1', '0', '1', '2'],
        category: 'Smearing',
        wikiPage: 'ISMEAR'
    },
    'SIGMA': {
        description: 'Width of smearing in eV.',
        type: 'float',
        default: '0.2',
        category: 'Smearing',
        wikiPage: 'SIGMA'
    },
    'KSPACING': {
        description: 'Distance between k-points in reciprocal space (1/Å). Automatically generates KPOINTS.',
        type: 'float',
        category: 'K-points',
        wikiPage: 'KSPACING'
    },
    'KGAMMA': {
        description: 'If .TRUE., the generated k-point grid is centered at the Gamma point.',
        type: 'boolean',
        default: '.TRUE.',
        category: 'K-points',
        wikiPage: 'KGAMMA'
    },

    // ============================================================
    // Ionic Relaxation
    // ============================================================
    'NSW': {
        description: 'Maximum number of ionic steps.',
        type: 'integer',
        default: '0',
        category: 'Ionic',
        wikiPage: 'NSW'
    },
    'IBRION': {
        description: 'Ionic relaxation algorithm. -1=no update, 0=MD, 1=quasi-Newton, 2=CG, 5-8=phonons.',
        type: 'integer',
        default: '-1 (NSW=0,1) or 0 (NSW>1)',
        values: ['-1', '0', '1', '2', '3', '5', '6', '7', '8', '44'],
        category: 'Ionic',
        wikiPage: 'IBRION'
    },
    'ISIF': {
        description: 'Stress tensor and what to relax. 0-1=atoms, 2=atoms+stress, 3=atoms+cell, 4-7=shape variants.',
        type: 'integer',
        default: '0 (IBRION=0) or 2 (IBRION>0)',
        values: ['0', '1', '2', '3', '4', '5', '6', '7'],
        category: 'Ionic',
        wikiPage: 'ISIF'
    },
    'EDIFFG': {
        description: 'Force convergence criterion. If negative, specifies max force in eV/Å.',
        type: 'float',
        default: 'EDIFF*10',
        category: 'Ionic',
        wikiPage: 'EDIFFG'
    },
    'POTIM': {
        description: 'Step width scaling for ionic movements (fs for MD, Å for relaxation).',
        type: 'float',
        default: '0.5',
        category: 'Ionic',
        wikiPage: 'POTIM'
    },
    'NFREE': {
        description: 'Number of displacements for finite difference phonon calculations.',
        type: 'integer',
        default: '2',
        values: ['1', '2', '4'],
        category: 'Ionic',
        wikiPage: 'NFREE'
    },
    'PSTRESS': {
        description: 'External pressure in kB. Used for enthalpy calculations (H = E + P*V).',
        type: 'float',
        default: '0',
        category: 'Ionic',
        wikiPage: 'PSTRESS'
    },

    // ============================================================
    // Spin and Magnetism
    // ============================================================
    'ISPIN': {
        description: 'Spin polarization: 1=non-spin-polarized, 2=spin-polarized.',
        type: 'integer',
        default: '1',
        values: ['1', '2'],
        category: 'Spin',
        wikiPage: 'ISPIN'
    },
    'MAGMOM': {
        description: 'Initial magnetic moment for each atom. For spin-polarized calculations.',
        type: 'float array',
        default: 'NIONS*1.0',
        category: 'Spin',
        wikiPage: 'MAGMOM'
    },
    'NUPDOWN': {
        description: 'Fix total magnetic moment to specified value (N_up - N_down).',
        type: 'integer',
        default: '-1 (not fixed)',
        category: 'Spin',
        wikiPage: 'NUPDOWN'
    },
    'LSORBIT': {
        description: 'Enable spin-orbit coupling. Requires non-collinear version of VASP.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Spin-Orbit',
        wikiPage: 'LSORBIT'
    },
    'LNONCOLLINEAR': {
        description: 'Enable non-collinear magnetism.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Spin-Orbit',
        wikiPage: 'LNONCOLLINEAR'
    },
    'SAXIS': {
        description: 'Quantization axis for spin (3D vector).',
        type: 'float array',
        default: '0 0 1',
        category: 'Spin-Orbit',
        wikiPage: 'SAXIS'
    },
    'LORBMOM': {
        description: 'Calculate and print orbital moments.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Spin-Orbit',
        wikiPage: 'LORBMOM'
    },

    // ============================================================
    // DFT+U
    // ============================================================
    'LDAU': {
        description: 'Enable DFT+U (LDA+U/GGA+U) calculations for strongly correlated systems.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'DFT+U',
        wikiPage: 'LDAU'
    },
    'LDAUTYPE': {
        description: 'Type of DFT+U approach: 1=Liechtenstein, 2=Dudarev (rotationally invariant).',
        type: 'integer',
        default: '2',
        values: ['1', '2', '4'],
        category: 'DFT+U',
        wikiPage: 'LDAUTYPE'
    },
    'LDAUL': {
        description: 'Angular momentum for each species (0=s, 1=p, 2=d, 3=f, -1=off).',
        type: 'integer array',
        category: 'DFT+U',
        wikiPage: 'LDAUL'
    },
    'LDAUU': {
        description: 'On-site Coulomb interaction U for each species (eV).',
        type: 'float array',
        category: 'DFT+U',
        wikiPage: 'LDAUU'
    },
    'LDAUJ': {
        description: 'On-site exchange interaction J for each species (eV).',
        type: 'float array',
        category: 'DFT+U',
        wikiPage: 'LDAUJ'
    },
    'LDAUPRINT': {
        description: 'Verbosity of DFT+U output: 0=silent, 1=occupancy matrix, 2=all.',
        type: 'integer',
        default: '0',
        values: ['0', '1', '2'],
        category: 'DFT+U',
        wikiPage: 'LDAUPRINT'
    },
    'LMAXMIX': {
        description: 'Maximum l-quantum number for one-center PAW charge densities passed to mixer.',
        type: 'integer',
        default: '2 (s,p,d) or 4 (for f-elements)',
        values: ['2', '4', '6'],
        category: 'DFT+U',
        wikiPage: 'LMAXMIX'
    },

    // ============================================================
    // Hybrid Functionals
    // ============================================================
    'LHFCALC': {
        description: 'Enable hybrid functional calculation with Hartree-Fock exchange.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Hybrid',
        wikiPage: 'LHFCALC'
    },
    'HFSCREEN': {
        description: 'Screening length for HSE-type hybrid functionals (1/Å). 0.2 for HSE06.',
        type: 'float',
        default: '0',
        category: 'Hybrid',
        wikiPage: 'HFSCREEN'
    },
    'AEXX': {
        description: 'Fraction of exact exchange. 0.25 for PBE0/HSE06.',
        type: 'float',
        default: '0 (GGA) or 0.25 (hybrid)',
        category: 'Hybrid',
        wikiPage: 'AEXX'
    },
    'AGGAX': {
        description: 'Fraction of GGA exchange.',
        type: 'float',
        default: '1 - AEXX',
        category: 'Hybrid',
        wikiPage: 'AGGAX'
    },
    'AGGAC': {
        description: 'Fraction of GGA correlation.',
        type: 'float',
        default: '1.0',
        category: 'Hybrid',
        wikiPage: 'AGGAC'
    },
    'ALDAC': {
        description: 'Fraction of LDA correlation.',
        type: 'float',
        default: '1.0',
        category: 'Hybrid',
        wikiPage: 'ALDAC'
    },
    'TIME': {
        description: 'Time step for ALGO=Damped. Typically 0.4 for hybrids.',
        type: 'float',
        default: '0.4',
        category: 'Hybrid',
        wikiPage: 'TIME'
    },
    'PRECFOCK': {
        description: 'FFT grid for exact exchange: Fast, Normal, Accurate.',
        type: 'string',
        default: 'Normal',
        values: ['Low', 'Medium', 'Fast', 'Normal', 'Accurate'],
        category: 'Hybrid',
        wikiPage: 'PRECFOCK'
    },
    'NKRED': {
        description: 'Reduce k-points for HF by this factor in all directions.',
        type: 'integer',
        default: '1',
        category: 'Hybrid',
        wikiPage: 'NKRED'
    },

    // ============================================================
    // Van der Waals
    // ============================================================
    'IVDW': {
        description: 'Van der Waals correction method. 0=none, 10=DFT-D2, 11=DFT-D3, 12=DFT-D3(BJ), 2=TS, 202=MBD.',
        type: 'integer',
        default: '0',
        values: ['0', '1', '2', '4', '10', '11', '12', '20', '21', '202', '263', '300'],
        category: 'Van der Waals',
        wikiPage: 'IVDW'
    },
    'LUSE_VDW': {
        description: 'Enable vdW-DF type non-local van der Waals functionals.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Van der Waals',
        wikiPage: 'LUSE_VDW'
    },
    'VDW_RADIUS': {
        description: 'Cutoff radius for vdW interactions (Å).',
        type: 'float',
        default: '50.0',
        category: 'Van der Waals',
        wikiPage: 'VDW_RADIUS'
    },
    'VDW_S6': {
        description: 'Global scaling parameter s6 for DFT-D2.',
        type: 'float',
        default: 'functional dependent',
        category: 'Van der Waals',
        wikiPage: 'VDW_S6'
    },
    'VDW_S8': {
        description: 'Scaling parameter s8 for DFT-D3.',
        type: 'float',
        default: 'functional dependent',
        category: 'Van der Waals',
        wikiPage: 'VDW_S8'
    },
    'VDW_A1': {
        description: 'Damping parameter a1 for DFT-D3(BJ).',
        type: 'float',
        category: 'Van der Waals',
        wikiPage: 'VDW_A1'
    },
    'VDW_A2': {
        description: 'Damping parameter a2 for DFT-D3(BJ).',
        type: 'float',
        category: 'Van der Waals',
        wikiPage: 'VDW_A2'
    },
    'PARAM1': {
        description: 'Parameter 1 for vdW-DF functionals.',
        type: 'float',
        category: 'Van der Waals',
        wikiPage: 'PARAM1'
    },
    'PARAM2': {
        description: 'Parameter 2 for vdW-DF functionals.',
        type: 'float',
        category: 'Van der Waals',
        wikiPage: 'PARAM2'
    },
    'ZAB_VDW': {
        description: 'vdW-DF kernel parameter. -1.8867 for vdW-DF2.',
        type: 'float',
        default: '-0.8491',
        category: 'Van der Waals',
        wikiPage: 'ZAB_VDW'
    },

    // ============================================================
    // Molecular Dynamics
    // ============================================================
    'MDALGO': {
        description: 'Molecular dynamics algorithm: 0=NVE, 1=Andersen, 2=Nose-Hoover, 3=Langevin.',
        type: 'integer',
        default: '0',
        values: ['0', '1', '2', '3', '11', '13', '21'],
        category: 'Molecular Dynamics',
        wikiPage: 'MDALGO'
    },
    'SMASS': {
        description: 'Nose mass parameter. -3=NVE (microcanonical), >=0 for Nose-Hoover thermostat.',
        type: 'float',
        default: '-3',
        category: 'Molecular Dynamics',
        wikiPage: 'SMASS'
    },
    'TEBEG': {
        description: 'Initial temperature for MD (Kelvin).',
        type: 'float',
        default: '0',
        category: 'Molecular Dynamics',
        wikiPage: 'TEBEG'
    },
    'TEEND': {
        description: 'Final temperature for MD simulated annealing (Kelvin).',
        type: 'float',
        default: 'TEBEG',
        category: 'Molecular Dynamics',
        wikiPage: 'TEEND'
    },
    'NBLOCK': {
        description: 'After NBLOCK steps, pair correlation function and DOS are calculated.',
        type: 'integer',
        default: '1',
        category: 'Molecular Dynamics',
        wikiPage: 'NBLOCK'
    },
    'KBLOCK': {
        description: 'After KBLOCK*NBLOCK steps, averages are written to XDATCAR.',
        type: 'integer',
        default: 'NSW',
        category: 'Molecular Dynamics',
        wikiPage: 'KBLOCK'
    },
    'LANGEVIN_GAMMA': {
        description: 'Friction coefficient for Langevin thermostat (ps^-1) per species.',
        type: 'float array',
        category: 'Molecular Dynamics',
        wikiPage: 'LANGEVIN_GAMMA'
    },
    'PMASS': {
        description: 'Mass of the lattice degree of freedom for NPT dynamics.',
        type: 'float',
        category: 'Molecular Dynamics',
        wikiPage: 'PMASS'
    },

    // ============================================================
    // Output Control
    // ============================================================
    'LWAVE': {
        description: 'Write WAVECAR file (wavefunctions).',
        type: 'boolean',
        default: '.TRUE.',
        category: 'Output',
        wikiPage: 'LWAVE'
    },
    'LCHARG': {
        description: 'Write CHGCAR file (charge density).',
        type: 'boolean',
        default: '.TRUE.',
        category: 'Output',
        wikiPage: 'LCHARG'
    },
    'LVTOT': {
        description: 'Write LOCPOT file (local potential).',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Output',
        wikiPage: 'LVTOT'
    },
    'LVHAR': {
        description: 'Write only Hartree potential to LOCPOT.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Output',
        wikiPage: 'LVHAR'
    },
    'LELF': {
        description: 'Write ELFCAR file (electron localization function).',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Output',
        wikiPage: 'LELF'
    },
    'LAECHG': {
        description: 'Write all-electron charge density (AECCAR0, AECCAR2).',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Output',
        wikiPage: 'LAECHG'
    },
    'LORBIT': {
        description: 'Write PROCAR/DOSCAR with projected DOS. 10-12 for atom-projected.',
        type: 'integer',
        default: '0',
        values: ['0', '1', '2', '5', '10', '11', '12', '13', '14'],
        category: 'Output',
        wikiPage: 'LORBIT'
    },
    'NEDOS': {
        description: 'Number of grid points in DOS.',
        type: 'integer',
        default: '301',
        category: 'Output',
        wikiPage: 'NEDOS'
    },
    'EMIN': {
        description: 'Minimum energy for DOS/band evaluation (eV relative to E_Fermi).',
        type: 'float',
        default: 'lowest eigenvalue',
        category: 'Output',
        wikiPage: 'EMIN'
    },
    'EMAX': {
        description: 'Maximum energy for DOS/band evaluation (eV relative to E_Fermi).',
        type: 'float',
        default: 'highest eigenvalue',
        category: 'Output',
        wikiPage: 'EMAX'
    },
    'NWRITE': {
        description: 'Verbosity of OUTCAR. 0=minimal, 1=less, 2=normal, 3=verbose, 4=debug.',
        type: 'integer',
        default: '2',
        values: ['0', '1', '2', '3', '4'],
        category: 'Output',
        wikiPage: 'NWRITE'
    },
    'NBANDS': {
        description: 'Number of bands included in the calculation.',
        type: 'integer',
        default: 'depends on NELECT',
        category: 'Electronic',
        wikiPage: 'NBANDS'
    },
    'NELECT': {
        description: 'Total number of electrons in the calculation.',
        type: 'float',
        default: 'from POTCAR',
        category: 'Electronic',
        wikiPage: 'NELECT'
    },

    // ============================================================
    // Symmetry
    // ============================================================
    'ISYM': {
        description: 'Symmetry mode: -1=off, 0=on (no symmetry of k), 1=auto, 2=force, 3=no time-reversal.',
        type: 'integer',
        default: '1 or 2',
        values: ['-1', '0', '1', '2', '3'],
        category: 'Symmetry',
        wikiPage: 'ISYM'
    },
    'SYMPREC': {
        description: 'Precision for finding symmetry (Å).',
        type: 'float',
        default: '1E-5',
        category: 'Symmetry',
        wikiPage: 'SYMPREC'
    },

    // ============================================================
    // Parallelization
    // ============================================================
    'NCORE': {
        description: 'Number of cores per orbital. NCORE * NPAR = total cores.',
        type: 'integer',
        default: '1',
        category: 'Parallel',
        wikiPage: 'NCORE'
    },
    'NPAR': {
        description: 'Number of groups for band parallelization.',
        type: 'integer',
        default: 'total cores / NCORE',
        category: 'Parallel',
        wikiPage: 'NPAR'
    },
    'KPAR': {
        description: 'Number of k-point groups for k-point parallelization.',
        type: 'integer',
        default: '1',
        category: 'Parallel',
        wikiPage: 'KPAR'
    },
    'LPLANE': {
        description: 'Use plane-wave parallelization (usually .TRUE.).',
        type: 'boolean',
        default: '.TRUE.',
        category: 'Parallel',
        wikiPage: 'LPLANE'
    },
    'LSCALU': {
        description: 'Use parallel LU decomposition.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Parallel',
        wikiPage: 'LSCALU'
    },
    'NSIM': {
        description: 'Number of bands handled simultaneously in RMM-DIIS.',
        type: 'integer',
        default: '4',
        category: 'Parallel',
        wikiPage: 'NSIM'
    },

    // ============================================================
    // Dipole Corrections
    // ============================================================
    'LDIPOL': {
        description: 'Enable dipole corrections for surface calculations.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Dipole',
        wikiPage: 'LDIPOL'
    },
    'IDIPOL': {
        description: 'Direction of dipole correction: 1=x, 2=y, 3=z, 4=all directions.',
        type: 'integer',
        values: ['1', '2', '3', '4'],
        category: 'Dipole',
        wikiPage: 'IDIPOL'
    },
    'DIPOL': {
        description: 'Center of the cell for dipole correction (fractional coordinates).',
        type: 'float array',
        default: 'center of ionic charges',
        category: 'Dipole',
        wikiPage: 'DIPOL'
    },
    'EPSILON': {
        description: 'Dielectric constant of bulk for dipole correction.',
        type: 'float',
        default: '1.0',
        category: 'Dipole',
        wikiPage: 'EPSILON'
    },
    'EFIELD': {
        description: 'Applied electric field in eV/Å.',
        type: 'float',
        default: '0',
        category: 'Dipole',
        wikiPage: 'EFIELD'
    },

    // ============================================================
    // GW and BSE
    // ============================================================
    'NOMEGA': {
        description: 'Number of frequency points for GW/dielectric function.',
        type: 'integer',
        default: '0',
        category: 'GW/BSE',
        wikiPage: 'NOMEGA'
    },
    'ENCUTGW': {
        description: 'Energy cutoff for response function in GW (eV).',
        type: 'float',
        default: '2/3 * ENCUT',
        category: 'GW/BSE',
        wikiPage: 'ENCUTGW'
    },
    'NBANDSO': {
        description: 'Number of occupied bands in BSE calculation.',
        type: 'integer',
        category: 'GW/BSE',
        wikiPage: 'NBANDSO'
    },
    'NBANDSV': {
        description: 'Number of virtual (unoccupied) bands in BSE calculation.',
        type: 'integer',
        category: 'GW/BSE',
        wikiPage: 'NBANDSV'
    },
    'ANTIRES': {
        description: 'Antiresonant contributions in BSE: 0=with, 1=without (Tamm-Dancoff).',
        type: 'integer',
        default: '0',
        values: ['0', '1', '2'],
        category: 'GW/BSE',
        wikiPage: 'ANTIRES'
    },
    'OMEGAMAX': {
        description: 'Maximum frequency for optical properties (eV).',
        type: 'float',
        default: '-30 to 30',
        category: 'GW/BSE',
        wikiPage: 'OMEGAMAX'
    },

    // ============================================================
    // Optical Properties
    // ============================================================
    'LOPTICS': {
        description: 'Calculate frequency-dependent dielectric function.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Optics',
        wikiPage: 'LOPTICS'
    },
    'CSHIFT': {
        description: 'Complex shift for dielectric function (eV). Broadening of spectra.',
        type: 'float',
        default: '0.1',
        category: 'Optics',
        wikiPage: 'CSHIFT'
    },
    'LEPSILON': {
        description: 'Calculate static dielectric tensor and Born effective charges.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Optics',
        wikiPage: 'LEPSILON'
    },
    'LRPA': {
        description: 'Include local field effects in RPA.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Optics',
        wikiPage: 'LRPA'
    },

    // ============================================================
    // Machine Learning Force Fields
    // ============================================================
    'ML_LMLFF': {
        description: 'Enable machine learning force field.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'MLFF',
        wikiPage: 'ML_LMLFF'
    },
    'ML_MODE': {
        description: 'MLFF operation mode: run, select, train, refit.',
        type: 'string',
        default: 'run',
        values: ['run', 'select', 'train', 'refit'],
        category: 'MLFF',
        wikiPage: 'ML_MODE'
    },
    'ML_ISTART': {
        description: 'When to start MLFF learning (ionic step).',
        type: 'integer',
        default: '0',
        category: 'MLFF',
        wikiPage: 'ML_ISTART'
    },
    'ML_RCUT1': {
        description: 'Cutoff radius for 2-body descriptors (Å).',
        type: 'float',
        default: '8.0',
        category: 'MLFF',
        wikiPage: 'ML_RCUT1'
    },
    'ML_RCUT2': {
        description: 'Cutoff radius for 3-body descriptors (Å).',
        type: 'float',
        default: '4.0',
        category: 'MLFF',
        wikiPage: 'ML_RCUT2'
    },
    'ML_MB': {
        description: 'Maximum number of 2-body basis functions.',
        type: 'integer',
        default: '8000',
        category: 'MLFF',
        wikiPage: 'ML_MB'
    },
    'ML_MB3': {
        description: 'Maximum number of 3-body basis functions.',
        type: 'integer',
        default: '8000',
        category: 'MLFF',
        wikiPage: 'ML_MB'
    },
    'ML_WTIFOR': {
        description: 'Weight for forces in MLFF training.',
        type: 'float',
        default: '10.0',
        category: 'MLFF',
        wikiPage: 'ML_WTIFOR'
    },
    'ML_WTOTEN': {
        description: 'Weight for total energy in MLFF training.',
        type: 'float',
        default: '1.0',
        category: 'MLFF',
        wikiPage: 'ML_WTOTEN'
    },
    'ML_WTSIF': {
        description: 'Weight for stress tensor in MLFF training.',
        type: 'float',
        default: '1.0',
        category: 'MLFF',
        wikiPage: 'ML_WTSIF'
    },

    // ============================================================
    // NEB (Nudged Elastic Band)
    // ============================================================
    'IMAGES': {
        description: 'Number of NEB images between initial and final states.',
        type: 'integer',
        category: 'NEB',
        wikiPage: 'IMAGES'
    },
    'SPRING': {
        description: 'Spring constant for NEB (eV/Å²).',
        type: 'float',
        default: '-5.0',
        category: 'NEB',
        wikiPage: 'SPRING'
    },
    'LCLIMB': {
        description: 'Enable climbing image NEB for better transition state.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'NEB',
        wikiPage: 'LCLIMB'
    },
    'LTANGENTOLD': {
        description: 'Use old tangent definition for NEB.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'NEB',
        wikiPage: 'LTANGENTOLD'
    },
    'LNEBCELL': {
        description: 'Allow cell shape to change in NEB.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'NEB',
        wikiPage: 'LNEBCELL'
    },
    'ICHAIN': {
        description: 'Chain type for NEB: 0=normal, 1=reversal.',
        type: 'integer',
        default: '0',
        values: ['0', '1', '2', '3'],
        category: 'NEB',
        wikiPage: 'ICHAIN'
    },

    // ============================================================
    // Phonons
    // ============================================================
    'DQ': {
        description: 'Displacement for finite difference phonon calculations (Å).',
        type: 'float',
        default: '0.015',
        category: 'Phonons',
        wikiPage: 'IBRION'
    },

    // ============================================================
    // Mixing Parameters
    // ============================================================
    'AMIX': {
        description: 'Linear mixing parameter for charge density.',
        type: 'float',
        default: '0.4',
        category: 'Mixing',
        wikiPage: 'AMIX'
    },
    'BMIX': {
        description: 'Cutoff wave vector for Kerker mixing scheme.',
        type: 'float',
        default: '1.0',
        category: 'Mixing',
        wikiPage: 'BMIX'
    },
    'AMIX_MAG': {
        description: 'Linear mixing parameter for magnetization density.',
        type: 'float',
        default: '1.6',
        category: 'Mixing',
        wikiPage: 'AMIX_MAG'
    },
    'AMIN': {
        description: 'Minimal mixing parameter.',
        type: 'float',
        default: '0.1',
        category: 'Mixing',
        wikiPage: 'AMIN'
    },
    'IMIX': {
        description: 'Mixing type: 0=no, 1=Kerker, 2=Tchebychev, 4=Broyden.',
        type: 'integer',
        default: '4',
        values: ['0', '1', '2', '4'],
        category: 'Mixing',
        wikiPage: 'IMIX'
    },
    'INIMIX': {
        description: 'Initial mixing type.',
        type: 'integer',
        default: '1',
        category: 'Mixing',
        wikiPage: 'INIMIX'
    },
    'MIXPRE': {
        description: 'Preconditioning for Broyden mixer.',
        type: 'integer',
        default: '1',
        category: 'Mixing',
        wikiPage: 'MIXPRE'
    },
    'MAXMIX': {
        description: 'Maximum number of steps stored in Broyden mixer.',
        type: 'integer',
        default: '-45',
        category: 'Mixing',
        wikiPage: 'MAXMIX'
    },
    'WC': {
        description: 'Weight factor for each step in Broyden mixing.',
        type: 'float',
        default: '1000',
        category: 'Mixing',
        wikiPage: 'WC'
    },

    // ============================================================
    // Algorithm Parameters
    // ============================================================
    'IALGO': {
        description: 'Algorithm for electronic minimization (deprecated, use ALGO).',
        type: 'integer',
        default: '38',
        category: 'Algorithm',
        wikiPage: 'IALGO'
    },
    'LDIAG': {
        description: 'Sub-space diagonalization in eigenvalue solver.',
        type: 'boolean',
        default: '.TRUE.',
        category: 'Algorithm',
        wikiPage: 'LDIAG'
    },
    'LSUBROT': {
        description: 'Optimize orbitals by subspace rotation.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Algorithm',
        wikiPage: 'LSUBROT'
    },
    'WEIMIN': {
        description: 'Maximum weight for eigenvalue minimization.',
        type: 'float',
        default: '0.001',
        category: 'Algorithm',
        wikiPage: 'WEIMIN'
    },
    'EBREAK': {
        description: 'Break condition for eigenvalue optimization.',
        type: 'float',
        category: 'Algorithm',
        wikiPage: 'EBREAK'
    },
    'DEPER': {
        description: 'Relative energy change for break condition.',
        type: 'float',
        default: '0.3',
        category: 'Algorithm',
        wikiPage: 'DEPER'
    },
    'ENINI': {
        description: 'Initial cutoff energy for wavefunctions.',
        type: 'float',
        category: 'Algorithm',
        wikiPage: 'ENINI'
    },

    // ============================================================
    // Grid Parameters
    // ============================================================
    'NGXF': {
        description: 'FFT grid points in x for augmentation charges.',
        type: 'integer',
        category: 'Grid',
        wikiPage: 'NGXF'
    },
    'NGYF': {
        description: 'FFT grid points in y for augmentation charges.',
        type: 'integer',
        category: 'Grid',
        wikiPage: 'NGYF'
    },
    'NGZF': {
        description: 'FFT grid points in z for augmentation charges.',
        type: 'integer',
        category: 'Grid',
        wikiPage: 'NGZF'
    },
    'NGX': {
        description: 'FFT grid points in x direction.',
        type: 'integer',
        category: 'Grid',
        wikiPage: 'NGX'
    },
    'NGY': {
        description: 'FFT grid points in y direction.',
        type: 'integer',
        category: 'Grid',
        wikiPage: 'NGY'
    },
    'NGZ': {
        description: 'FFT grid points in z direction.',
        type: 'integer',
        category: 'Grid',
        wikiPage: 'NGZ'
    },

    // ============================================================
    // Additional Common Parameters
    // ============================================================
    'APACO': {
        description: 'Distance for pair correlation function (Å).',
        type: 'float',
        default: '16.0',
        category: 'MD',
        wikiPage: 'APACO'
    },
    'NPACO': {
        description: 'Number of slots for pair correlation function.',
        type: 'integer',
        default: '256',
        category: 'MD',
        wikiPage: 'NPACO'
    },
    'NLSPLINE': {
        description: 'Spline interpolation for projectors.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Basic',
        wikiPage: 'NLSPLINE'
    },
    'LCOMPAT': {
        description: 'Compatibility mode with VASP 4.4.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Basic',
        wikiPage: 'LCOMPAT'
    },
    'GGA_COMPAT': {
        description: 'GGA compatibility with VASP 4.4-4.6.',
        type: 'boolean',
        default: '.TRUE.',
        category: 'Basic',
        wikiPage: 'GGA_COMPAT'
    },
    'LMAXPAW': {
        description: 'Maximum l for on-site density.',
        type: 'integer',
        default: '-100 (auto)',
        category: 'PAW',
        wikiPage: 'LMAXPAW'
    },
    'ROPT': {
        description: 'Real-space optimization parameters per species.',
        type: 'float array',
        category: 'Basic',
        wikiPage: 'ROPT'
    },
    'VOSKOWN': {
        description: 'Vosko-Wilk-Nusair interpolation: 0=Perdew, 1=VWN.',
        type: 'integer',
        default: '0',
        values: ['0', '1'],
        category: 'XC',
        wikiPage: 'VOSKOWN'
    },

    // ============================================================
    // POTCAR / Pseudopotential Parameters
    // ============================================================
    'POMASS': {
        description: 'Ionic mass for each species from POTCAR (in atomic mass units).',
        type: 'float array',
        category: 'POTCAR',
        wikiPage: 'POMASS'
    },
    'ZVAL': {
        description: 'Number of valence electrons for each species from POTCAR.',
        type: 'float array',
        category: 'POTCAR',
        wikiPage: 'ZVAL'
    },
    'RWIGS': {
        description: 'Wigner-Seitz radius for each species (Å). Used for LORBIT projected DOS.',
        type: 'float array',
        category: 'POTCAR',
        wikiPage: 'RWIGS'
    },
    'EATOM': {
        description: 'Atomic reference energy from POTCAR (eV).',
        type: 'float',
        category: 'POTCAR',
        wikiPage: 'POTCAR'
    },
    'RMAX': {
        description: 'Maximum radius for radial grids in POTCAR (a.u.).',
        type: 'float',
        category: 'POTCAR',
        wikiPage: 'POTCAR'
    },
    'RAUG': {
        description: 'Augmentation sphere radius from POTCAR (a.u.).',
        type: 'float',
        category: 'POTCAR',
        wikiPage: 'POTCAR'
    },
    'RDEP': {
        description: 'Core radius for depletion charge from POTCAR (a.u.).',
        type: 'float',
        category: 'POTCAR',
        wikiPage: 'POTCAR'
    },
    'RDEPT': {
        description: 'Core radius for depletion charge kinetic energy (a.u.).',
        type: 'float',
        category: 'POTCAR',
        wikiPage: 'POTCAR'
    },
    'QCUT': {
        description: 'Plane wave cutoff for POTCAR projectors.',
        type: 'float',
        category: 'POTCAR',
        wikiPage: 'POTCAR'
    },
    'QGAM': {
        description: 'Gamma for POTCAR optimization.',
        type: 'float',
        category: 'POTCAR',
        wikiPage: 'POTCAR'
    },
    'DEXC': {
        description: 'Exchange-correlation energy difference from POTCAR.',
        type: 'float',
        category: 'POTCAR',
        wikiPage: 'POTCAR'
    },
    'EREF': {
        description: 'Reference energy from POTCAR.',
        type: 'float',
        category: 'POTCAR',
        wikiPage: 'POTCAR'
    },
    'LEXCH': {
        description: 'Exchange-correlation type in POTCAR (e.g., PE=PBE, CA=LDA).',
        type: 'string',
        category: 'POTCAR',
        wikiPage: 'POTCAR'
    },

    // ============================================================
    // Hybrid Functional Additional Parameters
    // ============================================================
    'HFSCREENC': {
        description: 'Screening for correlation in hybrid functionals.',
        type: 'float',
        default: '0',
        category: 'Hybrid',
        wikiPage: 'HFSCREEN'
    },
    'HFRCUT': {
        description: 'Real-space cutoff for HF exchange (Å).',
        type: 'float',
        default: '0',
        category: 'Hybrid',
        wikiPage: 'HFRCUT'
    },
    'HFALPHA': {
        description: 'Yukawa potential screening for exchange.',
        type: 'float',
        default: '0',
        category: 'Hybrid',
        wikiPage: 'HFALPHA'
    },
    'HFKIDENT': {
        description: 'Use k-point identity for HF.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Hybrid',
        wikiPage: 'HFKIDENT'
    },
    'LRHFCALC': {
        description: 'Long-range Hartree-Fock calculation.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Hybrid',
        wikiPage: 'LRHFCALC'
    },
    'LHARTREE': {
        description: 'Include Hartree energy in HF.',
        type: 'boolean',
        default: '.TRUE.',
        category: 'Hybrid',
        wikiPage: 'LHFCALC'
    },
    'ENCUTFOCK': {
        description: 'FFT grid cutoff for exact exchange (eV).',
        type: 'float',
        category: 'Hybrid',
        wikiPage: 'ENCUTFOCK'
    },
    'LMAXFOCK': {
        description: 'Maximum L for charge augmentation in HF.',
        type: 'integer',
        category: 'Hybrid',
        wikiPage: 'LMAXFOCK'
    },
    'LMAXFOCKAE': {
        description: 'Maximum L for AE charge augmentation in HF.',
        type: 'integer',
        category: 'Hybrid',
        wikiPage: 'LMAXFOCKAE'
    },
    'NMAXFOCKAE': {
        description: 'Maximum number of AE augmentation channels.',
        type: 'integer',
        category: 'Hybrid',
        wikiPage: 'NMAXFOCKAE'
    },
    'ALDAX': {
        description: 'Fraction of LDA exchange (usually 0 for hybrids).',
        type: 'float',
        default: '0',
        category: 'Hybrid',
        wikiPage: 'ALDAX'
    },
    'EXXOEP': {
        description: 'Exact exchange OEP method: 0=standard, 1=local, 2=mixed.',
        type: 'integer',
        default: '0',
        values: ['0', '1', '2'],
        category: 'Hybrid',
        wikiPage: 'EXXOEP'
    },
    'LMODELHF': {
        description: 'Use model HF (Thomas-Fermi screening).',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Hybrid',
        wikiPage: 'LMODELHF'
    },

    // ============================================================
    // GW and BSE Additional Parameters
    // ============================================================
    'NBANDSGW': {
        description: 'Number of bands for GW calculation.',
        type: 'integer',
        category: 'GW/BSE',
        wikiPage: 'NBANDSGW'
    },
    'NBANDSGWLOW': {
        description: 'Lowest band for GW calculation.',
        type: 'integer',
        category: 'GW/BSE',
        wikiPage: 'NBANDSGWLOW'
    },
    'ENCUTGWSOFT': {
        description: 'Soft cutoff for response function (eV).',
        type: 'float',
        category: 'GW/BSE',
        wikiPage: 'ENCUTGWSOFT'
    },
    'ENCUTLF': {
        description: 'Energy cutoff for local field effects (eV).',
        type: 'float',
        category: 'GW/BSE',
        wikiPage: 'ENCUTLF'
    },
    'ENCUT4O': {
        description: 'Cutoff for 4-orbital integrals (eV).',
        type: 'float',
        category: 'GW/BSE',
        wikiPage: 'ENCUT4O'
    },
    'NOMEGAR': {
        description: 'Number of frequency points for real-axis.',
        type: 'integer',
        category: 'GW/BSE',
        wikiPage: 'NOMEGAR'
    },
    'OMEGAMIN': {
        description: 'Minimum frequency for dielectric function (eV).',
        type: 'float',
        category: 'GW/BSE',
        wikiPage: 'OMEGAMIN'
    },
    'OMEGAGRID': {
        description: 'Type of frequency grid for GW.',
        type: 'integer',
        category: 'GW/BSE',
        wikiPage: 'OMEGAGRID'
    },
    'OMEGATL': {
        description: 'Maximum frequency for time-dependent calculations (eV).',
        type: 'float',
        category: 'GW/BSE',
        wikiPage: 'OMEGATL'
    },
    'OMEGAPAR': {
        description: 'Parallelization over frequency points.',
        type: 'integer',
        category: 'GW/BSE',
        wikiPage: 'OMEGAPAR'
    },
    'LFERMIGW': {
        description: 'Update Fermi energy in GW.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'GW/BSE',
        wikiPage: 'LFERMIGW'
    },
    'LSPECTRAL': {
        description: 'Use spectral method for dielectric function.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'GW/BSE',
        wikiPage: 'LSPECTRAL'
    },
    'LSPECTRALGW': {
        description: 'Calculate spectral function in GW.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'GW/BSE',
        wikiPage: 'LSPECTRALGW'
    },
    'SELFENERGY': {
        description: 'Self-energy evaluation method.',
        type: 'integer',
        category: 'GW/BSE',
        wikiPage: 'SELFENERGY'
    },
    'NATURALO': {
        description: 'Natural orbital occupations.',
        type: 'integer',
        category: 'GW/BSE',
        wikiPage: 'NATURALO'
    },
    'LSINGLES': {
        description: 'Include single excitations.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'GW/BSE',
        wikiPage: 'LSINGLES'
    },
    'LTRIPLET': {
        description: 'Calculate triplet excitations in BSE.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'GW/BSE',
        wikiPage: 'LTRIPLET'
    },
    'LADDER': {
        description: 'Include ladder diagrams.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'GW/BSE',
        wikiPage: 'LADDER'
    },
    'LFXC': {
        description: 'Include exchange-correlation kernel.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'GW/BSE',
        wikiPage: 'LFXC'
    },
    'LFXCEPS': {
        description: 'Include fxc in dielectric function.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'GW/BSE',
        wikiPage: 'LFXCEPS'
    },
    'LFXHEG': {
        description: 'Use homogeneous electron gas fxc.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'GW/BSE',
        wikiPage: 'LFXHEG'
    },
    'IBSE': {
        description: 'BSE type: 0=direct, 1=iterative, 2=Haydock.',
        type: 'integer',
        default: '0',
        values: ['0', '1', '2'],
        category: 'GW/BSE',
        wikiPage: 'IBSE'
    },
    'LRSRPA': {
        description: 'Long-range RPA correlation.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'GW/BSE',
        wikiPage: 'LRSRPA'
    },
    'LRSCOR': {
        description: 'Long-range correlation.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'GW/BSE',
        wikiPage: 'LRSCOR'
    },
    'NELMHF': {
        description: 'Number of electronic steps for HF.',
        type: 'integer',
        default: '1',
        category: 'GW/BSE',
        wikiPage: 'NELMHF'
    },
    'SCISSOR': {
        description: 'Scissors operator shift (eV).',
        type: 'float',
        default: '0',
        category: 'GW/BSE',
        wikiPage: 'SCISSOR'
    },
    'L2ORDER': {
        description: 'Second order contribution in RPA.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'GW/BSE',
        wikiPage: 'L2ORDER'
    },

    // ============================================================
    // Van der Waals Additional Parameters
    // ============================================================
    'BPARAM': {
        description: 'B parameter for DFT-D corrections.',
        type: 'float',
        category: 'Van der Waals',
        wikiPage: 'DFT-D2'
    },
    'CPARAM': {
        description: 'C parameter for DFT-D corrections.',
        type: 'float',
        category: 'Van der Waals',
        wikiPage: 'DFT-D2'
    },
    'LVDWSCS': {
        description: 'Self-consistent screening for vdW.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Van der Waals',
        wikiPage: 'LVDWSCS'
    },

    // ============================================================
    // Berry Phase and Electric Field
    // ============================================================
    'LBERRY': {
        description: 'Calculate Berry phase for polarization.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Berry Phase',
        wikiPage: 'LBERRY'
    },
    'IGPAR': {
        description: 'Direction for Berry phase: 1=a, 2=b, 3=c.',
        type: 'integer',
        values: ['1', '2', '3'],
        category: 'Berry Phase',
        wikiPage: 'IGPAR'
    },
    'NPPSTR': {
        description: 'Number of k-points along Berry phase direction.',
        type: 'integer',
        category: 'Berry Phase',
        wikiPage: 'NPPSTR'
    },
    'LCALCEPS': {
        description: 'Calculate dielectric tensor.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Berry Phase',
        wikiPage: 'LCALCEPS'
    },
    'LCALCPOL': {
        description: 'Calculate polarization.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Berry Phase',
        wikiPage: 'LCALCPOL'
    },
    'LPEAD': {
        description: 'Derivative of orbitals w.r.t. k for optical properties.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Berry Phase',
        wikiPage: 'LPEAD'
    },

    // ============================================================
    // Solvation and Implicit Solvent
    // ============================================================
    'LSOL': {
        description: 'Enable implicit solvation model (VASPsol).',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Solvation',
        wikiPage: 'LSOL'
    },
    'EB_K': {
        description: 'Bulk dielectric constant of solvent.',
        type: 'float',
        default: '78.4',
        category: 'Solvation',
        wikiPage: 'EB_K'
    },
    'TAU': {
        description: 'Surface tension parameter for solvation (eV/Å²).',
        type: 'float',
        category: 'Solvation',
        wikiPage: 'TAU'
    },
    'LAMBDA': {
        description: 'Debye screening length for solvation.',
        type: 'float',
        category: 'Solvation',
        wikiPage: 'LAMBDA'
    },

    // ============================================================
    // Core Level and Spectroscopy
    // ============================================================
    'ICORELEVEL': {
        description: 'Core level shift calculation: 0=none, 1=initial, 2=final.',
        type: 'integer',
        default: '0',
        values: ['0', '1', '2'],
        category: 'Core Level',
        wikiPage: 'ICORELEVEL'
    },
    'CLNT': {
        description: 'Principal quantum number for core level.',
        type: 'integer',
        category: 'Core Level',
        wikiPage: 'CLNT'
    },
    'CLN': {
        description: 'Principal quantum number n for core hole.',
        type: 'integer',
        category: 'Core Level',
        wikiPage: 'CLN'
    },
    'CLL': {
        description: 'Angular momentum l for core hole.',
        type: 'integer',
        category: 'Core Level',
        wikiPage: 'CLL'
    },
    'CLZ': {
        description: 'Effective nuclear charge for core level.',
        type: 'float',
        category: 'Core Level',
        wikiPage: 'CLZ'
    },

    // ============================================================
    // NMR and Magnetic Response
    // ============================================================
    'LCHIMAG': {
        description: 'Calculate magnetic susceptibility.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'NMR',
        wikiPage: 'LCHIMAG'
    },
    'LNMR_SYM_RED': {
        description: 'Reduce symmetry for NMR calculations.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'NMR',
        wikiPage: 'LNMR_SYM_RED'
    },
    'ORBITALMAG': {
        description: 'Calculate orbital magnetization.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'NMR',
        wikiPage: 'ORBITALMAG'
    },
    'LNABLA': {
        description: 'Use nabla operator for momentum matrix.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'NMR',
        wikiPage: 'LNABLA'
    },
    'LEFG': {
        description: 'Calculate electric field gradient.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'NMR',
        wikiPage: 'LEFG'
    },
    'LLRAUG': {
        description: 'Include augmentation for linear response.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'NMR',
        wikiPage: 'LLRAUG'
    },
    'LRHOB': {
        description: 'Include B-field in density.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'NMR',
        wikiPage: 'LRHOB'
    },

    // ============================================================
    // Spiral and Non-Collinear Magnetism
    // ============================================================
    'LSPIRAL': {
        description: 'Generalized Bloch spin spiral calculation.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Spin Spiral',
        wikiPage: 'LSPIRAL'
    },
    'QSPIRAL': {
        description: 'Spin spiral propagation vector.',
        type: 'float array',
        category: 'Spin Spiral',
        wikiPage: 'QSPIRAL'
    },
    'LZEROZ': {
        description: 'Force zero z-component of magnetization.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Spin Spiral',
        wikiPage: 'LZEROZ'
    },
    'LMAGBLOCH': {
        description: 'Calculate Bloch representation of magnetization.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Spin Spiral',
        wikiPage: 'LMAGBLOCH'
    },

    // ============================================================
    // Parallelization Additional Parameters
    // ============================================================
    'NKREDX': {
        description: 'Reduce k-points for HF in x direction.',
        type: 'integer',
        default: '1',
        category: 'Parallel',
        wikiPage: 'NKREDX'
    },
    'NKREDY': {
        description: 'Reduce k-points for HF in y direction.',
        type: 'integer',
        default: '1',
        category: 'Parallel',
        wikiPage: 'NKREDY'
    },
    'NKREDZ': {
        description: 'Reduce k-points for HF in z direction.',
        type: 'integer',
        default: '1',
        category: 'Parallel',
        wikiPage: 'NKREDZ'
    },
    'NKREDLFX': {
        description: 'Reduce k-points for local field effects.',
        type: 'integer',
        default: '1',
        category: 'Parallel',
        wikiPage: 'NKREDLFX'
    },

    // ============================================================
    // Advanced Algorithm Parameters
    // ============================================================
    'LBFGS': {
        description: 'Use L-BFGS optimizer (requires IBRION=3).',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Algorithm',
        wikiPage: 'IBRION'
    },
    'LHF': {
        description: 'Enable Hartree-Fock (same as LHFCALC).',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Hybrid',
        wikiPage: 'LHFCALC'
    },
    'LIBXC': {
        description: 'Use LibXC library for XC functionals.',
        type: 'integer',
        category: 'XC',
        wikiPage: 'LIBXC'
    },
    'ICHIBARE': {
        description: 'Bare susceptibility calculation type.',
        type: 'integer',
        category: 'Response',
        wikiPage: 'ICHIBARE'
    },
    'SMIX': {
        description: 'Mixing parameter for Kerker scheme.',
        type: 'float',
        default: '0.4',
        category: 'Mixing',
        wikiPage: 'SMIX'
    },
    'IWAVPR': {
        description: 'Wavefunction prediction: 0=none, 1=charge, 2=wave, 3=both.',
        type: 'integer',
        default: '2',
        values: ['0', '1', '2', '3'],
        category: 'Algorithm',
        wikiPage: 'IWAVPR'
    },
    'IRESTART': {
        description: 'Restart mode for MD.',
        type: 'integer',
        default: '0',
        category: 'MD',
        wikiPage: 'IRESTART'
    },
    'NREBOOT': {
        description: 'Number of reboots for MD.',
        type: 'integer',
        default: '0',
        category: 'MD',
        wikiPage: 'NREBOOT'
    },
    'SCALEE': {
        description: 'Scaling for kinetic energy in MD.',
        type: 'float',
        default: '1.0',
        category: 'MD',
        wikiPage: 'SCALEE'
    },
    'LMONO': {
        description: 'Monopole corrections.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Dipole',
        wikiPage: 'LMONO'
    },
    'LASYNC': {
        description: 'Asynchronous communication.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Parallel',
        wikiPage: 'LASYNC'
    },
    'LINTERFAST': {
        description: 'Fast interpolation.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Algorithm',
        wikiPage: 'LINTERFAST'
    },

    // ============================================================
    // DFPT and Linear Response Additional
    // ============================================================
    'LTCTE': {
        description: 'Thermal conductivity (electron).',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Transport',
        wikiPage: 'LTCTE'
    },
    'LTETE': {
        description: 'Thermal transport (electron-electron).',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Transport',
        wikiPage: 'LTETE'
    },
    'LTHOMAS': {
        description: 'Thomas-Fermi screening.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Screening',
        wikiPage: 'LTHOMAS'
    },
    'LTDEP': {
        description: 'Temperature-dependent effective potential.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Phonons',
        wikiPage: 'LTDEP'
    },
    'LSYMGRAD': {
        description: 'Symmetrize gradient.',
        type: 'boolean',
        default: '.TRUE.',
        category: 'Symmetry',
        wikiPage: 'LSYMGRAD'
    },
    'LFOCKAEDFT': {
        description: 'Include Fock term in AE DFT.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Hybrid',
        wikiPage: 'LFOCKAEDFT'
    },
    'LVEL': {
        description: 'Velocity operator for optical properties.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Optics',
        wikiPage: 'LVEL'
    },
    'LDNEB': {
        description: 'Dynamical NEB.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'NEB',
        wikiPage: 'LDNEB'
    },
    'LUSEW': {
        description: 'Use wavefunctions for transition state.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'NEB',
        wikiPage: 'LUSEW'
    },
    'LUSEVDW': {
        description: 'Use vdW-corrected energies.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Van der Waals',
        wikiPage: 'LUSEVDW'
    },
    'LCORR': {
        description: 'Harris-Foulkes correction to forces.',
        type: 'boolean',
        default: '.TRUE.',
        category: 'Basic',
        wikiPage: 'LCORR'
    },

    // ============================================================
    // Output Energy Terms
    // ============================================================
    'TOTEN': {
        description: 'Total free energy including entropy term (eV).',
        type: 'float',
        category: 'Output',
        wikiPage: 'TOTEN'
    },
    'E0': {
        description: 'Energy sigma->0 extrapolated (eV).',
        type: 'float',
        category: 'Output',
        wikiPage: 'TOTEN'
    },
    'EBANDS': {
        description: 'Band structure energy (eV).',
        type: 'float',
        category: 'Output',
        wikiPage: 'OUTCAR'
    },
    'EFERMI': {
        description: 'Fermi energy (eV).',
        type: 'float',
        category: 'Output',
        wikiPage: 'EFERMI'
    },

    // ============================================================
    // System Information
    // ============================================================
    'NKPTS': {
        description: 'Total number of k-points.',
        type: 'integer',
        category: 'K-points',
        wikiPage: 'KPOINTS'
    },
    'NIONS': {
        description: 'Total number of ions (atoms).',
        type: 'integer',
        category: 'System',
        wikiPage: 'POSCAR'
    },
    'NPLWV': {
        description: 'Total number of plane waves.',
        type: 'integer',
        category: 'System',
        wikiPage: 'OUTCAR'
    },
    'ALAT': {
        description: 'Lattice constant (Å).',
        type: 'float',
        category: 'System',
        wikiPage: 'POSCAR'
    },
    'A1': {
        description: 'First lattice vector.',
        type: 'float array',
        category: 'System',
        wikiPage: 'POSCAR'
    },
    'A2': {
        description: 'Second lattice vector.',
        type: 'float array',
        category: 'System',
        wikiPage: 'POSCAR'
    },
    'A3': {
        description: 'Third lattice vector.',
        type: 'float array',
        category: 'System',
        wikiPage: 'POSCAR'
    },
    'DIM': {
        description: 'Dimensionality of the system.',
        type: 'integer',
        category: 'System',
        wikiPage: 'OUTCAR'
    },
    'KPOINT': {
        description: 'Current k-point index.',
        type: 'integer',
        category: 'K-points',
        wikiPage: 'KPOINTS'
    },
    'KINTER': {
        description: 'K-point interpolation scheme.',
        type: 'integer',
        category: 'K-points',
        wikiPage: 'KINTER'
    },
    'NQ': {
        description: 'Number of q-points for response.',
        type: 'integer',
        category: 'Response',
        wikiPage: 'OUTCAR'
    },
    'QGRID': {
        description: 'Q-point grid specification.',
        type: 'integer array',
        category: 'Response',
        wikiPage: 'QGRID'
    },
    'IXMIN': {
        description: 'Minimum FFT index in x.',
        type: 'integer',
        category: 'Grid',
        wikiPage: 'OUTCAR'
    },
    'IXMAX': {
        description: 'Maximum FFT index in x.',
        type: 'integer',
        category: 'Grid',
        wikiPage: 'OUTCAR'
    },
    'SHIFTRED': {
        description: 'Shift reduction for band structures.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'K-points',
        wikiPage: 'SHIFTRED'
    },
    'EVENONLY': {
        description: 'Use only even k-points.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'K-points',
        wikiPage: 'EVENONLY'
    },
    'ODDONLY': {
        description: 'Use only odd k-points.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'K-points',
        wikiPage: 'ODDONLY'
    },
    'EVENONLYGW': {
        description: 'Use only even k-points for GW.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'GW/BSE',
        wikiPage: 'EVENONLYGW'
    },
    'ODDONLYGW': {
        description: 'Use only odd k-points for GW.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'GW/BSE',
        wikiPage: 'ODDONLYGW'
    },
    'LMAXMP2': {
        description: 'Maximum L for MP2 calculations.',
        type: 'integer',
        category: 'GW/BSE',
        wikiPage: 'LMAXMP2'
    },

    // ============================================================
    // MD Additional Parameters
    // ============================================================
    'TAUPAR': {
        description: 'Parrinello-Rahman barostat time constant.',
        type: 'float',
        category: 'MD',
        wikiPage: 'TAUPAR'
    },
    'TEIN': {
        description: 'Initial kinetic energy for MD.',
        type: 'float',
        category: 'MD',
        wikiPage: 'TEIN'
    },
    'TELESCOPE': {
        description: 'Telescope sampling for response.',
        type: 'integer',
        category: 'Response',
        wikiPage: 'TELESCOPE'
    },
    'TURBO': {
        description: 'Turbo mode for faster calculations.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Algorithm',
        wikiPage: 'TURBO'
    },
    'VCA': {
        description: 'Virtual crystal approximation weights.',
        type: 'float array',
        category: 'VCA',
        wikiPage: 'VCA'
    },
    'RTIME': {
        description: 'Real-time TDDFT time step.',
        type: 'float',
        category: 'TDDFT',
        wikiPage: 'RTIME'
    },
    'MAXITERFT': {
        description: 'Maximum iterations for response function.',
        type: 'integer',
        category: 'Response',
        wikiPage: 'MAXITERFT'
    },
    'MCALPHA': {
        description: 'Alpha parameter for meta-GGA.',
        type: 'float',
        category: 'XC',
        wikiPage: 'MCALPHA'
    },
    'DEG_THRESHOLD': {
        description: 'Threshold for degeneracy detection.',
        type: 'float',
        category: 'Algorithm',
        wikiPage: 'DEG_THRESHOLD'
    },
};

// OUTCAR section patterns for navigation
const OUTCAR_SECTIONS = [
    { pattern: /^\s*Startparameter for this run/i, name: 'Startparameters' },
    { pattern: /^\s*Dimension of arrays/i, name: 'Dimension of arrays' },
    { pattern: /^\s*POSCAR:/i, name: 'POSCAR' },
    { pattern: /^\s*k-points\s+in\s+BZ/i, name: 'K-points' },
    { pattern: /^\s*Electronic Relaxation/i, name: 'Electronic Relaxation' },
    { pattern: /^\s*Ionic Relaxation/i, name: 'Ionic Relaxation' },
    { pattern: /^\s*DOS related values/i, name: 'DOS Parameters' },
    { pattern: /^\s*Iteration\s+(\d+)\s*\(/i, name: 'Iteration $1' },
    { pattern: /^\s*aborting loop/i, name: 'Loop Aborted' },
    { pattern: /^\s*FREE ENERGIE OF THE ION-ELECTRON SYSTEM/i, name: 'Free Energy' },
    { pattern: /^\s*VOLUME and BASIS/i, name: 'Volume and Basis' },
    { pattern: /^\s*TOTAL-FORCE/i, name: 'Total Force' },
    { pattern: /^\s*STRESS TENSOR/i, name: 'Stress Tensor' },
    { pattern: /^\s*EIGENVALUE/i, name: 'Eigenvalues' },
    { pattern: /^\s*E-fermi\s*:/i, name: 'Fermi Energy' },
    { pattern: /^\s*BRION:\s*g\(F\)/i, name: 'Ionic Convergence' },
    { pattern: /^\s*reached required accuracy/i, name: 'Converged!' },
    { pattern: /^\s*General timing and accounting/i, name: 'Timing Summary' },
    { pattern: /^\s*writing wavefunctions/i, name: 'Writing Wavefunctions' },
];

// Document symbol provider for OUTCAR navigation
class OutcarSymbolProvider implements vscode.DocumentSymbolProvider {
    public provideDocumentSymbols(
        document: vscode.TextDocument,
        _token: vscode.CancellationToken
    ): vscode.ProviderResult<vscode.DocumentSymbol[]> {
        const symbols: vscode.DocumentSymbol[] = [];
        const text = document.getText();
        const lines = text.split('\n');

        for (let i = 0; i < lines.length; i++) {
            const line = lines[i];

            for (const section of OUTCAR_SECTIONS) {
                const match = line.match(section.pattern);
                if (match) {
                    let name = section.name;
                    // Replace $1, $2, etc. with captured groups
                    if (match[1]) {
                        name = name.replace('$1', match[1]);
                    }

                    const range = new vscode.Range(i, 0, i, line.length);
                    const symbol = new vscode.DocumentSymbol(
                        name,
                        '',
                        vscode.SymbolKind.Function,
                        range,
                        range
                    );
                    symbols.push(symbol);
                    break; // Only match one pattern per line
                }
            }
        }

        return symbols;
    }
}

// CodeLens provider for VASP files
class VaspCodeLensProvider implements vscode.CodeLensProvider {
    public provideCodeLenses(
        document: vscode.TextDocument,
        _token: vscode.CancellationToken
    ): vscode.ProviderResult<vscode.CodeLens[]> {
        const codeLenses: vscode.CodeLens[] = [];
        const topRange = new vscode.Range(0, 0, 0, 0);
        const filePath = document.uri.fsPath;
        const dir = path.dirname(filePath);

        // Build description based on what's available
        const parts: string[] = [];

        // Check for POSCAR/CONTCAR
        const poscarPath = path.join(dir, 'POSCAR');
        const contcarPath = path.join(dir, 'CONTCAR');
        if (fs.existsSync(poscarPath) || fs.existsSync(contcarPath)) {
            parts.push('Structure');
        }

        // Check for vasprun.xml
        const vasprunPath = path.join(dir, 'vasprun.xml');
        if (fs.existsSync(vasprunPath)) {
            parts.push('Summary');
            // Check for relaxation
            try {
                const content = fs.readFileSync(vasprunPath, 'utf-8');
                const calcBlocks = content.match(/<calculation>/g);
                if (calcBlocks && calcBlocks.length > 1) {
                    parts.push(calcBlocks.length + ' steps');
                }
            } catch (e) {
                // Ignore errors
            }
        }

        if (parts.length > 0) {
            codeLenses.push(new vscode.CodeLens(topRange, {
                title: '$(beaker) VASPsum: ' + parts.join(' + '),
                command: 'vasp.vaspsum',
                tooltip: 'Show VASP calculation summary and structure'
            }));
        }

        return codeLenses;
    }
}

// Create hover provider for INCAR files
class VaspHoverProvider implements vscode.HoverProvider {
    public provideHover(
        document: vscode.TextDocument,
        position: vscode.Position,
        _token: vscode.CancellationToken
    ): vscode.ProviderResult<vscode.Hover> {
        const range = document.getWordRangeAtPosition(position, /[A-Za-z_][A-Za-z0-9_]*/);
        if (!range) {
            return null;
        }

        const word = document.getText(range).toUpperCase();
        const param = VASP_PARAMETERS[word];

        if (!param) {
            return null;
        }

        const markdown = new vscode.MarkdownString();
        markdown.isTrusted = true;

        // Parameter name and category
        markdown.appendMarkdown(`## ${word}\n\n`);
        markdown.appendMarkdown(`**Category:** ${param.category}\n\n`);

        // Description
        markdown.appendMarkdown(`${param.description}\n\n`);

        // Type and default
        markdown.appendMarkdown(`**Type:** \`${param.type}\`\n\n`);
        if (param.default) {
            markdown.appendMarkdown(`**Default:** \`${param.default}\`\n\n`);
        }

        // Possible values
        if (param.values && param.values.length > 0) {
            markdown.appendMarkdown(`**Values:** ${param.values.map(v => `\`${v}\``).join(', ')}\n\n`);
        }

        // Link to VASP wiki
        const wikiUrl = `${VASP_WIKI_BASE}${param.wikiPage}`;
        markdown.appendMarkdown(`---\n\n`);
        markdown.appendMarkdown(`[Open VASP Wiki documentation](${wikiUrl})\n`);

        return new vscode.Hover(markdown, range);
    }
}

// Provide completions for INCAR parameters
class VaspCompletionProvider implements vscode.CompletionItemProvider {
    public provideCompletionItems(
        document: vscode.TextDocument,
        position: vscode.Position,
        _token: vscode.CancellationToken,
        _context: vscode.CompletionContext
    ): vscode.ProviderResult<vscode.CompletionItem[]> {
        const linePrefix = document.lineAt(position).text.substring(0, position.character);

        // Only complete at the start of a line or after whitespace
        if (!/^[\s]*$/.test(linePrefix) && !/[\s=]$/.test(linePrefix)) {
            // Check if we're typing a parameter name
            if (!/^[\s]*[A-Za-z_][A-Za-z0-9_]*$/.test(linePrefix)) {
                return [];
            }
        }

        const completions: vscode.CompletionItem[] = [];

        for (const [name, param] of Object.entries(VASP_PARAMETERS)) {
            const item = new vscode.CompletionItem(name, vscode.CompletionItemKind.Property);
            item.detail = `[${param.category}] ${param.type}`;
            item.documentation = new vscode.MarkdownString(
                `${param.description}\n\n` +
                (param.default ? `**Default:** \`${param.default}\`\n\n` : '') +
                `[VASP Wiki](${VASP_WIKI_BASE}${param.wikiPage})`
            );

            // Add snippet for boolean parameters
            if (param.type === 'boolean') {
                item.insertText = new vscode.SnippetString(`${name} = \${1|.TRUE.,.FALSE.|}`);
            } else if (param.values && param.values.length > 0) {
                const choices = param.values.join(',');
                item.insertText = new vscode.SnippetString(`${name} = \${1|${choices}|}`);
            } else {
                item.insertText = new vscode.SnippetString(`${name} = $1`);
            }

            completions.push(item);
        }

        return completions;
    }
}

// System Info Panel for debugging and support
// POSCAR/CONTCAR structure data interface
interface PoscarData {
    title: string;
    scale: number;
    lattice: number[][];  // 3x3
    elements: string[];
    counts: number[];
    isDirect: boolean;
    positions: number[][];  // Nx3
    selectiveDynamics: boolean;
    constraints: boolean[][];  // Nx3, true = constrained (F), false = free (T)
}

// CPK color scheme for common elements
const CPK_COLORS: Record<string, string> = {
    'H': '#FFFFFF', 'He': '#D9FFFF', 'Li': '#CC80FF', 'Be': '#C2FF00',
    'B': '#FFB5B5', 'C': '#909090', 'N': '#3050F8', 'O': '#FF0D0D',
    'F': '#90E050', 'Ne': '#B3E3F5', 'Na': '#AB5CF2', 'Mg': '#8AFF00',
    'Al': '#BFA6A6', 'Si': '#F0C8A0', 'P': '#FF8000', 'S': '#FFFF30',
    'Cl': '#1FF01F', 'Ar': '#80D1E3', 'K': '#8F40D4', 'Ca': '#3DFF00',
    'Sc': '#E6E6E6', 'Ti': '#BFC2C7', 'V': '#A6A6AB', 'Cr': '#8A99C7',
    'Mn': '#9C7AC7', 'Fe': '#E06633', 'Co': '#F090A0', 'Ni': '#50D050',
    'Cu': '#C88033', 'Zn': '#7D80B0', 'Ga': '#C28F8F', 'Ge': '#668F8F',
    'As': '#BD80E3', 'Se': '#FFA100', 'Br': '#A62929', 'Kr': '#5CB8D1',
    'Rb': '#702EB0', 'Sr': '#00FF00', 'Y': '#94FFFF', 'Zr': '#94E0E0',
    'Nb': '#73C2C9', 'Mo': '#54B5B5', 'Tc': '#3B9E9E', 'Ru': '#248F8F',
    'Rh': '#0A7D8C', 'Pd': '#006985', 'Ag': '#C0C0C0', 'Cd': '#FFD98F',
    'In': '#A67573', 'Sn': '#668080', 'Sb': '#9E63B5', 'Te': '#D47A00',
    'I': '#940094', 'Xe': '#429EB0', 'Cs': '#57178F', 'Ba': '#00C900',
    'La': '#70D4FF', 'Ce': '#FFFFC7', 'Pr': '#D9FFC7', 'Nd': '#C7FFC7',
    'Pm': '#A3FFC7', 'Sm': '#8FFFC7', 'Eu': '#61FFC7', 'Gd': '#45FFC7',
    'Tb': '#30FFC7', 'Dy': '#1FFFC7', 'Ho': '#00FF9C', 'Er': '#00E675',
    'Tm': '#00D452', 'Yb': '#00BF38', 'Lu': '#00AB24', 'Hf': '#4DC2FF',
    'Ta': '#4DA6FF', 'W': '#2194D6', 'Re': '#267DAB', 'Os': '#266696',
    'Ir': '#175487', 'Pt': '#D0D0E0', 'Au': '#FFD123', 'Hg': '#B8B8D0',
    'Tl': '#A6544D', 'Pb': '#575961', 'Bi': '#9E4FB5', 'Po': '#AB5C00',
    'At': '#754F45', 'Rn': '#428296', 'Fr': '#420066', 'Ra': '#007D00',
    'Ac': '#70ABFA', 'Th': '#00BAFF', 'Pa': '#00A1FF', 'U': '#008FFF',
    'Np': '#0080FF', 'Pu': '#006BFF', 'Am': '#545CF2', 'Cm': '#785CE3'
};

// Parse POSCAR/CONTCAR file content
function parsePoscar(content: string): PoscarData | null {
    const lines = content.split('\n').map(l => l.trim()).filter(l => l.length > 0 && !l.startsWith('#'));

    if (lines.length < 7) {
        return null;
    }

    let lineIndex = 0;

    // Line 1: Title/comment
    const title = lines[lineIndex++];

    // Line 2: Scale factor
    const scale = parseFloat(lines[lineIndex++]);
    if (isNaN(scale)) {
        return null;
    }

    // Lines 3-5: Lattice vectors
    const lattice: number[][] = [];
    for (let i = 0; i < 3; i++) {
        const parts = lines[lineIndex++].split(/\s+/).map(parseFloat);
        if (parts.length < 3 || parts.some(isNaN)) {
            return null;
        }
        lattice.push(parts.slice(0, 3));
    }

    // Line 6: Element symbols (VASP 5+) or atom counts (VASP 4)
    let elements: string[] = [];
    let counts: number[] = [];

    const line6Parts = lines[lineIndex].split(/\s+/);
    const isVasp5 = isNaN(parseInt(line6Parts[0]));

    if (isVasp5) {
        // VASP 5+ format: element symbols on this line
        elements = line6Parts.filter(p => p.match(/^[A-Z][a-z]?$/));
        lineIndex++;
        // Next line has counts
        counts = lines[lineIndex++].split(/\s+/).map(n => parseInt(n)).filter(n => !isNaN(n));
    } else {
        // VASP 4 format: counts directly, try to infer elements from title
        counts = line6Parts.map(n => parseInt(n)).filter(n => !isNaN(n));
        lineIndex++;
        // Try to extract elements from title
        const titleMatch = title.match(/([A-Z][a-z]?\d*)+/g);
        if (titleMatch) {
            elements = titleMatch.join('').match(/[A-Z][a-z]?/g) || [];
        }
        // If still no elements, use placeholder symbols
        if (elements.length < counts.length) {
            elements = counts.map((_, i) => `X${i + 1}`);
        }
    }

    // Ensure elements and counts match
    elements = elements.slice(0, counts.length);
    while (elements.length < counts.length) {
        elements.push(`X${elements.length + 1}`);
    }

    // Check for Selective dynamics line
    let selectiveDynamics = false;
    const nextLine = lines[lineIndex];
    if (nextLine && (nextLine.toLowerCase().startsWith('s') && !nextLine.toLowerCase().startsWith('d') && !nextLine.toLowerCase().startsWith('c'))) {
        selectiveDynamics = true;
        lineIndex++; // Skip selective dynamics line
    }

    // Coordinate type line (Direct/Cartesian)
    const coordLine = lines[lineIndex++];
    const isDirect = coordLine.toLowerCase().startsWith('d');

    // Atomic positions and constraints
    const totalAtoms = counts.reduce((a, b) => a + b, 0);
    const positions: number[][] = [];
    const constraints: boolean[][] = [];

    for (let i = 0; i < totalAtoms && lineIndex < lines.length; i++) {
        const lineParts = lines[lineIndex++].split(/\s+/);
        const coords = lineParts.slice(0, 3).map(parseFloat);
        if (coords.length >= 3 && !coords.some(isNaN)) {
            positions.push(coords);

            // Parse selective dynamics flags (T = free, F = constrained)
            if (selectiveDynamics && lineParts.length >= 6) {
                const flags = lineParts.slice(3, 6).map(f => f.toUpperCase() === 'F');
                constraints.push(flags);
            } else {
                // No constraints - all free
                constraints.push([false, false, false]);
            }
        }
    }

    if (positions.length !== totalAtoms) {
        return null;
    }

    return {
        title,
        scale,
        lattice,
        elements,
        counts,
        isDirect,
        positions,
        selectiveDynamics,
        constraints
    };
}

// Convert POSCAR data to Cartesian coordinates
function toCartesian(data: PoscarData): number[][] {
    if (!data.isDirect) {
        // Already Cartesian, just apply scale
        return data.positions.map(pos => pos.map(x => x * data.scale));
    }

    // Convert from fractional to Cartesian
    const scaledLattice = data.lattice.map(row => row.map(x => x * data.scale));

    return data.positions.map(frac => {
        const cart = [0, 0, 0];
        for (let i = 0; i < 3; i++) {
            for (let j = 0; j < 3; j++) {
                cart[j] += frac[i] * scaledLattice[i][j];
            }
        }
        return cart;
    });
}

// Calculate cell volume
function cellVolume(data: PoscarData): number {
    const a = data.lattice[0].map(x => x * data.scale);
    const b = data.lattice[1].map(x => x * data.scale);
    const c = data.lattice[2].map(x => x * data.scale);

    // Volume = a · (b × c)
    const cross = [
        b[1] * c[2] - b[2] * c[1],
        b[2] * c[0] - b[0] * c[2],
        b[0] * c[1] - b[1] * c[0]
    ];
    return Math.abs(a[0] * cross[0] + a[1] * cross[1] + a[2] * cross[2]);
}

// Generate chemical formula
function getFormula(data: PoscarData): string {
    return data.elements.map((el, i) => {
        const count = data.counts[i];
        return count === 1 ? el : `${el}${count}`;
    }).join('');
}

// Structure Preview Panel for POSCAR/CONTCAR visualization
class StructurePreviewPanel {
    public static currentPanel: StructurePreviewPanel | undefined;
    private readonly _panel: vscode.WebviewPanel;
    private _disposables: vscode.Disposable[] = [];
    private _currentUri: vscode.Uri | undefined;

    public static createOrShow(uri: vscode.Uri) {
        const column = vscode.ViewColumn.Beside;

        if (StructurePreviewPanel.currentPanel) {
            StructurePreviewPanel.currentPanel._panel.reveal(column);
            StructurePreviewPanel.currentPanel._update(uri);
            return;
        }

        const panel = vscode.window.createWebviewPanel(
            'vaspStructurePreview',
            'Structure: ' + path.basename(uri.fsPath),
            column,
            {
                enableScripts: true,
                retainContextWhenHidden: true
            }
        );

        StructurePreviewPanel.currentPanel = new StructurePreviewPanel(panel, uri);
    }

    private constructor(panel: vscode.WebviewPanel, uri: vscode.Uri) {
        this._panel = panel;
        this._currentUri = uri;
        this._update(uri);

        this._panel.onDidDispose(() => this.dispose(), null, this._disposables);
    }

    public dispose() {
        StructurePreviewPanel.currentPanel = undefined;
        this._panel.dispose();
        while (this._disposables.length) {
            const x = this._disposables.pop();
            if (x) {
                x.dispose();
            }
        }
    }

    private async _update(uri: vscode.Uri) {
        this._currentUri = uri;
        this._panel.title = 'Structure: ' + path.basename(uri.fsPath);

        try {
            const content = fs.readFileSync(uri.fsPath, 'utf-8');
            const data = parsePoscar(content);

            if (data) {
                this._panel.webview.html = this._getHtmlForWebview(data);
            } else {
                this._panel.webview.html = this._getErrorHtml('Failed to parse POSCAR/CONTCAR file');
            }
        } catch (error) {
            this._panel.webview.html = this._getErrorHtml(`Error reading file: ${error}`);
        }
    }

    private _getHtmlForWebview(data: PoscarData): string {
        const cartesian = toCartesian(data);
        const volume = cellVolume(data);
        const formula = getFormula(data);
        const totalAtoms = data.counts.reduce((a, b) => a + b, 0);

        // Build atom data for 3Dmol.js
        let atomIndex = 0;
        const atoms: { elem: string; x: number; y: number; z: number; color: string; constrained: boolean }[] = [];

        for (let i = 0; i < data.elements.length; i++) {
            const element = data.elements[i];
            const count = data.counts[i];
            const color = CPK_COLORS[element] || '#FF00FF';

            for (let j = 0; j < count; j++) {
                const pos = cartesian[atomIndex];
                const constraint = data.constraints[atomIndex] || [false, false, false];
                // Atom is constrained if ANY direction is constrained
                const isConstrained = constraint.some(c => c);
                atomIndex++;
                atoms.push({
                    elem: element,
                    x: pos[0],
                    y: pos[1],
                    z: pos[2],
                    color: color,
                    constrained: isConstrained
                });
            }
        }

        // Scaled lattice vectors for unit cell
        const scaledLattice = data.lattice.map(row => row.map(x => x * data.scale));

        return `<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <meta http-equiv="Content-Security-Policy" content="default-src 'none'; script-src https://3dmol.org 'unsafe-inline' 'unsafe-eval'; style-src 'unsafe-inline'; img-src https: data:;">
    <title>Structure Preview</title>
    <script src="https://3dmol.org/build/3Dmol-min.js"></script>
    <style>
        body {
            font-family: var(--vscode-font-family);
            font-size: var(--vscode-font-size);
            color: var(--vscode-foreground);
            background-color: var(--vscode-editor-background);
            margin: 0;
            padding: 10px;
            display: flex;
            flex-direction: column;
            height: 100vh;
            box-sizing: border-box;
            overflow-y: auto;
        }
        .info-panel {
            background-color: var(--vscode-editor-inactiveSelectionBackground);
            padding: 10px 15px;
            border-radius: 5px;
            margin-bottom: 10px;
            flex-shrink: 0;
            position: relative;
            z-index: 10;
        }
        .info-panel h2 {
            margin: 0 0 8px 0;
            font-size: 1.1em;
            color: var(--vscode-textLink-foreground);
        }
        .info-grid {
            display: grid;
            grid-template-columns: auto 1fr;
            gap: 4px 15px;
        }
        .info-label {
            font-weight: bold;
            color: var(--vscode-symbolIcon-propertyForeground);
        }
        .info-value {
            font-family: var(--vscode-editor-font-family);
        }
        .controls {
            display: flex;
            gap: 10px;
            margin-bottom: 10px;
            flex-wrap: wrap;
            flex-shrink: 0;
            position: relative;
            z-index: 10;
        }
        button {
            background-color: var(--vscode-button-background);
            color: var(--vscode-button-foreground);
            border: none;
            padding: 6px 12px;
            cursor: pointer;
            border-radius: 3px;
            font-size: 12px;
        }
        button:hover {
            background-color: var(--vscode-button-hoverBackground);
        }
        button.active {
            background-color: var(--vscode-button-hoverBackground);
            outline: 1px solid var(--vscode-focusBorder);
        }
        #viewer {
            flex: 1;
            min-height: 300px;
            border: 1px solid var(--vscode-widget-border);
            border-radius: 5px;
            position: relative;
            z-index: 1;
        }
        .legend {
            display: flex;
            flex-wrap: wrap;
            gap: 10px;
            margin-top: 10px;
            padding: 8px;
            background-color: var(--vscode-editor-inactiveSelectionBackground);
            border-radius: 5px;
            flex-shrink: 0;
            position: relative;
            z-index: 10;
        }
        .legend-item {
            display: flex;
            align-items: center;
            gap: 5px;
        }
        .legend-color {
            width: 16px;
            height: 16px;
            border-radius: 50%;
            border: 1px solid var(--vscode-widget-border);
        }
        .legend-marker {
            color: red;
            font-size: 16px;
            font-weight: bold;
            width: 16px;
            text-align: center;
        }
        .legend-separator {
            margin-left: 10px;
            padding-left: 10px;
            border-left: 1px solid var(--vscode-widget-border);
        }
        .repeat-controls {
            display: flex;
            align-items: center;
            gap: 15px;
            margin-left: auto;
            padding-left: 15px;
            border-left: 1px solid var(--vscode-widget-border);
        }
        .repeat-controls label {
            font-size: 12px;
            color: var(--vscode-descriptionForeground);
        }
        .repeat-group {
            display: flex;
            align-items: center;
            gap: 3px;
        }
        .repeat-group span {
            font-size: 12px;
            font-weight: bold;
            min-width: 12px;
        }
        .spinner {
            display: flex;
            flex-direction: column;
        }
        .spinner button {
            padding: 0;
            width: 18px;
            height: 12px;
            font-size: 8px;
            line-height: 1;
            display: flex;
            align-items: center;
            justify-content: center;
            border-radius: 2px;
        }
        .spinner button:first-child {
            border-radius: 2px 2px 0 0;
        }
        .spinner button:last-child {
            border-radius: 0 0 2px 2px;
        }
        .repeat-input {
            width: 30px;
            height: 24px;
            text-align: center;
            background-color: var(--vscode-input-background);
            color: var(--vscode-input-foreground);
            border: 1px solid var(--vscode-input-border);
            border-radius: 2px;
            font-size: 12px;
        }
        .repeat-input:focus {
            outline: 1px solid var(--vscode-focusBorder);
        }
    </style>
</head>
<body>
    <div class="info-panel">
        <h2>${this._escapeHtml(data.title)}</h2>
        <div class="info-grid">
            <span class="info-label">Formula:</span>
            <span class="info-value">${this._escapeHtml(formula)}</span>
            <span class="info-label">Atoms:</span>
            <span class="info-value">${totalAtoms}</span>
            <span class="info-label">Volume:</span>
            <span class="info-value">${volume.toFixed(2)} &#8491;&sup3;</span>
            <span class="info-label">Coordinates:</span>
            <span class="info-value">${data.isDirect ? 'Direct (fractional)' : 'Cartesian'}</span>
            <span class="info-label">Constraints:</span>
            <span class="info-value">${data.selectiveDynamics ? atoms.filter(a => a.constrained).length + ' constrained (shown with X)' : 'None'}</span>
        </div>
    </div>

    <div class="controls">
        <button id="resetView">Reset View</button>
        <button id="toggleCell" class="active">Unit Cell</button>
        <button id="toggleLabels">Labels</button>
        <button id="styleStick">Stick</button>
        <button id="styleSphere" class="active">Sphere</button>
        <div class="repeat-controls">
            <label>Repeat:</label>
            <div class="repeat-group">
                <span>a</span>
                <input type="number" id="repeatA" class="repeat-input" value="1" min="1" max="10">
                <div class="spinner">
                    <button onclick="incrementRepeat('A')">&#9650;</button>
                    <button onclick="decrementRepeat('A')">&#9660;</button>
                </div>
            </div>
            <div class="repeat-group">
                <span>b</span>
                <input type="number" id="repeatB" class="repeat-input" value="1" min="1" max="10">
                <div class="spinner">
                    <button onclick="incrementRepeat('B')">&#9650;</button>
                    <button onclick="decrementRepeat('B')">&#9660;</button>
                </div>
            </div>
            <div class="repeat-group">
                <span>c</span>
                <input type="number" id="repeatC" class="repeat-input" value="1" min="1" max="10">
                <div class="spinner">
                    <button onclick="incrementRepeat('C')">&#9650;</button>
                    <button onclick="decrementRepeat('C')">&#9660;</button>
                </div>
            </div>
        </div>
    </div>

    <div id="viewer"></div>

    <div class="legend">
        ${data.elements.map((el, i) => `
            <div class="legend-item">
                <div class="legend-color" style="background-color: ${CPK_COLORS[el] || '#FF00FF'};"></div>
                <span>${el} (${data.counts[i]})</span>
            </div>
        `).join('')}
        ${data.selectiveDynamics ? `
            <div class="legend-item legend-separator">
                <div class="legend-marker">✱</div>
                <span>Constrained atom</span>
            </div>
        ` : ''}
        <div class="legend-item legend-separator">
            <button id="savePng" title="Save as PNG">Save PNG</button>
        </div>
    </div>

    <script>
        const baseAtoms = ${JSON.stringify(atoms)};
        const lattice = ${JSON.stringify(scaledLattice)};

        let viewer;
        let showCell = true;
        let showLabels = false;
        let styleMode = 'sphere';
        let repeatA = 1, repeatB = 1, repeatC = 1;

        // Generate repeated atoms based on current repeat values
        function getRepeatedAtoms() {
            const repeated = [];
            for (let ia = 0; ia < repeatA; ia++) {
                for (let ib = 0; ib < repeatB; ib++) {
                    for (let ic = 0; ic < repeatC; ic++) {
                        // Calculate translation vector
                        const tx = ia * lattice[0][0] + ib * lattice[1][0] + ic * lattice[2][0];
                        const ty = ia * lattice[0][1] + ib * lattice[1][1] + ic * lattice[2][1];
                        const tz = ia * lattice[0][2] + ib * lattice[1][2] + ic * lattice[2][2];

                        baseAtoms.forEach(atom => {
                            repeated.push({
                                elem: atom.elem,
                                x: atom.x + tx,
                                y: atom.y + ty,
                                z: atom.z + tz,
                                color: atom.color,
                                constrained: atom.constrained
                            });
                        });
                    }
                }
            }
            return repeated;
        }

        // Generate XYZ format string from atoms array
        function generateXyz(atomList) {
            return atomList.length + "\\n\\n" + atomList.map(a =>
                a.elem + " " + a.x.toFixed(6) + " " + a.y.toFixed(6) + " " + a.z.toFixed(6)
            ).join("\\n");
        }

        // Helper to determine if background is light or dark
        function isLightBackground(color) {
            // Parse hex or rgb color and calculate luminance
            let r, g, b;
            if (color.startsWith('#')) {
                const hex = color.slice(1);
                r = parseInt(hex.substr(0, 2), 16);
                g = parseInt(hex.substr(2, 2), 16);
                b = parseInt(hex.substr(4, 2), 16);
            } else if (color.startsWith('rgb')) {
                const match = color.match(/\\d+/g);
                if (match) {
                    r = parseInt(match[0]);
                    g = parseInt(match[1]);
                    b = parseInt(match[2]);
                }
            } else {
                return false; // Default to dark theme
            }
            // Calculate relative luminance
            const luminance = (0.299 * r + 0.587 * g + 0.114 * b) / 255;
            return luminance > 0.5;
        }

        let unitCellColor = 'yellow'; // Will be set in initViewer

        function initViewer() {
            const element = document.getElementById('viewer');
            // Get theme background color - this will adapt to light/dark themes
            const bgColor = getComputedStyle(document.body).getPropertyValue('--vscode-editor-background').trim() || 'black';
            // Set unit cell color based on background luminance
            unitCellColor = isLightBackground(bgColor) ? '#555555' : 'yellow';
            viewer = $3Dmol.createViewer(element, {
                backgroundColor: bgColor
            });

            addAtoms();
            // addUnitCell is called at the end of addAtoms
            viewer.zoomTo();
            viewer.render();
        }

        function addAtoms() {
            viewer.removeAllModels();
            viewer.removeAllLabels();
            viewer.removeAllShapes();

            const atoms = getRepeatedAtoms();

            // Add atoms as spheres with individual styling
            atoms.forEach((atom, i) => {
                const radius = styleMode === 'sphere' ? 0.5 : 0.2;
                const opacity = atom.constrained ? 0.9 : 1.0;

                viewer.addSphere({
                    center: {x: atom.x, y: atom.y, z: atom.z},
                    radius: radius,
                    color: atom.color,
                    opacity: opacity
                });

                // Add X marker for constrained atoms (X in XY plane + vertical line)
                if (atom.constrained) {
                    const size = radius * 1.2;
                    // Draw X using two crossed cylinders in XY plane
                    viewer.addCylinder({
                        start: {x: atom.x - size, y: atom.y - size, z: atom.z},
                        end: {x: atom.x + size, y: atom.y + size, z: atom.z},
                        radius: 0.05,
                        color: 'red',
                        fromCap: true,
                        toCap: true
                    });
                    viewer.addCylinder({
                        start: {x: atom.x + size, y: atom.y - size, z: atom.z},
                        end: {x: atom.x - size, y: atom.y + size, z: atom.z},
                        radius: 0.05,
                        color: 'red',
                        fromCap: true,
                        toCap: true
                    });
                    // Add vertical line (Z direction) for visibility from all angles
                    viewer.addCylinder({
                        start: {x: atom.x, y: atom.y, z: atom.z - size},
                        end: {x: atom.x, y: atom.y, z: atom.z + size},
                        radius: 0.05,
                        color: 'red',
                        fromCap: true,
                        toCap: true
                    });
                }

                if (showLabels) {
                    viewer.addLabel(atom.elem, {
                        position: {x: atom.x, y: atom.y, z: atom.z},
                        fontSize: 12,
                        fontColor: 'white',
                        backgroundOpacity: 0.5
                    });
                }
            });

            // Re-add unit cell since we cleared all shapes
            addUnitCell();
        }

        function addUnitCell() {
            // Note: Don't removeAllShapes here - atoms are also shapes now

            if (!showCell) return;

            const a = lattice[0];
            const b = lattice[1];
            const c = lattice[2];

            // Draw supercell based on repeat values
            // Scale lattice vectors by repeat factors
            const superA = [a[0] * repeatA, a[1] * repeatA, a[2] * repeatA];
            const superB = [b[0] * repeatB, b[1] * repeatB, b[2] * repeatB];
            const superC = [c[0] * repeatC, c[1] * repeatC, c[2] * repeatC];

            const origin = [0, 0, 0];

            // Define the 8 corners of the supercell
            const corners = [
                origin,
                superA,
                superB,
                superC,
                [superA[0] + superB[0], superA[1] + superB[1], superA[2] + superB[2]],
                [superA[0] + superC[0], superA[1] + superC[1], superA[2] + superC[2]],
                [superB[0] + superC[0], superB[1] + superC[1], superB[2] + superC[2]],
                [superA[0] + superB[0] + superC[0], superA[1] + superB[1] + superC[1], superA[2] + superB[2] + superC[2]]
            ];

            // Define the 12 edges
            const edges = [
                [0, 1], [0, 2], [0, 3],
                [1, 4], [1, 5],
                [2, 4], [2, 6],
                [3, 5], [3, 6],
                [4, 7], [5, 7], [6, 7]
            ];

            edges.forEach(([i, j]) => {
                viewer.addCylinder({
                    start: {x: corners[i][0], y: corners[i][1], z: corners[i][2]},
                    end: {x: corners[j][0], y: corners[j][1], z: corners[j][2]},
                    radius: 0.03,
                    color: unitCellColor,
                    fromCap: true,
                    toCap: true
                });
            });
        }

        function updateView() {
            addAtoms();
            // addUnitCell is called at the end of addAtoms
            viewer.render();
        }

        // Repeat control functions
        function incrementRepeat(dir) {
            const input = document.getElementById('repeat' + dir);
            const val = parseInt(input.value) || 1;
            if (val < 10) {
                input.value = val + 1;
                updateRepeatValue(dir);
            }
        }

        function decrementRepeat(dir) {
            const input = document.getElementById('repeat' + dir);
            const val = parseInt(input.value) || 1;
            if (val > 1) {
                input.value = val - 1;
                updateRepeatValue(dir);
            }
        }

        function updateRepeatValue(dir) {
            const input = document.getElementById('repeat' + dir);
            let val = parseInt(input.value) || 1;
            val = Math.max(1, Math.min(10, val));
            input.value = val;

            if (dir === 'A') repeatA = val;
            else if (dir === 'B') repeatB = val;
            else if (dir === 'C') repeatC = val;

            updateView();
            viewer.zoomTo();
            viewer.render();
        }

        document.getElementById('resetView').addEventListener('click', () => {
            viewer.zoomTo();
            viewer.render();
        });

        document.getElementById('toggleCell').addEventListener('click', (e) => {
            showCell = !showCell;
            e.target.classList.toggle('active', showCell);
            updateView();
        });

        document.getElementById('toggleLabels').addEventListener('click', (e) => {
            showLabels = !showLabels;
            e.target.classList.toggle('active', showLabels);
            viewer.removeAllLabels();
            updateView();
        });

        document.getElementById('styleStick').addEventListener('click', (e) => {
            styleMode = 'stick';
            document.getElementById('styleSphere').classList.remove('active');
            e.target.classList.add('active');
            updateView();
        });

        document.getElementById('styleSphere').addEventListener('click', (e) => {
            styleMode = 'sphere';
            document.getElementById('styleStick').classList.remove('active');
            e.target.classList.add('active');
            updateView();
        });

        // Add event listeners for repeat input fields
        ['A', 'B', 'C'].forEach(dir => {
            const input = document.getElementById('repeat' + dir);
            input.addEventListener('change', () => updateRepeatValue(dir));
            input.addEventListener('keyup', (e) => {
                if (e.key === 'Enter') updateRepeatValue(dir);
            });
        });

        // Save PNG button
        document.getElementById('savePng').addEventListener('click', () => {
            const pngUri = viewer.pngURI();
            const link = document.createElement('a');
            link.href = pngUri;
            link.download = 'structure.png';
            link.click();
        });

        // Debug: show base atoms and check if 3Dmol loaded
        console.log("Base atoms:", baseAtoms);
        console.log("Lattice:", lattice);
        console.log("3Dmol loaded:", typeof $3Dmol !== 'undefined');

        try {
            initViewer();
        } catch (e) {
            console.error("Error initializing viewer:", e);
            document.getElementById('viewer').innerHTML = '<p style="color:red;padding:20px;">Error: ' + e.message + '</p>';
        }
    </script>
</body>
</html>`;
    }

    private _getErrorHtml(message: string): string {
        return `<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Structure Preview Error</title>
    <style>
        body {
            font-family: var(--vscode-font-family);
            color: var(--vscode-foreground);
            background-color: var(--vscode-editor-background);
            padding: 20px;
        }
        .error {
            color: var(--vscode-errorForeground);
            padding: 20px;
            background-color: var(--vscode-inputValidation-errorBackground);
            border: 1px solid var(--vscode-inputValidation-errorBorder);
            border-radius: 5px;
        }
    </style>
</head>
<body>
    <div class="error">
        <h2>Error</h2>
        <p>${this._escapeHtml(message)}</p>
    </div>
</body>
</html>`;
    }

    private _escapeHtml(text: string): string {
        return text
            .replace(/&/g, '&amp;')
            .replace(/</g, '&lt;')
            .replace(/>/g, '&gt;')
            .replace(/"/g, '&quot;')
            .replace(/'/g, '&#039;');
    }
}

class SystemInfoPanel {
    public static currentPanel: SystemInfoPanel | undefined;
    private readonly _panel: vscode.WebviewPanel;
    private _disposables: vscode.Disposable[] = [];

    public static createOrShow(context: vscode.ExtensionContext) {
        const column = vscode.window.activeTextEditor
            ? vscode.window.activeTextEditor.viewColumn
            : undefined;

        if (SystemInfoPanel.currentPanel) {
            SystemInfoPanel.currentPanel._panel.reveal(column);
            SystemInfoPanel.currentPanel._update(context);
            return;
        }

        const panel = vscode.window.createWebviewPanel(
            'vaspSystemInfo',
            'VASP System Info',
            column || vscode.ViewColumn.One,
            {
                enableScripts: true,
                retainContextWhenHidden: true
            }
        );

        SystemInfoPanel.currentPanel = new SystemInfoPanel(panel, context);
    }

    private constructor(panel: vscode.WebviewPanel, context: vscode.ExtensionContext) {
        this._panel = panel;
        this._update(context);

        this._panel.onDidDispose(() => this.dispose(), null, this._disposables);
    }

    public dispose() {
        SystemInfoPanel.currentPanel = undefined;
        this._panel.dispose();
        while (this._disposables.length) {
            const x = this._disposables.pop();
            if (x) {
                x.dispose();
            }
        }
    }

    private _update(context: vscode.ExtensionContext) {
        const systemInfo = this._gatherSystemInfo(context);
        this._panel.webview.html = this._getHtmlForWebview(systemInfo);
    }

    private _gatherSystemInfo(context: vscode.ExtensionContext): Record<string, string | Record<string, string>> {
        const info: Record<string, string | Record<string, string>> = {};

        // Extension Information
        const extensionInfo: Record<string, string> = {};
        const packageJson = context.extension.packageJSON;
        extensionInfo['Extension Name'] = packageJson.displayName || packageJson.name;
        extensionInfo['Extension Version'] = packageJson.version;
        extensionInfo['Extension ID'] = context.extension.id;
        extensionInfo['Extension Path'] = context.extensionPath;
        info['Extension'] = extensionInfo;

        // VS Code Information
        const vscodeInfo: Record<string, string> = {};
        vscodeInfo['VS Code Version'] = vscode.version;
        vscodeInfo['UI Kind'] = vscode.env.uiKind === vscode.UIKind.Desktop ? 'Desktop' : 'Web';
        vscodeInfo['App Name'] = vscode.env.appName;
        vscodeInfo['App Root'] = vscode.env.appRoot;
        vscodeInfo['Language'] = vscode.env.language;
        vscodeInfo['Shell'] = vscode.env.shell;
        info['VS Code'] = vscodeInfo;

        // System Information
        const sysInfo: Record<string, string> = {};
        sysInfo['OS Platform'] = os.platform();
        sysInfo['OS Release'] = os.release();
        sysInfo['OS Type'] = os.type();
        sysInfo['Architecture'] = os.arch();
        sysInfo['Hostname'] = os.hostname();
        sysInfo['Home Directory'] = os.homedir();
        sysInfo['Temp Directory'] = os.tmpdir();
        sysInfo['Total Memory'] = `${Math.round(os.totalmem() / (1024 * 1024 * 1024))} GB`;
        sysInfo['Free Memory'] = `${Math.round(os.freemem() / (1024 * 1024 * 1024))} GB`;
        sysInfo['CPUs'] = `${os.cpus().length} cores (${os.cpus()[0]?.model || 'Unknown'})`;
        sysInfo['Node.js Version'] = process.version;
        info['System'] = sysInfo;

        // VASP Environment Variables
        const vaspEnv: Record<string, string> = {};
        const vaspEnvVars = [
            'VASP_PP_PATH', 'VASP_POTCAR_PATH', 'VASP_PSP_PATH',
            'VASP_ROOT', 'VASP_HOME', 'VASP_DIR', 'VASP_BIN',
            'VASP_COMMAND', 'ASE_VASP_COMMAND', 'VASP_SCRIPT',
            'VASP_GAMMA', 'VASP_STD', 'VASP_NCL', 'VASP_GAM',
            'POT_GGA_PAW', 'POT_LDA_PAW', 'POTCAR_PATH',
            'VASP_MPI_PROCS', 'NSLOTS', 'OMP_NUM_THREADS',
            'MKL_NUM_THREADS', 'MPI_PROCS'
        ];

        for (const envVar of vaspEnvVars) {
            const value = process.env[envVar];
            if (value) {
                vaspEnv[envVar] = value;
            }
        }

        // Check common PATH additions for VASP
        const pathDirs = process.env.PATH?.split(path.delimiter) || [];
        const vaspPathDirs = pathDirs.filter(p =>
            p.toLowerCase().includes('vasp') ||
            p.toLowerCase().includes('mpi') ||
            p.toLowerCase().includes('intel')
        );
        if (vaspPathDirs.length > 0) {
            vaspEnv['VASP-related PATH entries'] = vaspPathDirs.join('\n');
        }

        if (Object.keys(vaspEnv).length === 0) {
            vaspEnv['Status'] = 'No VASP environment variables detected';
        }
        info['VASP Environment'] = vaspEnv;

        // VASP Installation Detection
        const vaspInstall: Record<string, string> = {};

        // Try to find VASP executable
        const vaspCommands = ['vasp_std', 'vasp_gam', 'vasp_ncl', 'vasp', 'vasp6', 'vasp5'];
        for (const cmd of vaspCommands) {
            try {
                let result: string;
                if (os.platform() === 'win32') {
                    result = execSync(`where ${cmd}`, { encoding: 'utf-8', timeout: 5000 }).trim();
                } else {
                    result = execSync(`which ${cmd}`, { encoding: 'utf-8', timeout: 5000 }).trim();
                }
                if (result) {
                    vaspInstall[`${cmd} location`] = result;
                }
            } catch {
                // Command not found, continue
            }
        }

        // Try to get VASP version if available
        if (process.env.VASP_COMMAND || vaspInstall['vasp_std location']) {
            try {
                const vaspCmd = process.env.VASP_COMMAND || 'vasp_std';
                // Note: VASP doesn't have a --version flag, so we just note it's installed
                vaspInstall['VASP Command'] = vaspCmd;
            } catch {
                // Ignore errors
            }
        }

        // Check for common POTCAR directories
        const potcarPaths = [
            process.env.VASP_PP_PATH,
            process.env.VASP_POTCAR_PATH,
            process.env.POT_GGA_PAW,
            path.join(os.homedir(), 'vasp', 'potcars'),
            path.join(os.homedir(), '.vasp', 'potpaw_PBE'),
            '/opt/vasp/potcars',
            '/usr/local/vasp/potcars'
        ].filter(Boolean) as string[];

        for (const potPath of potcarPaths) {
            if (potPath && fs.existsSync(potPath)) {
                vaspInstall['POTCAR Directory'] = potPath;
                try {
                    const contents = fs.readdirSync(potPath).slice(0, 10);
                    vaspInstall['Available POTCARs (sample)'] = contents.join(', ') + (contents.length >= 10 ? '...' : '');
                } catch {
                    // Ignore read errors
                }
                break;
            }
        }

        if (Object.keys(vaspInstall).length === 0) {
            vaspInstall['Status'] = 'No VASP installation detected in PATH';
        }
        info['VASP Installation'] = vaspInstall;

        // Workspace Information
        const workspaceInfo: Record<string, string> = {};
        const workspaceFolders = vscode.workspace.workspaceFolders;
        if (workspaceFolders && workspaceFolders.length > 0) {
            workspaceInfo['Workspace Folders'] = workspaceFolders.map(f => f.uri.fsPath).join('\n');

            // Check for VASP files in workspace
            const vaspFiles: string[] = [];
            for (const folder of workspaceFolders) {
                const vaspFileNames = ['INCAR', 'POSCAR', 'CONTCAR', 'KPOINTS', 'OUTCAR', 'POTCAR', 'CHGCAR', 'WAVECAR', 'DOSCAR'];
                for (const fileName of vaspFileNames) {
                    const filePath = path.join(folder.uri.fsPath, fileName);
                    if (fs.existsSync(filePath)) {
                        vaspFiles.push(fileName);
                    }
                }
            }
            if (vaspFiles.length > 0) {
                workspaceInfo['VASP Files Found'] = vaspFiles.join(', ');
            }
        } else {
            workspaceInfo['Status'] = 'No workspace folder open';
        }
        info['Workspace'] = workspaceInfo;

        // Python/ASE Information (often used with VASP)
        const pythonInfo: Record<string, string> = {};
        try {
            const pythonVersion = execSync('python3 --version', { encoding: 'utf-8', timeout: 5000 }).trim();
            pythonInfo['Python Version'] = pythonVersion;
        } catch {
            try {
                const pythonVersion = execSync('python --version', { encoding: 'utf-8', timeout: 5000 }).trim();
                pythonInfo['Python Version'] = pythonVersion;
            } catch {
                pythonInfo['Python'] = 'Not found in PATH';
            }
        }

        // Check for ASE
        try {
            const aseVersion = execSync('python3 -c "import ase; print(ase.__version__)"', { encoding: 'utf-8', timeout: 5000 }).trim();
            pythonInfo['ASE Version'] = aseVersion;
        } catch {
            try {
                const aseVersion = execSync('python -c "import ase; print(ase.__version__)"', { encoding: 'utf-8', timeout: 5000 }).trim();
                pythonInfo['ASE Version'] = aseVersion;
            } catch {
                pythonInfo['ASE'] = 'Not installed';
            }
        }

        // Check for pymatgen
        try {
            const pmgVersion = execSync('python3 -c "import pymatgen; print(pymatgen.__version__)"', { encoding: 'utf-8', timeout: 5000 }).trim();
            pythonInfo['pymatgen Version'] = pmgVersion;
        } catch {
            try {
                const pmgVersion = execSync('python -c "import pymatgen; print(pymatgen.__version__)"', { encoding: 'utf-8', timeout: 5000 }).trim();
                pythonInfo['pymatgen Version'] = pmgVersion;
            } catch {
                // pymatgen not installed, don't report
            }
        }

        info['Python Environment'] = pythonInfo;

        return info;
    }

    private _getHtmlForWebview(info: Record<string, string | Record<string, string>>): string {
        // Build the content string for copying
        let textContent = 'VASP Extension System Information\n';
        textContent += '='.repeat(50) + '\n\n';

        for (const [section, data] of Object.entries(info)) {
            textContent += `## ${section}\n`;
            textContent += '-'.repeat(30) + '\n';
            if (typeof data === 'string') {
                textContent += `${data}\n`;
            } else {
                for (const [key, value] of Object.entries(data)) {
                    textContent += `${key}: ${value}\n`;
                }
            }
            textContent += '\n';
        }

        textContent += '\nGenerated: ' + new Date().toISOString() + '\n';

        // Build HTML sections
        let htmlSections = '';
        for (const [section, data] of Object.entries(info)) {
            htmlSections += `<div class="section">
                <h2>${this._escapeHtml(section)}</h2>
                <table>`;

            if (typeof data === 'string') {
                htmlSections += `<tr><td colspan="2">${this._escapeHtml(data)}</td></tr>`;
            } else {
                for (const [key, value] of Object.entries(data)) {
                    const displayValue = value.includes('\n')
                        ? `<pre>${this._escapeHtml(value)}</pre>`
                        : this._escapeHtml(value);
                    htmlSections += `<tr><td class="key">${this._escapeHtml(key)}</td><td class="value">${displayValue}</td></tr>`;
                }
            }

            htmlSections += `</table></div>`;
        }

        return `<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>VASP System Info</title>
    <style>
        body {
            font-family: var(--vscode-font-family);
            font-size: var(--vscode-font-size);
            color: var(--vscode-foreground);
            background-color: var(--vscode-editor-background);
            padding: 20px;
            line-height: 1.6;
        }
        h1 {
            color: var(--vscode-titleBar-activeForeground);
            border-bottom: 2px solid var(--vscode-titleBar-activeForeground);
            padding-bottom: 10px;
        }
        h2 {
            color: var(--vscode-textLink-foreground);
            margin-top: 20px;
            margin-bottom: 10px;
        }
        .section {
            margin-bottom: 25px;
            background-color: var(--vscode-editor-inactiveSelectionBackground);
            padding: 15px;
            border-radius: 5px;
        }
        table {
            width: 100%;
            border-collapse: collapse;
        }
        td {
            padding: 8px;
            border-bottom: 1px solid var(--vscode-widget-border);
            vertical-align: top;
        }
        .key {
            font-weight: bold;
            width: 35%;
            color: var(--vscode-symbolIcon-propertyForeground);
        }
        .value {
            word-break: break-all;
            font-family: var(--vscode-editor-font-family);
        }
        pre {
            margin: 0;
            white-space: pre-wrap;
            word-wrap: break-word;
            font-family: var(--vscode-editor-font-family);
            font-size: var(--vscode-editor-font-size);
        }
        .button-container {
            margin: 20px 0;
            display: flex;
            gap: 10px;
        }
        button {
            background-color: var(--vscode-button-background);
            color: var(--vscode-button-foreground);
            border: none;
            padding: 10px 20px;
            cursor: pointer;
            border-radius: 3px;
            font-size: 14px;
        }
        button:hover {
            background-color: var(--vscode-button-hoverBackground);
        }
        .copy-success {
            color: var(--vscode-testing-iconPassed);
            margin-left: 10px;
            display: none;
        }
        .timestamp {
            color: var(--vscode-descriptionForeground);
            font-size: 0.9em;
            margin-top: 20px;
        }
        textarea {
            position: absolute;
            left: -9999px;
        }
    </style>
</head>
<body>
    <h1>VASP Extension System Information</h1>

    <div class="button-container">
        <button onclick="copyToClipboard()">Copy All to Clipboard</button>
        <span class="copy-success" id="copySuccess">Copied!</span>
    </div>

    ${htmlSections}

    <p class="timestamp">Generated: ${new Date().toISOString()}</p>

    <textarea id="copyText" readonly>${this._escapeHtml(textContent)}</textarea>

    <script>
        function copyToClipboard() {
            const text = document.getElementById('copyText').value;
            navigator.clipboard.writeText(text).then(() => {
                const successMsg = document.getElementById('copySuccess');
                successMsg.style.display = 'inline';
                setTimeout(() => {
                    successMsg.style.display = 'none';
                }, 2000);
            }).catch(err => {
                // Fallback for older browsers
                const textarea = document.getElementById('copyText');
                textarea.style.position = 'static';
                textarea.select();
                document.execCommand('copy');
                textarea.style.position = 'absolute';
                textarea.style.left = '-9999px';

                const successMsg = document.getElementById('copySuccess');
                successMsg.style.display = 'inline';
                setTimeout(() => {
                    successMsg.style.display = 'none';
                }, 2000);
            });
        }
    </script>
</body>
</html>`;
    }

    private _escapeHtml(text: string): string {
        return text
            .replace(/&/g, '&amp;')
            .replace(/</g, '&lt;')
            .replace(/>/g, '&gt;')
            .replace(/"/g, '&quot;')
            .replace(/'/g, '&#039;');
    }
}

// Interface for VASP calculation summary data
interface VaspSummaryData {
    directory: string;
    formula: string;
    atomCount: number;
    energy: number | null;
    energyPerAtom: number | null;
    fermiEnergy: number | null;
    maxForce: number | null;
    forces: number[][] | null;
    stress: number[][] | null;
    pressure: number | null;
    converged: boolean | null;
    electronicSteps: number | null;
    ionicSteps: number | null;
    kpoints: string;
    encut: number | null;
    xc: string;
    ispin: number | null;
    ediff: number | null;
    ediffg: number | null;
    ibrion: number | null;
    nsw: number | null;
    elements: string[];
    atomCounts: number[];
    lattice: number[][] | null;
    volume: number | null;
}

// Simple XML tag extraction helper
function extractXmlTag(xml: string, tagName: string): string | null {
    const regex = new RegExp(`<${tagName}[^>]*>([\\s\\S]*?)</${tagName}>`, 'i');
    const match = xml.match(regex);
    return match ? match[1].trim() : null;
}

// Extract all matches for a tag
function extractAllXmlTags(xml: string, tagName: string): string[] {
    const regex = new RegExp(`<${tagName}[^>]*>([\\s\\S]*?)</${tagName}>`, 'gi');
    const matches: string[] = [];
    let match;
    while ((match = regex.exec(xml)) !== null) {
        matches.push(match[1].trim());
    }
    return matches;
}

// Extract attribute value from an XML tag
function extractXmlAttr(xml: string, attrName: string): string | null {
    const regex = new RegExp(`${attrName}="([^"]*)"`, 'i');
    const match = xml.match(regex);
    return match ? match[1] : null;
}

// Parse vasprun.xml to extract calculation summary
function parseVasprunXml(content: string): VaspSummaryData | null {
    try {
        // Extract basic structure info
        const atomInfoBlock = extractXmlTag(content, 'atominfo');
        let elements: string[] = [];
        let atomCounts: number[] = [];
        let atomCount = 0;

        if (atomInfoBlock) {
            // Extract number of atoms (handle whitespace around number)
            const atomsMatch = atomInfoBlock.match(/<atoms>\s*(\d+)\s*<\/atoms>/);
            if (atomsMatch) {
                atomCount = parseInt(atomsMatch[1]);
            }

            // Extract element types from atomtypes array
            const typesBlock = extractXmlTag(atomInfoBlock, 'array');
            if (typesBlock) {
                const typeMatches = typesBlock.match(/<c>([A-Z][a-z]?)\s*<\/c>/g);
                if (typeMatches) {
                    elements = typeMatches.map(m => {
                        const match = m.match(/<c>([A-Z][a-z]?)\s*<\/c>/);
                        return match ? match[1].trim() : '';
                    }).filter(e => e.length > 0);
                }
            }

            // Count atoms per element
            const atomsBlock = content.match(/<array name="atoms"[\s\S]*?<\/array>/);
            if (atomsBlock) {
                const atomTypes = atomsBlock[0].match(/<c>([A-Z][a-z]?)\s*<\/c>/g);
                if (atomTypes) {
                    const counts: { [key: string]: number } = {};
                    atomTypes.forEach(m => {
                        const match = m.match(/<c>([A-Z][a-z]?)\s*<\/c>/);
                        if (match) {
                            const el = match[1].trim();
                            counts[el] = (counts[el] || 0) + 1;
                        }
                    });
                    elements = Object.keys(counts);
                    atomCounts = elements.map(e => counts[e]);
                }
            }
        }

        // Generate formula
        const formula = elements.map((el, i) => {
            const count = atomCounts[i] || 1;
            return count === 1 ? el : `${el}${count}`;
        }).join('');

        // Extract final energy from last calculation block
        let energy: number | null = null;
        const energyMatches = content.match(/<i name="e_fr_energy">\s*([-\d.Ee+]+)\s*<\/i>/g);
        if (energyMatches && energyMatches.length > 0) {
            const lastMatch = energyMatches[energyMatches.length - 1].match(/<i name="e_fr_energy">\s*([-\d.Ee+]+)\s*<\/i>/);
            if (lastMatch) {
                energy = parseFloat(lastMatch[1]);
            }
        }

        // Extract Fermi energy
        let fermiEnergy: number | null = null;
        const fermiMatch = content.match(/<i name="efermi">\s*([-\d.Ee+]+)\s*<\/i>/);
        if (fermiMatch) {
            fermiEnergy = parseFloat(fermiMatch[1]);
        }

        // Extract forces from last ionic step
        let forces: number[][] | null = null;
        let maxForce: number | null = null;
        const forceBlocks = content.match(/<varray name="forces"[\s\S]*?<\/varray>/g);
        if (forceBlocks && forceBlocks.length > 0) {
            const lastForceBlock = forceBlocks[forceBlocks.length - 1];
            const forceVectors = lastForceBlock.match(/<v>\s*([-\d.Ee+\s]+)\s*<\/v>/g);
            if (forceVectors) {
                forces = forceVectors.map(v => {
                    const match = v.match(/<v>\s*([-\d.Ee+\s]+)\s*<\/v>/);
                    if (match) {
                        return match[1].trim().split(/\s+/).map(parseFloat);
                    }
                    return [0, 0, 0];
                });
                // Calculate max force magnitude
                maxForce = Math.max(...forces.map(f =>
                    Math.sqrt(f[0] * f[0] + f[1] * f[1] + f[2] * f[2])
                ));
            }
        }

        // Extract stress tensor from last ionic step
        let stress: number[][] | null = null;
        let pressure: number | null = null;
        const stressBlocks = content.match(/<varray name="stress"[\s\S]*?<\/varray>/g);
        if (stressBlocks && stressBlocks.length > 0) {
            const lastStressBlock = stressBlocks[stressBlocks.length - 1];
            const stressVectors = lastStressBlock.match(/<v>\s*([-\d.Ee+\s]+)\s*<\/v>/g);
            if (stressVectors) {
                stress = stressVectors.map(v => {
                    const match = v.match(/<v>\s*([-\d.Ee+\s]+)\s*<\/v>/);
                    if (match) {
                        return match[1].trim().split(/\s+/).map(parseFloat);
                    }
                    return [0, 0, 0];
                });
                // Calculate pressure as -1/3 * trace(stress) in kBar, convert to GPa
                if (stress.length >= 3) {
                    pressure = -(stress[0][0] + stress[1][1] + stress[2][2]) / 3 / 10; // kBar to GPa
                }
            }
        }

        // Extract lattice vectors
        let lattice: number[][] | null = null;
        let volume: number | null = null;
        const basisBlocks = content.match(/<varray name="basis"[\s\S]*?<\/varray>/g);
        if (basisBlocks && basisBlocks.length > 0) {
            const lastBasisBlock = basisBlocks[basisBlocks.length - 1];
            const basisVectors = lastBasisBlock.match(/<v>\s*([-\d.Ee+\s]+)\s*<\/v>/g);
            if (basisVectors && basisVectors.length >= 3) {
                lattice = basisVectors.slice(0, 3).map(v => {
                    const match = v.match(/<v>\s*([-\d.Ee+\s]+)\s*<\/v>/);
                    if (match) {
                        return match[1].trim().split(/\s+/).map(parseFloat);
                    }
                    return [0, 0, 0];
                });
                // Calculate volume
                const a = lattice[0];
                const b = lattice[1];
                const c = lattice[2];
                const cross = [
                    b[1] * c[2] - b[2] * c[1],
                    b[2] * c[0] - b[0] * c[2],
                    b[0] * c[1] - b[1] * c[0]
                ];
                volume = Math.abs(a[0] * cross[0] + a[1] * cross[1] + a[2] * cross[2]);
            }
        }

        // Extract parameters
        let encut: number | null = null;
        const encutMatch = content.match(/<i[^>]*name="ENCUT"[^>]*>\s*([-\d.Ee+]+)\s*<\/i>/i);
        if (encutMatch) {
            encut = parseFloat(encutMatch[1]);
        }

        let xc = 'Unknown';
        const xcMatch = content.match(/<i[^>]*name="GGA"[^>]*>\s*(\w+)\s*<\/i>/i);
        if (xcMatch) {
            xc = xcMatch[1].trim();
        }
        // Also check for METAGGA
        const metaGgaMatch = content.match(/<i[^>]*name="METAGGA"[^>]*>\s*(\w+)\s*<\/i>/i);
        if (metaGgaMatch && metaGgaMatch[1].trim() !== '--') {
            xc = metaGgaMatch[1].trim();
        }

        let ispin: number | null = null;
        const ispinMatch = content.match(/<i[^>]*name="ISPIN"[^>]*>\s*(\d+)\s*<\/i>/i);
        if (ispinMatch) {
            ispin = parseInt(ispinMatch[1]);
        }

        let ediff: number | null = null;
        const ediffMatch = content.match(/<i[^>]*name="EDIFF"[^>]*>\s*([-\d.Ee+]+)\s*<\/i>/i);
        if (ediffMatch) {
            ediff = parseFloat(ediffMatch[1]);
        }

        let ediffg: number | null = null;
        const ediffgMatch = content.match(/<i[^>]*name="EDIFFG"[^>]*>\s*([-\d.Ee+]+)\s*<\/i>/i);
        if (ediffgMatch) {
            ediffg = parseFloat(ediffgMatch[1]);
        }

        let ibrion: number | null = null;
        const ibrionMatch = content.match(/<i[^>]*name="IBRION"[^>]*>\s*([-\d]+)\s*<\/i>/i);
        if (ibrionMatch) {
            ibrion = parseInt(ibrionMatch[1]);
        }

        let nsw: number | null = null;
        const nswMatch = content.match(/<i[^>]*name="NSW"[^>]*>\s*(\d+)\s*<\/i>/i);
        if (nswMatch) {
            nsw = parseInt(nswMatch[1]);
        }

        // Extract k-points info
        let kpoints = 'Unknown';
        const kpointsMatch = content.match(/<varray name="kpointlist"[\s\S]*?<\/varray>/);
        if (kpointsMatch) {
            const kptCount = (kpointsMatch[0].match(/<v>/g) || []).length;
            kpoints = `${kptCount} k-points`;
        }

        // Count electronic and ionic steps
        const scfBlocks = content.match(/<scstep>/g);
        const electronicSteps = scfBlocks ? scfBlocks.length : null;

        const ionicStepBlocks = content.match(/<calculation>/g);
        const ionicSteps = ionicStepBlocks ? ionicStepBlocks.length : null;

        // Check convergence
        let converged: boolean | null = null;
        // Electronic convergence - check if NELM was reached
        const nelmMatch = content.match(/<i[^>]*name="NELM"[^>]*>\s*(\d+)\s*<\/i>/i);
        if (nelmMatch && scfBlocks && ionicStepBlocks) {
            const nelm = parseInt(nelmMatch[1]);
            const avgScfPerIonic = scfBlocks.length / ionicStepBlocks.length;
            converged = avgScfPerIonic < nelm;
        }
        // For relaxations (IBRION > 0, NSW > 0): if it stopped before NSW steps, it converged
        // This is more reliable than checking maxForce because VASP only considers forces on free atoms
        if (ibrion !== null && ibrion > 0 && nsw !== null && nsw > 0 && ionicSteps !== null) {
            converged = ionicSteps < nsw;
        }

        return {
            directory: '',
            formula: formula || 'Unknown',
            atomCount,
            energy,
            energyPerAtom: energy !== null && atomCount > 0 ? energy / atomCount : null,
            fermiEnergy,
            maxForce,
            forces,
            stress,
            pressure,
            converged,
            electronicSteps,
            ionicSteps,
            kpoints,
            encut,
            xc,
            ispin,
            ediff,
            ediffg,
            ibrion,
            nsw,
            elements,
            atomCounts,
            lattice,
            volume
        };
    } catch (error) {
        console.error('Error parsing vasprun.xml:', error);
        return null;
    }
}

// Find vasprun.xml in workspace or active file's directory
async function findVasprunXml(): Promise<string | null> {
    // First check if a vasprun.xml is open
    const editor = vscode.window.activeTextEditor;
    if (editor) {
        const currentFile = editor.document.uri.fsPath;
        if (path.basename(currentFile) === 'vasprun.xml') {
            return currentFile;
        }
        // Check in same directory as current file
        const dir = path.dirname(currentFile);
        const vasprunPath = path.join(dir, 'vasprun.xml');
        if (fs.existsSync(vasprunPath)) {
            return vasprunPath;
        }
    }

    // Check workspace folders
    const workspaceFolders = vscode.workspace.workspaceFolders;
    if (workspaceFolders) {
        for (const folder of workspaceFolders) {
            const vasprunPath = path.join(folder.uri.fsPath, 'vasprun.xml');
            if (fs.existsSync(vasprunPath)) {
                return vasprunPath;
            }
        }
    }

    return null;
}

// VASP Summary Panel for displaying calculation results
class VaspSummaryPanel {
    public static currentPanel: VaspSummaryPanel | undefined;
    private readonly _panel: vscode.WebviewPanel;
    private _disposables: vscode.Disposable[] = [];

    public static async createOrShow() {
        const column = vscode.ViewColumn.Beside;

        // Find vasprun.xml
        const vasprunPath = await findVasprunXml();
        if (!vasprunPath) {
            vscode.window.showErrorMessage('No vasprun.xml found in current directory or workspace');
            return;
        }

        if (VaspSummaryPanel.currentPanel) {
            VaspSummaryPanel.currentPanel._panel.reveal(column);
            VaspSummaryPanel.currentPanel._update(vasprunPath);
            return;
        }

        const panel = vscode.window.createWebviewPanel(
            'vaspSummary',
            'VASP Summary',
            column,
            {
                enableScripts: true,
                retainContextWhenHidden: true
            }
        );

        VaspSummaryPanel.currentPanel = new VaspSummaryPanel(panel, vasprunPath);
    }

    private constructor(panel: vscode.WebviewPanel, vasprunPath: string) {
        this._panel = panel;
        this._update(vasprunPath);

        this._panel.onDidDispose(() => this.dispose(), null, this._disposables);
    }

    public dispose() {
        VaspSummaryPanel.currentPanel = undefined;
        this._panel.dispose();
        while (this._disposables.length) {
            const x = this._disposables.pop();
            if (x) {
                x.dispose();
            }
        }
    }

    private _update(vasprunPath: string) {
        this._panel.title = 'VASP Summary: ' + path.basename(path.dirname(vasprunPath));

        try {
            const content = fs.readFileSync(vasprunPath, 'utf-8');
            const data = parseVasprunXml(content);
            if (!data) {
                this._panel.webview.html = this._getErrorHtml('Failed to parse vasprun.xml');
                return;
            }
            data.directory = path.dirname(vasprunPath);
            this._panel.webview.html = this._getHtmlForWebview(data);
        } catch (error) {
            this._panel.webview.html = this._getErrorHtml(`Error reading file: ${error}`);
        }
    }

    private _getErrorHtml(message: string): string {
        return `<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <style>
        body {
            font-family: var(--vscode-font-family);
            color: var(--vscode-foreground);
            background-color: var(--vscode-editor-background);
            padding: 20px;
        }
        .error { color: var(--vscode-errorForeground); }
    </style>
</head>
<body>
    <h2 class="error">Error</h2>
    <p>${message}</p>
</body>
</html>`;
    }

    private _getHtmlForWebview(data: VaspSummaryData): string {
        const convergenceStatus = data.converged === null ? 'Unknown' :
            (data.converged ? '✓ Converged' : '✗ Not Converged');
        const convergenceClass = data.converged === null ? 'unknown' :
            (data.converged ? 'converged' : 'not-converged');

        const formatEnergy = (e: number | null) => e !== null ? e.toFixed(6) : 'N/A';
        const formatForce = (f: number | null) => f !== null ? f.toFixed(6) : 'N/A';
        const formatPressure = (p: number | null) => p !== null ? p.toFixed(3) : 'N/A';

        // Build stress tensor display
        let stressHtml = 'N/A';
        if (data.stress && data.stress.length >= 3) {
            stressHtml = `<table class="matrix">
                ${data.stress.map(row => `<tr>${row.map(v => `<td>${v.toFixed(2)}</td>`).join('')}</tr>`).join('')}
            </table>`;
        }

        // Build forces summary
        let forcesHtml = '';
        if (data.forces && data.elements.length > 0) {
            let atomIdx = 0;
            for (let i = 0; i < data.elements.length; i++) {
                const count = data.atomCounts[i] || 1;
                for (let j = 0; j < count && atomIdx < data.forces.length; j++) {
                    const f = data.forces[atomIdx];
                    const mag = Math.sqrt(f[0] * f[0] + f[1] * f[1] + f[2] * f[2]);
                    forcesHtml += `<tr>
                        <td>${atomIdx + 1}</td>
                        <td>${data.elements[i]}</td>
                        <td>${f[0].toFixed(4)}</td>
                        <td>${f[1].toFixed(4)}</td>
                        <td>${f[2].toFixed(4)}</td>
                        <td class="${mag > 0.05 ? 'highlight' : ''}">${mag.toFixed(4)}</td>
                    </tr>`;
                    atomIdx++;
                }
            }
        }

        return `<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <meta http-equiv="Content-Security-Policy" content="default-src 'none'; style-src 'unsafe-inline'; script-src 'unsafe-inline';">
    <title>VASP Summary</title>
    <style>
        body {
            font-family: var(--vscode-font-family);
            font-size: var(--vscode-font-size);
            color: var(--vscode-foreground);
            background-color: var(--vscode-editor-background);
            margin: 0;
            padding: 15px;
            line-height: 1.5;
        }
        h1 {
            color: var(--vscode-textLink-foreground);
            font-size: 1.4em;
            margin-bottom: 5px;
            border-bottom: 1px solid var(--vscode-widget-border);
            padding-bottom: 10px;
        }
        h2 {
            color: var(--vscode-textLink-foreground);
            font-size: 1.1em;
            margin-top: 20px;
            margin-bottom: 10px;
        }
        .directory {
            font-size: 0.85em;
            color: var(--vscode-descriptionForeground);
            margin-bottom: 15px;
        }
        .section {
            background-color: var(--vscode-editor-inactiveSelectionBackground);
            padding: 12px 15px;
            border-radius: 5px;
            margin-bottom: 15px;
        }
        .grid {
            display: grid;
            grid-template-columns: 160px 1fr;
            gap: 8px 15px;
        }
        .label {
            font-weight: bold;
            color: var(--vscode-symbolIcon-propertyForeground);
        }
        .value {
            font-family: var(--vscode-editor-font-family);
        }
        .energy {
            font-size: 1.1em;
            font-weight: bold;
        }
        .converged { color: #4ec9b0; }
        .not-converged { color: #f14c4c; }
        .unknown { color: var(--vscode-descriptionForeground); }
        .highlight { color: #f14c4c; font-weight: bold; }
        table {
            border-collapse: collapse;
            font-family: var(--vscode-editor-font-family);
            font-size: 0.9em;
        }
        table.matrix td {
            padding: 3px 10px;
            text-align: right;
        }
        table.forces {
            width: 100%;
            margin-top: 10px;
        }
        table.forces th, table.forces td {
            padding: 5px 8px;
            text-align: right;
            border-bottom: 1px solid var(--vscode-widget-border);
        }
        table.forces th {
            background-color: var(--vscode-editor-background);
            font-weight: bold;
        }
        table.forces td:first-child, table.forces th:first-child {
            text-align: center;
        }
        table.forces td:nth-child(2), table.forces th:nth-child(2) {
            text-align: left;
        }
        .params-grid {
            display: grid;
            grid-template-columns: repeat(auto-fill, minmax(180px, 1fr));
            gap: 8px;
        }
        .param-item {
            display: flex;
            gap: 8px;
        }
        .collapsible {
            cursor: pointer;
            user-select: none;
        }
        .collapsible:before {
            content: '▶ ';
            font-size: 0.8em;
        }
        .collapsible.active:before {
            content: '▼ ';
        }
        .collapse-content {
            display: none;
            margin-top: 10px;
        }
        .collapse-content.show {
            display: block;
        }
    </style>
</head>
<body>
    <h1>${this._escapeHtml(data.formula)}</h1>
    <div class="directory">${this._escapeHtml(data.directory)}</div>

    <div class="section">
        <div class="grid">
            <span class="label">Total Energy:</span>
            <span class="value energy">${formatEnergy(data.energy)} eV</span>
            <span class="label">Energy/Atom:</span>
            <span class="value">${formatEnergy(data.energyPerAtom)} eV</span>
            <span class="label">Fermi Energy:</span>
            <span class="value">${formatEnergy(data.fermiEnergy)} eV</span>
            <span class="label">Convergence:</span>
            <span class="value ${convergenceClass}">${convergenceStatus}</span>
        </div>
    </div>

    <div class="section">
        <h2>Structure</h2>
        <div class="grid">
            <span class="label">Formula:</span>
            <span class="value">${this._escapeHtml(data.formula)}</span>
            <span class="label">Atoms:</span>
            <span class="value">${data.atomCount}</span>
            <span class="label">Volume:</span>
            <span class="value">${data.volume !== null ? data.volume.toFixed(3) : 'N/A'} Å³</span>
            <span class="label">Composition:</span>
            <span class="value">${data.elements.map((e, i) => `${e}: ${data.atomCounts[i]}`).join(', ')}</span>
        </div>
    </div>

    <div class="section">
        <h2>Forces & Stress</h2>
        <div class="grid">
            <span class="label">Max Force:</span>
            <span class="value ${data.maxForce !== null && data.maxForce > 0.05 ? 'highlight' : ''}">${formatForce(data.maxForce)} eV/Å</span>
            <span class="label">Pressure:</span>
            <span class="value">${formatPressure(data.pressure)} GPa</span>
        </div>
        ${data.stress ? `
        <h3 class="collapsible" onclick="toggleCollapse(this)">Stress Tensor (kBar)</h3>
        <div class="collapse-content">${stressHtml}</div>
        ` : ''}
        ${data.forces && data.forces.length > 0 ? `
        <h3 class="collapsible" onclick="toggleCollapse(this)">Atomic Forces (eV/Å)</h3>
        <div class="collapse-content">
            <table class="forces">
                <tr><th>#</th><th>Element</th><th>Fx</th><th>Fy</th><th>Fz</th><th>|F|</th></tr>
                ${forcesHtml}
            </table>
        </div>
        ` : ''}
    </div>

    <div class="section">
        <h2>Calculation Parameters</h2>
        <div class="params-grid">
            <div class="param-item"><span class="label">XC:</span><span class="value">${this._escapeHtml(data.xc)}</span></div>
            <div class="param-item"><span class="label">ENCUT:</span><span class="value">${data.encut !== null ? data.encut.toFixed(0) : 'N/A'} eV</span></div>
            <div class="param-item"><span class="label">K-points:</span><span class="value">${this._escapeHtml(data.kpoints)}</span></div>
            <div class="param-item"><span class="label">ISPIN:</span><span class="value">${data.ispin !== null ? data.ispin : 'N/A'}</span></div>
            <div class="param-item"><span class="label">IBRION:</span><span class="value">${data.ibrion !== null ? data.ibrion : 'N/A'}</span></div>
            <div class="param-item"><span class="label">NSW:</span><span class="value">${data.nsw !== null ? data.nsw : 'N/A'}</span></div>
            <div class="param-item"><span class="label">EDIFF:</span><span class="value">${data.ediff !== null ? data.ediff.toExponential(1) : 'N/A'}</span></div>
            <div class="param-item"><span class="label">EDIFFG:</span><span class="value">${data.ediffg !== null ? data.ediffg.toExponential(1) : 'N/A'}</span></div>
        </div>
    </div>

    <div class="section">
        <h2>Iteration Summary</h2>
        <div class="grid">
            <span class="label">Ionic Steps:</span>
            <span class="value">${data.ionicSteps !== null ? data.ionicSteps : 'N/A'}</span>
            <span class="label">Electronic Steps:</span>
            <span class="value">${data.electronicSteps !== null ? data.electronicSteps : 'N/A'}</span>
        </div>
    </div>

    <script>
        function toggleCollapse(element) {
            element.classList.toggle('active');
            const content = element.nextElementSibling;
            content.classList.toggle('show');
        }
    </script>
</body>
</html>`;
    }

    private _escapeHtml(text: string): string {
        return text
            .replace(/&/g, '&amp;')
            .replace(/</g, '&lt;')
            .replace(/>/g, '&gt;')
            .replace(/"/g, '&quot;')
            .replace(/'/g, '&#039;');
    }
}

// Interface for relaxation step data
interface RelaxationStep {
    energy: number;
    positions: number[][];  // Nx3 fractional coordinates
    forces: number[][] | null;
    maxForce: number | null;
}

// Interface for full relaxation data
interface RelaxationData {
    formula: string;
    elements: string[];
    atomCounts: number[];
    lattice: number[][];
    scale: number;
    steps: RelaxationStep[];
    constraints: boolean[];  // true = atom is constrained (any direction has F)
}

// Parse relaxation steps from vasprun.xml
function parseRelaxationSteps(content: string): RelaxationData | null {
    try {
        // Extract atom info
        const atomInfoBlock = content.match(/<atominfo>[\s\S]*?<\/atominfo>/);
        let elements: string[] = [];
        let atomCounts: number[] = [];

        if (atomInfoBlock) {
            const atomsBlock = content.match(/<array name="atoms"[\s\S]*?<\/array>/);
            if (atomsBlock) {
                const atomTypes = atomsBlock[0].match(/<c>([A-Z][a-z]?)\s*<\/c>/g);
                if (atomTypes) {
                    const counts: { [key: string]: number } = {};
                    atomTypes.forEach(m => {
                        const match = m.match(/<c>([A-Z][a-z]?)\s*<\/c>/);
                        if (match) {
                            const el = match[1].trim();
                            counts[el] = (counts[el] || 0) + 1;
                        }
                    });
                    elements = Object.keys(counts);
                    atomCounts = elements.map(e => counts[e]);
                }
            }
        }

        const formula = elements.map((el, i) => {
            const count = atomCounts[i] || 1;
            return count === 1 ? el : `${el}${count}`;
        }).join('');

        // Extract initial lattice (use first basis found)
        let lattice: number[][] = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
        const basisMatch = content.match(/<varray name="basis"[\s\S]*?<\/varray>/);
        if (basisMatch) {
            const basisVectors = basisMatch[0].match(/<v>\s*([-\d.Ee+\s]+)\s*<\/v>/g);
            if (basisVectors && basisVectors.length >= 3) {
                lattice = basisVectors.slice(0, 3).map(v => {
                    const match = v.match(/<v>\s*([-\d.Ee+\s]+)\s*<\/v>/);
                    if (match) {
                        return match[1].trim().split(/\s+/).map(parseFloat);
                    }
                    return [0, 0, 0];
                });
            }
        }

        // Extract selective dynamics (constraints) - atom is constrained if ANY direction has F
        let constraints: boolean[] = [];
        const selectiveMatch = content.match(/<varray name="selective"[^>]*>[\s\S]*?<\/varray>/);
        if (selectiveMatch) {
            const selectiveVectors = selectiveMatch[0].match(/<v[^>]*>\s*([TF]\s+[TF]\s+[TF])\s*<\/v>/g);
            if (selectiveVectors) {
                constraints = selectiveVectors.map(v => {
                    const match = v.match(/<v[^>]*>\s*([TF])\s+([TF])\s+([TF])\s*<\/v>/);
                    if (match) {
                        // Atom is constrained if ANY direction is F (fixed)
                        return match[1] === 'F' || match[2] === 'F' || match[3] === 'F';
                    }
                    return false;
                });
            }
        }

        // Extract all calculation blocks (ionic steps)
        const calculationBlocks = content.match(/<calculation>[\s\S]*?<\/calculation>/g);
        if (!calculationBlocks || calculationBlocks.length === 0) {
            return null;
        }

        const steps: RelaxationStep[] = [];

        for (const calcBlock of calculationBlocks) {
            // Extract energy
            let energy = 0;
            const energyMatch = calcBlock.match(/<i name="e_fr_energy">\s*([-\d.Ee+]+)\s*<\/i>/);
            if (energyMatch) {
                energy = parseFloat(energyMatch[1]);
            }

            // Extract positions
            let positions: number[][] = [];
            const posMatch = calcBlock.match(/<varray name="positions"[\s\S]*?<\/varray>/);
            if (posMatch) {
                const posVectors = posMatch[0].match(/<v>\s*([-\d.Ee+\s]+)\s*<\/v>/g);
                if (posVectors) {
                    positions = posVectors.map(v => {
                        const match = v.match(/<v>\s*([-\d.Ee+\s]+)\s*<\/v>/);
                        if (match) {
                            return match[1].trim().split(/\s+/).map(parseFloat);
                        }
                        return [0, 0, 0];
                    });
                }
            }

            // Extract forces
            let forces: number[][] | null = null;
            let maxForce: number | null = null;
            const forceMatch = calcBlock.match(/<varray name="forces"[\s\S]*?<\/varray>/);
            if (forceMatch) {
                const forceVectors = forceMatch[0].match(/<v>\s*([-\d.Ee+\s]+)\s*<\/v>/g);
                if (forceVectors) {
                    forces = forceVectors.map(v => {
                        const match = v.match(/<v>\s*([-\d.Ee+\s]+)\s*<\/v>/);
                        if (match) {
                            return match[1].trim().split(/\s+/).map(parseFloat);
                        }
                        return [0, 0, 0];
                    });
                    // Calculate max force only on UNCONSTRAINED atoms
                    const forceMagnitudes = forces.map((f, i) => {
                        // If atom is constrained, don't include its force in max calculation
                        if (constraints.length > 0 && constraints[i]) {
                            return 0;
                        }
                        return Math.sqrt(f[0] * f[0] + f[1] * f[1] + f[2] * f[2]);
                    });
                    maxForce = Math.max(...forceMagnitudes);
                }
            }

            steps.push({ energy, positions, forces, maxForce });
        }

        return {
            formula,
            elements,
            atomCounts,
            lattice,
            scale: 1.0,
            steps,
            constraints
        };
    } catch (error) {
        console.error('Error parsing relaxation steps:', error);
        return null;
    }
}

// Relaxation Viewer Panel
class RelaxationViewerPanel {
    public static currentPanel: RelaxationViewerPanel | undefined;
    private readonly _panel: vscode.WebviewPanel;
    private _disposables: vscode.Disposable[] = [];

    public static async createOrShow() {
        const column = vscode.ViewColumn.One;

        const vasprunPath = await findVasprunXml();
        if (!vasprunPath) {
            vscode.window.showErrorMessage('No vasprun.xml found in current directory or workspace');
            return;
        }

        if (RelaxationViewerPanel.currentPanel) {
            RelaxationViewerPanel.currentPanel._panel.reveal(column);
            RelaxationViewerPanel.currentPanel._update(vasprunPath);
            return;
        }

        const panel = vscode.window.createWebviewPanel(
            'vaspRelaxation',
            'Relaxation Progress',
            column,
            {
                enableScripts: true,
                retainContextWhenHidden: true
            }
        );

        RelaxationViewerPanel.currentPanel = new RelaxationViewerPanel(panel, vasprunPath);
    }

    private constructor(panel: vscode.WebviewPanel, vasprunPath: string) {
        this._panel = panel;
        this._update(vasprunPath);

        this._panel.onDidDispose(() => this.dispose(), null, this._disposables);
    }

    public dispose() {
        RelaxationViewerPanel.currentPanel = undefined;
        this._panel.dispose();
        while (this._disposables.length) {
            const x = this._disposables.pop();
            if (x) {
                x.dispose();
            }
        }
    }

    private _update(vasprunPath: string) {
        this._panel.title = 'Relaxation: ' + path.basename(path.dirname(vasprunPath));

        try {
            const content = fs.readFileSync(vasprunPath, 'utf-8');
            const data = parseRelaxationSteps(content);
            if (!data || data.steps.length === 0) {
                this._panel.webview.html = this._getErrorHtml('No relaxation steps found in vasprun.xml');
                return;
            }
            this._panel.webview.html = this._getHtmlForWebview(data);
        } catch (error) {
            this._panel.webview.html = this._getErrorHtml(`Error reading file: ${error}`);
        }
    }

    private _getErrorHtml(message: string): string {
        return `<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <style>
        body { font-family: var(--vscode-font-family); color: var(--vscode-foreground); background-color: var(--vscode-editor-background); padding: 20px; }
        .error { color: var(--vscode-errorForeground); }
    </style>
</head>
<body>
    <h2 class="error">Error</h2>
    <p>${message}</p>
</body>
</html>`;
    }

    private _getHtmlForWebview(data: RelaxationData): string {
        // Build atom data for each step
        const stepsJson = JSON.stringify(data.steps.map(step => ({
            energy: step.energy,
            maxForce: step.maxForce,
            positions: step.positions
        })));

        // Build element info for coloring
        let atomElements: string[] = [];
        for (let i = 0; i < data.elements.length; i++) {
            for (let j = 0; j < data.atomCounts[i]; j++) {
                atomElements.push(data.elements[i]);
            }
        }

        const energies = data.steps.map(s => s.energy);
        const minEnergy = Math.min(...energies);
        const maxEnergy = Math.max(...energies);

        return `<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <meta http-equiv="Content-Security-Policy" content="default-src 'none'; script-src https://cdn.jsdelivr.net https://3dmol.org 'unsafe-inline' 'unsafe-eval'; style-src 'unsafe-inline'; img-src https: data:;">
    <title>Relaxation Progress</title>
    <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
    <script src="https://3dmol.org/build/3Dmol-min.js"></script>
    <style>
        body {
            font-family: var(--vscode-font-family);
            color: var(--vscode-foreground);
            background-color: var(--vscode-editor-background);
            margin: 0;
            padding: 15px;
            display: flex;
            flex-direction: column;
            height: 100vh;
            box-sizing: border-box;
        }
        h1 {
            margin: 0 0 10px 0;
            font-size: 1.3em;
            color: var(--vscode-textLink-foreground);
        }
        .info-row {
            display: flex;
            gap: 20px;
            margin-bottom: 10px;
            font-size: 0.9em;
        }
        .info-item {
            display: flex;
            gap: 5px;
        }
        .info-label {
            font-weight: bold;
            color: var(--vscode-symbolIcon-propertyForeground);
        }
        .main-content {
            display: flex;
            flex: 1;
            gap: 15px;
            min-height: 0;
        }
        .left-panel {
            flex: 1;
            display: flex;
            flex-direction: column;
            min-width: 0;
        }
        .right-panel {
            flex: 1;
            display: flex;
            flex-direction: column;
            min-width: 0;
        }
        .chart-container {
            flex: 1;
            min-height: 200px;
            position: relative;
        }
        .slider-container {
            padding: 10px 0;
            display: flex;
            align-items: center;
            gap: 10px;
        }
        .slider-container label {
            font-weight: bold;
            min-width: 50px;
        }
        .slider-container input[type="range"] {
            flex: 1;
        }
        .slider-value {
            min-width: 80px;
            text-align: right;
            font-family: var(--vscode-editor-font-family);
        }
        #viewer {
            flex: 1;
            min-height: 300px;
            border: 1px solid var(--vscode-widget-border);
            border-radius: 5px;
        }
        .step-info {
            padding: 10px;
            background-color: var(--vscode-editor-inactiveSelectionBackground);
            border-radius: 5px;
            margin-top: 10px;
            font-family: var(--vscode-editor-font-family);
        }
        .step-info .value {
            font-weight: bold;
        }
    </style>
</head>
<body>
    <h1>${this._escapeHtml(data.formula)} - Relaxation Progress</h1>
    <div class="info-row">
        <div class="info-item">
            <span class="info-label">Steps:</span>
            <span>${data.steps.length}</span>
        </div>
        <div class="info-item">
            <span class="info-label">Energy range:</span>
            <span>${minEnergy.toFixed(4)} to ${maxEnergy.toFixed(4)} eV</span>
        </div>
        <div class="info-item">
            <span class="info-label">ΔE:</span>
            <span>${(maxEnergy - minEnergy).toFixed(4)} eV</span>
        </div>
    </div>

    <div class="slider-container">
        <label>Step:</label>
        <input type="range" id="stepSlider" min="0" max="${data.steps.length - 1}" value="${data.steps.length - 1}">
        <span class="slider-value" id="sliderValue">${data.steps.length} / ${data.steps.length}</span>
    </div>

    <div class="main-content">
        <div class="left-panel">
            <div class="chart-container">
                <canvas id="energyChart"></canvas>
            </div>
            <div class="step-info" id="stepInfo">
                <div><span class="info-label">Energy:</span> <span class="value" id="currentEnergy">${data.steps[data.steps.length - 1].energy.toFixed(6)}</span> eV</div>
                <div><span class="info-label">Max Force:</span> <span class="value" id="currentForce">${data.steps[data.steps.length - 1].maxForce?.toFixed(4) || 'N/A'}</span> eV/Å</div>
            </div>
        </div>
        <div class="right-panel">
            <div id="viewer"></div>
        </div>
    </div>

    <script>
        const steps = ${stepsJson};
        const lattice = ${JSON.stringify(data.lattice)};
        const atomElements = ${JSON.stringify(atomElements)};
        const cpkColors = ${JSON.stringify(CPK_COLORS)};

        let currentStep = steps.length - 1;
        let viewer;
        let chart;
        let highlightPoint = null;

        // Convert fractional to Cartesian
        function fracToCart(frac) {
            return [
                frac[0] * lattice[0][0] + frac[1] * lattice[1][0] + frac[2] * lattice[2][0],
                frac[0] * lattice[0][1] + frac[1] * lattice[1][1] + frac[2] * lattice[2][1],
                frac[0] * lattice[0][2] + frac[1] * lattice[1][2] + frac[2] * lattice[2][2]
            ];
        }

        // Initialize chart
        function initChart() {
            const ctx = document.getElementById('energyChart').getContext('2d');
            const energyData = steps.map((s, i) => ({ x: i + 1, y: s.energy }));

            chart = new Chart(ctx, {
                type: 'line',
                data: {
                    datasets: [{
                        label: 'Energy (eV)',
                        data: energyData,
                        borderColor: '#4ec9b0',
                        backgroundColor: 'rgba(78, 201, 176, 0.1)',
                        fill: true,
                        tension: 0.1,
                        pointRadius: 3,
                        pointHoverRadius: 6
                    }, {
                        label: 'Current',
                        data: [{ x: currentStep + 1, y: steps[currentStep].energy }],
                        borderColor: '#f14c4c',
                        backgroundColor: '#f14c4c',
                        pointRadius: 8,
                        pointHoverRadius: 10,
                        showLine: false
                    }]
                },
                options: {
                    responsive: true,
                    maintainAspectRatio: false,
                    scales: {
                        x: {
                            type: 'linear',
                            title: { display: true, text: 'Ionic Step', color: '#ccc' },
                            ticks: { color: '#ccc' },
                            grid: { color: 'rgba(255,255,255,0.1)' }
                        },
                        y: {
                            title: { display: true, text: 'Energy (eV)', color: '#ccc' },
                            ticks: { color: '#ccc' },
                            grid: { color: 'rgba(255,255,255,0.1)' }
                        }
                    },
                    plugins: {
                        legend: { display: false }
                    },
                    onClick: (e) => {
                        const points = chart.getElementsAtEventForMode(e, 'nearest', { intersect: true }, false);
                        if (points.length > 0 && points[0].datasetIndex === 0) {
                            const idx = points[0].index;
                            document.getElementById('stepSlider').value = idx;
                            updateStep(idx);
                        }
                    }
                }
            });
        }

        // Helper to determine if background is light or dark
        function isLightBackground(color) {
            let r, g, b;
            if (color.startsWith('#')) {
                const hex = color.slice(1);
                r = parseInt(hex.substr(0, 2), 16);
                g = parseInt(hex.substr(2, 2), 16);
                b = parseInt(hex.substr(4, 2), 16);
            } else if (color.startsWith('rgb')) {
                const match = color.match(/\\d+/g);
                if (match) {
                    r = parseInt(match[0]);
                    g = parseInt(match[1]);
                    b = parseInt(match[2]);
                }
            } else {
                return false;
            }
            const luminance = (0.299 * r + 0.587 * g + 0.114 * b) / 255;
            return luminance > 0.5;
        }

        let unitCellColor = 'yellow';

        // Initialize 3D viewer
        function initViewer() {
            const element = document.getElementById('viewer');
            // Get theme background color - this will adapt to light/dark themes
            const bgColor = getComputedStyle(document.body).getPropertyValue('--vscode-editor-background').trim() || 'black';
            unitCellColor = isLightBackground(bgColor) ? '#555555' : 'yellow';
            viewer = $3Dmol.createViewer(element, { backgroundColor: bgColor });
            updateStructure();
            viewer.zoomTo();
            viewer.render();
        }

        // Update structure for current step
        function updateStructure() {
            viewer.removeAllShapes();
            viewer.removeAllModels();

            const positions = steps[currentStep].positions;

            // Add atoms
            positions.forEach((frac, i) => {
                const cart = fracToCart(frac);
                const elem = atomElements[i];
                const color = cpkColors[elem] || '#FF00FF';

                viewer.addSphere({
                    center: { x: cart[0], y: cart[1], z: cart[2] },
                    radius: 0.5,
                    color: color
                });
            });

            // Add unit cell
            const a = lattice[0];
            const b = lattice[1];
            const c = lattice[2];
            const corners = [
                [0, 0, 0], a, b, c,
                [a[0] + b[0], a[1] + b[1], a[2] + b[2]],
                [a[0] + c[0], a[1] + c[1], a[2] + c[2]],
                [b[0] + c[0], b[1] + c[1], b[2] + c[2]],
                [a[0] + b[0] + c[0], a[1] + b[1] + c[1], a[2] + b[2] + c[2]]
            ];
            const edges = [[0,1],[0,2],[0,3],[1,4],[1,5],[2,4],[2,6],[3,5],[3,6],[4,7],[5,7],[6,7]];

            edges.forEach(([i, j]) => {
                viewer.addCylinder({
                    start: { x: corners[i][0], y: corners[i][1], z: corners[i][2] },
                    end: { x: corners[j][0], y: corners[j][1], z: corners[j][2] },
                    radius: 0.03,
                    color: unitCellColor,
                    fromCap: true,
                    toCap: true
                });
            });

            viewer.render();
        }

        // Update step
        function updateStep(idx) {
            currentStep = idx;

            // Update slider display
            document.getElementById('sliderValue').textContent = (idx + 1) + ' / ' + steps.length;

            // Update info
            document.getElementById('currentEnergy').textContent = steps[idx].energy.toFixed(6);
            document.getElementById('currentForce').textContent = steps[idx].maxForce?.toFixed(4) || 'N/A';

            // Update chart highlight point
            chart.data.datasets[1].data = [{ x: idx + 1, y: steps[idx].energy }];
            chart.update('none');

            // Update structure
            updateStructure();
        }

        // Slider event
        document.getElementById('stepSlider').addEventListener('input', (e) => {
            updateStep(parseInt(e.target.value));
        });

        // Initialize
        initChart();
        initViewer();
    </script>
</body>
</html>`;
    }

    private _escapeHtml(text: string): string {
        return text
            .replace(/&/g, '&amp;')
            .replace(/</g, '&lt;')
            .replace(/>/g, '&gt;')
            .replace(/"/g, '&quot;')
            .replace(/'/g, '&#039;');
    }
}

// Unified VaspSumPanel - combines structure preview and calculation summary
class VaspSumPanel {
    public static currentPanel: VaspSumPanel | undefined;
    private readonly _panel: vscode.WebviewPanel;
    private _disposables: vscode.Disposable[] = [];
    private _directory: string;

    public static async createOrShow(context: vscode.ExtensionContext) {
        const column = vscode.ViewColumn.Beside;

        // Get directory from active editor or workspace
        const editor = vscode.window.activeTextEditor;
        let directory: string | undefined;

        if (editor) {
            directory = path.dirname(editor.document.uri.fsPath);
        } else if (vscode.workspace.workspaceFolders && vscode.workspace.workspaceFolders.length > 0) {
            directory = vscode.workspace.workspaceFolders[0].uri.fsPath;
        }

        if (!directory) {
            vscode.window.showErrorMessage('No active file or workspace to determine VASP directory');
            return;
        }

        if (VaspSumPanel.currentPanel) {
            VaspSumPanel.currentPanel._panel.reveal(column);
            VaspSumPanel.currentPanel._update(directory);
            return;
        }

        const panel = vscode.window.createWebviewPanel(
            'vaspsum',
            'VASPsum',
            column,
            {
                enableScripts: true,
                retainContextWhenHidden: true
            }
        );

        VaspSumPanel.currentPanel = new VaspSumPanel(panel, directory);
    }

    private constructor(panel: vscode.WebviewPanel, directory: string) {
        this._panel = panel;
        this._directory = directory;
        this._update(directory);

        this._panel.onDidDispose(() => this.dispose(), null, this._disposables);
    }

    public dispose() {
        VaspSumPanel.currentPanel = undefined;
        this._panel.dispose();
        while (this._disposables.length) {
            const x = this._disposables.pop();
            if (x) {
                x.dispose();
            }
        }
    }

    private _update(directory: string) {
        this._directory = directory;
        this._panel.title = 'VASPsum: ' + path.basename(directory);

        // Check what files are available
        const poscarPath = fs.existsSync(path.join(directory, 'CONTCAR'))
            ? path.join(directory, 'CONTCAR')
            : fs.existsSync(path.join(directory, 'POSCAR'))
                ? path.join(directory, 'POSCAR')
                : null;
        const vasprunPath = fs.existsSync(path.join(directory, 'vasprun.xml'))
            ? path.join(directory, 'vasprun.xml')
            : null;

        // Parse structure if available
        let poscarData: PoscarData | null = null;
        if (poscarPath) {
            try {
                const content = fs.readFileSync(poscarPath, 'utf-8');
                poscarData = parsePoscar(content);
            } catch (error) {
                console.error('Error reading POSCAR:', error);
            }
        }

        // Parse vasprun.xml if available
        let summaryData: VaspSummaryData | null = null;
        let relaxationData: RelaxationData | null = null;
        if (vasprunPath) {
            try {
                const content = fs.readFileSync(vasprunPath, 'utf-8');
                summaryData = parseVasprunXml(content);
                if (summaryData) {
                    summaryData.directory = directory;
                }
                // Check for relaxation (multiple ionic steps)
                relaxationData = parseRelaxationSteps(content);
                if (relaxationData && relaxationData.steps.length <= 1) {
                    relaxationData = null; // Only show relaxation viewer for multi-step
                }
            } catch (error) {
                console.error('Error reading vasprun.xml:', error);
            }
        }

        if (!poscarData && !summaryData) {
            this._panel.webview.html = this._getErrorHtml('No POSCAR/CONTCAR or vasprun.xml found in ' + directory);
            return;
        }

        this._panel.webview.html = this._getHtmlForWebview(poscarData, summaryData, relaxationData, directory);
    }

    private _getErrorHtml(message: string): string {
        return `<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <style>
        body {
            font-family: var(--vscode-font-family);
            color: var(--vscode-foreground);
            background-color: var(--vscode-editor-background);
            padding: 20px;
        }
        .error { color: var(--vscode-errorForeground); }
    </style>
</head>
<body>
    <h2 class="error">Error</h2>
    <p>${message}</p>
</body>
</html>`;
    }

    private _getHtmlForWebview(poscarData: PoscarData | null, summaryData: VaspSummaryData | null, relaxationData: RelaxationData | null, directory: string): string {
        // Build structure section if available
        let structureHtml = '';
        let structureScript = '';
        let atoms: { elem: string; x: number; y: number; z: number; color: string; constrained: boolean }[] = [];
        let scaledLattice: number[][] = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];

        if (poscarData) {
            const cartesian = toCartesian(poscarData);
            const volume = cellVolume(poscarData);
            const formula = getFormula(poscarData);
            const totalAtoms = poscarData.counts.reduce((a, b) => a + b, 0);

            // Build atom data
            let atomIndex = 0;
            for (let i = 0; i < poscarData.elements.length; i++) {
                const element = poscarData.elements[i];
                const count = poscarData.counts[i];
                const color = CPK_COLORS[element] || '#FF00FF';

                for (let j = 0; j < count; j++) {
                    const pos = cartesian[atomIndex];
                    const constraint = poscarData.constraints[atomIndex] || [false, false, false];
                    const isConstrained = constraint.some(c => c);
                    atomIndex++;
                    atoms.push({
                        elem: element,
                        x: pos[0],
                        y: pos[1],
                        z: pos[2],
                        color: color,
                        constrained: isConstrained
                    });
                }
            }

            scaledLattice = poscarData.lattice.map(row => row.map(x => x * poscarData.scale));

            structureHtml = `
            <div class="structure-section">
                <div class="info-panel">
                    <h2>${this._escapeHtml(poscarData.title)}</h2>
                    <div class="info-grid">
                        <span class="info-label">Formula:</span>
                        <span class="info-value">${this._escapeHtml(formula)}</span>
                        <span class="info-label">Atoms:</span>
                        <span class="info-value">${totalAtoms}</span>
                        <span class="info-label">Volume:</span>
                        <span class="info-value">${volume.toFixed(2)} &#8491;&sup3;</span>
                        <span class="info-label">Coordinates:</span>
                        <span class="info-value">${poscarData.isDirect ? 'Direct (fractional)' : 'Cartesian'}</span>
                        <span class="info-label">Constraints:</span>
                        <span class="info-value">${poscarData.selectiveDynamics ? atoms.filter(a => a.constrained).length + ' constrained (shown with X)' : 'None'}</span>
                    </div>
                </div>

                <div class="controls">
                    <button id="resetView">Reset View</button>
                    <button id="toggleCell" class="active">Unit Cell</button>
                    <button id="toggleLabels">Labels</button>
                    <button id="styleStick">Stick</button>
                    <button id="styleSphere" class="active">Sphere</button>
                    <div class="repeat-controls">
                        <label>Repeat:</label>
                        <div class="repeat-group">
                            <span>a</span>
                            <input type="number" id="repeatA" class="repeat-input" value="1" min="1" max="10">
                            <div class="spinner">
                                <button onclick="incrementRepeat('A')">&#9650;</button>
                                <button onclick="decrementRepeat('A')">&#9660;</button>
                            </div>
                        </div>
                        <div class="repeat-group">
                            <span>b</span>
                            <input type="number" id="repeatB" class="repeat-input" value="1" min="1" max="10">
                            <div class="spinner">
                                <button onclick="incrementRepeat('B')">&#9650;</button>
                                <button onclick="decrementRepeat('B')">&#9660;</button>
                            </div>
                        </div>
                        <div class="repeat-group">
                            <span>c</span>
                            <input type="number" id="repeatC" class="repeat-input" value="1" min="1" max="10">
                            <div class="spinner">
                                <button onclick="incrementRepeat('C')">&#9650;</button>
                                <button onclick="decrementRepeat('C')">&#9660;</button>
                            </div>
                        </div>
                    </div>
                </div>

                <div id="viewer"></div>

                <div class="legend">
                    ${poscarData.elements.map((el, i) => `
                        <div class="legend-item">
                            <div class="legend-color" style="background-color: ${CPK_COLORS[el] || '#FF00FF'};"></div>
                            <span>${el} (${poscarData.counts[i]})</span>
                        </div>
                    `).join('')}
                    ${poscarData.selectiveDynamics ? `
                        <div class="legend-item legend-separator">
                            <div class="legend-marker">✱</div>
                            <span>Constrained atom</span>
                        </div>
                    ` : ''}
                    <div class="legend-item legend-separator">
                        <button id="savePng" title="Save as PNG">Save PNG</button>
                    </div>
                </div>
            </div>`;
        }

        // Build summary section if available
        let summaryHtml = '';
        if (summaryData) {
            const convergenceStatus = summaryData.converged === null ? 'Unknown' :
                (summaryData.converged ? '✓ Converged' : '✗ Not Converged');
            const convergenceClass = summaryData.converged === null ? 'unknown' :
                (summaryData.converged ? 'converged' : 'not-converged');

            const formatEnergy = (e: number | null) => e !== null ? e.toFixed(6) : 'N/A';
            const formatForce = (f: number | null) => f !== null ? f.toFixed(6) : 'N/A';
            const formatPressure = (p: number | null) => p !== null ? p.toFixed(3) : 'N/A';

            // Build stress tensor display
            let stressHtml = 'N/A';
            if (summaryData.stress && summaryData.stress.length >= 3) {
                stressHtml = `<table class="matrix">
                    ${summaryData.stress.map(row => `<tr>${row.map(v => `<td>${v.toFixed(2)}</td>`).join('')}</tr>`).join('')}
                </table>`;
            }

            // Build forces summary
            let forcesHtml = '';
            if (summaryData.forces && summaryData.elements.length > 0) {
                let atomIdx = 0;
                for (let i = 0; i < summaryData.elements.length; i++) {
                    const count = summaryData.atomCounts[i] || 1;
                    for (let j = 0; j < count && atomIdx < summaryData.forces.length; j++) {
                        const f = summaryData.forces[atomIdx];
                        const mag = Math.sqrt(f[0] * f[0] + f[1] * f[1] + f[2] * f[2]);
                        forcesHtml += `<tr>
                            <td>${atomIdx + 1}</td>
                            <td>${summaryData.elements[i]}</td>
                            <td>${f[0].toFixed(4)}</td>
                            <td>${f[1].toFixed(4)}</td>
                            <td>${f[2].toFixed(4)}</td>
                            <td class="${mag > 0.05 ? 'highlight' : ''}">${mag.toFixed(4)}</td>
                        </tr>`;
                        atomIdx++;
                    }
                }
            }

            summaryHtml = `
            <div class="summary-section">
                <h1>Calculation Summary</h1>
                <div class="directory">${this._escapeHtml(directory)}</div>

                <div class="section">
                    <div class="grid">
                        <span class="label">Total Energy:</span>
                        <span class="value energy">${formatEnergy(summaryData.energy)} eV</span>
                        <span class="label">Energy/Atom:</span>
                        <span class="value">${formatEnergy(summaryData.energyPerAtom)} eV</span>
                        <span class="label">Fermi Energy:</span>
                        <span class="value">${formatEnergy(summaryData.fermiEnergy)} eV</span>
                        <span class="label">Convergence:</span>
                        <span class="value ${convergenceClass}">${convergenceStatus}</span>
                    </div>
                </div>

                <div class="section">
                    <h2>Forces & Stress</h2>
                    <div class="grid">
                        <span class="label">Max Force:</span>
                        <span class="value ${summaryData.maxForce !== null && summaryData.maxForce > 0.05 ? 'highlight' : ''}">${formatForce(summaryData.maxForce)} eV/Å</span>
                        <span class="label">Pressure:</span>
                        <span class="value">${formatPressure(summaryData.pressure)} GPa</span>
                    </div>
                    ${summaryData.stress ? `
                    <h3 class="collapsible" onclick="toggleCollapse(this)">Stress Tensor (kBar)</h3>
                    <div class="collapse-content">${stressHtml}</div>
                    ` : ''}
                    ${summaryData.forces && summaryData.forces.length > 0 ? `
                    <h3 class="collapsible" onclick="toggleCollapse(this)">Atomic Forces (eV/Å)</h3>
                    <div class="collapse-content">
                        <table class="forces">
                            <tr><th>#</th><th>Element</th><th>Fx</th><th>Fy</th><th>Fz</th><th>|F|</th></tr>
                            ${forcesHtml}
                        </table>
                    </div>
                    ` : ''}
                </div>

                <div class="section">
                    <h2>Calculation Parameters</h2>
                    <div class="params-grid">
                        <div class="param-item"><span class="label">XC:</span><span class="value">${this._escapeHtml(summaryData.xc)}</span></div>
                        <div class="param-item"><span class="label">ENCUT:</span><span class="value">${summaryData.encut !== null ? summaryData.encut.toFixed(0) : 'N/A'} eV</span></div>
                        <div class="param-item"><span class="label">K-points:</span><span class="value">${this._escapeHtml(summaryData.kpoints)}</span></div>
                        <div class="param-item"><span class="label">ISPIN:</span><span class="value">${summaryData.ispin !== null ? summaryData.ispin : 'N/A'}</span></div>
                        <div class="param-item"><span class="label">IBRION:</span><span class="value">${summaryData.ibrion !== null ? summaryData.ibrion : 'N/A'}</span></div>
                        <div class="param-item"><span class="label">NSW:</span><span class="value">${summaryData.nsw !== null ? summaryData.nsw : 'N/A'}</span></div>
                        <div class="param-item"><span class="label">EDIFF:</span><span class="value">${summaryData.ediff !== null ? summaryData.ediff.toExponential(1) : 'N/A'}</span></div>
                        <div class="param-item"><span class="label">EDIFFG:</span><span class="value">${summaryData.ediffg !== null ? summaryData.ediffg.toExponential(1) : 'N/A'}</span></div>
                    </div>
                </div>

                <div class="section">
                    <h2>Iteration Summary</h2>
                    <div class="grid">
                        <span class="label">Ionic Steps:</span>
                        <span class="value">${summaryData.ionicSteps !== null ? summaryData.ionicSteps : 'N/A'}</span>
                        <span class="label">Electronic Steps:</span>
                        <span class="value">${summaryData.electronicSteps !== null ? summaryData.electronicSteps : 'N/A'}</span>
                    </div>
                </div>
            </div>`;
        }

        // Build relaxation section if available
        let relaxationHtml = '';
        let relaxationScript = '';
        if (relaxationData && relaxationData.steps.length > 1) {
            relaxationHtml = `
            <div class="relaxation-section">
                <h2>Relaxation Progress (${relaxationData.steps.length} ionic steps)</h2>
                <div class="chart-container">
                    <canvas id="energyChart"></canvas>
                </div>
                <div class="slider-container">
                    <label>Ionic Step:</label>
                    <input type="range" id="stepSlider" min="0" max="${relaxationData.steps.length - 1}" value="${relaxationData.steps.length - 1}">
                    <span id="sliderValue">${relaxationData.steps.length} / ${relaxationData.steps.length}</span>
                </div>
                <div class="step-info">
                    <span class="label">Energy:</span>
                    <span id="currentEnergy">${relaxationData.steps[relaxationData.steps.length - 1].energy.toFixed(6)}</span> eV
                    <span class="label" style="margin-left: 20px;">Max Force:</span>
                    <span id="currentForce">${relaxationData.steps[relaxationData.steps.length - 1].maxForce?.toFixed(4) || 'N/A'}</span> eV/Å
                </div>
            </div>`;

            relaxationScript = `
            // Relaxation data and chart
            const relaxSteps = ${JSON.stringify(relaxationData.steps.map(s => ({ energy: s.energy, maxForce: s.maxForce, positions: s.positions })))};
            const relaxLattice = ${JSON.stringify(relaxationData.lattice)};
            const relaxElements = ${JSON.stringify(relaxationData.elements)};
            const relaxAtomCounts = ${JSON.stringify(relaxationData.atomCounts)};
            let chart;
            let currentRelaxStep = relaxSteps.length - 1;

            // Build element list for each atom
            const relaxAtomElements = [];
            relaxElements.forEach((el, i) => {
                for (let j = 0; j < relaxAtomCounts[i]; j++) {
                    relaxAtomElements.push(el);
                }
            });

            // Convert fractional to Cartesian coordinates
            function fracToCart(frac) {
                return [
                    frac[0] * relaxLattice[0][0] + frac[1] * relaxLattice[1][0] + frac[2] * relaxLattice[2][0],
                    frac[0] * relaxLattice[0][1] + frac[1] * relaxLattice[1][1] + frac[2] * relaxLattice[2][1],
                    frac[0] * relaxLattice[0][2] + frac[1] * relaxLattice[1][2] + frac[2] * relaxLattice[2][2]
                ];
            }

            // Update structure viewer with positions from current step
            function updateStructureFromRelax() {
                if (typeof viewer === 'undefined' || !viewer) return;

                viewer.removeAllModels();
                viewer.removeAllLabels();
                viewer.removeAllShapes();

                const positions = relaxSteps[currentRelaxStep].positions;

                // Add atoms
                positions.forEach((frac, i) => {
                    const cart = fracToCart(frac);
                    const elem = relaxAtomElements[i];
                    const color = baseAtoms[i]?.color || '#FF00FF';
                    const constrained = baseAtoms[i]?.constrained || false;
                    const radius = styleMode === 'sphere' ? 0.5 : 0.2;
                    const opacity = constrained ? 0.9 : 1.0;

                    viewer.addSphere({
                        center: { x: cart[0], y: cart[1], z: cart[2] },
                        radius: radius,
                        color: color,
                        opacity: opacity
                    });

                    // Add X marker for constrained atoms (X in XY plane + vertical line)
                    if (constrained) {
                        const size = radius * 1.2;
                        viewer.addCylinder({
                            start: { x: cart[0] - size, y: cart[1] - size, z: cart[2] },
                            end: { x: cart[0] + size, y: cart[1] + size, z: cart[2] },
                            radius: 0.05,
                            color: 'red',
                            fromCap: true,
                            toCap: true
                        });
                        viewer.addCylinder({
                            start: { x: cart[0] + size, y: cart[1] - size, z: cart[2] },
                            end: { x: cart[0] - size, y: cart[1] + size, z: cart[2] },
                            radius: 0.05,
                            color: 'red',
                            fromCap: true,
                            toCap: true
                        });
                        // Add vertical line (Z direction) for visibility from all angles
                        viewer.addCylinder({
                            start: { x: cart[0], y: cart[1], z: cart[2] - size },
                            end: { x: cart[0], y: cart[1], z: cart[2] + size },
                            radius: 0.05,
                            color: 'red',
                            fromCap: true,
                            toCap: true
                        });
                    }

                    if (showLabels) {
                        viewer.addLabel(elem, {
                            position: { x: cart[0], y: cart[1], z: cart[2] },
                            fontSize: 12,
                            fontColor: 'white',
                            backgroundOpacity: 0.5
                        });
                    }
                });

                // Re-add unit cell
                addUnitCell();
                viewer.render();
            }

            function initRelaxChart() {
                const ctx = document.getElementById('energyChart').getContext('2d');
                const energyData = relaxSteps.map((s, i) => ({ x: i + 1, y: s.energy }));

                chart = new Chart(ctx, {
                    type: 'line',
                    data: {
                        datasets: [
                            {
                                label: 'Energy (eV)',
                                data: energyData,
                                borderColor: '#4ec9b0',
                                backgroundColor: 'rgba(78, 201, 176, 0.1)',
                                fill: true,
                                tension: 0.1,
                                pointRadius: 3
                            },
                            {
                                label: 'Current',
                                data: [{ x: relaxSteps.length, y: relaxSteps[relaxSteps.length - 1].energy }],
                                borderColor: 'red',
                                backgroundColor: 'red',
                                pointRadius: 8,
                                pointStyle: 'circle',
                                showLine: false
                            }
                        ]
                    },
                    options: {
                        responsive: true,
                        maintainAspectRatio: false,
                        plugins: {
                            legend: { display: false }
                        },
                        scales: {
                            x: {
                                type: 'linear',
                                title: { display: true, text: 'Ionic Step' },
                                ticks: { stepSize: 1 }
                            },
                            y: {
                                title: { display: true, text: 'Energy (eV)' }
                            }
                        }
                    }
                });
            }

            function updateRelaxStep(idx) {
                currentRelaxStep = idx;
                document.getElementById('sliderValue').textContent = (idx + 1) + ' / ' + relaxSteps.length;
                document.getElementById('currentEnergy').textContent = relaxSteps[idx].energy.toFixed(6);
                document.getElementById('currentForce').textContent = relaxSteps[idx].maxForce?.toFixed(4) || 'N/A';

                chart.data.datasets[1].data = [{ x: idx + 1, y: relaxSteps[idx].energy }];
                chart.update('none');

                // Update structure viewer with new positions
                updateStructureFromRelax();
            }

            document.getElementById('stepSlider').addEventListener('input', (e) => {
                updateRelaxStep(parseInt(e.target.value));
            });
            `;
        }

        return `<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <meta http-equiv="Content-Security-Policy" content="default-src 'none'; script-src https://3dmol.org https://cdn.jsdelivr.net 'unsafe-inline' 'unsafe-eval'; style-src 'unsafe-inline'; img-src https: data:;">
    <title>VASPsum</title>
    ${poscarData ? '<script src="https://3dmol.org/build/3Dmol-min.js"></script>' : ''}
    ${relaxationData ? '<script src="https://cdn.jsdelivr.net/npm/chart.js"></script>' : ''}
    <style>
        body {
            font-family: var(--vscode-font-family);
            font-size: var(--vscode-font-size);
            color: var(--vscode-foreground);
            background-color: var(--vscode-editor-background);
            margin: 0;
            padding: 10px;
            overflow-y: auto;
        }
        h1 {
            color: var(--vscode-textLink-foreground);
            font-size: 1.4em;
            margin-bottom: 5px;
            border-bottom: 1px solid var(--vscode-widget-border);
            padding-bottom: 10px;
        }
        h2 {
            color: var(--vscode-textLink-foreground);
            font-size: 1.1em;
            margin-top: 20px;
            margin-bottom: 10px;
        }
        h3 { font-size: 1em; margin: 10px 0 5px 0; }
        .directory {
            font-size: 0.85em;
            color: var(--vscode-descriptionForeground);
            margin-bottom: 15px;
        }
        .structure-section {
            margin-bottom: 20px;
        }
        .summary-section, .relaxation-section {
            border-top: 1px solid var(--vscode-widget-border);
            padding-top: 15px;
            margin-top: 15px;
        }
        .info-panel {
            background-color: var(--vscode-editor-inactiveSelectionBackground);
            padding: 10px 15px;
            border-radius: 5px;
            margin-bottom: 10px;
            position: relative;
            z-index: 10;
        }
        .info-panel h2 {
            margin: 0 0 8px 0;
            font-size: 1.1em;
            color: var(--vscode-textLink-foreground);
        }
        .info-grid {
            display: grid;
            grid-template-columns: auto 1fr;
            gap: 4px 15px;
        }
        .info-label {
            font-weight: bold;
            color: var(--vscode-symbolIcon-propertyForeground);
        }
        .info-value {
            font-family: var(--vscode-editor-font-family);
        }
        .controls {
            display: flex;
            gap: 10px;
            margin-bottom: 10px;
            flex-wrap: wrap;
            position: relative;
            z-index: 10;
        }
        button {
            background-color: var(--vscode-button-background);
            color: var(--vscode-button-foreground);
            border: none;
            padding: 6px 12px;
            cursor: pointer;
            border-radius: 3px;
            font-size: 12px;
        }
        button:hover {
            background-color: var(--vscode-button-hoverBackground);
        }
        button.active {
            background-color: var(--vscode-button-hoverBackground);
            outline: 1px solid var(--vscode-focusBorder);
        }
        #viewer {
            height: 350px;
            min-height: 300px;
            border: 1px solid var(--vscode-widget-border);
            border-radius: 5px;
            position: relative;
            z-index: 1;
        }
        .legend {
            display: flex;
            flex-wrap: wrap;
            gap: 10px;
            margin-top: 10px;
            padding: 8px;
            background-color: var(--vscode-editor-inactiveSelectionBackground);
            border-radius: 5px;
            position: relative;
            z-index: 10;
        }
        .legend-item {
            display: flex;
            align-items: center;
            gap: 5px;
        }
        .legend-color {
            width: 16px;
            height: 16px;
            border-radius: 50%;
            border: 1px solid var(--vscode-widget-border);
        }
        .legend-marker {
            color: red;
            font-size: 16px;
            font-weight: bold;
            width: 16px;
            text-align: center;
        }
        .legend-separator {
            margin-left: 10px;
            padding-left: 10px;
            border-left: 1px solid var(--vscode-widget-border);
        }
        .repeat-controls {
            display: flex;
            align-items: center;
            gap: 15px;
            margin-left: auto;
            padding-left: 15px;
            border-left: 1px solid var(--vscode-widget-border);
        }
        .repeat-controls label {
            font-size: 12px;
            color: var(--vscode-descriptionForeground);
        }
        .repeat-group {
            display: flex;
            align-items: center;
            gap: 3px;
        }
        .repeat-group span {
            font-size: 12px;
            font-weight: bold;
            min-width: 12px;
        }
        .spinner {
            display: flex;
            flex-direction: column;
        }
        .spinner button {
            padding: 0;
            width: 18px;
            height: 12px;
            font-size: 8px;
            line-height: 1;
            display: flex;
            align-items: center;
            justify-content: center;
            border-radius: 2px;
        }
        .repeat-input {
            width: 30px;
            height: 24px;
            text-align: center;
            background-color: var(--vscode-input-background);
            color: var(--vscode-input-foreground);
            border: 1px solid var(--vscode-input-border);
            border-radius: 2px;
            font-size: 12px;
        }
        .section {
            background-color: var(--vscode-editor-inactiveSelectionBackground);
            padding: 12px 15px;
            border-radius: 5px;
            margin-bottom: 15px;
        }
        .grid {
            display: grid;
            grid-template-columns: 160px 1fr;
            gap: 8px 15px;
        }
        .label {
            font-weight: bold;
            color: var(--vscode-symbolIcon-propertyForeground);
        }
        .value {
            font-family: var(--vscode-editor-font-family);
        }
        .energy {
            font-size: 1.1em;
            font-weight: bold;
        }
        .converged { color: #4ec9b0; }
        .not-converged { color: #f14c4c; }
        .unknown { color: var(--vscode-descriptionForeground); }
        .highlight { color: #f14c4c; font-weight: bold; }
        table {
            border-collapse: collapse;
            font-family: var(--vscode-editor-font-family);
            font-size: 0.9em;
        }
        table.matrix td {
            padding: 3px 10px;
            text-align: right;
        }
        table.forces {
            width: 100%;
            margin-top: 10px;
        }
        table.forces th, table.forces td {
            padding: 5px 8px;
            text-align: right;
            border-bottom: 1px solid var(--vscode-widget-border);
        }
        table.forces th {
            background-color: var(--vscode-editor-background);
            font-weight: bold;
        }
        .params-grid {
            display: grid;
            grid-template-columns: repeat(auto-fill, minmax(180px, 1fr));
            gap: 8px;
        }
        .param-item {
            display: flex;
            gap: 8px;
        }
        .collapsible {
            cursor: pointer;
            user-select: none;
        }
        .collapsible:before {
            content: '▶ ';
            font-size: 0.8em;
        }
        .collapsible.active:before {
            content: '▼ ';
        }
        .collapse-content {
            display: none;
            margin-top: 10px;
        }
        .collapse-content.show {
            display: block;
        }
        .chart-container {
            height: 200px;
            margin: 10px 0;
        }
        .slider-container {
            display: flex;
            align-items: center;
            gap: 10px;
            margin: 10px 0;
        }
        .slider-container input[type="range"] {
            flex: 1;
        }
        .step-info {
            font-family: var(--vscode-editor-font-family);
            padding: 8px;
            background-color: var(--vscode-editor-inactiveSelectionBackground);
            border-radius: 5px;
        }
    </style>
</head>
<body>
    ${structureHtml}
    ${relaxationHtml}
    ${summaryHtml}

    <script>
        function toggleCollapse(element) {
            element.classList.toggle('active');
            const content = element.nextElementSibling;
            content.classList.toggle('show');
        }

        ${poscarData ? `
        // Structure viewer code
        const baseAtoms = ${JSON.stringify(atoms)};
        const lattice = ${JSON.stringify(scaledLattice)};

        let viewer;
        let showCell = true;
        let showLabels = false;
        let styleMode = 'sphere';
        let repeatA = 1, repeatB = 1, repeatC = 1;

        function getRepeatedAtoms() {
            const repeated = [];
            for (let ia = 0; ia < repeatA; ia++) {
                for (let ib = 0; ib < repeatB; ib++) {
                    for (let ic = 0; ic < repeatC; ic++) {
                        const tx = ia * lattice[0][0] + ib * lattice[1][0] + ic * lattice[2][0];
                        const ty = ia * lattice[0][1] + ib * lattice[1][1] + ic * lattice[2][1];
                        const tz = ia * lattice[0][2] + ib * lattice[1][2] + ic * lattice[2][2];

                        baseAtoms.forEach(atom => {
                            repeated.push({
                                elem: atom.elem,
                                x: atom.x + tx,
                                y: atom.y + ty,
                                z: atom.z + tz,
                                color: atom.color,
                                constrained: atom.constrained
                            });
                        });
                    }
                }
            }
            return repeated;
        }

        // Helper to determine if background is light or dark
        function isLightBackground(color) {
            let r, g, b;
            if (color.startsWith('#')) {
                const hex = color.slice(1);
                r = parseInt(hex.substr(0, 2), 16);
                g = parseInt(hex.substr(2, 2), 16);
                b = parseInt(hex.substr(4, 2), 16);
            } else if (color.startsWith('rgb')) {
                const match = color.match(/\\d+/g);
                if (match) {
                    r = parseInt(match[0]);
                    g = parseInt(match[1]);
                    b = parseInt(match[2]);
                }
            } else {
                return false;
            }
            const luminance = (0.299 * r + 0.587 * g + 0.114 * b) / 255;
            return luminance > 0.5;
        }

        let unitCellColor = 'yellow';

        function initViewer() {
            const element = document.getElementById('viewer');
            // Get theme background color - this will adapt to light/dark themes
            const bgColor = getComputedStyle(document.body).getPropertyValue('--vscode-editor-background').trim() || 'black';
            unitCellColor = isLightBackground(bgColor) ? '#555555' : 'yellow';
            viewer = $3Dmol.createViewer(element, { backgroundColor: bgColor });
            addAtoms();
            viewer.zoomTo();
            viewer.render();
        }

        function addAtoms() {
            viewer.removeAllModels();
            viewer.removeAllLabels();
            viewer.removeAllShapes();

            const atoms = getRepeatedAtoms();

            atoms.forEach((atom, i) => {
                const radius = styleMode === 'sphere' ? 0.5 : 0.2;
                const opacity = atom.constrained ? 0.9 : 1.0;

                viewer.addSphere({
                    center: {x: atom.x, y: atom.y, z: atom.z},
                    radius: radius,
                    color: atom.color,
                    opacity: opacity
                });

                // Add X marker for constrained atoms (X in XY plane + vertical line)
                if (atom.constrained) {
                    const size = radius * 1.2;
                    viewer.addCylinder({
                        start: {x: atom.x - size, y: atom.y - size, z: atom.z},
                        end: {x: atom.x + size, y: atom.y + size, z: atom.z},
                        radius: 0.05,
                        color: 'red',
                        fromCap: true,
                        toCap: true
                    });
                    viewer.addCylinder({
                        start: {x: atom.x + size, y: atom.y - size, z: atom.z},
                        end: {x: atom.x - size, y: atom.y + size, z: atom.z},
                        radius: 0.05,
                        color: 'red',
                        fromCap: true,
                        toCap: true
                    });
                    // Add vertical line (Z direction) for visibility from all angles
                    viewer.addCylinder({
                        start: {x: atom.x, y: atom.y, z: atom.z - size},
                        end: {x: atom.x, y: atom.y, z: atom.z + size},
                        radius: 0.05,
                        color: 'red',
                        fromCap: true,
                        toCap: true
                    });
                }

                if (showLabels) {
                    viewer.addLabel(atom.elem, {
                        position: {x: atom.x, y: atom.y, z: atom.z},
                        fontSize: 12,
                        fontColor: 'white',
                        backgroundOpacity: 0.5
                    });
                }
            });

            addUnitCell();
        }

        function addUnitCell() {
            if (!showCell) return;

            const a = lattice[0];
            const b = lattice[1];
            const c = lattice[2];

            const superA = [a[0] * repeatA, a[1] * repeatA, a[2] * repeatA];
            const superB = [b[0] * repeatB, b[1] * repeatB, b[2] * repeatB];
            const superC = [c[0] * repeatC, c[1] * repeatC, c[2] * repeatC];

            const origin = [0, 0, 0];
            const corners = [
                origin,
                superA,
                superB,
                superC,
                [superA[0] + superB[0], superA[1] + superB[1], superA[2] + superB[2]],
                [superA[0] + superC[0], superA[1] + superC[1], superA[2] + superC[2]],
                [superB[0] + superC[0], superB[1] + superC[1], superB[2] + superC[2]],
                [superA[0] + superB[0] + superC[0], superA[1] + superB[1] + superC[1], superA[2] + superB[2] + superC[2]]
            ];

            const edges = [[0,1],[0,2],[0,3],[1,4],[1,5],[2,4],[2,6],[3,5],[3,6],[4,7],[5,7],[6,7]];

            edges.forEach(([i, j]) => {
                viewer.addCylinder({
                    start: {x: corners[i][0], y: corners[i][1], z: corners[i][2]},
                    end: {x: corners[j][0], y: corners[j][1], z: corners[j][2]},
                    radius: 0.03,
                    color: unitCellColor,
                    fromCap: true,
                    toCap: true
                });
            });
        }

        function updateView() {
            addAtoms();
            viewer.render();
        }

        function incrementRepeat(dir) {
            const input = document.getElementById('repeat' + dir);
            const val = parseInt(input.value) || 1;
            if (val < 10) {
                input.value = val + 1;
                updateRepeatValue(dir);
            }
        }

        function decrementRepeat(dir) {
            const input = document.getElementById('repeat' + dir);
            const val = parseInt(input.value) || 1;
            if (val > 1) {
                input.value = val - 1;
                updateRepeatValue(dir);
            }
        }

        function updateRepeatValue(dir) {
            const input = document.getElementById('repeat' + dir);
            let val = parseInt(input.value) || 1;
            val = Math.max(1, Math.min(10, val));
            input.value = val;

            if (dir === 'A') repeatA = val;
            else if (dir === 'B') repeatB = val;
            else if (dir === 'C') repeatC = val;

            updateView();
            viewer.zoomTo();
            viewer.render();
        }

        document.getElementById('resetView').addEventListener('click', () => {
            viewer.zoomTo();
            viewer.render();
        });

        document.getElementById('toggleCell').addEventListener('click', (e) => {
            showCell = !showCell;
            e.target.classList.toggle('active', showCell);
            updateView();
        });

        document.getElementById('toggleLabels').addEventListener('click', (e) => {
            showLabels = !showLabels;
            e.target.classList.toggle('active', showLabels);
            viewer.removeAllLabels();
            updateView();
        });

        document.getElementById('styleStick').addEventListener('click', (e) => {
            styleMode = 'stick';
            document.getElementById('styleSphere').classList.remove('active');
            e.target.classList.add('active');
            updateView();
        });

        document.getElementById('styleSphere').addEventListener('click', (e) => {
            styleMode = 'sphere';
            document.getElementById('styleStick').classList.remove('active');
            e.target.classList.add('active');
            updateView();
        });

        ['A', 'B', 'C'].forEach(dir => {
            const input = document.getElementById('repeat' + dir);
            input.addEventListener('change', () => updateRepeatValue(dir));
            input.addEventListener('keyup', (e) => {
                if (e.key === 'Enter') updateRepeatValue(dir);
            });
        });

        // Save PNG button
        document.getElementById('savePng').addEventListener('click', () => {
            const pngUri = viewer.pngURI();
            const link = document.createElement('a');
            link.href = pngUri;
            link.download = 'structure.png';
            link.click();
        });

        try {
            initViewer();
        } catch (e) {
            console.error("Error initializing viewer:", e);
            document.getElementById('viewer').innerHTML = '<p style="color:red;padding:20px;">Error: ' + e.message + '</p>';
        }
        ` : ''}

        ${relaxationScript}

        ${relaxationData ? 'initRelaxChart();' : ''}
    </script>
</body>
</html>`;
    }

    private _escapeHtml(text: string): string {
        return text
            .replace(/&/g, '&amp;')
            .replace(/</g, '&lt;')
            .replace(/>/g, '&gt;')
            .replace(/"/g, '&quot;')
            .replace(/'/g, '&#039;');
    }
}

export function activate(context: vscode.ExtensionContext) {
    console.log('VASPsum extension is now active');

    // Register hover provider for INCAR and OUTCAR files
    const hoverProvider = new VaspHoverProvider();
    context.subscriptions.push(
        vscode.languages.registerHoverProvider('incar', hoverProvider),
        vscode.languages.registerHoverProvider('outcar', hoverProvider)
    );

    // Register document symbol provider for OUTCAR navigation
    const outcarSymbolProvider = new OutcarSymbolProvider();
    context.subscriptions.push(
        vscode.languages.registerDocumentSymbolProvider('outcar', outcarSymbolProvider)
    );

    // Register completion provider for INCAR files
    const completionProvider = new VaspCompletionProvider();
    context.subscriptions.push(
        vscode.languages.registerCompletionItemProvider('incar', completionProvider)
    );

    // Register CodeLens provider for all VASP file types
    const codeLensProvider = new VaspCodeLensProvider();
    const vaspLanguages = ['incar', 'poscar', 'kpoints', 'outcar', 'chgcar', 'doscar', 'xdatcar', 'oszicar', 'eigenval', 'procar', 'ibzkpt', 'pcdat'];
    vaspLanguages.forEach(lang => {
        context.subscriptions.push(
            vscode.languages.registerCodeLensProvider(lang, codeLensProvider)
        );
    });

    // Register command to show system info
    context.subscriptions.push(
        vscode.commands.registerCommand('vasp.showSystemInfo', () => {
            SystemInfoPanel.createOrShow(context);
        })
    );

    // Register command to preview POSCAR/CONTCAR structure
    context.subscriptions.push(
        vscode.commands.registerCommand('vasp.previewStructure', () => {
            const editor = vscode.window.activeTextEditor;
            if (!editor) {
                vscode.window.showErrorMessage('No active editor');
                return;
            }

            const document = editor.document;
            const languageId = document.languageId;

            if (languageId !== 'poscar') {
                vscode.window.showErrorMessage('Structure preview is only available for POSCAR/CONTCAR files');
                return;
            }

            StructurePreviewPanel.createOrShow(document.uri);
        })
    );

    // Register a command to open VASP wiki for the word under cursor
    context.subscriptions.push(
        vscode.commands.registerCommand('vasp.openWiki', () => {
            const editor = vscode.window.activeTextEditor;
            if (!editor) {
                return;
            }

            const position = editor.selection.active;
            const range = editor.document.getWordRangeAtPosition(position, /[A-Za-z_][A-Za-z0-9_]*/);
            if (!range) {
                return;
            }

            const word = editor.document.getText(range).toUpperCase();
            const param = VASP_PARAMETERS[word];

            if (param) {
                const wikiUrl = `${VASP_WIKI_BASE}${param.wikiPage}`;
                vscode.env.openExternal(vscode.Uri.parse(wikiUrl));
            } else {
                vscode.window.showInformationMessage(`No VASP documentation found for '${word}'`);
            }
        })
    );

    // Register command to show VASP calculation summary
    context.subscriptions.push(
        vscode.commands.registerCommand('vasp.showSummary', () => {
            VaspSummaryPanel.createOrShow();
        })
    );

    // Register command to show relaxation progress
    context.subscriptions.push(
        vscode.commands.registerCommand('vasp.showRelaxation', () => {
            RelaxationViewerPanel.createOrShow();
        })
    );

    // Register unified VASPsum command
    context.subscriptions.push(
        vscode.commands.registerCommand('vasp.vaspsum', () => {
            VaspSumPanel.createOrShow(context);
        })
    );
}

export function deactivate() {}
