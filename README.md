# Membrane Bending Energy Calculation from Molecular Dynamics Trajectories

## Overview

This Python script calculates the bending energy of lipid membranes from GROMACS molecular dynamics (MD) trajectories.

The method reconstructs a two-dimensional membrane height surface from selected lipid beads and calculates membrane curvature using four different numerical approaches:

1. Small Gradient Approximation (SGA)
2. Extended Small Gradient Approximation (SGA2)
3. Full nonlinear curvature (Gamma)
4. Full nonlinear curvature with central-difference derivatives (Gamma2)

The script uses:

- **MDAnalysis** for reading GROMACS topology and trajectories
- **NumPy** for numerical calculations

The output provides frame-by-frame bending energies that can be used to analyze membrane deformation caused by proteins, ligands, phase separation, or other molecular interactions.

---

# Features

- Reads GROMACS `.gro` topology files
- Reads GROMACS `.xtc` trajectory files
- User-defined membrane bead selection
- User-defined grid resolution
- Periodic boundary condition treatment
- Calculates bending energy for every trajectory frame
- Calculates four curvature estimators
- Saves frame-wise bending energies

---

# Requirements

Python packages required:

```
numpy
MDAnalysis
```

Install dependencies:

```bash
pip install numpy MDAnalysis
```

Recommended environment:

```
Python >= 3.9
NumPy >= 1.22
MDAnalysis >= 2.4
```
---
# Repository Structure

The repository is organized as follows:

```text
membrane-bending-energy/
│
├── src/
│   └── bending_energy.py
│       # Main Python script for calculating membrane bending energy
│
├── examples/
│   ├── run_example.sh
│   │   # Example command for running the analysis
│   │
│   └── example_output.txt
│       # Example bending energy output
│
├── README.md
│   # Documentation, theory, installation, usage, and output description
│
├── requirements.txt
│   # Python package dependencies
│
├── LICENSE
│   # Software license information
│
└── .gitignore
    # Files excluded from Git version control
```

## File Description

| File | Description |
|------|-------------|
| `src/bending_energy.py` | Main Python script for calculating membrane bending energy from GROMACS molecular dynamics trajectories. |
| `examples/run_example.sh` | Example command showing how to run the analysis. |
| `examples/example_output.txt` | Example frame-by-frame bending energy output. |
| `README.md` | Complete documentation including theory, installation, usage, input/output description, and interpretation. |
| `requirements.txt` | Python dependencies required to run the script. |
| `LICENSE` | License information for using and distributing the software. |
| `.gitignore` | Prevents temporary files and large simulation files from being uploaded. |

## Data Management

Large molecular dynamics trajectory files are not included in this repository.

The following files should normally **not** be uploaded:

```text
*.xtc
*.trr
*.tpr
*.edr
```

Users should provide their own GROMACS topology and trajectory files.

---

# Files Required

The script requires:

## 1. GROMACS topology file

Example:

```
system.gro
```

This contains:

- atom names
- residue names
- coordinates
- box dimensions

---

## 2. GROMACS trajectory file

Example:

```
trajectory.xtc
```

The trajectory contains the membrane coordinates during simulation.

---

# Running the Script

The general command format is:

```bash
python3 bending_energy.py GRO_FILE XTC_FILE N_GRIDS "ATOM_SELECTION"
```

---

# Example

Example system:

```
raw_20.gro
raw_0.xtc
```

Run:

```bash
python3 bending_energy.py \
raw_20.gro \
raw_0.xtc \
8 \
"resname POPC CHOL DBSM"
```

---

# Input Parameters

## 1. GRO_FILE

Example:

```
raw_20.gro
```

The GROMACS topology file.

---

## 2. XTC_FILE

Example:

```
raw_0.xtc
```

The trajectory file.

The script checks whether the file exists and is not empty before loading.

---

## 3. N_GRIDS

Example:

```
8
```

Number of grid divisions along the x and y directions.

The membrane surface is divided into:

$$
N_{grid} \times N_{grid}
$$

cells.

Example:

```
N_GRIDS = 8
```

creates:

$$
8 \times 8 = 64
$$

surface elements.

Higher values:

Advantages:
- better spatial resolution

Disadvantages:
- more sensitive to noise


Lower values:

Advantages:
- smoother surface

Disadvantages:
- loss of local curvature information

---

## 4. ATOM_SELECTION

The membrane surface is reconstructed from selected atoms/beads using MDAnalysis selection syntax.

Example:

```bash
"resname POPC CHOL DBSM"
```

selects all beads belonging to:

- POPC
- CHOL
- DBSM

Other examples:

Only POPC:

```bash
"resname POPC"
```

Multiple lipid species:

```bash
"resname POPC CHOL POPS"
```

Specific bead:

```bash
"resname POPC and name PO4"
```

---

# How the Calculation Works

The workflow is:

```
GROMACS trajectory
        |
        |
        v
Select lipid beads
        |
        |
        v
Divide membrane into grid
        |
        |
        v
Calculate average z-height
        |
        |
        v
Generate membrane height surface h(x,y)
        |
        |
        v
Calculate curvature
        |
        |
        v
Calculate bending energy
        |
        |
        v
Save frame-wise energy
```

---

# Membrane Height Calculation

For every trajectory frame:

The selected lipid beads are placed into an XY grid.

For every grid cell:

$$
h(x,y)=\frac{1}{N}\sum z_i
$$

where:

- \(h(x,y)\) is membrane height
- \(z_i\) are bead coordinates

This produces a two-dimensional membrane surface.

---

# Curvature Calculations

The script calculates four bending-energy estimates.

---

# 1. SGA Energy

Output:

```
SGA
```

Small Gradient Approximation.

The membrane curvature is approximated as:

$$
H \approx \nabla^2 h
$$


The bending energy is:

$$
E =
\int
(\nabla^2h)^2 dxdy
$$


Advantages:

- Fast
- Stable
- Suitable for nearly flat membranes


Limitations:

- Less accurate for highly curved membranes

---

# 2. SGA2 Energy

Output:

```
SGA2
```

Uses the same approximation but calculates derivatives using a larger finite difference stencil.

Instead of:

$$
h(x+\Delta)-h(x-\Delta)
$$

it uses:

$$
h(x+2\Delta)-h(x-2\Delta)
$$


Advantages:

- Reduced numerical noise
- Smoother curvature calculation


Limitations:

- Lower spatial resolution

---

# 3. Gamma Energy

Output:

```
Gamma
```

Uses the full nonlinear membrane curvature:

$$
H=
\nabla \cdot
\left(
\frac{\nabla h}
{\sqrt{1+|\nabla h|^2}}
\right)
$$


The energy is:

$$
E=
\int
\sqrt{1+|\nabla h|^2}
H^2 dxdy
$$


Advantages:

- More physically accurate
- Includes membrane slope effects
- Suitable for highly deformed membranes

---

# 4. Gamma2 Energy

Output:

```
Gamma2
```

Same nonlinear curvature method as Gamma but uses central differences:

$$
\frac{f(x+\Delta)-f(x-\Delta)}
{2\Delta}
$$


Advantages:

- Higher numerical accuracy
- Reduced discretization error

Recommended for publication-quality analysis.

---

# Output

After completion the script generates:

```
energy_binding_results.txt
```

---

## Output File Format

Example:

```
Frame   SGA       SGA2      Gamma     Gamma2
0       12.345    11.876    13.521    13.210
1       12.512    12.021    13.782    13.450
2       12.734    12.213    14.012    13.821
```

---

# Output Columns

| Column | Description |
|---|---|
| Frame | MD trajectory frame number |
| SGA | Small-gradient bending energy |
| SGA2 | Extended small-gradient bending energy |
| Gamma | Nonlinear curvature bending energy |
| Gamma2 | Nonlinear curvature using central differences |

---

# Interpretation of Results

For a flat membrane:

```
SGA ≈ SGA2 ≈ Gamma ≈ Gamma2
```

For strongly curved membranes:

```
Gamma and Gamma2 > SGA and SGA2
```

because nonlinear curvature includes geometric corrections ignored by the small-gradient approximation.

---

# Example Analysis

Load results:

```python
import numpy as np

data = np.loadtxt(
    "energy_binding_results.txt",
    skiprows=1
)

frame = data[:,0]

SGA = data[:,1]

SGA2 = data[:,2]

Gamma = data[:,3]

Gamma2 = data[:,4]
```

Calculate average bending energy:

```python
print(np.mean(Gamma2))
```

---

# Notes

## Periodic Boundary Conditions

The script applies periodic indexing:

```python
(i+1)%N_grid
```

to handle membrane boundaries.

---

## Grid Selection

Recommended:

| Membrane size | Grid |
|-|-|
| Small membrane | 6-10 |
| Medium membrane | 10-20 |
| Large membrane | 20+ |

---

# Troubleshooting

## Error: No usable XTC files found

Check:

- trajectory path
- file size
- permissions


## Error: Empty output

Possible causes:

- incorrect atom selection
- wrong residue names
- missing membrane beads


Check selection:

```python
print(universe.residues.resnames)
```

---

# Citation and References

This script is an independent analysis tool developed for calculating membrane bending energy from molecular dynamics trajectories.

If you use this script in your research, please cite the relevant theoretical framework and software dependencies:

## Membrane curvature theory

The bending energy calculation is based on the Helfrich membrane model:

Helfrich, W.  
"Elastic properties of lipid bilayers: theory and possible experiments"  
*Zeitschrift für Naturforschung C*, 28(11–12), 693–703 (1973).

## Molecular dynamics trajectory analysis

Trajectory processing is performed using MDAnalysis:

Michaud-Agrawal, N., Denning, E. J., Woolf, T. B., & Beckstein, O.  
"MDAnalysis: a toolkit for the analysis of molecular dynamics simulations."  
*Journal of Computational Chemistry*, 32(10), 2319–2327 (2011).


# Author

Developed for membrane curvature and bending-energy analysis from molecular dynamics simulations.
## License

This project is released under the MIT License.

See the [LICENSE](LICENSE) file for details.

