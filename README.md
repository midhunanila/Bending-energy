# Membrane Bending Energy

Python implementation of a grid-based membrane bending energy calculation from molecular dynamics trajectories.

## Features

- MDAnalysis based
- Supports multiple XTC trajectories
- Adjustable grid resolution
- Custom lipid selections
- Periodic boundary conditions
- Computes four bending-energy estimators

---

## Installation

```bash
git clone https://github.com/username/membrane_bending_energy.git

cd membrane_bending_energy

pip install -r requirements.txt
```

---

## Requirements

- Python ≥3.10
- numpy
- scipy
- MDAnalysis

Install

```bash
pip install -r requirements.txt
```

---

## Usage

```python
from bending_energy import calculate_bending_energy

results = calculate_bending_energy(

    gro_file="system.gro",

    xtc_files=["traj.xtc"],

    lipid_selection="resname POPC CHOL",

    n_grids=10,

)
```

---

## Input

| Parameter | Type | Description |
|-----------|------|-------------|
| gro_file | str | GROMACS topology |
| xtc_files | list | One or more XTC trajectories |
| lipid_selection | str | MDAnalysis atom selection |
| n_grids | int | Number of XY grid divisions |
| use_periodic | bool | Apply periodic boundary conditions |

---

## Output

The function returns a dictionary.

```python
{
 'energy_sga': [...],
 'energy_sga2': [...],
 'energy_gamma': [...],
 'energy_gamma2': [...],
 'time': [...]
}
```

---

## Method

1. Select membrane beads.
2. Divide the membrane plane into an N×N grid.
3. Compute the mean height (z) in every grid.
4. Calculate membrane gradients.
5. Compute local curvature.
6. Integrate curvature over the membrane surface to estimate bending energy.

---

## Citation

If you use this software in published work, please cite

Author et al.

Journal

Year

DOI
