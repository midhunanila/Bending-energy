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
## Output

The `calculate_bending_energy()` function returns a dictionary containing the bending energy calculated for every trajectory frame.

```python
results = calculate_bending_energy(...)

results["energy_sga"]
results["energy_sga2"]
results["energy_gamma"]
results["energy_gamma2"]
results["time"]
```

Each energy array contains one value per frame in the trajectory.

Example:

```python
{
    "energy_sga":   [115.2, 116.5, 118.1, ...],
    "energy_sga2":  [113.8, 114.9, 117.0, ...],
    "energy_gamma": [121.7, 123.4, 122.8, ...],
    "energy_gamma2":[120.5, 122.1, 121.9, ...],
    "time":         [0, 1, 2, 3, ...]
}
```

where

| Key | Description |
|------|-------------|
| `energy_sga` | Membrane bending energy calculated using the Small Gradient Approximation (SGA). |
| `energy_sga2` | Bending energy using the Small Gradient Approximation with a wider finite-difference stencil. |
| `energy_gamma` | Bending energy computed from the full nonlinear curvature expression. |
| `energy_gamma2` | Bending energy computed from the nonlinear curvature using central-difference derivatives. |
| `time` | Trajectory frame index corresponding to each calculated energy value. |

---

# Description of the Four Curvature Methods

The package computes four numerical estimates of membrane bending energy. They are **different numerical approximations of the same physical quantity**, allowing users to compare numerical accuracy and assess the validity of the small-gradient approximation.

## 1. Small Gradient Approximation (`energy_sga`)

This method assumes that the membrane is nearly flat,

\[
|\nabla h| \ll 1,
\]

where \(h(x,y)\) is the membrane height.

Under this approximation, the mean curvature is approximated by the Laplacian of the height field,

\[
H \approx \nabla^2 h.
\]

The bending energy is then calculated as

\[
E \propto
\int (\nabla^2 h)^2\,dx\,dy.
\]

### Advantages

- Computationally efficient
- Numerically stable
- Appropriate for nearly planar membranes

### Limitations

- Ignores nonlinear geometric effects
- Underestimates curvature for strongly bent membranes

---

## 2. Small Gradient Approximation with Extended Stencil (`energy_sga2`)

This method is identical to `energy_sga` except that second derivatives are evaluated using a larger finite-difference stencil (two grid spacings instead of one).

Compared with `energy_sga`, this produces a smoother curvature field and reduces numerical noise.

### Advantages

- Reduced sensitivity to local fluctuations
- Better suited for noisy trajectories
- More stable on coarse grids

### Limitations

- Slight loss of spatial resolution
- May smooth out highly localized curvature

---

## 3. Nonlinear Curvature (`energy_gamma`)

This method evaluates the complete mean curvature without assuming small membrane slopes.

The curvature is calculated as

\[
H=
\nabla\cdot
\left(
\frac{\nabla h}
{\sqrt{1+|\nabla h|^2}}
\right).
\]

The bending energy becomes

\[
E=
\int
\sqrt{1+|\nabla h|^2}
H^2
\,dx\,dy,
\]

where the additional factor

\[
\sqrt{1+|\nabla h|^2}
\]

accounts for the true membrane surface area.

### Advantages

- Physically accurate for strongly curved membranes
- Includes nonlinear geometric corrections
- Closely follows the Helfrich bending-energy formalism

### Limitations

- Computationally more expensive
- More sensitive to numerical noise

---

## 4. Nonlinear Curvature with Central Differences (`energy_gamma2`)

This method computes the same nonlinear curvature as `energy_gamma` but evaluates spatial derivatives using central differences,

\[
\frac{f(x+\Delta)-f(x-\Delta)}{2\Delta},
\]

instead of forward differences.

Central differences generally provide second-order numerical accuracy and reduce directional bias.

### Advantages

- Improved derivative accuracy
- Reduced discretization error
- Recommended for high-quality analyses

### Limitations

- Requires a sufficiently smooth height field
- Can be slightly more sensitive to noisy data

---

# Interpretation of the Results

The four energy estimates should agree closely for nearly flat membranes.

As membrane curvature increases, the nonlinear methods (`energy_gamma` and `energy_gamma2`) generally predict larger bending energies because they account for nonlinear geometric effects neglected by the small-gradient approximation.

For example,

| Method | Bending Energy |
|---------|---------------:|
| `energy_sga` | 132.4 |
| `energy_sga2` | 129.8 |
| `energy_gamma` | 141.5 |
| `energy_gamma2` | 139.7 |

In this example, the nonlinear methods predict higher bending energies due to stronger membrane deformations, while the small-gradient approximations underestimate the curvature contribution.

For highly curved membranes (e.g., protein-induced deformation, budding, or tubulation), `energy_gamma2` is generally recommended as the primary estimator because it combines the full nonlinear curvature expression with a more accurate central-difference approximation.
## Citation

If you use this software in published work, please cite

Author et al.

Journal

Year

DOI
