# 3D-Lithofacies-Manifold-TPM

# 3D Lithofacies Manifold TPM Simulation

This repository provides MATLAB scripts for performing three-dimensional (3D) stochastic lithofacies simulations integrating **manifold embedding** with **transition probability modeling (TPM)**. The method explicitly incorporates structural geological constraints (pole-to-plane orientations), enabling realistic non-stationary modeling of subsurface lithological structures without the need for predefined Training Images (TIs).

The provided MATLAB script is specifically configured for **unconditional simulation** (related to Figure 4 in the associated manuscript).

---

## Main MATLAB script

- **`unconditional_simulation.m`**: Performs an unconditional lithofacies stochastic simulation using synthetic manifold embedding and TPM.

---

## Explanation of simulation modes

The script can operate in two modes:

- **Unconditional simulation (`con_unc=0`)**
  - Does not require external data files.
  - Uses internally defined mean lengths (MLs) and target lithological proportions.

- **Conditional simulation (`con_unc=1`)**
  - Requires external data files (`cdat.txt` and `TPM.txt`).
  - Conditioning information such as known lithology points and predefined transition probabilities must be provided explicitly by users.

Currently, the script is set to `con_unc=0` (unconditional simulation) and runs as a standalone code.

---

## Key Input Parameters

| Parameter   | Explanation                                                     | Default Setting (Unconditional mode) |
|-------------|-----------------------------------------------------------------|--------------------------------------|
| `nx, ny, nz`| Number of grid cells in the x, y, and z directions              | 100, 100, 100                        |
| `nl`        | Number of lithologies                                           | 6                                    |
| `MLengths`  | Mean lengths of lithological units                              | `[50, 30, 100, 5, 20, 40]`           |
| `p_tgt`     | Target lithological proportions                                 | `[0.2115, 0.1304, 0.3599, 0.0330, 0.0835, 0.1817]` |
| `seed`      | Random seed (for reproducibility)                               | 1234                                 |
| `bet`       | Pole-to-plane gradient scaling coefficient                      | 5                                    |
| `nabsc`     | Number of numerical integration points (abscissae)              | 5                                    |

These parameters explicitly control the generated lithofacies domains' geological structure and complexity. Users can adjust them according to their specific requirements.

---

## How to Run

1. Ensure you have MATLAB installed (preferably MATLAB R2021a or newer).
2. Clone or download this repository.
3. Open MATLAB and navigate to the directory containing `unconditional_simulation.m`.
4. Run the script directly:

```matlab
>> unconditional_simulation
```

---

## Conditional Simulation

To run a **conditional simulation**, set:

```matlab
con_unc=1;
```

Then, prepare these files:
- **`cdat.txt`**: Conditioning lithology data. Each row: `[x_coordinate, y_coordinate, z_coordinate, lithology_code]`
- **`TPM.txt`**: Transition Probability Matrix (pre-defined TPM). This file must be a square matrix (nl × nl), explicitly representing transition probabilities between lithologies.

Example format for `cdat.txt`:
```
10 15 20 2
11 15 20 1
12 15 20 3
...
```

Example format for `TPM.txt`:
```
0.70 0.10 0.20
0.20 0.60 0.20
0.10 0.30 0.60
```

---

## Output

- **Visualization**: Generated lithofacies distributions are explicitly visualized in MATLAB figure windows.
- **Data arrays**: Lithological domain is available as `dom_3D`, a 3D numerical array (`ny × nx × nz`) representing the simulated lithologies.

---

## Important MATLAB Functions Explained

| Function Name  | Description                                            |
|----------------|--------------------------------------------------------|
| `geodist`      | Computes geodesic distances considering manifold curvature |
| `lgwt`         | Computes Legendre-Gauss abscissae and weights for numerical integration |
| `lith_est`     | Determines lithology stochastically based on TPM and distance metrics |
| `TPM_gen1`     | Generates TPM based on arbitrary counts and mean lengths |
| `TPM_gen2`     | Generates TPM constrained by target lithological proportions |

These functions explicitly implement critical algorithmic steps detailed in the associated manuscript.

---

## Citation

If you use this code, please cite our manuscript explicitly:

```
Park, E., [Full Author list]. (Year). Title of the paper. Journal Name, Volume(Issue), Page Numbers. DOI.
```

---

## License

This repository and code are provided under the MIT License.

---

## Contact

For questions, issues, or collaboration opportunities, please contact:

- **Eungyu Park**
- **park.eungyu@gmail.com**
- **Dept. Geology, Kyungpook National University**
