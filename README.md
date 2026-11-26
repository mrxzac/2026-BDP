# CIV102 Bridge Design Project (BDP)

A Python-based structural analysis framework for the CIV102 Matboard Bridge Design Project, University of Toronto, Engineering Science, Fall 2026.

This project provides a computational approach to designing, simulating, visualizing, and evaluating matboard box-girder bridges under moving train loads using elastic beam theory and thin-plate buckling theory.

## Features

- Object-oriented bridge geometry modelling.
- Arbitrary polygonal matboard panels, including trapezoidal webs.
- Variable-depth and variable-cross-section bridge designs.
- Unit-aware numerical calculations using Pint.
- Moving six-axle train load simulation in both orientations.
- Shear Force Diagrams (SFD) and Bending Moment Diagrams (BMD).
- Shear Force Envelopes (SFE) and Bending Moment Envelopes (BME).
- Cross-sectional centroid, first moment of area, and second moment of area calculations.
- Material failure and plate buckling analysis.
- Position-dependent glue shear analysis.
- Factor of Safety (FOS) and theoretical failure load prediction.
- Interactive 3D bridge visualization and cross-section inspection.
- CSV export of structural analysis results.

## Installation

Python 3 is required.

Install the dependencies:

```bash
pip install numpy matplotlib pint
```

## Quick Start

Run both standard loading cases:

```bash
python run.py
```

Run the analysis without graphical output:

```bash
python run.py --case both --no-plot
```

### Loading Cases

| Load case | Description | Total train load |
|---|---|---:|
| 1 | Uniformly loaded train | 400 N |
| 2 | Base loading with a heavier locomotive | 452 N |
| later | Subsequent loading protocol | User-defined |

```bash
python run.py --case 1
python run.py --case 2
python run.py --case later --weight 1000
```

### Command-Line Options

| Argument | Description |
|---|---|
| `--case` | Select `1`, `2`, `both`, or `later` |
| `--weight` | Total train load for subsequent loading cases |
| `--design` | Select configuration module and beam object |
| `--stations` | Number of longitudinal analysis stations |
| `--csv-dir` | Export numerical results to a directory |
| `--no-plot` | Disable graphical output |

Example:

```bash
python run.py --case both --stations 2501 --csv-dir results
```

## Project Structure

```text
BDP/
├── Beam.py
├── Matboard.py
├── Simulation.py
├── compute.py
├── strength.py
├── config_0.py
├── bridge_diagram.py
├── run.py
├── test_bdp.py
└── README.md
```

| File | Description |
|---|---|
| `Matboard.py` | Polygonal matboard geometry and cross-section extraction |
| `Beam.py` | Bridge assembly and sectional properties |
| `Simulation.py` | Moving-load simulation and force envelopes |
| `compute.py` | Numerical geometry and unit-handling utilities |
| `strength.py` | Material failure, buckling, and FOS calculations |
| `config_0.py` | Default bridge configuration |
| `bridge_diagram.py` | Interactive 3D geometry and sectional visualization |
| `run.py` | Main structural analysis entry point |
| `test_bdp.py` | Numerical regression tests against Design 0 |

## Coordinate System

The bridge is defined using a Cartesian coordinate system:

- **X-axis:** longitudinal direction of the bridge.
- **Y-axis:** vertical direction, positive upwards.
- **Z-axis:** transverse direction of the bridge.

All geometric dimensions are expressed in millimetres by default. Forces are expressed in newtons and stresses in megapascals.

## Creating a Bridge Configuration

A matboard is represented by an ordered set of vertices defining its mid-plane polygon and a nominal thickness direction.

```python
from Beam import Beam
from Matboard import Matboard, u

web = Matboard([
    (0, 0, 30),
    (600, 0, 30),
    (600, 120, 30),
    (0, 80, 30),
], "z")
```

The above example defines a trapezoidal web panel.

Multiple matboard panels can be assembled into a beam:

```python
beam_0 = Beam(
    longitudinal_boards + diaphragms,
    glue_joints
)
```

Glue joints are represented by their interface elevation and total bonded width. Position-dependent interfaces can be defined using callable functions.

Additional configurations can be stored in separate Python modules and selected without modifying the main analysis program:

```bash
python run.py --design config_1:beam_0
```

For compatibility with the strength solver, the longitudinal cross-section must contain two distinct structural webs, and each matboard must produce a connected cross-section at the requested cutting plane.

## Material Properties

Based on the CIV102 Bridge Design Project handout.

| Property | Value |
|---|---:|
| Matboard dimensions | 813 × 1016 mm |
| Matboard thickness | 1.27 mm |
| Tensile strength | 30 MPa |
| Compressive strength | 6 MPa |
| Shear strength | 4 MPa |
| Young's modulus | 4000 MPa |
| Poisson's ratio | 0.20 |
| Contact cement shear strength | Up to 2 MPa |

## Structural Analysis

The analysis follows the following computational workflow:

1. Define the bridge geometry and material properties.
2. Simulate the train moving across the bridge.
3. Calculate support reactions, shear forces, and bending moments.
4. Construct signed shear force and bending moment envelopes.
5. Extract cross-sectional geometry at each analysis station.
6. Calculate sectional properties and induced stresses.
7. Evaluate material failure and plate buckling capacities.
8. Determine the minimum Factor of Safety and governing failure mechanism.
9. Predict the theoretical maximum train load for the specified axle-load ratios.

### Failure Mechanisms

The strength solver evaluates eight failure mechanisms:

- Matboard tensile failure.
- Matboard compressive failure.
- Matboard shear failure.
- Glue shear failure.
- Middle flange buckling.
- Side flange buckling.
- Web compression buckling.
- Web shear buckling.

## Visualization

Launch the interactive 3D bridge viewer:

```bash
python bridge_diagram.py
```

Inspect specific cross-sections:

```bash
python bridge_diagram.py --x 300 --reference 625
```

Export a bridge visualization:

```bash
python bridge_diagram.py --x 300 --save bridge.png
```

The main analysis program also generates plots of force envelopes, structural capacities, and centroidal second moments of area.

## Testing

Run the regression tests:

```bash
python -m unittest test_bdp
```

The original regression tests are based on the prescribed Design 0 geometry and its analytical solutions. Retain the original Design 0 configuration when running these tests, or update the expected values appropriately for other designs.

## Assumptions and Limitations

- Linear elastic material behaviour.
- Simply supported beam with end overhangs.
- Small-deformation beam theory.
- Concentrated axle loads.
- Idealized plate buckling boundary conditions.
- Perfect composite action between bonded layers.
- Nominal axis-directed matboard thickness.
- Two-web box-girder geometry for strength calculations.
- No explicit modelling of progressive damage, adhesive peeling, splice failure, or local wheel contact failure.

Transverse diaphragms determine the local web buckling length but are excluded from longitudinal strength calculations.

Inclined flange thickness and variable-depth geometry are approximated within the beam-theory model. The software does not replace experimental testing or detailed finite element analysis.

Predicted failure loads represent theoretical estimates rather than guaranteed physical capacities.

## Future Development

- Automated cross-sectional geometry optimization.
- Material usage and two-dimensional cutting optimization.
- Local web and flange reinforcement modelling.
- Improved splice connection analysis.
- Generalized multi-web cross-sections.
- Deflection and stiffness calculations.
- Finite element analysis integration.
- Experimental validation.

## References

University of Toronto, *CIV102 Matboard Bridge Design Project*, Fall 2026, Revision 3.

CIV102 course lecture notes on beam theory, shear stress, and thin-plate buckling.

---

Developed for the University of Toronto Engineering Science CIV102 Bridge Design Project.
