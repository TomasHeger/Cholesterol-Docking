# ChimeraX Cholesterol Docking

**Cholesterol Docking V7** is a plugin for [UCSF ChimeraX](https://www.cgl.ucsf.edu/chimerax/) designed for automated identification, docking, and evaluation of cholesterol molecules in cryo-electron microscopy (cryo-EM) density maps.

The plugin identifies membrane-associated density features potentially corresponding to cholesterol molecules and performs rigid-body docking of cholesterol (CLR) into candidate sites. It evaluates alternative molecular orientations, including head-to-tail arrangements, and jointly optimizes neighboring cholesterol candidates to improve their placement.

Candidate docking sites are classified into confidence tiers, and detailed diagnostic information is generated to support manual inspection and structural interpretation.

The tool is intended to assist cryo-EM researchers in identifying potential cholesterol-binding sites, examining membrane-associated density, and building cholesterol-containing structural models.

## Features

- **Automated candidate detection:** Identifies membrane-associated density features potentially corresponding to cholesterol molecules.
- **Rigid-body cholesterol docking:** Fits cholesterol molecules into candidate density regions using multiple orientations.
- **Orientation optimization:** Evaluates alternative cholesterol orientations, including head-to-tail arrangements.
- **Multi-site optimization:** Jointly optimizes neighboring cholesterol candidates.
- **Confidence classification:** Groups accepted cholesterol placements into high-, medium-, and low-confidence categories.
- **Rejection diagnostics:** Records information about rejected candidates to facilitate evaluation of docking results.
- **Structured output:** Generates PDB models and CSV reports for subsequent visualization and analysis.

## Applications

Cholesterol Docking V7 is intended to support cryo-EM structure determination and model building, particularly for membrane proteins with surrounding density features that may correspond to cholesterol molecules.

Potential applications include:

- Identifying candidate cholesterol-binding sites in membrane protein structures.
- Investigating unassigned membrane-associated cryo-EM density.
- Generating initial cholesterol placements for subsequent refinement.
- Comparing alternative cholesterol orientations within candidate binding sites.
- Evaluating neighboring cholesterol placements.
- Prioritizing candidate sites for manual inspection and structural interpretation.

Automated cholesterol docking provides candidate molecular assignments rather than definitive identification of cholesterol. Interpretation should consider map quality, local resolution, molecular environment, and independent experimental evidence where available.

## How It Works

The plugin processes a cryo-EM density map and an associated atomic model through the following steps:

1. **Candidate detection:** Identifies membrane-associated density regions that may accommodate cholesterol molecules.
2. **Cholesterol template preparation:** Obtains a cholesterol (CLR) molecular template from the Chemical Component Dictionary.
3. **Rigid-body docking:** Evaluates multiple orientations of cholesterol within candidate density regions.
4. **Orientation evaluation:** Considers alternative head-to-tail orientations to improve candidate placement.
5. **Neighbor optimization:** Jointly optimizes neighboring cholesterol candidates.
6. **Candidate evaluation:** Classifies accepted placements into confidence tiers and records rejection diagnostics.
7. **Output generation:** Exports docked cholesterol models and structured reports for further analysis.

## Installation

The plugin is provided as a ChimeraX bundle and can be built and installed from its source directory.

The supplied version is intended for **ChimeraX 1.12**.

### Installation steps

1. Download or clone the repository to a local directory.
2. Open UCSF ChimeraX.
3. Build the bundle using the ChimeraX command line:

   ```
   devel build /path/to/CholesterolDocking
   ```

4. Install the bundle:

   ```
   devel install /path/to/CholesterolDocking noDeps true user true
   ```

5. Restart ChimeraX.

Replace `/path/to/CholesterolDocking` with the actual location of the downloaded repository.

### Graphical user interface

After installation, the tool can be accessed through:

**Tools → Structure Analysis → Cholesterol Docking V7**

## Command-Line Usage

The plugin also provides the ChimeraX command `choldock7`.

Example:

```
choldock7 #1 proteinModels #2 detectionLevel 0.141 supportLevel 0.05 orientations 1500 maxsites 60 out /path/to/output
```

In this example:

- `#1` specifies the cryo-EM density map model.
- `#2` specifies the protein atomic model.
- `detectionLevel` specifies the density level used for candidate detection.
- `supportLevel` specifies the density level used for supporting density analysis.
- `orientations` specifies the number of orientations evaluated during docking.
- `maxsites` limits the number of candidate sites processed.
- `out` specifies the output directory.

The example values are illustrative and may require adjustment depending on the density map and structural model.

## Output Files

The plugin generates the following files:

| File | Description |
|---|---|
| `all_accepted_cholesterols.pdb` | All accepted cholesterol placements. |
| `high_confidence_cholesterols.pdb` | Cholesterol placements assigned to the high-confidence tier. |
| `medium_confidence_cholesterols.pdb` | Cholesterol placements assigned to the medium-confidence tier. |
| `low_confidence_cholesterols.pdb` | Cholesterol placements assigned to the low-confidence tier. |
| `final_selection.csv` | Summary of final candidate selection results. |
| `site_diagnostics.csv` | Site-level diagnostic information. |
| `neighbor_groups.csv` | Information about neighboring candidate groups. |
| `membrane_model.csv` | Information associated with the membrane model used during analysis. |

These outputs facilitate further examination of candidate cholesterol placements and their incorporation into subsequent structural analysis workflows.

## Requirements

- UCSF ChimeraX 1.12
- ChimeraX Core and UI components
- SciPy, including `scipy.ndimage`, for connected-component segmentation
- Internet access from ChimeraX to retrieve the cholesterol component `CCD:CLR`

## Limitations

The identification and docking of cholesterol molecules depend on the quality of the input cryo-EM map, the selected density thresholds, and the surrounding structural environment.

Density features attributed to cholesterol may also originate from other lipids, detergents, or molecular components. Consequently, a favorable docking result does not independently establish the chemical identity of a density feature.

The high-, medium-, and low-confidence classifications represent evidence tiers assigned by the docking procedure. They should not be interpreted as calibrated probabilities of correct molecular identification.

All automated cholesterol assignments should be inspected manually before biological interpretation, further refinement, or inclusion in a final structural model.

## License

This project is distributed under the MIT License. See the `LICENSE` file for details.
