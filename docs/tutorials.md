# Step-by-Step Tutorials

## Example 1: Basic Allosteric Network Analysis

### 1. Prepare Input
- Obtain a PDB file of your protein (e.g., `myprotein.pdb`).

### 2. Select Perturbation Residue
- Use the GUI to select a residue visually, or specify via CLI:
  ```bash
  mdpertool --input myprotein.pdb --source ALA123 --output results/
  ```

### 3. Run Analysis
- Launch the GUI:
  ```bash
  python -m mdpertool
  ```
- Or run from the command line as above.

### 4. View Results
- Inspect CSV and PNG outputs in the results folder.
- In the GUI, view PC plots and network diagrams interactively.

### 5. Interpret Results
- High PC values: Allosteric hubs or amplifiers.
- Low PC values: Relays, collectors, or isolated nodes.

---

## Example 2: Systematic Residue Scanning

1. For each residue of interest, repeat the analysis to map the allosteric landscape.
2. Compare PC plots to identify key communication pathways.

---

## Example 3: Troubleshooting

- **Problem:** PC values are all zero.
  - **Solution:** Check that the selected residue exists in the network and is not isolated.
- **Problem:** Network is disconnected.
  - **Solution:** Ensure the input structure is complete and residue numbering matches.

---

For more advanced tutorials, see the [docs/tutorials/](tutorials/) folder.
