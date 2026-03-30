# FAQ & Troubleshooting

## Frequently Asked Questions

**Q: Which residue should I select as the perturbation source?**
A: See the [Residue Selection Guide](index.md#residue-selection-guide). Prefer active/allosteric sites, conserved residues, or use the GUI for exploratory analysis.

**Q: My network is disconnected or PC is zero everywhere!**
A: Check your input structure, residue numbering, and ensure the selected residue is present in the network.

**Q: How do I interpret the PC plot?**
A: High PC values indicate residues that efficiently propagate signals; these are often allosteric hubs.

**Q: Can I analyze multiple residues systematically?**
A: Yes, you can run MDPerTool for each residue of interest and compare the results.

**Q: What file formats are supported?**
A: Standard PDB for structure, XTC/DCD for trajectories (for energy decomposition).

---

## Troubleshooting

- **Problem:** PC values are all zero.
  - **Solution:** Check that the selected residue exists in the network and is not isolated.
- **Problem:** Network is disconnected.
  - **Solution:** Ensure the input structure is complete and residue numbering matches.
- **Problem:** GUI does not display plots.
  - **Solution:** Check that all dependencies are installed and outputs are generated.

---

For more help, see the [Quickstart](index.md#quickstart) and [Theory](theory.md) sections.
