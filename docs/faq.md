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

---

## Advanced Scientific FAQ

**Q: What is the minimum input required for an analysis?**  
A: A high-quality PDB file (all residues present, no missing atoms in key regions) and at least one perturbation (allosteric) site.

**Q: How do I choose the best perturbation site?**  
A: Start with known allosteric or regulatory residues ([PMC3282753](https://pmc.ncbi.nlm.nih.gov/articles/PMC3282753/)). Use sequence conservation and literature data to guide selection. Systematic scanning can reveal unexpected hotspots.

**Q: Why do some residues have very low or zero propagation coefficient (PC)?**  
A: These residues are peripheral to the main signal transduction pathways, may be structurally or dynamically isolated, or not involved in allosteric communication.

**Q: What if my network has no clear super-hubs?**  
A: Not all proteins have highly connected super-hubs; this is system-dependent. Check your input structure and simulation parameters for errors.

**Q: How can I validate my computational predictions?**  
A: Compare predicted key residues with experimental mutagenesis data ([PMC3282753, Table 3](https://pmc.ncbi.nlm.nih.gov/articles/PMC3282753/table/pone-0031529-t003/)). High overlap indicates strong predictive power.

**Q: Can I analyze multimeric or heteromeric proteins?**  
A: Yes. MDPerTool supports analysis of both intra- and inter-subunit communication ([Tandfonline 2013](https://www.tandfonline.com/doi/full/10.1080/07391102.2013.855145?needAccess=true#d1e164)). Analyze each subunit and the full complex for a complete picture.

**Q: What are the limitations of the energy dissipation model?**  
A: Relies on the accuracy of the input structure and simulation parameters. May not capture allosteric effects mediated by long-range water or ligand interactions not present in the structure. The model assumes energy input is localized and does not account for all possible biological perturbations.

**Q: How do I interpret dynamic modules?**  
A: Modules group residues by their response time to perturbation. Early modules often contain key regulatory residues; later modules may be more structural.

**Q: What if my results differ from published data?**  
A: Check for differences in input structure, simulation settings, or perturbation sites. Biological variability and experimental conditions can also cause differences.

**Q: Where can I find more scientific background?**  
A: See the [theory_energy_dissipation.md](theory_energy_dissipation.md), [advanced_scientific_insights.md](advanced_scientific_insights.md), and the referenced literature for in-depth explanations.

---

If you have additional questions, please consult the original articles or contact the MDPerTool development team.
