
# MDPerTool Documentation

Welcome to the **MDPerTool** documentation! This guide is designed to help users of all backgrounds understand, install, and use MDPerTool for allosteric network analysis, with a focus on clarity, scientific context, and practical workflows.

---

## Table of Contents

- [Introduction](#introduction)
- [Quickstart](#quickstart)
- [Theory & Background](#theory--background)
- [Workflow Overview](#workflow-overview)
- [Residue Selection Guide](#residue-selection-guide)
- [Step-by-Step Tutorials](#step-by-step-tutorials)
- [FAQ & Troubleshooting](#faq--troubleshooting)
- [References](#references)

---

## Introduction

**MDPerTool** is a toolkit for analyzing allosteric communication and signal propagation in protein structures using network-based metrics such as the Propagation Coefficient (PC). It is inspired by key publications in the field and provides both command-line and GUI interfaces for flexible analysis.

---

## Quickstart

1. **Installation**
	```bash
	pip install mdpertool
	# or clone and install locally
	git clone https://github.com/Bio-Otto/MDPerTool_GUI.git
	cd MDPerTool_GUI
	pip install .
	```
2. **Basic Usage**
	```bash
	mdpertool --input myprotein.pdb --source ALA123 --output results/
	```
3. **GUI Launch**
	```bash
	python -m mdpertool
	```

---

## Theory & Background

### What is Allosteric Network Analysis?
Allosteric network analysis models a protein as a graph, where nodes are residues and edges represent communication or energy transfer. The goal is to identify key residues (hubs, relays, sources) that mediate signal propagation.

### Propagation Coefficient (PC)
The PC metric quantifies a residue's ability to propagate a perturbation through the network. It is calculated as:

$$
PC(i) = \frac{m_i \cdot n_i}{\sum_j m_j n_j + l}
$$

Where:
- $m_i$ = in-degree (number of incoming edges)
- $n_i$ = out-degree (number of outgoing edges)
- $l$ = number of edges bypassing the current layer

### Key References
- [PMC3195717](https://pmc.ncbi.nlm.nih.gov/articles/PMC3195717/)
- [PMC3282753](https://pmc.ncbi.nlm.nih.gov/articles/PMC3282753/)
- [Tandfonline 2013](https://www.tandfonline.com/doi/full/10.1080/07391102.2013.855145?needAccess=true#d1e164)

---

## Workflow Overview

```mermaid
graph TD;
	 A[Input Structure (PDB)] --> B[Network Construction];
	 B --> C[Residue Selection];
	 C --> D[PC Calculation];
	 D --> E[Visualization & Output];
```

---

## Residue Selection Guide

!!! tip "How to Select a Perturbation Residue?"
	 - **Active/Allosteric Sites:** Prefer residues in known functional or allosteric regions.
	 - **Conserved Residues:** Highly conserved residues are often functionally important.
	 - **Experimental Data:** Use literature or mutagenesis data if available.
	 - **Exploratory:** You can systematically try different residues to map the allosteric landscape.
	 - **No Prior Knowledge?** Start with the active site or use the GUI to visually select residues.

---

## Step-by-Step Tutorials

### Example: Allosteric Network Analysis

1. **Prepare Input**
	- Obtain a PDB file of your protein.
2. **Select Perturbation Residue**
	- Use the GUI or specify via CLI (e.g., `--source ALA123`).
3. **Run Analysis**
	- View results in the GUI or as CSV/PNG outputs.
4. **Interpret Results**
	- Identify super-hubs, relays, and key allosteric pathways.

---

## FAQ & Troubleshooting

**Q: Which residue should I select as the perturbation source?**
A: See the [Residue Selection Guide](#residue-selection-guide).

**Q: My network is disconnected or PC is zero everywhere!**
A: Check your input structure, residue numbering, and ensure the selected residue is present in the network.

**Q: How do I interpret the PC plot?**
A: High PC values indicate residues that efficiently propagate signals; these are often allosteric hubs.

---

## References
- [PMC3195717](https://pmc.ncbi.nlm.nih.gov/articles/PMC3195717/)
- [PMC3282753](https://pmc.ncbi.nlm.nih.gov/articles/PMC3282753/)
- [Tandfonline 2013](https://www.tandfonline.com/doi/full/10.1080/07391102.2013.855145?needAccess=true#d1e164)

---

For more, see the full documentation and tutorials in the `docs/` folder.