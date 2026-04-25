"""Residue Manager - Handle residue selection and conservation data."""

from __future__ import annotations

import csv
import os
from typing import Any, Dict, List, Sequence

from Bio.PDB import PDBParser


def _extract_residue_labels(pdb_path: str) -> List[str]:
    """Parse a PDB file and return labels formatted as `{ResName}{ResNum}{ChainID}`.

    Heteroatoms (waters, ligands) are skipped so only protein residues are returned,
    matching the original prody `select('protein')` behaviour the GUI depends on.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("s", pdb_path)

    labels: List[str] = []
    for model in structure:
        for chain in model:
            chain_id = chain.id
            for residue in chain:
                if residue.id[0] != " ":
                    continue
                labels.append(f"{residue.resname}{residue.id[1]}{chain_id}")
    return labels


class ResidueManager:
    """Manage residue selection, conservation scores, and related operations."""

    def fill_residue_combobox(self, pdb_path: str) -> List[str]:
        """Extract residue labels from PDB file and populate combobox."""
        return _extract_residue_labels(pdb_path)

    def _get_target_residues(
        self,
        all_residue_as_target: bool,
        target_res_comboBox: Any,
        selected_target_residues_listWidget: Any,
    ) -> List[str]:
        """Get list of target residues for analysis.
        
        Args:
            all_residue_as_target: Boolean flag to use all residues
            target_res_comboBox: Qt combobox with all available residues
            selected_target_residues_listWidget: Qt listwidget with selected residues
        
        Returns:
            List of target residue names
        """
        if all_residue_as_target:
            return [target_res_comboBox.itemText(i).strip()
                   for i in range(target_res_comboBox.count())]
        else:
            return [selected_target_residues_listWidget.item(x).text().strip()
                   for x in range(selected_target_residues_listWidget.count())]

    def _get_conservation_settings(
        self,
        use_conservation_checkBox: Any,
        conservation_PDB_ID_lineEdit: Any,
        conservation_pdb_chain_id_lineedit: Any,
        conserv_score_doubleSpinBox: Any,
    ) -> Dict[str, Any]:
        """Extract conservation calculation settings from UI controls.
        
        Returns:
            Dictionary with conservation configuration
        """
        return {
            'use_conservation': use_conservation_checkBox.isChecked(),
            'pdb_id': conservation_PDB_ID_lineEdit.text(),
            'chain': conservation_pdb_chain_id_lineedit.text(),
            'conservation_threshold': conserv_score_doubleSpinBox.value(),
            'save_conservation_scores': False
        }

    def _save_conservation_scores(
        self,
        res_IDs: Sequence[Any],
        con_scores: Sequence[Any],
        pdb_id: str,
        output_folder_directory: str,
    ) -> None:
        """Save conservation scores to CSV file.
        
        Args:
            res_IDs: List of residue identifiers
            con_scores: List of conservation scores
            pdb_id: PDB ID for filename
            output_folder_directory: Directory to save file in
        """
        rows = zip(res_IDs, con_scores)
        with open(os.path.join(output_folder_directory, f'conservation_{pdb_id}.csv'), "w", newline='') as f:
            writer = csv.writer(f)
            for row in rows:
                writer.writerow(row)

    def create_shortest_path_string(self, sp: Sequence[Any]) -> str:
        """Format shortest path data for display.
        
        Args:
            sp: Shortest path object/list
        
        Returns:
            Formatted string representation of shortest path
        """
        # Placeholder implementation - extend based on actual data structure
        shortest_str_form = ''
        return shortest_str_form
