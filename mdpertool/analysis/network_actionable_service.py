"""Build user-focused actionable residue insights from network outputs."""

from __future__ import annotations

from typing import Dict, List, Sequence, Tuple


ActionableRow = Tuple[str, float, int, str, str]


def _normalize_map(values: Dict[str, float]) -> Dict[str, float]:
    if not values:
        return {}
    max_value = max(values.values())
    if max_value <= 0:
        return {key: 0.0 for key in values}
    return {key: float(val) / float(max_value) for key, val in values.items()}


def build_actionable_residue_rows(
    critical_rows: Sequence[Tuple[str, int, float, float]],
    superhub_rows: Sequence[Tuple[str, int, float, float, str]],
    top_n: int = 12,
) -> List[ActionableRow]:
    """Rank residues that are most actionable for follow-up mutation or validation."""
    critical_composite: Dict[str, float] = {str(res): float(comp) for res, _hit, _bw, comp in critical_rows}
    critical_hits: Dict[str, int] = {str(res): int(hit) for res, hit, _bw, _comp in critical_rows}

    superhub_betweenness: Dict[str, float] = {str(node): float(bw) for node, _deg, bw, _cl, _nbr in superhub_rows}

    normalized_critical = _normalize_map(critical_composite)
    normalized_hub = _normalize_map(superhub_betweenness)

    candidate_residues = set(normalized_critical.keys()) | set(normalized_hub.keys())
    ranked: List[ActionableRow] = []

    for residue in sorted(candidate_residues):
        critical_score = normalized_critical.get(residue, 0.0)
        hub_score = normalized_hub.get(residue, 0.0)
        final_score = (0.7 * critical_score) + (0.3 * hub_score)
        path_hits = critical_hits.get(residue, 0)

        if critical_score > 0.0 and hub_score > 0.0:
            rationale = "High path impact + central hub"
            priority = "High"
        elif critical_score > 0.0:
            rationale = "High shortest-path involvement"
            priority = "Medium"
        else:
            rationale = "Hub-dominant signal relay"
            priority = "Medium"

        ranked.append((residue, final_score, path_hits, priority, rationale))

    ranked.sort(key=lambda row: (-row[1], -row[2], row[0]))
    return ranked[:top_n]
