import numpy as np
import pandas as pd


RESPONSE_THRESHOLD = 0.01


def get_residue_response_times(reference_name, perturbed_name, output_name='responseTimes.csv'):
    """Calculate the first responsive frame for each residue and write it to disk."""

    reference_energies = pd.read_csv(reference_name)
    perturbed_energies = pd.read_csv(perturbed_name)

    energy_diff = reference_energies.sub(perturbed_energies)
    energy_diff = energy_diff.mask(energy_diff.abs() < RESPONSE_THRESHOLD, 0.0)

    num_frames = len(energy_diff.index)
    residue_response_times = []

    for column_name in energy_diff.columns:
        residue_diff = energy_diff[column_name].to_numpy()
        non_zero_indices = np.flatnonzero(residue_diff)

        if non_zero_indices.size:
            residue_response_times.append(int(non_zero_indices[0]))
        else:
            residue_response_times.append(num_frames)

    response_time_array = np.asarray(residue_response_times)
    np.savetxt(output_name, response_time_array, delimiter=',')
    return response_time_array


def getResidueResponseTimes(referenceName, perturbedName, outputName='responseTimes.csv'):
    """Backward-compatible wrapper for legacy callers."""

    return get_residue_response_times(referenceName, perturbedName, outputName)


__all__ = ["RESPONSE_THRESHOLD", "get_residue_response_times", "getResidueResponseTimes"]


# getResidueResponseTimes('reference_energy_file.csv', 'modified_energy_file.csv')
