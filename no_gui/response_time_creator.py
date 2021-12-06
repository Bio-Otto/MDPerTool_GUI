import numpy as np
import pandas as pd


def getResidueResponseTimes(referenceName, perturbedName, outputName='responseTimes.csv'):
    """
    calculate residue response times and write them into a file
    """

    refEnergies = pd.read_csv(referenceName)
    perEnergies = pd.read_csv(perturbedName)

    numEnergies = len(refEnergies.columns)
    numFrames = len(refEnergies.values)

    energyDiff = []

    for i in range(0, numEnergies):
        columnName = refEnergies.columns[i]
        diff = refEnergies[columnName] - perEnergies[columnName]
        for k in range(0, numFrames):
            if abs(diff[k]) < 0.01:
                # print(diff[k])
                diff[k] = 0
        energyDiff.append(diff)

    residueResponseTimes = []

    for residueDiff in energyDiff:
        # For residues with no response, assign a value of length of trajectory.
        # For residues with response, assign the first frame with a non-zero element.
        if np.sum(residueDiff) != 0:
            responseTime = np.array(residueDiff).nonzero()[0][0]
            residueResponseTimes.append(responseTime)
        else:
            residueResponseTimes.append(0)

    residueResponseTimes = np.asarray(residueResponseTimes)
    print("NUMBER OF RESIDUES: %s" % numEnergies)
    print("NUMBER OF FRAME: %s" % numFrames)
    print(residueResponseTimes)
    np.savetxt(outputName, residueResponseTimes, delimiter=',')


# getResidueResponseTimes('reference_energy_file.csv', 'modified_energy_file.csv')
