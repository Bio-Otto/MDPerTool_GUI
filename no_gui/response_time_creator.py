import numpy as np
import pandas as pd


def getResidueResponseTimes(referenceName, perturbedName, outputName='responseTimes.csv'):
    # calculate residue response times and write them into a file

    refEnergies = pd.read_csv(referenceName)
    perEnergies = pd.read_csv(perturbedName)

    numEnergies = len(refEnergies.columns)
    print(numEnergies)
    numFrames = len(refEnergies.values)
    print(numFrames)

    energyDiff = list()

    for i in range(0, numEnergies):
        columnName = refEnergies.columns[i]


        diff = refEnergies[columnName] - perEnergies[columnName]
        for k in range(0, numFrames):
            if abs(diff[k]) < 0.01:
                # print(diff[k])
                diff[k] = 0
        energyDiff.append(diff)

    residueResponseTimes = list()

    np.savetxt('TIME_OMM.csv', energyDiff, delimiter=',')

    for residueDiff in energyDiff:
        # For residues with no response, assign a value of length of trajectory.
        # For residues with response, assign the first frame with a non-zero element.
        if int(np.sum(residueDiff)) != 0:
            # responseTime = np.to_numpy().nonzero(residueDiff)[0][0]
            responseTime = np.array(residueDiff).nonzero()[0][0]
            print(responseTime)

            residueResponseTimes.append(responseTime)
        else:
            residueResponseTimes.append(0)


    residueResponseTimes = np.asarray(residueResponseTimes)
    print(residueResponseTimes)
    np.savetxt(outputName, residueResponseTimes, delimiter=',')


# getResidueResponseTimes('reference_energy_file.csv', 'modified_energy_file.csv')

# getResidueResponseTimes('3IBF_G188L_2_50ns_200ps_ref_energies_Total.csv', '3IBF_G188L_2_50ns_200ps_290590_x4_perturb_energies_Total.csv')