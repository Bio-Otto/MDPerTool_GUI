import numpy as np

array1 = np.array([200., 300., 400., 500., 1., 2.])
array2 = np.array([200., 300., 400., 500., 501.])

# Pad the shorter array with zeros to make their shapes equal
length = max(len(array1), len(array2))
array1 = np.pad(array1, (0, length - len(array1)))
array2 = np.pad(array2, (0, length - len(array2)))

# Combine arrays using numpy.where
result = np.where(array2 > array1, array2, array1)
print(result)
