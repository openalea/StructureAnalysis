import example
import numpy as np


# Fonction add
a = example.add(1,5)
b = example.add()
print(a, b)  # 6, 3


# Fonction fois_deux_bis
data = np.array([[[1,2,3], [4,5,6]], [[1,2,3], [4,5,6]]], dtype=np.float64)
print(data.shape)      # (2,2,3)

x = example.fois_deux_bis(data)
print(data)            # data*2
print(x, type(x))      # 42.0  <class 'float'>




# Class Pet
