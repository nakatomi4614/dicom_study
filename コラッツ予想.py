from matplotlib import pyplot as plt
import numpy as np
from numba import jit

k = int(input("コラッツ予想したい整数を入力："))
calc_data = np.empty(0, dtype=int)
while k != 1:
    if k % 2 == 0:
        k = k / 2
    else:
        k = k * 3 + 1
        calc_data = np.append(calc_data, k)
print(calc_data)
x = np.arange(1, len(calc_data) + 1)
plt.plot(x, calc_data)
plt.show( )
