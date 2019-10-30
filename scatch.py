import numpy as np
import matplotlib.pyplot as plt


x = np.arange(100)

start = (0, 30)
stop = (50, 1)
offset = -stop([0])

A = start[1]

y = A*((x + offset)**2) + stop[1]

(y - stop[1])/((x - offset)**2)