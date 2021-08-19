import numpy as np

d15 = 69.2e-12
d22 = 20.8e-12
d31 = -0.85e-12
d33 = 6e-12
d = np.array([[0, 0, 0, 0, d15, -2*d22], [-d22, d22,
                                          0, d15, 0, 0], [d31, d31, d33, 0, 0, 0]], dtype=np.float64)