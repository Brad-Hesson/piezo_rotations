import numpy as np


# Piezoelectric constants [pm/V] and matrix (source: https://www.bostonpiezooptics.com/lithium-niobate)
d15: float = 69.2e-12
d22: float = 20.8e-12
d31: float = -0.85e-12
d33: float = 6e-12
d: np.ndarray = np.array(
    [
        [0, 0, 0, 0, d15, -2*d22],
        [-d22, d22, 0, d15, 0, 0],
        [d31, d31, d33, 0, 0, 0]
    ],
    dtype=np.float64)

# Lattice constants [m] (source: http://www.mt-berlin.com/frames_cryst/descriptions/lnb_lta.htm)
a: float = 5.148e-10
c: float = 13.863e-10
