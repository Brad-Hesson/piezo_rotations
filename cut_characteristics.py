import numpy as np

from lib.transformations import *
from lib.lithium_niobate import d


def main():
    axes = {"X": 0, "Y": 1, "Z": 2}
    trans = {"X": Az, "Y": Ax, "Z": Ay}
    cuts = [
        ("X", 0),
        ("Y", 0),
        ("Y", 36),
        ("Y", 41),
        ("Y", 64),
        ("Y", 128),
        ("Z", 0),
    ]
    for cut in cuts:
        T = cut[1]
        t = T / 180 * np.pi
        A = trans[cut[0]]
        dp = A(t) @ d @ np.transpose(N(A(t)))
        axis = axes[cut[0]]
        other1 = (axis+1) % 3
        other2 = (axis+2) % 3
        shear = np.sqrt(
            np.power(dp[axis, other1 + 3], 2) +
            np.power(dp[axis, other2 + 3], 2)
        )
        long_perp = [dp[axis, other1], dp[axis, other2]]
        long_perp_max = long_perp[np.argmax(np.abs(long_perp))]
        long_perp_min = long_perp[np.argmin(np.abs(long_perp))]
        print(f"{T}° {cut[0]}-Cut:")
        print(f"|    Shear Perp:%8.2f pm/V" % (shear * 1e12))
        print(f"|    Shear Parr:%8.2f pm/V" % (dp[axis, axis+3] * 1e12))
        print(f"| Long Max Perp:%8.2f pm/V" % (long_perp_max * 1e12))
        print(f"| Long Min Perp:%8.2f pm/V" % (long_perp_min * 1e12))
        print(f"|     Long Parr:%8.2f pm/V" % (dp[axis, axis] * 1e12))
        print()


if __name__ == "__main__":
    main()
