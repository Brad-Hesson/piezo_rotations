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
        long_perp1 = dp[axis, other1]
        long_perp2 = dp[axis, other2]
        print(f"{T}° {cut[0]}-Cut:")
        print(f"|    Shear Perp:  {shear}")
        print(f"|    Shear Parr:  {dp[axis, axis+3]}")
        print(f"| Long Max Perp:  {max(long_perp1, long_perp2)}")
        print(f"| Long Min Perp:  {min(long_perp1, long_perp2)}")
        print(f"|     Long Parr:  {dp[axis,axis]}")
        print()


if __name__ == "__main__":
    main()
