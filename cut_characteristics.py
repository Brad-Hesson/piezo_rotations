from typing import Tuple
import numpy as np

from lib.transformations import *
from lib.lithium_niobate import d


def cut_characteristics(cut: Tuple[str, float]):
    axes_ind = {"X": 0, "Y": 1, "Z": 2}
    axes_name = {v: k for k, v in axes_ind.items()}
    transform = {"X": Az, "Y": Ax, "Z": Ay}

    rot_angle = cut[1] / 180 * np.pi
    A = transform[cut[0]]
    dp = A(rot_angle) @ d @ np.transpose(N(A(rot_angle)))
    axis = axes_ind[cut[0]]
    other1 = (axis+1) % 3
    other2 = (axis+2) % 3
    shear = np.sqrt(
        np.power(dp[axis, other1 + 3], 2) +
        np.power(dp[axis, other2 + 3], 2)
    )
    shear_angle = np.arctan2(dp[axis, other1 + 3],
                             dp[axis, other2 + 3])
    shear_angle = shear_angle if abs(
        shear_angle) <= np.pi/2 else shear_angle % (2*np.pi)-np.pi
    long_perp = [dp[axis, other1], dp[axis, other2]]
    long_perp_max = long_perp[np.argmax(np.abs(long_perp))]
    long_perp_min = long_perp[np.argmin(np.abs(long_perp))]
    return {'shear_perp': shear,
            'sher_parr': dp[axis, axis+3],
            'shear_ang': (shear_angle, axes_name[other2],
                          axes_name[other1]),
            'long_perp_max': long_perp_max,
            'long_perp_min': long_perp_min,
            'long_parr': dp[axis, axis]}


def main():
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
        cut_chars = cut_characteristics(cut)
        print(f"{cut[1]}° {cut[0]}-Cut:")
        print(f"|    Shear Perp:%8.2f pm/V" % (cut_chars['shear_perp'] * 1e12))
        print(f"|    Shear Parr:%8.2f pm/V" % (cut_chars['sher_parr'] * 1e12))
        print(f"|   Shear Angle:%8.2f degrees from %s towards %s" %
              (cut_chars['shear_ang'][0]/np.pi*180, *cut_chars['shear_ang'][1:]))
        print(f"| Long Max Perp:%8.2f pm/V" %
              (cut_chars['long_perp_max'] * 1e12))
        print(f"| Long Min Perp:%8.2f pm/V" %
              (cut_chars['long_perp_min'] * 1e12))
        print(f"|     Long Parr:%8.2f pm/V" % (cut_chars['long_parr'] * 1e12))
        print()


if __name__ == "__main__":
    main()
