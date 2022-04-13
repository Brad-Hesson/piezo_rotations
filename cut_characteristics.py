import numpy as np
from matplotlib import pyplot as plt

from lib.transformations import *


def cut_characteristics(d_matrix: np.ndarray, transform: np.ndarray, axis: int) -> tuple:
    assert transform.shape == (3, 3)
    E = np.transpose(np.array([0, 0, 0]))
    E[axis] = 1
    # A': Cut Basis     -> Crystal Basis
    # A : Crystal Basis -> Cut Basis
    # strain = N(A) * d' * A' * E
    strain = N(transform) @ np.transpose(d_matrix) @ np.transpose(transform) @ E
    strain = np.array(strain, dtype=np.float64)
    shear_vec = strain[3:]
    long_parr = strain[axis]
    rotor = rotors[axis]
    def long_perp_rotor(x: float): return (N(rotor(x)) @ strain)[(axis-1) % 3]
    return (shear_vec, long_parr, long_perp_rotor)


def cut_params(axis: int, angle: float, shear_vec: np.ndarray, long_parr: float, long_perp_rotor) -> tuple:
    E = np.transpose(np.array([0, 0, 0]))
    E[axis] = 1

    shear_parr = np.dot(E, shear_vec)
    shear_perp = np.linalg.norm(shear_vec - shear_parr)
    indeces = [(axis-2) % 3, (axis-1) % 3]
    shear_perp_angle = np.arctan2(*shear_vec[indeces])
    shear_perp_angle = shear_perp_angle if abs(
        shear_perp_angle) <= np.pi/2 else shear_perp_angle % (2*np.pi)-np.pi

    Xs = np.linspace(0, 2*np.pi, 1000)
    long_perp_max_angle = Xs[np.argmax(np.abs(long_perp_rotor(Xs)))]
    long_perp_max = long_perp_rotor(long_perp_max_angle)
    long_perp_min = long_perp_rotor(long_perp_max_angle + np.pi/2)

    return (long_parr, long_perp_max, long_perp_min, long_perp_max_angle, shear_parr, shear_perp, shear_perp_angle)


def main() -> None:
    axes_ind = {"X": 0, "Y": 1, "Z": 2}
    # Add cuts to the list below to get their parameters.  The
    # final cut in the list will have it's longitudinal
    # perpendicular function plotted.
    cuts = [
        ("X", 0),
        ("Y", 0),
        ("Y", 36),
        ("Y", 41),
        ("Y", 64),
        ("Y", 128),
        ("Z", 0),
    ]
    from materials.lithium_niobate import d
    for (axis_name, angle) in cuts:
        axis = axes_ind[axis_name]
        rotor = rotors[(axis-1) % 3]
        angle *= np.pi/180

        print(f"%d° %s-Cut:" % (angle/np.pi*180, axis_name))

        (shear_vec,
         long_parr,
         long_perp_rotor) = cut_characteristics(d, rotor(angle), axis)
        (long_parr,
         long_perp_max,
         long_perp_min,
         long_perp_angle,
         shear_parr,
         shear_perp,
         shear_perp_angle,) = cut_params(axis, angle, shear_vec, long_parr, long_perp_rotor)

        print("%.4e" % shear_vec[0])
        print("%.4e" % shear_vec[1])
        print("%.4e" % shear_vec[2])
        print(f"|            Longitudinal Parallel:%8.2f pm/V" % (long_parr * 1e12))
        print(f"|   Longitudinal Perpendicular Max:%8.2f pm/V" % (long_perp_max * 1e12))
        print(f"|   Longitudinal Perpendicular Min:%8.2f pm/V" % (long_perp_min * 1e12))
        print(f"| Longitudinal Perpendicular Angle:%8.2f°" % (long_perp_angle * 180 / np.pi))
        print(f"|                   Shear Parallel:%8.2f pm/V" % (shear_parr * 1e12))
        print(f"|              Shear Perpendicular:%8.2f pm/V" % (shear_perp * 1e12))
        print(f"|        Shear Perpendicular Angle:%8.2f°" % (shear_perp_angle * 180 / np.pi))
        print()

    Xs = np.linspace(0, 2*np.pi, 1000)
    plt.plot(Xs, long_perp_rotor(Xs))
    plt.show()


if __name__ == "__main__":
    main()
