import numpy as np
from matplotlib import pyplot as plt

from lib.transformations import *
from lib.lithium_niobate import d


def cut_characteristics(d_matrix, transform, axis):
    assert transform.shape == (3, 3)
    E = np.transpose(np.array([0, 0, 0]))
    E[axis] = 1
    # strain = N(A) * d' * A' * E
    strain = N(transform) @ np.transpose(d_matrix) @ np.transpose(transform) @ E
    strain = np.array(strain, dtype=np.float64)
    shear_vec = strain[3:]
    long_parr = strain[axis]
    print(strain[:3]*1e12)
    print(strain[3:]*1e12)
    rotor = rotors[axis]
    def long_perp_rotor(x): return (N(rotor(x)) @ strain)[(axis-1) % 3]
    return (shear_vec, long_parr, long_perp_rotor)


def cut_params(axis, angle, shear_vec, long_parr, long_perp_rotor):
    E = np.transpose(np.array([0, 0, 0]))
    E[axis] = 1

    shear_parr = np.dot(E, shear_vec)
    shear_perp = np.linalg.norm(np.cross(E, shear_vec))
    indeces = [(axis-1) % 3, (axis-2) % 3]
    shear_perp_angle = np.arctan2(*np.cross(E, shear_vec)[indeces])
    shear_perp_angle = shear_perp_angle if abs(
        shear_perp_angle) <= 90 else shear_perp_angle % 360-180

    Xs = np.linspace(0, 2*np.pi, 1000)
    long_perp_max_angle = Xs[np.argmax(np.abs(long_perp_rotor(Xs)))]
    long_perp_max = long_perp_rotor(long_perp_max_angle)
    long_perp_min = long_perp_rotor(long_perp_max_angle + np.pi/2)

    return (long_perp_max, long_perp_min, long_perp_max_angle, long_parr, shear_parr, shear_perp, shear_perp_angle)


def main():
    axes_ind = {"X": 0, "Y": 1, "Z": 2}
    cuts = [
        ("X", 0),
        ("Y", 0),
        ("Y", 36),
        ("Y", 41),
        ("Y", 64),
        ("Y", 128),
        ("Z", 0),
    ]
    for (axis_name, angle) in cuts:
        axis = axes_ind[axis_name]
        rotor = rotors[(axis-1) % 3]
        angle *= np.pi/180

        print(f"%d° %s-Cut:" % (angle/np.pi*180, axis_name))

        (shear_vec,
         long_parr,
         long_perp_rotor) = cut_characteristics(d, rotor(angle), axis)
        (long_perp_max,
         long_perp_min,
         long_perp_angle,
         long_parr,
         shear_parr,
         shear_perp,
         shear_perp_angle,) = cut_params(axis, angle, shear_vec, long_parr, long_perp_rotor)

        print(f"|        Long Parr:%8.2f pm/V" % (long_parr * 1e12))
        print(f"|    Long Perp Max:%8.2f pm/V" % (long_perp_max * 1e12))
        print(f"|    Long Perp Min:%8.2f pm/V" % (long_perp_min * 1e12))
        print(f"|  Long Perp Angle:%8.2f°" % (long_perp_angle * 180 / np.pi))
        print(f"|       Shear Parr:%8.2f pm/V" % (shear_parr * 1e12))
        print(f"|       Shear Perp:%8.2f pm/V" % (shear_perp * 1e12))
        print(f"| Shear Perp Angle:%8.2f°" % (shear_perp_angle * 180 / np.pi))
        print()

    Xs = np.linspace(0, 2*np.pi, 1000)
    plt.plot(Xs, long_perp_rotor(Xs))
    plt.show()


if __name__ == "__main__":
    main()
