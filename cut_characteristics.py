import numpy as np

from lib.transformations import *
from lib.lithium_niobate import d


def cut_characteristics(transform, axis):
    assert transformation.shape == (3, 3)
    rotor = rotors[axis]
    E = np.transpose(np.array([0, 0, 0]))
    E[axis] = 1
    # strain = N(A) * d' * A' * E
    strain = N(transform) @ np.transpose(d) @ np.transpose(transform) @ E
    strain = np.array(strain, dtype=np.float64)
    shear_vec = strain[3:]
    long_parr = strain[axis]
    def long_perp_rotor(x): return (N(rotor(x)) @ strain)[(axis-1) % 3]
    return (shear_vec, long_parr, long_perp_rotor)


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

        (shear_vec, long_parr, long_perp_rotor) = cut_characteristics(rotor(angle), axis)

        E = np.transpose(np.array([0, 0, 0]))
        E[axis] = 1

        shear_parr = np.dot(E, shear_vec)
        shear_perp = np.linalg.norm(np.cross(E, shear_vec))

        Xs = np.linspace(0, 2*np.pi, 1000)
        long_perp_max_angle = Xs[np.argmax(np.abs(long_perp_rotor(Xs)))]
        long_perp_max = long_perp_rotor(long_perp_max_angle)
        long_perp_min = long_perp_rotor(long_perp_max_angle + np.pi/2)

        print(f"%d° %s-Cut:" % (angle/np.pi*180, axis_name))
        print(f"|    Shear Perp:%8.2f pm/V" % (shear_perp * 1e12))
        print(f"|    Shear Parr:%8.2f pm/V" % (shear_parr * 1e12))
        print(f"| Long Max Perp:%8.2f pm/V" % (long_perp_max * 1e12))
        print(f"| Long Min Perp:%8.2f pm/V" % (long_perp_min * 1e12))
        print(f"|     Long Parr:%8.2f pm/V" % (long_parr * 1e12))
        print()


if __name__ == "__main__":
    main()
