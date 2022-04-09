import numpy as np
from matplotlib import pyplot as plt

from lib.transformations import *

# This function takes
# `d_matrix`: the matrix of piezoelectric coefficients for the material in question
# `transform`: the 3d rotation that represents the cutting plane of the crystal
# `axis`: the axis across which the electric field would be applied
# 
# it returns
# `shear_vec`: a vector of the shear response
# `long_parr`: how much the crystal expands (or contracts if negative) in the axis parallel to the applied electric field 
# `long_perp_rotor`: a function that tells you how much the crystal expands (or contracts if negative) perpendicular to the
#   applied electric field at any given angle about that axis
# 
def cut_characteristics(d_matrix, transform, axis):
    # throw an error if the transform is not the correct shape
    assert transform.shape == (3, 3)
    # create the electric field vector based on the given `axis` variable
    E = np.transpose(np.array([0, 0, 0]))
    E[axis] = 1
    # here `A` is the given `transform` variable which is the rotation matrix
    # multiplying on the left by `A` converts from the crystal basis to the cut basis
    # multiplying on the left by `A'` (the transpose or inverse of `A`) converts from the cut basis to the crystal basis
    # `@` is matrix multiplication
    # d_t = A @ d @ N(A)' 
    # d_t' = N(A) @ d' @ A'
    # strain = d_t' @ E
    # strain = N(A) @ d' @ A' @ E
    # see `2. Theoretical Background` in `papers/LiNbO3 Orientation Dependence.pdf`
    strain = N(transform) @ np.transpose(d_matrix) @ np.transpose(transform) @ E
    strain = np.array(strain, dtype=np.float64)
    # `strain` is now a 6-vector representing the strain matrix of the cut crystal as follows
    # |s1 s6 s5|
    # |s6 s2 s4|
    # |s5 s4 s3|
    # `shear_vec` is |s4 s5 s6|
    # `long_parr` is s[axis]
    # `rotor` is a function that gives information on the perpendicular component of the longitudinal motion
    shear_vec = strain[3:]
    long_parr = strain[axis]
    rotor = rotors[axis]
    def long_perp_rotor(x): return (N(rotor(x)) @ strain)[(axis-1) % 3]
    return (shear_vec, long_parr, long_perp_rotor)


# this function takes in the information from `cut_characteristics` and returns more fine
# grained information.
#           `long_parr`: the longitudinal motion perpendicular to the cutting plane
#       `long_perp_max`: the maximum signed longitudinal motion in the cutting plane
#       `long_perp_min`: the minimum signed longitudinal motion in the cutting plane
# `long_perp_max_angle`: the angle of the maximum longitudial motion in the cutting plane
#          `shear_parr`: the shear motion around the axis perpendicular to the cutting plane
#          `shear_perp`: the shear motion about an axis in the cutting plane
#    `shear_perp_angle`: the angle of the perpendicular shear axis in the cutting plane
def cut_params(axis, angle, shear_vec, long_parr, long_perp_rotor):
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


# takes a bunch of common lithium niobate cuts and calculates their parameters of motion
# using the functions above
def main():
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
