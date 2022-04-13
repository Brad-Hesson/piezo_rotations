import numpy as np

# This is a function that takes a 3d rotation matrix and
# converts it into a 6-vector rotation matrix as detailed
# in `Section 7` of `papers/Crystal Properties.pdf`
def N(a):
    assert(a.shape == (3, 3))
    s0r = [[0, 0, 0], [1, 1, 1], [2, 2, 2]]
    s1r = [[2, 2, 2], [0, 0, 0], [1, 1, 1]]
    s2r = [[1, 1, 1], [2, 2, 2], [0, 0, 0]]
    s0c = np.transpose(s0r)
    s1c = np.transpose(s1r)
    s2c = np.transpose(s2r)
    na11 = a[s0r, s0c] * a[s0r, s0c]
    na12 = a[s0r, s1c] * a[s0r, s2c]
    na21 = a[s1r, s0c] * a[s2r, s0c] * 2
    na22 = a[s1r, s1c] * a[s2r, s2c] + a[s1r, s2c] * a[s2r, s1c]
    return np.vstack([np.hstack([na11, na12]), np.hstack([na21, na22])])


# Generates a matrix that will rotate about the x-axis by `t` radians
def Ax(t):
    return np.array([[1, 0, 0], [0, np.cos(t), np.sin(t)], [0, -np.sin(t), np.cos(t)]], dtype=object)


# Generates a matrix that will rotate about the y-axis by `t` radians
def Ay(t):
    return np.array([[np.cos(t), 0, np.sin(t)], [0, 1, 0], [-np.sin(t), 0, np.cos(t)]], dtype=object)


# Generates a matrix that will rotate about the z-axis by `t` radians
def Az(t):
    return np.array([[np.cos(t), np.sin(t), 0], [-np.sin(t), np.cos(t), 0], [0, 0, 1]], dtype=object)

# a list of `rotor` matrix generators, this is useful in other functions
# that need to programatically generate rotations about different axes
rotors = [Ax, Ay, Az]
