import numpy as np
from dataclasses import dataclass
from materials.lithium_niobate import a, c
from cut_characteristics import cut_characteristics


@dataclass
class CVec:
    x: float
    y: float
    z: float

    def vec(self) -> np.ndarray:
        return np.array([[self.x], [self.y], [self.z]])


@dataclass
class LVec:
    h: float
    k: float
    l: float

    def vec(self) -> np.ndarray:
        return np.array([[self.h], [self.k], [self.l]])
    
    def inverted(self):
        return LVec(1/self.h, 1/self.k, 1/self.l)


lattice_mat = np.array([
    [a,  -a/2,           0],
    [0,  np.sqrt(3)/2*a, 0],
    [0,  0,              c]])


def lattice_to_cart(vec: LVec) -> CVec:
    mat = lattice_mat
    out: np.ndarray = np.matmul(mat, vec.vec())
    return CVec(out[0, 0], out[1, 0], out[2, 0])


def cart_to_lattice(vec: CVec) -> LVec:
    mat = np.linalg.inv(lattice_mat)
    out: np.ndarray = np.matmul(mat, vec.vec())
    return LVec(out[0, 0], out[1, 0], out[2, 0])


def find_close_lattice_point(vec: LVec, search_len = 20):
    vec = vec.vec()
    vec /= np.max(np.abs(vec))
    min_err = np.Inf
    min_err_vec = None
    for i in range(1, search_len):
        err = np.sum(np.power(i*vec - np.round(i*vec), 2))
        if min_err_vec is None or err < min_err:
            min_err = err
            min_err_vec = np.round(i*vec)
    return LVec(*[v[0] for v in min_err_vec]), np.sqrt(min_err)


def main():
    from materials.lithium_niobate import d
    from lib.transformations import rotors

    axis = 0
    angle = 0

    rotor = rotors[(axis-1)%3]
    transform = rotor(angle/180*np.pi)
    shear_cutb, _, _ = cut_characteristics(d, transform, axis)

    E = np.transpose(np.array([0, 0, 0]))
    E[axis] = 1
    shear_parr_cutb = np.dot(E, shear_cutb) * E
    shear_perp_cutb = shear_cutb - shear_parr_cutb
    shear_perp_crystalb = np.transpose(transform) @ shear_perp_cutb
    shear_perp_crystalb = CVec(*shear_perp_crystalb)
    shear_perp_latticeb = cart_to_lattice(shear_perp_crystalb)

    lattice_vector, err = find_close_lattice_point(shear_perp_latticeb, 100)
    lattice_plane, _ = find_close_lattice_point(lattice_vector.inverted(), 1000)

    print(" Lattice Vector:", lattice_vector, "+-", err)
    print("  Lattice Plane:", lattice_plane)


if __name__ == "__main__":
    main()
