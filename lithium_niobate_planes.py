import numpy as np
from dataclasses import dataclass
from lib.lithium_niobate import a, c
from cut_characteristics import cut_characteristics


@dataclass
class VectorCartesian:
    x: float
    y: float
    z: float

    def vec(self) -> np.ndarray:
        return np.array([[self.x], [self.y], [self.z]])


@dataclass
class VectorLattice:
    h: float
    k: float
    l: float

    def vec(self) -> np.ndarray:
        return np.array([[self.h], [self.k], [self.l]])


lattice_mat = np.array([
    [a,  -a/2,           0],
    [0,  np.sqrt(3)/2*a, 0],
    [0,  0,              c]])


def lattice_to_cart(vec: VectorLattice) -> VectorCartesian:
    mat = lattice_mat
    out: np.ndarray = np.matmul(mat, vec.vec())
    return VectorCartesian(out[0, 0], out[1, 0], out[2, 0])


def cart_to_lattice(vec: VectorCartesian) -> VectorLattice:
    mat = np.linalg.inv(lattice_mat)
    out: np.ndarray = np.matmul(mat, vec.vec())
    return VectorLattice(out[0, 0], out[1, 0], out[2, 0])


def find_close_lattice_point(vec: VectorLattice):
    vec = vec.vec()
    vec /= np.max(np.abs(vec))
    min_err = np.Inf
    min_err_vec = None
    for i in range(1, 50):
        err = np.sum(np.power(i*vec - np.round(i*vec), 2))
        if min_err_vec is None or err < min_err:
            min_err = err
            min_err_vec = vec*i
    return VectorLattice(*[v[0] for v in min_err_vec])


def main():
    cut = cut_characteristics(("X", 0))
    angle = cut['shear_ang'][0]/np.pi*180
    print(f"Axis Angle: %0.2f degrees from %s towards %s" % (angle, *cut['shear_ang'][1:]))

    cart_vec = VectorCartesian(0, 0, 0)
    cart_vec.y = np.sin(cut['shear_ang'][0])
    cart_vec.z = np.cos(cut['shear_ang'][0])
    print("Starting Vector: ", cart_vec)

    lattice = cart_to_lattice(cart_vec)
    print("Lattice Vector: ", find_close_lattice_point(lattice))

    cvec = lattice_to_cart(lattice)
    print("Converted Back: ", cvec)
    angle = np.arctan2(cvec.y, cvec.z)/np.pi*180
    print(angle)
    angle = angle if abs(angle) <= 90 else angle % 360-180
    print(angle)


if __name__ == "__main__":
    main()
