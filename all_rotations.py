import numpy as np
from numpy import pi
from matplotlib import pyplot as plt
from skimage import measure
import plotly.graph_objects as go

from transformations import *
from lithium_niobate import d


def main():
    def tp(x): return np.transpose(x)

    def dp(a, b, y): return Az(
        a) @ Ax(b) @ Az(y) @ d @ tp(N(Az(y))) @ tp(N(Ax(b))) @ tp(N(Az(a)))

    a, b, y = pi*np.mgrid[0:2:100j, 0:1:100j, 0:2:100j]
    dp31 = dp(a, b, y)[2, 0]
    dp32 = dp(a, b, y)[2, 1]
    dp33 = dp(a, b, y)[2, 2]
    dp34 = dp(a, b, y)[2, 3]
    dp35 = dp(a, b, y)[2, 4]
    dp36 = dp(a, b, y)[2, 5]
    shear = np.sqrt(
        np.power(dp34, 2) + np.power(dp35, 2))
    print("Generated Matricies")

    iso1 = go.Isosurface(
        x=a.flatten()/pi*180,
        y=b.flatten()/pi*180,
        z=y.flatten()/pi*180,
        colorscale='Reds',
        value=shear.flatten(),
        opacity=0.9,
        surface_count=1,
        isomin=8e-11,
        isomax=8e-11, caps=dict(x_show=False, y_show=False,  z_show=False),
    )
    iso2 = go.Isosurface(
        x=a.flatten()/pi*180,
        y=b.flatten()/pi*180,
        z=y.flatten()/pi*180,
        colorscale='BlueRed',
        value=dp36.flatten(),
        opacity=0.6,
        surface_count=1,
        isomin=0,
        isomax=0, caps=dict(x_show=False, y_show=False,  z_show=False),
    )

    fig = go.Figure(data=[iso1, iso2])

    fig.show()

    # iso_val = 8.01e-11
    # verts, faces, _, _ = measure.marching_cubes(
    #     shear, iso_val, spacing=(360/60, 180/30, 360/60))

    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    # ax.plot_trisurf(verts[:, 0], verts[:, 1], faces, verts[:, 2],
    #                 cmap='Spectral', lw=1)
    # plt.show()


if __name__ == "__main__":
    main()
