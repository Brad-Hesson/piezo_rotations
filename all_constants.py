import numpy as np
from matplotlib import pyplot as plt
import itertools as it

from transformations import *
from lithium_niobate import d


def main():
    def plot(dp):
        xs = np.linspace(0, np.pi, 500)
        for r, c in it.product([0, 1, 2], [0, 1, 2, 3, 4, 5]):
            plt.subplot(3, 6, r*6+c+1)
            try:
                plt.plot(xs/np.pi*180, np.abs(dp(xs)[r, c]))
            except:
                plt.plot(xs, xs*0)
            plt.ylim(0, 9e-11)
            plt.title(f"dp{r+1}{c+1}")
    # -----------------------------------------------------------
    plt.figure()
    plot(lambda t: Az(t) @ d @ np.transpose(N(Az(t))))
    plt.suptitle("X rotated about Z")
    # -----------------------------------------------------------
    plt.figure()
    plot(lambda t: Ax(t) @ d @ np.transpose(N(Ax(t))))
    plt.suptitle("Y rotated about X")
    # -----------------------------------------------------------
    plt.figure()
    plot(lambda t: Ay(t) @ d @ np.transpose(N(Ay(t))))
    plt.suptitle("Z rotated about Y")
    # -----------------------------------------------------------
    plt.show()


if __name__ == "__main__":
    main()
