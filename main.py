import numpy as np
from matplotlib import pyplot as plt
import itertools as it


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


def Ax(t):
    return np.array([[1, 0, 0], [0, np.cos(t), np.sin(t)], [0, -np.sin(t), np.cos(t)]])


def Ay(t):
    return np.array([[np.cos(t), 0, np.sin(t)], [0, 1, 0], [-np.sin(t), 0, np.cos(t)]])


def Az(t):
    return np.array([[np.cos(t), np.sin(t), 0], [-np.sin(t), np.cos(t), 0], [0, 0, 1]])


def main():
    d15 = 69.2e-12
    d22 = 20.8e-12
    d31 = -0.85e-12
    d33 = 6e-12
    d = np.array([[0, 0, 0, 0, d15, -2*d22], [-d22, d22,
                 0, d15, 0, 0], [d31, d31, d33, 0, 0, 0]], dtype=np.float64)
    
    def plot(dp):
        xs = np.linspace(0, np.pi, 500)
        for r,c in it.product([0,1,2],[0,1,2,3,4,5]):
            plt.subplot(3,6, r*6+c+1)
            try:
                plt.plot(xs/np.pi*180,np.abs(dp(xs)[r,c]))
            except:
                plt.plot(xs, xs*0)
            plt.ylim(0,9e-11)
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
    T = 36
    t = T / 180 * np.pi
    A = Az
    dp = A(t) @ d @ np.transpose(N(A(t)))
    print(f"{T}° X-Cut: {np.sqrt(np.power(dp[0,4],2) + np.power(dp[0,5],2))}")
    A = Ax
    dp = A(t) @ d @ np.transpose(N(A(t)))
    print(f"{T}° Y-Cut: {np.sqrt(np.power(dp[1,3],2) + np.power(dp[1,5],2))}")
    A = Ay
    dp = A(t) @ d @ np.transpose(N(A(t)))
    print(f"{T}° Z-Cut: {np.sqrt(np.power(dp[2,3],2) + np.power(dp[2,4],2))}")



if __name__ == "__main__":
    main()
