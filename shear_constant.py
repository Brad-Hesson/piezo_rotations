import numpy as np

from lib.transformations import *
from lib.lithium_niobate import d


def main():
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
