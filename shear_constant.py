import numpy as np

from lib.transformations import *
from lib.lithium_niobate import d


def main():
    T = 0
    t = T / 180 * np.pi
    A = Az
    dp = A(t) @ d @ np.transpose(N(A(t)))
    print(f"{T}° X-Cut:")
    print(f"       Shear: {np.sqrt(np.power(dp[0,4],2) + np.power(dp[0,5],2))}")
    print(f"    Deform Y: {dp[0,1]}")
    print(f"    Deform Z: {dp[0,2]}")
    print()
    A = Ax
    dp = A(t) @ d @ np.transpose(N(A(t)))
    print(f"{T}° Y-Cut:")
    print(f"       Shear: {np.sqrt(np.power(dp[1,3],2) + np.power(dp[1,5],2))}")
    print(f"    Deform X: {dp[1,0]}")
    print(f"    Deform Z: {dp[1,2]}")
    print()
    A = Ay
    dp = A(t) @ d @ np.transpose(N(A(t)))
    print(f"{T}° Z-Cut:")
    print(f"       Shear: {np.sqrt(np.power(dp[2,3],2) + np.power(dp[2,4],2))}")
    print(f"    Deform X: {dp[2,0]}")
    print(f"    Deform Y: {dp[2,1]}")
    print()


if __name__ == "__main__":
    main()
