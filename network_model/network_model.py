import random

import numpy as np
import scipy as sp
import scipy.sparse as sparse
import scipy.sparse.linalg as sla
from scipy.special import kv
from typing import Tuple
import numpy as np
import random
from typing import Tuple, Optional

mu_B = 5.788e-5  # Bohr-magneton


class MicroMagnet:
    """
    A class for simulating a micro-magnetic system.

    Args:
        shape (Tuple[int, int]): The shape of the micro-magnetic system.
        seed (Optional[int]): Seed for random number generator (optional).

    Attributes:
        shape (Tuple[int, int]): The shape of the micro-magnetic system.
        config (np.ndarray): The configuration of the micro-magnetic system.
        up_idxs (list): List of indices representing 'up' spins.
        down_idxs (list): List of indices representing 'down' spins.

    Methods:
        set_M(new_M: float) -> None:
            Set the magnetization of the system to a new value.

        get_M() -> float:
            Get the current magnetization of the system.

    Usage:
    mm = MicroMagnet(shape=(5, 5), seed=42)
    mm.set_M(0.2)
    print(mm.get_M())
    """

    def __init__(self, shape: Tuple[int, int], seed: Optional[int] = None):
        """
        Initializes a MicroMagnet instance.

        Args:
            shape (Tuple[int, int]): The shape of the micro-magnetic system.
            seed (Optional[int]): Seed for random number generator (optional).
        """
        self.shape = shape
        self.config = np.ones(shape)
        self.up_idxs = list(range(shape[0] * shape[1]))

        if seed is not None:
            random.seed(seed)

        random.shuffle(self.up_idxs)
        self.down_idxs = []

    def set_M(self, new_M: float) -> None:
        """
        Set the magnetization of the system to a new value.

        Args:
            new_M (float): The new magnetization value.
        """
        assert abs(new_M) <= 1

        to_flip_num = int((new_M - self.get_M()) * (self.shape[0] * self.shape[1]) / 2)

        if to_flip_num > 0:
            for i in range(to_flip_num):
                idx = self.down_idxs.pop()
                self.up_idxs.append(idx)
                self.config[np.unravel_index(idx, self.shape)] = +1
        else:
            to_flip_num = -to_flip_num
            for i in range(to_flip_num):
                idx = self.up_idxs.pop()
                self.down_idxs.append(idx)
                self.config[np.unravel_index(idx, self.shape)] = -1

    def get_M(self) -> float:
        """
        Get the current magnetization of the system.

        Returns:
            float: The magnetization value.
        """
        return (len(self.up_idxs) - len(self.down_idxs)) / (
            self.shape[0] * self.shape[1]
        )


def block_sum(input_array: np.ndarray, M: int) -> np.ndarray:
    """
    Compute the sum of MxM subblocks in a 2D array and return the resulting NxN array.

    Args:
        input_array (np.ndarray): The input 2D array of dimension NMxNM.
        M (int): The size of the subblock.

    Returns:
        np.ndarray: The resulting NxN array where each value is the sum of an MxM subblock.

    Raises:
        ValueError: If the dimensions of the input array are not divisible by M.

    """
    NM = input_array.shape[0]
    N = NM // M
    result_array = np.empty((N, N), dtype=input_array.dtype)

    if NM % M != 0:
        raise ValueError("Dimensions of input_array are not divisible by M.")

    for i in range(N):
        for j in range(N):
            # Extract the MxM subblock
            block = input_array[i * M : (i + 1) * M, j * M : (j + 1) * M]

            # Calculate the sum of the subblock and store it in the result array
            result_array[i, j] = np.sum(block)

    return result_array


def generate_averaging_filter(xi: float) -> np.ndarray:
    """
    Generate averaging filter.

    - xi: float
        Spin-averaging length in domain size.

    Returns:
    - np.ndarray
        The averaging filter.
    """
    domains_N = 41
    points_per_domain = 40
    points_N = points_per_domain * domains_N

    # Generate the grid
    x_ax = np.linspace(-domains_N, domains_N, points_N)
    dx = x_ax[1] - x_ax[0]
    x, y = np.meshgrid(x_ax, x_ax)
    r = np.sqrt(x**2 + y**2)

    # Evaluate the Bessel function on the grid
    K0 = 1 / (2 * np.pi) * kv(0, r / xi)

    # Block integration to get coars-grained values
    averaging_filter = block_sum(K0 * dx**2, points_per_domain)
    averaging_filter /= averaging_filter.sum()

    return averaging_filter


def sgnd_pdf(x, alpha, beta):
    """
    PDF of symmetric generalized normal distribution.
    """

    return (
        beta
        / (2 * alpha * sp.special.gamma(1 / beta))
        * np.exp(-((abs(x) / alpha) ** beta))
    )


def solve_circuit(x_ax, y_ax, rho, deltaV):
    """
    Solve the Laplace equation for a ohmic material with inhomogeneous resistivity.
    """

    Nx = x_ax.shape[0]
    Ny = y_ax.shape[0]

    dx = x_ax[1] - x_ax[0]
    dy = y_ax[1] - y_ax[0]

    rho = rho.flatten()

    I_x = sparse.eye(Nx)
    I_y = sparse.eye(Ny)

    Shpx = sparse.kron(
        I_y,
        sparse.diags(
            [np.ones(Nx - 1)],
            [1],
            shape=(Nx, Nx),
            format="lil",
        ),
    )

    Shmx = sparse.kron(
        I_y,
        sparse.diags(
            [np.ones(Nx - 1)],
            [-1],
            shape=(Nx, Nx),
            format="lil",
        ),
    )

    Shpy = sparse.kron(
        sparse.diags(
            [np.ones(Ny - 1)],
            [1],
            shape=(Ny, Ny),
            format="lil",
        ),
        I_x,
    )

    Shmy = sparse.kron(
        sparse.diags(
            [np.ones(Ny - 1)],
            [-1],
            shape=(Ny, Ny),
            format="lil",
        ),
        I_x,
    )

    A = (
        sparse.diags(1 / rho / dx) @ Shmx
        - (sparse.diags(1 / rho / dx) + sparse.diags(Shpx @ (1 / rho / dx)))
        + sparse.diags(Shpx @ (1 / rho / dx)) @ Shpx
    ) + (
        sparse.diags(1 / rho / dy) @ Shmy
        - (sparse.diags(1 / rho / dy) + sparse.diags(Shpy @ (1 / rho / dy)))
        + sparse.diags(Shpy @ (1 / rho / dy)) @ Shpy
    )
    B = np.zeros((Ny, Nx))

    # Left BC
    for n in range(Ny):
        B[n, 0] = -deltaV / rho[n * Nx] / dx

    # Right BC
    for n in range(Ny):
        A[Nx - 1 + n * Nx, Nx - 1 + n * Nx] += -1 / rho[-1] / dx

    # Bottom BC
    for n in range(Nx):
        A[n, n] += +1 / rho[n] / dy

    B = B.flatten()
    phi = sla.spsolve(sparse.csr_matrix(A), B)
    E_x = -(sparse.eye(Nx * Ny) - Shmx) @ phi / dx
    E_y = -(sparse.eye(Nx * Ny) - Shmy) @ phi / dy
    j_x = E_x / rho
    j_y = E_y / rho
    charge = (sparse.eye(Nx * Ny) - Shmx) @ E_x + (sparse.eye(Nx * Ny) - Shmx) @ E_y

    phi = phi.reshape(Ny, Nx)
    E_x = E_x.reshape(Ny, Nx)
    E_y = E_y.reshape(Ny, Nx)
    j_x = j_x.reshape(Ny, Nx)
    j_y = j_y.reshape(Ny, Nx)
    rho = rho.reshape(Ny, Nx)
    charge = charge.reshape(Ny, Nx)

    return (phi, E_x, E_y, j_x, j_y, rho, charge)


def compute_resistance(x_ax, y_ax, rho):
    """
    Given a resistivity field, calculate the resistance for leads attached at the two ends.
    """

    dx = x_ax[1] - x_ax[0]
    dy = y_ax[1] - y_ax[0]

    (phi, E_x, E_y, j_x, j_y, rho, charge) = solve_circuit(x_ax, y_ax, rho, deltaV=1)
    I = np.sum(j_x, axis=0) * dy

    return 1 / np.mean(I[1:])
