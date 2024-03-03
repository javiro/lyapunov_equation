import numpy as np
import matplotlib.pyplot as plt
# import RungeKutta as rk

from copy import deepcopy
from scipy.integrate import odeint
from scipy.linalg import expm
from abc import ABC, abstractmethod


class DynamicalSystem(ABC):

    def __init__(self, ss_dimension: int, t_initial: float, num_t: int):
        self.ss_dimension = ss_dimension
        self.t_initial = t_initial
        self.num_t = num_t
    
    @property
    @abstractmethod
    def _name(self):
        """Implements the name of the class."""

    @abstractmethod
    def velocity(self, ssp, t):
        """Implements the velocity."""

    @abstractmethod
    def stability_matrix(self, ssp):
        """Implements the stability matrix."""

    def flow(self, ssp0, deltat):
        ssp_solution = odeint(self.velocity, ssp0, [0.0, deltat])
        ssp_deltat = ssp_solution[-1, :]
        return ssp_deltat

    def jacobian_velocity(self, ssp_jacobian, t):
        delta = self.D * np.identity(self.ss_dimension)
        ssp = ssp_jacobian[0:self.ss_dimension]
        jacobian = ssp_jacobian[self.ss_dimension:].reshape((self.ss_dimension, self.ss_dimension))
        vel_j = np.zeros(np.size(ssp_jacobian))
        vel_j[0:self.ss_dimension] = self.velocity(ssp, t)
        vel_tangent = self.stability_matrix(ssp) @ jacobian + jacobian @ np.transpose(self.stability_matrix(ssp)) + delta
        vel_j[self.ss_dimension:] = np.reshape(vel_tangent, self.ss_dimension ** 2)
        return vel_j

    def jacobian(self, ssp, t):
        period = t
        jacobian0 = np.identity(self.ss_dimension)
        ssp_jacobian0 = np.zeros(self.ss_dimension + self.ss_dimension ** 2)
        ssp_jacobian0[0:self.ss_dimension] = ssp
        ssp_jacobian0[self.ss_dimension:] = np.reshape(jacobian0, self.ss_dimension ** 2)
        t_final = period
        t_array = np.linspace(self.t_initial, t_final, self.num_t)
        ssp_jacobian_solution = odeint(self.jacobian_velocity, ssp_jacobian0, t_array)
        jacobian = ssp_jacobian_solution[-self.ss_dimension, self.ss_dimension:] \
            .reshape((self.ss_dimension, self.ss_dimension))
        return jacobian
        # return expm(self.stability_matrix(ssp) * t)
    
    def q_velocity(self, ssp_q, t):
        delta = self.D * np.identity(self.ss_dimension)
        ssp = ssp_q[0:self.ss_dimension]
        q = ssp_q[self.ss_dimension:].reshape((self.ss_dimension, self.ss_dimension))
        vel_j = np.zeros(np.size(ssp_q))
        vel_j[0:self.ss_dimension] = self.velocity(ssp, t)
        vel_tangent = self.stability_matrix(ssp) @ q + q @ np.transpose(self.stability_matrix(ssp)) + delta
        vel_j[self.ss_dimension:] = np.reshape(vel_tangent, self.ss_dimension ** 2)
        return vel_j

    def q(self, ssp, t):
        period = t
        q0 = np.identity(self.ss_dimension)
        ssp_q0 = np.zeros(self.ss_dimension + self.ss_dimension ** 2)
        ssp_q0[0:self.ss_dimension] = ssp
        ssp_q0[self.ss_dimension:] = np.reshape(q0, self.ss_dimension ** 2)
        t_final = period
        t_array = np.linspace(self.t_initial, t_final, self.num_t)
        ssp_q_solution = odeint(self.q_velocity, ssp_q0, t_array)
        q = ssp_q_solution[-1, self.ss_dimension:] \
            .reshape((self.ss_dimension, self.ss_dimension))
        return q
    
    def get_eigenvalues_of_q(self, ssp_solution, t_array):
        return np.array([np.linalg.eigvals(np.linalg.inv(self.q(ssp_solution[i], t))) for i, t in enumerate(t_array)])
    
    def get_eigenvectors_of_q(self, ssp_solution, t_array):
        output = [np.linalg.eig(np.linalg.inv(self.q(ssp_solution[i], t))) for i, t in enumerate(t_array)]
        return np.array([e[1].T[0] for e in output]), [e[1].T[1] for e in output]
    
    def get_norm_of_q(self, ssp_solution, t_array):
        output = [np.linalg.norm(self.q(ssp_solution[i], t)) for i, t in enumerate(t_array)]
        return np.array(output)

    # def q(self, t, d, t_array, ssp_solution):
    #     q0 = np.identity(self.ss_dimension)
    #     delta = d * np.identity(self.ss_dimension)
    #     index_t_array = np.where(t_array < t)[0]
    #     # TODO: clarify what is delta. Can q0 be the identity matrix?
    #     # q0 must to by symmetric
    #     integral = sum([self.jacobian(ssp_solution[i], t) @ delta @ np.transpose(self.jacobian(ssp_solution[i], t))
    #                     for i in index_t_array])
    #     return self.jacobian(ssp_solution[0], t) @ q0 @ np.transpose(self.jacobian(ssp_solution[0], t)) + integral


class DSNoisyCircle(DynamicalSystem):

    LAMBDA = 1.0
    R_C = 1.0
    OMEGA = 1.0
    D = 4.0
    # TODO: clarify how to implement the noise
    EPSILON = 0.0001

    def __init__(self, ss_dimension: int, t_initial: float, num_t: int):
        self.ss_dimension = ss_dimension
        self.t_initial = t_initial
        self.num_t = num_t
        super().__init__(self.ss_dimension, self.t_initial, self.num_t)
    
    @property
    def _name(self):
        """Implements the name of the class."""
        return "DSNoisyCircle"

    def velocity(self, ssp, t):
        """"Implements the velocity."""
        x, y = ssp
        epx = self.EPSILON
        epy = self.EPSILON
        # print(epx, epy)
        dxdt = self.LAMBDA * (self.R_C - np.sqrt(x ** 2 + y ** 2)) * x - self.OMEGA * y + \
            np.sqrt(2 * self.D) * epx
        dydt = self.LAMBDA * (self.R_C - np.sqrt(x ** 2 + y ** 2)) * y + self.OMEGA * x + \
            np.sqrt(2 * self.D) * epy
        vel = np.array([dxdt, dydt], float)
        return vel

    def stability_matrix(self, ssp):
        """"Implements the stability matrix."""
        x, y = ssp
        return np.array([[
            self.LAMBDA * (self.R_C - np.sqrt(x ** 2 + y ** 2) - x ** 2 / np.sqrt(x ** 2 + y ** 2)),
            -self.LAMBDA * x * y / np.sqrt(x ** 2 + y ** 2) - self.OMEGA],
            [-self.LAMBDA * x * y / np.sqrt(x ** 2 + y ** 2) + self.OMEGA,
            self.LAMBDA * (self.R_C - np.sqrt(x ** 2 + y ** 2) - y ** 2 / np.sqrt(x ** 2 + y ** 2))]], float)


class VdPOscillator(DynamicalSystem):

    LAMBDA = 3.0
    D = 0.1
    EPSILON = 0.01

    def __init__(self, ss_dimension: int, t_initial: float, num_t: int):
        self.ss_dimension = ss_dimension
        self.t_initial = t_initial
        self.num_t = num_t
        super().__init__(self.ss_dimension, self.t_initial, self.num_t)
    
    @property
    def _name(self):
        """Implements the name of the class."""
        return "VdPOscillator"

    def velocity(self, ssp, t):
        """"Implements the velocity."""
        x, y = ssp
        dxdt = self.LAMBDA * (x - 1/3 * x ** 3 - y) + np.sqrt(2 * self.D) * self.EPSILON
        dydt = 1 / self.LAMBDA * x + np.sqrt(2 * self.D) * self.EPSILON
        vel = np.array([dxdt, dydt], float)
        return vel

    def stability_matrix(self, ssp):
        """"Implements the stability matrix."""
        x, y = ssp
        return np.array([[
            self.LAMBDA * (1 - x ** 2),
            -self.LAMBDA],
            [1 / self.LAMBDA, 0]], float)


if __name__ == "__main__":

    tInitial = 0.0  # Initial time
    tFinal = 10.0  # Final time
    Nt = 100  # Number of time points to be used in the integration

    tArray = np.linspace(tInitial, tFinal, Nt)  # Time array for solution
    # ssp0 = np.array([1.0,
    #                  1.0,
    #                  1.0], float)  # Initial condition for the solution

    noisy_system = DSNoisyCircle(ss_dimension=2, t_initial=tInitial, num_t=Nt)
    # noisy_system = VdPOscillator(ss_dimension=2, t_initial=tInitial, num_t=Nt)
    ssp0 = np.array([1.2, 1.2])

    # sspSolution = rk.RK4(noisy_system.velocity, ssp0, tArray)
    sspSolution = odeint(noisy_system.velocity, ssp0, tArray)
    # sspSolution_aux = deepcopy(sspSolution)
    eigenvalues = noisy_system.get_eigenvalues_of_q(sspSolution, tArray)
    eigenvectors = noisy_system.get_eigenvectors_of_q(sspSolution, tArray)
    q_norm = noisy_system.get_norm_of_q(sspSolution, tArray)
    # eigenvalues, eigenvectors = noisy_system.get_eigenvectors_of_q(sspSolution, tArray)

    xt = sspSolution[:, 0]  # Read x(t)
    yt = sspSolution[:, 1]  # Read y(t)

    # print(sum([np.array([[1.2, 1.2], [1.2, 1.2]]), np.array([[1.2, 1.2], [1.2, 1.2]])]))
    # print([noisy_system.q(t, tArray, sspSolution) for t in tArray])  # Print final point
    # print(eigenvalues)
    # print(sspSolution)

    sigma = [np.sqrt(1/(2*e[0])) for e in eigenvalues]
    sigma1 = [np.sqrt(1/(2*e[1])) for e in eigenvalues]
    e0 = [e[0] for e in eigenvalues]
    e1 = [e[1] for e in eigenvalues]

    u0 = [x[0] for x in eigenvectors[0]]
    v0 = [x[1] for x in eigenvectors[0]]
    u1 = [x[0] for x in eigenvectors[1]]
    v1 = [x[1] for x in eigenvectors[1]]

    fig, axs = plt.subplots(3, 2, figsize=(15, 15))  # Create a figure instance
    # ax = fig.gca()  # Get current axes in 3D projection
    axs[0,0].set_title("Solution of the numerically integrated, noiseless")
    axs[0][0].plot(xt, yt, '.')  # Plot the solution

    axs[0][1].set_title("Eigenvectors")
    axs[0][1].quiver(xt, yt, u0, v0, color='red')
    axs[0][1].quiver(xt, yt, u1, v1, color='blue')

    axs[1][0].set_title("sigma 0")
    axs[1][0].plot(sigma)
    axs[1][1].set_title("sigma 1")
    axs[1][1].plot(sigma1)
    # axs[1][0].plot(e1, '.')
    axs[2][0].set_title("Eigenvalues mu")
    axs[2][0].plot(eigenvalues, '.')
    axs[2][1].set_title("Q norm")
    axs[2][1].plot(q_norm)
    # plt.plot(eigenvalues, '.', c='b')
    # ax.set_xlabel('x')  # Set x label
    # ax.set_ylabel('y')  # Set y label
    plt.savefig(noisy_system._name, bbox_inches='tight')  # Show the figure
