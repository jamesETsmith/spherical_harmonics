import numpy as np
from scipy.special import sph_harm
from matplotlib import cm, colors
from mpl_toolkits.mplot3d import Axes3D


class ylm:
    """
    A class for spherical harmonic objects.

    Attributes:
    -----------

    unit_sphere : np.ndarray
        The cartesian coordinates for the unit sphere.

    m : int

    l : int

    fcolors :  np.ndarray

    """

    def __init__(self, l, m):
        """
        """
        self.l = l
        self.m = m

        self.unit_sphere, self.phi, self.theta = make_unit_sphere_()

        # Calculate the spherical harmonic Y(l,m) and normalize to [0,1]
        self.fcolors = sph_harm(m, l, self.theta, self.phi)
        print(type(self.fcolors))

    def plot(self, ax, part):
        """
        Plot ``ylm``.

        Part can be "re", "im", "abs"
        """

        fcs = None

        if part == "abs":
            fcs = np.absolute(self.fcolors)
        elif part == "re":
            fcs = self.fcolors.real
        elif part == "im":
            fcs = self.fcolors.imag
        else:
            raise ValueError("Please enter a valid type 'abs', 're', or 'im'")

        fmax, fmin = fcs.max(), fcs.min()
        fcs = (fcs - fmin) / (fmax - fmin)

        plot = ax.plot_surface(
            self.unit_sphere[0],
            self.unit_sphere[1],
            self.unit_sphere[2],
            rstride=1,
            cstride=1,
            facecolors=cm.seismic(fcs))

        return plot

    def propagate(self, t):
        self.fcolors *= np.exp(-1j * t * self.l * (self.l + 1))


def make_unit_sphere_():
    """
    Function to make the unit sphere coordinates.
    """

    phi = np.linspace(0, np.pi, 100)
    theta = np.linspace(0, 2 * np.pi, 100)
    phi, theta = np.meshgrid(phi, theta)
    phi, theta = phi, theta

    # The Cartesian coordinates of the unit sphere
    x = np.sin(phi) * np.cos(theta)
    y = np.sin(phi) * np.sin(theta)
    z = np.cos(phi)

    unit_sphere = np.array([x, y, z])
    return unit_sphere, phi, theta
