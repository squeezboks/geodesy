#!/usr/bin/env python3
""" A geodesy package for Python3 - entity definitions """

# import statements
import numpy as np

# DEBUG CODE - remove on dist
DEBUG = False  # set to False to clear up console output


# -- SOME USEFUL CLASSES -- #
class RefEllipsoid(object):
    # instance init
    def __init__(self, datum='WGS84'):
        self.datum = datum

        if self.datum == 'WGS84':
            # defined parameters
            self.a = 6378.137e3  # meters, semi-major axis
            self.f_inv = 298.257223563  # flattening

            # derived parameters
            self.f = 1 / self.f_inv
            self.b = (1 - self.f) * self.a  # meters, semi-minor axis
            self.a2 = np.power(self.a, 2)
            self.b2 = np.power(self.b, 2)
            self.e2 = (self.a2 - self.b2) / self.a2
            self.e = np.sqrt(self.e2)  # first eccentricity
            self.mean_r = (2 * self.a + self.b) / 3  # mean radius to minimize square relative error

    @property
    def radius_x(self):
        return self.a

    @property
    def radius_y(self):
        return self.a

    @property
    def radius_z(self):
        return self.b

    def curvature_meridian(self, lat):
        w = np.sqrt(1 - self.e2 * np.power(np.sin(lat), 2))
        result = self.a * (1 - self.e2) / np.power(w, 3)
        return result

    def curvature_prime_vertical(self, lat):
        w = np.sqrt(1 - self.e2 * np.power(np.sin(lat), 2))
        result = self.a / w
        return result