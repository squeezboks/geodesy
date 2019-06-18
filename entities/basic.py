#!/usr/bin/env python3
""" A geodesy package for Python3 - entity definitions """

# import statements
import numpy as np
from ..entities.ellipsoids import RefEllipsoid
from ..functions import convert_coordinates

class GeoPoint(object):
    def __init__(self, lon=0, lat=0, alt=0, xyz=None, lla=None, name=''):
        """ enter lat/lon in deg, alt in m or xyx in m """

        # spheroid
        self.spheroid = RefEllipsoid(datum='WGS84')

        # initialize ecef
        self._x = 0
        self._y = 0
        self._z = 0

        # set geodetic position
        self._lon = np.deg2rad(lon) # rad
        self._lat = np.deg2rad(lat) # rad
        self._alt = alt # m

        if xyz is not None:
            self.set_ecef(xyz[0], xyz[1], xyz[2])
            self._update_geodetic()
        elif lla is not None:
            self.set_geodetic(lla[0], lla[1], lla[2])
            self._update_ecef()
        else:
            self._update_ecef()

        self.name = name

    @property
    def x(self):
        return self._x

    @x.setter
    def x(self, value):
        self._x = value
        self._update_geodetic()

    @property
    def y(self):
        return self._y

    @y.setter
    def y(self, value):
        self._y = value
        self._update_geodetic()

    @property
    def z(self):
        return self._z

    @z.setter
    def z(self, value):
        self._z = value
        self._update_geodetic()

    @property
    def xyz(self):
        return [self.x, self.y, self.z]

    @property
    def lon(self):
        return np.rad2deg(self._lon)

    @lon.setter
    def lon(self, value):
        self._lon = np.deg2rad(value)
        self._update_ecef()

    @property
    def lat(self):
        return np.rad2deg(self._lat)

    @lat.setter
    def lat(self, value):
        self._lat = np.deg2rad(value)
        self._update_ecef()

    @property
    def alt(self):
        return self._alt

    @alt.setter
    def alt(self, value):
        self._alt = value
        self._update_ecef()

    @property
    def lla(self):
        return self.lon, self.lat, self.alt

    def set_ecef(self, x, y, z):
        self._x = x
        self._y = y
        self._z = z
        self._update_geodetic()

    def set_geodetic(self, lon, lat, alt):
        self._lon = np.deg2rad(lon)
        self._lat = np.deg2rad(lat)
        self._alt = alt
        self._update_ecef()

    def _update_ecef(self):
        x, y, z = convert_coordinates.geodetic2ecef(self._lon, self._lat, self._alt, self.spheroid)

        self._x = x
        self._y = y
        self._z = z

    def _update_geodetic(self):
        lon, lat, alt = convert_coordinates.ecef2geodetic(self._x, self._y, self._z, self.spheroid)

        self._lon = lon
        self._lat = lat
        self._alt = alt
