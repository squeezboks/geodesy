#!/usr/bin/env python3
""" A geodesy package for Python3 - utility functions """

# import statements
import numpy as np
from numpy.linalg import norm
from .convert_coordinates import ecef2ned, ned2enu

# DEBUG CODE - TODO: remove on dist
DEBUG = True  # set to False to clear up console output


def get_nvector(lon, lat):
    """ Returns the local normal vector in ECEF.
        Inputs:
            lon   : [float]  - longitude, rad
            lat   : [float]  - latitude, rad
        Outputs:
            x,y,z : [float] - unit length
    """

    x = np.cos(lat) * np.cos(lon)
    y = np.cos(lat) * np.sin(lon)
    z = np.sin(lat)

    return x, y, z


def elevation(a, b):
    """ Returns the elevation angle to a target orbital position from a terrestrial position.
        Inputs:
            a : [GeoPoint] - terrestrial position
            b : [GeoPoint] - orbital position
        Outputs:
            elev_angle : [float] - elevation angle, in degrees
    """
    # determine satellites position in NED
    n, e, d = ecef2ned(b.x, b.y, b.z, a)

    # determine position in ENU
    e, n, u = ned2enu(n, e, d)

    # calculate elevation angle
    elev_angle = np.arctan2(u, (np.sqrt(np.power(e, 2) + np.power(n, 2))))
    return np.rad2deg(elev_angle)


def scan_angle(a, b):
    """ Returns the scan angle off satellites boresight for a target location.
        Inputs:
            a : [GeoPoint] - satellites location
            b : [GeoPoint] - target location
        Outputs:
            scan_angle : [float] - scan angle, in degrees
    """
    # convert target to NED coordinates
    tx, ty, tz = ecef2ned(b.x, b.y, b.z, a)
    target = np.array([tx, ty, tz])

    # in NED, boresight is z basis vector
    bore = np.array([0, 0, 1])

    # calculate direction cosine between the vectors
    gamma = np.dot(target, bore)/norm(target)

    # calculate angle
    angle = np.arccos(gamma)

    return np.rad2deg(angle)

