#!/usr/bin/env python3
""" A geodesy package for Python3 - distance calculation functions """

# import statements
import numpy as np
from numpy.linalg import norm

# DEBUG CODE - remove on dist
DEBUG = True  # set to False to clear up console output


def haversine(a, b):
    """ Returns the haversine distance between two GeoPoints.
        Inputs:
            a [GeoPoint] - first location
            b [GeoPoint] - second location
        Outputs:
            distance [float] - great circle distance, in meters
    """

    try:
        assert a.spheroid.datum == b.spheroid.datum
    except AssertionError as e:
        print("Reference spheroids must have same datum.")
        raise e
    else:
        a_lon = np.deg2rad(a.lon)
        a_lat = np.deg2rad(a.lat)
        b_lon = np.deg2rad(b.lon)
        b_lat = np.deg2rad(b.lat)
        distance = _haversine(a_lon, a_lat, b_lon, b_lat, a.spheroid)
    return distance


def euclidean(a, b):
    """ Returns the euclidean distance between two points.
        Inputs:
            a : [GeoPoint] - first location
            b : [GeoPoint] - second location
        Outputs:
            distance : [float] - euclidean distance, in meters
    """
    # check that the datums are the same
    try:
        assert a.spheroid.datum == b.spheroid.datum
    except AssertionError as e:
        print("Error: Inputs must have same reference spheroid.")
        raise e
    else:
        distance = _euclidean(a.x, a.y, a.z, b.x, b.y, b.z)
    return distance


# base case functions
def _haversine(lon1, lat1, lon2, lat2, ref_spheroid):
    """ Returns the haversine distance between two points based on mean earth radius.
        Inputs:
            lon1 [float] - longitude of first location, in radians
            lat1 [float] - latitude of first location, in radians
            lon2 [float] - longitude of second location, in radians
            lat2 [float] - latitude of second location, in radians
        Outputs:
            distance [float] - great circle distance, in meters
    """

    # get spheroid
    globe = ref_spheroid

    # calculate deltas
    dlon = lon2 - lon1
    dlat = lat2 - lat1

    # intermediate calculations
    a = np.square(np.sin(dlat / 2)) + np.cos(lat1) * np.cos(lat2) * np.square(np.sin(dlon / 2))
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))

    # calculate distance
    distance = globe.mean_r * c
    return distance


def _euclidean(x1, y1, z1, x2, y2, z2):
    """ Returns the euclidean distance between two points.
        Inputs:
            x1, y1, z1 - float, meters
            x2, y2, z2 - float, meters
        Outputs:
            distance - meters, float
    """
    # calculate deltas
    x_vec = np.array([x1, y1, z1])
    y_vec = np.array([x2, y2, z2])

    # calculate distance
    distance = norm(y_vec - x_vec)
    return distance

