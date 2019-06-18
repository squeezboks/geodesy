#!/usr/bin/env python3
""" A geodesy package for Python3 - conversion functions """

# import statements
import numpy as np
from numpy.linalg import norm

# DEBUG CODE - remove on dist
DEBUG = False  # set to False to clear up console output


def geodetic2ecef(lon, lat, alt, spheroid, units='rad'):
    """ Returns the ECEF position vector (x,y,z) corresponding to the provided geodetic position vector (lon,lat,alt).
        Inputs:
            lon   : [float]  - longitude, in degrees or radian
            lat   : [float]  - geodetic latitude, in degrees or radians
            alt   : [float]  - altitude (or height) above the reference datum, in meters
            units : [string] - if 'deg' will convert to radians
        Outputs:
            x : [float] - distance, in meters, along the x-axis of the ECEF coordinate system
            y : [float] - distance, in meters, along the x-axis of the ECEF coordinate system
            z : [float] - distance, in meters, along the x-axis of the ECEF coordinate system
    """
    # convert to radians
    if units == 'deg':
        lat = np.deg2rad(lat)
        lon = np.deg2rad(lon)

    # grab some spheroid params
    b = spheroid.b
    e2 = spheroid.e2

    # intermediate calculations
    n = spheroid.curvature_prime_vertical(lat)
    p = (n + alt) * np.cos(lat)

    if lat == 90:
        x = 0
        y = 0
        z = b + alt
    elif lat == -90:
        x = 0
        y = 0
        z = -(b + alt)
    else:
        x = p * np.cos(lon)
        y = p * np.sin(lon)
        z = (n * (1 - e2) + alt) * np.sin(lat)

    return x, y, z


def ecef2geodetic(x, y, z, spheroid, units='rad'):
    """ Returns the geodetic position vector (lon,lat,alt) corresponding to the provided ECEF position vector (x,y,z).
        Inputs:
            x : [float] - distance, in meters, along the x-axis of the ECEF coordinate system
            y : [float] - distance, in meters, along the x-axis of the ECEF coordinate system
            z : [float] - distance, in meters, along the x-axis of the ECEF coordinate system
        Outputs:
            lon : [float] - longitude, in degrees
            lat : [float] - geodetic latitude, in degrees
            alt : [float] - altitude (or height) above the reference datum, in meters
    """
    a = spheroid.a
    b = spheroid.b
    e2 = spheroid.e2
    ep = np.sqrt(1 - e2)

    # intermediate calculations
    p = np.hypot(x, y)

    if p > 0:
        # 'Fast transform from geocentric to geodetic coordinates'
        # T. Fukushima
        # Journal of Geodesy, 1999

        # gather parameters
        zp = ep * z
        c = a * e2
        u = 2 * (zp - c)
        v = 2 * (zp + c)

        # check Tm value
        tm = (c - zp) / p
        tm3 = np.power(tm, 3)
        tm4 = np.power(tm, 4)
        fm = p * tm4 + u * tm3 + v * tm - p

        # calculate bounds
        t0 = p / (zp + c)
        t1 = (p - c + zp) / (p - c + 2 * zp)

        # initial guess
        if tm <= 0:  # case 1
            if DEBUG: print(" - Case 1: Tm <= 0 (Tm = {:g})".format(tm))
            t = t1
        elif tm >= 1:  # case 2
            if DEBUG: print(" - Case 2: Tm >= 1 (Tm = {:g})".format(tm))
            t = t0
        elif fm >= 0:  # case 3a
            if DEBUG: print(" - Case 3a:  0 < Tm < 1 & Fm >= 0 (Tm = {:g}, Fm = {:g})".format(tm, fm))
            t = t0
        else:  # case 3b
            if DEBUG: print(" - Case 3a:  0 < Tm < 1 & Fm < 0 (Tm = {:g}, Fm = {:g})".format(tm, fm))
            t = t1

        # compute initial powers once
        t2 = np.power(t, 2)
        t3 = np.power(t, 3)
        t4 = np.power(t, 4)

        # set limits
        error_limit = 1e-9  # arbitrary - can be less or more
        loop_limit = 10

        # iterate
        count = 0
        exit_flag = False
        while not exit_flag:
            # calculate delta_t
            delta_t = (p - (p * t4 + u * t3 + v * t)) / (4 * p * t3 + 3 * u * t2 + v)
            # increment t
            t = t + delta_t
            # compute powers once
            t2 = np.power(t, 2)
            t3 = np.power(t, 3)
            t4 = np.power(t, 4)

            # check for exit conditions
            if abs(delta_t) < error_limit:
                if DEBUG: print("Exiting on error limit")
                exit_flag = True
            elif count >= loop_limit:
                print("Exiting on loop limit")
                exit_flag = True

            # increment counter
            count = count + 1

            # DEBUG CODE - remove on dist
            if DEBUG:
                print("loop_count = " + str(count))
                print("delta_t: " + str(delta_t))

        # solve longitude
        lon = np.arctan2(y, x)

        # solve latitude
        lat = np.arctan2((1 - t2), (2 * ep * t))

        # solve altitude
        alt = (2 * p * ep * t + z * (1 - t2) - a * ep * (1 + t2)) / np.sqrt(np.power((1 + t2), 2) - 4 * e2 * t2)

    else:
        if z >= 0:
            lat = np.pi
        else:
            lat = -1 * np.pi
        lon = 0
        alt = abs(z) - b

    # convert to radians
    if units == 'deg':
        lon = np.rad2deg(lon)
        lat = np.rad2deg(lat)

    return lon, lat, alt


def ecef2ned(x, y, z, origin):
    """ Returns the NED position vector corresponding to the provided ECEF position vector (x,y,z) and origin.
        Inputs:
            x : [float] - distance, in meters, along the x-axis of the ECEF coordinate system
            y : [float] - distance, in meters, along the x-axis of the ECEF coordinate system
            z : [float] - distance, in meters, along the x-axis of the ECEF coordinate system
            origin [GeoPoint] - instance of a GeoPoint object describing the NED origin
        Outputs:
            x : [float] - distance, in meters, along the x-axis of the local NED coordinate system
            y : [float] - distance, in meters, along the x-axis of the local NED coordinate system
            z : [float] - distance, in meters, along the x-axis of the local NED coordinate system
    """
    p_ecef = np.matrix([[x], [y], [z]])  # point in ECEF
    p0_ecef = np.matrix([[origin.x], [origin.y], [origin.z]])  # origin in ECEF

    # grab origin lon/lat in radians
    lon = np.deg2rad(origin.lon)
    lat = np.deg2rad(origin.lat)

    # construct rotation matrix ECEF -> NED
    rotation_matrix = np.matrix([[-1 * np.sin(lat) * np.cos(lon), -1 * np.sin(lat) * np.sin(lon), np.cos(lat)],
                                 [-1 * np.sin(lon), np.cos(lon), 0.],
                                 [-1 * np.cos(lat) * np.cos(lon), -1 * np.cos(lat) * np.sin(lon), -np.sin(lat)]])

    # solve with matrix multiplication
    p_ned = rotation_matrix * (p_ecef - p0_ecef)  # point in NED

    # assign results from vector
    x = p_ned.item(0)
    y = p_ned.item(1)
    z = p_ned.item(2)

    return x, y, z


def ned2ecef(x, y, z, origin):
    """ Returns the ECEF position vector corresponding to the provided NED position vector (x,y,z) and origin.
        Inputs:
            x      : [float]    - distance, in meters, along the x-axis of the local NED coordinate system
            y      : [float]    - distance, in meters, along the y-axis of the local NED coordinate system
            z      : [float]    - distance, in meters, along the z-axis of the local NED coordinate system
            origin : [GeoPoint] - instance of a GeoPoint object describing the NED origin
        Outputs:
            x : [float] - distance, in meters, along the x-axis of the ECEF coordinate system
            y : [float] - distance, in meters, along the y-axis of the ECEF coordinate system
            z : [float] - distance, in meters, along the z-axis of the ECEF coordinate system
    """
    v_ned = np.matrix([[x], [y], [z]])  # point in NED
    p0_ecef = np.matrix([[origin.x], [origin.y], [origin.z]])  # origin in ECEF

    # grab origin lon/lat in radians
    lon = np.deg2rad(origin.lon)
    lat = np.deg2rad(origin.lat)

    # construct rotation matrix NED -> ECEF
    rotation_matrix = np.matrix([[-1 * np.sin(lat) * np.cos(lon), -1 * np.sin(lat) * np.sin(lon), np.cos(lat)],
                                 [-1 * np.sin(lon), np.cos(lon), 0.],
                                 [-1 * np.cos(lat) * np.cos(lon), -1 * np.cos(lat) * np.sin(lon), -np.sin(lat)]])
    rotation_matrix = rotation_matrix.transpose()

    # solve with matrix multiplication
    p_ecef = rotation_matrix * v_ned + p0_ecef  # point in ECEF

    # assign results from vector
    x = p_ecef.item(0)
    y = p_ecef.item(1)
    z = p_ecef.item(2)

    return x, y, z


def ned2enu(x, y, z):
    """ Returns the ENU position vector corresponding to a provided NED position vector (x,y,z).
        Inputs:
            x : [float] - distance, in meters, along the x-axis of the local NED coordinate system
            y : [float] - distance, in meters, along the y-axis of the local NED coordinate system
            z : [float] - distance, in meters, along the z-axis of the local NED coordinate system
        Outputs:
            x : [float] - distance, in meters, along the x-axis of the local ENU coordinate system
            y : [float] - distance, in meters, along the x-axis of the local ENU coordinate system
            z : [float] - distance, in meters, along the x-axis of the local ENU coordinate system
    """

    v_ned = np.matrix([[x], [y], [z]])  # point in NED

    # construct rotation matrix NED -> ENU
    rotation_matrix = np.matrix([[0., 1., 0.],
                                 [1., 0., 0.],
                                 [0., 0.,-1.]])

    # solve with matrix multiplication
    p_enu = rotation_matrix * v_ned  # point in ENU

    # assign results from vector
    x = p_enu.item(0)
    y = p_enu.item(1)
    z = p_enu.item(2)

    return x, y, z


def enu2ned(x, y, z):
    """ Returns the NED position vector corresponding to a provided ENU position vector (x,y,z).
        Inputs:
            x : [float] - distance, in meters, along the x-axis of the local ENU coordinate system
            y : [float] - distance, in meters, along the y-axis of the local ENU coordinate system
            z : [float] - distance, in meters, along the z-axis of the local ENU coordinate system
        Outputs:
            x : [float] - distance, in meters, along the x-axis of the local NED coordinate system
            y : [float] - distance, in meters, along the x-axis of the local NED coordinate system
            z : [float] - distance, in meters, along the x-axis of the local NED coordinate system
    """

    v_enu = np.matrix([[x], [y], [z]])  # point in ENU

    # construct rotation matrix NED -> ENU
    rotation_matrix = np.matrix([[0., 1., 0.],
                                 [1., 0., 0.],
                                 [0., 0., -1.]])

    # solve with matrix multiplication
    p_ned = rotation_matrix * v_enu  # point in ENU

    # assign results from vector
    x = p_ned.item(0)
    y = p_ned.item(1)
    z = p_ned.item(2)

    return x, y, z


def aer2xyz(az, el, r, units='deg'):

    if units == 'deg':
        az = np.deg2rad(az)
        el = np.deg2rad(el)

    x = r * np.cos(el) * np.cos(az)
    y = r * np.cos(el) * np.sin(az)
    z = r * np.sin(el)

    return x, y, z


def xyz2aer(x, y, z):

    vec3 = np.matrix([[x], [y], [z]])
    vec2 = np.matrix([[x], [y]])

    r = norm(vec3)
    az = np.arctan2(y, x)
    el = np.arctan2(z, norm(vec2))

    return az, el, r
