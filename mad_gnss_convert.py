"""
Example command

python mad_gnss_convert.py ~/Downloads/los_20180506.001.h5.hdf5 2018-05-06T05:20:05 600 45 60 -130 -100
"""

from __future__ import annotations

from pathlib import Path
from datetime import datetime, timedelta
import argparse
import dateutil.parser

import h5py
import numpy as np


def mad_times(file: Path):
    """
    Load all times in a Madrigal GNSS TEC file
    Knowing the time is necesary to index the unordered data in the file
    """

    with h5py.File(file, "r") as f:
        return f["Data"]["Table Layout"]["ut1_unix"]


def mad_coord(file: Path) -> tuple:
    """
    Load lat, lon grid in a Madrigal GNSS TEC file
    These are used to index unordered data in the file
    """

    with h5py.File(file, "r") as f:
        gdlat = f["Data"]["Table Layout"]["gdlat"]
        glon = f["Data"]["Table Layout"]["glon"]

    return gdlat, glon


def get_time_index(times, time: datetime, tolerance: timedelta):
    """
    return boolean array corresponding to the single closest time index
    most closely matching the requested time within tolerance
    """

    i = abs(times - time.timestamp()).argmin()
    if abs(times[i] - time.timestamp()) > tolerance.total_seconds():
        raise ValueError(f"Time {time} not found in file within tolerance {tolerance}")

    return times == times[i]


def time_range(times, start: datetime, end: datetime):
    """
    return boolean array corresponding to the time range requested
    inclusive of endpoints.
    """

    i = (times >= start.timestamp()) & (times <= end.timestamp())
    if i.sum() == 0:
        raise ValueError(f"No data found for time range {start} to {end}")

    return i


def latlon_range(
    gdlat, glon, lat_bounds: tuple[float, float], lon_bounds: tuple[float, float]
):
    """
    return boolean array corresponding to the lat lon range requested
    inclusive of endpoints.
    """
    ilat = (gdlat >= lat_bounds[0]) & (gdlat <= lat_bounds[1])
    if ilat.sum() == 0:
        print(f"gdlat: {gdlat} min: {gdlat.min()}  max: {gdlat.max()}")
        raise ValueError(f"No observations found in the requested latitude range.")

    ilon = (glon >= lon_bounds[0]) & (glon <= lon_bounds[1])

    if ilon.sum() == 0:
        print(f"glon: {glon}  min: {glon.min()}  max: {glon.max()}")
        raise ValueError(f"No observations found in the requested longitude range.")

    return ilat & ilon


def get_index(file: Path, start: datetime, end: datetime):
    """
    return boolean vector for file of desired time and lat lon range
    """
    print(f"loading times from {file}, this may take a couple minutes...")
    time = mad_times(file)
    index = time_range(time, start, end)
    time = time[index]
    with h5py.File(output, "w") as f:
        f["/time"] = time
    # %% Filter by lat lon
    print(f"loading lat lon from {file}, this may take a couple minutes...")
    gdlat, glon = mad_coord(file)
    index &= latlon_range(gdlat, glon, P.lat_bounds, P.lon_bounds)
    gdlat = gdlat[index]
    glon = glon[index]
    with h5py.File(output, "a") as f:
        f["/gdlat"] = gdlat
        f["/glon"] = glon
        f["/index"] = index

    return index


if __name__ == "__main__":
    """
    if __name__ == "__main__": allows functions to be imported from this file,
    or for the file to be run as a script itself.
    """

    p = argparse.ArgumentParser()
    p.add_argument("file", help="Madrigal GNSS TEC file to convert")
    p.add_argument("time", help="time request (UTC)")
    p.add_argument("atol", help="time tolerance (seconds)", type=float)
    p.add_argument("lat_bounds", help="latitude bounds (degrees)", nargs=2, type=float)
    p.add_argument("lon_bounds", help="longitude bounds (degrees)", nargs=2, type=float)
    p.add_argument("-o", "--output", help="output file name")
    P = p.parse_args()

    file = Path(P.file).expanduser().resolve(strict=True)

    if not P.output:
        output = file.parent / (file.stem + "_slice.h5")
    print(f"output file: {output}")

    treq = dateutil.parser.parse(P.time)
    # treq = datetime(2018, 5, 6, 5, 20, 5)
    print("Requested time:", treq, "tolerance:", P.atol, "seconds")

    atol = timedelta(seconds=P.atol)
    # %% Filter by Time
    start_time = treq - atol
    end_time = treq + atol

    index = None
    if output.is_file():
        with h5py.File(output, "r") as f:
            if "index" in f:
                index = f["/index"][:]

    if index is None:
        index = get_index(file, start_time, end_time)

    print(f"Using {index.sum()} values from {file}")

    with h5py.File(file, "r") as f:
        working = f["Data"]["Table Layout"][index]

    unique_times = np.unique(working["ut1_unix"])
    # implicitly sorted by numpy.unique
    print(f"{unique_times.size} unique times found from {start_time} to {end_time}.")
