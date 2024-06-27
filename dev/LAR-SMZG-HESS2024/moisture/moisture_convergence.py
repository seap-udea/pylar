"""Computes the moisture convergence for a basin using ERA5 data."""

# Copyright: 2020, 2021, 2024 Ruben Molina <ruben.molina@udea.edu.co>
# License: BSD-2-clause

import collections
import os.path

import cdsapi
import geopandas
import matplotlib.pyplot as plt
import numpy as np
from cartopy.crs import PlateCarree
from cartopy.mpl.geoaxes import GeoAxes
from geographiclib.geodesic import Constants, Geodesic
from matplotlib.path import Path
from mpl_toolkits.axes_grid1 import AxesGrid
from netCDF4 import MFDataset  # pylint: disable=no-name-in-module # type: ignore

WGS84 = Geodesic(Constants.WGS84_a, Constants.WGS84_f)

NetcdfVariable = collections.namedtuple("NetcdfVariable", ["filename", "variable"])

CellGeometry = collections.namedtuple(
    "CellGeometry",
    [
        "eastern_side_length",
        "western_side_length",
        "northern_side_length",
        "southern_side_length",
        "cell_area",
    ],
)

BasinMask = collections.namedtuple(
    "BasinMask",
    ["eastward_in", "eastward_out", "northward_in", "northward_out", "basin"],
)

Coordinates = collections.namedtuple("Coordinates", ["latitude", "longitude"])


LONGITUDE_DELTA = 0.25  # cell size from ERA5 = 0.25 deg
LATITUDE_DELTA = -0.25  # cell size from ERA5 = 0.25 deg


GRDC_SHP = "shapefiles/%s.shp"

# Vertical integral of eastward water vapour flux
VIWVE = NetcdfVariable(
    "data/era5-vertical_integral_of_eastward_water_vapour_flux-*.nc", "p71.162"
)

# Vertical integral of northward water vapour flux
VIWVN = NetcdfVariable(
    "data/era5-vertical_integral_of_northward_water_vapour_flux-*.nc", "p72.162"
)

# Vertical integral of eastward cloud liquid water flux
VILWE = NetcdfVariable(
    "data/era5-vertical_integral_of_eastward_cloud_liquid_water_flux-*.nc", "p88.162"
)

# Vertical integral of northward cloud liquid water flux
VILWN = NetcdfVariable(
    "data/era5-vertical_integral_of_northward_cloud_liquid_water_flux-*.nc", "p89.162"
)

# Vertical integral of eastward cloud frozen water flux
VIIWE = NetcdfVariable(
    "data/era5-vertical_integral_of_eastward_cloud_frozen_water_flux-*.nc", "p90.162"
)

# Vertical integral of northward cloud frozen water flux
VIIWN = NetcdfVariable(
    "data/era5-vertical_integral_of_northward_cloud_frozen_water_flux-*.nc", "p91.162"
)


# 2 metre temperature
TEMP = NetcdfVariable("data/era5-2m_temperature-*.nc", "t2m")

# Evaporation
EVAP = NetcdfVariable("data/era5-evaporation-*.nc", "e")

# Total precipitation
PREC = NetcdfVariable("data/era5-total_precipitation-*.nc", "tp")


def wrap_lon360(lon):
    """Source: https://github.com/pyoceans/python-oceans/blob/master/oceans/ocfis/ocfis.py"""
    lon = np.atleast_1d(lon).copy()
    positive = lon > 0
    lon = lon % 360
    lon[np.logical_and(lon == 0, positive)] = 360
    return lon


def wrap_lon180(lon):
    """Source: https://github.com/pyoceans/python-oceans/blob/master/oceans/ocfis/ocfis.py"""
    lon = np.atleast_1d(lon).copy()
    angles = np.logical_or((lon < -180), (lon > 180))
    lon[angles] = wrap_lon360(lon[angles] + 180) - 180
    return lon


def era5_download(arg_variable, arg_year):
    """
    Retrieves monthly-averaged reanalysis data.

    Parameters
    ----------
    variable : str
        Variable name (e.g., '2m_temperature', 'total_precipitation', ...).
    year : int
        Year.

    Returns
    -------
    None.

    """

    client = cdsapi.Client()

    client.retrieve(
        "reanalysis-era5-single-levels-monthly-means",
        {
            "format": "netcdf",
            "product_type": "monthly_averaged_reanalysis",
            "variable": arg_variable,
            "year": arg_year,
            "month": [f"{(month + 1):02}" for month in range(12)],
            "time": "00:00",
        },
        f"data/era5-{arg_variable}-{arg_year}.nc",
    )


def download_all():
    """Download all the required files."""
    for year in range(1979, 2021):
        for variable in [
            "2m_temperature",
            "evaporation",
            "total_precipitation",
            "vertical_integral_of_eastward_cloud_frozen_water_flux",
            "vertical_integral_of_eastward_cloud_liquid_water_flux",
            "vertical_integral_of_eastward_water_vapour_flux",
            "vertical_integral_of_northward_cloud_frozen_water_flux",
            "vertical_integral_of_northward_cloud_liquid_water_flux",
            "vertical_integral_of_northward_water_vapour_flux",
        ]:
            if not os.path.exists(f"data/era5-{variable}-{year}.nc"):
                era5_download(variable, year)


def shapefile_reader(basin_name, station_name):
    """Read the shapefile."""
    shapefile = geopandas.read_file(GRDC_SHP % basin_name.split(" ")[0].lower())
    selection = np.isin(shapefile["station"], [station_name.upper()])
    selected = shapefile[selection]
    multipoly = selected.geometry[0]
    return list(multipoly.geoms)[0]  # first polygon in the multipolygon


def get_box_coordinates(basin):
    """Define the box coordinates."""
    # (bounding box + 1 degree buffer)

    north = np.ceil(max(basin.exterior.coords.xy[1])) + 1
    south = np.floor(min(basin.exterior.coords.xy[1])) - 1
    print("box.lats:", south, north)

    east = np.ceil(max(basin.exterior.coords.xy[0])) + 1
    west = np.floor(min(basin.exterior.coords.xy[0])) - 1
    print("box.lons:", west, east)

    box_lats = np.arange(north, south + LATITUDE_DELTA, LATITUDE_DELTA)
    box_lons = np.arange(west, east + LONGITUDE_DELTA, LONGITUDE_DELTA)

    return box_lats, box_lons


def basin_raster(basin, latitudes, longitudes):
    """
    Rasterize the basin polygon into a boolean mask with the shape of the data

    Parameters
    ----------
    basin : LinearRing
        Basin boundaries.

    Returns
    -------
    mask : bool[:, :]
        Rasterized polygon.

    """
    coords = np.array(basin.xy).T
    path = Path(coords)
    width, height = longitudes.size, latitudes.size
    xmg, ymg = np.meshgrid(longitudes, latitudes)
    points = np.vstack((xmg.ravel(), ymg.ravel())).T
    # https://stackoverflow.com/questions/3654289/scipy-create-2d-polygon-mask
    mask = path.contains_points(points).reshape(height, width)
    return mask


def eastward_masks(basin):
    """
    Identify the input and output cells for an eastward flux.

    Parameters
    ----------
    basin : bool[:, :]
        Boolean mask identifying the basin cells.

    Returns
    -------
    eastward_in : bool[:, :]
        Boolean mask identifying the input cells for an eastward flux.
    eastward_out : bool[:, :]
        Boolean mask identifying the output cells for an eastward flux.

    """
    height, width = basin.shape

    eastward_in = np.zeros_like(basin, dtype="bool")
    eastward_out = np.zeros_like(basin, dtype="bool")

    for row in range(height):
        if basin[row].any():
            target = True
            for col in range(width):
                if basin[row, col] == target:
                    if target:
                        eastward_in[row, col] = True
                    else:
                        eastward_out[row, col] = True
                    target = not target

    return eastward_in, eastward_out


def northward_masks(basin):
    """
    Identify the input and output cells for a northward flux.

    Parameters
    ----------
    basin : bool[:, :]
        Boolean mask identifying the basin cells.

    Returns
    -------
    northward_in : bool[:, :]
        Boolean mask identifying the input cells for a northward flux.
    northward_out : bool[:, :]
        Boolean mask identifying the output cells for a northward flux.

    """
    height, width = basin.shape

    northward_in = np.zeros_like(basin, dtype="bool")
    northward_out = np.zeros_like(basin, dtype="bool")

    for col in range(width):
        if basin[:, col].any():
            target = True
            for row in range(height):
                if basin[row, col] == target:
                    if target:
                        northward_out[row, col] = True
                    else:
                        northward_in[row, col] = True
                    target = not target

    return northward_in, northward_out


def make_masks(basin_poly, latitudes, longitudes):
    """
    Rasterizes the polygon and identifies the contour segments for the entry
    an exit cells of Eastward and Northward fluxes.

    Parameters
    ----------
    basin_poly : shapely.geometry.polygon.LinearRing
        Shapely geometry read from a Shapefile.

    Returns
    -------
    BasinMask
        Named tuple with eastward_in, eastward_out, northward_in,
        northward_out, and basin masks.

    """
    basin = basin_raster(basin_poly, latitudes, longitudes)
    print(basin.shape)

    eastward_in, eastward_out = eastward_masks(basin)
    northward_in, northward_out = northward_masks(basin)

    return BasinMask(eastward_in, eastward_out, northward_in, northward_out, basin)


def plot_masks(basin_name, box_lons, box_lats, masks):
    """Plot the masks."""

    projection = PlateCarree()

    axes_class = (GeoAxes, {"projection": projection})
    fig = plt.figure()
    grid = AxesGrid(
        fig,
        111,
        axes_class=axes_class,
        nrows_ncols=(3, 2),
        axes_pad=0.1,
    )

    grid[0].pcolor(box_lons, box_lats, masks.basin, shading="auto")
    grid[1].pcolor(
        box_lons,
        box_lats,
        np.logical_or(
            np.logical_or(masks.eastward_in, masks.eastward_out),
            np.logical_or(masks.northward_in, masks.northward_out),
        ),
        shading="auto",
    )

    grid[2].pcolor(box_lons, box_lats, masks.eastward_in, shading="auto")
    grid[3].pcolor(box_lons, box_lats, masks.eastward_out, shading="auto")

    grid[4].pcolor(box_lons, box_lats, masks.northward_in, shading="auto")
    grid[5].pcolor(box_lons, box_lats, masks.northward_out, shading="auto")

    print(f"Saving {basin_name}_masks.png")
    fig.savefig(
        f"output/{basin_name}_masks.png", dpi=600, bbox_inches="tight", pad_inches=0.1
    )


def cell_geometry(box_lons, box_lats):
    """
    Computes the area and side lenghts of each cell. Northern and Southern
    sides are different, but Eastern and Western sides are equal. Values vary
    with latitude but not longitude.

    Returns
    -------
    east_length : float[:, :]
        Cell eastern-side's length in meters.
        eastern-and western side are equal.
    north_length : float[:, :]
        Cell northern-side's length in meters.
    south_length : float[:, :]
        Cell southern-side's length in meters.
    cell_area : float[:, :]
        Cell area in square meters.

    """
    east_length = np.zeros(shape=(box_lats.size, box_lons.size))
    north_length = np.zeros(shape=(box_lats.size, box_lons.size))
    south_length = np.zeros(shape=(box_lats.size, box_lons.size))
    cell_area = np.zeros(shape=(box_lats.size, box_lons.size))

    for j in range(box_lats.size):

        # https://gis.stackexchange.com/a/245879
        # https://geographiclib.sourceforge.io/html/python/examples.html

        east_length[j] = WGS84.Inverse(
            box_lats[j] + np.abs(LATITUDE_DELTA / 2),
            box_lons[1],
            box_lats[j] - np.abs(LATITUDE_DELTA / 2),
            box_lons[1],
        )["s12"]

        north_length[j] = WGS84.Inverse(
            box_lats[j] + np.abs(LATITUDE_DELTA / 2),
            box_lons[1],
            box_lats[j] + np.abs(LATITUDE_DELTA / 2),
            box_lons[0],
        )["s12"]

        south_length[j] = WGS84.Inverse(
            box_lats[j] - np.abs(LATITUDE_DELTA / 2),
            box_lons[1],
            box_lats[j] - np.abs(LATITUDE_DELTA / 2),
            box_lons[0],
        )["s12"]

        cell_area[j] = np.abs(
            # N Segment
            Geodesic.Inverse(
                box_lats[j] + np.abs(LATITUDE_DELTA / 2),
                box_lons[1],
                box_lats[j] + np.abs(LATITUDE_DELTA / 2),
                box_lons[0],
                Geodesic.AREA,
            )["S12"]
            -
            # S Segment
            WGS84.Inverse(
                box_lats[j] - np.abs(LATITUDE_DELTA / 2),
                box_lons[1],
                box_lats[j] - np.abs(LATITUDE_DELTA / 2),
                box_lons[0],
                Geodesic.AREA,
            )["S12"]
        )

    return CellGeometry(east_length, east_length, north_length, south_length, cell_area)


def plot_geoms(basin_name, box_lons, box_lats, geoms, masks):
    """Plot the cell geometries."""

    projection = PlateCarree()

    axes_class = (GeoAxes, {"projection": projection})
    fig = plt.figure()
    grid = AxesGrid(
        fig,
        111,
        axes_class=axes_class,
        nrows_ncols=(2, 2),
        axes_pad=0.1,
        cbar_location="right",
        cbar_mode="edge",
        cbar_pad=0.1,
        cbar_size="3%",
    )
    vmin = min(
        np.ma.masked_where(~masks.eastward_in, geoms.western_side_length).min(),
        np.ma.masked_where(~masks.eastward_out, geoms.eastern_side_length).min(),
    )
    vmax = max(
        np.ma.masked_where(~masks.eastward_in, geoms.western_side_length).max(),
        np.ma.masked_where(~masks.eastward_out, geoms.eastern_side_length).max(),
    )

    west = grid[0].pcolor(
        box_lons,
        box_lats,
        np.ma.masked_where(~masks.eastward_in, geoms.western_side_length),
        shading="auto",
        vmin=vmin,
        vmax=vmax,
    )

    grid[1].pcolor(
        box_lons,
        box_lats,
        np.ma.masked_where(~masks.eastward_out, geoms.eastern_side_length),
        shading="auto",
        vmin=vmin,
        vmax=vmax,
    )

    vmin = min(
        np.ma.masked_where(~masks.northward_in, geoms.southern_side_length).min(),
        np.ma.masked_where(~masks.northward_out, geoms.northern_side_length).min(),
    )
    vmax = max(
        np.ma.masked_where(~masks.northward_in, geoms.southern_side_length).max(),
        np.ma.masked_where(~masks.northward_out, geoms.northern_side_length).max(),
    )

    south = grid[2].pcolor(
        box_lons,
        box_lats,
        np.ma.masked_where(~masks.northward_in, geoms.southern_side_length),
        shading="auto",
        vmin=vmin,
        vmax=vmax,
    )

    grid[3].pcolor(
        box_lons,
        box_lats,
        np.ma.masked_where(~masks.northward_out, geoms.northern_side_length),
        shading="auto",
        vmin=vmin,
        vmax=vmax,
    )

    grid.cbar_axes[0].colorbar(west)
    grid.cbar_axes[1].colorbar(south)

    print(f"Saving {basin_name}_geoms.png")
    fig.savefig(
        f"output/{basin_name}_geoms.png", dpi=600, bbox_inches="tight", pad_inches=0.1
    )


def get_latlon_indices(box_lats, box_lons):
    """Get the indices for the box coordinates."""
    west = box_lons[0]
    east = box_lons[-1]
    north = box_lats[0]
    south = box_lats[-1]

    # TIP: Any NC file will work here. We only want the lats/lons
    with MFDataset(VIWVN.filename, aggdim="time") as dataset:

        ix_west, ix_east = np.searchsorted(
            dataset.variables["longitude"][:], wrap_lon360([west, east])
        )
        ix_east += 1
        lons = wrap_lon180(dataset.variables["longitude"][int(ix_west) : int(ix_east)])

        print("ix_lons:", ix_west, ix_east)

        assert np.allclose(lons, box_lons)

        ix_south, ix_north = dataset.variables["latitude"][:].size - np.searchsorted(
            dataset.variables["latitude"][:][::-1], [south, north]
        )
        ix_north -= 1
        lats = dataset.variables["latitude"][int(ix_north) : int(ix_south)]

        print("ix_lats:", ix_north, ix_south)

        assert np.allclose(lats, box_lats)

    return Coordinates(slice(ix_north, ix_south), slice(ix_west, ix_east))


def repeat_mask_ntimes(basin_mask, ntimes):
    """
    Repeat a basin_mask ntimes for array masking of arrays including a time
    dimenson

    Parameters
    ----------
    basin_mask : BasinMask
        Named tuple of masks.
    ntimes : int
        Numer of times to repeat.

    Returns
    -------
    BasinMask
        Named tuple of repeated masks.

    """

    eastward_in = np.repeat(basin_mask.eastward_in[np.newaxis, ...], ntimes, axis=0)
    northward_in = np.repeat(basin_mask.northward_in[np.newaxis, ...], ntimes, axis=0)
    eastward_out = np.repeat(basin_mask.eastward_out[np.newaxis, ...], ntimes, axis=0)
    northward_out = np.repeat(basin_mask.northward_out[np.newaxis, ...], ntimes, axis=0)
    basin_area = np.repeat(basin_mask.basin[np.newaxis, ...], ntimes, axis=0)

    return BasinMask(eastward_in, eastward_out, northward_in, northward_out, basin_area)


def flux_balance_out(northward, eastward, slices, masks, geoms):
    """Balances the horizontal fluxes."""
    with MFDataset(eastward.filename, aggdim="time") as dataset:
        flux_e = dataset.variables[eastward.variable][
            :, slices.latitude, slices.longitude
        ]

    with MFDataset(northward.filename, aggdim="time") as dataset:
        flux_n = dataset.variables[northward.variable][
            :, slices.latitude, slices.longitude
        ]

    ntimes = flux_e.shape[0]
    many_masks = repeat_mask_ntimes(masks, ntimes)

    flux_ei = np.sum(
        geoms.western_side_length[many_masks.eastward_in[0]]
        * flux_e[many_masks.eastward_in].reshape(ntimes, -1),
        axis=1,
    )
    flux_eo = np.sum(
        geoms.eastern_side_length[many_masks.eastward_out[0]]
        * flux_e[many_masks.eastward_out].reshape(ntimes, -1),
        axis=1,
    )

    flux_ni = np.sum(
        geoms.southern_side_length[many_masks.northward_in[0]]
        * flux_n[many_masks.northward_in].reshape(ntimes, -1),
        axis=1,
    )
    flux_no = np.sum(
        geoms.northern_side_length[many_masks.northward_out[0]]
        * flux_n[many_masks.northward_out].reshape(ntimes, -1),
        axis=1,
    )

    return (flux_ei - flux_eo) + (flux_ni - flux_no)


def integrate_rate(ncdata, slices, basin_mask, cell_geom):
    """
    Integrates a variable for the BASIN_AREA as the sum of the product of the
    variable's value by the cell's area at each cell inside the BASIN_AREA.


    Parameters
    ----------
    dataset : netCDF4.Dataset
        File descriptor to an open dataset.
    variable : str
        Name of the variable inside the dataset.

    Returns
    -------
    integrated : float(num_times) [units: m**3 s**-1]
        Sum of the product of the variable's value by the cell's area at each
        cell inside the BASIN_AREA.

    """
    with MFDataset(ncdata.filename, aggdim="time") as dataset:
        print(dataset.variables[ncdata.variable].units)
        var = dataset.variables[ncdata.variable][
            :, slices.latitude, slices.longitude
        ]  # m
        ntimes = dataset.variables["time"][...].size

    many_masks = repeat_mask_ntimes(basin_mask, ntimes)

    integrated = np.sum(
        cell_geom.cell_area[many_masks.basin[0]]
        * var[many_masks.basin].reshape(ntimes, -1),
        axis=1,
    )

    # print(integrated.shape)

    return integrated  # m ** 3


def zonal_average(ncdata, slices, basin_mask):
    """
    Integrates a variable for the BASIN_AREA as the sum of the product of the
    variable's value by the cell's area at each cell inside the BASIN_AREA.


    Parameters
    ----------
    dataset : netCDF4.Dataset
        File descriptor to an open dataset.
    variable : str
        Name of the variable inside the dataset.

    Returns
    -------
    integrated : float(num_times) [units: kg s**-1]
        Sum of the product of the variable's value by the cell's area at each
        cell inside the BASIN_AREA.

    """
    with MFDataset(ncdata.filename, aggdim="time") as dataset:
        print(dataset.variables[ncdata.variable])
        var = dataset.variables[ncdata.variable][
            :, slices.latitude, slices.longitude
        ]  # m
        ntimes = dataset.variables["time"][...].size

    many_masks = repeat_mask_ntimes(basin_mask, ntimes)

    integrated = np.mean(var[many_masks.basin].reshape(ntimes, -1), axis=1)

    print(integrated.shape)

    return integrated  # m ** 3


def moisture_convergence(basin_name, station_name):
    """Computes moisture convergence for the basin."""
    basin = shapefile_reader(basin_name, station_name)
    box_lats, box_lons = get_box_coordinates(basin)

    masks = make_masks(basin.exterior, box_lats, box_lons)  # 6147 celdas
    plot_masks(basin_name, box_lons, box_lats, masks)

    cell_geoms = cell_geometry(box_lons, box_lats)
    plot_geoms(basin_name, box_lons, box_lats, cell_geoms, masks)
    slices = get_latlon_indices(box_lats, box_lons)

    vapor_balance = flux_balance_out(VIWVN, VIWVE, slices, masks, cell_geoms)
    liquid_balance = flux_balance_out(VILWN, VILWE, slices, masks, cell_geoms)
    ice_balance = flux_balance_out(VIIWN, VIIWE, slices, masks, cell_geoms)

    precipitation = integrate_rate(PREC, slices, masks, cell_geoms)
    evaporation = integrate_rate(EVAP, slices, masks, cell_geoms)

    temperature = zonal_average(TEMP, slices, masks)
    output = np.stack(
        [
            vapor_balance,
            liquid_balance,
            ice_balance,
            precipitation,
            evaporation,
            temperature,
        ]
    ).T
    print("output", output.shape)
    np.savetxt(
        f"output/{basin_name}_output.csv",
        output,
        delimiter=",",
        newline="\n",
        fmt="%.2f",
        comments="",
        header=",".join(
            [
                "VAPOR_WATER_FLUX_KG/S",
                "CLOUD_LIQUID_WATER_FLUX_KG/S",
                "CLOUD_FROZEN_WATER_FLUX_KG/S",
                "PRECIPITATION_KG/S",
                "EVAPORATION_KG/S",
                "2M_TEMPERATURE_K",
            ]
        ),
    )

    return output


def main():
    """Main function."""
    # download_all()  # Uncomment to download
    plt.style.use("ggplot")
    moisture_convergence("AMAZONAS", "OBIDOS")
    moisture_convergence("PARANA", "TIMBUES")
    moisture_convergence("MISSISSIPPI RIVER", "VICKSBURG, MS")
    moisture_convergence("OB", "SALEKHARD")
    moisture_convergence("YENISEY", "IGARKA")
    moisture_convergence("CONGO", "KINSHASA")


if __name__ == "__main__":
    main()
