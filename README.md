# microphysics_code
Standalone [RegCM](https://github.com/ICTP/RegCM/) microphysics code driven by
[ERA5](https://cds-beta.climate.copernicus.eu/datasets/reanalysis-era5-complete?tab=overview)
atmospheric dataset on model level and the original reduced gaussian
[N640](https://confluence.ecmwf.int/display/EMOS/N640) grid.

# Description
The aim of this software is to provide a reference baseline for optimization.

# Data
Data to run the application needs to be downloaded from ECMWF. In the *data*
directory, the *mlevel* and *surface* directories must be populated. The user
is provided with two python scripts using the [CDS API](https://cds.climate.copernicus.eu/api-how-to)
that can be used with a *python3* to download the dataset. Static data populate
the other directories.

# Compilation
To compile the code, a Fortran 2003 compiler is required (tested with GNU
*gfortran9* on Ubuntu Linux), along with the [netCDF](https://www.unidata.ucar.edu/software/netcdf/)
Fortran library (tested with version 4.6.1) and the *make* program
(tested with GNU *make* 4.2.1). Type

`make`

and the program *microphysics_code* is built and can be run. Only text output
is provided, any modification involves editing the Fortran soirce code and is
left to the user.

# Happy hacking.

Graziano.
