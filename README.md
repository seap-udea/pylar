# PyLAR
## Python utilities for the LAR model

<!-- This are visual tags that you may add to your package at the beginning with useful information on your package --> 
[![version](https://img.shields.io/pypi/v/ipylar?color=blue)](https://pypi.org/project/ipylar/)
[![downloads](https://img.shields.io/pypi/dw/ipylar)](https://pypi.org/project/ipylar/)
[![license](https://img.shields.io/pypi/l/ipylar)](https://pypi.org/project/ipylar/)
[![implementation](https://img.shields.io/pypi/implementation/ipylar)](https://pypi.org/project/ipylar/)
[![pythonver](https://img.shields.io/pypi/pyversions/ipylar)](https://pypi.org/project/ipylar/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.12167661.svg)](https://doi.org/10.5281/zenodo.12167661)
[![DOI](https://img.shields.io/badge/10.5194%2Fhess-2023-172)](https://doi.org/10.5194/hess-2023-172)
<!-- Static badge generator: https://shields.io/badges/static-badge: put the badge number and ready -->

`PyLAR` is a python package and datasets intended to test and use the LAR (Land-Atmospheric and Reservoir) model. The LAR model is intended to describe the changes in the water storage in large river basin around the world, including
atmospheric processes as a critical component of the basin
water budget. 

<p align="center"><img src="https://github.com/seap-udea/pylar/blob/main/tutorials/resources/LAR-Conceptual.png?raw=true" alt="Conceptual illustration of LAR"/></p>

For the science behind the LAR model please refer to the following paper:

> Juan F. Salazar, Rubén D. Molina, Jorge I. Zuluaga, and Jesus D. Gomez-Velez (2024), **Wetting and drying trends in the land–atmosphere reservoir of large basins around the world**, [Hydrology and Earth System Sciences, in publication (2024)](https://hess.copernicus.org/preprints/hess-2023-172/), [doi.org/10.5194/hess-2023-172](https://doi.org/10.5194/hess-2023-172).

All the notebooks and data required to reproduce the results of this paper, and other papers produced by our group, are available in the [`dev` directory in this repository](https://github.com/seap-udea/pylar/tree/main/dev).

## Downloading and Installing `PyLAR` 

`PyLAR` is available at the `Python` package index and can be installed using:

```bash
$ sudo pip install ipylar
```
as usual this command will install all dependencies and download some useful data, scripts and constants.

> **NOTE**: If you don't have access to `sudo`, you can install `PyLAR` in your local environmen (usually at `~/.local/`). In that case you need to add to your `PATH` environmental variable the location of the local python installation. Add to `~/.bashrc` the line `export PATH=$HOME/.local/bin:$PATH`

## Quickstart

To start using `PyLAR`, you should first obtain data for a large river basin. We have provided with the package a dataset especially prepared for the Amazonas Basin we will use in this quickstart. 

You must start by importing the package:

```python
import ipylar as lar
```

Create a basin:

```python
amazonas = lar.Basin(key='amazonas',name='Amazonas')
```

Once created, you should read the timeseries for the basin and load it into the pandas dataframe `amazonas.data`. The present version of `PyLAR` includes sample data. You may read the sample data using:

```python
amazonas.read_basin_data()
```

Once the data is loaded you can perform operations on the data, for instance, you can plot it:

```python
fig = amazonas.plot_basin_series()
```

<p align="center"><img src="https://github.com/seap-udea/pylar/blob/main/tutorials/resources/amazonas-lar-timeseries.png?raw=true" alt="Amazonas LAR time-series"/></p>

## Tutorials

We have prepared a set of [basic tutorials](https://github.com/seap-udea/pylar/tree/main/tutorials) for illustrating the usage of some of the tools including in `PyLAR`. The tutorials can be ran in `Google Colab`.

## What's new

For a detailed list of the newest characteristics of the code see the file [What's new](https://github.com/seap-udea/pylar/blob/master/WHATSNEW.md).

------------

This package has been designed and written by Jorge I. Zuluaga, Ruben D. Molina, Juan F. Salazar and Jesus D. Gomez-Velez (C) 2024