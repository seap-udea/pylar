# PyLAR
## Python utilities for the LAR model

<!-- This are visual tags that you may add to your package at the beginning with useful information on your package --> 
[![version](https://img.shields.io/pypi/v/ipylar?color=blue)](https://pypi.org/project/ipylar/)
[![downloads](https://img.shields.io/pypi/dw/ipylar)](https://pypi.org/project/ipylar/)
[![license](https://img.shields.io/pypi/l/ipylar)](https://pypi.org/project/ipylar/)
[![implementation](https://img.shields.io/pypi/implementation/ipylar)](https://pypi.org/project/ipylar/)
[![pythonver](https://img.shields.io/pypi/pyversions/ipylar)](https://pypi.org/project/ipylar/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.12167662.svg)](https://doi.org/10.5281/zenodo.12167662)

<!--
<a target="_blank" href="https://colab.research.google.com/github/seap-udea/fargopy/blob/main/README.ipynb">
  <img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/>
</a>
-->

`PyLAR` is a python package and datasets intended to test and use the LAR (Land Atmospheric and Reservoir) model. The LAR model is intended to describe the changes in the water storage in large river basin around the world, including
atmospheric processes as a critical component of the basin
water budget. 

For the science behind the LAR model please refer to the following paper:

> Juan F. Salazar, Rubén D. Molina, Jorge I. Zuluaga, and Jesus D. Gomez-Velez (2024), **Wetting and drying trends in the land–atmosphere reservoir of large basins around the world**, [Hydrology and Earth System Sciences, in publication (2024)](https://hess.copernicus.org/preprints/hess-2023-172/), [doi.org/10.5194/hess-2023-172](https://doi.org/10.5194/hess-2023-172).

All the notebooks and data required to reproduce the results of this paper, and other papers produced by our group, are available in the `dev` directory in this repository. [Link](https://github.com/seap-udea/pylar/tree/main/dev).

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
$ import ipylar as lar
```

You can load the data using:

```python
$ amazonas = lar.Basin('amazonas')
```

## What's new

For a detailed list of the newest characteristics of the code see the file [What's new](https://github.com/seap-udea/pylar/blob/master/WHATSNEW.md).

------------

This package has been designed and written by Jorge I. Zuluaga, Ruben D. Molina, Juan F. Salazar and Jesus D. Gomez-Velez (C) 2024