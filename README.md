# Overview

This repository contains tools and instructions for accessing and interacting with the data cubes (often referred to as "grids") from the Arecibo Legacy Fast ALFA (Arecibo L-band Feed Array) survey, or ALFALFA. Below are the basic instructions for getting started, and more details about accessing the data, the metadata, and the survey can be found in the [Wiki](https://github.com/jonesmg/ALFALFA_Legacy/wiki) associated with this repository. We strongly recommend reviewing this material if you plan to use these data.

# Getting Started

The ALFALFA grids are stored in the [NRAO data archive](https://data.nrao.edu/portal). The tools in this repository are intended to be used with data cubes downloaded from this archive. These tools are all based in Python (mostly Jupyter notebooks) and in order for them to run as intended, we recommend building a dedicated Python environment using the [environment.yml](environment.yml) file in this repository (instructions below). The notebooks folder contains these notebooks which show worked examples of standard tasks that user are likely to want to perform. For example, identifying which ALFALFA grid is needed to a particular position on the sky, extracting a 1-dimensional spectrum from a data cube, or producing a moment zero map and overlaying it on an optical image. The scripts folder includes a Python file defining a function to read in data cubes. We recommend using this function if interacting with the cubes in Python as this updates the FITS headers to the latest standards and applies an approximate correction to the world coordinate system in the header.

## I don't want to install anything.

Although it is strongly recommended to download this repository and build the Python environment if you wish to use or adapt the code hosted here, if you don't want to install or download anything then the simplest way to use this repository is to launch it in the cloud using [mybinder](https://mybinder.org/) by clicking this button: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/jonesmg/ALFALFA_Legacy/HEAD)

Note that when you close (or lose) your connection to mybinder, all of the open work that you haven't manually downloaded and saved **will be lost**. This interface is generally intended for testing things out, not for seriously working on code or data.

## Downloading this repository

If you use git then the easiest way to get this repo is to run the following command in a terminal in a suitable location:

```bash
git clone https://github.com/jonesmg/ALFALFA_Legacy.git
```

If you are new to git, installation instructions can be found [here](https://git-scm.com/downloads), and [here](https://www.geeksforgeeks.org/basic-git-commands-with-examples/) are some examples of basic git commands.

If you just want to download the repo, but don't want to use git, then you can simply download it as a zip file at the following url:

[https://github.com/jonesmg/ALFALFA_Legacy/archive/refs/heads/main.zip](https://github.com/jonesmg/ALFALFA_Legacy/archive/refs/heads/main.zip)

## Building the Python environment

The Python environment is defined in the file [environment.yml](environment.yml). This is a [YAML](https://yaml.org/spec/1.2.2/#chapter-1-introduction-to-yaml) file that lists all of the specific package versions such that the software environment used can be common and reproducible. 

If you already know what you're doing then you can use whatever tool you want to build and launch this environment. However, we recommend using conda. If you already have a version of conda installed (see below) then navigate to where you downloaded this repo and in a terminal run the command:

```bash
conda env create -f environment.yml
``` 
This will create a new conda environment called "AALegacy". Note that if running on a newer Mac, the environment may fail to build, in which case you should try instead using the command:

```bash
conda env create -f environment.yml --platform osx-64
```

Once you have successfully created the environment, you can then activate it with:

```bash
conda activate AALegacy
```

### Installing Conda

If you don't already have a version of conda installed then you can follow the instructions at the following links to install either [Miniconda](https://www.anaconda.com/docs/getting-started/miniconda/install) or [Anaconda](https://www.anaconda.com/download/success). Both will work the same, but Miniconda might save a little space, if that is a concern.


## Downloading example data

Most of the notebooks require an ALFALFA grid as input to run. To obtain the grid used in most of our examples navigate in your terminal to the scripts directory as run the following command:

```bash
python download_example_data.py
```

This script will automatically download a tar file (approximately 1.5 GB) contain a grid and extract it into a new directory where the notebooks expect to find the data.


## Launching a notebook

After building and activating the AALegacy environment, navigate to your local version of this repository in your terminal and run the command:

```bash
jupyter notebook
```

This will launch a notebook server in your default browser (if the browser fails to open, then try clicking the link for the notebook server in the terminal output), after which you should be able to open the notebooks folder and launch any of the example notebooks.