# Overview

This repository is intended to be a central location to collect and document code developed for interacting with the [ALFALFA](https://egg.astro.cornell.edu/alfalfa/data/index.php) data cubes hosted at [NRAO](data.nrao.edu).

# Getting Started

The basic function of this repo is to ensure that all the code for interacting with the ALFALFA cubes uses a single Python environment and is stored in one place.

## I don't want to install anything

Although it is strongly recommended to download this repository and build the Python environment if you wish to use or adapt the code hosted here, if you don't want to install or download anything then the simplest way to use this repository is to launch it in the cloud using [mybinder](https://mybinder.org/) by clicking this button: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/jonesmg/ALFALFA_Legcay/HEAD)

Note that when you close (or lose) your connection to mybinder all of the open work that you haven't manually downloaded and saved **will be lost**. This interface is generally intended for testing things out, not for seriously working on code.

## Downloading this repository

If you use git then the easiest way to get the repo is to run the following command in a terminal in a suitable location:

```bash
git clone git@github.com:jonesmg/ALFALFA_Legacy.git
```

If you are planning to work on developing code for this repo then we strongly recommend, using git to allow all the changes to be kept track of. For git installation instructions please click [here](https://git-scm.com/downloads).

If you just want to download the repo, but don't want to develop the code or use git, then you can simply download it as a zip file at the following url:

[https://github.com/jonesmg/ALFALFA_Legacy/archive/refs/heads/main.zip](https://github.com/jonesmg/ALFALFA_Legacy/archive/refs/heads/main.zip)

## Building the Python environment

This Python environment is defined in the file [environment.yml](environment.yml). This is a [YAML](https://yaml.org/spec/1.2.2/#chapter-1-introduction-to-yaml) file that lists all of the specific package versions such that the software environment used can be common and reproducible. 

If you already know what you're doing then you can use whatever tool you want to build and launch this environment. However, we recommend using conda. If you already have a version of conda installed (see below) then navigate to where you downloaded this repo and in a terminal run the command:

```bash
conda env create -f environment.yml
``` 
This will create a new conda environment called "AALegacy" which you can then activate with:

```bash
conda activate AALegacy
```

### Installing Conda

If you don't already have a version of conda installed then you can follow the instructions at the following links to install either [Miniconda](https://www.anaconda.com/docs/getting-started/miniconda/install) or [Anaconda](https://www.anaconda.com/download/success). Both will work the same, but Miniconda might save a little space, if that is a concern.

#Writing Documentation

Please document any code that you add/edit. We recommend having a separate readme file for each tool, and writing those file in Markdown. GitHub will automatically display these files correct if they have the .md extension. [Here](https://www.markdownguide.org/basic-syntax/) is a quick guide for writing in Markdown. If you've never used Markdown before, you can just write in plain text and it'll work just fine.

If you are planning on creating a new tool for interacting with the ALFALFA data then please create a new [branch](https://www.geeksforgeeks.org/how-to-create-a-new-branch-in-git/) of the repo using

```bash
git checkout -b new_of_new_tool
```

before starting to make changes to the code. This will be enormously helpful in keeping track of everything that is developed without leading to a lot of unnecessary conflicts.