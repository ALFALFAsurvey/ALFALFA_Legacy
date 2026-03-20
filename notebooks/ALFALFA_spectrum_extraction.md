---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.19.1
kernelspec:
  display_name: Python [conda env:AALegacy]
  language: python
  name: conda-env-AALegacy-py
---

```{code-cell} ipython3
import os, sys
import numpy, astropy
import astropy.io.fits as fits
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.table import Table
import astropy.units as u
import matplotlib.pyplot as plt
from photutils.aperture import SkyEllipticalAperture, EllipticalAperture, SkyRectangularAperture, RectangularAperture 

#Add the scripts directory path
scripts_path = os.path.abspath("../scripts")
sys.path.append(scripts_path)

from load_grid import load_grid
```

This notebook takes an ALFALFA cube and extract a 1-dimensional spectrum at a specified location. It does this first over an area equal to the beam shape and then over a larger square region.

**Note**: You will need the ALFALFA 1044+13 grid in that data directory in order for this notebook to run successfully. You can obtain this by running the [download_example_data.py](../scripts/download_example_data.py) script in the scripts directory. If you wish to run this with other grids then you can find instructions for accessing the grids in the [grid_access.md](../docs/grid_access.md) file in the docs folder or illustrated instructions on the [wiki](https://github.com/jonesmg/ALFALFA_Legacy/wiki/Grid-access-via-NRAO-archive).

```{code-cell} ipython3
#Define the path to the data
#You should already have run the data download script
cwd = os.getcwd()+'/'
data_path = cwd+'../data/A2010/pipeline.unknown_date/'
```

```{code-cell} ipython3
#Define the grid you are using
grid_ra = '1044'
grid_dec = '13'
freq_slice = 'a'

#Use the load_grid function to load the cube and WCS
cube, freq, vel, grid_wcs, header = load_grid(data_path,grid_ra,grid_dec,freq_slice)

#Average together the two polarizations
cube = numpy.mean(cube,axis=0)

#Extract beam parameters
bmaj, bmin, bpa = header['BMAJ']*u.deg, header['BMIN']*u.deg, 0*u.deg
pixsize = header['CDELT2']*u.deg

#Make the beam correction factor
beam_factor = (numpy.pi*bmaj*bmin/(pixsize**2.))/(4.*numpy.log(2.))

#Extract the channel width from the header and build frequency array
chan_df = header['CDELT3']*u.MHz
```

```{code-cell} ipython3
#Start by extracting a spectrum in the shape of the beam

#Define a position to extract a spectrum at
#In this case we will use NGC3338
spec_pos = SkyCoord(160.531332, 13.747050, unit=u.deg)

#Make an aperture that's at least twice the size of the beam
sky_aperture = SkyEllipticalAperture(spec_pos, a=bmaj, b=bmin, theta=bpa)

#Convert to a pixel aperture
pix_aperture = sky_aperture.to_pixel(grid_wcs.celestial)
```

```{code-cell} ipython3
#Plot the aperture on the grid

fig = plt.figure(figsize=(8,8))
ax = plt.subplot(projection=grid_wcs.celestial)
plt.imshow(cube.sum(axis=0), origin='lower', cmap='Greys', aspect='equal')

pix_aperture.plot(color='orange',lw=1.5,ls='-')

plt.xlabel(r'RA')
plt.ylabel(r'Dec')
```

```{code-cell} ipython3
#Now extract the spectrum
pix_aperture_mask = pix_aperture.to_mask(method='exact')

#Create empty spectrum
spec = numpy.zeros(len(freq))

#Step through cube channels and extract flux in each
for i in range(len(spec)):
    spec[i] = numpy.sum(pix_aperture_mask.multiply(cube[i]))/beam_factor

#Set units
spec = spec*u.mJy
```

```{code-cell} ipython3
plt.plot(vel,spec)

plt.axhline(0.,c='k',lw=1)
plt.ylim(-10.,600.)

plt.xlabel('Velocity [km/s] (Heliocentric)')
plt.ylabel('Flux Density [mJy]')
```

The spectral profile looks quite peculiar, this is likely because the source is extended and we are therefore not including all the emission.</br></br> Now let's switch to using a large 20'x20' square extraction region.

```{code-cell} ipython3
#Make a rectangular aperture that's 20x20 arcmin
sky_aperture = SkyRectangularAperture(spec_pos,w=20*u.arcmin,h=20*u.arcmin)

#Convert to a pixel aperture
pix_aperture = sky_aperture.to_pixel(grid_wcs.celestial)
```

```{code-cell} ipython3
#Plot the aperture on the grid

fig = plt.figure(figsize=(8,8))
ax = plt.subplot(projection=grid_wcs.celestial)
plt.imshow(cube.sum(axis=0), origin='lower', cmap='Greys', aspect='equal')

pix_aperture.plot(color='orange',lw=1.5,ls='-')

plt.xlabel(r'RA')
plt.ylabel(r'Dec')
```

```{code-cell} ipython3
#Now extract the spectrum
pix_aperture_mask = pix_aperture.to_mask(method='exact')

#Create empty spectrum
spec = numpy.zeros(len(freq))

#Step through cube channels and extract flux in each
for i in range(len(spec)):
    spec[i] = numpy.sum(pix_aperture_mask.multiply(cube[i]))/beam_factor

#Set units
spec = spec*u.mJy
```

```{code-cell} ipython3
plt.plot(vel,spec)

plt.axhline(0.,c='k',lw=1)
plt.ylim(-50.,1200.)

plt.xlabel('Velocity [km/s] (Heliocentric)')
plt.ylabel('Flux Density [mJy]')
```

Now let's compare with the spectrum catalogued in ALFALFA. You will see that a.100 extracted spectrum misses half of the flux. This is expected for the extended source as noted in Haynes+18. 

```{code-cell} ipython3
alf_spec = Table.read("https://vizier.cds.unistra.fr/viz-bin/ftp-index?J/ApJ/861/49/sp/A005826.fits")
```

```{code-cell} ipython3
plt.plot(vel,spec, label="Extracted")
plt.plot(alf_spec["VHELIO"], alf_spec["FLUX"], label="a.100")

plt.axhline(0.,c='k',lw=1)
plt.ylim(-50.,1200.)

plt.xlabel('Velocity [km/s] (Heliocentric)')
plt.ylabel('Flux Density [mJy]')

plt.legend()
```
