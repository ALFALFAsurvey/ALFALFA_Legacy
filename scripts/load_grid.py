import astropy.io.fits as fits
from astropy.wcs import WCS
import numpy as np
import astropy
import astropy.units as u

def load_grid(dir_path, grid_ra, grid_dec, freq_slice, include_weights=False):
    '''
    This function loads the requested ALFALFA grid and applies corrections to
    the header keywords and an approximate fix for the world coordinate system.
    The weights cubes is also loaded if requested.

    INPUTS:
    dir_path: (String) Path to the directory containing the data.
    grid_ra: (String) Right Ascension of the grid to be loaded, e.g., "1044".
    grid_dec: (String) Declination of the grid to be loaded, e.g., "13".
    freq_slice: (String) Frequency slice of the ALFALFA data, e.g., "a", "b", "c", or "d".
    include_weights: (Bool) Return weights cube as well as data (optional).

    OUTPUTS:
    cube: (Array) Data cube.
    freq: (Array) Frequency array for the cube (MHz).
    vel: (Array) Heliocentric velocity array for the cube (km/s).
    wcs: (WCS Object) World Coordinate System object for the grid (2-dimensional).
    header: (Header Object) FITS header of the cube.
    (Optional) weights_cube: (Array) Weights cube.
    (Optional) weights_header: (Header Object) FITS header of the weights cube.
    '''

    grid_filename = f'{dir_path}{grid_ra}+{grid_dec}{freq_slice}_spectral.fits'

    #Open the fits file
    hdu = fits.open(grid_filename)

    #Extract the data from fits file
    cube = hdu[0].data

    #Start by extracting the channel width from the header and build frequency array
    chan_df = hdu[0].header['CDELT3']
    freq = (hdu[0].header['CRVAL3'] + hdu[0].header['CDELT3']*np.arange(0,1024,1))*u.MHz
    
    #Now build velocity array
    rest_freq = (hdu[0].header['RESTFREQ']/1E6)*u.MHz
    vel = ((astropy.constants.c.value/1000)*(rest_freq-freq)/freq)*u.km/u.s

    #Build new WCS
    #Note that the ALFALFA source positions were corrected as described in section 3.4.2
    #of Brian Kent's PhD thesis: https://ecommons.cornell.edu/items/a278d737-b863-4e07-b637-88587dff669f
    #However, these corrections were NOT applied to the grids themselves.
    #Here we construct an approximate WCS that should be suitable for most applications.
    wcs_new = WCS(naxis = 2)
    wcs_new.wcs.crpix = [72.5,72.5]
    wcs_new.wcs.cdelt = [-1/60, 1/60]
    cen_ra = (float(grid_ra[:2]) + float(grid_ra[2:])/60)*15
    cen_dec = float(grid_dec)
    wcs_new.wcs.crval = [cen_ra, cen_dec]
    wcs_new.wcs.ctype = ["RA---CAR", "DEC--CAR"]

    #Update the header with the same WCS
    header_new = hdu["PRIMARY"].header.copy()
    header_new['CRVAL1'] = cen_ra
    header_new['CDELT1'] = -1/60
    header_new['CRPIX1'] = 72.5
    header_new['CRVAL2'] = cen_dec
    header_new['CDELT2'] = 1/60
    header_new['CRPIX2'] = 72.5

    #Update various header keywords to latest standards
    header_new["CTYPE3"]  = "FREQ"
    header_new["CUNIT3"]  = "MHz"
    header_new["SPECSYS"] = "HELIOCEN"
    header_new["RESTFRQ"] = 1420.405751768e6  # (Hz)
    header_new["CNAME3"]  = "FREQ-HEL"
    header_new["BMAJ"]  = 3.8/60 #arcmin
    header_new["BMIN"]  = 3.3/60 #arcmin
    header_new["BPA"]   = 0
    header_new["CRVAL4"] = -2  # LL and RR
    
    if include_weights:
        wgts_filename = f'{dir_path}{grid_ra}+{grid_dec}{freq_slice}_spectralweights.fits'
        
        #Open the fits file
        wgts_hdu = fits.open(wgts_filename)
    
        #Extract the data from fits file
        wgts_cube = wgts_hdu[0].data
    
        #Update the header with the same WCS
        wgts_header_new = wgts_hdu["PRIMARY"].header.copy()
        wgts_header_new['CRVAL1'] = cen_ra
        wgts_header_new['CDELT1'] = -1/60
        wgts_header_new['CRPIX1'] = 72.5
        wgts_header_new['CRVAL2'] = cen_dec
        wgts_header_new['CDELT2'] = 1/60
        wgts_header_new['CRPIX2'] = 72.5
    
        #Update various header keywords to latest standards
        wgts_header_new["CTYPE3"]  = "FREQ"
        wgts_header_new["CUNIT3"]  = "MHz"
        wgts_header_new["SPECSYS"] = "HELIOCEN"
        wgts_header_new["RESTFRQ"] = 1420.405751768e6  # (Hz)
        wgts_header_new["CNAME3"]  = "FREQ-HEL"
        wgts_header_new["BMAJ"]  = 3.8/60 #arcmin
        wgts_header_new["BMIN"]  = 3.3/60 #arcmin
        wgts_header_new["BPA"]   = 0
        wgts_header_new["CRVAL4"] = -2  # LL and RR

        return cube, freq, vel, wcs_new, header_new, wgts_cube, wgts_header_new
    else:
        return cube, freq, vel, wcs_new, header_new