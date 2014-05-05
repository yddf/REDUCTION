#!/usr/bin/env python

from astropy.io import fits

def addObjectName(filenm, object_name):
	#filenm = '/raw/yddf/140306/yddf140306.3000.fits'
	hdulist = fits.open(filenm, mode='update')

	hdr = hdulist[0].header
	hdr.set('object', object_name, 'name of object observed')

	#hdulist.writeto(filenm)
	hdulist.close()

import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Add OBJECT keyword to YDDF-taken FITS files.')
    parser.add_argument('filename',
        help='Name of input FITS file', type=str)
    parser.add_argument('object_name',
        help='The name of the object to add to the FITS header', type=int)
    args = parser.parse_args()

    addObjectName(args.filename,args.object_name)
