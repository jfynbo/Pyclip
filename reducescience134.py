# This is a python program to make a Flat frame
print('Script running')
nframes = 7
xsize = 2048
ysize = 2064
xoutstart = 155
youtstart = 130
xoutend = 1970
youtend = 1955
badsig = 1.e13

from astropy.io import fits
import numpy as np
from astropy.io.fits import getheader

Outnames = ['data1340001.fits','data1340002.fits','data1340003.fits','data1340004.fits','data1340005.fits','data1340006.fits','data1340007.fits']

#Read in the raw flat frames and subtact mean of overscan region and BIAS and
#divide with the flat field
list = open('scienceNB.list')
BIASframe = fits.open('BIAS.fits')
BIAS = np.array(BIASframe[0].data)
FLATframe = fits.open('Flat134.fits')
FLAT = np.array(FLATframe[0].data)
for i in range(0,nframes):
   outframe = np.zeros((ysize,xsize))
   print('Image number:', i)
   filenm = list.readline()
   imagename = filenm.strip("\n")
   rawimage = fits.open(imagename)
   hdr = rawimage[1].header
   print('Info on file:')
   print(rawimage.info())
   data = np.array(rawimage[1].data)
   median = np.mean(data[1:2048,1:50])
   data = data - median
   print('Subtracted the median value of the overscan :',median)
   outframe[0:ysize-1,0:xsize-1] = data[0:ysize-1,51:51+xsize-1]
   outframe = outframe-BIAS
   outframe = outframe/FLAT
   badsigpix = np.where(FLAT < 0.5)
   outframe[badsigpix] = badsig
   outframe = outframe[youtstart:youtend,xoutstart:xoutend]
   print(Outnames[i])
   fits.writeto(Outnames[i],outframe,hdr,overwrite=True)
list.close()


#Also save the trimmed flat-frame for combination code
FLAT = FLAT[youtstart:youtend,xoutstart:xoutend]
fits.writeto('Flat134_trim.fits',FLAT,overwrite=True)
