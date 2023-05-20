# This is a python program to make a Flat frame
print('Script running')
nframes = 19
ysize = 2064
xsize = 2048

from astropy.io import fits
import numpy

#Read in the raw flat frames and subtact mean of overscan region and BIAS and
#normalise
list = open('Flats134.list')
bigflat = numpy.zeros((nframes,ysize,xsize),float)
BIASframe = fits.open('BIAS.fits')
BIAS = numpy.array(BIASframe[0].data)
#bigflat = numpy.zeros((nframes,3,3))
for i in range(0,nframes):
   print('Image number:', i)
   filenm = list.readline()
   rawflat = fits.open('../raw/'+filenm.strip("\n"))
   print('Info on file:')
   print(rawflat.info())
   data = numpy.array(rawflat[1].data)
   median = numpy.mean(data[1:2064,1:50])
   data = data - median
   print('Subtracted the median value of the overscan :',median)
   bigflat[i-1,0:ysize-1,0:xsize-1] = data[0:ysize-1,51:51+xsize-1]
   bigflat[i-1,:,:] = bigflat[i-1,:,:]-BIAS
   norm = numpy.median(bigflat[i-1,500:1500,500:1500])
   print('Normalised with the median of the frame :',norm)
   bigflat[i-1,:,:] = bigflat[i-1,:,:]/norm
list.close()

#Calculate flat as median at each pixel
medianflat = numpy.median(bigflat,axis=0)

#Write out result to fitsfile. Import header from one of the raw files
from astropy.io.fits import getheader
hdr = getheader('../raw/ALGc190033.fits',0)
fits.writeto('Flat134.fits',medianflat,hdr,overwrite=True)
