# This is a python program to make a BIAS frame
print('Script running')
nframes = 22
ysize = 2064
xsize = 2048

from astropy.io import fits
import numpy

#Read in the raw bias frames and subtact mean of overscan region 
list = open('bias.list')
bigbias = numpy.zeros((nframes,ysize,xsize),float)
#bigbias = numpy.zeros((nframes,3,3))
for i in range(0,nframes):
   print('Image number:', i)
   filenm = list.readline()
   rawbias = fits.open('../raw/'+filenm.strip("\n"))
   print('Info on file:')
   print(rawbias.info())
   data = numpy.array(rawbias[1].data)
   median = numpy.mean(data[1:2064,1:50])
   data = data - median
   print('Subtracted the median value of the overscan :',median)
   bigbias[i-1,0:ysize-1,0:xsize-1] = data[0:ysize-1,51:51+xsize-1]
list.close()

#Calculate bias is median at each pixel
medianbias = numpy.median(bigbias,axis=0)

#Write out result to fitsfile
from astropy.io.fits import getheader
hdr = getheader('../raw/ALGc200256.fits',0)
fits.writeto('BIAS.fits',medianbias,hdr,overwrite=True)
