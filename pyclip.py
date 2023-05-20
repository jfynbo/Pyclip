#le -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.patches import Ellipse
from matplotlib import axes
import numpy as np
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)
from astropy.utils.data import get_pkg_data_filename
from astropy.io import fits
from astropy.io.fits import getheader
from scipy import optimize
from scipy import stats
from astroscrappy import detect_cosmics
from scipy.stats import sigmaclip

#matplotlib.rcParams['backend'] = 'Qt4Agg' # before 'TkAgg'


# The program should take as input N fits images of the same size and scale,
# taken with the same filter. It should then produce an optimally combined 
# output image using weighted average and sigmaclipping. The program should 
# be a python version of Palle Møller's original clip.f written in Fortran.
#
# Input images:
#   1. Data frames, flat fielded but NOT rebinned to match.
#   2. Flat field used.
#
# Output images:
#   1. Combined data.
#   2. Corresponding variance.
#   3. A "rejection frame" containing in each pixel the number of
#      values rejected at that position.
#
# Output plot (optional):
#   Chi distribution versus surface brightness.
#
#   JF 18.09.01, 04.04.23


# Declarations

axlen = np.arange(2,dtype=int)
start = np.arange(2,dtype=int)
step = np.arange(2,dtype=int)

# Maximum good counts (saturation in ADUs).
badsig=1.e7

# Number of loops for the sigma clipping, maximum number of rejections
# in a single pixel.

nclip=4
maxrej=2

# Read the configuration file
file = open('config.txt', 'r')
# Read info in the top
print('Required info in the config file:')
for i in range(0,4):
    a = file.readline()
fltnm = file.readline().strip("\n")
nimags = int(file.readline())
nall = nimags+5
iarray=1
imname = ["" for i in range(nall)] 
xoff = np.arange(nimags,dtype=int)
yoff = np.arange(nimags,dtype=int)
sky = np.zeros(nimags,dtype=float)
scal = np.zeros(nimags,dtype=float)
scal2 = np.zeros(nimags,dtype=float)
weight = np.zeros(nimags,dtype=float)
weigsq = np.zeros(nimags,dtype=float)
irw = np.arange(nimags)
skip = np.zeros(nimags)
skipx = np.zeros(nimags)
phots = np.zeros(nimags,dtype=float)*0.
vphots = np.zeros(nimags,dtype=float)*0.
ww = np.zeros(nimags,dtype=float)*0.
ww2 = np.zeros(nimags,dtype=float)*0.
xphots = np.zeros(nimags,dtype=float)*0.
xvphot = np.zeros(nimags,dtype=float)*0.
xww = np.zeros(nimags,dtype=float)*0.
xww2 = np.zeros(nimags,dtype=float)*0.
chiarray = np.zeros(nimags,dtype=float)*0.
sigmaarray = np.zeros(nimags,dtype=float)*0.

# Read the config file
a=next(file).split()
chiinp=float(a[0])
csky=float(a[1])
csscal=float(a[2])
a=next(file).split()
ron=float(a[0])
phpadu=float(a[1])
ipix1=float(a[2])
ipix2=float(a[3])
irow1=float(a[4])
irow2=float(a[5])
file.close()

# Read the offsets, sky and scaling parameters
file = open('params.data', 'r')
for nima in range(0,nimags):
   a=next(file).split()
   xoff[nima]=float(a[0])
   yoff[nima]=float(a[1])
   sky[nima]=float(a[2])
   scal[nima]=float(a[3])
file.close()

#    Since we combine images with different seeing, and since we only
#    transform with integral pixels, strong sources will cause large chi
#    values. In an attempt to prevent the rejection of these pixels, we
#    here use the following scheme:
#
#   We define a function chilim(signal). "signal" is here the expected
#   signal in a pixel after sky subtraction (estimated as the weighted
#   mean of non-rejected c measurements). The function chilim() is
#   defined by three parameters; chiinp, areg, and c.
#   chilim(signal) = chiinp                   for    signal < c * sky
#   chilim(signal) = areg*signal/sky + breg   for    signal >= c * sky
#   breg can be found from the other three parameters given that the
#   function chilim() must be continuous in c * sky.
#   breg = chiinp - a * c
#   For the sky value, we pick a typical value of the sky: skytyp
#   All values measured above the function chilim() are rejected.

areg = csscal * 1.
breg = chiinp - areg * csky
skytyp = sky[3] # arbitraary choice
skysig = np.sqrt(skytyp*phpadu + ron*ron)
aregnm = areg/skysig

answer = input('Do you wish to see a chi-plot (will make the combination last longer) (y/n)? ')
if answer == 'n': plot = 0
if answer == 'y':
    plot = 1
    chiupp = float(input('Enter upper value of chi: '))
    if chiupp < 0: chiupp=8 
    photup = float(input('Upper signal (units of sky sigma): '))
    if photup > 0: photup=10
    photup = photup * skysig
    photup = (chiupp - breg) / aregnm
    xmin = -4.*skysig
#    plt.ion()
    fig, ax = plt.subplots()
    x = []
    y = []
    x.append(0)
    y.append(0)
    line1, = ax.plot([xmin,csky*skysig],[chiinp,chiinp],color='r')
    line2, = ax.plot([csky*skysig,photup],[chiinp,aregnm*photup+breg],color='r')
    line3, = ax.plot([x],[y],'go')
    ax.set_xlim([xmin,photup])
    ax.set_ylim([0,chiupp])
    ax.set_xlabel('sigma')
    ax.set_ylabel('chi')

for i in range(0, nimags): scal2[i] = scal[i] * scal[i]

ceilin = csky * skytyp * phpadu
# For now, set ceiling high.
ceilin = 1.e15
phpad2 = phpadu * phpadu

# Calculate weights to give maximum S/N combination for 0 source counts.
for i in range(0, nimags):
    weight[i] = 1. / ((ron * ron + phpadu * sky[i]) * scal2[i])
    weigsq[i] = weight[i] * weight[i]

# Read image names
file = open('files.txt', 'r')
imname[0] = 'rawfiles/Flat'+fltnm+'_trim.fits'
imname[1] = 'comb'+fltnm+'.fits'
imname[2] = 'var'+fltnm+'.fits'
imname[3] = 'sig'+fltnm+'.fits'
imname[4] = 'rej'+fltnm+'.fits'
for i in range(5, nimags+5):
    a = file.readline()
    imname[i]=a.strip("\n")
file.close()

# Size of the images
file = fits.open(imname[0])
flat = file[0].data
hdr = file[0].header
nx = hdr['NAXIS1']
ny = hdr['NAXIS2']

# Open the images
im = np.zeros((nimags,ny,nx),float)
for i in range(0,nimags):
    file = fits.open(imname[i+5])  
    im[i,:,:] = file[0].data[:,:] 

# Calculate the size of the resulting frames.
xmin = 1.e6
xmax =-1.e6
ymin = 1.e6
ymax =-1.e6
for i in range(0, nimags):
    if xoff[i] < xmin: xmin = xoff[i]
    if xoff[i] > xmax: xmax = xoff[i]
    if yoff[i] < ymin: ymin = yoff[i]
    if yoff[i] > ymax: ymax = yoff[i]
      
axlen[0] = ny + round(ymax - ymin)
axlen[1] = nx + round(xmax - xmin)

# Make the arrays for the result frames
sig = np.zeros((nx,nimags),float)
var = np.zeros((nx,nimags),float)
finsig = np.zeros(axlen[1], dtype=float)
finvar = np.zeros(axlen[1], dtype=float)
combout = np.zeros((axlen[0],axlen[1]),float)
varout = np.zeros((axlen[0],axlen[1]),float)
rejout = np.zeros((axlen[0],axlen[1]),int)

#-----Images are now all opened---------------------------------------
#=====================================================================

nreje1 = 0
nreje2 = 0
nreje3 = 0
nreje4 = 0
lorej  = 0
hirej  = 0
nused = 0

#=====================================================================
#-------Here the real work starts-------------------------------------
#     
#     WHAT WE DO...
#     What we do here is first to shift all images to nearest integer
#     pixel, and then combine information in that pixel in all frames. The
#     indivitual pixels can either contain good information, or no
#     information. If it contains no information it is either "skipped", or
#     "rejected". If a pixel is KNOWN a priori not to contain information
#     (hot pixel, bad column, saturated etc.) it has already been marked
#     as such in the input file (set to very large number). Such pixels,
#     and those which are outside of the used area on the CCD, are being
#     "skipped" by the programme. Pixels which are not skipped, but which
#     are found to diverge strongly from the expected value (cosmics etc.)
#     are "rejected".
#     All non-skipped, non-rejected pixel values (and corresponding
#     variances) are kept in the data areas. All skipped and rejected pixels
#     are substituted with the weighted mean of the rest, and the variance
#     is set to 1.e15.
# 
#     Loop over all rows, read the corresponding rows of each image
#     (using the y-offsets) and the row of the flat field.
#     irow in the final image corresponds to irw[i] = irow - nint[yoff(i)]
#     of individual images.
# 
#=====================================================================

# This for loop ends at mark 210
for irow in range(0,axlen[0]):
    print('Working on row:',irow)
    for i in range(0,nimags):
            irw[i] = irow - round(yoff[i] - ymin)
# Read rows, read corresponding flatfield row, calculate variance.
    for i in range(0,nimags):
        if irw[i] < irow1 or irw[i] > irow2: 
            skip[i] = 1 
        if irw[i] >= irow1 and irw[i] <= irow2: 
            ffrow = flat[irw[i],:]
            sig[:,i] = im[i,irw[i],:]
# Data has already been flat fielded, just subtract the sky and scale.
# Also put everything on real photon counts, and calculate the
# corresponding variance of each pixel.
            for ipix in range(0,nx):
               tsig = sig[ipix,i]
               if tsig > badsig: 
                   tsig = badsig
#               if np.isnan(tsig): tsig = badsig
               if ffrow[ipix] <= 0.: ffrow[ipix] = 1.e-6
               var[ipix,i] = scal2[i]*phpadu*abs(tsig/ffrow[ipix])
               var[ipix,i] = var[ipix,i]+(scal[i]*ron/ffrow[ipix])**2
#               if var[ipix,i] == 0.: 
#                  sig[ipix,i] = badsig+1
#                  var[ipix,i] = 1.e13+1.
               sig[ipix,i] = (sig[ipix,i]-sky[i])*scal[i]*phpadu
               skip[i] = 0
               if np.isnan(sig[ipix,i]): 
                  sig[ipix,i] = badsig+1
                  var[ipix,i] = 1.e13+1.

# Now do the comparison. Do that by applying the x-offsets when getting
# the values. Skip the ones that should be skipped, and count how many
# to use.
    
    for ipix in range(0, axlen[1]): 
        chilim = chiinp
        nrejec = 0
        lowrej = 0
        higrej = 0
        nskipx = 0.
        nnext = 1
#
        for i in range(0, nimags):
            skipx[i] = 0
# First do the "skipping" part.
            if skip[i] == 1: 
               skipx[i] = 1
            jpix = ipix - round(xoff[i] - xmin)
            if jpix <= ipix1 or jpix >= ipix2:
               skipx[i] = 1
               nskipx = nskipx + 1

            if skipx[i] == 0 and jpix >= ipix1 and jpix <= ipix2:
               phots[nnext-1]  = sig[jpix, i]
               vphots[nnext-1] = var[jpix, i]
               ww[nnext-1]     = weight[i]
               ww2[nnext-1]    = weigsq[i]
               nnext         = nnext + 1
               nused = nnext - 1 
               nrejec = nimags - nused
    
# The above initial list of non-skipped, non-rejected pixels shall be
# needed again later, hence remember it.

        nxused = nused

        if nused > 0:
            for i in range(0, nused):
                xphots[i]  = phots[i]
                xvphot[i]  = vphots[i]                                                                                           
                xww[i]     = ww[i]
                xww2[i]    = ww2[i]

            nused = nxused

# This is where the clipping should be done. First find a preliminary estimate of the value of this pixel.

        done = False
        while done == False:
# 0 acceptable pixels:
           if nused == 0:
               finsig[ipix] = 1.e9 * phpadu
               finvar[ipix] = 1.e9 * phpad2
               done = True
               break
           if nused == 1:
# 1 acceptable pixel
               finsig[ipix] = phots[nused-1]
               finvar[ipix] = vphots[nused-1]
               done = True
               break
           if nused == 2:
# 2 acceptable pixels: Calculate weighted mean, and chi of fit.
               wsum = ww[0] + ww[1]
               sigmn = (ww[0] * phots[0] + ww[1] * phots[1])/wsum
               varmn = (ww2[0]*vphots[0] + ww2[1]*vphots[1])/(wsum*wsum)
               chi1sq = ((sigmn - phots[0])**2) / vphots[0]
               chi2sq = ((sigmn - phots[1])**2) / vphots[1]
               chi = np.sqrt(chi1sq + chi2sq)
# If both are very big, dont care about the chi just accept it.
               if phots[0] > ceilin and phots[1] > ceilin:
                   finsig[ipix] = sigmn
                   finvar[ipix] = varmn
                   done = True
                   break
               if chi < chilim:
                   finsig[ipix] = sigmn
                   finvar[ipix] = varmn
                   done = True
                   break
               if chi > chilim:
# Reject the biggest one.
                   if phots[0] > phots[1]:
                         phots[0]  = phots[1]
                         vphots[0] = vphots[1]
                         ww[0]     = ww[1]
                         ww2[0]    = ww2[1]
                   nused = 1
                   higrej += 1
                   nrejec += 1
           if nused > 2:
# 3 or more acceptable pixels: Calculate weighted mean, and chi of fit.
               wsum   = 0.
               sigsum = 0.
               varsum = 0.
               for i in range(0, nused):
                  wsum   = wsum + ww[i]
                  sigsum = sigsum + ww[i] * phots[i]
                  varsum = varsum + ww2[i] * vphots[i]
               sigmn  = sigsum / wsum
               varmn  = varsum / (wsum * wsum)
# If all values are very big, dont care about the chi just accept it.
               allhigh = True
               for i in range (0, nused):
                   if (phots[i] < ceilin): allhigh = False
               if allhigh == True:
                   finsig[ipix] = sigmn
                   finvar[ipix] = varmn
                   nrejec = maxrej + 1

               chisum = 0.
               for i in range(0,nused): 
                   chisum = chisum + ((sigmn - phots[i])**2) / vphots[i]
               chi = np.sqrt(chisum / float(nused - 1))

# Here we need to handle Palles jump 231. Search for an interrim value
# Accept this if chi < given limit.

               interrim = False
               while interrim == False:
                   chilim = max([chiinp, finsig[ipix]*aregnm + breg])
                   if chi < chilim:
                       finsig[ipix] = sigmn
                       finvar[ipix] = varmn
                       interrim = True
                       break
# Reject the highest and the lowest, one by one, and calculate the
# chi of the remaining pixels. Choose the one with the smallest chi.

#-----------------------Reject highest value----------------------------
                   phtmax = -1.e20
                   jmax = 0
                   for j in range(0, nused):
                       if phots[j] > phtmax:
                         phtmax = phots[j]
                         jmax   = j
                   wsum   = 0.
                   sigsum = 0.
                   varsum = 0.
                   for i in range(0, nused):
                       if i != jmax:
                           wsum   = wsum + ww[i]
                           sigsum = sigsum + ww[i] * phots[i]
                           varsum = varsum + ww2[i] * vphots[i]
                   if wsum == 0: wsum = 1
                   sigmax  = sigsum / wsum
                   varmax  = varsum / (wsum * wsum)
                   chisum = 0.
                   for i in range(0, nused):
                       if i != jmax:
                           chisum = chisum + ((sigmax - phots[i])**2) / vphots[i]
                   chimax = np.sqrt(chisum / float(nused - 2))

#-----------------------Reject lowest value----------------------------
                   phtmin = 1.e20
                   jmin = 0
                   for j in range(0, nused):
                       if phots[j] < phtmin:
                         phtmin = phots[j]
                         jmin   = j
                   wsum   = 0.
                   sigsum = 0.
                   varsum = 0.
                   for i in range(0, nused):
                       if i != jmin:
                           wsum   = wsum + ww[i]
                           sigsum = sigsum + ww[i] * phots[i]
                           varsum = varsum + ww2[i] * vphots[i]
                   if wsum == 0: wsum = 1
                   sigmin  = sigsum / wsum
                   varmin  = varsum / (wsum * wsum)
                   chisum = 0.
                   for i in range(0, nused):
                       if i != jmax:
                           chisum = chisum + ((sigmin - phots[i])**2) / vphots[i]
                   chimin = np.sqrt(chisum / float(nused - 2))

# Now determine which is best
                   if chimin < chimax:
                       sigmn  = sigmin
                       varmn  = varmin
                       chi    = chimin
                       jrej   = jmax
                       higrej += 1
                   if chimin > chimax:
                       sigmn  = sigmax
                       varmn  = varmax
                       chi    = chimax
                       jrej   = jmin
                       lowrej += 1
                   nrejec = nrejec + 1

# Accept this if chi < given limit.

                   if chi <= chilim:
                       finsig[ipix] = sigmn
                       finvar[ipix] = varmn
                       interrim = True
                       break
                   if jrej != nused-1:
                       phots[jrej]  = phots[nused-1]
                       vphots[jrej] = vphots[nused-1]
                       ww[jrej]     = ww[nused-1]
                       ww2[jrej]    = ww2[nused-1]                                                                                                 
                   interrim = True

## A good guess for the real value is found, now use sigma clipping.
## Do the clipping in a loop. Clip nclip times.
## As before apply the x-offsets when getting the values. Skip the ones
## that should be skipped, and count how many to use.
## this loop ends at 400

               #This loop ends at Mark 400
               for iclip in range(0, nclip):
                    chilim = max([chiinp, finsig[ipix]*aregnm + breg])
                    nrejec = 0
                    lowrej = 0
                    highrej = 0

# Go back to original list of pixel values.

                    nused = nxused
                    nrejec = nimags-nxused

                    for i in range(0, nused):
                        phots[i]  = xphots[i]
                        vphots[i] = xvphot[i]
                        ww[i]     = xww[i]
                        ww2[i]    = xww2[i]

#   If a plot was requested, and this is the first run, plot the chi
#     distribution.
                    if plot == 1 and iclip == 0: 
                        for i in range(0,nused):
                           chi=abs((finsig[ipix]-phots[i])/np.sqrt(vphots[i]))
                           sigmaarray[iarray-1] = finsig[ipix]
                           chiarray[iarray-1] = chi
                           iarray = iarray + 1
                        line3.set_data(sigmaarray,chiarray)
                        ax.scatter(sigmaarray,chiarray)
                        plt.pause(0.00001)
                        fig.canvas.update()
                        fig.canvas.draw()
                        iarray = 1

#     First find the chi of all individual values, reject them if they fall
#     outside chilim.


                    done2 = False
                    while done2 == False:
                        ngood = 0
                        for i in range(0, nused):
                            chi= abs((finsig[ipix] - phots[i])/np.sqrt(vphots[i]))
                            if chi > chilim:
                                 if finsig[ipix] < phots[i]: lowrej += 1
                                 if finsig[ipix] >= phots[i]: higrej += 1
                                 phots[i] = phots[nused-1]
                                 vphots[i] = vphots[nused-1]
                                 ww[i]     = ww[nused-1]
                                 ww2[i]    = ww2[nused-1]
                                 nrejec    = nrejec + 1
                                 nused     = nused-1
                            if chi < chilim: ngood = ngood + 1
                        if ngood == nused and nused > 0: 
                           done = True
                           break
                        if nused <= 0:
                           nrejec = nimags
                           finsig[ipix] = 1.e9 * phpadu
                           finvar[ipix] = 1.e9 * phpad2
                           done = True
                           break

# Now calculate the weighted mean of the rest.

                    if nrejec < maxrej:
                       wsum = 0.
                       sigsum = 0.
                       varsum = 0.
                       for i in range(0,nused):
                           wsum = wsum + ww[i]
                           sigsum = sigsum + ww[i]*phots[i]
                           varsum = varsum + ww2[i]*vphots[i]
                       finsig[ipix] = sigsum / wsum
                       finvar[ipix] = varsum / (wsum * wsum)
   
# Mark 400
                    done = True

# Here the clipping should be complete


        rejout[irow,ipix] = nrejec
        combout[irow,ipix] = finsig[ipix] / phpadu
        varout[irow,ipix] = finvar[ipix] / phpad2
        if nused == 0:
            rejout[irow,ipix] = nimags 
            combout[irow,ipix] = 0.
            varout[irow,ipix] = 0.

# 230
# 210

# That should then be it for this row foevks. Now write it into final
# signal and variance image.


file = fits.open(imname[5])
hdr = file[0].header

#Set very high numbers to 0
highpix = np.where(combout > 1.e7)
combout[highpix] = 0.

fits.writeto(imname[1],combout,hdr,overwrite=True)
fits.writeto(imname[2],varout,hdr,overwrite=True)
fits.writeto(imname[3],np.sqrt(varout),hdr,overwrite=True)
fits.writeto(imname[4],rejout,hdr,overwrite=True)
   
