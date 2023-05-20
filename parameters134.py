# -*- coding: utf-8 -*-
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
from astroscrappy import detect_cosmics
from astropy.stats import SigmaClip
from photutils.background import ModeEstimatorBackground
from photutils import CircularAperture
from photutils import CircularAnnulus
from photutils import aperture_photometry
import astroscrappy


#Input images and parameters
nframes = 7

gain = 0.16 
ron = 4.3 
frac = 0.01
objlim = 15
sigclip = 5
niter = 2

vmin = 100
vmax = 10000

aperture_radius = 8.0
r_in = 15
r_out = 20

Imagelist = open('files.txt')

#Parameter arrays
names = []
xoff = []
yoff = []
sky = []
counts = []

#Definitions
fws = 20 #fit-window-size

def gaussian(height, center_x, center_y, width_x, width_y):
    """Returns a gaussian function with the given parameters"""
    width_x = float(width_x)
    width_y = float(width_y)
    return lambda x,y: height*np.exp(
                -(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)

def moments(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution by calculating its
    moments """
    total = data.sum()
    X, Y = np.indices(data.shape)
    x = (X*data).sum()/total
    y = (Y*data).sum()/total
    col = data[:, int(y)]
    width_x = np.sqrt(np.abs((np.arange(col.size)-x)**2*col).sum()/col.sum())
    row = data[int(x), :]
    width_y = np.sqrt(np.abs((np.arange(row.size)-y)**2*row).sum()/row.sum())
    height = data.max()
    return height, x, y, width_x, width_y

def fitgaussian(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution found by a fit"""
    params = moments(data)
    errorfunction = lambda p: np.ravel(gaussian(*p)(*np.indices(data.shape)) -
                                 data)
    p, success = optimize.leastsq(errorfunction, params)
    return p

#Loop over list of images

for i in range(0,nframes):
#Load image
   image = Imagelist.readline()
   names.append(image)
   print('Image:',image)
   image_file = get_pkg_data_filename(image.strip("\n"))
   fits.info(image_file)
   image_data = fits.getdata(image_file, ext=0)
   crmask, clean_arr = astroscrappy.detect_cosmics(image_data, sigclip=sigclip, sigfrac=frac, objlim=objlim, cleantype='medmask', niter=niter, sepmed=True, verbose=True)
   image_data = clean_arr 

##############################################################################
#Background
   sigma_clip = SigmaClip(sigma=3.0)
   bkg = ModeEstimatorBackground(median_factor=3.0, mean_factor=2.0,
                                sigma_clip=sigma_clip) 
   bkg_value = bkg.calc_background(image_data[150:1700,150:1700])
   sky.append(bkg_value)
   print(image,bkg_value)

##############################################################################
#Select star

   plt.figure()
   plt.imshow(image_data, cmap='gray_r', vmin=vmin, vmax=vmax, origin='lower')
   plt.colorbar()

   print('Click on star with the cursor. End with q')
   tpoints = plt.ginput(n=1, timeout=30, show_clicks=True, mouse_add=1, mouse_stop=2)
   xstar = int(tpoints[0][0])
   ystar = int(tpoints[0][1])
   plt.show()

#Now fit a 2d-gaussian to a region with size fws x fws around this position
   subimage = image_data[ystar-fws:ystar+fws,xstar-fws:xstar+fws]
   plt.matshow(subimage, cmap=plt.cm.gist_earth_r, origin='lower')
   params = fitgaussian(subimage)
   fit = gaussian(*params)
   plt.contour(fit(*np.indices(subimage.shape)), cmap=plt.cm.copper)
   ax = plt.gca()
   (height, x, y, width_x, width_y) = params

   plt.text(0.95, 0.05, """
   x : %.1f
   y : %.1f
   width_x : %.1f
   width_y : %.1f""" %(y+xstar-fws, x+ystar-fws, width_x, width_y),
        fontsize=16, horizontalalignment='right',
        verticalalignment='bottom', transform=ax.transAxes)
   plt.show()

   xoff.append(y+xstar-fws)
   yoff.append(x+ystar-fws)

   aperture = CircularAperture((y+xstar-fws,x+ystar-fws), r=aperture_radius)
   annulus = CircularAnnulus((y+xstar-fws,x+ystar-fws), r_in=r_in, r_out=r_out)
   apers = [aperture, annulus]
   phot_table_local_bkg = aperture_photometry(image_data, apers)
   bkg_mean = phot_table_local_bkg['aperture_sum_1'] / annulus.area
   bkg_sum = bkg_mean * aperture.area
   final_sum = phot_table_local_bkg['aperture_sum_0'] - bkg_sum
   phot_table_local_bkg['residual_aperture_sum'] = final_sum
   counts.append(final_sum)

   plt.figure()
   plt.imshow(image_data, cmap='gray_r', vmin=vmin, vmax=vmax, origin='lower')
   plt.ylim(x+ystar-fws-100,x+ystar-fws+100)
   plt.xlim(y+xstar-fws-100,y+xstar-fws+100)
   ap_patches = aperture.plot(color='white', lw=2,
                           label='Photometry aperture')
   ann_patches = annulus.plot(color='red', lw=2,
                                    label='Background annulus')
   handles = (ap_patches[0], ann_patches[0])
   plt.legend(loc=(0.17, 0.05), facecolor='#458989', labelcolor='white',
           handles=handles, prop={'weight': 'bold', 'size': 11})

   plt.colorbar()
   plt.show()

#Scaling from aperture photometry (important that the star is not saturated)

print(counts,  max(counts))

#Now work out offsets
for n in range(0,nframes):
   shiftx = int(np.mean(xoff)-xoff[n])
   shifty = int(np.mean(yoff)-yoff[n])
   scale = counts[n] / max(counts) 
   print(names[n],shiftx,shifty, sky[n], scale)

