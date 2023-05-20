# Pyclip
A python script for co-adding fits-images using sign-clipping.

pyclip.py is a python script to co-add a set of fits-files using sigma-clipping. The script is an implementation of the algoritm described in Sect. 2.1 of Møller & Warren 1993 (https://ui.adsabs.harvard.edu/abs/1993A%26A...270...43M/abstract). The algorith is further detailed in this PhD-thesis: https://ui.adsabs.harvard.edu/abs/2000PhDT.........4F/abstract

Here is provided a test-dataset. It is a set of narrow-band images of the quasar GQ1218+0832. Included here are bias-frames, flat-frames and seven raw science frames. The raw data can be reduced using:

mkbias.py
mkflat134.py
reducescience134.py

The script parameters134.py will determine the offsets between the frames, the flux-scaling and the background level of each frame.

pyclip.py is the actual combination script. It takes as input config.txt, files.txt, and params.txt in addtion to the reduced science frames and the trimmed Flat-field frame. 

The script is rather slow. It needs to be optimized.  
