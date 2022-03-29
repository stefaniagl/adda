'''
# This code represents the comparison of scattering intensities calculated
# with GLMT (reference data) and DDA (with ADDA) for the scattering of
# Bessel beam (CS type) by a coated sphere.
# Reference data were kindly provided by Zhuyang Chen
# (DOI:10.1088/2040-8978/16/5/055701).
'''

import os
import matplotlib.pyplot as plt
import bb_module as bb

option = 1 # - scattering by a coated sphere

if not os.path.exists('dda'):
    os.makedirs('dda')
bb.adda_run(option)

# extraction of ADDA results
dda_theta,dda_iper,dda_ipar = bb.extractData('dda',option)
# extraction of reference data (GLMT)
ref_thetaper,ref_iper,ref_thetapar,ref_ipar = bb.extractData('glmt',option)

# data visualisation
fig = plt.figure(figsize=(9,4),constrained_layout=True)
fig.set_constrained_layout_pads(wspace=0.1)

# Parallel scattering intensity (1)
ax = fig.add_subplot(121)
bb.plotData(dda_theta,dda_ipar,ref_thetapar,ref_ipar,1,ax)

#image of a scattering particle
img = plt.imread('particles/shape_image_'+str(option)+'.png')
newax = fig.add_axes([0.12,0.4,0.3,0.3], anchor='NE', zorder=1)
newax.imshow(img)
newax.axis('off')

# Perpendicular scattering intensity (2)
ax = fig.add_subplot(122)
bb.plotData(dda_theta,dda_iper,ref_thetaper,ref_iper,2,ax)

# data save
os.makedirs('saved', exist_ok=True)
plt.savefig('saved/fig13_coatedsphere.pdf', bbox_inches='tight')

plt.show()
