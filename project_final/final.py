'''
created: dec 2020
last modified: dec 2020
owner: bill

bio: astro 533 final
note: 
'''


import illustris_python_true as il
import numpy as np
from astropy.cosmology import WMAP9 as cosmo

# scp final.py billchen@login.rc.fas.harvard.edu:astro533/ 
# scp billchen@login.rc.fas.harvard.edu:astro533/hmf.txt .
# mnt/c/Users/Bill/Documents/Semester 20-21 Fall/ASTRO-533_galaxy/Project_final


def kde_gauss(data, x, h):
    k = sum(np.exp(-((x-data[:,None])/h)**2/2.))
    return k/(np.sqrt(2.*np.pi)*h)

base = '/n/holylfs/LABS/hernquist_lab/IllustrisTNG/Runs/L75n1820TNG/output/' # tng100-1
# base = '/n/holylfs/LABS/hernquist_lab/IllustrisTNG/Runs/L205n2500TNG/output/' # tng300-1

fields = ['GroupMass', 'Group_M_Crit200', 'GroupMassType']
halos = il.groupcat.loadHalos(base, 99, fields=fields)

fields = ['SubhaloMass', 'SubhaloMassType', 'SubhaloStellarPhotometrics', 'SubhaloVelDisp']
subhalos = il.groupcat.loadSubhalos(base, 99, fields=fields)


########## mass #########

# log_bins = np.linspace(6, 16, 100+1)
# bins = 10**log_bins

# hist, bins = np.histogram(halos['GroupMassType'][:,4]*1e10, bins=bins)
# hist_kde = kde_gauss(np.log10(halos['GroupMassType'][:,4])+10, log_bins, h=0.1)

# np.savetxt('hsmf.txt', hist)
# np.savetxt('hsmf_kde.txt', hist_kde)

# hist, bins = np.histogram(subhalos['SubhaloMassType'][:,4]*1e10, bins=bins)
# hist_kde = kde_gauss(np.log10(subhalos['SubhaloMassType'][:,4])+10, log_bins, h=0.1)

# np.savetxt('sub_hsmf.txt', hist)
# np.savetxt('sub_hsmf_kde.txt', hist_kde)


########## photometry ###########

bins = np.linspace(-23, -13, 100+1)

# hist_kde = kde_gauss(subhalos['SubhaloStellarPhotometrics'][:,5]-5*np.log10(0.7), bins, h=0.1)

# np.savetxt('sub_hlf.txt', hist_kde)


# gr = subhalos['SubhaloStellarPhotometrics'][:,4] - subhalos['SubhaloStellarPhotometrics'][:,5]

# idx_red = np.where(gr>0.6)[0]
# idx_blue = np.where(gr<0.6)[0]

# hist_kde_red = kde_gauss(subhalos['SubhaloStellarPhotometrics'][idx_red,5]-5*np.log10(0.7), bins, h=0.1)
# hist_kde_blue = kde_gauss(subhalos['SubhaloStellarPhotometrics'][idx_blue,5]-5*np.log10(0.7), bins, h=0.1)

# np.savetxt('sub_hlf_red_tng300.txt', hist_kde_red)
# np.savetxt('sub_hlf_blue_tng300.txt', hist_kde_blue)


########## velocity dispersion ##########

idx_high_sigma = np.where(subhalos['SubhaloVelDisp']>70)[0]
idx_low_sigma = np.where(subhalos['SubhaloVelDisp']<70)[0]

hist_kde_high_sigma = kde_gauss(subhalos['SubhaloStellarPhotometrics'][idx_high_sigma,5]-5*np.log10(0.7), 
	bins, h=0.1)
hist_kde_low_sigma = kde_gauss(subhalos['SubhaloStellarPhotometrics'][idx_low_sigma,5]-5*np.log10(0.7), 
	bins, h=0.1)

np.savetxt('sub_hlf_high_sigma.txt', hist_kde_high_sigma)
np.savetxt('sub_hlf_low_sigma.txt', hist_kde_low_sigma)

