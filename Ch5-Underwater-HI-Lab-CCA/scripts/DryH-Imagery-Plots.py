# -*- coding: utf-8 -*-
"""
Created on Wed Apr 26 18:38:40 2023

@author: jcmontes

UHI - Tank

"""

from pathlib import Path
import numpy as np

import spectral as sp
import scipy.signal as ss
import pysptools.spectro as spectro

# import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from cycler import cycler

# Files path
files_path = Path(r"C:\Users\jcmontes\OneDrive - University of Tasmania\01_Projects_Drive\UHI_Tank\data\JCMH-FC-Reflectance")
out_p = r"C:\Users\jcmontes\OneDrive - University of Tasmania\01_Projects_Drive\UHI_Tank\results\UHI-Imagery\Dry-Experiment-TwoIndex/"

#Files Live-Epibiota-Shells
f1_path = files_path / 'T4-Dry-Epibiota-Shells'

# Reflectance images
hdr_list = list(f1_path.glob('*.hdr'))
bin_list = list(f1_path.glob('*.img'))

# List of files
targets = list(zip(hdr_list,bin_list))

# Several tests to get this dimensions. Can be fixed with improvements to Translation stage
sample_dims = [
[(50,190),(90,460)],
[(55,205),(70,430)],
[(35,150),(40,420)],
[(15,110),(60,430)],
[(15,120),(60,450)],
[(40,145),(30,450)],
[(45,155),(25,440)],
]


sample_key = [] # Empty list for sample names (strings)
sample_dict = {} # Empty dict for Original Reflectance data cubes (np.ndarray)
dd_dict = {}
kmeans_dict = {} # Empty dict for k-means results
map_dict = {}

#K mean number of clusters & colors
# n_clusters = 5

## Color maps
# Set the colormap for plt.imshow
cmap = plt.cm.Dark2
# Set the color cycle for plt.plot
# colors = cmap(np.linspace(0, 1, n_clusters))
# plt.rc('axes', prop_cycle=(cycler('color', colors)))

norm = Normalize(vmin=0.03, vmax=0.1, clip=True)


# Begin
for i, sample in enumerate(targets):
    sample_key.append(sample[0].stem[20:-12])
    print(sample_key[i])
    target_open = sp.envi.open(sample[0], sample[1])
    target_data = target_open.load()
    
    # Build dictionaries - Reflectance and Derivatives
    sample_dict[sample_key[i]] = ss.savgol_filter(target_data[:,:,3:238], 3, 1, axis=2).copy()
    dd_dict[sample_key[i]] = ss.savgol_filter(target_data[:,:,3:238], 15, 3, axis=2, deriv=2).copy()
    
    #Adjust the arrays to specific dimensions with only the targets 
    sample_dict[sample_key[i]] = sample_dict[sample_key[i]][sample_dims[i][0][0]:sample_dims[i][0][1], 
                                                            sample_dims[i][1][0]:sample_dims[i][1][1], :].copy()
    # dd_dict[sample_key[i]] = dd_dict[sample_key[i]][sample_dims[i][0][0]:sample_dims[i][0][1], sample_dims[i][1][0]:sample_dims[i][1][1], :].copy()
    map_dict[sample_key[i]] = np.empty((sample_dict[sample_key[i]].shape[0], sample_dict[sample_key[i]].shape[1], 1), dtype=np.float64)
    
    # Masking
    mask = sample_dict[sample_key[i]][:,:,70] < 0.01 # 500 nm threshold
    sample_m = sample_dict[sample_key[i]].copy() # Copy array
    sample_m[mask] = 0 # Apply mask
    sample_dict[sample_key[i]] = sample_m.copy() # Save
    
    # sample_mdd = dd_dict[sample_key[i]].copy()
    # sample_mdd[mask] = 0.01 # Apply mask on dd array
    # dd_dict[sample_key[i]] = sample_mdd.copy()
    
    ## K means
    # (m,c) = sp.kmeans(sample_dict[sample_key[i]],n_clusters,8)
    # kmeans_dict[sample_key[i]] = (m.copy(),c.copy()) # Store k-means results in the kmeans_dict dictionary
    
print('Done :)')

# Save wavelengths
wvl = [float(i) for i in target_open.metadata['wavelength'][3:238]]

i = 0

while i < len(targets):
    
    #RGB
    plt.imshow((sample_dict[sample_key[i]][:, :, (145,92,12)]**0.45), interpolation='spline36')
    plt.title(sample_key[i] + ' RGB - 647, 555, 419 nm', fontsize=12)
    plt.tight_layout()
    plt.savefig(str(out_p) + sample_key[i] + '-Dry-RGB.png', dpi=600)
    plt.show()
    plt.close()
    
    #NDVI
    # plt.imshow((sample_dict[sample_key[i]][:, :, 179] - sample_dict[sample_key[i]][:, :, 159]) /
    #                      (sample_dict[sample_key[i]][:, :, 179] + sample_dict[sample_key[i]][:, :, 159]), vmin=0.1, vmax=1, interpolation='spline36')
    # plt.title(sample_key[i] + ' NDVI 700 nm', fontsize=12)
    # # cbar = plt.colorbar(shrink=0.5)
    # # cbar.ax.get_yaxis().labelpad = 12
    # plt.tight_layout()
    # plt.savefig(str(out_p) + sample_key[i] + '-Dry-NDVI.png', dpi=600)
    # plt.show()
    # plt.close()
    
    # ## dd562.6
    # plt.imshow((dd_dict[sample_key[i]][:, :, 96]), 
    #                      norm=mpl.colors.Normalize(vmin=0.0002, vmax=dd_dict[sample_key[i]][:, :, 96].max()), interpolation='spline36')
    # plt.title(sample_key[i] + ' δδ 563 nm', fontsize=12)
    # # cbar = plt.colorbar(shrink=0.5)
    # # cbar.ax.get_yaxis().labelpad = 12
    # plt.tight_layout()
    # plt.savefig(str(out_p) + sample_key[i] + '-Dry-dd563.png', dpi=600)
    # plt.show()
    # plt.close()
    
    # ## K-Means
    # plt.imshow(kmeans_dict[sample_key[i]][0], cmap=cmap)
    # plt.title(sample_key[i] + ' K-Mean clusters', fontsize=12)
    # plt.tight_layout()
    # plt.savefig(str(out_p) + sample_key[i] + '-Dry-k-means-img.png', dpi=600)
    # plt.show()
    # plt.close()
    
    # fig, ax = plt.subplots()
    # ax.set_ylim([0,0.8])
    # for j in range(kmeans_dict[sample_key[i]][1].shape[0]):
    #     ax.plot(np.linspace(399,800,235), kmeans_dict[sample_key[i]][1][j,:])
    # plt.title(sample_key[i] + ' K-Mean clusters spectral signatures', fontsize=12)
    # plt.tight_layout()
    # plt.savefig(str(out_p) + sample_key[i] + '-Dry-K-means-spectra.png', dpi=600)
    # plt.show()
    
    ## dd563
    # pigment_map = np.empty((sample_dict[sample_key[i]].shape[0], sample_dict[sample_key[i]].shape[1], 1), dtype=np.float64)
    
    # for y in range(0, sample_dict[sample_key[i]].shape[0]):
    #     for x in range(0, sample_dict[sample_key[i]].shape[1]):
    #         if sample_dict[sample_key[i]][y,x].any() == 0:
    #             pigment_map[y,x] = 0
    #         else:
    #             second_d = ss.savgol_filter(sample_dict[sample_key[i]][y,x], 15, 3, deriv=2)
    #             pigment_map[y,x] = 0.02675 + (57.24781 * second_d[96]) # dd 563 nm
                
    # dd563_map = pigment_map.copy()
    
    # plt.imshow(dd563_map, interpolation=None, norm=norm, cmap='magma')
    # plt.title(sample_key[i] + 'R-Phycoerythrin - δδ 563 nm', fontsize=12)
    # # cbar = plt.colorbar(shrink=0.5)
    # # cbar.ax.get_yaxis().labelpad = 18
    # # cbar.ax.set_ylabel('$\mathregular{mg/mm^{2}}$', rotation=270, fontsize=12)
    # plt.tight_layout()
    # plt.savefig(str(out_p) + sample_key[i] + '-dd563_dry_map.png', dpi=600)
    # plt.show()
    # plt.close()
    
    ## ANMB 563_40
    pigment_map = np.empty((sample_dict[sample_key[i]].shape[0], sample_dict[sample_key[i]].shape[1], 1), dtype=np.float64)
    
    for y in range(0, sample_dict[sample_key[i]].shape[0]):
        for x in range(0, sample_dict[sample_key[i]].shape[1]):
            if sample_dict[sample_key[i]][y,x].any() == 0:
                pigment_map[y,x] = 0
            else:
                # sav_gol_filter = ss.savgol_filter(sample[y,x][88:113], 3, 1)
                region = sample_dict[sample_key[i]][y,x][85:108]
                # pre_continuum = sav_gol_filter.tolist()
                pre_continuum = region.tolist()
                bands = list(wvl[85:108]) # Select bands (PE563, 88:113), CHL [12:36]
                conv_h_q = spectro.FeaturesConvexHullQuotient(pre_continuum, bands, baseline=0.9999)
                auc = conv_h_q.get_area(feat_no=1)
                mbd = conv_h_q.get_absorbtion_depth(feat_no=1)
                anmb = auc/mbd
                
                pigment_map[y,x] = 0.03139 + (0.00678 * anmb) # ANMB 563 nm 
    
    anmb563_40map = pigment_map.copy()
    
    plt.imshow(anmb563_40map, interpolation=None, norm=norm, cmap='magma')
    plt.title(sample_key[i] + ' R-Phycoerythrin - ANMBD [543 to 583 nm]', fontsize=12)
    # cbar = plt.colorbar(shrink=0.5)
    # cbar.ax.get_yaxis().labelpad = 16
    # cbar.ax.set_ylabel('$\mathregular{mg/mm^{2}}$', rotation=270, fontsize=12)
    plt.tight_layout()
    plt.savefig(str(out_p) + sample_key[i] + '-anmb563_40_dry_map.png', dpi=600)
    plt.show()
    plt.close()
    
    i += 1
