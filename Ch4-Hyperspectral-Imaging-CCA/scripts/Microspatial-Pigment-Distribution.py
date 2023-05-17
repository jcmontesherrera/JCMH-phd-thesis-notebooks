# -*- coding: utf-8 -*-
"""
Created on Sun Jan 29 15:50:30 2023

@author: jcmontes
"""

from pathlib import Path
import numpy as np

import spectral as sp
import pysptools.spectro as spectro
import scipy.signal as ss

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
# import matplotlib.patches as patches
from matplotlib.colors import Normalize

plt.rcParams['figure.figsize'] = [20, 9]

# Files path
files_path = Path(r"C:\Users\jcmontes\OneDrive - University of Tasmania\01_Projects_Drive\Imaging_spectroscopy\Phenotyping_macroalgae\data\NIWA-Antarctic-CCA\Kestrel-Reflectance")

# Reflectance images
hdr_list = list(files_path.glob('*.hdr'))
bin_list = list(files_path.glob('*.img'))

# List of files
targets = list(zip(hdr_list,bin_list))

# Several tests to get this dimensions. Can be fixed with improvements to Translation stage
sample_dims = [
    [(5,150), (5,265)],  # ['CE-01_']
    [(20,145), (25,220)], # ['CE-02_']
    [(190,300), (40,250)],  # ['CE-03_']
    [(190,300), (52,223)], # ['CE-04_']
    [(103,195), (66,245)], # ['CE-05_']
    [(106,180), (55,170)], # ['GHN-01']
    [(95,165), (100,203)], # 02
    [(130,225), (65,215)], # 03
    [(37,127), (65,235)], # 04
    [(73,177), (37,187)], # 05
    [(155,225), (52,215)], # 06
    [(75,170), (75,210)], # 07
    [(85,140), (98,205)], # 08
    [(108,192), (55,220)], # GHS01
    [(290,445), (60,245)], # GHS02
    [(180,320), (40,220)], # GHS03
    [(0,90), (40,210)], # WHite1
    [(0,90), (40,210)] # White2
   ]

# Sample Dictionary
sample_key = [] # Empty list for sample names (strings)
sample_dict = {} # Empty dict for Original Reflectance data cubes (np.ndarray)
map_dict = {} # Pigment map dictionary (2D arrays with same Y & X dims as samples) - The canvas

# Begin
for i, sample in enumerate(targets):
    sample_key.append(sample[0].stem[9:15])
    # print(sample_key[i])
    target_open = sp.envi.open(sample[0], sample[1])
    target_data = target_open.load()
    
    # Build dictionaries - Reflectance and Derivatives
    sample_dict[sample_key[i]] = ss.savgol_filter(target_data[100:500,80:360,3:238], 3, 1, axis=2).copy() # Read Target Image w/ from band 3 to 238
    # dd_dict[sample_key[i]] = ss.savgol_filter(target_data[100:500,80:360,:], 15, 3, axis=2, deriv=2).copy()
    
    # Adjust the arrays to specific dimensions with only the targets 
    sample_dict[sample_key[i]] = sample_dict[sample_key[i]][sample_dims[i][0][0]:sample_dims[i][0][1], 
                                                            sample_dims[i][1][0]:sample_dims[i][1][1], :].copy()
    # dd_dict[sample_key[i]] = dd_dict[sample_key[i]][sample_dims[i][0][0]:sample_dims[i][0][1], sample_dims[i][1][0]:sample_dims[i][1][1], :].copy()
    map_dict[sample_key[i]] = np.empty((sample_dict[sample_key[i]].shape[0], sample_dict[sample_key[i]].shape[1], 1), dtype=np.float64)
    
    # Masking
    mask = sample_dict[sample_key[i]][:,:,234] < 0.07 # 800 nm threshold
    sample_m = sample_dict[sample_key[i]].copy() # Copy array
    sample_m[mask] = 0 # Apply mask
    sample_dict[sample_key[i]] = sample_m.copy() # Save

# Save wavelengths
wvl = [float(i) for i in target_open.metadata['wavelength'][3:238]]

## ~~~~~~ Figure ~~~~~~~~~~~ ##
# Set color range
norm = Normalize(vmin=0.03, vmax=0.1, clip=True)
norm_chl = Normalize(vmin=10, vmax=20, clip=True)

i = 0

while i < len(targets):

    #Output path
    out_p = r"C:\Users\jcmontes\OneDrive - University of Tasmania\01_Projects_Drive\Imaging_spectroscopy\Phenotyping_macroalgae\results\Hyperspectral_Imaging\Sample-Pigment-Maps/"
    
    # Read sample
    sample = sample_dict[sample_key[i]]
    
    # Display RGB image + save
    plt.imshow(sample[:,:,(175,92,12)])
    plt.title(sample_key[i] + ' RGB - Red (698 nm) Green (555 nm) Blue (419 nm)', fontsize=15)
    plt.tight_layout()
    # plt.savefig(out_p + sample_key[i] + '-RGB.png', dpi=600)
    plt.show()
    plt.close()
    
    # #d563 ~~~~~~~~~~~~~~~~~~~~~
    # pigment_map = np.empty((sample.shape[0], sample.shape[1], 1), dtype=np.float64)
    
    # for y in range(0, sample.shape[0]):
    #     for x in range(0, sample.shape[1]):
    #         if sample[y,x].any() == 0:
    #             pigment_map[y,x] = 0
    #         else:
    #             first_d = ss.savgol_filter(sample[y,x], 30, 1, deriv=1)
    #             pigment_map[y,x] = 0.04088 + (33.33275 * first_d[96]) # d 563 nm
    
    # d563_map = pigment_map.copy()
    
    # plt.imshow(d563_map, interpolation='spline16', norm=norm, cmap='magma')
    # plt.title(sample_key[i] + ' R-Phycoerythrin - δ 563 nm', fontsize=15)
    # cbar = plt.colorbar()
    # cbar.ax.get_yaxis().labelpad = 16
    # cbar.ax.set_ylabel('$\mathregular{mg/mm^{2}}$', rotation=270, fontsize=16)
    # plt.tight_layout()
    # # plt.savefig(out_p + sample_key[i] + '-d563_map.png', dpi=600)
    # plt.show()
    # plt.close()

    # #dd563 ~~~~~~~~~~~~~~~~~~~~~
    # pigment_map = np.empty((sample.shape[0], sample.shape[1], 1), dtype=np.float64)
    
    # for y in range(0, sample.shape[0]):
    #     for x in range(0, sample.shape[1]):
    #         if sample[y,x].any() == 0:
    #             pigment_map[y,x] = 0
    #         else:
    #             second_d = ss.savgol_filter(sample[y,x], 15, 3, deriv=2)
    #             pigment_map[y,x] = 0.02675 + (57.24781 * second_d[96]) # dd 563 nm
                
    # dd563_map = pigment_map.copy()
    
    # plt.imshow(dd563_map, interpolation=None, norm=norm, cmap='magma')
    # plt.title(sample_key[i] + 'R-Phycoerythrin - δδ 563 nm', fontsize=15)
    # cbar = plt.colorbar()
    # cbar.ax.get_yaxis().labelpad = 18
    # cbar.ax.set_ylabel('$\mathregular{mg/mm^{2}}$', rotation=270, fontsize=16)
    # plt.tight_layout()
    # # plt.savefig(out_p + sample_key[i] + '-dd563_map_nointerp.png', dpi=600)
    # plt.show()
    # plt.close()

    #dd569 ~~~~~~~~~~~~~~~~~~~~~
    pigment_map = np.empty((sample.shape[0], sample.shape[1], 1), dtype=np.float64)
    
    for y in range(0, sample.shape[0]):
        for x in range(0, sample.shape[1]):
            if sample[y,x].any() == 0:
                pigment_map[y,x] = 0
            else:
                second_d = ss.savgol_filter(sample[y,x], 15, 3, deriv=2)
                pigment_map[y,x] = 0.02038 + (37.94613689 * second_d[100]) # dd 569 nm
                
    dd569_map = pigment_map.copy()
    
    plt.imshow(dd569_map, interpolation=None, norm=norm, cmap='magma')
    plt.title(sample_key[i] + 'R-Phycoerythrin - δδ 569 nm', fontsize=15)
    cbar = plt.colorbar()
    cbar.ax.get_yaxis().labelpad = 16
    cbar.ax.set_ylabel('$\mathregular{mg/mm^{2}}$', rotation=270, fontsize=16)
    plt.tight_layout()
    # plt.savefig(out_p + sample_key[i] + '-dd569_map_nointerp.png', dpi=600)
    plt.show()
    plt.close()
    
    ## ANMB 563_20 ~~~~~~~~~~~~~~~~~~~~~
    pigment_map = np.empty((sample.shape[0], sample.shape[1], 1), dtype=np.float64)
    
    for y in range(0, sample.shape[0]):
        for x in range(0, sample.shape[1]):
            if sample[y,x].any() == 0:
                pigment_map[y,x] = 0
            else:
                # sav_gol_filter = ss.savgol_filter(sample[y,x][88:113], 3, 1)
                region = sample[y,x][90:103]
                # pre_continuum = sav_gol_filter.tolist()
                pre_continuum = region.tolist()
                bands = list(wvl[90:103]) # Select bands (PE563, 88:113), CHL [12:36]
                conv_h_q = spectro.FeaturesConvexHullQuotient(pre_continuum, bands, baseline=0.9999)
                auc = conv_h_q.get_area(feat_no=1)
                mbd = conv_h_q.get_absorbtion_depth(feat_no=1)
                anmb = auc/mbd
                
                pigment_map[y,x] = 0.02692 + (0.06631 * anmb) # ANMB 563 nm 
    
    anmb563_map = pigment_map.copy()
    
    plt.imshow(anmb563_map, interpolation='spline16', norm=norm, cmap='magma')
    plt.title(sample_key[i] + ' R-Phycoerythrin - ANMBD 563 nm Δ 20 nm', fontsize=15)
    cbar = plt.colorbar()
    cbar.ax.get_yaxis().labelpad = 16
    cbar.ax.set_ylabel('$\mathregular{mg/mm^{2}}$', rotation=270, fontsize=16)
    plt.tight_layout()
    # plt.savefig(out_p + sample_key[i] + '-anmb563_map.png', dpi=600)
    plt.show()
    plt.close()
    
    ## ANMB 563_40 ~~~~~~~~~~~~~~~~~~~~~
    pigment_map = np.empty((sample.shape[0], sample.shape[1], 1), dtype=np.float64)
    
    for y in range(0, sample.shape[0]):
        for x in range(0, sample.shape[1]):
            if sample[y,x].any() == 0:
                pigment_map[y,x] = 0
            else:
                # sav_gol_filter = ss.savgol_filter(sample[y,x][88:113], 3, 1)
                region = sample[y,x][85:108]
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
    plt.title(sample_key[i] + ' R-Phycoerythrin - ANMBD 563 nm Δ 40 nm', fontsize=15)
    cbar = plt.colorbar()
    cbar.ax.get_yaxis().labelpad = 16
    cbar.ax.set_ylabel('$\mathregular{mg/mm^{2}}$', rotation=270, fontsize=16)
    plt.tight_layout()
    # plt.savefig(out_p + sample_key[i] + '-anmb563_40_map_nointerp.png', dpi=600)
    plt.show()
    plt.close()
    
    #ND 563 nm ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    pigment_map = np.empty((sample.shape[0], sample.shape[1], 1), dtype=np.float64)
    
    for y in range(0, sample.shape[0]):
        for x in range(0, sample.shape[1]):
            if sample[y,x].any() == 0:
                pigment_map[y,x] = 0
            else:
                pigment_map[y,x] = 0.01933 + (0.16628 * (
                    (sample[y,x,173] - sample[y,x,96]) / (sample[y,x,173] + sample[y,x,96])))
                
    ND563_map = pigment_map.copy()
    
    plt.imshow(ND563_map, interpolation=None, norm=norm, cmap='magma')
    plt.title(sample_key[i] + 'R-Phycoerythrin - ND563', fontsize=15)
    cbar = plt.colorbar()
    cbar.ax.get_yaxis().labelpad = 16
    cbar.ax.set_ylabel('$\mathregular{mg/mm^{2}}$', rotation=270, fontsize=16)
    plt.tight_layout()
    # plt.savefig(out_p + sample_key[i] + '-ND563_map_nointerp.png', dpi=600)
    plt.show()
    plt.close()
    
    # #R665 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # pigment_map = np.empty((sample.shape[0], sample.shape[1], 1), dtype=np.float64)
    
    # for y in range(0, sample.shape[0]):
    #     for x in range(0, sample.shape[1]):
    #         if sample[y,x].any() == 0:
    #             pigment_map[y,x] = 0
    #         else:
    #             pigment_map[y,x] = -2.62970 + (47.60657 * sample[y,x,156]) # R665 nm
    
    # R665_map = pigment_map.copy()
    
    # plt.imshow(R665_map, interpolation=None, norm=norm_chl, cmap='viridis')
    # plt.title(sample_key[i] + 'Chlorophyll / Phaeophytin - R 665 nm', fontsize=15)
    # cbar = plt.colorbar()
    # cbar.ax.get_yaxis().labelpad = 16
    # cbar.ax.set_ylabel('$\mathregular{mg/mm^{2}}$', rotation=270, fontsize=16)
    # plt.tight_layout()
    # plt.savefig(out_p + sample_key[i] + '-R665_map_nointerp.png', dpi=600)
    # plt.show()
    # plt.close()
    
    
    # # Subplot R-PE Maps ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    out_p = r"C:\Users\jcmontes\OneDrive - University of Tasmania\01_Projects_Drive\Imaging_spectroscopy\Phenotyping_macroalgae\results\Hyperspectral_Imaging\R-PE-Maps-Subplot/"

    fig, ((ax0,ax1), (ax2,ax3)) = plt.subplots(2,2, figsize=(13,7), squeeze=True, sharey=False)
    plt.subplots_adjust(wspace=0, hspace=0)
    
    ax0.imshow(sample[:,:,(178,95,15)])
    ax0.title.set_text(sample_key[i] + ' RGB (703, 560, 424 nm)')
    
    # im1 = ax1.imshow(d563_map, interpolation=None, norm=norm, cmap='magma')
    # divider = make_axes_locatable(ax1)
    # cax = divider.append_axes('right', size='5%', pad=0.05)
    # cbar = fig.colorbar(im1, cax=cax, orientation='vertical')
    # cbar.ax.get_yaxis().labelpad = 16
    # cbar.ax.set_ylabel('$\mathregular{mg/mm^{2}}$', rotation=270, fontsize=13)
    # ax1.title.set_text(sample_key[i] + ' R-Phycoerythrin - δ 563 nm')
    
    im1 = ax1.imshow(ND563_map, interpolation=None, norm=norm, cmap='magma')
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    cbar = fig.colorbar(im1, cax=cax, orientation='vertical')
    cbar.ax.get_yaxis().labelpad = 16
    cbar.ax.set_ylabel('$\mathregular{mg/mm^{2}}$', rotation=270, fontsize=13)
    ax1.title.set_text(sample_key[i] + 'R-Phycoerythrin - ND 563 nm')
    
    im2 = ax2.imshow(dd569_map, interpolation=None, norm=norm, cmap='magma')
    divider = make_axes_locatable(ax2)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    cbar = fig.colorbar(im2, cax=cax, orientation='vertical')
    cbar.ax.get_yaxis().labelpad = 16
    cbar.ax.set_ylabel('$\mathregular{mg/mm^{2}}$', rotation=270, fontsize=13)
    ax2.title.set_text(sample_key[i] + 'R-Phycoerythrin - δδ 569 nm')
    
    im3 = ax3.imshow(anmb563_map, interpolation=None, norm=norm, cmap='magma')
    divider = make_axes_locatable(ax3)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    cbar = fig.colorbar(im3, cax=cax, orientation='vertical')
    cbar.ax.get_yaxis().labelpad = 16
    cbar.ax.set_ylabel('$\mathregular{mg/mm^{2}}$', rotation=270, fontsize=13)
    ax3.title.set_text(sample_key[i] + ' R-Phycoerythrin - ANMBD [553 to 573 nm]')
    
    # im5 = ax5.imshow(anmb563_40map, interpolation=None, norm=norm, cmap='magma')
    # divider = make_axes_locatable(ax5)
    # cax = divider.append_axes('right', size='5%', pad=0.05)
    # cbar = fig.colorbar(im3, cax=cax, orientation='vertical')
    # cbar.ax.get_yaxis().labelpad = 16
    # cbar.ax.set_ylabel('$\mathregular{mg/mm^{2}}$', rotation=270, fontsize=13)
    # ax5.title.set_text(sample_key[i] + ' R-Phycoerythrin - ANMBD [543 to 583 nm]')
    
    plt.tight_layout(pad=0, h_pad=1, w_pad=0,)
    plt.savefig(out_p + sample_key[i] + '-R-PE_maps.png', dpi=600)
    plt.show()
    
    
    # # R-PE Zoom Subplot figure ~~~~~~~~~~~~~~~~~~~~~
    # # Change output path
    # out_p = r"C:\Users\jcmontes\OneDrive - University of Tasmania\01_Projects_Drive\Imaging_spectroscopy\Phenotyping_macroalgae\results\Hyperspectral_Imaging\R-PE-Maps-Zoom/"
    
    # # Zoom rectangle
    # rect = patches.Rectangle((0, 30), 60, 50, linewidth=1, edgecolor='r', facecolor='none')
    # # im0 = ax0.imshow(_map[20:70,20:85](Different Zoom region)
    
    # # Figure
    # fig, ((ax0,ax1), (ax2,ax3)) = plt.subplots(2,2, figsize=(13,7), squeeze=True, sharey=False)
    # # im0 = ax0.imshow(anmb563_map, interpolation='none', norm=None, cmap='magma')
    # plt.subplots_adjust(wspace=0, hspace=0)
    
    # ax0.imshow(sample[:,:,(178,95,15)])
    # ax0.title.set_text(sample_key[i] + ' RGB - Red (698 nm) Green (555 nm) Blue (419 nm)')
    # ax0.add_patch(rect)
    
    # im1 = ax1.imshow(d563_map[30:80, 0:60], cmap='magma')
    # divider = make_axes_locatable(ax1)
    # cax = divider.append_axes('right', size='5%', pad=0.05)
    # cbar = fig.colorbar(im1, cax=cax, orientation='vertical')
    # cbar.ax.get_yaxis().labelpad = 16
    # cbar.ax.set_ylabel('$\mathregular{mg/mm^{2}}$', rotation=270, fontsize=16)
    # ax1.title.set_text(sample_key[i] + ' R-Phycoerythrin - δ 563 nm')
    
    # im2 = ax2.imshow(dd569_map[30:80, 0:60], cmap='magma')
    # divider = make_axes_locatable(ax2)
    # cax = divider.append_axes('right', size='5%', pad=0.05)
    # cbar = fig.colorbar(im2, cax=cax, orientation='vertical')
    # cbar.ax.get_yaxis().labelpad = 16
    # cbar.ax.set_ylabel('$\mathregular{mg/mm^{2}}$', rotation=270, fontsize=16)
    # ax2.title.set_text(sample_key[i] + 'R-Phycoerythrin - δδ 569 nm')
    
    # im3 = ax3.imshow(anmb563_map[30:80, 0:60], cmap='magma')
    # divider = make_axes_locatable(ax3)
    # cax = divider.append_axes('right', size='5%', pad=0.05)
    # cbar = fig.colorbar(im3, cax=cax, orientation='vertical')
    # cbar.ax.get_yaxis().labelpad = 16
    # cbar.ax.set_ylabel('$\mathregular{mg/mm^{2}}$', rotation=270, fontsize=16)
    # ax3.title.set_text(sample_key[i] + ' R-Phycoerythrin - ANMBD 563 nm Δ 20 nm')
    
    # plt.tight_layout(pad=0, h_pad=1, w_pad=0,)
    # # plt.savefig(out_p + sample_key[i] + '-R-PE-Zoom.png', dpi=600)
    # plt.show()
    
    # Zoom R-PE vs Chlorophyl/Phaeophytin ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Change output path
    # out_p = r"C:\Users\jcmontes\OneDrive - University of Tasmania\01_Projects_Drive\Imaging_spectroscopy\Phenotyping_macroalgae\results\Hyperspectral_Imaging\R-PEvsChl-Maps/"
    
    # # Zoom rectangle
    # # rect = patches.Rectangle((20, 30), 60, 50, linewidth=1, edgecolor='r', facecolor='none')
    
    # # Figure
    # fig, ((ax0,ax1), (ax2,ax3)) = plt.subplots(2,2, figsize=(13,7), squeeze=True)
    # # im0 = ax0.imshow(anmb563_map, interpolation='none', norm=None, cmap='magma')
    # plt.subplots_adjust(wspace=0, hspace=0)
    
    # ax0.imshow(sample[:,:,(178,95,15)])
    # ax0.title.set_text(sample_key[i] + ' RGB - Red (698 nm) Green (555 nm) Blue (419 nm)')
    # # ax0.add_patch(rect)
    
    # im1 = ax1.imshow(R665_map, cmap='viridis', interpolation='spline16', norm=norm_chl)
    # divider = make_axes_locatable(ax1)
    # cax = divider.append_axes('right', size='5%', pad=0.1)
    # cbar = fig.colorbar(im1, cax=cax, orientation='vertical')
    # cbar.ax.get_yaxis().labelpad = 16
    # cbar.ax.set_ylabel('$\mathregular{mg/mm^{2}}$', rotation=270, fontsize=16)
    # ax1.title.set_text(sample_key[i] + ' Chlorophyll / Phaeophytin - R 665 nm')
    
    # im2 = ax2.imshow(dd569_map, cmap='magma', interpolation='spline16', norm=norm)
    # divider = make_axes_locatable(ax2)
    # cax = divider.append_axes('right', size='5%', pad=0.05)
    # cbar = fig.colorbar(im2, cax=cax, orientation='vertical')
    # cbar.ax.set_ylabel('$\mathregular{mg/mm^{2}}$', rotation=270, fontsize=16)
    # cbar.ax.get_yaxis().labelpad = 16
    # ax2.title.set_text(sample_key[i] + 'R-Phycoerythrin - δδ 569 nm')
    
    # im3 = ax3.imshow(sample[:,:,173])
    # ax3.title.set_text(sample_key[i] + '695 nm - Red edge')
    
    # plt.tight_layout(pad=0, h_pad=1, w_pad=0,)
    # plt.savefig(out_p + sample_key[i] + '-R-PEvsChl.png', dpi=600)
    # plt.show()
    
    i += 1

print('All Done')
