# -*- coding: utf-8 -*-
"""
Created on Fri Jun  3 11:44:28 2022

@author: jcmontes
"""

import os
from pathlib import Path
import numpy as np
import spectral as sp
import scipy.signal as ss

## Function: Radiance to Reflectance
def radiance2reflectance(img_array, bands, spectralon):
    unfold_array = img_array.reshape(-1, len(bands))
    _reflectance = np.divide(unfold_array, spectralon)
    img_reflectance = np.reshape(_reflectance, img_array.shape)
    print('Completed Radiance to Reflectance')
    return img_reflectance

## Directories - Plastic
# white_path = Path(r"R:\SET\Spatial Science\UAV\_data\IMAS\Marine_plastics_project\20210324_FX17_plastics_scans\calib_reflec_targets_corrected")
# files_path = Path(r"R:\SET\Spatial Science\UAV\_data\IMAS\Marine_plastics_project\20210324_FX17_plastics_scans\bird_gut_corrected")
# out_path = Path(r"R:\SET\Spatial Science\UAV\_data\IMAS\Marine_plastics_project\20220425_FX17_plastic_scans_Reflectance\Bird_guts_reflectance")
# sg_out_path = Path(r"R:\SET\Spatial Science\UAV\_data\IMAS\Marine_plastics_project\20220425_FX17_plastic_scans_Reflectance\Bird_guts_reflectance_SG")

## Directories - CCA
files_path = Path(r"C:\Users\jcmontes\OneDrive - University of Tasmania\01_Projects_Drive\Imaging_spectroscopy\Phenotyping_macroalgae\data\Kestrel_CCA_Radiance")
hdr_list = list(files_path.glob('*.hdr'))
bin_list = list(files_path.glob('*.dat'))
out_path = Path(r"C:\Users\jcmontes\OneDrive - University of Tasmania\01_Projects_Drive\Imaging_spectroscopy\Phenotyping_macroalgae\data\Reflectance")
sg_out_path = Path(r"C:\Users\jcmontes\OneDrive - University of Tasmania\01_Projects_Drive\Imaging_spectroscopy\Phenotyping_macroalgae\data\SG_reflectance")

### Spectralon panel
## Plastics
# whiteref_hdr = white_path / 'FX17_spectral_evo_WR_2021-03-24_02-12-32_radiance.HDR'
# whiteref_dat = white_path / 'FX17_spectral_evo_WR_2021-03-24_02-12-32_radiance.DAT'

## CCA
whiteref_open = sp.envi.open(hdr_list[11], bin_list[11])

## Band metadata
bands_vector = whiteref_open.metadata['wavelength']

# Load spectralon
whiteref_data = whiteref_open.load()
whiteref_data = whiteref_data[360:450,890:1170,:] ## Y1:Y2,X1:X2,All bands

## Mean per band
spectralon = whiteref_data.mean(axis=(0, 1))

## Can save this
#np.savetxt(out_path / 'spectralon.csv', spectralon, delimiter=',')

# Create a list of HDR files
os.chdir(files_path)
#hdr_list = list(files_path.glob('*.hdr'))
#bin_list = list(files_path.glob('*.dat'))

for i, file in enumerate(hdr_list):
    print(hdr_list[i])
    target_open = sp.envi.open(hdr_list[i], bin_list[i])
    target_data = target_open.load()
    target_reflectance = radiance2reflectance(target_data, bands_vector, spectralon)
    reflectance = target_reflectance[:,500:1500,:]
    sg_reflectance = ss.savgol_filter(reflectance, 3, 1, axis=2)
    os.chdir(out_path)
    sp.envi.save_image(hdr_list[i].stem.replace('radiance','reflectance')+'.hdr', reflectance, interleave='bil', metadata=target_open.metadata)
    os.chdir(sg_out_path)
    sp.envi.save_image(hdr_list[i].stem.replace('radiance','reflectance_SG')+'.hdr', sg_reflectance, interleave='bil', metadata=target_open.metadata)

print('All done :)')