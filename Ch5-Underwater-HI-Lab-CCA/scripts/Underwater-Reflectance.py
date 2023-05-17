import os
from pathlib import Path
import numpy as np
import spectral as sp
# import scipy.signal as ss
# import matplotlib.pyplot as plt

## Function: Radiance to Reflectance - Underwater
def radiance2reflectance(img_array, bands, dry_plastic, wet_plastic):
    unfold_array = img_array.reshape(-1, len(bands))
    _reflectance = (np.divide(unfold_array, wet_plastic) * dry_plastic)
    img_reflectance = np.reshape(_reflectance, img_array.shape)
    print('Completed Radiance to Reflectance')
    return img_reflectance

## Directories
files_path = Path(r"C:\Users\jcmontes\OneDrive - University of Tasmania\01_Projects_Drive\UHI_Tank\data\JCMH-FC-RadCorr")
hdr_list = list(files_path.glob('*.HDR'))
rad_list = list(files_path.glob('*.DAT'))
out_path = Path(r"C:\Users\jcmontes\OneDrive - University of Tasmania\01_Projects_Drive\UHI_Tank\data\JCMH-FC-Reflectance")
# sg_out_path = Path(r"C:\Users\jcmontes\OneDrive - University of Tasmania\01_Projects_Drive\UHI_Tank\data\JCMH-FC-SG-Reflectance")

## Reference plaques
open_refl_spectralon = sp.envi.open(r"C:\Users\jcmontes\OneDrive - University of Tasmania\01_Projects_Drive\UHI_Tank\data\JCMH-FC-Refl-Panels\Refl_spectralon.hdr",
                               r"C:\Users\jcmontes\OneDrive - University of Tasmania\01_Projects_Drive\UHI_Tank\data\JCMH-FC-Refl-Panels\Refl_spectralon.img")

refl_spectralon = open_refl_spectralon.load()

r_spectralon = refl_spectralon.mean(axis=(0,1))
lu_spectralon = np.load(r"C:\Users\jcmontes\OneDrive - University of Tasmania\01_Projects_Drive\UHI_Tank\data\JCMH-FC-Numpy\Spectralon.npy")
lu_dry_pvc = np.load(r"C:\Users\jcmontes\OneDrive - University of Tasmania\01_Projects_Drive\UHI_Tank\data\JCMH-FC-Numpy\White-Dry-PVC.npy")
lu_wet_pvc = np.load(r"C:\Users\jcmontes\OneDrive - University of Tasmania\01_Projects_Drive\UHI_Tank\data\JCMH-FC-Numpy\White-UW-PVC.npy")

r_dry_pvc = (np.divide(lu_dry_pvc, lu_spectralon) * r_spectralon)

## Reflectance loop
for i, file in enumerate(hdr_list):
    print(hdr_list[i])
    target_open = sp.envi.open(hdr_list[i], rad_list[i])
    bands_vector = target_open.metadata['wavelength']
    target_data = target_open.load()
    target_reflectance = radiance2reflectance(target_data, bands_vector, r_dry_pvc, lu_wet_pvc)
    reflectance = target_reflectance[:,190:700,:]
    # sg_reflectance = ss.savgol_filter(reflectance, 3, 1, axis=2)
    os.chdir(out_path)
    sp.envi.save_image(hdr_list[i].stem.replace('radiance','reflectance')+'.hdr', reflectance, interleave='bil', metadata=target_open.metadata)
    sp.save_rgb(hdr_list[i].stem.replace('radiance','rgb')+'.tiff', reflectance, [144,92,42])
    # os.chdir(sg_out_path)
    # sp.envi.save_image(hdr_list[i].stem.replace('radiance','reflectance_SG')+'.hdr', sg_reflectance, interleave='bil', metadata=target_open.metadata)
    # sp.save_rgb(hdr_list[i].stem.replace('radiance','rgb')+'.tiff', sg_reflectance, [144,92,42])

print('All done :)')


