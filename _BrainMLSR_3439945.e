#######################################################################################
after load module ,plz source /public/software/apps/freesurfer_infant/freesurfer7.3.2/freesurfer/7.3.2-1/SetUpFreeSurfer.sh
#######################################################################################
mri_convert: /public/software/compiler/gnu/10.1.0/7//lib64/libgomp.so.1: no version information available (required by mri_convert)
mri_convert: /public/software/compiler/gnu/10.1.0/7//lib64/libgomp.so.1: no version information available (required by mri_convert)
WARNING ==================++++++++++++++++++++++++=======================================
The physical sizes are (176.00 mm, 219.87 mm, 256.16 mm), which cannot fit in 256^3 mm^3 volume.
The resulting volume will have 513 slices.
If you find problems, please let us know (freesurfer@nmr.mgh.harvard.edu).
==================================================++++++++++++++++++++++++===============

mri_convert: /public/software/compiler/gnu/10.1.0/7//lib64/libgomp.so.1: no version information available (required by mri_convert)
mri_convert: /public/software/compiler/gnu/10.1.0/7//lib64/libgomp.so.1: no version information available (required by mri_convert)
WARNING ==================++++++++++++++++++++++++=======================================
The physical sizes are (176.00 mm, 219.87 mm, 256.16 mm), which cannot fit in 256^3 mm^3 volume.
The resulting volume will have 513 slices.
If you find problems, please let us know (freesurfer@nmr.mgh.harvard.edu).
==================================================++++++++++++++++++++++++===============

Traceback (most recent call last):
  File "/home_data/home/caoshui2024/MutiSurf_pipeline/BrainMLSR/BrainMLSR/Step00_Register.py", line 2, in <module>
    import ants
ModuleNotFoundError: No module named 'ants'
