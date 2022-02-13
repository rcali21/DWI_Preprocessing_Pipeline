from nipype import Node, Workflow
import nipype.interfaces.mrtrix3 as mrt3
import nipype.interfaces.mrtrix as mrt
from nipype.interfaces import fsl
import os
import glob
import nibabel as nb
import matplotlib.pyplot as plt
from pathlib import Path
import time
import shutil
import argparse

path = os.getcwd()

parser = argparse.ArgumentParser()
parser.add_argument("--pe1", "-pe1")
parser.add_argument("--pe2", "-pe2")
parser.add_argument("--bvec", "-bvec", required=False)
parser.add_argument("--bval", "-bval", required=False)
args = parser.parse_args()

PE_1 = args.pe1
PE_2 = args.pe2


data_dir = ["orig", "mrtrix_files", "dtifit", "bedpostx"]


def files(path):
    for file in os.listdir(path):
        if os.path.isfile(os.path.join(path, file)):
            yield file


def organize_dir(dir_name):
    if not os.path.isdir(dir_name):
        print(dir_name + " directory created, all related files will be moved here...")
        os.mkdir(dir_name)
    elif os.path.isdir(dir_name):
        print(
            dir_name
            + " directory found, moving all related outputs to this directory..."
        )


organize_dir(data_dir[0])


ext = [".nii.gz", ".nii ", ".bvec", ".bval", ".json"]

orig_dir = os.path.join(path, "orig")

file_list = []

for file in files(path):
    if file.endswith(tuple(ext)):
        file_list.append(file)
        shutil.copy(file, orig_dir)

print("-----File Manifest-----")
print("")
for xx in file_list:
    print("")
    print(xx + " \u2713")

direction_1_string = Path(PE_1).stem.split(".")[0]
direction_2_string = Path(PE_2).stem.split(".")[0]

dir_1_list = []
dir_2_list = []


for file in files(path):
    if file.startswith(direction_1_string):
        # print(os.path.join(path, file))
        dir_1_list.append(file)

for file in files(path):
    if file.startswith(direction_2_string):
        # print(os.path.join(path, file))
        dir_2_list.append(file)

for a in dir_1_list:
    if a.endswith(".nii.gz" or ".nii"):
        os.rename(a, "dwi_PE_1.nii.gz")

    elif a.endswith(".bvec"):
        os.rename(a, "dwi_PE_1.bvec")

    elif a.endswith(".bval"):
        os.rename(a, "dwi_PE_1.bval")

    elif a.endswith(".json"):
        os.rename(a, "dwi_PE_1.json")

# IMPORTANT - All files belonging to a dwi volume MUST CONTAIN the same basename (i.e., dwi.nii, dwi.bvec, dwi.bval, dwi.json)

for t in dir_2_list:
    if t.endswith(".nii.gz" or ".nii"):
        os.rename(t, "dwi_PE_2.nii.gz")

    elif t.endswith(".bvec"):
        os.rename(t, "dwi_PE_2.bvec")

    elif t.endswith(".bval"):
        os.rename(t, "dwi_PE_2.bval")

    elif t.endswith(".json"):
        os.rename(t, "dwi_PE_2.json")


mrconvert = mrt3.MRConvert()
mrconvert.inputs.in_file = "dwi_PE_1.nii.gz"
mrconvert.inputs.out_file = "dwi_PE_1.mif"
mrconvert.inputs.args = (
    "-fslgrad dwi_PE_1.bvec dwi_PE_1.bval -json_import dwi_PE_1.json -force"
)
mrconvert.run()

mrconvert = mrt3.MRConvert()
mrconvert.inputs.in_file = "dwi_PE_2.nii.gz"
mrconvert.inputs.out_file = "dwi_PE_2.mif"
mrconvert.inputs.args = (
    "-fslgrad dwi_PE_2.bvec dwi_PE_2.bval -json_import dwi_PE_2.json -force"
)
mrconvert.run()

mrcat = mrt3.MRCat()
mrcat.inputs.in_files = ["dwi_PE_1.mif", "dwi_PE_2.mif"]
mrcat.inputs.out_file = "dwi_merged_PE.mif"
mrcat.run()

mrconvert = mrt3.MRConvert()
mrconvert.inputs.in_file = "dwi_merged_PE.mif"
mrconvert.inputs.out_file = "dwi_merged_PE.mif"
mrconvert.inputs.args = (
    "-export_pe_eddy index.txt config.txt -export_grad_fsl bvecs bvals -force"
)
mrconvert.run()


dwidenoise = mrt3.DWIDenoise()
dwidenoise.inputs.in_file = "dwi_merged_PE.mif"
dwidenoise.inputs.out_file = "denoised_merged_dwi.mif"
dwidenoise.inputs.noise = "noise.mif"
dwidenoise.inputs.args = "-extent 5 -force"
dwidenoise.run()

dwiextract = mrt3.DWIExtract()
dwiextract.inputs.in_file = "denoised_merged_dwi.mif"
dwiextract.inputs.bzero = True
dwiextract.inputs.out_file = "b0_vols.mif"
dwiextract.inputs.grad_fsl = ("bvecs", "bvals")
dwiextract.run()

preproc = mrt3.DWIPreproc()
preproc.inputs.in_file = "denoised_merged_dwi.mif"
preproc.inputs.rpe_options = "all"
preproc.inputs.out_file = "preproc.mif"
preproc.inputs.eddy_options = (
    "--slm=linear --cnr_maps"  # linear second level model and replace outliers
)
preproc.inputs.args = "-nthreads 5 -eddyqc_all eddy -force"
preproc.inputs.export_grad_mrtrix = True
preproc.inputs.pe_dir = "AP"
preproc.run()

mrconvert = mrt3.MRConvert()
mrconvert.inputs.in_file = "preproc.mif"
mrconvert.inputs.out_file = "preproc.mif"
mrconvert.inputs.args = "-export_grad_fsl eddy_bvecs eddy_bvals -force"
mrconvert.run()

mrconvert = mrt3.MRConvert()
mrconvert.inputs.in_file = "dwi_merged_PE.mif"
mrconvert.inputs.out_file = "dwi_merged_PE.nii.gz"
mrconvert.run()

mrconvert = mrt3.MRConvert()
mrconvert.inputs.in_file = "noise.mif"
mrconvert.inputs.out_file = "noise.nii.gz"
mrconvert.run()

# Compute residuals by subtracting noise from raw image
maths = fsl.ImageMaths()
maths.inputs.in_file = "dwi_merged_PE.nii.gz"
maths.inputs.args = "-sub noise.nii.gz"
maths.inputs.out_file = "residuals.nii.gz"
maths.run()

dwiextract = mrt3.DWIExtract()
dwiextract.inputs.in_file = "preproc.mif"
dwiextract.inputs.bzero = True
dwiextract.inputs.out_file = "preproc_b0_vols.mif"
dwiextract.inputs.grad_fsl = ("eddy_bvecs", "eddy_bvals")
dwiextract.run()

mrconvert = mrt3.MRConvert()
mrconvert.inputs.in_file = "preproc_b0_vols.mif"
mrconvert.inputs.out_file = "preproc_b0_vols.nii.gz"
mrconvert.run()


mrconvert = mrt3.MRConvert()
mrconvert.inputs.in_file = "preproc.mif"
mrconvert.inputs.out_file = "preproc.nii.gz"
mrconvert.run()

bet = fsl.BET()
bet.inputs.in_file = "preproc_b0_vols.nii.gz"
bet.inputs.mask = True
bet.inputs.frac = 0.15
bet.run()

dti = fsl.DTIFit()
dti.inputs.dwi = "preproc.nii.gz"
dti.inputs.bvecs = "eddy_bvecs"
dti.inputs.bvals = "eddy_bvals"
dti.inputs.base_name = "dtifit"
dti.inputs.mask = "preproc_b0_vols_brain_mask.nii.gz"
dti.run()

mrtrix_dir = os.path.join(path, "mrtrix_files")

organize_dir(data_dir[1])

for s in files(path):
    if s.endswith(".mif"):
        b = os.path.join(path, s)
        c = os.path.join(mrtrix_dir, s)
        os.replace(b, c)


dti_fit_dir = os.path.join(path, "dtifit")

organize_dir(data_dir[2])

for o in files(path):
    if o.startswith("dtifit"):
        q = os.path.join(path, o)
        y = os.path.join(dti_fit_dir, o)
        os.replace(q, y)

organize_dir(data_dir[3])
