from nipype import Node, Workflow
import nipype.interfaces.mrtrix3 as mrt3
import nipype.interfaces.fsl as fsl
import os
import matplotlib.pyplot as plt
from pathlib import Path
import shutil
import argparse
import fnmatch
import sys


#TO-DO Add subject ID flags




parser = argparse.ArgumentParser()
parser.add_argument("--pe1", "-pe1", help="dMRI file using first phase-encoding direction.", required=False)
parser.add_argument("--pe2", "-pe2", help="dMRI file using second phase-encoding direction.", required=False)
parser.add_argument("--subject", "-s", help="Path to subject directory.", required=True)
parser.add_argument("--config", "-c", help="Path to text eddy-style config file containing phase-encoding information. If not provided, default will be sourced from image header. IMPORTANT: PHILLIPS Scanners will require their own config and index files for eddy/topup.", required=False)
parser.add_argument("--index", "-i", help="Path to eddy-style index file containing phase-encoding information. If not provided, default will be sourced from image header. IMPORTANT: PHILLIPS Scanners will require their own config and index files for eddy/topup.", required=False)
parser.add_argument("--single_pe", help="Run diffusion pipeline on a single phase-encoding direction acquired dMRI sequence.", required=False) # TODO: Add this step!!!!!
parser.add_argument("--derivatives", "-d", help="Path to derivatives directory. If not provided, default is sourcedata directory.", required=False)
args = parser.parse_args()


path = args.subject
PE_1 = args.pe1
PE_2 = args.pe2

os.chdir(path)
sys.stdout = open('dwi_log.txt', 'w')

# acqp_path = '/autofs/space/nicc_001/users/rcali/test_environment/acqp.txt'
# index_path = '/autofs/space/nicc_001/users/rcali/test_environment/index.txt'





pe_dir_1 = ['AP', 'RL', 'SI']
pe_dir_2 = ['PA', 'LR', 'IS']

data_dir = ["orig", "mrtrix_files", "dtifit"]


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

def mrzip(in_file):
    if in_file.endswith('.nii'):
        mrconvert = mrt3.MRConvert()
        mrconvert.inputs.in_file = in_file
        mrconvert.inputs.out_file = in_file + '.gz'
        mrconvert.inputs.args = '-force'
        mrconvert.run()
        os.remove(in_file)
        return

def throw_error(output_file):
    if os.path.exists(output_file):
        print(output_file + " generated, continuing processing pipeline...")
    else:
        print("Error generating [" + output_file + "], check last step for potential issues.")
        sys.exit()

organize_dir(data_dir[0])

ext = [".nii.gz", ".nii", ".bvec", ".bval", ".json"]

orig_dir = os.path.join(path, "orig")

file_list = []

for file in files(path):
    if file.endswith(tuple(ext)):
        file_list.append(file)
        
        shutil.copy(file, orig_dir)

        

print("-----File Manifest-----")
print("")
for xx in file_list:
    print(xx)
    print("")


direction_1_string = Path(PE_1).stem.split(".")[0]
direction_2_string = Path(PE_2).stem.split(".")[0]



dir_1_list = []
dir_2_list = []



for file in files(path):
    for pe_dir in pe_dir_1:
        if pe_dir[:] in file:
            dir_1_list.append(file)  

                

for file2 in files(path):
    for pe_dir in pe_dir_2:
        if pe_dir[:] in file2:
            dir_2_list.append(file2)



# if filename contains 'AP', 'RL', 'SI'....

for a in dir_1_list:
    if a.endswith(".nii.gz"):
        os.rename(a, "dwi_PE_1.nii.gz")
    elif a.endswith(".nii"):
        os.rename(a, "dwi_PE_1.nii")
    elif a.endswith(".bvec"):
        os.rename(a, "dwi_PE_1.bvec")
    elif a.endswith(".bval"):
        os.rename(a, "dwi_PE_1.bval")
    elif a.endswith(".json"):
        os.rename(a, "dwi_PE_1.json")

# # IMPORTANT - All files belonging to a dwi volume MUST CONTAIN the same basename (i.e., dwi.nii, dwi.bvec, dwi.bval, dwi.json)

for t in dir_2_list:
    if t.endswith(".nii.gz"):
        os.rename(t, "dwi_PE_2.nii.gz")
    elif t.endswith(".nii"):
        os.rename(t, "dwi_PE_2.nii")
    elif t.endswith(".bvec"):
        os.rename(t, "dwi_PE_2.bvec")
    elif t.endswith(".bval"):
        os.rename(t, "dwi_PE_2.bval")
    elif t.endswith(".json"):
        os.rename(t, "dwi_PE_2.json")




for qq in os.listdir(path):
    if qq.endswith('.nii'):
        mrzip(qq)



# mrconvert = mrt3.MRConvert()
# mrconvert.inputs.in_file = "dwi_PE_1.nii.gz"
# mrconvert.inputs.out_file = "dwi_PE_1.mif"
# mrconvert.inputs.args = (
#     "-fslgrad dwi_PE_1.bvec dwi_PE_1.bval -json_import dwi_PE_1.json -import_pe_eddy /autofs/space/nicc_001/users/rcali/test_environment/acqp.txt /autofs/space/nicc_001/users/rcali/test_environment/index.txt -force"
# )
# mrconvert.run()

# throw_error("dwi_PE_1.mif")


mrconvert = mrt3.MRConvert()
mrconvert.inputs.in_file = "dwi_PE_1.nii.gz"
mrconvert.inputs.out_file = "dwi_PE_1.mif"
mrconvert.inputs.args = (
    "-fslgrad dwi_PE_1.bvec dwi_PE_1.bval -json_import dwi_PE_1.json -force"
)
mrconvert.run()

throw_error("dwi_PE_1.mif")



mrconvert = mrt3.MRConvert()
mrconvert.inputs.in_file = "dwi_PE_2.nii.gz"
mrconvert.inputs.out_file = "dwi_PE_2.mif"
mrconvert.inputs.args = (
    "-fslgrad dwi_PE_2.bvec dwi_PE_2.bval -json_import dwi_PE_2.json -force"
)
mrconvert.run()

throw_error("dwi_PE_2.mif")

mrcat = mrt3.MRCat()
mrcat.inputs.in_files = ["dwi_PE_1.mif", "dwi_PE_2.mif"]
mrcat.inputs.out_file = "dwi_merged_PE.mif"
mrcat.inputs.args = "-force"
mrcat.run()

throw_error("dwi_merged_PE.mif")

mrconvert = mrt3.MRConvert()
mrconvert.inputs.in_file = "dwi_merged_PE.mif"
mrconvert.inputs.out_file = "dwi_merged_PE.mif"
mrconvert.inputs.args = (
    "-export_pe_eddy index.txt config.txt -export_grad_fsl bvecs bvals -force"
)
mrconvert.run()

throw_error("dwi_merged_PE.mif")


dwidenoise = mrt3.DWIDenoise()
dwidenoise.inputs.in_file = "dwi_merged_PE.mif"
dwidenoise.inputs.out_file = "denoised_merged_dwi.mif"
dwidenoise.inputs.noise = "noise.mif"
dwidenoise.inputs.args = "-extent 5 -force"
dwidenoise.run()

throw_error("denoised_merged_dwi.mif")

dwiextract = mrt3.DWIExtract()
dwiextract.inputs.in_file = "denoised_merged_dwi.mif"
dwiextract.inputs.bzero = True
dwiextract.inputs.out_file = "b0_vols.mif"
dwiextract.inputs.grad_fsl = ("bvecs", "bvals")
dwiextract.run()

throw_error("b0_vols.mif")

# preproc = mrt3.DWIPreproc()
# preproc.inputs.in_file = "denoised_merged_dwi.mif"
# preproc.inputs.rpe_options = "all"
# preproc.inputs.out_file = "preproc.mif"
# preproc.inputs.eddy_options = (
#     "--slm=linear --cnr_maps"  # linear second level model and replace outliers
# )
# preproc.inputs.args = "-nthreads 5 -eddyqc_all eddy -force"
# preproc.inputs.export_grad_mrtrix = True
# preproc.inputs.pe_dir = "LR"
# preproc.run()



preproc = mrt3.DWIPreproc()
preproc.inputs.in_file = "denoised_merged_dwi.mif"
preproc.inputs.rpe_options = "all"
preproc.inputs.out_file = "preproc.mif"
preproc.inputs.eddy_options = (
    "--slm=linear --cnr_maps"  # linear second level model and replace outliers
)
preproc.inputs.args = "-nthreads 5 -eddyqc_all eddy -force"
preproc.inputs.export_grad_mrtrix = True
preproc.inputs.pe_dir = "LR"
preproc.run()




throw_error("preproc.mif")

mrconvert = mrt3.MRConvert()
mrconvert.inputs.in_file = "preproc.mif"
mrconvert.inputs.out_file = "preproc.mif"
mrconvert.inputs.args = "-export_grad_fsl eddy_bvecs eddy_bvals -force"
mrconvert.run()

throw_error("preproc.mif")

mrconvert = mrt3.MRConvert()
mrconvert.inputs.in_file = "dwi_merged_PE.mif"
mrconvert.inputs.out_file = "dwi_merged_PE.nii.gz"
mrconvert.run()

throw_error("dwi_merged_PE.nii.gz")

mrconvert = mrt3.MRConvert()
mrconvert.inputs.in_file = "noise.mif"
mrconvert.inputs.out_file = "noise.nii.gz"
mrconvert.run()

throw_error("noise.nii.gz")

# Compute residuals by subtracting noise from raw image
maths = fsl.ImageMaths()
maths.inputs.in_file = "dwi_merged_PE.nii.gz"
maths.inputs.args = "-sub noise.nii.gz"
maths.inputs.out_file = "residuals.nii.gz"
maths.run()

throw_error("residuals.nii.gz")

dwiextract = mrt3.DWIExtract()
dwiextract.inputs.in_file = "preproc.mif"
dwiextract.inputs.bzero = True
dwiextract.inputs.out_file = "preproc_b0_vols.mif"
dwiextract.inputs.grad_fsl = ("eddy_bvecs", "eddy_bvals")
dwiextract.run()

throw_error("preproc_b0_vols.mif")

mrconvert = mrt3.MRConvert()
mrconvert.inputs.in_file = "preproc_b0_vols.mif"
mrconvert.inputs.out_file = "preproc_b0_vols.nii.gz"
mrconvert.run()

throw_error("preproc_b0_vols.nii.gz")

mrconvert = mrt3.MRConvert()
mrconvert.inputs.in_file = "preproc.mif"
mrconvert.inputs.out_file = "preproc.nii.gz"
mrconvert.run()

throw_error("preproc.nii.gz")

bet = fsl.BET()
bet.inputs.in_file = "preproc_b0_vols.nii.gz"
bet.inputs.mask = True
bet.inputs.frac = 0.15
bet.run()

throw_error("preproc_b0_vols.nii.gz")

dti = fsl.DTIFit()
dti.inputs.dwi = "preproc.nii.gz"
dti.inputs.bvecs = "eddy_bvecs"
dti.inputs.bvals = "eddy_bvals"
dti.inputs.base_name = "dtifit"
dti.inputs.mask = "preproc_b0_vols_brain_mask.nii.gz"
dti.run()


mrtrix_dir = os.path.join(path, "mrtrix_files")

organize_dir(data_dir[1])

eddy_dir = os.path.join(path, 'eddy')
trix_files = [".mif", "grad.b"]
qc_files = ["residuals.nii.gz", "noise.nii.gz"]
dti_fit_files = ['dtifit_FA.nii.gz', 'dtifit_L1.nii.gz', 'dtifit_L2.nii.gz', 'dtifit_L3.nii.gz', 'dtifit_MD.nii.gz', 'dtifit_MO.nii.gz', 'dtifit_S0.nii.gz', 'dtifit_V1.nii.gz', 'dtifit_V2.nii.gz', 'dtifit_V3.nii.gz']
for s in files(path):
    if s.endswith(tuple(trix_files[:])):
        b = os.path.join(path, s)
        c = os.path.join(mrtrix_dir, s)
        os.replace(b, c)
    elif s.startswith(tuple(qc_files[:])):
        b = os.path.join(path, s)
        c = os.path.join(eddy_dir, s)
        os.replace(b, c)

for dti_fit_file in dti_fit_files:
    throw_error(dti_fit_file)
        
dti_fit_dir = os.path.join(path, "dtifit")

organize_dir(data_dir[2])

for o in files(path):
    if o.startswith("dtifit"):
        q = os.path.join(path, o)
        y = os.path.join(dti_fit_dir, o)
        os.replace(q, y)

for w in files(path):
    if w.startswith('dwi_PE'):
        os.remove(w)






sys.stdout.close()