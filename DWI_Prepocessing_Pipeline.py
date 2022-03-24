import nipype.interfaces.mrtrix3 as mrt3
import nipype.interfaces.fsl as fsl
import os
import matplotlib.pyplot as plt
from pathlib import Path
import shutil
import argparse
import fnmatch
import sys


# TO-DO Add subject ID flags

parser = argparse.ArgumentParser()
parser.add_argument(
    "--pe1",
    "-pe1",
    help="dMRI file using first phase-encoding direction.",
    required=False,
)
parser.add_argument(
    "--pe2",
    "-pe2",
    help="dMRI file using second phase-encoding direction.",
    required=False,
)
parser.add_argument("--subject", "-s", help="Path to subject directory.", required=True)
parser.add_argument(
    "--config",
    "-c",
    help="Path to text eddy-style config file containing phase-encoding information. If not provided, default will be sourced from image header. IMPORTANT: PHILLIPS Scanners will require their own config and index files for eddy/topup.",
    required=False,
)
parser.add_argument(
    "--index",
    "-i",
    help="Path to eddy-style index file containing phase-encoding information. If not provided, default will be sourced from image header. IMPORTANT: PHILLIPS Scanners will require their own config and index files for eddy/topup.",
    required=False,
)
parser.add_argument(
    "--single_pe",
    help="Run diffusion pipeline on a single phase-encoding direction acquired dMRI sequence.",
    required=False) 
parser.add_argument(
    "--derivatives",
    "-d",
    help="Path to derivatives directory. If not provided, default is sourcedata directory.",
    required=False,
)
args = parser.parse_args()

# Assign vars to args
path = args.subject
PE_1 = args.pe1
PE_2 = args.pe2
# if args.single_pe:
#     single_pe_file = args.single_pe
# elif args.pe1 and args.pe2:
#     print('pe 2 if statement works')
# else:
#     print(args)




os.chdir(path)
# Initiate logging
# sys.stdout = open("dwi_log.txt", "w")

# acqp_path = '/autofs/space/nicc_001/users/rcali/test_environment/acqp.txt'
# index_path = '/autofs/space/nicc_001/users/rcali/test_environment/index.txt'

# Set file lists
pe_dir_1 = ["AP", "RL", "SI", "BlipA"]
pe_dir_2 = ["PA", "LR", "IS", "BlipP"]
ext = [".nii.gz", ".nii", ".bvec", ".bval", ".json"]
data_dir = ["orig", "mrtrix_files", "dtifit"]
trix_files = [".mif", "grad.b"]
qc_files = ["residuals.nii.gz", "noise.nii.gz"]
mif_files = ['preproc.mif', 'dwi_merged_PE.mif', 'noise.mif', 'preproc_b0_vols.mif']
dti_fit_files = [
    "dtifit_FA.nii.gz",
    "dtifit_L1.nii.gz",
    "dtifit_L2.nii.gz",
    "dtifit_L3.nii.gz",
    "dtifit_MD.nii.gz",
    "dtifit_MO.nii.gz",
    "dtifit_S0.nii.gz",
    "dtifit_V1.nii.gz",
    "dtifit_V2.nii.gz",
    "dtifit_V3.nii.gz",
]
file_list = []
dir_1_list = []
dir_2_list = []


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
    if in_file.endswith(".nii"):
        mrconvert = mrt3.MRConvert()
        mrconvert.inputs.in_file = in_file
        mrconvert.inputs.out_file = in_file + ".gz"
        mrconvert.inputs.args = "-force"
        mrconvert.run()
        os.remove(in_file)
        return


def file_check(output_file):
    if os.path.exists(output_file):
        print(output_file + " generated, continuing processing pipeline...")
    else:
        print(
            "Error generating ["
            + output_file
            + "], check last step for potential issues."
        )
        sys.exit()

#------------------------Initiate reverse phase-encoded data processing pipeline---------------------------
def reverse_pe_wf(direction_1, direction_2):

    mrconvert = mrt3.MRConvert()
    mrconvert.inputs.in_file = direction_1
    mrconvert.inputs.out_file = "dwi_PE_1.mif"
    mrconvert.inputs.args = (
        "-fslgrad dwi_PE_1.bvec dwi_PE_1.bval -json_import dwi_PE_1.json -force"
    )
    mrconvert.run()

    file_check("dwi_PE_1.mif")


    mrconvert = mrt3.MRConvert()
    mrconvert.inputs.in_file = direction_2
    mrconvert.inputs.out_file = "dwi_PE_2.mif"
    mrconvert.inputs.args = (
        "-fslgrad dwi_PE_2.bvec dwi_PE_2.bval -json_import dwi_PE_2.json -force"
    )
    mrconvert.run()



    file_check("dwi_PE_2.mif")

    mrcat = mrt3.MRCat()
    mrcat.inputs.in_files = ["dwi_PE_1.mif", "dwi_PE_2.mif"]
    mrcat.inputs.out_file = "dwi_merged_PE.mif"
    mrcat.inputs.args = "-force"
    mrcat.run()

    file_check("dwi_merged_PE.mif")

    mrconvert = mrt3.MRConvert()
    mrconvert.inputs.in_file = "dwi_merged_PE.mif"
    mrconvert.inputs.out_file = "dwi_merged_PE.mif"
    mrconvert.inputs.args = (
        "-export_pe_eddy index.txt config.txt -export_grad_fsl bvecs bvals -force"
    )
    mrconvert.run()

    file_check("dwi_merged_PE.mif")


    dwidenoise = mrt3.DWIDenoise()
    dwidenoise.inputs.in_file = "dwi_merged_PE.mif"
    dwidenoise.inputs.out_file = "denoised_merged_dwi.mif"
    dwidenoise.inputs.noise = "noise.mif"
    dwidenoise.inputs.args = "-extent 5 -force"
    dwidenoise.run()

    file_check("denoised_merged_dwi.mif")

    dwiextract = mrt3.DWIExtract()
    dwiextract.inputs.in_file = "denoised_merged_dwi.mif"
    dwiextract.inputs.bzero = True
    dwiextract.inputs.out_file = "b0_vols.mif"
    dwiextract.inputs.grad_fsl = ("bvecs", "bvals")
    dwiextract.run()

    file_check("b0_vols.mif")

    preproc = mrt3.DWIPreproc()
    preproc.inputs.in_file = "denoised_merged_dwi.mif"
    preproc.inputs.rpe_options = "header"
    preproc.inputs.out_file = "preproc.mif"
    preproc.inputs.eddy_options = (
        "--slm=linear --cnr_maps"  # linear second level model and replace outliers
    )
    preproc.inputs.args = "-nthreads 5 -eddyqc_all eddy -force"
    preproc.inputs.export_grad_mrtrix = True
    preproc.inputs.pe_dir = "AP"
    preproc.run()



def import_tables(config, index):
    
    mrconvert = mrt3.MRConvert()
    mrconvert.inputs.in_file = 'dwi_Pe_1.nii.gz'
    mrconvert.inputs.out_file = "dwi_PE_1.mif"
    mrconvert.inputs.args = (
        "-fslgrad dwi_PE_2.bvec dwi_PE_2.bval import_pe_eddy" + config + index + "-force"
    )
    mrconvert.run()

    file_check("dwi_PE_1.mif")


    mrconvert = mrt3.MRConvert()
    mrconvert.inputs.in_file = 'dwi_Pe_1.nii.gz'
    mrconvert.inputs.out_file = "dwi_PE_2.mif"
    mrconvert.inputs.args = (
        "-fslgrad dwi_PE_2.bvec dwi_PE_2.bval import_pe_eddy" + config + index + "-force"
    )
    mrconvert.run()



    file_check("dwi_PE_2.mif")



if args.config and args.index:
    config = args.config
    index = args.index
    import_tables(config, index)
    


#-------------------------Intitiate single phase-encoding direction processing stream----------------------------------

def single_pe_wf(subject_dir, in_file):
    for xx in os.listdir(subject_dir):
        if ext[0] in xx:
            mrconvert = mrt3.MRConvert()
            mrconvert.inputs.in_file = in_file
            mrconvert.inputs.out_file = "dwi_PE_1.mif"
            mrconvert.inputs.args = (
            "-fslgrad dwi_PE_1.bvec dwi_PE_1.bval -json_import dwi_PE_1.json -force"
            )
            mrconvert.run()
            
            dwidenoise = mrt3.DWIDenoise()
            dwidenoise.inputs.in_file = "dwi_PE_1.mif"
            dwidenoise.inputs.out_file = "denoise_dwi.mif"
            dwidenoise.inputs.noise = "noise.mif"
            dwidenoise.inputs.args = "-extent 5 -force"
            dwidenoise.run()

            degibbs = mrt3.MRDeGibbs()
            degibbs.inputs_in_file = 'denoise_dwi.mif'
            degibbs.inputs_out_file = 'denoise_degib_dwi.mif'
            degibbs.inputs_args = '-axes 0,1 -maxW 3 -minW 1 -nshifts 20'


# if single_pe_file.endswith(ext[1]):
#     print('Found unzipped, single phase-encoded image [ ' + single_pe_file + ' ], compressing image...')
#     mrzip(single_pe_file)
# elif single_pe_file.endswith(ext[0]):
#     print('Found zipped, single phase-encoded image [ ' + single_pe_file + ' ], proceeding with single phase-encoding processing pipeline...')

# single_pe_wf(path, single_pe_file)


def mif_to_nifti(mif_file_list):
    
    mrconvert = mrt3.MRConvert()
    mrconvert.inputs.in_file = mif_file_list[0]
    mrconvert.inputs.out_file = "preproc.mif"
    mrconvert.inputs.args = "-export_grad_fsl eddy_bvecs eddy_bvals -force"
    mrconvert.run()


    mrconvert = mrt3.MRConvert()
    mrconvert.inputs.in_file = mif_file_list[0]
    mrconvert.inputs.out_file = "preproc.nii.gz"
    mrconvert.run()


    mrconvert = mrt3.MRConvert()
    mrconvert.inputs.in_file = mif_file_list[1]
    mrconvert.inputs.out_file = "dwi_merged_PE.nii.gz"
    mrconvert.run()


    mrconvert = mrt3.MRConvert()
    mrconvert.inputs.in_file = mif_file_list[2]
    mrconvert.inputs.out_file = "noise.nii.gz"
    mrconvert.run()


    mrconvert = mrt3.MRConvert()
    mrconvert.inputs.in_file = mif_file_list[3]
    mrconvert.inputs.out_file = "preproc_b0_vols.nii.gz"
    mrconvert.run()


def tensor_fit(preproc_file):

    dwiextract = mrt3.DWIExtract()
    dwiextract.inputs.in_file = preproc_file
    dwiextract.inputs.bzero = True
    dwiextract.inputs.out_file = "preproc_b0_vols.mif"
    dwiextract.inputs.grad_fsl = ("eddy_bvecs", "eddy_bvals")
    dwiextract.run()


    bet = fsl.BET()
    bet.inputs.in_file = "preproc_b0_vols.nii.gz"
    bet.inputs.mask = True
    bet.inputs.frac = 0.15
    bet.ignore_exception = True
    bet.run()


    dti = fsl.DTIFit()
    dti.inputs.dwi = "preproc.nii.gz"
    dti.inputs.bvecs = "eddy_bvecs"
    dti.inputs.bvals = "eddy_bvals"
    dti.inputs.base_name = "dtifit"
    dti.inputs.mask = "preproc_b0_vols_brain_mask.nii.gz"
    dti.run()


def calc_residuals(in_file):
    maths = fsl.ImageMaths()
    maths.inputs.in_file = in_file
    maths.inputs.args = "-sub noise.nii.gz"
    maths.inputs.out_file = "residuals.nii.gz" 
    maths.run()

#denoise/degibb/mask/eddy/mask/dwigradcheck


organize_dir(data_dir[0])


orig_dir = os.path.join(path, "orig")


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

reverse_pe_wf('dwi_PE_1.nii.gz', 'dwi_PE_2.nii.gz')


mrtrix_dir = os.path.join(path, "mrtrix_files")

organize_dir(data_dir[1])

eddy_dir = os.path.join(path, 'eddy')

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
    file_check(dti_fit_file)

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

# sys.stdout.close()


