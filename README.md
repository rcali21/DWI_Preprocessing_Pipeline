## Diffusion Processing Pipeline ##

This pipeline leverages functions from the nipype library to process both single phase encoding and reverse phase encoding acquired diffusion weighted MRI images.

Command line tool dependencies:

* FSL: https://fsl.fmrib.ox.ac.uk/fsl/fslwiki
* Mrtrix: https://www.mrtrix.org/

Nipype: https://nipype.readthedocs.io/en/latest/


### Important: ###
Some of the nipype Mrtrix functions (MRCat, etc.) are awaiting push to master so they will have to be manually implemented until they can be pulled from the nipype repo.

