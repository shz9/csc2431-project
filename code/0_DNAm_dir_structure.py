"""
Author: Shadi Zabad
Date: March 31th, 2018

"""

import os
import glob
from shutil import copyfile
import errno

# Constants and global variables:
# -------------------------------
IDAT_FILES_DIR = "/valr/tahmid/csc2431/GSE97362"
OUTPUT_DIR = "/valr/szabad/csc2431_project/data"

# Move .idat files:
for idat_f in glob.glob(os.path.join(IDAT_FILES_DIR, "*.idat")):

    idat_basename = os.path.basename(idat_f)
    sub_dir = idat_basename.split("_")[1]

    try:
        os.makedirs(os.path.join(OUTPUT_DIR, sub_dir))
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    copyfile(idat_f, os.path.join(OUTPUT_DIR, sub_dir, idat_basename))

