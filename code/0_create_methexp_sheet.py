"""
Author: Shadi Zabad
Date: March 30th, 2018

Patient annotations table was extracted from:
https://www.ncbi.nlm.nih.gov/geo/geo2r/?acc=GSE97362

Access date: 2018-03-30
"""

import pandas as pd
import os
import glob

# Constants and global variables:
# -------------------------------
PATIENT_ANNOTATIONS_FILE = "C:/Users/Moses Lab/Code/csc2431/patient_profiles_GSE97362.csv"
METHEXP_SHEET_FILE = "C:/Users/Moses Lab/Code/csc2431/GSE97362_methexp_sheet.csv"
IDAT_FILES_DIR = "C:/Users/Moses Lab/Documents/Data/DNAm_csc2431/idat/5723646052"

OUTPUT_COLUMNS = ["Sample_Name", "Sample_Well", "Sample_Plate",
                  "Sample_Group", "Pool_ID", "Sentrix_ID", "Sentrix_Position",
                  "GSM_ID", "Age", "Gender", "Disease_State"]

# Auxiliary functions:
# -------------------------------


def parse_age(age_str):
    try:
        return float(age_str.split(":")[-1].strip())
    except ValueError:
        return None


def parse_title(title_str):
    return title_str.split()[0].strip()


def get_sentrix_id(gsm_id):
    return idat_dict[gsm_id]["Sentrix_ID"]


def get_sentrix_pos(gsm_id):
    return idat_dict[gsm_id]["Sentrix_Position"]


# Retrieve the .idat files and extract Sentrix ID and position:
# -------------------------------

idat_dict = dict()

for idat_f in glob.glob(os.path.join(IDAT_FILES_DIR, "*.idat")):

    idat_filename = os.path.basename(idat_f)
    gsm_id, sentrix_id, sentrix_pos = idat_filename.split("_")[:3]
    idat_dict[gsm_id] = {
        "Sentrix_ID": sentrix_id,
        "Sentrix_Position": sentrix_pos
    }


# Parse the annotations file:
# -------------------------------

# Load patient annotations file:
patient_annotations = pd.read_csv(PATIENT_ANNOTATIONS_FILE, header=0)

# Transforming the dataframe:

patient_annotations["Age"] = patient_annotations["Characteristics"].apply(parse_age)
patient_annotations["Sample_Name"] = patient_annotations["Title"].apply(parse_title)

patient_annotations["Sample_Group"] = patient_annotations["Sample type"]
patient_annotations["GSM_ID"] = patient_annotations["Accession"]
patient_annotations["Disease_State"] = patient_annotations["Disease state"]

patient_annotations["Sample_Well"] = None
patient_annotations["Sample_Plate"] = None
patient_annotations["Pool_ID"] = None

patient_annotations["Sentrix_ID"] = patient_annotations["Accession"].apply(get_sentrix_id)
patient_annotations["Sentrix_Position"] = patient_annotations["Accession"].apply(get_sentrix_pos)

# Reordering and excluding unnecessary columns:
# -------------------------------
patient_annotations = patient_annotations[OUTPUT_COLUMNS]

# Outputting sample sheet
# -------------------------------
patient_annotations.to_csv(METHEXP_SHEET_FILE, index=None)
