# Author: Shadi Zabad
# Date: Mar. 29th, 2018
# ---------------------
# This file is largely based on the following tutorials:
#
# https://f1000research.com/articles/5-1281
# https://www.bioconductor.org/help/course-materials/2014/BioC2014/minfi_BioC2014.pdf
#
# The tutorials will be quoted here to justify some of the steps in the workflow.
# ---------------------

# ---------------------
# Load relevant libraries:
# ---------------------

library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(plyr)
library(limma)

# ---------------------
# Define global variables and constants:
# ---------------------

# Constants
LR_PVAL_CUTOFF = 0.01 # LR -> Limma Regression
MWUT_PVAL_CUTOFF = 0.01 # MWUT -> Mann-Whitney U Test
ABS_BETA_DIFF = 0.1 

# Input paths

TARGETS_FILE <- "../output/0_processed_DNAm_data/m_targets.csv"
M_MATRIX_FILE <- "../output/0_processed_DNAm_data/m_m_matrix.csv"
BETA_MATRIX_FILE <- "../output/0_processed_DNAm_data/m_beta_matrix.csv"

# Output paths

OUTPUT_DIR <- "../output/1_differential_methylation_analysis/"

# ---------------------
# Workflow steps:
# ---------------------

# Load modified targets file:

print("=> Reading the targets file")

targets <- read.csv(TARGETS_FILE)

# Extract subject IDs for each category:

KMT2D_subjects <- targets[targets$Sample_Group == "KMT2D_LOF_discovery_cohort",
                          "Sample_Name"]
KMT2D_controls <- targets[targets$Sample_Group == "Control_for_KMT2D_LOF_discovery_cohort",
                          "Sample_Name"]

CHD7_subjects <- targets[targets$Sample_Group == "CHD7_LOF_discovery_cohort",
                         "Sample_Name"]
CHD7_controls <- targets[targets$Sample_Group == "Control_for_CHD7_LOF_discovery_cohort",
                         "Sample_Name"]

# -----
# Load M and beta matrices:

print ("=> Loading the M and beta matrices") 

fm_vals <- read.csv(M_MATRIX_FILE, row.names=1)
fbeta_vals <- read.csv(BETA_MATRIX_FILE, row.names=1)


##########################
# Probe-wise differential methylation analysis (limma regression):

print("==> Probe-wise differential methylation analysis (limma regression)...")

# Factor of interest
sample_group <- factor(targets$Sample_Group)

# Factors to account for: Age and Gender

patient_age <- factor(targets$Age)
patient_gender <- factor(targets$Gender)

design_mat <- model.matrix(~0+sample_group+patient_age+patient_gender, data=targets)

# Drop 2 individuals with missing age information (if considered):
fm_vals <- fm_vals[, !colnames(fm_vals) %in% c("control-74", "control-80")]

fit <- lmFit(fm_vals, design_mat)

contMatrix <- makeContrasts(sample_groupKMT2D_LOF_discovery_cohort-sample_groupControl_for_KMT2D_LOF_discovery_cohort,
                            sample_groupCHD7_LOF_discovery_cohort-sample_groupControl_for_CHD7_LOF_discovery_cohort,
                            levels=design_mat)

fit_contrast <- contrasts.fit(fit, contMatrix)
ebayes_fit <- eBayes(fit_contrast)

ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann450kSub <- ann450k[match(rownames(fm_vals),ann450k$Name),
                      c(1:4,12:19,24:ncol(ann450k))]

KMT2D_dmp <- topTable(ebayes_fit,
                      num=Inf,
                      p.value=LR_PVAL_CUTOFF,
                      coef=1,
                      genelist=ann450kSub,
                      adjust.method="fdr")

CHD7_dmp <- topTable(ebayes_fit,
                     num=Inf,
                     p.value=LR_PVAL_CUTOFF,
                     coef=2,
                     genelist=ann450kSub,
                     adjust.method="fdr")


##########################
# Non-parametric probe-wise differential methylation analysis (Mann-Whitney U test):

print("==> Non-parametric probe-wise differential methylation analysis (Mann-Whitney U test)...")

nonparam_test <- apply(fbeta_vals, 1,
                       function(x) {
                         pairwise.wilcox.test(x, targets$Sample_Group, p.adj="fdr")
                         })

nonparam_CHD7 <- data.frame(sapply(nonparam_test,
                                   function(x) {
                                     x$p.value["Control_for_CHD7_LOF_discovery_cohort",
                                               "CHD7_LOF_discovery_cohort"]
                                     }))

colnames(nonparam_CHD7) <- c("MW_U_pval")

# Filter to only keep p-values less than defined threshold
nonparam_CHD7 <- nonparam_CHD7[nonparam_CHD7$MW_U_pval < MWUT_PVAL_CUTOFF, , F]

nonparam_KMT2D <- data.frame(sapply(nonparam_test,
                                   function(x) {
                                     x$p.value["KMT2D_LOF_discovery_cohort",
                                               "Control_for_KMT2D_LOF_discovery_cohort"]
                                   }))

colnames(nonparam_KMT2D) <- c("MW_U_pval")

nonparam_KMT2D <- nonparam_KMT2D[nonparam_KMT2D$MW_U_pval < MWUT_PVAL_CUTOFF, , F]


##########################
# Probe-wise differential methylation analysis using absolute difference:

print("==> Probe-wise differential methylation analysis using absolute difference...")

KMT2D_beta_diff <- data.frame(rowMeans(fbeta_vals[, KMT2D_subjects], na.rm=TRUE) -
  rowMeans(fbeta_vals[, KMT2D_controls], na.rm=TRUE))

colnames(KMT2D_beta_diff) <- c("beta_diff")

KMT2D_beta_diff <- KMT2D_beta_diff[abs(KMT2D_beta_diff$beta_diff) > ABS_BETA_DIFF, , F]

CHD7_beta_diff <- data.frame(rowMeans(fbeta_vals[, CHD7_subjects], na.rm=TRUE) -
  rowMeans(fbeta_vals[, CHD7_controls], na.rm=TRUE))

colnames(CHD7_beta_diff) <- c("beta_diff")

CHD7_beta_diff <- CHD7_beta_diff[abs(CHD7_beta_diff$beta_diff) > ABS_BETA_DIFF, , F]


##########################
# Matching significant probes from the different criteria:

print("==> Matching significant probes from the different criteria...")

CHD7_dmp$probe_id <- rownames(CHD7_dmp)
nonparam_CHD7$probe_id <- rownames(nonparam_CHD7)
CHD7_beta_diff$probe_id <- rownames(CHD7_beta_diff)

CHD7_sig_probes_beta <- join_all(list(CHD7_dmp,
                                      nonparam_CHD7,
                                      CHD7_beta_diff),
                                 by="probe_id",
                                 type = "inner")

KMT2D_dmp$probe_id <- rownames(KMT2D_dmp)
nonparam_KMT2D$probe_id <- rownames(nonparam_KMT2D)
KMT2D_beta_diff$probe_id <- rownames(KMT2D_beta_diff)

KMT2D_sig_probes_beta <- join_all(list(KMT2D_dmp,
                                      nonparam_KMT2D,
                                      KMT2D_beta_diff),
                                 by="probe_id",
                                 type = "inner")


##########################
# Output differentially methylated probes:

print("=> Outputting files for differentially methylated probes:")

# For CHD7:

write.csv(CHD7_sig_probes_beta, paste0(OUTPUT_DIR, "CHD7_sig_probes.csv"))
write.csv(fbeta_vals[rownames(fbeta_vals) %in% CHD7_sig_probes_beta$probe_id, ], paste0(OUTPUT_DIR, "CHD7_beta_matrix.csv"))
write.csv(fbeta_vals[rownames(fm_vals) %in% CHD7_sig_probes_beta$probe_id, ], paste0(OUTPUT_DIR, "CHD7_M_matrix.csv"))

# For KMT2D:

write.csv(KMT2D_sig_probes_beta, paste0(OUTPUT_DIR, "KMT2D_sig_probes.csv"))
write.csv(fbeta_vals[rownames(fbeta_vals) %in% KMT2D_sig_probes_beta$probe_id, ], paste0(OUTPUT_DIR, "KMT2D_beta_matrix.csv"))
write.csv(fbeta_vals[rownames(fm_vals) %in% KMT2D_sig_probes_beta$probe_id, ], paste0(OUTPUT_DIR, "KMT2D_M_matrix.csv"))


##########################
# Checking results against those published in Butcher et al. 2017

KMT2D_ref_probes <- read.csv("../metadata/Reference_KMT2D_Diff_Methylated.csv")

CHD7_ref_probes <- read.csv("../metadata/Reference_CHD7_Diff_Methylated.csv")

print("=> Total number of reference KMT2D probes:")

print(nrow(KMT2D_ref_probes))

print("=> Number of matching KMT2D probes:")

print(sum(!is.na(match(KMT2D_ref_probes$DIFF_METH_PROBES,
                       KMT2D_sig_probes_beta$probe_id))))

#############

print("=> Total number of reference CHD7 probes:")

print(nrow(CHD7_ref_probes))

print("=> Number of matching CHD7 probes:")

print(sum(!is.na(match(CHD7_ref_probes$DIFF_METH_PROBES,
                       CHD7_sig_probes_beta$probe_id))))



