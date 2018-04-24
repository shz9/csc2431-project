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
library(limma)
library(DMRcate)
library(RColorBrewer)

# ---------------------
# Define global variables and constants:
# ---------------------

# Constants

EXCLUDE_VALIDATION_SUBJECTS <- FALSE # If you change this, remember to change the OUTPUT directories!
EXCLUDE_POOR_SAMPLES <- TRUE
GENERATE_QUALITY_REPORT <- TRUE
NORMALIZATION_METHOD <- "Illumina"
DECILE_DIFF_CUTOFF <- 0.1

# Input files and directories
METHEXP_FILES_DIR <- "../data"
METHEXP_SHEET <- "GSE97362_methexp_sheet.csv"

NON_SPECIFIC_PROBES_FILE <- "../metadata/48639-non-specific-probes-Illumina450k.csv"
MULTI_MAP_PROBES_FILE <- "../metadata/HumanMethylation450_15017482_v.1.1_hg19_bowtie_multimap.txt"

# Output files and directories

QUALITY_REPORT_PATH <- "../quality_report/GSE97362_qc.pdf"
OUTPUT_TARGETS_PATH <- "../output/0_processed_DNAm_data/validation/m_targets.csv"
OUTPUT_BETA_MATRIX <- "../output/0_processed_DNAm_data/validation/m_beta_matrix.csv"
OUTPUT_M_MATRIX <- "../output/0_processed_DNAm_data/validation/m_m_matrix.csv"
OUTPUT_DIR <- "../output/"
PLOTS_DIR <- "../plots/"

# Color Pallette

pal <- brewer.pal(8,"Dark2")

# ---------------------
# Plotting functions:
# ---------------------

generate_density_plot <- function(file_name, raw_data, normalized_data) {
  
  png(paste0(PLOTS_DIR, file_name))
  
  par(mfrow=c(1,2))
  densityPlot(raw_data, sampGroups=targets$Sample_Group, main="Raw", legend=FALSE)
  legend("top", legend = levels(factor(targets$Sample_Group)),
         text.col=brewer.pal(8, "Dark2"))
  densityPlot(getBeta(normalized_data), sampGroups=targets$Sample_Group,
              main="Normalized", legend=FALSE)
  legend("top", legend = levels(factor(targets$Sample_Group)),
         text.col=brewer.pal(8, "Dark2"))
  
  dev.off()
  
}

generate_mds_plot <- function(file_name, grs_data) {
  
  png(paste0(PLOTS_DIR, file_name))
  
  plotMDS(getM(grs_data), top=1000, gene.selection="common",
          col=pal[factor(targets$Sample_Group)])
  legend("top", legend=levels(factor(targets$Sample_Group)), text.col=pal,
         bg="white")
  
  dev.off()
  
}

# ---------------------
# Auxiliary functions:
# ---------------------

decile_range <- function (x) { 
  d <- quantile(x, c(.1, .9))
  d["90%"] - d["10%"] 
}

# ---------------------
# Workflow steps:
# ---------------------

##########################
# Read metharray data:

print("==> Reading data...")

targets <- read.metharray.sheet(METHEXP_FILES_DIR, 
                                pattern=METHEXP_SHEET)

# Remove spaces from names:
targets$Sample_Group <- gsub(" ", "_", targets$Sample_Group)
targets$Disease_State <- gsub(" ", "_", targets$Disease_State)

if (EXCLUDE_VALIDATION_SUBJECTS) {
    # Filter targets to keep discovery cohort only:
    targets <- targets[grepl("discovery", targets$Sample_Group),]
}

metharray_dat <- read.metharray.exp(targets=targets)

sampleNames(metharray_dat) <- targets$Sample_Name

##########################
# Quality control:
print("==> Quality control...")

if (GENERATE_QUALITY_REPORT) {
  qcReport(metharray_dat, pdf=QUALITY_REPORT_PATH)
}

detP <- detectionP(metharray_dat)

if (EXCLUDE_POOR_SAMPLES) {
  keep <- colMeans(detP) < 0.05
  metharray_dat <- metharray_dat[, keep]
}

##########################
# Preprocessing and normalization:
print("==> Preprocessing and normalization...")

if (NORMALIZATION_METHOD == "quantile") {
  # "preprocessQuantile function is more suited 
  # for datasets where you do not expect global differences between your samples, 
  # for example a single tissue."
  grs <- preprocessQuantile(metharray_dat)
} else if (NORMALIZATION_METHOD == "raw") {
  grs <- preprocessRaw(metharray_dat)
} else {
  # Illumina preprocessing is the default
  mset <- preprocessIllumina(metharray_dat, normalize = "controls")
  rset <- ratioConvert(mset, what = "both", keepCN = T)
  grs <- mapToGenome(rset)
  grs <- addSnpInfo(grs)
}

##########################
# Generate density & MDS plots:
print("==> Generating density plot...")
generate_density_plot("beta_vals_density.png", metharray_dat, grs)

print("==> Generating MDS plot...")
generate_mds_plot("MDS_plot.png", grs)

##########################
# Filtering:
print("==> Filtering...")

# Poor performing probes:

detP <- detP[match(featureNames(grs), rownames(detP)), ]
#poor_perf_probes <- rowSums(detP < 0.01) != ncol(grs)

poor_perf_probes <- c()

# Cross-reactive/non-specific probes
cross_react <- read.csv(NON_SPECIFIC_PROBES_FILE, head = T, as.is = T)
cross_react_probes <- as.character(cross_react$TargetID)

# BOWTIE2 multi-mapped probes
multi_map <- read.csv(MULTI_MAP_PROBES_FILE, head = F, as.is = T)
multi_map_probes <- as.character(multi_map$V1)

# X/Y probes
grs_ann <- getAnnotation(grs)
XY_probes <- rownames(grs_ann)[which(grs_ann$chr == "chrX" | 
                                       grs_ann$chr == "chrY")]

# Determine unique probes
filter_probes <- unique(c(poor_perf_probes, 
                          cross_react_probes, 
                          multi_map_probes, 
                          XY_probes))

grs <- grs[!rownames(grs) %in% filter_probes, ]

grs <- dropLociWithSnps(grs, maf=0.01)

##########################
# Generate MDS plots after filtering:

print("==> Generating MDS plot...")
generate_mds_plot("MDS_plot_filt.png", grs)

##########################
# Obtain M and beta values:

print("==> Generating Beta and M matrices...")

# "[D]ue to their distributional properties, M-values are more 
# appropriate for statistical testing"

m_vals <- getM(grs)

# "Beta values are generally preferable for describing the 
# level of methylation at a locus or for graphical presentation 
# because percentage methylation is easily interpretable."

beta_vals <- getBeta(grs)

##########################
# Filter beta and M matrices:

print("==> Filtering Beta and M matrices...")

# Beta matrix
nbeta_vals <- rmSNPandCH(beta_vals, 
                         dist = 5, 
                         mafcut = 0.01, 
                         and = T, 
                         rmcrosshyb = T, 
                         rmXY = T)

fbeta_vals <- nbeta_vals[complete.cases(nbeta_vals) & is.finite(rowSums(nbeta_vals)), ]

# Remove probes where the top deciles (10% and 90% differ by less than 0.1):
beta_decile_range <- apply(fbeta_vals, 1, decile_range)
fbeta_vals <- fbeta_vals[beta_decile_range > DECILE_DIFF_CUTOFF, ]

# M matrix
nm_vals <- rmSNPandCH(m_vals,
                      dist = 5,
                      mafcut = 0.01,
                      and = T,
                      rmcrosshyb = T,
                      rmXY = T)

fm_vals <- nm_vals[complete.cases(nm_vals) & is.finite(rowSums(nm_vals)), ]

# Ensure that the same probes are in the M and beta matrices:
# (Some probes will have Inf M values, and these will be removed.
# So, we need to make sure they're removed from both matrices.)

fbeta_vals <- fbeta_vals[rownames(fbeta_vals) %in% rownames(fm_vals), ]
fm_vals <- fm_vals[rownames(fm_vals) %in% rownames(fbeta_vals), ]

print("==> Dimensions of probes in Beta/M matrices...")
print(dim(fbeta_vals))

# Output of preprocessing:

# modified targets file:

write.csv(targets, OUTPUT_TARGETS_PATH, row.names=F)

# beta matrix:

write.csv(fbeta_vals, OUTPUT_BETA_MATRIX)

# M matrix:

write.csv(fm_vals, OUTPUT_M_MATRIX)

print("=> Done!")

