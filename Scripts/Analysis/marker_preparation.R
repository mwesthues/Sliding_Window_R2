#!/usr/bin/Rscript --vanilla --slave
source("./Scripts/Functions/Marker_Functions.R")
source("./Scripts/Functions/Misc_Functions.R")
pkgTest("reshape2")
pkgTest("data.table")
pkgTest("Matrix")
pkgTest("dplyr")
pkgTest("ggplot2")


# PREPARATION OF PHENOTYPIC DATA ------------------------------------------
load(file = "./Data/genotypes")
init_hybrids <- gsub(pattern = "DF_", replacement = "", x = genotypes)
# Remove hybrids whose genotypes were not genotyped.
no_marker <- c("812", "1309", "1310", "1311", "1312")
ind_dent <- which(unlist(lapply(X = strsplit(x = init_hybrids, split = "_"), 
                                FUN = "[[", 1)) %in% no_marker)
ind_flint <- which(unlist(lapply(X = strsplit(x = init_hybrids, split = "_"), 
                                 FUN = "[[", 2)) %in% no_marker)
ind <- union(ind_dent, ind_flint)
kernel_df <- kernel_df[-ind, ]

# Replace the old genotype (prior to "=" sign) with the new genotype (after
# the "=" sign). The latter genotypes are the subsequent generation of the 
# former genotypes (S5 generation). (Email from TAS on 2015-01-09)
repl_marker <- c("1307" = "387", "1308" = "394", "1313" = "430", 
                 "1314" = "429", "4508" = "548", "4509" = "545",
                 "4510" = "724", "4511" = "725", "4512" = "727", 
                 "4513" = "431")
dent <- unlist(lapply(X = strsplit(as.character(kernel_df$Hybrid),
                                   split = "_"), FUN = "[[", 1))
dent_repl <- unname(repl_marker[dent[dent %in% names(repl_marker)]])
dent[dent %in% names(repl_marker)] <- dent_repl
flint <- unlist(lapply(X = strsplit(as.character(kernel_df$Hybrid), 
                                    split = "_"), FUN = "[[", 2))
flint_repl <- unname(repl_marker[flint[flint %in% names(repl_marker)]])
flint[flint %in% names(repl_marker)] <- flint_repl
hybrid <- paste(dent, flint, sep = "_")
kernel_df$Hybrid <- as.factor(hybrid)
kernel_df$Dent <- as.factor(dent)
kernel_df$Flint <- as.factor(flint)
kernel_df <- kernel_df[c("Hybrid", "Dent", "Flint", "KTM", "KTS")]
dent_parents <- unique(dent)
flint_parents <- unique(flint)

# Matrix with phenotypic data.
kernel_matrix <- as.matrix(kernel_df[, c("KTM", "KTS")])
rownames(kernel_matrix) <- as.character(kernel_df$Hybrid)


# MARKER MATRIX PREPARATION -----------------------------------------------
# Load data frame with meta information on marker names and their positions in 
# the genome.
marker_map <- read.table("./Server/marker_map.txt", header = TRUE, sep = "\t")
marker_map <- marker_map[complete.cases(marker_map), ]
marker_map <- marker_map[marker_map$chromosome != 0, ]
map_marker_names <- as.character(marker_map[, "markername"])
map_marker_names <- gsub(pattern= "[.]", replacement= "_", x= map_marker_names)
map_marker_names <- gsub(pattern= "[-]", replacement= "_", x= map_marker_names)
marker_map$markername <- map_marker_names

######################### PREPARING MARKER DATA ###############################
# Load the marker data.
x_marker <- read.table(file = "./Data/marker_sample_matrix_FWD.txt", 
                       header = TRUE, sep = "\t")

# Load the meta data on samples from the UHOH database.
x_sample <- read.table(file = "./Data/collect_samples.txt",
                       header = TRUE, sep = "\t", stringsAsFactors = FALSE)

x_sample$id.GTP <- as.factor(x_sample$id.GTP)
x_sample$GTP.DXM.id <- as.factor(x_sample$GTP.DXM.id)
names(x_sample)[names(x_sample) == 'id.GTP'] <- 'GTP_ID'

# Names for sample ID conversion.
smp_conv <- read.table(file = "./Data/uhoh_smp_id.txt", header = TRUE, 
                       sep = "\t", stringsAsFactors = FALSE)
colnames(x_marker)[2:ncol(x_marker)] <- smp_conv[match(colnames(x_marker)[
    2:ncol(x_marker)], smp_conv$GTP.DXM.id), "id.GTP"]

unique_lines <- unique(c(dent, flint))
present_genos <- colnames(x_marker)[!colnames(x_marker) %in% "markername"]
# Stop if not any parental inbred line is genotyped.
stopifnot(sum(!unique_lines %in% present_genos) == 0)



#################### SPECIFYING MARKER SELECTION CRITERIA #####################
marker_origin_excl <- NULL
callfreq <- 0.95 # remove if callfreq <= value
rating <- 3
smp_callrate <- 0
maf <- 0.05 # remove if maf <= value
snp_het <- 0.05 # remove if snp_het > 0.05


############################## SNP QUALITY CHECKS ############################
pool_pass <- droplevels(with(x_sample, 
                             x_sample[GTP_ID %in% c(flint_parents,
                                                    dent_parents), ]))

## Keep all samples equal to or smaller than the rating threshold.
rating_pass <- as.character(pool_pass[pool_pass$rating <= rating, 
                                      "GTP_ID"])

# All marker origins in the UHOH database.
marker_origin <- c('PZE', 'SYN', 'PTU', 'PZA', 'ZM0', 'PZD', 'PZB', 'PHM',
                   'MISC')

if(is.null(marker_origin_excl))
{
    # All markers are included.
    marker_origin_incl <- marker_origin
    
} else
{
    ### Remove markers from sets, which shall not be included.
    unwanted <- unique(grep(pattern = paste(marker_origin_excl,
                                            collapse = "|"),
                            x = x_marker[, "markername"], value = TRUE))
    
    # Markers from origins, which shall be kept.
    x_marker <- x_marker[! as.character(x_marker$markername) %in% unwanted, ]
    
    # Which marker origins are included in the output?
    marker_origin_incl <- marker_origin[!marker_origin %in% marker_origin_excl]
}


# Convert DF with markers to matrix.
init_markmat <- as.matrix(x_marker[, colnames(x_marker) %in% rating_pass])

# Extract marker names and convert special character to '_'.
mnames <- as.character(x_marker[, 1])
mnames <- gsub(pattern = "[.]", replacement = "_", x = mnames)
mnames <- gsub(pattern = "[-]", replacement = "_", x = mnames)

# Add marker names
rownames(init_markmat) <- mnames

# Change all NAs to '?'.
init_markmat[is.na(init_markmat)] <- "??"
init_markmat[init_markmat == "XX"] <- "??"
init_markmat <- t(init_markmat)
dent_init_markmat <- init_markmat[rownames(init_markmat) %in% dent_parents, ]
flint_init_markmat <- init_markmat[rownames(init_markmat) %in% flint_parents, ]

# Call frequency
call_dent_names <- callfreq.fun(x = dent_init_markmat, output = "markerNames",
                                callThresh = callfreq)
call_flint_names <- callfreq.fun(x = flint_init_markmat, 
                                 output = "markerNames", callThresh = callfreq)
dent_cfreq_mat <- dent_init_markmat[, call_dent_names]
flint_cfreq_mat <- flint_init_markmat[, call_flint_names]


# Remove markers without any information as well as markers that are 
# monomorphic and determine the reference allele for each marker.
poly_dent_names <- maf.fun(x = dent_cfreq_mat, output = "markerNames",
                           mafThresh = maf)
poly_flint_names <- maf.fun(x = flint_cfreq_mat, output = "markerNames",
                            mafThresh = maf)

# Remove markers with high heterozygosity levels.
het_dent_names <- heterozygosity.fun(x = dent_cfreq_mat, 
                                     output = "markerNames",
                                     hetThresh = snp_het)
het_flint_names <- heterozygosity.fun(x = flint_cfreq_mat, 
                                      output = "markerNames",
                                      hetThresh = snp_het)


# Combine the markers, which fullfill the quality check criteria in each group.
crit_dent_names <- intersect(poly_dent_names, het_dent_names)
crit_flint_names <- intersect(poly_flint_names, het_flint_names)
crit_names <- union(crit_dent_names, crit_flint_names)
poly_dent <- init_markmat[dent_parents, crit_dent_names]
poly_flint <- init_markmat[flint_parents, crit_flint_names]

# Get mutual major and minor genotypes for each locus in the combined Dent
# and Flint data.
geno_dent_list <- maf.fun(x = poly_dent, output = "genoList")
geno_flint_list <- maf.fun(x = poly_flint, output = "genoList")

# Recode dent genotypes.
dent_reco <- recode.fun(x = poly_dent, 
                        major = geno_dent_list[["major_allele"]],
                        minor = geno_dent_list[["minor_allele"]],
                        major_coding = "XX", minor_coding = "YY", 
                        het_coding = "XY", na_coding = "??")
# Recode flint genotypes.
flint_reco <- recode.fun(poly_flint, 
                         major = geno_flint_list[["major_allele"]],
                         minor = geno_flint_list[["minor_allele"]],
                         major_coding = "XX", minor_coding = "YY", 
                         het_coding = "XY", na_coding = "??")


######################## IMPUTATION (PREPARATION) ############################
### 1) Preparation of map data
dent_map <- marker_map[marker_map$markername %in% colnames(dent_reco), ]
flint_map <- marker_map[marker_map$markername %in% colnames(flint_reco), ]

### 2) Preparation of marker data
# Initial marker matrix for heterotic groups.
# Sort matrices and check congruency.
dent_marker <- t(dent_reco)
dent_marker <- dent_marker[dent_map$markername, ]
identical(rownames(dent_marker), dent_map$markername)

flint_marker <- t(flint_reco)
flint_marker <- flint_marker[flint_map$markername, ]
identical(rownames(flint_marker), flint_map$markername)

# Store all map and marker matrices in one list.
m_list <- list()
m_list[["Dent"]][["marker"]] <- dent_marker 
m_list[["Dent"]][["map"]] <- dent_map
m_list[["Flint"]][["marker"]] <- flint_marker 
m_list[["Flint"]][["map"]] <- flint_map




################################ BEAGLE INPUT FILES ###########################
type <- c("Dent", "Flint")
chromosome <- seq_len(10)

for(i in type){
    
    dta_raw <- m_list[[i]][["marker"]]
    map_raw <- m_list[[i]][["map"]]
    map_raw <- map_raw[order(map_raw$position), ]
    dta_raw <- dta_raw[map_raw$markername, ]    
    names(map_raw)[c(3, 4)] <- c("chr", "pos")
    rownames(map_raw) <- map_raw$markername
    map_raw <- map_raw[, c("chr", "pos")]
    
    for(j in chromosome){
        
        # Marker file of chromosome j
        dta_chr <- dta_raw[map_raw$chr == j, ]    
        # Map file of chromosome j
        map_chr <- map_raw[map_raw$chr == j, ]        
        # Data translation for each chromosome
        y <- c()
        
        for(l in 1:ncol(dta_chr)){
            
            x <- na.omit(dta_chr[, l])
            x1 <- unlist(lapply(X = strsplit(x, ""), FUN = "[[", 1))
            x2 <- unlist(lapply(X = strsplit(x, ""), FUN = "[[", 2))
            y <- cbind(y, x1, x2)
            stopifnot(!any(is.na(y)))
        }
        
        # Generate the first row of the matrix required for Beagle.
        line_names <- rep(colnames(dta_chr), each = 2)
        ifelse(identical(line_names[(1:length(line_names)) %% 2 == 0], 
                         colnames(dta_chr)),
               yes = y1 <- rbind(line_names, y),
               no = warning("Lack of congruency"))
        # Generate the first column required for Beagle.
        I <- c("I", rep("M", nrow(dta_chr)))
        # Generate the second columnn required for Beagle.
        mk_names <- c("id", rownames(dta_chr))
        # Produce the final Beagle input file.
        y2 <- cbind(I, mk_names, y1)        
        map_chr$alleleA <- rep("X", times = nrow(map_chr))
        map_chr$alleleB <- rep("Y", times = nrow(map_chr))
        map_chr$chr <- NULL
        map_chr <- as.matrix(map_chr)
        colnames(map_chr) <- NULL
        # Output files (== Beagle input files)
        if(identical(rownames(map_chr), as.vector(y2[-1, "mk_names"])))
        {        
            # Marker file
            write.table(y2, paste0("./Derived/SNP_Data/kernel_data/all_lines/", 
                                   "Imputed_SNPs/input",
                                   "/", get("i"), "_chr", 
                                   get("j"), ".txt"),
                        sep = "\t", quote = FALSE, col.names = FALSE, 
                        row.names = FALSE)
            
            # Map file
            write.table(map_chr, 
                        paste0("./Derived/SNP_Data/kernel_data/all_lines/", 
                               "Imputed_SNPs/input",
                               "/", get("i"), "_map_chr",
                               get("j"), ".txt"), sep = "\t", quote = FALSE,
                        col.names = FALSE, row.names = TRUE)
        }
        message(paste(i, paste0("chr_", j)))
    }
} 




############################## BEAGLE OUTPUT FILES ############################
# Extract the names of all Beagle input files.
all_files <- list.files(path = paste0("./Derived/SNP_Data/kernel_data/", 
                                      "all_lines/Imputed_SNPs/input"))
type <- c("Dent", "Flint")
chr <- paste0("chr", seq_len(10))

# Make sure to adjust the 'Xmx...m' number, depending on available memory.
# Create Beagle output.
for(i in type){
    for(j in chr){
        
        chr_files <- all_files[!grepl(pattern = "map", x = all_files)]
        map_files <- all_files[grepl(pattern = "map", x = all_files)]
        
        # Browse the list of geno-files and map-files with matching names,
        # then select the matching geno and map files and impute missing
        # values.
        if(j == "chr1"){
            chr_files <- chr_files[!grepl(pattern = "chr10", x = chr_files)]
            map_files <- map_files[!grepl(pattern = "chr10", x = map_files)]
        }
        
        chr_matches <- sapply(X = c(get("i"), get("j")),
                              FUN = grep, fixed = TRUE, chr_files)
        
        perf_chr_match <- intersect(chr_matches[[1]], chr_matches[[2]])
        
        
        map_matches <- sapply(X = c(get("i"), get("j")),
                              FUN = grep, map_files)
        perf_map_match <- intersect(map_matches[[1]], map_matches[[2]])
        
        system(paste("java -Xmx15000m -jar", 
                     paste0(getwd(), "/Data/", "beagle.jar"), 
                     paste0("unphased=", getwd(),
                            "/Derived/SNP_Data/kernel_data/all_lines/", 
                            "Imputed_SNPs/input/",
                            chr_files[perf_chr_match]),
                     paste0("markers=", getwd(),
                            "/Derived/SNP_Data/kernel_data/all_lines/", 
                            "Imputed_SNPs/input/",
                            map_files[perf_map_match]),
                     "missing=?", 
                     paste0("out=", getwd(),
                            "/Derived/SNP_Data/kernel_data/all_lines/", 
                            "Imputed_SNPs/output/out"),
                     "niterations=25 nsamples=20"))
    }
}


########################## HETEROZYGOUS BEAGLE OUTPUT ########################
# Merge alleles from Beagle output to genotypes for each marker in order to be
# able of determining the reference alleles.
phased_gz <- list.files(path = paste0("./Derived/SNP_Data/kernel_data/", 
                                      "all_lines/Imputed_SNPs/output"))
phased_gz <- phased_gz[grepl(pattern = "phased", x = phased_gz)]

for(i in phased_gz){
    system(paste("gzip -d -f",
                 paste0(getwd(), "/Derived/SNP_Data/kernel_data/all_lines/",
                        "Imputed_SNPs/output/", i)))
}

phased_unz <- list.files(path = paste0("./Derived/SNP_Data/kernel_data/", 
                                       "all_lines/Imputed_SNPs/output"))
phased_unz <- phased_unz[grepl(pattern = "phased", x = phased_unz)]

# Merge the two alleles at each locus, thereby converting them back to 
# genotypes.
for(i in phased_unz){
    phased <- read.table(paste0(getwd(), "/Derived/SNP_Data/kernel_data/", 
                                "all_lines/Imputed_SNPs/output/", i),
                         check.names = FALSE, header = TRUE,
                         stringsAsFactors = FALSE)
    rownames(phased) <- phased$id
    phased$I <- NULL
    phased$id <- NULL
    phased <- as.matrix(phased)
    oddvals <- seq(from = 1, to = ncol(phased), by = 2)
    evenvals <- seq(from = 2, to = ncol(phased), by = 2)
    oddMat <- phased[, oddvals] 
    evenMat <- phased[, evenvals]
    mPhased <- matrix(paste0(oddMat, evenMat), 
                      nrow = nrow(phased), ncol = (0.5 * ncol(phased)),
                      dimnames = list(rownames(phased), colnames(evenMat)))
    outname <- get("i")
    outname <- gsub(pattern = ".txt.phased", replacement = "", x = outname)
    outname <- gsub(pattern = "out.", replacement = "", x = outname)
    outname <- paste0("phase.", outname)    
    write.table(x = mPhased, 
                file = paste0(getwd(), "/Derived/SNP_Data/kernel_data/", 
                              "all_lines/Imputed_SNPs/", "output/", outname))   
}

# Dent (chromosome concatenation)
chromosome <- seq_len(10)
dentList <- list()
for(i in chromosome){
    impDf <- read.table(file = paste0("./Derived/SNP_Data/kernel_data/", 
                                      "all_lines/Imputed_SNPs/output/", 
                                      "phase.Dent_chr", i),
                        stringsAsFactors = FALSE, header = TRUE)
    dentList[[i]]  <- impDf
}
dentMat <- as.matrix(do.call(rbind, dentList))
colnames(dentMat) <- gsub(pattern = "X", replacement = "",
                          x = colnames(dentMat))
dent_mat <- t(dentMat[, dent_parents])

# Preparation for 'maf.fun' function.
dent_mat[dent_mat == "XX"] <- "AA"
dent_mat[dent_mat == "XY"] <- "AB"
dent_mat[dent_mat == "YY"] <- "BB"

# Flint (chromosome concatenation)
chromosome <- seq_len(10)
flintList <- list()
for(i in chromosome){
    impDf <- read.table(file = paste0("./Derived/SNP_Data/kernel_data/", 
                                      "all_lines/Imputed_SNPs/output/", 
                                      "phase.Flint_chr", i),
                        stringsAsFactors = FALSE, header = TRUE)
    flintList[[i]]  <- impDf
}
flintMat <- as.matrix(do.call(rbind, flintList))
colnames(flintMat) <- gsub(pattern = "X", replacement = "",
                           x = colnames(flintMat))
flint_mat <- t(flintMat[, flint_parents])

# Preparation for 'maf.fun' function.
flint_mat[flint_mat == "XX"] <- "AA"
flint_mat[flint_mat == "XY"] <- "AB"
flint_mat[flint_mat == "YY"] <- "BB"


# Marker matrices with markers that pass the quality check criteria in both
# heterotic groups.
mut_dent <- dent_mat[, intersect(colnames(dent_mat), colnames(flint_mat))]
mut_flint <- flint_mat[, intersect(colnames(dent_mat), colnames(flint_mat))]


#################### DUPLICATES (1st OPTION) ################################
# Unique markers in Dent material for chromosome 1.
cfMarkDent <- callfreq.fun(x = mut_dent, output = "markerNames",
                           callThresh = 0.95)
mafMarkDent <- maf.fun(x = mut_dent, output = "markerNames", mafThresh = 0.05)
hetMarkDent <- heterozygosity.fun(x = mut_dent, output = "markerNames", 
                                  hetThresh = 0.05)
cfmafMarkDent <- intersect(intersect(cfMarkDent, mafMarkDent), hetMarkDent)
curedMarkDent <- dent_mat[, cfmafMarkDent]
dent_map <- marker_map[marker_map$markername %in% cfmafMarkDent, ]


# Unique markers in Flint material for chromosome 1.
cfMarkFlint <- callfreq.fun(x = mut_flint, output = "markerNames",
                            callThresh = 0.95)
mafMarkFlint <- maf.fun(x = mut_flint, output = "markerNames", 
                        mafThresh = 0.05)
hetMarkFlint <- heterozygosity.fun(x = mut_flint, output = "markerNames",
                                   hetThresh = 0.05)
cfmafMarkFlint <- intersect(cfMarkFlint, mafMarkFlint)
curedMarkFlint <- flint_mat[, cfmafMarkFlint]
flint_map <- marker_map[marker_map$markername %in% cfmafMarkFlint, ]


mList <- lapply(X = seq_len(10), FUN = function(i){
    markMatFlint <- curedMarkFlint[, flint_map[flint_map$chromosome == i,
                                              "markername"]]
    duplMarkFlint <- dupfreqm.fun(x = markMatFlint, output = "names")
    uniqueMarkFlint <- colnames(markMatFlint)[!colnames(markMatFlint) %in%
                                                  duplMarkFlint]
    uniqueMatFlint <- markMatFlint[, uniqueMarkFlint]
    if(any(duplicated(uniqueMatFlint, MARGIN = 2))){
        stop("Duplicated markers present")}
    
    markMatDent <- curedMarkDent[, dent_map[dent_map$chromosome == i,
                                           "markername"]]
    duplMarkDent <- dupfreqm.fun(x = markMatDent, output = "names")
    uniqueMarkDent <- colnames(markMatDent)[!colnames(markMatDent) %in% 
                                                duplMarkDent]
    uniqueMatDent <- markMatDent[, uniqueMarkDent]
    if(any(duplicated(uniqueMatDent, MARGIN = 2))){
        stop("Duplicated markers present")}
    
    # Combine the selected markers for the Dent and the Flint set, respectively.
    # Check the union of markers for duplicates in either pool, then keep all
    # markers that are still unique. Finally, keep the intersect of markers, 
    # which are unique in either heterotic groups.
    augMatDent <- markMatDent[, colnames(markMatDent) %in%
                                  union(uniqueMarkDent, uniqueMarkFlint)]
    updUnqDent <- setdiff(x = colnames(augMatDent), 
                          y = dupfreqm.fun(augMatDent, output = "names"))
    augMatFlint <- markMatFlint[, colnames(markMatFlint) %in% 
                                    union(uniqueMarkDent, uniqueMarkFlint)]
    updUnqFlint <- setdiff(x = colnames(augMatFlint),
                           y = dupfreqm.fun(augMatFlint, output = "names"))
    return(intersect(updUnqDent, updUnqFlint))    
})
sum(unlist(lapply(mList, length)))


# Plot of unique markers per chromosome.
centromeres = c(134400000, 92900000, 99700000, 105300000, 102300000, 49600000,
                54600000, 49000000, 72200000, 50100000)
pdf(file = paste0("./Derived/Plots/kernel_scanOption1Markers.pdf"), 
    height = 12, width = 14)
opar <- par(no.readonly = TRUE)
par(mfrow = c(5, 2), oma = c(0, 0, 3, 0))
for(i in seq_len(10)){
    plot(1, type = "n", yaxt = "n", 
         xlab = "Chromosome Position [Mbp]", ylab = "",
         xlim = c(0, max(marker_map[marker_map$chromosome == i, "position"])),
         ylim = c(0, 1), 
         main = paste("Chromosome", get("i"), "with",
                      length(mList[[i]]), "SNPs"))
    abline(v = marker_map[marker_map$markername %in% mList[[i]], "position"], 
           col = "#377eb8", lwd = 0.05)
    abline(v = centromeres[i], col = "#d95f02", lwd = 3)
}
# ... add the metabolite name in green color to the plot.
mtext("1st Option", cex = 1.5, outer = TRUE, col = "darkblue", font = 2)
par(opar)
dev.off()


# MARKER MATRIX PREPARATION -----------------------------------------------
# Concatenate the imputed marker matrices for dent and flint lines, 
# respectively.
imp_mat <- rbind(mut_dent, mut_flint)
ref_geno <- maf.fun(imp_mat[, unq_marker], output = "genoList")
rec_dent <- recode.fun(x = imp_mat[dent_parents, unq_marker],
                       major = ref_geno[["major_allele"]],
                       minor = ref_geno[["minor_allele"]],
                       major_coding = 2, minor_coding = 0,
                       het_coding = 1, na_coding = 999)
rec_dent <- matrix(as.numeric(rec_dent), 
                   nrow = nrow(rec_dent), ncol = ncol(rec_dent),
                   dimnames = dimnames(rec_dent))
if(any(rec_dent == 999)) stop("Missing values despite imputation")

rec_flint <- recode.fun(x = imp_mat[flint_parents, unq_marker],
                       major = ref_geno[["major_allele"]],
                       minor = ref_geno[["minor_allele"]],
                       major_coding = 2, minor_coding = 0,
                       het_coding = 1, na_coding = 999)
rec_flint <- matrix(as.numeric(rec_flint), 
                   nrow = nrow(rec_flint), ncol = ncol(rec_flint),
                   dimnames = dimnames(rec_flint))
if(any(rec_flint == 999)) stop("Missing values despite imputation")

# KINSHIP MATRIX STABILITY ------------------------------------------------
# Rationale: Selected markers shall not be confounded when building the kinship
# matrices. Therefore, the method of leaving out one chromosome
# (Rincent et al., 2014) is combined with picking distal loci at equi-distance.
# The following procedure aims at comparing the congruency of kinship matrices
# where all markers were picked at equi-distance but from different starting
# positions on the chromosome. If the correlation between the different 
# kinship matrices is high, any of them can be used for our QTL analyses since
# they reflect the same kinship among genotypes.
unq_marker <- unlist(mList)
unq_map <- droplevels(marker_map[marker_map$markername %in% unq_marker, ])
# Select equi-distant markers for each chromosome.
chr_list <- split(x = unq_map, f = list(unq_map$chromosome), drop = TRUE)

# Select equi-distant loci whose neighbors are ~ 1 Mbp apart and repeat this
# procedure 20 times.
equi_marker <- equi.dist.marker(markerMapList = chr_list, distance = 1e6,
                                rounds = 20)

dent_sim_kinship <- rounds.kinship(roundsMarkerDf = equi_marker, 
                                   genoMat = rec_dent)
flint_sim_kinship <- rounds.kinship(roundsMarkerDf = equi_marker, 
                                    genoMat = rec_flint)

library("reshape2")
dent_sim_list <- split(dent_sim_kinship, f = list(dent_sim_kinship$CHR))
dent_kin_cor <- vapply(X = dent_sim_list, FUN = function(x){
    tmpCast <- dcast(x, ID ~ Round, value.var = "Value")
    tmpCast$ID <- NULL
    corOut <- cor(tmpCast)[t(upper.tri(matrix(nrow = ncol(tmpCast),
                                              ncol = ncol(tmpCast))))]
    min(corOut)
}, FUN.VALUE = double(1))
if(any(dent_kin_cor < 0.95)) warning("Correlation among rounds < 0.95")

flint_sim_list <- split(flint_sim_kinship, f = list(flint_sim_kinship$CHR))
flint_kin_cor <- vapply(X = flint_sim_list, FUN = function(x){
    tmpCast <- dcast(x, ID ~ Round, value.var = "Value")
    tmpCast$ID <- NULL
    corOut <- cor(tmpCast)[t(upper.tri(matrix(nrow = ncol(tmpCast),
                                              ncol = ncol(tmpCast))))]
    min(corOut)
}, FUN.VALUE = double(1))
if(any(flint_kin_cor < 0.95)) warning("Correlation among rounds < 0.95")


# FINAL KINSHIP MATRIX DEFINITION -----------------------------------------
# Use marker from the first round of kinship marker selection.
equi_dent_kin <- droplevels(equi_marker[equi_marker$Round == 1, ])
equi_flint_kin <- droplevels(equi_marker[equi_marker$Round == 1, ])

# Generate lists of kinship matrices (one matrix for each chromosome) 
# according to Rincent et al. (2014).
dent_rincent <- rincent14.fun(markerMap = equi_dent_kin, genoMat = rec_dent)
names(dent_rincent) <- paste0("Chr.", seq_len(10))
flint_rincent <- rincent14.fun(markerMap = equi_flint_kin, genoMat = rec_flint)
names(flint_rincent) <- paste0("Chr.", seq_len(10))
K_chr <- c(dent_rincent, flint_rincent)
saveRDS(object = K_chr, file = "./Derived/SNP_Data/K_Chr_Rincent2014")


# Divide each chromosome into bins with a range of 5Mb ('5e6'), each. Bins can
# then be used for the computation of linkage disequilibria.
library("dplyr")
ord_unq_map <- unq_map[order(unq_map$chromosome, unq_map$position), ]
ord_unq_map$chromosome <- as.factor(as.character(ord_unq_map$chromosome))
ord_unq_map <- ord_unq_map %>%
    group_by(chromosome) %>%
    mutate(Bin = findInterval(position, seq(from = min(range(position)), 
                                            to = max(range(position)),
                                            by = 5e6)))
ord_unq_map <- as.data.frame(ord_unq_map)
ord_unq_map$Bin <- as.factor(as.character(ord_unq_map$Bin))
detach("package:dplyr")



library("reshape2")
geno_mat <- rbind(rec_dent, rec_flint)
geno_df <- data.frame(geno_mat)
geno_df$Pool <- NA
geno_df[rownames(geno_df) %in% rownames(rec_dent), "Pool"] <- "Dent"
geno_df[rownames(geno_df) %in% rownames(rec_flint), "Pool"] <- "Flint"
geno_df$Pool <- as.factor(geno_df$Pool)
geno_df$Genotype <- as.factor(rownames(geno_df))
geno_rec <- melt(data = geno_df,
                id.vars = c("Genotype", "Pool"),
                variable.name = "SNP",
                value.name = "allele") 

geno_rec <- merge(x = geno_rec,
                 y = ord_unq_map[, c("markername", "chromosome", 
                                   "position", "Bin")], 
                 by.x= "SNP", by.y= "markername", all.x= TRUE)
saveRDS(object= geno_rec, file= "./Derived/SNP_Data/mutMarkerWithPosition")
rm(list = ls())


