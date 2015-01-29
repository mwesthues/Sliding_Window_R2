#!/usr/bin/Rscript --vanilla --slave
# REQUIRED PACKAGES -------------------------------------------------------
# install.packages("reshape2")


# REQUIRED OBJECTS --------------------------------------------------------
# mutMarkerWithPosition
# K_Chr_Rincent2014


# Select data from chromosome 1.
cmbGeno <- readRDS("./Derived/SNP_Data/mutMarkerWithPosition")
chr1Geno <- droplevels(cmbGeno[cmbGeno$chromosome == 1, ])
chr1List <- split(x = chr1Geno, f = list(chr1Geno$Pool), drop = TRUE)
names(chr1List) <- c("Dent.1", "Flint.1")
rm(list = c("cmbGeno", "chr1Geno"))

# Load the kinship matrices.
K_chr <- readRDS("./Derived/SNP_Data/K_Chr_Rincent2014")
K_chr <- K_chr[c("Dent.1", "Flint.1")]


expand.grid.unique <- function(x, y, include.equals = FALSE){
    x <- unique(x)
    y <- unique(y)
    g <- function(i){
        z <- setdiff(y, x[seq_len(i - include.equals)])
        if(length(z)) cbind(x[i], z, deparse.level = 0)
    }    
    do.call(rbind, lapply(seq_along(x), g))
}

# Calculate kinship-corrected R^2 (R2V)
res <- lapply(X = names(chr1List), FUN = function(x){
    library("reshape2")
    V_inv <- K_chr[[x]][["KinMatrix"]]
    i <- chr1List[[x]]
    pool <- strsplit(x, split = "[.]")[[1]][1]
    chromosome <- as.numeric(strsplit(x, split = "[.]")[[1]][2])
    castDf <- dcast(i, Genotype ~ SNP, value.var = "allele")
    castDf$Genotype <- NULL
    snpPos <- unique(i[, c("SNP", "position")])
    position <- snpPos[, "position"]
    snpNames <- as.vector(snpPos[, "SNP"])
    names(position) <- snpNames
    castDf <- castDf[, snpNames]
    combs <- data.frame(expand.grid.unique(x = snpNames, y = snpNames))
    dmat <- as.matrix(dist(position))
    combs[, "SNPdist"] <- dmat[t(upper.tri(matrix(nrow = ncol(dmat), 
                                                  ncol = ncol(dmat))))]
    ### Kinship-correction 
    # Mangin et al., 2012, "Novel measures of linkage disequilibrium that 
    # correct the bias due to population structure and relatedness". 
    # Heredity,page 287
    fract <- (rep(1, nrow(castDf)) / rep(1, nrow(castDf)) %*% V_inv %*% 
                  rep(1, nrow(castDf))) %*% t(rep(1, nrow(castDf))) %*% V_inv
    loci <- as.matrix(castDf)
    block <- loci - fract %*% loci
    sig <- t(block) %*% V_inv %*% block
    rv <- cov2cor(sig)
    r2v <- cov2cor(sig)^2
    combs[, "R2V"] <- r2v[t(upper.tri(matrix(nrow = ncol(r2v),
                                             ncol = ncol(r2v))))]
    combs[, "RV"] <- rv[t(upper.tri(matrix(nrow = ncol(rv),
                                           ncol = ncol(rv))))]
    ###
    combs <- combs[, !colnames(combs) %in%  c("X1", "X2", "RV")]
    combs$Pool <- pool
    return(combs)
})
genoDf <- do.call("rbind", res)
genoDf$Pool <- as.factor(genoDf$Pool)
rm(list = c("res", "chr1List"))


dent <- droplevels(genoDf[genoDf$Pool == "Dent", c("SNPdist", "R2V")])
saveRDS(object = dent, file = "/home/westhues/Dropbox/chromosome1/dent_chr1")
flint <- droplevels(genoDf[genoDf$Pool == "Flint", c("SNPdist", "R2V")])
saveRDS(object = flint, file = "/home/westhues/Dropbox/chromosome1/flint_chr1")
