#!/usr/bin/Rscript --vanilla --slave
rm(list = ls())

# LD DECAY ----------------------------------------------------------------
cmbGeno <- readRDS("./Derived/SNP_Data/mutMarkerWithPosition")
genoList <- split(x = cmbGeno, f = list(cmbGeno$Pool, cmbGeno$chromosome),
                  drop = TRUE)
rm(cmbGeno)


# Generate sliding windows. 
# halfwin: size of the window (in base pairs) to each side of the central locus
# slide: step size of the sliding windows (in base pairs)
binList <- lapply(X = names(genoList), FUN = function(x, halfwin = 1.25e6,
                                                      slide = 2.5e5){
    library("reshape2")
    library("data.table")
    i <- genoList[[x]]
    Geno <- dcast(i, Genotype ~ SNP, value.var = "allele")
    Geno$Genotype <- NULL
    genoMat <- as.matrix(Geno)
    snpPos <- unique(i[, c("SNP", "position")])
    snpPos <- snpPos[order(snpPos$position), ]
    genoMat <- genoMat[, as.character(snpPos$SNP)]
    colnames(genoMat) <- snpPos$position
    # Ideal values (central position for each sliding window)
    sequ <- seq(from = halfwin, to = (max(snpPos[, 2]) - halfwin), by = slide)
    # Actual positions
    posi <- snpPos$position
    # Determine the SNP-marker position that is closest to the central 
    # position of the intended sliding window.
    dt <- data.table(posi, val = posi) 
    setattr(dt, "sorted", "posi")  
    dtval <- dt[J(sequ), roll = "nearest"]
    winCentr <- dtval$val
    Bins <- lapply(X = winCentr, FUN = function(y){
        binPos <- with(snpPos, snpPos[position > (y - halfwin) &
                                          position < (y + halfwin), 
                                      "position", drop = TRUE])
        binMat <- genoMat[, as.character(binPos)]
        binMat
        })
    names(Bins) <- winCentr
    # Remove bins/windows with only one locus.
    noMatrix <- names(Filter(function(x) class(x) != "matrix", x = Bins))
    Bins <- Bins[!names(Bins) %in% noMatrix]
    return(Bins)
    })
names(binList) <- names(genoList)
rm(genoList)

### Inverse of the relationship matrix using the formula by 
# Yang et al., 2010, "Common SNPs explain a large proportion of the 
# heritability for human height". Nature Genetics, Online Methods
### using the kinship correction suggested by
# Rincent et al., 2014, "Recovering Power in Association Mapping Panels
# with Variable Levels of Linkage Disequilibrium". Genetics, page 377    
K_chr <- readRDS("./Derived/SNP_Data/K_Chr_Rincent2014")


expand.grid.unique <- function(x, y, include.equals = FALSE){
    x <- unique(x)
    y <- unique(y)
    g <- function(i){
        z <- setdiff(y, x[seq_len(i - include.equals)])
        if(length(z)) cbind(x[i], z, deparse.level = 0)
    }    
    do.call(rbind, lapply(seq_along(x), g))
}

library("reshape2")
windList <- lapply(X = names(binList), FUN = function(i){
    V_inv <- K_chr[[i]][["KinMatrix"]]
    y <- binList[[i]]
    pool <- strsplit(i, split = "[.]")[[1]][1]
    chromosome <- as.numeric(strsplit(i, split = "[.]")[[1]][2])
    
    wins <- lapply(X = names(y), FUN = function(center){
        Geno <- y[[center]]
        position <- colnames(Geno)
        if(length(position) != length(unique(position))){
            position[duplicated(position)] <- as.numeric(
                position[duplicated(position)]) + 1
        }
        combs <- data.frame(expand.grid.unique(x = position, y = position))
        dmat <- as.matrix(dist(position))
        combs[, "SNPdist"] <- dmat[t(upper.tri(matrix(nrow = ncol(dmat), 
                                                      ncol = ncol(dmat))))]
        ### Kinship-correction 
        # Mangin et al., 2012, "Novel measures of linkage disequilibrium that 
        # correct the bias due to population structure and relatedness". 
        # Heredity,page 287
        fract <- (rep(1, nrow(Geno)) / rep(1, nrow(Geno)) %*% V_inv %*% 
                      rep(1, nrow(Geno))) %*% t(rep(1, nrow(Geno))) %*% V_inv
        loci <- as.matrix(Geno)
        block <- loci - fract %*% loci
        sig <- t(block) %*% V_inv %*% block
        rv <- cov2cor(sig)
        r2v <- cov2cor(sig)^2
        combs[, "R2V"] <- r2v[t(upper.tri(matrix(nrow = ncol(r2v),
                                                 ncol = ncol(r2v))))]
        combs[, "RV"] <- rv[t(upper.tri(matrix(nrow = ncol(rv),
                                               ncol = ncol(rv))))]
        ###
        combs <- combs[, !colnames(combs) %in%  c("X1", "X2")]
        combs$chromosome <- chromosome
        combs$Pool <- pool
        combs$Center <- center
        combs
    })
    windowDf <- do.call("rbind", wins)
    return(windowDf)
})
windowDf <- do.call("rbind", windList)
rm(list = setdiff(ls(), "windowDf"))
windowDf$Pool <- as.factor(windowDf$Pool)
windowDf$chromosome <- as.factor(as.character(windowDf$chromosome))
windowDf$Center <- as.numeric(windowDf$Center)
saveRDS(object = windowDf, file = "./Derived/SNP_Data/slidingWindowR2V")



# PLOTS -------------------------------------------------------------------
library("dplyr")
windowDf$Center <- as.factor(as.character(windowDf$Center))
windCount <- windowDf %>%
    group_by(Center) %>%
    summarize(count = n())
range(windCount$count)

library("ggplot2")
pdf(file = "./Derived/Plots/window_members.pdf", height = 12, width = 14)
(p1 <- ggplot(data = windCount, aes(x = count)) +
    geom_histogram(binwidth = 300, aes(fill = ..count..)) +
    scale_fill_gradient("Count", low = "green", high = "red") +
    theme_bw() +
    xlim(c(0, 5000)) +
    xlab("Pairwise Comparisons per Window") +
    ylab("Number of Windows"))
dev.off()
