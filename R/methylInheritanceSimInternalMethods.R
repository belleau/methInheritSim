#' @title Estimate the alpha parameter of a Beta distribution
#'
#' @description Estimate the alpha parameter from the mean and the variance
#' of a Beta distribution.
#'
#' @param valCtrl a \code{vector} with 2 entries, the first value is the mean
#' and the second value is the variance of the controls (CTRL) at a specific 
#' CpG site.
#'
#' @param minVal a \code{double}, the minimum value accepted for the mean
#' value. If the first entry of \code{valCrtl} is smaller than 
#' \code{minVal}, then \code{minVal} is used in the calculation of the alpha
#' parameter. 
#' Default: \code{1e-06}.
#'
#' @return a \code{double}, the alpha parameter of a Beta distribution.
#'
#' @examples
#'
#' ## Estimate alpha parameters with mean = 0.5 and variance = 0.1
#' methylInheritanceSim:::estBetaAlpha(c(0.5,0.1))
#'
#' @author Pascal Belleau
#' @keywords internal
estBetaAlpha <- function(valCtrl, minVal = 1e-06){
    
    mu <- max(valCtrl[1], minVal)
    
    sigma2 <- max(valCtrl[2], ifelse(mu < 0.01, min(minVal, mu/1000), minVal))
    
    # mu must be smaller than 1 
    mu <- min(mu, 1 - min(0.001, sigma2 * 10^(-log10(minVal)/2)))
    
    return(max(0, -mu * (sigma2 + mu^2 - mu) / sigma2))
}


#' @title Estimate the beta paramater of a beta distribution
#'
#' @description Estimate the beta paramater from the mean and the variance
#' of a beta distribution.
#'
#' @param valCtrl  a \code{vector} with 2 entries, the first value is the mean
#' and the second value is the variance of the controls (CTRL) at a specific 
#' CpG site.
#'
#' @param minVal a \code{double}, the minimum value accepted for the mean
#' value. If the first entry of \code{valCrtl} is smaller than 
#' \code{minVal}, then \code{minVal} is used in the calculation of the beta
#' paramter. 
#' Default: \code{1e-06}.
#'
#' @return a \code{double}, the beta parameter of a Beta distribution.
#'
#' @examples
#'
#' ## Estimate beta parameters with mean = 0.5, variance = 0.1
#' methylInheritanceSim:::estBetaBeta(c(0.5,0.1))
#'
#' @author Pascal Belleau
#' @keywords internal
estBetaBeta <- function(valCtrl, minVal = 1e-06){
    
    mu <- max(valCtrl[1], minVal)
    
    # variance is at least minVal or mu / 1000
    sigma2 <- max(valCtrl[2], ifelse(mu < 0.01, min(minVal, mu/1000), minVal))
    
    # mu must be smaller than 1. 
    mu <- min(mu, 1 - min(0.001, sigma2 * 10^(-log10(minVal)/2)))
    
    return(max(0, (sigma2 + mu^2 - mu) * (mu -1) / sigma2))
}


#' @title Create a synthetic chromosome with the CTRL genome
#'
#' @description Create a synthetic chromosome with the sampling of a specified 
#' number of blocks and a specified number of consecutive CpG.
#'
#' @param methInfo is object of class \code{methylBase}, the CpG information
#' from controls (CTRL) that will be used to create the sythetic chromosome. 
#' The object can also contain information from cases but onl the controls will
#' be used.
#'
#' @param nbBlock \code{integer}, the number of blocks used for sampling.
#'
#' @param nbCpG a \code{integer}, the number of consecutive CpG positions used
#' for sampling from \code{methInfo}.
#'
#' @return a \code{GRanges} object, the synthetic chromosome
#'
#' @examples
#'
#' ## Load methyl information
#' data(samplesForChrSynthetic)
#' 
#' ## Ensure results are reproducible
#' set.seed(32)
#' 
#' ## Create synthetic chromosome
#' methylInheritanceSim:::getSyntheticChr(methInfo = samplesForChrSynthetic, 
#' nbBlock = 10, nbCpG = 20)
#'
#' @author Pascal Belleau
#' @importFrom stats runif var
#' @importFrom methylKit getData
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @keywords internal
getSyntheticChr <- function(methInfo, nbBlock, nbCpG) {
    
    seqChr <- unique(methInfo$chr) # list of Chr
    
    # Sample m chromosomes, the same can be sample more than once
    chr <- seqChr[sample(1:length(seqChr), nbBlock, replace=TRUE)]
    
    # The position of the CTRL
    posCTRL <- which(methInfo@treatment == 0)
    
    # Init the data.frame
    total <- nbBlock * nbCpG
    res <- data.frame(chr = rep("S", total), start=rep(0, total), 
                        end = rep(0, total),
                        meanCTRL = rep(0, total), 
                        varCTRL = rep(0, total),
                        alphaCTRL = rep(0, total), 
                        betaCTRL = rep(0, total),
                        chrOri = rep(0, total), 
                        startOri = rep(0, total))
    
    # First position
    l <- 1000
    
    for(i in 1:nbBlock){ # For each block of CpG
        
        # Select a position in the block
        v <- round(runif(1, 0, 1) *
                    (length(methInfo$start[methInfo$chr == chr[i]]) - nbCpG))
        
        methBlock <- getData(methInfo[methInfo$chr == 
                                            chr[i]][(v + 1):(v + nbCpG)])
        matProp <- sapply(posCTRL, function(x, methCur) {
                unname(unlist(methCur[3*(x - 1) + 2]/methCur[3*(x - 1) + 1]))},
                methCur = methBlock[5:length(methBlock)])
        
        res$chrOri[((i - 1) * nbCpG + 1):(i * nbCpG)] <- rep(chr[i], nbCpG)
        res$startOri[((i - 1) * nbCpG + 1):(i * nbCpG)] <- 
                unname(unlist(methBlock[2]))
        res$start[((i - 1) * nbCpG + 1):(i * nbCpG)] <- 
                unname(unlist(methBlock[2])) - 
                unname(unlist(methBlock[2]))[1] + l
        res$end[((i - 1) * nbCpG + 1):(i * nbCpG)] <- 
                res$start[((i - 1) * nbCpG + 1):(i * nbCpG)]
        res$meanCTRL[((i - 1) * nbCpG + 1):(i * nbCpG)] <- rowMeans(matProp)
        res$varCTRL[((i - 1) * nbCpG + 1):(i * nbCpG)] <- 
                apply(matProp, 1, var)
        l <- res$start[i * nbCpG] + 10000
    }
    
    res$alphaCTRL <- apply(res[, c(4,5)], 1, estBetaAlpha)
    res$betaCTRL  <- apply(res[, c(4,5)], 1, estBetaBeta)
    res <- GRanges(seqnames = res$chr, 
                    ranges = IRanges(start = res$start, end = res$end),
                    strand = rep("+", total), chrOri = res$chrOri,
                    startOri = res$startOri, meanCTRL = res$meanCTRL,
                    varCTRL = res$varCTRL)
    
    return(res)
}


#' @title get a proportion c/t for a case at a differentially 
#' 
#'
#' @description TODO
#' 
#'
#' @param x TODO
#'
#' @param nb TODO
#'
#' @param sDiff TODO
#'
#' @param diffCase TODO
#'
#' @param propDiffsd
#' -c * exp(b * log(distance between two CpG))
#'
#' @param b constant for distance
#'
#' @param endLength TODO
#'
#' @return TODO
#'
#' @examples
#'
#' ## TODO
#'
#'
#' @author Pascal Belleau
#' @keywords internal
getDiffCase <- function(x, nb, sDiff, diffCase, propDiffsd){
    
    meanDiff <- 0
    
    if(x[3] == 0){
        val <- rbeta(nb, estBetaAlpha(x[1:2]), estBetaBeta(x[1:2]))
        meanDiff <- x[1]
        partitionDiff <- c(0, nb)
    } else{
        meanDiff <- ifelse(x[1] < 0.5,
                            min(1, x[1]+sDiff),
                            max(0, x[1]-sDiff))
        
        partitionDiff <- c(diffCase, nb - diffCase)
        
        val <- c(rbeta(partitionDiff[1], estBetaAlpha(c(meanDiff, x[2])),
                    estBetaBeta(c(meanDiff, x[2]))),
                    rbeta(partitionDiff[2], estBetaAlpha(x), estBetaBeta(x)))
    }
    
    return(c(meanDiff, partitionDiff, val))
}


#' @title TODO
#'
#' @description TODO
#'
#' @param nbCtrl a non-negative \code{integer}, the number of controls
#'
#' @param nbCase a non-negative \code{integer}, the number of cases
#'
#' @param generation TODO
#'
#' @param stateInfo TODO
#'
#' @param stateDiff
#'
#' @param diffValue
#'
#' @param propDiff
#'
#' @param propDiffsd Default: \code{0.1}
#'
#' @param propInheritance
#'
#' @param propHetero
#'
#' @return TODO
#'
#' @examples
#'
#' ## TODO
#'
#'
#' @author Pascal Belleau
#' @importFrom msm rtnorm
#' @importFrom GenomicRanges GRangesList GRanges
#' @importFrom GenomeInfoDb seqnames
#' @importFrom IRanges ranges
#' @importFrom BiocGenerics strand
#' @importFrom S4Vectors mcols
#' @keywords internal
getSim <- function(nbCtrl, nbCase, generation, stateInfo, stateDiff, 
                    diffValue, propDiff, propDiffsd = 0.1, propInheritance, 
                    propHetero) {
    inR <- propDiff
    #res<-list()
    res <- GRangesList()
    if(propDiffsd < 0.0000001){
        diffCase <- round(nbCase * inR)
    } else{
        diffCase <- round(nbCase * rtnorm(1, mean = inR, sd = propDiffsd, 
                                            lower = 0, upper = 1))
    }
    
    ctrl <- t(apply(mcols(stateInfo)[3:4], 1, function(x, nb){
            rbeta(nb, estBetaAlpha(x), estBetaBeta(x))}, nb = nbCtrl))
    
    case <- t(apply(cbind(matrix(unlist(mcols(stateInfo)[3:4]) , ncol = 2), 
                            stateDiff$stateDiff), 1, getDiffCase, nb=nbCase, 
                            sDiff = diffValue, diffCase = diffCase, 
                            propDiffsd = propDiffsd))
    
    res[[1]] <- GRanges(seqnames = seqnames(stateInfo),
                        ranges = ranges(stateInfo),
                        strand =  strand(stateInfo),
                        meanDiff = case[, 1], meanCTRL = mcols(stateInfo)[3],
                        partitionCase = case[, 2], partitionCtrl = case[, 3],
                        ctrl = ctrl, case = case[, 4:length(case[1,])])

    for(i in 2:generation)
    {
        rm(case, ctrl)
        inR <- propDiff * propInheritance^(i - 2) # One generation move
        diffCur <- diffValue * propHetero       # Change diffValue
        if(propDiffsd < 0.0000001) {
            diffCase <- round(nbCase * inR)
        } else{
            diffCase <- round(nbCase * rtnorm(1, mean = inR, 
                                        sd = propDiffsd, lower = 0, upper = 1))
        }

        # Note mcols(stateInfo)[3:4] is a matrix with foreach position a row 
        # meanCTRL, varianceCTRL
        ctrl <- t(apply(mcols(stateInfo)[3:4], 1, function(x, nb) {
            rbeta(nb, estBetaAlpha(x), estBetaBeta(x))}, nb = nbCtrl))
        
        # matrix(unlist(mcols(stateInfo)[3:4]), nc = 2) is a matrix with 
        # foreach position a row with meanCTRL, varianceCTRL 
        case <- t(apply(cbind(matrix(unlist(mcols(stateInfo)[3:4]), ncol = 2),
                        stateDiff$stateInherite), 1,
                        getDiffCase, nb = nbCase, sDiff = diffCur, 
                        diffCase = diffCase, propDiffsd = propDiffsd))
        
        res[[i]] <- GRanges(seqnames = seqnames(stateInfo),
                        ranges = ranges(stateInfo),
                        strand =  strand(stateInfo),
                        meanDiff = case[, 1], meanCTRL = mcols(stateInfo)[3],
                        partitionCase = case[, 2], partitionCtrl = case[, 3],
                        ctrl = ctrl, case = case[,4:length(case[1,])])
    }
    
    return(res)
}

#' @title TODO
#'
#' @description Identify the pos where the case are Diff meth and which one are heritable
#'
#' @param stateInfo
#'
#' @param rateDiff
#'
#' @param minRate
#'
#' @param propInherite
#'
#' @param c
#'
#' @param b
#'
#' @param endLength
#'
#' @return TODO
#'
#' @examples
#'
#' ## TODO
#' 
#' @author Pascal Belleau
#' @importFrom stats rbeta rexp runif rpois
#' @keywords internal
getDiffMeth <- function(stateInfo, rateDiff, minRate, propInherite, 
                            c = 1.0, b = -1e-01, endLength=1000) {
    
    nbPos <- length(stateInfo)
    nbTry <- 1
    flag  <-  TRUE
    
    while(nbTry < 1000 & flag) {
        stateDiff <- rep(0, nbPos)
        stateInherite <- rep(0, nbPos)
        vExp <- rexp(nbPos, rateDiff)
        
        i<-1
        m<- max(1, round(vExp[i]))
        
        while(m <= nbPos){
            flagInherite <- ifelse(runif(1,0,1) < propInherite, TRUE,FALSE)
            stateDiff[m] <- 1
            if(flagInherite){
                stateInherite[m] <- 1
            }
            m <- m+1
            while(m <= nbPos && 
                  (start(stateInfo)[m] - start(stateInfo)[m-1]) <= endLength){
                cutOff <- c * 
                    exp(b*log(start(stateInfo)[m] - start(stateInfo)[m-1]))
                u <- runif(1,0,1)
                if(u < cutOff){
                    stateDiff[m] <- 1
                    if(flagInherite){
                        stateInherite[m] <- 1
                    }
                }
                m <- m + 1
            }
            i <- i+1
            m <- m + round(vExp[i])
        }
        length(which(stateDiff == 1))
        if(length(which(stateDiff == 1))>= minRate * nbPos){
            flag <- FALSE
        } else{
            warning(paste0("Nb: ", length(which(stateDiff == 1))))
        }
        nbTry <- nbTry + 1
    }
    if(flag){
        stateDiff <- NULL
        warning("Enable to generate the differentially methyyleted proportion fin\n")
    }
    
    return(list(stateDiff = stateDiff, stateInherite = stateInherite))
}


#' @title TODO
#'
#' @description TODO
#'
#' @param pathOut
#'
#' @param pref
#'
#' @param k
#'
#' @param nbCtrl
#'
#' @param nbCase
#'
#' @param treatment
#'
#' @param sample.id
#'
#' @param generation
#'
#' @param res
#'
#' @param rateDiff
#'
#' @param minRate
#'
#' @param propInherite
#'
#' @param diffValue
#'
#' @param propDiff
#'
#' @param propDiffsd
#'
#' @param propInheritance
#'
#' @param propHetero
#'
#' @param diffRes (default \code{NULL})
#'
#' @return TODO
#'
#' @examples
#'
#' ## TODO
#' 
#' @author Pascal Belleau
#' @importFrom methylKit read filterByCoverage normalizeCoverage unite calculateDiffMeth get.methylDiff getData tileMethylCounts methRead
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom stats rpois
#' @importFrom methods new
#' @importFrom stats start end
#' @keywords internal
simInheritance <- function(pathOut, pref, k, nbCtrl, nbCase, treatment, 
                    sample.id, generation, stateInfo,
                    rateDiff ,minRate, propInherite, diffValue, propDiff,
                    propDiffsd, propInheritance, propHetero, 
                    minReads = 10, maxPercReads = 99.9,
                    context = "CpG", assembly="Rnor_5.0",
                    meanCov = 80,
                    diffRes = NULL, saveGRanges = TRUE,
                    saveMethylKit = TRUE, anaMethylKit = TRUE){
    
    # Test if the simulation was done before
    # if just a part of the simulation is done it do it again
    
    alreadyDone <- TRUE
    if(! (file.exists(paste0(pathOut, "/stateDiff_", pref , "_", k,".rds"))) 
       || ! (file.exists(paste0(pathOut, "/simV0.1_", pref , "_", k,".rds")))){
        alreadyDone <- FALSE
    }
    if(saveGRanges && 
       ! (file.exists(paste0(pathOut, "/methylGR_", pref , "_", k,".rds")))) {
        alreadyDone <- FALSE
       }
    if(saveMethylKit &&
       ! (file.exists(paste0(pathOut, "/methylObj_", pref , "_", k,".rds")))) {
        alreadyDone <- FALSE
    }
    if(anaMethylKit &&
       ( ! (file.exists(paste0(pathOut, "/meth_", pref , "_", k,".rds")))
         || ! (file.exists(paste0(pathOut, "/methDiff_", pref , "_", k,".rds")))
       )) {
        alreadyDone <- FALSE
    }
    
    if(!(alreadyDone)){
        if(is.null(diffRes)){
            diffRes <- getDiffMeth(stateInfo=stateInfo,
                                   rateDiff=rateDiff, minRate=minRate,
                                   propInherite=propInherite)
        }
        simV0.1 <- getSim(nbCtrl=nbCtrl, nbCase=nbCase, generation=generation,
                          stateInfo=stateInfo, stateDiff=diffRes , diffValue=diffValue,
                          propDiff=propDiff, propDiffsd=propDiffsd,
                          propInheritance=propInheritance, propHetero = propHetero)
        
        saveRDS(diffRes,file = paste0(pathOut, "/stateDiff_", pref , "_", k,".rds"))
        saveRDS(simV0.1,file = paste0(pathOut, "/simV0.1_", pref , "_", k,".rds"))
        
        
        
        myobj <- list()
        myGR <- list()
        myMat <- list()
        myTr <- list()
        meth <- list()
        myDiff <- list()
        
        for(i in 1:generation){
            outList <- list()
            outGR <- GRangesList() 
            for(j in 1:(nbCtrl+nbCase)){
                coverage <- rpois( length(stateInfo), meanCov) + 1
                
                testM <- GRanges(seqnames = seqnames( stateInfo), 
                                    ranges = ranges(stateInfo),
                                    strand = strand(stateInfo),
                                    coverage=coverage,
                                    numCs = round(coverage * 
                                            unlist(mcols(simV0.1[[i]])[4+j])))
                if(saveMethylKit){
                    obj<-new("methylRaw", data.frame(chr = seqnames( testM), 
                                                    start = start(testM), 
                                                    end = end(testM),
                                                    strand = strand(testM),
                                                    coverage = testM$coverage, 
                                                    numCs = testM$numCs,
                                                    numTs = testM$coverage - testM$numCs)
                             ,sample.id=sample.id[[i]][[j]],assembly=assembly,
                             context=context,resolution='base')
                    
                    outList[[j]]<-obj
                }
                if(saveGRanges){
                    outGR[[j]]<-testM
                }
            }
            myMat[[i]] <- outList
            myTr[[i]] <- treatment
            
            if(saveMethylKit){
                myobj[[i]] <- new("methylRawList",outList,treatment=treatment)
            }
            
            if(saveGRanges){
                myGR[[i]] <- outList
            }
            
            if(anaMethylKit){
                filtered.myobj <- filterByCoverage(myobj[[i]],
                                                    lo.count = minReads, 
                                                    lo.perc=NULL, hi.count=NULL,
                                                    hi.perc= maxPercReads)
                filtered.myobj <- normalizeCoverage(filtered.myobj, "median")
                meth[[i]] <- unite(filtered.myobj, destrand=FALSE)
                myDiff[[i]] <- calculateDiffMeth(meth[[i]])
            }
        }
        
        if(saveGRanges){
            saveRDS(myGR,file = paste0(pathOut, "/methylGR_", pref , "_", k,".rds"))
        }
        
        if(saveMethylKit){
            saveRDS(myobj,file = paste0(pathOut, "/methylObj_", pref , "_", k,".rds"))
        }
        
        if(anaMethylKit){
            saveRDS(meth,file = paste0(pathOut, "/meth_", pref , "_", k,".rds"))
            saveRDS(myDiff,file = paste0(pathOut, "/methDiff_", pref , "_", k,".rds"))
        }
    }
    
}



