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
#' The object can also contain information from cases but only the controls 
#' will be used.
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
#' @param nbCtrl a positive \code{integer}, the number of controls.
#'
#' @param nbCase a positive \code{integer}, the number of cases.
#'
#' @param generation a positive \code{integer}, the number of generations.
#'
#' @param stateInfo a \code{GRanges} object, the synthetic chromosome 
#' generated by \code{getSyntheticChr} function. 
#' TODO (ajouter les champs de metadata ?)
#'
#' @param stateDiff a \code{list} with 2 entries. The first entry is called 
#' \code{stateDiff} and contains a \code{vector} of \code{integer} (\code{0} 
#' and \code{1}) with 
#' a length corresponding the length of \code{stateInfo}. The \code{statDiff}
#' indicates, using a \code{1}, the positions where the CpG sites are
#' differentially methylated. The second entry is
#' called \code{statInherite} and contains a \code{vector} of \code{integer} 
#' (\code{0} and \code{1})
#' with a length corresponding the length of \code{stateInfo}. The 
#' \code{statInherite}
#' indicates, using a \code{1}, the positions where the CpG values are
#' inherited.
#'
#' @param diffValue a non-negative \code{double} between between [0,1], the 
#' proportion of C/T for a case differentially methylated following a 
#' beta distribution 
#' where the mean is shifted of \code{diffValue} from the CTRL distribution.
#'
#' @param propDiff a \code{double} superior to \code{0} and inferior or equal 
#' to \code{1}, the mean value for the proportion of samples that will have,
#' for a specific position, differentially methylated values. It can be 
#' interpreted as the penetrance.
#'
#' @param propDiffsd a non-negative \code{double}, the standard deviation 
#' associated to the \code{propDiff}.
#'
#' @param propInheritance a non-negative \code{double} between [0,1], the 
#' proportion of case that inherite differentially methylated sites.
#'
#' @param propHetero TODO a non-negative \code{double} between [0,1], the 
#' reduction of \code{diffValue} for the second and following generations.
#'
#' @return TODO
#'
#' @examples
#'
#' ## Fix seed to have reproducible results
#' set.seed(312)
#' 
#' ## Load dataset
#' data("samplesForChrSynthetic")
#' 
#' ## Generate a stateInfo object using samples
#' stateInfo <- methylInheritanceSim:::getSyntheticChr(methInfo = 
#' samplesForChrSynthetic, nbBlock = 1, nbCpG = 3)
#' 
#' ## Generate a stateDiff object with length corresponding to
#' ## nbBlock * nbCpG from stateInfo
#' stateDiff <- list()
#' stateDiff[["stateDiff"]] <- c(1, 0, 1)
#' stateDiff[["stateInherite"]] <- c(1, 0, 0)
#' 
#' ## Create a simulation using stateInfo and stateDiff
#' methylInheritanceSim:::getSim(nbCtrl = 3, nbCase = 2, generation = 3, 
#' stateInfo = stateInfo, stateDiff = stateDiff, diffValue = 10, 
#' propDiff = 0.8, propDiffsd = 0.2, propInheritance = 0.8, propHetero = 0.1)
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
#' @param c a \code{double}, TODO Default: \code{1.0}.
#'
#' @param b a \code{double}, TODO . Default: \code{-1e-01}.
#'
#' @param endLength TODO . Default: \code{1000}.
#'
#' @return a \code{list} containing the following elements:
#' \itemize{
#' \item \code{stateDiff} a \code{vector} TODO. The length of the 
#' \code{vector} corresponds to the length of 
#' the \code{stateInfo} parameter.
#' \item \code{stateDiff} a \code{vector} TODO. The length of the 
#' \code{vector} corresponds to the length of 
#' the \code{stateInfo} parameter.
#' }
#'
#' @examples
#'
#' ## TODO
#' 
#' @author Pascal Belleau
#' @importFrom BiocGenerics start
#' @importFrom stats rbeta rexp runif rpois
#' @keywords internal
getDiffMeth <- function(stateInfo, rateDiff, minRate, propInherite, 
                            c = 1.0, b = -1e-01, endLength = 1000) {
    
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
            flagInherite <- ifelse(runif(1, 0, 1) < propInherite, TRUE, FALSE)
            stateDiff[m] <- 1
            if(flagInherite){
                stateInherite[m] <- 1
            }
            m <- m+1
            while(m <= nbPos && 
                (start(stateInfo)[m] - start(stateInfo)[m - 1]) <= endLength){
                print(paste0("Aye ", m))
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
        
        
        if (length(which(stateDiff == 1)) >= minRate * nbPos) {
            flag <- FALSE
        }
        
        nbTry <- nbTry + 1
    }
    
    if (flag) {
        stateDiff <- NULL
        stop("Enable to generate the differentially methyyleted proportion fin\n")
    }
    
    return(list(stateDiff = stateDiff, stateInherite = stateInherite))
}


#' @title TODO
#'
#' @description TODO
#'
#' @param pathOut TODO
#'
#' @param pref TODO
#'
#' @param k TODO
#'
#' @param nbCtrl a positive \code{integer}, the number of controls.
#'
#' @param nbCase a positive \code{integer}, the number of cases.
#'
#' @param treatment TODO
#'
#' @param sample.id TODO
#'
#' @param generation a positive \code{integer}, the number of generations.
#'
#' @param stateInfo TODO
#'
#' @param rateDiff TODO
#'
#' @param minRate TODO
#'
#' @param propInherite TODO
#'
#' @param diffValue TODO
#'
#' @param propDiff TODO
#'
#' @param propDiffsd TODO
#'
#' @param propInheritance TODO
#'
#' @param propHetero TODO
#'
#' @param minReads TODO
#' 
#' @param maxPercReads TODO
#' 
#' @param context TODO. Default: \code{"CpG"}.
#' 
#' @param assembly TODO. Default: \code{"Rnor_5.0"}.
#' 
#' @param meanCov TODO. Default: \code{80}.
#' 
#' @param diffRes TODO. Default: \code{NULL}.
#'
#' @param saveGRanges TODO. Default: \code{TRUE}.
#' 
#' @param saveMethylKit TODO. Default: \code{TRUE}.
#' 
#' @param anaMethylKit TODO. Default: \code{TRUE}.
#' 
#' @return \code{0} indicating that the function has been successful.
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
#' @importFrom BiocGenerics start end strand
#' @keywords internal
simInheritance <- function(pathOut, pref, k, nbCtrl, nbCase, treatment, 
                    sample.id, generation, stateInfo,
                    rateDiff ,minRate, propInherite, diffValue, propDiff,
                    propDiffsd, propInheritance, propHetero, 
                    minReads = 10, maxPercReads = 99.9,
                    context = "CpG", assembly="Rnor_5.0",
                    meanCov = 80, diffRes = NULL, saveGRanges = TRUE,
                    saveMethylKit = TRUE, anaMethylKit = TRUE) {
    
    # Test if the simulation was done before
    # if just a part of the simulation is done it do it again
    
    alreadyDone <- TRUE
    if(! (file.exists(paste0(pathOut, "/stateDiff_", pref, "_", k, ".rds"))) 
       || ! (file.exists(paste0(pathOut, 
                                    "/simV0.1_", pref, "_", k, ".rds")))) {
        alreadyDone <- FALSE
    }
    if(saveGRanges && 
       ! (file.exists(paste0(pathOut, "/methylGR_", pref, "_", k,".rds")))) {
        alreadyDone <- FALSE
       }
    if(saveMethylKit &&
       ! (file.exists(paste0(pathOut, 
                                "/methylObj_", pref, "_", k, ".rds")))) {
        alreadyDone <- FALSE
    }
    if(anaMethylKit &&
       ( ! (file.exists(paste0(pathOut, "/meth_", pref, "_", k,".rds")))
         || ! (file.exists(paste0(pathOut, 
                                    "/methDiff_", pref, "_", k, ".rds"))))) {
        alreadyDone <- FALSE
    }
    
    if (!(alreadyDone)) {
        if (is.null(diffRes)) {
            diffRes <- getDiffMeth(stateInfo = stateInfo,
                                    rateDiff = rateDiff, minRate = minRate,
                                    propInherite = propInherite)
        }
        
        simV0.1 <- getSim(nbCtrl = nbCtrl, nbCase = nbCase, 
                            generation = generation, stateInfo = stateInfo, 
                            stateDiff = diffRes, diffValue = diffValue,
                            propDiff = propDiff, propDiffsd = propDiffsd,
                            propInheritance = propInheritance, 
                            propHetero = propHetero)
        
        saveRDS(diffRes,file = paste0(pathOut, "/stateDiff_", pref , 
                                        "_", k,".rds"))
        saveRDS(simV0.1,file = paste0(pathOut, "/simV0.1_", pref , 
                                        "_", k,".rds"))
        
        myobj <- list()
        myGR <- list()
        myMat <- list()
        #myTr <- list()
        meth <- list()
        myDiff <- list()
        
        for(i in 1:generation) {
            outList <- list()
            outGR <- GRangesList() 
            for(j in 1:(nbCtrl+nbCase)){
                coverage <- rpois( length(stateInfo), meanCov) + 1
                
                testM <- GRanges(seqnames = seqnames( stateInfo), 
                                    ranges = ranges(stateInfo),
                                    strand = strand(stateInfo),
                                    coverage = coverage,
                                    numCs = round(coverage * 
                                            unlist(mcols(simV0.1[[i]])[4+j])))
                if(saveMethylKit){
                    obj<-new("methylRaw", data.frame(chr = seqnames( testM), 
                                        start = start(testM), 
                                        end = end(testM),
                                        strand = strand(testM),
                                        coverage = testM$coverage, 
                                        numCs = testM$numCs,
                                        numTs = testM$coverage - testM$numCs),
                            sample.id = sample.id[[i]][[j]], 
                            assembly = assembly,
                            context = context, resolution = 'base')
                    
                    outList[[j]]<-obj
                }
                if(saveGRanges){
                    outGR[[j]]<-testM
                }
            }
            myMat[[i]] <- outList
            #myTr[[i]] <- treatment
            
            if(saveMethylKit){
                myobj[[i]] <- new("methylRawList", outList,
                                    treatment = treatment)
            }
            
            if(saveGRanges){
                myGR[[i]] <- outList
            }
            
            if(anaMethylKit){
                filtered.myobj <- filterByCoverage(myobj[[i]],
                                        lo.count = minReads, 
                                        lo.perc = NULL, hi.count = NULL,
                                        hi.perc = maxPercReads)
                filtered.myobj <- normalizeCoverage(filtered.myobj, "median")
                meth[[i]] <- unite(filtered.myobj, destrand = FALSE)
                myDiff[[i]] <- calculateDiffMeth(meth[[i]])
            }
        }
        
        if(saveGRanges){
            saveRDS(myGR,file = paste0(pathOut, "/methylGR_", pref, 
                                            "_", k, ".rds"))
        }
        
        if(saveMethylKit){
            saveRDS(myobj,file = paste0(pathOut, "/methylObj_", pref, 
                                            "_", k, ".rds"))
        }
        
        if(anaMethylKit){
            saveRDS(meth,file = paste0(pathOut, "/meth_", pref, 
                                            "_", k,".rds"))
            saveRDS(myDiff,file = paste0(pathOut, "/methDiff_", pref, 
                                            "_", k, ".rds"))
        }
    }
    
    return(0)
}


#' @title Parameters validation for the \code{\link{runSim}} function
#'
#' @description Validation of all parameters needed by the public
#' \code{\link{runSim}} function.
#'
#' @param outputDir a string of \code{character}, the path where the 
#' files created by the function will be saved or \code{NULL}. 
#'
#' @param fileGen a string of \code{character}, TODO 
#' include each output file. Each output 
#' file are 
#' composed with a type name (methylGR, methylObj, ...), _, fileGen (ex F1),
#' parameters of the simulation and ".rds". 
#'
#' @param nbSynCHR a positive \code{integer}, the number of distinct synthetic 
#' chromosomes. generate. TODO
#'
#' @param methData an object of class \code{methylBase}, the CpG information
#' from controls (CTRL) that will be used to create the sythetic chromosome. 
#' The \code{methData} object can also contain information from cases but 
#' only the controls will be used.
#'
#' @param nbBlock \code{integer}, the number of blocks used for sampling.
#'
#' @param lBlock a \code{integer}, the number of consecutive CpG positions used
#' for sampling from \code{methInfo}.
#'
#' @param vNbSample a \code{vector} of positive \code{integer}, the number of 
#' methData (CTRL) and cases in the the simulation dataset. In 
#' the simulated dataset, the number of CTRL equals the number of Case. 
#' The number of CTRL do not need to be equal to the number of Case in
#' the real dataset.
#'
#' @param nbGeneration a positive \code{integer}, the number of generations.
#'
#' @param vpDiff a \code{double} superior to \code{0} and inferior or equal 
#' to \code{1}, the mean value for the proportion of samples that will have,
#' for a specific position, differentially methylated values. It can be 
#' interpreted as the penetrance.
#' 
#' @param vpDiffsd a non-negative \code{double}, the standard deviation 
#' associated to the \code{propDiff}.
#'
#' @param vDiff a positive \code{double} between [0,1], the proportion of 
#' C/T for a case differentially methylated follow a beta distribution 
#' where the mean is shifted of \code{vDiff} from the CTRL distribution
#'
#' @param vInheritance a positive \code{double} between [0,1], the 
#' proportion of cases that inherited differentially sites.
#' 
#' @param propInherite a non-negative \code{double} inferior or equal to \code{1}, 
#' proportion of differentially methylated site
#' are inherated
#'
#' @param rateDiff a positive \code{double} inferior to \code{1}, the mean of 
#' the chance that a site is differentially 
#' methylated.
#'
#' @param minRate a non-negative \code{double} inferior to \code{1}, the minimum 
#' rate of differentially methylated sites.
#'
#' @param propHetero a positive \code{double} between [0,1], the 
#' reduction of vDiff for the second and following generations.
#' 
#' @param minReads a positive \code{integer} Bases and regions having lower
#' coverage than this count are discarded. The parameter
#' correspond to the \code{lo.count} parameter in the \code{methylKit} package.
#' 
#' @param maxPercReads a \code{double} between [0,100], the percentile of read
#' counts that is going to be used as upper cutoff. Bases ore regions
#' having higher
#' coverage than this percentile are discarded. Parameter used for both CpG
#' sites and tiles analysis. The parameter
#' correspond to the \code{hi.perc} parameter in the  \code{methylKit} package.
#' 
#' @param context a string of \code{character}, the methylation context 
#' string, ex: CpG,CpH,CHH, etc. 
#' Default: \code{"CpG"}.
#' 
#' @param assembly a string of \code{character}, the short description of the 
#' genome assembly. Ex: mm9,hg18 etc.
#' Default: \code{"Rnor_5.0"}
#' 
#' @param meanCov a positive \code{integer} represent the mean of the coverage
#' at the CpG site Default: \code{80}.
#' 
#' @param n a positive \code{integer}, the number of simulation for each 
#' parameters (\code{vNbSample}, \code{vpDiff}, \code{vDiff} and
#' \code{vInheritance}).
#' 
#' @param keepDiff \code{logical} if true, the differentially methyled sites
#' will be the same for each parameter (\code{vpDiff}, 
#' \code{vDiff} and \code{vInheritance}). Default: \code{FALSE}.
#' 
#' @param saveGRanges \code{logical} if true, save a list of \code{nbGeneration}
#' GRangesList. Each GRangeList contain \code{vNbSample} GRanges which contain
#' for each CpG site: the position, the coverage and the proportion of the C/T.
#' A files treatment_... which contains the position of the CTRL and Case are 
#' save too. Default: \code{TRUE}.
#' 
#' @param saveMethylKit a \code{logical}, TODO. Default: \code{TRUE}.
#' 
#' @param anaMethylKit a \code{logical}, TODO. Default: \code{TRUE}
#' 
#' @param nbCores a positive \code{integer}, the number of cores to use when
#' processing the analysis. Default: \code{1}.
#' 
#' @param vSeed a \code{integer}, a seed used when reproducible results are
#' needed. When a value inferior or equal to zero is given, a random integer
#' is used. Default: \code{-1}.
#' 
#' @return \code{0} indicating that the function has been successful.
#'
#' @examples
#'
#' ## TODO
#' 
#' @author Pascal Belleau, Astrid Deschenes
#' @importFrom S4Vectors isSingleInteger isSingleNumber
#' @keywords internal
validateRunSimParameters <-function(outputDir, fileGen, nbSynCHR, methData, 
                                    nbBlock, lBlock,
                                    vNbSample, nbGeneration, vpDiff, 
                                    vpDiffsd, vDiff, 
                                    vInheritance,
                                    propInherite, rateDiff, minRate, 
                                    propHetero, 
                                    minReads, 
                                    maxPercReads, context, assembly,
                                    meanCov, n, keepDiff,
                                    saveGRanges, saveMethylKit,
                                    anaMethylKit,
                                    nbCores, vSeed) {
    
    ## Validate that the outputDir is an not empty string
    if (!is.null(outputDir) && !is.character(outputDir)) {
        stop("outputDir must be a character string or NULL")
    }
    
    ## Validate that the fileGen is an not empty string
    if (!is.null(fileGen) && !is.character(fileGen)) {
        stop("fileGen must be a character string or NULL")
    }
    
    ## Validate that nbSynCHR is an positive integer
    if (!(isSingleInteger(nbSynCHR) || isSingleNumber(nbSynCHR)) ||
        as.integer(nbSynCHR) < 1) {
        stop("nbSynCHR must be a positive integer or numeric")
    }
    
    ## Validate that methData is methylBase class from methylKit
    if(!(class(methData) == "methylBase") ){
        stop("methylBase must be methylBase class from methylKit")
    }
    
    ## Validate that nbBlock is an positive integer
    if (!(isSingleInteger(nbBlock) || isSingleNumber(nbBlock)) ||
        as.integer(nbBlock) < 1) {
        stop("nbBlock must be a positive integer or numeric")
    }
    
    ## Validate that lBlock is an positive integer
    if (!(isSingleInteger(lBlock) || isSingleNumber(lBlock)) ||
        as.integer(lBlock) < 1) {
        stop("lBlock must be a positive integer or numeric")
    }
    
    ## Validate that vNbSample is an positive integer
    if (!(isSingleInteger(vNbSample) || isSingleNumber(vNbSample)) ||
        as.integer(vNbSample) < 1) {
        stop("vNbSample must be a positive integer or numeric")
    }
    
    ## Validate that nbGeneration is an positive integer
    if (!(isSingleInteger(nbGeneration) || isSingleNumber(nbGeneration)) ||
        as.integer(nbGeneration) < 1) {
        stop("nbGeneration must be a positive integer or numeric")
    }
    
    ## Validate that vpDiff is an positive double between (0,1]
    if (!(isSingleNumber(vpDiff)) ||
        vpDiff <= 0.00 || vpDiff > 1.00) {
        stop("vpDiff must be a positive double between (0,1]")
    }
    
    ## Validate that vpDiffsd is an non-negative double 
    if (!(isSingleNumber(vpDiffsd)) ||
        vpDiff < 0.00) {
        stop("vpDiffsd is an non-negative double")
    }
    
    ## Validate that vDiff is an positive double between [0,1]
    if (!(isSingleNumber(vDiff)) ||
        vDiff < 0.00 || vDiff > 1.00) {
        stop("vDiff must be a positive double between [0,1]")
    }
    
    ## Validate that vInheritance is an positive double between [0,1]
    if (!(isSingleNumber(vInheritance)) ||
        vInheritance < 0.00 || vInheritance > 1.00) {
        stop("vInheritance must be a positive double between [0,1]")
    }
    
    ## Validate that propInherite is an positive double between [0,1]
    if (!(isSingleNumber(propInherite)) ||
        propInherite < 0.00 || propInherite > 1.00) {
        stop("propInherite must be a positive double between [0,1]")
    }
    
    ## Validate that rateDiff is an positive double between (0,1)
    if (!(isSingleNumber(rateDiff)) ||
        rateDiff <= 0.00 || rateDiff >= 1.00) {
        stop("rateDiff must be a positive double between (0,1)")
    }
    
    ## Validate that minRate is an positive double between [0,1)
    if (!(isSingleNumber(minRate)) ||
        minRate < 0.00 || minRate >= 1.00) {
        stop("minRate must be a positive double between [0,1)")
    }
    
    ## Validate that propHetero is an positive double between [0,1]
    if (!(isSingleNumber(propHetero)) ||
        propHetero < 0.00 || propHetero > 1.00) {
        stop("propHetero must be a positive double between [0,1]")
    }
    
    ## Validate that minReads is an positive integer
    if (!(isSingleInteger(minReads) || isSingleNumber(minReads)) ||
        as.integer(minReads) < 1) {
        stop("minReads must be a positive integer or numeric")
    }
    
    ## Validate that maxPercReads is an positive double between [0,100]
    if (!(isSingleNumber(maxPercReads)) ||
        maxPercReads < 0.00 || maxPercReads > 100.00) {
        stop("maxPercReads must be a positive double between [0,100]")
    }
    
#    context, assembly,
#    meanCov, n,
    ## Validate that keepDiff is a logical
    if (!is.logical(keepDiff)) {
        stop("keepDiff must be a logical")
    }
    
#    saveGRanges, saveMethylKit,
#    anaMethylKit,
#    nbCores, vSeed
}
