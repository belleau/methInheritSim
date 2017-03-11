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


#' @title Estimate the beta parameter of a beta distribution
#'
#' @description Estimate the beta parameter from the mean and the variance
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
estBetaBeta <- function(valCtrl, minVal = 1e-06) {
    
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
#' @param nbBlock a positive \code{integer}, the number of blocks used 
#' for sampling.
#'
#' @param nbCpG a \code{integer}, the number of consecutive CpG positions used
#' for sampling from \code{methInfo}.
#'
#' @return a \code{GRanges} object, the synthetic chromosome.
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
    
    # Sample nbBlock chromosomes, the same can be sample more than once
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
    
    ## Create returned value
    res <- GRanges(seqnames = res$chr, 
                    ranges = IRanges(start = res$start, end = res$end),
                    strand = rep("+", total), chrOri = res$chrOri,
                    startOri = res$startOri, meanCTRL = res$meanCTRL,
                    varCTRL = res$varCTRL)
    
    return(res)
}


#' @title Get a proportion C/T for each case for the sites selected as
#' differentially methylated or not
#' 
#'
#' @description Simulate the proportion of C/T for each case at the sites 
#' selected as differentially methylated or not.
#' 
#'
#' @param x a \code{vector} of \code{double} containing 3 entries: 
#' \itemize{
#' \item the mean of the CTRL at this sites
#' \item the variance of the CTRL at this sites
#' \item \code{1} if the site is selected as differentially methylated, 
#' otherwise \code{0}
#' }
#'
#' @param nb an \code{integer}, the number of cases.
#'
#' @param sDiff a non-negative \code{double} 
#' included in [0,1], the proportion of C/T for a case differentially 
#' methylated that follows 
#' a beta distribution where the mean is shifted of \code{vDiff} 
#' from the CTRL distribution.
#'
#' @param diffCase an \code{integer}, the number of cases differentially at 
#' the selected as differentially methylated site.
#'
#' @return a \code{vector} containing 3 + nb entries:
#' \itemize{
#' \item mean of proportion of C/T of the differentially methylated case
#' \item The number of case simulate with shifted distribution
#' \item The number of case simulate with the control distribution
#' \item the proportion of C/T for each case
#' }
#'
#' @examples
#'
#' ## Create vector containing 3 + nb entries:
#' ## 1 - The mean of the CTRL at this sites
#' ## 2 - The variance of the CTRL at this sites
#' ## 3 - 1 when DMS; otherwise 0
#' x <- c(0.9814562, 0.0003607153, 0)
#' 
#' ## Get the proportion of C/T for each case at a specific site.
#' methylInheritanceSim:::getDiffCase(x=x, nb=6, sDiff = 0.8, 
#' diffCase = round(6 * 0.9))
#' 
#' @author Pascal Belleau, Astrid Deschenes
#' @importFrom stats rbeta
#' @keywords internal
getDiffCase <- function(x, nb, sDiff, diffCase) {
    
    meanDiff <- 0
    
    if(x[3] == 0) {
        val <- rbeta(nb, estBetaAlpha(x[1:2]), estBetaBeta(x[1:2]))
        meanDiff <- x[1]
        partitionDiff <- c(0, nb)
    } else {
        meanDiff <- ifelse(x[1] < 0.5,
                            min(1, x[1] + sDiff),
                            max(0, x[1] - sDiff))
        
        partitionDiff <- c(diffCase, nb - diffCase)
        
        val <- c(rbeta(partitionDiff[1], estBetaAlpha(c(meanDiff, x[2])),
                    estBetaBeta(c(meanDiff, x[2]))),
                    rbeta(partitionDiff[2], estBetaAlpha(x), estBetaBeta(x)))
    }
    
    return(c(meanDiff, partitionDiff, val))
}


#' @title Simulate the proportion of C/T at each site of synthetic CHR for 
#' each control and case
#'
#' @description For each control and case, generate the proportion of C/T at 
#' each of the synthetic CHR.
#'
#' @param nbCtrl a positive \code{integer}, the number of controls.
#'
#' @param nbCase a positive \code{integer}, the number of cases.
#'
#' @param generation a positive \code{integer}, the number of generations.
#'
#' @param stateInfo a \code{GRanges} object, the synthetic chromosome 
#' generated by \code{getSyntheticChr} function. 
#'
#' @param stateDiff a \code{list} with 2 entries:
#' \itemize{
#' \item \code{stateDiff} a \code{vector} of \code{integer} (\code{0} 
#' and \code{1}) with length corresponding the length of \code{stateInfo}.
#' The \code{vector}
#' indicates, using a \code{1}, the positions where the CpG sites are
#' differentially methylated.
#' \item \code{stateInherite} a \code{vector} of \code{integer} (\code{0} and 
#' \code{1})
#' with length corresponding the length of \code{stateInfo}. The 
#' \code{vector}
#' indicates, using a \code{1}, the positions where the CpG values are
#' inherited.
#' }
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
#' @param propHetero a non-negative \code{double} between [0,1], the 
#' reduction of \code{vDiff} for the second and following generations.
#'
#' @return a \code{GRangesList} object contains information about the 
#' simulation. The file have four metadata related to real dataset:
#' \itemize{
#' \item meanDiff, the means of the shifted distribution
#' \item meanCTRL, the means of the control distribution
#' \item partitionCase, the number of cases simulated with the shifted 
#' distribution
#' \item partitionCtrl, the number of cases simulated with the control 
#' distribution and a metadata for each cases and controls 
#' the proportion of C/T.
#' }
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
#' stateInformation <- methylInheritanceSim:::getSyntheticChr(methInfo = 
#'     samplesForChrSynthetic, nbBlock = 1, nbCpG = 3)
#' 
#' ## Generate a stateDiff object with length corresponding to
#' ## nbBlock * nbCpG from stateInformation
#' stateDiff <- list()
#' stateDiff[["stateDiff"]] <- c(1, 0, 1)
#' stateDiff[["stateInherite"]] <- c(1, 0, 0)
#' 
#' ## Create a simulation using stateInfo and stateDiff
#' methylInheritanceSim:::getSim(nbCtrl = 3, nbCase = 2, generation = 3, 
#'     stateInfo = stateInformation, stateDiff = stateDiff, diffValue = 10, 
#'     propDiff = 0.8, propDiffsd = 0.2, propInheritance = 0.8, 
#'     propHetero = 0.1)
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
    
    if (propDiffsd < 0.0000001) {
        diffCase <- round(nbCase * inR)
    } else{
        diffCase <- round(nbCase * rtnorm(1, mean = inR, sd = propDiffsd, 
                                            lower = 0, upper = 1))
    }
    
    ctrl <- t(apply(mcols(stateInfo)[3:4], 1, function(x, nb){
            rbeta(nb, estBetaAlpha(x), estBetaBeta(x))}, nb = nbCtrl))
    
    case <- t(apply(cbind(matrix(unlist(mcols(stateInfo)[3:4]) , ncol = 2), 
                            stateDiff$stateDiff), 1, getDiffCase, nb=nbCase, 
                            sDiff = diffValue, diffCase = diffCase))
    
    # TODO change meanCTRL.meanCTRL in meanCTRL
    #tmpCol <- matrix(mcols(stateInfo)[3]$meanCTRL, nc = 1)
    res[[1]] <- GRanges(seqnames = seqnames(stateInfo),
                        ranges = ranges(stateInfo),
                        strand =  strand(stateInfo),
                        meanDiff = case[, 1], 
                        meanCTRL = mcols(stateInfo)[3],
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
                        diffCase = diffCase))
        
        res[[i]] <- GRanges(seqnames = seqnames(stateInfo),
                        ranges = ranges(stateInfo),
                        strand =  strand(stateInfo),
                        meanDiff = case[, 1], 
                        meanCTRL = mcols(stateInfo)[3],
                        partitionCase = case[, 2], partitionCtrl = case[, 3],
                        ctrl = ctrl, case = case[,4:length(case[1,])])
    }
    
    return(res)
}

#' @title Identify differentially methylated sites and among those, the ones
#' that are inherited.
#'
#' @description Identify the site positions where the cases are differentially 
#' methylated and, among those, the one that are inherited.
#'
#' @param stateInfo a \code{GRanges} that contains the CpG (or methylated 
#' sites).
#' The \code{GRanges} have four metadata from the real dataset:
#' \itemize{
#' \item chrOri, the chromosome from the real dataset
#' \item startOri, the position of the site in the real dataset
#' \item meanCTRL, the mean of the control in the real dataset
#' \item varCTRL, the variance of the control in the real dataset
#' }
#'
#' @param rateDiff a positive \code{double} inferior to \code{1}, the mean of 
#' the chance that a site is differentially 
#' methylated.
#'
#' @param minRate a non-negative \code{double} inferior to \code{1}, the 
#' minimum rate for differentially methylated sites.
#'
#' @param propInherite a non-negative \code{double} inferior or equal 
#' to \code{1}, 
#' the proportion of differentially methylated regions that 
#' are inherated.
#'
#' @param c a positive \code{double}, a factor in the formula to compute the 
#' probabylity of site to be diffentially methylated in a differentially
#' methylated region. 
#' The probability formula of site in differentially methylated region is
#' \code{c} * exp(\code{b} * log(distance with the preceding sites))
#' Default: \code{1.0}.
#'
#' @param b a negative \code{double}, a factor in the formula to compute the 
#' probabylity of site to be diffentially methylated in a differentially
#' methylated region. 
#' The probability formula of site in differentially methylated region is
#' \code{c} * exp(\code{b} * log(distance with the preceding sites)). 
#' Default: \code{-1e-01}.
#'
#' @param endLength a positive \code{integer}, when the distance with the 
#' preceding sites in a differentially
#' methylated region is larger than \code{endLength}, the differentially
#' methylated region is finished. Default: \code{1000}.
#'
#' @return a \code{list} containing the 2 following elements:
#' \itemize{
#' \item \code{stateDiff} a \code{vector} of \code{integer} (\code{0} 
#' and \code{1}) with length corresponding the length of \code{stateInfo}.
#' The \code{vector}
#' indicates, using \code{1}, the positions where the CpG sites are
#' differentially methylated.
#' \item \code{stateInherite} a \code{vector} of \code{integer} (\code{0} and 
#' \code{1})
#' with length corresponding the length of \code{stateInfo}. The 
#' \code{vector}
#' indicates, using \code{1}, the positions where the CpG values are
#' inherited.
#' }
#'
#' @examples
#' 
#' ## Load dataset containing a list of objects used by 
#' ## methylInheritanceSim internal functions
#' data(dataSimExample)
#' 
#' ## Identify differentially methylated sites and among those, the ones
#' ## that are inherited
#' methylInheritanceSim:::getDiffMeth(stateInfo = 
#'     dataSimExample$stateInfo, rateDiff = 0.3, minRate = 0.3,
#'     propInherite = 0.3)
#' 
#' @author Pascal Belleau, Astrid Deschenes
#' @importFrom BiocGenerics start
#' @importFrom stats rbeta rexp runif rpois
#' @keywords internal
getDiffMeth <- function(stateInfo, rateDiff, minRate, propInherite, 
                            c = 1.0, b = -1e-01, endLength = 1000) {
    
    nbPos <- length(stateInfo)
    nbTry <- 1
    flag  <-  TRUE
    
    while (nbTry < 1000 & flag) {
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
            m <- m + 1
            while(m <= nbPos && 
                (start(stateInfo)[m] - start(stateInfo)[m - 1]) <= endLength) {
                cutOff <- c * 
                    exp(b*log(start(stateInfo)[m] - start(stateInfo)[m-1]))
                
                u <- runif(1,0,1)
                if(u < cutOff) {
                    stateDiff[m] <- 1
                    if(flagInherite) {
                        stateInherite[m] <- 1
                    }
                }
                m <- m + 1
            }
            i <- i + 1
            m <- m + round(vExp[i])
        }
        
        
        if (length(which(stateDiff == 1)) >= minRate * nbPos) {
            flag <- FALSE
        }
        
        nbTry <- nbTry + 1
    }
    
    if (flag) {
        stateDiff <- NULL
        stop(paste0("Enable to generate the differentially methyyleted ", 
                        "proportion fin\n"))
    }
    
    return(list(stateDiff = stateDiff, stateInherite = stateInherite))
}


#' @title Simulate a multigeneration methylation experiment with inheritance
#'
#' @description Simulate a multigeneration methylation case versus control 
#' experiment 
#' with inheritance relation using a real control dataset. 
#' 
#' The simulation can 
#' be parametrized to fit different models. The number of cases and controls, 
#' the proportion of the case affected 
#' by the treatment (penetrance), the effect of the treatment on the mean of 
#' the distribution, the proportion of sites inherited, the proportion of the 
#' differentially methylated sites from the precedent generation inherited, 
#' etc..
#' 
#' The function simulates a multigeneration dataset like a bisulfite 
#' sequencing experiment. The simulation includes the information about 
#' control and case for each generation.
#'
#' @param pathOut a string of \code{character} or \code{NULL}, the path 
#' where the 
#' files created by the function will be saved. When \code{NULL}, the files
#' are saved in the current directory.
#'
#' @param pref a string of \code{character} representing the parameters of
#' specific simulation
#' the string is 
#' composed of those elements, separated by "_":
#' \itemize{ 
#' \item a \code{fileID}
#' \item the chromosome number, a number between 1 and \code{nbSynCHR}
#' \item the number of samples, a number in the \code{vNbSample} \code{vector}
#' \item the mean proportion of samples that has,
#' for a specific position, differentially methylated values, a 
#' number in the \code{vpDiff} \code{vector}
#' \item the proportion of 
#' C/T for a case differentially methylated that follows a shifted beta 
#' distribution, a
#' number in the \code{vDiff} \code{vector}
#' \item the 
#' proportion of cases that inherits differentially sites, a number in the
#' \code{vInheritance} \code{vector}
#' }
#'
#' @param k a positive \code{integer}, a Id for the current simulation.
#'
#' @param nbCtrl a positive \code{integer}, the number of controls.
#'
#' @param nbCase a positive \code{integer}, the number of cases.
#'
#' @param treatment a numeric vector denoting controls and cases
#'
#' @param sample.id a matrix the name of each samples for each generation (row)
#' and each case and control (column).
#'
#' @param generation a positive \code{integer}, the number of generations
#' simulated.
#'
#' @param stateInfo a GRanges that contains the CpG (or methylated sites).
#' The GRamges have four metadata from the real dataset:
#' chrOri the chromosome from the real dataset
#' startOri the position of the site in the real dataset
#' meanCTRL the mean of the control in the real dataset
#' varCTRL the variance of the control in the real dataset.
#'
#' @param rateDiff a positive \code{double} inferior to \code{1}, the mean of 
#' the chance that a site is differentially methylated.
#'
#' @param minRate a non-negative \code{double} inferior to \code{1}, the 
#' minimum rate for differentially methylated sites.
#' Default: \code{0.01}.
#'
#' @param propInherite a non-negative \code{double} inferior or equal 
#' to \code{1}, 
#' the proportion of differentially methylated regions that 
#' are inherated.
#'
#' @param diffValue a non-negative \code{double} 
#' included in [0,1], the proportion of C/T for a case differentially 
#' methylated that follows 
#' a beta distribution where the mean is shifted by \code{vDiff} 
#' from the CTRL distribution.
#'
#' @param propDiff a \code{double} superior to 
#' \code{0} and inferior or equal 
#' to \code{1}, the mean value for the proportion of samples that will have,
#' for a specific position, differentially methylated values. It can be 
#' interpreted as the penetrance.
#'
#' @param propDiffsd a non-negative \code{double}, the 
#' standard deviation associated to the \code{vpDiff}. Note that 
#' \code{vpDiff} and \code{vpDiffsd} must be the same length.
#'
#' @param propInheritance a non-negative \code{double} 
#' included in [0,1], the proportion of cases 
#' that inherits differentially methylated sites.
#'
#' @param propHetero a non-negative \code{double} between [0,1], the 
#' reduction of \code{vDiff} for the second and following generations.
#'
#' @param minReads a positive \code{integer}, sites and regions having lower
#' coverage than this count are discarded. The parameter
#' corresponds to the \code{lo.count} parameter in 
#' the \code{methylKit} package.
#' 
#' @param maxPercReads a \code{double} between [0,100], the percentile of read
#' counts that is going to be used as upper cutoff. Sites and regions
#' having higher
#' coverage than \code{maxPercReads} are discarded. This parameter is used for 
#' both CpG sites and tiles analysis. The parameter
#' correspond to the \code{hi.perc} parameter in the \code{methylKit} package.
#' 
#' @param context a string of \code{character}, the short description of the 
#' methylation context, such as "CpG", "CpH", "CHH", etc.. Default: "CpG"
#' 
#' @param assembly a string of \code{character}, the short description of the 
#' genome assembly, such as "mm9", "hg18", etc..
#' 
#' @param meanCov a positive \code{integer}, the mean of the coverage
#' at the simulated CpG sites.
#' 
#' @param diffRes a \code{list} with 2 entries:
#' \itemize{
#' \item \code{stateDiff} a \code{vector} of \code{integer} (\code{0} 
#' and \code{1}) with length corresponding the length of \code{stateInfo}.
#' The \code{vector}
#' indicates, using a \code{1}, the positions where the CpG sites are
#' differentially methylated.
#' \item \code{stateInherite} a \code{vector} of \code{integer} (\code{0} and 
#' \code{1})
#' with length corresponding the length of \code{stateInfo}. The 
#' \code{vector}
#' indicates, using a \code{1}, the positions where the CpG values are
#' inherited.
#' } when is \code{NULL} generate a new ones with \code{getDiffMeth}.
#'
#' @param saveGRanges a \code{logical}, when \code{true}, the package save two 
#' files type. The first generate for each simulation contains a \code{list}. 
#' The length of the \code{list} corresponds to the number of generation. 
#' The generation are stored in order (first entry = first generation, 
#' second entry = second generation, etc..). All samples related to one 
#' generations are contained in a \code{GRangesList}. 
#' The \code{GRangeaList} store a \code{list} of \code{GRanges}. Each 
#' \code{GRanges} stores the raw mehylation data of one sample.
#' The second file a numeric \code{vector} denoting controls and cases 
#' (a file is generates by entry in the \code{vector} parameters 
#' \code{vNbSample}).
#' 
#' @param saveMethylKit a \code{logical}, when \code{TRUE}, the package save 
#' a file contains a \code{list}. The length of the 
#' \code{list} corresponds to the number of generation. The generation are 
#' stored in order (first entry = first generation, 
#' second entry = second generation, etc..). All samples related to one 
#' generations are contained in a S4 \code{methylRawList} object. The 
#' \code{methylRawList} object contains two Slots:
#' 1. treatment: A numeric \code{vector} denoting controls and cases.
#' 2. .Data: A \code{list} of \code{methylRaw} objects. Each object stores the 
#' raw methylation data of one sample.
#' 
#' @param runAnalysis a \code{logical}, if \code{TRUE}, two files are saved :
#' \itemize{
#' \item 1. The first file is the methylObj... file formated 
#' with the \code{methylkit} package in a S4 \code{methylBase} 
#' object (with the \code{methylKit} 
#' functions: \code{filterByCoverage}, \code{normalizeCoverage} and 
#' \code{unite}).
#' \item 2. The second file contains a S4 \code{calculateDiffMeth} object 
#' generated with the \code{methylKit} functions \code{calculateDiffMeth} 
#' using the first file.
#' }
#' 
#' @return \code{0} indicating that the function has been successful.
#'
#' @examples
#'
#' temp_dir <- "test_simInheritance"
#' data(dataSimExample)
#' 
#' \dontrun{methylInheritanceSim:::simInheritance(pathOut = temp_dir,
#' pref = paste0("S1_", "6_0.9_0.8_0.5"),
#' k = 1, nbCtrl = 6, nbCase = 6, 
#' treatment = dataSimExample$treatment, 
#' sample.id = dataSimExample$sample.id,
#' generation = 3, 
#' stateInfo = dataSimExample$stateInfo,
#' propDiff = 0.9, propDiffsd = 0.1,
#' diffValue = 0.8,
#' propInheritance = 0.5,
#' rateDiff = 0.3, minRate = 0.3,
#' propInherite = 0.3, 
#' propHetero = 0.5,
#' saveGRanges = FALSE,
#' saveMethylKit = FALSE,
#' runAnalysis = FALSE
#' )}
#' 
#' ## Delete temp_dir
#' \dontrun{if (dir.exists(temp_dir)) {
#' unlink(temp_dir, recursive = TRUE, force = FALSE)
#' }}
#' 
#' @author Pascal Belleau, Astrid Deschenes
#' @importFrom methylKit read filterByCoverage normalizeCoverage unite 
#' calculateDiffMeth get.methylDiff getData tileMethylCounts methRead
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom stats rpois
#' @importFrom methods new
#' @importFrom BiocGenerics start end strand
#' @keywords internal
simInheritance <- function(pathOut, pref, k, nbCtrl, nbCase, treatment, 
                        sample.id, generation, stateInfo,
                        propDiff, propDiffsd, diffValue, propInheritance, 
                        rateDiff , minRate, propInherite, 
                        propHetero, minReads, maxPercReads,
                        context = "CpG", assembly, meanCov, diffRes, 
                        saveGRanges,
                        saveMethylKit, runAnalysis) {
    
    # Test if the simulation was done before
    # if just a part of the simulation is done it do it again
    if (!is.null(pathOut) && !dir.exists(pathOut)) {
            dir.create(pathOut, showWarnings = TRUE)
    }
    
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
    if(runAnalysis &&
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
                if (saveGRanges) {
                    outGR[[j]]<-testM
                }
            }
            myMat[[i]] <- outList
            
            if (saveMethylKit) {
                myobj[[i]] <- new("methylRawList", outList,
                                    treatment = treatment)
            }
            
            if (saveGRanges) {
                myGR[[i]] <- outList
            }
            
            if (runAnalysis) {
                filtered.myobj <- filterByCoverage(myobj[[i]],
                                        lo.count = minReads, 
                                        lo.perc = NULL, hi.count = NULL,
                                        hi.perc = maxPercReads)
                filtered.myobj <- normalizeCoverage(filtered.myobj, "median")
                meth[[i]]   <- unite(filtered.myobj, destrand = FALSE)
                myDiff[[i]] <- calculateDiffMeth(meth[[i]])
            }
        }
        
        if (saveGRanges) {
            saveRDS(myGR, file = paste0(pathOut, "/methylGR_", pref, 
                                            "_", k, ".rds"))
        }
        
        if (saveMethylKit) {
            saveRDS(myobj, file = paste0(pathOut, "/methylObj_", pref, 
                                            "_", k, ".rds"))
        }
        
        if (runAnalysis) {
            saveRDS(meth, file = paste0(pathOut, "/meth_", pref, 
                                            "_", k,".rds"))
            saveRDS(myDiff, file = paste0(pathOut, "/methDiff_", pref, 
                                            "_", k, ".rds"))
        }
    }
    
    return(0)
}


#' @title Parameters validation for the \code{\link{runSim}} function. Only
#' integer parameters are validated.
#'
#' @description Validation of all parameters needed by the public
#' \code{\link{runSim}} function. Only integer parameters are validated.
#'
#' @param nbSynCHR a positive \code{integer}, the number of distinct 
#' synthetic chromosomes that will be generated.
#' 
#' @param nbSimulation a positive \code{integer}, the number of simulations 
#' for each parameter (\code{vNbSample}, \code{vpDiff}, \code{vDiff} and
#' \code{vInheritance}).
#'
#' @param nbBlock a positive \code{integer}, the number of blocks used 
#' for sampling.
#'
#' @param nbCpG a positive \code{integer}, the number of consecutive CpG 
#' positions used for sampling from \code{methInfo}.
#'
#' @param vNbSample a \code{vector} of positive \code{integer}, the number of 
#' methData (CTRL) and cases in the the simulation dataset. In 
#' the simulated dataset, the number of CTRL equals the number of Case. 
#' The number of CTRL do not need to be equal to the number of Case in
#' the real dataset.
#'
#' @param nbGeneration a positive \code{integer}, the number of generations.
#' 
#' @param minReads a positive \code{integer} Bases and regions having lower
#' coverage than this count are discarded. The parameter
#' correspond to the \code{lo.count} parameter in the \code{methylKit} package.
#' 
#' @param meanCov a positive \code{integer} represent the mean of the coverage
#' at the CpG site Default: \code{80}.
#' 
#' @param nbCores a positive \code{integer}, the number of cores to use when
#' creating the simulated datasets. Default: \code{1} and always 
#' \code{1} for Windows.
#' 
#' @param vSeed a \code{integer}, a seed used when reproducible results are
#' needed. When a value inferior or equal to zero is given, a random integer
#' is used. Default: \code{-1}.
#' 
#' @return \code{0} indicating that the function has been successful.
#'
#' @examples
#' 
#' ## The function returns 0 when all paramaters are valid
#' methylInheritanceSim:::validateRunSimIntegerParameters(nbSynCHR = 1, 
#' nbSimulation = 2, nbBlock = 10, nbCpG = 4, vNbSample = 10, 
#' nbGeneration = 3, minReads = 10, meanCov = 80, 
#' nbCores = 1, vSeed = -1)
#' 
#' @author Pascal Belleau, Astrid Deschenes
#' @importFrom S4Vectors isSingleInteger isSingleNumber
#' @keywords internal
validateRunSimIntegerParameters <-function(nbSynCHR, nbSimulation, nbBlock, 
                                    nbCpG, vNbSample, nbGeneration, minReads, 
                                    meanCov, nbCores, vSeed) {

    ## Validate that nbSynCHR is an positive integer
    if (!(isSingleInteger(nbSynCHR) || isSingleNumber(nbSynCHR)) ||
        as.integer(nbSynCHR) < 1) {
        stop("nbSynCHR must be a positive integer or numeric")
    }
    
    ## Validate that nbSimulation is an positive integer
    if (!(isSingleInteger(nbSimulation) || isSingleNumber(nbSimulation)) ||
        as.integer(nbSimulation) < 1) {
        stop("nbSimulation must be a positive integer or numeric")
    }
    
    ## Validate that nbBlock is an positive integer
    if (!(isSingleInteger(nbBlock) || isSingleNumber(nbBlock)) ||
        as.integer(nbBlock) < 1) {
        stop("nbBlock must be a positive integer or numeric")
    }
    
    ## Validate that nbCpG is an positive integer
    if (!(isSingleInteger(nbCpG) || isSingleNumber(nbCpG)) ||
        as.integer(nbCpG) < 1) {
        stop("nbCpG must be a positive integer or numeric")
    }
    
    ## Validate that nbGeneration is an positive integer
    if (!(isSingleInteger(nbGeneration) || isSingleNumber(nbGeneration)) ||
        as.integer(nbGeneration) < 1) {
        stop("nbGeneration must be a positive integer or numeric")
    }
    
    ## Validate that vNbSample is a vector of distinct positive integer
    if (! is.numeric(vNbSample) || anyDuplicated(vNbSample) > 0 ||
        any(vNbSample < 1) || ! all(as.integer(vNbSample) == vNbSample)) {
        stop("vNbSample must be a vector of distinct positive integer")
    }
    
    ## Validate that minReads is an positive integer
    if (!(isSingleInteger(minReads) || isSingleNumber(minReads)) ||
        as.integer(minReads) < 1) {
        stop("minReads must be a positive integer or numeric")
    }
    
    ## Validate that meanCov is an positive integer
    if (!(isSingleInteger(meanCov) || isSingleNumber(meanCov)) ||
        as.integer(meanCov) < 1) {
        stop("meanCov must be a positive integer or numeric")
    }
    
    ## Validate that nbCores is an positive integer
    if (!(isSingleInteger(nbCores) || isSingleNumber(nbCores)) ||
        as.integer(nbCores) < 1) {
        stop("nbCores must be a positive integer or numeric")
    }
    
    ## Validate that nbCores is set to 1 on Windows system
    if (Sys.info()["sysname"] == "Windows" && as.integer(nbCores) != 1) {
        stop("nbCores must be 1 on a Windows system")
    }
    
    ## Validate that vSeed is an integer
    if (!(isSingleInteger(vSeed) || isSingleNumber(vSeed))) {
        stop("vSeed must be an integer or numeric")
    }
    
    return(0)
}


#' @title Parameters validation for the \code{\link{runSim}} function. Only
#' double parameters are validated.
#'
#' @description Validation of all parameters needed by the public
#' \code{\link{runSim}} function. Only double parameters are validated.
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
#' @param propInherite a non-negative \code{double} inferior or equal to 
#' \code{1}, the proportion of differentially methylated site
#' are inherated
#'
#' @param rateDiff a positive \code{double} inferior to \code{1}, the mean of 
#' the chance that a site is differentially methylated.
#'
#' @param minRate a non-negative \code{double} inferior to \code{1}, the 
#' minimum rate of differentially methylated sites.
#'
#' @param propHetero a positive \code{double} between [0,1], the 
#' reduction of vDiff for the second and following generations.
#' 
#' @param maxPercReads a \code{double} between [0,100], the percentile of read
#' counts that is going to be used as upper cutoff. Bases ore regions
#' having higher
#' coverage than this percentile are discarded. Parameter used for both CpG
#' sites and tiles analysis. The parameter
#' correspond to the \code{hi.perc} parameter in the  \code{methylKit} package.
#' 
#' @return \code{0} indicating that the function has been successful.
#'
#' @examples
#'
#' ## The function returns 0 when all paramaters are valid
#' methylInheritanceSim:::validateRunSimDoubleParameters(vpDiff =0.2, 
#' vpDiffsd = 0.3, vDiff = 0.4, vInheritance = 0.2, propInherite = 0.5, 
#' rateDiff = 0.2, minRate = 0.1, propHetero = 0.2, maxPercReads = 99.1)
#' 
#' @author Pascal Belleau, Astrid Deschenes
#' @importFrom S4Vectors isSingleInteger isSingleNumber
#' @keywords internal
validateRunSimDoubleParameters <-function(vpDiff, vpDiffsd, vDiff, 
                                        vInheritance,propInherite, rateDiff, 
                                        minRate, propHetero, minReads, 
                                        maxPercReads) {
    
    ## Validate that vpDiff is an positive double include in (0,1]
    if (! is.numeric(vpDiff) || anyDuplicated(vpDiff) > 0 ||
        any(vpDiff <= 0.00) || any(vpDiff > 1.00)) {
        stop(paste0("vpDiff must be a vector of distinct positive double ", 
                    "include in (0,1]"))
    }
    
    ## Validate that vpDiffsd is an non-negative double
    if (! is.numeric(vpDiffsd) || any(vpDiffsd < 0.00) ) {
        stop("vpDiffsd must be a vector of non-negative double")
    }
    
    ## Validate that vpDiff and vpDiffsd must be the same length
    if (length(vpDiff) != length(vpDiffsd)) {
        stop("vpDiff and vpDiffsd must be the same length")
    }
    
    ## Validate that vDiff is a vector of distinct non-negative double 
    ## include in  [0,1]
    if (! is.numeric(vDiff) || anyDuplicated(vDiff) > 0 ||
        any(vDiff < 0.00) || any(vDiff > 1.00)) {
        stop(paste0("vDiff must be a vector of distinct non-negative double", 
                    " include in [0,1]"))
    }
    
    ## Validate that vInheritance is a vector of distinct non-negative double 
    ## include in [0,1]
    if (! is.numeric(vInheritance) || anyDuplicated(vInheritance) > 0 ||
        any(vInheritance < 0.00) || any(vInheritance > 1.00)) {
        stop(paste0("vInheritance must be a vector of distinct non-negative ", 
                    "double include in [0,1]"))
    }
    
    ## TODO finish unit test
    ## Validate that rateDiff is an positive double include in (0,1)
    if (!(isSingleNumber(rateDiff)) || rateDiff <= 0.00 || rateDiff >= 1.00) {
        stop("rateDiff must be a positive double include in (0,1)")
    }
    
    ## Validate that minRate is an non-negative double include in [0,1)
    if (!(isSingleNumber(minRate)) || minRate < 0.00 || minRate >= 1.00) {
        stop("minRate must be a non-negative double include in [0,1)")
    }
    
    ## Validate that propInherite is a non-negative double include in [0,1]
    if (!(isSingleNumber(propInherite)) ||
        propInherite < 0.00 || propInherite > 1.00) {
        stop("propInherite must be a non-negative double include in [0,1]")
    }
    
    
    ## Validate that propHetero is an non-negative double include in [0,1]
    if (!(isSingleNumber(propHetero)) ||
        propHetero < 0.00 || propHetero > 1.00) {
        stop("propHetero must be a non-negative double include in [0,1]")
    }
    
    ## Validate that maxPercReads is an positive double between [0,100]
    if (!(isSingleNumber(maxPercReads)) ||
        maxPercReads < 0.00 || maxPercReads > 100.00) {
        stop("maxPercReads must be a positive double between [0,100]")
    }
    
    return(0)
}


#' @title Parameters validation for the \code{\link{runSim}} function. Only 
#' logical parameters are validated.
#'
#' @description Validation of all parameters needed by the public
#' \code{\link{runSim}} function. Only logical parameters are validated.
#' 
#' @param keepDiff \code{logical} if true, the differentially methyled sites
#' will be the same for each parameter (\code{vpDiff}, 
#' \code{vDiff} and \code{vInheritance}). Default: \code{FALSE}.
#' 
#' @param saveGRanges a \code{logical}, when \code{true}, the package save two 
#' files type. The first generate for each simulation contains a \code{list}. 
#' The length of the \code{list} corresponds to the number of generation. 
#' The generation are stored in order (first entry = first generation, 
#' second entry = second generation, etc..). All samples related to one 
#' generations are contained in a \code{GRangesList}. 
#' The \code{GRangeaList} store a \code{list} of \code{GRanges}. Each 
#' \code{GRanges} stores the raw mehylation data of one sample.
#' The second file a numeric \code{vector} denoting controls and cases 
#' (a file is generates by entry in the \code{vector} parameters 
#' \code{vNbSample}).
#' 
#' @param saveMethylKit a \code{logical}, when \code{TRUE}, for each 
#' simulations save a file contains a \code{list}. The length of the 
#' \code{list} corresponds to the number of generation. The generation are 
#' stored in order (first entry = first generation, 
#' second entry = second generation, etc..). All samples related to one 
#' generations are contained in a S4 \code{methylRawList} object. The 
#' \code{methylRawList} object contains two Slots:
#' 1. treatment: A numeric \code{vector} denoting controls and cases.
#' 2. .Data: A \code{list} of \code{methylRaw} objects. Each object stores the 
#' raw methylation data of one sample.
#' 
#' @param runAnalysis a \code{logical}, if \code{TRUE}, two files are saved 
#' for each simulation:
#' \itemize{
#' \item 1. The first file is the methylObj... file formated with 
#' the \code{methylkit} package in a S4 \code{methylBase} object 
#' (with the \code{methylKit} functions: \code{filterByCoverage}, 
#' \code{normalizeCoverage} and \code{unite}).
#' \item 2. The second file contains a S4 \code{calculateDiffMeth} object 
#' generated with the \code{methylKit} functions \code{calculateDiffMeth} on 
#' the first file.
#' }
#' 
#' @return \code{0} indicating that the function has been successful.
#'
#' @examples
#'
#' ## Load dataset
#' data("samplesForChrSynthetic")
#'
#' ## The function returns 0 when all paramaters are valid
#' methylInheritanceSim:::validateRunSimLogicalParameters(keepDiff = FALSE, 
#' saveGRanges = TRUE, saveMethylKit = FALSE, runAnalysis = FALSE)
#' 
#' @author Pascal Belleau, Astrid Deschenes
#' @importFrom S4Vectors isSingleInteger isSingleNumber
#' @keywords internal
validateRunSimLogicalParameters <-function(keepDiff, saveGRanges, 
                                            saveMethylKit, runAnalysis) {
    
    ## Validate that keepDiff is a logical
    if (!is.logical(keepDiff)) {
        stop("keepDiff must be a logical")
    }
    
    ## Validate that saveGRanges is a logical
    if (!is.logical(saveGRanges)) {
        stop("saveGRanges must be a logical")
    }
    
    ## Validate that saveMethylKit is a logical
    if (!is.logical(saveMethylKit)) {
        stop("saveMethylKit must be a logical")
    }
    
    ## Validate that runAnalysis is a logical
    if (!is.logical(runAnalysis)) {
        stop("runAnalysis must be a logical")
    }
    
    return(0)
}


#' @title Parameters validation for the \code{\link{runSim}} function
#'
#' @description Validation of all parameters needed by the public
#' \code{\link{runSim}} function.
#'
#' @param outputDir a string of \code{character} or \code{NULL}, the path 
#' where the 
#' files created by the function will be saved. When \code{NULL}, the files
#' are saved in the current directory. Default: \code{NULL}.
#'
#' @param fileID a string of \code{character}, a identifiant that will be 
#' included in each output file name. Each output 
#' file name is 
#' composed of those elements, separated by "_":
#' \itemize{ 
#' \item a type name, ex: methylGR, methylObj, etc..
#' \item a \code{fileID}
#' \item the chromosome number, a number between 1 and \code{nbSynCHR}
#' \item the number of samples, a number in the \code{vNbSample} \code{vector}
#' \item the mean proportion of samples that has,
#' for a specific position, differentially methylated values, a 
#' number in the \code{vpDiff} \code{vector}
#' \item the proportion of 
#' C/T for a case differentially methylated that follows a shifted beta 
#' distribution, a
#' number in the \code{vDiff} \code{vector}
#' \item the 
#' proportion of cases that inherits differentially sites, a number in the
#' \code{vInheritance} \code{vector}
#' \item the identifiant for the simulation, a number 
#' between 1 and \code{nbSimulation}
#' \item the file extension ".rds"
#' }
#'
#' @param methData an object of class \code{methylBase}, the CpG information
#' from controls (CTRL) that will be used to create the sythetic chromosome. 
#' The \code{methData} object can also contain information from cases but 
#' only the controls will be used.
#' 
#' @param context a string of \code{character}, the methylation context 
#' string, ex: CpG,CpH,CHH, etc. 
#' Default: \code{"CpG"}.
#' 
#' @param assembly a string of \code{character}, the short description of the 
#' genome assembly. Ex: mm9,hg18 etc.
#' Default: \code{"Rnor_5.0"}
#' 
#' @return \code{0} indicating that the function has been successful.
#'
#' @examples
#'
#' ## Load dataset
#' data("samplesForChrSynthetic")
#'
#' ## The function returns 0 when all paramaters are valid
#' methylInheritanceSim:::validateRunSimOtherParameters(
#' outputDir = "test", fileID = "test", methData = samplesForChrSynthetic, 
#' context = "CpG", assembly = "Rnor_5.0")
#' 
#' @author Pascal Belleau, Astrid Deschenes
#' @importFrom S4Vectors isSingleInteger isSingleNumber
#' @keywords internal
validateRunSimOtherParameters <-function(outputDir, fileID, methData, 
                                    context, assembly) {
    
    ## Validate that the outputDir is an not empty string
    if (!is.null(outputDir) && !is.character(outputDir)) {
        stop("outputDir must be a character string or NULL")
    }
    
    ## Validate that the fileID is an not empty string
    if (!is.null(fileID) && !is.character(fileID)) {
        stop("fileID must be a character string or NULL")
    }
    
    ## Validate that the methData is a methylBase object
    if (!"methylBase" %in% class(methData)) {
        stop("methData must be an object of class \"methyBase\"")
    }
    
    ## Validate that context is a character string
    if(!is.character(context)) {
        stop("context must be a character string")
    }
    
    ## Validate that assembly is a character string
    if (!is.character(assembly)) {
        stop("assembly must be a character string")
    }
    
    return(0)
}
