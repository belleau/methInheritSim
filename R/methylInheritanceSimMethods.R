#' @title TODO
#'
#' @description TODO
#'
#' @param pathOut the path where the the out put file is seved
#'
#' @param fileGen a \code{string} include each output file. Each output 
#' file are 
#' composed with a type name (methylGR, methylObj, ...), _, fileGen (ex F1),
#' parameters of the simulation and ".rds". 
#'
#' @param nbSynCHR a \code{integer}, the number of distinct synthetics 
#' chromosomes. generate.
#'
#' @param methData a \code{methylBase} object from \code{methylKit}, a real 
#' dataset used to sample block for the synthetic chromosome
#'
#' @param nbBlock a \code{integer}, the number of blocks sample from methData 
#' genome to create the synthetic chromosome.
#'
#' @param lBlock a \code{integer}, the number of the CpG in each block sampled
#'
#' @param vNbSample a \code{vector} of the number of CTRL and case in the the 
#' simulation dataset. In 
#' the simulation dataset the number of CTRL equal the number of Case. 
#' The number of CTRL do not need to be equal to the number of Case in
#' the real dataset
#'
#' @param nbGeneration a \code{integer}, the number of generation simulate
#'
#' @param vpDiff a positive \code{double} inferior ot \code{1}, the mean 
#' proportion of case simulate as differentially methylated when the site 
#' is selected as differentially methylated (this is like the penetrance).
#' 
#' @param vpDiffsd a \code{double} the variance of the parameter \code{vpDIff}
#'
#' @param vDiff a positive \code{double} inferior ot \code{1}, in the case of 
#' a differentially methylated site a case 
#' differentially methylated the proportion of C/T follow a beta distribution 
#' where the mean is shifted of vDiff from the CTRL distribution
#'
#' @param vInheritance a positive \code{double} inferior ot \code{1}, the 
#' proportion of case inherite the inherited sites.
#' 
#' @param propInherite a positive \code{double} inferior ot \code{1}, 
#' proportion of differentially methylated site
#' are inherated
#'
#' @param rateDiff a positive \code{double} inferior ot \code{1}, the mean of 
#' the chance that a site is differentially 
#' methylated
#'
#' @param minRate a positive \code{double} inferior ot \code{1}, the minimum 
#' number of diferentially methylated site
#'
#' @param propHetero a positive \code{double} inferior ot \code{1}, the 
#' reduction vDiff for the intergeneration
#' 
#' @param minReads a positive \code{integer} Bases and regions having lower
#' coverage than this count are discarded. The parameter
#' correspond to the \code{lo.count} parameter in the  \code{methylKit} package.
#' 
#' @param maxPercReads a \code{double} between [0,100], the percentile of read
#' counts that is going to be used as upper cutoff. Bases ore regions
#' having higher
#' coverage than this percentile are discarded. Parameter used for both CpG
#' sites and tiles analysis. The parameter
#' correspond to the \code{hi.perc} parameter in the  \code{methylKit} package.
#' 
#' @param context Default: \code{"CpG"}
#' 
#' @param assembly Default: \code{"Rnor_5.0"}
#' 
#' @param meanCov Default: \code{80}
#' 
#' @param keepDiff Default: \code{FALSE}
#' 
#' @param saveGRanges Default: \code{TRUE}
#' 
#' @param saveMethylKit Default: \code{TRUE}
#' 
#' @param anaMethylKit Default: \code{TRUE}
#' 
#' @param nbCores a positive \code{integer}, the number of cores to use when
#' processing the analysis. Default: \code{1}.
#' 
#' @param vSeed a \code{integer}, a seed used when reproducible results are
#' needed. When a value inferior or equal to zero is given, a random integer
#' is used. Default: \code{-1}.
#' 
#' @param vSeed a \code{integer}, a seed used when reproducible results are
#' needed. When a value inferior or equal to zero is given, a random integer
#' is used. Default: \code{-1}.
#'
#' @return TODO
#'
#' @examples
#'
#' ## TODO
#'
#'
#' @author Pascal Belleau
#' @importFrom parallel mclapply
#' @export



runSim <- function(pathOut, fileGen, nbSynCHR, methData, nbBlock, lBlock,
                    nbC, nbGeneration, vpDiff, vpDiffsd, vDiff, vInheritance,
                    propInherite ,rateDiff, minRate, propHetero, minReads = 10, 
                    maxPercReads = 99.9, context = "CpG", assembly="Rnor_5.0",
                    meanCov = 80, n, keepDiff = FALSE,
                    saveGRanges = TRUE, saveMethylKit = TRUE,
                    anaMethylKit = TRUE,
                    nbCores = 1, vSeed = -1){
    
    if (vSeed <= -1) {
        tSeed <- as.numeric(Sys.time())
        vSeed <- 1e8 * (tSeed - floor(tSeed))
    }
    set.seed(vSeed)
    
    for(s in 1:nbSynCHR){
        
        # get the synthetic chr
        res <- getSyntheticChr(methInfo=methData, nbBlock=nbBlock, nbCpG=lBlock)
        adPref <- paste0(fileGen, "_", s)
        saveRDS(res, file=paste0(pathOut, "/stateInfo_", adPref, ".rds"))
        
        for(nbSample in nbC){
            
            adPrefSample <- paste0(adPref, "_", nbSample)
            nbCtrl <- nbSample
            nbCase <- nbSample
            
            # Define tretment and sample.id 
            treatment=c(rep(0,nbSample),rep(1,nbSample))
            sample.id <- list()
            
            for(i in 1:(2*nbSample)){
                if(i <= nbSample){
                    for(j in 1:nbGeneration){
                        if(i == 1){
                            sample.id[[j]] <- list()
                        }
                        sample.id[[j]][[i]] <- paste0("F", j,"_",i,"_C")
                    }
                } else{
                    for(j in 1:nbGeneration){
                        sample.id[[j]][[i]] <- paste0("F", j,"_",i,"_OC")
                    }
                }
            }
            if(keepDiff == TRUE){
                diffRes <- getDiffMeth(stateInfo = res,
                                       rateDiff = rateDiff, minRate = minRate,
                                       propInherite = propInherite)
            } else{
                diffRes <- NULL
            }
            for(i in 1:length(vInheritance)){
                for(j in 1:length(vpDiff)){
                    propDiff <- vpDiff[j]
                    propDiffsd <- vpDiffsd[j]
                    for(k in 1:length(vDiff)){
                        diffValue <- vDiff[k]
                        propInheritance <- ifelse(vInheritance[i]>=0, vInheritance[i], vpDiff[j])
                        prefBase <- paste0(adPrefSample , "_", propDiff, "_", diffValue, "_", propInheritance)
                            
                        a <- mclapply(1:n, FUN = simInheritance, 
                                      pathOut = pathOut, 
                                      pref = prefBase, nbCtrl = nbCtrl,
                                      nbCase = nbCase, treatment = treatment, 
                                      sample.id = sample.id, 
                                      generation = nbGeneration,
                                      stateInfo = res, rateDiff = rateDiff,
                                      minRate = minRate, 
                                      propInherite = propInherite,
                                      diffValue = diffValue, 
                                      propDiff = propDiff,
                                      propDiffsd = propDiffsd, 
                                      propInheritance = propInheritance,
                                      propHetero = propHetero, 
                                      minReads = minReads, 
                                      maxPercReads = maxPercReads,
                                      context = context, 
                                      assembly = assembly,
                                      meanCov = meanCov, diffRes = diffRes,
                                      anaMethylKit = FALSE, 
                                      mc.cores = nbCores,
                                      mc.preschedule = FALSE)
                    }
                }
            }
        }
    }
}



