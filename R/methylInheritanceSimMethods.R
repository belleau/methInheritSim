#' @title TODO
#'
#' @description TODO
#'
#' @param pathOut
#'
#' @param fileGen
#'
#' @param nbSynCHR
#'
#' @param methData
#'
#' @param nbBlock
#'
#' @param lBlock
#'
#' @param nbC
#'
#' @param nbGeneration
#'
#' @param vpDiff
#'
#' @param vpDiffsd
#'
#' @param vDiff
#'
#' @param vInheritance
#' 
#' @param propInherite
#'
#' @param rateDiff
#'
#' @param minRate
#'
#' @param propHetero
#' 
#' @param nbCores
#' 
#' @param keepDiff
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
                   meanCov = 80, n, 
                   nbCores = 1, keepDiff = FALSE, vSeed = -1){
    
    if (vSeed <= -1) {
        tSeed <- as.numeric(Sys.time())
        vSeed <- 1e8 * (tSeed - floor(tSeed))
    }
    set.seed(vSeed)
    
    for(s in 1:nbSynCHR){
        
        # get the synthetic chr
        res <- getSyntheticChr(methInfo=methData, m=nbBlock, nbCpG=lBlock)
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



