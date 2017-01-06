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
                   propInherite ,rateDiff, minRate, propHetero, n, 
                   nbCores, keepDiff = FALSE){
    for(s in 1:nbSynCHR){
        
        # get the synthetic chr
        res <- getStateData(methInfo=methData, m=nbBlock, block=lBlock)
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
                diffRes <- getDiffMeth(stateInfo=res,
                                       rateDiff=rateDiff, minRate=minRate,
                                       propInherite=propInherite)
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
                            
                        a <- mclapply(1:n, FUN = simPed, pathOut=pathOut, 
                                      pref=prefBase, nbCtrl=nbCtrl,
                                      nbCase=nbCase, treatment=treatment, 
                                      sample.id=sample.id, 
                                      generation=nbGeneration,
                                      res=res, rateDiff=rateDiff,
                                      minRate=minRate, 
                                      propInherite=propInherite,
                                      diffValue=diffValue, propDiff=propDiff,
                                      propDiffsd=propDiffsd, 
                                      propInheritance=propInheritance,
                                      propHetero=propHetero, diffRes=diffRes,
                                      mc.cores=nbCores,
                                      mc.preschedule = FALSE)
                    }
                }
            }
        }
    }
}


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

reRunSim <- function(pathOut, fileGen, nbSynCHR, methData, nbBlock, lBlock,
                   nbC, nbGeneration, vpDiff, vpDiffsd, vDiff, vInheritance,
                   propInherite ,rateDiff, minRate, propHetero, n, 
                   nbCores, keepDiff = FALSE){
    
    
    
    for(s in 1:nbSynCHR){
        
        res <- NULL
        diffRes <- NULL
        # get the synthetic chr
        
        if( file.exists(paste0(pathOut, "/stateInfo_", adPref, ".rds"))){
            res <- readRDS(paste0(pathOut, "/stateInfo_", adPref, ".rds"))
        } else{
            res <- getStateData(methInfo=methData, m=nbBlock, block=lBlock)
            adPref <- paste0(fileGen, "_", s)
            saveRDS(res, file=paste0(pathOut, "/stateInfo_", adPref, ".rds"))
        }
        
        
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
                diffRes <- getDiffMeth(stateInfo=res,
                                       rateDiff=rateDiff, minRate=minRate,
                                       propInherite=propInherite)
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
                        
                        a <- mclapply(1:n, FUN = restartSim, pathOut=pathOut, 
                                      pref=prefBase, nbCtrl=nbCtrl,
                                      nbCase=nbCase, treatment=treatment, 
                                      sample.id=sample.id, 
                                      generation=nbGeneration,
                                      res=res, rateDiff=rateDiff,
                                      minRate=minRate, 
                                      propInherite=propInherite,
                                      diffValue=diffValue, propDiff=propDiff,
                                      propDiffsd=propDiffsd, 
                                      propInheritance=propInheritance,
                                      propHetero=propHetero, diffRes=diffRes,
                                      mc.cores=nbCores,
                                      mc.preschedule = FALSE)
                    }
                }
            }
        }
    }
}

