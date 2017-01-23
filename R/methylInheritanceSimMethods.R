#' @title TODO
#'
#' @description TODO
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
#' \item the mean proportion of samples that will have,
#' for a specific position, differentially methylated values, a 
#' number in the \code{vpDiff} \code{vector}
#' \item the proportion of 
#' C/T for a case differentially methylated follow a beta distribution, a
#' number in the \code{vDiff} \code{vector}
#' \item the proportion of 
#' C/T for a case differentially methylated follow a beta distribution, a
#' number in the \code{vDiff} \code{vector}
#' \item the 
#' proportion of cases that inherited differentially sites, a number in the
#' \code{vInheritance} \code{vector}
#' \item TODO
#' \item the file extension ".rds"
#' }
#'
#' @param nbSynCHR a positive \code{integer}, the number of distinct synthetic 
#' chromosomes that will be generated. Default: \code{1}.
#'
#' @param methData an object of class \code{methylBase}, the CpG information
#' from controls (CTRL) that will be used to create the sythetic chromosome. 
#' The \code{methData} object can also contain information from cases but 
#' only the controls will be used.
#'
#' @param nbBlock a positive \code{integer}, the number of blocks used 
#' for sampling.
#'
#' @param lBlock a positive \code{integer}, the number of consecutive CpG 
#' positions used for sampling from \code{methInfo}.
#'
#' @param vNbSample a \code{vector} of positive \code{integer}, the number of 
#' controls (CTRL) and cases in the simulated dataset. In 
#' the simulated dataset, the number of CTRL equals the number of cases. 
#' The number of CTRL do not need to be equal to the number of Case in
#' the real dataset \code{methData}.
#'
#' @param nbGeneration a positive \code{integer}, the number of generations.
#'
#' @param vpDiff a \code{vector} of \code{double} superior to \code{0} and 
#' inferior or equal 
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
#' @param vInheritance a \code{vector} of positive \code{double} between [0,1], 
#' the proportion of cases that inherited differentially sites.
#' 
#' @param propInherite a non-negative \code{double} inferior or equal 
#' to \code{1}, 
#' proportion of differentially methylated site
#' are inherated
#'
#' @param rateDiff a positive \code{double} inferior to \code{1}, the mean of 
#' the chance that a site is differentially 
#' methylated.
#'
#' @param minRate a non-negative \code{double} inferior to \code{1}, the 
#' minimum rate for differentially methylated sites.
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
#' genome assembly. Ex: "mm9", "hg18", etc.
#' Default: \code{"Rnor_5.0"}
#' 
#' @param meanCov a positive \code{integer}, the mean of the coverage
#' at the CpG sites. Default: \code{80}.
#' 
#' @param n a positive \code{integer}, the number of simulation for each 
#' parameter (\code{vNbSample}, \code{vpDiff}, \code{vDiff} and
#' \code{vInheritance}).
#' 
#' @param keepDiff a \code{logical}, when \code{TRUE}, the 
#' differentially methyled sites
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
#' creating the simulated datasets. Default: \code{1} and always 
#' \code{1} for Windows.
#' 
#' @param vSeed a \code{integer}, a seed used when reproducible results are
#' needed. When a value inferior or equal to zero is given, a random integer
#' is used. Default: \code{-1}.
#'
#' @return TODO
#'
#' @examples
#'
#' ## TODO outputDir where to write file ?
#' ## Load methyl information
#' data(samplesForChrSynthetic)
#' 
#' \dontrun{runSim(outputDir = "testData", fileID = "F1", nbSynCHR = 1, 
#' methData = samplesForChrSynthetic, nbBlock = 10, lBlock = 20,
#' vNbSample = c(6), nbGeneration = 3, vpDiff = c(0.9), 
#' vpDiffsd = c(0.1), vDiff = c(0.8), 
#' vInheritance = c(0.5), propInherite = 0.3,
#' rateDiff = 0.3, minRate = 0.2, propHetero = 0.5, n = 5, 
#' nbCores= 1, vSeed = 32)}
#' 
#' @author Pascal Belleau
#' @importFrom parallel mclapply
#' @export
runSim <- function(outputDir = NULL, fileID, nbSynCHR = 1, methData, 
                    nbBlock, lBlock,
                    vNbSample, nbGeneration, vpDiff, vpDiffsd, vDiff, 
                    vInheritance,
                    propInherite, rateDiff, minRate, propHetero, 
                    minReads = 10, 
                    maxPercReads = 99.9, context = "CpG", assembly="Rnor_5.0",
                    meanCov = 80, n, keepDiff = FALSE,
                    saveGRanges = TRUE, saveMethylKit = TRUE,
                    anaMethylKit = TRUE,
                    nbCores = 1, vSeed = -1) {
    
    validateRunSimParameters(outputDir = outputDir, fileID = fileID, 
                                nbSynCHR = nbSynCHR, methData = methData, 
                                nbBlock = nbBlock, lBlock  = lBlock,
                                vNbSample = vNbSample, 
                                nbGeneration = nbGeneration, 
                                vpDiff = vpDiff, vpDiffsd = vpDiffsd, 
                                vDiff = vDiff, 
                                vInheritance = vInheritance,
                                propInherite = propInherite, 
                                rateDiff, minRate, propHetero, 
                                minReads = minReads, 
                                maxPercReads = maxPercReads, 
                                context = context, assembly = assembly,
                                meanCov = meanCov, n, keepDiff = keepDiff,
                                saveGRanges = saveGRanges, 
                                saveMethylKit = saveMethylKit,
                                anaMethylKit = anaMethylKit,
                                nbCores = nbCores, vSeed = vSeed)
    
    ## Fix seed
    if (vSeed <= -1) {
        tSeed <- as.numeric(Sys.time())
        vSeed <- 1e8 * (tSeed - floor(tSeed))
    }
    set.seed(vSeed)
    
    if (!is.null(outputDir) && !dir.exists(outputDir)) {
        dir.create(outputDir, showWarnings = TRUE)
    }
    
    for(s in 1:nbSynCHR) {
        
        # Create a synthetic chromosome
        res <- getSyntheticChr(methInfo = methData, nbBlock = nbBlock, 
                                nbCpG = lBlock)
        
        adPref <- paste0(fileID, "_", s)
        saveRDS(res, file = paste0(outputDir, "/stateInfo_", adPref, ".rds"))
        
        for(nbSample in vNbSample) {
            
            adPrefSample <- paste0(adPref, "_", nbSample)
            nbCtrl <- nbSample
            nbCase <- nbSample
            
            # Define tretment and sample.id 
            treatment <- c(rep(0,nbSample), rep(1,nbSample))
            if(saveGRanges){
                saveRDS(treatment, file = paste0(outputDir, "/treatment_", 
                                            adPrefSample, ".rds"))
            }
            
            sample.id <- list()
            
            for (i in 1:(2*nbSample)) {
                if (i <= nbSample) {
                    for(j in 1:nbGeneration) {
                        if(i == 1){
                            sample.id[[j]] <- list()
                        }
                        sample.id[[j]][[i]] <- paste0("F", j, "_", i, "_C")
                    }
                } else {
                    for(j in 1:nbGeneration){
                        sample.id[[j]][[i]] <- paste0("F", j, "_", i, "_OC")
                    }
                }
            }
            if (keepDiff == TRUE) {
                diffRes <- getDiffMeth(stateInfo = res,
                                        rateDiff = rateDiff, minRate = minRate,
                                        propInherite = propInherite)
            } else {
                diffRes <- NULL
            }
            
            for (i in 1:length(vInheritance)) {
                for(j in 1:length(vpDiff)) {
                    propDiff <- vpDiff[j]
                    propDiffsd <- vpDiffsd[j]
                    for (k in 1:length(vDiff)) {
                        diffValue <- vDiff[k]
                        
                        propInheritance <- ifelse(vInheritance[i] >= 0, 
                                                vInheritance[i], vpDiff[j])
                        
                        prefBase <- paste0(adPrefSample , "_", propDiff, 
                                            "_", diffValue, "_", 
                                            propInheritance)
                            
                        a <- mclapply(1:n, FUN = simInheritance, 
                                        pathOut = outputDir, 
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
                                        saveGRanges = saveGRanges,
                                        saveMethylKit = saveMethylKit,
                                        anaMethylKit = anaMethylKit, 
                                        mc.cores = nbCores,
                                        mc.preschedule = FALSE)
                    }
                }
            }
        }
    }
    
    return(0)
}
