#' @title Simulate a multigeneration methylation experiment with inheritance 
#' 
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
#' Default: \code{"s"}.
#' 
#' @param nbSynCHR a positive \code{integer}, the number of distinct synthetic 
#' chromosomes that will be generated. Default: \code{1}.
#'
#' @param methData an object of class \code{methylBase}, the CpG information
#' from controls (CTRL) that will be used to create the synthetic chromosome. 
#' The \code{methData} object can also contain information from cases but 
#' only the controls are used.
#' 
#' @param nbSimulation a positive \code{integer}, the number of simulations 
#' generated 
#' for each parameter (\code{vNbSample}, \code{vpDiff}, \code{vDiff} and
#' \code{vInheritance}). 
#' The total number of simulation is 
#' nbSimulation * \code{length(vNbSample)} * \code{length(vpDiff)} *
#' \code{length(vInheritance)})
#' Default: \code{10}.
#'
#' @param nbBlock a positive \code{integer}, the number of blocks used 
#' for sampling.
#' Default: \code{100}.
#'
#' @param nbCpG a positive \code{integer}, the number of consecutive CpG 
#' positions used for sampling from \code{methInfo}.
#' Default: \code{50}.
#'
#' @param nbGeneration a positive \code{integer}, the number of generations
#' simulated.
#' Default: \code{3}.
#'
#' @param vNbSample a \code{vector} of distinct positive \code{integer}, 
#' the number of controls (CTRL) and cases in the simulated dataset. In 
#' the simulated dataset, the number of CTRL equals the number of cases. 
#' The number of CTRL do not need to be equal to the number of Case in
#' the real \code{methData} dataset.
#' Default: \code{c(3, 6)}.
#'
#' @param vpDiff a \code{vector} of distinct \code{double} superior to 
#' \code{0} and inferior or equal 
#' to \code{1}, the mean value for the proportion of samples that will have,
#' for a specific position, differentially methylated values. It can be 
#' interpreted as the penetrance. Note that \code{vpDiff} and \code{vpDiffsd}
#'  must be the same length.
#' Default: \code{c(0.9)}.
#' 
#' @param vpDiffsd a \code{vector} of a non-negative \code{double}, the standard
#' deviation associated to the \code{vpDiff}. Note vpDiff and vpDiffsd must be
#' the same length.
#' Default: \code{c(0.1)}
#'
#' @param vDiff a \code{vector} of distinct non-negative \code{double} include 
#' in [0,1], the proportion of C/T for a case differentially methylated follow 
#' a beta distribution where the mean is shifted of \code{vDiff} 
#' from the CTRL distribution
#' Default: \code{c(0.8)}
#'
#' @param vInheritance a \code{vector} of distinct non-negative \code{double} 
#' include in [0,1], the proportion of cases 
#' that inherited differentially sites.
#' Default: \code{c(0.5)}
#' 
#' @param rateDiff a positive \code{double} inferior to \code{1}, the mean of 
#' the chance that a site is differentially 
#' methylated.
#' Default: \code{0.01}
#'
#' @param minRate a non-negative \code{double} inferior to \code{1}, the 
#' minimum rate for differentially methylated sites.
#' Default: \code{0.01}
#'
#' @param propInherite a non-negative \code{double} inferior or equal 
#' to \code{1}, 
#' proportion of differentially methylated region
#' are inherated
#' Default: \code{0.3}
#'
#' @param propHetero a non-negative \code{double} between [0,1], the 
#' reduction of vDiff for the second and following generations.
#' Default: \code{0.5}
#' 
#' @param minReads a positive \code{integer} Bases and regions having lower
#' coverage than this count are discarded. The parameter
#' correspond to the \code{lo.count} parameter in the \code{methylKit} package.
#' Default: \code{10}
#' 
#' @param maxPercReads a \code{double} between [0,100], the percentile of read
#' counts that is going to be used as upper cutoff. Bases ore regions
#' having higher
#' coverage than this percentile are discarded. Parameter used for both CpG
#' sites and tiles analysis. The parameter
#' correspond to the \code{hi.perc} parameter in the  \code{methylKit} package.
#' Default: \code{99.9}
#' 
#' 
#' @param meanCov a positive \code{integer}, the mean of the coverage
#' at the CpG sites.
#' Default: \code{80}
#' 
#' 
#' @param context a string of \code{character}, the methylation context 
#'  one of the CpG,CpH,CHH, none. 
#' Default: \code{"CpG"}.
#' 
#' @param assembly a string of \code{character}, the short description of the 
#' genome assembly. Ex: "mm9", "hg18", etc.
#' Default: \code{"Rnor_5.0"}
#' 
#' @param keepDiff a \code{logical}, if \code{TRUE}, the 
#' differentially methyled sites
#' will be the same for each parameter (\code{vpDiff}, 
#' \code{vDiff} and \code{vInheritance}).
#' Default: \code{FALSE}
#' 
#' @param saveGRanges a \code{logical}, if \code{true}, for each simulation, 
#' save a file that contains a 
#' \code{list} of \code{GRangesList} (each \code{GRangesList} corresponding to 
#' a generation). Each \code{GRangeList} contain \code{2 * vNbSample} 
#' \code{GRanges} that contains for each CpG site: the position, the coverage 
#' and the proportion of the C/T.
#' 
#' A files treatment which contains the position of the CTRL and Case in 
#' the GRangesList are 
#' saved too. TODO Default: \code{TRUE}.
#' 
#' @param saveMethylKit a \code{logical}, if \code{TRUE}, for each simulation, 
#' save a file that contains a list of \code{methylRawList} 
#' (each \code{methylRawList} corresponding to a generation) representing
#' the simulation.
#' The file name is describe in the parameter \code{fileID},
#' in this case the file name start with: methylObj
#' Default: \code{TRUE}.
#' 
#' @param runAnalysis a \code{logical}, if \code{TRUE}, for each simulation run
#' a differentially methylation sites analysis with methylKit on each 
#' generation and save two files The first contains
#' a list of normalised \code{mehylBase} and the second a list of
#' \code{methylDiff} (each \code{methylDiff} is corresponding ) that contains 
#' the comparaison between case and control.
#' The file name is describe in the parameter \code{fileID},
#' in this case the file name start with: meth and methDiff
#' Default: \code{FALSE}
#' 
#' @param nbCores a positive \code{integer}, the number of cores to use when
#' creating the simulated datasets. Default: \code{1} and always 
#' \code{1} for Windows.
#' 
#' @param vSeed a \code{integer}, a seed used when reproducible results are
#' needed. When a value inferior or equal to zero is given, a random integer
#' is used. Default: \code{-1}.
#'
#' @return \code{0} indicating that all parameters validations have been
#' successful
#' 
#' @examples
#'
#' ## Load dataset containing methyl information
#' data(samplesForChrSynthetic)
#' 
#' ## Set the output directory where files will be created
#' temp_dir <- "test_runSim"
#' 
#' ## Create 2 simulated dataset (nbSimulation) 
#' ## over 3 generations (nbGenration = 3) with
#' ## 6 cases and 6 controls (nNbsample = 6) using only one set
#' ## of parameters (vpDiff = 0.9, vpDiffsd = 0.1, vDiff = 0.8)
#' \dontrun{runSim(outputDir = temp_dir, fileID = "F1", nbSynCHR = 1, 
#' methData = samplesForChrSynthetic, nbSimulation = 2, 
#' nbBlock = 10, nbCpG = 20,
#' nbGeneration = 3, vNbSample = c(6), vpDiff = c(0.9), 
#' vpDiffsd = c(0.1), vDiff = c(0.8), 
#' vInheritance = c(0.5), propInherite = 0.3,
#' rateDiff = 0.3, minRate = 0.2, propHetero = 0.5, 
#' nbCores = 1, vSeed = 32)}
#' 
#' ## To delete the output directory and its content
#' \dontrun{if (dir.exists(temp_dir)) {
#' unlink(temp_dir, recursive = TRUE, force = FALSE)
#' }}
#' 
#' @author Pascal Belleau
#' @importFrom parallel mclapply
#' @export
runSim <- function(outputDir = NULL, fileID = "s", 
                    nbSynCHR = 1, methData, 
                    nbSimulation = 10, nbBlock = 100, nbCpG = 50,
                    nbGeneration = 3, vNbSample = c(3, 6), vpDiff = c(0.9), 
                    vpDiffsd = c(0.1), vDiff = c(0.8), 
                    vInheritance = c(0.5),
                    rateDiff = 0.01, minRate = 0.01, propInherite = 0.3,
                    propHetero = 0.5, 
                    minReads = 10, 
                    maxPercReads = 99.9, meanCov = 80,
                    context = "CpG", assembly="Rnor_5.0",
                    keepDiff = FALSE,
                    saveGRanges = TRUE, saveMethylKit = TRUE,
                    runAnalysis = FALSE,
                    nbCores = 1, vSeed = -1) {
    
    validateRunSimParameters(outputDir = outputDir, fileID = fileID, 
                                nbSynCHR = nbSynCHR, methData = methData,
                                nbSimulation = nbSimulation,
                                nbBlock = nbBlock, nbCpG  = nbCpG,
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
                                meanCov = meanCov, keepDiff = keepDiff,
                                saveGRanges = saveGRanges, 
                                saveMethylKit = saveMethylKit,
                                runAnalysis = runAnalysis,
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
        
        # Create synthetic chromosome
        res <- getSyntheticChr(methInfo = methData, nbBlock = nbBlock, 
                                nbCpG = nbCpG)
        
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
                            
                        a <- mclapply(1:nbSimulation, FUN = simInheritance, 
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
                                        runAnalysis = runAnalysis, 
                                        mc.cores = nbCores,
                                        mc.preschedule = FALSE)
                    }
                }
            }
        }
    }
    
    return(0)
}
