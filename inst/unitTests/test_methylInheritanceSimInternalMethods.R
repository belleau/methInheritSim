###################################################
##
## Test the methylInheritanceSimInternalMethods functions
##
###################################################

###  Test setup

library( "GenomicRanges" )

data(samplesForChrSynthetic)
data(dataSimExample)

###################################################
## estBetaAlpha() function
###################################################

test.estBetaAlpha_good_01 <- function() {
    
    obs <- methylInheritanceSim:::estBetaAlpha(c(0.5, 0.2), 0.00001)
    
    exp <- 0.125
    
    message <- paste0("test.estBetaAlpha_good_01() ",
                    "- Valid parameters did not generated expected results.")
    
    checkEquals(obs, exp, message)
}


###################################################
## estBetaBeta() function
###################################################

test.estBetaBeta_good_01 <- function() {
    
    obs <- methylInheritanceSim:::estBetaBeta(c(0.3, 0.2), 0.00001)
    
    exp <- 0.035
    
    message <- paste0("test.estBetaBeta_good_01() ",
                      "- Valid parameters did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

###################################################
## getDiffCase() function
###################################################

test.getDiffCase_good_01 <- function() {
    set.seed(322)
    
    x <- c(0.14562, 0.0003607153, 1)
    obs <- methylInheritanceSim:::getDiffCase(x = x, nb = 4, sDiff = 0.8, 
                                              diffCase = 3)
    exp <- c(0.945620000000, 3.000000000000, 1.000000000000, 0.947694615429, 
             0.965193968711, 0.906084052941, 0.122224066759)
    
    message <- paste0("test.getDiffCase_good_01() ",
                    "- Valid parameters did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.getDiffCase_good_no_DMS <- function() {
    set.seed(22)
    
    x <- c(0.14562, 0.0003607153, 0)
    obs <- methylInheritanceSim:::getDiffCase(x = x, nb = 4, sDiff = 0.8, 
                                              diffCase = 3)
    exp <- c(0.1456200000, 0.000000000000, 4.000000000000, 0.1348148850, 
             0.1698118247, 0.1520623197, 0.1411707002)
    
    message <- paste0("test.getDiffCase_good_no_DMS() ",
                    "- Valid parameters with no DMS did not generated expected results.")
    
    checkEquals(obs, exp, message)
}


###################################################
## getSim() function
###################################################

test.getSim_good_01 <- function() {
    
    set.seed(22212)
    
    stateInformation <- methylInheritanceSim:::getSyntheticChr(methInfo = samplesForChrSynthetic, 
                                                                nbBlock = 1, nbCpG = 3)
    
    stateDiff <- list()
    stateDiff[["stateDiff"]] <- c(1, 0, 1)
    stateDiff[["stateInherite"]] <- c(1, 0, 0)

    obs <- methylInheritanceSim:::getSim(nbCtrl = 2, nbCase = 3, generation = 2, 
                stateInfo = stateInformation, stateDiff = stateDiff, diffValue = 10, 
                propDiff = 0.8, propDiffsd = 0.2, propInheritance = 0.8, propHetero = 0.1)
    
    exp <- GRangesList()
    
    exp[[1]] <- GRanges(seqnames = seqnames(stateInformation),
                        ranges = ranges(stateInformation),
                        strand =  strand(stateInformation),
                        meanDiff = c(0.000000000000000, 0.288519632432523, 1.000000000000000), 
                        meanCTRL = mcols(stateInformation)[3],
                        partitionCase = c(2, 0, 2), partitionCtrl = c(1, 3, 1),
                        ctrl.V1 = c(0.599248212271003, 0.104673052925196, 0.0163259007084415),
                        ctrl.V2 = c(0.660878027594768, 0.0737704350938676, 0.00412068310774805),
                        case.V1 = c(0.000000000000000, 0.297550189486176, 1.000000000000000),
                        case.V2 = c(0.000000000000000, 0.152428811636903, 1.000000000000000),
                        case.V3 = c(0.140108118897419, 0.22397752927231, 0.0242355851971473))
    
    exp[[2]] <- GRanges(seqnames = seqnames(stateInformation),
                        ranges = ranges(stateInformation),
                        strand =  strand(stateInformation),
                        meanDiff = c(0.000000000000000, 0.288519632432523, 0.0103452693986541), 
                        meanCTRL = mcols(stateInformation)[3],
                        partitionCase = c(2, 0, 0), partitionCtrl = c(1, 3, 3),
                        ctrl.V1 = c(0.0851227837846553, 0.281923948948365, 0.0198953426761693),
                        ctrl.V2 = c(0.393339608860298, 0.118756249033386, 0.019181448769631),
                        case.V1 = c(1.000000000000000, 0.274490475976465, 0.0345473331486896),
                        case.V2 = c(0.000000000000000, 0.140418843034454, 0.000824070856929141),
                        case.V3 = c(0.935659881913647, 0.282094425404783, 0.0143114682779736))
    
    message <- paste0("test.getSim_good_01() ",
                      "- Valid parameters for getSim() did not generated expected results.")
    
    checkEquals(obs, exp, message)
}


###################################################
## getDiffMeth() function
###################################################

test.getDiffMeth_good_01 <- function() {
    set.seed(3222)
    
    t<-dataSimExample$stateInfo[1:10,]
    
    obs <- methylInheritanceSim:::getDiffMeth(stateInfo = 
                    t, rateDiff = 0.3, minRate = 0.1,
                    propInherite = 0.2)
                                                            
    exp <- list()
    exp$stateDiff <- c(0, 0, 0, 0, 0, 1, 1, 0, 1, 1)
    exp$stateInherite <- c(0, 0, 0, 0, 0, 1, 0, 0, 0, 0)
    
    message <- paste0("test.getDiffMeth_good_01() ",
                      "- Valid parameters did not generated expected results.")
    
    checkEquals(obs, exp, message)
}


###################################################
## getSyntheticChr() function
###################################################

test.getSyntheticChr_good_01 <- function() {
    
    set.seed(322)
    
    obs <- methylInheritanceSim:::getSyntheticChr(methInfo = 
                                samplesForChrSynthetic, nbBlock = 1, nbCpG = 2)
    
    exp <- GenomicRanges::GRanges(seqnames = rep("S", 2), 
                   ranges = IRanges::IRanges(start = c(1000, 2514), 
                                                end = c(1000, 2514)),
                   strand = rep("+", 2), chrOri = rep(1, 2),
                   startOri = c(11690624, 11692138), 
                   meanCTRL = c(0.934017763674095, 0.957297904589756),
                   varCTRL = c(0.00293610808482296, 0.00038750651540637))
    
    message <- paste0("test.getSyntheticChr_good_01() ",
                      "- Valid parameters did not generated expected results.")
    
    checkEquals(obs, exp, message)
}


###################################################
## validateRunSimOtherParameters() function
###################################################

test.validateRunSimOtherParameters_outputDir_number <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimOtherParameters(outputDir = 33,
                            fileID = "F1", methData = samplesForChrSynthetic,
                            context = "CpG", assembly = "hg19"),
                    error=conditionMessage)
    
    exp <- "outputDir must be a character string or NULL"
    
    message <- paste0("test.validateRunSimOtherParameters_outputDir_number() ",
                      "- Number as outputDir parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimOtherParameters_fileID_number <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimOtherParameters(outputDir = "test",
                            fileID = 2, methData = samplesForChrSynthetic,
                            context = "CpG", assembly = "hg19"),
                    error=conditionMessage)
    
    exp <- "fileID must be a character string or NULL"
    
    message <- paste0("test.validateRunSimOtherParameters_fileID_number() ",
                      "- Number as fileID parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimOtherParameters_methBase_number <- function() {
    obs <- tryCatch(
        methylInheritanceSim:::validateRunSimOtherParameters(outputDir = "test",
                        fileID = "F1", methData = 33, context = "CpG", assembly = "hg19"),
        error=conditionMessage)
    
    exp <- "methData must be an object of class \"methyBase\""
    
    message <- paste0("test.validateRunSimOtherParameters_methBase_number() ",
                      "- Number as methBase parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimOtherParameters_assembly_double <- function() {
    obs <- tryCatch(
        methylInheritanceSim:::validateRunSimOtherParameters(outputDir = "test",
                    fileID = "F1", methData = samplesForChrSynthetic,
                    context = "CpG", assembly = 0.44),
        error=conditionMessage)
    
    exp <- "assembly must be a character string"
    
    message <- paste0("test.validateRunSimOtherParameters_assembly_double() ",
                    "- Double as assembly parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimOtherParameters_context_double <- function() {
    obs <- tryCatch(
        methylInheritanceSim:::validateRunSimOtherParameters(outputDir = "test",
                                                        fileID = "F1", 
                                                        methData = samplesForChrSynthetic,
                                                        context = 0.332, assembly = "hg19"),
        error=conditionMessage)
    
    exp <- "context must be a character string"
    
    message <- paste0("test.validateRunSimOtherParameters_context_double() ",
                      "- Double as context parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}


###################################################
## validateRunSimLogicalParameters() function
###################################################


test.validateRunSimLogicalParameters_keepDiff_double <- function() {
    obs <- tryCatch(
        methylInheritanceSim:::validateRunSimLogicalParameters(keepDiff = 0.22, 
                                saveGRanges = FALSE, saveMethylKit = FALSE,
                                runAnalysis = FALSE),
        error=conditionMessage)
    
    exp <- "keepDiff must be a logical"
    
    message <- paste0("test.validateRunSimLogicalParameters_keepDiff_double() ",
                      "- Double as keepDiff parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimLogicalParameters_saveGRanges_double <- function() {
    obs <- tryCatch(
        methylInheritanceSim:::validateRunSimLogicalParameters(keepDiff = TRUE, 
                                    saveGRanges = 0.2, saveMethylKit = FALSE,
                                    runAnalysis = FALSE),
        error=conditionMessage)
    
    exp <- "saveGRanges must be a logical"
    
    message <- paste0("test.validateRunSimLogicalParameters_saveGRanges_double() ",
                      "- Double as saveGRanges parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}


test.validateRunSimLogicalParameters_runAnalysis_double <- function() {
    obs <- tryCatch(
        methylInheritanceSim:::validateRunSimLogicalParameters(keepDiff = TRUE, 
                                saveGRanges = FALSE, saveMethylKit = FALSE,
                                runAnalysis = 0.01),
        error=conditionMessage)
    
    exp <- "runAnalysis must be a logical"
    
    message <- paste0("test.validateRunSimLogicalParameters_runAnalysis_double() ",
                      "- Double as runAnalysis parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}


test.validateRunSimLogicalParameters_saveMethylKit_double <- function() {
    obs <- tryCatch(
        methylInheritanceSim:::validateRunSimLogicalParameters(keepDiff = TRUE, 
                saveGRanges = FALSE, saveMethylKit = 0.2, runAnalysis = FALSE),
        error=conditionMessage)
    
    exp <- "saveMethylKit must be a logical"
    
    message <- paste0("test.validateRunSimLogicalParameters_saveMethylKit_double() ",
                      "- Double as saveMethylKit parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimLogicalParameters_runAnalysis_double <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimLogicalParameters(keepDiff = TRUE, 
                saveGRanges = FALSE, saveMethylKit = FALSE, runAnalysis = 0.33),
            error=conditionMessage)
    
    exp <- "runAnalysis must be a logical"
    
    message <- paste0("test.validateRunSimLogicalParameters_runAnalysis_double() ",
                      "- Double as runAnalysis parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimLogicalParameters_keepDiff_double <- function() {
    obs <- tryCatch(
        methylInheritanceSim:::validateRunSimLogicalParameters(keepDiff = 0.22, 
                saveGRanges = FALSE, saveMethylKit = FALSE, runAnalysis = FALSE),
        error=conditionMessage)
    
    exp <- "keepDiff must be a logical"
    
    message <- paste0("test.validateRunSimLogicalParameters_keepDiff_double() ",
                      "- Double as keepDiff parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}


###################################################
## validateRunSimNumberParameters() function
###################################################

test.validateRunSimNumberParameters_nbSynCHR_not_number <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimNumberParameters(nbSynCHR = "hi", 
                nbSimulation = 10, nbBlock = 3, nbCpG  = 2, vNbSample = 2, 
                nbGeneration = 3, vpDiff = 0.2, vpDiffsd = 0.1, vDiff = 0.8, 
                vInheritance = 0.5, propInherite = 0.8, rateDiff = 2, 
                minRate = 1, propHetero = 0.4, minReads = 2, maxPercReads = 99.9, 
                meanCov = 10, nbCores = 1, vSeed = -1),
            error=conditionMessage)
    
    exp <- "nbSynCHR must be a positive integer or numeric"
    
    message <- paste0("test.validateRunSimNumberParameters_nbSynCHR_not_number() ",
                      "- Not a number as nbSynCHR parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimNumberParameters_nbSynCHR_vector_number <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimNumberParameters(nbSynCHR = c(1,2), 
                nbSimulation = 10, nbBlock = 3, nbCpG  = 2,
                vNbSample = 2, nbGeneration = 3, vpDiff = 0.2, vpDiffsd = 0.1, vDiff = 0.8, 
                vInheritance = 0.5, propInherite = 0.8, rateDiff = 2, minRate = 1, 
                propHetero = 0.4, minReads = 2, maxPercReads = 99.9, meanCov = 10, 
                nbCores = 1, vSeed = -1),
            error=conditionMessage)
    
    exp <- "nbSynCHR must be a positive integer or numeric"
    
    message <- paste0("test.validateRunSimNumberParameters_nbSynCHR_vector_number() ",
                      "- Vector of numbers as nbSynCHR parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimNumberParameters_nbSynCHR_zero <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimNumberParameters(nbSynCHR = 0, 
                nbSimulation = 10, nbBlock = 3, nbCpG  = 2, vNbSample = 2, nbGeneration = 3, 
                vpDiff = 0.22, vpDiffsd = 0.1, vDiff = 0.8, vInheritance = 0.5, 
                propInherite = 0.8, rateDiff = 2, minRate = 1, propHetero = 0.4, 
                minReads = 2, maxPercReads = 99.9, 
                meanCov = 10, nbCores = 1, vSeed = -1),
            error=conditionMessage)
    
    exp <- "nbSynCHR must be a positive integer or numeric"
    
    message <- paste0("test.validateRunSimNumberParameters_nbSynCHR_zero() ",
                      "- Zero as nbSynCHR parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimNumberParameters_nbSimulation_not_number <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimNumberParameters(nbSynCHR = 1, 
                nbSimulation = "Hi", nbBlock = 3, nbCpG  = 2, vNbSample = 2, 
                nbGeneration = 3, vpDiff = 0.2, vpDiffsd = 0.1, vDiff = 0.8, 
                vInheritance = 0.5, propInherite = 0.8, rateDiff = 2, minRate = 1, 
                propHetero = 0.4, minReads = 2, maxPercReads = 99.9, meanCov = 10, 
                nbCores = 1, vSeed = -1),
            error=conditionMessage)
    
    exp <- "nbSimulation must be a positive integer or numeric"
    
    message <- paste0("test.validateRunSimNumberParameters_nbSimulation_not_number() ",
                      "- Not a number as nbSimulation parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimNumberParameters_nbSimulation_not_number <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimNumberParameters(nbSynCHR = 1, 
                nbSimulation = "Hi", nbBlock = 3, nbCpG  = 2, vNbSample = 2, 
                nbGeneration = 3, vpDiff = 0.2, vpDiffsd = 0.1, vDiff = 0.8, 
                vInheritance = 0.5, propInherite = 0.8, rateDiff = 2,  minRate = 1, 
                propHetero = 0.4, minReads = 2, maxPercReads = 99.9, meanCov = 10, 
                nbCores = 1, vSeed = -1),
            error=conditionMessage)
    
    exp <- "nbSimulation must be a positive integer or numeric"
    
    message <- paste0("test.validateRunSimNumberParameters_nbSimulation_not_number() ",
                      "- Not a number as nbSimulation parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimNumberParameters_nbSimulation_vector_number <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimNumberParameters(nbSynCHR = 1, 
                        nbSimulation = c(1,2), nbBlock = 3, nbCpG  = 2, vNbSample = 2, 
                        nbGeneration = 3, vpDiff = 0.2, vpDiffsd = 0.1, vDiff = 0.8, 
                        vInheritance = 0.5, propInherite = 0.8, rateDiff = 2, 
                        minRate = 1, propHetero = 0.4, minReads = 2, maxPercReads = 99.9, 
                        meanCov = 10, nbCores = 1, vSeed = -1),
            error=conditionMessage)
    
    exp <- "nbSimulation must be a positive integer or numeric"
    
    message <- paste0("test.validateRunSimNumberParameters_nbSimulation_vector_number() ",
                      "- Vector of numbers as nbSimulation parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimNumberParameters_nbSimulation_zero <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimNumberParameters(nbSynCHR = 1, 
                    nbSimulation = 0, nbBlock = 3, nbCpG  = 2, vNbSample = 2, 
                    nbGeneration = 3, vpDiff = 0.22, vpDiffsd = 0.1, vDiff = 0.8, 
                    vInheritance = 0.5, propInherite = 0.8, rateDiff = 2, minRate = 1, 
                    propHetero = 0.4, minReads = 2, maxPercReads = 99.9, 
                    meanCov = 10, nbCores = 1, vSeed = -1),
            error=conditionMessage)
    
    exp <- "nbSimulation must be a positive integer or numeric"
    
    message <- paste0("test.validateRunSimNumberParameters_nbSimulation_zero() ",
                      "- Zero as nbSimulation parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimNumberParameters_nbBlock_number_vector <- function() {
    obs <- tryCatch(
        methylInheritanceSim:::validateRunSimNumberParameters(nbSynCHR = 1, 
                    nbSimulation = 10, nbBlock = c(3, 1), nbCpG  = 2, vNbSample = 2, 
                    nbGeneration = 3, vpDiff = 0.22, vpDiffsd = 0.1, vDiff = 0.8, 
                    vInheritance = 0.5, propInherite = 0.8, rateDiff = 2, 
                    minRate = 1, propHetero = 0.4, minReads = 2, maxPercReads = 99.9, 
                    meanCov = 10, nbCores = 1, vSeed = -1),
        error=conditionMessage)
    
    exp <- "nbBlock must be a positive integer or numeric"
    
    message <- paste0("test.validateRunSimNumberParameters_nbBlock_number_vector() ",
                      "- Number vector as nbBlock parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimNumberParameters_nbBlock_string <- function() {
    obs <- tryCatch(
        methylInheritanceSim:::validateRunSimNumberParameters(nbSynCHR = 1, 
                    nbSimulation = 10, nbBlock = "hi", nbCpG  = 2, vNbSample = 2, 
                    nbGeneration = 3, vpDiff = 0.22, vpDiffsd = 0.1, vDiff = 0.8, 
                    vInheritance = 0.5, propInherite = 0.8, rateDiff = 2, minRate = 1, 
                    propHetero = 0.4, minReads = 2, maxPercReads = 99.9, meanCov = 10, 
                    nbCores = 1, vSeed = -1),
        error=conditionMessage)
    
    exp <- "nbBlock must be a positive integer or numeric"
    
    message <- paste0("test.validateRunSimNumberParameters_nbBlock_string() ",
                      "- String as nbBlock parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimNumberParameters_nbBlock_zero <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimNumberParameters(nbSynCHR = 1, 
                        nbSimulation = 10, nbBlock = 0, nbCpG  = 2, vNbSample = 2, 
                        nbGeneration = 3, vpDiff = 0.22, vpDiffsd = 0.1, vDiff = 0.8, 
                        vInheritance = 0.5, propInherite = 0.8, rateDiff = 2, minRate = 1, 
                        propHetero = 0.4, minReads = 2, maxPercReads = 99.9, meanCov = 10,
                        nbCores = 1, vSeed = -1),
                    error=conditionMessage)
    
    exp <- "nbBlock must be a positive integer or numeric"
    
    message <- paste0("test.validateRunSimNumberParameters_nbBlock_zero() ",
                      "- Zero as nbBlock parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimNumberParameters_nbCpG_number_vector <- function() {
    obs <- tryCatch(
        methylInheritanceSim:::validateRunSimNumberParameters(nbSynCHR = 1, 
                        nbSimulation = 10, nbBlock = 2, nbCpG = c(3,2), vNbSample = 2, 
                        nbGeneration = 3, vpDiff = 0.2, vpDiffsd = 0.1, vDiff = 0.8, 
                        vInheritance = 0.5, propInherite = 0.8, rateDiff = 2, minRate = 1, 
                        propHetero = 0.4, minReads = 2, maxPercReads = 99.9, meanCov = 10, 
                        nbCores = 1, vSeed = -1),
        error=conditionMessage)
    
    exp <- "nbCpG must be a positive integer or numeric"
    
    message <- paste0("test.validateRunSimNumberParameters_nbCpG_number_vector() ",
                      "- Number vector as nbCpG parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimNumberParameters_nbCpG_string <- function() {
    obs <- tryCatch(
        methylInheritanceSim:::validateRunSimNumberParameters(nbSynCHR = 1, 
                nbSimulation = 10, nbBlock = 2, nbCpG  = "hi", vNbSample = 2, 
                nbGeneration = 3, vpDiff = 0.2, vpDiffsd = 0.1, vDiff = 0.8, 
                vInheritance = 0.5, propInherite = 0.8, rateDiff = 2, minRate = 1, 
                propHetero = 0.4, minReads = 2, maxPercReads = 99.9, 
                meanCov = 10, nbCores = 1, vSeed = -1),
        error=conditionMessage)
    
    exp <- "nbCpG must be a positive integer or numeric"
    
    message <- paste0("test.validateRunSimNumberParameters_nbCpG_string() ",
                      "- String as nbCpG parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimNumberParameters_nbCpG_zero <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimNumberParameters(nbSynCHR = 1, 
                nbSimulation = 10, nbBlock = 2, nbCpG  = 0, vNbSample = 2, 
                nbGeneration = 3, vpDiff = 0.22, vpDiffsd = 0.1, vDiff = 0.8, 
                vInheritance = 0.5, propInherite = 0.8, rateDiff = 2, minRate = 1, 
                propHetero = 0.4, minReads = 2, maxPercReads = 99.9, meanCov = 10,
                nbCores = 1, vSeed = -1),
            error=conditionMessage)
    
    exp <- "nbCpG must be a positive integer or numeric"
    
    message <- paste0("test.validateRunSimParameters_nbCpG_zero() ",
                      "- Zero as nbCpG parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimNumberParameters_nbGeneration_number_vector <- function() {
    obs <- tryCatch(
        methylInheritanceSim:::validateRunSimNumberParameters(nbSynCHR = 1, 
            nbSimulation = 10, nbBlock = 2, nbCpG  = 2, nbGeneration = c(3,2), 
            vNbSample = 2, vpDiff = 0.2, vpDiffsd = 0.1, vDiff = 0.8, 
            vInheritance = 0.5, propInherite = 0.8, rateDiff = 2, minRate = 1, 
            propHetero = 0.4, minReads = 2, maxPercReads = 99.9,  meanCov = 10, 
            nbCores = 1, vSeed = -1),
        error=conditionMessage)
    
    exp <- "nbGeneration must be a positive integer or numeric"
    
    message <- paste0("test.validateRunSimNumberParameters_nbGeneration_number_vector() ",
                      "- Number vector as nbGeneration parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimNumberParameters_nbGeneration_string <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimNumberParameters(nbSynCHR = 1, 
            nbSimulation = 10, nbBlock = 2, nbCpG  = 2, nbGeneration = "Hi",
            vNbSample = 2, vpDiff = 0.2, vpDiffsd = 0.1, vDiff = 0.8, 
            vInheritance = 0.5, propInherite = 0.8, rateDiff = 2, minRate = 1, 
            propHetero = 0.4, minReads = 2, maxPercReads = 99.9,  meanCov = 10, 
            nbCores = 1, vSeed = -1),
        error=conditionMessage)
    
    exp <- "nbGeneration must be a positive integer or numeric"
    
    message <- paste0("test.validateRunSimNumberParameters_nbGeneration_string() ",
                      "- String as nbGeneration parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimNumberParameters_nbGeneration_zero <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimNumberParameters(nbSynCHR = 1, 
                nbSimulation = 10, nbBlock = 2, nbCpG  = 2, nbGeneration = 0,
                vNbSample = 2, vpDiff = 0.22, vpDiffsd = 0.1, vDiff = 0.8, 
                vInheritance = 0.5, propInherite = 0.8, rateDiff = 2, minRate = 1, 
                propHetero = 0.4, minReads = 2, maxPercReads = 99.9, meanCov = 10,
                nbCores = 1, vSeed = -1),
            error=conditionMessage)
    
    exp <- "nbGeneration must be a positive integer or numeric"
    
    message <- paste0("test.validateRunSimNumberParameters_nbGeneration_zero() ",
                      "- Zero as nbGeneration parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimNumberParameters_vNbSample_string <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimNumberParameters(nbSynCHR = 1, 
                nbSimulation = 10, nbBlock = 2, nbCpG  = 2, nbGeneration = 3,
                vNbSample = "Hi", vpDiff = 0.22, vpDiffsd = 0.1, vDiff = 0.8, 
                vInheritance = 0.5, propInherite = 0.8, rateDiff = 2, minRate = 1, 
                propHetero = 0.4, minReads = 2, maxPercReads = 99.9, meanCov = 10,
                nbCores = 1, vSeed = -1),
            error=conditionMessage)
    
    exp <- "vNbSample must be a vector of distinct positive integer"
    
    message <- paste0("test.validateRunSimNumberParameters_vNbSample_string() ",
                      "- String as vNbSample parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimNumberParameters_vNbSample_duplicate <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimNumberParameters(nbSynCHR = 1, 
                nbSimulation = 10, nbBlock = 2, nbCpG  = 2, nbGeneration = 3,
                vNbSample = c(1, 2, 2), vpDiff = 0.22, vpDiffsd = 0.1, vDiff = 0.8, 
                vInheritance = 0.5, propInherite = 0.8,  rateDiff = 2, minRate = 1, 
                propHetero = 0.4, minReads = 2, maxPercReads = 99.9, meanCov = 10,
                nbCores = 1, vSeed = -1),
            error=conditionMessage)
    
    exp <- "vNbSample must be a vector of distinct positive integer"
    
    message <- paste0("test.validateRunSimNumberParameters_vNbSample_duplicate() ",
                      "- Duplicate elements in vNbSample parameter did not 
                      generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimNumberParameters_vNbSample_float <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimNumberParameters(nbSynCHR = 1, 
                nbSimulation = 10, nbBlock = 2, nbCpG  = 2, nbGeneration = 3,
                vNbSample = c(1, 2, 2.1), vpDiff = 0.22, vpDiffsd = 0.1, 
                vDiff = 0.8, vInheritance = 0.5, propInherite = 0.8, rateDiff = 2, 
                minRate = 1, propHetero = 0.4, minReads = 2, maxPercReads = 99.9, 
                meanCov = 10, nbCores = 1, vSeed = -1),
            error=conditionMessage)
    
    exp <- "vNbSample must be a vector of distinct positive integer"
    
    message <- paste0("test.validateRunSimNumberParameters_vNbSample_float() ",
                      "- Float as vNbSample parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}


test.validateRunSimNumberParameters_vpDiff_string <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimNumberParameters(nbSynCHR = 1, 
                nbSimulation = 10, nbBlock = 2, nbCpG  = 2, nbGeneration = 3,
                vNbSample = 3, vpDiff = "Hi", vpDiffsd = 0.1, vDiff = 0.8, 
                vInheritance = 0.5, propInherite = 0.8, rateDiff = 2, minRate = 1, 
                propHetero = 0.4, minReads = 2, maxPercReads = 99.9, meanCov = 10,
                nbCores = 1, vSeed = -1),
        error=conditionMessage)
    
    exp <- "vpDiff must be a vector of distinct positive double include in (0,1]"
    
    message <- paste0("test.validateRunSimNumberParameters_vpDiff_string() ",
                    "- String as vpDiff parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimNumberParameters_vpDiff_duplicate <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimNumberParameters(nbSynCHR = 1, 
                nbSimulation = 10, nbBlock = 2, nbCpG  = 2, nbGeneration = 3,
                vNbSample = 3, vpDiff = c(0.5, 0.5, 0.8), vpDiffsd = 0.1, 
                vDiff = 0.9, vInheritance = 0.5, propInherite = 0.8, rateDiff = 2, 
                minRate = 1, propHetero = 0.4, minReads = 2, maxPercReads = 99.9, 
                meanCov = 10, nbCores = 1, vSeed = -1),
            error=conditionMessage)
    
    exp <- "vpDiff must be a vector of distinct positive double include in (0,1]"
    
    message <- paste0("test.validateRunSimNumberParameters_vpDiff_duplicate() ",
                      "- Duplicate elements in vpDiff parameter did not 
                      generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimNumberParameters_vpDiff_zero <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimNumberParameters(nbSynCHR = 1, 
                nbSimulation = 10, nbBlock = 2, nbCpG  = 2, nbGeneration = 3,
                vNbSample = 3, vpDiff = c(0.5, 0, 0.8), vpDiffsd = 0.1, vDiff = 0.8, 
                vInheritance = 0.5,propInherite = 0.8, rateDiff = 2, minRate = 1, 
                propHetero = 0.4, minReads = 2, maxPercReads = 99.9, meanCov = 10,
                nbCores = 1, vSeed = -1),
            error=conditionMessage)
    
    exp <- "vpDiff must be a vector of distinct positive double include in (0,1]"
    
    message <- paste0("test.validateRunSimNumberParameters_vpDiff_zero() ",
                      "- Zero in vpDiff parameter did not 
                      generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimNumberParameters_vpDiff_two <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimNumberParameters(nbSynCHR = 1, 
                nbSimulation = 10, nbBlock = 2, nbCpG  = 2, nbGeneration = 3,
                vNbSample = 3, vpDiff = c(0.5, 2, 0.8), vpDiffsd = 0.1, vDiff = 0.8, 
                vInheritance = 0.5, propInherite = 0.8, rateDiff = 2, minRate = 1, 
                propHetero = 0.4, minReads = 2, maxPercReads = 99.9, meanCov = 10,
                nbCores = 1, vSeed = -1),
            error=conditionMessage)
    
    exp <- "vpDiff must be a vector of distinct positive double include in (0,1]"
    
    message <- paste0("test.validateRunSimNumberParameters_vpDiff_two() ",
                      "- Two in vpDiff parameter did not 
                      generated expected results.")
    
    checkEquals(obs, exp, message)
}


test.validateRunSimNumberParameters_vpDiffsd_string <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimNumberParameters(nbSynCHR = 1, 
                nbSimulation = 10, nbBlock = 2, nbCpG  = 2, nbGeneration = 3,
                vNbSample = 3, vpDiff = 0.9, vpDiffsd = "Hi", vDiff = 0.8, 
                vInheritance = 0.5, propInherite = 0.8, rateDiff = 2, 
                minRate = 1, propHetero = 0.4, minReads = 2, maxPercReads = 99.9, 
                meanCov = 10, nbCores = 1, vSeed = -1),
            error=conditionMessage)
    
    exp <- "vpDiffsd must be a vector of non-negative double"
    
    message <- paste0("test.validateRunSimNumberParameters_vpDiffsd_string() ",
                      "- String as vpDiffsd parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimNumberParameters_vpDiffsd_negative <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimNumberParameters(nbSynCHR = 1, 
                nbSimulation = 10, nbBlock = 2, nbCpG  = 2, nbGeneration = 3,
                vNbSample = 3, vpDiff = 0.9, vpDiffsd = -1, vDiff = 0.8, 
                vInheritance = 0.5, propInherite = 0.8, rateDiff = 2, minRate = 1, 
                propHetero = 0.4, minReads = 2, maxPercReads = 99.9, meanCov = 10,
                nbCores = 1, vSeed = -1),
            error=conditionMessage)
    
    exp <- "vpDiffsd must be a vector of non-negative double"
    
    message <- paste0("test.validateRunSimNumberParameters_vpDiffsd_negative() ",
                      "- Negative as vpDiffsd parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimNumberParameters_vpDiffsd_vpDiff <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimNumberParameters(nbSynCHR = 1, 
                nbSimulation = 10, nbBlock = 2, nbCpG  = 2, nbGeneration = 3,
                vNbSample = 3, vpDiff = 0.9, vpDiffsd = c(0.1, 0.1), vDiff = 0.8, 
                vInheritance = 0.5, propInherite = 0.8, rateDiff = 2, minRate = 1, 
                propHetero = 0.4, minReads = 2, maxPercReads = 99.9, meanCov = 10,
                nbCores = 1, vSeed = -1),
            error=conditionMessage)
    
    exp <- "vpDiff and vpDiffsd must be the same length"
    
    message <- paste0("test.validateRunSimNumberParameters_vpDiffsd_vpDiff() ",
                      "- vpDiffsd vpDiff not same length did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimNumberParameters_vDiff_string <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimNumberParameters(nbSynCHR = 1, 
                nbSimulation = 10,nbBlock = 2, nbCpG  = 2, nbGeneration = 3,
                vNbSample = 3, vpDiff = 0.9, vpDiffsd = 0.1,  vDiff = "Hi", 
                vInheritance = 0.5, propInherite = 0.8, rateDiff = 2, minRate = 1, 
                propHetero = 0.4, minReads = 2, maxPercReads = 99.9, meanCov = 10,
                nbCores = 1, vSeed = -1),
            error=conditionMessage)
    
    exp <- "vDiff must be a vector of distinct non-negative double include in [0,1]"
    
    message <- paste0("test.validateRunSimNumberParameters_vDiff_string() ",
                      "- String as vDiff parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimNumberParameters_vDiff_duplicate <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimNumberParameters(nbSynCHR = 1, 
            nbSimulation = 10, nbBlock = 2, nbCpG  = 2, nbGeneration = 3,
            vNbSample = 3, vpDiff = c(0.5), vpDiffsd = 0.1, vDiff = c(0.5, 0.5, 0.8), 
            vInheritance = 0.5, propInherite = 0.8, rateDiff = 2, minRate = 1, 
            propHetero = 0.4, minReads = 2, maxPercReads = 99.9, meanCov = 10,
            nbCores = 1, vSeed = -1),
        error=conditionMessage)
    
    exp <- "vDiff must be a vector of distinct non-negative double include in [0,1]"
    
    message <- paste0("test.validateRunSimNumberParameters_vDiff_duplicate() ",
                    "- Duplicate elements in vDiff parameter did not ", 
                    "generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimNumberParameters_vDiff_negative <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimNumberParameters(nbSynCHR = 1, 
                nbSimulation = 10, nbBlock = 2, nbCpG  = 2, nbGeneration = 3,
                vNbSample = 3, vpDiff = c(0.5), vpDiffsd = 0.1, vDiff = -1, 
                vInheritance = 0.5, propInherite = 0.8, rateDiff = 2, minRate = 1, 
                propHetero = 0.4, minReads = 2, maxPercReads = 99.9, meanCov = 10,
                nbCores = 1, vSeed = -1),
            error=conditionMessage)
    
    exp <- "vDiff must be a vector of distinct non-negative double include in [0,1]"
    
    message <- paste0("test.validateRunSimNumberParameters_vDiff_negative() ",
                      "- Negative in vDiff parameter did not ",
                      "generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimNumberParameters_vDiff_two <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimNumberParameters(nbSynCHR = 1, 
                    nbSimulation = 10, nbBlock = 2, nbCpG  = 2, nbGeneration = 3,
                    vNbSample = 3, vpDiff = c(0.5), vpDiffsd = 0.1, vDiff = 2, 
                    vInheritance = 0.5, propInherite = 0.8, rateDiff = 2, 
                    minRate = 1, propHetero = 0.4, minReads = 2, maxPercReads = 99.9, 
                    meanCov = 10, nbCores = 1, vSeed = -1),
            error=conditionMessage)
    
    exp <- "vDiff must be a vector of distinct non-negative double include in [0,1]"
    
    message <- paste0("test.validateRunSimNumberParameters_vDiff_two() ",
                      "- Two in vDiff parameter did not 
                      generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimNumberParameters_vInheritance_string <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimNumberParameters(nbSynCHR = 1, 
                    nbSimulation = 10, nbBlock = 2, nbCpG  = 2, nbGeneration = 3,
                    vNbSample = 3, vpDiff = 0.9, vpDiffsd = 0.1, vDiff = 0.8, 
                    vInheritance = "Hi", propInherite = 0.8, rateDiff = 2, 
                    minRate = 1, propHetero = 0.4, minReads = 2, maxPercReads = 99.9, 
                    meanCov = 10, nbCores = 1, vSeed = -1),
            error=conditionMessage)
    
    exp <- paste0("vInheritance must be a vector of distinct non-negative ", 
                  "double include in [0,1]")
    
    message <- paste0("test.validateRunSimNumberParameters_vInheritance_string() ",
                      "- String as vInheritance parameter did not generated ", 
                      "expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimNumberParameters_vInheritance_duplicate <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimNumberParameters(nbSynCHR = 1, 
                        nbSimulation = 10, nbBlock = 2, nbCpG  = 2, nbGeneration = 3,
                        vNbSample = 3, vpDiff = c(0.5), vpDiffsd = 0.1, 
                        vDiff = c(0.5, 0.7, 0.8), vInheritance = c(0.5, 0.5, 0.7),
                        propInherite = 0.8, rateDiff = 2, minRate = 1, propHetero = 0.4, 
                        minReads = 2, maxPercReads = 99.9, meanCov = 10,
                        nbCores = 1, vSeed = -1),
            error=conditionMessage)
    
    exp <- paste0("vInheritance must be a vector of distinct non-negative ", 
                  "double include in [0,1]")
    
    message <- paste0("test.validateRunSimNumberParameters_vInheritance_duplicate() ",
                      "- Duplicate elements in vInheritance parameter did not 
                      generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimNumberParameters_vInheritance_negative <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimNumberParameters(nbSynCHR = 1, 
                nbSimulation = 10, nbBlock = 2, nbCpG  = 2, nbGeneration = 3,
                vNbSample = 3, vpDiff = c(0.5), vpDiffsd = 0.1, vDiff = 0.8, 
                vInheritance = -1, propInherite = 0.8, rateDiff = 2, minRate = 1, 
                propHetero = 0.4, minReads = 2, maxPercReads = 99.9, meanCov = 10,
                nbCores = 1, vSeed = -1),
            error=conditionMessage)
    
    exp <- paste0("vInheritance must be a vector of distinct non-negative ", 
                  "double include in [0,1]")
    
    message <- paste0("test.validateRunSimNumberParameters_vInheritance_negative() ",
                      "- Negative in vInheritance parameter did not ",
                      "generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimNumberParameters_vInheritance_two <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimNumberParameters(nbSynCHR = 1, 
                        nbSimulation = 10,  nbBlock = 2, nbCpG  = 2, nbGeneration = 3,
                        vNbSample = 3, vpDiff = c(0.5), vpDiffsd = 0.1, vDiff = 0.8, 
                        vInheritance = 2, propInherite = 0.8, rateDiff = 2, 
                        minRate = 1, propHetero = 0.4, minReads = 2, maxPercReads = 99.9, 
                        meanCov = 10, nbCores = 1, vSeed = -1),
            error=conditionMessage)
    
    exp <- paste0("vInheritance must be a vector of distinct non-negative ", 
                  "double include in [0,1]")
    
    message <- paste0("test.validateRunSimNumberParameters_vInheritance_two() ",
                      "- Two in vInheritance parameter did not 
                      generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimNumberParameters_nbCores_string <- function() {
    obs <- tryCatch(
        methylInheritanceSim:::validateRunSimNumberParameters(nbSynCHR = 1, 
                nbBlock = 2, nbCpG  = 2, vNbSample = 2, nbGeneration = 3, 
                vpDiff = 0.2, vpDiffsd = 0.1, vDiff = 0.2, vInheritance = 0.2,
                propInherite = 0.8, rateDiff = 0.8, minRate = 0.1, propHetero = 0.4, 
                minReads = 2, maxPercReads = 99.9, meanCov = 10, n = 3, 
                nbCores = "hi", vSeed = -1),
        error=conditionMessage)
    
    exp <- "nbCores must be a positive integer or numeric"
    
    message <- paste0("test.validateRunSimNumberParameters_nbCores_string() ",
                      "- String as nbCores parameter did not generated expected
                      results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimNumberParameters_nbCores_zero <- function() {
    obs <- tryCatch(
        methylInheritanceSim:::validateRunSimNumberParameters(nbSynCHR = 1, 
                nbBlock = 2, nbCpG = 2, vNbSample = 2, nbGeneration = 3, 
                vpDiff = 0.2, vpDiffsd = 0.1,  vDiff = 0.2, vInheritance = 0.2,
                propInherite = 0.8, rateDiff = 0.8, minRate = 0.1, propHetero = 0.4, 
                minReads = 2, maxPercReads = 99.9, meanCov = 10, n = 3, 
                nbCores = 0, vSeed = -1),
        error=conditionMessage)
    
    exp <- "nbCores must be a positive integer or numeric"
    
    message <- paste0("test.validateRunSimNumberParameters_nbCores_zero() ",
                      "- Zero as nbCores parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimNumberParameters_minrate_1 <- function() {
    obs <- tryCatch(
        methylInheritanceSim:::validateRunSimNumberParameters(nbSynCHR = 1, 
                    nbBlock = 2, nbCpG = 2, vNbSample = 2, nbGeneration = 3, 
                    vpDiff = 0.2, vpDiffsd = 0.1, vDiff = 0.2, vInheritance = 0.2,
                    propInherite = 0.8, rateDiff = 0.8, minRate = 1, propHetero = 0.4, 
                    minReads = 2, maxPercReads = 99.9, meanCov = 10, n = 3, 
                    nbCores = 1, vSeed = -1),
        error=conditionMessage)
    
    exp <- "minRate must be a non-negative double include in [0,1)"
    
    message <- paste0("test.validateRunSimNumberParameters_minrate_sup_to_1() ",
                      "- 1 as minRate parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}


test.validateRunSimNumberParameters_minrate_vector <- function() {
    obs <- tryCatch(
        methylInheritanceSim:::validateRunSimNumberParameters(nbSynCHR = 1, 
                nbBlock = 2, nbCpG = 2, vNbSample = 2, nbGeneration = 3, 
                vpDiff = 0.2, vpDiffsd = 0.1, vDiff = 0.2, vInheritance = 0.2,
                propInherite = 0.8, rateDiff = 0.8, minRate = c(0.01, 0.002), 
                propHetero = 0.4, minReads = 2, maxPercReads = 99.9, 
                meanCov = 10, n = 3, nbCores = 1, vSeed = -1),
        error=conditionMessage)
    
    exp <- "minRate must be a non-negative double include in [0,1)"
    
    message <- paste0("test.validateRunSimNumberParameters_minrate_vector() ",
                      "- Vector as minRate parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimNumberParameters_propHetero_vector <- function() {
    obs <- tryCatch(
        methylInheritanceSim:::validateRunSimNumberParameters(nbSynCHR = 1, 
            nbBlock = 2, nbCpG = 2, vNbSample = 2, nbGeneration = 3, 
            vpDiff = 0.2, vpDiffsd = 0.1, vDiff = 0.2, vInheritance = 0.2,
            propInherite = 0.8, rateDiff = 0.8, minRate = 0.1, propHetero = c(0.4, 0.5), 
            minReads = 2, maxPercReads = 99.9, meanCov = 10, n = 3, 
            nbCores = 1, vSeed = -1),
        error=conditionMessage)
    
    exp <- "propHetero must be a non-negative double include in [0,1]"
    
    message <- paste0("test.validateRunSimNumberParameters_propHetero_vector() ",
                      "- Vector as propHetero parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimNumberParameters_propHetero_sup_to_1 <- function() {
    obs <- tryCatch(
        methylInheritanceSim:::validateRunSimNumberParameters(nbSynCHR = 1, 
                nbBlock = 2, nbCpG = 2, vNbSample = 2, nbGeneration = 3, 
                vpDiff = 0.2, vpDiffsd = 0.1, vDiff = 0.2, vInheritance = 0.2,
                propInherite = 0.8, rateDiff = 0.8, minRate = 0.1, 
                propHetero = 1.001, minReads = 2, maxPercReads = 99.9, 
                meanCov = 10, n = 3, nbCores = 1, vSeed = -1),
        error=conditionMessage)
    
    exp <- "propHetero must be a non-negative double include in [0,1]"
    
    message <- paste0("test.validateRunSimNumberParameters_propHetero_sup_to_1() ",
                      "- Superior to 1 as propHetero parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimNumberParameters_propHetero_inf_zero <- function() {
    obs <- tryCatch(
        methylInheritanceSim:::validateRunSimNumberParameters(nbSynCHR = 1, 
                    nbBlock = 2, nbCpG = 2, vNbSample = 2, nbGeneration = 3, 
                    vpDiff = 0.2, vpDiffsd = 0.1, vDiff = 0.2, vInheritance = 0.2,
                    propInherite = 0.8, rateDiff = 0.8, minRate = 0.1, 
                    propHetero = -0.001, minReads = 2, maxPercReads = 99.9, 
                    meanCov = 10, n = 3, nbCores = 1, vSeed = -1),
        error=conditionMessage)
    
    exp <- "propHetero must be a non-negative double include in [0,1]"
    
    message <- paste0("test.validateRunSimNumberParameters_propHetero_inf_zero() ",
                      "- Inferior to zero as propHetero parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}


test.validateRunSimNumberParameters_propInherite_inf_zero <- function() {
    obs <- tryCatch(
        methylInheritanceSim:::validateRunSimNumberParameters(nbSynCHR = 1, 
                nbBlock = 2, nbCpG = 2, vNbSample = 2, nbGeneration = 3, 
                vpDiff = 0.2, vpDiffsd = 0.1, vDiff = 0.2, vInheritance = 0.2,
                propInherite = -0.001, rateDiff = 0.8, minRate = 0.1, 
                propHetero = 0.03, minReads = 2, maxPercReads = 99.9, 
                meanCov = 10, n = 3, nbCores = 1, vSeed = -1),
        error=conditionMessage)
    
    exp <- "propInherite must be a non-negative double include in [0,1]"
    
    message <- paste0("test.validateRunSimNumberParameters_propInherite_inf_zero() ",
                      "- Inferior to zero as propInherite parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimNumberParameters_propInherite_sup_to_1 <- function() {
    obs <- tryCatch(
        methylInheritanceSim:::validateRunSimNumberParameters(nbSynCHR = 1, 
            nbBlock = 2, nbCpG = 2, vNbSample = 2, nbGeneration = 3, 
            vpDiff = 0.2, vpDiffsd = 0.1, vDiff = 0.2, vInheritance = 0.2,
            propInherite = 1.001, rateDiff = 0.8, minRate = 0.1, propHetero = 0.3, 
            minReads = 2, maxPercReads = 99.9, meanCov = 10, n = 3, 
            nbCores = 1, vSeed = -1),
        error=conditionMessage)
    
    exp <- "propInherite must be a non-negative double include in [0,1]"
    
    message <- paste0("test.validateRunSimNumberParameters_propInherite_sup_to_1() ",
                      "- Superior to 1 as propInherite parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimNumberParameters_propInherite_vector <- function() {
    obs <- tryCatch(
        methylInheritanceSim:::validateRunSimNumberParameters(nbSynCHR = 1, 
                nbBlock = 2, nbCpG = 2, vNbSample = 2, nbGeneration = 3, 
                vpDiff = 0.2, vpDiffsd = 0.1, vDiff = 0.2, vInheritance = 0.2,
                propInherite = c(0.4, 0.5), rateDiff = 0.8, minRate = 0.1, 
                propHetero = 0.3, minReads = 2, maxPercReads = 99.9, 
                meanCov = 10, n = 3, nbCores = 1, vSeed = -1),
        error=conditionMessage)
    
    exp <- "propInherite must be a non-negative double include in [0,1]"
    
    message <- paste0("test.validateRunSimNumberParameters_propInherite_vector() ",
                      "- Vector as propInherite parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimNumberParameters_minReads_zero <- function() {
    obs <- tryCatch(
        methylInheritanceSim:::validateRunSimNumberParameters(nbSynCHR = 1, 
                nbBlock = 2, nbCpG = 2, vNbSample = 2, nbGeneration = 3, 
                vpDiff = 0.2, vpDiffsd = 0.1, vDiff = 0.2, vInheritance = 0.2,
                propInherite = 0.6, rateDiff = 0.8, minRate = 0.1, 
                propHetero = 0.3, minReads = 0, maxPercReads = 99.9, 
                meanCov = 10, n = 3, nbCores = 1, vSeed = -1),
        error=conditionMessage)
    
    exp <- "minReads must be a positive integer or numeric"
    
    message <- paste0("test.validateRunSimNumberParameters_minReads_zero() ",
                      "- Zero as minReads parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimNumberParameters_minReads_vector <- function() {
    obs <- tryCatch(
        methylInheritanceSim:::validateRunSimNumberParameters(nbSynCHR = 1, 
                nbBlock = 2, nbCpG = 2, vNbSample = 2, nbGeneration = 3, 
                vpDiff = 0.2, vpDiffsd = 0.1, vDiff = 0.2, vInheritance = 0.2,
                propInherite = 0.6, rateDiff = 0.8, minRate = 0.1, 
                propHetero = 0.3, minReads = c(1, 3), maxPercReads = 99.9, 
                meanCov = 10, n = 3, nbCores = 1, vSeed = -1),
        error=conditionMessage)
    
    exp <- "minReads must be a positive integer or numeric"
    
    message <- paste0("test.validateRunSimNumberParameters_minReads_vector() ",
                      "- Vector as minReads parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimNumberParameters_maxPercReads_inf_to_zero <- function() {
    obs <- tryCatch(
        methylInheritanceSim:::validateRunSimNumberParameters(nbSynCHR = 1, 
                    nbBlock = 2, nbCpG = 2, vNbSample = 2, nbGeneration = 3, 
                    vpDiff = 0.2, vpDiffsd = 0.1, vDiff = 0.2, vInheritance = 0.2,
                    propInherite = 0.6,  rateDiff = 0.8, minRate = 0.1, 
                    propHetero = 0.3, minReads = 4, maxPercReads = -0.01, 
                    meanCov = 10, n = 3, nbCores = 1, vSeed = -1),
        error=conditionMessage)
    
    exp <- "maxPercReads must be a positive double between [0,100]"
    
    message <- paste0("test.validateRunSimNumberParameters_maxPercReads_inf_to_zero() ",
                      "- Inferior to zero as maxPercReads parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimNumberParameters_maxPercReads_sup_to_100 <- function() {
    obs <- tryCatch(
        methylInheritanceSim:::validateRunSimNumberParameters(nbSynCHR = 1, 
                nbBlock = 2, nbCpG = 2, vNbSample = 2, nbGeneration = 3, 
                vpDiff = 0.2, vpDiffsd = 0.1, vDiff = 0.2, vInheritance = 0.2,
                propInherite = 0.6, rateDiff = 0.8, minRate = 0.1, 
                propHetero = 0.3, minReads = 4, maxPercReads = 100.001, 
                meanCov = 10, n = 3, nbCores = 1, vSeed = -1),
        error=conditionMessage)
    
    exp <- "maxPercReads must be a positive double between [0,100]"
    
    message <- paste0("test.validateRunSimNumberParameters_maxPercReads_sup_to_100() ",
                      "- Superior to 100 as maxPercReads parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimNumberParameters_maxPercReads_vector <- function() {
    obs <- tryCatch(
        methylInheritanceSim:::validateRunSimNumberParameters(nbSynCHR = 1, 
                nbBlock = 2, nbCpG = 2, vNbSample = 2, nbGeneration = 3, 
                vpDiff = 0.2, vpDiffsd = 0.1, vDiff = 0.2, vInheritance = 0.2,
                propInherite = 0.6, rateDiff = 0.8, minRate = 0.1, propHetero = 0.3, 
                minReads = 4, maxPercReads = c(99.8, 99.6), meanCov = 10, n = 3, 
                nbCores = 1, vSeed = -1),
        error=conditionMessage)
    
    exp <- "maxPercReads must be a positive double between [0,100]"
    
    message <- paste0("test.validateRunSimNumberParameters_maxPercReads_vector() ",
                      "- Vector as maxPercReads parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimNumberParameters_rateDiff_1 <- function() {
    obs <- tryCatch(
        methylInheritanceSim:::validateRunSimNumberParameters(nbSynCHR = 1, 
                nbBlock = 2, nbCpG = 2, vNbSample = 2, nbGeneration = 3, 
                vpDiff = 0.2, vpDiffsd = 0.1, vDiff = 0.2, vInheritance = 0.2,
                propInherite = 0.6, rateDiff = 1, minRate = 0.1, propHetero = 0.3, 
                minReads = 4, maxPercReads = 99.9, meanCov = 10, n = 3, 
                nbCores = 1, vSeed = -1),
        error=conditionMessage)
    
    exp <- "rateDiff must be a positive double include in (0,1)"
    
    message <- paste0("test.validateRunSimNumberParameters_rateDiff_1() ",
                      "- 1 as rateDiff parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimNumberParameters_rateDiff_zero <- function() {
    obs <- tryCatch(
        methylInheritanceSim:::validateRunSimNumberParameters(nbSynCHR = 1, 
                nbBlock = 2, nbCpG = 2, vNbSample = 2, nbGeneration = 3, 
                vpDiff = 0.2, vpDiffsd = 0.1, vDiff = 0.2, vInheritance = 0.2,
                propInherite = 0.6, rateDiff = 0, minRate = 0.1, propHetero = 0.3, 
                minReads = 4, maxPercReads = 99.9, meanCov = 10, n = 3, 
                nbCores = 1, vSeed = -1),
        error=conditionMessage)
    
    exp <- "rateDiff must be a positive double include in (0,1)"
    
    message <- paste0("test.validateRunSimNumberParameters_rateDiff_zero() ",
                      "- Zero as rateDiff parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimNumberParameters_rateDiff_vector <- function() {
    obs <- tryCatch(
        methylInheritanceSim:::validateRunSimNumberParameters(nbSynCHR = 1, 
            nbBlock = 2, nbCpG = 2, vNbSample = 2, nbGeneration = 3, 
            vpDiff = 0.2, vpDiffsd = 0.1, vDiff = 0.2, vInheritance = 0.2,
            propInherite = 0.6, rateDiff = c(0.02, 0.3), minRate = 0.1, 
            propHetero = 0.3, minReads = 4, maxPercReads = 99.9, 
            meanCov = 10, n = 3, nbCores = 1, vSeed = -1),
        error=conditionMessage)
    
    exp <- "rateDiff must be a positive double include in (0,1)"
    
    message <- paste0("test.validateRunSimNumberParameters_rateDiff_vector() ",
                      "- Vector as rateDiff parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimNumberParameters_meanCov_string <- function() {
    obs <- tryCatch(
        methylInheritanceSim:::validateRunSimNumberParameters(nbSynCHR = 1, 
            nbBlock = 2, nbCpG = 2, vNbSample = 2, nbGeneration = 3, 
            vpDiff = 0.2, vpDiffsd = 0.1, vDiff = 0.2, vInheritance = 0.2,
            propInherite = 0.6, rateDiff = 0.4, minRate = 0.1, propHetero = 0.3, 
            minReads = 4, maxPercReads = 99.9, meanCov = "allo", n = 3, 
            nbCores = 1, vSeed = -1),
        error=conditionMessage)
    
    exp <- "meanCov must be a positive integer or numeric"
    
    message <- paste0("test.validateRunSimNumberParameters_meanCov_string() ",
                      "- String as meanCov parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimNumberParameters_meanCov_zero <- function() {
    obs <- tryCatch(
        methylInheritanceSim:::validateRunSimNumberParameters(nbSynCHR = 1, 
            nbBlock = 2, nbCpG = 2, vNbSample = 2, nbGeneration = 3, 
            vpDiff = 0.2, vpDiffsd = 0.1, vDiff = 0.2, vInheritance = 0.2,
            propInherite = 0.6, rateDiff = 0.4, minRate = 0.1, 
            propHetero = 0.3, minReads = 4, maxPercReads = 99.9, 
            meanCov = 0, n = 3, nbCores = 1, vSeed = -1),
        error=conditionMessage)
    
    exp <- "meanCov must be a positive integer or numeric"
    
    message <- paste0("test.validateRunSimNumberParameters_meanCov_zero() ",
                      "- Zero as meanCov parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimNumberParameters_vSeed_string <- function() {
    obs <- tryCatch(
        methylInheritanceSim:::validateRunSimNumberParameters(nbSynCHR = 1, 
            nbBlock = 2, nbCpG = 2, vNbSample = 2, nbGeneration = 3, 
            vpDiff = 0.2, vpDiffsd = 0.1, vDiff = 0.2, vInheritance = 0.2,
            propInherite = 0.6,  rateDiff = 0.4, minRate = 0.1, propHetero = 0.3, 
            minReads = 4, maxPercReads = 99.9, meanCov = 10, n = 3, 
            nbCores = 1, vSeed = "test"),
        error=conditionMessage)
    
    exp <- "vSeed must be an integer or numeric"
    
    message <- paste0("test.validateRunSimNumberParameters_vSeed_string() ",
                      "- String as vSeed parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimNumberParameters_good_01 <- function() {
    obs <- methylInheritanceSim:::validateRunSimNumberParameters(nbSynCHR = 1, 
                nbBlock = 2, nbCpG = 2, vNbSample = 2, nbGeneration = 3, 
                vpDiff = 0.2, vpDiffsd = 0.1, vDiff = 0.2, vInheritance = 0.2,
                propInherite = 0.6, rateDiff = 0.4, minRate = 0.1, propHetero = 0.3, 
                minReads = 4, maxPercReads = 99.9, meanCov = 10, n = 3, 
                nbCores = 1, vSeed = -1)
    
    exp <- 0
    
    message <- paste0("test.validateRunSimNumberParameters_good_01() ",
                      "- All valid parameters did not generated expected results.")
    
    checkEquals(obs, exp, message)
}