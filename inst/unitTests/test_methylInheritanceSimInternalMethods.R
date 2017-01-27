###################################################
##
## Test the methylInheritanceSimInternalMethods functions
##
###################################################

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
##validateRunSimParameters() function
###################################################


###################################################
## validateRunSimParameters() function
###################################################

test.validateRunSimParameters_outputDir_number <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimParameters(outputDir = 33,
                                    fileID = "F1", 
                                    nbSynCHR = 1, 
                                    methData = samplesForChrSynthetic,
                                    nbSimulation = 10,
                                    nbBlock = 3, nbCpG  = 2,
                                    vNbSample = 2, 
                                    nbGeneration = 3, 
                                    vpDiff = 2, vpDiffsd = 0.1, 
                                    vDiff = 0.22, 
                                    vInheritance = 0.5,
                                    propInherite = 0.8, 
                                    rateDiff = 2, 
                                    minRate = 1, 
                                    propHetero = 0.4, 
                                    minReads = 2, 
                                    maxPercReads = 99.9, 
                                    context = "CpG", assembly = "hg19",
                                    meanCov = 10, 
                                    keepDiff = TRUE, saveGRanges = FALSE, 
                                    saveMethylKit = FALSE,
                                    runAnalysis = FALSE,
                                    nbCores = 1, vSeed = -1),
                            error=conditionMessage)

    exp <- "outputDir must be a character string or NULL"
    
    message <- paste0("test.validateRunSimParameters_outputDir_number() ",
                      "- Number as outputDir parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimParameters_fileID_number <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimParameters(outputDir = "test",
                                                                    fileID = 2, 
                                                                    nbSynCHR = 1, 
                                                                    methData = samplesForChrSynthetic,
                                                                    nbSimulation = 10,
                                                                    nbBlock = 3, nbCpG  = 2,
                                                                    vNbSample = 2, 
                                                                    nbGeneration = 3, 
                                                                    vpDiff = 0.2, vpDiffsd = 0.1, 
                                                                    vDiff = 0.8, 
                                                                    vInheritance = 0.5,
                                                                    propInherite = 0.8, 
                                                                    rateDiff = 2, 
                                                                    minRate = 1, 
                                                                    propHetero = 0.4, 
                                                                    minReads = 2, 
                                                                    maxPercReads = 99.9, 
                                                                    context = "CpG", assembly = "hg19",
                                                                    meanCov = 10, 
                                                                    keepDiff = TRUE, saveGRanges = FALSE, 
                                                                    saveMethylKit = FALSE,
                                                                    runAnalysis = FALSE,
                                                                    nbCores = 1, vSeed = -1),
                    error=conditionMessage)
    
    exp <- "fileID must be a character string or NULL"
    
    message <- paste0("test.validateRunSimParameters_fileID_number() ",
                      "- Number as fileID parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimParameters_nbSynCHR_not_number <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimParameters(outputDir = "test",
                                                        fileID = "F1", 
                                                        nbSynCHR = "hi", 
                                                        methData = samplesForChrSynthetic,
                                                        nbSimulation = 10,
                                                        nbBlock = 3, nbCpG  = 2,
                                                        vNbSample = 2, 
                                                        nbGeneration = 3, 
                                                        vpDiff = 0.2, vpDiffsd = 0.1, 
                                                        vDiff = 0.8, 
                                                        vInheritance = 0.5,
                                                        propInherite = 0.8, 
                                                        rateDiff = 2, 
                                                        minRate = 1, 
                                                        propHetero = 0.4, 
                                                        minReads = 2, 
                                                        maxPercReads = 99.9, 
                                                        context = "CpG", assembly = "hg19",
                                                        meanCov = 10, 
                                                        keepDiff = TRUE, saveGRanges = FALSE, 
                                                        saveMethylKit = FALSE,
                                                        runAnalysis = FALSE,
                                                        nbCores = 1, vSeed = -1),
                    error=conditionMessage)
    
    exp <- "nbSynCHR must be a positive integer or numeric"
    
    message <- paste0("test.validateRunSimParameters_nbSynCHR_not_number() ",
                      "- Not a number as nbSynCHR parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimParameters_nbSynCHR_vector_number <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimParameters(outputDir = "test",
                                                                    fileID = "F1", 
                                                                    nbSynCHR = c(1,2), 
                                                                    methData = samplesForChrSynthetic,
                                                                    nbSimulation = 10,
                                                                    nbBlock = 3, nbCpG  = 2,
                                                                    vNbSample = 2, 
                                                                    nbGeneration = 3, 
                                                                    vpDiff = 0.2, vpDiffsd = 0.1, 
                                                                    vDiff = 0.8, 
                                                                    vInheritance = 0.5,
                                                                    propInherite = 0.8, 
                                                                    rateDiff = 2, 
                                                                    minRate = 1, 
                                                                    propHetero = 0.4, 
                                                                    minReads = 2, 
                                                                    maxPercReads = 99.9, 
                                                                    context = "CpG", assembly = "hg19",
                                                                    meanCov = 10, 
                                                                    keepDiff = TRUE, saveGRanges = FALSE, 
                                                                    saveMethylKit = FALSE,
                                                                    runAnalysis = FALSE,
                                                                    nbCores = 1, vSeed = -1),
                    error=conditionMessage)
    
    exp <- "nbSynCHR must be a positive integer or numeric"
    
    message <- paste0("test.validateRunSimParameters_nbSynCHR_vector_number() ",
                      "- Vector of numbers as nbSynCHR parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimParameters_nbSynCHR_zero <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimParameters(outputDir = "test",
                                                                    fileID = "F1", 
                                                                    nbSynCHR = 0, 
                                                                    methData = samplesForChrSynthetic,
                                                                    nbSimulation = 10,
                                                                    nbBlock = 3, nbCpG  = 2,
                                                                    vNbSample = 2, 
                                                                    nbGeneration = 3, 
                                                                    vpDiff = 0.22, vpDiffsd = 0.1, 
                                                                    vDiff = 0.8, 
                                                                    vInheritance = 0.5,
                                                                    propInherite = 0.8, 
                                                                    rateDiff = 2, 
                                                                    minRate = 1, 
                                                                    propHetero = 0.4, 
                                                                    minReads = 2, 
                                                                    maxPercReads = 99.9, 
                                                                    context = "CpG", assembly = "hg19",
                                                                    meanCov = 10, 
                                                                    keepDiff = TRUE, saveGRanges = FALSE, 
                                                                    saveMethylKit = FALSE,
                                                                    runAnalysis = FALSE,
                                                                    nbCores = 1, vSeed = -1),
                    error=conditionMessage)
    
    exp <- "nbSynCHR must be a positive integer or numeric"
    
    message <- paste0("test.validateRunSimParameters_nbSynCHR_zero() ",
                      "- Zero as nbSynCHR parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimParameters_methBase_number <- function() {
    obs <- tryCatch(
            methylInheritanceSim:::validateRunSimParameters(outputDir = "test",
                                                            fileID = "F1", 
                                                            nbSynCHR = 1, 
                                                            methData = 33, 
                                                            nbSimulation = 10,
                                                            nbBlock = 3, nbCpG  = 2,
                                                            vNbSample = 2, 
                                                            nbGeneration = 3, 
                                                            vpDiff = 0.22, vpDiffsd = 0.1, 
                                                            vDiff = 0.8, 
                                                            vInheritance = 0.5,
                                                            propInherite = 0.8, 
                                                            rateDiff = 2, 
                                                            minRate = 1, 
                                                            propHetero = 0.4, 
                                                            minReads = 2, 
                                                            maxPercReads = 99.9, 
                                                            context = "CpG", assembly = "hg19",
                                                            meanCov = 10, 
                                                            keepDiff = TRUE, saveGRanges = FALSE, 
                                                            saveMethylKit = FALSE,
                                                            runAnalysis = FALSE,
                                                            nbCores = 1, vSeed = -1),
                    error=conditionMessage)
    
    
    exp <- "methData must be an object of class \"methyBase\""
    
    message <- paste0("test.validateRunSimParameters_methBase_number() ",
                      "- Number as methBase parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimParameters_nbSimulation_not_number <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimParameters(outputDir = "test",
                                                                    fileID = "F1", 
                                                                    nbSynCHR = 1, 
                                                                    methData = samplesForChrSynthetic,
                                                                    nbSimulation = "Hi",
                                                                    nbBlock = 3, nbCpG  = 2,
                                                                    vNbSample = 2, 
                                                                    nbGeneration = 3, 
                                                                    vpDiff = 0.2, vpDiffsd = 0.1, 
                                                                    vDiff = 0.8, 
                                                                    vInheritance = 0.5,
                                                                    propInherite = 0.8, 
                                                                    rateDiff = 2, 
                                                                    minRate = 1, 
                                                                    propHetero = 0.4, 
                                                                    minReads = 2, 
                                                                    maxPercReads = 99.9, 
                                                                    context = "CpG", assembly = "hg19",
                                                                    meanCov = 10, 
                                                                    keepDiff = TRUE, saveGRanges = FALSE, 
                                                                    saveMethylKit = FALSE,
                                                                    runAnalysis = FALSE,
                                                                    nbCores = 1, vSeed = -1),
                    error=conditionMessage)
    
    exp <- "nbSimulation must be a positive integer or numeric"
    
    message <- paste0("test.validateRunSimParameters_nbSimulation_not_number() ",
                      "- Not a number as nbSimulation parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimParameters_nbSimulation_vector_number <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimParameters(outputDir = "test",
                                                                    fileID = "F1", 
                                                                    nbSynCHR = 1, 
                                                                    methData = samplesForChrSynthetic,
                                                                    nbSimulation = c(1,2),
                                                                    nbBlock = 3, nbCpG  = 2,
                                                                    vNbSample = 2, 
                                                                    nbGeneration = 3, 
                                                                    vpDiff = 0.2, vpDiffsd = 0.1, 
                                                                    vDiff = 0.8, 
                                                                    vInheritance = 0.5,
                                                                    propInherite = 0.8, 
                                                                    rateDiff = 2, 
                                                                    minRate = 1, 
                                                                    propHetero = 0.4, 
                                                                    minReads = 2, 
                                                                    maxPercReads = 99.9, 
                                                                    context = "CpG", assembly = "hg19",
                                                                    meanCov = 10, 
                                                                    keepDiff = TRUE, saveGRanges = FALSE, 
                                                                    saveMethylKit = FALSE,
                                                                    runAnalysis = FALSE,
                                                                    nbCores = 1, vSeed = -1),
                    error=conditionMessage)
    
    exp <- "nbSimulation must be a positive integer or numeric"
    
    message <- paste0("test.validateRunSimParameters_nbSimulation_vector_number() ",
                      "- Vector of numbers as nbSimulation parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimParameters_nbSimulation_zero <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimParameters(outputDir = "test",
                                                                    fileID = "F1", 
                                                                    nbSynCHR = 1, 
                                                                    methData = samplesForChrSynthetic,
                                                                    nbSimulation = 0,
                                                                    nbBlock = 3, nbCpG  = 2,
                                                                    vNbSample = 2, 
                                                                    nbGeneration = 3, 
                                                                    vpDiff = 0.22, vpDiffsd = 0.1, 
                                                                    vDiff = 0.8, 
                                                                    vInheritance = 0.5,
                                                                    propInherite = 0.8, 
                                                                    rateDiff = 2, 
                                                                    minRate = 1, 
                                                                    propHetero = 0.4, 
                                                                    minReads = 2, 
                                                                    maxPercReads = 99.9, 
                                                                    context = "CpG", assembly = "hg19",
                                                                    meanCov = 10,
                                                                    keepDiff = TRUE, saveGRanges = FALSE, 
                                                                    saveMethylKit = FALSE,
                                                                    runAnalysis = FALSE,
                                                                    nbCores = 1, vSeed = -1),
                    error=conditionMessage)
    
    exp <- "nbSimulation must be a positive integer or numeric"
    
    message <- paste0("test.validateRunSimParameters_nbSimulation_zero() ",
                      "- Zero as nbSimulation parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}


test.validateRunSimParameters_nbBlock_number_vector <- function() {
    obs <- tryCatch(
        methylInheritanceSim:::validateRunSimParameters(outputDir = "test",
                                                        fileID = "F1", 
                                                        nbSynCHR = 1, 
                                                        methData = samplesForChrSynthetic,
                                                        nbSimulation = 10,
                                                        nbBlock = c(3, 1), nbCpG  = 2,
                                                        vNbSample = 2, 
                                                        nbGeneration = 3, 
                                                        vpDiff = 0.22, vpDiffsd = 0.1, 
                                                        vDiff = 0.8, 
                                                        vInheritance = 0.5,
                                                        propInherite = 0.8, 
                                                        rateDiff = 2, 
                                                        minRate = 1, 
                                                        propHetero = 0.4, 
                                                        minReads = 2, 
                                                        maxPercReads = 99.9, 
                                                        context = "CpG", assembly = "hg19",
                                                        meanCov = 10, 
                                                        keepDiff = TRUE, saveGRanges = FALSE, 
                                                        saveMethylKit = FALSE,
                                                        runAnalysis = FALSE,
                                                        nbCores = 1, vSeed = -1),
        error=conditionMessage)
    
    exp <- "nbBlock must be a positive integer or numeric"
    
    message <- paste0("test.validateRunSimParameters_nbBlock_number_vector() ",
                      "- Number vector as nbBlock parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimParameters_nbBlock_string <- function() {
    obs <- tryCatch(
        methylInheritanceSim:::validateRunSimParameters(outputDir = "test",
                                                        fileID = "F1", 
                                                        nbSynCHR = 1, 
                                                        methData = samplesForChrSynthetic,
                                                        nbSimulation = 10,
                                                        nbBlock = "hi", nbCpG  = 2,
                                                        vNbSample = 2, 
                                                        nbGeneration = 3, 
                                                        vpDiff = 0.22, vpDiffsd = 0.1, 
                                                        vDiff = 0.8, 
                                                        vInheritance = 0.5,
                                                        propInherite = 0.8, 
                                                        rateDiff = 2, 
                                                        minRate = 1, 
                                                        propHetero = 0.4, 
                                                        minReads = 2, 
                                                        maxPercReads = 99.9, 
                                                        context = "CpG", assembly = "hg19",
                                                        meanCov = 10, 
                                                        keepDiff = TRUE, saveGRanges = FALSE, 
                                                        saveMethylKit = FALSE,
                                                        runAnalysis = FALSE,
                                                        nbCores = 1, vSeed = -1),
        error=conditionMessage)
    
    exp <- "nbBlock must be a positive integer or numeric"
    
    message <- paste0("test.validateRunSimParameters_nbBlock_string() ",
                      "- String as nbBlock parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimParameters_nbBlock_zero <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimParameters(outputDir = "test",
                                                                    fileID = "F1", 
                                                                    nbSynCHR = 1, 
                                                                    methData = samplesForChrSynthetic,
                                                                    nbSimulation = 10,
                                                                    nbBlock = 0, nbCpG  = 2,
                                                                    vNbSample = 2, 
                                                                    nbGeneration = 3, 
                                                                    vpDiff = 0.22, vpDiffsd = 0.1, 
                                                                    vDiff = 0.8, 
                                                                    vInheritance = 0.5,
                                                                    propInherite = 0.8, 
                                                                    rateDiff = 2, 
                                                                    minRate = 1, 
                                                                    propHetero = 0.4, 
                                                                    minReads = 2, 
                                                                    maxPercReads = 99.9, 
                                                                    context = "CpG", assembly = "hg19",
                                                                    meanCov = 10,
                                                                    keepDiff = TRUE, saveGRanges = FALSE, 
                                                                    saveMethylKit = FALSE,
                                                                    runAnalysis = FALSE,
                                                                    nbCores = 1, vSeed = -1),
                    error=conditionMessage)
    
    exp <- "nbBlock must be a positive integer or numeric"
    
    message <- paste0("test.validateRunSimParameters_nbBlock_zero() ",
                      "- Zero as nbBlock parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}


test.validateRunSimParameters_nbCpG_number_vector <- function() {
    obs <- tryCatch(
        methylInheritanceSim:::validateRunSimParameters(outputDir = "test",
                                                        fileID = "F1", 
                                                        nbSynCHR = 1, 
                                                        methData = samplesForChrSynthetic,
                                                        nbSimulation = 10,
                                                        nbBlock = 2, nbCpG  = c(3,2),
                                                        vNbSample = 2, 
                                                        nbGeneration = 3, 
                                                        vpDiff = 0.2, vpDiffsd = 0.1, 
                                                        vDiff = 0.8, 
                                                        vInheritance = 0.5,
                                                        propInherite = 0.8, 
                                                        rateDiff = 2, 
                                                        minRate = 1, 
                                                        propHetero = 0.4, 
                                                        minReads = 2, 
                                                        maxPercReads = 99.9, 
                                                        context = "CpG", assembly = "hg19",
                                                        meanCov = 10, 
                                                        keepDiff = TRUE, saveGRanges = FALSE, 
                                                        saveMethylKit = FALSE,
                                                        runAnalysis = FALSE,
                                                        nbCores = 1, vSeed = -1),
        error=conditionMessage)
    
    exp <- "nbCpG must be a positive integer or numeric"
    
    message <- paste0("test.validateRunSimParameters_nbCpG_number_vector() ",
                      "- Number vector as nbCpG parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimParameters_nbCpG_string <- function() {
    obs <- tryCatch(
        methylInheritanceSim:::validateRunSimParameters(outputDir = "test",
                                                        fileID = "F1", 
                                                        nbSynCHR = 1, 
                                                        methData = samplesForChrSynthetic,
                                                        nbSimulation = 10,
                                                        nbBlock = 2, nbCpG  = "hi",
                                                        vNbSample = 2, 
                                                        nbGeneration = 3, 
                                                        vpDiff = 0.2, vpDiffsd = 0.1, 
                                                        vDiff = 0.8, 
                                                        vInheritance = 0.5,
                                                        propInherite = 0.8, 
                                                        rateDiff = 2, 
                                                        minRate = 1, 
                                                        propHetero = 0.4, 
                                                        minReads = 2, 
                                                        maxPercReads = 99.9, 
                                                        context = "CpG", assembly = "hg19",
                                                        meanCov = 10, 
                                                        keepDiff = TRUE, saveGRanges = FALSE, 
                                                        saveMethylKit = FALSE,
                                                        runAnalysis = FALSE,
                                                        nbCores = 1, vSeed = -1),
        error=conditionMessage)
    
    exp <- "nbCpG must be a positive integer or numeric"
    
    message <- paste0("test.validateRunSimParameters_nbCpG_string() ",
                      "- String as nbCpG parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimParameters_nbCpG_zero <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimParameters(outputDir = "test",
                                                                    fileID = "F1", 
                                                                    nbSynCHR = 1, 
                                                                    methData = samplesForChrSynthetic,
                                                                    nbSimulation = 10,
                                                                    nbBlock = 2, nbCpG  = 0,
                                                                    vNbSample = 2, 
                                                                    nbGeneration = 3, 
                                                                    vpDiff = 0.22, vpDiffsd = 0.1, 
                                                                    vDiff = 0.8, 
                                                                    vInheritance = 0.5,
                                                                    propInherite = 0.8, 
                                                                    rateDiff = 2, 
                                                                    minRate = 1, 
                                                                    propHetero = 0.4, 
                                                                    minReads = 2, 
                                                                    maxPercReads = 99.9, 
                                                                    context = "CpG", assembly = "hg19",
                                                                    meanCov = 10,
                                                                    keepDiff = TRUE, saveGRanges = FALSE, 
                                                                    saveMethylKit = FALSE,
                                                                    runAnalysis = FALSE,
                                                                    nbCores = 1, vSeed = -1),
                    error=conditionMessage)
    
    exp <- "nbCpG must be a positive integer or numeric"
    
    message <- paste0("test.validateRunSimParameters_nbCpG_zero() ",
                      "- Zero as nbCpG parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimParameters_nbGeneration_number_vector <- function() {
    obs <- tryCatch(
        methylInheritanceSim:::validateRunSimParameters(outputDir = "test",
                                                        fileID = "F1", 
                                                        nbSynCHR = 1, 
                                                        methData = samplesForChrSynthetic,
                                                        nbSimulation = 10,
                                                        nbBlock = 2, nbCpG  = 2,
                                                        nbGeneration = c(3,2), 
                                                        vNbSample = 2,
                                                        vpDiff = 0.2, vpDiffsd = 0.1, 
                                                        vDiff = 0.8, 
                                                        vInheritance = 0.5,
                                                        propInherite = 0.8, 
                                                        rateDiff = 2, 
                                                        minRate = 1, 
                                                        propHetero = 0.4, 
                                                        minReads = 2, 
                                                        maxPercReads = 99.9, 
                                                        context = "CpG", assembly = "hg19",
                                                        meanCov = 10, 
                                                        keepDiff = TRUE, saveGRanges = FALSE, 
                                                        saveMethylKit = FALSE,
                                                        runAnalysis = FALSE,
                                                        nbCores = 1, vSeed = -1),
        error=conditionMessage)
    
    exp <- "nbGeneration must be a positive integer or numeric"
    
    message <- paste0("test.validateRunSimParameters_nbGeneration_number_vector() ",
                      "- Number vector as nbGeneration parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimParameters_nbGeneration_string <- function() {
    obs <- tryCatch(
        methylInheritanceSim:::validateRunSimParameters(outputDir = "test",
                                                        fileID = "F1", 
                                                        nbSynCHR = 1, 
                                                        methData = samplesForChrSynthetic,
                                                        nbSimulation = 10,
                                                        nbBlock = 2, nbCpG  = 2,
                                                        nbGeneration = "Hi",
                                                        vNbSample = 2, 
                                                        vpDiff = 0.2, vpDiffsd = 0.1, 
                                                        vDiff = 0.8, 
                                                        vInheritance = 0.5,
                                                        propInherite = 0.8, 
                                                        rateDiff = 2, 
                                                        minRate = 1, 
                                                        propHetero = 0.4, 
                                                        minReads = 2, 
                                                        maxPercReads = 99.9, 
                                                        context = "CpG", assembly = "hg19",
                                                        meanCov = 10, 
                                                        keepDiff = TRUE, saveGRanges = FALSE, 
                                                        saveMethylKit = FALSE,
                                                        runAnalysis = FALSE,
                                                        nbCores = 1, vSeed = -1),
        error=conditionMessage)
    
    exp <- "nbGeneration must be a positive integer or numeric"
    
    message <- paste0("test.validateRunSimParameters_nbGeneration_string() ",
                      "- String as nbGeneration parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimParameters_nbGeneration_zero <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimParameters(outputDir = "test",
                                                                    fileID = "F1", 
                                                                    nbSynCHR = 1, 
                                                                    methData = samplesForChrSynthetic,
                                                                    nbSimulation = 10,
                                                                    nbBlock = 2, nbCpG  = 2,
                                                                    nbGeneration = 0,
                                                                    vNbSample = 2,
                                                                    vpDiff = 0.22, vpDiffsd = 0.1, 
                                                                    vDiff = 0.8, 
                                                                    vInheritance = 0.5,
                                                                    propInherite = 0.8, 
                                                                    rateDiff = 2, 
                                                                    minRate = 1, 
                                                                    propHetero = 0.4, 
                                                                    minReads = 2, 
                                                                    maxPercReads = 99.9, 
                                                                    context = "CpG", assembly = "hg19",
                                                                    meanCov = 10,
                                                                    keepDiff = TRUE, saveGRanges = FALSE, 
                                                                    saveMethylKit = FALSE,
                                                                    runAnalysis = FALSE,
                                                                    nbCores = 1, vSeed = -1),
                    error=conditionMessage)
    
    exp <- "nbGeneration must be a positive integer or numeric"
    
    message <- paste0("test.validateRunSimParameters_nbGeneration_zero() ",
                      "- Zero as nbGeneration parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimParameters_vNbSample_string <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimParameters(outputDir = "test",
                                                                    fileID = "F1", 
                                                                    nbSynCHR = 1, 
                                                                    methData = samplesForChrSynthetic,
                                                                    nbSimulation = 10,
                                                                    nbBlock = 2, nbCpG  = 2,
                                                                    nbGeneration = 3,
                                                                    vNbSample = "Hi",
                                                                    vpDiff = 0.22, vpDiffsd = 0.1, 
                                                                    vDiff = 0.8, 
                                                                    vInheritance = 0.5,
                                                                    propInherite = 0.8, 
                                                                    rateDiff = 2, 
                                                                    minRate = 1, 
                                                                    propHetero = 0.4, 
                                                                    minReads = 2, 
                                                                    maxPercReads = 99.9, 
                                                                    context = "CpG", assembly = "hg19",
                                                                    meanCov = 10,
                                                                    keepDiff = TRUE, saveGRanges = FALSE, 
                                                                    saveMethylKit = FALSE,
                                                                    runAnalysis = FALSE,
                                                                    nbCores = 1, vSeed = -1),
                    error=conditionMessage)
    
    exp <- "vNbSample must be a vector of distinct positive integer"
    
    message <- paste0("test.validateRunSimParameters_vNbSample_string() ",
                      "- String as vNbSample parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimParameters_vNbSample_duplicate <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimParameters(outputDir = "test",
                                                                    fileID = "F1", 
                                                                    nbSynCHR = 1, 
                                                                    methData = samplesForChrSynthetic,
                                                                    nbSimulation = 10,
                                                                    nbBlock = 2, nbCpG  = 2,
                                                                    nbGeneration = 3,
                                                                    vNbSample = c(1, 2, 2),
                                                                    vpDiff = 0.22, vpDiffsd = 0.1, 
                                                                    vDiff = 0.8, 
                                                                    vInheritance = 0.5,
                                                                    propInherite = 0.8, 
                                                                    rateDiff = 2, 
                                                                    minRate = 1, 
                                                                    propHetero = 0.4, 
                                                                    minReads = 2, 
                                                                    maxPercReads = 99.9, 
                                                                    context = "CpG", assembly = "hg19",
                                                                    meanCov = 10,
                                                                    keepDiff = TRUE, saveGRanges = FALSE, 
                                                                    saveMethylKit = FALSE,
                                                                    runAnalysis = FALSE,
                                                                    nbCores = 1, vSeed = -1),
                    error=conditionMessage)
    
    exp <- "vNbSample must be a vector of distinct positive integer"
    
    message <- paste0("test.validateRunSimParameters_vNbSample_duplicate() ",
                    "- Duplicate elements in vNbSample parameter did not 
                    generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimParameters_vNbSample_float <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimParameters(outputDir = "test",
                                                                    fileID = "F1", 
                                                                    nbSynCHR = 1, 
                                                                    methData = samplesForChrSynthetic,
                                                                    nbSimulation = 10,
                                                                    nbBlock = 2, nbCpG  = 2,
                                                                    nbGeneration = 3,
                                                                    vNbSample = c(1, 2, 2.1),
                                                                    vpDiff = 0.22, vpDiffsd = 0.1, 
                                                                    vDiff = 0.8, 
                                                                    vInheritance = 0.5,
                                                                    propInherite = 0.8, 
                                                                    rateDiff = 2, 
                                                                    minRate = 1, 
                                                                    propHetero = 0.4, 
                                                                    minReads = 2, 
                                                                    maxPercReads = 99.9, 
                                                                    context = "CpG", assembly = "hg19",
                                                                    meanCov = 10,
                                                                    keepDiff = TRUE, saveGRanges = FALSE, 
                                                                    saveMethylKit = FALSE,
                                                                    runAnalysis = FALSE,
                                                                    nbCores = 1, vSeed = -1),
                    error=conditionMessage)
    
    exp <- "vNbSample must be a vector of distinct positive integer"
    
    message <- paste0("test.validateRunSimParameters_vNbSample_float() ",
                      "- Float as vNbSample parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimParameters_vpDiff_string <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimParameters(outputDir = "test",
                                                                    fileID = "F1", 
                                                                    nbSynCHR = 1, 
                                                                    methData = samplesForChrSynthetic,
                                                                    nbSimulation = 10,
                                                                    nbBlock = 2, nbCpG  = 2,
                                                                    nbGeneration = 3,
                                                                    vNbSample = 3,
                                                                    vpDiff = "Hi", vpDiffsd = 0.1, 
                                                                    vDiff = 0.8, 
                                                                    vInheritance = 0.5,
                                                                    propInherite = 0.8, 
                                                                    rateDiff = 2, 
                                                                    minRate = 1, 
                                                                    propHetero = 0.4, 
                                                                    minReads = 2, 
                                                                    maxPercReads = 99.9, 
                                                                    context = "CpG", assembly = "hg19",
                                                                    meanCov = 10,
                                                                    keepDiff = TRUE, saveGRanges = FALSE, 
                                                                    saveMethylKit = FALSE,
                                                                    runAnalysis = FALSE,
                                                                    nbCores = 1, vSeed = -1),
                    error=conditionMessage)
    
    exp <- "vpDiff must be a vector of distinct positive double include in (0,1]"
    
    message <- paste0("test.validateRunSimParameters_vpDiff_string() ",
                      "- String as vpDiff parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimParameters_vpDiff_duplicate <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimParameters(outputDir = "test",
                                                                    fileID = "F1", 
                                                                    nbSynCHR = 1, 
                                                                    methData = samplesForChrSynthetic,
                                                                    nbSimulation = 10,
                                                                    nbBlock = 2, nbCpG  = 2,
                                                                    nbGeneration = 3,
                                                                    vNbSample = 3,
                                                                    vpDiff = c(0.5, 0.5, 0.8), vpDiffsd = 0.1, 
                                                                    vDiff = 0.9, 
                                                                    vInheritance = 0.5,
                                                                    propInherite = 0.8, 
                                                                    rateDiff = 2, 
                                                                    minRate = 1, 
                                                                    propHetero = 0.4, 
                                                                    minReads = 2, 
                                                                    maxPercReads = 99.9, 
                                                                    context = "CpG", assembly = "hg19",
                                                                    meanCov = 10,
                                                                    keepDiff = TRUE, saveGRanges = FALSE, 
                                                                    saveMethylKit = FALSE,
                                                                    runAnalysis = FALSE,
                                                                    nbCores = 1, vSeed = -1),
                    error=conditionMessage)
    
    exp <- "vpDiff must be a vector of distinct positive double include in (0,1]"
    
    message <- paste0("test.validateRunSimParameters_vpDiff_duplicate() ",
                      "- Duplicate elements in vpDiff parameter did not 
                      generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimParameters_vpDiff_zero <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimParameters(outputDir = "test",
                                                                    fileID = "F1", 
                                                                    nbSynCHR = 1, 
                                                                    methData = samplesForChrSynthetic,
                                                                    nbSimulation = 10,
                                                                    nbBlock = 2, nbCpG  = 2,
                                                                    nbGeneration = 3,
                                                                    vNbSample = 3,
                                                                    vpDiff = c(0.5, 0, 0.8), vpDiffsd = 0.1, 
                                                                    vDiff = 0.8, 
                                                                    vInheritance = 0.5,
                                                                    propInherite = 0.8, 
                                                                    rateDiff = 2, 
                                                                    minRate = 1, 
                                                                    propHetero = 0.4, 
                                                                    minReads = 2, 
                                                                    maxPercReads = 99.9, 
                                                                    context = "CpG", assembly = "hg19",
                                                                    meanCov = 10,
                                                                    keepDiff = TRUE, saveGRanges = FALSE, 
                                                                    saveMethylKit = FALSE,
                                                                    runAnalysis = FALSE,
                                                                    nbCores = 1, vSeed = -1),
                    error=conditionMessage)
    
    exp <- "vpDiff must be a vector of distinct positive double include in (0,1]"
    
    message <- paste0("test.validateRunSimParameters_vpDiff_zero() ",
                      "- Zero in vpDiff parameter did not 
                      generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimParameters_vpDiff_two <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimParameters(outputDir = "test",
                                                                    fileID = "F1", 
                                                                    nbSynCHR = 1, 
                                                                    methData = samplesForChrSynthetic,
                                                                    nbSimulation = 10,
                                                                    nbBlock = 2, nbCpG  = 2,
                                                                    nbGeneration = 3,
                                                                    vNbSample = 3,
                                                                    vpDiff = c(0.5, 2, 0.8), vpDiffsd = 0.1, 
                                                                    vDiff = 0.8, 
                                                                    vInheritance = 0.5,
                                                                    propInherite = 0.8, 
                                                                    rateDiff = 2, 
                                                                    minRate = 1, 
                                                                    propHetero = 0.4, 
                                                                    minReads = 2, 
                                                                    maxPercReads = 99.9, 
                                                                    context = "CpG", assembly = "hg19",
                                                                    meanCov = 10,
                                                                    keepDiff = TRUE, saveGRanges = FALSE, 
                                                                    saveMethylKit = FALSE,
                                                                    runAnalysis = FALSE,
                                                                    nbCores = 1, vSeed = -1),
                    error=conditionMessage)
    
    exp <- "vpDiff must be a vector of distinct positive double include in (0,1]"
    
    message <- paste0("test.validateRunSimParameters_vpDiff_two() ",
                      "- Two in vpDiff parameter did not 
                      generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimParameters_vpDiffsd_string <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimParameters(outputDir = "test",
                                                                    fileID = "F1", 
                                                                    nbSynCHR = 1, 
                                                                    methData = samplesForChrSynthetic,
                                                                    nbSimulation = 10,
                                                                    nbBlock = 2, nbCpG  = 2,
                                                                    nbGeneration = 3,
                                                                    vNbSample = 3,
                                                                    vpDiff = 0.9, vpDiffsd = "Hi", 
                                                                    vDiff = 0.8, 
                                                                    vInheritance = 0.5,
                                                                    propInherite = 0.8, 
                                                                    rateDiff = 2, 
                                                                    minRate = 1, 
                                                                    propHetero = 0.4, 
                                                                    minReads = 2, 
                                                                    maxPercReads = 99.9, 
                                                                    context = "CpG", assembly = "hg19",
                                                                    meanCov = 10,
                                                                    keepDiff = TRUE, saveGRanges = FALSE, 
                                                                    saveMethylKit = FALSE,
                                                                    runAnalysis = FALSE,
                                                                    nbCores = 1, vSeed = -1),
                    error=conditionMessage)
    
    exp <- "vpDiffsd must be a vector of non-negative double"
    
    message <- paste0("test.validateRunSimParameters_vpDiffsd_string() ",
                      "- String as vpDiffsd parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimParameters_vpDiffsd_negative <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimParameters(outputDir = "test",
                                                                    fileID = "F1", 
                                                                    nbSynCHR = 1, 
                                                                    methData = samplesForChrSynthetic,
                                                                    nbSimulation = 10,
                                                                    nbBlock = 2, nbCpG  = 2,
                                                                    nbGeneration = 3,
                                                                    vNbSample = 3,
                                                                    vpDiff = 0.9, vpDiffsd = -1, 
                                                                    vDiff = 0.8, 
                                                                    vInheritance = 0.5,
                                                                    propInherite = 0.8, 
                                                                    rateDiff = 2, 
                                                                    minRate = 1, 
                                                                    propHetero = 0.4, 
                                                                    minReads = 2, 
                                                                    maxPercReads = 99.9, 
                                                                    context = "CpG", assembly = "hg19",
                                                                    meanCov = 10,
                                                                    keepDiff = TRUE, saveGRanges = FALSE, 
                                                                    saveMethylKit = FALSE,
                                                                    runAnalysis = FALSE,
                                                                    nbCores = 1, vSeed = -1),
                    error=conditionMessage)
    
    exp <- "vpDiffsd must be a vector of non-negative double"
    
    message <- paste0("test.validateRunSimParameters_vpDiffsd_negative() ",
                      "- Negative as vpDiffsd parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimParameters_vpDiffsd_vpDiff <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimParameters(outputDir = "test",
                                                                    fileID = "F1", 
                                                                    nbSynCHR = 1, 
                                                                    methData = samplesForChrSynthetic,
                                                                    nbSimulation = 10,
                                                                    nbBlock = 2, nbCpG  = 2,
                                                                    nbGeneration = 3,
                                                                    vNbSample = 3,
                                                                    vpDiff = 0.9, vpDiffsd = c(0.1, 0.1), 
                                                                    vDiff = 0.8, 
                                                                    vInheritance = 0.5,
                                                                    propInherite = 0.8, 
                                                                    rateDiff = 2, 
                                                                    minRate = 1, 
                                                                    propHetero = 0.4, 
                                                                    minReads = 2, 
                                                                    maxPercReads = 99.9, 
                                                                    context = "CpG", assembly = "hg19",
                                                                    meanCov = 10,
                                                                    keepDiff = TRUE, saveGRanges = FALSE, 
                                                                    saveMethylKit = FALSE,
                                                                    runAnalysis = FALSE,
                                                                    nbCores = 1, vSeed = -1),
                    error=conditionMessage)
    
    exp <- "vpDiff and vpDiffsd must be the same length"
    
    message <- paste0("test.validateRunSimParameters_vpDiffsd_vpDiff() ",
                      "- vpDiffsd vpDiff not same length did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimParameters_vDiff_string <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimParameters(outputDir = "test",
                                                                    fileID = "F1", 
                                                                    nbSynCHR = 1, 
                                                                    methData = samplesForChrSynthetic,
                                                                    nbSimulation = 10,
                                                                    nbBlock = 2, nbCpG  = 2,
                                                                    nbGeneration = 3,
                                                                    vNbSample = 3,
                                                                    vpDiff = 0.9, vpDiffsd = 0.1, 
                                                                    vDiff = "Hi", 
                                                                    vInheritance = 0.5,
                                                                    propInherite = 0.8, 
                                                                    rateDiff = 2, 
                                                                    minRate = 1, 
                                                                    propHetero = 0.4, 
                                                                    minReads = 2, 
                                                                    maxPercReads = 99.9, 
                                                                    context = "CpG", assembly = "hg19",
                                                                    meanCov = 10,
                                                                    keepDiff = TRUE, saveGRanges = FALSE, 
                                                                    saveMethylKit = FALSE,
                                                                    runAnalysis = FALSE,
                                                                    nbCores = 1, vSeed = -1),
                    error=conditionMessage)
    
    exp <- "vDiff must be a vector of distinct non-negative double include in [0,1]"
    
    message <- paste0("test.validateRunSimParameters_vDiff_string() ",
                      "- String as vDiff parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimParameters_vDiff_duplicate <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimParameters(outputDir = "test",
                                                                    fileID = "F1", 
                                                                    nbSynCHR = 1, 
                                                                    methData = samplesForChrSynthetic,
                                                                    nbSimulation = 10,
                                                                    nbBlock = 2, nbCpG  = 2,
                                                                    nbGeneration = 3,
                                                                    vNbSample = 3,
                                                                    vpDiff = c(0.5), vpDiffsd = 0.1, 
                                                                    vDiff = c(0.5, 0.5, 0.8), 
                                                                    vInheritance = 0.5,
                                                                    propInherite = 0.8, 
                                                                    rateDiff = 2, 
                                                                    minRate = 1, 
                                                                    propHetero = 0.4, 
                                                                    minReads = 2, 
                                                                    maxPercReads = 99.9, 
                                                                    context = "CpG", assembly = "hg19",
                                                                    meanCov = 10,
                                                                    keepDiff = TRUE, saveGRanges = FALSE, 
                                                                    saveMethylKit = FALSE,
                                                                    runAnalysis = FALSE,
                                                                    nbCores = 1, vSeed = -1),
                    error=conditionMessage)
    
    exp <- "vDiff must be a vector of distinct non-negative double include in [0,1]"
    
    message <- paste0("test.validateRunSimParameters_vDiff_duplicate() ",
                      "- Duplicate elements in vDiff parameter did not 
                      generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimParameters_vDiff_negative <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimParameters(outputDir = "test",
                                                                    fileID = "F1", 
                                                                    nbSynCHR = 1, 
                                                                    methData = samplesForChrSynthetic,
                                                                    nbSimulation = 10,
                                                                    nbBlock = 2, nbCpG  = 2,
                                                                    nbGeneration = 3,
                                                                    vNbSample = 3,
                                                                    vpDiff = c(0.5), vpDiffsd = 0.1, 
                                                                    vDiff = -1, 
                                                                    vInheritance = 0.5,
                                                                    propInherite = 0.8, 
                                                                    rateDiff = 2, 
                                                                    minRate = 1, 
                                                                    propHetero = 0.4, 
                                                                    minReads = 2, 
                                                                    maxPercReads = 99.9, 
                                                                    context = "CpG", assembly = "hg19",
                                                                    meanCov = 10,
                                                                    keepDiff = TRUE, saveGRanges = FALSE, 
                                                                    saveMethylKit = FALSE,
                                                                    runAnalysis = FALSE,
                                                                    nbCores = 1, vSeed = -1),
                    error=conditionMessage)
    
    exp <- "vDiff must be a vector of distinct non-negative double include in [0,1]"
    
    message <- paste0("test.validateRunSimParameters_vDiff_negative() ",
                      "- Negative in vDiff parameter did not 
                      generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimParameters_vDiff_two <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimParameters(outputDir = "test",
                                                                    fileID = "F1", 
                                                                    nbSynCHR = 1, 
                                                                    methData = samplesForChrSynthetic,
                                                                    nbSimulation = 10,
                                                                    nbBlock = 2, nbCpG  = 2,
                                                                    nbGeneration = 3,
                                                                    vNbSample = 3,
                                                                    vpDiff = c(0.5), vpDiffsd = 0.1, 
                                                                    vDiff = 2, 
                                                                    vInheritance = 0.5,
                                                                    propInherite = 0.8, 
                                                                    rateDiff = 2, 
                                                                    minRate = 1, 
                                                                    propHetero = 0.4, 
                                                                    minReads = 2, 
                                                                    maxPercReads = 99.9, 
                                                                    context = "CpG", assembly = "hg19",
                                                                    meanCov = 10,
                                                                    keepDiff = TRUE, saveGRanges = FALSE, 
                                                                    saveMethylKit = FALSE,
                                                                    runAnalysis = FALSE,
                                                                    nbCores = 1, vSeed = -1),
                    error=conditionMessage)
    
    exp <- "vDiff must be a vector of distinct non-negative double include in [0,1]"
    
    message <- paste0("test.validateRunSimParameters_vDiff_two() ",
                      "- Two in vDiff parameter did not 
                      generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimParameters_vInheritance_string <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimParameters(outputDir = "test",
                                                                    fileID = "F1", 
                                                                    nbSynCHR = 1, 
                                                                    methData = samplesForChrSynthetic,
                                                                    nbSimulation = 10,
                                                                    nbBlock = 2, nbCpG  = 2,
                                                                    nbGeneration = 3,
                                                                    vNbSample = 3,
                                                                    vpDiff = 0.9, vpDiffsd = 0.1, 
                                                                    vDiff = 0.8, 
                                                                    vInheritance = "Hi",
                                                                    propInherite = 0.8, 
                                                                    rateDiff = 2, 
                                                                    minRate = 1, 
                                                                    propHetero = 0.4, 
                                                                    minReads = 2, 
                                                                    maxPercReads = 99.9, 
                                                                    context = "CpG", assembly = "hg19",
                                                                    meanCov = 10,
                                                                    keepDiff = TRUE, saveGRanges = FALSE, 
                                                                    saveMethylKit = FALSE,
                                                                    runAnalysis = FALSE,
                                                                    nbCores = 1, vSeed = -1),
                    error=conditionMessage)
    
    exp <- paste0("vInheritance must be a vector of distinct non-negative ", 
                  "double include in [0,1]")
    
    message <- paste0("test.validateRunSimParameters_vInheritance_string() ",
                    "- String as vInheritance parameter did not generated ", 
                    "expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimParameters_vInheritance_duplicate <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimParameters(outputDir = "test",
                                                                    fileID = "F1", 
                                                                    nbSynCHR = 1, 
                                                                    methData = samplesForChrSynthetic,
                                                                    nbSimulation = 10,
                                                                    nbBlock = 2, nbCpG  = 2,
                                                                    nbGeneration = 3,
                                                                    vNbSample = 3,
                                                                    vpDiff = c(0.5), vpDiffsd = 0.1, 
                                                                    vDiff = c(0.5, 0.7, 0.8), 
                                                                    vInheritance = c(0.5, 0.5, 0.7),
                                                                    propInherite = 0.8, 
                                                                    rateDiff = 2, 
                                                                    minRate = 1, 
                                                                    propHetero = 0.4, 
                                                                    minReads = 2, 
                                                                    maxPercReads = 99.9, 
                                                                    context = "CpG", assembly = "hg19",
                                                                    meanCov = 10,
                                                                    keepDiff = TRUE, saveGRanges = FALSE, 
                                                                    saveMethylKit = FALSE,
                                                                    runAnalysis = FALSE,
                                                                    nbCores = 1, vSeed = -1),
                    error=conditionMessage)
    
    exp <- paste0("vInheritance must be a vector of distinct non-negative ", 
                  "double include in [0,1]")
    
    message <- paste0("test.validateRunSimParameters_vInheritance_duplicate() ",
                      "- Duplicate elements in vInheritance parameter did not 
                      generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimParameters_vInheritance_negative <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimParameters(outputDir = "test",
                                                                    fileID = "F1", 
                                                                    nbSynCHR = 1, 
                                                                    methData = samplesForChrSynthetic,
                                                                    nbSimulation = 10,
                                                                    nbBlock = 2, nbCpG  = 2,
                                                                    nbGeneration = 3,
                                                                    vNbSample = 3,
                                                                    vpDiff = c(0.5), vpDiffsd = 0.1, 
                                                                    vDiff = 0.8, 
                                                                    vInheritance = -1,
                                                                    propInherite = 0.8, 
                                                                    rateDiff = 2, 
                                                                    minRate = 1, 
                                                                    propHetero = 0.4, 
                                                                    minReads = 2, 
                                                                    maxPercReads = 99.9, 
                                                                    context = "CpG", assembly = "hg19",
                                                                    meanCov = 10,
                                                                    keepDiff = TRUE, saveGRanges = FALSE, 
                                                                    saveMethylKit = FALSE,
                                                                    runAnalysis = FALSE,
                                                                    nbCores = 1, vSeed = -1),
                    error=conditionMessage)
    
    exp <- paste0("vInheritance must be a vector of distinct non-negative ", 
                    "double include in [0,1]")
    
    message <- paste0("test.validateRunSimParameters_vInheritance_negative() ",
                      "- Negative in vInheritance parameter did not 
                      generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimParameters_vInheritance_two <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimParameters(outputDir = "test",
                                                                    fileID = "F1", 
                                                                    nbSynCHR = 1, 
                                                                    methData = samplesForChrSynthetic,
                                                                    nbSimulation = 10,
                                                                    nbBlock = 2, nbCpG  = 2,
                                                                    nbGeneration = 3,
                                                                    vNbSample = 3,
                                                                    vpDiff = c(0.5), vpDiffsd = 0.1, 
                                                                    vDiff = 0.8, 
                                                                    vInheritance = 2,
                                                                    propInherite = 0.8, 
                                                                    rateDiff = 2, 
                                                                    minRate = 1, 
                                                                    propHetero = 0.4, 
                                                                    minReads = 2, 
                                                                    maxPercReads = 99.9, 
                                                                    context = "CpG", assembly = "hg19",
                                                                    meanCov = 10,
                                                                    keepDiff = TRUE, saveGRanges = FALSE, 
                                                                    saveMethylKit = FALSE,
                                                                    runAnalysis = FALSE,
                                                                    nbCores = 1, vSeed = -1),
                    error=conditionMessage)
    
    exp <- paste0("vInheritance must be a vector of distinct non-negative ", 
                  "double include in [0,1]")
    
    message <- paste0("test.validateRunSimParameters_vInheritance_two() ",
                      "- Two in vInheritance parameter did not 
                      generated expected results.")
    
    checkEquals(obs, exp, message)
}


test.validateRunSimParameters_nbCores_string <- function() {
    obs <- tryCatch(
        methylInheritanceSim:::validateRunSimParameters(outputDir = "test",
                                                        fileID = "F1", 
                                                        nbSynCHR = 1, 
                                                        methData = samplesForChrSynthetic,
                                                        nbBlock = 2, nbCpG  = 2,
                                                        vNbSample = 2, 
                                                        nbGeneration = 3, 
                                                        vpDiff = 0.2, vpDiffsd = 0.1, 
                                                        vDiff = 0.2, 
                                                        vInheritance = 0.2,
                                                        propInherite = 0.8, 
                                                        rateDiff = 0.8, 
                                                        minRate = 0.1, 
                                                        propHetero = 0.4, 
                                                        minReads = 2, 
                                                        maxPercReads = 99.9, 
                                                        context = "CpG", assembly = "hg19",
                                                        meanCov = 10, n = 3, 
                                                        keepDiff = TRUE, saveGRanges = FALSE, 
                                                        saveMethylKit = FALSE,
                                                        runAnalysis = FALSE,
                                                        nbCores = "hi", vSeed = -1),
        error=conditionMessage)
    
    exp <- "nbCores must be a positive integer or numeric"
    
    message <- paste0("test.validateRunSimParameters_nbCores_string() ",
                    "- String as nbCores parameter did not generated expected
                    results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimParameters_nbCores_zero <- function() {
    obs <- tryCatch(
        methylInheritanceSim:::validateRunSimParameters(outputDir = "test",
                                                        fileID = "F1", 
                                                        nbSynCHR = 1, 
                                                        methData = samplesForChrSynthetic,
                                                        nbBlock = 2, nbCpG = 2,
                                                        vNbSample = 2, 
                                                        nbGeneration = 3, 
                                                        vpDiff = 0.2, vpDiffsd = 0.1, 
                                                        vDiff = 0.2, 
                                                        vInheritance = 0.2,
                                                        propInherite = 0.8, 
                                                        rateDiff = 0.8, 
                                                        minRate = 0.1, 
                                                        propHetero = 0.4, 
                                                        minReads = 2, 
                                                        maxPercReads = 99.9, 
                                                        context = "CpG", assembly = "hg19",
                                                        meanCov = 10, n = 3, 
                                                        keepDiff = TRUE, saveGRanges = FALSE, 
                                                        saveMethylKit = FALSE,
                                                        runAnalysis = FALSE,
                                                        nbCores = 0, vSeed = -1),
        error=conditionMessage)
    
    exp <- "nbCores must be a positive integer or numeric"
    
    message <- paste0("test.validateRunSimParameters_nbCores_zero() ",
                      "- Zero as nbCores parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}
