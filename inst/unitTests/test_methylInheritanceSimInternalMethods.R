###################################################
##
## Test the methylInheritanceSimInternalMethods functions
##
###################################################

data(samplesForChrSynthetic)


###################################################
## estBetaAlpha() function
###################################################

test.estBetaAlpha_good_01 <- function() {
    
    obs <- methylInheritanceSim:::estBetaAlpha(c(0.5, 0.2), 0.00001)
    
    exp <- 0.125
    
    message <- paste0("test.estBetaAlpha_good_01() ",
                    "- Valid paramters did not generated expected results.")
    
    checkEquals(obs, exp, message)
}


###################################################
## estBetaBeta() function
###################################################

test.estBetaBeta_good_01 <- function() {
    
    obs <- methylInheritanceSim:::estBetaBeta(c(0.3, 0.2), 0.00001)
    
    exp <- 0.035
    
    message <- paste0("test.estBetaBeta_good_01() ",
                      "- Valid paramters did not generated expected results.")
    
    checkEquals(obs, exp, message)
}


###################################################
##getSyntheticChr() function
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
                      "- Valid paramters did not generated expected results.")
    
    checkEquals(obs, exp, message)
}


###################################################
##validateRunSimParameters() function
###################################################


###################################################
## getSim() function
###################################################

test.validateRunSimParameters_outputDir_number <- function() {
    
    obs <- tryCatch(methylInheritanceSim:::validateRunSimParameters(outputDir = 33,
                                    fileID = "F1", 
                                    nbSynCHR = 1, 
                                    methData = samplesForChrSynthetic, 
                                    nbBlock = 3, lBlock  = 2,
                                    vNbSample = 2, 
                                    nbGeneration = 3, 
                                    vpDiff = 2, vpDiffsd = 1, 
                                    vDiff = 2, 
                                    vInheritance = 2,
                                    propInherite = 0.8, 
                                    rateDiff = 2, 
                                    minRate = 1, 
                                    propHetero = 0.4, 
                                    minReads = 2, 
                                    maxPercReads = 99.9, 
                                    context = "CpG", assembly = "hg19",
                                    meanCov = 10, n = 3, 
                                    keepDiff = TRUE, saveGRanges = FALSE, 
                                    saveMethylKit = FALSE,
                                    anaMethylKit = FALSE,
                                    nbCores = 1, vSeed = -1),
                            error=conditionMessage)

    
    exp <- "outputDir must be a character string or NULL"
    
    
    message <- paste0("test.validateRunSimParameters_outputDir_number() ",
                      "- Number as outputDir paramter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimParameters_fileID_number <- function() {
    
    obs <- tryCatch(methylInheritanceSim:::validateRunSimParameters(outputDir = "test",
                                                                    fileID = 2, 
                                                                    nbSynCHR = 1, 
                                                                    methData = samplesForChrSynthetic, 
                                                                    nbBlock = 3, lBlock  = 2,
                                                                    vNbSample = 2, 
                                                                    nbGeneration = 3, 
                                                                    vpDiff = 2, vpDiffsd = 1, 
                                                                    vDiff = 2, 
                                                                    vInheritance = 2,
                                                                    propInherite = 0.8, 
                                                                    rateDiff = 2, 
                                                                    minRate = 1, 
                                                                    propHetero = 0.4, 
                                                                    minReads = 2, 
                                                                    maxPercReads = 99.9, 
                                                                    context = "CpG", assembly = "hg19",
                                                                    meanCov = 10, n = 3, 
                                                                    keepDiff = TRUE, saveGRanges = FALSE, 
                                                                    saveMethylKit = FALSE,
                                                                    anaMethylKit = FALSE,
                                                                    nbCores = 1, vSeed = -1),
                    error=conditionMessage)
    
    
    exp <- "fileID must be a character string or NULL"
    
    
    message <- paste0("test.validateRunSimParameters_fileID_number() ",
                      "- Number as fileID paramter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimParameters_nbSynCHR_not_number <- function() {
    
    obs <- tryCatch(methylInheritanceSim:::validateRunSimParameters(outputDir = "test",
                                                        fileID = "F1", 
                                                        nbSynCHR = "hi", 
                                                        methData = samplesForChrSynthetic, 
                                                        nbBlock = 3, lBlock  = 2,
                                                        vNbSample = 2, 
                                                        nbGeneration = 3, 
                                                        vpDiff = 2, vpDiffsd = 1, 
                                                        vDiff = 2, 
                                                        vInheritance = 2,
                                                        propInherite = 0.8, 
                                                        rateDiff = 2, 
                                                        minRate = 1, 
                                                        propHetero = 0.4, 
                                                        minReads = 2, 
                                                        maxPercReads = 99.9, 
                                                        context = "CpG", assembly = "hg19",
                                                        meanCov = 10, n = 3, 
                                                        keepDiff = TRUE, saveGRanges = FALSE, 
                                                        saveMethylKit = FALSE,
                                                        anaMethylKit = FALSE,
                                                        nbCores = 1, vSeed = -1),
                    error=conditionMessage)
    
    
    exp <- "nbSynCHR must be a positive integer or numeric"
    
    
    message <- paste0("test.validateRunSimParameters_nbSynCHR_not_number() ",
                      "- Not a number as nbSynCHR paramter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimParameters_nbSynCHR_vector_number <- function() {
    
    obs <- tryCatch(methylInheritanceSim:::validateRunSimParameters(outputDir = "test",
                                                                    fileID = "F1", 
                                                                    nbSynCHR = c(1,2), 
                                                                    methData = samplesForChrSynthetic, 
                                                                    nbBlock = 3, lBlock  = 2,
                                                                    vNbSample = 2, 
                                                                    nbGeneration = 3, 
                                                                    vpDiff = 2, vpDiffsd = 1, 
                                                                    vDiff = 2, 
                                                                    vInheritance = 2,
                                                                    propInherite = 0.8, 
                                                                    rateDiff = 2, 
                                                                    minRate = 1, 
                                                                    propHetero = 0.4, 
                                                                    minReads = 2, 
                                                                    maxPercReads = 99.9, 
                                                                    context = "CpG", assembly = "hg19",
                                                                    meanCov = 10, n = 3, 
                                                                    keepDiff = TRUE, saveGRanges = FALSE, 
                                                                    saveMethylKit = FALSE,
                                                                    anaMethylKit = FALSE,
                                                                    nbCores = 1, vSeed = -1),
                    error=conditionMessage)
    
    
    exp <- "nbSynCHR must be a positive integer or numeric"
    
    
    message <- paste0("test.validateRunSimParameters_nbSynCHR_vector_number() ",
                      "- Vector of numbers as nbSynCHR paramter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimParameters_nbSynCHR_zero <- function() {
    
    obs <- tryCatch(methylInheritanceSim:::validateRunSimParameters(outputDir = "test",
                                                                    fileID = "F1", 
                                                                    nbSynCHR = 0, 
                                                                    methData = samplesForChrSynthetic, 
                                                                    nbBlock = 3, lBlock  = 2,
                                                                    vNbSample = 2, 
                                                                    nbGeneration = 3, 
                                                                    vpDiff = 2, vpDiffsd = 1, 
                                                                    vDiff = 2, 
                                                                    vInheritance = 2,
                                                                    propInherite = 0.8, 
                                                                    rateDiff = 2, 
                                                                    minRate = 1, 
                                                                    propHetero = 0.4, 
                                                                    minReads = 2, 
                                                                    maxPercReads = 99.9, 
                                                                    context = "CpG", assembly = "hg19",
                                                                    meanCov = 10, n = 3, 
                                                                    keepDiff = TRUE, saveGRanges = FALSE, 
                                                                    saveMethylKit = FALSE,
                                                                    anaMethylKit = FALSE,
                                                                    nbCores = 1, vSeed = -1),
                    error=conditionMessage)
    
    
    exp <- "nbSynCHR must be a positive integer or numeric"
    
    
    message <- paste0("test.validateRunSimParameters_nbSynCHR_zero() ",
                      "- Zero as nbSynCHR paramter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}
