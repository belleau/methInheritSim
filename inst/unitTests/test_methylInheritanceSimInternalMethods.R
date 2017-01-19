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