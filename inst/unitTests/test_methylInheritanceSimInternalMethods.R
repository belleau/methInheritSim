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
## getSim() function
###################################################

test.getSim_good_01 <- function() {
    
    set.seed(312)
    
    stateInfo <- methylInheritanceSim:::getSyntheticChr(methInfo = 
                    samplesForChrSynthetic, nbBlock = 1, nbCpG = 2)
    
    stateDiff <- list()
    stateDiff[["stateDiff"]] <- c(0,1)
    stateDiff[["stateInherite"]] <- c(0,1)
    
    obs <- methylInheritanceSim:::getSim(2, 2, 2, stateInfo, stateDiff, 
                                            10, 0.8, 0.2, 0.8, 0.1)
    
    entry_01 <- GenomicRanges::GRanges(seqnames = rep("S", 2), 
                ranges = IRanges::IRanges(start = c(1000, 1007), 
                                        end = c(1000, 1007)),
                strand = rep("+", 2), meanDiff = c(0.987112859251413, 1),
                meanCTRL.meanCTRL = c(0.987112859251413, 0.0711020304345812), 
                partitionCase = c(0, 2), partitionCtrl = c(2,0),
                ctrl.V1 = c(0.984744151048604, 0.0751012481595866),
                ctrl.V2 = c(0.994061018191887, 0.0462071869386155),
                case.V1 = c(0.952275081359202, 0),
                case.V2 = c(0.998424919118535, 1))
    
    entry_02 <- GenomicRanges::GRanges(seqnames = rep("S", 2), 
                ranges = IRanges::IRanges(start = c(1000, 1007), 
                                end = c(1000, 1007)),
                strand = rep("+", 2), meanDiff = c(0.987112859251413, 1),
                meanCTRL.meanCTRL = c(0.987112859251413, 0.0711020304345812), 
                partitionCase = c(0, 1), partitionCtrl = c(2, 1),
                ctrl.V1 = c(0.996109673466832, 0.00572813753446487),
                ctrl.V2 = c(0.997316622965915, 0.0355629065128061),
                case.V1 = c(0.997846790687883, 1),
                case.V2 = c(0.974685083913519, 0.150238331165474))
    
    exp <- GenomicRanges::GRangesList(entry_01, entry_02)
    
    
    message <- paste0("test.getSim_good_01() ",
                      "- Valid paramters did not generated expected results.")
    
    checkEquals(obs, exp, message)
}