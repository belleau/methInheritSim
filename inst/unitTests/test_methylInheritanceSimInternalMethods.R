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

test.estBetaAlpha_good_02 <- function() {
    
    obs <- methylInheritanceSim:::estBetaAlpha(c(0.5, 0.2), 0.4)
    
    exp <- 0
    
    message <- paste0("test.estBetaAlpha_good_02() ",
                      "- Valid parameters did not generated expected results.")
    
    checkEquals(obs, exp, message)
}



###################################################
## estBetaAlphaNew() function
###################################################

test.estBetaAlphaNew_good_01 <- function() {
    
    obs <- methylInheritanceSim:::estBetaAlphaNew(0.5, 0.2, 0.00001)
    
    exp <- 0.125
    
    message <- paste0("test.estBetaAlphaNew_good_01() ",
                      "- Valid parameters did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.estBetaAlphaNew_good_02 <- function() {
    
    obs <- methylInheritanceSim:::estBetaAlphaNew(0.5, 0.2, 0.4)
    
    exp <- 0
    
    message <- paste0("test.estBetaAlphaNew_good_02() ",
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

test.estBetaBeta_good_02 <- function() {
    
    obs <- methylInheritanceSim:::estBetaBeta(c(0.3, 0.2), 0.35)
    
    exp <- 0
    
    message <- paste0("test.estBetaBeta_good_02() ",
                      "- Valid parameters did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.estBetaBeta_good_03 <- function() {
    
    obs <- methylInheritanceSim:::estBetaBeta(c(0.2, 0.002), 0.00001)
    
    exp <- 63.2
    
    message <- paste0("test.estBetaBeta_good_03() ",
                      "- Valid parameters did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

###################################################
## estBetaBetaNew() function
###################################################

test.estBetaBetaNew_good_01 <- function() {
    
    obs <- methylInheritanceSim:::estBetaBetaNew(meanCtrl=0.3, varCtrl=0.2, 0.00001)
    
    exp <- 0.035
    
    message <- paste0("test.estBetaBetaNew_good_01() ",
                      "- Valid parameters did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.estBetaBetaNew_good_02 <- function() {
    
    obs <- methylInheritanceSim:::estBetaBetaNew(meanCtrl=0.3, varCtrl=0.2, 0.35)
    
    exp <- 0
    
    message <- paste0("test.estBetaBetaNew_good_02() ",
                      "- Valid parameters did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.estBetaBeta_good_03 <- function() {
    
    obs <- methylInheritanceSim:::estBetaBetaNew(meanCtrl=0.2, varCtrl=0.002, 0.00001)
    
    exp <- 63.2
    
    message <- paste0("test.estBetaBeta_good_03() ",
                      "- Valid parameters did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

###################################################
## getDiffCaseNew() function
###################################################

test.getDiffCaseNew_good_01 <- function() {
    set.seed(322)
    
    obs <- methylInheritanceSim:::getDiffCaseNew(ctrlMean = 0.14562,
                ctrlVar = 0.0003607153, selectedAsDM = 1, nbCase = 4, sDiff = 0.8,
                nbDiffCase = 3)
    exp <- c(0.945620000000, 3.000000000000, 1.000000000000, 0.947694615429,
             0.965193968711, 0.906084052941, 0.122224066759)
    
    message <- paste0("test.getDiffCaseNew_good_01() ",
                      "- Valid parameters did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.getDiffCaseNew_good_no_DMS <- function() {
    set.seed(22)
    
    obs <- methylInheritanceSim:::getDiffCaseNew(ctrlMean = 0.14562,
                ctrlVar = 0.0003607153, selectedAsDM = 0, nbCase = 4, sDiff = 0.8,
                nbDiffCase = 3)
    exp <- c(0.1456200000, 0.000000000000, 4.000000000000, 0.1348148850,
             0.1698118247, 0.1520623197, 0.1411707002)
    
    message <- paste0("test.getDiffCaseNew_good_no_DMS() ",
                      "- Valid parameters with no DMS did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.getDiffCaseNew_good_02 <- function() {
    
    set.seed(322)
    
    obs <- methylInheritanceSim:::getDiffCaseNew(ctrlMean = 0.14562,
                ctrlVar = 0.0003607153, selectedAsDM = 1, nbCase = 4, sDiff = 0.1,
                nbDiffCase = 3)
    exp <- c(0.2456200000, 3.000000000000, 1.000000000000, 0.2435759950,
             0.2230879838, 0.2761777813, 0.1222240668)
    
    message <- paste0("test.getDiffCaseNew_good_02() ",
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

test.getDiffCase_good_02 <- function() {
    
    set.seed(322)
    
    x <- c(0.14562, 0.0003607153, 1)
    obs <- methylInheritanceSim:::getDiffCase(x = x, nb = 4, sDiff = 0.1,
                    diffCase = 3)
    exp <- c(0.2456200000, 3.000000000000, 1.000000000000, 0.2435759950,
             0.2230879838, 0.2761777813, 0.1222240668)
    
    message <- paste0("test.getDiffCase_good_02() ",
                      "- Valid parameters did not generated expected results.")
    
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


test.getSim_good_02 <- function() {
    
    set.seed(22212)
    
    stateInformation <- methylInheritanceSim:::getSyntheticChr(methInfo = samplesForChrSynthetic, 
                                                               nbBlock = 1, nbCpG = 3)
    
    stateDiff <- list()
    stateDiff[["stateDiff"]] <- c(1, 1, 1)
    stateDiff[["stateInherite"]] <- c(0, 0, 1)
    
    obs <- methylInheritanceSim:::getSim(nbCtrl = 2, nbCase = 3, generation = 2, 
                                         stateInfo = stateInformation, stateDiff = stateDiff, diffValue = 10, 
                                         propDiff = 0.6, propDiffsd = 0.3, propInheritance = 0.7, propHetero = 0.2)
    
    exp <- GRangesList()
    
    exp[[1]] <- GRanges(seqnames = seqnames(stateInformation),
                        ranges = ranges(stateInformation),
                        strand =  strand(stateInformation),
                        meanDiff = c(0.000000000000000, 1.000000000000000, 1.000000000000000), 
                        meanCTRL = mcols(stateInformation)[3],
                        partitionCase = c(2, 2, 2), partitionCtrl = c(1, 1, 1),
                        ctrl.V1 = c(0.599248212271003, 0.104673052925196, 0.0163259007084415),
                        ctrl.V2 = c(0.660878027594768, 0.0737704350938676, 0.00412068310774805),
                        case.V1 = c(0.000000000000000, 1.000000000000000, 1.000000000000000),
                        case.V2 = c(0.000000000000000, 1.000000000000000, 1.000000000000000),
                        case.V3 = c(0.140108118897419, 0.152428811636903, 0.0242355851971473))
    
    exp[[2]] <- GRanges(seqnames = seqnames(stateInformation),
                        ranges = ranges(stateInformation),
                        strand =  strand(stateInformation),
                        meanDiff = c(0.534404290296076, 0.288519632432523, 1.000000000000000), 
                        meanCTRL = mcols(stateInformation)[3],
                        partitionCase = c(0, 0, 1), partitionCtrl = c(3, 3, 2),
                        ctrl.V1 = c(0.0851227837846553, 0.281923948948365, 0.0198953426761693),
                        ctrl.V2 = c(0.393339608860298, 0.118756249033386, 0.019181448769631),
                        case.V1 = c(0.960713455049359, 0.140418843034454, 1.000000000000000),
                        case.V2 = c(0.935659881913647, 0.282094425404783, 0.0137838586380334),
                        case.V3 = c(0.495670278966184, 0.476562011603015, 0.00593249270592643))
    
    message <- paste0("test.getSim_good_02() ",
                      "- Valid parameters for getSim() did not generated expected results.")
    
    checkEquals(obs, exp, message)
}


test.getSim_good_03 <- function() {
    
    set.seed(25212)
    
    stateInformation <- methylInheritanceSim:::getSyntheticChr(methInfo = samplesForChrSynthetic, 
                            nbBlock = 1, nbCpG = 3)
    
    stateDiff <- list()
    stateDiff[["stateDiff"]] <- c(1, 1, 1)
    stateDiff[["stateInherite"]] <- c(1, 1, 1)
    
    obs <- methylInheritanceSim:::getSim(nbCtrl = 3, nbCase = 3, generation = 2, 
                stateInfo = stateInformation, stateDiff = stateDiff, diffValue = 0.4, 
                propDiff = 0.8, propDiffsd = 0.3, propInheritance = 0.7, propHetero = 0.2)
    
    exp <- GRangesList()
    
    exp[[1]] <- GRanges(seqnames = seqnames(stateInformation),
                        ranges = ranges(stateInformation),
                        strand =  strand(stateInformation),
                        meanDiff = c(0.558913110086267, 0.449962615641849, 0.479661899006178), 
                        meanCTRL = mcols(stateInformation)[3],
                        partitionCase = c(2, 2, 2), partitionCtrl = c(1, 1, 1),
                        ctrl.V1 = c(0.999890625604128, 0.0224875348296994, 0.145194054662506),
                        ctrl.V2 = c(0.972502722616637, 0.186942812560022, 0.0673588416556606),
                        ctrl.V3 = c(0.987586164771293, 0.020569016934098, 0.131640777446695),
                        case.V1 = c(0.421402396211426, 0.448889310453412, 0.473769208188253),
                        case.V2 = c(0.447204191206416, 0.448558536584455, 0.543957006411924),
                        case.V3 = c(0.999997880496921, 0.0584002453973412, 0.0385719160422928))
    
    exp[[2]] <- GRanges(seqnames = seqnames(stateInformation),
                        ranges = ranges(stateInformation),
                        strand =  strand(stateInformation),
                        meanDiff = c(0.878913110086267, 0.129962615641849, 0.159661899006178), 
                        meanCTRL = mcols(stateInformation)[3],
                        partitionCase = c(2, 2, 2), partitionCtrl = c(1, 1, 1),
                        ctrl.V1 = c(0.99998519489511, 0.0315260421172793, 0.079646895269669),
                        ctrl.V2 = c(0.986784892583002, 0.0131363053820546, 0.0409011080082275),
                        ctrl.V3 = c(0.998699787293772, 0.0549282859188099, 0.12809374105219),
                        case.V1 = c(0.749498444969622, 0.146170576424504, 0.16625757987502),
                        case.V2 = c(0.935273128088281, 0.204586530901366, 0.137162736236601),
                        case.V3 = c(0.958785951984364, 0.0640823013094492, 0.122777372842008))
    
    message <- paste0("test.getSim_good_03() ",
                      "- Valid parameters for getSim() did not generated expected results.")
    
    checkEquals(obs, exp, message)
}


###################################################
## getSimNew() function
###################################################

test.getSimNew_good_01 <- function() {
    
    set.seed(22212)
    
    stateInformation <- methylInheritanceSim:::getSyntheticChr(methInfo = samplesForChrSynthetic, 
                                                    nbBlock = 1, nbCpG = 3)
    
    stateDiff <- c(1, 0, 1)
    stateInherite <- c(1, 0, 0)
    
    obs <- methylInheritanceSim:::getSimNew(nbCtrl = 2, nbCase = 3, generation = 2, 
                                        stateInfo = stateInformation, stateDiff = stateDiff, 
                                        stateInherite = stateInherite, diffValue = 10, 
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
    
    message <- paste0("test.getSimNew_good_01() ",
                      "- Valid parameters for getSim() did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.getSimNew_good_02 <- function() {
    
    set.seed(22212)
    
    stateInformation <- methylInheritanceSim:::getSyntheticChr(methInfo = samplesForChrSynthetic, 
                                                               nbBlock = 1, nbCpG = 3)
    
    stateDiff <- c(1, 1, 1)
    stateInherite <- c(0, 0, 1)
    
    obs <- methylInheritanceSim:::getSimNew(nbCtrl = 2, nbCase = 3, generation = 2, 
                                    stateInfo = stateInformation, stateDiff = stateDiff, 
                                    stateInherite = stateInherite, diffValue = 10, 
                                    propDiff = 0.6, propDiffsd = 0.3, propInheritance = 0.7, propHetero = 0.2)
    
    exp <- GRangesList()
    
    exp[[1]] <- GRanges(seqnames = seqnames(stateInformation),
                        ranges = ranges(stateInformation),
                        strand =  strand(stateInformation),
                        meanDiff = c(0.000000000000000, 1.000000000000000, 1.000000000000000), 
                        meanCTRL = mcols(stateInformation)[3],
                        partitionCase = c(2, 2, 2), partitionCtrl = c(1, 1, 1),
                        ctrl.V1 = c(0.599248212271003, 0.104673052925196, 0.0163259007084415),
                        ctrl.V2 = c(0.660878027594768, 0.0737704350938676, 0.00412068310774805),
                        case.V1 = c(0.000000000000000, 1.000000000000000, 1.000000000000000),
                        case.V2 = c(0.000000000000000, 1.000000000000000, 1.000000000000000),
                        case.V3 = c(0.140108118897419, 0.152428811636903, 0.0242355851971473))
    
    exp[[2]] <- GRanges(seqnames = seqnames(stateInformation),
                        ranges = ranges(stateInformation),
                        strand =  strand(stateInformation),
                        meanDiff = c(0.534404290296076, 0.288519632432523, 1.000000000000000), 
                        meanCTRL = mcols(stateInformation)[3],
                        partitionCase = c(0, 0, 1), partitionCtrl = c(3, 3, 2),
                        ctrl.V1 = c(0.0851227837846553, 0.281923948948365, 0.0198953426761693),
                        ctrl.V2 = c(0.393339608860298, 0.118756249033386, 0.019181448769631),
                        case.V1 = c(0.960713455049359, 0.140418843034454, 1.000000000000000),
                        case.V2 = c(0.935659881913647, 0.282094425404783, 0.0137838586380334),
                        case.V3 = c(0.495670278966184, 0.476562011603015, 0.00593249270592643))
    
    message <- paste0("test.getSimNew_good_02() ",
                      "- Valid parameters for getSim() did not generated expected results.")
    
    checkEquals(obs, exp, message)
}


test.getSimNew_good_03 <- function() {
    
    set.seed(25212)
    
    stateInformation <- methylInheritanceSim:::getSyntheticChr(methInfo = samplesForChrSynthetic, 
                                                               nbBlock = 1, nbCpG = 3)
    
    stateDiff <- list()
    stateDiff <- c(1, 1, 1)
    stateInherite <- c(1, 1, 1)
    
    obs <- methylInheritanceSim:::getSimNew(nbCtrl = 3, nbCase = 3, generation = 2, 
                stateInfo = stateInformation, stateDiff = stateDiff, stateInherite = stateInherite, diffValue = 0.4, 
                propDiff = 0.8, propDiffsd = 0.3, propInheritance = 0.7, propHetero = 0.2)
    
    exp <- GRangesList()
    
    exp[[1]] <- GRanges(seqnames = seqnames(stateInformation),
                        ranges = ranges(stateInformation),
                        strand =  strand(stateInformation),
                        meanDiff = c(0.558913110086267, 0.449962615641849, 0.479661899006178), 
                        meanCTRL = mcols(stateInformation)[3],
                        partitionCase = c(2, 2, 2), partitionCtrl = c(1, 1, 1),
                        ctrl.V1 = c(0.999890625604128, 0.0224875348296994, 0.145194054662506),
                        ctrl.V2 = c(0.972502722616637, 0.186942812560022, 0.0673588416556606),
                        ctrl.V3 = c(0.987586164771293, 0.020569016934098, 0.131640777446695),
                        case.V1 = c(0.421402396211426, 0.448889310453412, 0.473769208188253),
                        case.V2 = c(0.447204191206416, 0.448558536584455, 0.543957006411924),
                        case.V3 = c(0.999997880496921, 0.0584002453973412, 0.0385719160422928))
    
    exp[[2]] <- GRanges(seqnames = seqnames(stateInformation),
                        ranges = ranges(stateInformation),
                        strand =  strand(stateInformation),
                        meanDiff = c(0.878913110086267, 0.129962615641849, 0.159661899006178), 
                        meanCTRL = mcols(stateInformation)[3],
                        partitionCase = c(2, 2, 2), partitionCtrl = c(1, 1, 1),
                        ctrl.V1 = c(0.99998519489511, 0.0315260421172793, 0.079646895269669),
                        ctrl.V2 = c(0.986784892583002, 0.0131363053820546, 0.0409011080082275),
                        ctrl.V3 = c(0.998699787293772, 0.0549282859188099, 0.12809374105219),
                        case.V1 = c(0.749498444969622, 0.146170576424504, 0.16625757987502),
                        case.V2 = c(0.935273128088281, 0.204586530901366, 0.137162736236601),
                        case.V3 = c(0.958785951984364, 0.0640823013094492, 0.122777372842008))
    
    message <- paste0("test.getSimNew_good_03() ",
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

test.validateRunSimOtherParameters_good_01 <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimOtherParameters(outputDir = "test",
                                fileID = "F1", methData = samplesForChrSynthetic,
                                context = "hi", assembly = "hg19"),
        error=conditionMessage)
    
    exp <- 0
    
    message <- paste0("test.validateRunSimOtherParameters_good_01() ",
                      "- All good parameters did not generated expected results.")
    
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

test.validateRunSimLogicalParameters_good_01 <- function() {
    obs <- tryCatch(
        methylInheritanceSim:::validateRunSimLogicalParameters(keepDiff = TRUE, 
            saveGRanges = FALSE, saveMethylKit = FALSE, runAnalysis = FALSE),
        error=conditionMessage)
    
    exp <- 0
    
    message <- paste0("test.validateRunSimLogicalParameters_good_01() ",
                      "- All good parameters did not generated expected results.")
    
    checkEquals(obs, exp, message)
}


###################################################
## validateRunSimIntegerParameters() function
###################################################

test.validateRunSimIntegerParameters_nbSynCHR_not_number <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimIntegerParameters(nbSynCHR = "hi", 
            nbSimulation = 10, nbBlock = 3, nbCpG  = 2, vNbSample = 2, 
            nbGeneration = 3, minReads = 2, meanCov = 10, nbCores = 1, vSeed = -1),
        error=conditionMessage)
    
    exp <- "nbSynCHR must be a positive integer or numeric"
    
    message <- paste0("test.validateRunSimIntegerParameters_nbSynCHR_not_number() ",
                      "- Not a number as nbSynCHR parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}


test.validateRunSimIntegerParameters_nbSynCHR_vector_number <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimIntegerParameters(nbSynCHR = c(1,2), 
            nbSimulation = 10, nbBlock = 3, nbCpG  = 2,
                vNbSample = 2, nbGeneration = 3, minReads = 2, meanCov = 10, 
                nbCores = 1, vSeed = -1),
            error=conditionMessage)
    
    exp <- "nbSynCHR must be a positive integer or numeric"
    
    message <- paste0("test.validateRunSimIntegerParametersnbSynCHR_vector_number() ",
                      "- Vector of numbers as nbSynCHR parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimIntegerParameters_nbSynCHR_zero <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimIntegerParameters(nbSynCHR = 0, 
                nbSimulation = 10, nbBlock = 3, nbCpG  = 2, vNbSample = 2, nbGeneration = 3, 
                meanCov = 10, nbCores = 1, vSeed = -1),
            error=conditionMessage)
    
    exp <- "nbSynCHR must be a positive integer or numeric"
    
    message <- paste0("test.validateRunSimIntegerParameters_nbSynCHR_zero() ",
                      "- Zero as nbSynCHR parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimIntegerParameters_nbSimulation_not_number <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimIntegerParameters(nbSynCHR = 1, 
                nbSimulation = "Hi", nbBlock = 3, nbCpG  = 2, vNbSample = 2, 
                nbGeneration = 3, minReads = 2, maxPercReads = 99.9, meanCov = 10, 
                nbCores = 1, vSeed = -1),
            error=conditionMessage)
    
    exp <- "nbSimulation must be a positive integer or numeric"
    
    message <- paste0("test.validateRunSimIntegerParameters_nbSimulation_not_number() ",
                    "- Not a number as nbSimulation parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}


test.validateRunSimIntegerParameters_nbSimulation_not_number <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimIntegerParameters(nbSynCHR = 1, 
                nbSimulation = "Hi", nbBlock = 3, nbCpG  = 2, vNbSample = 2, 
                nbGeneration = 3, minReads = 2, meanCov = 10, nbCores = 1, vSeed = -1),
            error=conditionMessage)
    
    exp <- "nbSimulation must be a positive integer or numeric"
    
    message <- paste0("test.validateRunSimIntegerParameters_nbSimulation_not_number() ",
                    "- Not a number as nbSimulation parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimIntegerParameters_nbSimulation_vector_number <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimIntegerParameters(nbSynCHR = 1, 
                        nbSimulation = c(1,2), nbBlock = 3, nbCpG  = 2, vNbSample = 2, 
                        nbGeneration = 3,  minReads = 2, 
                        meanCov = 10, nbCores = 1, vSeed = -1),
            error=conditionMessage)
    
    exp <- "nbSimulation must be a positive integer or numeric"
    
    message <- paste0("test.validateRunSimIntegerParameters_nbSimulation_vector_number() ",
                "- Vector of numbers as nbSimulation parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimIntegerParameters_nbSimulation_zero <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimIntegerParameters(nbSynCHR = 1, 
                    nbSimulation = 0, nbBlock = 3, nbCpG = 2, vNbSample = 2, 
                    nbGeneration = 3, minReads = 2, 
                    meanCov = 10, nbCores = 1, vSeed = -1),
            error=conditionMessage)
    
    exp <- "nbSimulation must be a positive integer or numeric"
    
    message <- paste0("test.validateRunSimIntegerParameters_nbSimulation_zero() ",
                "- Zero as nbSimulation parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimIntegerParameters_nbBlock_number_vector <- function() {
    obs <- tryCatch(
        methylInheritanceSim:::validateRunSimIntegerParameters(nbSynCHR = 1, 
                    nbSimulation = 10, nbBlock = c(3, 1), nbCpG  = 2, vNbSample = 2, 
                    nbGeneration = 3, minReads = 2,
                    meanCov = 10, nbCores = 1, vSeed = -1),
        error=conditionMessage)
    
    exp <- "nbBlock must be a positive integer or numeric"
    
    message <- paste0("test.validateRunSimIntegerParameters_nbBlock_number_vector() ",
                    "- Number vector as nbBlock parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimIntegerParameters_nbBlock_string <- function() {
    obs <- tryCatch(
        methylInheritanceSim:::validateRunSimIntegerParameters(nbSynCHR = 1, 
                    nbSimulation = 10, nbBlock = "hi", nbCpG  = 2, vNbSample = 2, 
                    nbGeneration = 3, minReads = 2, meanCov = 10, 
                    nbCores = 1, vSeed = -1),
        error=conditionMessage)
    
    exp <- "nbBlock must be a positive integer or numeric"
    
    message <- paste0("test.validateRunSimIntegerParameters_nbBlock_string() ",
                      "- String as nbBlock parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimIntegerParameters_nbBlock_zero <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimIntegerParameters(nbSynCHR = 1, 
                        nbSimulation = 10, nbBlock = 0, nbCpG  = 2, vNbSample = 2, 
                        nbGeneration = 3, minReads = 2, meanCov = 10,
                        nbCores = 1, vSeed = -1),
                    error=conditionMessage)
    
    exp <- "nbBlock must be a positive integer or numeric"
    
    message <- paste0("test.validateRunSimIntegerParameters_nbBlock_zero() ",
                      "- Zero as nbBlock parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimIntegerParameters_nbCpG_number_vector <- function() {
    obs <- tryCatch(
        methylInheritanceSim:::validateRunSimIntegerParameters(nbSynCHR = 1, 
                        nbSimulation = 10, nbBlock = 2, nbCpG = c(3,2), vNbSample = 2, 
                        nbGeneration = 3, minReads = 2, meanCov = 10, 
                        nbCores = 1, vSeed = -1),
        error=conditionMessage)
    
    exp <- "nbCpG must be a positive integer or numeric"
    
    message <- paste0("test.validateRunSimIntegerParameters_nbCpG_number_vector() ",
                      "- Number vector as nbCpG parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimIntegerParameters_nbCpG_string <- function() {
    obs <- tryCatch(
        methylInheritanceSim:::validateRunSimIntegerParameters(nbSynCHR = 1, 
                nbSimulation = 10, nbBlock = 2, nbCpG  = "hi", vNbSample = 2, 
                nbGeneration = 3, minReads = 2,
                meanCov = 10, nbCores = 1, vSeed = -1),
        error=conditionMessage)
    
    exp <- "nbCpG must be a positive integer or numeric"
    
    message <- paste0("test.validateRunSimIntegerParameters_nbCpG_string() ",
                      "- String as nbCpG parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimIntegerParameters_nbCpG_zero <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimIntegerParameters(nbSynCHR = 1, 
                nbSimulation = 10, nbBlock = 2, nbCpG  = 0, vNbSample = 2, 
                nbGeneration = 3, minReads = 2, meanCov = 10,
                nbCores = 1, vSeed = -1),
            error=conditionMessage)
    
    exp <- "nbCpG must be a positive integer or numeric"
    
    message <- paste0("test.validateRunSimIntegerParameters_nbCpG_zero() ",
                      "- Zero as nbCpG parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimIntegerParameters_nbGeneration_number_vector <- function() {
    obs <- tryCatch(
        methylInheritanceSim:::validateRunSimIntegerParameters(nbSynCHR = 1, 
            nbSimulation = 10, nbBlock = 2, nbCpG  = 2, nbGeneration = c(3,2), 
            vNbSample = 2, minReads = 2, meanCov = 10, 
            nbCores = 1, vSeed = -1),
        error=conditionMessage)
    
    exp <- "nbGeneration must be a positive integer or numeric"
    
    message <- paste0("test.validateRunSimIntegerParameters_nbGeneration_number_vector() ",
                      "- Number vector as nbGeneration parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimIntegerParameters_nbGeneration_string <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimIntegerParameters(nbSynCHR = 1, 
            nbSimulation = 10, nbBlock = 2, nbCpG  = 2, nbGeneration = "Hi",
            vNbSample = 2, minReads = 2, meanCov = 10, 
            nbCores = 1, vSeed = -1),
        error=conditionMessage)
    
    exp <- "nbGeneration must be a positive integer or numeric"
    
    message <- paste0("test.validateRunSimIntegerParameters_nbGeneration_string() ",
                      "- String as nbGeneration parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimIntegerParameters_nbGeneration_zero <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimIntegerParameters(nbSynCHR = 1, 
                nbSimulation = 10, nbBlock = 2, nbCpG  = 2, nbGeneration = 0,
                vNbSample = 2, minReads = 2, meanCov = 10,
                nbCores = 1, vSeed = -1),
            error=conditionMessage)
    
    exp <- "nbGeneration must be a positive integer or numeric"
    
    message <- paste0("test.validateRunSimIntegerParameters_nbGeneration_zero() ",
                      "- Zero as nbGeneration parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimIntegerParameters_vNbSample_string <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimIntegerParameters(nbSynCHR = 1, 
                nbSimulation = 10, nbBlock = 2, nbCpG  = 2, nbGeneration = 3,
                vNbSample = "Hi", minReads = 2, meanCov = 10,
                nbCores = 1, vSeed = -1),
            error=conditionMessage)
    
    exp <- "vNbSample must be a vector of distinct positive integer"
    
    message <- paste0("test.validateRunSimIntegerParameterss_vNbSample_string() ",
                      "- String as vNbSample parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimIntegerParametersvNbSample_duplicate <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimIntegerParameters(nbSynCHR = 1, 
                nbSimulation = 10, nbBlock = 2, nbCpG  = 2, nbGeneration = 3,
                vNbSample = c(1, 2, 2), minReads = 2, meanCov = 10,
                nbCores = 1, vSeed = -1),
            error=conditionMessage)
    
    exp <- "vNbSample must be a vector of distinct positive integer"
    
    message <- paste0("test.validateRunSimIntegerParameters_vNbSample_duplicate() ",
                      "- Duplicate elements in vNbSample parameter did not 
                      generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimIntegerParameters_vNbSample_float <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimIntegerParameters(nbSynCHR = 1, 
                nbSimulation = 10, nbBlock = 2, nbCpG  = 2, nbGeneration = 3,
                vNbSample = c(1, 2, 2.1), minReads = 2, 
                meanCov = 10, nbCores = 1, vSeed = -1),
            error=conditionMessage)
    
    exp <- "vNbSample must be a vector of distinct positive integer"
    
    message <- paste0("test.validateRunSimIntegerParameters_vNbSample_float() ",
                      "- Float as vNbSample parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimIntegerParameters_minReads_zero <- function() {
    obs <- tryCatch(
        methylInheritanceSim:::validateRunSimIntegerParameters(nbSynCHR = 1, 
                    nbBlock = 2, nbCpG = 2, vNbSample = 2, nbGeneration = 3, 
                    minReads = 0, meanCov = 10, n = 3, nbCores = 1, vSeed = -1),
        error=conditionMessage)
    
    exp <- "minReads must be a positive integer or numeric"
    
    message <- paste0("test.validateRunSimIntegerParameters_minReads_zero() ",
                      "- Zero as minReads parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}



test.validateRunSimIntegerParameters_nbCores_string <- function() {
    obs <- tryCatch(
        methylInheritanceSim:::validateRunSimIntegerParameters(nbSynCHR = 1, 
                nbBlock = 2, nbCpG  = 2, vNbSample = 2, nbGeneration = 3,  
                minReads = 2, meanCov = 10, n = 3, nbCores = "hi", vSeed = -1),
        error=conditionMessage)
    
    exp <- "nbCores must be a positive integer or numeric"
    
    message <- paste0("test.validateRunSimIntegerParameters_nbCores_string() ",
                      "- String as nbCores parameter did not generated expected
                      results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimIntegerParameters_nbCores_zero <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimIntegerParameters(nbSynCHR = 1, 
                    nbBlock = 2, nbCpG = 2, vNbSample = 2, nbGeneration = 3, 
                    minReads = 2, meanCov = 10, n = 3, nbCores = 0, vSeed = -1),
        error=conditionMessage)
    
    exp <- "nbCores must be a positive integer or numeric"
    
    message <- paste0("test.validateRunSimIntegerParameters_nbCores_zero() ",
                      "- Zero as nbCores parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimIntegerParameters_minReads_vector <- function() {
    obs <- tryCatch(
        methylInheritanceSim:::validateRunSimIntegerParameters(nbSynCHR = 1, 
                nbBlock = 2, nbCpG = 2, vNbSample = 2, nbGeneration = 3, 
                minReads = c(1, 3), meanCov = 10, n = 3, nbCores = 1, vSeed = -1),
        error=conditionMessage)
    
    exp <- "minReads must be a positive integer or numeric"
    
    message <- paste0("test.validateRunSimIntegerParameters_minReads_vector() ",
                      "- Vector as minReads parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}


test.validateRunSimIntegerParameters_meanCov_zero <- function() {
    obs <- tryCatch(
        methylInheritanceSim:::validateRunSimIntegerParameters(nbSynCHR = 1, 
            nbBlock = 2, nbCpG = 2, vNbSample = 2, nbGeneration = 3, 
            minReads = 4, meanCov = 0, n = 3, nbCores = 1, vSeed = -1),
        error=conditionMessage)
    
    exp <- "meanCov must be a positive integer or numeric"
    
    message <- paste0("test.validateRunSimIntegerParameters_meanCov_zero() ",
                      "- Zero as meanCov parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimIntegerParameters_vSeed_string <- function() {
    obs <- tryCatch(
        methylInheritanceSim:::validateRunSimIntegerParameters(nbSynCHR = 1, 
            nbBlock = 2, nbCpG = 2, vNbSample = 2, nbGeneration = 3, 
            minReads = 4, meanCov = 10, n = 3, nbCores = 1, vSeed = "test"),
        error=conditionMessage)
    
    exp <- "vSeed must be an integer or numeric"
    
    message <- paste0("test.validateRunSimIntegerParameters_vSeed_string() ",
                      "- String as vSeed parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimIntegerParameters_meanCov_string <- function() {
    obs <- tryCatch(
        methylInheritanceSim:::validateRunSimIntegerParameters(nbSynCHR = 1, 
                nbBlock = 2, nbCpG = 2, vNbSample = 2, nbGeneration = 3, 
                minReads = 4, meanCov = "allo", n = 3, nbCores = 1, vSeed = -1),
        error=conditionMessage)
    
    exp <- "meanCov must be a positive integer or numeric"
    
    message <- paste0("test.validateRunSimIntegerParameters_meanCov_string() ",
                      "- String as meanCov parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimIntegerParameters_good_01 <- function() {
    obs <- methylInheritanceSim:::validateRunSimIntegerParameters(nbSynCHR = 1, 
                nbBlock = 2, nbCpG = 2, vNbSample = 2, nbGeneration = 3, 
                minReads = 4, meanCov = 10, n = 3, nbCores = 1, vSeed = -1)
    
    exp <- 0
    
    message <- paste0("test.validateRunSimIntegerParameters_good_01() ",
                      "- All valid parameters did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

###################################################
## validateRunSimDoubleParameters() function
###################################################

test.validateRunSimDoubleParameters_vpDiff_string <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimDoubleParameters(
                vpDiff = "Hi", vpDiffsd = 0.1, vDiff = 0.8, 
                vInheritance = 0.5, propInherite = 0.8, rateDiff = 2, minRate = 1, 
                propHetero = 0.4, maxPercReads = 99.9),
        error=conditionMessage)
    
    exp <- "vpDiff must be a vector of distinct positive double include in (0,1]"
    
    message <- paste0("test.validateRunSimDoubleParameters_vpDiff_string() ",
                    "- String as vpDiff parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimDoubleParameters_vpDiff_duplicate <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimDoubleParameters(
                vpDiff = c(0.5, 0.5, 0.8), vpDiffsd = 0.1, 
                vDiff = 0.9, vInheritance = 0.5, propInherite = 0.8, rateDiff = 2, 
                minRate = 1, propHetero = 0.4, maxPercReads = 99.9),
            error=conditionMessage)
    
    exp <- "vpDiff must be a vector of distinct positive double include in (0,1]"
    
    message <- paste0("test.validateRunSimDoubleParameters_vpDiff_duplicate() ",
                      "- Duplicate elements in vpDiff parameter did not 
                      generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimDoubleParameters_vpDiff_zero <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimDoubleParameters(
                vpDiff = c(0.5, 0, 0.8), vpDiffsd = 0.1, vDiff = 0.8, 
                vInheritance = 0.5,propInherite = 0.8, rateDiff = 2, minRate = 1, 
                propHetero = 0.4, maxPercReads = 99.9),
            error=conditionMessage)
    
    exp <- "vpDiff must be a vector of distinct positive double include in (0,1]"
    
    message <- paste0("test.validateRunSimDoubleParameters_vpDiff_zero() ",
                      "- Zero in vpDiff parameter did not 
                      generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimDoubleParameters_vpDiff_two <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimDoubleParameters(
                vpDiff = c(0.5, 2, 0.8), vpDiffsd = 0.1, vDiff = 0.8, 
                vInheritance = 0.5, propInherite = 0.8, rateDiff = 2, minRate = 1, 
                propHetero = 0.4, maxPercReads = 99.9),
            error=conditionMessage)
    
    exp <- "vpDiff must be a vector of distinct positive double include in (0,1]"
    
    message <- paste0("test.validateRunSimDoubleParameters_vpDiff_two() ",
                      "- Two in vpDiff parameter did not 
                      generated expected results.")
    
    checkEquals(obs, exp, message)
}


test.validateRunSimDoubleParameters_vpDiffsd_string <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimDoubleParameters(
                vpDiff = 0.9, vpDiffsd = "Hi", vDiff = 0.8, 
                vInheritance = 0.5, propInherite = 0.8, rateDiff = 2, 
                minRate = 1, propHetero = 0.4, maxPercReads = 99.9),
            error=conditionMessage)
    
    exp <- "vpDiffsd must be a vector of non-negative double"
    
    message <- paste0("test.validateRunSimDoubleParameters_vpDiffsd_string() ",
                      "- String as vpDiffsd parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimDoubleParameters_vpDiffsd_negative <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimDoubleParameters(
                vpDiff = 0.9, vpDiffsd = -1, vDiff = 0.8, 
                vInheritance = 0.5, propInherite = 0.8, rateDiff = 2, minRate = 1, 
                propHetero = 0.4, maxPercReads = 99.9),
            error=conditionMessage)
    
    exp <- "vpDiffsd must be a vector of non-negative double"
    
    message <- paste0("test.validateRunSimDoubleParameters_vpDiffsd_negative() ",
                      "- Negative as vpDiffsd parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimDoubleParameters_vpDiffsd_vpDiff <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimDoubleParameters(
                vpDiff = 0.9, vpDiffsd = c(0.1, 0.1), vDiff = 0.8, 
                vInheritance = 0.5, propInherite = 0.8, rateDiff = 2, minRate = 1, 
                propHetero = 0.4, maxPercReads = 99.9),
            error=conditionMessage)
    
    exp <- "vpDiff and vpDiffsd must be the same length"
    
    message <- paste0("test.validateRunSimDoubleParameters_vpDiffsd_vpDiff() ",
                      "- vpDiffsd vpDiff not same length did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimDoubleParameters_vDiff_string <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimDoubleParameters(
                vpDiff = 0.9, vpDiffsd = 0.1,  vDiff = "Hi", 
                vInheritance = 0.5, propInherite = 0.8, rateDiff = 2, minRate = 1, 
                propHetero = 0.4, maxPercReads = 99.9),
            error=conditionMessage)
    
    exp <- "vDiff must be a vector of distinct non-negative double include in [0,1]"
    
    message <- paste0("test.validateRunSimDoubleParameters_vDiff_string() ",
                      "- String as vDiff parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimDoubleParameters_vDiff_duplicate <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimDoubleParameters(
            vpDiff = c(0.5), vpDiffsd = 0.1, vDiff = c(0.5, 0.5, 0.8), 
            vInheritance = 0.5, propInherite = 0.8, rateDiff = 2, minRate = 1, 
            propHetero = 0.4, maxPercReads = 99.9),
        error=conditionMessage)
    
    exp <- "vDiff must be a vector of distinct non-negative double include in [0,1]"
    
    message <- paste0("test.validateRunSimDoubleParameters_vDiff_duplicate() ",
                    "- Duplicate elements in vDiff parameter did not ", 
                    "generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimDoubleParameters_vDiff_negative <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimDoubleParameters(
                vpDiff = c(0.5), vpDiffsd = 0.1, vDiff = -1, 
                vInheritance = 0.5, propInherite = 0.8, rateDiff = 2, minRate = 1, 
                propHetero = 0.4, maxPercReads = 99.9),
            error=conditionMessage)
    
    exp <- "vDiff must be a vector of distinct non-negative double include in [0,1]"
    
    message <- paste0("test.validateRunSimDoubleParameters_vDiff_negative() ",
                      "- Negative in vDiff parameter did not ",
                      "generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimDoubleParameters_vDiff_two <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimDoubleParameters(
                    vpDiff = c(0.5), vpDiffsd = 0.1, vDiff = 2, 
                    vInheritance = 0.5, propInherite = 0.8, rateDiff = 2, 
                    minRate = 1, propHetero = 0.4, maxPercReads = 99.9),
            error=conditionMessage)
    
    exp <- "vDiff must be a vector of distinct non-negative double include in [0,1]"
    
    message <- paste0("test.validateRunSimDoubleParameters_vDiff_two() ",
                      "- Two in vDiff parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimDoubleParameters_vInheritance_string <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimDoubleParameters(
                    vpDiff = 0.9, vpDiffsd = 0.1, vDiff = 0.8, 
                    vInheritance = "Hi", propInherite = 0.8, rateDiff = 2, 
                    minRate = 1, propHetero = 0.4, maxPercReads = 99.9),
            error=conditionMessage)
    
    exp <- paste0("vInheritance must be a vector of distinct non-negative ", 
                  "double include in [0,1]")
    
    message <- paste0("test.validateRunSimDoubleParameters_vInheritance_string() ",
                      "- String as vInheritance parameter did not generated ", 
                      "expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimDoubleParameters_vInheritance_duplicate <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimDoubleParameters(
                        vpDiff = c(0.5), vpDiffsd = 0.1, 
                        vDiff = c(0.5, 0.7, 0.8), vInheritance = c(0.5, 0.5, 0.7),
                        propInherite = 0.8, rateDiff = 2, minRate = 1, propHetero = 0.4, 
                        maxPercReads = 99.9),
            error=conditionMessage)
    
    exp <- paste0("vInheritance must be a vector of distinct non-negative ", 
                  "double include in [0,1]")
    
    message <- paste0("test.validateRunSimDoubleParameters_vInheritance_duplicate() ",
                      "- Duplicate elements in vInheritance parameter did not 
                      generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimDoubleParameters_vInheritance_negative <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimDoubleParameters(
                vpDiff = c(0.5), vpDiffsd = 0.1, vDiff = 0.8, 
                vInheritance = -1, propInherite = 0.8, rateDiff = 2, minRate = 1, 
                propHetero = 0.4, maxPercReads = 99.9),
            error=conditionMessage)
    
    exp <- paste0("vInheritance must be a vector of distinct non-negative ", 
                  "double include in [0,1]")
    
    message <- paste0("test.validateRunSimDoubleParameters_vInheritance_negative() ",
                      "- Negative in vInheritance parameter did not ",
                      "generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimDoubleParameters_vInheritance_two <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimDoubleParameters(
                        vpDiff = c(0.5), vpDiffsd = 0.1, vDiff = 0.8, 
                        vInheritance = 2, propInherite = 0.8, rateDiff = 2, 
                        minRate = 1, propHetero = 0.4, maxPercReads = 99.9),
            error=conditionMessage)
    
    exp <- paste0("vInheritance must be a vector of distinct non-negative ", 
                  "double include in [0,1]")
    
    message <- paste0("test.validateRunSimDoubleParameters_vInheritance_two() ",
                      "- Two in vInheritance parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimDoubleParameters_minrate_1 <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimDoubleParameters(
                    vpDiff = 0.2, vpDiffsd = 0.1, vDiff = 0.2, vInheritance = 0.2,
                    propInherite = 0.8, rateDiff = 0.8, minRate = 1, propHetero = 0.4, 
                    maxPercReads = 99.9),
        error=conditionMessage)
    
    exp <- "minRate must be a non-negative double include in [0,1)"
    
    message <- paste0("test.validateRunSimDoubleParameters_minrate_sup_to_1() ",
                      "- 1 as minRate parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}


test.validateRunSimDoubleParameters_minrate_vector <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimDoubleParameters(vpDiff = 0.2, 
                vpDiffsd = 0.1, vDiff = 0.2, vInheritance = 0.2,
                propInherite = 0.8, rateDiff = 0.8, minRate = c(0.01, 0.002), 
                propHetero = 0.4, maxPercReads = 99.9),
        error=conditionMessage)
    
    exp <- "minRate must be a non-negative double include in [0,1)"
    
    message <- paste0("test.validateRunSimDoubleParameters_minrate_vector() ",
                      "- Vector as minRate parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimDoubleParameters_propHetero_vector <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimDoubleParameters(
            vpDiff = 0.2, vpDiffsd = 0.1, vDiff = 0.2, vInheritance = 0.2,
            propInherite = 0.8, rateDiff = 0.8, minRate = 0.1, propHetero = c(0.4, 0.5), 
            maxPercReads = 99.9),
        error=conditionMessage)
    
    exp <- "propHetero must be a non-negative double include in [0,1]"
    
    message <- paste0("test.validateRunSimDoubleParameters_propHetero_vector() ",
                      "- Vector as propHetero parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimDoubleParameters_propHetero_sup_to_1 <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimDoubleParameters(
                vpDiff = 0.2, vpDiffsd = 0.1, vDiff = 0.2, vInheritance = 0.2,
                propInherite = 0.8, rateDiff = 0.8, minRate = 0.1, 
                propHetero = 1.001, maxPercReads = 99.9),
        error=conditionMessage)
    
    exp <- "propHetero must be a non-negative double include in [0,1]"
    
    message <- paste0("test.validateRunSimDoubleParameters_propHetero_sup_to_1() ",
                      "- Superior to 1 as propHetero parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimDoubleParameters_propHetero_inf_zero <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimDoubleParameters(
                    vpDiff = 0.2, vpDiffsd = 0.1, vDiff = 0.2, vInheritance = 0.2,
                    propInherite = 0.8, rateDiff = 0.8, minRate = 0.1, 
                    propHetero = -0.001, maxPercReads = 99.9),
        error=conditionMessage)
    
    exp <- "propHetero must be a non-negative double include in [0,1]"
    
    message <- paste0("test.validateRunSimDoubleParameters_propHetero_inf_zero() ",
                      "- Inferior to zero as propHetero parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}


test.validateRunSimDoubleParameters_propInherite_inf_zero <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimDoubleParameters(
                vpDiff = 0.2, vpDiffsd = 0.1, vDiff = 0.2, vInheritance = 0.2,
                propInherite = -0.001, rateDiff = 0.8, minRate = 0.1, 
                propHetero = 0.03, maxPercReads = 99.9),
        error=conditionMessage)
    
    exp <- "propInherite must be a non-negative double include in [0,1]"
    
    message <- paste0("test.validateRunSimDoubleParameters_propInherite_inf_zero() ",
                      "- Inferior to zero as propInherite parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimDoubleParameters_propInherite_sup_to_1 <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimDoubleParameters( 
            vpDiff = 0.2, vpDiffsd = 0.1, vDiff = 0.2, vInheritance = 0.2,
            propInherite = 1.001, rateDiff = 0.8, minRate = 0.1, propHetero = 0.3, 
            maxPercReads = 99.9),
        error=conditionMessage)
    
    exp <- "propInherite must be a non-negative double include in [0,1]"
    
    message <- paste0("test.validateRunSimDoubleParameters_propInherite_sup_to_1() ",
                      "- Superior to 1 as propInherite parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimDoubleParameters_propInherite_vector <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimDoubleParameters(
                vpDiff = 0.2, vpDiffsd = 0.1, vDiff = 0.2, vInheritance = 0.2,
                propInherite = c(0.4, 0.5), rateDiff = 0.8, minRate = 0.1, 
                propHetero = 0.3, maxPercReads = 99.9),
        error=conditionMessage)
    
    exp <- "propInherite must be a non-negative double include in [0,1]"
    
    message <- paste0("test.validateRunSimDoubleParameters_propInherite_vector() ",
                      "- Vector as propInherite parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}


test.validateRunSimDoubleParameters_maxPercReads_inf_to_zero <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimDoubleParameters(
                vpDiff = 0.2, vpDiffsd = 0.1, vDiff = 0.2, vInheritance = 0.2,
                    propInherite = 0.6,  rateDiff = 0.8, minRate = 0.1, 
                    propHetero = 0.3, maxPercReads = -0.01),
        error=conditionMessage)
    
    exp <- "maxPercReads must be a positive double between [0,100]"
    
    message <- paste0("test.validateRunSimDoubleParameters_maxPercReads_inf_to_zero() ",
                      "- Inferior to zero as maxPercReads parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimDoubleParameters_maxPercReads_sup_to_100 <- function() {
    obs <- tryCatch(methylInheritanceSim:::validateRunSimDoubleParameters(
                vpDiff = 0.2, vpDiffsd = 0.1, vDiff = 0.2, vInheritance = 0.2,
                propInherite = 0.6, rateDiff = 0.8, minRate = 0.1, 
                propHetero = 0.3, maxPercReads = 100.001),
        error=conditionMessage)
    
    exp <- "maxPercReads must be a positive double between [0,100]"
    
    message <- paste0("test.validateRunSimDoubleParameters_maxPercReads_sup_to_100() ",
                      "- Superior to 100 as maxPercReads parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimDoubleParameters_maxPercReads_vector <- function() {
    obs <- tryCatch(
        methylInheritanceSim:::validateRunSimDoubleParameters(
                vpDiff = 0.2, vpDiffsd = 0.1, vDiff = 0.2, vInheritance = 0.2,
                propInherite = 0.6, rateDiff = 0.8, minRate = 0.1, propHetero = 0.3, 
                maxPercReads = c(99.8, 99.6)),
        error=conditionMessage)
    
    exp <- "maxPercReads must be a positive double between [0,100]"
    
    message <- paste0("test.validateRunSimDoubleParameters_maxPercReads_vector() ",
                      "- Vector as maxPercReads parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimDoubleParameters_rateDiff_1 <- function() {
    obs <- tryCatch(
        methylInheritanceSim:::validateRunSimDoubleParameters(
                vpDiff = 0.2, vpDiffsd = 0.1, vDiff = 0.2, vInheritance = 0.2,
                propInherite = 0.6, rateDiff = 1, minRate = 0.1, propHetero = 0.3, 
                maxPercReads = 99.9),
        error=conditionMessage)
    
    exp <- "rateDiff must be a positive double include in (0,1)"
    
    message <- paste0("test.validateRunSimDoubleParameters_rateDiff_1() ",
                      "- 1 as rateDiff parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimDoubleParameters_rateDiff_zero <- function() {
    obs <- tryCatch(
        methylInheritanceSim:::validateRunSimDoubleParameters(
                vpDiff = 0.2, vpDiffsd = 0.1, vDiff = 0.2, vInheritance = 0.2,
                propInherite = 0.6, rateDiff = 0, minRate = 0.1, propHetero = 0.3, 
                maxPercReads = 99.9),
        error=conditionMessage)
    
    exp <- "rateDiff must be a positive double include in (0,1)"
    
    message <- paste0("test.validateRunSimDoubleParameters_rateDiff_zero() ",
                      "- Zero as rateDiff parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimDoubleParameters_rateDiff_vector <- function() {
    obs <- tryCatch(
        methylInheritanceSim:::validateRunSimDoubleParameters(
            vpDiff = 0.2, vpDiffsd = 0.1, vDiff = 0.2, vInheritance = 0.2,
            propInherite = 0.6, rateDiff = c(0.02, 0.3), minRate = 0.1, 
            propHetero = 0.3, maxPercReads = 99.9),
        error=conditionMessage)
    
    exp <- "rateDiff must be a positive double include in (0,1)"
    
    message <- paste0("test.validateRunSimDoubleParameters_rateDiff_vector() ",
                      "- Vector as rateDiff parameter did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimDoubleParameters_good_01 <- function() {
    obs <- tryCatch(
        methylInheritanceSim:::validateRunSimDoubleParameters(
            vpDiff = 0.2, vpDiffsd = 0.1, vDiff = 0.2, vInheritance = 0.2,
            propInherite = 0.6, rateDiff = 0.2, minRate = 0.1, 
            propHetero = 0.3, maxPercReads = 99.9),
        error=conditionMessage)
    
    exp <- 0
    
    message <- paste0("test.validateRunSimDoubleParameters_good_01() ",
                      "- All good parameters did not generated expected results.")
    
    checkEquals(obs, exp, message)
}


###################################################
## validateRunSimParameters() function
###################################################

test.validateRunSimParameters_good_01 <- function() {
    obs <- tryCatch(
        methylInheritanceSim:::validateRunSimParameters (
            vpDiff = 0.2, vpDiffsd = 0.1, vDiff = 0.2, vInheritance = 0.2,
            propInherite = 0.6, rateDiff = 0.2, minRate = 0.1, 
            propHetero = 0.3, maxPercReads = 99.9, 
            nbSynCHR = 1, nbSimulation = 1, nbBlock = 1, nbCpG = 3, 
            vNbSample = 6, nbGeneration = 3, minReads = 10, meanCov = 40, 
            nbCores = 1, vSeed = 1010, keepDiff = TRUE, saveGRanges = FALSE, 
            saveMethylKit = FALSE, runAnalysis = FALSE, 
            outputDir = "toto", fileID = "4_10_3",
            methData = samplesForChrSynthetic, 
            context = "Cpg", assembly = "RNOR_5.0"),
        error=conditionMessage)
    
    exp <- 0
    
    message <- paste0("test.validateRunSimParameters_good_01() ",
                      "- All good parameters did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimParameters_double_validation <- function() {
    obs <- tryCatch(
        methylInheritanceSim:::validateRunSimParameters (
            vpDiff = 0.2, vpDiffsd = 0.1, vDiff = 0.2, vInheritance = 0.2,
            propInherite = 0.6, rateDiff = "hi", minRate = 0.1, 
            propHetero = 0.3, maxPercReads = 99.9, 
            nbSynCHR = 1, nbSimulation = 1, nbBlock = 1, nbCpG = 3, 
            vNbSample = 6, nbGeneration = 3, minReads = 10, meanCov = 40, 
            nbCores = 1, vSeed = 1010, keepDiff = TRUE, saveGRanges = FALSE, 
            saveMethylKit = FALSE, runAnalysis = FALSE, 
            outputDir = "toto", fileID = "4_10_3",
            methData = samplesForChrSynthetic, 
            context = "Cpg", assembly = "RNOR_5.0"),
        error=conditionMessage)
    
    exp <- "rateDiff must be a positive double include in (0,1)"
    
    message <- paste0("test.validateRunSimParameters_double_validation() ",
                      "- A string instead of a double did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimParameters_logical_validation <- function() {
    obs <- tryCatch(
        methylInheritanceSim:::validateRunSimParameters (
            vpDiff = 0.2, vpDiffsd = 0.1, vDiff = 0.2, vInheritance = 0.2,
            propInherite = 0.6, rateDiff = 0.2, minRate = 0.1, 
            propHetero = 0.3, maxPercReads = 99.9, 
            nbSynCHR = 1, nbSimulation = 1, nbBlock = 1, nbCpG = 3, 
            vNbSample = 6, nbGeneration = 3, minReads = 10, meanCov = 40, 
            nbCores = 1, vSeed = 1010, keepDiff = TRUE, saveGRanges = 0.22, 
            saveMethylKit = FALSE, runAnalysis = FALSE, 
            outputDir = "toto", fileID = "4_10_3",
            methData = samplesForChrSynthetic, 
            context = "Cpg", assembly = "RNOR_5.0"),
        error=conditionMessage)
    
    exp <- "saveGRanges must be a logical"
    
    message <- paste0("test.validateRunSimParameters_logical_validation() ",
                      "- All double instead of a logical did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimParameters_other_validation <- function() {
    obs <- tryCatch(
        methylInheritanceSim:::validateRunSimParameters (
            vpDiff = 0.2, vpDiffsd = 0.1, vDiff = 0.2, vInheritance = 0.2,
            propInherite = 0.6, rateDiff = 0.2, minRate = 0.1, 
            propHetero = 0.3, maxPercReads = 99.9, 
            nbSynCHR = 1, nbSimulation = 1, nbBlock = 1, nbCpG = 3, 
            vNbSample = 6, nbGeneration = 3, minReads = 10, meanCov = 40, 
            nbCores = 1, vSeed = 1010, keepDiff = TRUE, saveGRanges = TRUE, 
            saveMethylKit = FALSE, runAnalysis = FALSE, 
            outputDir = "toto", fileID = 0.33,
            methData = samplesForChrSynthetic, 
            context = "Cpg", assembly = "RNOR_5.0"),
        error=conditionMessage)
    
    exp <- "fileID must be a character string or NULL"
    
    message <- paste0("test.validateRunSimParameters_other_validation() ",
                      "- All double instead of a string did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.validateRunSimParameters_integer_validation <- function() {
    obs <- tryCatch(
        methylInheritanceSim:::validateRunSimParameters (
            vpDiff = 0.2, vpDiffsd = 0.1, vDiff = 0.2, vInheritance = 0.2,
            propInherite = 0.6, rateDiff = 0.2, minRate = 0.1, 
            propHetero = 0.3, maxPercReads = 99.9, 
            nbSynCHR = 1, nbSimulation = "home", nbBlock = 1, nbCpG = 3, 
            vNbSample = 6, nbGeneration = 3, minReads = 10, meanCov = 40, 
            nbCores = 1, vSeed = 1010, keepDiff = TRUE, saveGRanges = TRUE, 
            saveMethylKit = FALSE, runAnalysis = FALSE, 
            outputDir = "toto", fileID = "5_3",
            methData = samplesForChrSynthetic, 
            context = "Cpg", assembly = "RNOR_5.0"),
        error=conditionMessage)
    
    exp <- "nbSimulation must be a positive integer or numeric"
    
    message <- paste0("test.validateRunSimParameters_integer_validation() ",
                      "- All string instead of a integer did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

###################################################
## fixSeed() function
###################################################

test.fixSeed_value_not_minus_one <- function() {
    
    set.seed(1010)
    obs <- tryCatch(
        methylInheritanceSim:::fixSeed(vSeed = 101),
        error=conditionMessage)
    
    exp <- 101
    
    message <- paste0("test.fixSeed_value_minus_one() ",
                "- vSeed not equal -1 did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

###################################################
## calculateNbDiffCase() function
###################################################

test.calculateNbDiffCase_good_01 <- function() {
    
    set.seed(1010)
    obs <- tryCatch(
            methylInheritanceSim:::calculateNbDiffCase(nbCase = 12, 
                        propDiff = 0.8, propDiffSd = 0.6), 
            error=conditionMessage)
    
    exp <- 7
    
    message <- paste0("test.calculateNbDiffCase_good_01() ",
                      "- Parameters did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.calculateNbDiffCase_good_sd_small <- function() {
    
    obs <- tryCatch({set.seed(1010)
                methylInheritanceSim:::calculateNbDiffCase(nbCase = 6, 
                    propDiff = 0.8, propDiffSd  = 0.000000001)}, 
            error=conditionMessage)
    
    exp <- 5
    
    message <- paste0("test.calculateNbDiffCase_good_sd_small() ",
                      "- Parameters did not generated expected results.")
    
    checkEquals(obs, exp, message)
}


###################################################
## testIfAlreadyDone() function
###################################################

test.testIfAlreadyDone_false_001 <- function() {
    
    obs <- tryCatch({methylInheritanceSim:::testIfAlreadyDone(pathOut = ".",
                                preference = "AEI", id = 10200011,
                                saveGRanges = FALSE, saveMethylKit = FALSE,
                                runAnalysis = FALSE)}, 
        error=conditionMessage)
    
    exp <- FALSE
    
    message <- paste0("test.testIfAlreadyDone_false_001() ",
                      "- Parameters did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.testIfAlreadyDone_false_002 <- function() {
    
    obs <- tryCatch({methylInheritanceSim:::testIfAlreadyDone(pathOut = ".",
                                preference = "AEI", id = 10200021,
                                saveGRanges = TRUE, saveMethylKit = FALSE,
                                runAnalysis = FALSE)}, 
        error=conditionMessage)
    
    exp <- FALSE
    
    message <- paste0("test.testIfAlreadyDone_false_002() ",
                      "- Parameters did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.testIfAlreadyDone_false_003 <- function() {
    
    obs <- tryCatch({methylInheritanceSim:::testIfAlreadyDone(pathOut = ".",
                                    preference = "AEI", id = 10202021,
                                    saveGRanges = FALSE, saveMethylKit = TRUE,
                                    runAnalysis = FALSE)}, 
                    error=conditionMessage)
    
    exp <- FALSE
    
    message <- paste0("test.testIfAlreadyDone_false_003() ",
                      "- Parameters did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

test.testIfAlreadyDone_false_003 <- function() {
    
    obs <- tryCatch({methylInheritanceSim:::testIfAlreadyDone(pathOut = ".",
                            preference = "AEI", id = 10222021,
                            saveGRanges = FALSE, saveMethylKit = FALSE,
                            runAnalysis = TRUE)}, 
                    error=conditionMessage)
    
    exp <- FALSE
    
    message <- paste0("test.testIfAlreadyDone_false_003() ",
                      "- Parameters did not generated expected results.")
    
    checkEquals(obs, exp, message)
}

###################################################
## simEachGeneration() function
###################################################

test.getSimNew_good_02 <- function() {
    
    set.seed(22212)
    
    stateInformation <- methylInheritanceSim:::getSyntheticChr(methInfo = samplesForChrSynthetic, 
                                                               nbBlock = 1, nbCpG = 3)
    
    stateDiff <- c(1, 1, 1)
    stateInherite <- c(0, 0, 1)
    
    obs <- methylInheritanceSim:::getSimNew(nbCtrl = 2, nbCase = 3, generation = 2, 
                                            stateInfo = stateInformation, stateDiff = stateDiff, 
                                            stateInherite = stateInherite, diffValue = 10, 
                                            propDiff = 0.6, propDiffsd = 0.3, propInheritance = 0.7, propHetero = 0.2)
    
    exp <- GRangesList()
    
    exp[[1]] <- GRanges(seqnames = seqnames(stateInformation),
                        ranges = ranges(stateInformation),
                        strand =  strand(stateInformation),
                        meanDiff = c(0.000000000000000, 1.000000000000000, 1.000000000000000), 
                        meanCTRL = mcols(stateInformation)[3],
                        partitionCase = c(2, 2, 2), partitionCtrl = c(1, 1, 1),
                        ctrl.V1 = c(0.599248212271003, 0.104673052925196, 0.0163259007084415),
                        ctrl.V2 = c(0.660878027594768, 0.0737704350938676, 0.00412068310774805),
                        case.V1 = c(0.000000000000000, 1.000000000000000, 1.000000000000000),
                        case.V2 = c(0.000000000000000, 1.000000000000000, 1.000000000000000),
                        case.V3 = c(0.140108118897419, 0.152428811636903, 0.0242355851971473))
    
    exp[[2]] <- GRanges(seqnames = seqnames(stateInformation),
                        ranges = ranges(stateInformation),
                        strand =  strand(stateInformation),
                        meanDiff = c(0.534404290296076, 0.288519632432523, 1.000000000000000), 
                        meanCTRL = mcols(stateInformation)[3],
                        partitionCase = c(0, 0, 1), partitionCtrl = c(3, 3, 2),
                        ctrl.V1 = c(0.0851227837846553, 0.281923948948365, 0.0198953426761693),
                        ctrl.V2 = c(0.393339608860298, 0.118756249033386, 0.019181448769631),
                        case.V1 = c(0.960713455049359, 0.140418843034454, 1.000000000000000),
                        case.V2 = c(0.935659881913647, 0.282094425404783, 0.0137838586380334),
                        case.V3 = c(0.495670278966184, 0.476562011603015, 0.00593249270592643))
    
    message <- paste0("test.getSimNew_good_02() ",
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
## createSampleID() function
###################################################

test.createSampleID_good_01() {
    
    obs <- methylInheritanceSim:::createSampleID(nbGeneration = 3, 
                                                nbSample = 5)
    
    exp <- list()
    exp[[1]] <- list("F1_1_C", "F1_2_C", "F1_3_C", "F1_4_C", "F1_5_C",
                     "F1_1_C", "F1_2_OC", "F1_3_OC", "F1_4_OC", "F1_5_OC")
    exp[[2]] <- list("F2_1_C", "F2_2_C", "F2_3_C", "F2_4_C", "F2_5_C",
                     "F2_1_C", "F2_2_OC", "F2_3_OC", "F2_4_OC", "F2_5_OC")
    exp[[3]] <- list("F3_1_C", "F3_2_C", "F3_3_C", "F3_4_C", "F3_5_C",
                     "F3_1_C", "F3_2_OC", "F3_3_OC", "F3_4_OC", "F3_5_OC")
    
    message <- paste0("test.createSampleID_good_01() ",
                      "- Valid parameters did not generated expected results.")
    
    checkEquals(obs, exp, message)
}


###################################################
## simInheritance() function
###################################################

test.simInheritance_001 <- function() {
    
    temp_dir <- "simInheritance_001"
    
    stateDiff <- list()
    stateDiff[["stateDiff"]] <- c(1, 0, 1)
    stateDiff[["stateInherite"]] <- c(1, 0, 0)
    
    pref = "S1_6_0.9_0.8_0.5"
    
    set.seed(1022211)
    
    methylInheritanceSim:::simInheritance(pathOut = temp_dir,
                                             pref = pref, k = 1, nbCtrl = 2, nbCase = 2, 
                                             treatment = dataSimExample$treatment, sample.id = dataSimExample$sample.id,
                                             generation = 3, stateInfo = dataSimExample$stateInfo[1:3],
                                             propDiff = 0.3, propDiffsd = 0.1, diffValue = 0.4, 
                                             propInheritance = 0.5, rateDiff = 0.3, minRate = 0.3,
                                             propInherite = 0.3, propHetero = 0.5, minReads = 10, maxPercReads = 99, 
                                             assembly="RNOR_5.0", context="Cpg", meanCov = 40, diffRes = stateDiff,
                                             saveGRanges = FALSE, saveMethylKit = FALSE, runAnalysis = FALSE)
    
    checkTrue(file.exists(paste0(temp_dir, "/simV0.1_", pref, "_1.rds")))
    checkTrue(file.exists(paste0(temp_dir, "/stateDiff_", pref, "_1.rds")))
    
    obsA <- readRDS(paste0(temp_dir, "/simV0.1_", pref, "_1.rds"))
    
    obsB <- readRDS(paste0(temp_dir, "/stateDiff_", pref, "_1.rds"))
    
    expA_01 <- GenomicRanges::GRanges(seqnames = rep("S", 3), 
                                      ranges = IRanges::IRanges(start = c(1000, 1038, 1061), 
                                                                end = c(1000, 1038, 1061)),
                                      strand = rep("+", 3), 
                                      meanDiff = c(0.581456213959732, 0.98720191508078, 0.591465190869137), 
                                      meanCTRL.meanCTRL = c(0.981456213959732, 0.98720191508078, 0.991465190869137),
                                      partitionCase = c(1, 0, 1), partitionCtrl = c(1, 2, 1),
                                      ctrl.V1 = c(0.965762614001235, 0.989802034851579, 0.99493764010319),
                                      ctrl.V2 = c(0.998157064048423, 0.990140423002109, 0.998725359354295),
                                      case.V1 = c(0.565030506723245, 0.984007596913212, 0.596109858589778),
                                      case.V2 = c(0.971297235048743, 0.995951445476427, 0.983912507813916)
    )
    
    expA_02 <- GenomicRanges::GRanges(seqnames = rep("S", 3), 
                                      ranges = IRanges::IRanges(start = c(1000, 1038, 1061), 
                                                                end = c(1000, 1038, 1061)),
                                      strand = rep("+", 3), 
                                      meanDiff = c(0.781456213959732, 0.98720191508078, 0.991465190869137), 
                                      meanCTRL.meanCTRL = c(0.981456213959732, 0.98720191508078, 0.991465190869137),
                                      partitionCase = c(0, 0, 0), partitionCtrl = c(2, 2, 2),
                                      ctrl.V1 = c(0.993284935389844, 0.994683456177534, 0.991836135973971),
                                      ctrl.V2 = c(0.996819782343076, 0.992060001214114, 0.95774752683318),
                                      case.V1 = c(0.987046268678369, 0.996191533238895, 0.999325606332785),
                                      case.V2 = c(0.998751494085842, 0.99679989523906, 0.999875895212403)
    )
    
    expA_03 <- GenomicRanges::GRanges(seqnames = rep("S", 3), 
                                      ranges = IRanges::IRanges(start = c(1000, 1038, 1061), 
                                                                end = c(1000, 1038, 1061)),
                                      strand = rep("+", 3), 
                                      meanDiff = c(0.781456213959732, 0.98720191508078, 0.991465190869137), 
                                      meanCTRL.meanCTRL = c(0.981456213959732, 0.98720191508078, 0.991465190869137),
                                      partitionCase = c(0, 0, 0), partitionCtrl = c(2, 2, 2),
                                      ctrl.V1 = c(0.982955148466405, 0.99375655446083, 0.997292420254736),
                                      ctrl.V2 = c(0.976433855312209, 0.970577990158609, 0.999528942183618),
                                      case.V1 = c(0.965176751427782, 0.934997728804191, 0.970239614281013),
                                      case.V2 = c(0.961234834897086, 0.929809760271278, 0.994619422251482)
    )
    
    
    message <- paste0("test.simInheritance_001() ",
                      "- Valid parameters did not generated expected results.")
    
    expA <- GRangesList(list(expA_01, expA_02, expA_03))
    
    checkEquals(obsA, expA, message)
    checkEquals(obsB$stateDiff, c(1,0,1), message)
    checkEquals(obsB$stateInherite, c(1,0,0), message)
    
    if (dir.exists(temp_dir)) {
        unlink(temp_dir, recursive = TRUE, force = FALSE)
    }
}

test.simInheritance_diffRes_NULL <- function() {
    
    temp_dir <- "simInheritance_diffRes_NULL"
    
    pref = "S1_6_0.9_0.8_0.3"
    
    set.seed(10211211)
    
    methylInheritanceSim:::simInheritance(pathOut = temp_dir,
                                             pref = pref, k = 1, nbCtrl = 2, nbCase = 2, 
                                             treatment = dataSimExample$treatment, sample.id = dataSimExample$sample.id,
                                             generation = 3, stateInfo = dataSimExample$stateInfo[1:3],
                                             propDiff = 0.4, propDiffsd = 0.1, diffValue = 0.3, 
                                             propInheritance = 0.5, rateDiff = 0.3, minRate = 0.3,
                                             propInherite = 0.6, propHetero = 0.6, minReads = 10, maxPercReads = 99, 
                                             assembly="RNOR_5.0", context="Cpg", meanCov = 40, diffRes = NULL,
                                             saveGRanges = FALSE, saveMethylKit = FALSE, runAnalysis = FALSE)
    
    checkTrue(file.exists(paste0(temp_dir, "/simV0.1_", pref, "_1.rds")))
    checkTrue(file.exists(paste0(temp_dir, "/stateDiff_", pref, "_1.rds")))
    
    obsA <- readRDS(paste0(temp_dir, "/simV0.1_", pref, "_1.rds"))
    
    obsB <- readRDS(paste0(temp_dir, "/stateDiff_", pref, "_1.rds"))
    
    expA_01 <- GenomicRanges::GRanges(seqnames = rep("S", 3), 
                                      ranges = IRanges::IRanges(start = c(1000, 1038, 1061), 
                                                                end = c(1000, 1038, 1061)),
                                      strand = rep("+", 3), 
                                      meanDiff = c(0.981456213959732, 0.68720191508078 , 0.691465190869137), 
                                      meanCTRL.meanCTRL = c(0.981456213959732, 0.98720191508078, 0.991465190869137),
                                      partitionCase = c(0, 1, 1), partitionCtrl = c(2, 1, 1),
                                      ctrl.V1 = c(0.980405366890465, 0.986409829956374, 0.904795963240746),
                                      ctrl.V2 = c(0.929665588115341, 0.981403800451126, 0.99696121383443),
                                      case.V1 = c(0.967830792739589, 0.702550883295416, 0.683578328213083),
                                      case.V2 = c(0.991596262178303, 0.994776962242735, 0.957878906176306)
    )
    
    expA_02 <- GenomicRanges::GRanges(seqnames = rep("S", 3), 
                                      ranges = IRanges::IRanges(start = c(1000, 1038, 1061), 
                                                                end = c(1000, 1038, 1061)),
                                      strand = rep("+", 3), 
                                      meanDiff = c(0.981456213959732, 0.80720191508078, 0.811465190869137), 
                                      meanCTRL.meanCTRL = c(0.981456213959732, 0.98720191508078, 0.991465190869137),
                                      partitionCase = c(0, 1, 1), partitionCtrl = c(2, 1, 1),
                                      ctrl.V1 = c(0.98822985553992, 0.979963548035665, 0.99980239380256),
                                      ctrl.V2 = c(0.996994140414723, 0.954115772022289, 0.98886421783845),
                                      case.V1 = c(0.977768449003788, 0.808780405265338, 0.821550385656033),
                                      case.V2 = c(0.984957410039115, 0.99531216553522, 0.973207250124185)
    )
    
    expA_03 <- GenomicRanges::GRanges(seqnames = rep("S", 3), 
                                      ranges = IRanges::IRanges(start = c(1000, 1038, 1061), 
                                                                end = c(1000, 1038, 1061)),
                                      strand = rep("+", 3), 
                                      meanDiff = c(0.981456213959732, 0.80720191508078, 0.811465190869137), 
                                      meanCTRL.meanCTRL = c(0.981456213959732, 0.98720191508078, 0.991465190869137),
                                      partitionCase = c(0, 0, 0), partitionCtrl = c(2, 2, 2),
                                      ctrl.V1 = c(0.995045148055923, 0.986799274399891, 0.996865263430278),
                                      ctrl.V2 = c(0.992783863201566, 0.975460428136066, 0.992805282909862),
                                      case.V1 = c(0.94472919453394, 0.996687514947743, 0.999929013471419),
                                      case.V2 = c(0.972333959759893, 0.991803664258718, 0.998985460163771)
    )
    
    
    message <- paste0("test.simInheritanceNew_diffRes_NULL() ",
                      "- Valid parameters did not generated expected results.")
    
    expA <- GRangesList(list(expA_01, expA_02, expA_03))
    
    checkEquals(obsA, expA, message)
    checkEquals(obsB$stateDiff, c(0,1,1), message)
    checkEquals(obsB$stateInherite, c(0,1,1), message)
    
    if (dir.exists(temp_dir)) {
        unlink(temp_dir, recursive = TRUE, force = FALSE)
    }
}

###################################################
## simInheritanceNew() function
###################################################

test.simInheritanceNew_001 <- function() {
    
    temp_dir <- "simInheritanceNew_001"
    
    stateDiff <- list()
    stateDiff[["stateDiff"]] <- c(1, 0, 1)
    stateDiff[["stateInherite"]] <- c(1, 0, 0)
    
    pref = "S1_6_0.9_0.8_0.5"
    
    set.seed(1022211)
    
    methylInheritanceSim:::simInheritanceNew(pathOut = temp_dir,
        pref = pref, k = 1, nbCtrl = 2, nbCase = 2, 
        treatment = dataSimExample$treatment, sample.id = dataSimExample$sample.id,
                        generation = 3, stateInfo = dataSimExample$stateInfo[1:3],
                        propDiff = 0.3, propDiffsd = 0.1, diffValue = 0.4, 
                        propInheritance = 0.5, rateDiff = 0.3, minRate = 0.3,
                        propInherite = 0.3, propHetero = 0.5, minReads = 10, maxPercReads = 99, 
                        assembly="RNOR_5.0", context="Cpg", meanCov = 40, diffRes = stateDiff,
                        saveGRanges = FALSE, saveMethylKit = FALSE, runAnalysis = FALSE)
    
    checkTrue(file.exists(paste0(temp_dir, "/simV0.1_", pref, "_1.rds")))
    checkTrue(file.exists(paste0(temp_dir, "/stateDiff_", pref, "_1.rds")))
    
    obsA <- readRDS(paste0(temp_dir, "/simV0.1_", pref, "_1.rds"))
    
    obsB <- readRDS(paste0(temp_dir, "/stateDiff_", pref, "_1.rds"))
    
    expA_01 <- GenomicRanges::GRanges(seqnames = rep("S", 3), 
                                  ranges = IRanges::IRanges(start = c(1000, 1038, 1061), 
                                                            end = c(1000, 1038, 1061)),
                                  strand = rep("+", 3), 
                                  meanDiff = c(0.581456213959732, 0.98720191508078, 0.591465190869137), 
                                  meanCTRL.meanCTRL = c(0.981456213959732, 0.98720191508078, 0.991465190869137),
                                  partitionCase = c(1, 0, 1), partitionCtrl = c(1, 2, 1),
                                  ctrl.V1 = c(0.965762614001235, 0.989802034851579, 0.99493764010319),
                                  ctrl.V2 = c(0.998157064048423, 0.990140423002109, 0.998725359354295),
                                  case.V1 = c(0.565030506723245, 0.984007596913212, 0.596109858589778),
                                  case.V2 = c(0.971297235048743, 0.995951445476427, 0.983912507813916)
                                  )
    
    expA_02 <- GenomicRanges::GRanges(seqnames = rep("S", 3), 
                                      ranges = IRanges::IRanges(start = c(1000, 1038, 1061), 
                                                                end = c(1000, 1038, 1061)),
                                      strand = rep("+", 3), 
                                      meanDiff = c(0.781456213959732, 0.98720191508078, 0.991465190869137), 
                                      meanCTRL.meanCTRL = c(0.981456213959732, 0.98720191508078, 0.991465190869137),
                                      partitionCase = c(0, 0, 0), partitionCtrl = c(2, 2, 2),
                                      ctrl.V1 = c(0.993284935389844, 0.994683456177534, 0.991836135973971),
                                      ctrl.V2 = c(0.996819782343076, 0.992060001214114, 0.95774752683318),
                                      case.V1 = c(0.987046268678369, 0.996191533238895, 0.999325606332785),
                                      case.V2 = c(0.998751494085842, 0.99679989523906, 0.999875895212403)
    )
    
    expA_03 <- GenomicRanges::GRanges(seqnames = rep("S", 3), 
                                      ranges = IRanges::IRanges(start = c(1000, 1038, 1061), 
                                                                end = c(1000, 1038, 1061)),
                                      strand = rep("+", 3), 
                                      meanDiff = c(0.781456213959732, 0.98720191508078, 0.991465190869137), 
                                      meanCTRL.meanCTRL = c(0.981456213959732, 0.98720191508078, 0.991465190869137),
                                      partitionCase = c(0, 0, 0), partitionCtrl = c(2, 2, 2),
                                      ctrl.V1 = c(0.982955148466405, 0.99375655446083, 0.997292420254736),
                                      ctrl.V2 = c(0.976433855312209, 0.970577990158609, 0.999528942183618),
                                      case.V1 = c(0.965176751427782, 0.934997728804191, 0.970239614281013),
                                      case.V2 = c(0.961234834897086, 0.929809760271278, 0.994619422251482)
    )
    
    
    message <- paste0("test.simInheritanceNew_001() ",
                      "- Valid parameters did not generated expected results.")
    
    expA <- GRangesList(list(expA_01, expA_02, expA_03))
    
    checkEquals(obsA, expA, message)
    checkEquals(obsB$stateDiff, c(1,0,1), message)
    checkEquals(obsB$stateInherite, c(1,0,0), message)
    
    if (dir.exists(temp_dir)) {
         unlink(temp_dir, recursive = TRUE, force = FALSE)
    }
}

test.simInheritanceNew_diffRes_NULL <- function() {
    
    temp_dir <- "simInheritanceNew_diffRes_NULL"
    
    pref = "S1_6_0.9_0.8_0.3"
    
    set.seed(10211211)
    
    methylInheritanceSim:::simInheritanceNew(pathOut = temp_dir,
                                             pref = pref, k = 1, nbCtrl = 2, nbCase = 2, 
                                             treatment = dataSimExample$treatment, sample.id = dataSimExample$sample.id,
                                             generation = 3, stateInfo = dataSimExample$stateInfo[1:3],
                                             propDiff = 0.4, propDiffsd = 0.1, diffValue = 0.3, 
                                             propInheritance = 0.5, rateDiff = 0.3, minRate = 0.3,
                                             propInherite = 0.6, propHetero = 0.6, minReads = 10, maxPercReads = 99, 
                                             assembly="RNOR_5.0", context="Cpg", meanCov = 40, diffRes = NULL,
                                             saveGRanges = FALSE, saveMethylKit = FALSE, runAnalysis = FALSE)
    
    checkTrue(file.exists(paste0(temp_dir, "/simV0.1_", pref, "_1.rds")))
    checkTrue(file.exists(paste0(temp_dir, "/stateDiff_", pref, "_1.rds")))
    
    obsA <- readRDS(paste0(temp_dir, "/simV0.1_", pref, "_1.rds"))
    
    obsB <- readRDS(paste0(temp_dir, "/stateDiff_", pref, "_1.rds"))
    
    expA_01 <- GenomicRanges::GRanges(seqnames = rep("S", 3), 
                                      ranges = IRanges::IRanges(start = c(1000, 1038, 1061), 
                                                                end = c(1000, 1038, 1061)),
                                      strand = rep("+", 3), 
                                      meanDiff = c(0.981456213959732, 0.68720191508078 , 0.691465190869137), 
                                      meanCTRL.meanCTRL = c(0.981456213959732, 0.98720191508078, 0.991465190869137),
                                      partitionCase = c(0, 1, 1), partitionCtrl = c(2, 1, 1),
                                      ctrl.V1 = c(0.980405366890465, 0.986409829956374, 0.904795963240746),
                                      ctrl.V2 = c(0.929665588115341, 0.981403800451126, 0.99696121383443),
                                      case.V1 = c(0.967830792739589, 0.702550883295416, 0.683578328213083),
                                      case.V2 = c(0.991596262178303, 0.994776962242735, 0.957878906176306)
    )
    
    expA_02 <- GenomicRanges::GRanges(seqnames = rep("S", 3), 
                                      ranges = IRanges::IRanges(start = c(1000, 1038, 1061), 
                                                                end = c(1000, 1038, 1061)),
                                      strand = rep("+", 3), 
                                      meanDiff = c(0.981456213959732, 0.80720191508078, 0.811465190869137), 
                                      meanCTRL.meanCTRL = c(0.981456213959732, 0.98720191508078, 0.991465190869137),
                                      partitionCase = c(0, 1, 1), partitionCtrl = c(2, 1, 1),
                                      ctrl.V1 = c(0.98822985553992, 0.979963548035665, 0.99980239380256),
                                      ctrl.V2 = c(0.996994140414723, 0.954115772022289, 0.98886421783845),
                                      case.V1 = c(0.977768449003788, 0.808780405265338, 0.821550385656033),
                                      case.V2 = c(0.984957410039115, 0.99531216553522, 0.973207250124185)
    )
    
    expA_03 <- GenomicRanges::GRanges(seqnames = rep("S", 3), 
                                      ranges = IRanges::IRanges(start = c(1000, 1038, 1061), 
                                                                end = c(1000, 1038, 1061)),
                                      strand = rep("+", 3), 
                                      meanDiff = c(0.981456213959732, 0.80720191508078, 0.811465190869137), 
                                      meanCTRL.meanCTRL = c(0.981456213959732, 0.98720191508078, 0.991465190869137),
                                      partitionCase = c(0, 0, 0), partitionCtrl = c(2, 2, 2),
                                      ctrl.V1 = c(0.995045148055923, 0.986799274399891, 0.996865263430278),
                                      ctrl.V2 = c(0.992783863201566, 0.975460428136066, 0.992805282909862),
                                      case.V1 = c(0.94472919453394, 0.996687514947743, 0.999929013471419),
                                      case.V2 = c(0.972333959759893, 0.991803664258718, 0.998985460163771)
    )
    
    
    message <- paste0("test.simInheritanceNew_diffRes_NULL() ",
                      "- Valid parameters did not generated expected results.")
    
    expA <- GRangesList(list(expA_01, expA_02, expA_03))
    
    checkEquals(obsA, expA, message)
    checkEquals(obsB$stateDiff, c(0,1,1), message)
    checkEquals(obsB$stateInherite, c(0,1,1), message)
    
    if (dir.exists(temp_dir)) {
        unlink(temp_dir, recursive = TRUE, force = FALSE)
    }
}
 
test.simInheritanceNew_saveGRanges_TRUE <- function() {
    
    temp_dir <- "simInheritanceNew_saveGRanges_TRUE"
    
    pref = "S1_6_0.9_0.8_0.3"
    
    set.seed(123)
    
    sampleID <- list()
    sampleID[[1]] <- list("1_1", "1_2", "1_3", "1_4")
    sampleID[[2]] <- list("2_1", "2_2", "2_3", "2_4")
    sampleID[[3]] <- list("3_1", "3_2", "3_3", "3_4")
    
    methylInheritanceSim:::simInheritanceNew(pathOut = temp_dir,
                                pref = pref, k = 1, nbCtrl = 2, nbCase = 2, 
                                treatment = c(0,0,1,1), sample.id = sampleID,
                                generation = 3, stateInfo = dataSimExample$stateInfo[1:3],
                                propDiff = 0.3, propDiffsd = 0.1, diffValue = 0.32, 
                                propInheritance = 0.7, rateDiff = 0.8, minRate = 0.2,
                                propInherite = 0.6, propHetero = 0.6, minReads = 10, maxPercReads = 99, 
                                assembly="RNOR_5.0", context="Cpg", meanCov = 50, diffRes = NULL,
                                saveGRanges = TRUE, saveMethylKit = FALSE, runAnalysis = FALSE)
    
    checkTrue(file.exists(paste0(temp_dir, "/simV0.1_", pref, "_1.rds")))
    checkTrue(file.exists(paste0(temp_dir, "/stateDiff_", pref, "_1.rds")))
    
    obsA <- readRDS(paste0(temp_dir, "/simV0.1_", pref, "_1.rds"))
    
    obsB <- readRDS(paste0(temp_dir, "/stateDiff_", pref, "_1.rds"))
    


    expA_01 <- GenomicRanges::GRanges(seqnames = rep("S", 3),
                                      ranges = IRanges::IRanges(start = c(1000, 1038, 1061),
                                                                end = c(1000, 1038, 1061)),
                                      strand = rep("+", 3),
                                      meanDiff = c(0.661456213959732, 0.98720191508078 , 0.671465190869137),
                                      meanCTRL.meanCTRL = c(0.981456213959732, 0.98720191508078, 0.991465190869137),
                                      partitionCase = c(1, 0, 1), partitionCtrl = c(1, 2, 1),
                                      ctrl.V1 = c(0.985138973501651, 0.994352676705209, 0.997152045841383),
                                      ctrl.V2 = c(0.999356286203539, 0.956946767638153, 0.993679915803624),
                                      case.V1 = c(0.673440460748842, 0.996378572252127, 0.685307020644639),
                                      case.V2 = c(0.992193492009494, 0.970646786573561, 0.984602636334363)
    )

    expA_02 <- GenomicRanges::GRanges(seqnames = rep("S", 3),
                                      ranges = IRanges::IRanges(start = c(1000, 1038, 1061),
                                                                end = c(1000, 1038, 1061)),
                                      strand = rep("+", 3),
                                      meanDiff = c(0.981456213959732, 0.98720191508078, 0.991465190869137),
                                      meanCTRL.meanCTRL = c(0.981456213959732, 0.98720191508078, 0.991465190869137),
                                      partitionCase = c(0, 0, 0), partitionCtrl = c(2, 2, 2),
                                      ctrl.V1 = c(0.995819865253582, 0.998804046362085, 0.981024294172973),
                                      ctrl.V2 = c(0.985741837588672, 0.993989191533387, 0.999318420358833),
                                      case.V1 = c(0.995689488001639, 0.995234901080396, 0.980455170735616),
                                      case.V2 = c(0.994475519231408, 0.985737562045015, 0.976069838301402)
    )

    expA_03 <- GenomicRanges::GRanges(seqnames = rep("S", 3),
                                      ranges = IRanges::IRanges(start = c(1000, 1038, 1061),
                                                                end = c(1000, 1038, 1061)),
                                      strand = rep("+", 3),
                                      meanDiff = c(0.981456213959732, 0.98720191508078, 0.991465190869137),
                                      meanCTRL.meanCTRL = c(0.981456213959732, 0.98720191508078, 0.991465190869137),
                                      partitionCase = c(0, 0, 0), partitionCtrl = c(2, 2, 2),
                                      ctrl.V1 = c(0.973568319095094, 0.974062630653083, 0.999038035583776),
                                      ctrl.V2 = c(0.997994564848677, 0.970964498276891, 0.989303882856854),
                                      case.V1 = c(0.988002379181952, 0.986204793093944, 0.984274060618534),
                                      case.V2 = c(0.980530184730523, 0.999069145947362, 0.999911242227282)
    )


    message <- paste0("test.simInheritanceNew_saveGRanges_TRUE() ",
                      "- Valid parameters did not generated expected results.")

    expA <- GRangesList(list(expA_01, expA_02, expA_03))

    checkEquals(obsA, expA, message)
    checkEquals(obsB$stateDiff, c(1,0,1), message)
    checkEquals(obsB$stateInherite, c(0,0,0), message)

    if (dir.exists(temp_dir)) {
        unlink(temp_dir, recursive = TRUE, force = FALSE)
    }
}                                                   
                       

###################################################
## simEachGeneration() function
###################################################
                             
test.simEachGeneration_all_save_false <- function() {
    
    stateInformation <- methylInheritanceSim:::getSyntheticChr(methInfo = samplesForChrSynthetic, 
                                                               nbBlock = 1, nbCpG = 3)
    
    stateDiff  <- c(1, 0, 1)
    stateInherite <- c(1, 0, 0)
    
    sim <- methylInheritanceSim:::getSimNew(nbCtrl = 3, nbCase = 1, 
            generation = 3, stateInfo = stateInformation, stateDiff = stateDiff, 
            stateInherite = stateInherite, diffValue = 10, propDiff = 0.8, 
            propDiffsd = 0.2, propInheritance = 0.8, propHetero = 0.1)
    
    obs <- methylInheritanceSim:::simEachGeneration(simulation = sim, nbCtrl = 3, nbCase = 1, treatment = c(0,0,0,1), 
                sample.id = dataSimExample$sample.id, generation = 3, stateInfo = stateInformation, minReads = 10, 
                maxPercReads = 99, context = "Cpg", assembly = "RNOR_5.0", meanCov = 80, 
                saveGRanges = FALSE, saveMethylKit = FALSE, runAnalysis = FALSE)
    
    message <- paste0("test.simEachGeneration_all_save_false() - ",
                      "All logicial parameters to FALSE did not generate expected results.")
    checkEquals(obs$myObj, list())
    checkEquals(obs$myGR, list())
    checkEquals(obs$meth, list())
    checkEquals(obs$myDiff, list())
}

test.simEachGeneration_all_saveGRanges_true <- function() {
    
    stateInformation <- methylInheritanceSim:::getSyntheticChr(methInfo = samplesForChrSynthetic, 
                                                               nbBlock = 1, nbCpG = 3)
    
    stateDiff  <- c(1, 0, 1)
    stateInherite <- c(1, 0, 0)
    
    sim <- methylInheritanceSim:::getSimNew(nbCtrl = 3, nbCase = 1, 
                                            generation = 3, stateInfo = stateInformation, stateDiff = stateDiff, 
                                            stateInherite = stateInherite, diffValue = 10, propDiff = 0.8, 
                                            propDiffsd = 0.2, propInheritance = 0.8, propHetero = 0.1)
    
    obs <- methylInheritanceSim:::simEachGeneration(simulation = sim, nbCtrl = 3, nbCase = 1, treatment = c(0,0,0,1), 
                                                    sample.id = dataSimExample$sample.id, generation = 3, stateInfo = stateInformation, minReads = 10, 
                                                    maxPercReads = 99, context = "Cpg", assembly = "RNOR_5.0", meanCov = 80, 
                                                    saveGRanges = TRUE, saveMethylKit = FALSE, runAnalysis = FALSE)
    
    message <- paste0("test.simEachGeneration_all_saveGRanges_true() - ",
                      "saveGRanges to TRUE did not generate expected results.")
    checkEquals(obs$myObj, list())
    ## TODO
    ##checkEquals(obs$myGR, list())
    checkEquals(obs$meth, list())
    checkEquals(obs$myDiff, list())
}

test.simEachGeneration_all_saveMethylKit_true <- function() {
    
    set.seed(10112)
    
    stateInformation <- methylInheritanceSim:::getSyntheticChr(methInfo = samplesForChrSynthetic, 
                                                               nbBlock = 1, nbCpG = 3)
    
    stateDiff  <- c(1, 0, 1)
    stateInherite <- c(1, 0, 0)
    
    sampleID <- list()
    sampleID[[1]] <- list("F1_1_C", "F1_2_C", "F1_3_OC")
    sampleID[[2]] <- list("F2_1_C", "F2_2_C", "F2_3_OC")
    sampleID[[3]] <- list("F3_1_C", "F3_2_C", "F3_3_OC")
    
    sim <- methylInheritanceSim:::getSimNew(nbCtrl = 2, nbCase = 1, 
                generation = 3, stateInfo = stateInformation, stateDiff = stateDiff, 
                stateInherite = stateInherite, diffValue = 10, propDiff = 0.8, 
                propDiffsd = 0.2, propInheritance = 0.8, propHetero = 0.1)
    
    obs <- methylInheritanceSim:::simEachGeneration(simulation = sim, nbCtrl = 2, nbCase = 1, treatment = c(0,0,1), 
                sample.id = sampleID, generation = 3, stateInfo = stateInformation, minReads = 10, 
                maxPercReads = 99, context = "Cpg", assembly = "RNOR_5.0", meanCov = 80, 
                saveGRanges = FALSE, saveMethylKit = TRUE, runAnalysis = FALSE)
    
    expGR_1 <- list()
    expGR_1[[1]] <- new("methylRaw", data.frame(chr = rep("S", 3), 
                    start = c(1000, 1011, 1017), end = c(1000, 1011, 1017), 
                    strand = strand(rep("+", 3)), coverage = c(71, 90, 95), 
                    numCs = c(0, 1, 1), numTs = c(71, 89, 94)), sample.id = "F1_1_C", 
                    assembly = "RNOR_5.0", context = "Cpg", resolution = 'base')
    expGR_1[[2]] <- new("methylRaw", data.frame(chr = rep("S", 3),  
                    start = c(1000, 1011, 1017), end = c(1000, 1011, 1017), 
                    strand = strand(rep("+", 3)), coverage = c(93, 92, 93), 
                    numCs = c(1, 4, 0), numTs = c(92, 88, 93)), sample.id = "F1_2_C", 
                    assembly = "RNOR_5.0", context = "Cpg", resolution = 'base')
    expGR_1[[3]] <- new("methylRaw", data.frame(chr = rep("S", 3), 
                    start = c(1000, 1011, 1017), end = c(1000, 1011, 1017), 
                    strand = strand(rep("+", 3)), coverage = c(89, 87, 92), 
                    numCs = c(80, 2, 92), numTs = c(9, 85, 0)), sample.id = "F1_3_OC", 
                    assembly = "RNOR_5.0", context = "Cpg", resolution = 'base')
    expGR_2 <- list()
    expGR_2[[1]] <- new("methylRaw", data.frame(chr = rep("S", 3), 
                        start = c(1000, 1011, 1017), end = c(1000, 1011, 1017), 
                        strand = strand(rep("+", 3)), coverage = c(90, 85, 79), 
                        numCs = c(0, 0, 1), numTs = c(90, 85, 78)), sample.id = "F2_1_C", 
                        assembly = "RNOR_5.0", context = "Cpg", resolution = 'base')
    expGR_2[[2]] <- new("methylRaw", data.frame(chr = rep("S", 3),  
                        start = c(1000, 1011, 1017), end = c(1000, 1011, 1017), 
                        strand = strand(rep("+", 3)), coverage = c(73, 93, 78), 
                        numCs = c(0, 2, 0), numTs = c(73, 91, 78)), sample.id = "F2_2_C", 
                        assembly = "RNOR_5.0", context = "Cpg", resolution = 'base')
    expGR_2[[3]] <- new("methylRaw", data.frame(chr = rep("S", 3), 
                        start = c(1000, 1011, 1017), end = c(1000, 1011, 1017), 
                        strand = strand(rep("+", 3)), coverage = c(83, 79, 78), 
                        numCs = c(83, 1, 0), numTs = c(0, 78, 78)), sample.id = "F2_3_OC", 
                        assembly = "RNOR_5.0", context = "Cpg", resolution = 'base')
    expGR_3 <- list()
    expGR_3[[1]] <- new("methylRaw", data.frame(chr = rep("S", 3), 
                        start = c(1000, 1011, 1017), end = c(1000, 1011, 1017), 
                        strand = strand(rep("+", 3)), coverage = c(80, 73, 84), 
                        numCs = c(0, 0, 1), numTs = c(80, 73, 83)), sample.id = "F3_1_C", 
                        assembly = "RNOR_5.0", context = "Cpg", resolution = 'base')
    expGR_3[[2]] <- new("methylRaw", data.frame(chr = rep("S", 3),  
                        start = c(1000, 1011, 1017), end = c(1000, 1011, 1017), 
                        strand = strand(rep("+", 3)), coverage = c(77, 80, 94), 
                        numCs = c(0, 2, 2), numTs = c(77, 78, 92)), sample.id = "F3_2_C", 
                        assembly = "RNOR_5.0", context = "Cpg", resolution = 'base')
    expGR_3[[3]] <- new("methylRaw", data.frame(chr = rep("S", 3), 
                        start = c(1000, 1011, 1017), end = c(1000, 1011, 1017), 
                        strand = strand(rep("+", 3)), coverage = c(86, 79, 79), 
                        numCs = c(86, 2, 0), numTs = c(0, 77, 79)), sample.id = "F3_3_OC", 
                        assembly = "RNOR_5.0", context = "Cpg", resolution = 'base')
    expGR <- list()
    expGR[[1]] <- new("methylRawList", expGR_1, treatment = c(0, 0, 1))
    expGR[[2]] <- new("methylRawList", expGR_2, treatment = c(0, 0, 1))
    expGR[[3]] <- new("methylRawList", expGR_3, treatment = c(0, 0, 1))
    
    message <- paste0("test.simEachGeneration_all_saveMethylKit_true() - ",
                      "saveMethylKit to TRUE did not generate expected results.")
    
    checkEquals(obs$myObj, expGR, message)
    checkEquals(obs$myGR, list(), message)
    checkEquals(obs$meth, list(), message)
    checkEquals(obs$myDiff, list(), message)
}


test.simEachGeneration_all_runAnalysis_true <- function() {
    
    set.seed(1222122)
    
    stateInformation <- methylInheritanceSim:::getSyntheticChr(methInfo = samplesForChrSynthetic, 
                                                               nbBlock = 1, nbCpG = 4)
    
    stateDiff  <- c(1, 0, 1, 1)
    stateInherite <- c(1, 0, 0, 1)
    
    sampleID <- list()
    sampleID[[1]] <- list("F1_1_C", "F1_2_C", "F1_3_C", "F1_1_OC", "F1_2_OC", "F1_3_OC")
    sampleID[[2]] <- list("F2_1_C", "F2_2_C", "F2_3_C", "F2_1_OC", "F2_2_OC", "F2_3_OC")
    sampleID[[3]] <- list("F3_1_C", "F3_2_C", "F3_3_C", "F3_1_OC", "F3_2_OC", "F3_3_OC")
    
    sim <- methylInheritanceSim:::getSimNew(nbCtrl = 3, nbCase = 3, 
                                            generation = 3, stateInfo = stateInformation, stateDiff = stateDiff, 
                                            stateInherite = stateInherite, diffValue = 10, propDiff = 0.8, 
                                            propDiffsd = 0.2, propInheritance = 0.8, propHetero = 0.1)
    
    obs <- methylInheritanceSim:::simEachGeneration(simulation = sim, nbCtrl = 3, nbCase = 3, treatment = c(0,0,0, 1, 1,1), 
                            sample.id = sampleID, generation = 3, stateInfo = stateInformation, minReads = 3, 
                            maxPercReads = 99, context = "Cpg", assembly = "RNOR_5.0", meanCov = 80, 
                            saveGRanges = FALSE, saveMethylKit = FALSE, runAnalysis = TRUE)
    
    expGR_1 <- list()
    expGR_1[[1]] <- new("methylRaw", data.frame(chr = rep("S", 3), 
                                start = c(1000, 1011, 1017), end = c(1000, 1011, 1017), 
                                strand = strand(rep("+", 3)), coverage = c(71, 90, 95), 
                                numCs = c(0, 1, 1), numTs = c(71, 89, 94)), sample.id = "F1_1_C", 
                        assembly = "RNOR_5.0", context = "Cpg", resolution = 'base')
    expGR_1[[2]] <- new("methylRaw", data.frame(chr = rep("S", 3),  
                                start = c(1000, 1011, 1017), end = c(1000, 1011, 1017), 
                                strand = strand(rep("+", 3)), coverage = c(93, 92, 93), 
                                numCs = c(1, 4, 0), numTs = c(92, 88, 93)), sample.id = "F1_2_C", 
                        assembly = "RNOR_5.0", context = "Cpg", resolution = 'base')
    expGR_1[[3]] <- new("methylRaw", data.frame(chr = rep("S", 3), 
                                start = c(1000, 1011, 1017), end = c(1000, 1011, 1017), 
                                strand = strand(rep("+", 3)), coverage = c(89, 87, 92), 
                                numCs = c(80, 2, 92), numTs = c(9, 85, 0)), sample.id = "F1_3_OC", 
                        assembly = "RNOR_5.0", context = "Cpg", resolution = 'base')
    expGR_2 <- list()
    expGR_2[[1]] <- new("methylRaw", data.frame(chr = rep("S", 3), 
                                start = c(1000, 1011, 1017), end = c(1000, 1011, 1017), 
                                strand = strand(rep("+", 3)), coverage = c(90, 85, 79), 
                                numCs = c(0, 0, 1), numTs = c(90, 85, 78)), sample.id = "F2_1_C", 
                        assembly = "RNOR_5.0", context = "Cpg", resolution = 'base')
    expGR_2[[2]] <- new("methylRaw", data.frame(chr = rep("S", 3),  
                                start = c(1000, 1011, 1017), end = c(1000, 1011, 1017), 
                                strand = strand(rep("+", 3)), coverage = c(73, 93, 78), 
                                numCs = c(0, 2, 0), numTs = c(73, 91, 78)), sample.id = "F2_2_C", 
                        assembly = "RNOR_5.0", context = "Cpg", resolution = 'base')
    expGR_2[[3]] <- new("methylRaw", data.frame(chr = rep("S", 3), 
                                start = c(1000, 1011, 1017), end = c(1000, 1011, 1017), 
                                strand = strand(rep("+", 3)), coverage = c(83, 79, 78), 
                                numCs = c(83, 1, 0), numTs = c(0, 78, 78)), sample.id = "F2_3_OC", 
                        assembly = "RNOR_5.0", context = "Cpg", resolution = 'base')
    expGR_3 <- list()
    expGR_3[[1]] <- new("methylRaw", data.frame(chr = rep("S", 3), 
                                start = c(1000, 1011, 1017), end = c(1000, 1011, 1017), 
                                strand = strand(rep("+", 3)), coverage = c(80, 73, 84), 
                                numCs = c(0, 0, 1), numTs = c(80, 73, 83)), sample.id = "F3_1_C", 
                        assembly = "RNOR_5.0", context = "Cpg", resolution = 'base')
    expGR_3[[2]] <- new("methylRaw", data.frame(chr = rep("S", 3),  
                                start = c(1000, 1011, 1017), end = c(1000, 1011, 1017), 
                                strand = strand(rep("+", 3)), coverage = c(77, 80, 94), 
                                numCs = c(0, 2, 2), numTs = c(77, 78, 92)), sample.id = "F3_2_C", 
                        assembly = "RNOR_5.0", context = "Cpg", resolution = 'base')
    expGR_3[[3]] <- new("methylRaw", data.frame(chr = rep("S", 3), 
                                start = c(1000, 1011, 1017), end = c(1000, 1011, 1017), 
                                strand = strand(rep("+", 3)), coverage = c(86, 79, 79), 
                                numCs = c(86, 2, 0), numTs = c(0, 77, 79)), sample.id = "F3_3_OC", 
                        assembly = "RNOR_5.0", context = "Cpg", resolution = 'base')
    expGR <- list()
    expGR[[1]] <- new("methylRawList", expGR_1, treatment = c(0, 0, 1))
    expGR[[2]] <- new("methylRawList", expGR_2, treatment = c(0, 0, 1))
    expGR[[3]] <- new("methylRawList", expGR_3, treatment = c(0, 0, 1))
    
    
    expDiff <- list()
    expDiff[[1]] <- new("methylDiff", data.frame(chr = c("S"), start = c(1000), end = c(1000), 
                        strand = strand(c("+")), pvalue=c(1.079053269533417e-72),
                        qvalue = c(0), meth.diff=c(63.066202090592341)), sample.ids = unlist(sampleID[[1]]), destranded = FALSE,
                        assembly = "RNOR_5.0", context = "Cpg", treatment = c(0,0,0,1,1,1), resolution = 'base')
    expDiff[[2]] <- new("methylDiff", data.frame(chr = c("S"), start = c(1021), end = c(1021), 
                        strand = strand(c("+")), pvalue=c(0.710492384389323),
                        qvalue = c(0), meth.diff=c(0.421067968237779)), sample.ids = unlist(sampleID[[2]]), destranded = FALSE,
                        assembly = "RNOR_5.0", context = "Cpg", treatment = c(0,0,0,1,1,1), resolution = 'base')
    expDiff[[3]] <- new("methylDiff", data.frame(chr = c("S"), start = c(1037), end = c(1037), 
                        strand = strand(c("+")), pvalue=c(2.227413948761407e-57),
                        qvalue = c(-3.681152020132998e-73), meth.diff=c(63.994307400379505)), sample.ids = unlist(sampleID[[3]]), destranded = FALSE,
                        assembly = "RNOR_5.0", context = "Cpg", treatment = c(0,0,0,1,1,1), resolution = 'base')
    
    message <- paste0("test.simEachGeneration_all_runAnalysis_true() - ",
                      "runAnalysis to TRUE did not generate expected results.")
    
    checkEquals(length(obs$myObj), 3, message)
    checkEquals(obs$myGR, list(), message)
    checkEquals(length(obs$meth), 3, message)
    checkEquals(obs$myDiff, expDiff, message)
}