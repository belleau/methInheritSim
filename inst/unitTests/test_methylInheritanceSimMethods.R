###################################################
##
## Test the methylInheritanceSimMethods functions
##
###################################################


###  Test setup

library("GenomicRanges")

data(samplesForChrSynthetic)
data(dataSimExample)



###################################################
## runSim() function
###################################################

test.runSim_good_001 <- function() {
    
    ## Set the output directory where files will be created
    temp_dir <- "test_runSim_good_01"
    
    ## Create 2 simulated dataset (nbSimulation = 2)
    ## over 3 generations (nbGeneration = 3) with
    ## 6 cases and 6 controls (vNbsample = 6) using only one set
    ## of parameters (vpDiff = 0.9, vpDiffsd = 0.1, vDiff = 0.8)
    result <- runSim(outputDir = temp_dir, fileID = "F1", nbSynCHR = 1,
        methData = samplesForChrSynthetic, nbSimulation = 2,
        nbBlock = 2, nbCpG = 4, nbGeneration = 3, vNbSample = c(6), vpDiff = c(0.9),
        vpDiffsd = c(0.1), vDiff = c(0.8), vInheritance = c(0.5), propInherite = 0.3,
        rateDiff = 0.3, minRate = 0.2, propHetero = 0.5, nbCores = 1, vSeed = 232, 
        saveGRanges = T, saveMethylKit = T, runAnalysis = T)
    
    
    message <- paste0("test.simInheritanceNew_saveGRanges_TRUE() ",
                      "- Valid parameters did not generated expected results.")
    
    checkEquals(result, 0, message)
    checkTrue(file.exists(paste0(temp_dir, "/meth_F1_1_6_0.9_0.8_0.5_1.rds")))
    checkTrue(file.exists(paste0(temp_dir, "/meth_F1_1_6_0.9_0.8_0.5_2.rds")))
    checkTrue(file.exists(paste0(temp_dir, "/methDiff_F1_1_6_0.9_0.8_0.5_1.rds")))
    checkTrue(file.exists(paste0(temp_dir, "/methDiff_F1_1_6_0.9_0.8_0.5_2.rds")))
    checkTrue(file.exists(paste0(temp_dir, "/methylGR_F1_1_6_0.9_0.8_0.5_1.rds")))
    checkTrue(file.exists(paste0(temp_dir, "/methylGR_F1_1_6_0.9_0.8_0.5_2.rds")))
    checkTrue(file.exists(paste0(temp_dir, "/methylObj_F1_1_6_0.9_0.8_0.5_1.rds")))
    checkTrue(file.exists(paste0(temp_dir, "/methylObj_F1_1_6_0.9_0.8_0.5_2.rds")))
    checkTrue(file.exists(paste0(temp_dir, "/simV0.1_F1_1_6_0.9_0.8_0.5_1.rds")))
    checkTrue(file.exists(paste0(temp_dir, "/simV0.1_F1_1_6_0.9_0.8_0.5_2.rds")))
    checkTrue(file.exists(paste0(temp_dir, "/stateDiff_F1_1_6_0.9_0.8_0.5_1.rds")))
    checkTrue(file.exists(paste0(temp_dir, "/stateDiff_F1_1_6_0.9_0.8_0.5_2.rds")))
    checkTrue(file.exists(paste0(temp_dir, "/stateInfo_F1_1.rds")))
    checkTrue(file.exists(paste0(temp_dir, "/treatment_F1_1_6.rds")))
    
    expTreatment <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
    obsTreatment <- readRDS(paste0(temp_dir, "/treatment_F1_1_6.rds"))
    
    checkEquals(obsTreatment, expTreatment)
    
    expStateInfo <- GenomicRanges::GRanges(seqnames = rep("S", 8), 
                                ranges = IRanges::IRanges(start = c(1000, 1003, 1022, 1029, 11029, 11048, 11052, 11058), 
                                                              end = c(1000, 1003, 1022, 1029, 11029, 11048, 11052, 11058)),
                                strand = rep("+", 8), 
                                chrOri = rep(1, 8),
                                startOri = c(3966339, 3966342, 3966361, 3966368, 17813961, 17813980, 17813984, 17813990),
                                meanCTRL = c(0.0189487186549124, 0.019742922537726, 0.00968460541244091,
                                             0.00785548164521422, 0.00886015325670498, 0.00762063636561119,
                                             0.00627637130801688, 0.0155857779775954),
                                varCTRL = c(5.74940270224117e-05, 6.0841609497317e-05, 6.45952177849217e-05,
                                            3.71976397677139e-05, 0.000279443013167746, 3.5149491082417e-05,
                                            0.000109774742295572, 8.53562292547186e-05))
    obsStateInfo <- readRDS(paste0(temp_dir, "/stateInfo_F1_1.rds"))
    
    checkEquals(obsTreatment, expTreatment)
    
    expStateDiff_1 <- list(stateDiff = c(0, 1, 0, 0, 0, 1, 1, 1), stateInherite = rep(0, 8))
    expStateDiff_2 <- list(stateDiff = c(1, 1, 1, 0, 1, 0, 1, 1), stateInherite = c(0, 0, 0, 0, 1, 0, 1, 1))
    obsStateDiff_1 <- readRDS(paste0(temp_dir, "/stateDiff_F1_1_6_0.9_0.8_0.5_1.rds"))
    obsStateDiff_2 <- readRDS(paste0(temp_dir, "/stateDiff_F1_1_6_0.9_0.8_0.5_2.rds"))
    
    checkEquals(obsStateDiff_1, expStateDiff_1)
    checkEquals(obsStateDiff_2, expStateDiff_2)
    
    expStateDiff_1 <- list(stateDiff = c(0, 1, 0, 0, 0, 1, 1, 1), stateInherite = rep(0, 8))
    expStateDiff_2 <- list(stateDiff = c(1, 1, 1, 0, 1, 0, 1, 1), stateInherite = c(0, 0, 0, 0, 1, 0, 1, 1))
    obsStateDiff_1 <- readRDS(paste0(temp_dir, "/stateDiff_F1_1_6_0.9_0.8_0.5_1.rds"))
    obsStateDiff_2 <- readRDS(paste0(temp_dir, "/stateDiff_F1_1_6_0.9_0.8_0.5_2.rds"))
    
    checkEquals(obsStateDiff_1, expStateDiff_1)
    checkEquals(obsStateDiff_2, expStateDiff_2)
    
    obsMethDiff_1 <- readRDS(paste0(temp_dir, "/methDiff_F1_1_6_0.9_0.8_0.5_1.rds"))
    obsMethDiff_2 <- readRDS(paste0(temp_dir, "/methDiff_F1_1_6_0.9_0.8_0.5_2.rds"))
    
    checkEquals(obsMethDiff_1[[1]]@sample.ids, c("F1_1_C", "F1_2_C", "F1_3_C", "F1_4_C", "F1_5_C", 
                        "F1_6_C", "F1_7_OC", "F1_8_OC", "F1_9_OC", "F1_10_OC", "F1_11_OC", "F1_12_OC"))
    checkEquals(obsMethDiff_1[[2]]@sample.ids, c("F2_1_C", "F2_2_C", "F2_3_C", "F2_4_C", "F2_5_C", 
                        "F2_6_C", "F2_7_OC", "F2_8_OC", "F2_9_OC", "F2_10_OC", "F2_11_OC", "F2_12_OC"))
    checkEquals(obsMethDiff_1[[3]]@sample.ids, c("F3_1_C", "F3_2_C", "F3_3_C", "F3_4_C", "F3_5_C", 
                        "F3_6_C", "F3_7_OC", "F3_8_OC", "F3_9_OC", "F3_10_OC", "F3_11_OC", "F3_12_OC"))
    
    checkEquals(obsMethDiff_2[[1]]@sample.ids, c("F1_1_C", "F1_2_C", "F1_3_C", "F1_4_C", "F1_5_C", 
                        "F1_6_C", "F1_7_OC", "F1_8_OC", "F1_9_OC", "F1_10_OC", "F1_11_OC", "F1_12_OC"))
    checkEquals(obsMethDiff_2[[2]]@sample.ids, c("F2_1_C", "F2_2_C", "F2_3_C", "F2_4_C", "F2_5_C", 
                        "F2_6_C", "F2_7_OC", "F2_8_OC", "F2_9_OC", "F2_10_OC", "F2_11_OC", "F2_12_OC"))
    checkEquals(obsMethDiff_2[[3]]@sample.ids, c("F3_1_C", "F3_2_C", "F3_3_C", "F3_4_C", "F3_5_C", 
                        "F3_6_C", "F3_7_OC", "F3_8_OC", "F3_9_OC", "F3_10_OC", "F3_11_OC", "F3_12_OC"))
    
    checkEquals(obsMethDiff_1[[1]]@treatment, c(rep(0, 6), rep(1, 6)))
    checkEquals(obsMethDiff_1[[2]]@treatment, c(rep(0, 6), rep(1, 6)))
    checkEquals(obsMethDiff_1[[3]]@treatment, c(rep(0, 6), rep(1, 6)))
    
    checkEquals(obsMethDiff_2[[1]]@treatment, c(rep(0, 6), rep(1, 6)))
    checkEquals(obsMethDiff_2[[2]]@treatment, c(rep(0, 6), rep(1, 6)))
    checkEquals(obsMethDiff_2[[3]]@treatment, c(rep(0, 6), rep(1, 6)))
    
    checkEquals(obsMethDiff_1[[1]]@destranded, FALSE, message)
    checkEquals(obsMethDiff_1[[2]]@destranded, FALSE, message)
    checkEquals(obsMethDiff_1[[3]]@destranded, FALSE, message)
    
    checkEquals(obsMethDiff_2[[1]]@destranded, FALSE, message)
    checkEquals(obsMethDiff_2[[2]]@destranded, FALSE, message)
    checkEquals(obsMethDiff_2[[3]]@destranded, FALSE, message)
    
    obsMethylGR_1 <- readRDS(paste0(temp_dir, "/methylGR_F1_1_6_0.9_0.8_0.5_1.rds"))
    methylGR_1_1_1 <- new("methylRaw", data.frame(chr = rep("S", 8), 
                                       start = c(1000, 1009, 1021, 42325, 52325, 52340, 52342, 52344), 
                                       end = c(1000, 1009, 1021, 42325, 52325, 52340, 52342, 52344),
                                       strand = strand("+"),
                                       coverage = c(89, 77, 71, 72, 83, 94, 92, 70), 
                                       numCs = c(88, 75,  0,  2,  0,  0,  0,  1),
                                       numTs = c(1, 2, 71, 70, 83, 94, 92, 69)),
               sample.id = "F1_1_C", assembly = "Rnor_5.0",
               context = "CpG", resolution = 'base')
    
    methylGR_1_1_4 <- new("methylRaw", data.frame(chr = rep("S", 8), 
                                                  start = c(1000, 1009, 1021, 42325, 52325, 52340, 52342, 52344), 
                                                  end = c(1000, 1009, 1021, 42325, 52325, 52340, 52342, 52344),
                                                  strand = strand("+"),
                                                  coverage = c(90, 87, 79, 78, 90, 82, 78, 87), 
                                                  numCs = c(84, 82, 0, 1, 0, 1, 2, 1),
                                                  numTs = c(6, 5, 79, 77, 90, 81, 76, 86)),
                          sample.id = "F1_4_C", assembly = "Rnor_5.0",
                          context = "CpG", resolution = 'base')
    
    methylGR_1_1_7 <- new("methylRaw", data.frame(chr = rep("S", 8), 
                                                  start = c(1000, 1009, 1021, 42325, 52325, 52340, 52342, 52344), 
                                                  end = c(1000, 1009, 1021, 42325, 52325, 52340, 52342, 52344),
                                                  strand = strand("+"),
                                                  coverage = c(74, 79, 80, 92, 82, 106, 87, 84), 
                                                  numCs = c(73, 15, 1, 0, 0, 84, 71, 68),
                                                  numTs = c(1, 64, 79, 92, 82, 22, 16, 16)),
                          sample.id = "F1_7_OC", assembly = "Rnor_5.0",
                          context = "CpG", resolution = 'base')
    
    methylGR_1_1_9 <- new("methylRaw", data.frame(chr = rep("S", 8), 
                                                  start = c(1000, 1009, 1021, 42325, 52325, 52340, 52342, 52344), 
                                                  end = c(1000, 1009, 1021, 42325, 52325, 52340, 52342, 52344),
                                                  strand = strand("+"),
                                                  coverage = c(81, 84, 86, 84, 72, 72, 86, 72), 
                                                  numCs = c(76, 15, 0, 1, 0, 58, 70, 58),
                                                  numTs = c(5, 69, 86, 83, 72, 14, 16, 14)),
                          sample.id = "F1_9_OC", assembly = "Rnor_5.0",
                          context = "CpG", resolution = 'base')
    
    checkEquals(obsMethylGR_1[[1]][[1]], methylGR_1_1_1)
    checkEquals(obsMethylGR_1[[1]][[4]], methylGR_1_1_4)
    checkEquals(obsMethylGR_1[[1]][[7]], methylGR_1_1_7)
    checkEquals(obsMethylGR_1[[1]][[9]], methylGR_1_1_9)
    
    if (dir.exists(temp_dir)) {
        unlink(temp_dir, recursive = TRUE, force = FALSE)
    }
}

test.runSim_NULL_outputDir <- function() {
    
    ## Create 1simulated dataset (nbSimulation = 1)
    ## over 2 generations (nbGeneration = 2) with
    ## 2 cases and 2 controls (vNbsample = 2) using only one set
    ## of parameters (vpDiff = 0.9, vpDiffsd = 0.1, vDiff = 0.8)
    result <- runSim(outputDir = NULL, fileID = "F1", nbSynCHR = 1,
                     methData = samplesForChrSynthetic, nbSimulation = 1,
                     nbBlock = 2, nbCpG = 4, nbGeneration = 2, vNbSample = c(2), vpDiff = c(0.9),
                     vpDiffsd = c(0.1), vDiff = c(0.8), vInheritance = c(0.5), propInherite = 0.3,
                     rateDiff = 0.3, minRate = 0.2, propHetero = 0.5, nbCores = 1, vSeed = 22132, 
                     saveGRanges = T, saveMethylKit = T, runAnalysis = T)
    
    
    message <- paste0("test.runSim_NULL_outputDir() ",
                      "- Null outputDir did not generated expected results.")
    
    checkEquals(result, 0, message)
    checkTrue(file.exists("outputDir"))
    
    if (dir.exists("outputDir")) {
        unlink("outputDir", recursive = TRUE, force = FALSE)
    }
}


