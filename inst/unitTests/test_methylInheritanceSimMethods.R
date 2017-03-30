###################################################
##
## Test the methylInheritanceSimMethods functions
##
###################################################

###  Test setup

library("GenomicRanges")
library("IRanges")

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
    
    message <- paste0("test.runSim_good_001() ",
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
    checkTrue(file.exists(paste0(temp_dir, "/simData_F1_1_6_0.9_0.8_0.5_1.rds")))
    checkTrue(file.exists(paste0(temp_dir, "/simData_F1_1_6_0.9_0.8_0.5_2.rds")))
    checkTrue(file.exists(paste0(temp_dir, "/stateDiff_F1_1_6_0.9_0.8_0.5_1.rds")))
    checkTrue(file.exists(paste0(temp_dir, "/stateDiff_F1_1_6_0.9_0.8_0.5_2.rds")))
    checkTrue(file.exists(paste0(temp_dir, "/syntheticChr_F1_1.rds")))
    checkTrue(file.exists(paste0(temp_dir, "/treatment_F1_1_6.rds")))
    
    ## Check treatment file
    expTreatment <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
    obsTreatment <- readRDS(paste0(temp_dir, "/treatment_F1_1_6.rds"))
    
    checkEquals(obsTreatment, expTreatment)
    
    ## Check stateInfo file
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
    obsStateInfo <- readRDS(paste0(temp_dir, "/syntheticChr_F1_1.rds"))
    
    checkEquals(obsTreatment, expTreatment)
    
    ## Check stateDiff files
    expStateDiff_1 <- list(stateDiff = c(0, 1, 0, 0, 0, 1, 1, 1), stateInherite = rep(0, 8))
    expStateDiff_2 <- list(stateDiff = c(1, 1, 1, 0, 1, 0, 1, 1), stateInherite = c(0, 0, 0, 0, 1, 0, 1, 1))
    obsStateDiff_1 <- readRDS(paste0(temp_dir, "/stateDiff_F1_1_6_0.9_0.8_0.5_1.rds"))
    obsStateDiff_2 <- readRDS(paste0(temp_dir, "/stateDiff_F1_1_6_0.9_0.8_0.5_2.rds"))
    
    checkEquals(obsStateDiff_1, expStateDiff_1)
    checkEquals(obsStateDiff_2, expStateDiff_2)
    
    ## Check methDiff files
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
    
    ## check MethylGR files
    obsMethylGR_1 <- readRDS(paste0(temp_dir, "/methylGR_F1_1_6_0.9_0.8_0.5_1.rds"))
    obsMethylGR_2 <- readRDS(paste0(temp_dir, "/methylGR_F1_1_6_0.9_0.8_0.5_2.rds"))
    
    methylGR_1_1_1 <- GenomicRanges::GRanges(seqnames = rep("S", 8),
                            ranges = IRanges::IRanges(start = c(1000, 1009, 1021, 42325, 52325, 52340, 52342, 52344), 
                                        end = c(1000, 1009, 1021, 42325, 52325, 52340, 52342, 52344)),
                            strand = strand("+"), coverage = c(89, 77, 71, 72, 83, 94, 92, 70), 
                            numCs = unlist(list(ctrl.V1=88, ctrl.V1=75, ctrl.V1=0, ctrl.V1=2, ctrl.V1=0, ctrl.V1=0, 
                                                ctrl.V1=0, ctrl.V1=1)))
    
    methylGR_1_1_4 <- GenomicRanges::GRanges(seqnames = rep("S", 8),
                            ranges = IRanges::IRanges(start = c(1000, 1009, 1021, 42325, 52325, 52340, 52342, 52344), 
                                    end = c(1000, 1009, 1021, 42325, 52325, 52340, 52342, 52344)),
                            strand = strand("+"), coverage = c(90, 87, 79, 78, 90, 82, 78, 87), 
                            numCs = unlist(list(ctrl.V4=84, ctrl.V4=82, ctrl.V4=0, ctrl.V4=1, ctrl.V4=0, ctrl.V4=1, 
                                                ctrl.V4=2, ctrl.V4=1)))
    
    methylGR_1_1_7 <- GenomicRanges::GRanges(seqnames = rep("S", 8),
                            ranges = IRanges::IRanges(start = c(1000, 1009, 1021, 42325, 52325, 52340, 52342, 52344), 
                                    end = c(1000, 1009, 1021, 42325, 52325, 52340, 52342, 52344)),
                            strand = strand("+"), coverage = c(74, 79, 80, 92, 82, 106, 87, 84), 
                            numCs = unlist(list(case.V1=73, case.V1=15, case.V1=1, case.V1=0, case.V1=0, case.V1=84, 
                                                case.V1=71, case.V1=68)))
    
    methylGR_1_1_9 <- GenomicRanges::GRanges(seqnames = rep("S", 8),
                            ranges = IRanges::IRanges(start = c(1000, 1009, 1021, 42325, 52325, 52340, 52342, 52344), 
                                    end = c(1000, 1009, 1021, 42325, 52325, 52340, 52342, 52344)),
                            strand = strand("+"), coverage = c(81, 84, 86, 84, 72, 72, 86, 72), 
                            numCs = unlist(list(case.V3=76, case.V3=15, case.V3=0, case.V3=1, case.V3=0, case.V3=58, 
                                        case.V3=70, case.V3=58)))
    
    methylGR_1_2_2 <- GenomicRanges::GRanges(seqnames = rep("S", 8),
                            ranges = IRanges::IRanges(start = c(1000, 1009, 1021, 42325, 52325, 52340, 52342, 52344), 
                                                  end = c(1000, 1009, 1021, 42325, 52325, 52340, 52342, 52344)),
                            strand = strand("+"), coverage = c(80, 91, 77, 92, 79, 74, 72, 83), 
                            numCs = unlist(list(ctrl.V2=79, ctrl.V2=84, ctrl.V2=22, ctrl.V2=5, ctrl.V2=0, ctrl.V2=0, 
                                        ctrl.V2=1, ctrl.V2=2)))
    
    methylGR_1_2_5 <- GenomicRanges::GRanges(seqnames = rep("S", 8),
                            ranges = IRanges::IRanges(start = c(1000, 1009, 1021, 42325, 52325, 52340, 52342, 52344), 
                                                  end = c(1000, 1009, 1021, 42325, 52325, 52340, 52342, 52344)),
                            strand = strand("+"), coverage = c(85, 69, 72, 105, 67, 61, 74, 96), 
                            numCs = unlist(list(ctrl.V5=73, ctrl.V5=67, ctrl.V5=8, ctrl.V5=1, ctrl.V5=0, ctrl.V5=1, 
                                                ctrl.V5=4, ctrl.V5=1)))
    
    methylGR_1_2_8 <- GenomicRanges::GRanges(seqnames = rep("S", 8),
                            ranges = IRanges::IRanges(start = c(1000, 1009, 1021, 42325, 52325, 52340, 52342, 52344), 
                                                  end = c(1000, 1009, 1021, 42325, 52325, 52340, 52342, 52344)),
                            strand = strand("+"), coverage = c(84, 98, 79, 70, 93, 86, 92, 77), 
                            numCs = unlist(list(case.V2=84, case.V2=91, case.V2=2, case.V2=8, case.V2=1, case.V2=0, 
                                                case.V2=6, case.V2=1)))
    
    methylGR_1_3_3 <- GenomicRanges::GRanges(seqnames = rep("S", 8),
                            ranges = IRanges::IRanges(start = c(1000, 1009, 1021, 42325, 52325, 52340, 52342, 52344), 
                                                  end = c(1000, 1009, 1021, 42325, 52325, 52340, 52342, 52344)),
                            strand = strand("+"), coverage = c(90, 90, 97, 83, 76, 82, 87, 85), 
                            numCs = unlist(list(ctrl.V3=86, ctrl.V3=87, ctrl.V3=10, ctrl.V3=4, ctrl.V3=4, ctrl.V3=0, 
                                                ctrl.V3=1, ctrl.V3=1)))
    
    methylGR_1_3_11 <- GenomicRanges::GRanges(seqnames = rep("S", 8),
                            ranges = IRanges::IRanges(start = c(1000, 1009, 1021, 42325, 52325, 52340, 52342, 52344), 
                                                   end = c(1000, 1009, 1021, 42325, 52325, 52340, 52342, 52344)),
                            strand = strand("+"), coverage = c(92, 81, 74, 79, 80, 98, 65, 75), 
                            numCs = unlist(list(case.V5=89, case.V5=79, case.V5=8, case.V5=2, case.V5=1, case.V5=0, 
                                                case.V5=1, case.V5=1)))
    
    methylGR_2_1_3 <- GenomicRanges::GRanges(seqnames = rep("S", 8),
                            ranges = IRanges::IRanges(start = c(1000, 1009, 1021, 42325, 52325, 52340, 52342, 52344), 
                                            end = c(1000, 1009, 1021, 42325, 52325, 52340, 52342, 52344)),
                            strand = strand("+"), coverage = c(75, 79, 80, 81, 87, 83, 67, 66), 
                            numCs = unlist(list(ctrl.V3=63, ctrl.V3=77, ctrl.V3=14, ctrl.V3=0, ctrl.V3=0, ctrl.V3=0, 
                                                 ctrl.V3=2, ctrl.V3=1)))
    
    methylGR_2_1_10 <- GenomicRanges::GRanges(seqnames = rep("S", 8),
                            ranges = IRanges::IRanges(start = c(1000, 1009, 1021, 42325, 52325, 52340, 52342, 52344), 
                                                   end = c(1000, 1009, 1021, 42325, 52325, 52340, 52342, 52344)),
                            strand = strand("+"), coverage = c(81, 86, 88, 87, 91, 95, 79, 81), 
                            numCs = unlist(list(case.V4=14, case.V4=16, case.V4=74, case.V4=5, case.V4=73, 
                                                    case.V4=0, case.V4=65, case.V4=67)))
    
    methylGR_2_2_2 <- GenomicRanges::GRanges(seqnames = rep("S", 8),
                            ranges = IRanges::IRanges(start = c(1000, 1009, 1021, 42325, 52325, 52340, 52342, 52344), 
                                                end = c(1000, 1009, 1021, 42325, 52325, 52340, 52342, 52344)),
                            strand = strand("+"), coverage = c(77, 93, 75, 83, 102, 72, 89, 88), 
                            numCs = unlist(list(ctrl.V2=77, ctrl.V2=92, ctrl.V2=3, ctrl.V2=4, ctrl.V2=1, ctrl.V2=0, 
                                                ctrl.V2=2, ctrl.V2=1)))
    
    methylGR_2_2_7 <- GenomicRanges::GRanges(seqnames = rep("S", 8),
                            ranges = IRanges::IRanges(start = c(1000, 1009, 1021, 42325, 52325, 52340, 52342, 52344), 
                                        end = c(1000, 1009, 1021, 42325, 52325, 52340, 52342, 52344)),
                            strand = strand("+"), coverage = c(90, 78, 82, 79, 82, 80, 89, 73), 
                            numCs = unlist(list(case.V1=79, case.V1=76, case.V1=3, case.V1=1, case.V1=33, case.V1=3, 
                                                case.V1=38, case.V1=30)))
    
    methylGR_2_3_4 <- GenomicRanges::GRanges(seqnames = rep("S", 8),
                            ranges = IRanges::IRanges(start = c(1000, 1009, 1021, 42325, 52325, 52340, 52342, 52344), 
                                        end = c(1000, 1009, 1021, 42325, 52325, 52340, 52342, 52344)),
                            strand = strand("+"), coverage = c(90, 80, 87, 71, 89, 73, 73, 66), 
                            numCs = unlist(list(ctrl.V4=88, ctrl.V4=74, ctrl.V4=1, ctrl.V4=5, ctrl.V4=0, ctrl.V4=1, 
                                                ctrl.V4=5, ctrl.V4=0)))
    
    methylGR_2_3_12 <- GenomicRanges::GRanges(seqnames = rep("S", 8),
                            ranges = IRanges::IRanges(start = c(1000, 1009, 1021, 42325, 52325, 52340, 52342, 52344), 
                                                   end = c(1000, 1009, 1021, 42325, 52325, 52340, 52342, 52344)),
                            strand = strand("+"), coverage = c(84, 79, 77, 83, 65, 85, 81, 77), 
                            numCs = unlist(list(case.V6=83, case.V6=74, case.V6=9, case.V6=1, case.V6=0, case.V6=0, 
                                                case.V6=2, case.V6=1)))
    
    checkEquals(obsMethylGR_1[[1]][[1]], methylGR_1_1_1, message)
    checkEquals(obsMethylGR_1[[1]][[4]], methylGR_1_1_4, message)
    checkEquals(obsMethylGR_1[[1]][[7]], methylGR_1_1_7, message)
    checkEquals(obsMethylGR_1[[1]][[9]], methylGR_1_1_9, message)
    checkEquals(obsMethylGR_1[[2]][[2]], methylGR_1_2_2, message)
    checkEquals(obsMethylGR_1[[2]][[5]], methylGR_1_2_5, message)
    checkEquals(obsMethylGR_1[[2]][[8]], methylGR_1_2_8, message)
    checkEquals(obsMethylGR_1[[3]][[3]], methylGR_1_3_3, message)
    checkEquals(obsMethylGR_1[[3]][[11]], methylGR_1_3_11, message)
    checkEquals(obsMethylGR_2[[1]][[3]], methylGR_2_1_3, message)
    checkEquals(obsMethylGR_2[[1]][[10]], methylGR_2_1_10, message)
    checkEquals(obsMethylGR_2[[2]][[2]], methylGR_2_2_2, message)
    checkEquals(obsMethylGR_2[[2]][[7]], methylGR_2_2_7, message)
    checkEquals(obsMethylGR_2[[3]][[4]], methylGR_2_3_4, message)
    checkEquals(obsMethylGR_2[[3]][[12]], methylGR_2_3_12, message)
    
    if (dir.exists(temp_dir)) {
        unlink(temp_dir, recursive = TRUE, force = FALSE)
    }
}

test.runSim_NULL_outputDir <- function() {
    
    ## Create 1 simulated dataset (nbSimulation = 1)
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

test.runSim_keepDiff_true <- function() {
    
    temp_dir = "runSim_keepDiff_true"
    
    ## Create 1simulated dataset (nbSimulation = 1)
    ## over 2 generations (nbGeneration = 2) with
    ## 2 cases and 2 controls (vNbsample = 2) using only one set
    ## of parameters (vpDiff = 0.85, vpDiffsd = 0.1, vDiff = 0.8)
    result <- runSim(outputDir = paste0(temp_dir, "/"), fileID = "F1", nbSynCHR = 2,
                     methData = samplesForChrSynthetic, nbSimulation = 1, keepDiff = TRUE,
                     nbBlock = 2, nbCpG = 4, nbGeneration = 2, vNbSample = c(2), vpDiff = c(0.85),
                     vpDiffsd = c(0.1), vDiff = c(0.8), vInheritance = c(0.5), propInherite = 0.6,
                     rateDiff = 0.3, minRate = 0.2, propHetero = 0.5, nbCores = 1, vSeed = 21132, 
                     saveGRanges = FALSE, saveMethylKit = FALSE, runAnalysis = FALSE)
    
    
    message <- paste0("test.runSim_keepDiff_true() ",
                      "- Valid parameters did not generated expected results.")
    
    checkEquals(result, 0, message)
    checkTrue(file.exists(temp_dir))
    checkTrue(file.exists(paste0(temp_dir, "/simData_F1_1_2_0.85_0.8_0.5_1.rds")))
    checkTrue(file.exists(paste0(temp_dir, "/simData_F1_2_2_0.85_0.8_0.5_1.rds")))
    checkTrue(file.exists(paste0(temp_dir, "/syntheticChr_F1_1.rds")))
    checkTrue(file.exists(paste0(temp_dir, "/syntheticChr_F1_2.rds")))
    checkTrue(file.exists(paste0(temp_dir, "/stateDiff_F1_1_2_0.85_0.8_0.5_1.rds")))
    checkTrue(file.exists(paste0(temp_dir, "/stateDiff_F1_2_2_0.85_0.8_0.5_1.rds")))
    
    ## check MethylGR files
    sim_1 <- readRDS(paste0(temp_dir, "/simData_F1_1_2_0.85_0.8_0.5_1.rds"))
    sim_2 <- readRDS(paste0(temp_dir, "/simData_F1_2_2_0.85_0.8_0.5_1.rds"))
    
    checkTrue(any(IRanges::start(sim_1[[1]]) != IRanges::start(sim_2[[1]])))
    
    if (dir.exists(temp_dir)) {
        unlink(temp_dir, recursive = TRUE, force = FALSE)
    }
}


