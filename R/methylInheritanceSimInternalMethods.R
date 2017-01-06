#' @title TODO
#'
#' @description Estime the paramater alpha from the mean and the variance
#' of a beta distribution
#'
#' @param valCtrl array with first value mean
#' and second
#'
#' @param minVal float (default 1e-06)
#'
#' @return The alpha parameter of a Beta distribution
#'
#' @examples
#'
#' ## TODO
#'
#'
#' @author Pascal Belleau
#' @keywords internal
estBetaAlpha <- function(valCtrl, minVal=1e-06){
    mu <- max(valCtrl[1], minVal)
    
    sigma2 <- max(valCtrl[2], ifelse(mu < 0.01, min(minVal, mu/1000), minVal))
    mu <- min(mu, 1-min(0.001, sigma2 * 10^(-log10(minVal)/2)))
    
    max(0, -mu * (sigma2 + mu^2 - mu) / sigma2)
}


#' @title TODO
#'
#' @description Estime the paramater beta from the mean and the variance
#' of a beta distribution
#'
#' @param valCtrl array with first value mean
#' and second
#'
#' @param minVal float (default 1e-06)
#'
#' @return The beta parameter of a Beta distribution
#'
#' @examples
#'
#' ## TODO
#'
#'
#' @author Pascal Belleau
#' @keywords internal
estBetaBeta <- function(valCtrl, minVal=1e-06){
    mu <- max(valCtrl[1], minVal)
    
    sigma2 <- max(valCtrl[2], ifelse(mu < 0.01,min(minVal, mu/1000), minVal))
    #mu <- min(mu, 1-min(0.05, 100*sigma2))
    mu <- min(mu, 1-min(0.001, sigma2 * 10^(-log10(minVal)/2)))
    max(0, (sigma2 + mu^2 - mu) * (mu -1) / sigma2)
}


#' @title TODO
#'
#' @description TODO
#'
#' @param methInfo is methylKitList from CTRL data.
#'
#' @param m blocks
#'
#' @param block a \code{integer} consecutive position from methInfo
#'
#' @return TODO
#'
#' @examples
#'
#' ## TODO
#'
#'
#' @author Pascal Belleau
#' @keywords internal

getStateData <- function(methInfo, m, block){
    
    seqChr <- unique(methInfo$chr) # list of Chr
    # Sample m chromosome could be the same
    chr <- seqChr[sample(1:length(seqChr),m,replace=TRUE)]
    
    # The position of the CTRL
    posCTRL <- which(methInfo@treatment==0)
    
    # Init the data.frame
    res <- data.frame(chr=rep("S",m*block), start=rep(0,m*block), end=rep(0,m*block),
                      meanCTRL=rep(0,m*block), varCTRL=rep(0,m*block),
                      alphaCTRL=rep(0,m*block), betaCTRL=rep(0,m*block),
                      chrOri=rep(0, m*block), startOri=rep(0,m*block))
    
    # First position
    l <- 1000
    
    for(i in 1:m){ # For each block of CpG
        
        # Select a position of block
        v <- round(runif(1,0,1) *
                       (length(methInfo$start[methInfo$chr==chr[i]]) - block))
        
        methBlock <- getData(methInfo[methInfo$chr==chr[i]][(v+1):(v+block)])
        matProp <- sapply(posCTRL, function(x,methCur){unname(unlist(methCur[3*(x-1)+2]/methCur[3*(x-1)+1]))}, methCur=methBlock[5:length(methBlock)])
        res$chrOri[((i-1)*block + 1):(i*block)] <- rep(chr[i], block)
        res$startOri[((i-1)*block + 1):(i*block)] <- unname(unlist(methBlock[2]))
        res$start[((i-1)*block + 1):(i*block)] <- unname(unlist(methBlock[2])) - unname(unlist(methBlock[2]))[1] + l
        res$end[((i-1)*block + 1):(i*block)] <- res$start[((i-1)*block + 1):(i*block)]
        res$meanCTRL[((i-1)*block + 1):(i*block)] <- rowMeans(matProp)
        res$varCTRL[((i-1)*block + 1):(i*block)] <- apply(matProp, 1, var)
        l <- res$start[i*block] + 10000
    }
    res$alphaCTRL <- apply(res[,c(4,5)], 1, estBetaAlpha)
    res$betaCTRL <- apply(res[,c(4,5)], 1, estBetaBeta)
    res
}


#' @title TODO
#'
#' @description TODO
#'
#' @param x
#'
#' @param nb
#'
#' @param sDiff
#'
#' @param diffCase
#'
#' @param propDiffsd
#' -c * exp(b * log(distance between two CpG))
#'
#' @param b constant for distance
#'
#' @param endLength TODO
#'
#' @return TODO
#'
#' @examples
#'
#' ## TODO
#'
#'
#' @author Pascal Belleau
#' @keywords internal
# Identify the pos where the case are Diff meth and which one are heritable
# stateInfo contient en start les position des CpG (peut-etre le data.frame de getStateData)


getDiffCase <- function(x,nb, sDiff, diffCase, propDiffsd){
    meanDiff <- 0
    if(x[3] == 0){
        val <- rbeta(nb,estBetaAlpha(x[1:2]), estBetaBeta(x[1:2]))
        meanDiff <- x[1]
        partitionDiff <- c(0,nb)
    }
    else{
        meanDiff <- ifelse(x[1]<0.5,
                           min(1, x[1]+sDiff),
                           max(0, x[1]-sDiff))
        
        partitionDiff <- c(diffCase, nb -diffCase)
        
        
        val <- c(rbeta(partitionDiff[1], estBetaAlpha(c(meanDiff,x[2])),
                       estBetaBeta(c(meanDiff,x[2]))),
                 rbeta(partitionDiff[2],estBetaAlpha(x), estBetaBeta(x)))
    }
    #list(sDiff=sDiff, partitionDiff=partitionDiff, val=val)
    c(meanDiff, partitionDiff, val)
}


#' @title TODO
#'
#' @description TODO
#'
#' @param nbCtrl
#'
#' @param nbCase
#'
#' @param generation
#'
#' @param stateInfo
#'
#' @param stateDiff
#'
#' @param diffValue
#'
#' @param propDiff
#'
#' @param propDiffsd=0.1
#'
#' @param propInheritance
#'
#' @param propHetero
#'
#' @return TODO
#'
#' @examples
#'
#' ## TODO
#'
#'
#' @author Pascal Belleau
#' @importFrom msm rtnorm
#' @keywords internal
# Identify the pos where the case are Diff meth and which one are heritable
# stateInfo contient en start les position des CpG (peut-etre le data.frame de getStateData)


getSim <- function(nbCtrl,nbCase, generation, stateInfo, stateDiff, diffValue, propDiff, propDiffsd=0.1, propInheritance, propHetero){
    inR <- propDiff
    res<-list()
    if(propDiffsd < 0.0000001){
        diffCase <- round(nbCase * inR)
    } else{
        diffCase <- round(nbCase * rtnorm(1,mean = inR,sd = propDiffsd, lower=0,upper=1))
    }
    ctrl <- t(apply(stateInfo[,c(4,5)], 1, function(x,nb){rbeta(nb,estBetaAlpha(x), estBetaBeta(x))}, nb=nbCtrl))
    case <- t(apply(cbind(stateInfo[,c(4,5)], stateDiff$stateDiff), 1, getDiffCase, nb=nbCase, sDiff=diffValue, diffCase=diffCase, propDiffsd=propDiffsd))
    
    res[[1]] <- cbind(case[,c(1:3)], ctrl, case[,4:length(case[1,])])
    for(i in 2:generation)
    {
        rm(case,ctrl)
        inR <- propDiff * propInheritance^(i-2) # Ici decale d'une generation
        diffCur <- diffValue * propHetero       # Change diffValue
        if(propDiffsd < 0.0000001){
            diffCase <- round(nbCase * inR)
        } else{
            diffCase <- round(nbCase * rtnorm(1,mean = inR,sd = propDiffsd, lower=0,upper=1))
        }
        
        
        ctrl <- t(apply(stateInfo[,c(4,5)], 1, function(x,nb){
            rbeta(nb,estBetaAlpha(x), estBetaBeta(x))}, nb=nbCtrl))
        case <- t(apply(cbind(stateInfo[,c(4,5)], stateDiff$stateInherite), 1,
                        getDiffCase, nb=nbCase, sDiff=diffCur, diffCase=diffCase,
                        propDiffsd=propDiffsd))
        
        res[[i]] <- cbind(case[,c(1:3)], ctrl, case[,4:length(case[1,])])
    }
    res
}

#' @title TODO
#'
#' @description Identify the pos where the case are Diff meth and which one are heritable
#'
#' @param stateInfo
#'
#' @param rateDiff
#'
#' @param minRate
#'
#' @param propInherite
#'
#' @param c
#'
#' @param b
#'
#' @param endLength
#'
#' @return TODO
#'
#' @examples
#'
#' ## TODO
#' 
#' @author Pascal Belleau
#' @keywords internal

getDiffMeth <- function(stateInfo,rateDiff, minRate, propInherite, c=0.2, b=-2e-01, endLength=1000){
    
    nbPos <- length(stateInfo[,1])
    nbTry <- 1
    flag <-  TRUE
    
    while(nbTry < 1000 & flag){
        stateDiff <- rep(0, nbPos)
        stateInherite <- rep(0, nbPos)
        vExp <- rexp(nbPos, rateDiff)
        
        i<-1
        m<- max(1, round(vExp[i]))
        
        while(m <= nbPos){
            flagInherite <- ifelse(runif(1,0,1) < propInherite, TRUE,FALSE)
            stateDiff[m] <- 1
            if(flagInherite){
                stateInherite[m] <- 1
            }
            m <- m+1
            while(m <= nbPos && 
                  (stateInfo$start[m] - stateInfo$start[m-1]) <= endLength){
                cutOff <- -c * 
                    exp(b*log(stateInfo$start[m] - stateInfo$start[m-1]))
                u <- runif(1,0,1)
                if(u < cutOff){
                    stateDiff[m] <- 1
                    if(flagInherite){
                        stateInherite[m] <- 1
                    }
                }
                m <- m + 1
            }
            m <- m + round(vExp[i])
        }
        if(length(which(stateDiff == 1))>= minRate * nbPos){
            flag <- FALSE
        }
        nbTry <- nbTry + 1
    }
    if(flag){
        stateDiff <- NULL
    }
    list(stateDiff=stateDiff, stateInherite=stateInherite)
}



#' @title TODO
#'
#' @description TODO
#'
#' @param pathOut
#'
#' @param pref
#'
#' @param k
#'
#' @param nbCtrl
#'
#' @param nbCase
#'
#' @param treatment
#'
#' @param sample.id
#'
#' @param generation
#'
#' @param res
#'
#' @param rateDiff
#'
#' @param minRate
#'
#' @param propInherite
#'
#' @param diffValue
#'
#' @param propDiff
#'
#' @param propDiffsd
#'
#' @param propInheritance
#'
#' @param propHetero
#'
#' @param diffRes (default \code{NULL})
#'
#' @return TODO
#'
#' @examples
#'
#' ## TODO
#' 
#' @author Pascal Belleau
#' @importFrom methylKit read filterByCoverage normalizeCoverage unite calculateDiffMeth get.methylDiff getData tileMethylCounts methRead
#' @keywords internal
# Identify the pos where the case are Diff meth and which one are heritable
# stateInfo contient en start les position des CpG (peut-etre le data.frame de getStateData)

simInheritance <- function(pathOut, pref, k, nbCtrl, nbCase, treatment, sample.id, generation, res,
                   rateDiff ,minRate, propInherite, diffValue, propDiff,
                   propDiffsd, propInheritance, propHetero, diffRes = NULL){
    
    # Hard coded for the moment
    context="CpG"
    assembly="Rnor_5.0"
    meanCov <- 80
    if(is.null(diffRes)){
        diffRes <- getDiffMeth(stateInfo=res,
                               rateDiff=rateDiff, minRate=minRate,
                               propInherite=propInherite)
    }
    simV0.1 <- getSim(nbCtrl=nbCtrl, nbCase=nbCase, generation=generation,
                      stateInfo=res, stateDiff=diffRes , diffValue=diffValue,
                      propDiff=propDiff, propDiffsd=propDiffsd,
                      propInheritance=propInheritance, propHetero = propHetero)
    
    saveRDS(diffRes,file = paste0(pathOut, "/stateDiff_", pref , "_", k,".rds"))
    saveRDS(simV0.1,file = paste0(pathOut, "/simV0.1_", pref , "_", k,".rds"))
    
    
    
    myobj <- list()
    myMat <- list()
    myTr <- list()
    meth <- list()
    myDiff <- list()
    
    for(i in 1:generation){
        outList=list()
        for(j in 1:(nbCtrl+nbCase)){
            coverage <- rpois(length(res[,3]), meanCov)+1
            testM <- data.frame(chr=res[,1],start=res[,2],end=res[,3]
                                ,strand=rep("+",length(res[,3]))
                                ,coverage=coverage,numCs=round(coverage * simV0.1[[i]][,3+j]),numTs=round(coverage * (1-simV0.1[[i]][,3+j])))
            obj<-new("methylRaw",testM,sample.id=sample.id[[i]][[j]],assembly=assembly,
                     context=context,resolution='base')
            
            outList[[j]]<-obj
        }
        myMat[[i]] <- outList
        myTr[[i]] <- treatment
        myobj[[i]] <- new("methylRawList",outList,treatment=treatment)
        
        #meth=unite(myobj, destrand=FALSE)
        filtered.myobj <- filterByCoverage(myobj[[i]],lo.count=10,lo.perc=NULL,
                                           hi.count=NULL,hi.perc=99.9)
        filtered.myobj <- normalizeCoverage(filtered.myobj, "median")
        meth[[i]] <- unite(filtered.myobj, destrand=FALSE)
        myDiff[[i]] <- calculateDiffMeth(meth[[i]])
    }
    
    saveRDS(myobj,file = paste0(pathOut, "/methylObj_", pref , "_", k,".rds"))
    saveRDS(meth,file = paste0(pathOut, "/meth_", pref , "_", k,".rds"))
    saveRDS(myDiff,file = paste0(pathOut, "/methDiff_", pref , "_", k,".rds"))
    
}

#' @title TODO
#'
#' @description TODO
#'
#' @param pathOut
#'
#' @param pref
#'
#' @param k
#'
#' @param nbCtrl
#'
#' @param nbCase
#'
#' @param treatment
#'
#' @param sample.id
#'
#' @param generation
#'
#' @param res
#'
#' @param rateDiff
#'
#' @param minRate
#'
#' @param propInherite
#'
#' @param diffValue
#'
#' @param propDiff
#'
#' @param propDiffsd
#'
#' @param propInheritance
#'
#' @param propHetero
#'
#' @param diffRes (default \code{NULL})
#'
#' @return TODO
#'
#' @examples
#'
#' ## TODO
#' 
#' @author Pascal Belleau
#' @importFrom methylKit read filterByCoverage normalizeCoverage unite calculateDiffMeth get.methylDiff getData tileMethylCounts methRead
#' @keywords internal
# Identify the pos where the case are Diff meth and which one are heritable
# stateInfo contient en start les position des CpG (peut-etre le data.frame de getStateData)

restartSim <- function(pathOut, pref, k, nbCtrl, nbCase, treatment, sample.id, generation, res,
                   rateDiff ,minRate, propInherite, diffValue, propDiff,
                   propDiffsd, propInheritance, propHetero, diffRes = NULL){
    
    # Hard coded for the moment
    context="CpG"
    assembly="Rnor_5.0"
    meanCov <- 80
    if( file.exists(paste0(pathOut, "/stateDiff_", pref , "_", k,".rds"))){
        if( is.null(diffRes)){
            diffRes <- readRDS(paste0(pathOut, "/stateDiff_", pref , "_", k,".rds"))
        }
    } else{
        if( is.null(diffRes)){
            diffRes <- getDiffMeth(stateInfo=res,
                                   rateDiff=rateDiff, minRate=minRate,
                                   propInherite=propInherite)
        }
        saveRDS(diffRes,file = paste0(pathOut, "/stateDiff_", pref , "_", k,".rds"))
    }
    simV0.1 <- list()
    if( file.exists(paste0(pathOut, "/simV0.1_", pref , "_", k,".rds"))){
        simV0.1 <- readRDS(paste0(pathOut, "/simV0.1_", pref , "_", k,".rds"))
    } else{
        simV0.1 <- getSim(nbCtrl=nbCtrl, nbCase=nbCase, generation=generation,
                          stateInfo=res, stateDiff=diffRes , diffValue=diffValue,
                          propDiff=propDiff, propDiffsd=propDiffsd,
                          propInheritance=propInheritance, propHetero = propHetero)
        
        
        saveRDS(simV0.1,file = paste0(pathOut, "/simV0.1_", pref , "_", k,".rds"))
    }
    
    
    
    if(!(file.exists( paste0(pathOut, "/methylObj_", pref , "_", k,".rds")))
       | !(file.exists( paste0(pathOut, "/meth_", pref , "_", k,".rds"))) 
       | !(file.exists( paste0(pathOut, "/methDiff_", pref , "_", k,".rds"))) ){
        
        myobj <- list()
        myMat <- list()
        myTr <- list()
        meth <- list()
        myDiff <- list()
        
        for(i in 1:generation){
            outList=list()
            for(j in 1:(nbCtrl+nbCase)){
                coverage <- rpois(length(res[,3]), meanCov)+1
                testM <- data.frame(chr=res[,1],start=res[,2],end=res[,3]
                                    ,strand=rep("+",length(res[,3]))
                                    ,coverage=coverage,numCs=round(coverage * simV0.1[[i]][,3+j]),numTs=round(coverage * (1-simV0.1[[i]][,3+j])))
                obj<-new("methylRaw",testM,sample.id=sample.id[[i]][[j]],assembly=assembly,
                         context=context,resolution='base')
                
                outList[[j]]<-obj
            }
            myMat[[i]] <- outList
            myTr[[i]] <- treatment
            myobj[[i]] <- new("methylRawList",outList,treatment=treatment)
            
            #meth=unite(myobj, destrand=FALSE)
            filtered.myobj <- filterByCoverage(myobj[[i]],lo.count=10,lo.perc=NULL,
                                               hi.count=NULL,hi.perc=99.9)
            filtered.myobj <- normalizeCoverage(filtered.myobj, "median")
            meth[[i]] <- unite(filtered.myobj, destrand=FALSE)
            myDiff[[i]] <- calculateDiffMeth(meth[[i]])
        }
        
        saveRDS(myobj,file = paste0(pathOut, "/methylObj_", pref , "_", k,".rds"))
        saveRDS(meth,file = paste0(pathOut, "/meth_", pref , "_", k,".rds"))
        saveRDS(myDiff,file = paste0(pathOut, "/methDiff_", pref , "_", k,".rds"))
    }
}



