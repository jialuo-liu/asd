library(mvtnorm)

# ---- Part 0: Data Generation and Utilities ----

# ---- Numerical Integration Constants ----
# Constants for mvtnorm::pmvnorm and related functions
MAXPTS <- 1e5    # Maximum number of function evaluations
ABSEPS <- 1e-5   # Absolute error tolerance
RELEPS <- 0      # Relative error tolerance (0 = absolute only)
MVN_SEED <- 2025 # Random seed for reproducibility

# ---- GenerateDataFixedAllocation: Simulate data with fixed stage 2 sample sizes ----
#' Generates simulation data for pick-the-winner designs where stage 2 sample sizes
#' are predetermined and fixed across all treatment arms.
#'
#' @param nSim number of simulations
#' @param vN1 vector of stage 1 sample sizes (control, dose 1, dose 2, ...), 
#' can be a single number if all arms have the same sample size
#' @param vN2 vector of stage 2 sample sizes (control, dose 1, dose 2, ...), 
#' can be a single number if all arms have the same sample size
#' @param vMu vector of means (control, dose 1, dose 2, ...)
#' @param dSD standard deviation (assume same and known standard deviation for all arms)
#' @return A list containing: mZ1 (stage 1 Z-statistics matrix, nSim x nTreatments), 
#'   mZ2 (stage 2 Z-statistics matrix), mZ (pooled Z-statistics matrix), 
#'   vN1 (stage 1 sample sizes vector), mN2 (stage 2 sample sizes matrix)

GenerateDataFixedAllocation <- function(nSim, vN1, vN2, vMu, dSD){
  
  # Input validation
  if(nSim <= 0 || !is.finite(nSim)) stop("nSim must be a positive finite number")
  if(dSD <= 0 || !is.finite(dSD)) stop("dSD must be a positive finite number")
  if(any(vN1 <= 0) || any(!is.finite(vN1))) stop("All vN1 values must be positive and finite")
  if(any(vN2 <= 0) || any(!is.finite(vN2))) stop("All vN2 values must be positive and finite")
  if(any(!is.finite(vMu))) stop("All vMu values must be finite")
  
  if(length(vN1) == 1) vN1 <- rep(vN1, length(vMu))
  if(length(vN2) == 1) vN2 <- rep(vN2, length(vMu))
  
  if(length(vN1) != length(vMu)) stop("vN1 and vMu should have the same length.")
  if(length(vN2) != length(vMu)) stop("vN2 and vMu should have the same length.")
  
  # needed for naive pooling (need to recalculate the pooled Z stats)
  mXBar1 <- mvtnorm::rmvnorm( nSim, mean = vMu/dSD, sigma = diag(1/vN1))
  mXBar2 <- mvtnorm::rmvnorm( nSim, mean = vMu/dSD, sigma = diag(1/vN2))
  
  mZ1 <- t(t(mXBar1[,-1,drop=F] - mXBar1[,1])/sqrt(1/vN1[1]+1/vN1[-1]))
  mZ2 <- t(t(mXBar2[,-1,drop=F] - mXBar2[,1])/sqrt(1/vN2[1]+1/vN2[-1]))
  
  vN <- vN1 + vN2
  
  vMeanCtl <- (mXBar1[, 1]*vN1[ 1] + mXBar2[, 1]*vN2[ 1])/vN[ 1]
  mMeanTrt <- t((t(mXBar1[,-1,drop=F])*vN1[-1] + t(mXBar2[,-1,drop=F])*vN2[-1])/vN[-1] )
  
  mZ <- t(t(mMeanTrt - vMeanCtl)/sqrt(1/vN[1] + 1/vN[-1]))
  mN2 <- matrix( rep( vN2, nSim ), nSim, length(vN2), byrow = T)
  return( list(mZ1 = mZ1, mZ2 = mZ2, mZ = mZ, vN1=vN1, mN2=mN2) )
}

# ---- GenerateDataAdaptiveAllocation: Simulate with dynamic sample size reallocation ----
#' Generates simulation data for drop-the-loser designs where stage 2 sample sizes
#' are dynamically reallocated based on interim treatment selection.
#' Arms with Z-statistics below dZCut are dropped, and their sample sizes are
#' redistributed among continuing arms.
#'
#' @param nSim number of simulations
#' @param vN1 vector of stage 1 sample sizes
#' @param vN2 vector of original stage 2 sample sizes (before reallocation)
#' @param vMu vector of means (control, dose 1, dose 2, ...)
#' @param dSD standard deviation (assume same and known standard deviation for all arms)
#' @param dZCut Z-statistic cutoff for treatment selection (default 0)
#' @return A list containing: mZ1 (stage 1 Z-statistics matrix), mZ2 (stage 2 Z-statistics matrix, -Inf for dropped arms),
#'   mZ (pooled Z-statistics matrix, -Inf for dropped arms), vN1 (stage 1 sample sizes vector), 
#'   mN2 (reallocated stage 2 sample sizes matrix), mSelArm (logical matrix indicating which arms were selected)

GenerateDataAdaptiveAllocation <- function(nSim, vN1, vN2, vMu, dSD, dZCut = 0){
  
  # Input validation
  if(nSim <= 0 || !is.finite(nSim)) stop("nSim must be a positive finite number")
  if(dSD <= 0 || !is.finite(dSD)) stop("dSD must be a positive finite number")
  if(!is.finite(dZCut)) stop("dZCut must be finite")
  if(any(vN1 <= 0) || any(!is.finite(vN1))) stop("All vN1 values must be positive and finite")
  if(any(vN2 <= 0) || any(!is.finite(vN2))) stop("All vN2 values must be positive and finite")
  if(any(!is.finite(vMu))) stop("All vMu values must be finite")
  
  if(length(vN1) == 1) vN1 <- rep(vN1, length(vMu))
  if(length(vN2) == 1) vN2 <- rep(vN2, length(vMu))
  
  if(length(vN1) != length(vMu)) stop("vN1 and vMu should have the same length.")
  if(length(vN2) != length(vMu)) stop("vN2 and vMu should have the same length.")
  
  # Generate stage 1 data
  mXBar1 <- mvtnorm::rmvnorm(nSim, mean = vMu/dSD, sigma = diag(1/vN1))
  mZ1 <- (mXBar1[,-1,drop=F] - mXBar1[,1])/sqrt(1/vN1[1] + 1/vN1[-1])
  
  mSelArm <- (mZ1 >= dZCut)
  vN2New <- sum(vN2)/(rowSums(mSelArm) + 1)
  vN2New[rowSums(mSelArm) == 0] <- 0
  mN2 <- cbind(1, ifelse(mSelArm, 1, 0))* vN2New
  
  # Generate stage 2 data
  mMuS2 <- matrix(rep(vMu/dSD, nSim), nrow = nSim, byrow = TRUE)
  # Avoid division by zero: set sigma to Inf for dropped arms (mN2 == 0)
  mSigmaS2 <- sqrt(ifelse(mN2 == 0, 0, 1/mN2))
  mXBar2 <- mMuS2 + mSigmaS2 * matrix(rnorm(nSim * length(vMu)), nrow = nSim)
  
  # Calculate Z statistics for stage 1 and stage 2
  mZ2 <- (mXBar2[,-1,drop=F] - mXBar2[,1])/sqrt(1/mN2[,1] + 1/mN2[,-1])
  
  # Calculate pooled Z statistics
  mN <- matrix(rep(vN1, nSim), nrow = nSim, byrow = TRUE) + mN2
  vMeanCtl <- (mXBar1[,1] * vN1[1] + mXBar2[,1] * mN2[,1]) / mN[,1]
  mMeanTrt <- (mXBar1[,-1,drop=F] * vN1[-1] + mXBar2[,-1,drop=F] * mN2[,-1])/mN[,-1]
  
  mZ <- (mMeanTrt - vMeanCtl)/sqrt(1/mN[,1] + 1/mN[,-1])
  mZ2[!mSelArm] <- -Inf
  mZ[!mSelArm] <- -Inf
  return(list(mZ1 = mZ1, mZ2 = mZ2, mZ = mZ, vN1 = vN1, mN2 = mN2,
              mSelArm = mSelArm))
}

# ---- SubsetSimulationData: Extract subset of simulation iterations ----
#' Subsets simulation data to specific iterations (e.g., for block processing)
#'
#' @param lData list containing simulation data matrices (mZ1, mZ2, mZ, vN1, mN2, optionally mSelArm)
#' @param vInd vector of iteration indices to extract
#' @return A list with the same structure as lData but only containing the specified iterations

SubsetSimulationData <- function(lData, vInd){
  # Enhanced version that handles mSelArm for drop-the-loser scenarios
  if("mSelArm" %in% names(lData)){
    return( list(mZ1 = lData$mZ1[vInd,,drop=F],
                 mZ2 = lData$mZ2[vInd,,drop=F],
                 mZ  = lData$mZ[vInd,,drop=F],
                 vN1 = lData$vN1,
                 mN2 = lData$mN2[vInd,,drop=F],
                 mSelArm = lData$mSelArm[vInd,,drop=F])
    )
  }else{
    return( list(mZ1 = lData$mZ1[vInd,,drop=F],
                 mZ2 = lData$mZ2[vInd,,drop=F],
                 mZ  = lData$mZ[vInd,,drop=F],
                 vN1 = lData$vN1,
                 mN2 = lData$mN2[vInd,,drop=F]) )
  }
}

# ---- Part I: Methods ----
# ---- IsVector: check if an object is a vector ----
#' @param vZ an object
#' @return Logical value indicating if the object is a vector (TRUE) or not (FALSE)

IsVector <- function( vZ ){
  if( is.vector(vZ) ) {
    return(TRUE)
  }else if( is.matrix(vZ)){
    if(nrow(vZ) == 1){
      return(TRUE)
    }else{
      return(FALSE)
    }
  }else{
    return(FALSE)
  }
}

# ---- BonfPvalue: Get Bonferroni P-value  ----
#' @param vZ a vector of Z statistics
#' @return Adjusted p-value using Bonferroni correction (min of nArm * min(p-values), 1)

BonfPvalue <- function( vZ ){
  
  if( !IsVector(vZ) ) stop("Please provide a vector.")
  
  if( length(vZ) == 0 ) {
    dPval <- 1
  }else{
    nArm <- length(vZ)
    dPval <- nArm * min( pnorm(vZ, lower.tail = F) )
  }
  
  return( min(dPval,1) )
}

# ---- DunnPvalue: Get Dunnett P-value  ----
#' @param vZ a vector, Z test statistics for each active arm
#' @param vN a vector, if not provided or length = 1, all arms have equal sample
#'  size; if length > 1, sample sizes (or relative proportion) for each arm, 
#'  where the first element represents the control arm.
#' @return Adjusted p-value using Dunnett's method for multiple comparisons with control

DunnPvalue <- function( vZ, vN = NULL ){
  
  if( !IsVector(vZ) ) stop("Please provide a vector.")
  
  if( length(vZ) == 0 ) {
    
    dPval <- 1
    
  }else{
    
    nArm <- length(vZ)
    dZMax <- max( vZ )
    
    if( is.null(vN) ) 
      message("Dunnett Test: sample size per arm not provided, assumed equal sample size.")
    if( is.null(vN) || length( unique(vN) ) == 1 ){
      
      Pfun <- function(x) ( pnorm( dZMax * sqrt( 2 ) + x ) )^nArm * dnorm(x)
      dPval <- 1 - stats::integrate(Pfun, -Inf, Inf)$value
      
    }else{
      if(length(vN) - 1 != length(vZ) )
        stop("The length of vN should equal the total number of arms, including control.")
      vRR <- vN[-1]/(vN[-1] + vN[1])
      
      if( length( unique(vRR) ) == 1 ){
        Pfun <- function(x) ( pnorm( ( dZMax + sqrt(vRR[1]) * x )/sqrt(1-vRR[1]) ) )^nArm * dnorm(x)
        dPval <- 1 - stats::integrate(Pfun, -Inf, Inf)$value
      }else{
        mCov <- sqrt ( outer( vRR, vRR ) )
        diag(mCov) <- 1
        dPval <- 1 - mvtnorm::pmvnorm( upper = rep(dZMax, nArm),
                                       sigma = mCov,
                                       algorithm = GenzBretz(maxpts = MAXPTS, 
                                                             abseps = ABSEPS,
                                                             releps = RELEPS),
                                       seed = MVN_SEED)
      }
      
    }
    
  }
  
  return( min(dPval,1) )
}

# ---- SimesPvalue: Get Simes P-value  ----
#' @param vZ a vector of Z statistics
#' @return Adjusted p-value using Simes' method for testing intersection hypotheses

SimesPvalue <- function( vZ ){
  if( !IsVector(vZ) ) stop("Please provide a vector.")
  
  if( length(vZ) == 0 ) {
    dPval <- 1
  }else{
    nArm <- length(vZ)
    dPval <- nArm * min( sort( pnorm(vZ, lower.tail = F) ) / (1:nArm) )
  }
  
  return( min(dPval,1) )
}

# ---- GetPval: Wrapper function to get p-value  ----
#' @param vZ a vector of Z statistics
#' @param vN a vector, if not provided or length = 1, all arms have equal sample 
#' size; if length > 1, sample sizes (or relative proportion) for each arm, 
#' where the first element represents the control arm. Only required for Dunnett 
#' test
#' @param vMethod a vector of multiple testing procedures, "Bonf" = Bonferroni,
#'  "Simes" = Simes, "Dunn" = Dunnett method.
#' @return Data frame containing p-values for each requested method

GetPval <- function(vZ, vN = NULL, vMethod = c("Bonf", "Simes", "Dunn")){
  if( !IsVector(vZ) ) stop("Please provide a vector.")
  
  dfPval <- NULL
  if("Bonf" %in% vMethod){
    dfPval$dPvalBonf  <-  BonfPvalue(vZ)
  }
  
  if("Simes" %in% vMethod){
    dfPval$dPvalSimes <- SimesPvalue(vZ)
  }
  
  if("Dunn" %in% vMethod){
    dfPval$dPvalDunn  <-  DunnPvalue(vZ, vN = vN)
  }
  
  return( as.data.frame( dfPval ) )
}

# ---- CTPTest: Closed Testing Procedure----
#' Performs closed testing procedure for multiple comparison tests
#'
#' @param vZ a vector, Z values for each treatment arm.
#' @param vSelArm a logical vector
#' @param vN a vector, if not provided or length = 1, all arms have equal sample
#'  size; if length > 1, sample sizes (or relative proportion) for each arm, 
#'  where the first element represents the control arm. Only required for 
#'  Dunnett test
#' @param vMethod a vector of multiple testing procedures, "Bonf" = Bonferroni,
#'  "Simes" = Simes, "Dunn" = Dunnett method.
#' @return A list containing: vPval (vector of minimum p-values across methods),
#'   dfTestK (data frame with results for all intersection hypotheses), 
#'   dfTest (data frame with adjusted p-values for each treatment)

CTPTest <- function( vZ, vSelArm, vN=NULL, vMethod = c("Bonf","Simes","Dunn")){
  
  nArm <- length(vZ)
  
  if( !IsVector(vZ) ) stop("Please provide a vector.")
  if( missing(vSelArm) || is.null(vSelArm) ) vSelArm <- rep(T, nArm)
  if( length(vSelArm) != nArm ) stop("vSelArm must have the same length as vZ.")
  
  if( length(vN) == 1) vN <- rep(vN, nArm + 1)
  
  if( !is.null(vN) & length(vN) != nArm + 1) stop("vN should have the same length as the total number of arms, including control")
  
  if( nArm == 1 ) {
    dfTestK <- GetPval(vZ = vZ, vN = vN, vMethod = vMethod)
    vPval <- sapply( dfTestK , min )
    dfTestK$nTrt <- 1
    dfTest <- dfTestK
  }else{
    # each hypothesis either included in the intersection hypothesis or not
    # There are a total of 2^(nArm) - 1 intersection hypotheses
    dfHypo <- expand.grid( rep( list(c(TRUE,FALSE)) , nArm ) )
    dfHypo <- dfHypo[rowSums(dfHypo) != 0, ] # remove H_K where K is empty
    colnames(dfHypo) <- paste0("H", 1:nArm)
    
    dfHypo <- as.data.frame(dfHypo)
    
    dfTestK <- NULL
    # for each intersection hypothesis
    for(i in 1:nrow(dfHypo)){
      vTrt <- unlist(dfHypo[i, 1:nArm])
      # get the pvalue for the H_{K\cap K_j}
      dfTestKTmp <- GetPval(vZ = vZ[vTrt&vSelArm],
                            vN = vN[c(T,vTrt&vSelArm)],
                            vMethod = vMethod)
      dfTestK <- rbind(dfTestK, dfTestKTmp)
    }
    dfTestK <- as.data.frame(dfTestK)
    
    dfTest <- NULL
    for(j in 1:nArm){
      dfContain <- dfTestK[dfHypo[,paste0("H",j)],,drop=FALSE]
      if(nrow(dfContain) == 0) next
      dfTestTmp <- apply(dfContain, 2, max)
      dfTest <- rbind(dfTest, dfTestTmp)
    }
    dfTest <- as.data.frame(dfTest)
    vPval <- sapply( dfTest , min )
    
    dfTestK <- cbind(dfHypo, dfTestK)
    
    dfTest$nTrt <- 1:nArm
    
    rownames(dfTestK) <- rownames(dfTest) <- NULL
  }
  
  return(
    list(
      vPval    = vPval,
      dfTestK = dfTestK,
      dfTest  = dfTest)
  )
  
}

# ---- Method 1: Second Stage Only ----
# ---- S2Only: 2nd Stage only analysis ----
#' Analyzes only the second stage data using closed testing procedure
#'
#' @param lData a list, data by stage and pooled
#' @param mSelArm a logical matrix, each row representing whether an arm is 
#' selected or not, e.g., c(T, F) if dose 1 is selected and dose 2 is not
#' @param vMethod a vector of multiple testing procedures, "Bonf" = Bonferroni,
#'  "Simes" = Simes, "Dunn" = Dunnett method.
#' @return A list containing: dfPval (data frame with minimum p-values per iteration),
#'   dfTest (data frame with detailed test results for all arms and iterations)

S2Only <- function( lData, mSelArm,
                    vMethod = c("Bonf","Simes","Dunn")){
  
  if( ! identical( dim(lData$mZ2), dim(mSelArm)) )
    stop("The dimension of mZ2 and mSelArm do not match!")
  
  mZ2 <- lData$mZ2
  nK2 <- ncol(mZ2)
  nSim <- nrow(mZ2)
  mN2 <- lData$mN2
  
  # for drop the loser rule, some trials may not even have stage 2 data
  # for those trials, the p-values are set as 1.
  # only summarize trials with stage 2 data
  vCont <- which(rowSums(mSelArm) != 0)
  lTest <- lapply(seq_len(nSim), function(i) {
    if(sum(mSelArm[i,]) == 0){
      dfTestK <- as.data.frame( matrix(1, nK2, length(vMethod)) )
      colnames(dfTestK) <- paste0("dPval", vMethod)
      dfTestK$nTrt <- seq_len(nK2)
      dfTestK$nIter <- i
      dfTestK
    }else{
      dfTestK <- CTPTest(vZ = mZ2[ i , mSelArm[i,] ],
                         # only use Z stats/sample size of selected arms
                         vN = mN2[ i , c(T,mSelArm[i,]) ],
                         vMethod = vMethod)$dfTest
      TrtMap <- stats::setNames( (1:nK2)[mSelArm[i,]], 1:sum(mSelArm[i,]) )
      dfTestK$nTrt <- TrtMap[dfTestK$nTrt]
      dfTestKU <- as.data.frame( matrix( 1, sum(!mSelArm[i,] ), length(vMethod) + 1) )
      colnames( dfTestKU ) <- colnames( dfTestK )
      dfTestKU$nTrt <- (1:nK2)[!mSelArm[i,]]
      dfTestK <- rbind(dfTestK, dfTestKU)
      dfTestK <- dfTestK[order(dfTestK$nTrt),]
      dfTestK$nIter <- i
      dfTestK
    }
  } )
  dfTest <- as.data.frame( do.call(rbind,lTest) )
  dfTmp <- dfTest[,c("nIter",paste0("dPval",vMethod))]
  dfPval <- stats::aggregate(. ~ nIter, data = dfTmp, FUN = min)
  return( list( dfPval = dfPval,
                dfTest = dfTest) )
}

# ---- Method 2: Combination Function ----
# ---- PvalComb: Combine P value across stages ----
#' Combines p-values from two stages using inverse normal method with weights
#' This function also apply to matrices, treating each matrix as a vector and 
#' return a matrix with same dimension
#' @param vPval1 a vector, first stage p-values.
#' @param vPval2 a vector, second stage p-values.
#' @param vW     a vector, information fractions at first stage; 2nd stage is 1- vW.
#' @return A vector of combined p-values

PvalComb <- function( vPval1, vPval2, vW){
  if( length(vW) != 1 & var( c( length(vPval1), length(vPval2), length(vW)) ) != 0 )
    stop("The dimensions of vW, vPval1, and vPval2 do not match.")
  vPval <- 1 - pnorm( sqrt(vW) * qnorm( 1 - vPval1 ) + sqrt(1-vW) * qnorm( 1 - vPval2 ) )
  
  return( vPval )
}

# ---- CombFunInvN: Combination Function Approach by Weighted InvNorm ----
#' Combines p-values across stages using inverse normal combination with sample size weighting
#'
#' @param lData a list, data by stage and pooled
#' @param mSelArm a logical matrix, each row representing whether an arm is 
#' selected or not, e.g., c(T, F) if dose 1 is selected and dose 2 is not
#' @param vMethod a vector of multiple testing procedures, "Bonf" = Bonferroni,
#'  "Simes" = Simes, "Dunn" = Dunnett method.
#' @param strW a character, "ctrl" = by control sample size per stage, 
#' "total" = by total sample size per stage
#' @return A list containing: dfPval (data frame with minimum p-values per iteration),
#'   dfTest (data frame with detailed test results), dfTestK (all intersection hypotheses results)

CombFunInvN <- function( lData, mSelArm,
                         vMethod = c("Bonf","Simes","Dunn"),
                         strW){
  
  if( ! identical( dim(lData$mZ), dim(mSelArm)) )
    stop("The dimensions of mZ and mSelArm do not match!")
  
  mZ1 <- lData$mZ1
  mZ2 <- lData$mZ2
  vN1 <- lData$vN1
  
  nSim <- nrow(mZ1)
  nK1 <- ncol(mZ1)
  
  # Adjust sample size
  mN1 <- matrix( rep( vN1, nSim ), nSim, length(vN1), byrow = T)
  mN2 <- lData$mN2
  mN2[!cbind(TRUE,mSelArm)] <- 0
  mN <- mN1 + mN2
  
  if(strW == "ctrl"){
    vW <- mN1[,1] / ( mN[,1] )
  }else if(strW == "total"){
    vW <- rowSums(mN1) / ( rowSums(mN) )
  }else if(strW == "obs"){ 
    vW <- rowSums(mN1[,cbind(T,mSelArm),drop = F])/rowSums(mN[,cbind(T,mSelArm),drop = F])
  }else{
    stop("Please provide a valid strW, e.g., 'ctrl', 'total' or 'obs'.")
  }
  
  lTestK1 <- lapply(1:nSim, function(i) {
    dfTestK <- CTPTest(mZ1[i,],
                       vN = mN1[i,])$dfTestK
    dfTestK$nIter <- i
    dfTestK$dW <- vW[i]
    dfTestK
  } )
  dfTestK1 <- do.call(rbind, lTestK1)
  
  lTestK2 <- lapply(1:nSim, function(i) {
    dfTestK <- CTPTest(mZ2[i,],
                       vSelArm = mSelArm[i,], 
                       vN = mN2[i,])$dfTestK
    dfTestK$nIter <- i
    dfTestK
  } )
  
  dfTestK2 <- do.call(rbind, lTestK2)
  
  dfTestK <- merge(dfTestK1,
                   dfTestK2, all = T, by = c(paste0("H",1:nK1),"nIter"),
                   suffixes = c(1,2))
  
  # combine p-values across stages
  for(strMethod in vMethod){
    dfTestK[,paste0("dPval",strMethod)] <- PvalComb(
      vPval1 = dfTestK[,paste0("dPval",strMethod, 1)],
      vPval2 = dfTestK[,paste0("dPval",strMethod, 2)],
      vW     = dfTestK[,"dW"]
    )
    dfTestK[,paste0("dPval",strMethod, 1)] <- NULL
    dfTestK[,paste0("dPval",strMethod, 2)] <- NULL
  }
  
  # Handle drop-the-loser: set p-values to 1 when all arms dropped
  vStop <- which( rowSums(mSelArm) == 0 )
  dfTestK[dfTestK$nIter %in% vStop, paste0("dPval", vMethod)] <- 1
  
  dfTest <- NULL
  for(j in 1:nK1){
    
    dfContain <- dfTestK[dfTestK[,paste0("H",j)],,drop=FALSE]
    if(nrow(dfContain) == 0) next
    dfContain <- dfContain[,c("nIter",paste0("dPval",vMethod))]
    dfTestTmp <- stats::aggregate(. ~ nIter, data = dfContain, FUN = max)
    dfTestTmp$nTrt <- j
    dfTest <- rbind(dfTest, dfTestTmp)
    
  }
  
  dfTmp <- dfTest[,c("nIter",paste0("dPval",vMethod))]
  dfPval <- stats::aggregate(. ~ nIter, data = dfTmp, FUN = min)
  
  return( list(dfPval = dfPval,
               dfTestK = dfTestK,
               dfTest = dfTest) )
  
}

# ---- Method 3: Marginal Combination Function ----
# ---- MargPInvN: Marginal P-value combination function approach ----
#' Combines marginal p-values across stages using inverse normal method for each arm separately
#'
#' @param lData a list, data by stage and pooled
#' @param mSelArm a logical matrix, each row representing whether an arm is 
#' selected or not, e.g., c(T, F) if dose 1 is selected and dose 2 is not
#' @param vMethod a vector of multiple testing procedures, "Bonf" = Bonferroni
#' @param strW weighting scheme: "ctrl", "total", or "obs"
#' @return A list containing: dfPval (minimum p-values per iteration), 
#'   dfTestK (all intersection hypothesis results), dfTest (per-treatment results)

MargPInvN <- function( lData, mSelArm,
                       vMethod = c("Bonf"),
                       strW = "total"){
  
  if( ! identical( dim(lData$mZ), dim(mSelArm)) )
    stop("The dimensions of mZ and mSelArm do not match!")
  
  mZ1 <- lData$mZ1
  mZ2 <- lData$mZ2
  vN1 <- lData$vN1
  mZ2[!mSelArm] <- -Inf
  nSim <- nrow(mZ1)
  nK1 <- ncol(mZ1)
  
  # Adjust sample size
  mN1 <- matrix( rep( vN1, nSim ), nSim, length(vN1), byrow = T)
  mN2 <- lData$mN2
  mN2[!cbind(TRUE,mSelArm)] <- 0
  mN <- mN1 + mN2
  
  mPval1 <- 1 - pnorm(mZ1)
  mPval2 <- 1 - pnorm(mZ2)
  
  if(strW == "ctrl"){
    mW <- matrix( rep( mN1[,1] / mN[,1], nK1 ),
                  nrow = nSim, ncol = nK1, byrow = F)
  }else if(strW == "total"){
    mW <- matrix( rep( rowSums(mN1) / rowSums(mN), nK1 ),
                  nrow = nSim, ncol = nK1, byrow = F)
  }else if(strW == "obs"){
    vWTmp <- rowSums(mN1[,cbind(T,mSelArm),drop = F])/rowSums(mN[,cbind(T,mSelArm),drop = F])
    mW <- matrix( rep( vWTmp, nK1 ),
                  nrow = nSim, ncol = nK1, byrow = F)
  }else{
    stop("Please provide a valid strW, e.g., 'ctrl', 'total' or 'obs'.")
  }
  
  mPval <- PvalComb(
    vPval1 = mPval1,
    vPval2 = mPval2,
    vW     = mW
  )
  
  # Manually set p-value to 1 for non-selected arms
  mPval[!mSelArm] <- 1
  lTest <- lapply( 1:nSim, function(i) {
    dfTestK <- CTPTest(vZ = qnorm(1-mPval)[i, ], 
                       vN = mN[i, ], 
                       vMethod = vMethod )$dfTest
    dfTestK$nIter <- i
    dfTestK
  })
  
  dfTest <- as.data.frame( do.call(rbind, lTest) )
  
  dfTmp <- dfTest[,c("nIter",paste0("dPval",vMethod))]
  dfPval <- stats::aggregate(. ~ nIter, data = dfTmp, FUN = min)
  
  return( list( dfPval = dfPval,
                dfTest = dfTest ) )
  
}

# ---- Method 4: Conditional Error Function ----
# ---- CalculateDunnettCriticalValue: Critical value of the classical Dunnett test ----
#' Calculates the critical value for Dunnett's test at a given alpha level
#'
#' @param dAlpha alpha
#' @param vN sample sizes per arm, where the first element represents the 
#' control arm.
#' @return Critical value (dDs) for Dunnett's test

CalculateDunnettCriticalValue <- function( dAlpha, vN){
  
  # Input validation
  if(dAlpha <= 0 || dAlpha >= 1) stop("dAlpha must be between 0 and 1")
  if(length(vN) < 2) stop("vN must have at least 2 elements (control + 1 treatment)")
  if(any(vN <= 0) || any(!is.finite(vN))) stop("All vN values must be positive and finite")
  
  # make sure random seed won't change after this function
  if (exists(".Random.seed", envir = .GlobalEnv)) {
    old_seed <- .Random.seed
    # Restore the seed when function exits (normal or error)
    on.exit(.Random.seed <<- old_seed)
  } else {
    # If no seed existed before, remove it on exit
    on.exit({
      if (exists(".Random.seed", envir = .GlobalEnv)) {
        rm(.Random.seed, envir = .GlobalEnv)
      }
    })
  }
  
  vRR <- vN[-1]/(vN[-1] + vN[1])
  
  if( length( unique(vRR) ) == 1 ){
    
    Afun <- function(z) {
      Pfun <- function(x) ( pnorm( ( z + sqrt(vRR[1]) * x )/sqrt(1-vRR[1]) ) )^length(vRR) * dnorm(x)
      1 - stats::integrate(Pfun, -Inf, Inf)$value - dAlpha
    }
    dDs <- stats::uniroot(Afun, c(-1e5,1e5))$root
    
  }else{
    mCorr <- sqrt ( outer( vRR, vRR ) )
    diag(mCorr) <- 1
    dDs <- mvtnorm::qmvnorm( 1-dAlpha, sigma = mCorr, tail = "lower",
                             algorithm = GenzBretz(maxpts = MAXPTS, 
                                                   abseps = ABSEPS, 
                                                   releps = RELEPS),
                             seed = MVN_SEED)$quantile
    
  }
  
  attributes(dDs) <- NULL
  
  return( dDs )
}

# ---- DunnCondError: Calculate Conditional Error ----
#' Calculates conditional error function for Dunnett's test given stage 1 data
#'
#' @param vN1 a vector of total sample sizes per arm (control, dose 1, dose 2, ...)
#' @param vN2 a vector of stage 2 sample sizes per arm (control, dose 1, dose 2, ...)
#' @param dDs critical value
#' @param vZ1 a vector of stage 1 Z test statistics (dose 1, dose 2, ...)
#' @return Conditional error (probability of rejecting at least one hypothesis)

DunnCondError <- function( vN1, vN2, dDs, vZ1 ){
  
  if( var( c(length(vN1)-1, length(vN2)-1, length(vZ1)) ) !=0 )
    stop("The length of vN1, vN2, vZ1, vZ2 should match.")
  
  if( length( unique( round( vN1 / vN2, 6) ) ) != 1 )
    stop("Constant randomization ratio is required across stages.")
  
  vN <- vN1 + vN2
  
  dInfoFrac <- sum(vN1)/sum(vN)
  
  nArm <- length(vZ1)
  
  vRR <- vN[-1]/(vN[-1] + vN[1])
  
  if( length( unique(vRR) ) == 1){
    
    Pfun <- function(x) {
      dInt <- 1
      for(j in 1:length(vZ1)){
        dInt <- dInt * pnorm( ( dDs - sqrt(dInfoFrac) * vZ1[j] )/sqrt( (1-vRR[1])*(1-dInfoFrac) ) + sqrt(vRR[1]/(1-vRR[1])) * x)
      }
      dInt * dnorm(x)
    }
    dCondErr <- 1 - integrate(Pfun, -Inf, Inf)$value
    
  }else{
    
    mCorr <- sqrt ( outer( vRR, vRR ) )
    diag(mCorr) <- 1
    mCov <- mCorr * (1-dInfoFrac)
    
    dCondErr <- 1 - mvtnorm::pmvnorm( upper = rep(dDs, nArm),
                                      mean =  sqrt(dInfoFrac)*vZ1,
                                      sigma = mCov,
                                      algorithm = GenzBretz(maxpts = MAXPTS, 
                                                            abseps = ABSEPS, 
                                                            releps = RELEPS),
                                      seed = MVN_SEED)
    attributes(dCondErr) <- NULL
    
  }
  
  return( dCondErr )
  
}

# ---- CondDunnPvalue: Conditional second-stage Dunnett P-value  ----
#' Calculates conditional Dunnett p-value combining stage 1 and stage 2 data
#'
#' @param vZ1 a vector of stage 1 Z test statistics (dose 1, dose 2, ...)
#' @param vZ2 a vector of stage 2 Z test statistics (dose 1, dose 2, ...)
#' @param vN1 a vector of total sample sizes per arm (control, dose 1, dose 2, ...)
#' @param vN2 a vector of stage 2 sample sizes per arm (control, dose 1, dose 2, ...)
#' @return Conditional Dunnett p-value

CondDunnPvalue <- function( vZ1, vZ2, vN1, vN2 ){
  if( var( c(length(vN1)-1, length(vN2)-1, length(vZ1), length(vZ2)) ) !=0 )
    stop("The length of vN1, vN2, vZ1, vZ2 should match.")
  
  if( length( unique( round( vN1 / vN2, 6) ) ) != 1 )
    stop("Constant randomization ratio is required across stages.")
  
  if(length(vZ1) == 0){
    dPvalue <- 1
  }else{
    vN <- vN1 + vN2
    dInfoFrac <- vN1[1]/vN[1] # should be equivalent to sum(vN1)/sum(vN)
    vZ <- sqrt(dInfoFrac)*vZ1 + sqrt(1-dInfoFrac)*vZ2
    
    dZMax <- max(vZ)
    
    nArm <- length(vZ)
    
    vRR <- vN[-1]/(vN[-1] + vN[1])
    if( length(unique(vRR)) == 1 ){
      Pfun <- function(x) {
        dInt <- 1
        for(j in 1:length(vZ1)){
          dInt <- dInt * pnorm( ( dZMax - sqrt(dInfoFrac) * vZ1[j] )/sqrt( (1-vRR[1])*(1-dInfoFrac) ) + sqrt(vRR[1]/(1-vRR[1])) * x)
        }
        dInt * dnorm(x)
      }
      dPvalue <- 1 - integrate(Pfun, -Inf, Inf)$value
    }else{
      mCorr <- sqrt ( outer( vRR, vRR ) )
      diag(mCorr) <- 1
      mCov <- mCorr * (1-dInfoFrac)
      dPvalue <- 1 - mvtnorm::pmvnorm( upper = rep(dZMax, nArm),
                                       mean =  sqrt(dInfoFrac)*vZ1,
                                       sigma = mCov,
                                       algorithm = GenzBretz(maxpts = MAXPTS, 
                                                             abseps = ABSEPS, 
                                                             releps = RELEPS),
                                       seed = MVN_SEED)
      attributes(dPvalue) <- NULL
    }
    
  }
  
  return( dPvalue )
}

# ---- AdaptDunnTest: Adaptive Dunnett Tests, Closure Principle ----
#' Performs adaptive Dunnett test using conditional error function with closed testing
#'
#' @param lData a list, data by stage and pooled
#' @param mSelArm a logical matrix, each row representing whether an arm is 
#' selected or not, e.g., c(T, F) if dose 1 is selected and dose 2 is not
#' @param dAlpha alpha
#' @param dEps a small value for finding root
#' @return A list containing: dfPval (minimum conditional p-values per iteration),
#'   dfTest (per-treatment test results), dfTestK (all intersection hypothesis results)

AdaptDunnTest <- function( lData, mSelArm, dAlpha, dEps = 1e-6){
  
  if( ! identical( dim(lData$mZ), dim(mSelArm)) )
    stop("The dimensions of mZ and mSelArm do not match!")
  
  mZ1 <- lData$mZ1
  mZ2 <- lData$mZ2
  vN1 <- lData$vN1
  mZ2[!mSelArm] <- -Inf
  nSim <- nrow(mZ1)
  nK1 <- ncol(mZ1)
  
  # Adjust sample size
  mN1 <- matrix( rep( vN1, nSim ), nSim, length(vN1), byrow = T)
  mN2 <- lData$mN2
  mN2[!cbind(TRUE,mSelArm)] <- 0
  mN <- mN1 + mN2
  
  dfHypo <- expand.grid( rep( list(c(TRUE,FALSE)) , nK1 ) )
  colnames(dfHypo) <- paste0("H", 1:nK1)
  dfHypo <- dfHypo[rowSums(dfHypo) != 0, ]
  dfHypo <- as.data.frame(dfHypo)
  
  dfTestK <- NULL
  for(i in 1:nrow(dfHypo)){
    
    vTrt <- unlist(dfHypo[i, 1:nK1])
    
    # get Ds outside the loop
    vN2Total <- rowSums(mN2) 
    vN2U <- unique(vN2Total)
    dfDs <- NULL
    for( j in 1:length(vN2U)){
      # add exception for those trials without stage 2
      if(vN2U[j] == 0){
        dfTmp <- data.frame(dDs = Inf, vN2Total = 0)
      }else{
        vN2Ori <- vN2U[j] * vN1/sum(vN1)
        dDs <- CalculateDunnettCriticalValue( dAlpha, vN1[c(T,vTrt)] + vN2Ori[c(T,vTrt)])
        dfTmp <- data.frame(dDs = dDs, vN2Total = vN2U[j])
      }
      dfDs <- rbind(dfDs, dfTmp)
    }
    
    dfDs <- merge(data.frame(nIter = seq_len(nSim), vN2Total = vN2Total), dfDs,
                  by = "vN2Total",all.x = T)
    dfDs <- dfDs[order(dfDs$nIter), ]
    
    lTestKTmp <- lapply(seq_len(nSim), function(nIter) {
      vTrtSub <- vTrt & mSelArm[nIter,]
      if(sum(vTrtSub) == 0) {
        dCondPval <- 1
        dCondErr <- 0
      }else{
        # using selected arm only to get 2nd stage p-value
        dCondPval <- CondDunnPvalue( vZ1 = mZ1[nIter,vTrtSub],
                                     vZ2 = mZ2[nIter,vTrtSub],
                                     vN1 = vN1[c(T,vTrtSub)],
                                     vN2 = mN2[nIter, c(T,vTrtSub)])
        # conditional error function based on all first stage data,
        # including those from unselected arms
        
        dCondErr <- DunnCondError( vN1 = vN1[c(T,vTrt)],
                                   vN2 = vN2Ori[c(T,vTrt)],
                                   dDs = dfDs[nIter, "dDs"],
                                   vZ1 = mZ1[nIter,vTrt])
        
      }
      bCondRej <- ifelse(dCondPval <= dCondErr, TRUE, FALSE)
      data.frame(nIter     = nIter,
                 vCondPval = dCondPval,
                 vCondErr  = dCondErr ,
                 vCondRej  = bCondRej)
      
    } )
    dfTestKTmp <- do.call(rbind, lTestKTmp)
    mSelSub <- mSelArm
    mSelSub[,!vTrt] <- FALSE
    colnames(mSelSub) <- paste0("H",1:nK1)
    
    dfTestKTmp <- cbind( dfTestKTmp, mSelSub )
    
    dfTestK <- rbind(dfTestK, dfTestKTmp)
  }
  
  dfTest <- NULL
  for(j in seq_len(nK1)){
    dfContain <- dfTestK[dfTestK[,paste0("H",j)],,drop=FALSE]
    if(nrow(dfContain) == 0) next
    dfTestTmp <- stats::aggregate(vCondRej ~ nIter, data = dfContain, FUN = all)
    vDrop <- base::setdiff( seq_len(nSim), dfTestTmp$nIter)
    if(length(vDrop) > 0){
      dfTestTmpU <- data.frame( nIter = vDrop, vCondRej = FALSE)
      dfTestTmp <- rbind(dfTestTmp, dfTestTmpU)
    }
    dfTestTmp <- dfTestTmp[order(dfTestTmp$nIter),]
    dfTestTmp$nTrt <- j
    dfTest <- rbind(dfTest, dfTestTmp)
  }
  
  dfRej <- stats::aggregate(vCondRej ~ nIter, data = dfTest, FUN = any)
  
  return( list(dfRej = dfRej,
               dfTestK = dfTestK,
               dfTest = dfTest) )
}

# ---- Method 5: Partial Conditional Error Function ----
# ---- InvNCondErr: Conditional Error based on Weighted Inverse normal Combination Function ----
#' Calculates conditional error using inverse normal combination function
#' This function also apply to matrices, treating each matrix as a vector and 
#' return a matrix with same dimension.
#' @param dAlpha alpha
#' @param vW     a vector, information fractions.
#' @param vZ1 a vector of stage 1 Z test statistics
#' @return Conditional error probabilities

InvNCondErr <- function(dAlpha, vW, vZ1){
  if( length(vW) != 1 & var( c( length(vZ1), length(vW)) ) != 0 )
    stop("The dimensions of vW, vPval1, and vPval2 do not match.")
  
  1-pnorm((qnorm(1-dAlpha) - sqrt(vW)*vZ1)/sqrt(1-vW))
}


# ---- BonfPartCondError: Partial Conditional Error for Bonferroni ----
#' Performs Bonferroni test using partial conditional error function
#'
#' @param lData a list, data by stage and pooled
#' @param mSelArm a logical matrix, each row representing whether an arm is
#'  selected or not, e.g., c(T, F) if dose 1 is selected and dose 2 is not
#' @param dAlpha alpha
#' @param strW a character, "ctrl" = by control sample size per stage, 
#' "total" = by total sample size per stage
#' @param dEps a small value for finding root
#' @return A list containing: dfPval (minimum conditional error/p-values per iteration),
#'   dfTest (per-treatment test results), dfTestK (all intersection hypothesis results)

BonfPartCondError <- function( lData, mSelArm, dAlpha, strW, dEps = 1e-6){
  
  if( ! identical( dim(lData$mZ), dim(mSelArm)) )
    stop("The dimensions of mZ and mSelArm do not match!")
  
  mZ1 <- lData$mZ1
  mZ2 <- lData$mZ2
  vN1 <- lData$vN1
  mZ2[!mSelArm] <- -Inf
  nSim <- nrow(mZ1)
  nK1 <- ncol(mZ1)
  
  # Adjust sample size
  mN1 <- matrix( rep( vN1, nSim ), nSim, length(vN1), byrow = T)
  mN2 <- lData$mN2
  mN2[!cbind(TRUE,mSelArm)] <- 0
  mN <- mN1 + mN2
  
  if(strW == "ctrl"){
    vW <- mN1[,1] / mN[,1]
  }else if(strW == "total"){
    vW <- rowSums(mN1) / rowSums(mN)
  }else if(strW == "obs"){
    vW <- rowSums(mN1[,cbind(T,mSelArm),drop = F])/rowSums(mN[,cbind(T,mSelArm),drop = F])
  }else{
    stop("Please provide a valid strW, e.g., 'ctrl', 'total' or 'obs'.")
  }
  
  mPval2 <- 1 - pnorm(mZ2)
  
  dfHypo <- expand.grid( rep( list(c(TRUE,FALSE)) , nK1 ) )
  dfHypo <- dfHypo[rowSums(dfHypo) != 0, ] # remove H_K where K is empty
  colnames(dfHypo) <- paste0("H", 1:nK1)
  
  dfHypo <- as.data.frame(dfHypo)
  
  dfTestK <- NULL
  # for each intersection hypothesis
  for(i in 1:nrow(dfHypo)){
    vTrt <- unlist(dfHypo[i, 1:nK1])
    
    lTestTmp <- lapply(1:nSim, function(nIter) {
      vCondErrTmp <- InvNCondErr( dAlpha = dAlpha/sum(vTrt),
                                  vW = vW[nIter],
                                  vZ1 = mZ1[nIter,vTrt])
      dCondErr <- sum(vCondErrTmp)
      
      vCondErr <- numeric(nK1)
      vCondErr[vTrt] <- vCondErrTmp
      
      vTrtSub <- vTrt&mSelArm[nIter,]
      
      vCondErrUp <- numeric(nK1)
      if(sum(vTrtSub) != 0){
        if( sum(vTrtSub) != sum(vTrt) ){
          f1 <- function(alpha) {
            vCondErrTmp <- InvNCondErr( dAlpha = alpha/sum(vTrtSub),
                                        vW = vW[nIter],
                                        vZ1 = mZ1[nIter,vTrtSub])
            (sum(vCondErrTmp) - dCondErr)^2
          }
          
          dAlphaAdapt <- stats::optimize(f1, c(0, 1))$minimum
          vCondErrTmp <- InvNCondErr( dAlpha = dAlphaAdapt/sum(vTrtSub),
                                      vW = vW[nIter],
                                      vZ1 = mZ1[nIter,vTrtSub])
        }
        vCondErrUp[vTrtSub] <- vCondErrTmp
      }
      
      if(dCondErr >= 1){
        bRejKlin <- TRUE
      }else{
        bRejKlin <- any(mPval2[nIter,vTrtSub] <= vCondErrUp[vTrtSub])
      }
      
      bRejPosch <- any(mPval2[nIter,vTrtSub] <= vCondErrUp[vTrtSub]/dCondErr*min(dCondErr,1))
      bRejMid <- any(mPval2[nIter,vTrtSub] <= vCondErrUp[vTrtSub])
      
      data.frame(nIter,dCondErr,bRejKlin,bRejPosch,bRejMid)
      
    } )
    dfTestKTmp <- do.call(rbind, lTestTmp)
    mSelSub <- mSelArm
    mSelSub[,!vTrt] <- FALSE
    colnames(mSelSub) <- paste0("H",1:nK1)
    
    dfTestKTmp <- cbind( dfTestKTmp, mSelSub )
    dfTestK <- rbind(dfTestK, dfTestKTmp)
  }
  
  dfTest <- NULL
  for(j in seq_len(nK1)){
    
    dfContain <- dfTestK[dfTestK[,paste0("H",j)],,drop=FALSE]
    if(nrow(dfContain) == 0) next
    dfTestTmp <- stats::aggregate(cbind(bRejKlin, bRejPosch, bRejMid) ~ nIter, 
                                  data = dfContain, FUN = all)
    
    vDrop <- base::setdiff( 1:nSim, dfTestTmp$nIter)
    if(length(vDrop) > 0){
      dfTestTmpU <- data.frame( nIter = vDrop,
                                bRejKlin  = FALSE,
                                bRejPosch = FALSE,
                                bRejMid   = FALSE
      )
      dfTestTmp <- rbind(dfTestTmp, dfTestTmpU)
    }
    dfTestTmp <- dfTestTmp[order(dfTestTmp$nIter),]
    dfTestTmp$nTrt <- j
    dfTest <- rbind(dfTest, dfTestTmp)
  }
  
  dfRej <- stats::aggregate(cbind(bRejKlin, bRejPosch, bRejMid) ~ nIter, 
                            data = dfTest, FUN = any)
  
  return( list(dfRej = dfRej,
               dfTestK = dfTestK,
               dfTest = dfTest) )
}
