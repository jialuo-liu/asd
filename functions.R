# ---- Part I: Methods ----
# ---- IsVector: check if an object is a vector ----
#' @param vZ an object

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
#' @param vN a vector, if not provided or length = 1, all arms have equal sample size; if length > 1, sample sizes (or relative proportion) for each arm, where the first element represents the control arm.

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
      dPval <- 1 - integrate(Pfun, -Inf, Inf)$value
      
    }else{
      if(length(vN) - 1 != length(vZ) )
        stop("The length of vN should equal the total number of arms, including control.")
      vRR <- vN[-1]/(vN[-1] + vN[1])
      
      if( length( unique(vRR) ) == 1 ){
        Pfun <- function(x) ( pnorm( ( dZMax + sqrt(vRR[1]) * x )/sqrt(1-vRR[1]) ) )^nArm * dnorm(x)
        dPval <- 1 - integrate(Pfun, -Inf, Inf)$value
      }else{
        mCov <- sqrt ( outer( vRR, vRR ) )
        diag(mCov) <- 1
        dPval <- 1 - mvtnorm::pmvnorm( upper = rep(dZMax, nArm),
                                       sigma = mCov,
                                       algorithm = GenzBretz(maxpts = 1e5, abseps = 1e-5, releps = 0))
      }
    }
  }
  
  return( min(dPval,1) )
}

# ---- SimesPvalue: Get Simes P-value  ----
#' @param vZ a vector of Z statistics

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
#' @param vN a vector, if not provided or length = 1, all arms have equal sample size; if length > 1, sample sizes (or relative proportion) for each arm, where the first element represents the control arm. Only required for Dunnett test
#' @param vMethod a vector of multiple testing procedures, "Bonf" = Bonferroni, "Simes" = Simes, "Dunn" = Dunnett method.

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
#' @param vZ a vector, Z values for each treatment arm.
#' @param vSelArm a logical vector
#' @param vN a vector, if not provided or length = 1, all arms have equal sample size; if length > 1, sample sizes (or relative proportion) for each arm, where the first element represents the control arm. Only required for Dunnett test
#' @param vMethod a vector of multiple testing procedures, "Bonf" = Bonferroni, "Simes" = Simes, "Dunn" = Dunnett method.

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
#' @param lData a list, data by stage and pooled
#' @param mSelArm a logical matrix, each row representing whether an arm is selected or not, e.g., c(T, F) if dose 1 is selected and dose 2 is not
#' @param vMethod a vector of multiple testing procedures, "Bonf" = Bonferroni, "Simes" = Simes, "Dunn" = Dunnett method.

S2Only <- function( lData, mSelArm,
                    vMethod = c("Bonf","Simes","Dunn")){
  
  if( ! identical( dim(lData$mZ2), dim(mSelArm)) )
    stop("The dimensions of mZ2 and mSelArm do not match!")
  
  mZ2 <- lData$mZ2
  nK2 <- ncol(mZ2)
  nSim <- nrow(mZ2)
  mN2 <- lData$mN2
  lTest <- lapply(1:nSim, function(i) {
    dfTestK <- CTPTest(vZ = mZ2[ i , mSelArm[i,] ],
                       # only use Z stats/sample size of selected arms
                       vN = mN2[ i , c(T,mSelArm[i,]) ],
                       vMethod = vMethod)$dfTest
    TrtMap <- setNames( (1:nK2)[mSelArm[i,]], 1:sum(mSelArm[i,]) )
    dfTestK$nTrt <- TrtMap[dfTestK$nTrt]
    dfTestKU <- as.data.frame( matrix( 1, sum(!mSelArm[i,] ), length(vMethod) + 1) )
    colnames( dfTestKU ) <- colnames( dfTestK )
    dfTestKU$nTrt <- (1:nK2)[!mSelArm[i,]]
    dfTestK <- rbind(dfTestK, dfTestKU)
    dfTestK <- dfTestK[order(dfTestK$nTrt),]
    dfTestK$nIter <- i
    dfTestK
  } )
  dfTest <- as.data.frame( do.call(rbind,lTest) )
  dfTmp <- dfTest[,c("nIter",paste0("dPval",vMethod))]
  dfPval <- aggregate(. ~ nIter, data = dfTmp, FUN = min)
  return( list( dfPval = dfPval,
                dfTest = dfTest) )
}

# ---- Method 2: Combination Function ----
# ---- PvalComb: Combine P value across stages ----
# This function also apply to matrices, treating each matrix as a vector and return a matrix with same dimension
#' @param vPval1 a vector, first stage p-values.
#' @param vPval2 a vector, second stage p-values.
#' @param vW     a vector, information fractions at first stage; 2nd stage is 1- vW.

PvalComb <- function( vPval1, vPval2, vW){
  if( length(vW) != 1 & var( c( length(vPval1), length(vPval2), length(vW)) ) != 0 )
    stop("The dimension of vW, vPval1, and vPval2 do not match.")
  vPval <- 1 - pnorm( sqrt(vW) * qnorm( 1 - vPval1 ) + sqrt(1-vW) * qnorm( 1 - vPval2 ) )
  
  return( vPval )
}

# ---- CombFunInvN: Combination Function Approach by Weighted InvNorm ----
#' @param lData a list, data by stage and pooled
#' @param mSelArm a logical matrix, each row representing whether an arm is selected or not, e.g., c(T, F) if dose 1 is selected and dose 2 is not
#' @param vMethod a vector of multiple testing procedures, "Bonf" = Bonferroni, "Simes" = Simes, "Dunn" = Dunnett method.
#' @param strW a character, "ctrl" = by control sample size per stage, "total" = by total sample size per stage

CombFunInvN <- function( lData, mSelArm,
                         vMethod = c("Bonf","Simes","Dunn"),
                         strW){
  
  if( ! identical( dim(lData$mZ), dim(mSelArm)) )
    stop("The dimensions of mZ and mSelArm do not match!")
  
  mZ1 <- lData$mZ1
  mZ2 <- lData$mZ2
  vN1 <- lData$vN1
  #vN2 <- lData$vN2
  
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
                       vSelArm = mSelArm[i,], # either add vSelArm or subsetting mZ2 and vN.
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
  
  dfTest <- NULL
  for(j in 1:nK1){
    
    dfContain <- dfTestK[dfTestK[,paste0("H",j)],,drop=FALSE]
    if(nrow(dfContain) == 0) next
    dfContain <- dfContain[,c("nIter",paste0("dPval",vMethod))]
    dfTestTmp <- aggregate(. ~ nIter, data = dfContain, FUN = max)
    dfTestTmp$nTrt <- j
    dfTest <- rbind(dfTest, dfTestTmp)
    
  }
  
  dfTmp <- dfTest[,c("nIter",paste0("dPval",vMethod))]
  dfPval <- aggregate(. ~ nIter, data = dfTmp, FUN = min)
  
  return( list(dfPval = dfPval,
               dfTestK = dfTestK,
               dfTest = dfTest) )
}

# ---- Method 4: Marginal Combination Function ----
# ---- MargPInvN: Marginal P-value combination function approach ----
#' @param lData a list, data by stage and pooled
#' @param mSelArm a logical matrix, each row representing whether an arm is selected or not, e.g., c(T, F) if dose 1 is selected and dose 2 is not
#' @param vMethod a vector of multiple testing procedures, "Bonf" = Bonferroni, "Simes" = Simes, "Dunn" = Dunnett method.

MargPInvN <- function( lData, mSelArm,
                       vMethod = c("Bonf","Simes","Dunn"),
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
    stop("Please provide a valid strW, e.g., 'ctrl' or 'total'.")
  }
  
  mPval <- PvalComb(
    vPval1 = mPval1,
    vPval2 = mPval2,
    vW     = mW
  )
  
  lTest <- lapply( 1:nSim, function(i) {
    dfTestK <- CTPTest(vZ = qnorm(1-mPval)[i, ], vN = mN[i, ], vMethod = vMethod )$dfTest
    dfTestK$nIter <- i
    dfTestK
  })
  
  dfTest <- as.data.frame( do.call(rbind, lTest) )
  
  dfTmp <- dfTest[,c("nIter",paste0("dPval",vMethod))]
  dfPval <- aggregate(. ~ nIter, data = dfTmp, FUN = min)
  
  return( list( dfPval = dfPval,
                dfTest = dfTest ) )
}

# ---- Method 5: Conditional Error Function ----
# ---- getDs: ds is the critical value of the classical Dunnett test ----
#' @param dAlpha alpha
#' @param vN sample sizes per arm, where the first element represents the control arm.

getDs <- function( dAlpha, vN){
  
  vRR <- vN[-1]/(vN[-1] + vN[1])
  
  if( length( unique(vRR) ) == 1 ){
    
    Afun <- function(z) {
      Pfun <- function(x) ( pnorm( ( z + sqrt(vRR[1]) * x )/sqrt(1-vRR[1]) ) )^length(vRR) * dnorm(x)
      1 - integrate(Pfun, -Inf, Inf)$value - dAlpha
    }
    dDs <- uniroot(Afun, c(-1e5,1e5))$root
    
  }else{
    mCorr <- sqrt ( outer( vRR, vRR ) )
    diag(mCorr) <- 1
    dDs <- mvtnorm::qmvnorm( 1-dAlpha, sigma = mCorr, tail = "lower",
                             algorithm = GenzBretz(maxpts = 1e5, abseps = 1e-5, releps = 0))$quantile
    
  }
  
  attributes(dDs) <- NULL
  
  return( dDs )
}

# ---- DunnCondError: Calculate Conditional Error ----
#' @param vN1 a vector of total sample sizes per arm (control, dose 1, dose 2, ...)
#' @param vN2 a vector of stage 2 sample sizes per arm (control, dose 1, dose 2, ...)
#' @param dDs critical value
#' @param vZ1 a vector of stage 1 Z test statistics (dose 1, dose 2, ...)

DunnCondError <- function( vN1, vN2, dDs, vZ1 ){
  
  if( var( c(length(vN1)-1, length(vN2)-1, length(vZ1)) ) !=0 )
    stop("The length of vN1, vN2, vZ1, vZ2 should match.")
  
  vN <- vN1 + vN2
  nArm <- length(vZ1)
  vRR <- vN[-1]/(vN[-1] + vN[1])
  
  mR11 <- sqrt( outer( vN1[-1]/(vN1[1]+vN1[-1]), vN1[-1]/(vN1[1]+vN1[-1]) ) )
  diag(mR11) <- 1
  
  mR22 <- sqrt( outer( vN[-1]/(vN[1]+vN[-1]), vN[-1]/(vN[1]+vN[-1]) ) )
  diag(mR22) <- 1
  
  mR12 <- sqrt( vN1[-1]/vN[-1] * outer( vN1[-1]/(vN1[1]+vN1[-1]), vN[-1]/(vN[1]+vN[-1]) ) )
  diag(mR12) <- sqrt( (vN1[1]+vN1[-1])/(vN[1]+vN[-1]) )
  
  vMean <- c( t(mR12) %*% solve( mR11, vZ1) )
  mCov  <- mR22 - t(mR12) %*% solve( mR11, mR12 )
  
  dCondErr <- 1 - mvtnorm::pmvnorm( upper = rep(dDs, nArm),
                                    mean =  vMean,
                                    sigma = mCov,
                                    algorithm = GenzBretz(maxpts = 1e5, abseps = 1e-5, releps = 0))
  attributes(dCondErr) <- NULL
  return( dCondErr )
}

# ---- CondDunnPvalue: Conditional second-stage Dunnett P-value  ----
#' @param vZ1 a vector of stage 1 Z test statistics (dose 1, dose 2, ...)
#' @param vZ2 a vector of stage 2 Z test statistics (dose 1, dose 2, ...)
#' @param vN1 a vector of total sample sizes per arm (control, dose 1, dose 2, ...)
#' @param vN2 a vector of stage 2 sample sizes per arm (control, dose 1, dose 2, ...)

CondDunnPvalue <- function( vZ1, vZ2, vN1, vN2 ){
  
  if( var( c(length(vN1)-1, length(vN2)-1, length(vZ1), length(vZ2)) ) !=0 )
    stop("The length of vN1, vN2, vZ1, vZ2 should match.")

  if(length(vZ1) == 0){
    dPvalue <- 1
  }else{
    vN <- vN1 + vN2

    dInfoFrac <- sum(vN1)/sum(vN) # should be equivalent to sum(vN1)/sum(vN)
    vZ <- sqrt(dInfoFrac)*vZ1 + sqrt(1-dInfoFrac)*vZ2
    
    dZMax <- max(vZ)
    
    nArm <- length(vZ)
    mR11 <- sqrt( outer( vN1[-1]/(vN1[1]+vN1[-1]), vN1[-1]/(vN1[1]+vN1[-1]) ) )
    diag(mR11) <- 1
    
    mR22 <- sqrt( outer( vN[-1]/(vN[1]+vN[-1]), vN[-1]/(vN[1]+vN[-1]) ) )
    diag(mR22) <- 1
    
    mR12 <- sqrt( vN1[-1]/vN[-1] * outer( vN1[-1]/(vN1[1]+vN1[-1]), vN[-1]/(vN[1]+vN[-1]) ) )
    diag(mR12) <- sqrt( (vN1[1]+vN1[-1])/(vN[1]+vN[-1]) )
    
    vMean <- c( t(mR12) %*% solve( mR11, vZ1) )
    mCov  <- mR22 - t(mR12) %*% solve( mR11, mR12 )
    dPvalue <- 1 - mvtnorm::pmvnorm( upper = rep(dZMax, nArm),
                                     mean =  vMean,
                                     sigma = mCov,
                                     algorithm = GenzBretz(maxpts = 1e5, abseps = 1e-5, releps = 0))
    attributes(dPvalue) <- NULL
  }
  
  return( dPvalue )
}

# ---- AdaptDunnTest: Adaptive Dunnett Tests, Closure Principle ----
#' @param lData a list, data by stage and pooled
#' @param mSelArm a logical matrix, each row representing whether an arm is selected or not, e.g., c(T, F) if dose 1 is selected and dose 2 is not
#' @param dAlpha alpha
#' @param dEps a small value for finding root

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
  
  vHypo <- paste0("H",apply(dfHypo,1,function(x) paste((1:3)[x],collapse = "")))
  dfTestK <- NULL
  for(i in 1:nrow(dfHypo)){
    
    vTrt <- unlist(dfHypo[i, 1:nK1])
    
    # get Ds outside the loop
    vN2Total <- rowSums(mN2) 
    vN2U <- unique(vN2Total)
    dfDs <- NULL
    for( j in 1:length(vN2U)){
      vN2Ori <- vN2U[j] * vN1/sum(vN1)
      dDs <- getDs( dAlpha, vN1[c(T,vTrt)] + vN2Ori[c(T,vTrt)])
      dfTmp <- data.frame(dDs = dDs, vN2Total = vN2U[j])
      dfDs <- rbind(dfDs, dfTmp)
    }
    
    dfDs <- merge(data.frame(vN2Total = vN2Total), dfDs, by = "vN2Total",all.x = T)
    
    lTestKTmp <- lapply(1:nSim, function(nIter) {
      vTrtSub <- vTrt & mSelArm[nIter,]
      if(sum(vTrtSub) == 0) {
        dCondPval <- 1
        dCondErr <- 0
      }else{
        # using selected arm only to get 2nd stage p-value
        dCondPval <- CondDunnPvalue(vZ1 = mZ1[nIter,vTrtSub],
                                     vZ2 = mZ2[nIter,vTrtSub],
                                     vN1 = vN1[c(T,vTrtSub)],
                                     vN2 = mN2[nIter, c(T,vTrtSub)])
        
        # conditional error function based on all first stage data, including those from unselected arms
        dCondErr <- DunnCondError( vN1 = vN1[c(T,vTrt)],
                                   vN2 = vN2Ori[c(T,vTrt)],
                                   dDs = dfDs[nIter, "dDs"],
                                   vZ1 = mZ1[nIter,vTrt])
        
      }
      bCondRej <- ifelse(dCondPval <= dCondErr, TRUE, FALSE)
      data.frame(nIter     = nIter,
                 strHypo   = vHypo[i],
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
  for(j in 1:nK1){
    dfContain <- dfTestK[dfTestK[,paste0("H",j)],,drop=FALSE]
    if(nrow(dfContain) == 0) next
    dfTestTmp <- aggregate(vCondRej ~ nIter, data = dfContain, FUN = all)
    vDrop <- setdiff( 1:nSim, dfTestTmp$nIter)
    if(length(vDrop) > 0){
      dfTestTmpU <- data.frame( nIter = vDrop, vCondRej = FALSE)
      dfTestTmp <- rbind(dfTestTmp, dfTestTmpU)
    }
    dfTestTmp <- dfTestTmp[order(dfTestTmp$nIter),]
    dfTestTmp$nTrt <- j
    dfTest <- rbind(dfTest, dfTestTmp)
  }
  
  dfRej <- aggregate(vCondRej ~ nIter, data = dfTest, FUN = any)
  
  return( list(dfRej = dfRej,
               dfTestK = dfTestK,
               dfTest = dfTest) )
}

# ---- Method 6: Partial Conditional Error Function ----
# ---- InvNCondErr: Conditional Error based on Weighted Inverse normal Combination Function ----
# This function also apply to matrices, treating each matrix as a vector and return a matrix with same dimension.
#' @param dAlpha alpha
#' @param vW     a vector, information fractions.
#' @param vZ1 a vector of stage 1 Z test statistics

InvNCondErr <- function(dAlpha, vW, vZ1){
  if( length(vW) != 1 & var( c( length(vZ1), length(vW)) ) != 0 )
    stop("The dimension of vW, vPval1, and vPval2 do not match.")
  
  1-pnorm((qnorm(1-dAlpha) - sqrt(vW)*vZ1)/sqrt(1-vW))
}


InvNCondErr1 <- function(dAlpha, vW, vZ1){
  if( length(vW) != 1 & var( c( length(vZ1), length(vW)) ) != 0 )
    stop("The dimension of vW, vPval1, and vPval2 do not match.")
  
  vMean <- vZ1*sqrt(vW)
  mCov <- diag( length(vZ1) )
  diag(mCov) <- 1 - vW
  dCondErr <- 1 - mvtnorm::pmvnorm(upper = rep(qnorm(1-dAlpha),length(vZ1)),
                                   mean =  vMean,
                                   sigma = mCov,
                                   algorithm = GenzBretz(maxpts = 1e5, abseps = 1e-5, releps = 0))
  attributes(dCondErr) <- NULL
  return(dCondErr)
}
# ---- BonfPartCondError: Partial Conditional Error for Bonferroni ----
#' @param lData a list, data by stage and pooled
#' @param mSelArm a logical matrix, each row representing whether an arm is selected or not, e.g., c(T, F) if dose 1 is selected and dose 2 is not
#' @param dAlpha alpha
#' @param strW a character, "ctrl" = by control sample size per stage, "total" = by total sample size per stage

BonfPartCondError <- function( lData, mSelArm, dAlpha, strW ){
  
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
    stop("Please provide a valid strW, e.g., 'ctrl' or 'total'.")
  }
  
  mPval2 <- 1 - pnorm(mZ2)
  
  dfHypo <- expand.grid( rep( list(c(TRUE,FALSE)) , nK1 ) )
  dfHypo <- dfHypo[rowSums(dfHypo) != 0, ] # remove H_K where K is empty
  colnames(dfHypo) <- paste0("H", 1:nK1)
  
  dfHypo <- as.data.frame(dfHypo)
  vHypo <- paste0("H",apply(dfHypo,1,function(x) paste((1:3)[x],collapse = "")))
  
  dfTestK <- NULL
  # for each intersection hypothesis
  for(i in 1:nrow(dfHypo)){
    vTrt <- unlist(dfHypo[i, 1:nK1])
    
    lTestTmp <- lapply(1:nSim, function(nIter) {
      dCondErr <- sum(InvNCondErr( dAlpha = dAlpha/sum(vTrt),
                                   vW = vW[nIter],
                                   vZ1 = mZ1[nIter,vTrt]))
      vTrtSub <- vTrt&mSelArm[nIter,]
      if(dCondErr >= 1){
        bRejKlin <- TRUE
      }else{
        bRejKlin <- any(mPval2[nIter,vTrtSub] <= dCondErr/sum(vTrtSub))
      }
      bRejPosch <- any(mPval2[nIter,vTrtSub] <= min(dCondErr,1)/sum(vTrtSub))
      bRejMid <- any(mPval2[nIter,vTrtSub] <= dCondErr/sum(vTrtSub))
      strHypo <- vHypo[i]
      data.frame(nIter,strHypo,dCondErr,bRejKlin,bRejPosch,bRejMid)
    } )
    dfTestKTmp <- do.call(rbind, lTestTmp)
    mSelSub <- mSelArm
    mSelSub[,!vTrt] <- FALSE
    colnames(mSelSub) <- paste0("H",1:nK1)
    
    dfTestKTmp <- cbind( dfTestKTmp, mSelSub )
    dfTestK <- rbind(dfTestK, dfTestKTmp)
  }
  
  dfTest <- NULL
  for(j in 1:nK1){
    
    dfContain <- dfTestK[dfTestK[,paste0("H",j)],,drop=FALSE]
    if(nrow(dfContain) == 0) next
    dfTestTmp <- aggregate(cbind(bRejKlin, bRejPosch, bRejMid) ~ nIter, data = dfContain, FUN = all)
    
    vDrop <- setdiff( 1:nSim, dfTestTmp$nIter)
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
  
  dfRej <- aggregate(cbind(bRejKlin, bRejPosch, bRejMid) ~ nIter, data = dfTest, FUN = any)
  
  return( list(dfRej = dfRej,
               dfTestK = dfTestK,
               dfTest = dfTest) )
}

# ---- Part II: Simulation ----
# ---- GenerateDataAll: Simulate Xbar and Z statistics ----
#' @param nSim number of simulations
#' @param vN1 vector of stage 1 sample sizes (control, dose 1, dose 2, ...), can be a single number if all arms have the same sample size
#' @param vN2 vector of stage 2 sample sizes (control, dose 1, dose 2, ...), can be a single number if all arms have the same sample size
#' @param vMu vector of means        (control, dose 1, dose 2, ...)
#' @param dSD standard deviation     (assume same and known standard deviation for all arms)
GenerateDataAll <- function(nSim, vN1, vN2, vMu, dSD){
  
  if(length(vN1) == 1) vN1 <- rep(vN1, length(vMu))
  if(length(vN2) == 1) vN2 <- rep(vN2, length(vMu))
  
  if(length(vN1) != length(vMu)) stop("vN1 and vMu should have the same length.")
  if(length(vN2) != length(vMu)) stop("vN2 and vMu should have the same length.")
  
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

# ---- RunSimFixTotal: Run simulations with Total sample size (across stages) fixed  ----
#' @param i the ith simulation setting
#' @param arg a data frame of specifications of all simulation scenarios 

RunSimFixTotal  <- function(i, arg)
{
  lPara  <- unlist(arg[i, ], recursive = FALSE, use.names = TRUE)
  nSim <- lPara$nSim
  
  tic <- proc.time()
  
  vN1 <- lPara$vRR * lPara$vSS[1]
  vN2Ori <- lPara$vRR * lPara$vSS[2]
  nArm <- length(vN1) - 1
  
  vDenom <- unique(switch(lPara$strRule,
                          SelectAll = apply(combn(nArm,nArm),2,function(x) sum(lPara$vRR[c(1,x+1)])),
                          SelectTwo = apply(combn(nArm,2),2,function(x) sum(lPara$vRR[c(1,x+1)])),
                          SelectOne = apply(combn(nArm,1),2,function(x) sum(lPara$vRR[c(1,x+1)]))))
  
  lDataList <- list()
  for(k in 1:length(vDenom)){
    set.seed( lPara$nSeed )
    lDataList[[k]] <- GenerateDataAll(nSim = lPara$nSim,
                                      vN1 = vN1,
                                      vN2 = sum(vN2Ori)*lPara$vRR/vDenom[k],
                                      vMu = lPara$vMu,
                                      dSD = lPara$dSD)
    
  }
  
  nK1 <- ncol(lDataList[[1]]$mZ1)
  if(lPara$strRule == "SelectTwo"){
    mSelArm <- t(apply(lDataList[[1]]$mZ1,1,function(x) {
      tmp <- rep(TRUE, length(x))
      tmp[which.min(x)] <- FALSE
      tmp}))
  }
  if(lPara$strRule == "SelectOne"){
    mSelArm <- t(apply(lDataList[[1]]$mZ1,1,function(x) {
      tmp <- rep(FALSE, length(x))
      tmp[which.max(x)] <- TRUE
      tmp}))
  }
  if(lPara$strRule == "SelectAll"){
    mSelArm <- t(apply(lDataList[[1]]$mZ1,1,function(x) {
      tmp <- rep(TRUE, length(x))
      tmp}))
  }
  
  mZ <- mZ2 <- matrix(0, nrow = nSim, ncol = nArm)
  mN2 <- matrix(0, nrow = nSim, ncol = nArm+1)
  for(k in 1:length(vDenom)){
    mTmp <- lDataList[[k]]$mZ2
    mTmp[!(cbind(T,mSelArm) %*% lPara$vRR == vDenom[k]),] <- 0
    mZ2 <- mZ2 + mTmp
    
    mTmp <- lDataList[[k]]$mZ
    mTmp[!(cbind(T,mSelArm) %*% lPara$vRR == vDenom[k]),] <- 0
    mZ <- mZ + mTmp
    
    mTmp <- lDataList[[k]]$mN2
    mTmp[!(cbind(T,mSelArm) %*% lPara$vRR == vDenom[k]),] <- 0
    mN2 <- mN2 + mTmp
  }
  
  lData <- list(mZ1 = lDataList[[1]]$mZ1,
                mZ2 = mZ2, mZ = mZ, vN1=lDataList[[1]]$vN1, mN2=mN2)
  
  lS2 <- S2Only( lData, mSelArm,
                 vMethod = c(c("Bonf","Simes","Dunn")))
  
  lAdaptDunnTest <- AdaptDunnTest(lData, mSelArm,
                                  dAlpha = 0.025
  )
  
  lMargP <- MargPInvN( lData, mSelArm,
                       vMethod = c("Bonf","Simes","Dunn"),
                       strW = "total")
  
  lCombFnWTotal <- CombFunInvN( lData, mSelArm,
                                vMethod = c(c("Bonf","Simes","Dunn")),
                                strW = "total")
  
  lBonfPCEWTotal <- BonfPartCondError( lData, mSelArm, dAlpha=0.025,strW = "total" )
  
  lRej <- list()
  lTest <- list()
  
  for(k in (1:nK1)[apply(mSelArm, 2, any)]){
    lTest[[k]] <- data.frame(
      nTrt          = k,
      nIter         = 1:nSim,
      marP_Dunn     = with( subset(lMargP$dfTest,       nTrt==k),  dPvalDunn[ order(nIter)] ),
      combFnT_Dunn  = with( subset(lCombFnWTotal$dfTest,nTrt==k),  dPvalDunn[ order(nIter)] ),
      S2_Dunn       = with( subset(lS2$dfTest          ,nTrt==k),  dPvalDunn[ order(nIter)] ),
      marP_Simes    = with( subset(lMargP$dfTest       ,nTrt==k), dPvalSimes[ order(nIter)] ),
      combFnT_Simes = with( subset(lCombFnWTotal$dfTest,nTrt==k), dPvalSimes[ order(nIter)] ),
      S2_Simes      = with( subset(lS2$dfTest          ,nTrt==k), dPvalSimes[ order(nIter)] ),
      marP_Bonf     = with( subset(lMargP$dfTest       ,nTrt==k),  dPvalBonf[ order(nIter)] ),
      combFnT_Bonf  = with( subset(lCombFnWTotal$dfTest,nTrt==k),  dPvalBonf[ order(nIter)] ),
      S2_Bonf       = with( subset(lS2$dfTest          ,nTrt==k),  dPvalBonf[ order(nIter)] ))
    lRej[[k]] <- data.frame(
      nTrt          = k,
      nIter         = 1:nSim,
      CEF_Dunn      = with( subset(lAdaptDunnTest$dfTest,nTrt==k), vCondRej [ order(nIter) ] ) ,
      KlinT_Bonf    = with( subset(lBonfPCEWTotal$dfTest,nTrt==k), bRejKlin [ order(nIter) ] ) ,
      PoschT_Bonf   = with( subset(lBonfPCEWTotal$dfTest,nTrt==k), bRejPosch[ order(nIter) ] ) ,
      MidT_Bonf     = with( subset(lBonfPCEWTotal$dfTest,nTrt==k), bRejMid  [ order(nIter) ] ) )
  }
  
  dfRej <- do.call(rbind, lRej)
  dfTest <- do.call(rbind, lTest)
  
  toc <- proc.time() - tic
  
  lRes <- list(lPara = lPara, dfRej = dfRej, dfTest = dfTest,
               pid = Sys.getpid(), time = toc[3], mSelArm = mSelArm)
  
  saveRDS(lRes, file=paste0("FixTotal_Res",i,".rds"))
  return(lRes)
  
}

