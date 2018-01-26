
#############################################################################
###
###  Supporting Functions
###
#############################################################################

###turn a vector into a matrix
mat.vec    <- function(the.vec,nr)
{
    matrix(the.vec,byrow=TRUE,ncol=length(the.vec),nrow=nr)
}

###center a vector
col.center <- function(the.mat) {
    n <- nrow(the.mat); c <- ncol(the.mat);
    cMeans <- .colMeans(the.mat, m = n, n = c, na.rm=TRUE);
    if(sum(cMeans^2) < 10^-20) the.mat else
    the.mat - mat.vec(cMeans, nr=n)}

###calculate the sd for columns in a matrix
col.sd     <- function(the.mat)
{
    n <- nrow(the.mat); c <- ncol(the.mat);
    cMeans <- .colMeans(the.mat, m = n, n = c, na.rm=TRUE);
    if(sum(cMeans^2) < 10^-20) {
        sqrt((n/(n-1))*.colMeans(the.mat^2, m = n, n = c,na.rm=TRUE))
    } else {
        sqrt((n/(n-1))*.colMeans((col.center(the.mat))^2, m = n, n = c,
                                 na.rm=TRUE))
    }
}

###calculate the sd for columns in a matrix
col.norm   <- function(the.mat) col.center(the.mat)/mat.vec(col.sd(the.mat),
                                                            nr=nrow(the.mat))

##Bonferroni approach for medTest.SBMH
proc.intersection.adaptivebonf = function(pv1,pv2, t1=0.025, t2=0.025, lambda=0){
  if (lambda>0){
    t1 = min(t1,lambda)
    t2= min(t2,lambda)
  }
  
  selected =which((pv1<=t1) & (pv2<=t2))
  R1 = sum(pv1<=t1)
  R2 = sum(pv2<=t2)
  
  if(length(selected)==0){
    return(list(r.value = rep(NA,length(pv1)), R1=R1, R2=R2,  selected = selected))
  }
  
  S1 = (pv1<=t1) #the selected from study 1
  pi2 = (1+sum(pv2[S1]>lambda))/(R1*(1-lambda)) #the fraction of nulls in study 2 among S1
  
  S2 = (pv2<=t2) #the selected from study 2
  pi1 = (1+sum(pv1[S2]>lambda))/(R2*(1-lambda)) #the fraction of nulls in study 1 among S2
  
  if (lambda==0){ #ie nonadaptive
    pi2=1
    pi1=1
  }
  
  if((pi1 > 1) | (pi2 > 1))
  {
    warning("At least one of the estimated fraction of nulls is > 1: Using lambda=0 to set them to 1.")
    pi2=1
    pi1=1
  }
  
  r.value = rep(NA,length(pv1))
  r.value[selected] = pmax(pi1*R2*pv1[selected]/0.5, pi2*R1*pv2[selected]/0.5)
  
  return(list(r.value = r.value, R1=R1, R2=R2,selected = selected,pi1=pi1, pi2=pi2))
}

##FDR approach for medTest.SBMH
proc.intersection.adaptiveFDR = function(pv1,pv2, t1=0.025, t2=0.025,lambda=0.0){
  #if adaptive=TRUE, estimates the fraction of zeros of study i in selection j, for (i,j) = (1,2) or (2,1)
  #incorporate it into procedure by multiplying the p-value of study i by this fraction.   
  #print(c(t1,t2))
  #print(summary(pv1))
  #print(summary(pv2))
  #print(lambda)
  if (lambda>0){
    t1 = min(t1,lambda)
    t2= min(t2,lambda)
  }
  selected =which((pv1<=t1) & (pv2<=t2))
  R1 = sum(pv1<=t1)
  R2 = sum(pv2<=t2)
  if(length(selected)==0){
    return(list(r.value = rep(NA,length(pv1)), R1=R1, R2=R2, selected = selected))
  }
  
  S1 = (pv1<=t1) #the selected from study 1
  pi2 = (1+sum(pv2[S1]>lambda))/(R1*(1-lambda)) #the fraction of nulls in study 2 among S1
  
  S2 = (pv2<=t2) #the selected from study 2
  pi1 = (1+sum(pv1[S2]>lambda))/(R2*(1-lambda)) #the fraction of nulls in study 1 among S2
  
  if (lambda==0){ #ie nonadaptive
    pi2=1
    pi1=1
  }
  
  if((pi1 > 1) | (pi2 > 1))
  {
    warning("At least one of the estimated fraction of nulls is > 1: Using lambda=0 to set them to 1.")
    pi2=1
    pi1=1
  }
  
  Z.selected <- pmax(pi1*R2*pv1/0.5, pi2*R1*pv2/0.5)[selected]
  oz <- order(Z.selected, decreasing =TRUE)
  ozr <- order(oz)
  r.value.selected <- cummin( (Z.selected/rank(Z.selected, ties.method= "max"))[oz] )[ozr]
  r.value.selected <- pmin(r.value.selected,1)
  r.value  = rep(NA,length(pv1))
  r.value[selected] = r.value.selected
  return(list(r.value = r.value, R1=R1, R2=R2, selected = selected, pi1=pi1, pi2=pi2))
}

#############################################################################
###
###  medTest
###
#############################################################################

###INPUT:
###Y: 		    Outcome  (n x 1 matrix)
###E: 		    Exposure (n x 1 matrix)
###M: 	     	Mediator (n x p matrix, where p is the number of mediators)
##Z:             Additional covariates (n x c matrix, where c is the number of
##additional covariates)
###nperm:  	number of permutations for estimating p-value
###w:        weight assigned to each subject
###(numeric vector of either length 1 (in which case,
###the same weight is assigned to each subject) or length n
###
###OUTPUT:
###Sp:		Statistic (absolute value of the product of correlations) and
###p-value in a p x 2 matrix

medTest <- function(E,M,Y,Z=NULL,useWeightsZ=TRUE,nperm=100,w=1){

  n <- length(E)
  ##if M is a vector or a data frame, change into a matrix
  if(is.null(dim(M)))
  {
      M <- as.matrix(M, ncol=1)
  }
  M <- as.matrix(M)
  p <- ncol(M)

  ##make sure w is correct
  if(length(w)!=1 & length(w)!=n)
  {
    stop("The length of w must be either 1 or the length of E.")
  }
  if(length(w)==1)
  {
    w <- rep(w, n)
  }

  ##standardize w
  w <- w/sum(w)

  ##if E and Y are vectors, change them into matrices
  E <- matrix(E, nrow=n, ncol=1)
  Y <- matrix(Y, nrow=n, ncol=1)

  ##if Z is not null, take the residuals of E, M, and Y on Z before
  ##performing the remaining analyses
  if(!is.null(Z))
  {
      ##in case Z is a data frame, change into matrix
      Z <- as.matrix(Z)
      ##add intercept
      Z1 <- cbind(1, Z)

      Y <- lm.fit(x=Z1, y=Y)$residuals
      if(useWeightsZ)
      {
          for(m in 1:p)
          {
              M[,m] <- lm.wfit(x=Z1, y=M[,m], w=w)$residuals
          }
          E <- lm.wfit(x=Z1, y=E, w=w)$residuals
      } else {
          for(m in 1:p)
          {
              M[,m] <- lm.fit(x=Z1, y=M[,m])$residuals
          }
          E <- lm.fit(x=Z1, y=E)$residuals
      }

      ##residuals will be returned as vectors, so change back to matrices
      E <- matrix(E, nrow=n, ncol=1)
      Y <- matrix(Y, nrow=n, ncol=1)
  }

  ##normalization
  En		<- E-sum(E)/n
  Mn		<- col.center(M)
  Yn		<- Y-sum(Y)/n
  sdE <- sqrt(sum(En^2)/(n-1))
  sdY <- sqrt(sum(Yn^2)/(n-1))

  #getting the residuals
  tEn <- t.default(En)
  invCrossEn <- 1/sum(En^2)

  B1.obs 	<- invCrossEn*sum(En*Yn)
  B2.obs	<- invCrossEn*tEn%*%Mn
  rY		<- Yn-B1.obs*En
  rM		<- Mn-En%*%B2.obs

  ##calculate these only once
  col.norm.Mn <- col.norm(Mn)
  col.norm.rM <- col.norm(rM)

  ##getting the observed value of the statistic
  ##use weighted E-M correlation
  En.weight <- sqrt(w)*(E-sum(w*E))
  Mn.weight <- sqrt(w)*(M-mat.vec(.colSums(w*M, m=n, n=p, FALSE),
                                  nr=n))
  sdE.weight  <- sqrt(sum(En.weight^2))
  sdMn.weight <- mat.vec(sqrt(.colSums(Mn.weight^2, m=n, n=p, na.rm=TRUE)),
                         nr=n)

  cEM	<- t.default(En.weight/sdE.weight)%*%(Mn.weight/sdMn.weight)
  cYM	<- t.default(col.norm(rY))%*%col.norm.rM/(n-1)

  S     <- abs(cEM*cYM)
  nmed	<- ncol(rM)

  ##the remainder is dedicated to getting p-values

  ##identify the metabolites in Group A, where want to randomize Y
  groupA <- abs(cEM) >= abs(cYM)
  groupB <- !(groupA)

  ##max(S) from nperm permutations
  max.S.mat	<- matrix(nrow=nperm,ncol=1)

  ##only want to calculate this value once
  EEiE     <- invCrossEn*En

  for (the.perm in 1:nperm){
    ##Randomize Y
    if(the.perm == nperm) ##for the last permutation, use observed data
    {
      rYt <- rY
    } else {
      rYt	  	<- rY[sample.int(n), , drop=FALSE]
    }
    B1 		<- c(sum(EEiE*rYt))
    rY2		<- rYt-B1*En
    cYM.A	<- t.default(col.norm(rY2))%*%col.norm.rM/(n-1)

    ##Randomize E
    if(the.perm == nperm) ##for the last permutation, use observed data
    {
      Ent <- En
      Et		<- E
    } else {
      permInd <- sample.int(n)
      Ent   <- En[permInd]
      Et		<- E[permInd, , drop=FALSE]
    }
    tEnt <- t.default(Ent)

    B2	 	<- invCrossEn*tEnt%*%Mn[,groupB]
    rM3		<- rM
    rM3[,groupB]  	<- Mn[,groupB]-Ent%*%B2
    B3	 	<- invCrossEn*sum(Ent*Y)
    rY3		<- rY-B3*Ent

    ##use weighted E-M correlation
    Ent.weight    <- sqrt(w)*(Et-sum(w*Et))
    sdEt.weight  <- sqrt(sum(Ent.weight^2))

    cEM.B	<- t.default(Ent.weight/sdEt.weight)%*%(Mn.weight/sdMn.weight)
    cYM.B <- cYM

    max.S.mat[the.perm] <- max(abs(cEM.B*cYM.B)[groupB],abs(cEM*cYM.A)[groupA])
  }

  pval	<- matrix(nrow=nmed,ncol=1)
  for (i in 1:nmed) pval[i] <- sum(max.S.mat > S[i],na.rm=TRUE)/nperm

  Sp    <- cbind(c(S),c(pval))
  colnames(Sp) <- c("S","p")

  Sp

}

#############################################################################
###
###  medTest.SBMH (Mediator Test based on Bogomolov & Heller)
###
#############################################################################

###INPUT: 
###pEM:     a vector of size m (where m = number of mediators). Entries are the p-values for the E,M_j relationship 
###pMY:     a vector of size m (where m = number of mediators). Entries are the p-values for the M_j,Y|E relationship
###MCP.type:    multiple comparison procedure - either "FWER" or "FDR"
###t1:   threshold for determining the cutoff to be one of the top S_1 E/M_j relationships 
###t2:   threshold for determining the cutoff to be one of the top S_2 M_j/Y relationships 
###adaptive:  FALSE/TRUE depending on whether an adaptive threshold should be used
###
###OUTPUT:
###m x 1 matrix - either p-values (if MCP.type = "FWER") or q-values (if MCP.type = "FDR")

medTest.SBMH <- function(pEM,pMY,MCP.type="FWER",t1=0.05,t2=0.05,lambda=0){
  if (MCP.type=="FWER")  possVal <- proc.intersection.adaptivebonf(pEM,pMY, t1=t1, t2=t2,lambda=lambda)$r.value
  if (MCP.type=="FDR")   possVal <- proc.intersection.adaptiveFDR( pEM,pMY, t1=t1, t2=t2,lambda=lambda)$r.value
  possVal <- ifelse(is.na(possVal),1,possVal)
  ##threshold values at 1
  ifelse(possVal > 1, 1, possVal)
}
