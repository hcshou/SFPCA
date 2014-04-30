########################################################################################
###         This file contains R codes to perform HD-MFPCA with unbalanced design in R
###         combined with the low-dimensional R code
###         Edited 4/30/2014 by Haochang                          
########################################################################################


########################################################################################
###
### Input of the function "HD_MFPCA_full"
###
### y:      a data matrix of dimension n by p containing functional responses. Each row contains measurements 
###         from a function at a set of grid points, and each column contains measurements
###         of all functions at a particular grid point.
### id:     a vector contains the id information that is used to identify clusters. 
### visit:  optional argument if functions from the same cluster have natural orders 
###         such as visit time.
### I:      number of clusters  
### J:      cluster size
### p:      number of grid points per function. 
### tstart: the start point of the grid points
### tlength:the length of the interval that the functions are evaluated at.
### hd:     a logical argument indicates whether if TRUE, use the high-dimensional projection developed by Zipunnikov et al. (2014); FALSE, Di (2009)
### twoway: a logical argument indicates whether oneway or twoway functional ANOVA (analysis
###         of variance) decomposition is more proper for the problem. "twoway=TRUE" will carry
###         out twoway ANOVA, and the default, "twoway=TRUE" will carry out twoway ANOVA.
### smoothing: a logical argument that controls whether or not smoothing is conducted. When 
###         the functional data is measured without noise, smoothing is not needed. On the other
###         hand, if measurement error is present, specifying "smoothing=TRUE" will incorporate 
###         smoothing steps. As a result, the estimated eigenvalues will be less biased and the
###         estimated eigenfunctions will be more smooth. The default is "smoothing=FALSE".
### min.K1: the minimum number of components to keep in level 1 space. The actual number of 
###         components will be the maximum of min.K1 and the number determined by the data.
### min.K2: the minimum number of components to keep in level 2 space.
### 
########################################################################################



########################################################################################
###
### Output of the function "HD_MFPCA_full"
### 
###     The output of the function is a list that contains the following elements. 
###
### N1:         the number of components at level 1
### N2:         the number of components at level 2
### lambda:     estimated reliability ratio
### eigen1:     a vector containing all level 1 eigenvalues in non-increasing order
### eigen2:     a vector containing all level 2 eigenvalues in non-increasing order
### phi1:       a matrix containing all level 1 eigenfunctions. Each row contains an 
###             eigenfunction evaluated at the same set of grid points as the input data.
###             The eigenfunctions are in the same order as the corresponding eigenvalues.
### phi2:       a matrix containing all level 2 eigenfunctions. Each row contains an 
###             eigenfunction evaluated at the same set of grid points as the input data.
###             The eigenfunctions are in the same order as the corresponding eigenvalues.
### scores1:    a matrix containing estimated level 1 principal component scores. Each row
###             corresponding to the level 1 scores for a particular subject in a cluster.
###             The number of rows is the same as that of the input matrix y. Each column 
###             contains the scores for a level 1 component for all subjects.
### scores2:    a matrix containing estimated level 2 principal component scores. Each row
###             corresponding to the level 2 scores for a particular subject in a cluster.
###             The number of rows is the same as that of the input matrix y. Each column 
###             contains the scores for a level 2 component for all subjects.
### mu:         a vector containing the overall mean function.
### eta:        a matrix containing the deviation from overall mean function to visit specific 
###             mean function. The number of rows is the number of visits.
###
############################################################################################
HD_MFPCA_full<- function(y, id=NULL, visit=NULL, J=NULL,I=NULL, p=NULL, tstart=0, tlength=1, twoway=TRUE, hd=TRUE,
                     smoothing=FALSE, L=0.9, min.K1=4, min.K2=4) {
  
  ###     check the missingness of arguments
  ###     1. Either (id,visit) or (I,J) is needed for identifying clusters. 
  ###         If both are missing, give error message "not enough information";
  ###         if only (I,J) is provided, the function would generate (id,visit);
  ###         if only (id,visit) is provided, the function would generate (I,J);
  ###         
  ###     2. If the number of grid points p is missing, the function would generate it.
  
  if( ( is.null(id)|is.null(visit) ) && ( is.null(I)|is.null(J) )  )
     stop("Not enough information! Please provide (id, visit) or (I,J) !")
  if( !( is.null(I)|is.null(J) ) && ( is.null(id) | is.null(visit) ) ) {
     id <- rep(1:I,each=J)
     visit <- rep(1:J, I)
  }
  if( !( is.null(id) | is.null(visit) ) && (is.null(I) | is.null(J) ) )  {
     I <- length(table(id))
     J <- length(table(visit))        
  }
  if(is.null(p)) 
     p <- dim(y)[2]
  
  ###     Generate the set of grid points that are equally spaced.
    t <- seq(tstart, tstart + tlength, length=p)
  
  ## step 1.  Calculate overall mean function 
  ###     If twoway functional ANOVA is needed ("twoway=TRUE"),  the visit specific mean function
  ###     and the deviation from the overall mean to visit specific mean functions are also
  ###     computed. 
    mu   <- apply(y, 2, mean)
    resd <- t( t(y) - mu ) 
  
    if( twoway == TRUE ) {
      T = visit
      eta <- matrix(0, length(unique(T)), p)
      for( j in unique(T) ) {
        if( sum( T == j ) == 0 ) next
        if( sum( T == j ) == 1 ) { 
          eta[ which( unique(T) == j ), ] <- y[ T == j,  ] - mu
        }
        else 
          eta[ which( unique(T) == j ), ] <- apply( y[ T == j,  ], 2, mean ) - mu
      }
      
      ### Calculate residuals by subtracting visit-specific mean from original functions for 
      ### 'twoway == TRUE', or subtracting overall mean function for 'twoway == FALSE'. 
      
      for( j in unique(T) ) { 
        if( sum( T == j) == 0 ) next
        resd[T == j, ] <- t ( t( y[ T == j,  ] ) - (mu + eta[ which( unique(T) == j ), ] ) )
       }				
     }
    
      X <- resd 
 
  
    ## step 2: calculate the covariance operator in high or low dimensional case
      if(hd==TRUE){
        system.time(UUt <- X %*% t(X))
        SVD = svd(UUt)
        U  = SVD$u
        S  = SVD$d
      }else{U <- X}
      
    
      ### Create index for subjects and visits 
      tid  <- table(id)
      n_I0 <- tid[match(names(tid),unique(id))]
      k2   <- sum(n_I0^2);  n    <- dim(y)[1]
      
      U_I  <- rowsum(U,id)
      
      ##step 1 obtain covariance matrix in the reduced dimension of U
      Ku_W = (t(U)%*%diag(n_I0[rep(1:I,n_I0)])%*%U-t(U_I)%*%(U_I))/(k2-n)
      Ku_T = t(U)%*% U/n
      Ku_B = Ku_T - Ku_W
      
      if(hd==TRUE){
        Ku_B<-t(t(sqrt(S)*Ku_B)*sqrt(S))
        Ku_W<-t(t(sqrt(S)*Ku_W)*sqrt(S))
      }   
      
      #get the eigen values and eigenvectors
      eB<-eigen(Ku_B); Sigma_B<- eB$values; A_B<- eB$vectors
      eW<-eigen(Ku_W); Sigma_W<- eW$values; A_W<- eW$vectors      
      # end of Step 2 
      
      # Step 2s
      ###  get the eigenvalues and determine the number of PCs to keep
      fpca1.value <- Sigma_B *sqrt(p/tlength)
      fpca2.value <- Sigma_W *sqrt(p/tlength)
  
      ###     Keep only non-negative eigenvalues
      fpca1.value <- ifelse(fpca1.value>=0, fpca1.value, 0)
      fpca2.value <- ifelse(fpca2.value>=0, fpca2.value, 0)
      ###     Calculate the percentage of variance that are explained by the components
      percent1 <- (fpca1.value)/sum(fpca1.value)
      percent2 <- (fpca2.value)/sum(fpca2.value)
      
      ###     Decide the number of components that are kept at level 1 and 2. The general
      ###     rule is to stop at the component where the cumulative percentage of variance 
      ###     explained is greater than 90% and the variance explained by any single component
      ###     after is less than 1/p. The number of components are also no less than the 
      ###     pre-determined minimum values min.K1 or min.K2.
      N1 <- max( which(cumsum(percent1) < L | percent1 > 1/p ) + 1, min.K1 )
      N2 <- max( which(cumsum(percent2) < L | percent2 > 1/p ) + 1, min.K2 )
      # end of Step 2s
  
      # output reliability ratio
      lambda<-sum(fpca1.value)/sum(fpca1.value+fpca2.value)
  
      # end of Step 2 
  
      # Step 3 
      if(hd==TRUE){
        Cs = t(UT_1)%*%X1+t(UT_2)%*%X2
        # record on a hard-drive. otherwise the memory will be exhaused.
        Phi_B = t(A_B/sqrt(S))%*%Cs
        Phi_W = t(A_W/sqrt(S))%*%Cs
      } 
      # end of Step 3 
        
  
     AN1_B = A_B[, 1:N1]; AN2_W = A_W[, 1:N2]
     lambda1e = fpca1.value[1:N1]*(fpca1.value[1:N1]>0);
     lambda2e = fpca2.value[1:N2]*(fpca2.value[1:N2]>0)
     phi1e = t(AN1_B); phi2e = t(AN2_W)
  
      # Step 4 obtain PC scores in low-dimensional case
      ###     First, calculate the inner product (the cosine of the angles) between 
      ###     level 1 eigenfunctions and level 2 eigenfunctions
      cross.integral <- phi1e
      for(i in 1:N1)
        for(j in 1:N2) 
          cross.integral[i,j] <- sum(fpca1.vectors[,i]* fpca2.vectors[,j]) *tlength/p
  
      ###     Next, calculate the inner product of each centered function with the 
      ###     level 1 or level 2 eigenfunctions
  
      int1 <- phi1e %*% t(resd) * tlength /p
      int2 <- phi2e %*% t(resd) * tlength /p
      
      
  ###     Finally, calculate the principal component scores based on the formulas
  ###     given in the paper.
  s1 <- matrix(0, M*J, K1)
  s2 <- matrix(0, M*J, K2)
  library(MASS)
  design.xi <- ginv( diag(rep(1,K1)) - cross.integral %*% t(cross.integral) )
  for(m in 1:M) {
    resid <- rep(0, K1)
    for(j in 1:J) {
      index <- (m-1) * J + j
      resid <- resid + ( int1[index,] - drop(cross.integral %*% int2[index,]) )/J
    }
    index.m <- ( (m-1) * J + 1 ) : (m*J)
    xi.temp <- design.xi %*% resid
    s1[index.m,] <- matrix(rep(xi.temp, each=J), nrow=J)
    for(j in 1:J) {
      index <- (m-1) * J + j
      s2[index,] <- int2[index,] - drop( t(cross.integral) %*% xi.temp )
    }
  }
  
  
  
  
      # Step 4 ## Obtain PC scores
      dj = rep(1,J)
      C_BW = t(AN1_B)%*% AN2_W
      
      D11 = J*diag(rep(1,N1))
      D12 = kronecker(t(dj),C_BW)
      D21 = t(D12)
      D22 = diag(rep(1,J*N2))
      
      D = rbind(cbind(D11, D12),cbind(D21, D22))
      
      # get by subject singular vectors
      
      U_PC = matrix(1,I*J,I*J)
      for (i in 1:I)
      {
        U_PC[2*(i-1)+1,] = UT_1[i,]
        U_PC[2*(i-1)+2,] = UT_2[i,]
      }
      
      ps_B = matrix(0,N1,I)
      ps_W = matrix(0,J*N2,I)
      
      
      for (i in 1:I)
      {
        rs_i = numeric(0)
        rsb_i = t(AN1_B) %*% sqrt(diag(S))%*% t(U_PC[(2*(i-1)+1):(2*i),])%*%dj
        rs_i = rbind(rs_i, rsb_i)
        for (j in 1:J)
        {
          rsw_ij = t(AN2_W) %*% sqrt(diag(S)) %*% U_PC[(2*(i-1)+j),]
          rs_i = rbind(rs_i, rsw_ij) 
        } 
        
        ps_i = solve(D) %*% rs_i
        ps_B[,i] = ps_i[1:N1,]
        ps_W[,i] = ps_i[(N1+1):(N1+N2*J),]
      }
  
  
  ###     Return the results from multilevel FPCA as a list
  return( list(K1=N1, K2=N2, eigen1=fpca1.value, eigen2=fpca2.value, phi1=phi1e, 
               phi2=phi2e, scores1=ps_B, scores2=ps_W, mu=mu, eta=t(eta) ) )
} 
