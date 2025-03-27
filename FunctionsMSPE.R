# THESE FUNCTIONS ARE USED FOR CONSTRUCTING MULTISTRATUM RESPONSE SURFACE
# DESIGNS BY APPLYING THE METHODS IN XXXXXXX.

# THEY IMPROVE AND EXTEND THE APPROACHES OF TRINCA AND GILMOUR (2001, 2015, 2017)
# AND GILMOUR AND TRINCA (2012). 

# BY APPLYING THE CODE IN STAGES (STRATA), DESIGNS IN COMPLEX LAYOUTS CAN BE BUILD 
# AS MIXING STRUCTURES OF NESTING AND CROSSING UNITS.

#################################################################################################
# FUNCTIONS FOR UNBLOCKED DESIGNS OR FIRST PHASE UNBLOCKED MULTISTRATUM
# FUNCTION FOR GENERATING INITIAL DESIGN
dinicial <- function(cand,Ncand,N,Npar)
{
   ind <- sample(seq(1:Ncand),N,replace=TRUE)
   ind <- ind[order(ind)]
   di <- cand[ind,]
   d <- exp(sum(log(round(eigen(crossprod(di[,-1]), symmetric=TRUE, 
                               only.values=TRUE)$values, 6)))/Npar)
   d[is.nan(d)] <- 0

#   d <- prod(round(eigen(t(di[,-1])%*%di[,-1], symmetric=TRUE, 
#                         only.values=TRUE)$values, 6))^(1/Npar)
   list(di=di,d=d)
}

# CONTROLS FOR NO-SINGULAR INITIAL DESIGN
SampleD <- function(cand,Ncand,N,Npar)
{
   d <- 0
   while(d <= 10^(-6))
     {
      dini <- dinicial(cand,Ncand,N,Npar)
      d <- dini$d
     }
   X <- dini$di
   return(X)
}
  
# FUNCTION TO EXCHANGE ROWS
swap <- function(cand,D,crita,dfPE,W,Ncand,N,Npar,Kappa,probDP, probLP)
{
   KapD <- Kappa[1]; KapDP <- Kappa[2]; KapLP <- Kappa[3]; KapL <- Kappa[4]; KapDF <- Kappa[5]
   df2 <- dfPE
   improve <- 0
   
   for(i in 1:Ncand)
     {
      if(D[1,1]!=cand[i,1])
        {  
         Xout <- D[1,]
         D[1,] <- cand[i,]
         Q <- diag(1,nc=N, nr=N)-matrix(1,nr=N,nc=1)%*%matrix(1,nc=N, nr=1)/N
         M <- t(D[,-c(1,2)])%*%Q%*%D[,-c(1,2)]
         criteC <- criterion(M,W,Npar,KapL, KapLP)
         critD <- criteC$critD
         critL <- criteC$critL
         Critcomp <- 0
         if(critD>10^(-6))
           {                   
            critDPK <- 1
            critLPK <- 1
            critDF <- 1
            if(max(KapDP,KapLP,KapDF)>0)
             {
              critDPK <- 0
              critLPK <- 0
              df2 <- N-nlevels(as.factor(D[,1]))
              critDF <- (N-df2)   # OK, result=v in the paper
              if(df2>0)
                {
                 critDPK <- (critD/qf(probDP,Npar,df2))^KapDP
                 critLPK <- (critL/qf(probLP,1,df2))^KapLP
                }
             }
            Critcomp <- (critD^KapD)*(critDPK)*(critL^KapL)*(critLPK)*(critDF^KapDF)
           }
         if(Critcomp>crita)
           {
            crita <- Critcomp
            improve <- 1
            dfPE <- df2
           } else D[1,] <- Xout
        }
     }
   for(l in 2:N) 
     {
      if(D[l,1]!=D[(l-1),1])
        {
         for(i in 1:Ncand)
           {
            if(D[l,1]!=cand[i,1])
              {
               Xout <- D[l,] 
               D[l,] <- cand[i,]   
               Q <- diag(1,nc=N, nr=N)-matrix(1,nr=N,nc=1)%*%matrix(1,nc=N, nr=1)/N
               M <- t(D[,-c(1,2)])%*%Q%*%D[,-c(1,2)]
               criteC <- criterion(M,W,Npar,KapL, KapLP)
               critD <- criteC$critD
               critL <- criteC$critL
               Critcomp <- 0
               if(critD>10^(-6))
               { 
                critDPK <- 1
                critLPK <- 1
                critDF <- 1
                if(max(KapDP,KapLP,KapDF)>0)
                {
                 df2 <- N-nlevels(as.factor(D[,1]))
                 critDF <- (N-df2)    # OK, result=v in the paper
                 critDPK <- 0
                 critLPK <- 0
                 if(df2>0)
                 {
                   critDPK <- (critD/qf(probDP,Npar,df2))^KapDP
                   critLPK <- (critL/qf(probLP,1,df2))^KapLP
                 }
               }
                Critcomp <- (critD^KapD)*(critDPK)*(critL^KapL)*(critLPK)*(critDF^KapDF)
               }
               if(Critcomp>crita)
                 {
                  crita <- Critcomp
                  improve <- 1
                  dfPE <- df2
                 }else D[l,] <- Xout
              }
           }
        }
     }
  list(D=D,crita=crita,df2=dfPE,improve=improve)
}

# FUNCTION TO DRIVE THE SEARCH
SearchTreat <- function(N, Npar, W, Kappa, Ntries, K, Nlev, Levels=NULL, candExternal=NULL, model, 
                        Fnames, probDP, probLP, MC,
                        parallelTries=FALSE, update.info=FALSE)
{
  Kappa <- Kappa/sum(Kappa)
  KapD <- Kappa[1]; KapDP <- Kappa[2]; KapLP <- Kappa[3]; KapL <- Kappa[4]; KapDF <- Kappa[5]
  if(!is.null(MC)) {if(MC=="Y") probLP <- probLP^(1/Npar)}
  TimeB <- Sys.time()
  if(is.null(Levels)) 
    {
     Levels <- list(1:Nlev[1])
     if(K>1) {for(i in 2:K) Levels <- c(Levels,list(1:Nlev[i]))}
  }
  cand <- candExternal
  if(is.null(candExternal))
  {      
    cand <-as.matrix(expand.grid(Levels))
    cand <- apply(cand,2,CodeX)
  } 
  Ncand <- nrow(cand)
  if(Cubic=='N') cand <- as.matrix(Sphcand(cand)$cand) 
  cand <- TreatLabels(cand,Ncand)$Treat
  dcand <- data.frame(cand)
  colnames(dcand) <- c("Tr", Fnames)
  cand <- cbind(dcand[,1],model.matrix(model, dcand[,-1,drop=FALSE],keep.order=TRUE))
  ########
  mx <- 1:Ntries
  progressr::with_progress({
    p <- progressr::progressor(along = mx)
    if (isTRUE(parallelTries)) {
      designs <- foreach::foreach(m = mx, .options.future = list(seed = TRUE)) %dofuture% 
        {
          if (isTRUE(update.info)) 
            p(message = sprintf("Current iteration: %i out of %i", 
                                m, Ntries))
          X <- SampleD(cand,Ncand,N,Npar)
          Q <- diag(1,nc=N, nr=N)-matrix(1,nr=N,nc=1)%*%matrix(1,nc=N, nr=1)/N
          M <- t(X[,-c(1,2)])%*%Q%*%X[,-c(1,2)]
          criteC <- criterion(M,W,Npar,KapL, KapLP)
          critD <- criteC$critD
          critL <- criteC$critL
          Critcomp <- 0
          df2 <- NULL
          if(critD > 10^(-6))
          {
            critDPK <- 1
            critLPK <- 1
            critDF <- 1
            if(max(KapDP,KapLP,KapDF)>0)
            {
              df2 <- N-nlevels(as.factor(X[,1]))
              critDF <- (N-df2)  # OK, result=v in the paper
              critDPK <- 0
              critLPK <- 0
              if(df2>0)
              {
                critDPK <- (critD/qf(probDP,Npar,df2))^KapDP
                critLPK <- (critL/qf(probLP,1,df2))^KapLP
              }
            }
            Critcomp <- (critD^KapD)*(critDPK)*(critL^KapL)*(critLPK)*(critDF^KapDF)
          }
          improve <- 1
          while(improve==1)
          {
            Xs <- swap(cand,X,Critcomp,df2,W,Ncand,N,Npar,Kappa,probDP, probLP)
            X <- Xs$D
            Critcomp <- Xs$crita
            improve <- Xs$improve
            df2 <- Xs$df2
          }
          list(X = X, value = Critcomp, df2=df2)
        }
      Critall <- sapply(designs, function(x) x$value)
      final <- designs[[which.max(Critall)]]
      Xopt <- final$X
      dfPE <- final$df2
      Critopt <- final$value
    } else {  
            Critall <- numeric() 
            for(m in 1:Ntries) 
            {
             X <- SampleD(cand,Ncand,N,Npar)
             Q <- diag(1,nc=N, nr=N)-matrix(1,nr=N,nc=1)%*%matrix(1,nc=N, nr=1)/N
             M <- t(X[,-c(1,2)])%*%Q%*%X[,-c(1,2)]
             criteC <- criterion(M,W,Npar,KapL, KapLP)
             critD <- criteC$critD
             critL <- criteC$critL
             Critcomp <- 0
             df2 <- NULL
             if(critD > 10^(-6))
             {
              critDPK <- 1
              critLPK <- 1
              critDF <- 1
              if(max(KapDP,KapLP,KapDF)>0)
              {
               df2 <- N-nlevels(as.factor(X[,1]))
               critDF <- (N-df2)  # OK, result=v in the paper
               critDPK <- 0
               critLPK <- 0
               if(df2>0)
               {
                critDPK <- (critD/qf(probDP,Npar,df2))^KapDP
                critLPK <- (critL/qf(probLP,1,df2))^KapLP
               }
              }
             Critcomp <- (critD^KapD)*(critDPK)*(critL^KapL)*(critLPK)*(critDF^KapDF)
             }
            improve <- 1
            while(improve==1)
            {
             Xs <- swap(cand,X,Critcomp,df2,W,Ncand,N,Npar,Kappa,probDP, probLP)
             X <- Xs$D
             Critcomp <- Xs$crita
             improve <- Xs$improve
             df2 <- Xs$df2
            }
            Critall[m] <- Critcomp
            if(m==1)
            {
             Critopt <- Critcomp
             Xopt <- X
             dfPE <- df2
            } else {
                    if(Critcomp > Critopt)
                    {
                     Xopt <- X
                     Critopt <- Critcomp
                     dfPE <- df2
                    }
                   }
            }
           }
    }) 
  dfPE.true <- N-nlevels(as.factor(Xopt[,1]))
  if(is.null(dfPE)&max(KapDP,KapLP,KapDF)<=0) dfPE <- dfPE.true
  DF <- N-dfPE.true
  Q <- diag(1,nc=N, nr=N)-matrix(1,nr=N,nc=1)%*%matrix(1,nc=N, nr=1)/N
  M <- t(Xopt[,-c(1,2)])%*%Q%*%Xopt[,-c(1,2)]
  D <- criterion(M,W,Npar,KapL, KapLP)$critD
  V <- diag(solve(M))
  L <- 1/sum(W*V)
  DP <-  0
  LP <- 0
  if(dfPE.true>0)
  {
    DP <- D/qf(probDP,Npar,dfPE.true)
    LP <- L/qf(probLP,1,dfPE.true)
  }
  
  CritCheck <- c((D^KapD)*(DP^KapDP)*(L^KapL)*(LP^KapLP)*(DF^KapDF),dfPE.true=dfPE.true)
  Xopt <- Xopt[order(Xopt[,1]),]
  Xopt <- Xopt[,3:(K+3-1), drop=FALSE]
  
  Runtime <- difftime(Sys.time(), TimeB, units = "secs")
  CritParts <- c(D=D, DP=DP, LP=LP, L=L, dfPE=dfPE)
  list(Xopt=Xopt,Critopt=Critopt,CritCheck=CritCheck,CritParts=CritParts, Kappa=Kappa,
       Critall=Critall,V=V, Runtime=Runtime)
}

#############################################################################################
# THESE FUNCTIONS ARE COMMON FOR ALL LAYOUTS
# FUNCTION FOR CODING THE LEVELS
CodeX <- function(x)
{
  xc <- (x-(max(x)+min(x))/2)/((max(x)-min(x))/2)
}

TreatLabels <- function(X,n)
 {
  Treat <- matrix(c(1,rep(0,n-1)),n,1)
  Label <- 1
  for (i in 2:n)
    {
     Label <- Label+1
     Treat[i] <- Label
     for (j in 1:(i-1))
      {
       if(min(as.numeric(X[i,]==X[j,]))==1) Treat[i] <- Treat[j]
      }
     if(Treat[i] < Label) Label <- Label-1
    }
   Treat <- as.matrix(cbind(Treat,X))
   list(Treat=Treat)
}


# THIS DOES NOT MAKE SENSE FOR MULTISTRATUM
Sphcand <- function(X)
{
  candc <- t[apply(X^2,1,sum)==K,]
  for(k in (K-1):1)
  {
    candc <- rbind(candc,sqrt(K/k)*t[apply(X^2,1,sum)==k,])
  }
  candc <- rbind(candc,t[apply(X^2,1,sum)==0,])
  list(candc=candc)
}

# FUNCTION FOR CALCULATING THE CRITERION
criterion <- function(M,W,Npar,KapL, KapLP)
{
  critD <- exp(sum(log(round(eigen(M, symmetric=TRUE, 
                                   only.values=TRUE)$values, 6)))/Npar)
  critD[is.nan(critD)] <- 0  
  critL <- 0
  if(critD > 10^(-6)) critL <- ifelse(max(KapL,KapLP)>0,1/sum(W*diag(solve(M))),1)
  list(critD=critD,critL=critL)
}


###################################################################################################
###################################################################################################
# FUNCTIONS FOR BLOCKING (BLOCKED DESIGNS OR FIRST PHASE BLOCKED MULTISTRATUM DESIGNS)

# FUNCTION FOR GENERATING INITIAL DESIGN
dinicial_b <- function(cand,Q,Ncand,N,Npar)
{
  ind <- sample(seq(1:Ncand),N,replace=TRUE)
  di <- cand[ind,]
  V <- t(di[,-1])%*%Q%*%di[,-1] 
  d <- exp(sum(log(round(eigen(V, symmetric=TRUE, 
                               only.values=TRUE)$values, 6)))/Npar)
  d[is.nan(d)] <- 0
  list(di=di,d=d)
}

# CONTROLS FOR NO-SINGULAR INITIAL DESIGN
SampleD_b <- function(cand,Q,Ncand,N,Npar)
{
  d <- 0
  while(d < 10^(-6))
    {
     dini <- dinicial_b(cand,Q,Ncand,N,Npar)
     d <- dini$d
    }
  X <- as.matrix(dini$di)
  list(X=X,d=d)
}


# FUNCTION TO EXCHANGE ROWS
swap_b <- function(cand,D,crita,dfPE,Q,b,Nbloc,W,Ncand,N,Npar,Blocsize, Kappa, probDP, probLP)
{
  KapD <- Kappa[1]; KapDP <- Kappa[2]; KapLP <- Kappa[3]; KapL <- Kappa[4]; KapDF <- Kappa[5]
  df2 <- dfPE
  improve <- 0
  whichblock <- numeric()
  l <- 1
  whichblock <- c(whichblock,1)
  for(i in 1:Ncand)
    {
    if(D[1,1]!=cand[i,1])
      {  
       Xout <- D[1,]    
       D[1,] <- cand[i,]
       M <- t(D[,-1])%*%Q%*%D[,-1]
       criteC <- criterion(M,W,Npar,KapL, KapLP)
       critD <- criteC$critD
       critL <- criteC$critL
       Critcomp <- 0
       if(critD>10^(-6))
         {
          critDPK <- 1
          critLPK <- 1
          critDF <- 1
          if(max(KapDP,KapLP,KapDF)>0)
            {
             MM <- model.matrix(~b+as.factor(D[,1]))
             df2 <- N-rankMatrix(MM)[1]
             critDF <- (N-Nbloc-df2) + 1  # OK, "v" for blocks
             critDPK <- 0
             critLPK <- 0
             if(df2>0)
               {
                critDPK <- (critD/qf(probDP,Npar,df2))^KapDP
                critLPK <- (critL/qf(probLP,1,df2))^KapLP
               }
            }
          Critcomp <- (critD^KapD)*(critDPK)*(critL^KapL)*(critLPK)*(critDF^KapDF) 
         }
       if(Critcomp>crita)
         {
          crita <- Critcomp
          improve <- 1
          dfPE <- df2
         } else {D[1,] <- Xout}
      }
    }
  for(l in 2:N) 
    {
    j <- (l+Blocsize-1)%/%Blocsize
    whichblock <- c(whichblock,j)
    if(D[l,1]!=D[(l-1),1]|whichblock[l]!=whichblock[(l-1)])
    {
     for(i in 1:Ncand)
       {
        if(D[l,1]!=cand[i,1])
          {
           Xout <- D[l,]
           D[l,] <- cand[i,] 
           Q <- diag(1,nc=N, nr=N)-matrix(1,nr=N,nc=1)%*%matrix(1,nc=N, nr=1)/N
           M <- t(D[,-1])%*%Q%*%D[,-1]
           criteC <- criterion(M,W,Npar,KapL, KapLP)
           critD <- criteC$critD
           critL <- criteC$critL
           Critcomp <- 0
           if(critD>10^(-6))
             {
              critDPK <- 1
              critLPK <- 1
              critDF <- 1
              if(max(KapDP,KapLP,KapDF)>0)
              {
               MM <- model.matrix(~b+as.factor(D[,1]))
               df2 <- N-rankMatrix(MM)[1]
               critDF <- (N-Nbloc-df2) + 1  # OK, "v" for blocks
               critDPK <- 0
               critLPK <- 0
               if(df2>0)
               {
                 critDPK <- (critD/qf(probDP,Npar,df2))^KapDP
                 critLPK <- (critL/qf(probLP,1,df2))^KapLP
               }
             }
             Critcomp <- (critD^KapD)*(critDPK)*(critL^KapL)*(critLPK)*(critDF^KapDF) 
             }
           if(Critcomp>crita)
             {
              crita <- Critcomp
              improve <- 1
              dfPE <- df2} else {D[l,] <- Xout
             }
          }
       }
    }
  }
  list(D=D,crita=crita,improve=improve, dfPE=dfPE)
}

# FUNCTION TO DRIVE THE SEARCH FOR BLOCKED DESIGNS
SearchTreat_b <- function(Nbloc, Blocsize, Npar, W, Kappa, Ntries, K, Nlev, Levels=NULL, candExternal=NULL, model, 
                          Fnames, probDP, probLP, MC, parallelTries=FALSE, update.info=FALSE)
{
  Kappa <- Kappa/sum(Kappa)
  KapD <- Kappa[1]; KapDP <- Kappa[2]; KapLP <- Kappa[3]; KapL <- Kappa[4]; KapDF <- Kappa[5]
  if(!is.null(MC)) {if(MC=="Y") probLP <- probLP^(1/Npar)}
  TimeB <- Sys.time()
  N <- Nbloc*Blocsize
  if(is.null(Levels)) 
  {
    Levels <- list(1:Nlev[1])
    if(K>1) {for(i in 2:K) Levels <- c(Levels,list(1:Nlev[i]))}
  }
  cand <- candExternal
  if(is.null(candExternal))
  {      
    cand <-as.matrix(expand.grid(Levels))
    cand <- apply(cand,2,CodeX)
  } 
  Ncand <- nrow(cand)
  if(Cubic=='N') cand <- as.matrix(Sphcand(cand)$cand)
  cand <- TreatLabels(cand,Ncand)$Treat
  
  dcand <- data.frame(cand)
  colnames(dcand) <- Fnames
  cand <- cbind(cand[,1],model.matrix(model, dcand[,-1], keep.order=TRUE))
  b <- as.factor(rep(1:Nbloc,Blocsize))
  B <- as.matrix(model.matrix(~-1+b))
  Q <- diag(rep(1,N))-B%*%solve(crossprod(B))%*%t(B)
  mx <- 1:Ntries
  progressr::with_progress({
    p <- progressr::progressor(along = mx)
    if (isTRUE(parallelTries)) {
      designs <- foreach::foreach(m = mx, .options.future = list(seed = TRUE)) %dofuture% 
        {
          if (isTRUE(update.info)) 
            p(message = sprintf("Current iteration: %i out of %i", 
                                m, Ntries))
          X <- SampleD_b(cand,Q,Ncand,N,Npar)$X
          M <- t(X[,-1])%*%Q%*%X[,-1]
          criteC <- criterion(M,W,Npar,KapL, KapLP)
          critD <- criteC$critD
          critL <- criteC$critL
          Critcomp <- 0
          df2 <- NULL
          if(critD>10^(-6))
          {
            critDPK <- 1
            critLPK <- 1
            critDF <- 1
            if(max(KapDP,KapLP,KapDF)>0)
            {
              MM <- model.matrix(~b+as.factor(X[,1]))
              df2 <- N-rankMatrix(MM)[1]
              critDF <- (N-Nbloc-df2) + 1  # OK, "v" for blocks
              critDPK <- 0
              critLPK <- 0
              if(df2>0)
              {
                critDPK <- (critD/qf(probDP,Npar,df2))^KapDP
                critLPK <- (critL/qf(probLP,1,df2))^KapLP
              } 
            } 
            Critcomp <- (critD^KapD)*(critDPK)*(critL^KapL)*(critLPK)*(critDF^KapDF) 
          }         
          improve <- 1
          while(improve==1)
          {
            Xs <- swap_b(cand,X,crite,df2,Q,b,Nbloc,W,Ncand,N,Npar,Blocsize,Kappa, probDP, probLP)
            X <- Xs$D
            Critcomp <- Xs$crita
            improve <- Xs$improve
            df2 <- Xs$df2
          }
          list(X = X, value = Critcomp, df2=df2)
        }
      Critall <- sapply(designs, function(x) x$value)
      final <- designs[[which.max(Critall)]]
      Xopt <- final$X
      dfPE <- final$df2
      Critopt <- final$value
    } 
    else { 
          Critall <- numeric() 
          for(m in 1:Ntries) 
          {
           X <- SampleD_b(cand,Q,Ncand,N,Npar)$X
           M <- t(X[,-1])%*%Q%*%X[,-1]
           criteC <- criterion(M,W,Npar,KapL, KapLP)
           critD <- criteC$critD
           critL <- criteC$critL
           Critcomp <- 0
           df2 <- NULL
           if(critD>10^(-6))
           {
            critDPK <- 1
            critLPK <- 1
            critDF <- 1
            if(max(KapDP,KapLP,KapDF)>0)
            {
             MM <- model.matrix(~b+as.factor(X[,1]))
             df2 <- N-rankMatrix(MM)[1]
             critDF <- (N-Nbloc-df2) + 1  # OK, "v" for blocks
             critDPK <- 0
             critLPK <- 0
             if(df2>0)
             {
              critDPK <- (critD/qf(probDP,Npar,df2))^KapDP
              critLPK <- (critL/qf(probLP,1,df2))^KapLP
             } 
            } 
            Critcomp <- (critD^KapD)*(critDPK)*(critL^KapL)*(critLPK)*(critDF^KapDF) 
           }         
           improve <- 1
           while(improve==1)
           {
            Xs <- swap_b(cand,X,crite,Q,b,Nbloc,W,Ncand,N,Npar,Blocsize,Kappa,probDP, probLP)
            X <- Xs$D
            Critcomp <- Xs$crita
            improve <- Xs$improve
            df2 <- Xs$df2
           }
           Critall[m] <- Critcomp
           if(m==1)
           {
            Critopt <- Critcomp
            Xopt <- X
            dfPE <- df2
           } else {
                   if(Critcomp > Critopt)
                   {
                    Xopt <- X
                    Critopt <- Critcomp
                    dfPE <- df2
                   }
                  }
          }
    }
  })
  MM <- model.matrix(~b+as.factor(Xopt[,1]))
  dfPE.true <- N-rankMatrix(MM)[1]
  if(is.null(dfPE)&max(KapDP,KapLP,KapDF)<=0) dfPE <- dfPE.true
  Xopt <- Xopt[,-1]
  M <- t(Xopt)%*%Q%*%Xopt
  V <- diag(solve(M))
  D <- criterion(M,W,Npar,KapL, KapLP)$critD
  L <- 1/sum(W*V)
  if(dfPE.true>0)
  {
    DP <- D/qf(probDP,Npar,dfPE.true)
    LP <- L/qf(probLP,1,dfPE.true)
  } 
  DF <- (N-Nbloc-dfPE.true) + 1  # OK, "v" for blocks
  CritCheck <- c((D^KapD)*(DP^KapDP)*(L^KapL)*(LP^KapLP)*(DF^KapDF),dfPE.true)
  Runtime <- difftime(Sys.time(), TimeB, units = "secs")
  CritParts <- c(D=D, DP=DP, LP=LP, L=L, dfPE=dfPE)
  Xopt <- Xopt[,1:K]
  list(Xopt=Xopt,Critopt=Critopt,CritCheck=CritCheck,CritParts=CritParts, Kappa=Kappa,
       Critall=Critall,V=V,Runtime=Runtime)
}

#################################################################################################
#################################################################################################
# BEGIN FUNCTIONS FOR STRATUM BY STRATUM DESIGN

# CONTROLS FOR NO-SINGULAR INITIAL DESIGN (POINT EXCHANGE)
SS_SampleD <- function(cand,Ncand,Q,model,K,Xw,N,Fnames)
{
  logd <- -14
  while(logd <= log(10^(-6)))
    {
     ind <- sample(seq(1:Ncand),N,replace=TRUE)
     X0 <- cbind(Xw,cand[ind,-1])
     dX0 <- data.frame(X0)
     colnames(dX0) <- Fnames
     X <- model.matrix(model,dX0,keep.orde=TRUE)[,-1]
     M <- t(X)%*%Q%*%X
     T2 <- cand[ind,1]
     logd <- sum(log(round(eigen(M, symmetric=TRUE, 
                                  only.values=TRUE)$values,6)))
     logd[is.nan(logd)] <- -14 
    }
  list(X=X,M=M,logd=logd,T2=T2)
}

# FUNCTIONS TO EXPAND WHOLE PLOTS
Expand <- function(x,Swhole)
{
  rep(x,Swhole)
}

# Update M and V    
SSCritUp <- function(i, j, Kappa, improve, Xout, T2=NULL, tout=NULL, W, X, tB, Bsize, V, M, logDeterm, Crit, Npar, dfPE, N, Nbloc, MM, Q=NULL, probDP, probLP)
{
  KapD <- Kappa[1]; KapDP <- Kappa[2]; KapLP <- Kappa[3]; KapL <- Kappa[4]; KapDF <- Kappa[5]
  Critup <- 0
  df2 <- dfPE
  dij <- matrix((X[i,] - Xout),nc=1)
  ftil <- c(tB[j,]%*%X/Bsize)
  aij <- X[i,] - ftil + t(dij)/Bsize
  bij <- Xout - ftil
  Aij <- cbind(t(aij),dij)
  VA <- V%*%Aij
  Bij <- t(cbind(dij,bij))      
  Mup <- M + Aij%*%Bij
  Min <- diag(rep(1,2))+Bij%*%VA
  detMin <- (Min[1,1]*Min[2,2]-Min[1,2]*Min[2,1])
  detMin[is.nan(detMin)] <- 0 
  
  if(detMin>10^(-6))
  {
    # determ <- Determ*detMin
    logdeterm <- logDeterm+log(detMin)
    # Mf <- t(X)%*%Q%*%X
    # logd <- sum(log(eigen(Mf, symmetric=TRUE, 
    #                       only.values=TRUE)$values))
    # print(c(logdeterm-logd)>10^-6)
    Min_inv <- solve(Min)
    Vup <- V-VA%*%Min_inv%*%Bij%*%V
     # dif <- c(Mup-Mf)
     # dif <- dif[abs(dif)>10^-6]
     # print(print(dif))
     # 
    critDPK <- 1
    critLPK <- 1
    critDF <- 1
    if(logdeterm > log(10^(-6)))
    {
      critL <- 1/sum(W*diag(Vup))
      critD <- exp(logdeterm/Npar)
      if(max(KapDP,KapLP,KapDF) > 0)
      {
        df2 <- N-rankMatrix(MM)[1]
        critDF <- (N-Nbloc-df2) + 1  # OK, "v" for blocks
        critDPK <- 0
        critLPK <- 0
        if(df2>0)
          {
           critDPK <- (critD/qf(probDP,Npar,df2))^KapDP
           critLPK <- (critL/qf(probLP,1,df2))^KapLP
          }
      }
      Critup <- (critD^KapD)*(critDPK)*(critL^KapL)*(critLPK)*(critDF^KapDF)
    }
  }
  if(Critup > Crit)
  {
    improve <- 1
    Crit <- Critup
    logDeterm <- logdeterm
    M <- Mup
    V <- Vup
    dfPE <- df2
  }else {
         T2[i] <- tout
         X[i,] <- Xout
  }
  list(X=X,V=V,M=M,logDeterm=logDeterm,Crit=Crit,T2=T2,improve=improve,dfPE=dfPE)
} 

# No updating     
SSCrit <- function(i, j=NULL, Kappa, improve, Xout, T2=NULL, tout=NULL, W, X, tB, Bsize, V, M, logDeterm, Crit, Npar, dfPE, N, Nbloc, MM, Q, probDP, probLP)
{
  KapD <- Kappa[1]; KapDP <- Kappa[2]; KapLP <- Kappa[3]; KapL <- Kappa[4]; KapDF <- Kappa[5]
  df2 <- dfPE
  Mup <- t(X)%*%Q%*%X
  logdeterm <- sum(log(round(eigen(Mup, symmetric=TRUE, 
                                   only.values=TRUE)$values,6)))
  critD <- exp(logdeterm/Npar)
  critD[is.nan(critD)] <- 0 
  Critup <- 0
  if(critD > 10^(-6))
  {
    Vup <- solve(Mup)
    critDPK <- 1
    critLPK <- 1
    critDF <- 1
    critL <- 1/sum(W*diag(Vup))
    if(max(KapDP,KapLP,KapDF) > 0)
    {
      df2 <- N-rankMatrix(MM)[1]
      critDF <- (N-Nbloc-df2) + 1  # OK, "v" for blocks
      critDPK <- 0
      critLPK <- 0
      if(df2 > 0)
      {
        critDPK <- (critD/qf(probDP,Npar,df2))^KapDP
        critLPK <- (critL/qf(probLP,1,df2))^KapLP
      }
    }
    Critup <- (critD^KapD)*(critDPK)*(critL^KapL)*(critLPK)*(critDF^KapDF)
  }
  if(Critup > Crit)
  {
    improve <- 1
    Crit <- Critup
    logDeterm <- logdeterm
    M <- Mup
    V <- Vup
    dfPE <- df2
  }else {
         T2[i] <- tout
         X[i,] <- Xout
  }
  list(X=X,V=V,M=M,logDeterm=logDeterm,Crit=Crit,T2=T2,improve=improve,dfPE=dfPE)
} 

# FUNCTION TO EXCHANGE COORDINATES 
# THIS FUNCTION RELABELS TREATMENTS AFTER EACH COORDINATE CHANGE
SS_Swap_Coord_re <- function(Kappa, Xw, X, Crit, T2=NULL, K2, W, Npar, M, V, tB, b,
                          N, Nbloc, Bsize,logDeterm, M1T, df2, Fnames, Nlev, 
                          Levels, treat.restriction, Q, probDP, probLP, crit.fun.name)
{ 
  improve <- 0
  whichblock <- numeric()
  for(i in 1:N)
  {
    j <- (i+Bsize-1)%/%Bsize
    whichblock <- c(whichblock,j)
      for(k in 1:K2)
      {
        match <- match(X[i,k],Levels[[k]])
        set.lev <- setdiff(1:Nlev[[k]],match)
        for(l in set.lev)
        {
          dif <- l-match
          if(abs(dif)>0)
          { 
            if(sum(X[i,(1:K2)[treat.restriction[k,]==1]])<1)
            {
              Xout <- X[i,]
              X[i,k] <- Levels[[k]][l]
              M2T <- model.matrix(~as.factor(TreatLabels(X[,1:K2],N)$Treat[,1]))[,-1]
              din <- data.frame(matrix(c(Xw[i,],X[i,1:K2]),nr=1))
              X[i,] <- model.matrix(model,din, keep.order=TRUE)[,-1]
              MM <- model.matrix(~ b+M2T+M1T:M2T)
              XUp <- crit.fun.name(i, j, Kappa, improve, Xout, T2=NULL, tout=NULL, W, X, tB, Bsize, V, M, logDeterm, Crit, Npar, df2, N, Nbloc, MM, Q, probDP, probLP)
              X <- XUp$X
              Crit <- XUp$Crit
              improve <- XUp$improve
              M <- XUp$M
              V <- XUp$V
              logDeterm <- XUp$logDeterm
              df2 <- XUp$dfPE
              T2 <- XUp$T2
            }
          }
        }
      }
  }
  list(X=X,Crit=Crit,improve=improve,M=M,V=V,logDeterm=logDeterm,df2=df2,T2=T2)
}

# FUNCTION TO EXCHANGE COORDINATES 
# THIS FUNCTION UPDATES TREATMENT LABEL AFTER EACH COORDINATE CHANGE
SS_Swap_Coord <- function(Kappa, Xw, X, Crit, T2, K2, W, Npar, M, V, tB, b,
                          N, Nbloc, Bsize,logDeterm, M1T, df2, Fnames, Nlev, 
                          Levels, treat.restriction, Q, probDP, probLP, crit.fun.name)
{ 
  improve <- 0
  whichblock <- numeric()
  i <- 1
  whichblock <- c(whichblock,1)
  ###########
  ndif <- 1
  for(k in 1:K2)
  {
    T2.aux <- T2
    ndif <- ifelse(k>1, ndif*Nlev[k-1], ndif)
    match <- match(X[i,k],Levels[[k]])
    set.lev <- setdiff(c(1:Nlev[k]),match)
    for(l in set.lev)
    {
      dif <- l-match
      if(abs(dif)>0)
      { 
        if(sum(X[i,(1:K2)[treat.restriction[k,]==1]])<1)
        {
          tout <- T2[i]
          Xout <- X[i,]
          T2[i] <- T2.aux[i]+dif*ndif
          X[i,k] <- Levels[[k]][l]
          din <- data.frame(matrix(c(Xw[i,],X[i,1:K2]),nr=1))
          X[i,] <- model.matrix(model,din, keep.order=TRUE)[,-1]
          M2T <- model.matrix(~as.factor(T2))[,-1]
          MM <- model.matrix(~ b+M2T+M1T:M2T)
          XUp <- crit.fun.name(i, 1, Kappa, improve, Xout, T2, tout, W, X, tB, Bsize, V, M, logDeterm, Crit, Npar, df2, N, Nbloc, MM, Q, probDP, probLP)
          X <- XUp$X
          T2 <- XUp$T2
          Crit <- XUp$Crit
          improve <- XUp$improve
          M <- XUp$M
          V <- XUp$V
          logDeterm <- XUp$logDeterm
          df2 <- XUp$dfPE
          #                 print(c(tout,improve,T2[i],k,X[i,1:K2]))
        }
      }
    }
  }
  for(i in 2:N)
  {
    ndif <- 1
    j <- (i+Bsize-1)%/%Bsize
    whichblock <- c(whichblock,j)
    if(T2[i]!=T2[(i-1)]|whichblock[i]!=whichblock[(i-1)])
    {
    for(k in 1:K2)
    {
      T2.aux <- T2
      #        improve <- 0
      ndif <- ifelse(k>1, ndif*Nlev[k-1], ndif)
      match <- match(X[i,k],Levels[[k]])
      set.lev <- setdiff(c(1:Nlev[k]),match)
      for(l in set.lev)
      {
        dif <- l-match
        if(abs(dif)>0)
        { 
          if(sum(X[i,(1:K2)[treat.restriction[k,]==1]])<1)
          {
            tout <- T2[i]
            Xout <- X[i,]
            T2[i] <- T2.aux[i]+dif*ndif
            X[i,k] <- Levels[[k]][l]
            din <- data.frame(matrix(c(Xw[i,],X[i,1:K2]),nr=1))
            X[i,] <- model.matrix(model,din, keep.order=TRUE)[,-1]
            M2T <- model.matrix(~as.factor(T2))[,-1]
            MM <- model.matrix(~ b+M2T+M1T:M2T)
            XUp <- crit.fun.name(i, j, Kappa, improve, Xout, T2, tout, W, X, tB, Bsize, V, M, logDeterm, Crit, Npar, df2, N, Nbloc, MM, Q, probDP, probLP)
            X <- XUp$X
            T2 <- XUp$T2
            Crit <- XUp$Crit
            improve <- XUp$improve
            M <- XUp$M
            V <- XUp$V
            logDeterm <- XUp$logDeterm
            df2 <- XUp$dfPE
            #                 print(c(tout,improve,T2[i],k,X[i,1:K2]))
          }
        }
      }
    }
   }
  }
  list(X=X,T2=T2,Crit=Crit,improve=improve,M=M,V=V,logDeterm=logDeterm,df2=df2)
}

# FUNCTION TO DRIVE THE CE SEARCH 
SS_CoordExch <- function(inverse.updating, Kappa, Xwhole, Nbloc, Bsize, K1, K2, W, Npar, Nlev, Levels, Ntries, probDP, probLP,
                         model, Fnames, candExternal, treat.restriction, parallelTries, update.info)
{
  Kappa <- Kappa/sum(Kappa)
  KapD <- Kappa[1]; KapDP <- Kappa[2]; KapLP <- Kappa[3]; KapL <- Kappa[4]; KapDF <- Kappa[5]
  TimeB <- Sys.time()
  N <- Bsize*Nbloc
  K <- K1+K2  
  crit.fun.name <- ifelse(inverse.updating=="Y", SSCritUp, SSCrit)
  if(is.null(Levels)) 
  {
    Levels <- list(CodeX(1:Nlev[1]))
    if(K2>1) {for(i in 2:K2) Levels <- c(Levels,list(CodeX(1:Nlev[i])))}
  }
  cand <- candExternal
  if(is.null(candExternal))
  {      
    cand <-as.matrix(expand.grid(Levels))
#    cand <- apply(cand,2,CodeX)
  } 
  Ncand <- nrow(cand)
  cand.aux <- TreatLabels(cand,Ncand)$Treat

  if(any(treat.restriction>0)) 
    {
    for(k in 1:K2)
      {
       if(sum(treat.restriction[k,])>0)
       {  
        select <- apply(cand.aux[,(2:(K2+1))[treat.restriction[k,]==1]],1,sum)<=1
        cand.aux <- cand.aux[select,]
       }
    }
    Ncand <- nrow(cand.aux)
  }
  cand <- cand.aux
  Xw <- matrix(apply(Xwhole,1,Expand,Bsize),nc=K1,byrow=TRUE)
  T1 <- as.factor(TreatLabels(Xw,N)$Treat[,1])
  M1T <- model.matrix(~ T1)[,-1]
  b <- as.factor(rep(1:Nbloc,rep(Bsize,Nbloc)))
  B <- as.matrix(model.matrix(~-1+b))
  tB <- t(B)
  Q <- diag(rep(1,N))-B%*%solve(tB%*%B)%*%tB
  ####
  mx <- 1:Ntries
  progressr::with_progress({
    p <- progressr::progressor(along = mx)
    if (isTRUE(parallelTries)) {
      designs <- foreach::foreach(m = mx, .options.future = list(seed = TRUE)) %dofuture% 
        {
          if (isTRUE(update.info)) 
            p(message = sprintf("Current iteration: %i out of %i", 
                                m, Ntries))
          Xi <- SS_SampleD(cand,Ncand,Q,model,K,Xw,N, Fnames)
          X <- Xi$X
          M <- Xi$M
          logDeterm <- Xi$logd
          T2 <- Xi$T2
          M2T <- model.matrix(~ as.factor(T2))[,-1]
          V <- solve(M)
          critD <- exp(logDeterm/Npar)
          Critcomp <- 0
          df2 <- NULL
          if(critD>10^(-6))
          {
            critDPK <- 1
            critLPK <- 1
            critDF <- 1
            critL <- 1/sum(W*diag(V))
            if(max(KapDP,KapLP,KapDF)>0)
            {
              MM <- model.matrix(~ b+M2T+M1T:M2T)
              df2 <- N-rankMatrix(MM)[1]
              critDF <- (N-Nbloc-df2) + 1  # OK, "v" for blocks
              critDPK <- 0
              critLPK <- 0
              if(df2>0)
              {
                critDPK <- (critD/qf(probDP,Npar,df2))^KapDP
                critLPK <- (critL/qf(probLP,1,df2))^KapLP
              }
            } 
            Critcomp <- (critD^KapD)*(critDPK)*(critL^KapL)*(critLPK)*(critDF^KapDF)
          }
          improve <- 1
          while(improve==1)
          {      
            Xs <- SS_Swap_Coord(Kappa, Xw, X, Critcomp, T2, K2, W, Npar, M, V, tB, b,
                                N, Nbloc, Bsize, logDeterm, M1T, df2, Fnames, Nlev, Levels, 
                                treat.restriction, Q, probDP, probLP, crit.fun.name)             
            X <- Xs$X
            Critcomp <- Xs$Crit
            logDeterm <- Xs$logDeterm
            M <- Xs$M
            V <- Xs$V
            df2 <- Xs$df2
            T2 <- Xs$T2
            improve <- Xs$improve
          }
          list(X = X, value = Critcomp, df2=df2)
        }
      Critall <- sapply(designs, function(x) x$value)
      final <- designs[[which.max(Critall)]]
      Xopt <- final$X
      dfPE <- final$df2
      Critopt <- final$value
    }        
    else {  
      Critall <- numeric() 
      for(m in 1:Ntries) 
      { 
        Xi <- SS_SampleD(cand,Ncand,Q,model,K,Xw,N, Fnames)
        X <- Xi$X
        M <- Xi$M
        logDeterm <- Xi$logd
        T2 <- Xi$T2
        M2T <- model.matrix(~ as.factor(T2))[,-1]
        V <- solve(M)
        critD <- exp(logDeterm/Npar)
        df2 <- NULL
        Critcomp <- 0
        if(critD>10^(-6))
        {
          critDPK <- 1
          critLPK <- 1
          critDF <- 1
          critL <- 1/sum(W*diag(V))
          if(max(KapDP,KapLP,KapDF)>0)
          {
            MM <- model.matrix(~ b+M2T+M1T:M2T)
            df2 <- N-rankMatrix(MM)[1]
            critDF <- (N-Nbloc-df2) + 1  # OK, "v" for blocks
            critDPK <- 0
            critLPK <- 0
            if(df2>0)
            {
              critDPK <- (critD/qf(probDP,Npar,df2))^KapDP
              critLPK <- (critL/qf(probLP,1,df2))^KapLP
            }
          } 
          Critcomp <- (critD^KapD)*(critDPK)*(critL^KapL)*(critLPK)*(critDF^KapDF)
        }
        improve <- 1
        while(improve==1)
        {      
          Xs <- SS_Swap_Coord(Kappa, Xw, X, Critcomp, T2, K2, W, Npar, M, V, tB, b,
                              N, Nbloc, Bsize, logDeterm, M1T, df2, Fnames, Nlev, Levels, 
                              treat.restriction, Q, probDP, probLP, crit.fun.name)             
          X <- Xs$X
          Critcomp <- Xs$Crit
          logDeterm <- Xs$logDeterm
          M <- Xs$M
          V <- Xs$V
          df2 <- Xs$df2
          T2 <- Xs$T2
          improve <- Xs$improve
        }
        Critall[m] <- Critcomp
        if(m==1)
        {
          Critopt <- Critcomp
          Xopt <- X
          Mopt <- M
          Vopt <- V
          dfPE <- df2
        }else {
          if(Critcomp > Critopt)
          {
            Xopt <- X
            Critopt <- Critcomp
            Mopt <- M
            Vopt <- V
            dfPE <- df2
          }
        }   
      }
    }
  })
  M <- t(Xopt)%*%Q%*%Xopt
  V <- diag(solve(M))
  D <- exp(sum(log(round(eigen(M, symmetric=TRUE, 
                               only.values=TRUE)$values, 6)))/Npar)
  T2 <- TreatLabels(Xopt[,1:K2],N)$Treat[,1]
  M2T <- model.matrix(~as.factor(T2))[,-1]
  MM <- model.matrix(~ b+M2T+M1T:M2T)
  dfPE.true <- N-rankMatrix(MM)[1]
  if(is.null(dfPE)&max(KapDP,KapLP,KapDF)<=0) dfPE <- dfPE.true
  L <- 1/sum(W*V)
  DP <- 1
  LP <- 1
  if(dfPE.true>0)
  {
    DP <- D/qf(probDP,Npar,dfPE.true)
    LP <- L/qf(probLP,1,dfPE.true)
  } 
  DF <- (N-Nbloc-dfPE.true) + 1  # OK, "v" for blocks
  CritCheck <- c((D^KapD)*(DP^KapDP)*(L^KapL)*(LP^KapLP)*(DF^KapDF), dfPE.true=dfPE.true)
  Runtime <- difftime(Sys.time(), TimeB, units = "secs")
  CritParts <- c(D=D, DP=DP, LP=LP, L=L, dfPE=dfPE)
  Xopt <- cbind(Xw,Xopt)[,1:K]
  
  list(Xopt=Xopt,Critopt=Critopt,CritCheck=CritCheck,CritParts=CritParts, Kappa=Kappa,
       Critall=Critall,V=V, Runtime=Runtime)
}

# FUNCTION TO EXCHANGE ROWS 
SS_Swap_Point <- function(Kappa, Xw, cand, X, Crit, W, Npar, M, V, tB, b, Ncand, N, Nbloc, Bsize,
                          logDeterm, M1T, T2, df2, Fnames, Q, probDP, probLP, crit.fun.name)
{
    KapD <- Kappa[1]; KapDP <- Kappa[2]; KapLP <- Kappa[3]; KapL <- Kappa[4]; KapDF <- Kappa[5]
    improve <- 0
    whichblock <- numeric()
    i <- 1
    whichblock <- c(whichblock,1)
    for(jjj in 1:Ncand)
    { 
      if(T2[i]!=cand[jjj,1])
      {
        tout <- T2[i]
        Xout <- X[i,]
        din <- data.frame(matrix(c(Xw[i,],cand[jjj,-1]),nr=1))
        colnames(din) <- Fnames
        X[i,] <- model.matrix(model,din, keep.order=TRUE)[,-1] 
        T2[i] <- cand[jjj,1]
        M2T <- model.matrix(~as.factor(T2))[,-1]
        MM <- model.matrix(~ b+M2T+M1T:M2T) 
        XUp <- crit.fun.name(i,1,Kappa,improve,Xout,T2,tout,W,X,tB,Bsize,V,M,logDeterm,Crit,Npar,df2,N,Nbloc,MM,Q, probDP, probLP)
        X <- XUp$X
        T2 <- XUp$T2
        Crit <- XUp$Crit
        improve <- XUp$improve
        M <- XUp$M
        V <- XUp$V
        logDeterm <- XUp$logDeterm
        df2 <- XUp$dfPE
      }
    }
    for(i in 2:N)
    {
      j <- (i+Bsize-1)%/%Bsize
      whichblock <- c(whichblock,j)
      if(T2[i]!=T2[(i-1)]|whichblock[i]!=whichblock[(i-1)])
      {
        for(jjj in 1:Ncand)
        {
          if(T2[i]!=cand[jjj,1])
          {
            tout <- T2[i]
            Xout <- X[i,]
            din <- data.frame(matrix(c(Xw[i,],cand[jjj,-1]),nr=1))
            colnames(din) <- Fnames
            X[i,] <- model.matrix(model,din, keep.order=TRUE)[,-1] 
            T2[i] <- cand[jjj,1]
            M2T <- model.matrix(~as.factor(T2))[,-1]
            MM <- model.matrix(~ b+M2T+M1T:M2T) 
            XUp <- crit.fun.name(i,j,Kappa,improve,Xout,T2,tout,W,X,tB,Bsize,V,M,logDeterm,Crit,Npar,df2,N,Nbloc,MM,Q,probDP, probLP)
            X <- XUp$X
            T2 <- XUp$T2
            Crit <- XUp$Crit
            improve <- XUp$improve
            M <- XUp$M
            V <- XUp$V
            logDeterm <- XUp$logDeterm
            df2 <- XUp$dfPE
          }
        } 
      }
    }
    list(X=X,T2=T2,Crit=Crit,improve=improve,M=M,V=V,logDeterm=logDeterm,df2=df2)
  }


# FUNCTION TO DRIVE THE SEARCH FOR POINT EXCHANGE
SS_PointExch <- function(inverse.updating, Kappa, Xwhole, Nbloc, Bsize, K1, K2, W, Npar, Nlev, Levels, Ntries,  probDP, probLP, 
                         model, Fnames, candExternal, treat.restriction, parallelTries, update.info)
{
  KapD <- Kappa[1]; KapDP <- Kappa[2]; KapLP <- Kappa[3]; KapL <- Kappa[4]; KapDF <- Kappa[5]
  crit.fun.name <- ifelse(inverse.updating=="Y", SSCritUp, SSCrit)
  TimeB <- Sys.time()
  N <- Bsize*Nbloc
  K <- K1+K2  
  if(is.null(Levels)) 
  {
    Levels <- list(CodeX(1:Nlev[1]))
    if(K2>1) {for(i in 2:K2) Levels <- c(Levels,list(CodeX(1:Nlev[i])))}
  }
  cand <- candExternal
  if(is.null(candExternal))
  {      
    cand <-as.matrix(expand.grid(Levels))
    cand <- apply(cand,2,CodeX)
  } 
  Ncand <- nrow(cand)
  Xw <- matrix(apply(Xwhole,1,Expand,Bsize),nc=K1,byrow=TRUE)
  T1 <- as.factor(TreatLabels(Xw,N)$Treat[,1])
  M1T <- model.matrix(~ T1)[,-1]
  cand <- TreatLabels(cand,Ncand)$Treat
  b <- as.factor(rep(1:Nbloc,rep(Bsize,Nbloc)))
  B <- as.matrix(model.matrix(~-1+b))
  tB <- t(B)
  Q <- diag(rep(1,N))-B%*%solve(tB%*%B)%*%tB
  ########
  mx <- 1:Ntries
  progressr::with_progress({
    p <- progressr::progressor(along = mx)
    if (isTRUE(parallelTries)) {
      designs <- foreach::foreach(m = mx, .options.future = list(seed = TRUE)) %dofuture% 
        {
          if (isTRUE(update.info)) 
            p(message = sprintf("Current iteration: %i out of %i", 
                                m, Ntries))
            Xi <- SS_SampleD(cand,Ncand,Q,model,K,Xw,N, Fnames)
            X <- Xi$X
            M <- Xi$M
            logDeterm <- Xi$logd
            T2 <- Xi$T2
            M2T <- model.matrix(~ as.factor(T2))[,-1]
            V <- solve(M)
            critD <- exp(logDeterm/Npar)
            df2 <- NULL
            Critcomp <- 0
            if(critD>10^(-6))
            {
             critDPK <- 1
             critLPK <- 1
             critDF <- 1
             critL <- 1/sum(W*diag(V))
             if(max(KapDP,KapLP,KapDF)>0)
              {
               MM <- model.matrix(~ b+M2T+M1T:M2T)
               df2 <- N-rankMatrix(MM)[1]
               critDF <- (N-Nbloc-df2) + 1  # OK, "v" for blocks
               critDPK <- 0
               critLPK <- 0
               if(df2>0)
                {
                 critDPK <- (critD/qf(probDP,Npar,df2))^KapDP
                 critLPK <- (critL/qf(probLP,1,df2))^KapLP
                }
              } 
              Critcomp <- (critD^KapD)*(critDPK)*(critL^KapL)*(critLPK)*(critDF^KapDF)
            }
           improve <- 1
           while(improve==1)
           {      
             Xs <- SS_Swap_Point(Kappa,Xw,cand,X,Critcomp,W,Npar,M,V,tB,b,Ncand,N,Nbloc,Bsize,
                                logDeterm,M1T,T2,df2, Fnames, Q, probDP, probLP, crit.fun.name)             
             X <- Xs$X
             T2 <- Xs$T2
             Critcomp <- Xs$Crit
             logDeterm <- Xs$logDeterm
             M <- Xs$M
             V <- Xs$V
             df2 <- Xs$df2
             improve <- Xs$improve
           }
           list(X = X, value = Critcomp, df2=df2)
        }
        Critall <- sapply(designs, function(x) x$value)
        final <- designs[[which.max(Critall)]]
        Xopt <- final$X
        dfPE <- final$df2
        Critopt <- final$value
     }        
    else {  
          Critall <- numeric() 
          for(m in 1:Ntries) 
             { 
               Xi <- SS_SampleD(cand,Ncand,Q,model,K,Xw,N, Fnames)
               X <- Xi$X
               M <- Xi$M
               logDeterm <- Xi$logd
               T2 <- Xi$T2
               M2T <- model.matrix(~ as.factor(T2))[,-1]
               V <- solve(M)
               critD <- exp(logDeterm/Npar)
               Critcomp <- 0
               df2 <- NULL
              if(critD>10^(-6))
              {
               critDPK <- 1
               critLPK <- 1
               critDF <- 1
               critL <- 1/sum(W*diag(V))
               if(max(KapDP,KapLP,KapDF)>0)
               {
                MM <- model.matrix(~ b+M2T+M1T:M2T)
                df2 <- N-rankMatrix(MM)[1]
                critDF <- (N-Nbloc-df2) + 1  # OK, "v" for blocks
                critDPK <- 0
                critLPK <- 0
                if(df2>0)
                {
                 critDPK <- (critD/qf(probDP,Npar,df2))^KapDP
                 critLPK <- (critL/qf(probLP,1,df2))^KapLP
                }
               } 
               Critcomp <- (critD^KapD)*(critDPK)*(critL^KapL)*(critLPK)*(critDF^KapDF)
              }
              improve <- 1
              while(improve==1)
              {  
               Xs <- SS_Swap_Point(Kappa, Xw,cand,X,Critcomp,W,Npar,M,V,tB,b,Ncand,N,Nbloc,Bsize,
                                   logDeterm,M1T,T2,df2, Fnames, Q, probDP, probLP, crit.fun.name)  
               X <- Xs$X
               T2 <- Xs$T2
               Critcomp <- Xs$Crit
               logDeterm <- Xs$logDeterm
               M <- Xs$M
               V <- Xs$V
               df2 <- Xs$df2
               improve <- Xs$improve
              }
             Critall[m] <- Critcomp
             if(m==1)
             {
              Critopt <- Critcomp
              Xopt <- X
              Mopt <- M
              Vopt <- V
              dfPE <- df2
             }else {
                    if(Critcomp > Critopt)
                    {
                     Xopt <- X
                     Critopt <- Critcomp
                     Mopt <- M
                     Vopt <- V
                     dfPE <- df2
                    }
                   }   
             }
         }
  }) 
  M <- t(Xopt)%*%Q%*%Xopt
  V <- diag(solve(M))
  D <- exp(sum(log(round(eigen(M, symmetric=TRUE, 
                               only.values=TRUE)$values, 6)))/Npar)
  T2 <- TreatLabels(Xopt[,1:K2],N)$Treat[,1]
  M2T <- model.matrix(~as.factor(T2))[,-1]
  MM <- model.matrix(~ b+M2T+M1T:M2T)
  dfPE.true <- N-rankMatrix(MM)[1]
  if(is.null(dfPE)&max(KapDP,KapLP,KapDF)<=0) dfPE <- dfPE.true
  L <- 1/sum(W*V)
  DP <- 1
  LP <- 1
  if(dfPE.true>0)
  {
    DP <- D/qf(probDP,Npar,dfPE.true)
    LP <- L/qf(probLP,1,dfPE.true)
  } 
  DF <- (N-Nbloc-dfPE.true) + 1  # OK, "v" for blocks
  CritCheck <- c((D^KapD)*(DP^KapDP)*(L^KapL)*(LP^KapLP)*(DF^KapDF), dfPE.true=dfPE.true)
  Runtime <- difftime(Sys.time(), TimeB, units = "secs")
  CritParts <- c(D=D, DP=DP, LP=LP, L=L, dfPE=dfPE)
  Xopt <- cbind(Xw,Xopt)[,1:K]
  
  list(Xopt=Xopt,Critopt=Critopt,CritCheck=CritCheck,CritParts=CritParts, Kappa=Kappa,
       Critall=Critall,V=V, Runtime=Runtime)
}


#################################################################################
# CALL ALL FUNCTIONS
SS_Search <- function(algorithm, inverse.updating, Kappa, Xwhole, Nbloc, Bsize, K1, K2, W,
                      Npar, Nlev, Levels, Ntries, probDP, probLP, MC, 
                      model, Fnames, candExternal, treat.restriction, parallelTries=FALSE,
                      update.info = FALSE)
{
  algorithm.name <- ifelse(algorithm=="Point", SS_PointExch, SS_CoordExch)
  Kappa <- Kappa/sum(Kappa)
  if(!is.null(MC)) {if(MC=="Y") probLP <- probLP^(1/Npar)}
  result <- algorithm.name(inverse.updating, Kappa, Xwhole, Nbloc, Bsize, K1, K2, W,
                           Npar, Nlev, Levels, Ntries, probDP, probLP, 
                           model, Fnames, candExternal, treat.restriction, 
                           parallelTries, update.info)
  return(result)
}

################################################################################
# FUNCTIONS FOR CROSSED IN PHASE 2 OR HIGHER
# FUNCTION TO EXCHANGE ROWS
SS_Swap_RC <- function(cand,Ncand,model,X,Xw,N,Crit,Q,Rfac,Cfac,T2,M1T,NRow,NCol,Npar,Kappa,W,probDP,probLP)
{
  KapD <- Kappa[1]; KapDP <- Kappa[2]; KapLP <- Kappa[3]; KapL <- Kappa[4]; KapDF <- Kappa[5]
  improve <- 0
  for(i in 1:N)
  {
    for(j in 1:Ncand)
    {
      if(T2[i]!=cand[j,1])
      {  
        tout <- T2[i]
        Xout <- X[i,]
        din <- data.frame(matrix(c(Xw[i,],cand[j,-1]),nr=1))
        colnames(din) <- Fnames
        X[i,] <- model.matrix(model,din, keep.order=TRUE)[,-1] 
        T2[i] <- cand[j,1]
        M <- t(X)%*%Q%*%X
        crite <- criterion(M,W,Npar,KapL,KapLP)
        critD <- crite$critD
        critL <- crite$critL
        Critcomp <- 0
        if(critD>10^(-6))
        {
          critDPK <- 1
          critLPK <- 1
          critDF <- 1
          if(max(KapDP,KapLP,KapDF)>0)
          {
            critDPK <- 0
            critLPK <- 0
            M2T <- model.matrix(~as.factor(T2))[,-1]
            MM <- model.matrix(~Rfac+Cfac+M2T+M1T:M2T) 
            df2 <- N-rankMatrix(MM)[1]
            critDF <- N-NRow-NCol-df2 + 1 
            if(df2 > 0)
            {
              critDPK <- (critD/qf(probDP,Npar,df2))^KapDP
              critLPK <- (critL/qf(probLP,1,df2))^KapLP
            }
            Critcomp <- (critD^KapD)*(critDPK)*(critL^KapL)*(critLPK)*(critDF^KapDF)
          }
          if(Critcomp>Crit)
          {
            Crit <- Critcomp
            improve <- 1
          } else {X[i,] <- Xout
          T2[i] <- tout}
        }
      }
    }
  }
  list(X=X,Crit=Crit,improve=improve,T2=T2,M=M,df2=df2)
}


# FUNCTION TO DRIVE THE SEARCH
# change name SS_Search_RC
SS_Search_RC <- function(Kappa, Xwhole, NRow, NCol, Bsize, K1, K2, W, Npar, Nlev, 
                         Levels, Ntries,  probDP, probLP, model, Fnames, candExternal, 
                         treat.restriction, parallelTries, update.info)
{
  KapD <- Kappa[1]; KapDP <- Kappa[2]; KapLP <- Kappa[3]; KapL <- Kappa[4]; KapDF <- Kappa[5]
  TimeB <- Sys.time()
  K <- K1+K2  
  N <- NRow*NCol*Bsize
  if(!is.null(MC)) {if(MC=="Y") probLP <- probLP^(1/Npar)}
  if(is.null(Levels)) 
  {
    Levels <- list(CodeX(1:Nlev[1]))
    if(K2>1) {for(i in 2:K2) Levels <- c(Levels,list(CodeX(1:Nlev[i])))}
  }
  cand <- candExternal
  if(is.null(candExternal))
  {      
    cand <-as.matrix(expand.grid(Levels))
    cand <- apply(cand,2,CodeX)
  } 
  Ncand <- nrow(cand)
  Xw <- matrix(apply(Xwhole,1,Expand,(NCol*Bsize)),nc=K1,byrow=T)
  T1 <- as.factor(TreatLabels(Xw,N)$Treat[,1])
  M1T <- model.matrix(~ T1)[,-1]
  cand <- TreatLabels(cand,Ncand)$Treat
  R <- rbind(kronecker(diag(rep(1,NRow)),matrix(1,nr=NCol,nc=1)))
  C <- rbind(kronecker(matrix(1,nr=NRow,nc=1),diag(rep(1,NCol))))
  Con <- matrix(c(rep(0,NRow),rep(1,NCol)),nr=1)
  Rfac <- as.factor(rep(1:NRow,rep(NCol,NRow)))
  Cfac <- as.factor(rep(1:NCol,NRow))
  B <- cbind(R,C)
  BtB <- crossprod(B) + crossprod(Con)
  Q <- diag(rep(1,N))-B%*%solve(BtB)%*%t(B)
  ########
  mx <- 1:Ntries
  progressr::with_progress({
    p <- progressr::progressor(along = mx)
    if (isTRUE(parallelTries)) {
      designs <- foreach::foreach(m = mx, .options.future = list(seed = TRUE)) %dofuture% 
        {
          if (isTRUE(update.info)) 
            p(message = sprintf("Current iteration: %i out of %i", 
                                m, Ntries))
          Xi <- SS_SampleD(cand,Ncand,Q,model,K,Xw,N,Fnames)
          X <- Xi$X
          M <- Xi$M
          T2 <- Xi$T2
          M2T <- model.matrix(~as.factor(T2))[,-1]
          crite <- criterion(M,W,Npar,KapL,KapLP)
          critD <- crite$critD
          critL <- crite$critL
          df2 <- NULL
          Critcomp <- 0
          if(critD>10^(-6))
          {
            critDPK <- 1
            critLPK <- 1
            critDF <- 1
            if(max(KapDP,KapLP,KapDF)>0)
            {
              M2T <- model.matrix(~as.factor(T2))[,-1]
              MM <- model.matrix(~Rfac+Cfac+M2T+M1T:M2T)
              df2 <- N-rankMatrix(MM)[1]
              critDF <- N-NRow-NCol-df2 + 1 
              critDPK <- 0
              critLPK <- 0
              if(df2>0)
              {
                critDPK <- (critD/qf(probDP,Npar,df2))^KapDP
                critLPK <- (critL/qf(probLP,1,df2))^KapLP
              }
            }
            Critcomp <- (critD^KapD)*(critDPK)*(critL^KapL)*(critLPK)*(critDF^KapDF) 
          }         
          improve <- 1
          while(improve==1)
          {           
            Xs <- SS_Swap_RC(cand,Ncand,model,X,Xw,N,Critcomp,Q,Rfac,Cfac,T2,M1T,
                             NRow,NCol,Npar,Kappa,W,probDP,probLP)
            X <- Xs$X
            T2 <- Xs$T2
            Crit <- Xs$Crit
            Critcomp <- Xs$Crit
            improve <- Xs$improve
            M <- Xs$M
            df2 <- Xs$df2
          }
          list(X = X, value = Critcomp, df2=df2)
        }
      Critall <- sapply(designs, function(x) x$value)
      final <- designs[[which.max(Critall)]]
      Xopt <- final$X
      dfPE <- final$df2
      Critopt <- final$value
    }        
    else {  
      Critall <- numeric() 
      for(m in 1:Ntries) 
      {
        Xi <- SS_SampleD(cand,Ncand,Q,model,K,Xw,N,Fnames)
        X <- Xi$X
        M <- Xi$M
        T2 <- Xi$T2
        crite <- criterion(M,W,Npar,KapL,KapLP)
        critD <- crite$critD
        critL <- crite$critL
        Critcomp <- 0
        if(critD>10^(-6))
        {
          critDPK <- 1
          critLPK <- 1
          critDF <- 1
          if(max(KapDP,KapLP,KapDF)>0)
          {
            M2T <- model.matrix(~as.factor(T2))[,-1]
            MM <- model.matrix(~Rfac+Cfac+M2T+M1T:M2T)
            df2 <- N-rankMatrix(MM)[1]
            critDF <- N-NRow-NCol-df2 + 1 
            critDPK <- 0
            critLPK <- 0
            if(df2>0)
            {
              critDPK <- (critD/qf(probDP,Npar,df2))^KapDP
              critLPK <- (critL/qf(probLP,1,df2))^KapLP
            }
          }
          Critcomp <- (critD^KapD)*(critDPK)*(critL^KapL)*(critLPK)*(critDF^KapDF) 
        }         
        improve <- 1
        while(improve==1)
        {           
          Xs <- SS_Swap_RC(cand,Ncand,model,X,Xw,N,Critcomp,Q,Rfac,Cfac,T2,M1T,NRow,
                           NCol,Npar,Kappa,W,probDP,probLP)
          X <- Xs$X
          T2 <- Xs$T2
          Crit <- Xs$Crit
          Critcomp <- Xs$Crit
          improve <- Xs$improve
          M <- Xs$M
          df2 <- Xs$df2
        }
        Critall[m] <- Critcomp
        if(m==1)
        {
          Critopt <- Critcomp
          Xopt <- X
          Mopt <- M
          dfPE <- df2
        } else {
          if(Critcomp > Critopt)
          {
            Critopt <- Critcomp
            Xopt <- X
            Mopt <- M
            dfPE <- df2
          }
        }
      }
    }
  })  
  M <- t(Xopt)%*%Q%*%Xopt
  T2 <- TreatLabels(Xopt[,1:K2],N)$Treat[,1]
  M2T <- model.matrix(~as.factor(T2))[,-1]
  MM <- model.matrix(~Rfac+Cfac+M2T+M1T:M2T)
  df2 <- N-rankMatrix(MM)[1]
  # if(is.null(dfPE)) dfPE <- df2
  V <- diag(solve(M))
  D <- exp(sum(log(round(eigen(M, symmetric=TRUE, 
                               only.values=TRUE)$values, 6)))/Npar)
  L <- 1/sum(W*V)
  DP <- 1
  LP <- 1
  if(df2>0)
  {
    DP <- D/qf(probDP,Npar,df2)
    LP <- L/qf(probLP,1,df2)
  } 
  DF <- (N-NRow-NCol-df2) + 1  # OK, "v" for blocks
  CritCheck <- c(crit=(D^KapD)*(DP^KapDP)*(L^KapL)*(LP^KapLP)*(DF^KapDF), df2=df2)
  if(df2==0) {DP <- 0; LP <-0}
  CritParts <- c(D=D, DP=DP, LP=LP, L=L, df2=df2, dfPE=dfPE)
  Xopt <- cbind(Xw,Xopt)[,1:K]
  Runtime <- difftime(Sys.time(), TimeB, units = "secs")
  list(Xopt=Xopt,Critopt=Critopt,CritCheck=CritCheck,CritParts=CritParts, Kappa=Kappa,
       Critall=c(Critall),V=V, Runtime=Runtime)
}

#################################################################################
# FUNCTION FOR FINAL ADJUSTMENT
# 
FinalAdjust <- function(Xfull, N, NRow, Rowsize, NCol, Colsize, PLOTsize, K1, K2, model, 
                        Npar, W, Fnames, Kappa, probDP, probLP)
{
  Kappa <- Kappa/sum(Kappa)
  if(!is.null(MC)) {if(MC=="Y") probLP <- probLP^(1/Npar)}
  KapD <- Kappa[1];  KapDP <- Kappa[2];  KapLP <- Kappa[3];  KapL <- Kappa[4]; KapDF <- Kappa[5]
  Cols <- factor(Xfull[,2])
  Rows <- factor(Xfull[,1])
  Xfull <- Xfull[,-c(1:2)]
  Xw <- Xfull[,1:K1]
  T1 <- TreatLabels(Xw,N)$Treat[,1]
  Z <- model.matrix(~-1+Cols+Rows)
  ZtZ <- crossprod(Z)
  Q <- diag(rep(1,N))-Z%*%solve(ZtZ)%*%t(Z)
  Xfull <- data.frame(Xfull)
  colnames(Xfull) <- Fnames
  X <- model.matrix(model, Xfull, keep.order=TRUE)[,-1]
  M <- t(X)%*%Q%*%X
  logDeterm <- sum(log(round(eigen(M, symmetric = TRUE, only.values = TRUE)$values,6)))
  logDeterm[is.nan(logDeterm)] <- -18
  critD <- exp(logDeterm/Npar)
  Critcomp <- 0
  if(critD>10^(-6))
  { 
    critDPK <- 1
    critLPK <- 1
    critDF <- 1 
    critL <- 1
    if(max(KapL,KapDP,KapLP,KapDF)>0)
    { 
      critL <- criterionAdj(M,W,Npar)$critL
      if(max(KapDP,KapLP,KapDF))
      {
        T2 <- as.factor(TreatLabels(X[,(K1+1):(K1+K2)],N)$Treat[,1])
        MM <- model.matrix(~Z+T1+factor(T2)+factor(T1):T2)
        df2 <- N-rankMatrix(MM)[1]
        critDF <- (N-NCol-NRow-df2) + 1  # OK, "v" for blocks
        critDPK <- 0
        critLPK <- 0
        if(df2>0)
        {
          critDPK <- (critD/qf(probDP,Npar,df2))^KapDP
          critLPK <- (critL/qf(probLP,1,df2))^KapLP
        }
      }
    }
    Critcomp <- (critD^KapD)*(critDPK)*(critLPK)*(critL^KapL)*(critDF^KapDF) 
  }
  Crit_i <- Critcomp
  improve <- 1
  while(improve==1)
  {         
    Xnew <- SwapBlocksPE(N,Z,NRow,Rowsize,NCol,Colsize,PLOTsize,X,Q,W,Npar,Critcomp,T1,
                         K1,K2)  
    X <- Xnew$X
    Critcomp <- Xnew$Crit
    improve <- Xnew$improve
  } 
  XF <- cbind(Rows,Cols,Xw,X[,1:K2]) 
  Xopt <- XF[order(XF[,1],XF[,2]),]
  list(Xopt=Xopt,Critopt=Critcomp, Crit_i=Crit_i) 
}

SwapBlocksPE <- function(N,Z,NRow,Rowsize,NCol,Colsize,PLOTsize,X,Q,W,Npar,Crit,T1,K1,K2)  
{
  improve <- 0
  for(j1 in seq(1,N-Colsize,by=PLOTsize))
  {
    b1 <- (j1%/%Colsize)+1
    for(j2 in seq((N-2*Colsize+1),(N-PLOTsize+1),by=PLOTsize))
    {
      b2 <- j2%/%Colsize+1
      if(b1!=b2)
      {
        if(T1[j1]==T1[j2])
        {
          Xnew <- X
          Xnew[j1:(j1+PLOTsize-1),] <- X[j2:(j2+PLOTsize-1),]
          Xnew[j2:(j2+PLOTsize-1),] <- X[j1:(j1+PLOTsize-1),]
          M <- t(Xnew)%*%Q%*%Xnew
          logDeterm <- sum(log(round(eigen(M, symmetric = TRUE, only.values = TRUE)$values,6)))
          logDeterm[is.nan(logDeterm)] <- -18
          critD <- exp(logDeterm/Npar)
          Critcomp <- 0
          if(critD>10^(-6))
          { 
            critDPK <- 1
            critLPK <- 1
            critDF <- 1 
            critL <- 1
            if(max(KapL,KapDP,KapLP,KapDF)>0)
            { 
              critL <- criterionAdj(M,W,Npar)$critL
              if(max(KapDP,KapLP,KapDF))
              {
                T2 <- as.factor(TreatLabels(X[,(K1+1):(K1+K2)],N)$Treat[,1])
                MM <- model.matrix(~Z+factor(T1)+T2+factor(T1):T2)
                df2 <- N-rankMatrix(MM)[1]
                critDF <- (N-NCol-NRow-df2) + 1  # OK, "v" for blocks
                critDPK <- 0
                critLPK <- 0
                if(df2>0)
                {
                  critDPK <- (critD/qf(probDP,Npar,df2))^KapDP
                  critLPK <- (critL/qf(probLP,1,df2))^KapLP
                }
              }
            }
            Critcomp <- (critD^KapD)*(critDPK)*(critLPK)*(critDF^KapDF)*(critL^KapL)
          }
          if(Critcomp>Crit)
          {
            X <- Xnew
            Crit <- Critcomp
            improve <- 1
          }
        }
      }
    }
  }
  list(X=X,improve=improve,Crit=Crit)
}

criterionAdj <- function(M,W,Npar)
{
  critL <- 1
  if(max(KapL,KapLP)>0)
  {
    Minv <- solve(M)
    critL <- 1/sum(W*diag(Minv))
  }
  list(critL=critL)
}