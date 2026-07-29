#### Basic function ####

# code for MSAR


rpldis = function (n, xmin, alpha, discrete_max = 10000) 
{
  if (length(n) > 1L) 
    n = length(n)
  xmin = floor(xmin)
  u = runif(n)
  if (discrete_max > 0.5) {
    constant = zeta(alpha)
    if (xmin > 1) 
      constant = constant - sum((1:(xmin - 1))^(-alpha))
    cdf = c(0, 1 - (constant - cumsum((xmin:discrete_max)^(-alpha)))/constant)
    dups = duplicated(cdf, fromLast = TRUE)
    if (any(dups)) 
      cdf = cdf[1:which.min(!dups)]
    rngs = as.numeric(cut(u, cdf)) + xmin - 1
    is_na = is.na(rngs)
    if (any(is_na)) 
      rngs[is_na] = floor((xmin - 0.5) * (1 - u[is_na])^(-1/(alpha - 
                                                               1)) + 0.5)
  }
  else {
    rngs = floor((xmin - 0.5) * (1 - u)^(-1/(alpha - 1)) + 
                   0.5)
  }
  rngs
}

getYmat<-function(D, W, N, X, B, Sige, is_nor)
{
  p = nrow(D)
  In = Diagonal(n = N, x = 1)
  if (is_nor)
    eps = rnorm(N*p)
  else
    eps = rt(N*p, df = 5)
  ee = eigen(Sige)
  eps = kronecker(t(sqrt(ee$values)*t(ee$vectors)), In)%*%eps
  Xbeta = as.vector(X%*%B)
  eps1 = eps+Xbeta
  DkW = kronecker(t(D),(W))
  #InvDkW = solve(Diagonal(N*p)-DkW)
  #Yvec = InvDkW%*%eps
  DkW2 = DkW%*%DkW; DkW3 = DkW2%*%DkW; DkW4 = DkW3%*%DkW
  
  #D2 = D%*%D; D3 = D2%*%D; D4 = D3%*%D
  #W2 = W%*%W; W3 = W2%*%W; W4 = W3%*%W
  #DkW2 = kronecker(t(D2),(W2)); DkW3 = kronecker(t(D3),(W3)); DkW4 = kronecker(t(D4),(W4))
  Yvec = eps1 + DkW%*%eps1 + DkW2%*%eps1+ DkW3%*%eps1+ DkW4%*%eps1
  Ymat = matrix(Yvec, ncol = p) 
  return(Ymat)
}



getYmat.gamma<-function(D, W, N, X, B, Sige, is_nor)
{
  p = nrow(D)
  In = Diagonal(n = N, x = 1)
  eps = rgamma(N*p, shape = 1.5, scale = 2)
  eps = eps - mean(eps)
  ee = eigen(Sige)
  eps = kronecker(t(sqrt(ee$values)*t(ee$vectors)), In)%*%eps
  Xbeta = as.vector(X%*%B)
  eps1 = eps+Xbeta
  DkW = kronecker(t(D),(W))
  #InvDkW = solve(Diagonal(N*p)-DkW)
  #Yvec = InvDkW%*%eps
  DkW2 = DkW%*%DkW; DkW3 = DkW2%*%DkW; DkW4 = DkW3%*%DkW
  
  #D2 = D%*%D; D3 = D2%*%D; D4 = D3%*%D
  #W2 = W%*%W; W3 = W2%*%W; W4 = W3%*%W
  #DkW2 = kronecker(t(D2),(W2)); DkW3 = kronecker(t(D3),(W3)); DkW4 = kronecker(t(D4),(W4))
  Yvec = eps1 + DkW%*%eps1 + DkW2%*%eps1+ DkW3%*%eps1+ DkW4%*%eps1
  Ymat = matrix(Yvec, ncol = p) 
  return(Ymat)
}


getYmat.heter<-function(D, W, N, X, B, Sige, is_nor)
{
  p = nrow(D)
  In = Diagonal(n = N, x = 1)
  if (is_nor)
    eps = rnorm(N*p)
  else
    eps = rt(N*p, df = 5)
  d = (rowSums(W>0))
  eps = eps*rep(d, p) #rnorm(N*p, sd = 3)
  ee = eigen(Sige)
  eps = kronecker(t(sqrt(ee$values)*t(ee$vectors)), In)%*%eps
  Xbeta = as.vector(X%*%B)
  eps1 = eps + Xbeta
  DkW = kronecker(t(D),(W))
  #InvDkW = solve(Diagonal(N*p)-DkW)
  #Yvec = InvDkW%*%eps
  DkW2 = DkW%*%DkW; DkW3 = DkW2%*%DkW; DkW4 = DkW3%*%DkW
  
  #D2 = D%*%D; D3 = D2%*%D; D4 = D3%*%D
  #W2 = W%*%W; W3 = W2%*%W; W4 = W3%*%W
  #DkW2 = kronecker(t(D2),(W2)); DkW3 = kronecker(t(D3),(W3)); DkW4 = kronecker(t(D4),(W4))
  Yvec = eps1 + DkW%*%eps1 + DkW2%*%eps1+ DkW3%*%eps1+ DkW4%*%eps1
  Ymat = matrix(Yvec, ncol = p) 
  return(Ymat)
}



getDyadW<-function(N, N1 = 10, delta = 1.2, normalize = T)                                             ### simulate Dyad network W; N1: number of mutual pairs 2N1*N, delta: P((0,1)) = P((1,0)) = 0.5*N^{delta}
{
  A = matrix(0, nrow = N, ncol = N)                                                                    ### use A to store network structure
  
  ######################################### mutual follower ###########################################
  ind = which(upper.tri(A), arr.ind = T)                                                               ### record all the index of upper triangular matrix
  indM = ind[sample(1:nrow(ind), N*N1),]                                                               ### sample N*N1 as mutual pairs in the upper triangular matrix
  A[indM] = 1                                                                                          ### the corresponding links are set to be 1
  A[indM[,2:1]] = 1                                                                                    ### the following matrix is set to be symmetric, as a result, mutual pairs are 2N1*N
  
  ######################################### single relationship #######################################
  ind1 = which(A==0&upper.tri(A), arr.ind = T)                                                         ### record all the zero pairs in the upper triangular matrix
  indS = ind1[sample(1:nrow(ind1), N^delta),]                                                          ### choose N^delta index as single relations
  tmp = sample(1:nrow(indS), floor(N^delta/2))                                                         ### randomly choose 0.5*N^delta to be in the lower triangular matrix
  indS[tmp,] = indS[tmp, 2:1]                                                                          ### change the corresponding index to be inverse
  A[indS] = 1                                                                                          ### the single following relation is set to be 
  diag(A) = 0                                                                                          ### aii = 0
  if (!normalize)
    return(A)
  W = A/rowSums(A)                                                                                     ### W is row-normalized
  W = as(W, "dgCMatrix")
  return(W)
}

getPowerLawW<-function(N, alpha, normalize = T)                                                        ### get power-law network W
{
  Nfollowers = rpldis(N, 1, alpha)                                                                     ### generate N random numbers following power-law(1, alpha): k1-kN
  A = sapply(Nfollowers, function(n) {                                                                 ### for node i, randomly select ki nodes to follow it
    vec = rep(0, N)
    vec[sample(1:N, min(n,N))] = 1
    return(vec)
  })
  diag(A) = 0
  ind = which(rowSums(A)==0)                                                                           ### in case some row sums are zero
  for (i in ind)
  {
    A[i, sample(setdiff(1:N,i), 3)] = 1                                                                ### for those node, randomly select 3 followees
  }
  if (!normalize)
    return(A)
  W = A/rowSums(A)
  W = as(W, "dgCMatrix")
  return(W)
}


getBlockW<-function(N, Nblock, normalize = T)                                                          ### get block network
{
  if (N%%Nblock==0){                                                                                   ### if N mod Nblock is integer
    isDiagList = rep(list(matrix(1, nrow = N/Nblock, ncol = N/Nblock)), Nblock)                        ### obtain the diagnal block list
    mList = rep(list(matrix(rbinom((N/Nblock)^2, size = 1, prob = 0.9*N^{-1}),                       ### generate following relations within the blocks
                            nrow = N/Nblock, ncol = N/Nblock)), Nblock)
  }
  else
  {
    isDiagList = rep(list(matrix(1, nrow = floor(N/Nblock), ncol = floor(N/Nblock))), Nblock-1)        ### if N mod Nblock is not integer
    isDiagList[[length(Nblock)]] = matrix(1, nrow = N%%Nblock, ncol = N%%Nblock)
    
    mList = rep(list(matrix(rbinom(floor(N/Nblock)^2, size = 1, prob = 0.9*N^{-1}),                  ### generate following relations within the blocks
                            nrow = floor(N/Nblock), ncol = floor(N/Nblock))), Nblock-1)
    mList[[Nblock]] = matrix(rbinom(floor(N/Nblock)^2, size = 1, prob = 0.9*N^{-1}),                 ### generate following relations within the blocks
                             nrow = floor(N/Nblock), ncol = floor(N/Nblock))
  }
  isDiag = bdiag(isDiagList)                                                                           ### combine the blocks in matrix
  offDiag = which(isDiag == 0, arr.ind = T)                                                            ### to calculate the index of the off digonal indexes
  mList = lapply(mList, function(M){
    ind = which(rowSums(M)==0)
    if (length(ind)>0)
      M[cbind(ind, sample(1:nrow(M), length(ind)))] = 1
    return(M)
  })
  bA = bdiag(mList)
  bA[offDiag] = rbinom(nrow(offDiag), size = 1, prob = 0.3/N)                                          ### people between blocks have 0.3 prob to follow
  bA = as.matrix(bA)
  upperInd = which(upper.tri(bA), arr.ind = T)
  
  ################ transform bA to be a symmetric matrix ##############################################
  bA[upperInd[,2:1]] = bA[upper.tri(bA)]
  diag(bA) = 0
  
  
  ind = which(rowSums(bA)==0)                                                                          ### in case some row sums are zero
  for (i in ind)
  {
    bA[i, sample(setdiff(1:N,i), 3)] = 1                                                               ### for those node, randomly select 3 followees
  }
  
  if (!normalize)
    return(bA)
  W = bA/rowSums(bA)                                                                                   ### row normalize bA
  W = as(W, "dgCMatrix")
  return(W)
}


MSAR.Lse.Sig1<-function(vecY, W, ww, N, p, q, D, Sige, Omee, 
                        S, tS, S1, m, OmeeS, tSe, SYX, IX, IXbeta)
{
  ### define frequently used matrices 
  In = Diagonal(N, x = 1)
  Inp = Diagonal(N*p, x = 1)
  
  
  ### by eigenvalue decomposition, we have \Sigma_e = Q Lambda t(Q)
  ### here Q = ee$vectors
  ee = eigen(Sige)
  Sige_half = Matrix(t(sqrt(ee$values)*t(ee$vectors))) # QLambda^{1/2}
  ISige_half = kronecker(Sige_half, In)
  tISige_half = kronecker(t(Sige_half), In)
  IOmee_half = kronecker(Omee%*%Sige_half, In)
  
  ### (I-t(D)\otimes W)^{-1}(X\beta)
  SIXbeta = S1%*%IXbeta
  
  ### matrices used
  mtSe = m*tSe
  m2tSe = m^2*tSe
  Sem2tSe = OmeeS%*%(m2tSe)
  
  ### first order: E(Q_{j1j2}^d Q_{k1k2}^d)
  A1 = list()
  A2 = list()
  A3 = list()
  A4 = list()
  G = matrix(0, nrow = N*p, ncol = p^2)
  dM = matrix(0, nrow = N*p, ncol = p^2)
  
  # mS1 = as.matrix(S1)
  tmp2 = tISige_half %*% tSe
  tmp3 = tISige_half %*% OmeeS
  #tmp4 = as.matrix(tmp3 %*% m2tSe)
  tmp4 = tmp3 %*% m2tSe
  
  right = S1%*%ISige_half
  
  # mS1 = as.matrix(S1)
  ### first order: E(Q_{j1j2}^d Q_beta)
  Sig1dx = matrix(0, p^2, p*q)
  SSX = -as.matrix(Sem2tSe) %*% S %*% as.matrix(m2tSe) %*% IX
  
  mOmeeS = as.matrix(OmeeS)
  mSem2tS = as.matrix(Sem2tSe)
  gc()
  #cat("now start computational\n")
  for (j1 in 1:p)
  {
    for (j2 in 1:p)
    {
      #cat(j1, j2, "\n")
      jj = (j2-1)*p+j1 # jj is the row
      Ij2j1 = Matrix(0, nrow = p, ncol = p)
      Ij2j1[j2, j1] = 1
      Ij1j2 = t(Ij2j1)
      
      ### matrices gradients
      Se_g = -kronecker(Omee%*%Ij2j1, W)
      S_g = -kronecker(Ij2j1, W)
      V_g = kronecker(Ij1j2%*%Omee%*%t(D)+D%*%Omee%*%Ij2j1, ww)
      m_g = -m^2*diag(V_g)
      
      ### calculate tr(Mj1j2 Mk1k2)
      A1[[jj]] = m*m_g*tS + m^2*t(S_g)
      A2[[jj]] = as.matrix(OmeeS%*%A1[[jj]])
      A3[[jj]] = m2tSe%*%S_g
      A4[[jj]] = as.matrix(mSem2tS %*% S_g)
      A1[[jj]] = as.matrix(A1[[jj]])
      A3[[jj]] = as.matrix(A3[[jj]])
      
      G1j = tmp3 %*% ((m*m_g)*tSe) %*% SYX
      G2j = tmp3 %*% ((m^2)*t(Se_g)) %*% SYX
      G3j = tISige_half %*% A4[[jj]] %*% vecY
      G[,jj] = G1j[,1] + G2j[,1] + G3j[,1]
      Sig1dx[jj,] = as.numeric(t(SIXbeta)%*%t(S_g)%*%SSX)
      #Sig1dx[jj,] = as.numeric(t(SSX)%*%S_g%*%vecY)
      
      ### calculate diag(Mj1j2)
      left = tmp4 %*% S_g
      dM[,jj] = rowSums(left*t(right)) + rowSums((tISige_half%*%A2[[jj]])* t(IOmee_half))
      gc()
    }
  }
  wdel = as.vector(kronecker((solve(Sige_half)), In)%*%SYX)
  del4 = mean(wdel^4) -  3*mean(wdel^2)^2
  GG = crossprod(G)
  dM2 = crossprod(dM)
  #cat("now no computational\n")
  ### calculate 
  Sig1d1 = matrix(0, p^2, p^2)
  # Sig1d2 = matrix(0, p^2, p^2)
  # Sig1d3 = matrix(0, p^2, p^2)
  # Sig1d4 = matrix(0, p^2, p^2)
  for (j1 in 1:p)
  {
    for (j2 in 1:p)
    {
      for (k1 in 1:p)
      {
        for (k2 in 1:p)
        {
          #cat(j1,j2,k1,k2,'\n')
          jj = (j2-1)*p+j1 # jj is the row
          kk = (k2-1)*p+k1 # kk is the column
          if (kk>=jj)
          {
            #cat(kk, jj, "\n")
            ### E(Q_{j1j2}^dQ_{k1k2}^d)
            ### the same as Sig1d[jj,kk] = as.numeric(tr(M_g[[jj]]%*%t(M_g[[kk]])) + tr(M_g[[jj]]%*%M_g[[kk]])+ t(U_g[[jj]])%*%U_g[[kk]])
            
            ### calculate tr(Mj1j2 Mk1k2)
            Sig1d1[jj,kk] = sum(A2[[jj]]*t(A2[[kk]]))+sum(A4[[jj]]*t(A1[[kk]]))+
              sum(A4[[kk]]*t(A1[[jj]])) + sum(A3[[jj]]*t(A3[[kk]]))
            
            # Sig1d2[jj,kk] = sum(M_g[[jj]] * (t(M_g[[kk]])))
            # 
            # Sig1d3[jj,kk] = sum(M_g[[jj]] * ( M_g[[kk]])) +
            #   sum(U_g[[jj]]*U_g[[kk]]) 
            
            # Sig1d4[jj,kk] = del4*sum(diag(M_g[[jj]])*diag(M_g[[kk]]))
            
            # Sig1d[jj,kk] = sum(M_g[[jj]] * (t(M_g[[kk]]) + M_g[[kk]])) +
            #   sum(U_g[[jj]]*U_g[[kk]]) +
            #   del4*sum(diag(M_g[[jj]])*diag(M_g[[kk]]))
            
          }
          gc()
        }
      }
    }
  }
  Sig1d1[lower.tri(Sig1d1)] = t(Sig1d1)[lower.tri(Sig1d1)]
  Sig1d = Sig1d1 + GG + dM2*del4
  
  ### first order: E(Q_beta Q_beta^\top)
  Sig1x = t(IX)%*%Sem2tSe%*%S%*%(m2tSe)%*%IX
  
  
  ### the exact form of Sig1 = E{Q(theta)Q(theta)^\top}
  Sig1 = matrix(0, nrow = p^2+p*q, ncol = p^2+p*q)
  Sig1[1:p^2, 1:p^2] = Sig1d
  Sig1[1:p^2, (p^2+1):(p^2+p*q)] = Sig1dx
  Sig1[(p^2+1):(p^2+p*q), 1:p^2] = t(Sig1dx)
  Sig1[(p^2+1):(p^2+p*q),(p^2+1):(p^2+p*q)] = as.matrix(Sig1x)
  Sig1 = 4*Sig1
  gc()
  
  return(Sig1)
}

ifelse<-function(x, val1, val0)
{
  if (x)
  {
    return(val1)
  }else
  {
    return(val0)
  }
}

### calculate the gradient and hessian matrix of Q(theta)
lse.grad.hessian<-function(D, Ymat, vecY, W, ww, N, p, X, B, Sige, infer = F, old = F)
{
  p = ncol(Ymat)
  q = ncol(X)
  if (is.null(q)){
    q <- 1
  }
  ############################################################
  #### matrices
  
  tW = t(W) # ww = tW%*%W
  
  ### identity matrix
  I = Diagonal(N*p, x = 1)
  Ip = Diagonal(p, x = 1)
  In = Diagonal(N, x = 1)
  Inp = Diagonal(N*p, x = 1)
  
  
  IX = kronecker(Ip, X)
  IXbeta = as.vector(X%*%B)
  
  Omee = solve(Sige)
  IOmee = kronecker(Omee, In)
  
  S = I - kronecker(t(D), W); tS = t(S)
  SYX = S%*%vecY - IXbeta # This \mE
  OmeeS = IOmee%*%S
  tSe = t(OmeeS)
  OmeSYX = tSe%*%SYX
  
  Ome = tSe%*%S
  m = 1/diag(Ome)
  
  MY = m*OmeSYX # this is F
  
  dww = diag(ww)
  vMY = as.vector(MY)
  mtSe = m*tSe
  
  # left1 = matrix(-m^2*vMY, ncol = p)
  # left2 = matrix(m*vMY, ncol = p)
  # left3 = matrix(OmeeS%*%(m*vMY), ncol = p)
  
  ### calculate gradient
  
  right1 = matrix(tSe%*%SYX, ncol = p)
  right2 = matrix(SYX, ncol = p)
  right3 = Ymat
  
  ###
  left1b = matrix(m^2*vMY, ncol = p)
  left2b = matrix(-m*vMY, ncol = p)
  
  
  #grad_D = rep(0, p^2)
  F_d = matrix(0, nrow = N*p, ncol = p^2)
  Q_db2 = Matrix(0, nrow = p^2, ncol = p*q)
  
  for (j1 in 1:p)
  {
    for (j2 in 1:p)
    {
      #cat(j1, j2, "\n")
      Ij2j1 = Matrix(0, nrow = p, ncol = p)
      Ij2j1[j2, j1] = 1
      Ij1j2 = t(Ij2j1)
      
      Vj1j21 = Ij1j2%*%Omee%*%t(D)+D%*%Omee%*%Ij2j1
      # right1j = (dww*right1)%*%diag(diag(Vj1j21))
      
      F_d1j = ifelse(p==1, -m^2*as.vector((dww*right1)%*%(Vj1j21)),
                     -m^2*as.vector((dww*right1)%*%diag(diag(Vj1j21))))
      
      # right2j = -t(W)%*%right2%*%t(Ij1j2%*%Omee)
      F_d2j = -m*as.vector(t(W)%*%right2%*%t(Ij1j2%*%Omee))
      
      # right3j = -W%*%right3%*%t(Ij2j1)
      F_d3j = -mtSe%*%as.vector(W%*%right3%*%t(Ij2j1))
      
      
      ii = (j2-1)*p+j1 # j2 is the column
      #grad_D[ii] = sum(left1*right1j) + sum(left2*right2j)+ sum(left3*right3j)
      F_d[,ii] = F_d1j + F_d2j + F_d3j[,1]
      
      ### for the calculation of F_db
      ### 1st term: 
      ### 1.1: diag(Vj1j21)*Omee \otimes (dww*X)
      ### 1.2: -diag(Vj1j21)*(Omee%*%t(D))\otimes (dww*W%*%X)
      Q_db1j = t(dww*X)%*%left1b%*%(diag(Vj1j21)*Omee)-
        t(dww*t(W)%*%X)%*%left1b%*%(diag(Vj1j21)*(D%*%Omee))
      ### 2nd term: -t(Omee%*%Ij2j1) \otimes t(W)%%X
      Q_db2j = -t(tW%*%X)%*%left2b%*%t(Omee%*%Ij2j1)
      Q_db2[ii,] = as.vector(Q_db1j + Q_db2j)*2
    }
  }
  
  F_b = -m*tSe%*%IX 
  
  ### calculate Q_{j1j2}^d and Q_{\beta} by (2.1)
  #grad_D = 2*grad_D
  grad_D = 2*colSums(vMY*F_d)
  grad_beta = 2*colSums(vMY*F_b)
  grad = c(grad_D, grad_beta)
  
  
  #############################################################
  #### 2nd order: Hessian
  
  ### modify the hessian matrix computation
  
  ### first term
  left11 = matrix(2*m^3*vMY, ncol = p)
  left12 = matrix(-m^2*vMY, ncol = p)
  right1 = matrix(OmeSYX, ncol = p)
  
  ### second term
  #### for Omega_k1k2
  left2 = matrix(-m^2*vMY, ncol = p)
  right21 = Ymat
  right22 = X%*%B
  
  ### 3rd term: just as the 2nd one
  ### 4th term
  left4 = matrix(m*vMY, ncol = p)
  right4 = Ymat
  
  Q_dd2 = matrix(0, p^2, p^2)
  
  system.time({
    for (j1 in 1:p)
    {
      for (j2 in 1:p)
      {
        for (k1 in 1:p)
        {
          for (k2 in 1:p)
          {
            jj = (j2-1)*p+j1 # jj is the row
            kk = (k2-1)*p+k1 # kk is the column
            if (jj <= kk)
            {
              ii = (kk-1)*p^2 + jj
              
              
              #j1 = 1;j2 = 1; k1= 1;k2 = 1
              ### get Ij1j2 Ij2j1 and Ik1k2 Ik2k1
              Ij2j1 = Matrix(0, nrow = p, ncol = p)
              Ij2j1[j2, j1] = 1
              Ij1j2 = t(Ij2j1)
              Ik2k1 = Matrix(0, nrow = p, ncol = p)
              Ik2k1[k2, k1] = 1
              Ik1k2 = t(Ik2k1)
              IIk2k1 = kronecker(Ik2k1, In)
              
              Vj1j21 = Ij1j2%*%Omee%*%t(D)+D%*%Omee%*%Ij2j1
              Vk1k21 = Ik1k2%*%Omee%*%t(D)+D%*%Omee%*%Ik2k1
              Vjk = Ij1j2%*%Omee%*%Ik2k1+Ik1k2%*%Omee%*%Ij2j1
              
              ### 1st term
              right11jk = ifelse(p==1, (dww^2*right1)%*%(Vj1j21)*diag(Vk1k21),
                                 (dww^2*right1)%*%diag(diag(Vj1j21)*diag(Vk1k21)))
              right12jk = ifelse(p==1, (dww*right1)%*%(Vjk), (dww*right1)%*%diag(diag(Vjk)))
              
              ### 2nd term
              right211jk = -(dww*tW)%*%right21%*%t(diag(Vj1j21)*(Ik1k2%*%Omee))
              right212jk = -(dww*W)%*%right21%*%t(diag(Vj1j21)*(Omee%*%Ik2k1))
              right213jk = (dww*ww)%*%right21%*%t(diag(Vj1j21)*Vk1k21)
              right22jk = (dww*tW)%*%right22%*%t(diag(Vj1j21)*t(Omee%*%Ik2k1))
              
              ### 3rd term
              right211kj = -(dww*tW)%*%right21%*%t(diag(Vk1k21)*(Ij1j2%*%Omee))
              right212kj = -(dww*W)%*%right21%*%t(diag(Vk1k21)*(Omee%*%Ij2j1))
              right213kj = (dww*ww)%*%right21%*%t(diag(Vk1k21)*Vj1j21)
              right22kj = (dww*tW)%*%right22%*%t(diag(Vk1k21)*t(Omee%*%Ij2j1))
              
              ### 4th term
              right4jk = ww%*%right4%*%t(Vjk)
              
              Q_dd2[jj,kk] = sum(left11*right11jk)+sum(left12*right12jk)+
                sum(left2*right211jk) + sum(left2*right212jk) + sum(left2*right213jk)+ sum(left2*right22jk)+
                sum(left2*right211kj) + sum(left2*right212kj) + sum(left2*right213kj)+ sum(left2*right22kj)+
                sum(left4*right4jk)
              
            }
          }
        }
        
      }
    }
  })
  
  Q_dd2[lower.tri(Q_dd2)] = t(Q_dd2)[lower.tri(Q_dd2)]
  
  ### calculate (2.8)--(2.10)
  Q_dd = 2*t(F_d)%*%F_d + Q_dd2*2
  Q_db = 2*t(F_d)%*%F_b+Q_db2
  Q_bb = 2*t(F_b)%*%F_b
  
  
  ### Hessian Matrix
  Q_h = Matrix(0, nrow = p^2+p*q+p^2, ncol = p^2+p*q+p^2)
  ind_d = 1:p^2; ind_b = (p^2+1):(p^2+p*q)
  
  Q_h[ind_d, ind_d] = Q_dd
  Q_h[ind_b, ind_b] = Q_bb
  
  Q_h[ind_d, ind_b] = Q_db
  Q_h[ind_b, ind_d] = t(Q_db)
  
  ind0 = c(ind_d, ind_b)
  Q_h = Q_h[ind0, ind0]
  
  if (infer) # if inference is required
  {
    tS = t(S)
    ### fast calculation of (I-t(D)\otimes W)^{-1}\approx I + dw+dw^2+dw^3...
    dw = kronecker(t(D), W)
    W2 = (W%*%W)
    W3 = W2%*%(W)
    dw2 = kronecker(t(D%*%D), W2)
    dw3 = kronecker(t(D%*%D%*%D),W3)
    S1 = Inp + dw + dw2 + dw3
    
    if (old)
    {
      Sig1 = MSAR.Lse.Sig1.old(vecY, W, ww, N, p, q, D, Sige, S, m, OmeeS, tSe, Ome, Omee, IX, IXbeta, SYX)
    }else{
      Sig1 = MSAR.Lse.Sig1(vecY, W, ww, N, p, q, D, Sige, Omee,
                           S, tS, S1, m, OmeeS, tSe, SYX, IX, IXbeta)
    }
    
    inv_Qh = solve(Q_h)
    covl = inv_Qh%*%Sig1%*%inv_Qh
    return(list(grad = grad, hessian = Q_h, cov = covl))
  }
  
  gc()
  return(list(grad = grad, hessian = Q_h, obj = mean(vMY^2)))
}

#### This file includes functions which calculate the newton-raphson algorithm of LSE


### calculate the trace of a matrix
tr<-function(M)
{
  return(sum(diag(M)))
}



MSAR.Lse<-function(Ymat, X, W, D, B, Sige, verbose = F, infer = T, old = F, tol = 10^{-3})  ### LSE algorithm
{
  tim = system.time({
    N = nrow(W); p = ncol(Ymat); vecY = as.vector(Ymat); ww = crossprod(W); q = nrow(B)
    Omee = solve(Sige)
    Del1 = 1
    n.iter = 1
    flag = FALSE
    ind_d = 1:p^2; ind_b = (p^2+1):(p^2+p*q); ind_e0 = which(lower.tri(Sige,diag = T))
    obj0 = -1
    while (mean(abs(Del1))>tol & n.iter<=100)
    {
      if (verbose)
        cat(mean(abs(Del1)), " ")
      if (any(abs(diag(D))>1))
      {
        D = Diagonal(n = p, x = runif(p, -0.5, 0.5))
        B = Matrix(runif(p*q), nrow = q, ncol = p)
      }
      vecTheta = c(as.vector(D), as.vector(B))
      
      
      system.time({grad_hessian = lse.grad.hessian(D, Ymat, vecY, W, ww, N, p, X, B, Sige)})
      
      hessian = grad_hessian$hessian
      Del = solve(hessian)%*%grad_hessian$grad
      if (any(abs(diag(D))>1))
        vecTheta = vecTheta - Del*0.1
      else
        vecTheta = vecTheta - Del
      D = Matrix(matrix(vecTheta[ind_d], nrow = p))
      B = Matrix(matrix(vecTheta[ind_b], nrow = q))
      E = Ymat - W%*%Ymat%*%D - X%*%B
      Sige = t(E)%*%E/N
      Del1 = abs(grad_hessian$obj-obj0)
      obj0 = grad_hessian$obj
      n.iter = n.iter+1
    }
  })
  if (infer)
  {
    grad_hessian_cov = lse.grad.hessian(D, Ymat, vecY, W, ww, N, p, X, B, Sige, infer = infer, old = old)
  }else{
    grad_hessian_cov = NULL
  }
  
  return(list(D = D, B = B, Sige = Sige, theta = vecTheta, cov = grad_hessian_cov$cov,
              iter = n.iter, time = tim[3]))
}



# Function to generate random skew-symmetric matrix
rSkewSym <- function(d,kappa=5*d^(1.5)) {
  m2 <- rep(1,d) / sqrt(sum((rep(1,d))^2))
  m1 <- rmovMF(1, kappa*m2)
  opt_trans <- log_rotation(m2,m1) 
  M <- opt_trans$L
  return(M)
}



## optimal transport
unit_norm <- function(x) x / sqrt(sum(x^2))

geo_dist <- function(x, y) {
  z <- sum(x * y)
  z <- max(min(z, 1), -1)
  acos(z)
}

rand_on_sphere <- function(d) {
  z <- rnorm(d)
  unit_norm(z)
}

# Orthonormal basis of span{x,y}
span_basis_xy <- function(x, y) {
  u1 <- unit_norm(x)
  proj <- sum(y * u1) * u1
  v <- y - proj
  nv <- sqrt(sum(v^2))
  if (nv < 1e-12) {
    tmp <- rnorm(length(x))
    tmp <- tmp - sum(tmp * u1) * u1
    u2 <- unit_norm(tmp)
  } else {
    u2 <- v / nv
  }
  list(u1 = u1, u2 = u2)
}



log_rotation <- function(x, y, tol = 1e-12) {
  x <- as.numeric(x); y <- as.numeric(y)
  stopifnot(length(x) == length(y))
  
  # ensure unit length
  x <- x / sqrt(sum(x^2))
  y <- y / sqrt(sum(y^2))
  
  cxy <- sum(x * y)
  cxy <- max(min(cxy, 1), -1)
  
  # identical or antipodal cases
  if (1 - abs(cxy) < tol) {
    if (cxy > 0) {
      return(list(L = matrix(0, length(x), length(x)), theta = 0, u1 = x, u2 = NA))
    } else {
      # antipodal: choose any unit u2 perpendicular to x
      j <- which.min(abs(x))
      e <- rep(0, length(x)); e[j] <- 1
      v <- e - sum(e * x) * x
      vn <- sqrt(sum(v^2))
      if (vn < tol) {
        if (length(x) < 2) stop("Need dimension >= 2 for antipodal case")
        v <- c(-x[2], x[1], if (length(x) > 2) rep(0, length(x) - 2) else NULL)
        vn <- sqrt(sum(v^2))
      }
      u2 <- v / vn
      theta <- pi
      Q <- tcrossprod(u2, x) - tcrossprod(x, u2)
      return(list(L = theta * Q, theta = theta, u1 = x, u2 = u2))
    }
  }
  
  # general case
  u1 <- x
  v <- y - cxy * u1
  u2 <- v / sqrt(sum(v^2))
  theta <- acos(cxy)
  Q <- tcrossprod(u2, u1) - tcrossprod(u1, u2)
  
  list(L = theta * Q, theta = theta, u1 = u1, u2 = u2)
}


exp_from_log_rotation_predict <- function(q,theta) {
  th <- theta
  L <- q
  if (is.null(th)) {
    eig <- eigen(L, symmetric = FALSE)
    V <- eig$vectors; D <- diag(exp(eig$values))
    Re(V %*% D %*% solve(V))
  } else {
    I <- diag(nrow(L))
    K <- L / max(th, 1e-16)
    I + sin(th) * K + (1 - cos(th)) * (K %*% K)
  }
}


exp_from_log_rotation <- function(L,a) {
  th <- L$theta
  L <- L$L
  if (is.null(th)) {
    eig <- eigen(L, symmetric = FALSE)
    V <- eig$vectors; D <- diag(exp(eig$values))
    Re(V %*% D %*% solve(V))
  } else {
    I <- diag(nrow(L))
    K <- a * L / max(th, 1e-16)
    I + sin(th) * K + (1 - cos(th)) * (K %*% K)
  }
}


#### intrinsic mean
intrinsic_mean_sphere <- function(y, tol = 1e-8, max_iter = 1000) {
  n <- nrow(y)
  d <- ncol(y)
  
  # Initialize with extrinsic mean
  mu <- colSums(y)
  mu <- mu / sqrt(sum(mu^2))
  
  for (iter in 1:max_iter) {
    tangent_sum <- rep(0, d)
    for (i in 1:n) {
      theta <- acos(sum(mu * y[i,]))
      if (theta > 1e-12) {
        tangent <- (theta / sin(theta)) * (y[i,] - cos(theta) * mu)
        tangent_sum <- tangent_sum + tangent
      }
    }
    tangent_mean <- tangent_sum / n
    norm_t <- sqrt(sum(tangent_mean^2))
    
    if (norm_t < tol) break
    mu <- cos(norm_t) * mu + sin(norm_t) * (tangent_mean / norm_t)
  }
  
  return(mu)
}

inv_clr <- function(clr_row) {
  x <- exp(clr_row)
  x / sum(x)
}

PSSAR_direct_sim_function <- function(seed_set,sample_size,dim_set,neighbor_set, rho0,dim_red='PCA'){
  set.seed(seed_set)
  # Parameters
  n <- sample_size      # number of observations
  d <- dim_set          # dimension of skew-symmetric matrices
  
  k <- neighbor_set         # number of non-zero neighbors per row
  
  W <- matrix(0, n, n)
  
  for (i in 1:n) {
    # choose k neighbors for node i
    neighbors <- sample(setdiff(1:n, i), k)
    
    for (j in neighbors) {
      # assign a random weight
      w <- runif(1)
      W[i, j] <- w
      W[j, i] <- w   # enforce symmetry
    }
  }
  
  # Row-standardize
  row_sums <- rowSums(W)
  W <- W / row_sums
  
  # Mean matrix (skew-symmetric)
  q_bar <- 0 #rSkewSym(d)
  
  # Noise terms
  epsilon <- lapply(1:n, function(i) rSkewSym(d))
  epsilon_mean <- Reduce("+", epsilon) / n
  
  # Solve SAR-like system for q_i
  I_n <- diag(n)
  
  E <- epsilon
  
  
  ## faster version
  d2 <- d * d
  
  # Stack E into 3D → matrix
  E_arr <- simplify2array(E)      # d x d x n
  E_center <- sweep(E_arr, c(1,2), epsilon_mean)
  E_mat <- matrix(aperm(E_center, c(3,1,2)), n, d2)  # n × d²
  
  # Q_resid
  I_n <- diag(n)
  Q_resid <- solve(I_n - rho0 * W)
  
  # Q_i (all at once)
  Q_mat <- Q_resid %*% E_mat
  
  # Add q_bar
  Q_mat <- sweep(Q_mat, 2, as.vector(q_bar), "+")
  
  # Gram matrix
  q_hat_vec <- colMeans(Q_mat)
  Q_centered <- sweep(Q_mat, 2, q_hat_vec)
  
  G_hat <- Q_centered %*% t(Q_centered)
  
  ## cost function
  objective <- function(rho, W, G_hat) {
    n <- nrow(W)
    S <- diag(n) - rho * W
    val <- sum(diag(t(S) %*% W %*% S %*% G_hat))
    return(val^2)
  }
  
  rho_hat <- optimize(objective, interval=c(-0.99, 0.99), W=W, G_hat=G_hat)$minimum
  
  
  ##hypothesis test
  
  if (dim_red=='PCA'){
    
    # Convert back to list of d × d matrices
    Q <- lapply(1:n, function(i) matrix(Q_mat[i, ], d, d))
    q_hat <- Reduce("+", Q) / n
    
    E_est_list <- list()
    
    for (i in 1:n) {
      E_est <- Q[[i]]-q_hat
      for (j in 1:n) {
        E_est <- E_est - rho_hat * W[i,j] * (Q[[j]]-q_hat)
      }
      E_est_list[[i]] <- E_est
    }
    
    
    # Step 1. Vectorize each epsilon_i
    U <- t(sapply(E_est_list, function(E) as.vector(E)))  # n x d^2
    
    # Step 2. Center columns
    Uc <- scale(U, center = TRUE, scale = FALSE)
    
    # Step 3. Truncated SVD (to avoid forming huge covariance)
    if (dim_set<=10){
      Sigma_hat <- cov(Uc)  # d^2 x d^2 matrix
      eig <- eigen(Sigma_hat, symmetric = TRUE)
      eigs <- eig$values
    } else {
      ##max_m = dim_set
      max_m = d
      m <- min(max_m, n - 1)
      sv <- irlba(Uc, nv = m, nu = 0)
      eigs <- (sv$d^2) / n   # eigenvalues of sample covariance
    }
    
    
    
    # Step 4. Choose k by variance explained
    cumvar <- cumsum(eigs) / sum(eigs)
    k <- which(cumvar >= 0.9)[1]
    
    
    # Step 5. Compute varpi’s using top-k eigenvalues
    lambda <- eigs[1:k]
    
    # Step 5. Compute varpi_1, varpi_3, varpi_4
    varpi_1 <- (sum(lambda))^2
    varpi_2 <- sum(lambda^2)
    sigma_cst <- varpi_2 / varpi_1
    sigma_est <- n * sigma_cst * sum(diag(W %*% (W+t(W)))) / (sum(diag((W+t(W)) %*% W %*% solve(I_n - rho_hat * W))))^2
  } else if (dim_red=='bootstrap'){
    
    Q_resid <- I_n
    
    ## faster
    ## ---- Precompute E_est_list ----
    # Q_centered_mat = Q_mat - q_hat_vec
    Q_centered_mat <- sweep(Q_mat, 2, q_hat_vec)
    
    # Compute E_est_mat under the null
    E_est_mat <- Q_centered_mat ##- rho_hat * (W %*% Q_centered_mat)
    
    # Convert back to list of d × d matrices
    E_est_list <- lapply(1:n, function(i) matrix(E_est_mat[i, ], d, d))
    
    
    
    rho_est_boots <- rep(0,500)
    for (boots_ind in 1:500) {
      print(boots_ind)
      
      ## faster 
      # 2. Bootstrap residual sample
      idx <- sample.int(n, n, TRUE)
      eps_boots_mat  <- E_est_mat[idx, , drop=FALSE]
      eps_boots_mean <- colMeans(eps_boots_mat)
      
      # 3. Compute Q_i
      eps_center_mat <- sweep(eps_boots_mat, 2, eps_boots_mean)
      Q_mat <- sweep(Q_resid %*% eps_center_mat, 2, q_hat_vec, "+")
      
      # 4. Compute G_hat
      q_hat_boots_vec <- colMeans(Q_mat)
      Q_centered_mat <- sweep(Q_mat, 2, q_hat_boots_vec)
      G_hat <- Q_centered_mat %*% t(Q_centered_mat)
      
      ## cost function
      objective <- function(rho, W, G_hat) {
        n <- nrow(W)
        S <- diag(n) - rho * W
        val <- sum(diag(t(S) %*% W %*% S %*% G_hat))
        return(val^2)
      }
      
      rho_est_boots[boots_ind] <- optimize(objective, interval=c(-0.99, 0.99), W=W, G_hat=G_hat)$minimum
    }
    
    ##sigma_est <- var(rho_est_boots) * n
    ##plot(density(rho_est_boots))
    ci <- quantile(rho_est_boots, c(0.025, 0.975))
  } else{
    
    # Convert back to list of d × d matrices
    Q <- lapply(1:n, function(i) matrix(Q_mat[i, ], d, d))
    q_hat <- Reduce("+", Q) / n
    
    varpi_1 <- 0
    
    E_est_list <- list()
    
    for (i in 1:n) {
      E_est <- Q[[i]]-q_hat
      for (j in 1:n) {
        E_est <- E_est - rho_hat * W[i,j] * (Q[[j]]-q_hat)
      }
      E_est_list[[i]] <- E_est
      cov_est <- sum(E_est*E_est) 
      varpi_1 <- varpi_1 + cov_est
    }
    varpi_1 <- (varpi_1/n)^2
    
    
    varpi_2 <- 0
    
    for (i in 1:n) {
      for (j in setdiff(1:n, i)) {       
        ip <- sum(E_est_list[[i]] * E_est_list[[j]])  
        varpi_2 <- varpi_2 + ip^2
      }
    }
    
    # divide by n*(n-1) for unbiased estimator
    varpi_2 <- varpi_2 / (n * (n - 1))
    sigma_cst <- varpi_2 / varpi_1
    sigma_est <- n * sigma_cst * sum(diag(W %*% (W+t(W)))) / (sum(diag((W+t(W)) %*% W %*% solve(I_n - rho_hat * W))))^2
  }
  
  
  if (dim_red=='bootstrap'){
    rej_result <- (rho_hat < ci[1] || rho_hat > ci[2])
  } else{
    threshold_set <- qchisq(0.05, df=1, lower.tail=FALSE)
    wald_test <- n * rho_hat^2 / sigma_est
    rej_result <- (wald_test>=threshold_set)
  }
  
  
  
  return(c(rho_hat,rej_result))
}


inv_clr <- function(clr_row) {
  x <- exp(clr_row)
  x / sum(x)
}


PSSAR_trans_sim_prediction_function2 <- function(seed_set,sample_size, rho0=0.8){
  set.seed(seed_set)
  n <- sample_size
  d <- 4           # dimension of spherical vectors, previous setting 3
  alpha <- rho0    # neighbor correlation strength (0 < alpha < 1)
  
  mu <- rep(1,d)   # shared mean vector
  mu <- mu / sqrt(sum(mu^2))  # normalize
  
  # ----------------------
  # Generate Spatial Coordinates and Distance Matrix
  # ----------------------
  coords <- matrix(runif(n * 2), ncol = 2)
  dist_matrix <- as.matrix(dist(coords))
  
  # ----------------------
  # Define neighbor adjacency matrix W (kNN based on distance)
  # ----------------------
  k <- 20         # number of non-zero neighbors per row
  
  # Initialize W as all zeros
  W <- matrix(0, n, n)
  
  for (i in 1:n) {
    # Order distances, drop the 1st (distance to self is 0)
    dists_i <- dist_matrix[i, ]
    ordered_idx <- order(dists_i)
    neighbors <- ordered_idx[ordered_idx != i][1:k]
    
    for (j in neighbors) {
      # assign weight based on inverse distance
      w_val <- 1 #1 / dist_matrix[i, j]
      W[i, j] <- w_val
      W[j, i] <- w_val   # enforce symmetry before row-standardization
    }
  }
  
  # Row-standardize
  row_sums <- rowSums(W)
  W <- W / row_sums
  
  q_bar <- 0 
  
  # Noise terms
  epsilon <- lapply(1:n, function(i) rSkewSym(d,kappa=5*d^(1.5)))
  epsilon_mean <- Reduce("+", epsilon) / n
  

  I_n <- diag(n)
  
  # Stack epsilon as vector of skew-symmetric matrices
  E <- epsilon
  
  # Compute residuals
  Q_resid <- solve(I_n - rho0 * W)  # operator on residuals
  # Apply to each element in E
  Q <- vector("list", n)
  for (i in 1:n) {
    Qi <- matrix(0, d, d)
    for (j in 1:n) {
      Qi <- Qi + Q_resid[i,j] * E[[j]] ##(E[[j]] - epsilon_mean) ##E[[j]] ##
    }
    Q[[i]] <- q_bar + Qi
  }
  
  y <- NULL
  for (i in 1:n) {
    #angle_n <- sqrt(sum((Q[[i]] %*% mu)^2))
    #y_i <- exp_from_log_rotation_predict(Q[[i]] ,angle_n) %*% mu
    
    y_i <- expm::expm(Q[[i]]) %*% mu
    
    y <- rbind(y,t(y_i))
  }
  
  
  ## consider use 1,...,n-1 to predict n

  n_train <- n-1
  y_train <- y[1:n_train,]
  W_train <- W[1:n_train,1:n_train]
  
  mu_intrinsic <- intrinsic_mean_sphere(y_train)
  ## compute G_n
  # Compute Gram matrix G_hat
  Q <- list()
  for (i in 1:n) {
    opt_tran <- log_rotation(mu_intrinsic,y[i,])
    Q[[i]] <- opt_tran$L
  }
  
  
  
  Q_train <- Q[1:n_train]
  
  q_hat <- Reduce("+", Q_train) / n_train
  G_hat <- matrix(0, n_train, n_train)
  for (j in 1:n_train) {
    for (k in 1:n_train) {
      A <- Q[[j]] - q_hat 
      B <- Q[[k]] - q_hat 
      G_hat[j,k] <- sum(diag(t(A) %*% B))  # Frobenius inner product
    }
  }
  
  ## cost function
  objective <- function(rho, W, G_hat) {
    n <- nrow(W)
    S <- diag(n) - rho * W
    val <- sum(diag(t(S) %*% W %*% S %*% G_hat))
    return(val^2)
  }
  
  rho_hat <- optimize(objective, interval=c(-0.99, 0.99), W=W_train, G_hat=G_hat)$minimum
  
  ##prediction
  # Weighted sum of deviations
  sum_dev <- Reduce("+", lapply(1:n_train, function(j) {
    W[n, j] * (Q[[j]] - q_hat)
  }))
  
  W_row <- W[n, ]        # n-th row
  w_nn <- W[n, n]
  q_n_pred <- q_hat + rho_hat * sum_dev / (1 - rho_hat * w_nn)
  
  ## project back
  y_pred <- expm::expm(q_n_pred) %*% mu_intrinsic
  
  y_true <- y[n,]
  pred_dist <- as.numeric(acos(t(y_pred) %*% y_true))
  
  ##MSAR
  p <- ncol(y_train)
  q <- 1 ##ncol(covariate)
  ### Initial value setup
  D0 = Matrix(0, p, p)
  B0 = Matrix(0, q, p)
  Sige0 = Diagonal(n = p, x = rep(0.4,p))
  Sige0[1,2] = 0.15
  Sige0[2,1] = 0.15
  
  Ymat <- y_train
  X <- rep(1,nrow(y_train))##covariate
  
  lse_res = MSAR.Lse(Ymat, X, W_train, D = D0, B = B0, Sige = Sige0) # LSE estimation
  
  w_n_vec <- W[n,1:n_train]
  y_pred_MSAR <- as.vector(w_n_vec %*% y_train %*% lse_res$D + lse_res$B)
  y_pred_MSAR <- y_pred_MSAR / sqrt(sum(y_pred_MSAR^2))
  
  pred_dist_ignore <- as.numeric(acos(t(y_pred_MSAR) %*% y_true))
  
  return(c(pred_dist,pred_dist_ignore))
}




PSSAR_sim_prediction_int_function <- function(seed_set,sample_size,d=6,neigh_set,rho0=0.8,alpha){
  set.seed(seed_set)
  n <- sample_size
  mu <- rep(1,d)   # shared mean vector
  mu <- mu / sqrt(sum(mu^2))  # normalize
  
  # ----------------------
  # Generate Spatial Coordinates and Distance Matrix
  # ----------------------
  coords <- matrix(runif(n * 2), ncol = 2)
  dist_matrix <- as.matrix(dist(coords))
  
  # ----------------------
  # Define neighbor adjacency matrix W (kNN based on distance)
  # ----------------------
  k <- neigh_set         # number of non-zero neighbors per row
  
  # Initialize W as all zeros
  W <- matrix(0, n, n)
  
  for (i in 1:n) {
    # Order distances, drop the 1st (distance to self is 0)
    dists_i <- dist_matrix[i, ]
    ordered_idx <- order(dists_i)
    neighbors <- ordered_idx[ordered_idx != i][1:k]
    
    for (j in neighbors) {
      # assign weight based on inverse distance
      w_val <- 1 #1 / dist_matrix[i, j]
      W[i, j] <- w_val
      W[j, i] <- w_val   
    }
  }
  
  # Row-standardize
  row_sums <- rowSums(W)
  W <- W / row_sums
  
  q_bar <- 0 
  
  # Noise terms
  epsilon <- lapply(1:n, function(i) rSkewSym(d,kappa=5*d^(1.5)))
  epsilon_mean <- Reduce("+", epsilon) / n
  
  # Solve SAR-like system for q_i
  I_n <- diag(n)
  
  # Stack epsilon as vector of skew-symmetric matrices
  E <- epsilon
  
  # Compute residuals
  Q_resid <- solve(I_n - rho0 * W)  # operator on residuals
  # Apply to each element in E
  Q <- vector("list", n)
  for (i in 1:n) {
    Qi <- matrix(0, d, d)
    for (j in 1:n) {
      Qi <- Qi + Q_resid[i,j] * E[[j]] 
    }
    Q[[i]] <- q_bar + Qi
  }
  
  y <- NULL
  for (i in 1:n) {
    y_i <- expm::expm(Q[[i]]) %*% mu
    
    y <- rbind(y,t(y_i))
  }
  
  ## consider use 1,...,n-1 to predict n
  
  n_train <- n-1
  y_train <- y[1:n_train,]
  W_train <- W[1:n_train,1:n_train]
  
  mu_intrinsic <- intrinsic_mean_sphere(y_train)
  ## compute G_n
  # Compute Gram matrix G_hat
  Q <- list()
  for (i in 1:n) {
    opt_tran <- log_rotation(mu_intrinsic,y[i,])
    Q[[i]] <- opt_tran$L
  }
  
  Q_train <- Q[1:n_train]
  
  q_hat <- Reduce("+", Q_train) / n_train
  G_hat <- matrix(0, n_train, n_train)
  for (j in 1:n_train) {
    for (k in 1:n_train) {
      A <- Q[[j]] - q_hat 
      B <- Q[[k]] - q_hat 
      G_hat[j,k] <- sum(diag(t(A) %*% B))  
    }
  }
  
  ## cost function
  objective <- function(rho, W, G_hat) {
    n <- nrow(W)
    S <- diag(n) - rho * W
    val <- sum(diag(t(S) %*% W %*% S %*% G_hat))
    return(val^2)
  }
  
  rho_hat <- optimize(objective, interval=c(-0.99, 0.99), W=W_train, G_hat=G_hat)$minimum
  
  ##prediction
  # Weighted sum of deviations
  sum_dev <- Reduce("+", lapply(1:n_train, function(j) {
    W[n, j] * (Q[[j]] - q_hat)
  }))
  
  W_row <- W[n, ]        # n-th row
  w_nn <- W[n, n]
  q_n_pred <- q_hat + rho_hat * sum_dev / (1 - rho_hat * w_nn)
  
  ## project back
  y_pred <- expm::expm(q_n_pred) %*% mu_intrinsic
  
  y_true <- y[n,]
  pred_dist <- as.numeric(acos(t(y_pred) %*% y_true))
  print(sum(y_pred^2))
  
  
  ## prediction interval
  # Split indices
  idx <- sample(1:(n-1), (n-1))
  train <- idx[1:round((n-1)/2)]
  calib <- idx[(round((n-1)/2) + 1):(n-1)]
  
  # 3. Fit SAR on training
  W_train <- W[train, train]
  y_train <- y[train,]
  mu_intrinsic <- intrinsic_mean_sphere(y_train)
  
  Q <- list()
  for (k in 1:(n-1)) {
    opt_tran <- log_rotation(mu_intrinsic,y[k,])
    Q[[k]] <- opt_tran$L
  }
  
  Q_train <- Q[train]
  train_len <- length(train)
  
  q_hat <- Reduce("+", Q_train) / train_len
  G_hat <- matrix(0, train_len, train_len)
  for (j in 1:train_len) {
    for (k in 1:train_len) {
      A <- Q_train[[j]] - q_hat 
      B <- Q_train[[k]] - q_hat 
      G_hat[j,k] <- sum(diag(t(A) %*% B))  # Frobenius inner product
    }
  }
  ## cost function
  objective <- function(rho, W, G_hat) {
    S <- diag(nrow(W)) - rho * W
    val <- sum(diag(t(S) %*% W %*% S %*% G_hat))
    return(val^2)
  }
  
  rho_hat <- optimize(objective, interval=c(-0.99, 0.99), W=W_train, G_hat=G_hat)$minimum
  
  
  # 4. Compute calibration residuals strictly OUT-OF-SAMPLE using only training data
  R_vec_calib <- rep(0, length(calib))
  
  for (idx_c in 1:length(calib)) {
    i_cal <- calib[idx_c] # The actual index of the calibration node
    
    # Predict calibration node using ONLY the training set
    sum_dev_cal <- Reduce("+", lapply(train, function(j) {
      W[i_cal, j] * (Q[[j]] - q_hat)
    }))
    
    q_i_pred <- q_hat + rho_hat * sum_dev_cal
    
    # Calculate the absolute residual (Frobenius norm)
    R_vec_calib[idx_c] <- sqrt(sum((Q[[i_cal]] - q_i_pred)^2))
  }
  
  scores <- R_vec_calib
  
  calib_distances <- dist_matrix[n, calib]
  
  # Use a Gaussian kernel to compute spatial weights
  eta <- median(dist_matrix[dist_matrix > 0])
  w_raw_initial <- exp(-(calib_distances^2) / (2 * eta^2))
  
  w_raw <- w_raw_initial
  sum_w <- sum(w_raw) + 1
  w_tilde <- w_raw / sum_w
  w_tilde_n <- 1 / sum_w
  
  # Create weighted empirical distribution
  ord <- order(scores)
  sort_scores <- scores[ord]
  sort_w <- w_tilde[ord]
  
  vals <- c(sort_scores, Inf)
  cum_probs <- cumsum(c(sort_w, w_tilde_n))
  # ---------------------------------------------
  
  if (length(alpha)==1){
    # Find the weighted quantile
    q_alpha <- vals[which(cum_probs >= (1 - alpha))[1]]
    
    ##prediction
    # Weighted sum of deviations
    sum_dev <- Reduce("+", lapply(train, function(j) {
      W[n, j] * (Q[[j]] - q_hat)
    }))
    
    q_n_pred <- q_hat + rho_hat * sum_dev 
    
    R_new <- log_rotation(mu_intrinsic,y[n,])
    predict_score <- sqrt(sum((R_new$L-q_n_pred)^2))
    
    # Check coverage
    covered <- predict_score <= q_alpha
    
    return(c(pred_dist,covered,q_alpha))
  } else{
    final_results <- pred_dist
    for (PI_level in 1:length(alpha)) {
      alpha_current <- alpha[PI_level]
      
      # Find the weighted quantile
      q_alpha <- vals[which(cum_probs >= (1 - alpha_current))[1]]
      
      
      ##prediction
      # Weighted sum of deviations
      sum_dev <- Reduce("+", lapply(train, function(j) {
        W[n, j] * (Q[[j]] - q_hat)
      }))
      
      q_n_pred <- q_hat + rho_hat * sum_dev 
      
      R_new <- log_rotation(mu_intrinsic,y[n,])
      predict_score <- sqrt(sum((R_new$L-q_n_pred)^2))
      
      # Check coverage
      covered <- predict_score <= q_alpha
      
      final_results <- c(final_results,covered,q_alpha)
    }
    return(final_results)
  }
}

#### Frechet regression ####

covariate_function <- function(period,sample_size){
  x_matrix <- rep(0,period)
  x_matrix[1] <- 1
  for (i in 2:sample_size) {
    x_vec <- rep(0,period)
    current_ind <- i + period - floor((i + period -1)/period) * period
    x_vec[current_ind] <- 1
    x_matrix <- rbind(x_matrix,x_vec)
  }
  return(x_matrix)
}


l2norm <- function(x){
  #sqrt(sum(x^2))
  as.numeric(sqrt(crossprod(x)))
}

SpheGeoDist <- function(y1,y2) {
  if (abs(length(y1) - length(y2)) > 0) {
    stop("y1 and y2 should be of the same length.")
  }
  if ( !isTRUE( all.equal(l2norm(y1),1) ) ) {
    stop("y1 is not a unit vector.")
  }
  if ( !isTRUE( all.equal(l2norm(y2),1) ) ) {
    stop("y2 is not a unit vector.")
  }
  y1 = y1 / l2norm(y1)
  y2 = y2 / l2norm(y2)
  if (sum(y1 * y2) > 1){
    return(0)
  } else if (sum(y1*y2) < -1){
    return(pi)
  } else return(acos(sum(y1 * y2)))
}


SpheGeoGrad <- function(x,y) { 
  tmp <- 1 - sum(x * y) ^ 2
  return(- (tmp) ^ (-0.5) * x)
  # if (tmp < tol) {
  #   return(- Inf * x)
  # } else {
  #   return(- (tmp) ^ (-0.5) * x)
  # }
}

SpheGeoHess <- function(x,y) { #,tol = 1e-10){
  return(- sum(x * y) * (1 - sum(x * y) ^ 2) ^ (-1.5) * x %*% t(x))
}

GloSpheReg <- function(xin=NULL, yin=NULL, xout=NULL){
  
  if (is.null(xout))
    xout <- xin
  
  
  yout <- GloSpheGeoReg(xin = xin, yin = yin, xout = xout)
  
  res <- list(xout = xout, yout = yout, xin = xin, yin = yin)
  class(res) <- "spheReg"
  return(res)
}




GloSpheGeoReg <- function(xin, yin, xout) {
  
  k = nrow(xout)
  n = nrow(xin)
  m = ncol(yin)
  
  ##xbar <- colMeans(xin)
  ##Sigma <- cov(xin) * (n-1) / n
  ##invSigma <- solve(Sigma)
  
  ## period setting
  xbar <- colMeans(xin)
  invSigma <- diag(1/xbar) + rep(1,length(xbar)) %*% t(rep(1,length(xbar))) / length(xbar)
  
  yout = sapply(1:k, function(j){
    s <- 1 + t(t(xin) - xbar) %*% invSigma %*% (xout[j,] - xbar)
    s <- as.vector(s)
    
    # initial guess
    y0 = colMeans(yin*s)
    y0 = y0 / l2norm(y0)
    if ( any( sapply( 1:n, function(i) isTRUE( all.equal( sum(yin[i,]*y0), 1 ) ) ) ) ){
      # if (sum(sapply(1:n, function(i) sum(yin[i,]*y0)) > 1-1e-8)){
      #if (sum( is.infinite (sapply(1:n, function(i) (1 - sum(yin[i,]*y0)^2)^(-0.5) )[ker((xout[j] - xin) / bw)>0] ) ) + 
      #   sum(sapply(1:n, function(i) 1 - sum(yin[i,] * y0)^2 < 0)) > 0){
      y0[1] = y0[1] + 1e-3
      y0 = y0 / l2norm(y0)
    }
    
    objFctn = function(y){
      if ( ! isTRUE( all.equal(l2norm(y),1) ) ) {
        return(list(value = Inf))
      }
      f = mean(s * sapply(1:n, function(i) SpheGeoDist(yin[i,], y)^2))
      g = 2 * colMeans(t(sapply(1:n, function(i) SpheGeoDist(yin[i,], y) * SpheGeoGrad(yin[i,], y))) * s)
      res = sapply(1:n, function(i){
        grad_i = SpheGeoGrad(yin[i,], y)
        return((grad_i %*% t(grad_i) + SpheGeoDist(yin[i,], y) * SpheGeoHess(yin[i,], y)) * s[i])
      }, simplify = "array")
      h = 2 * apply(res, 1:2, mean)
      return(list(value=f, gradient=g, hessian=h))
    }
    res = trust::trust(objFctn, y0, 0.1, 1)
    return(res$argument)
  })
  return(t(yout))
}


GloSpheRegNew <- function(xin=NULL, yin=NULL, xout=NULL, add_matrix){
  
  if (is.null(xout))
    xout <- xin
  
  yout <- GloSpheGeoRegNew(xin = xin, yin = yin, xout = xout, add_matrix)
  res <- list(xout = xout, yout = yout, xin = xin, yin = yin)
  class(res) <- "spheReg"
  return(res)
}


GloSpheGeoRegNew <- function(xin, yin, xout, add_matrix) {
  
  k = nrow(xout)
  n = nrow(xin)
  m = ncol(yin)
  
  xin <- cbind(1,xin)
  xout <- cbind(1,xout)
  
  yout = sapply(1:k, function(j){
    s <- t(xout[j,]) %*% solve(t(xin) %*% t(add_matrix) %*% add_matrix %*% xin) %*% t(xin) %*% t(add_matrix) %*% add_matrix
    
    s <- as.vector(s)
    
    # initial guess
    y0 = colMeans(yin*s)
    y0 = y0 / l2norm(y0)
    if ( any( sapply( 1:n, function(i) isTRUE( all.equal( sum(yin[i,]*y0), 1 ) ) ) ) ){
      # if (sum(sapply(1:n, function(i) sum(yin[i,]*y0)) > 1-1e-8)){
      #if (sum( is.infinite (sapply(1:n, function(i) (1 - sum(yin[i,]*y0)^2)^(-0.5) )[ker((xout[j] - xin) / bw)>0] ) ) + 
      #   sum(sapply(1:n, function(i) 1 - sum(yin[i,] * y0)^2 < 0)) > 0){
      y0[1] = y0[1] + 1e-3
      y0 = y0 / l2norm(y0)
    }
    
    objFctn = function(y){
      if ( ! isTRUE( all.equal(l2norm(y),1) ) ) {
        return(list(value = Inf))
      }
      f = mean(s * sapply(1:n, function(i) SpheGeoDist(yin[i,], y)^2))
      g = 2 * colMeans(t(sapply(1:n, function(i) SpheGeoDist(yin[i,], y) * SpheGeoGrad(yin[i,], y))) * s)
      res = sapply(1:n, function(i){
        grad_i = SpheGeoGrad(yin[i,], y)
        return((grad_i %*% t(grad_i) + SpheGeoDist(yin[i,], y) * SpheGeoHess(yin[i,], y)) * s[i])
      }, simplify = "array")
      h = 2 * apply(res, 1:2, mean)
      return(list(value=f, gradient=g, hessian=h))
    }
    res = trust::trust(objFctn, y0, 0.1, 1)
    return(res$argument)
  })
  return(t(yout))
}




SSARR_tangent_space_sim_function <- function(seed_set,sample_size,dim_set,neighbor_set, rho0,dim_red='PCA'){
  set.seed(seed_set)
  # Parameters
  n <- sample_size      # number of observations
  d <- dim_set          # dimension of skew-symmetric matrices
  concern_par <- floor(0.9*d)
  # Generate row-standardized weight matrix W:: sparse
  k <- neighbor_set         # number of non-zero neighbors per row
  
  # Initialize W as all zeros
  W <- matrix(0, n, n)
  
  for (i in 1:n) {
    # choose k neighbors for node i
    neighbors <- sample(setdiff(1:n, i), k)
    
    for (j in neighbors) {
      # assign a random weight
      w <- runif(1)
      W[i, j] <- w
      W[j, i] <- w   # enforce symmetry
    }
  }
  
  # Row-standardize
  row_sums <- rowSums(W)
  W <- W / row_sums
  
  mu_set <- list()
  for (mu_ind in 1:6) {
    m2 <- rep(1,d) / sqrt(sum((rep(1,d))^2))
    ##mu_current <- movMF::rmovMF(1, 5*m2)
    mu_current <- movMF::rmovMF(1, concern_par*m2) ## perform better
    mu_set[[mu_ind]] <- mu_current
  }
  mu <- matrix(0, n, d)
  for (i in 1:n) {
    mu_chosen <- ((i-1) %% 6) + 1
    mu[i,] <- mu_set[[mu_chosen]]
  }
  
  # ----------------------------
  # Generate SAR in tangent space
  # ----------------------------
  # Step 1: generate noise in tangent space
  epsilon <- movMF::rmovMF(n, concern_par*m2)

  
  # Step 2: solve SAR for xi_i
  I_n <- diag(n)
  xi_mat <- solve(I_n - rho0 * W) %*% epsilon  # n x d
  
  y <- mu + xi_mat
  row_norms <- sqrt(rowSums(y^2))
  y <- y / row_norms
  
  
  x <- covariate_function(6,n)
  mean_est <- GloSpheReg(x,y,x)$yout
  # 1. Compute all log rotations and store in 3D array
  Q <- array(0, dim = c(d, d, n))
  Q_LIST <- list()
  for (i in 1:n) {
    qsave <- log_rotation(mean_est[i,], y[i,])$L
    Q[,,i] <- qsave
    Q_LIST[[i]] <- qsave
  }
  # 2. q_hat
  q_hat <- apply(Q, c(1,2), mean)
  
  # 3. Center matrices
  Qc <- sweep(Q, c(1,2), q_hat)
  
  # 4. Vectorize
  d2 <- d*d
  Qc_mat <- matrix(0, d2, n)
  for (i in 1:n) {
    Qc_mat[,i] <- as.vector(Qc[,,i])
  }
  
  # 5. Gram matrix (fast)
  G_hat <- crossprod(Qc_mat)
  
  ## cost function
  objective <- function(rho, W, G_hat) {
    n <- nrow(W)
    S <- diag(n) - rho * W
    val <- sum(diag(t(S) %*% W %*% S %*% G_hat))
    return(val^2)
  }
  
  rho_hat <- optimize(objective, interval=c(-0.99, 0.99), W=W, G_hat=G_hat)$minimum
  
  
  
  ##hypothesis test
  
  if (dim_red=='PCA'){
    
    # Convert back to list of d × d matrices
    ##Q <- lapply(1:n, function(i) matrix(Q_mat[i, ], d, d))
    Q <- Q_LIST
    q_hat <- Reduce("+", Q) / n
    
    E_est_list <- list()
    
    for (i in 1:n) {
      E_est <- Q[[i]]-q_hat
      for (j in 1:n) {
        E_est <- E_est - rho_hat * W[i,j] * (Q[[j]]-q_hat)
      }
      E_est_list[[i]] <- E_est
    }
    
    
    # Step 1. Vectorize each epsilon_i
    U <- t(sapply(E_est_list, function(E) as.vector(E)))  # n x d^2
    
    # Step 2. Center columns
    Uc <- scale(U, center = TRUE, scale = FALSE)
    
    # Step 3. Truncated SVD (to avoid forming huge covariance)
    if (dim_set<=10){
      Sigma_hat <- cov(Uc)  # d^2 x d^2 matrix
      eig <- eigen(Sigma_hat, symmetric = TRUE)
      eigs <- eig$values
    } else {
      ##max_m = dim_set
      max_m = d
      m <- min(max_m, n - 1)
      sv <- irlba(Uc, nv = m, nu = 0)
      eigs <- (sv$d^2) / n   # eigenvalues of sample covariance
    }
    
    
    
    # Step 4. Choose k by variance explained
    cumvar <- cumsum(eigs) / sum(eigs)
    k <- which(cumvar >= 0.9)[1]
    
    
    # Step 5. Compute varpi’s using top-k eigenvalues
    lambda <- eigs[1:k]
    
    # Step 5. Compute varpi_1, varpi_3, varpi_4
    varpi_1 <- (sum(lambda))^2
    varpi_2 <- sum(lambda^2)
    sigma_cst <- varpi_2 / varpi_1
    sigma_est <- n * sigma_cst * sum(diag(W %*% (W+t(W)))) / (sum(diag((W+t(W)) %*% W %*% solve(I_n - rho_hat * W))))^2
  } else if (dim_red=='bootstrap'){

    Q_resid <- I_n
    
    ## faster
    ## ---- Precompute E_est_list ----

    Q_centered_mat <- t(Qc_mat)
    q_hat_vec <- as.vector(q_hat)
    # Compute E_est_mat under the null
    E_est_mat <- Q_centered_mat 
    
    # Convert back to list of d × d matrices
    E_est_list <- lapply(1:n, function(i) matrix(E_est_mat[i, ], d, d))
    
    
    
    rho_est_boots <- rep(0,500)
    for (boots_ind in 1:500) {
      print(boots_ind)
      
      ## faster 
      # 2. Bootstrap residual sample
      idx <- sample.int(n, n, TRUE)
      eps_boots_mat  <- E_est_mat[idx, , drop=FALSE]
      eps_boots_mean <- colMeans(eps_boots_mat)
      
      # 3. Compute Q_i
      eps_center_mat <- sweep(eps_boots_mat, 2, eps_boots_mean)
      Q_mat <- sweep(Q_resid %*% eps_center_mat, 2, q_hat_vec, "+")
      
      # 4. Compute G_hat
      q_hat_boots_vec <- colMeans(Q_mat)
      Q_centered_mat <- sweep(Q_mat, 2, q_hat_boots_vec)
      G_hat <- Q_centered_mat %*% t(Q_centered_mat)
      
      ## cost function
      objective <- function(rho, W, G_hat) {
        n <- nrow(W)
        S <- diag(n) - rho * W
        val <- sum(diag(t(S) %*% W %*% S %*% G_hat))
        return(val^2)
      }
      
      rho_est_boots[boots_ind] <- optimize(objective, interval=c(-0.99, 0.99), W=W, G_hat=G_hat)$minimum
    }
    

    ci <- quantile(rho_est_boots, c(0.025, 0.975))
  } else{
    

    Q <- Q_LIST
    q_hat <- Reduce("+", Q) / n
    
    varpi_1 <- 0
    
    E_est_list <- list()
    
    for (i in 1:n) {
      E_est <- Q[[i]]-q_hat
      for (j in 1:n) {
        E_est <- E_est - rho_hat * W[i,j] * (Q[[j]]-q_hat)
      }
      E_est_list[[i]] <- E_est
      cov_est <- sum(E_est*E_est) 
      varpi_1 <- varpi_1 + cov_est
    }
    varpi_1 <- (varpi_1/n)^2
    
    
    varpi_2 <- 0
    
    for (i in 1:n) {
      for (j in setdiff(1:n, i)) {     
        ip <- sum(E_est_list[[i]] * E_est_list[[j]])  
        varpi_2 <- varpi_2 + ip^2
      }
    }
    
    # divide by n*(n-1) for unbiased estimator
    varpi_2 <- varpi_2 / (n * (n - 1))
    sigma_cst <- varpi_2 / varpi_1
    sigma_est <- n * sigma_cst * sum(diag(W %*% (W+t(W)))) / (sum(diag((W+t(W)) %*% W %*% solve(I_n - rho_hat * W))))^2
  }
  
  
  if (dim_red=='bootstrap'){
    rej_result <- (rho_hat < ci[1] || rho_hat > ci[2])
  } else{
    threshold_set <- qchisq(0.05, df=1, lower.tail=FALSE)
    wald_test <- n * rho_hat^2 / sigma_est
    rej_result <- (wald_test>=threshold_set)
  }
  
  
  
  return(c(rho_hat,rej_result))
}



SSARR_sim_prediction_int_function <- function(seed_set,sample_size,d=6,neigh_set,rho0=0.8,alpha){
  set.seed(seed_set)
  n <- sample_size
  
  concern_par <- floor(1.1*d)
  # ----------------------
  # Define neighbor adjacency matrix W
  # ----------------------
  k <- neigh_set         # number of non-zero neighbors per row
  
  # ----------------------
  # Generate Spatial Coordinates and Distance Matrix
  # ----------------------
  coords <- matrix(runif(n * 2), ncol = 2)
  dist_matrix <- as.matrix(dist(coords))
  
  # ----------------------
  # Define neighbor adjacency matrix W (kNN based on distance)
  # ----------------------
  k <- neigh_set         # number of non-zero neighbors per row
  
  # Initialize W as all zeros
  W <- matrix(0, n, n)
  
  for (i in 1:n) {
    dists_i <- dist_matrix[i, ]
    ordered_idx <- order(dists_i)
    neighbors <- ordered_idx[ordered_idx != i][1:k]
    
    for (j in neighbors) {
      w_val <- 1 
      W[i, j] <- w_val
      W[j, i] <- w_val   # 
    }
  }
  
  # Row-standardize
  row_sums <- rowSums(W)
  W <- W / row_sums
  
  mu_set <- list()
  for (mu_ind in 1:6) {
    m2 <- rep(1,d) / sqrt(sum((rep(1,d))^2))
    mu_current <- movMF::rmovMF(1, concern_par*m2) 
    mu_set[[mu_ind]] <- mu_current
  }
  mu <- matrix(0, n, d)
  for (i in 1:n) {
    mu_chosen <- ((i-1) %% 6) + 1
    mu[i,] <- mu_set[[mu_chosen]]
  }
  
  # ----------------------------
  # Generate SAR in tangent space
  # ----------------------------
  # Step 1: generate noise in tangent space
  epsilon <- movMF::rmovMF(n, concern_par*m2)

  # Step 2: solve SAR for xi_i
  I_n <- diag(n)
  xi_mat <- solve(I_n - rho0 * W) %*% epsilon  # n x d
  
  y <- mu + xi_mat
  row_norms <- sqrt(rowSums(y^2))
  y <- y / row_norms
  
  ## consider use 1,...,n-1 to predict n
  
  
  
  n_train <- n-1
  y_train <- y[1:n_train,]
  W_train <- W[1:n_train,1:n_train]
  
  
  x <- covariate_function(6,n)
  x_train <- x[1:n_train,]
  mean_est <- GloSpheReg(x_train,y_train,x)$yout
  ## compute G_n
  # Compute Gram matrix G_hat
  Q <- list()
  for (i in 1:n) {
    opt_tran <- log_rotation(mean_est[i,],y[i,])
    Q[[i]] <- opt_tran$L
  }
  
  
  Q_train <- Q[1:n_train]
  
  q_hat <- Reduce("+", Q_train) / n_train
  G_hat <- matrix(0, n_train, n_train)
  for (j in 1:n_train) {
    for (k in 1:n_train) {
      A <- Q[[j]] - q_hat 
      B <- Q[[k]] - q_hat 
      G_hat[j,k] <- sum(diag(t(A) %*% B))  # Frobenius inner product
    }
  }
  
  ## cost function
  objective <- function(rho, W, G_hat) {
    n <- nrow(W)
    S <- diag(n) - rho * W
    val <- sum(diag(t(S) %*% W %*% S %*% G_hat))
    return(val^2)
  }
  
  rho_hat <- optimize(objective, interval=c(-0.99, 0.99), W=W_train, G_hat=G_hat)$minimum
  
  ##prediction
  # Weighted sum of deviations
  sum_dev <- Reduce("+", lapply(1:n_train, function(j) {
    W[n, j] * (Q[[j]] - q_hat)
  }))
  
  W_row <- W[n, ]        # n-th row
  w_nn <- W[n, n]
  q_n_pred <- q_hat + rho_hat * sum_dev / (1 - rho_hat * w_nn)
  
  ## project back
  angle_n <- sqrt(sum((q_n_pred %*% mean_est[n,])^2))
  y_pred <- exp_from_log_rotation_predict(q_n_pred,angle_n) %*% mean_est[n,]
  
  
  y_true <- y[n,]
  pred_dist <- as.numeric(acos(t(y_pred) %*% y_true))
  print(sum(y_pred^2))
  ## prediction interval
  # Split indices
  idx <- sample(1:(n-1), (n-1))
  train <- idx[1:round((n-1)/2)]
  calib <- idx[(round((n-1)/2) + 1):(n-1)]
  
  # 3. Fit SAR on training
  W_train <- W[train, train]
  y_train <- y[train,]
  
  x_train <- x[train,]
  mean_est <- GloSpheReg(x_train,y_train,x)$yout
  ## compute G_n
  # Compute Gram matrix G_hat
  Q <- list()
  for (i in 1:(n-1)) {
    opt_tran <- log_rotation(mean_est[i,],y[i,])
    Q[[i]] <- opt_tran$L
  }
  
  
  Q_train <- Q[train]
  train_len <- length(train)
  
  q_hat <- Reduce("+", Q_train) / train_len
  G_hat <- matrix(0, train_len, train_len)
  for (j in 1:train_len) {
    for (k in 1:train_len) {
      A <- Q_train[[j]] - q_hat 
      B <- Q_train[[k]] - q_hat 
      G_hat[j,k] <- sum(diag(t(A) %*% B))  # Frobenius inner product
    }
  }
  ## cost function
  objective <- function(rho, W, G_hat) {
    S <- diag(nrow(W)) - rho * W
    val <- sum(diag(t(S) %*% W %*% S %*% G_hat))
    return(val^2)
  }
  
  rho_hat <- optimize(objective, interval=c(-0.99, 0.99), W=W_train, G_hat=G_hat)$minimum
  

  R_vec_calib <- rep(0, length(calib))
  
  for (idx_c in 1:length(calib)) {
    i_cal <- calib[idx_c] # The actual index of the calibration node
    
    # Predict calibration node using ONLY the training set
    sum_dev_cal <- Reduce("+", lapply(train, function(j) {
      W[i_cal, j] * (Q[[j]] - q_hat)
    }))
    
    q_i_pred <- q_hat + rho_hat * sum_dev_cal
    
    # Calculate the absolute residual (Frobenius norm)
    R_vec_calib[idx_c] <- sqrt(sum((Q[[i_cal]] - q_i_pred)^2))
  }
  
  scores <- R_vec_calib
  
  # --- NONEXCHANGEABLE SPLIT CONFORMAL ---

  calib_distances <- dist_matrix[n, calib]
  
  # Use a Gaussian kernel to compute spatial weights

  eta <- median(dist_matrix[dist_matrix > 0])
  w_raw_initial <- exp(-(calib_distances^2) / (2 * eta^2))
  
  w_raw <- w_raw_initial
  sum_w <- sum(w_raw) + 1
  w_tilde <- w_raw / sum_w
  w_tilde_n <- 1 / sum_w
  
  # Create weighted empirical distribution
  ord <- order(scores)
  sort_scores <- scores[ord]
  sort_w <- w_tilde[ord]
  
  vals <- c(sort_scores, Inf)
  cum_probs <- cumsum(c(sort_w, w_tilde_n))
  # ---------------------------------------------
  
  if (length(alpha)==1){
    # Find the weighted quantile
    q_alpha <- vals[which(cum_probs >= (1 - alpha))[1]]
    
    # 5. Pick one test node (outside train/calib)
    
    ##prediction
    # Weighted sum of deviations
    sum_dev <- Reduce("+", lapply(train, function(j) {
      W[n, j] * (Q[[j]] - q_hat)
    }))
    
    q_n_pred <- q_hat + rho_hat * sum_dev 
    
    R_new <- log_rotation(mean_est[n,],y[n,])
    predict_score <- sqrt(sum((R_new$L-q_n_pred)^2))
    
    # Check coverage
    covered <- predict_score <= q_alpha
    
    return(c(pred_dist,covered,q_alpha))
  } else{
    final_results <- pred_dist
    for (PI_level in 1:length(alpha)) {
      alpha_current <- alpha[PI_level]
      q_alpha <- vals[which(cum_probs >= (1 - alpha_current))[1]]

      ##prediction
      # Weighted sum of deviations
      sum_dev <- Reduce("+", lapply(train, function(j) {
        W[n, j] * (Q[[j]] - q_hat)
      }))
      
      q_n_pred <- q_hat + rho_hat * sum_dev 
      
      R_new <- log_rotation(mean_est[n,],y[n,])
      predict_score <- sqrt(sum((R_new$L-q_n_pred)^2))
      
      # Check coverage
      covered <- predict_score <= q_alpha
      
      final_results <- c(final_results,covered,q_alpha)
    }
    return(final_results)
  }
  
}


SSARR_trans_sim_prediction_function2 <- function(seed_set,sample_size, rho0=0.8){
  set.seed(seed_set)
  n <- sample_size
  d <- 4           # dimension of spherical vectors, previous setting 3
  alpha <- rho0    # neighbor correlation strength (0 < alpha < 1)
  
  concern_par <- floor(1.1*d)
  # ----------------------
  # Define neighbor adjacency matrix W
  # ----------------------
  k <- 20         # number of non-zero neighbors per row
  
  
  
  # Initialize W as all zeros
  W <- matrix(0, n, n)
  
  for (i in 1:n) {
    # choose k neighbors for node i
    neighbors <- sample(setdiff(1:n, i), k)
    
    for (j in neighbors) {
      # assign a random weight
      w <- runif(1)
      W[i, j] <- w
      W[j, i] <- w   # enforce symmetry
    }
  }
  
  # Row-standardize
  row_sums <- rowSums(W)
  W <- W / row_sums
  
  mu_set <- list()
  for (mu_ind in 1:6) {
    m2 <- rep(1,d) / sqrt(sum((rep(1,d))^2))
    ##mu_current <- movMF::rmovMF(1, 5*m2)
    mu_current <- movMF::rmovMF(1, concern_par*m2) ## perform better
    mu_set[[mu_ind]] <- mu_current
  }
  mu <- matrix(0, n, d)
  for (i in 1:n) {
    mu_chosen <- ((i-1) %% 6) + 1
    mu[i,] <- mu_set[[mu_chosen]]
  }
  
  # ----------------------------
  # Generate SAR in tangent space
  # ----------------------------
  # Step 1: generate noise in tangent space
  epsilon <- movMF::rmovMF(n, concern_par*m2)
  ##for (i in 1:n) {
  # Project to tangent space at mu_i
  ##  epsilon[i,] <- epsilon[i,] - sum(epsilon[i,] * mu[i,]) * mu[i,]
  ##}
  
  # Step 2: solve SAR for xi_i
  I_n <- diag(n)
  xi_mat <- solve(I_n - rho0 * W) %*% epsilon  # n x d
  
  y <- mu + xi_mat
  row_norms <- sqrt(rowSums(y^2))
  y <- y / row_norms
  
  ## consider use 1,...,n-1 to predict n
  
  
  
  n_train <- n-1
  y_train <- y[1:n_train,]
  W_train <- W[1:n_train,1:n_train]
  
  x <- covariate_function(6,n)
  x_train <- x[1:n_train,]
  mean_est <- GloSpheReg(x_train,y_train,x)$yout
  ## compute G_n
  # Compute Gram matrix G_hat
  Q <- list()
  for (i in 1:n) {
    opt_tran <- log_rotation(mean_est[i,],y[i,])
    Q[[i]] <- opt_tran$L
  }
  
  
  
  Q_train <- Q[1:n_train]
  
  q_hat <- Reduce("+", Q_train) / n_train
  G_hat <- matrix(0, n_train, n_train)
  for (j in 1:n_train) {
    for (k in 1:n_train) {
      A <- Q[[j]] - q_hat ##q_bar
      B <- Q[[k]] - q_hat ##q_bar
      G_hat[j,k] <- sum(diag(t(A) %*% B))  # Frobenius inner product
    }
  }
  
  ## cost function
  objective <- function(rho, W, G_hat) {
    n <- nrow(W)
    S <- diag(n) - rho * W
    val <- sum(diag(t(S) %*% W %*% S %*% G_hat))
    return(val^2)
  }
  
  rho_hat <- optimize(objective, interval=c(-0.99, 0.99), W=W_train, G_hat=G_hat)$minimum
  
  ##prediction
  # Weighted sum of deviations
  sum_dev <- Reduce("+", lapply(1:n_train, function(j) {
    W[n, j] * (Q[[j]] - q_hat)
  }))
  
  W_row <- W[n, ]        # n-th row
  w_nn <- W[n, n]
  q_n_pred <- q_hat + rho_hat * sum_dev / (1 - rho_hat * w_nn)
  
  ## project back
  angle_n <- sqrt(sum((q_n_pred %*% mean_est[n,])^2))
  y_pred <- exp_from_log_rotation_predict(q_n_pred,angle_n) %*% mean_est[n,]
  
  
  y_true <- y[n,]
  pred_dist <- as.numeric(acos(t(y_pred) %*% y_true))
  
  ##MSAR
  p <- ncol(y_train)
  q <- 6 ##ncol(covariate)
  ### Initial value setup
  D0 = Matrix(0, p, p)
  B0 = Matrix(0, q, p)
  Sige0 = Diagonal(n = p, x = rep(0.4,p))
  Sige0[1,2] = 0.15
  Sige0[2,1] = 0.15
  
  Ymat <- y_train
  ##X <- rep(1,nrow(y_train))##covariate
  X <- x_train
  
  x_train <- x[1:n_train,]
  
  lse_res = MSAR.Lse(Ymat, X, W_train, D = D0, B = B0, Sige = Sige0) # LSE estimation
  
  w_n_vec <- W[n,1:n_train]
  y_pred_MSAR <- as.vector(w_n_vec %*% y_train %*% lse_res$D + x[n,] %*% lse_res$B)
  y_pred_MSAR <- y_pred_MSAR / sqrt(sum(y_pred_MSAR^2))
  
  pred_dist_ignore <- as.numeric(acos(t(y_pred_MSAR) %*% y_true))
  
  return(c(pred_dist,pred_dist_ignore))
}


