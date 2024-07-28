# required packages
library(e1071)
# required c-codes
Rcpp::sourceCpp('cfun-PIBP.cpp')
# load data
load("TRUTH.RData")

# calculate the mode
getmode = function(v) {
  uniquev = unique(v)
  uniquev[which.max(tabulate(match(v, uniquev)))]
}
# generate latent matrices A and B
dataAB = function(seed, P = 300, K = 5, r0 = - log(1 - 0.3), t0 = 0.2) {
  set.seed(seed)
  
  A0 = matrix(0, 48, K)
  for(k in 1 : K) {
    # specific tree structure
    n = as.numeric(runif(1) > exp(- r0 * t0))
    if(n == 1) A0[1 : 24, k] = 1 else {
      n = as.numeric(runif(1) > exp(- r0 * t0))
      if(n == 1) A0[1 : 12, k] = 1 else {
        n = as.numeric(runif(1) > exp(- r0 * t0))
        if(n == 1) A0[1 : 6, k] = 1 else {
          n = as.numeric(runif(1) > exp(- r0 * t0))
          if(n == 1) A0[1 : 3, k] = 1 else A0[1 : 3, k] = as.numeric(runif(3) > exp(- r0 * t0))
          n = as.numeric(runif(1) > exp(- r0 * t0))
          if(n == 1) A0[4 : 6, k] = 1 else A0[4 : 6, k] = as.numeric(runif(3) > exp(- r0 * t0))
        }
        n = as.numeric(runif(1) > exp(- r0 * t0))
        if(n == 1) A0[7 : 12, k] = 1 else {
          n = as.numeric(runif(1) > exp(- r0 * t0))
          if(n == 1) A0[7 : 9, k] = 1 else A0[7 : 9, k] = as.numeric(runif(3) > exp(- r0 * t0))
          n = as.numeric(runif(1) > exp(- r0 * t0))
          if(n == 1) A0[10 : 12, k] = 1 else A0[10 : 12, k] = as.numeric(runif(3) > exp(- r0 * t0))
        }
      }
      n = as.numeric(runif(1) > exp(- r0 * t0))
      if(n == 1) A0[13 : 24, k] = 1 else {
        n = as.numeric(runif(1) > exp(- r0 * t0))
        if(n == 1) A0[13 : 18, k] = 1 else {
          n = as.numeric(runif(1) > exp(- r0 * t0))
          if(n == 1) A0[13 : 15, k] = 1 else A0[13 : 15, k] = as.numeric(runif(3) > exp(- r0 * t0))
          n = as.numeric(runif(1) > exp(- r0 * t0))
          if(n == 1) A0[16 : 18, k] = 1 else A0[16 : 18, k] = as.numeric(runif(3) > exp(- r0 * t0))
        }
        n = as.numeric(runif(1) > exp(- r0 * t0))
        if(n == 1) A0[19 : 24, k] = 1 else {
          n = as.numeric(runif(1) > exp(- r0 * t0))
          if(n == 1) A0[19 : 21, k] = 1 else A0[19 : 21, k] = as.numeric(runif(3) > exp(- r0 * t0))
          n = as.numeric(runif(1) > exp(- r0 * t0))
          if(n == 1) A0[22 : 24, k] = 1 else A0[22 : 24, k] = as.numeric(runif(3) > exp(- r0 * t0))
        }
      }
    }
    n = as.numeric(runif(1) > exp(- r0 * t0))
    if(n == 1) A0[25 : 48, k] = 1 else {
      n = as.numeric(runif(1) > exp(- r0 * t0))
      if(n == 1) A0[25 : 36, k] = 1 else {
        n = as.numeric(runif(1) > exp(- r0 * t0))
        if(n == 1) A0[25 : 30, k] = 1 else {
          n = as.numeric(runif(1) > exp(- r0 * t0))
          if(n == 1) A0[25 : 27, k] = 1 else A0[25 : 27, k] = as.numeric(runif(3) > exp(- r0 * t0))
          n = as.numeric(runif(1) > exp(- r0 * t0))
          if(n == 1) A0[28 : 30, k] = 1 else A0[28 : 30, k] = as.numeric(runif(3) > exp(- r0 * t0))
        }
        n = as.numeric(runif(1) > exp(- r0 * t0))
        if(n == 1) A0[31 : 36, k] = 1 else {
          n = as.numeric(runif(1) > exp(- r0 * t0))
          if(n == 1) A0[31 : 33, k] = 1 else A0[31 : 33, k] = as.numeric(runif(3) > exp(- r0 * t0))
          n = as.numeric(runif(1) > exp(- r0 * t0))
          if(n == 1) A0[34 : 36, k] = 1 else A0[34 : 36, k] = as.numeric(runif(3) > exp(- r0 * t0))
        }
      }
      n = as.numeric(runif(1) > exp(- r0 * t0))
      if(n == 1) A0[37 : 48, k] = 1 else {
        n = as.numeric(runif(1) > exp(- r0 * t0))
        if(n == 1) A0[37 : 42, k] = 1 else {
          n = as.numeric(runif(1) > exp(- r0 * t0))
          if(n == 1) A0[37 : 39, k] = 1 else A0[37 : 39, k] = as.numeric(runif(3) > exp(- r0 * t0))
          n = as.numeric(runif(1) > exp(- r0 * t0))
          if(n == 1) A0[40 : 42, k] = 1 else A0[40 : 42, k] = as.numeric(runif(3) > exp(- r0 * t0))
        }
        n = as.numeric(runif(1) > exp(- r0 * t0))
        if(n == 1) A0[43 : 48, k] = 1 else {
          n = as.numeric(runif(1) > exp(- r0 * t0))
          if(n == 1) A0[43 : 45, k] = 1 else A0[43 : 45, k] = as.numeric(runif(3) > exp(- r0 * t0))
          n = as.numeric(runif(1) > exp(- r0 * t0))
          if(n == 1) A0[46 : 48, k] = 1 else A0[46 : 48, k] = as.numeric(runif(3) > exp(- r0 * t0))
        }
      }
    }
  }
  
  B0 = matrix(0, P, K)
  for(k in 1 : K) B0[(P / K * (k - 1) + 1) : (P / K * k), k] = 1
  # random perturbation
  B0[B0 == 0] = as.numeric(runif(sum(B0 == 0)) < 0.1)
  
  return(list(A0 = A0, B0 = B0))
}
# generate indicator Z and data X
dataZX = function(a0, b0, A0, B0, w0 = c(2.0, 2.5, 3.0, 3.5, 4.0), z0 = log(0.5), lower = 5e1, upper = 5e2) {
  N = nrow(A0); P = nrow(B0)
  
  Z0 = matrix(0, N, P)
  odds = exp(A0 %*% diag(w0) %*% t(B0) + z0)
  Z0[runif(N * P) < odds / (1 + odds)] = 1
  
  G0 = matrix(0, N, P)
  G0[Z0 == 1] = rgamma(sum(Z0 == 1), a0)
  G0[Z0 == 0] = rgamma(sum(Z0 == 0), b0)
  
  X0 = matrix(0, N, P)
  for(j in 1 : P) X0[, j] = rmultinom(1, runif(1, lower, upper), prob = G0[, j])
  
  return(X0)
}

# pIBP parents and neighbors
PARENT = c(0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14, 14, 15, 
           15, 16, 16, 16, 17, 17, 17, 18, 18, 18, 19, 19, 19, 20, 20, 20, 21, 21, 21, 22, 22, 22, 23, 23, 23,
           24, 24, 24, 25, 25, 25, 26, 26, 26, 27, 27, 27, 28, 28, 28, 29, 29, 29, 30, 30, 30, 31, 31, 31)
NEIGHBOR = vector("list", 79)
for(i in 1 : 79) {
  NEIGHBOR[[i]] = c(PARENT[i][PARENT[i] != 0], which(PARENT == i))
}
ROOT = 1 # arbitrary root

# potentials of dual cliques
potentC = function(m) {
  potA = matrix(c(exp(-m), 1 - exp(-m), 0, 1), 2)
  potB = matrix(c(exp(-m), 1 - exp(-m), 0, 0), 2)
  
  potC = array(0, dim = c(79, 79, 2, 2))
  for(i in 2 : 79) {
    if(PARENT[i] != 1) {
      potC[PARENT[i], i, , ] = potA; potC[i, PARENT[i], , ] = t(potA)
    } else {
      potC[1, i, , ] = potB; potC[i, 1, , ] = t(potB)
    }
  }
  
  return(potC)
}
sumproduct = function(evidence, nonevidence, z, potS, potC) {
  # evidence nodes
  potE = potS
  potE[, (32 : 79)[evidence]] = rbind(1 - z[evidence], z[evidence]) * potE[, (32 : 79)[evidence]]
  
  mess = array(0, dim = c(79, 79, 2)) # messages
  
  # send messages
  send = function(j, i) {
    product = c(1, 1)
    for(k in NEIGHBOR[[j]]) {
      if(k != i) product = product * mess[k, j, ]
    }
    
    mes = c(0, 0)
    for(k in 1 : 2) {
      mes = mes + potE[k, j] * product[k] * potC[j, i, , k]
    }
    
    mess[j, i, ] <<- mes
  }
  
  # collect message
  collect = function(i, j) {
    for(k in NEIGHBOR[[j]]) {
      if(k != i) collect(j, k)
    }
    
    send(j, i)
  }
  
  # distribute message
  distribute = function(i, j) {
    send(i, j)
    
    for(k in NEIGHBOR[[j]]) {
      if(k != i) distribute(j, k)
    }
  }
  
  # sum product algorithm
  for(e in NEIGHBOR[[ROOT]]) {
    collect(ROOT, e)
  }
  
  for(e in NEIGHBOR[[ROOT]]) {
    distribute(ROOT, e)
  }
  
  # compute marginal
  product = 1
  for(j in NEIGHBOR[[(32 : 79)[nonevidence]]]) {
    product = product * mess[j, (32 : 79)[nonevidence], ]
  }
    
  MARGIN = product * potS[, (32 : 79)[nonevidence]]
  MARGIN = MARGIN / sum(MARGIN)
    
  return(MARGIN)
}
marginSP = function(m, z) {
  # single potentials
  potS = matrix(rep(1, 79 * 2), 2)
  potS[, 1] = c(1, 0)
  
  # clique potentials
  potC = potentC(m)
  
  MARGIN = matrix(0, 48, 2)  # marginal distribution
  for(s in 1 : 48) MARGIN[s, ] = sumproduct((1 : 48)[-s], s, z, potS, potC)
    
  return(MARGIN[, 2])
}
jointSP = function(m, z) {
  # single potentials
  potS = matrix(rep(1, 79 * 2), 2)
  potS[, 1] = c(1, 0)
  
  # clique potentials
  potC = potentC(m)
  
  evidence = which(z == 1)[1]
  nonevidence = (1 : 48)[-evidence]
  
  JOINT = 1 # joint distribution
  for(s in nonevidence) {
    MARGIN = sumproduct(evidence, s, z, potS, potC)
    JOINT = JOINT * sum(MARGIN * rbind(1 - z[s], z[s]))
    
    evidence = c(evidence, s)
  }
  
  return(JOINT)
}
updateA = function(Z, A, B, W, MARGIN, JOINT, z, p, r, alpha, tA = 0.2, tB = 15.4, iter = 100, c = 0.06, d = 0.8, s = 0.1) {
  N = nrow(A); P = nrow(B)
  
  ANEW = A # save new matrix A
  
  for(i in 1 : N) {
    for(k in 1 : ncol(A)) {
      ANEW[i, k] = 1 - A[i, k]
      
      oddsnew = B %*% (ANEW[i, ] * W[i, ]) + z[i]
      oddsold = B %*% (A[i, ] * W[i, ]) + z[i]
      
      postnew = sum(oddsnew * Z[i, ] - log(1 + exp(oddsnew)))
      postold = sum(oddsold * Z[i, ] - log(1 + exp(oddsold)))
      
      if(A[i, k] == 1) {
        postnew = postnew + log(1 - MARGIN[i, k])
        postold = postold + log(MARGIN[i, k])
      } else {
        postnew = postnew + log(MARGIN[i, k])
        postold = postold + log(1 - MARGIN[i, k])
      }
      
      if(postold - postnew <= log(1 / runif(1) - 1)) {
        A[i, k] = ANEW[i, k]
        
        MARGIN[, k] = marginSP(- log(1 - p[k]) * tA, A[, k])
        if(sum(A[, k]) != 0) {
          JOINT[k] = jointSP(- log(1 - p[k]) * tA, A[, k])
        } else JOINT[k] = 0
      } else ANEW[i, k] = A[i, k]
    }
  
    kstar = rpois(1, alpha * (digamma(tA + tB + 1) - digamma(tB + 1)))
    
    if(kstar != 0) {
      AK = matrix(0, N, kstar); AK[i, ] = rep(1, kstar); AK = cbind(A, AK)
      BK = matrix(as.numeric(runif(P * kstar) < r), P, kstar); BK = cbind(B, BK)
      WK = matrix(rgamma(N * kstar, 1, s), N, kstar); WK = cbind(W, WK)
      
      oddsnum = BK %*% (AK[i, ] * WK[i, ]) + z[i]
      oddsden = B %*% (A[i, ] * W[i, ]) + z[i]
      
      num = sum(oddsnum * Z[i, ] - log(1 + exp(oddsnum)))
      den = sum(oddsden * Z[i, ] - log(1 + exp(oddsden)))
      
      if(log(runif(1)) <= num - den) {
        ANEW = A = AK; B = BK; W = WK
        
        pnew = JOINTNEW = rep(0, kstar); MARGINEW = matrix(0, N, kstar)
        
        for(k in 1 : kstar) {
          pg = 0.5
          
          for(t in 1 : iter) {
            pu = rnorm(1, pg, sqrt(c * pg * (1 - pg) + d))
            while(pu > 1 | pu < 0) {
              pu = rnorm(1, pg, sqrt(c * pg * (1 - pg) + d))
            }
            
            num = dnorm(pg, pu, sqrt(c * pu * (1 - pu) + d), log = TRUE) + log(1 - (1 - pu) ^ tA) + tB * log(1 - pu) - log(pu)
            den = dnorm(pu, pg, sqrt(c * pg * (1 - pg) + d), log = TRUE) + log(1 - (1 - pg) ^ tA) + tB * log(1 - pg) - log(pg)
            
            if(log(runif(1)) <= num - den) pg = pu
          }
          
          pnew[k] = pg
          
          JOINTNEW[k] = jointSP(- log(1 - pnew[k]) * tA, A[, ncol(A) - kstar + k])
          MARGINEW[, k] = marginSP(- log(1 - pnew[k]) * tA, A[, ncol(A) - kstar + k])
        }
        
        p = c(p, pnew); JOINT = c(JOINT, JOINTNEW); MARGIN = cbind(MARGIN, MARGINEW)
      }
    }
  }
  
  index = NULL
  
  for(k in 1 : ncol(A)) {
    if(sum(A[, k]) == 0) index = c(index, k)
  }
  
  if(length(index) != 0 & length(index) != ncol(A)) {
    A = as.matrix(A[, -index]); B = as.matrix(B[, -index])
    W = as.matrix(W[, -index]); MARGIN = as.matrix(MARGIN[, -index])
    p = p[-index]; JOINT = JOINT[-index]
  }
  
  if(ncol(A) == 1 & sum(A) == 0) A[which.max(MARGIN)] = 1
  
  return(list(A = A, B = B, W = W, MARGIN = MARGIN, JOINT = JOINT, p = p))
}
updatealpha = function(A, tA = 0.2, tB = 15.4) {
  a = ncol(A) + 1
  b = digamma(tA + tB + 1) - digamma(1) + 0.1
  
  return(rgamma(1, a, b))
}
updatep = function(A, MARGIN, JOINT, p, tA = 0.2, c = 0.06, d = 0.8) {
  for(k in 1 : length(p)) {
    pnew = rnorm(1, p[k], sqrt(c * p[k] * (1 - p[k]) + d))
    while(pnew > 1 | pnew < 0) {
      pnew = rnorm(1, p[k], sqrt(c * p[k] * (1 - p[k]) + d))
    }
    
    JOINTNEW = jointSP(- log(1 - pnew) * tA, A[, k])
    
    num = dnorm(p[k], pnew, sqrt(c * pnew * (1 - pnew) + d), log = TRUE) + log(JOINTNEW)
    den = dnorm(pnew, p[k], sqrt(c * p[k] * (1 - p[k]) + d), log = TRUE) + log(JOINT[k])
    
    if(log(runif(1)) <= num - den) {
      p[k] = pnew; JOINT[k] = JOINTNEW
      
      MARGIN[, k] = marginSP(- log(1 - p[k]) * tA, A[, k])
    }
  }
  
  return(list(MARGIN = MARGIN, JOINT = JOINT, p = p))
}

iterPIBP = function(X0, alpha0, beta0, maxit, burnin, interval, r = 0.5, alpha = 1, K0 = 10, tA = 0.2, tB = 15.4) {
  N = nrow(X0); P = ncol(X0)
  
  # initialize parameters
  p = rep(0.5, K0)
  A = matrix(as.numeric(runif(N * K0) < 0.5), N, K0)
  
  MARGIN = matrix(0, N, K0)
  for(k in 1 : K0) {
    MARGIN[, k] = marginSP(- log(1 - p[k]) * tA, A[, k])
  }
  
  JOINT = rep(0, K0)
  for(k in 1 : K0) {
    JOINT[k] = jointSP(- log(1 - p[k]) * tA, A[, k])
  }
  
  B = matrix(as.numeric(runif(P * K0) < r), P, K0)
  W = matrix(rgamma(N * K0, 1), N, K0)
  z = rnorm(N)
  a = b = rep(1, N)
  Z = matrix(as.numeric(runif(N * P) < 0.5), N, P)
  
  recordA = list(A)
  recordB = list(B)
  
  iter = 1
  
  while(iter < maxit) {
    Z = updateZ(X0, Z, A, B, W, z, a, b)
    a = updatea(X0, Z, a, b, alpha0, beta0)
    b = updateb(X0, Z, a, b, alpha0, beta0)
    
    alpha = updatealpha(A)
    
    RESULT = updateA(Z, A, B, W, MARGIN, JOINT, z, p, r, alpha)
    A = RESULT$A
    B = RESULT$B
    W = RESULT$W
    p = RESULT$p
    JOINT = RESULT$JOINT
    MARGIN = RESULT$MARGIN
    
    RESULT = updatep(A, MARGIN, JOINT, p)
    p = RESULT$p
    JOINT = RESULT$JOINT
    MARGIN = RESULT$MARGIN
    
    B = updateB(Z, A, B, W, z, r)
    r = updater(B)
    W = updateW(Z, A, B, W, z)
    z = updatez(Z, A, B, W, z)
    
    ## record the result
    recordA = c(recordA, list(A))
    recordB = c(recordB, list(B))
    
    iter = iter + 1
  }
  
  Khat = NULL
  for(i in seq(burnin + 1, maxit, interval)) {
    Khat = c(Khat, ncol(recordA[[i]]))
  }  # posterior mode
  
  Khat = getmode(Khat)
  
  num = 0; Ahat = Bhat = NULL
  for(i in seq(burnin + 1, maxit, interval)) {
    if(ncol(recordA[[i]]) == Khat) {
      num = num + 1
      
      Ahat = c(Ahat, list(recordA[[i]]))
      Bhat = c(Bhat, list(recordB[[i]]))
    }
  }
  
  Ahat = Ahat[[num]]
  Bhat = Bhat[[num]]
  
  return(list(Ahat = Ahat, Bhat = Bhat))
}
iterD = function(X0, maxit, burnin, interval, r = 0.5, K0 = 10, tA = 0.2, tB = 15.4, alpha = 1, s = 0.25) {
  N = nrow(X0); P = ncol(X0)
  
  # discretization
  Y0 = t(t(X0) / colSums(X0))
  
  Z0 = (Y0 > 1e-5)
  for(i in 1 : N) {
    if(sum(Z0[i, ] > 0) > (1 - s) * P) Z0[i, ] = (Y0[i, ] > quantile(Y0[i, ], s))
  }
  
  # initialize parameters
  p = rep(0.5, K0)
  A = matrix(as.numeric(runif(N * K0) < 0.5), N, K0)
  
  MARGIN = matrix(0, N, K0)
  for(k in 1 : K0) {
    MARGIN[, k] = marginSP(- log(1 - p[k]) * tA, A[, k])
  }
  
  JOINT = rep(0, K0)
  for(k in 1 : K0) {
    JOINT[k] = jointSP(- log(1 - p[k]) * tA, A[, k])
  }
  
  B = matrix(as.numeric(runif(P * K0) < r), P, K0)
  W = matrix(rgamma(N * K0, 1), N, K0)
  z = rnorm(N)
  a = b = rep(1, N)
  
  recordA = list(A)
  recordB = list(B)
  
  iter = 1
  
  while(iter < maxit) {
    alpha = updatealpha(A)
    
    RESULT = updateA(Z0, A, B, W, MARGIN, JOINT, z, p, r, alpha)
    A = RESULT$A
    B = RESULT$B
    W = RESULT$W
    p = RESULT$p
    JOINT = RESULT$JOINT
    MARGIN = RESULT$MARGIN
    
    RESULT = updatep(A, MARGIN, JOINT, p)
    p = RESULT$p
    JOINT = RESULT$JOINT
    MARGIN = RESULT$MARGIN
    
    B = updateB(Z0, A, B, W, z, r)
    r = updater(B)
    W = updateW(Z0, A, B, W, z)
    z = updatez(Z0, A, B, W, z)
    
    ## record the result
    recordA = c(recordA, list(A))
    recordB = c(recordB, list(B))
    
    iter = iter + 1
  }
  
  Khat = NULL
  for(i in seq(burnin + 1, maxit, interval)) {
    Khat = c(Khat, ncol(recordA[[i]]))
  }  # posterior mode
  
  Khat = getmode(Khat)
  
  num = 0; Ahat = Bhat = NULL
  for(i in seq(burnin + 1, maxit, interval)) {
    if(ncol(recordA[[i]]) == Khat) {
      num = num + 1
      
      Ahat = c(Ahat, list(recordA[[i]]))
      Bhat = c(Bhat, list(recordB[[i]]))
    }
  }
  
  Ahat = Ahat[[num]]
  Bhat = Bhat[[num]]
  
  return(list(Ahat = Ahat, Bhat = Bhat))
}
simulationPIBP = function(a0, b0, A0, B0, alpha0 = 1, beta0 = 0.1, maxit = 5000, burnin = 2500, interval = 5) {
  X0 = dataZX(a0, b0, A0, B0)
  
  RPIBP = iterPIBP(X0, alpha0, beta0, maxit, burnin, interval)
  RD = iterD(X0, maxit, burnin, interval)
  
  APIBP = RPIBP$Ahat
  BPIBP = RPIBP$Bhat
  
  AD = RD$Ahat
  BD = RD$Bhat
  
  return(list(APIBP, BPIBP, AD, BD))
}
evaluateR = function(A0, Ahat, B0, Bhat, K0 = 5, Khat) {
  if(ncol(A0) >= ncol(Ahat)) {
    ## list all permutations
    permutation = permutations(ncol(A0))
    permutation = permutation[, 1 : ncol(Ahat)]
    permutation = unique(permutation)
        
    DHamming = rep(0, nrow(permutation))
    for(l in 1 : nrow(permutation)) {
      Ah = A0[, permutation[l, ]]
      DHamming[l] = sum(abs(Ah - Ahat))
    }
        
    orderc = permutation[which.min(DHamming), ]
    errorA = sum(abs(Ahat - A0[, orderc])) / nrow(Ahat) / ncol(Ahat)
    errorB = sum(abs(Bhat - B0[, orderc])) / nrow(Bhat) / ncol(Bhat)
  }
    
  if(ncol(A0) < ncol(Ahat)) {
    ## list all permutations
    permutation = permutations(ncol(Ahat))
    permutation = permutation[, 1 : ncol(A0)]
    permutation = unique(permutation)
        
    DHamming = rep(0, nrow(permutation))
    for(l in 1 : nrow(permutation)) {
      Ah = Ahat[, permutation[l, ]]
      DHamming[l] = sum(abs(A0 - Ah))
    }
        
    orderc = permutation[which.min(DHamming), ]
    errorA = sum(abs(A0 - Ahat[, orderc])) / nrow(A0) / ncol(A0)
    errorB = sum(abs(B0 - Bhat[, orderc])) / nrow(B0) / ncol(B0)
  }
    
  errorK = as.numeric(Khat != K0)
    
  return(c(errorA, errorB, errorK))
}