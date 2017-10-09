ciEMSdR=function(x,alpha,eps=1e-04,expmethod="PadeRBS"){
  JVector = function(h,i){
    J = rep(0,h)
    J[i] = 1
    J
  }
  
  JMatrix = function(h,i,j){
    J = matrix(0,h,h)
    J[i,j] = 1
    J
  }
  
  signif_level = alpha
  Q = x$par
  h = nrow(Q)
  te = x$te
  tmabs = x$tm
  
  QUseable = matrix(0,h,h)
  
  for(i in 1:h){
    for(j in 1:h){
      if(Q[i,j]>eps && j!=i)  QUseable[i,j] = 1
    }
  }
  
  VQ = matrix(0,sum(sum(QUseable)),2) 
  Points = length(VQ[,1])
  
  counter = 1
  for(i in 1:h){
    for(j in 1:h){
      if(QUseable[i,j]==1){
        VQ[counter,1] = i
        VQ[counter,2] = j
        counter = counter+1
      }
    }
  }
  
  Q = Q*(Q>eps)
  diag(Q) = -rowSums(Q)
  
  Hess = matrix(0,Points,Points)
  
  LowerPart1 = matrix(0,h,h)
  LowerPart2 = matrix(0,2*h,2*h)
  
  for(i in 1:Points){
    for(j in 1:Points){
      alpha = VQ[i,1]
      beta = VQ[i,2]
      mu = VQ[j,1]
      nu = VQ[j,2]
      
      UpperPartgamma = Q[mu,nu]*JVector(h,mu)%*%t(JVector(h,nu))
      
      #Build the matrix to be exponentiated.
      Cgamma = rbind(cbind(Q,UpperPartgamma),cbind(LowerPart1,Q))
      MatrixExponentialgamma = expm(Cgamma*te,method=expmethod)
      
      #eta
      UpperParteta = JVector(h,alpha)%*%t(JVector(h,beta))-JVector(h,alpha)%*%t(JVector(h,alpha)) #JVector is a row vector
      #Build the matrix to be exponentiated.
      Ceta = rbind(cbind(Q,UpperParteta),cbind(LowerPart1,Q))
      MatrixExponentialeta = expm(Ceta*te,method=expmethod)
      
      #phi
      UpperPartphi = JVector(h,mu)%*%t(JVector(h,mu)) #JVector is a row vector
      #Build the matrix to be exponentiated.
      Cphi = rbind(cbind(Q,UpperPartphi),cbind(LowerPart1, Q))
      MatrixExponentialphi = expm(Cphi*te,method=expmethod)
      
      #Let us now calculate the bigger matrices
      
      #psi
      #derivative of Cgamma
      DCgammaUpper = JVector(h,mu)%*%t(JVector(h,nu))*(mu==alpha && nu==beta)
      UpperPartpsi = rbind(cbind(JVector(h,alpha)%*%t(JVector(h,beta))-JVector(h,alpha)%*%t(JVector(h,alpha)),DCgammaUpper),
                         cbind(LowerPart1, JVector(h,alpha)%*%t(JVector(h,beta))-JVector(h,alpha)%*%t(JVector(h,alpha))))
      
      #Build the matrix to be exponentiated.
      Cpsi = rbind(cbind(Cgamma,UpperPartpsi),cbind(LowerPart2, Cgamma))
      MatrixExponentialpsi = expm(Cpsi*te,method=expmethod)
      
      #omega
      #derivative of Cphi
      UpperPartomega = rbind(cbind(JVector(h,alpha)%*%t(JVector(h,beta))-JVector(h,alpha)%*%t(JVector(h,alpha)),LowerPart1),
                           cbind(LowerPart1, JVector(h,alpha)%*%t(JVector(h,beta))-JVector(h,alpha)%*%t(JVector(h,alpha))))
      
      #Build the matrix to be exponentiated.
      Comega = rbind(cbind(Cphi,UpperPartomega),
                   cbind(LowerPart2, Cphi))
      MatrixExponentialomega = expm(Comega*te,method=expmethod)
      
      
      ExpQ = expm(Q*te,method=expmethod)
      
      #Now we can calculate the entries in the Hessian.
      #The formula for this is written down explicitly in the paper
      for(s in 1:(h-1)){
        for(r in 1:h){
          if(tmabs[s,r]>0){
            
            Hess[i,j] = Hess[i,j]+tmabs[s,r]*(-(1/Q[mu,nu]^2)*(ExpQ[s,r])^(-1)*(MatrixExponentialgamma[s,r+h])*(mu==alpha && nu==beta)
                                            -(1/Q[mu,nu])*(ExpQ[s,r])^(-2)*(MatrixExponentialeta[s,r+h])*(MatrixExponentialgamma[s,r+h]) 
                                            +(1/Q[mu,nu])*(ExpQ[s,r])^(-1)*(MatrixExponentialpsi[s,r+3*h])
                                            + (ExpQ[s,r])^(-2)*(MatrixExponentialeta[s,r+h])*(MatrixExponentialphi[s,r+h]) 
                                            -(ExpQ[s,r])^(-1)*(MatrixExponentialomega[s,r+3*h]))
          }
        }
      }
      
    }
  }
  
  
  #In order to find the information matrix we must take the negative of the
  #Hessian and invert it. So the information matrix is,
  
  CVmat= -solve((Hess+t(Hess))/2)
  #The estimates of the variance of q is the diagonal elements of Fisher
  #Recall for the normal distribution it is 1.96 standard deviations from the
  #mean.
  
  SEvec = sqrt(diag(CVmat))
  SEmat = matrix(0, nrow(tmabs), nrow(tmabs))
  for (k in 1:length(VQ[,1])) {
    SEmat[VQ[k,1], VQ[k,2]] = SEvec[k]
  }
  diagse = vector(length = nrow(tmabs))
  for (i in unique(VQ[,1])) {
    elem = VQ[which(VQ[,1] == i),2]
    if (length(elem) == 1) {
      diagse[i] = SEmat[i, elem]
    }
    else {
      combs = combn(elem, 2)
      CVsum = 0
      for (k in 1:ncol(combs)) {
        par1 = intersect(which(VQ[,1] == i), which(VQ[,2] == 
                                                     combs[1, k]))
        par2 = intersect(which(VQ[,1] == i), which(VQ[,2] == 
                                                     combs[2, k]))
        CVsum = CVsum + CVmat[par1, par2]
      }
      diagse[i] = sqrt(sum(SEmat[i, elem]^2) + 2 * CVsum)
    }
  }
  diag(SEmat) = diagse
  lowermat = Q - qnorm(1 - signif_level/2) * SEmat
  lowermat[which(SEmat == 0)] = NA
  lowermat[which(diag(Q) == 0), ] = 0
  uppermat = Q + qnorm(1 - signif_level/2) * SEmat
  uppermat[which(SEmat == 0)] = NA
  uppermat[which(diag(Q) == 0), ] = 0
  limits = list(lower = lowermat, upper = uppermat)
  limits
}

