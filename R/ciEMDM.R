
ciEMDM=function (x, alpha, te, eps = 1e-04, expmethod = "PadeRBS") 
{
  JVector = function(h, i) {
    J = rep(0, h)
    J[i] = 1
    J
  }
  JMatrix = function(h, i, j) {
    J = matrix(0, h, h)
    J[i, j] = 1
    J
  }
  signif_level = alpha
  Q = x$par
  h = nrow(Q)
  te0 = x$te
  tmabs = x$tm
  QUseable = matrix(0, h, h)
  for (i in 1:h) {
    for (j in 1:h) {
      if (Q[i, j] > eps && j != i) 
        QUseable[i, j] = 1
    }
  }
  VQ = matrix(0, sum(sum(QUseable)), 2)
  Points = length(VQ[, 1])
  counter = 1
  for (i in 1:h) {
    for (j in 1:h) {
      if (QUseable[i, j] == 1) {
        VQ[counter, 1] = i
        VQ[counter, 2] = j
        counter = counter + 1
      }
    }
  }
  Hess = matrix(0, Points, Points)
  Jacobian=matrix(0,h^2,Points)
  LowerPart1 = matrix(0, h, h)
  LowerPart2 = matrix(0, 2 * h, 2 * h)
  ExpQ = expm(Q * te0, method = expmethod)
  Ceta=list()
  MatrixExponentialeta=list()
  MatrixExponentialetate=list()
  for(i in 1:Points){
    alpha = VQ[i, 1]
    beta = VQ[i, 2]
    UpperParteta = JVector(h, alpha) %*% t(JVector(h,beta)) - JVector(h, alpha) %*% t(JVector(h, alpha))
    Ceta[[i]] = rbind(cbind(Q, UpperParteta), cbind(LowerPart1, Q))
    MatrixExponentialeta[[i]] = expm(Ceta[[i]] * te0, method = expmethod)
    MatrixExponentialetate[[i]] = expm(Ceta[[i]] * te, method = expmethod)
  }
  for (i in 1:Points) {
    alpha = VQ[i, 1]
    beta = VQ[i, 2]
    for (j in 1:Points) {
      mu = VQ[j, 1]
      nu = VQ[j, 2]
      UpperPartxi = rbind(cbind(JVector(h, mu) %*% t(JVector(h, nu)) - JVector(h, mu) %*% t(JVector(h,mu)),
                                LowerPart1), cbind(LowerPart1, JVector(h,mu) %*% t(JVector(h, nu)) - JVector(h, mu) %*% t(JVector(h, mu))))
      Cxi = rbind(cbind(Ceta[[i]], UpperPartxi), cbind(LowerPart2, Ceta[[i]]))
      MatrixExponentialxi = expm(Cxi * te0, method = expmethod)
      for (s in 1:h) {
        for (r in 1:h) {
          Jacobian[(s-1)*h+r,i]=MatrixExponentialetate[[i]][s,h+r]
          if (tmabs[s, r] > 0) {
            Hess[i, j] = Hess[i, j] + tmabs[s, r] * ExpQ[s, r]^(-1) *( MatrixExponentialxi[s,r+3*h]
                                                                       -ExpQ[s, r]^(-1) * MatrixExponentialeta[[i]][s, h+ r]*MatrixExponentialeta[[j]][s, h+ r])
          }
        }
      }
    }
  }
  CVmat=solve(-Hess)
  SEmat=matrix(0,h,h)
  for(s0 in 1:h){
    for(sT in 1:h){
      SEmat[s0,sT]=t(Jacobian[(s0-1)*h+sT,])%*%CVmat%*%t(t(Jacobian[(s0-1)*h+sT,]))
    }
  }
  SEmat[which(diag(Q)==0),]=0
  SEmat=sqrt(SEmat)
  P=expm(Q*te)
  lowermat = P - qnorm(1 - signif_level/2) * SEmat
  lowermat[which(SEmat == 0)] = NA
  lowermat[which(diag(Q) == 0), ] = 0
  uppermat = P + qnorm(1 - signif_level/2) * SEmat
  uppermat[which(SEmat == 0)] = NA
  uppermat[which(diag(Q) == 0), ] = 0
  limits = list(lower = lowermat, upper = uppermat, SE=SEmat)
  limits
}