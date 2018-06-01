ciEMoFI=function (x, alpha, eps = 1e-04, expmethod = "PadeRBS") 
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
  te = x$te
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
  LowerPart1 = matrix(0, h, h)
  LowerPart2 = matrix(0, 2 * h, 2 * h)
  ExpQ = expm(Q * te, method = expmethod)
  Ceta=list()
  MatrixExponentialeta=list()
  for(i in 1:Points){
    alpha = VQ[i, 1]
    beta = VQ[i, 2]
    UpperParteta = JVector(h, alpha) %*% t(JVector(h,beta)) - JVector(h, alpha) %*% t(JVector(h, alpha))
    Ceta[[i]] = rbind(cbind(Q, UpperParteta), cbind(LowerPart1, Q))
    MatrixExponentialeta[[i]] = expm(Ceta[[i]] * te, method = expmethod)
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
      MatrixExponentialxi = expm(Cxi * te, method = expmethod)
      for (s in 1:h) {
        for (r in 1:h) {
          if (tmabs[s, r] > 0) {
            Hess[i, j] = Hess[i, j] + tmabs[s, r] * ExpQ[s, r]^(-1) *( MatrixExponentialxi[s,r+3*h]
                                                                       -ExpQ[s, r]^(-1) * MatrixExponentialeta[[i]][s, h+ r])*MatrixExponentialeta[[j]][s,h+r]
          }
        }
      }
    }
  }
  CVmat = -solve((Hess + t(Hess))/2)
  SEvec = sqrt(diag(CVmat))
  SEmat = matrix(0, nrow(tmabs), nrow(tmabs))
  for (k in 1:length(VQ[, 1])) {
    SEmat[VQ[k, 1], VQ[k, 2]] = SEvec[k]
  }
  diagse = vector(length = nrow(tmabs))
  for (i in unique(VQ[, 1])) {
    elem = VQ[which(VQ[, 1] == i), 2]
    if (length(elem) == 1) {
      diagse[i] = SEmat[i, elem]
    }
    else {
      combs = combn(elem, 2)
      CVsum = 0
      for (k in 1:ncol(combs)) {
        par1 = intersect(which(VQ[, 1] == i), which(VQ[,2] == combs[1, k]))
        par2 = intersect(which(VQ[, 1] == i), which(VQ[,2] == combs[2, k]))
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
  limits = list(lower = lowermat, upper = uppermat, FI=Hess, CVmat=CVmat)
  limits
}
