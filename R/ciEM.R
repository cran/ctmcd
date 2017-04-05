ciEM <-
function (x, alpha, eps = 1e-04, expmethod = "PadeRBS") 
{
    tmabs = x$tm
    te = x$te
    ERiT_func = function(x, Q, i0, j0, i, j, ...) {
        n = nrow(Q)
        Q[i0, j0] = x
        diag(Q) = 0
        diag(Q) = -rowSums(Q)
        dtt = expm(Q * te, method = expmethod)
        uitui = matrix(0, n, n)
        uitui[i, i] = 1
        augmat = rbind(cbind(Q, uitui), cbind(matrix(0, n, n), 
            Q))
        caseRiT = expm(te * augmat, method = expmethod)[1:n, 
            (n + 1):(2 * n)]
        ERiT = sum(tmabs * caseRiT/dtt, na.rm = TRUE)
        ERiT
    }
    ENijT_func = function(x, Q, i0, j0, i, j, ...) {
        n = nrow(Q)
        Q[i0, j0] = x
        diag(Q) = 0
        diag(Q) = -rowSums(Q)
        dtt = expm(Q * te, method = expmethod)
        uituj = matrix(0, n, n)
        uituj[i, j] = 1
        augmat = rbind(cbind(Q, uituj), cbind(matrix(0, n, n), 
            Q))
        caseNijT = Q[i, j] * expm(te * augmat, method = expmethod)[1:n, 
            (n + 1):(2 * n)]
        ENijT = sum(tmabs * caseNijT/dtt, na.rm = TRUE)
        ENijT
    }
    ParMat = x$par
    ParMat[which(x$par < eps)] = 0
    npar = sum(ParMat > 0)
    ivec = as.vector(row(ParMat) * (ParMat != 0))
    jvec = as.vector(col(ParMat) * (ParMat != 0))
    ivec = ivec[which(ivec != 0)]
    jvec = jvec[which(jvec != 0)]
    diag(ParMat) = -rowSums(ParMat)
    FI = matrix(0, length(ivec), length(jvec))
    for (k in 1:length(ivec)) {
        for (l in 1:length(ivec)) {
            if (l == k) {
                ERiTder = grad(ERiT_func, x = x$par[ivec[k], 
                  jvec[k]], Q = x$par, i0 = ivec[k], j0 = jvec[k], 
                  i = ivec[k], j = jvec[k])
                ENijTder = grad(ENijT_func, x = x$par[ivec[k], 
                  jvec[k]], Q = x$par, i0 = ivec[k], j0 = jvec[k], 
                  i = ivec[k], j = jvec[k])
                FI[k, k] = 1/x$par[ivec[k], jvec[k]]^2 * x$ENijT[ivec[k], 
                  jvec[k]] - 1/x$par[ivec[k], jvec[k]] * ENijTder + 
                  ERiTder
            }
            else {
                ERiTder = grad(ERiT_func, x = x$par[ivec[l], 
                  jvec[l]], Q = x$par, i0 = ivec[l], j0 = jvec[l], 
                  i = ivec[k], j = jvec[k])
                ENijTder = grad(ENijT_func, x = x$par[ivec[l], 
                  jvec[l]], Q = x$par, i0 = ivec[l], j0 = jvec[l], 
                  i = ivec[k], j = jvec[k])
                FI[k, l] = -1/x$par[ivec[k], jvec[k]] * ENijTder + 
                  ERiTder
            }
        }
    }
    CVmat = solve((FI + t(FI))/2)
    SEvec = sqrt(diag(CVmat))
    SEmat = matrix(0, nrow(tmabs), nrow(tmabs))
    for (k in 1:length(ivec)) {
        SEmat[ivec[k], jvec[k]] = SEvec[k]
    }
    diagse = vector(length = nrow(tmabs))
    for (i in unique(ivec)) {
        elem = jvec[which(ivec == i)]
        if (length(elem) == 1) {
            diagse[i] = SEmat[i, elem]
        }
        else {
            combs = combn(elem, 2)
            CVsum = 0
            for (k in 1:ncol(combs)) {
                par1 = intersect(which(ivec == i), which(jvec == 
                  combs[1, k]))
                par2 = intersect(which(ivec == i), which(jvec == 
                  combs[2, k]))
                CVsum = CVsum + CVmat[par1, par2]
            }
            diagse[i] = sqrt(sum(SEmat[i, elem]^2) + 2 * CVsum)
        }
    }
    diag(SEmat) = diagse
    lowermat = x$par - qnorm(1 - alpha/2) * SEmat
    lowermat[which(SEmat == 0)] = NA
    lowermat[which(diag(x$par) == 0), ] = 0
    uppermat = x$par + qnorm(1 - alpha/2) * SEmat
    uppermat[which(SEmat == 0)] = NA
    uppermat[which(diag(x$par) == 0), ] = 0
    limits = list(lower = lowermat, upper = uppermat)
    limits
}
