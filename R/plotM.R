plotM <-
function(mat,
                 mattext,
                 col = c("grey", "red"),
                 main,
                 las = 1,
                 xlab = "To",
                 ylab = "From",
                 xnames,
                 ynames,
                 cex = min(1, nrow(mat) / 8),
                 fig = 3,
                 opacity_factor) {
  mat = as.matrix(mat)
  
  if (missing(main)) {
    main = ""
  }
  
  if (missing(mattext)) {
    mattext = round(mat, fig)
  }
  
  if (missing(xnames)) {
    xnames = dimnames(mat)[[2]]
  }
  
  if (missing(ynames)) {
    ynames = dimnames(mat)[[1]]
  }
  
  nc = ncol(mat)
  nr = nrow(mat)
  
  posmat = mat
  posmat[which(posmat <= 0)] = NA
  negmat = mat
  negmat[which(mat >= 0)] = NA
  
  if (missing(opacity_factor)) {
    opacity_factor = vector(length = 2)
    if (prod(is.na(posmat)) == 0) {
      opacity_factor[1] = max(posmat[which(posmat > 0)]) / quantile(posmat[which(posmat >
                                                                                   0)], .75)[[1]]
    }
    if (prod(is.na(negmat)) == 0) {
      opacity_factor[2] = max(abs(negmat)[which(negmat < 0)]) / quantile(abs(negmat)[which(negmat <
                                                                                             0)], .75)[[1]]
    }
  }
  specp = rev(1 - ((0:(nr * nc)) / (nr * nc)) ^ opacity_factor[1])
  specn = rev(1 - ((0:(nr * nc)) / (nr * nc)) ^ opacity_factor[2])
  
  if (prod(is.na(posmat)) == 0) {
    image(
      t(apply(posmat, 2, rev)),
      col = rgb(t(col2rgb(col[1])) / 255, alpha = specp),
      main = main,
      axes = F,
      zlim = c(0, max(mat)),
      xlab = xlab,
      ylab = ylab
    )
  }
  if (prod(is.na(negmat)) == 0) {
    image(
      t(apply(abs(negmat), 2, rev)),
      col = rgb(t(col2rgb(col[2])) / 255, alpha = specn),
      main = main,
      axes = F,
      zlim = c(0, abs(min(mat))),
      xlab = xlab,
      ylab = ylab,
      add = T
    )
  }
  axis(1, (0:(nc - 1)) / (nc - 1), xnames, las = las)
  axis(2, (0:(nr - 1)) / (nr - 1), rev(ynames), las = las)
  
  rvec = (0:(nr - 1) / (nr - 1))
  cvec = (0:(nc - 1) / (nc - 1))
  for (j in 1:nr) {
    text(cvec, 1 - rvec[j], mattext[j,], cex = cex)
  }
}
