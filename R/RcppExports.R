# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

rNijTRiT_ModRej <- function(tmabs, te, gm) {
    .Call(`_ctmcd_rNijTRiT_ModRej`, tmabs, te, gm)
}

rNijTRiT_Unif <- function(tmabs, te, gm, tpm) {
    .Call(`_ctmcd_rNijTRiT_Unif`, tmabs, te, gm, tpm)
}

# Register entry points for exported C++ functions
methods::setLoadAction(function(ns) {
    .Call(`_ctmcd_RcppExport_registerCCallable`)
})
