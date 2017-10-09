expmMC =
  function(gm, t, method="PadeRBS",order=8) {
    if (!requireNamespace("markovchain", quietly = TRUE)) {
      stop("markovchain is needed for this function to work. Please install it.",
           call. = FALSE)
    }
    if(!(class(gm) %in% c("gm","matrix","Matrix")))
      stop("Error! Expecting either a matrix or a gm object")
    if ( class(gm) %in% c("matrix","Matrix")) generator_matrix = gm else generator_matrix = as.matrix(gm[["par"]])
    transitionMatrix = expm(generator_matrix*t,method=method,order=order)
    out = list()
    attr(out,"states") = letters[1:nrow(transitionMatrix)]
    attr(out,"byrow") = TRUE
    rownames(transitionMatrix) = letters[1:nrow(transitionMatrix)]
    colnames(transitionMatrix) = letters[1:ncol(transitionMatrix)]
    attr(out,"transitionMatrix") = transitionMatrix
    attr(out,"name") = "Unnamed Markov chain"
    class(out) = "markovchain"
    attr(attr(out,"class"),"package") = "markovchain"
    out = asS4(out)
    return(out)
  }



