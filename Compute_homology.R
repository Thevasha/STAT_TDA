newdiag2landscape <- function(homology, dimension, k=0, nseq=1000){
  ## PREPROCESSING
  #  homology
  if (!inherits(homology,"homology")){
    stop("* diag2landscape : input 'homology' is not a valid homology object. Please use an output from 'diagRips' or other construction algorithms.")
  }
  #  dimension
  dimension = round(dimension) # target dimension : cannot be multiple
  if (!(dimension %in% homology$Dimension)){
    stop("* diag2landscape : input 'dimension' does not have corresponding information in the given 'homology'.")
  }
  idin      = which(homology$Dimension==dimension)
  dat.dim   = round(homology$Dimension[idin])
  dat.birth = homology$Birth[idin]
  dat.death = homology$Death[idin]
  #  others
  myk = ifelse((round(k)<1), length(dat.birth), round(k))               # number of landscape functions
  myt = seq(from=0, to=10, length.out=max(10, round(nseq))) # time sequence part  
  
  ## MAIN COMPUTATION
  #  if too many features are expected, reduce the number
  if (length(dat.birth) < myk){
    myk = length(dat.birth)
  }
  #  main compute !
  lambdas = compute_lambdas(myt, dat.birth, dat.death, myk)
  
  ## WRAP AND RETURN
  res = list(lambda=lambdas, tseq=myt, dimension=dimension)
  class(res) = "landscape"
  return(res)
}

#' @keywords internal
#' @noRd
compute_lambdas <- function(tseq, births, deaths, maxK){
  ntest = length(births)
  ntime = length(tseq)
  
  output = array(0,c(ntime,ntest))
  for (i in 1:ntest){
    b = births[i]
    d = deaths[i]
    output[,i] = base::pmax(base::pmin(tseq-b, d-tseq), rep(0,ntime))
  }
  
  newout = array(0,c(ntime,ntest))
  for (i in 1:ntime){
    newout[i,] = sort(output[i,], decreasing = TRUE)
  }
  return((newout[,1:maxK]))
}