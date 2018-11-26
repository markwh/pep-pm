# Flow-law and mass-conservation functional programming


#' "params" argument is as follows:
#' - "manning_bam": n, then A0_median
#' - "metroman": a (vector), then b vector, then A0_median
fl <- function(swotlist, 
               method = c("bam_man", "metroman", "bam_amhg", "omniscient")) {
  method <- match.arg(method)
  
  if (method == "bam_man") {
    out <- fl_bamman(swotlist)
  }
  
  out
}

#' Manning BAM variant flow law
#' 
#' Returns a function of paramaters, a vector of length n_s + 1, with the first 
#' element being logn, the remaining elements as log(A_0) (median-referenced)
#' 
fl_bamman <- function(swotlist) {
  logW <- log(swotlist$W)
  logS <- log(swotlist$S)
  dA <- swotlist$dA
  
  out <- function(params) {
    logn <- params[1]
    A0_med <- exp(params[-1])
    Amat <- swot_A(A0vec = A0_med, dAmat = dA, zero = "median")
    logA <- log(Amat)
    
    logQ <- -2/3 * logW + 1/2 * logS + 5/3 * logA - logn
    
    outlist <- swotlist
    outlist$Ahat <- Amat
    outlist$Qhat <- exp(logQ)
    return(outlist)
  }
  
  qgrad <- function(params) {
    Qhat <- out(params)$Qhat
    Ahat <- out(params)$Ahat
    dqdlogn <- -as.vector(Qhat)
    logA0 <- params[-1]
    dqdlogA <- 5/3 * Qhat / Ahat * swot_vec2mat(exp(logA0), Qhat)
    
    # need a permutation matrix.
    diagmats <- replicate(diag(1, nr = length(logA0)), n = ncol(Qhat), 
                          simplify = FALSE)
    permmat <- Reduce(diagmats, f = cbind)
    dlogAmat <- permmat %*% diag(as.vector(dqdlogA))
    gradmat <- rbind(dqdlogn, dlogAmat)
    gradmat
  }
  
  attr(out, "gradient") <- qgrad
  out
}


# Mass-conservation functionals -------------------------------------------


#' Apply a mass-conservation operator to a flow-law function. 
#' 
#' Returns a function of the same parameters of which fl is a function.
#' 
#' @param flfun A "loaded" flow-law function of params
#' 
#' @importFrom Matrix Diagonal Matrix
mc <- function(flfun, 
               method = c("bam", "metroman", "median", "mean", "omniscient")) {
  method <- match.arg(method)
  
  if (method == "omniscient") {
    mcflfun <- mc_omniscient(flfun)
  } else if (method == "bam") {
    mcflfun <- mc_bam(flfun)
  # } else if (method == "median") {
  #   Qhat_new <- swot_vec2mat(apply(Qhat, 2, median), Qhat)
  } else if (method == "mean") {
    mcflfun <- mc_mean(flfun)
  } else if (method == "metroman") {
    stop("metroman mc not ready yet")
  }
  
  mcflfun
}

mc_omniscient <- function(flfun) {
  
  out <- function(params) {
    swotlist <- flfun(params)
    return(swotlist)
  }
  
  if (is.null(attr(flfun, "gradient"))) {
    message("flow-law function contains no derivative info. Gradient not used.")
    return(out)
  }
  
  gradfun <- function(params) {
    swotlist <- flfun(params)
    ns <- nrow(swotlist[["W"]])
    nt <- ncol(swotlist[["W"]])
    
    flgrad <- attr(flfun, "gradient")(params) # really a Jacobian matrix
    # Jacobian of omniscient mc function is identity matrix, so return flgrad
    flgrad
  }
  
  attr(out, "gradient") <- gradfun
  out
}


mc_mean <- function(flfun) {
  
  out <- function(params) {
    swotlist <- flfun(params)
    stopifnot(!is.null(swotlist[["Qhat"]]))
    Qhat <- swotlist[["Qhat"]]
    Qhat_new <- swot_vec2mat(apply(Qhat, 2, mean), Qhat)
    swotlist$Qhat <- Qhat_new
    swotlist
  }
  
  if (is.null(attr(flfun, "gradient"))) {
    message("flow-law function contains no derivative info. Gradient not used.")
    return(out)
  }
  
  gradfun <- function(params) {
    swotlist <- flfun(params)
    Qhat <- swotlist[["Qhat"]]
    
    # Begin Jacobian calculation (see mcflo-math document)
    ns <- nrow(Qhat)
    nt <- ncol(Qhat)
    tindmat <- matrix(1:nt, nrow = ns, ncol = nt, byrow = TRUE)
    tindvec <- as.vector(tindmat)
    # TODO use sparse matrices from Matrix pkg 
    # http://r.789695.n4.nabble.com/sparse-matrix-from-vector-outer-product-td4701795.html
    makezero <- outer(tindvec, tindvec, function(x, y) x != y)
    mcgrad <- matrix(1 / ns, nrow = ns * nt, ncol = ns * nt)
    mcgrad[makezero] <- 0
    
    flgrad <- attr(flfun, "gradient")(params)
    
    outmat <- flgrad %*% mcgrad
    outmat
  }
  
  attr(out, "gradient") <- gradfun
  
  out
}

mc_bam <- function(flfun) {

  out <- function(params) {
    swotlist <- flfun(params)
    stopifnot(!is.null(swotlist[["Qhat"]]))
    Qhat <- swotlist[["Qhat"]]
    Qhat_new <- swot_vec2mat(apply(Qhat, 2, geomMean), Qhat)
    swotlist$Qhat <- Qhat_new
    swotlist
  }
  
  if (is.null(attr(flfun, "gradient"))) {
    message("flow-law function contains no derivative info. Gradient not used.")
    return(out)
  }
  
  gradfun <- function(params) {
    swotlist <- flfun(params)
    Qhat <- swotlist[["Qhat"]]
    Qhat_new <- swot_vec2mat(apply(Qhat, 2, geomMean), Qhat)
    
    # Begin Jacobian calculation (see mcflo-math document)
    coefpiece <- log(Qhat) / (nrow(Qhat) * Qhat)
    coefvec <- as.vector(coefpiece)
    
    ns <- nrow(Qhat)
    nt <- ncol(Qhat)
    # browser()
    tindmat <- matrix(1:nt, nrow = ns, ncol = nt, byrow = TRUE)
    tindvec <- as.vector(tindmat)
    # TODO use sparse matrices from Matrix pkg 
    # http://r.789695.n4.nabble.com/sparse-matrix-from-vector-outer-product-td4701795.html
    makezero <- outer(tindvec, tindvec, function(x, y) x != y)
    mcgrad <- outer(as.vector(Qhat_new), coefvec)
    mcgrad[makezero] <- 0
    
    # browser()
    
    flgrad <- attr(flfun, "gradient")(params)
    
    outmat <- flgrad %*% mcgrad
    outmat
  }
  
  attr(out, "gradient") <- gradfun
  
  out
}

#' Units are important here. x: meters, t: days
mc_mm <- function(Qhat, Amat, tmat, xmat) {
  
  dAdt <- findif_t(Amat) / findif_t(tmat)
  dx <- findif_x(xmat)
  timeconv <- 1 / (3600 * 24)
  dQ <- dAdt * dx * timeconv # convert days to seconds
  deltaQ <- apply(dQ, 2, cumsum)
  
  qpiece <- match.arg(qpiece)
  qmat <- get_qmat(swotlist, qpiece)
  
  
  
  meanQ <- swot_vec2mat(apply(Qhat, 2, mean), Qhat)
  
}


# Performance metrics -----------------------------------------------------

#' Apply a performance measure to a mcfl function
#' 
#' @param mcflfun a "loaded" mass-conserved flow law function 
#' @param method Which metric to apply
#' 
metric <- function(mcflfun, method = c("rrmse", "nrmse")) {
  method <- match.arg(method)
  
  if (method == "rrmse") {
    out <- metric_rrmse(mcflfun)
  } else {
    stop(sprintf("%s not implemented", method))
  }
  
  out
}



# Explicitly joining mc to fl, because piping isn't working ---------------

#' objective, mass-conservation, flow law wrapper
#' 
#' Compose flow-law, mass-conservation, and objective functions, and carry through
#' gradient via Jacobians
#' 
mcflob <- function(swotlist,
                 mc = c("bam", "metroman", "median", "mean", "omniscient"),
                 fl = c("bam_man", "metroman", "bam_amhg", "omniscient"),
                 ob = c("rrmse", "nrmse")) {
  fl_method <- match.arg(fl)
  mc_method <- match.arg(mc)
  ob <- match.arg(ob)
  
  flfun <- fl(swotlist, method = fl_method)  # flow law
  mcflfun <- mc(flfun, method = mc_method) # mass-conserved flow law
  omcflfun <- metric(mcflfun, method = ob)
  
  omcflfun
}

