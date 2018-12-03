# Flow-law and mass-conservation functional programming


#' "params" argument is as follows:
#' - "manning_bam": n, then A0_median
#' - "metroman": a (vector), then b vector, then A0_median
fl <- function(swotlist, 
               method = c("bam_man", "metroman", "bam_amhg", "omniscient")) {
  method <- match.arg(method)
  
  if (method == "bam_man") {
    out <- fl_bamman(swotlist)
  } else if (method == "metroman") {
    out <- fl_mm(swotlist)
  } else if (method == "bam_amhg") {
    stop("AMHG flow law not implemented")
  } else if (method == "omniscient") {
    out <- fl_omni(swotlist)
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
    Qhat_new <- exp(logQ)
    outlist$Qhat <- Qhat_new
    
    # Gradient-------------------------------------------------------------
    dqdlogn <- -as.vector(Qhat_new)
    logA0 <- params[-1]
    A0mat <- swot_vec2mat(exp(logA0), Qhat_new)
    dqdlogA <- 5/3 * Qhat_new / Amat * A0mat
    
    # need a permutation matrix.
    diagmats <- replicate(diag(1, nr = length(logA0)), n = ncol(Amat), 
                          simplify = FALSE)
    permmat <- Reduce(diagmats, f = cbind)
    dlogAmat <- permmat %*% diag(as.vector(dqdlogA))
    gradmat <- rbind(dqdlogn, dlogAmat)
    
    attr(outlist, "gradient") <- gradmat
    return(outlist)
  }
  out
}

#' MetroMan flow law
#' 
#' Implements a power-law for Manning's n based on depth (A / W)

fl_mm <- function(swotlist) {
  logW <- log(swotlist$W)
  logS <- log(swotlist$S)
  dA <- swotlist$dA
  
  
  out <- function(params) {
    ns <- nrow(logW)
    a <- params[1:ns]
    b <- params[ns + (1:ns)]
    bmat <- swot_vec2mat(b, logW)
    logA0 <- params[-1:(-2 * ns)]
    A0_med <- exp(logA0)
    Amat <- swot_A(A0vec = A0_med, dAmat = dA, zero = "median")
    logA <- log(Amat)
    logD <- logA - logW
    logn <- swot_vec2mat(a, logW) + bmat * logD
    
    logQ <- -2/3 * logW + 1/2 * logS + 5/3 * logA - logn
    
    outlist <- swotlist
    outlist$Ahat <- Amat
    Qhat_new <- exp(logQ)
    outlist$Qhat <- Qhat_new
    
    # Gradient-------------------------------------------------------------
    dqda <- -as.vector(Qhat_new)
    dqdb <- dqda * as.vector((logA - logW))
    A0mat <- swot_vec2mat(A0_med, logW)
    dqdlogA <- (5 - 3 * bmat)/(3 * Amat) * Qhat_new * A0mat
    
    # need a permutation matrix.
    diagmats <- replicate(diag(1, nr = length(logA0)), n = ncol(Amat), 
                          simplify = FALSE)
    permmat <- Reduce(diagmats, f = cbind)
    dlogAmat <- permmat %*% diag(as.vector(dqdlogA))
    damat <- permmat %*% diag(dqda)
    dbmat <- permmat %*% diag(dqdb)
    gradmat <- rbind(damat, dbmat, dlogAmat)
    
    attr(outlist, "gradient") <- gradmat
    return(outlist)
  }
  out
}

#' Omniscient "flow law"
#' 
#' Benchmark case--returns true flow, takes no parameters.

fl_omni <- function(swotlist) {
  out <- function(params = numeric(0)) {
    if (length(params > 0)) 
      warning("Nonzero-length param argument supplied.\nOmniscient flow law takes no parameters.")
    swotlist$Qhat <- swotlist$Q
    swotlist$Ahat <- swotlist$A
    
    swotlist
  }
  
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
    mcflfun <- mc_mm(flfun)
  }
  
  mcflfun
}

mc_omniscient <- function(flfun) {
  
  gradfun <- function(params) {
    if (is.null(attr(flfun(params), "gradient"))) return(NULL)
    
    swotlist <- flfun(params)
    ns <- nrow(swotlist[["W"]])
    nt <- ncol(swotlist[["W"]])
    
    flgrad <- attr(flfun(params), "gradient") # really a Jacobian matrix
    # Jacobian of omniscient mc function is identity matrix, so return flgrad
    flgrad
  }
  
  out <- function(params) {
    swotlist <- flfun(params)
    attr(swotlist, "gradient") <- gradfun(params)
    return(swotlist)
  }
  
  out
}


mc_mean <- function(flfun) {
  
  gradfun <- function(params) {
    if (is.null(attr(flfun(params), "gradient"))) return(NULL)
    
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
    
    flgrad <- attr(flfun(params), "gradient")
    
    outmat <- flgrad %*% mcgrad
    outmat
  }
  
  out <- function(params) {
    swotlist <- flfun(params)
    stopifnot(!is.null(swotlist[["Qhat"]]))
    Qhat <- swotlist[["Qhat"]]
    Qhat_new <- swot_vec2mat(apply(Qhat, 2, mean), Qhat)
    swotlist$Qhat <- Qhat_new
    attr(swotlist, "gradient") <- gradfun(params)
    swotlist
  }
  
  attr(out, "gradient") <- gradfun
  
  out
}

mc_bam <- function(flfun) {

  gradfun <- function(params) {
    if (is.null(attr(flfun(params), "gradient"))) return(NULL)
    swotlist <- flfun(params)
    Qhat <- swotlist[["Qhat"]]
    
    ns <- nrow(Qhat)
    nt <- ncol(Qhat)
    
    Qhat_new <- swot_vec2mat(apply(Qhat, 2, geomMean), Qhat)
    
    # Begin Jacobian calculation (see mcflo-math document)
    # TODO use sparse matrices from Matrix pkg 
    # http://r.789695.n4.nabble.com/sparse-matrix-from-vector-outer-product-td4701795.html
    mcgradmat <- Qhat_new / (ns * Qhat)
    mcgradvec <- as.vector(mcgradmat)
    mcgrad <- matrix(mcgradvec, nrow = ns * nt, ncol = ns * nt, byrow = FALSE)
    
    # Jacobian is zero when t != t'
    tindmat <- matrix(1:nt, nrow = ns, ncol = nt, byrow = TRUE)
    tindvec <- as.vector(tindmat)
    
    makezero <- outer(tindvec, tindvec, function(x, y) x != y)
    mcgrad[makezero] <- 0

    flgrad <- attr(flfun(params), "gradient")
    
    outmat <- flgrad %*% mcgrad
    outmat
  }
  
  out <- function(params) {
    # browser()
    swotlist <- flfun(params)
    stopifnot(!is.null(swotlist[["Qhat"]]))
    Qhat <- swotlist[["Qhat"]]
    Qhat_new <- swot_vec2mat(apply(Qhat, 2, geomMean), Qhat)
    swotlist$Qhat <- Qhat_new
    attr(swotlist, "gradient") <- gradfun(params)
    swotlist
  }
  
  out
}

#' Units are important here. x: meters, t: days
mc_mm <- function(flfun) {
  
  out <- function(params) {
    swotlist <- mc_mean(flfun)(params)
    
    Ahat <- swotlist$Ahat
    tmat <- swotlist$t * 3600 * 24 # convert from days to seconds
    dAdt <- findif_t(Ahat) / findif_t(tmat)
    xmat <- swotlist$x - mean(swotlist$x)

    Qadj <- -dAdt * xmat
    Qadj <- Qadj - mean(Qadj)
    
    swotlist$Qhat <- swotlist$Qhat + Qadj
    swotlist
  }
  out
}


# Performance metrics -----------------------------------------------------

#' Apply a performance measure to a mcfl function
#' 
#' @param mcflfun a "loaded" mass-conserved flow law function 
#' @param method Which metric to apply
#' 
metric <- function(mcflfun, method = c("rrmse", "nrmse"), gradient = TRUE) {
  method <- match.arg(method)
  
  if (method == "rrmse") {
    out <- metric_rrmse(mcflfun, gradient = gradient)
  } else {
    out <- metric_nrmse(mcflfun, gradient = gradient)
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
                   ob = c("rrmse", "nrmse"),
                   gradient = TRUE) {
  fl_method <- match.arg(fl)
  mc_method <- match.arg(mc)
  ob <- match.arg(ob)
  
  flfun <- fl(swotlist, method = fl_method)  # flow law
  mcflfun <- mc(flfun, method = mc_method) # mass-conserved flow law
  omcflfun <- metric(mcflfun, method = ob, gradient = gradient)
  
  omcflfun
}

