# Functions to aid in optimization




# objective functions with gradients --------------------------------------

#' These return a function of parameters suitable to bue run in nlm().
#' 
#' @param mcflfun a "loaded" mass-conserved flow law function, ideally with a 
#'  gradient attribute

metric_rrmse <- function(mcflfun, gradient = TRUE) {
  out <- function(params) {
    if (missing(params)) {  # for omniscient flow law
      params <- numeric(0) 
      gradient <- FALSE
    }
    swotlist_post <- mcflfun(params)
    Qhat <- swotlist_post$Qhat
    Qobs <- swotlist_post[["Q"]]
    relres <- as.vector((Qhat - Qobs) / Qobs) # vector of relative resids
    
    relmse <- as.vector(relres %*% relres) / length(relres)
    rrmse <- sqrt(relmse)
    
    if (!gradient) return(rrmse)
    
    #-- begin gradient calculation-----
    # Jacobian matrix of Q predictor (mcfl)
    pgrad <- attr(mcflfun(params), "gradient")

    # objective gradient wrt Qbar (mcfl prediction)
    ograd1 <- as.vector(sqrt(1 / (length(relres) * relres %*% relres)))
    ograd2 <- relres / as.vector(Qobs)
    ograd <- ograd1 * ograd2
    
    # mcflo gradient (using Jacobian from mcfl)
    grad <- as.vector(pgrad %*% ograd)
    
    #--put it all together---------
    attr(rrmse, "gradient") <- grad
    rrmse
  }
  
  out
}


#' These return a function of parameters suitable to bue run in nlm().
#' 
#' @param mcflfun a "loaded" mass-conserved flow law function, ideally with a 
#'  gradient attribute

metric_nrmse <- function(mcflfun, gradient = TRUE) {
  out <- function(params) {
    
    if (missing(params)) {  # for omniscient flow law
      params <- numeric(0) 
      gradient <- FALSE
    }
    
    swotlist_post <- mcflfun(params)
    Qhat <- swotlist_post$Qhat
    Qobs <- swotlist_post[["Q"]]
    resids <- as.vector((Qhat - Qobs)) # vector of relative resids
    
    mse <- as.vector(resids %*% resids) / length(resids)
    nrmse <- sqrt(mse) / mean(Qobs)

    if (!gradient) return(nrmse)
    
    #-- begin gradient calculation-----
    # Jacobian matrix of Q predictor (mcfl)
    pgrad <- attr(mcflfun(params), "gradient")
    
    # objective gradient wrt Qbar (mcfl prediction)
    ograd1 <- resids / mean(Qobs)
    ograd <- ograd1 / sqrt(mse) / length(resids)
    
    # mcflo gradient (using Jacobian from mcfl)
    grad <- as.vector(pgrad %*% ograd)
    
    #--put it all together---------
    attr(nrmse, "gradient") <- grad
    nrmse
  }
  
  out
}

#' Manning n parameters for MetroMan
npars_metroman <- function(swotlist, mc = TRUE, area = c("stat", "true"),
                           qagfun = geomMean) {
  area <- match.arg(area)
  
  Wmat <- swotlist$W
  Smat <- swotlist$S
  
  if (area == "stat") {
    swotlist$dA <- rezero_dA(swotlist$dA, "median")
    A0vec <- apply(swotlist$A, 1, median)
    Amat <- swotlist$dA + swot_vec2mat(A0vec, Wmat)
  } else if (area == "true") {
    Amat <- swotlist$A
  }

  if (mc) {
    Qvec <- apply(swotlist$Q, 2, qagfun)
    Qmat <- swot_vec2mat(Qvec, Wmat)
  } else {
    Qmat <- swotlist$Q
  }
  Qdotmat <- Wmat^(-2/3) * Amat^(5/3) * Smat^(1/2)
  nmat <- Qdotmat / Qmat
  
  # log(n) = a + b * log(d) = a + b * log(A / W)
  logd <- log(Amat) - log(Wmat) # log-depth
  
  # simple linear regression of logn on logd.
  ddf <- as.data.frame(t(logd))
  ndf <- as.data.frame(t(log(nmat)))
  dncor <- map2_dbl(ddf, ndf, cor)
  dsd <- map_dbl(ddf, sd)
  nsd <- map_dbl(ndf, sd)
  dmean <- map_dbl(ddf, mean)
  nmean <- map_dbl(ndf, mean)
  bvec <- dncor * nsd / dsd
  avec <- nmean - (bvec * dmean)
  
  amat <- swot_vec2mat(avec, Wmat)
  bmat <- swot_vec2mat(bvec, Wmat)
  nmat_pred <- exp(amat + bmat * logd)
  
  out <- list(a = avec, b = bvec, n = nmat_pred)
  out
}

#' Get Manning parameters and n from stats (and possibly regression).
#' Also compute Q based on qagfun.
#' 
#' Even in the variable-n case, Q is the same everywhere, (test is for steady-state mc).
#' 
#' @importFrom swotr swot_vec2mat
#' @importFrom purrr map2_dbl map_dbl
manning_peek <- function(swotlist, man_n = c("single", "spatial", "metroman"),
                         qagfun = geomMean, nagfun = geomMean, 
                         A0fun = median) {
  man_n <- match.arg(man_n)
  swotlist$dA <- rezero_dA(swotlist$dA, "median")
  
  Qvec <- apply(swotlist$Q, 2, qagfun)
  A0vec <- apply(swotlist$A, 1, A0fun)
  
  Wmat <- swotlist$W
  Smat <- swotlist$S
  Amat <- swotlist$dA + swot_vec2mat(A0vec, Wmat)
  
  Qmat <- swot_vec2mat(Qvec, Wmat)
  Qdotmat <- Wmat^(-2/3) * Amat^(5/3) * Smat^(1/2)
  nmat <- Qdotmat / Qmat
  
  params <- list(A0 = A0vec)
  
  if (man_n == "single") {
    n <- nagfun(nmat)
    nmat_pred <- matrix(n, nrow = nrow(Wmat), ncol = ncol(Wmat))
    params$n <- n
  } else if (man_n == "spatial") {
    n <- apply(nmat, 1, nagfun)
    nmat_pred <- swot_vec2mat(n, Wmat)
    params$n <- n
  } else if (man_n == "metroman") {
    npars <- npars_metroman(swotlist, mc = TRUE, area = "stat")
    nmat_pred <- npars$n
    params$n_a <- npars$a
    params$n_b <- npars$b
  }

  Qpredmat <- Qdotmat / nmat_pred
  Qpred <- apply(Qpredmat, 2, qagfun)
  
  params$Q <- Qpred
  params
}

#' Get flow-law parameters using closure stats
#' 
#' 

fl_peek <- function(swotlist, 
                    fl = c("bam_man", "metroman", "bam_amhg", "omniscient"),
                    qagfun = geomMean,
                    nagfun = geomMean) {
  
  fl <- match.arg(fl)
  swotlist$dA <- rezero_dA(swotlist$dA, "median")
  
  Qvec <- apply(swotlist$Q, 2, qagfun)
  A0vec <- apply(swotlist$A, 1, median)
  
  Wmat <- swotlist$W
  Smat <- swotlist$S
  Amat <- swotlist$dA + swot_vec2mat(A0vec, Wmat)
  
  Qmat <- swot_vec2mat(Qvec, Wmat)
  Qdotmat <- Wmat^(-2/3) * Amat^(5/3) * Smat^(1/2)
  nmat <- Qdotmat / Qmat
  
  ns <- nrow(Wmat)
  nt <- ncol(Wmat)
  
  params <- list()
  
  if (fl == "bam_man") {
    n <- nagfun(nmat)
    nmat_pred <- matrix(n, nrow = nrow(Wmat), ncol = ncol(Wmat))
    
    pnames <- c("logn", paste0("logA0_", 1:ns))
    params <- c(log(n), log(A0vec))
    names(params) <- pnames
    
  } else if (fl == "metroman") {
    npars <- npars_metroman(swotlist, mc = TRUE, area = "stat")
    nmat_pred <- npars$n
    
    pnames <- c(paste0("a_", 1:ns), 
                paste0("b_", 1:ns),
                paste0("logA0_", 1:ns))
    params <- c(npars$a, npars$b, log(A0vec))
    names(params) <- pnames
  } else if (fl == "omniscient") {
    params <- NULL
  }
  
  params
}


# Composing a more thoughtful McFLO optimization. -------------------------

# Start with the same inputs:
# case, 
# fl = c("bam_man", "metroman", "omniscient"), 
# mc = c("bam", "metroman", "omniscient"), 
# metric = c("rrmse", "nrmse"), 
# method = c("stats", "optim") -- NOT NEEDED

# First generate a list of inputs to optimizaiton funciton, plus useful auxilary info:
# loaded function
# loaded gradient
# initial parameters
# parameter bounds
# data used for loading


# Next make a wrapper around an optimization function, with options:
# timeout
# method

#' Set up the optimization
mcflo_inps <- function(swotlist, 
                       fl = c("bam_man", "metroman", "omniscient"), 
                       mc = c("bam", "metroman", "omniscient"), 
                       metric = c("rrmse", "nrmse"),
                       msg = NA, startparams = NULL) {
  fl <- match.arg(fl)
  mc <- match.arg(mc)
  metric <- match.arg(metric)
  
  # Purge NA's from swotlist
  ndat1 <- nrow(swotlist$W) * ncol(swotlist$W)
  swotlist <- swot_purge_nas(swotlist)
  swotlist$dA <- rezero_dA(swotlist$dA)
  ndat2 <- nrow(swotlist$W) * ncol(swotlist$W)
  if (ndat2 < ndat1) {
    message(sprintf("Purged %.1f percent of data", (1 - ndat2 / ndat1) * 100))
  }
  
  # Optional message, useful for using in a loop / apply statement
  if (!is.na(msg)) {
    cat(msg, " ")
  }
  
  # Initialization
  statparams <- fl_peek(swotlist, fl = fl)
  if (is.null(startparams)) startparams <- statparams
  
  mcflo <- mcflob(swotlist, mc = mc, fl = fl, ob = metric, 
                  gradient = TRUE)
  
  # Parameter bounds (just on A0)
  # Requirement: A0 parameters are last n_s elements in vector
  minA0 <- -apply(swotlist[["dA"]], 1, min, na.rm = TRUE)
  ns <- length(minA0)
  
  minA0 <- -apply(swotlist[["dA"]], 1, min, na.rm = TRUE) + 0.1 # add a little tolerance
  ns <- length(minA0)
  mins0 <- rep(-Inf, length(statparams) - ns)
  mins <- c(mins0, log(minA0))
  
  # Compile into a list to return
  res <- list()
  res$f <- mcflo
  res$g <- function(params) attr(mcflo(params), "gradient")
  res$p_init <- startparams
  res$p_mins <- mins
  res$data <- swotlist
  
  res
}

#' Do the optimization
#' @param inplist a list of inputs, generated with \code{mcflo_inps}
mcflo_optim <- function(inplist, timeout = 600, ...) {
  
  wt <- R.utils::withTimeout
  
  # browser()
  resopt <- wt(optim(par = inplist$p_init, fn = inplist$f, 
                     lower = inplist$p_mins,
                     gr = inplist$g, method = "L-BFGS-B", ...), 
               timeout = timeout)
  
  res$metric <- mcflo(resopt$par)
  res$params <- resopt$par
  res$optim.info <- resopt
  
  res
}


#' Wrapper function to permit same operations as before
master_mcflo <- function(case, 
                         fl = c("bam_man", "metroman", "omniscient"), 
                         mc = c("bam", "metroman", "omniscient"), 
                         metric = c("rrmse", "nrmse"), 
                         method = c("stats", "optim"),
                         msg = NA, startparams = NULL,
                         timeout = 600,
                         ...) {
  fl <- match.arg(fl)
  mc <- match.arg(mc)
  metric <- match.arg(metric)
  method <- match.arg(method)
  swotlist <- reachdata[[case]]
  
  optim_inps <- mcflo_inps(swotlist = swotlist, fl = fl, mc = mc, 
                           metric = metric, msg = msg, 
                           startparams = startparams)
  
  res <- list()
  
  if (fl == "omniscient") {
    res$metric <- as.numeric(optim_inps$f())
    res$params <- numeric()
  } else if (method == "stats") {
    if (!is.null(startparams)) stop("for 'stats' method, startparams must be left NULL")
    statparams <- optim_inps$p_init
    res$metric <- as.numeric(optim_inps$f(statparams))
    res$params <- statparams
  } else if (method == "optim") {
    res <- mcflo_optim(inplist = optim_inps, timeout = timeout, ...)
  }
  
  res
}


