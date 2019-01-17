# Functions to aid in optimization




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
                       method = c("stats", "optim", "part_optim"),
                       metric = c("rrmse", "nrmse"),
                       msg = NA) {
  fl <- match.arg(fl)
  mc <- match.arg(mc)
  metric <- match.arg(metric)

  # Purge NA's from swotlist
  ndat1 <- nrow(swotlist$W) * ncol(swotlist$W)
  swotlist <- swot_purge_nas(swotlist)
  swotlist$dA <- rezero_dA(swotlist$dA, zero = "median")
  ndat2 <- nrow(swotlist$W) * ncol(swotlist$W)
  if (ndat2 < ndat1) {
    message(sprintf("Purged %.1f percent of data", (1 - ndat2 / ndat1) * 100))
  }
  
  # Optional message, useful for using in a loop / apply statement
  if (!is.na(msg)) {
    cat(msg, " ")
  }

  # Initialize list to return
  res <- list()

  # Initialization
  statparams <- fl_peek(swotlist, fl = fl)
  res$p_stat <- statparams
  
  if (method == "stats" || fl == "omniscient") {
    res$f <- mcflob(swotlist, mc = mc, fl = fl, ob = metric, 
                    gradient = FALSE)
    res$g <- NULL
  } else if (method == "part_optim") {
    
    mcflo <- mcflob(swotlist, mc = mc, fl = fl, ob = metric, 
                    real_A = TRUE, gradient = TRUE)
    mcflo_mem <- memoise(mcflo)
    
    keeppars <- which(!grepl("^logA0", names(statparams)))
    res$p_stat <- res$p_stat[keeppars]
    res$f <- mcflo_mem
    res$g <- function(params) attr(mcflo_mem(params), "gradient")
    
    mins <- rep(-Inf, length(keeppars))
    
    res$p_mins <- mins    
    
  } else if (method == "optim") {
    mcflo <- mcflob(swotlist, mc = mc, fl = fl, ob = metric, 
                    real_A = FALSE, gradient = TRUE)
    mcflo_mem <- memoise(mcflo)
    res$f <- mcflo_mem
    res$g <- function(params) attr(mcflo_mem(params), "gradient")
    
    # Parameter bounds (just on A0)
    # Requirement: A0 parameters are last n_s elements in vector
    minA0 <- -apply(swotlist[["dA"]], 1, min, na.rm = TRUE)
    ns <- length(minA0)
    
    minA0 <- -apply(swotlist[["dA"]], 1, min, na.rm = TRUE) + 1 # add a little tolerance
    ns <- length(minA0)
    mins0 <- rep(-Inf, length(statparams) - ns)
    mins <- c(mins0, log(minA0) + 0.01) # More tolerance added here
    
    res$p_mins <- mins    
  }

  res$data <- swotlist
  
  res
}

#' Do the optimization
#' @param inplist a list of inputs, generated with \code{mcflo_inps}
mcflo_optim <- function(inplist, timeout = 600, startparams = NULL, ...) {
  
  wt <- R.utils::withTimeout
  
  if (is.null(startparams)) startparams <- inplist$p_stat
  
  # browser()
  resopt <- wt(optim(par = startparams, fn = inplist$f, 
                     lower = inplist$p_mins,
                     gr = inplist$g, method = "L-BFGS-B", ...), 
               timeout = timeout)
  res <- list()
  res$metric <- inplist$f(resopt$par)
  res$params <- resopt$par
  res$optim.info <- resopt
  
  if (is.memoised(inplist$f)) {
    forget(inplist$f)
  }
  
  res
}


#' Wrapper function to permit same operations as before
master_mcflo <- function(case, 
                         fl = c("bam_man", "metroman", "omniscient"), 
                         mc = c("bam", "metroman", "omniscient"), 
                         metric = c("rrmse", "nrmse"), 
                         method = c("stats", "optim", "part_optim"),
                         area = c("unknown", "known"),
                         msg = NA, startparams = NULL,
                         timeout = 600,
                         ...) {

  fl <- match.arg(fl)
  mc <- match.arg(mc)
  metric <- match.arg(metric)
  method <- match.arg(method)
  area <- match.arg(area)
  swotlist <- reachdata[[case]]
  
  optim_inps <- mcflo_inps(swotlist = swotlist, fl = fl, mc = mc, 
                           method = method, metric = metric, msg = msg)

  res <- list()
  
  if (fl == "omniscient") {
    res$metric <- as.numeric(optim_inps$f())
    res$params <- numeric()
  } else if (method == "stats") {
    if (!is.null(startparams)) stop("for 'stats' method, startparams must be left NULL")
    statparams <- optim_inps$p_stat
    res$metric <- as.numeric(optim_inps$f(statparams))
    res$params <- statparams
  } else if (method %in% c("optim", "part_optim")) {
    res <- mcflo_optim(inplist = optim_inps, timeout = timeout, 
                       startparams = startparams, ...)
  } 
  
  res
}


