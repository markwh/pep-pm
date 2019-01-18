# Flow-law and mass-conservation functional programming


#' "params" argument is as follows:
#' - "manning_bam": n, then A0_median
#' - "metroman": a (vector), then b vector, then A0_median
#' @param swotlist a list of SWOT observables
#' @param method which flow-law to use
#' @param real_A Use real bathymetry or take this as a parameter?
#' @param gradient Return gradient with function value? 
fl <- function(swotlist, 
               method = c("bam_man", "metroman", "bam_amhg", "omniscient"),
               real_A = FALSE,
               gradient = TRUE) {
  method <- match.arg(method)
  
  if (method == "bam_man") {
    out <- fl_bamman(swotlist, real_A = real_A, gradient = gradient)
  } else if (method == "metroman") {
    out <- fl_mm(swotlist, real_A = real_A, gradient = gradient)
  } else if (method == "bam_amhg") {
    stop("AMHG flow law not implemented")
  } else if (method == "omniscient") {
    out <- fl_omni(swotlist, gradient = gradient) # Never should use gradient with omniscient
  }
  
  out
}

#' Manning BAM variant flow law
#' 
#' Returns a function of paramaters, a vector of length n_s + 1, with the first 
#' element being logn, the remaining elements as log(A_0) (median-referenced)
#' 
fl_bamman <- function(swotlist, real_A = FALSE, gradient = TRUE) {
  logW <- log(swotlist$W)
  logS <- log(swotlist$S)
  sum(is.na(logS))
  dA <- swotlist$dA
  
  ns <- nrow(dA)
  nt <- ncol(dA)

    
  out <- function(params) {
    logn <- params[1]
    
    if (real_A) {
      Amat <- swotlist$A
    } else {
      A0_med <- exp(params[-1])
      Amat <- swot_A(A0vec = A0_med, dAmat = dA, zero = "median")
    }
    
    logA <- log(Amat)
    logQ <- -2/3 * logW + 1/2 * logS + 5/3 * logA - logn
    
    outlist <- swotlist
    outlist$Ahat <- Amat # Could be misleading: this is exactly swotlist$A when real_A == TRUE
    Qhat_new <- exp(logQ)
    outlist$Qhat <- Qhat_new
    
    if (gradient) {
      # Gradient-------------------------------------------------------------
      dqdlogn <- -as.vector(Qhat_new)
      
      if (real_A) {
        gradmat <- matrix(dqdlogn, nrow = 1)
        
      } else {
        logA0 <- params[-1]
        A0mat <- swot_vec2mat(exp(logA0), Qhat_new)
        dqdlogA <- 5/3 * Qhat_new / Amat * A0mat
        
        # need a permutation matrix.
        diagmats <- replicate(diag(1, nr = ns), n = nt, 
                              simplify = FALSE)
        permmat <- Matrix(Reduce(diagmats, f = cbind))
        dlogAmat <- permmat %*% Diagonal(x = as.vector(dqdlogA))
        gradmat <- rbind(dqdlogn, dlogAmat)
      }
      
      attr(outlist, "gradient") <- gradmat
    }
    
    return(outlist)
  }
  out
}

#' MetroMan flow law
#' 
#' Implements a power-law for Manning's n based on depth (A / W)

fl_mm <- function(swotlist, real_A = FALSE, gradient = TRUE) {
  logW <- log(swotlist$W)
  logS <- log(swotlist$S)
  dA <- swotlist$dA
  ns <- nrow(logW)
  nt <- ncol(logW)
  
  out <- function(params) {
    a <- params[1:ns] # intercept in log space
    b <- params[ns + (1:ns)] # slope in log space
    bmat <- swot_vec2mat(b, logW)
    
    
    if (real_A) {
      Amat <- swotlist$A
    } else {
      logA0 <- params[-1:(-2 * ns)]
      A0_med <- exp(logA0)
      Amat <- swot_A(A0vec = A0_med, dAmat = dA, zero = "median")
    }
    
    logA <- log(Amat)
    logD <- logA - logW
    logn <- swot_vec2mat(a, logW) + bmat * logD
    
    logQ <- -2/3 * logW + 1/2 * logS + 5/3 * logA - logn
    
    outlist <- swotlist
    outlist$Ahat <- Amat
    Qhat_new <- exp(logQ)
    outlist$Qhat_fl <- Qhat_new
    outlist$Qhat <- Qhat_new
    
    if (gradient) {
      # Gradient-------------------------------------------------------------
      # need a permutation matrix.
      diagmats <- replicate(diag(1, nr = ns), n = nt, 
                            simplify = FALSE)
      permmat <- Matrix(Reduce(diagmats, f = cbind))
      
      # MetroMan parameters
      dqda <- -as.vector(Qhat_new)
      dqdb <- dqda * as.vector((logA - logW))
      
      damat <- permmat %*% Diagonal(x = dqda)
      dbmat <- permmat %*% Diagonal(x = dqdb)
      
      # A0 parameters (potentially)
      if (real_A) {
        gradmat <- rbind(damat, dbmat)
        
      } else {
        A0mat <- swot_vec2mat(A0_med, logW)
        dqdlogA <- (5 - 3 * bmat)/(3 * Amat) * Qhat_new * A0mat
        dlogAmat <- permmat %*% Diagonal(x = as.vector(dqdlogA))
        
        gradmat <- rbind(damat, dbmat, dlogAmat)
      }
      
      attr(outlist, "gradient") <- gradmat
    }

    return(outlist)
  }
  out
}

#' Omniscient "flow law"
#' 
#' Benchmark case--returns true flow, takes no parameters.

fl_omni <- function(swotlist, gradient = FALSE) {
  if (gradient) stop("nonsensical to use gradient for omniscient flow law")
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
               method = c("bam", "metroman", "median", "mean", "omniscient"),
               gradient = TRUE) {
  method <- match.arg(method)
  
  if (method == "omniscient") {
    mcflfun <- mc_omniscient(flfun)
  } else if (method == "bam") {
    mcflfun <- mc_bam(flfun, gradient = gradient)
  # } else if (method == "median") {
  #   Qhat_new <- swot_vec2mat(apply(Qhat, 2, median), Qhat)
  } else if (method == "mean") {
    mcflfun <- mc_mean(flfun, gradient = gradient)
  } else if (method == "metroman") {
    mcflfun <- mc_mm(flfun, gradient = gradient)
  }
  
  mcflfun
}

mc_omniscient <- function(flfun, gradient = TRUE) {
  
  out <- function(params) {
    swotlist <- flfun(params)
    if (!gradient) attr(swotlist, "gradient") <- NULL
    
    # Jacobian of omniscient mc function is identity matrix, so don't touch fl gradient
    return(swotlist)
  }
  
  out
}


mc_mean <- function(flfun, gradient = TRUE) {
  
  ff_mem <- memoise(flfun)
  
  out <- function(params) {
    swotlist <- ff_mem(params)
    stopifnot(!is.null(swotlist[["Qhat"]]))
    Qhat <- swotlist[["Qhat"]]
    Qhat_new <- swot_vec2mat(apply(Qhat, 2, mean), Qhat)
    swotlist$Qhat <- Qhat_new
    
    if (gradient) {
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
      
      flgrad <- attr(swotlist, "gradient")
      jacmat <- flgrad %*% mcgrad
      
      attr(swotlist, "gradient") <- jacmat
    }
    
    swotlist
  }
  
  out
}

mc_bam <- function(flfun, gradient = TRUE) {

  ff_mem <- memoise(flfun)
  
  out <- function(params) {
    # browser()
    swotlist <- ff_mem(params)
    stopifnot(!is.null(swotlist[["Qhat"]]))
    Qhat <- swotlist[["Qhat"]]
    Qhat_new <- swot_vec2mat(apply(Qhat, 2, geomMean), Qhat)
    swotlist$Qhat <- Qhat_new
    
    if (gradient) {
      ns <- nrow(Qhat)
      nt <- ncol(Qhat)
      
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
      
      flgrad <- attr(swotlist, "gradient")
      jacmat <- flgrad %*% mcgrad

      attr(swotlist, "gradient") <- jacmat
    }
    
    swotlist
  }
  
  out
}

#' Units are important here. x: meters, t: days
mc_mm <- function(flfun, gradient = TRUE) {
  
  out <- function(params) {
    swotlist <- mc_mean(flfun, gradient = gradient)(params)
    
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
                   real_A = FALSE,
                   gradient = TRUE) {
  fl_method <- match.arg(fl)
  mc_method <- match.arg(mc)
  ob <- match.arg(ob)
  
  flfun <- fl(swotlist, method = fl_method, 
              real_A = real_A, gradient = gradient)  # flow law
  mcflfun <- mc(flfun, method = mc_method, gradient = gradient) # mass-conserved flow law
  omcflfun <- metric(mcflfun, method = ob, gradient = gradient)
  
  omcflfun
}




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
    
    if (gradient) {
      #-- begin gradient calculation-----
      # Jacobian matrix of Q predictor (mcfl)
      # browser()
      pgrad <- attr(swotlist_post, "gradient")
      
      # objective gradient wrt Qbar (mcfl prediction)
      ograd1 <- as.vector(sqrt(1 / (length(relres) * relres %*% relres)))
      ograd2 <- relres / as.vector(Qobs)
      ograd <- ograd1 * ograd2
      
      # mcflo gradient (using Jacobian from mcfl)
      grad <- as.vector(pgrad %*% ograd)
      
      #--put it all together---------
      attr(rrmse, "gradient") <- grad
    }
    
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
    if (!is.finite(nrmse)) dput(params)
    
    if (gradient) {
      #-- begin gradient calculation-----
      # Jacobian matrix of Q predictor (mcfl)
      pgrad <- attr(swotlist_post, "gradient")
      
      # objective gradient wrt Qbar (mcfl prediction)
      ograd1 <- resids / mean(Qobs)
      ograd <- ograd1 / sqrt(mse) / length(resids)
      
      # mcflo gradient (using Jacobian from mcfl)
      grad <- as.vector(pgrad %*% ograd)
      
      #--put it all together---------
      attr(nrmse, "gradient") <- grad
    }

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
