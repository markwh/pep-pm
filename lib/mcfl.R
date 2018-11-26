# Flow-law and mass-conservation functional programming



#' smartly get q matrix (either Qhat or Q) from a swotlist
get_qmat <- function(swotlist, qpiece = c("choose", "Q", "Qhat")) {
  qpiece <- match.arg(qpiece)
  if (qpiece == "choose") {
    qpiece <- ifelse(is.null(swotlist[["Qhat"]]), "Q", "Qhat")
  }
  qmat <- swotlist[[qpiece]]
  if (is.null(qmat)) stop(sprintf("%s is missing", qpicece))  
  qmat
}

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
#' @param fl A flow-law function of params
mc <- function(oper, 
               method = c("bam", "metroman", "median", "mean", "omniscient")) {
  method <- match.arg(method)
  
  out <- function(params) {
    if (is.list(oper) && names(oper) == c("value", "visible")) {
      oper <- oper[["value"]]
    }
    stopifnot(is.function(oper))
    outlist <- oper(params)
    Qhat <- outlist$Qhat
    
    if (method == "omniscient") {
      Qhat_new <- Qhat
    } else if (method == "bam") {
      Qhat_new <- swot_vec2mat(apply(Qhat, 2, geomMean), Qhat)
    } else if (method == "median") {
      Qhat_new <- swot_vec2mat(apply(Qhat, 2, median), Qhat)
    } else if (method == "mean") {
      Qhat_new <- swot_vec2mat(apply(Qhat, 2, mean), Qhat)
    } else if (method == "metroman") {
      stop("metroman mc not ready yet")
    }
    outlist$Qhat <- Qhat_new
    outlist
  }
  
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

metric <- function(oper, method = c("rrmse", "nrmse")) {
  method <- match.arg(method)
  
  evalfun <- if (method == "rrmse") RRMSE else 
    if (method == "nrmse") NRMSE else
      stop("metric not implemented")
  
  out <- function(params) {
    outlist <- oper(params)
    
    out2 <- evalfun(pred = as.vector(outlist[["Qhat"]]), 
                  meas = as.vector(outlist[["Q"]]))
    out2
  }
}



# Explicitly joining mc to fl, because piping isn't working ---------------

mcflob <- function(swotlist,
                 mc = c("bam", "metroman", "median", "mean", "omniscient"),
                 fl = c("bam_man", "metroman", "bam_amhg", "omniscient"),
                 ob = c("rrmse", "nrmse")) {
  fl_method <- match.arg(fl)
  mc_method <- match.arg(mc)
  ob <- match.arg(ob)
  out1 <- fl(swotlist, method = fl_method)
  out2 <- mc(out1, method = mc_method)
  out3 <- metric(out2, method = ob)
  out3
}

