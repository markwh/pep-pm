# mcflo.R
# Mark Hagemann
# 12/20/2018
# Transferring mcflo analysis from notebooks 20181129, 20181204


# Wrapper function (sloppy) -----------------------------------------------

master_mcflo <- function(case, 
                         fl = c("bam_man", "metroman", "omniscient"), 
                         mc = c("bam", "metroman", "omniscient"), 
                         metric = c("rrmse", "nrmse"), 
                         method = c("stats", "optim"),
                         ...) {
  fl <- match.arg(fl)
  mc <- match.arg(mc)
  metric <- match.arg(metric)
  method <- match.arg(method)
  swotlist <- reachdata[[case]]
  statparams <- fl_peek(swotlist, fl = fl)
  # browser()
  mcflo <- mcflob(swotlist, mc = mc, fl = fl, ob = metric, 
                  gradient = FALSE)
  res <- list()
  
  if (fl == "omniscient") {
    res$metric <- mcflo()
    res$params <- numeric()
  } else if (method == "stats") {
    res$metric <- mcflo(statparams)
    res$params <- statparams
  } else if (method == "optim") {
    minA0 <- -apply(swotlist[["dA"]], 1, min, na.rm = TRUE)
    ns <- length(minA0)
    
    ui0 <- matrix(0, nrow = ns, ncol = length(statparams) - ns)
    ui <- cbind(ui0, diag(1, ns))
    ci <- log(minA0)
    wt <- R.utils::withTimeout
    
    mcflo_gr <- mcflob(swotlist, mc = mc, fl = fl, ob = metric, 
                       gradient = TRUE)
    gradfun <- function(params) {
      out <- attr(mcflo_gr(params), "gradient")
      out
    }
    
    resopt <- wt(constrOptim(theta = statparams, f = mcflo, ui = ui, ci = ci,
                             grad = gradfun, outer.iterations = 500, ...), 
                 # timeout = 3600 * 1.5)
                 timeout = 360)

    # res$metric <- mcflo(resopt$estimate)
    # res$params <- resopt$estimate
    res$optim.info <- resopt
  }
  
  res$call <- sys.call()
  res
}


# Matrix of inputs to run through master function -------------------------

cases <- names(reachdata)
flow_laws <- c("bam_man", "metroman", "omniscient")
mc_types <- c("bam", "metroman", "omniscient")
metrics <- c("rrmse", "nrmse")
methods <- c("stats", "optim")

inputMatrix <- expand.grid(cases, flow_laws, mc_types, metrics, methods, 
                           stringsAsFactors = FALSE, 
                           KEEP.OUT.ATTRS = FALSE) %>% 
  setNames(c("case", "fl", "mc", "metric", "method"))



# Run all inputs (takes a long time!) -------------------------------------

# resultList <- pmap(.l = inputMatrix,
#                    .f = possibly(master_mcflo, otherwise = NA))

# cache("resultList")



# Attach results to input matrix ------------------------------------------

resultVec <- resultList %>% 
  map_dbl(possibly(~.$metric, otherwise = NA))
failureVec <- is.na(resultList)

resultMatrix <- within(inputMatrix, {
  result = resultVec
  failed = failureVec
  })

cache("resultMatrix")

