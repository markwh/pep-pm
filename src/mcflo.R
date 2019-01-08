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
                         msg = NA, startparams = NULL,
                         timeout = 3600,
                         ...) {
  fl <- match.arg(fl)
  mc <- match.arg(mc)
  metric <- match.arg(metric)
  method <- match.arg(method)
  swotlist <- reachdata[[case]]
  
  ndat1 <- nrow(swotlist$W) * ncol(swotlist$W)
  swotlist <- swot_purge_nas(swotlist)
  swotlist$dA <- rezero_dA(swotlist$dA)
  ndat2 <- nrow(swotlist$W) * ncol(swotlist$W)
  if (ndat2 < ndat1) {
    message(sprintf("Purged %.1f percent of data", (1 - ndat2 / ndat1) * 100))
  }
  
  if (!is.na(msg)) {
    cat(msg, " ")
  }
  
  statparams <- fl_peek(swotlist, fl = fl)
  if (is.null(startparams)) startparams <- statparams
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
    
    minA0 <- -apply(swotlist[["dA"]], 1, min, na.rm = TRUE) + 0.1 # add a little tolerance
    ns <- length(minA0)
    mins0 <- rep(-Inf, length(statparams) - ns)
    mins <- c(mins0, log(minA0))
    
    wt <- R.utils::withTimeout
    
    mcflo_gr <- mcflob(swotlist, mc = mc, fl = fl, ob = metric, 
                       gradient = TRUE)
    gradfun <- function(params) {
      out <- attr(mcflo_gr(params), "gradient")
      out
    }
    print(mcflo(startparams))
    # browser()
    resopt <- wt(optim(par = startparams, fn = mcflo, lower = mins,
                             gr = gradfun, method = "L-BFGS-B", ...), 
                 # timeout = 3600 * 1.5)
                 timeout = timeout)

    res$metric <- mcflo(resopt$par)
    res$params <- resopt$par
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
  setNames(c("case", "fl", "mc", "metric", "method")) %>% 
  mutate(msg = 1:nrow(.))



# Run all inputs (takes a long time!) -------------------------------------

resultList <- pmap(.l = inputMatrix, #[-1:-342, ],
                   # .f = possibly(master_mcflo, otherwise = NA))
                   .f = possibly(mmcflo2, otherwise = NA))
                   # control = list(trace = NULL))
cache("resultList")

# # Go back and fix my stupid mistakes.
# fillin_missing <- function(prt, inputrow) {
#   if (is.na(prt)) return(prt)
#   
#   print(inputrow[["rowid"]])
#   
#   if (is.null(prt[["metric"]]) && ! is.null(prt[["optim.info"]])) {
#     case <- inputrow[["case"]]
#     flfun <- inputrow[["fl"]]
#     mcfun <- inputrow[["mc"]]
#     metricfun <- inputrow[["metric"]]
#     
#     swotlist <- reachdata[[case]]
#     mcflo <- mcflob(swotlist, mc = mcfun, fl = flfun, ob = metricfun, 
#                     gradient = FALSE)
#     
#     prt$params <- prt[["optim.info"]][["par"]]
#     prt$metric <- mcflo(prt[["params"]])
#   }
#   prt
# }
# 
# reslist2 <- inputMatrix %>% 
#   mutate(rowid = 1:nrow(.)) %>% 
#   split(f = 1:nrow(.)) %>% 
#   map2(resultList, ~fillin_missing(prt = .y, inputrow = .x))


# manually redo bad optimizations

# toredo <- which(convgvec == 1)
toredo <- 602

startparams <- c(-4.28992, -4.25515, -0.275068, -2.1414, -4.65791, -2.63876, 
                 -4.40692, 0.444417, -0.216653, -4.21376, -1.42375, 0.779172, 
                 -1.42048, 0.464116, 4.14221, 4.98235, 4.0535, 5.30748, 4.4135, 
                 4.42257, 6.1248)

# foo_sl <- within(swot_purge_nas(reachdata$Cumberland), {dA = rezero_dA(dA, "median")})

inputMatrix[toredo, ]

inps1 <- mcflo_inps(swotlist = reachdata$SacramentoUpstream, fl = "metroman", 
           mc = "metroman", metric = "nrmse")
opt1 <- mcflo_optim(inps1, control = list(trace = 6, maxit = 500))



mmfargs <- append(as.list(inputMatrix[toredo, ]), 
                  values = list(startparams = startparams,
                                control = list(trace = 6, maxit = 500)))
resultList[[toredo]] <- do.call(mmcflo2, 
                          args = mmfargs)


for (ind in toredo) {
  print(inputMatrix[ind, ])
  resultList[[ind]] <- do.call(master_mcflo, 
                            args = append(as.list(inputMatrix[ind, ]), 
                              values = list(control = list(trace = 6, maxit = 500))))
}


# Attach results to input matrix ------------------------------------------

resultVec <- resultList %>% 
  map_dbl(possibly(~.$metric, otherwise = NA))
failureVec <- is.na(resultList)

convgfun <- function(prt) {
  if (is.na(prt) || is.null(prt[["optim.info"]])) return(NA_integer_)
  convgcode <- prt[["optim.info"]][["convergence"]]
  convgcode
}
convgvec <- map_int(resultList, convgfun) # zeros indicate successful convergence

resultMatrix <- within(inputMatrix, {
  result = resultVec
  failed = failureVec
  })

cache("resultMatrix")

