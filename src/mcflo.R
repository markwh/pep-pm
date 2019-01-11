# mcflo.R
# Mark Hagemann
# 12/20/2018
# Transferring mcflo analysis from notebooks 20181129, 20181204

# Matrix of inputs to run through master function -------------------------

# Note: master function is now in optim.R.

cases <- names(reachdata)
cases <- cases[!grepl("^StLaw", cases)] # St Lawrence doesn't have A data
cases <- cases[!grepl("^Tanana", cases)] # Tanana doesn't have A data either

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

# torun <- which(failureVec)
torun <- 1:nrow(inputMatrix)

resultList[torun] <- pmap(.l = inputMatrix[torun, ], #[-1:-342, ],
                   .f = possibly(master_mcflo, otherwise = NA))
                   # .f = possibly(mmcflo2, otherwise = NA))
                   # control = list(trace = NULL))
cache("resultList")


# manually redo bad optimizations

toredo <- which(convgvec == 1)
# toredo <- which(failureVec)
# toredo <- 579

inputMatrix[toredo, ]

inps1 <- mcflo_inps(swotlist = reachdata$Connecticut, fl = "bam_man", 
           mc = "omniscient", metric = "rrmse", method = "optim")
opt1 <- mcflo_optim(inps1, control = list(trace = 3, maxit = 500))
# opt2 <- mcflo_optim(inps1, startparams = opt1$params, control = list(trace = 3, maxit = 500))

mmfargs <- append(as.list(inputMatrix[toredo, ]), 
                  values = list(control = list(trace = 3, 
                                               maxit = 500)))#,
                                # startparams = opt1$params))
resultList[[toredo]] <- do.call(master_mcflo, 
                                # args = inputMatrix[toredo, ])
                                args = mmfargs)


for (ind in toredo) {
  print(inputMatrix[ind, ])
  startpars <- resultList[[ind]]$params
  resultList[[ind]] <- do.call(master_mcflo, 
                            args = append(as.list(inputMatrix[ind, ]), 
                              values = list(
                                startparams = startpars,
                                control = list(trace = 3, maxit = 500))))
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

