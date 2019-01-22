# mcflo-omni.R
# Parameters obtained using omniscient mc, testing using all variants
# Because optimized mcflo using different mc's led to crazy fl results,
#  even though mcfl results were awesome.
# modified from notebook20190118.Rmd



# list of objective functions

inp_wrapper <- function(case, 
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
  
  optim_inps
}


funlist <- list()
for (i in 1:nrow(resultMatrix)) {
  # cat(i)
  argrow <- resultMatrix[i, ]
  funlist[[i]] <- with(argrow, inp_wrapper(
    case = case, fl = fl, mc = mc, metric = metric, 
    method = method, msg = msg
  ))$f
}


# new results using omnisicent-mc parameters ------------------------------

# I need a table for matching non-omnisicient parameters to omniscient parameters

omnimat <- resultMatrix %>% 
  filter(mc == "omniscient") %>% 
  mutate(omnirow = msg) %>% 
  select(case:method, -mc, omnirow)

matchmat <- resultMatrix %>% 
  mutate(origrow = msg) %>% 
  select(case:method, origrow) %>% 
  left_join(omnimat)

# Get result from funlist[[origrow]], using parameters from resultList[[omnirow]]

resultList_omni <- list()
for (i in 1:nrow(matchmat)) {
  cat(i, " ")
  origi <- matchmat[i, "origrow"]
  omnii <- matchmat[i, "omnirow"]
  funi <- funlist[[origi]]
  paramsi <- resultList[[omnii]]$params
  resultList_omni[[i]] <- funi(paramsi)
}

resultVec_omni <- unlist(resultList_omni)


resultMatrix_omni <- resultMatrix %>% 
  mutate(result = resultVec_omni)


matchmat_omni <- matchmat
cache("matchmat_omni")
cache("resultMatrix_omni")
