

# Load Pepsi datasets -----------------------------------------------------
# 
# ncfiles <- list.files("../SWOT/data/ncdata/", pattern = "\\.nc$", 
#                       full.names = TRUE)
# casenames <- list.files("../SWOT/data/ncdata/", pattern = "\\.nc$") %>% 
#   gsub("\\.nc$", "", .)
# 
# reachdata <- map(ncfiles, ~nc_reach(.)) %>% 
#   setNames(casenames)


# Reachdata, see munge script in mcfli-swotr project ----------------------


load("../mcfli-swotr/cache/reachdata.RData")

cache("reachdata")
