

## build graph 

# load the necessary graph function
source("../build_graph.R")
library(org.Hs.eg.db)

# select_data_ge: the dataset processed by feature screening method
# sel_no        : number of features 
# data_ge_probe : gene probe dataset
# return value  : edge information for the features processed by feature screening
graph_extract = function(sel_data_ge,sel_no,data_ge_probe){
  diff.ids           = data.frame(matrix(ncol = 3, nrow = sel_no))
  colnames(diff.ids) = c("probe","symbol","egid")
  diff.ids$probe  = colnames(sel_data_ge)[1:sel_no]
  diff.ids$symbol = data_ge_probe[data_ge_probe$ProbeSet %in% diff.ids$probe, 3]
  # extract egid and keggid
  sym = diff.ids$symbol
  EG_IDs = mget(sym, revmap(org.Hs.egSYMBOL),ifnotfound=NA)
  # build graph 
  graph= ENT2Graph(EG_IDs,edge=TRUE) 
  return(graph$E)
}


