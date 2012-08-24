err.msg <- function(err.func = "edge",msg) {
  cat('\n')
  cat('\t')
  cat('ERROR in the',err.func,'function: ','\n')
  cat('\t',msg,'\n\n')
}