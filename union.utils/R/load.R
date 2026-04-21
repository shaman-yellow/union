# ==========================================================================
# load union
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

list_unions <- function(path = dirname(find.package("union.utils"))) {
  base <- file.path(path, c("union.publish", "union.project"))
  if (!all(file.exists(base))) {
    stop('!all(file.exists(base)).')
  }
  series <- list.files(path, "union\\.series\\.", full.names = TRUE)
  c(base, series)
}

load_unions <- function(path = dirname(find.package("union.utils")), 
  pkgs = list_unions(path))
{
  for (i in pkgs) {
    devtools::load_all(i)
  }
}



