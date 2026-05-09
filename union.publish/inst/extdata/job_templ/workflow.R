# ==========================================================================
# workflow of .{{{name}}}.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_.{{{name}}}. <- setClass("job_.{{{name}}}.", 
  contains = c("job"),
  prototype = prototype(
    pg = ".{{{name}}}.",
    info = c(""),
    cite = "",
    method = "",
    tag = ".{{{name}}}.",
    analysis = ""
    ))

job_.{{{name}}}. <- function()
{
  .job_.{{{name}}}.()
}

setMethod("step0", signature = c(x = "job_.{{{name}}}."),
  function(x){
    step_message("Prepare your data with function `job_.{{{name}}}.`.")
  })

setMethod("step1", signature = c(x = "job_.{{{name}}}."),
  function(x){
    step_message("Quality control (QC).")
    return(x)
  })

setMethod("set_remote", signature = c(x = "job_.{{{name}}}."),
  function(x, wd)
  {
    x$wd <- wd
    return(x)
  })
