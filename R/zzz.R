.onLoad <- function(libname, pkgname){

  flogger <- logger::layout_glue_generator(format =
                                             paste('{time} {node} {pid} ',
                                                   '{namespace} {fn} ',
                                                   '{level}:  {msg}',
                                                   sep = ""))

  logger::log_layout(flogger, index = 1)
}

.onAttach <- function(libname, pkgname){
  packageStartupMessage(paste0(
    "              ****                 \n"
    ,"The CausalGPS version 0.5.x will undergo significant design changes. "
    ,"Please be prepared to update your code accordingly, or consider "
    ,"continuing with version 0.4.x until you are ready to update.\n"
    ,"              ****                 \n"
  ))
}
