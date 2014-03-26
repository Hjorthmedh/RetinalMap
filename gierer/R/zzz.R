## zzz.R --- general functions for loading/unloading package.

.onUnload <- function(libpath) {
  library.dynam.unload("gierer", libpath)
  cat("Detaching gierer package with shared objects\n")
}

#.onUnload is preferred to .last.Lib, see help file.
# use unloadNamespace("gierer") to detach the package.
