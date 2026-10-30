
.onLoad <- function(libname, pkgname) {
	tsrm_set_eigen_threads(1L)
}