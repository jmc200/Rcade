
.onAttach <- function(libname, pkgname) {

	packageStartupMessage("Welcome to Rcade - version ", packageDescription("Rcade", fields="Version"))
	packageStartupMessage('If you are new to Rcade, please consider reading the vignette through the command: vignette("Rcade").')
}

