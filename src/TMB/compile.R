
tmb_name <- "tsrm_TMBExports"

if (file.exists(paste0(tmb_name, ".cpp"))) {

    # eigene Flags von der Kommandozeile (optional)
    tmb_flags <- paste(commandArgs(trailingOnly = TRUE), collapse = " ")

    # Windows: grosse Objektdateien erlauben
    if (Sys.info()[["sysname"]] == "Windows") {
        tmb_flags <- paste(tmb_flags, "-Wa,-mbig-obj")
    }

    TMB::compile(file = paste0(tmb_name, ".cpp"),
                 flags = trimws(tmb_flags),
                 openmp = FALSE,
                 safebounds = FALSE, safeunload = FALSE)

    file.copy(from = paste0(tmb_name, .Platform$dynlib.ext),
              to = "..", overwrite = TRUE)
}