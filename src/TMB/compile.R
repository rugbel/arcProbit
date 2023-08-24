# compile tmb models
invisible(sapply(Sys.glob("*.cpp"),
                 TMB::compile,
                 safebounds = TRUE, safeunload = TRUE))
# copy dynlibs to src
invisible(file.copy(from = Sys.glob(paste0("*", .Platform$dynlib.ext)),
                    to = "..", overwrite = TRUE))
# cleanup done in ../Makevars[.win]

## args <- commandArgs(trailingOnly = TRUE)
## TMB::compile(args, safebounds = FALSE, safeunload = FALSE)
