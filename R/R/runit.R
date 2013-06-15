
runit.dir <- function()
{
    pkg <- "iproc"
    if (Sys.getenv("RCMDCHECK") == "FALSE") {
        path <- file.path(getwd(), "..", "inst", "unitTests")
    } else {
        path <- system.file(package=pkg, "unitTests")
    }
    path
}
