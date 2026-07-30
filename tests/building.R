#' script_path <- rstudioapi::getActiveDocumentContext()$path
#' print(script_path)
#' working_directory <- dirname(dirname(script_path))
#' setwd(working_directory)
#'
#' pkgdown::build_site()
#'
#' pkgload::load_all()
#' rstantools::use_rstan()
#' devtools::document()
#'
#' devtools::build_manual()
#'
#' devtools::build()
#'
#' detach("package:ForceChoice", unload = TRUE, character.only = TRUE)
#' .rs.restartR()
#' devtools::install()
#' devtools::check(manual = TRUE)
#'
#' a <- devtools::spell_check(vignettes = TRUE, use_wordlist = TRUE)
#' a
#' a[1]
#'
#' devtools::check_rhub()
#'
#' devtools::check_win_devel()
#'
#' devtools::release()
#' devtools::submit_cran()
#'
#'
#' NEWS.md
#' DESCRIPTION
#' cran-comments.md
#'
#'
#' pack <- available.packages()
#'
#' which(rownames(pack) == "Qval")
#'
#'
# sum <- 0
# library(cranlogs)
# downloads <- cran_downloads(packages = "LCAP", from = "2026-01-22", to = "last-day")
# sum <- sum + sum(downloads$count)
# downloads <- cran_downloads(packages = "Qval", from = "2014-06-30", to = "last-day")
# sum <- sum + sum(downloads$count)
# downloads <- cran_downloads(packages = "EFAfactors", from = "2024-02-04", to = "last-day")
# sum <- sum + sum(downloads$count)
# downloads <- cran_downloads(packages = "LSTMfactors", from = "2024-09-25", to = "last-day")
# sum <- sum + sum(downloads$count)
# sum
#'
#'
