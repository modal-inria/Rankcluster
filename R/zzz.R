#' @export
.onAttach <- function(lib, pkg) {
  packageStartupMessage(
    paste(
      "WARNING: Since Rancluster 0.92",
      "ranks have to be given in the ranking notation (see convertRank function), with the following convention:",
      "\n- missing positions are replaced by 0",
      "\n- tied positions are replaced by the lowest position they share"
    )
  )
}
