#' Calculate IBR curves at alpha and gamma scales and return them in a long format
#'
#' Please load tidyverse before using this.
#'
#' @param x a site by species abundance matrix
#'
#' @return a data.frame
#' @import tidyverse
#'

rarefy_long <- function(x) {
    require(tidyverse)
    if(is.matrix(x)==F) x=matrix(x,nrow = 1, byrow =T, dimnames= list("x", names(x)))
    alphas <-
        lapply(row.names(x), function(i)
            return(as.numeric(vegan::rarefy(
                x[i, ], sample = 1:sum(x[i, ])
            )))) %>%
        lapply(function(x)
            return(data.frame(
                S_n = as.numeric(x), N = 1:length(x)
            )))
    names(alphas) <- rownames(x)
    alphas <- alphas %>% plyr::ldply(.id = "Curve")
    alphas$type = "minor"
    mean_alpha <-
        data.frame(
            Curve = "mean_alpha",
            S_n = colMeans(as.matrix(vegan::rarefy(
                x, 1:min(rowSums(x))
            ))),
            N = 1:min(rowSums(x)),
            type = "major"
        )
    gamma <-
        data.frame(
            Curve = "gamma",
            S_n = as.numeric(vegan::rarefy(colSums(x), 1:sum(x))),
            N = 1:sum(x),
            type = "major"
        )
    out = alphas %>% full_join(mean_alpha, by = c("Curve", "S_n", "N",  "type")) %>% full_join(gamma, by = c("Curve", "S_n", "N", "type"))

}
