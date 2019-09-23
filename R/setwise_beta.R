#' Calculate pairwise/ setwise beta diversities
#'
#' This is an internal helper function for beta_stand(). Beta-diversities are calculated for the number of samples
#' in x and the the given index. Therefore, x must have exactly N_samples rows. This is a security feature to ensure that users don't
#' try to calculate another set size than given by the number of rows in x. Use beta_stand()!
#'
#' If extrapolate = F, observed S will be returned if N_stand exceeds sample size.
#'
#' @param x A site by species matrix.
#' @param index One of the following indices: "S_n", "S_PIE", "S_cov"
#' @param N_stand Number of individuals used if index = "S_n.
#' @param C_stand Coverage value used if index = "S_cov.
#' @param extrapolate Logical. Use extrapolation?
#' @param N_samples Number of samples in the set. Has to match the number of rows in x.
#'
#' @return A single numeric value for the beta diversity index.
#'
#'
#' @examples
setwise_beta <- function(x,
                         setsize = 2,
                         func,
                         ...) {
    if (dim(x)[1] != setsize)
        stop("x has to be a sites by species matrix with exactly N_samples sites!")
    beta=as.numeric(match.fun(func)(x, ...))
    return(beta)
}



#' Calculate pairwise/ setwise beta diversities
#'
#' This function makes all pairwise/ setwise combinations of sites in the community matrix x and calculates pairwise/ setwise beta diversities,
#' respectively. "Set" refers to a user-defined number of samples that will be used calculate gamma diversities. This is essentially
#' the number of plots in a sample-based rarefaction. It's default is 2 (i.e. pairs). The function accepts the following indices: "S_PIE", "S_n", "S_cov".
#' The latter two require a value for the standardisation of number of individuals (N_stand) and coverage (S_cov).
#' If these arguments are not supplied, the function will standardise to the number of individuals of the smallest sample (if extrapolate = FALSE),
#' or twice this number (if extrapolate = TRUE). For index = "S_cov" the corresponding expected coverages will be used if C_stand= NULL.
#'
#' If the number of possible sets/ pairwise comparisons exceeds the value of max_combn. A random subsample of all combinations will
#' be drawn. The number of random samples can be adjusted using the argument resamples. It is not recommended increase the
#' value of max_combn. The function choose() "n over k" can be used to manually compute the number of possible subsets where
#' n is the number of samples/rows in x and k is setsize. If summarise = T, mean and variance of all setwise
#' comparisons are returned. Otherwise, the function returns all individual beta values.
#'
#' The default value for setsize is 2 which means that pairwise beta diversities are returned.
#' Pairwise beta_Sn values are independent of the gamma_scale diversity and sampling effort.
#' Therefore, they allow for meaningfull comparisons between groups with different numbers of
#' samples (see Marion et al, 2017).
#'
#'
#' @param x A site by species matrix.
#' @inheritParams setwise_beta
#' @param setsize Number of samples per subset.
#' @param summarise Return mean and variance of all betas?
#' @param resamples Number of samples used if possible combinations exceeds max_combn
#' @param max_combn Number of combinations allowed before resampling is used instead.
#' @param verbose Print notifications?
#'
#' @return A data frame.
#' @export
#'
#' @examples
beta_stand <- function(x,
                       func = list("beta_true", "beta_SN"),
                       setsize= 2,
                       summarise = T,
                       resamples = 1000,
                       max_combn = 1000,
                       verbose= T,
                       args) {

    x <- as.matrix(x)
    rownames(x) = NULL

    N <- rowSums(x)

    if (choose(nrow(x), setsize) > max_combn) {
        if(verbose) print("The number of possible subsets is too large. They will be resampled")
        combinations <-
            sapply(1:resamples, function(y, nrows = nrow(x))
                sample(1:nrows, setsize))
    } else{
        combinations <- combn(1:nrow(x), setsize)
    }


    pairs <- vector("list", length = ncol(combinations))
    for (i in 1:ncol(combinations)) {
        pairs[[i]] <- x[combinations[, i],]
        pairs[[i]] <- pairs[[i]][, colSums(pairs[[i]]) > 0, drop = F]
    }
    all_betas<- as.list(rep(NA, length(func)))
    out<- as.list(rep(NA, length(func)))

    for (i in 1:length(func)){
        all_betas[[i]] <-
            sapply(
                pairs,
                FUN=function(x,args,...)  R.utils::doCall(match.fun(func[[i]]),x=x,args=args),
                args=args
            )

        if (summarise == T) {
            out[[i]] <- data.frame(mean=mean(all_betas[[i]], na.rm = T),
                                   var=var(all_betas[[i]]),
                                   func=i,
                                   setsize=setsize)

        } else {
            out[[i]] = all_betas[[i]]
        }
    }

    if (summarise ==T) out <- do.call(rbind, out)
    else{
        out = do.call(cbind, out)
    }
    return(out)
}
