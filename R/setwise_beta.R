#' Calculate pairwise/ setwise beta diversities
#'
#' This is an internal helper function for beta_stand(). Beta-diversity is calculated for the number of samples
#' in x and the the given index. Therefore, x must have exactly N_samples rows. This is a security feature to ensure that users don't
#' try to calculate another set size than given by the number of rows in x. Use beta_stand()!
#'
#' If extrapolate = F, observed S will be returned if N_stand exceeds sample size.
#'
#' @param setsize number of sites in the matrix
#' @param func function returning the desired beta-diversity metric
#' @param ... additional argumets passed on to func
#' @param x A site by species matrix.
#'
#' @return A single numeric value for the beta diversity index.
#'

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
#' This function makes pairwise/ setwise combinations of sites in the community matrix x and calculates pairwise/ setwise beta diversities,
#' respectively. "Set" refers to a user-defined number of samples that will be used. This is essentially
#' the number of plots in a sample-based rarefaction. It's default is 2 (i.e. pairs).
#'
#' The calulation of the beta-diversity metric is done according to the function in "func". This should be a function that returns a
#' (dis-)similarity or beta diversity metric for a site-by species matrix like x. "func" is applied to every subset. Addidtional arguments will be passed on to "func".
#'
#' If the number of possible sets/ pairwise comparisons exceeds the value of max_combn. A random subsample of all combinations will
#' be drawn. The number of random samples can be adjusted using the argument resamples. It is not recommended increase the
#' value of max_combn. The function choose() "n over k" can be used to manually compute the number of possible subsets where
#' n is the number of samples/rows in x and k is setsize. If summarise = T, mean and variance of all setwise
#' comparisons are returned. Otherwise, the function returns all individual beta values.
#'
#'
#' @param x A site by species matrix.
#' @param setsize Number of samples per subset.
#' @param func a list of function names to be used for the metric calculation
#' @param args a list of additional arguments to be passed on to the functions in func
#' @param summarise Return mean and variance of all betas?
#' @param max_combn Number of combinations allowed before resampling is used instead
#' @param resamples Number of samples used if possible combinations exceeds max_combn
#' @param verbose Print notifications?
#' @import utils
#'
#' @return A data frame
#' @export
#'
#'
#' @examples
#' \donttest{
#' library(vegan)
#' data(BCI)
#' beta_pairwise<-
#' beta_stand(BCI, func = list("beta_true", "beta_SN", "beta_C"), setsize=2,
#' args = list(N=150, C=0.5))
#' }
#'
beta_stand <- function(x,
                       setsize= 2,
                       func = list("beta_true"),
                       args = NULL,
                       summarise = T,
                       max_combn = 10000,
                       resamples = 1000,
                       verbose= T
                       ) {

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
        pairs[[i]] <- x[combinations[, i],, drop=F]
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
                                   func=func[[i]],
                                   setsize=setsize)

        } else {
            out[[i]] = all_betas[[i]]
        }
    }

    if (summarise ==T) out <- do.call(rbind, out)
    else{
        out = do.call(cbind, out)
        colnames(out)=unlist(func)

    }
    return(out)
}
