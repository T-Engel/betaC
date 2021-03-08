##################################################################
# Supplementary material S4: Novel code
##################################################################

# Below you find the new code for the calculation of betaC. Our new functions
# "beta_C" and "C_target" use some functions originally written by Hsieh et al (2016).
#
# Reference:
# Hsieh, T.C., Ma, K.H. and Chao, A. (2016),
# iNEXT: an R package for rarefaction and extrapolation of species diversity (Hill numbers).
# Methods Ecol Evol, 7: 1451-1456. https://doi.org/10.1111/2041-210X.12613


#' Calculate expected sample coverage C_hat
#'
#' Returns expected sample coverage  of a sample x for a smaller than observed sample size m (Chao & Jost, 2012).
#' This code was copied from INEXT's internal function Chat.Ind() (Hsieh et al 2016).
#'
#' @param x integer vector (species abundances)
#' @param m integer. (smaller than observed sample size)
#'
#' @return a numeric value.
#' @export
#'
#' @examples
#' \donttest{
#' library(vegan)
#' data(BCI)
#'
#' # What is the expected coverage corresponding to a sample size of 50 at the gamma scale?
#' Chat(colSums(BCI), 50)
#' }
Chat <- function (x, m)
{
    x <- x[x > 0]
    n <- sum(x)
    f1 <- sum(x == 1)
    f2 <- sum(x == 2)
    f0.hat <- ifelse(f2 == 0, (n - 1) / n * f1 * (f1 - 1) / 2, (n -
                                                                    1) / n * f1 ^
                         2 / 2 / f2)
    A <- ifelse(f1 > 0, n * f0.hat / (n * f0.hat + f1), 1)
    Sub <- function(m) {
        if (m < n) {
            xx <- x[(n - x) >= m]
            out <- 1 - sum(xx / n * exp(
                lgamma(n - xx + 1) - lgamma(n -
                                                xx - m + 1) - lgamma(n) + lgamma(n - m)
            ))
        }
        if (m == n)
            out <- 1 - f1 / n * A
        if (m > n)
            out <- 1 - f1 / n * A ^ (m - n + 1)
        out
    }
    sapply(m, Sub)
}

#' Number of individuals corresponding to a desired coverage (inverse C_hat)
#'
#' If you wanted to resample a vector to a certain expected sample coverage, how many individuals would you have to draw?
#' This is C_hat solved for the number of individuals.This code is a modification of INEXT's internal function invChat.Ind() (Hsieh et al 2016).
#'
#' @param x integer vector.
#' @param C numeric. between 0 and 1
#'
#' @return a numeric value
#' @export
#' @import stats
#' @examples
#' \donttest{
#' library(vegan)
#' data(BCI)
#'
#' # What sample size corresponds to an expected sample coverage of 55%?
#' invChat(colSums(BCI), 0.55)
#' }
#'
invChat <- function (x, C)
{
    m <- NULL
    n <- sum(x)
    refC <- Chat(x, n)
    f <- function(m, C)
        abs(Chat(x, m) - C)
    # for interpolation
    if (refC > C) {
        opt <- optimize(f,
                        C = C,
                        lower = 0,
                        upper = sum(x))
        mm <- opt$minimum
    }
    # for extrapolation
    if (refC <= C) {
        f1 <- sum(x == 1)
        f2 <- sum(x == 2)
        if (f1 > 0 & f2 > 0) {
            A <- (n - 1) * f1 / ((n - 1) * f1 + 2 * f2)
        }
        if (f1 > 1 & f2 == 0) {
            A <- (n - 1) * (f1 - 1) / ((n - 1) * (f1 - 1) + 2)
        }
        if (f1 == 1 & f2 == 0) {
            A <- 1
        }
        if (f1 == 0 & f2 == 0) {
            A <- 1
        }
        mm <- (log(n / f1) + log(1 - C)) / log(A) - 1
        mm <- n + mm

    }
    if (mm > 2 * n)
        warning(
            "The maximum size of the extrapolation exceeds double reference sample size, the results for q = 0 may be subject to large prediction bias."
        )
    return(mm)
}

#' Rarefaction curve (inter- and extrapolation)
#'
#' Returns the expected number of species for a samplesize m. Interpolation and extrapolation is possible.
#' This code was copied from the rPackage iNEXT (Hsieh et al 2016)
#' @param x integer vector of abundances
#' @param m sample size. can be a vector
#'
#' @return numerical (vector)
#' @export
#'
#' @examples
D0.hat <- function(x, m) {
    x <- x[x > 0]
    n <- sum(x)
    Sub <- function(m) {
        if (m <= n) {
            Fun <- function(x) {
                if (x <= (n - m))
                    exp(lgamma(n - x + 1) + lgamma(n - m + 1) -
                            lgamma(n - x - m + 1) - lgamma(n + 1))
                else
                    0
            }
            sum(1 - sapply(x, Fun))
        }
        else {
            Sobs <- sum(x > 0)
            f1 <- sum(x == 1)
            f2 <- sum(x == 2)
            f0.hat <-
                ifelse(f2 == 0, (n - 1) / n * f1 * (f1 - 1) / 2, (n - 1) / n * f1 ^ 2 /
                           2 / f2)
            A <- n * f0.hat / (n * f0.hat + f1)
            ifelse(f1 == 0, Sobs, Sobs + f0.hat * (1 - A ^ (m - n)))
        }
    }
    sapply(m, Sub)
}


#' Calculate beta_C
#'
#' Beta_C uses coverage-based rarefaction to standardize beta-diversity.
#'
#' @param x a site by species matrix
#' @param C target coverage. value between 0 and 1.
#' @param extrapolation logical. should extrapolation be used?
#' @param interrupt logical. Should the function throw an error when C exceeds the maximum recommendable coverage?
#'
#' @return a numeric value
#' @export
#'
#' @examples
#' \donttest{
#' library(vegan)
#' data(BCI)
#'
#' # What is beta_C for a coverage value of 60%?
#' beta_C(BCI,C = 0.6)
#' }
beta_C <- function(x,
                   C,
                   extrapolation = T,
                   interrupt = T) {
    x <- as.matrix(x)
    total <- colSums(x)
    N <- round(invChat(total, C))
    C_max = C_target(x, factor = ifelse(extrapolation, 2, 1))
    if (C > C_max & interrupt == T) {
        if (extrapolation == F) {
            stop(
                paste0(
                    "Coverage exceeds the maximum possible value for interpolation (i.e. C_target = ",
                    round(C_max, 4),
                    "). Use extrapolation or reduce the value of C."
                )
            )
        } else{
            stop(
                paste0(
                    "Coverage exceeds the maximum possible value recommendable for extrapolation (i.e. C_target = ",
                    round(C_max, 4),
                    "). Reduce the value of C."
                )
            )
        }
    }
    if (N > 1) {
        gamma_value = D0.hat(total, N)
        alpha_value = mean(apply(x, 1, D0.hat, m = N))
        beta = gamma_value / alpha_value
    } else {
        beta = NA
    }
    attr(beta, "C") = C
    attr(beta, "N") = N
    return(beta)

}

#' Calculate the recommended maximum coverage value for the computation of beta_C from a site by species matrix
#'
#' This returns the coverage of x at the gamma scale that corresponds to the smallest observed sample size at the alpha scale times an extrapolation factor.
#' The default (factor = 2) allows for extrapolation up to 2 times the observed sample size of the smallest alpha sample. For factor= 1, only interpolation is applied.
#' Its not recommendable to use factors largers than 2.
#'
#' @param x a site by specie matrix
#' @param factor numeric. how far
#'
#' @return numeric value
#' @export
#'
#' @examples
#' \donttest{
#' library(vegan)
#' data(BCI)
#'
#' # What is the largest possible C that I can use to calculate beta_C for my site by species matrix?
#' C_target(BCI)
#' }
#'

C_target <- function(x, factor = 2) {
    x <- as.matrix(x)
    n = min(factor * rowSums(x))
    out <- Chat(colSums(x), n)
    return(out)
}

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
    beta = as.numeric(match.fun(func)(x, ...))
    return(beta)
}



#' Calculate pairwise/ setwise beta diversities
#'
#' This function makes pairwise/ setwise combinations of sites in the community matrix x and calculates pairwise/ setwise beta diversities,
#' respectively. "Set" refers to a user-defined number of samples that will be used. This is essentially
#' the number of plots in a sample-based rarefaction. It's default is 2 (i.e. pairs).
#'
#' The calculation of the beta-diversity metric is done according to the function in "func". This should be a function that returns a
#' (dis-)similarity or beta diversity metric for a site-by species matrix like x. "func" is applied to every subset. Additional arguments will be passed on to "func".
#'
#' If the number of possible sets/ pairwise comparisons exceeds the value of max_combn. A random subsample of all combinations will
#' be drawn. The number of random samples can be adjusted using the argument resamples. It is not recommended to increase the
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
#' beta_stand(BCI, func = list( "beta_C"), setsize=2, summarise=F,
#' args = list(C=0.5))
#' }
#'
beta_stand <- function(x,
                       setsize = 2,
                       func = list("beta_true"),
                       args = NULL,
                       summarise = T,
                       max_combn = 10000,
                       resamples = 1000,
                       verbose = T) {
    x <- as.matrix(x)
    rownames(x) = NULL

    N <- rowSums(x)

    if (choose(nrow(x), setsize) > max_combn) {
        if (verbose)
            print("The number of possible subsets is too large. They will be resampled")
        combinations <-
            sapply(1:resamples, function(y, nrows = nrow(x))
                sample(1:nrows, setsize))
    } else{
        combinations <- combn(1:nrow(x), setsize)
    }


    pairs <- vector("list", length = ncol(combinations))
    for (i in 1:ncol(combinations)) {
        pairs[[i]] <- x[combinations[, i], , drop = F]
        pairs[[i]] <-
            pairs[[i]][, colSums(pairs[[i]]) > 0, drop = F]
    }
    all_betas <- as.list(rep(NA, length(func)))
    out <- as.list(rep(NA, length(func)))

    for (i in 1:length(func)) {
        all_betas[[i]] <-
            sapply(
                pairs,
                FUN = function(x, args, ...)
                    R.utils::doCall(match.fun(func[[i]]), x = x, args = args),
                args = args
            )

        if (summarise == T) {
            out[[i]] <- data.frame(
                mean = mean(all_betas[[i]], na.rm = T),
                var = var(all_betas[[i]]),
                func = func[[i]],
                setsize = setsize
            )

        } else {
            out[[i]] = all_betas[[i]]
        }
    }

    if (summarise == T)
        out <- do.call(rbind, out)
    else{
        out = do.call(cbind, out)
        colnames(out) = unlist(func)

    }
    return(out)
}
