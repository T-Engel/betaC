#' Calculate sample coverage
#'
#' Sample coverage is calculated according to Chao & Jost (2012).
#'
#' @param m site by species matrix
#'
#' @return a vector
#' @export
#'
#' @examples
#' \donttest{
#' library(vegan)
#' data(BCI)
#'
#' # at the alpha scale
#' coverage(BCI)
#'
#' # at the gamma scale
#' coverage(colSums(BCI))
#' }
#'
coverage <- function(m) {
  if (is.null(dim(m)))
    m = matrix(m, ncol = length(m))
  n <- apply(m, 1, sum)
  f1 <- apply(m, 1, function(x)
    length(x[x == 1]))
  f2 <- apply(m, 1, function(x)
    length(x[x == 2]))
  cover = 1 - (f1 / n) * (((n - 1) * f1) / ((n - 1) * f1 + 2 * f2))
  return(cover)
}


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
#' This is C_hat solved for the number of individuals. This code is a modification INEXT's internal function invChat.Ind() (Hsieh et al 2016).
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
#' This function was taken from the rPackage iNEXT (https://github.com/JohnsonHsieh/iNEXT).
#'
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
        else 0
      }
      sum(1 - sapply(x, Fun))
    }
    else {
      Sobs <- sum(x > 0)
      f1 <- sum(x == 1)
      f2 <- sum(x == 2)
      f0.hat <- ifelse(f2 == 0, (n - 1)/n * f1 * (f1 - 1)/2, (n - 1)/n * f1^2/2/f2)
      A <- n * f0.hat/(n * f0.hat + f1)
      ifelse(f1 == 0, Sobs, Sobs + f0.hat * (1 - A^(m - n)))
    }
  }
  sapply(m, Sub)
}





#' Calculate beta_C
#'
#' Beta_C uses coverage-based rarefaction to quantify the non-random component in beta-diversity.
#' The partitioning is done at the number of individuals that corresponds to a sample coverage of C at the gamma scale.
#'
#' @param x a site by species matrix
#' @param C target coverage. value between 0 and 1.
#' @param extrapolation logical. should extrapolation be used?
#' @param interrupt logical. SHould the function throw an error when C exceeds the maximum recommendable coverage?
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
beta_C <- function(x, C, extrapolation= T, interrupt=T) {
  x <- as.matrix(x)
  total <- colSums(x)
  N <- round(invChat(total, C))
  C_max=C_target(x, factor=ifelse(extrapolation,2,1))
  if(C>C_max& interrupt==T){
    if(extrapolation==F){
      stop(paste0("Coverage exceeds the maximum possible value for interpolation (i.e. C_target = ",round(C_max,4),"). Use extrapolation or reduce the value of C.")
      )
    }else{
      stop(paste0("Coverage exceeds the maximum possible value recommendable for extrapolation (i.e. C_target = ",round(C_max,4),"). Reduce the value of C.")
      )
    }
  }
  if(N>1){
    gamma_value = D0.hat(total, N)
    alpha_value = mean(apply(x,1, D0.hat, m=N))
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
#' This returns the coverage of x at the gamma scale that corresponds to the smalles observed sample size at the alpha scale times an extrapolation factor.
#' The default (factor = 2) allows for extrapolation up to 2 times the observed sample size of the smallest alpha sample. For factor= 1, only interpolation is applied.
#' its not recommendable to use factors largers than 2.
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

C_target<- function(x, factor=2) {
  x <- as.matrix(x)
  n = min(factor*rowSums(x))
  out <- Chat(colSums(x), n)
  return(out)
}


#' Calculate beta_Sn
#'
#' @param x a site by species matrix
#' @param N an integer value
#'
#' @return a numeric value
#' @export
#'
#' @examples
#' \donttest{
#' library(vegan)
#' data(BCI)
#' beta_SN(BCI, 50)
#' }
beta_SN<-function(x, N){
  x<-as.matrix(x)
  total<-colSums(x)
  C= Chat(total,N)

  if(N>1){
    gamma_value = as.numeric(vegan::rarefy(total, N))
    alpha_value = mean(vegan::rarefy(x, N))
    beta = gamma_value / alpha_value
  } else {
    beta = NA
  }

  attr(beta, "C") = C
  attr(beta, "N") = N
  return(beta)

}


#' Calculate Beta diversity
#'
#' @param x a site by species matrix
#' @param transformation apply Jost's transformation?
#'
#' @return numeric.
#' @export
#'
#' @examples
#' \donttest{
#' library(vegan)
#' data(BCI)
#' beta_true(BCI)
#' }
beta_true= function(x, transformation= F){
    x= as.matrix(x)
    alpha= mean(vegan::specnumber(x))
    gamma=vegan::specnumber(colSums(x))
    beta =gamma/alpha
    if(transformation) beta= 1-(1/beta)
    return(beta)
}


#' Calculate PIE and Effective number of species (S_PIE)
#'
#' This function was copied from  the package \href{https://github.com/MoBiodiv/mobr/}{mobr}.
#'
#' @param x site by species matrix
#' @param ENS return effective number of species?
#'
#' @return
#' @export
#'

calc_PIE = function(x, ENS=FALSE) {
  x = drop(as.matrix(x))
  if (any(x < 0, na.rm = TRUE))
    stop("input data must be non-negative")
  if (length(dim(x)) > 1) {
    total = apply(x, 1, sum)
    S = apply(x, 1, function(x) return(sum(x > 0)))
    x = sweep(x, 1, total, "/")
  } else {
    total = sum(x)
    S = sum(x > 0)
    x = x / total
  }
  x = x * x
  if (length(dim(x)) > 1) {
    H = rowSums(x, na.rm = TRUE)
  } else {
    H = sum(x, na.rm = TRUE)
  }
  # calculate PIE without replacement (for total >= 2)
  H = ifelse(total < 2, NA, (total / (total - 1) * (1 - H)))
  if (ENS) {
    # convert to effective number of species (except for PIE == 1)
    H = ifelse(H==1| S == total, NA, (1/ (1-H)))
  }
  return(H)
}


