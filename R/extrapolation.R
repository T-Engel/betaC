#' Rarefaction curve (inter- and extrapolation)
#'
#' Returns the expected number of species for a samplesize m. Interpolation and extrapolation is possible.
#' This function was taken from the rPackage iNEXT (https://github.com/JohnsonHsieh/iNEXT).
#'
#' @param x integer vector of abundances
#' @param m sample size. can be a vector
#'
#' @return
#' @export numerical (vector)
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


beta_C_extra <- function(x, C, extrapolation= T, interrupt=T) {
    x <- as.matrix(x)
    total <- colSums(x)
    N <- round(invChat(total, C))
    C_max=C_target_extra(x, extrapolation=extrapolation)
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

C_target_extra <- function(x, extrapolation=T) {
    factor=ifelse(extrapolation,2,1)
    x <- as.matrix(x)
    n = min(factor*rowSums(x))
    out <- Chat(colSums(x), n)
    return(out)
}

C_target_factor<- function(x, factor=1) {
    x <- as.matrix(x)
    n = min(factor*rowSums(x))
    out <- Chat(colSums(x), n)
    return(out)
}



calc_chao<-function(x){
    if (!is.numeric(x) & !is.matrix(x) & !is.data.frame(x))
        stop("invalid data structure")
    if (is.matrix(x) | is.data.frame(x)) {
        S_Chao1 = apply(x, 1, calc_chao)
    }
    else {
        n = sum(x)
        D = sum(x > 0)
        f1 = sum(x == 1)
        f2 = sum(x == 2)
        if (f1 > 0 & f2 > 0)
            S_Chao1 = D + (n - 1)/n * f1^2/(2 * f2)
        else if (f1 > 1 & f2 == 0)
            S_Chao1 = D + (n - 1)/n * f1 * (f1 - 1)/(2 * (f2 +
                                                              1))
        else S_Chao1 = NA
    }
    return(S_Chao1)

}
beta_asymptote <- function(x) {
    x <- as.matrix(x)
    total <- colSums(x)
    gamma_value = as.numeric(calc_chao(total))
    alpha_value = mean(calc_chao(x))
    beta = gamma_value / alpha_value


    return(beta)

}
