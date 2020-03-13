#load packages
library(tidyverse)
library(mobsim)
library(furrr)
library(betaC)

# Seed
set.seed(123)

# parameters for simulation
n_sim = 4000 #number of indivduals
cv = 1 #cv of lognormal SAD
sigma_common = 0.01 #width of cluster
s_pool = 10:500 #species pool sizes
motherpoints_common = c(1, 4, 10, 20, n_sim) #number of clusters, i.e. degree of aggregation
replicates = 3 #number of replicates of each combination
replicate = 1:replicates
n_quadrats = 4 #number of sampling quadrats
quadrat_area = 0.08 #area of sampling quadrats

## cross parameters
dat <- expand.grid(
    s_pool = s_pool,
    motherpoints_common = motherpoints_common,
    replicate = replicate
)
dat <- dat %>% mutate(
    n_sim = n_sim,
    cv =  cv,
    sigma_common = sigma_common,
    n_quadrats = n_quadrats,
    quadrat_area = quadrat_area,
    sigma_rare = sigma_common,
    motherpoints_rare =  motherpoints_common
)

plan(multiprocess)
## create total species pools
dat <- mutate(dat,
              total_com = future_pmap(list(s_pool = s_pool, n_sim = n_sim, cv = cv),
                                      function(s_pool, n_sim, cv) {
                                          sad_coef = list("cv_abund" = cv)
                                          return(sim_sad(
                                              s_pool = s_pool,
                                              n_sim = n_sim,
                                              sad_type = "lnorm",
                                              sad_coef = sad_coef
                                          ))
                                      }, .progress = T))

# calculate metrics of total species pools
dat <- dat %>% mutate(
    ENS_true = map_dbl(total_com, calc_PIE, ENS = T),
    N = map_dbl(total_com, sum),
    S_true = map_dbl(total_com, vegan::specnumber),
    S_common = S_true
)

# assign sigma values and mother points to the species.  here, rare and common species have the same values but this could be changed for other simulations
dat <- dat %>% mutate(sigma_vector = future_pmap(list(
    sigma_common = sigma_common,
    sigma_rare = sigma_rare,
    S_common = S_common,
    S = S_true
),
function(sigma_common, sigma_rare, S_common, S) {
    vect <- c(rep(sigma_common, S_common),
              rep(sigma_rare, S - S_common))
    return(vect)
}, .progress = T))

dat <- dat %>% mutate(motherpoint_vector = future_pmap(list(
    motherpoints_common = motherpoints_common,
    motherpoints_rare = motherpoints_rare,
    S_common = S_common,
    S = S_true
),
function(motherpoints_common,
         motherpoints_rare,
         S_common,
         S) {
    vect <- c(rep(motherpoints_common, S_common),
              rep(motherpoints_rare, S - S_common))
    return(vect)
}, .progress = T))

# make the communites spatially explicit
dat <- dat %>% mutate(spat_com = future_pmap(
    list(
        abund_vec = total_com,
        sigma = sigma_vector,
        mother_points = motherpoint_vector
    ),
    mobsim::sim_thomas_coords,
    .progress = T
))


# place sample quadrats in each simulated world
dat <-
    dat %>% mutate(
        alpha_samples = future_pmap(
            list(
                comm = spat_com,
                n_quadrats = n_quadrats,
                quadrat_area = quadrat_area
            ),
            sample_quadrats,
            plot = F,
            method = "grid",
            delta_x = 0.5,
            delta_y = 0.5,
            x0 = 0.1,
            y0 = 0.1,
            .progress = T
        )
    )
dat <-
    dat %>% mutate(alpha_samples =  map(alpha_samples, magrittr::extract2, 1))

# remove all the spatial point patterns to make the dataframe lighter
dat <- dat %>% select(-spat_com)
saveRDS(dat, "Simulations/simulation_data.rds")
