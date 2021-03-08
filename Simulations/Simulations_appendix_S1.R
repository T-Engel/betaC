#load packages
library(tidyverse)
library(mobsim)
library(furrr)
library(betaC)

# Seed
set.seed(123)

# parameters for simulation
n_sim = 4000 #number of indivduals
alpha= 1:100 #cv of lognormal SAD
sigma_common = c(0.01, 0.1,0.2,0.4, 0.8,1) #width of cluster
s_pool = NA #species pool sizes
motherpoints_common = 1 #number of clusters, i.e. degree of aggregation
replicates = 7 #number of replicates of each combination
replicate = 1:replicates
n_quadrats = 4 #number of sampling quadrats
quadrat_area = 0.08 #area of sampling quadrats

## cross parameters
dat <- expand.grid(
    alpha = alpha,
    sigma_common = sigma_common,
    replicate = replicate
)
dat <- dat %>% mutate(
    n_sim = n_sim,
    motherpoints_common = motherpoints_common,
    n_quadrats = n_quadrats,
    quadrat_area = quadrat_area,
    sigma_rare = sigma_common,
    motherpoints_rare =  motherpoints_common
)

plan(multiprocess)
## create total species pools
dat <- mutate(dat,
              total_com = future_pmap(list( n_sim = n_sim, alpha = alpha),
                                      function(n_sim, alpha) {
                                          sad_coef = list("N" = n_sim, "alpha" = alpha)
                                          return(sim_sad(
                                              s_pool = NULL,
                                              n_sim = n_sim,
                                              sad_type = "ls",
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

# remove 0s
dat<- dat %>%
    mutate(alpha_samples = map(alpha_samples, function(x) return(x[,colSums(x)>0])))
# remove all the spatial point patterns to make the dataframe lighter
dat <- dat %>% select(-spat_com)



# Whittaker's beta
plan(multiprocess)
dat<- dat %>% mutate(beta_true= map_dbl(alpha_samples, beta_true))

# remove cases with too few individuals
dat<- dat %>% mutate(N_min_alpha= map_dbl(alpha_samples, function(x) return(min(rowSums(x))))) # removing cases with empty samples
N_stand=50 ##
dat <- dat %>%
    filter(N_min_alpha>= (N_stand))

dat<- dat %>% mutate(N_obs=
                         map_dbl(alpha_samples,
                                 function(x) return(sum(x)))) # removing cases with empty samples



# beta_Sn
dat<-dat %>% mutate(beta_Sn= future_map_dbl(alpha_samples, beta_SN, N=N_stand, .progress=T))

# beta_C
dat<- dat %>% mutate(C_target= future_map_dbl(alpha_samples, C_target, .progress= T))
C= min(dat$C_target)

dat<-dat %>% mutate(beta_CC= future_map_dbl(alpha_samples, beta_C, C=C, .progress=T))

dat<- dat %>% mutate(sigma_common= as.factor(sigma_common))

dat %>%
    ggplot( aes(alpha,beta_true, group= sigma_common, col= sigma_common)) +
    geom_point(alpha=0.2)+
    geom_smooth()+
    labs(title= paste0("C = ", round(C, 3)))

dat %>%
    filter(N_min_alpha/N_obs >0.05) %>%  # exclude cases where quadrats fall into empty patches
    pivot_longer(c(beta_CC, beta_Sn, beta_true)) %>%
    ggplot( aes(alpha,value, group= sigma_common, col= sigma_common)) +
    geom_point(alpha=0.1)+
    geom_smooth(se=F,method="gam", formula = "value~s(alpha, by =  sigma_common)+ sigma_common" ) +
    facet_grid(~name)


## model

library(mgcv)
dat<-dat %>%
    filter(N_min_alpha/N_obs >0.05)
gam_Sn<- gam(beta_Sn~s(alpha, by= sigma_common)+ sigma_common, method = "REML" ,data= dat)
summary(gam_Sn)

gam_C<- gam(beta_CC~s(alpha, by= sigma_common)+ sigma_common, method = "REML" ,data= dat)
summary(gam_C)
gam.check(gam_C)

gambeta_true<- gam(beta_true~s(alpha, by =  sigma_common)+ sigma_common, method = "REML" ,data= dat)
summary(gambeta_true)


x_new <- expand.grid(
    alpha=alpha,
    sigma_common = sigma_common)

Sn_pred <- as_tibble(predict(gam_Sn,  x_new,se.fit = T)) %>% rename(gam_Sn = fit, se_Sn=se.fit)
c_pred <- as_tibble(predict(gam_C,  x_new,se.fit = T)) %>% rename(gam_C = fit, se_C=se.fit)
beta_true_pred <- as_tibble(predict(gambeta_true,  x_new,se.fit = T))%>% rename(gambeta_true = fit, se_beta_true=se.fit)
x_new<-bind_cols(x_new, Sn_pred, c_pred,beta_true_pred)

x_new_long <- x_new %>%
    gather(key="index", value= "estimate", starts_with("gam")) %>%
    gather(key="index.se", value="se", starts_with("se") )


x_new_long<- x_new_long %>% mutate(index= recode(index,
                                                 gam_C= "beta[C]",
                                                 gam_Sn= "beta[S[n]]",
                                                 gambeta_true = "beta"),
                                   upper= estimate+ 1.96*se,
                                   lower= estimate- 1.96*se)

index_labels= c("beta[C]", "beta[S[n]]", "beta")
names(index_labels)<-c("gam_C", "gam_Sn", "gambeta_true")
x_new_long$index <- factor(x_new_long$index, levels = c("beta","beta[S[n]]", "beta[C]"))
x_new_long$sigma_common <- factor(x_new_long$sigma_common)

dat_long=  dat %>%
    pivot_longer(c(beta_CC, beta_Sn, beta_true),
                 names_to = "index",
                 values_to = "estimate" ) %>%
    mutate(index=
               recode_factor(index,
                                beta_true= "beta",
                                beta_Sn = "beta[S[n]]",
                                beta_CC = "beta[C]"
    ))

bottom_row <-  ggplot()+
    geom_point( aes(x=alpha,y=estimate, col= sigma_common), data = dat_long, alpha= 0.1)+
    geom_line( aes(x=alpha,y=estimate, col= sigma_common), data = x_new_long, size=1)+
    labs(y="beta-diversity")+
    facet_grid(.~index, labeller =label_parsed) +
    guides(col=guide_legend(nrow=2,byrow=TRUE, title="sigma"))+
    theme_cowplot(12)+
    theme(legend.position = "bottom")

save_plot("Simulations/FigureS1.jpg", plot = bottom_row,base_asp = 1,
          ncol = 3)




#saveRDS(dat, "Simulations/simulation_data_appendix.rds")
