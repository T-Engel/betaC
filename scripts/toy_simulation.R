library(betaC)
library(vegan)
library(tidyverse)
library(mobsim)

# simulate communities (same aggregation, different species pool)
low=sim_thomas_community(s_pool = 50,n_sim = 2000,sigma = 0.1, mother_points = 2 )
high=sim_thomas_community(s_pool = 300,n_sim = 2000,sigma = 0.1, mother_points = 2 )

# sample
samples_low<-sample_quadrats(comm = low,quadrat_area = 0.1, avoid_overlap = T, n_quadrats = 4, method = "grid", delta_x = 0.6, delta_y = 0.6)[[1]]
samples_high<-sample_quadrats(comm = high,quadrat_area = 0.1, avoid_overlap = T, n_quadrats = 4, method = "grid", delta_x = 0.6, delta_y = 0.6)[[1]]

#rarefaction curves
curve_low<-betaC:::rarefy_long(as.matrix(samples_low))
curve_low %>% filter(type== "major") %>% ggplot(aes(N,S_n, col= Curve))+ geom_line(size=1)
curve_high<-betaC:::rarefy_long(as.matrix(samples_high))
curve_high %>% filter(type== "major") %>% ggplot(aes(N,S_n, col= Curve))+ geom_line(size=1)

# normal beta
specnumber(colSums(samples_low))/mean(specnumber(samples_low))
specnumber(colSums(samples_high))/mean(specnumber(samples_high))

# calculate beta_Sn
N= min(rowSums(samples_high), rowSums(samples_low))
beta_SN(samples_low, N)
beta_SN(samples_high, N)

# determine C
low_C <- C_target(samples_low)
high_C<- C_target(samples_high)
C=min(low_C, high_C)

# calculate beta_C
beta_C(samples_low,C)
beta_C(samples_high,C)


# calculate full curve of beta_C /beta_Sn
high_curve<-beta_C_curve(samples_high)
low_curve<- beta_C_curve(samples_low)
high_curve$scenario= "high"
low_curve$scenario= "low"
both_curves= full_join(high_curve,low_curve)

# plot
both_curves %>%  ggplot(aes(N, beta_Sn, col= scenario)) + geom_line(size = 1) + geom_rug()
both_curves %>%  ggplot(aes(C, beta_Sn, col= scenario)) + geom_line(size = 1) + geom_rug()

# repeat a couple of times



######


# Kraft null model

# beta_true= function(x, transformation= F){
#     alpha= mean(vegan::specnumber(x))
#     gamma=vegan::specnumber(colSums(x))
#     beta =gamma/alpha
#     if(transformation) beta= 1-(1/beta)
#     return(beta)
# }
#
# null_model<-function(x, permutations, func=beta_true,...){
#     x=as.matrix(x)
#     observed=match.fun(func)(x, ...)
#     null=sapply(1:permutations, function(i) {
#         n=rowSums(x)
#         gamma=colSums(x)
#         random<-tibble(species=sample(rep(1:length(gamma), gamma)),site=rep(1:length(n), n)) %>%
#             group_by(species, site) %>%
#             count %>%
#             spread(key=species, value= n, fill= 0) %>% column_to_rownames("site") %>%
#             as.matrix()
#         return(match.fun(func)(random,...))
#     })
#
#     mean_null=mean(null)
#     sd_null= sd(null)
#     beta_dev=(observed-mean_null)/sd_null
#     return(beta_dev)
# }
#
# beta_dev_high<-null_model(samples_high, 500)
#
# beta_dev_low<-null_model(samples_low, 500)
