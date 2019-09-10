library(betaC)
library(vegan)
library(tidyverse)
library(mobsim)

# simulate communities
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

low_C <- C_target(samples_low)
high_C<- C_target(samples_high)
C=min(low_C, high_C)

beta_C(samples_low,C)
beta_C(samples_high,C)

N= min(rowSums(samples_high), rowSums(samples_low))

beta_SN(samples_low, N)
beta_SN(samples_high, N)

specnumber(colSums(samples_low))/mean(specnumber(samples_low))
specnumber(colSums(samples_high))/mean(specnumber(samples_high))


high_curve<-beta_C_curve(samples_high)
low_curve<- beta_C_curve(samples_low)

high_curve$scenario= "high"
low_curve$scenario= "low"
both_curves= full_join(high_curve,low_curve)

both_curves %>%  ggplot(aes(N, beta_Sn, col= scenario)) + geom_line(size = 1) + geom_rug()

both_curves %>%  ggplot(aes(C, beta_Sn, col= scenario)) + geom_line(size = 1) + geom_rug()


