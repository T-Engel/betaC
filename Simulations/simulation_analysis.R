library(tidyverse)
library(furrr)
library(cowplot)
library(betaC)
theme_set(theme_cowplot())

dat<-readRDS( "Simulations/simulation_data.rds")

# Calculate metrics
# Richness
dat<- dat %>% mutate(S_gamma= map_dbl(alpha_samples, function(x)vegan::specnumber(colSums(x))))

# S_PIE
dat<- dat %>% mutate(S_PIE= map_dbl(total_com, calc_PIE, ENS =T))


# Whittaker's beta
dat<- dat %>% mutate(beta_true= map_dbl(alpha_samples, beta_true))

# remove cases with too few individuals
dat<- dat %>% mutate(N_min_alpha= map_dbl(alpha_samples, function(x) return(min(rowSums(x))))) # removing cases with empty samples
N_stand=50 ##
dat <- dat %>%
    filter(N_min_alpha>= (N_stand))

plan(multiprocess)

# beta_Sn
dat<-dat %>% mutate(beta_Sn= future_map_dbl(alpha_samples, beta_SN, N=N_stand, .progress=T))

# beta_C
dat<- dat %>% mutate(C_target= future_map_dbl(alpha_samples, C_target, .progress= T))
C= min(dat$C_target)

dat<-dat %>% mutate(beta_CC= future_map_dbl(alpha_samples, beta_C, C=C, .progress=T))

dat<- dat %>% mutate(motherpoints_common= as.factor(motherpoints_common))

beta_C_plot<-dat %>% ggplot( aes(alpha,beta_CC, group= motherpoints_common, col= motherpoints_common)) + geom_smooth()+
    labs(title= paste0("C = ", round(C, 3)))+ theme(legend.position = "none")

dat %>%
    ggplot( aes(s_pool,beta_CC, group= motherpoints_common, col= motherpoints_common)) +
    geom_point(alpha=0.2)+
    geom_smooth()+
    labs(title= paste0("C = ", round(C, 3)))+ theme(legend.position = "none")



##################
# Model and plot simulation results
library(mgcv)
gam_Sn<- gam(beta_Sn~s(s_pool, by= motherpoints_common)+ motherpoints_common, method = "REML" ,data= dat)
summary(gam_Sn)

gam_C<- gam(beta_CC~s(s_pool, by= motherpoints_common)+ motherpoints_common, method = "REML" ,data= dat)
summary(gam_C)
gam.check(gam_C)

gambeta_true<- gam(beta_true~s(s_pool, by =  motherpoints_common)+ motherpoints_common, method = "REML" ,data= dat)
summary(gambeta_true)

# new data for prediction
x_new <- expand.grid(s_pool=seq(10, max(dat$s_pool)), motherpoints_common =as.factor(c(1,4,10,20,max(dat$n_sim))))
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
  geom_point( aes(x=s_pool,y=estimate, col= motherpoints_common), data = dat_long, alpha= 0.05)+
  geom_line( aes(x=s_pool,y=estimate, col= motherpoints_common), data = x_new_long, size=1)+
    facet_grid(.~index, labeller =label_parsed) +
    labs(x ="Species pool size", y= expression(beta-"diversity"))+
    scale_colour_viridis_d(name  ="intraspecific aggregation\n(# conspecific clusters)",
                           breaks=as.factor(c(1,4,10,20, 4000)),
                           labels=c("1 cluster", "4 clusters", "10 clusters", "20 clusters", "random")) +

    guides(col=guide_legend(nrow=2,byrow=TRUE))+
    theme_cowplot(12)+
    theme(legend.position = "bottom")




# example
library(mobsim)
set.seed(1)
com<-sim_sad(100, 1000)
aggr_com<- sim_thomas_coords(com, sigma= 0.01, mother_points = 1)

x0 = 0.1
y0 = 0.1
quadrat_area= 0.08
quadrat_length= sqrt(quadrat_area)
delta_y = 0.5
delta_x = 0.5
xpos=c(x0,
       x0+  delta_x,
       x0,
       x0+  delta_x
)
ypos= c(y0,
        y0,
        y0+  delta_y,
        y0+  delta_y
)


p_aggr <- function() {
    par(
        mar = c(3, 3, 1, 1),
        mgp = c(2, 1, 0)
    )

    plot(aggr_com, main= "1 cluster", col=viridisLite::magma(100),xaxt="none",yaxt="none", xlab=NA, ylab=NA)
    graphics::rect(xpos, ypos, xpos + quadrat_length, ypos+
                       quadrat_length, lwd = 2, col = grDevices::adjustcolor("white",
                                                                             alpha.f = 0.6))
}

rand_com<- sim_thomas_coords(com, sigma= 0.02, mother_points = 4000)
p_rand <- function() {
    par(
        mar = c(3, 3, 1, 1),
        mgp = c(2, 1, 0)
    )

    plot(rand_com, main= "random", col=viridisLite::magma(100),xaxt="none",yaxt="none", xlab=NA, ylab=NA)
    graphics::rect(xpos, ypos, xpos + quadrat_length, ypos+
                       quadrat_length, lwd = 2, col = grDevices::adjustcolor("white",
                                                                             alpha.f = 0.6))
}
aggr_plot<-ggdraw(p_aggr)+coord_equal(xlim = c(0,1),ylim = c(0,1))
rand_plot<-ggdraw(p_rand)+coord_equal(xlim = c(0,1),ylim = c(0,1))
r1<-plot_grid(aggr_plot,rand_plot, labels = c("A", "B"))

fig4_new<-plot_grid(r1, bottom_row, ncol = 1, labels= c("", "C"))



ggsave("Figures/Figure4.pdf",fig4_new, width = 15.6, height = 15.6, units="cm")



