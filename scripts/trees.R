library(tidyverse)
library(vegan)
library(betaC)
hf<-read.csv("data/hf253-03-trees-2014.csv")

hf  %>% select(gx, gy, quadrat) %>% distinct() %>%
    ggplot(aes(gx,gy, col=quadrat))+
    geom_point()

hf<-hf %>% filter(dbh>=10, census.id==1) %>% #dbh filter
    mutate(ha_plot= round(quadrat/100)) #

hf_m<-hf %>% group_by(ha_plot, sp) %>% summarise(count=n()) %>%
    pivot_wider(names_from = sp, values_from = count, values_fill =list(count=0) ) %>%
    column_to_rownames("ha_plot") %>% as.matrix()

data(BCI)
bci_m<-BCI[1:35,] %>% as.matrix()

#betaC
C_hf<-min(betaC:::beta_stand(hf_m, setsize = 10, list("C_target"),summarise = F,resamples = 1000))

C_bci<-min(betaC:::beta_stand(bci_m, setsize = 10, list("C_target"),summarise = F,resamples = 1000))

C_stand=min(C_hf, C_bci)-0.01
hf_betaC<-betaC:::beta_stand(hf_m, setsize = 10, list("beta_C"),args = list(C= C_stand),summarise = F)
bci_betaC<-betaC:::beta_stand(bci_m, setsize = 10, list("beta_C"),args = list(C= C_stand),summarise = F)
results_betaC=cbind(hf_betaC, bci_betaC )
colnames(results_betaC)=c("hf", "bci")
colMeans(results_betaC)

betaC_plot10<-results_betaC %>% as.data.frame() %>% pivot_longer(everything(), names_to = "Site", values_to = "beta_C") %>%
    ggplot(aes(x= Site, y= beta_C))+geom_boxplot()+labs(y=expression(beta[C]))


# beta Whittaker

hf_betaTrue<-betaC:::beta_stand(hf_m, setsize = 10, list("beta_true"),summarise = F)
bci_betaTrue<-betaC:::beta_stand(bci_m, setsize = 10, list("beta_true"),summarise = F)
results_betaTrue=cbind(hf_betaTrue,bci_betaTrue )
colnames(results_betaTrue)=c("hf", "bci")


betaTrue_plot10<-results_betaTrue %>% as.data.frame() %>% pivot_longer(everything(), names_to = "Site", values_to = "Whittaker") %>%
    ggplot(aes(x= Site, y= Whittaker))+geom_boxplot()+labs(y=expression(beta))

library(cowplot)
theme_set(theme_cowplot())
plot_grid(betaTrue_plot10 , betaC_plot10, labels = "AUTO")

# beta Sn

N=min(c(rowSums(bci_m),rowSums(hf_m)))


hf_betaSN<-betaC:::beta_stand(hf_m, setsize = 10, list("beta_SN"),args = list(N= N),summarise = F)
bci_betaSN<-betaC:::beta_stand(bci_m, setsize = 10, list("beta_SN"),args = list(N= N),summarise = F)
results_betaSN=cbind(hf_betaSN, bci_betaSN )
colnames(results_betaSN)=c("hf", "bci")
colMeans(results_betaSN)

betaSN_plot10<-results_betaSN %>% as.data.frame() %>% pivot_longer(everything(), names_to = "Site", values_to = "beta_SN") %>%
    ggplot(aes(x= Site, y= beta_SN))+geom_boxplot()+labs(y=expression(beta[S[n]]))


library(cowplot)
plot_grid(betaTrue_plot10, betaC_plot10, nrow = 1, labels="AUTO")
ggsave("bci_hf.jpg", width = 21, height = 10, units="cm")



#####
#extrapolation

#betaC
C_extra_hf<-min(betaC:::beta_stand(hf_m, setsize = 10, list("C_target_extra"),summarise = F,resamples = 1000))

C_extra_bci<-min(betaC:::beta_stand(bci_m, setsize = 10, list("C_target_extra"),summarise = F,resamples = 1000))

C_stand_extra=min(C_extra_hf, C_extra_bci)-0.01
hf_betaC_extra<-betaC:::beta_stand(hf_m, setsize = 10, list("beta_C_extra"),args = list(C= C_stand_extra),summarise = F)
bci_betaC_extra<-betaC:::beta_stand(bci_m, setsize = 10, list("beta_C_extra"),args = list(C= C_stand_extra),summarise = F)
results_betaC_extra=cbind(hf_betaC_extra, bci_betaC_extra )
colnames(results_betaC_extra)=c("hf", "bci")
colMeans(results_betaC_extra)



betaC_plot10_extra<-results_betaC_extra %>% as.data.frame() %>% pivot_longer(everything(), names_to = "Site", values_to = "beta_C_extra") %>%
    ggplot(aes(x= Site, y= beta_C_extra))+geom_boxplot()+labs(y=expression(beta[C]))


plot_grid(betaC_plot10+labs(title= "interpolated C = 0.86"),betaC_plot10_extra+labs(title= "extrapolated C = 0.93"))
