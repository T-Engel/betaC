# The conceptual figures use random draws from SADs.
# Therefore they can look slightly different every time you run this script

##########################
# load packages and define some function

library(tidyverse)
library(vegan)
library(cowplot)
library(mobsim)
library(betaC)

rarefy_long <- function(x) {
    if(is.matrix(x)==F) x=matrix(x,nrow = 1, byrow =T, dimnames= list("x", names(x)))
    alphas <-
        lapply(row.names(x), function(i)
            return(as.numeric(vegan::rarefy(
                x[i, ], sample = 1:sum(x[i, ])
            )))) %>%
        lapply(function(x)
            return(data.frame(
                S_n = as.numeric(x), N = 1:length(x)
            )))
    names(alphas) <- rownames(x)
    alphas <- alphas %>% plyr::ldply(.id = "Curve")
    alphas$type = "minor"
    mean_alpha <-
        data.frame(
            Curve = "mean_alpha",
            S_n = colMeans(as.matrix(vegan::rarefy(
                x, 1:min(rowSums(x))
            ))),
            N = 1:min(rowSums(x)),
            type = "major"
        )
    gamma <-
        data.frame(
            Curve = "gamma",
            S_n = as.numeric(vegan::rarefy(colSums(x), 1:sum(x))),
            N = 1:sum(x),
            type = "major"
        )
    out = alphas %>% full_join(mean_alpha, by = c("Curve", "S_n", "N",  "type")) %>% full_join(gamma, by = c("Curve", "S_n", "N", "type"))
    return(out)

}
splitgamma <-
    function(x,
             type = c("distinct", "random", "on_curve"),
             n = round(sum(x) / 2),
             iter = 150) {
        if (type == "distinct") {
            alpha1 = x
            alpha2 = x

            #index1=sample(1:length(x),length(x)/2)
            #index2=setdiff(1:length(x),index1)
            #alpha1[index1] = 0
            #alpha2[index2] = 0

            # alpha1[seq(1, length(x), 2)] = 0
            # alpha2[seq(2, length(x), 2)] = 0
            alpha2[1]=0
            for (i in 2: length(x)){
                if(sum(alpha1) > sum(alpha2)){
                    alpha1[i]=0
                }else{
                    alpha2[i]=0
                }
            }
        }

        if (type == "random") {
            alpha1 = sample_sad_N(x,N =  n, replace = F)
            alpha2 = x - alpha1
        }

        if (type == "on_curve") {
            cases = lapply(1:iter, function(i)
                sample_sad_N(x = x, N = n))
            curves = lapply(cases, function(m)
                rarefy(m, 1:n))
            gamma_curve = rarefy(x, 1:n)
            SS = sapply(curves, function(Sn) {
                return(sum((gamma_curve - Sn) ^ 2))
            })
            alpha1 = cases[[order(SS)[1]]]
            alpha2 = x - alpha1
        }




        return(rbind(alpha1, alpha2))
    }

# Take subsamples of (Meta-)Community abundance vectors (individual based)

sample_sad_N<-function(x,N, replace=F){
    sample_fun<-function(x,N, replace){
        index=1:length(x)
        y=rep(index,x)
        samp<-data.frame(Species=sample(y, size = N, replace = replace)) %>%
            group_by_all() %>%
            count()
        missing=data.frame(Species=setdiff(index, samp$Species))
        samp=samp %>% full_join(missing, by = "Species") %>% arrange(Species) %>% pull(var = n)
        samp[is.na(samp)]<-0
        return(samp)
    }

    if(is.data.frame(x))x<- as.matrix(x)

    if(is.vector(x)){
        names=names(x)
        x<-matrix(x, byrow = T, nrow = 1)
    } else{
        names<- dimnames(x)
    }
    if(any(rowSums(x)==0)) stop("Remove sites without individuals!")
    out<-apply(x,1,sample_fun, replace = replace, N= N)
    out<-t(out)
    if(dim(out)[1]==1){
        out= out[1,,drop=T]
        names(out)= names
    }else{
        dimnames(out)<-names
    }
    return(out)
}
########################################################################################
# Styling
theme_set(theme_cowplot())
mytheme= theme(legend.position = "none",
               axis.text=element_text(size=8),
               axis.title=element_text(size=10),
               plot.title = element_text(size=10,face = "bold"))


#########################################################################

# Figure 1
# reference meta-community
# color palette
pal2<-viridisLite::magma(5)[c(1,4)]
base = as.integer(sim_sad(s_pool = 450, n_sim = 1000, sad_coef = list(cv_abund =2)) )
base_m = splitgamma(base, type = "on_curve",iter =300 )
base_curve <- rarefy_long(base_m)
base_plot <-
    base_curve %>% filter(type == "major") %>% ggplot(aes(N, S_n, col = Curve)) +
    geom_abline(intercept = specnumber(base), slope = 0, linetype=5, col=pal2[1])+
    geom_abline(intercept = mean(specnumber(base_m)), slope = 0, linetype=5, col=pal2[2])+
    geom_vline(xintercept = 250, linetype= "dashed", color ="grey")+
    geom_line(size = 1) +
    geom_text(aes(x= 1000, y= 55),col=1, label=paste0("beta == ", round(specnumber(base)/mean(specnumber(base_m)),2)), nudge_y = -30,parse = T ,hjust="right",vjust="center")+
    geom_text(aes(x= 600, y= 55),col=1, label="beta[S[n]] == 1", nudge_y = -30,parse = T ,hjust="right", vjust="center")+

    coord_cartesian(xlim = c(0, 1050), ylim = c(0, 300),expand = F)+
    labs(title = "reference", x= "Individuals", y= "Rarefied richness")+
    mytheme + scale_color_manual(values = pal2)

# fewer individuals
individuals = splitgamma(base, type = "on_curve", iter = 300,n = 500)[1, ]
individuals_m = splitgamma(individuals, type = "on_curve", iter = 200)
individuals_curve <- rarefy_long(individuals_m)
individuals_plot <-
    individuals_curve %>% filter(type == "major")  %>% ggplot(aes(N, S_n, col = Curve)) +
    geom_abline(intercept = specnumber(individuals), slope = 0, linetype=5, col=pal2[1])+
    geom_abline(intercept = mean(specnumber(individuals_m)), slope = 0, linetype=5,, col=pal2[2])+
    geom_vline(xintercept = 250, linetype= "dashed", color ="grey")+
    geom_line(data=base_curve %>% filter(type=="major", Curve=="gamma"), linetype= "dotted", size=1, col= "grey")  +
    geom_line(size = 1) +
    geom_text(aes(x= 1000, y= 55),col=1, label=paste0("beta == ", round(specnumber(individuals)/mean(specnumber(individuals_m)),2)), nudge_y = -30,parse = T ,hjust="right", vjust="center")+
    geom_text(aes(x= 600, y= 55),col=1, label="beta[S[n]] == 1", nudge_y = -30,parse = T ,hjust="right",vjust="center")+
    coord_cartesian(xlim = c(0, 1050), ylim = c(0, 300),expand = F) +
    labs(title = "fewer individuals", x= "Individuals", y= "Rarefied richness")+
    mytheme + scale_color_manual(values = pal2)

# SAD change
pool =   as.integer(sim_sad(s_pool = 80, n_sim = 1000, sad_coef = list(cv_abund = 2))  )# sim_ENS(30, 85, 1000)
pool_m = splitgamma(pool, type = "on_curve")
pool_curve <- rarefy_long(pool_m)
pool_plot <-
    pool_curve %>% filter(type == "major") %>% ggplot(aes(N, S_n, col = Curve)) +
    geom_abline(intercept = specnumber(pool), slope = 0, linetype=5, col=pal2[1])+
    geom_abline(intercept = mean(specnumber(pool_m)), slope = 0, linetype=5, col=pal2[2])+
    geom_vline(xintercept = 250, linetype= "dashed", color ="grey")+
    geom_line(size = 1) +
    geom_text(aes(x= 1000, y=55),col=1, label=paste0("beta == ", round(specnumber(pool)/mean(specnumber(pool_m)),2)), nudge_y = -30,parse = T ,hjust="right", vjust="center")+
    geom_text(aes(x= 600, y= 55),col=1, label="beta[S[n]] == 1", nudge_y = -30,parse = T ,hjust="right", vjust="center")+
    coord_cartesian(xlim = c(0, 1050), ylim = c(0, 300),expand = F) +
    labs(title = "smaller species pool", x= "Individuals", y= "Rarefied richness")+
    mytheme + scale_color_manual(values = pal2)

# aggregation
space_m = splitgamma(base, type = "distinct")
space_curve <- rarefy_long(space_m)
space_plot <-
    space_curve %>% filter(type == "major") %>% ggplot(aes(N, S_n, col = Curve)) +
    geom_abline(intercept = specnumber(base), slope = 0, linetype=5, col=pal2[1])+
    geom_abline(intercept = mean(rarefy(space_m, min(rowSums(space_m)))), slope = 0, linetype=5, col=pal2[2])+
    geom_vline(xintercept = 250, linetype= "dashed", color ="grey")+
    geom_line(size = 1) +
    geom_text(aes(x= 1000, y= 55),col=1, label=paste0("beta == ", round(specnumber(base)/mean(specnumber(space_m)),2)), nudge_y = -30,parse = T ,hjust="right", vjust="center")+
    geom_text(aes(x= 690, y= 55),col=1, label="beta[S[n]] == 1.38", nudge_y = -30,parse = T ,hjust="right",vjust="center")+

    labs(title = "intraspecific aggregation", x= "Individuals", y= "Rarefied richness")+
    coord_cartesian(xlim = c(0, 1050), ylim = c(0, 300),expand = F) +
    mytheme + scale_color_manual(values = pal2)


Figure1 <- plot_grid(
    NULL,
    base_plot,
    NULL,
    pool_plot,
    individuals_plot,
    space_plot,
    ncol = 3,
    labels = c(NA, "A", NA, "B", "C", "D"),
    align= "hv"
)
ggsave("conceptual figures/Figure1.jpg",Figure1, width = 18, height = 12, units="cm")

########################################################################################################

# Figures 2 and 3
library(mobsim)

# color palette
pal <- viridisLite::viridis(10)[c(1,8)]
names(pal) <- c("large", "small")

base = sim_sad(s_pool = 100, n_sim = 1000, sad_coef = list(cv_abund = 2))
space_m = splitgamma(base, type = "distinct")
space_curve <- rarefy_long(space_m)

pool2 = sim_sad(s_pool = 500, n_sim = 1000, sad_coef = list(cv_abund =2))
space2_m = splitgamma(pool2, type = "distinct")
space2_curve <- rarefy_long(space2_m)

N1<- min(rowSums(space_m))
gamma_Sn1<-rarefy(base,N1)
alpha_Sn1<- mean(rarefy(space_m,N1))

cov_value= Chat(pool2,min(rowSums(space2_m)))
cov_value_small= Chat(base,min(rowSums(space_m)))
N_low<-round(invChat(base, cov_value))
SnC_gamma <-D0.hat(base, N_low)
SnC_alpha <- mean(apply(space_m,1,D0.hat,m=N_low))
betaC = SnC_gamma/SnC_alpha

beta_C_small<-beta_C(space_m, cov_value)
beta_C_large<-beta_C(space2_m, cov_value)


space_curve$Curve<-relevel(space_curve$Curve, "gamma")
space2_curve$Curve<-relevel(space2_curve$Curve, "gamma")
small_plot_C <-
    ggplot() +
    geom_line(size = 1) +
    #geom_hline(yintercept = specnumber(base), linetype=5)+
    geom_hline(yintercept =SnC_gamma, linetype=5, col= "darkgrey")+
    geom_vline(xintercept = N_low, linetype = "dashed" , col= "darkgrey")+
    geom_hline(yintercept =SnC_alpha, linetype=5, col= "darkgrey")+
    geom_abline(slope = 1-cov_value, intercept = SnC_gamma - ((1-cov_value)*N_low), size=1, col= "darkgrey")+
    geom_line(aes(N, S_n, linetype = Curve), data= filter(space_curve,type == "major"), size = 1,col= pal[2]) +
    #geom_text(aes(x= N_low, y= 0), label=paste0("n = ", round(N_low,2)),nudge_x = 20, nudge_y = 10,parse = F ,hjust="left")+
    geom_text(aes(x= 1000, y= SnC_alpha+4), label=paste0("beta[C] == ", round(beta_C_small,3)),nudge_x = , nudge_y = -25,parse = T ,hjust="right", vjust="bottom")+
    labs(title="Small species pool\n(100 spp.)", x= "Individuals", y= "Rarefied richness")+ #"small pool - beta\nstandardised by coverage\nof large pool"
    coord_cartesian(xlim = c(0, 1050), ylim = c(0, 300),expand = F) +
    mytheme + theme(plot.title = element_text(colour = pal[2]))

small_plot_N <-
    ggplot(data =NULL) +
    geom_hline(yintercept =alpha_Sn1, linetype=5, col= "darkgrey")+
    geom_vline(xintercept = min(rowSums(space_m)), linetype = "dashed", col= "darkgrey")+
    geom_hline(yintercept = gamma_Sn1, linetype=5, col= "darkgrey")+
    geom_abline(slope = 1-cov_value_small, intercept = gamma_Sn1 - ((1-cov_value_small)*N1), size=1, col= "darkgrey")+
    geom_line(aes(N, S_n, linetype = Curve), data= filter(space_curve,type == "major"), size = 1, col= pal[2]) +
    labs(title="Small species pool\n(100 spp.)", x= "Individuals", y= "Rarefied richness")+
    #geom_text(aes(x= N1,y= 0), label=paste0("n = ", round(N1,2)),nudge_x = 20, nudge_y = 10, parse = F, hjust="left" )+
    geom_text(aes(x= 1000, y= alpha_Sn1), label=paste0("beta[s[n]] == ", round(gamma_Sn1/alpha_Sn1,3)),nudge_x = , nudge_y = -25,parse = T ,hjust="right", vjust="bottom")+

    coord_cartesian(xlim = c(0, 1050), ylim = c(0, 300), expand = F) +
    mytheme + theme(plot.title = element_text(colour = pal[2]))


N_low2 <-min(rowSums(space2_m))
gamma_Sn<-D0.hat(pool2,N_low2)
alpha_Sn<- mean(apply(space2_m,1, D0.hat,N_low2))
large_plot_N <-
    ggplot(data =NULL) +
    geom_hline(yintercept =alpha_Sn, linetype=5, col= "darkgrey")+
    geom_vline(xintercept = min(rowSums(space2_m)), linetype = "dashed", col= "darkgrey")+
    geom_hline(yintercept = gamma_Sn, linetype=5, col= "darkgrey")+
    geom_abline(slope = 1-cov_value, intercept = gamma_Sn - ((1-cov_value)*N_low2), size=1, col= "darkgrey")+
    geom_line(aes(N, S_n, linetype = Curve), data= filter(space2_curve,type == "major"), size = 1, col= pal[1]) +

    labs(title="Large species pool\n(500 spp.)", x= "Individuals", y= "Rarefied richness")+
    #geom_text(aes(x= N_low2,y= 0), label=paste0("n = ", round(N_low2,2)),nudge_x = 20, nudge_y = 10, parse = F, hjust="left" )+
    geom_text(aes(x= 1000, y= alpha_Sn), label=paste0("beta[s[n]] == ", round(gamma_Sn/alpha_Sn,3)), nudge_y = -25,parse = T ,hjust="right", vjust="bottom")+
    coord_cartesian(xlim = c(0, 1050), ylim = c(0, 300), expand = F) +
    mytheme + theme(plot.title = element_text(colour = pal[1]))

large_plot_C <-
    ggplot() +
    geom_hline(yintercept =alpha_Sn, linetype=5, col= "darkgrey")+
    geom_vline(xintercept = min(rowSums(space2_m)), linetype = "dashed", col= "darkgrey")+
    geom_hline(yintercept = gamma_Sn, linetype=5, col= "darkgrey")+
    geom_abline(slope = 1-cov_value, intercept = gamma_Sn - ((1-cov_value)*N_low2), size=1, col= "darkgrey")+
    geom_line(aes(N, S_n, linetype = Curve), data= filter(space2_curve,type == "major"), size = 1, col= pal[1]) +
    labs(title="Large species pool\n(500 spp.)", x= "Individuals", y= "Rarefied richness")+
    geom_text(aes(x= 1000, y= alpha_Sn), label=paste0("beta[C] == ", round(beta_C_large,3)), nudge_y = -25,parse = T ,hjust="right", vjust="bottom")+
    coord_cartesian(xlim = c(0, 1050), ylim = c(0, 300), expand = F) +
    mytheme +
    theme(plot.title = element_text(colour = pal[1]))



# scaling relationship

dat1=tibble(N=1:(2*min(rowSums(space_m))),
            C=map_dbl(N,function(N) Chat(colSums(space_m), N)),
            beta_Sn = map_dbl(N,function(N)beta_SN(space_m, N)),
            Species_pool= "small"
)

dat2=tibble(N=1:(2*min(rowSums(space2_m))),
            C=map_dbl(N,function(N) Chat(colSums(space2_m), N)),
            beta_Sn = map_dbl(N,function(N)beta_SN(space2_m, N)),
            Species_pool= "large"
)

dat=bind_rows(dat1,dat2)

N_plot<-ggplot(data=dat)+
    geom_line(aes(x=N, y=beta_Sn, col= Species_pool), size=1)+
    theme(legend.position = "bottom")+
    labs( x= "Individuals", y= expression(beta[s[n]]))+scale_color_manual(values = pal)+mytheme

C_plot<-ggplot(data=dat)+
    geom_line(aes(x=C, y=beta_Sn, col= Species_pool),size=1)+
    theme(legend.position = "bottom")+
    labs( x= "Estimated coverage", y= expression(beta[C]))+
    scale_color_manual(values = pal)+
    geom_vline(xintercept = cov_value, linetype = "dashed", col= "darkgrey")+
    mytheme


Figure2<-plot_grid(large_plot_N,small_plot_N, N_plot,ncol = 3, labels = "AUTO")
Figure3<-plot_grid(large_plot_C,small_plot_C, C_plot, ncol = 3, labels = "AUTO")

Figure2
Figure3
ggsave("conceptual figures/Figure2.jpg",Figure2, width = 18, height = 10, units="cm")
ggsave("conceptual figures/Figure3.jpg",Figure3, width = 18, height = 10, units="cm")
