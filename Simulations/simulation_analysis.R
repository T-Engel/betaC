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

ggsave("Simulations/Figure4.jpg",fig4_new, width = 21, height = 21, units="cm")

##
# Supplementary material


# beta_C_extra

# factors
factor= c(1,2,3,4,5,10,20, 100)

plotlist=list()
C_targ= numeric()
results=dat %>% select(s_pool, motherpoints_common)

for(i in 1:length(factor)){
    C_targ[i] <- min(future_map_dbl(dat$alpha_samples, C_target,factor=factor[i], .progress= T))
    results[,(i+2)]<-future_map_dbl(dat$alpha_samples, beta_C, C=C_targ[i],interrupt=F, .progress=T)
}
colnames(results)[3:(length(factor)+2)]= paste0(
    "C_target = ", round(C_targ,3),", Factor: ", factor
)

factor_plot<-results %>%
    pivot_longer(-c(s_pool, motherpoints_common),names_to = "Factor",values_to = "beta_C") %>%
    ggplot( aes(s_pool,beta_C, group= motherpoints_common, col= motherpoints_common))+
    geom_smooth()+
    scale_colour_viridis_d(name  ="intraspecific aggregation\n(# conspecific clusters)",
                           breaks=as.factor(c(1,4,10,20, 4000)),
                           labels=c("1 cluster", "4 clusters", "10 clusters", "20 clusters", "random")) +
    facet_wrap(~Factor,ncol = 2)+
    theme(legend.position = "bottom")


ggsave("Simulations/Supplement.jpg", factor_plot, width = 21, height = 27, units="cm")

# beta asymptotic


C_asymptotic=0.999999

dat<-dat %>% mutate(beta_CC_asymptotic= future_map_dbl(alpha_samples, beta_C_extra, C=C_asymptotic,interrupt=F, .progress=T))


beta_C_asymptotic_plot<-dat %>% ggplot( aes(s_pool,beta_CC_asymptotic, group= motherpoints_common, col= motherpoints_common)) + geom_smooth()+
    labs(title= "Asymptotic")+theme(legend.position = "none")

plot_grid(beta_C_plot, beta_C_extra_plot, beta_C_asymptotic_plot, nrow = 1)


#############
# beta_deviation


# Null model code from http://jsebastiantello.weebly.com/r-code.html

### FUNCTION CODE ##############################################################
assemble.from.pool.randA <- function(compo, rand.N=999, fix.local.abund=TRUE,
                                     fix.rSAD=TRUE, pool.type="data.defined", spp.pool=NULL, save.output=FALSE,
                                     save.format="matrices", path.to.save, print.progress=TRUE)
{
    ## Makes sure there is a path to save files if needed
    if(save.output==TRUE & missing(path.to.save)) path.to.save <- getwd()

    ## Makes sure 'compo' is in the right format
    original.dim.compo <- length(dim(compo))
    if(original.dim.compo==0)
        stop("'compo' does not have more than 1 dimension")
    if(original.dim.compo>2)
        stop("This function does not know what to do when 'compo' has more than 2
      dimensions")

    compo <- as.matrix(compo, ncol=ncol(compo))

    ## If missing, creates names for the composition matrix
    if(is.null(rownames(compo)))
        rownames(compo) <- paste("site_", 1:nrow(compo), sep="")
    if(is.null(colnames(compo)))
        colnames(compo) <- paste("sp_", 1:nrow(compo), sep="")


    ## Checks that values in the cells are all positive integers
    if(sum((compo - round(compo, 0)) != 0) > 0)
        stop("This function cannot randomize an array with species abundances that
      are not integers")
    if(sum(compo<0) > 0)
        stop("This function cannot randomize an array with negative species
      abundances")
    if(sum(is.na(compo)) > 0)
        stop("This function cannot randomize an array with NAs")


    ## Calculates abundances per species in the species pool
    if(pool.type=="user.defined")
    {
        if(length(spp.pool)<ncol(compo))
            stop("Species pool is smaller than the number of species in 'compo'")

        if(is.null(names(spp.pool)))
            names(spp.pool) <- paste("sp_", 1:length(spp.pool), sep="")

        if(max(table(names(spp.pool)))>1)
            stop("Species pool contains repeated species names")

        if(length(setdiff(colnames(compo), names(spp.pool)))>0)
            warning("At least some species in 'compo' are not part of the species pool")

        spp.in.pool <- names(spp.pool)
        pool.spp.abund <- spp.pool
    }

    if(pool.type=="data.defined")
    {
        spp.in.pool <- colnames(compo)
        pool.spp.abund <- colSums(compo)
    }

    if(min(pool.spp.abund)<=0)
        warning("Some species in the species pool have abundances of zero")


    ## Calculates total number of occurrences
    pool.abundance <- sum(pool.spp.abund)
    matrix.abundance <- sum(compo)
    if(pool.abundance < matrix.abundance & pool.type=="user.defined")
        stop("if 'pool.type=user.defined', then the number of individual in
      'spp.pool' cannot be less than the number of individuals in 'compo'")


    ## Calculates occurrences per site
    site.densities <- rowSums(compo)
    if(min(site.densities)<=0)
        warning("Some empirical sites seem to be empty of individuals")


    ## Finds the species with at least one individual
    spp.more.than.0.indices <- which(pool.spp.abund>0)
    spp.more.than.0.names <- spp.in.pool[spp.more.than.0.indices]


    ## Defines the site of occurrences in the dataset. Randomized species IDs will
    ## be assigned to each. This maintains the number of occurrences per site and
    ## per census constant in the null model.
    numbers.per.cell <- as.numeric((apply(compo, 1, function(x) as.numeric((x)))))
    emp.site.assignation <- rep(x=rep(rownames(compo), each=ncol(compo)),
                                times=numbers.per.cell)

    if(pool.type=="data.defined")
        emp.species.assignation <- rep(rep(spp.in.pool, times=nrow(compo)),
                                       times=numbers.per.cell)
    if(pool.type=="user.defined")
        emp.species.assignation <- rep(spp.in.pool, times=pool.spp.abund)


    ## Produces an empty list to save null composition matrices
    rand.datasets <- sapply(rep(NA, rand.N), list)
    names(rand.datasets) <- paste("RandDataset_", 1:rand.N, sep="")

    require(R.utils)
    pb <- txtProgressBar(min = 0, max = rand.N, style = 3)
    for(i in 1:rand.N)
    {
        if(print.progress==TRUE)
            setTxtProgressBar(pb, i)

        ## When the regional SAD is fixed and the pool is defined by the data, makes
        ## the null SAD identical to the empirical SAD
        if(fix.rSAD==TRUE & pool.type=="data.defined")
            null.species.assignation <- emp.species.assignation

        ## When the regional SAD is fixed and the pool is defined by the user, makes
        ## the null SAD a sample of 'matrix.abundance' individuals from the species
        ## pool
        if(fix.rSAD==TRUE & pool.type=="user.defined")
            null.species.assignation <- sample(emp.species.assignation,
                                               matrix.abundance)

        ## When the regional SAD is NOT fixed, produces a null SAD to be used in
        ## analyses by randomly assigning individuals to species
        if(fix.rSAD==FALSE)
            null.species.assignation <- sample(spp.in.pool, matrix.abundance,
                                               replace=TRUE)


        ## Assigns individuals to local sites at random
        if(fix.local.abund==TRUE)
            null.site.assignation <- sample(emp.site.assignation, replace=FALSE)

        if(fix.local.abund==FALSE)
            null.site.assignation <- sample(rownames(compo), matrix.abundance,
                                            replace=TRUE)


        ## Creates a null species composition matrix using the null assignation of
        ## individuals to sites
        null.compo <- table(null.site.assignation, null.species.assignation)


        ## Matches names in rows of the null and empirical composition matrices.
        ## This also matches names of the columns in the null composition matrix
        ## with the species in the species pool, creating empty columns (NAs) for
        ## species not sampled during the randomization
        null.empty.spp <- setdiff(spp.in.pool, colnames(null.compo))
        if(length(null.empty.spp)>0)
            warning("Some species with zero abundances have been produced")

        if(nrow(compo)>1)
        {
            null.compo <- null.compo[match(rownames(compo), rownames(null.compo)),]
            rownames(null.compo) <- rownames(compo)
        }

        if(ncol(compo)>1)
        {
            null.compo <- null.compo[,match(spp.in.pool, colnames(null.compo))]
            colnames(null.compo) <- spp.in.pool
        }


        ## Replaces NA values with 0. NAs can be generated because of species in the
        ## species pool that were not sampled during the randomization
        which.na.values <- which(is.na(null.compo))
        if(length(which.na.values)>0)
            null.compo[which.na.values] <- 0


        ## Adds the null composition matrix to the list that compiles the results
        if(save.output==FALSE | save.format=="list")
            rand.datasets[[i]] <- null.compo

        ## If requested, the null composition matrix is saved as a file
        if(save.output==TRUE & save.format=="matrices")
            write.table(null.compo, file=paste(path.to.save, "\\", "RandDataset_", i,
                                               ".txt", sep=""), quote=FALSE, sep="\t", na="NA", dec=".",
                        row.names=TRUE, col.names=TRUE)
    }
    close(pb)


    ## If saving a list is requested, the full list of null composition matrices
    ## is saved as a file
    if(save.format=="list")
        save(rand.datasets, file=paste(path.to.save, "\\", "RandDatasets", sep=""))

    ## Makes a list of the parameters used in the randomization
    rand.parameters <- c(fix.local.abund, fix.rSAD, rand.N)
    names(rand.parameters) <- c("fix.local.abund", "fix.rSAD", "rand.N")

    ## Creates the output to return, depending on whether the main output was
    ## saved or not
    if(save.output==TRUE)
        output <- list(rand.parameters, path.to.save)
    else
        output <- list(rand.parameters, rand.datasets)

    names(output) <- c("rand.parameters", "rand.datasets")

    output
}



beta_func= function(x){
    assemble.from.pool.randA(x,rand.N = 400) %>%
        magrittr::extract2(2) %>%
        sapply( beta_true)
}

plan(multiprocess)
dat<- dat %>% mutate(
    null_betas= future_map(alpha_samples,beta_func, .progress = T)
)



dat <- dat %>%  mutate(
    null_mean= map_dbl(null_betas, mean),
    null_var= map_dbl(null_betas, var),
    beta_dev= (beta_true- null_mean)/sqrt(null_var)

)

dat %>% ggplot( aes(s_pool,beta_dev, group= motherpoints_common, col= motherpoints_common)) +
    geom_point(alpha=0.1)+
    geom_smooth()

plot1<-
dat %>% ggplot( aes(s_pool,beta_dev, group= motherpoints_common, col= motherpoints_common)) +
    geom_point(alpha=0.05)+
    geom_smooth()+
    labs(y=expression(beta["dev"]), x="Species pool size")+
    scale_colour_viridis_d(name  ="intraspecific aggregation\n(# conspecific clusters)",
                           breaks=as.factor(c(1,4,10,20, 4000)),
                           labels=c("1 cluster", "4 clusters", "10 clusters", "20 clusters", "random"))

plot2<-
dat %>% ggplot( aes(s_pool,beta_CC, group= motherpoints_common, col= motherpoints_common)) +
    geom_point(alpha=0.05)+
    geom_smooth()+
    labs(y=expression(beta["C"]), x="Species pool size")+
    scale_colour_viridis_d(name  ="intraspecific aggregation\n(# conspecific clusters)",
                           breaks=as.factor(c(1,4,10,20, 4000)),
                           labels=c("1 cluster", "4 clusters", "10 clusters", "20 clusters", "random"))


save_plot( "Simulations/supplement/beta_dev_plot.jpg",
           plot_grid(plot1,plot2,nrow = 1, ncol=2, labels = "AUTO"),ncol = 2 )
dat %>% ggplot(aes(beta_CC, beta_dev))+geom_point()
dat %>% select(beta_n, beta_dev) %>%
    filter(is.finite(beta_dev)) %>%
    cor(method = "spearman")

dat %>%
    ggplot(aes(beta_CC, beta_dev, col= motherpoints_common)) +
    geom_point()

saveRDS(dat,"Simulations/simulation_data_beta_dev.rds")
