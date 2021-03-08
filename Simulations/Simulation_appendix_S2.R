#############
# Supplementary material S2: Beta-deviation
# Run "simulation_analysis" script first

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

##########################################################


# wrapper function
beta_func= function(x){
    assemble.from.pool.randA(x,rand.N = 400) %>%
        magrittr::extract2(2) %>%
        sapply( beta_true)
}

# Calculate Beta deviation
plan(multiprocess)
dat<- dat %>% mutate(
    null_betas= future_map(alpha_samples,beta_func, .progress = T)
)

dat <- dat %>%  mutate(
    null_mean= map_dbl(null_betas, mean),
    null_var= map_dbl(null_betas, var),
    beta_dev= (beta_true- null_mean)/sqrt(null_var)

)

# plot

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


save_plot( "Simulations/FigureS2.jpg",
           plot_grid(plot1,plot2,nrow = 1, ncol=2, labels = "AUTO"),ncol = 2 )


# Correlation
dat %>% ggplot(aes(beta_CC, beta_dev))+geom_point()
dat %>% select(beta_n, beta_dev) %>%
    filter(is.finite(beta_dev)) %>%
    cor(method = "spearman")

dat %>%
    ggplot(aes(beta_CC, beta_dev, col= motherpoints_common)) +
    geom_point()

#saveRDS(dat,"Simulations/simulation_data_beta_dev.rds")
