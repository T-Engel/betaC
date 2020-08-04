library(tidyverse)
library(vegan)
library(betaC)
library(cowplot)
theme_set(theme_cowplot())

######################################################
## First case study
# gentry dataset
library(BIEN)
#gentry<-BIEN_plot_dataset("Gentry Transect Dataset", all.metadata = T,only.new.world = T, all.taxonomy = T)
#saveRDS(gentry, "datasets/gentry.rds")
gentry<-readRDS("datasets/gentry.rds")

# Exclude occurences without individual counts
gentry<-subset(gentry, !is.na(individual_count))
# Capitalize coordinates
colnames(gentry)[match(c("latitude","longitude"),colnames(gentry))]<-c("Latitude","Longitude")

#extract plot level information
plots_gentry<-aggregate(gentry[,c("elevation_m","Latitude","Longitude")], by=list(plot_name=gentry$plot_name),function(x) return(mean(na.omit(x))) )

# Update coordinates in dataframe to make sure they are consistent for a plot
gentry[,"Latitude"]<-sapply(gentry$plot_name, function(x) return(plots_gentry[match(x,plots_gentry$plot_name),"Latitude"]),simplify = T)
gentry[,"Longitude"]<-sapply(gentry$plot_name, function(x) return(plots_gentry[match(x,plots_gentry$plot_name),"Longitude"]),simplify = T)

gentry=unite(gentry, "alpha_id", c("plot_name", "subplot"),remove = F)
gentry= gentry %>% mutate(abslat= abs(Latitude))

gentry_df<- gentry %>%
    group_by(alpha_id,plot_name, verbatim_scientific_name) %>%
    summarise(total=sum(individual_count)) %>%
    spread("verbatim_scientific_name", "total",fill= 0) %>%
    arrange(plot_name) %>%
    group_by(plot_name) %>%
    group_nest(.key = "SAD") %>%
    mutate(SAD = map(SAD, column_to_rownames,var= "alpha_id")) %>%
    mutate(SAD = map(SAD, function(x) x[,colSums(x)>0]))

plot_attr<- gentry %>% select(-verbatim_scientific_name, -individual_count)#
#determine variables to keep in mob_in
vars<- plot_attr %>% colnames()
index<-logical()
for( i in vars){
    index[i]=!any(rowSums(table(plot_attr$alpha_id, plot_attr[, i])>0)>1)
}


plot_attr<-plot_attr %>% select(names(index)[index]) %>% distinct
vars<- 1: ncol(plot_attr)
index<-logical()
for( i in vars){
    index[i]=plot_attr %>% group_by(plot_name) %>% select(i, plot_name) %>% n_distinct()>n_distinct(plot_attr$plot_name)

}
plot_attr<- plot_attr %>% group_by(plot_name) %>% nest(sample_attributes=vars[index])


gentry_df<-gentry_df %>% full_join(plot_attr,by = "plot_name")






gentry_df<- gentry_df %>% mutate(n_samples= map_dbl(SAD,nrow),
                                 n_min=map_dbl(SAD, function(x) min(rowSums(x))))
gentry_df<-gentry_df %>% filter(n_samples==10,
                                n_min > 5)

gentry_df<- gentry_df %>% mutate(beta_true= map_dbl(SAD, function(x) beta_true(x)))
gentry_df<- gentry_df %>% mutate(S_gamma= map_dbl(SAD, function(x) specnumber(colSums(x))))

gentry_df<- gentry_df %>% mutate(C_target= map_dbl(SAD, function(x) C_target(x)))
hist(gentry_df$C_target)



C_tar=min(gentry_df$C_target)
gentry_df<- gentry_df %>% mutate(betaC= map_dbl(SAD,beta_C, C=C_tar))
gentry_df=gentry_df %>% filter(!is.na(betaC))

beta_C_plot<-gentry_df %>%  ggplot(aes(abslat, betaC))+
    geom_point()+
    geom_smooth( col="black", method= "lm") +
    labs(x= "Absolute latitude", y= expression(beta[c]) )



beta_true_plot<-gentry_df  %>% ggplot(aes(abslat, beta_true))+
    geom_point()+
    geom_smooth( col="black", method= "lm") +
    labs(x= "Absolute latitude", y= expression(beta))

gentry_plot<-plot_grid(beta_true_plot , beta_C_plot, labels = c("A", "B"))

title <- ggdraw() +
    draw_label(
        expression('Case study 1: Gentry plots (C'[target]*' = '*0.10*")"),
        fontface = 'bold',
        x = 0,
        hjust = 0
    ) +
    theme(
        # add margin on the left of the drawing canvas,
        # so title is aligned with left edge of first plot
        plot.margin = margin(0, 0, 0, 7)
    )
first_row=plot_grid(
    title, gentry_plot,
    ncol = 1,
    # rel_heights values control vertical title margins
    rel_heights = c(0.1, 1)
)


######################################################
## Second case study
# HF dataset
HF<-read.csv("datasets/HF253-03-trees-2014.csv")

HF<-HF %>% filter(dbh>=10, census.id==1) %>% #dbh filter
    mutate(ha_plot= round(quadrat/100)) #

HF_m<-HF %>% group_by(ha_plot, sp) %>% summarise(count=n()) %>%
    pivot_wider(names_from = sp, values_from = count, values_fill =list(count=0) ) %>%
    column_to_rownames("ha_plot") %>% as.matrix()

#BCI dataset (fromn vegan package)
data(BCI)
BCI_m<-BCI[1:35,] %>% as.matrix()

# calc beta_C for subset of 10 plots (using beta_stand)
C_HF<-min(beta_stand(HF_m, setsize = 10, list("C_target"),summarise = F,resamples = 1000))

C_BCI<-min(beta_stand(BCI_m, setsize = 10, list("C_target"),summarise = F,resamples = 1000))

C_stand=round(min(C_HF, C_BCI)- 0.01, 2)
HF_betaC<-beta_stand(HF_m, setsize = 10, list("beta_C"),args = list(C= C_stand),summarise = F)
BCI_betaC<-beta_stand(BCI_m, setsize = 10, list("beta_C"),args = list(C= C_stand),summarise = F)
results_betaC=cbind(HF_betaC, BCI_betaC )
colnames(results_betaC)=c("HF", "BCI")
colMeans(results_betaC)

betaC_plot10<-results_betaC %>% as.data.frame() %>% pivot_longer(everything(), names_to = "Site", values_to = "beta_C") %>%
    ggplot(aes(x= Site, y= beta_C))+geom_boxplot()+labs(y=expression(beta[C]))


# beta Whittaker

HF_betaTrue<-betaC:::beta_stand(HF_m, setsize = 10, list("beta_true"),summarise = F)
BCI_betaTrue<-betaC:::beta_stand(BCI_m, setsize = 10, list("beta_true"),summarise = F)
results_betaTrue=cbind(HF_betaTrue,BCI_betaTrue )
colnames(results_betaTrue)=c("HF", "BCI")


betaTrue_plot10<-results_betaTrue %>% as.data.frame() %>% pivot_longer(everything(), names_to = "Site", values_to = "Whittaker") %>%
    ggplot(aes(x= Site, y= Whittaker))+geom_boxplot()+labs(y=expression(beta))

library(cowplot)
theme_set(theme_cowplot())
HF_BCI<-plot_grid(betaTrue_plot10 , betaC_plot10, labels = c("C", "D"))

title2 <- ggdraw() +
    draw_label(
        expression('Case study 2: BCI vs. HF (C'[target]*' = '*0.93*")"),
        fontface = 'bold',
        x = 0,
        hjust = 0
    ) +
    theme(
        # add margin on the left of the drawing canvas,
        # so title is aligned with left edge of first plot
        plot.margin = margin(0, 0, 0, 7)
    )
second_row=plot_grid(
    title2, HF_BCI,
    ncol = 1,
    # rel_heights values control vertical title margins
    rel_heights = c(0.1, 1)
)

# combine
fig_5<-plot_grid(first_row , second_row, nrow = 2)
ggsave("datasets/Figure5.jpg", fig_5, width = 21, height = 21, units="cm")



# coverage deficit BCI
coverage(colSums(BCI_m))- Chat(colSums(BCI_m),mean(rowSums(BCI_m)))

# coverage deficit HF
coverage(colSums(HF_m))- Chat(colSums(HF_m),mean(rowSums(HF_m)))

