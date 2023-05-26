#######################################################################################
###################### Code for Hagen & Mason 2023 meta-analysis ######################
#######################################################################################

#Load in necessary packages
library(metafor)
library(metagear)
library(dplyr)
library(ape)
library(ggplot2)
library(ggpubr)
library(emmeans)
library(xlsx)
library(orchaRd)
library(rlang)

#Load in effect sizes & phylogeny
ES_all <- read.xlsx("/Path/to/data/file", 1)
tree <- read.tree("/Path/to/tree/file")

#Change columns to numeric
ES_all[,c(22:27)] <- sapply(ES_all[,c(22:27)], as.numeric)

#Separate data
ES_in_search <- ES_all[which(ES_all$Search_group == 0),]
ES_not_in_search <- ES_all[which(ES_all$Search_group == 1),]

#Impute missing standard deviation values
impute_sd <- function(dat){
  sd_sum <- sum(na.omit(dat$sd1i))
  mean_sum <- sum(na.omit(dat$m1i))
  for(i in which(is.na(dat$sd1i))){
    dat$sd1i[i] <- (dat$m1i[i] * (sd_sum / mean_sum))
  }
  sd_sum <- sum(na.omit(dat$sd2i))
  mean_sum <- sum(na.omit(dat$m2i))
  for(i in which(is.na(dat$sd2i))){
    dat$sd2i[i] <- (dat$m2i[i] * (sd_sum / mean_sum))
  }
  dat$n1i <- as.numeric(dat$n1i) ; dat$n2i <- as.numeric(dat$n2i)
  return(dat)
}
ES_all <- impute_sd(ES_all)
ES_in_search <- impute_sd(ES_in_search)
ES_not_in_search <- impute_sd(ES_not_in_search)

#Calculate SMDs & add them to the datasets
get_smd <- function(dat){
  SMD <- escalc(measure = "SMD", n1i = dat$n1i, n2i = dat$n2i, m1i = dat$m1i, m2i = dat$m2i, sd1i = dat$sd1i, sd2i = dat$sd2i)
  dat <- bind_cols(dat, SMD)
  for(i in 1:nrow(dat)){
    if(is.na(dat$Effect_Direction[i])){
      next
    }else if(dat$Effect_Direction[i] == 1){
      dat$yi[i] <- dat$yi[i] * -1
    }else{ next }
  }
  return(dat)
}
ES_all <- get_smd(ES_all)
ES_in_search <- get_smd(ES_in_search)
ES_not_in_search <- get_smd(ES_not_in_search)

#Labeling NA values in descriptive columns as "unknown" to search for patterns
convert_na <- function(dat){
  dat$Polyploid_type[which(is.na(dat$Polyploid_type))] <- "Unknown"
  dat$Auto.allo[which(is.na(dat$Auto.allo))] <- "Unknown"
  dat$Cultivation_status[which(is.na(dat$Cultivation_status))] <- "Unknown"
  dat$Synth.est[which(is.na(dat$Synth.est))] <- "Unknown"
  dat$Pathogen_type[which(is.na(dat$Pathogen_type))] <- "Unknown"
  dat$Pathogen_lifestyle[which(is.na(dat$Pathogen_lifestyle))] <- "Unknown"
  return(dat)
}
ES_all <- convert_na(ES_all)
ES_in_search <- convert_na(ES_in_search)
ES_not_in_search <- convert_na(ES_not_in_search)

#Accounting for phylogeny in a multi-level meta-analytic model:
get_phylo_model <- function(dat){
  specs <- unique(dat$Family)
  phylo <- drop.tip(tree, tree$tip.label[which(!tree$tip.label %in% specs)])
  specs <- which(specs %in% phylo$tip.label)
  phy <- compute.brlen(phylo)
  cor <- vcv(phy, cor = T)
  comput_dat <- dat[which(dat$Family %in% phylo$tip.label),]
  names(comput_dat)[which(names(comput_dat) == "Family")] <- "phylo"
  phylo_m <- rma.mv(yi = yi, V = vi, W=1/vi, mods = ~ Polyploid_type + Cultivation_status + 
             Synth.est + Pathogen_type + Pathogen_lifestyle + Auto.allo, random = list(~1 | phylo, 
             ~1 | Study2, ~1 | ES2, ~1 | Diff_dips1), R = list(phylo = cor), method = "REML", 
             data = comput_dat)
  return(phylo_m)
}
ES_all_phymod <- get_phylo_model(ES_all)
summary(ES_all_phymod)
i2_ml(ES_all_phymod) #No effect of phylogeny

#Robust multilevel model
robust_m <- rma.mv(yi = yi, V = vi, W = 1/vi, mods = ~ Polyploid_type + Cultivation_status + Synth.est +
                   + Pathogen_type + Pathogen_lifestyle + Auto.allo, random = list(~1 | Study2, 
                   ~1 | ES2, ~1 | Diff_dips1), method = "REML", data = ES_all)
summary(robust_m)
i2_ml(robust_m)

#Extract mean effect size
mean(robust_m$b[,1]) #Effect size

#Removing largest families
pruned_down <- ES_all[-which(ES_all$Family == "Musaceae"),]
pruned_down <- pruned_down[-which(ES_all$Family == "Poaceae"),]
pruned_mod <- rma.mv(yi = yi, V = vi, W = 1/vi, mods = ~ factor(Diff_dips1) + Family,
                     random = list(~1 | Study2, ~1 | ES2), method = "REML", data = pruned_down)
summary(pruned_mod)

#Modeling just auto.allo
for_mod <- ES_all[-which(ES_all$Auto.allo == "Both"),]
for_mod <- for_mod[-which(for_mod$Auto.allo == "NA"),]
for_mod <- for_mod[,c(2,4,5,14,30,31)]
aa_mod <- rma.mv(yi = yi, V = vi, W = 1/vi, mods = ~ factor(Diff_dips1) + Auto.allo,
                 random = list(~1 | Study2, ~1 | ES2), method = "REML", data = for_mod)
summary(aa_mod)
i2_ml(aa_mod)
#Simple rma model; removing intercept so all moderators appear
aa_simple <- rma(yi = yi, vi = vi, mods = ~ Auto.allo - 1, method = "REML", data = for_mod)
summary(aa_simple)

#Modeling just pathogen lifestyle
pl_mod <- rma.mv(yi = yi, V = vi, W = 1/vi, mods = ~ factor(Diff_dips1) + Pathogen_lifestyle,
                 random = list(~1 | Study2, ~1 | ES2), method = "REML", data = ES_all)
summary(pl_mod)

#Egger's regression + influence of publication year (both of these require simpler rma models):
eggers <- rma(yi = yi, vi = vi, method = "REML", data = ES_all)
regtest(eggers, model='rma')
metareg.pubyear.all <- rma(yi, vi, mod= ~Year, method="REML", data=ES_all)
summary(metareg.pubyear.all)

#Fail-safe N:
#### Rosenthal fail-safe N
fsn(ES_all$yi, ES_all$vi, type="Rosenthal") #3408
fsn(ES_in_search$yi, ES_in_search$vi, type="Rosenthal") #666
fsn(ES_not_in_search$yi, ES_not_in_search$vi, type="Rosenthal") #9433
#### Orwin fail-safe N
fsn(ES_all$yi, ES_all$vi, type="Orwin", target=0.05) #0
fsn(ES_in_search$yi, ES_in_search$vi, type="Orwin", target=0.05) #237
fsn(ES_not_in_search$yi, ES_not_in_search$vi, type="Orwin", target=0.05) #172
#### Rosenberg fail-safe N
fsn(ES_all$yi, ES_all$vi, type="Rosenberg") #33141
fsn(ES_in_search$yi, ES_in_search$vi, type="Rosenberg") #194
fsn(ES_not_in_search$yi, ES_not_in_search$vi, type="Rosenberg") #24038


#######################################################################################
####################################### FIGURES #######################################
#######################################################################################

#Set the working directory for saving figures
setwd("/path/to/figure/directory/")

#Figure 1: PRISMA Diagram
phases <- c("START_PHASE: # of studies identified through database searching\n\n1602",
            "START_PHASE: # of additional studies identified through other sources\n\n55",
            "# of studies after duplicates removed\n\n1545",
            "# of studies with title and abstract screened\n\n1545",
            "EXCLUDE_PHASE: # of studies excluded\n\n1390",
            "# of full-text articles assessed for eligibility\n\n155",
            "EXCLUDE_PHASE: # of full-text excluded, not fitting eligibility criteria\n\n27",
            "final # of studies included in quantitative synthesis (meta-analysis)\n\n128")
jpeg("Fig1.jpeg", width = 10, height = 10, units = 'in', res = 300)
plot_PRISMA(phases, design="classic")
dev.off()

#Figure 2: Orchard plot of families
for_mod <- ES_all[-which(ES_all$Family == "Actinidiaceae"),] #Plot won't work with Actinidiaceae
pf_mod <- rma.mv(yi = yi, V = vi, W = 1/vi, mods = ~ factor(Diff_dips1) + Family,
                 random = list(~1 | Study2, ~1 | ES2), method = "REML", data = for_mod)
summary(pf_mod)
i2_ml(pf_mod)
g <- orchard_plot(pf_mod, mod="Family", data=for_mod, group="Study2", xlab="Standardized mean difference", transfm="none")
jpeg("Fig2.jpeg", width = 10, height = 10, units = 'in', res = 300)
g + theme(axis.text.x=element_text(size=8), axis.text.y=element_text(angle=0, size=8),
          axis.title.x=element_text(size=10), axis.title.y=element_text(size=10))
dev.off()

#Figure 3: Combined orchard plot of auto/allo and pathogen lifestyle
newcol <- paste0(ES_all$Auto.allo, "_", ES_all$Pathogen_lifestyle)
w <- 1/ES_all$vi
orch_dat <- as.data.frame(cbind(ES_all$ES2, ES_all$Study2, ES_all$Diff_dips1, newcol, ES_all$yi, ES_all$vi, w))
colnames(orch_dat) <- c("ES2", "Study2", "Diff_dips1", "Combo", "yi", "vi", "w")
greenlit <- c("Auto_Hemibiotrophic", "Allo_Necrotrophic", "Auto_Biotrophic", "Allo_Biotrophic", "Allo_Hemibiotrophic", "Auto_Necrotrophic")
orch_dat <- orch_dat[which(orch_dat$Combo %in% greenlit),]
combo_modfit <- rma.mv(yi = as.numeric(yi), V = as.numeric(vi), W = as.numeric(w), 
                  mods = ~ Combo, random = list(~1 | as.numeric(Study2), ~1 | as.numeric(ES2), 
                  ~1 | as.numeric(Diff_dips1)), method = "REML", data = orch_dat)
summary(combo_modfit)
i2_ml(combo_modfit)
g <- orchard_plot(combo_modfit, mod="Combo", data=orch_dat, group="Study2", xlab="Standardized mean difference", transfm="none")
jpeg("Fig3.jpeg", width = 7, height = 7, units = 'in', res = 300)
g + theme(axis.text.x=element_text(size=8), axis.text.y=element_text(angle=0, size=8),
          axis.title.x=element_text(size=10), axis.title.y=element_text(size=10))
dev.off()

#Figure 4: Combined orchard plot of auto/allo and synth/est
newcol <- paste0(ES_all$Auto.allo, "_", ES_all$Synth.est)
w <- 1/ES_all$vi
orch_dat <- as.data.frame(cbind(ES_all$ES2, ES_all$Study2, ES_all$Diff_dips1, newcol, ES_all$yi, ES_all$vi, w))
colnames(orch_dat) <- c("ES2", "Study2", "Diff_dips1", "Combo", "yi", "vi", "w")
allowed <- c("Auto_Unsure", "Allo_Unsure", "Both_Unsure", "Both_Synth", "Both_Est", "NA_Unsure", "NA_Synth", "NA_Est")
orch_dat$Combo[which(orch_dat$Combo %in% allowed)] <- "Other"
combo_simple <- rma(yi = as.numeric(yi), vi = as.numeric(vi), mods = ~ Combo - 1, method = "REML", data = orch_dat)
summary(combo_simple)
g <- orchard_plot(combo_simple, mod="Combo", data=orch_dat, group="Study2", xlab="Standardized mean difference", transfm="none")
jpeg("Fig4.jpeg", width = 7, height = 7, units = 'in', res = 300)
g + theme(axis.text.x=element_text(size=8), axis.text.y=element_text(angle=0, size=8),
          axis.title.x=element_text(size=10), axis.title.y=element_text(size=10))
dev.off()
#Just synth vs. unsure
synth_simple <- rma(yi = as.numeric(yi), vi = as.numeric(vi), mods = ~ Synth.est - 1, method = "REML", data = ES_all)
summary(synth_simple)

#Figure 5: Orchard plot of pathogen type
pt_mod <- rma.mv(yi = yi, V = vi, W = 1/vi, mods = ~ factor(Diff_dips1) + Pathogen_type,
                 random = list(~1 | Study2, ~1 | ES2), method = "REML", data = ES_all)
summary(pt_mod)
i2_ml(pt_mod)
g <- orchard_plot(pt_mod, mod="Pathogen_type", data=ES_all, group="Study2", xlab="Standardized mean difference", transfm="none")
jpeg("Fig5.jpeg", width = 7, height = 7, units = 'in', res = 300)
g + theme(axis.text.x=element_text(size=8), axis.text.y=element_text(angle=0, size=8),
          axis.title.x=element_text(size=10), axis.title.y=element_text(size=10))
dev.off()

#######################################################################################
################################ SUPPLEMENTARY FIGURES ################################
#######################################################################################

#Figure A1: Histogram of publication over time + pie charts
tb <- ES_all$Year
tb <- tb[order(tb)]
tbh <- data.frame(year=as.vector(tb))
p1 <- ggplot(tbh, aes(x=year)) +
  geom_histogram(color="black", fill="gray", binwidth=8, center=1) +
  scale_x_continuous(breaks=seq(1957,2021,8)) +
  scale_y_continuous(breaks=seq(0,70,10), lim=c(0,70)) +
  labs(y = "Number of Publications", x = "Year Published") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title=element_text(size=10), axis.text=element_text(size=8))
#Pie chart of pathogen types
tbp <- table(ES_all$Pathogen_type)
tbp <- tbp[order(tbp)]
tbp <- data.frame(type=names(tbp), num=as.vector(tbp))
p2 <- ggplot(tbp, aes(x="", y=num, fill=type)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() +
  theme(legend.position = 'left', legend.title=element_blank())
#Pie chart of plant families
tbf <- table(ES_all$Family)
tbf <- tbf[order(tbf)]
tbf <- data.frame(family=names(tbf)[14:21], num=as.vector(tbf)[14:21])
tbf_fixed <- data.frame(family="3 or fewer effect sizes", num=19)
tbf_fixed <- rbind(tbf_fixed, tbf)
tbf_fixed <- tbf_fixed[order(tbf_fixed$num),]
p3 <- ggplot(tbf_fixed, aes(x="", y=num, fill=family)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() +
  theme(legend.title=element_blank())
#Plot all 3 at once
jpeg("FigA1.jpeg", width = 9, height = 9, units = 'in', res = 300)
ggpubr::ggarrange(p1, ggarrange(p2, p3, ncol = 2, labels = c("B", "C"), font.label=list(size = 18, color = "black")), nrow = 2, labels = "A", font.label=list(size = 18, color = "black")) 
dev.off()

#Figure A2: Caterpillars plot
make_caterpillar <- function(model, modz=NULL, dat){
  I2 <- i2_ml(model)
  disp <- round(I2[1],1)
  if(is.null(modz)){
    model_results <- mod_results(model, mod="1", group="ES1", data=dat) 
    #1 denotes intercept (the default)
  }else{
    model_results <- mod_results(model, mod = modz, group="ES1", data=dat)
  }
  p <- caterpillars(model_results, mod = "Int", xlab = "Standardized mean difference") 
  q <- p + theme(axis.text.x=element_text(size=8), axis.title.x=element_text(size=10)) + 
       annotate(geom = "text", x = -3.5, y = 20, label = paste0("italic(I)^{2} == ", 
       disp), color = "black", parse = TRUE, size = 8)
  q + annotate(geom = "text", x = -2.35, y = 19, label = "%", color = "black", parse = FALSE, 
               size = 8)
}
jpeg("FigA2.jpeg", width = 9, height = 9, units = 'in', res = 300)
make_caterpillar(robust_m, dat=ES_all)
dev.off()

#Figure A3: Orchard plot of auto vs. allo
g <- orchard_plot(aa_mod, mod="Auto.allo", data=for_mod, group="Study2", xlab="Standardized mean difference", transfm="none")
jpeg("FigA3.jpeg", width = 5, height = 5, units = 'in', res = 300)
g + theme(axis.text.x=element_text(size=8), axis.text.y=element_text(angle=0, size=8),
          axis.title.x=element_text(size=10), axis.title.y=element_text(size=10))
dev.off()

#Figure A4: Funnel plot + trim-and-fill analysis (TAF requires simpler rma models)
simple_all <- rma(yi = yi, vi = vi, method = "REML", data = ES_all)
taf <- trimfill(simple_all)
jpeg("FigA4.jpeg", width = 6, height = 5, units = 'in', res = 300)
funnel(taf)
dev.off()
summary(taf)
simple_in_search <- rma(yi = yi, vi = vi, method = "REML", data = ES_in_search)
taf2 <- trimfill(simple_in_search)
funnel(taf2)
summary(taf2)
simple_not_in_search <- rma(yi = yi, vi = vi, method = "REML", data = ES_not_in_search)
taf3 <- trimfill(simple_not_in_search)
funnel(taf3)
summary(taf3)

#Figure A5: Orchard plot of effect sizes in vs. not in search results
for_mod <- ES_all
for_mod$Search_group[which(for_mod$Search_group == 0)] <- "In"
for_mod$Search_group[which(for_mod$Search_group == 1)] <- "Not"
in_mod <- rma.mv(yi = yi, V = vi, W = 1/vi, mods = ~ factor(Diff_dips1) + Search_group,
                 random = list(~1 | Study2, ~1 | ES2), method = "REML", data = for_mod)
summary(in_mod)
i2_ml(in_mod)
g <- orchard_plot(in_mod, mod="Search_group", data=for_mod, group="Study2", xlab="Standardized mean difference", transfm="none")
jpeg("FigA5.jpeg", width = 5, height = 5, units = 'in', res = 300)
g + theme(axis.text.x=element_text(size=8), axis.text.y=element_text(angle=0, size=8),
          axis.title.x=element_text(size=10), axis.title.y=element_text(size=10))
dev.off()


