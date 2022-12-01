#######################################################################################
###################### Code for Hagen & Mason 2023 meta-analysis ######################
#######################################################################################

#Load in necessary packages
library(metafor)
library(metagear)
library(dplyr)
library(ape)
library(ggplot2)
library(emmeans)
library(xlsx)
library(orchaRd)
library(rlang)

#Load in effect sizes & phylogeny
ES_all <- read.xlsx("/Path/to/data/file", 1)
tree <- read.tree("/Path/to/tree/file")

#Change columns to numeric
ES_all[,c(21:26)] <- sapply(ES_all[,c(21:26)],as.numeric)

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
  dat$Pathogen_type[which(is.na(dat$Pathogen_type))] <- "Unknown"
  dat$Pathogen_lifestyle[which(is.na(dat$Pathogen_lifestyle))] <- "Unknown"
  return(dat)
}
ES_all <- convert_na(ES_all)
ES_in_search <- convert_na(ES_in_search)
ES_not_in_search <- convert_na(ES_not_in_search)

#Accounting for phylogeny in a multi-level meta-analytic model:
get_phylo_model <- function(dat){
  specs <- dat$Family ; specs <- unique(specs)
  phylo <- drop.tip(tree, tree$tip.label[which(!tree$tip.label %in% specs)])
  specs <- which(specs %in% phylo$tip.label)
  phy <- compute.brlen(phylo)
  cor <- vcv(phy, cor = T)
  comput_dat <- dat[which(dat$Family %in% phylo$tip.label),]
  names(comput_dat)[which(names(comput_dat) == "Family")] <- "phylo"
  phylo_m <- rma.mv(yi = yi, V = vi, W=1/vi, mods = ~ Polyploid_type + Cultivation_status + 
             Pathogen_type + Pathogen_lifestyle + Auto.allo, random = list(~1 | phylo, 
             ~1 | Study2, ~1 | ES2, ~1 | Diff_dips2), R = list(phylo = cor), method = "REML", 
             data = comput_dat)
  return(phylo_m)
}
ES_all_phymod <- get_phylo_model(ES_all)
summary(ES_all_phymod)
i2_ml(ES_all_phymod)
#No effect of phylogeny

#Robust multilevel model
robust_m <- rma.mv(yi = yi, V = vi, W = 1/vi, mods = ~ Polyploid_type + Cultivation_status 
                   + Pathogen_type + Pathogen_lifestyle + Auto.allo, random = list(~1 | Study2, 
                   ~1 | ES2, ~1 | Diff_dips2), method = "REML", data = ES_all)
summary(robust_m)
i2_ml(robust_m)

#Figure 1: PRISMA Diagram
phases <- c("START_PHASE: # of studies identified through database searching\n\n1602",
            "START_PHASE: # of additional studies identified through other sources\n\n55",
            "# of studies after duplicates removed\n\n1545",
            "# of studies with title and abstract screened\n\n1545",
            "EXCLUDE_PHASE: # of studies excluded\n\n1390",
            "# of full-text articles assessed for eligibility\n\n155",
            "EXCLUDE_PHASE: # of full-text excluded, not fitting eligibility criteria\n\n27",
            "final # of studies included in quantitative synthesis (meta-analysis)\n\n128")
plot_PRISMA(phases, design="classic")

#Figure 2: Caterpillars plot
make_caterpillar <- function(model, modz=NULL, dat){
  I2 <- i2_ml(model)
  disp <- round(I2[1],2)
  if(is.null(modz)){
    model_results <- mod_results(model, mod="1", group="ES1", data=dat) 
    #1 denotes intercept (the default)
  }else{
    model_results <- mod_results(model, mod = modz, group="ES1", data=dat)
  }
  p <- caterpillars(model_results, mod = "Int", xlab = "Standardized mean difference") 
  q <- p + annotate(geom = "text", x = -3.5, y = 20, label = paste0("italic(I)^{2} == ", 
                          disp), color = "black", parse = TRUE, size = 8)
  q + annotate(geom = "text", x = -2.45, y = 18, label = "%", color = "black", parse = FALSE, 
               size = 8)
}
make_caterpillar(robust_m, dat=ES_all)

#Figure 3: Orchard plot of families
for_mod <- ES_all[-which(ES_all$Family == "Actinidiaceae"),] #Plot won't work with Actinidiaceae
pf_mod <- rma.mv(yi = yi, V = vi, W = 1/vi, mods = ~ factor(Diff_dips1) + Family,
                 random = list(~1 | Study2, ~1 | ES2), method = "REML", data = for_mod)
summary(pf_mod)
i2_ml(pf_mod)
orchard_plot(pf_mod, mod="Family", data=for_mod, group="Study2", xlab="Standardized mean difference", transfm="none")

#Removing largest families
pruned_down <- for_mod[-which(for_mod$Family == "Musaceae"),]
pruned_down <- pruned_down[-which(for_mod$Family == "Poaceae"),]
pruned_mod <- rma.mv(yi = yi, V = vi, W = 1/vi, mods = ~ factor(Diff_dips1) + Family,
              random = list(~1 | Study2, ~1 | ES2), method = "REML", data = pruned_down)
summary(pruned_mod)

#Figure 4: Orchard plot of auto vs. allo
for_mod <- ES_all[-which(ES_all$Auto.allo == "Both"),]
for_mod <- for_mod[-which(for_mod$Auto.allo == "NA"),]
for_mod <- for_mod[,c(2,4,5,14,29,30)]
aa_mod <- rma.mv(yi = yi, V = vi, W = 1/vi, mods = ~ factor(Diff_dips1) + Auto.allo,
               random = list(~1 | Study2, ~1 | ES2), method = "REML", data = for_mod)
summary(aa_mod)
i2_ml(aa_mod)
orchard_plot(aa_mod, mod="Auto.allo", data=for_mod, group="Study2", xlab="Standardized mean difference", transfm="none")
#Simple rma model; removing intercept so all moderators appear
aa_simple <- rma(yi = yi, vi = vi, mods = ~ Auto.allo - 1, method = "REML", data = for_mod)
summary(aa_simple)

#Figure 5: Combined orchard plot of auto/allo and pathogen lifestyle
newcol <- paste0(ES_all$Auto.allo, "_", ES_all$Pathogen_lifestyle)
w <- 1/ES_all$vi
orch_dat <- as.data.frame(cbind(ES_all$ES2, ES_all$Study2, ES_all$Diff_dips2, newcol, ES_all$yi, ES_all$vi, w))
colnames(orch_dat) <- c("ES2", "Study2", "Diff_dips2", "Combo", "yi", "vi", "w")
greenlit <- c("Auto_Hemibiotrophic", "Allo_Necrotrophic", "Auto_Biotrophic", "Allo_Biotrophic", "Allo_Hemibiotrophic", "Auto_Necrotrophic")
orch_dat <- orch_dat[which(orch_dat$Combo %in% greenlit),]
combo_modfit <- rma.mv(yi = as.numeric(yi), V = as.numeric(vi), W = as.numeric(w), 
                  mods = ~ Combo, random = list(~1 | as.numeric(Study2), ~1 | as.numeric(ES2), 
                  ~1 | as.numeric(Diff_dips2)), method = "REML", data = orch_dat)
summary(combo_modfit)
i2_ml(combo_modfit)
orchard_plot(combo_modfit, mod="Combo", data=orch_dat, group="Study2", xlab="Standardized mean difference", transfm="none")

#Modeling just pathogen lifestyle
pl_mod <- rma.mv(yi = yi, V = vi, W = 1/vi, mods = ~ factor(Diff_dips1) + Pathogen_lifestyle,
                     random = list(~1 | Study2, ~1 | ES2), method = "REML", data = ES_all)
summary(pl_mod)

#Figure 6: Orchard plot of pathogen type
pt_mod <- rma.mv(yi = yi, V = vi, W = 1/vi, mods = ~ factor(Diff_dips1) + Pathogen_type,
                 random = list(~1 | Study2, ~1 | ES2), method = "REML", data = ES_all)
summary(pt_mod)
i2_ml(pt_mod)
orchard_plot(pt_mod, mod="Pathogen_type", data=ES_all, group="Study2", xlab="Standardized mean difference", transfm="none")

#Figure 7: Funnel plot + trim-and-fill analysis (TAF requires simpler rma models)
simple_all <- rma(yi = yi, vi = vi, method = "REML", data = ES_all)
taf <- trimfill(simple_all)
funnel(taf)
summary(taf)
simple_in_search <- rma(yi = yi, vi = vi, method = "REML", data = ES_in_search)
taf2 <- trimfill(simple_in_search)
funnel(taf2)
summary(taf2)
simple_not_in_search <- rma(yi = yi, vi = vi, method = "REML", data = ES_not_in_search)
taf3 <- trimfill(simple_not_in_search)
funnel(taf3)
summary(taf3)

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

#Figure 8: Orchard plot of effect sizes in vs. not in search results
for_mod <- ES_all
for_mod$Search_group[which(for_mod$Search_group == 0)] <- "In"
for_mod$Search_group[which(for_mod$Search_group == 1)] <- "Not"
in_mod <- rma.mv(yi = yi, V = vi, W = 1/vi, mods = ~ factor(Diff_dips1) + Search_group,
                 random = list(~1 | Study2, ~1 | ES2), method = "REML", data = for_mod)
summary(in_mod)
i2_ml(in_mod)
orchard_plot(in_mod, mod="Search_group", data=for_mod, group="Study2", xlab="Standardized mean difference", transfm="none")

