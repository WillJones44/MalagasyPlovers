# http://www.mpcm-evolution.com/practice/online-practical-material-chapter-11/chapter-11-1-simple-model-mcmcglmm

library(MCMCglmm)
library(ape)
library(readxl)
library(phytools)
library(tidyverse)
library(geiger)
library(Rmisc)
library(fGarch)
library(effects)
library(reshape2)

#### data import and wrangle ----
## get the tree
plover <- read_excel("data/Species.xlsx")
tree <- read.newick("data/phylogeny.nwk")
data <- column_to_rownames(plover, var = "Taxa")

# Sanity check to see if the tree looks ok
plot(tree)

# We need to prune the tree to remove all of the species that we do not have 
# survival estimates for
obj <- name.check(tree, data)
obj

shorttree <- drop.tip(tree, obj$tree_not_data)

# Sanity checks to see whether all the species in the dataset match those for 
# the tree file and that the tree looks ok
name.check(shorttree, data)
plot(shorttree)

# function to approximate the sd given the 95% CIs
approx_sd <- function(x1, x2){
  (x2-x1) / (qnorm(0.975) - qnorm(0.025))
}

# import data
mcmc_data <- 
  read_excel("data/mcmcglmm.xlsx") %>% 
  column_to_rownames(var = "Pop") %>%
  mutate(Survival = as.numeric(Survival),
         LowerCI = as.numeric(LowerCI),
         UpperCI = as.numeric(UpperCI),
         Population = row.names(.)) %>% 
  mutate(sd_est = ifelse(!is.na(LowerCI), 
                         approx_sd(x1 = LowerCI, x2 = UpperCI),
                         LowerCI))

# prepare the phylogeny for the model
phylo_short <- shorttree
inv.phylo <- inverseA(phylo_short, nodes = "TIPS", scale = TRUE)

# prepare the priors for the random effects (assume a noninformative prior for
# the method)
phylo_method_prior <- 
  list(G = list(G1 = list(V = 1, nu = 0.02),
                G2 = list(V = 1, nu = 0)), 
       R = list(V = 1, nu = 0.02))

phylo_prior <- 
  list(G = list(G1 = list(V = 1, nu = 0.02)), 
       R = list(V = 1, nu = 0.02))

# visually inspect how the uncertainty varies by population and method
ggplot(data = mcmc_data) +
  geom_point(aes(x = Population, y = Survival)) +
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI, x = Population)) +
  facet_grid(Method ~ .)

#### Survival model ----
survival_model <- MCMCglmm(Survival ~ scale(Clutch) * scale(EquatorDist) + 
                             scale(LogMass), 
                           random = ~phylo + Method,
                           prior = phylo_method_prior, family = "gaussian", 
                           ginverse = list(phylo = inv.phylo$Ainv),
                           data = mcmc_data, 
                           nitt = 5000000, burnin = 10000, thin = 5000)

# forest plot of results
summary(survival_model)$solutions %>%
  as.data.frame(.) %>% 
  rownames_to_column(., var = "variable") %>% 
  dplyr::rename(lower95 = 'l-95% CI',
                upper95 = 'u-95% CI') %>% 
  ggplot(.) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  geom_errorbarh(aes(xmin = lower95,
                     xmax = upper95,
                     y = as.factor(variable)),
                 alpha = 1, color = "grey30", 
                 size = 0.5,
                 height = 0) +
  geom_point(aes(y = variable, x = post.mean),
             size = 3, shape = 21, 
             fill = "#ECEFF4", col = "grey30", 
             alpha = 1, stroke = 0.5) +
  theme_bw() +
  theme(
    text = element_text(family = "Franklin Gothic Book"),
    axis.title.x = element_text(size = 10),
    axis.text.x  = element_text(size = 10), 
    axis.title.y = element_text(size = 12),
    axis.text.y = element_text(size = 10),
    strip.text = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(size = 0.5, colour = "grey40"),
    axis.ticks.length = unit(0.2, "cm"),
    panel.border = element_rect(linetype = "solid", colour = "grey")
  ) +
  ylab("Predictors of clutch size") +
  xlab(expression(paste("Standardized effect size (", beta,")" %+-% "95% CI", sep = "")))

#### Clutch model ----
# model to check colinearity between clutch size and the other predictors
clutch_model <- MCMCglmm(Clutch ~ scale(EquatorDist) + scale(LogMass), 
                         random = ~phylo,
                         prior = phylo_prior, family = "gaussian", 
                         ginverse = list(phylo = inv.phylo$Ainv),
                         data = mcmc_data, 
                         nitt = 5000000, burnin = 10000, thin = 5000)

# forest plot of results
summary(clutch_model)$solutions %>%
  as.data.frame(.) %>% 
  rownames_to_column(., var = "variable") %>% 
  dplyr::rename(lower95 = 'l-95% CI',
         upper95 = 'u-95% CI') %>% 
  ggplot(.) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  geom_errorbarh(aes(xmin = lower95,
                     xmax = upper95,
                     y = as.factor(variable)),
                 alpha = 1, color = "grey30", 
                 size = 0.5,
                 height = 0) +
  geom_point(aes(y = variable, x = post.mean),
             size = 3, shape = 21, 
             fill = "#ECEFF4", col = "grey30", 
             alpha = 1, stroke = 0.5) +
  theme_bw() +
  theme(
    text = element_text(family = "Franklin Gothic Book"),
    axis.title.x = element_text(size = 10),
    axis.text.x  = element_text(size = 10), 
    axis.title.y = element_text(size = 12),
    axis.text.y = element_text(size = 10),
    strip.text = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(size = 0.5, colour = "grey40"),
    axis.ticks.length = unit(0.2, "cm"),
    panel.border = element_rect(linetype = "solid", colour = "grey")
  ) +
  ylab("Predictors of clutch size") +
  xlab(expression(paste("Standardized effect size (", beta,")" %+-% "95% CI", sep = "")))

#### Survival bootstrap function ----
# Function to randomly draw a survival estimate from the uncertainty distributions
# provided in the source and run the MCMCglmm model
plover_surv_boot <- 
  function(nreps = 1000, df,
           prior, inv.phylo, 
           nitt = 5000000, burnin = 10000, thin = 5000){
  
  coef_storage_matrix <- matrix(numeric(5 * nreps), nreps)
  lambda_storage_matrix <- matrix(numeric(4 * nreps), nreps)
  model_storage_list <- vector(mode = "list", length = nreps)
  
  for(i in 1:nreps){
    
    df$surv_boot <- 
      suppressWarnings(apply(df, 1, 
            function(x) rnorm(1, mean = as.numeric(x["Survival"]), 
                              sd = as.numeric(x["sd_est"]))))
    
    df$surv_boot <- 
      ifelse(!is.na(df$surv_boot), df$surv_boot, df$Survival)
    
    model_simple <- MCMCglmm(surv_boot ~ scale(Clutch) * scale(EquatorDist) + scale(LogMass), 
                             random = ~phylo + Method, prior = prior, 
                             family = "gaussian", 
                             ginverse = list(phylo = inv.phylo$Ainv), 
                             data = df, nitt = nitt, 
                             burnin = burnin, thin = thin)
  
    coef_storage_matrix[i, 1] <- summary(model_simple)$solutions[1, 1]
    coef_storage_matrix[i, 2] <- summary(model_simple)$solutions[2, 1]
    coef_storage_matrix[i, 3] <- summary(model_simple)$solutions[3, 1]
    coef_storage_matrix[i, 4] <- summary(model_simple)$solutions[4, 1]
    coef_storage_matrix[i, 5] <- summary(model_simple)$solutions[5, 1]
    
    lambda <- 
      model_simple$VCV[,'phylo'] /
      (model_simple$VCV[,'phylo'] + model_simple$VCV[,'units'])
    
    lambda_storage_matrix[i, 1] <- mean(lambda)
    lambda_storage_matrix[i, 2] <- HPDinterval(lambda)[1]
    lambda_storage_matrix[i, 3] <- HPDinterval(lambda)[2]
    lambda_storage_matrix[i, 4] <- posterior.mode(lambda)
  
    model_storage_list[[i]] <- model_simple
  }
  
  lambda_out <- 
    as.data.frame(lambda_storage_matrix) %>% 
    dplyr::rename(mean = V1,
                  lower95 = V2,
                  upper95 = V3, 
                  posterior_mode = V4) %>% 
    mutate(iteration = c(1:nreps))
  
  coefficients_out <- 
    as.data.frame(coef_storage_matrix) %>% 
    dplyr::rename(Intercept = V1,
                  Clutch = V2,
                  EquatorDist = V3, 
                  LogMass = V4,
                  ClutchxEquatorDist = V5) %>% 
    t(.) %>% 
    melt() %>% 
    dplyr::rename(variable = Var1,
                  iteration = Var2,
                  estimate = value)
  
  results_list <- list(lambda_out = lambda_out, 
                       coefficients_out = coefficients_out,
                       models = model_storage_list)
  
  return(results_list)
}

#### Run survival bootstrap ----
# NOTE: the defaults are "nitt = 5000000, burnin = 10000, thin = 5000"
# if you want to do a test run then, change these specs manually and the nreps
boot_out <- 
  plover_surv_boot(nreps = 1000, df = mcmc_data, 
                   prior = phylo_method_prior, 
                   inv.phylo = inv.phylo)

# plot of standardized effect size distribution of bootstrap results
boot_out$coefficients_out %>% 
  dplyr::filter(variable != "Intercept") %>% 
  ggplot(.) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  geom_violin(aes(factor(variable), estimate), draw_quantiles = c(0.025, 0.5, 0.975)) +
  theme_bw() +
  theme(
    text = element_text(family = "Franklin Gothic Book"),
    axis.title.x = element_text(size = 12),
    axis.text.x  = element_text(size = 10), 
    axis.title.y = element_text(size = 12),
    axis.text.y = element_text(size = 10),
    strip.text = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(size = 0.5, colour = "grey40"),
    axis.ticks.length = unit(0.2, "cm"),
    panel.border = element_rect(linetype = "solid", colour = "grey")
  ) +
  coord_flip() +
  ylab(expression(paste("Standardized effect size (", beta,") median" %+-% "95% CI", sep = ""))) +
  xlab("Predictors of apparent survival") +
  scale_x_discrete(labels = c("Clutch" = "clutch size",
                              "EquatorDist" = "distance from equator",
                              "LogMass" = "body mass (log)",
                              "ClutchxEquatorDist" = "clutch size : dist. from equator"))