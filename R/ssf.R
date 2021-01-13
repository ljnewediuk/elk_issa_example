# -------------------------------------------#

# Example code for elk iSSA examining effect
# of elk density on habitat selection via
# interaction between density at the start of
# the step and habitat at the end of a step.
# Using method from Muff et al. 2019 
# 10.1111/1365-2656.13087.

# -------------------------------------------#

########################
##### Load packages ----
libs <- c('data.table', 'dplyr', 'ggplot2', 'glmmTMB')
lapply(libs, require, character.only=TRUE)

####################
#### Load data ----
# Cleaned track data
dat <- readRDS('input/ssf_dat.rds')
# Data specifying whether or not habitats are used
used_ppns <- readRDS('input/used_proportions.rds')

#############################################
####   1 - Fit model using poisson method ###
#############################################
# Habitat covariates: agriculture, conifer, mixedwood
# Fixed effects: 
#     - Habitats (proportion): agriculture, mixedwood, coniferous
#     - Density at start of step and habitat at end interaction
# Random effects: 
#     - Intercept for step cluster
#     - Intercepts for elkyear
#     - Conditional effects for habitats and interactions

# Set max iterations to 10,000 to aid convergence
glmmTMBControl(optCtrl=list(iter.max=1e4,eval.max=1e4))

# Set up model but don't fit (some Matrix warnings)
pop_model <- glmmTMB(case_ ~ agricultural + mixedwood + coniferous + # Selection for habs at end of step
                       density_start:agricultural + density_start:mixedwood + density_start:coniferous + # Step-level effect of density on selection
                       # agricultural:avg_dens + mixedwood:avg_dens + coniferous:avg_dens + # Home range-level effect of density on selection
                       # agricultural:rel_pop_count + mixedwood:rel_pop_count + coniferous:rel_pop_count + # Population-level effect of density on selection
                       (1|step_id) + # Individual intercepts
                       (0 + agricultural| elkyear) + (0 + mixedwood| elkyear) + (0 + coniferous| elkyear) + # Random effects on selection for habs at end of step
                       # (0 + agricultural:avg_dens|id) + (0 + mixedwood:avg_dens|id) + (0 + coniferous:avg_dens|id) +  # Random effects of home range-level density on selection
                       (0 + density_start:agricultural| elkyear) + (0 + density_start:mixedwood| elkyear) + (0 + density_start:coniferous| elkyear), # Random effects of step-level density on selection
                     family=poisson(), 
                     data = dat, doFit=FALSE)

# Set variance of random intercept to 10^6
pop_model$parameters$theta[1] <- log(1e6)
nvar_parm <- length(pop_model$parameters$theta)
pop_model$mapArg <- list(theta = factor(c(NA, 1:(nvar_parm - 1))))

# Fit model using large fixed variance
pop_model <- glmmTMB:::fitTMB(pop_model)

########################################################
####   2 - Tidy and organize the output for plotting ###
########################################################

tidy_pop <- as.data.table(broom.mixed::tidy(pop_model, effects='fixed'))

# Add confidence intervals
tidy_pop[, c('lower', 'upper') := .(
  estimate-(std.error*1.96),
  estimate+(std.error*1.96)
)]

# Extract random effects
pop_ranefs <- as.data.table(ranef(pop_model)$cond$elkyear, keep.rownames=T)
colnames(pop_ranefs)[1] <- 'elkyear'

# Melt so estimates are measured by model term (data.table warning)
pop_ranefs <- melt(pop_ranefs, variable.name='term', value.name='ind_estimate')

# Remove any individuals without an estimate
pop_ranefs <- pop_ranefs[!ind_estimate==0]

# Factor fixed effects
tidy_pop$term <- factor(tidy_pop$term)
levels(tidy_pop$term) <- c('(Intercept)', 'Ag', 'Ag:D_St', 'Conifer', 'Conifer:D_St',  'MW', 'MW:D_St')

# Factor random effects
pop_ranefs$term <- factor(pop_ranefs$term)
levels(pop_ranefs$term) <- c('Ag', 'MW', 'Conifer', 'Ag:D_St', 'MW:D_St', 'Conifer:D_St')

# Add column for fixed + random effect
plot_dat <- merge(pop_ranefs, tidy_pop[,3:4], by='term', all=F)
plot_dat[, 'blup' := .(estimate+ind_estimate)]

# Add column for used proportions for each term
plot_dat <- merge(plot_dat, used_ppns, by=c('term', 'elkyear'))

#############################################
####   3 - Plot fixed and random effects ###
############################################
# Habitat covariates:
#   - "Ag" = agriculture; "Conifer" = coniferous; "MW" = mixedwood
#   - "Habitat:D_St" = habitat at end of step interacted with density at start of step
# Random effects (jittered points) are coloured by "ppn" (proportion). If the individual 
#   has > 0 used steps in a habitat, its ppn value is "used". If none, its ppn value is "not_used". 
#   If one of either the habitat or density layer is used but not both (applies to Habitat:D_St 
#   covariates), the pn value is "one used".
# Black points are fixed effects, bars are 95% confidence intervals

# Plot
# tiff('figures/population_model.tiff', width = 12, height = 8, units = 'in', res = 300)
ggplot(tidy_pop[!term %in% '(Intercept)'], aes(x=term, y=estimate)) + geom_point(size=3) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2) +
  geom_point(data=plot_dat, aes(x=term, y=blup, colour=ppn), pch=20, alpha=0.5, position='jitter') +
  geom_hline(yintercept = 0, alpha=0.5, linetype='dashed') +
  theme(panel.background = element_rect(fill='white', colour='black'), panel.grid = element_blank(),
        axis.text=element_text(size=15, colour='black'), axis.title.x=element_text(size=15, colour='black', vjust=-2),
        axis.title.y=element_text(size=15, colour='black', vjust=2), plot.margin=unit(c(1,0.5,1,0.5), 'cm')) +
  ylab('Selection coefficient') + xlab('')

