################################################################
################################################################
######                                                  ########
######       SCRIPT 3 -  POPULATION LEVEL MODELS        ########    
######                                                  ########
################################################################
################################################################

# First sets up model without fitting, then sets random 
# intercept variance to large fixed value (i.e. Muff method), 
# other variance components to 0, and fixed intercepts for each 
# step-id. Finally fits model using Poisson likelihood

########################
##### Load packages ----
libs <- c('lubridate', 'amt', 'data.table', 'tidyverse', 'glmmTMB')
lapply(libs, require, character.only=TRUE)

#################################
#### Load cleaned track data ----

# Select whether to use quadrat or spatial smoothing kernel density data
method <- 'kernel_d'
method <- 'quadrat'

# Select whether to use proportion or binary rasters
type <- 'proportion'
type <- 'binary'

dat <- readRDS(paste('output/tracks/stps_elk', paste(type, paste(method, 'cleaned.rds', sep='_'), sep='_'), sep='_'))

################################################################
#### Subset data around survey date, add col for HR density ----
dat_sub <- data.table()
for(i in unique(dat$id)){
  dat_id <- dat[id==i]
  # Check if > 2.5% steps in agriculture
  # if(nrow(dat_id[case_=='TRUE' & agricultural > 0])==0){
  #   next
  # }
  # Specify number of days around median survey date to buffer
  date_buffer <- 42
  # Check if fewer than 80% of expected days in date buffer period
  # if(nrow(dat_id[case_==TRUE])<((date_buffer*12))*.8){
  #   next
  # }
  # Subset rows
  dat_id <- dat_id[j_day >= dat_id[1]$survey_date-date_buffer & j_day <= dat_id[1]$survey_date+date_buffer]
  # Make a quick HR density value
  dat_id[, 'avg_dens' := mean(density_end, na.rm=T)]
  # Bind individuals back together
  dat_sub <- rbind(dat_sub, dat_id)
}

##########################################
####   1 - Fit model using Muff method ###
##########################################

# Habitat covariates: agriculture, conifer, mixedwood
# Interaction covariates: 
#     - Density at start of step and habitat at end, 
#     - Home range density and habitat at end,
#     - Population density and habitat at end
# Random covariates: 
#     - Intercept for each step cluster,
#     - Intercepts for each individual,
#     - Random slopes for all covariates

# Set max iterations to 10,000 to aid convergence
glmmTMBControl(optCtrl=list(iter.max=1e4,eval.max=1e4))

# Set up model but don't fit
pop_model <- glmmTMB(case_ ~ agricultural + mixedwood + coniferous + # Selection for habs at end of step
                       density_start:agricultural + density_start:mixedwood + density_start:coniferous + # Step-level effect of density on selection
                       # agricultural:avg_dens + mixedwood:avg_dens + coniferous:avg_dens + # Home range-level effect of density on selection
                       # agricultural:rel_pop_count + mixedwood:rel_pop_count + coniferous:rel_pop_count + # Population-level effect of density on selection
                       (1|step_id) + # Individual intercepts
                       (0 + agricultural|id) + (0 + mixedwood| id) + (0 + coniferous| id) + # Random effects on selection for habs at end of step
                       # (0 + agricultural:avg_dens|id) + (0 + mixedwood:avg_dens|id) + (0 + coniferous:avg_dens|id) +  # Random effects of home range-level density on selection
                       (0 + density_start:agricultural|id) + (0 + density_start:mixedwood| id) + (0 + density_start:coniferous| id), # Random effects of step-level density on selection
                     family=poisson(), 
                     data = dat_sub, doFit=FALSE)

# Set variance of random intercept to 10^6
pop_model$parameters$theta[1] <- log(1e6)
nvar_parm <- length(pop_model$parameters$theta)
pop_model$mapArg <- list(theta = factor(c(NA, 1:(nvar_parm - 1))))

# Fit model using large fixed variance
pop_model <- glmmTMB:::fitTMB(pop_model)

#######################################
####   2 - Tidy and save the output ###
#######################################

# Tidy output
tidy_pop <- as.data.table(broom.mixed::tidy(pop_model, effects='fixed'))

# Add confidence intervals
tidy_pop[, c('lower', 'upper') := .(
  estimate-(std.error*1.96),
  estimate+(std.error*1.96)
)]

# Extract random effects
pop_ranefs <- as.data.table(ranef(pop_model)$cond$id)
# Melt
pop_ranefs <- melt(pop_ranefs, variable.name='term', value.name='ind_estimate')
# Remove any individuals without a random effect
pop_ranefs <- pop_ranefs[!ind_estimate==0]

# Save 
saveRDS(tidy_pop, file = "output/models/population/muffmod_fixedeffs.rds")
saveRDS(pop_ranefs, file = "output/models/population/muffmod_ranefs.rds")

#########################################
####   3 - Factor and plot the output ###
#########################################

# Factor population model
tidy_pop$term <- factor(tidy_pop$term)
levels(tidy_pop$term) <- c('(Intercept)', 'Ag', 'Ag:HR_D', 'Ag:D_St', 'Conifer', 'Conifer:HR_D', 'Conifer:D_St',  'MW', 'MW:HR_D', 'MW:D_St')
levels(tidy_pop$term) <- c('(Intercept)', 'Ag', 'Ag:D_St', 'Conifer', 'Conifer:D_St',  'MW', 'MW:D_St')

# Factor random effects
pop_ranefs$term <- factor(pop_ranefs$term)
levels(pop_ranefs$term) <- c('Ag', 'MW', 'Conifer', 'Ag:HR_D', 'MW:HR_D', 'Conifer:HR_D', 'Ag:D_St', 'MW:D_St', 'Conifer:D_St')
levels(pop_ranefs$term) <- c('Ag', 'MW', 'Conifer', 'Ag:D_St', 'MW:D_St', 'Conifer:D_St')

# Add column for fixed + random effect
pop_ranefs <- merge(pop_ranefs, tidy_pop[,3:4], by='term', all=F)
pop_ranefs[, 'blup' := .(estimate+ind_estimate)]

# Plot
# tiff('figures/population_model.tiff', width = 12, height = 8, units = 'in', res = 300)
ggplot(tidy_pop[!term %in% '(Intercept)'], aes(x=term, y=estimate)) + geom_point(size=3) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2) +
  geom_point(data=pop_ranefs, aes(x=term, y=blup), pch=20, alpha=0.5, position='jitter') +
  geom_hline(yintercept = 0, alpha=0.5, linetype='dashed') +
  theme(panel.background = element_rect(fill='white', colour='black'), panel.grid = element_blank(),
        axis.text=element_text(size=15, colour='black'), axis.title.x=element_text(size=15, colour='black', vjust=-2),
        axis.title.y=element_text(size=15, colour='black', vjust=2), plot.margin=unit(c(1,0.5,1,0.5), 'cm')) +
  ylab('Selection coefficient') + xlab('')

