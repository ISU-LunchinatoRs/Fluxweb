rm(list=ls())
## Fluxweb trials ##
install.packages('fluxweb')
library(fluxweb)

# Set Working Directory 
setwd("C:/Users/tjbutts/Desktop/R/Clear Lake Stability/Fluxweb material")

# Load data lists 
load("species.level.RData")
load("groups.level.RData")

# attach a list 
attach(species.level)

#Parameters 
mat 
biomasses 
bodymasses
efficiencies
names

detach(species.level)

#Can also do groups 
attach(groups.level)
mat 
biomasses 
bodymasses 
efficiencies 
species.tgs

detach(groups.level)

#### Fluxing ####
attach(species.level)
# Define metabolic losses
# General allometric equation (Brown et al. 2004_"Toward a Metabolic Theory of Ecology")
met.rates = 0.71*bodymasses^-0.25

# Efficiencies are already defined in the efficiencies vector 
# Then the netowrk of fluxes is obtained with fluxes function 
mat.fluxes = fluxing(mat, biomasses, met.rates, efficiencies)

# Can now estimate ecosystem functions such as herbivory, 
## detritivory, or carnivory 
# First, need to define basal species (without prey)

# basal species are species without prey 
basals = colSums(mat.fluxes) == 0
names[basals]

# plants are basal species that are not organic matter or exudates 
plants = basals 
plants[which(names =='dead organic matter' | names == 'root exudates')] = FALSE 

# Herbivory is defined as the sum of fluxes outgoing from plants
## consumers 
herbivory = sum(rowSums(mat.fluxes[plants, ]))
# Carnivory is defined as the sum of fluxes outgoing from animals 
carnivory = sum(rowSums(mat.fluxes[!basals, ]))
# Detritivory is defined as the sum of fluxes outgoing from
## detritus consumers 
detritivory = sum(mat.fluxes[names == 'dead organic matter', ])
# total fluxes 
total = sum(mat.fluxes) 

# plot 
fluxes <- c(herbivory, carnivory, detritivory, total)
xlab <- c('herbivory', 'carnivory', 'detritivory', 'total')

#Create a window and margins for the plot
windows(height=6, width=7)
par(mai=c(1,1.25,0.25,0.25)) 
barplot(fluxes, ylab = 'Fluxes (J/yr)', names.arg = xlab)

detach(species.level)

#### Stability ####
species <- species.level

# Matrix of fluxes is computed through the call to the fluxing 
## function as before 
# But we can modify the function for various purposes 
species.fluxes <- fluxing(species$mat, species$biomasses, met.rates, species$efficiencies, bioms.prefs = TRUE, bioms.losses = TRUE, ef.level = "prey")

# bioms.prefs specifies that species preferences depend on prey
## abundances

# bioms.losses is set to TRUE to compute metabolic losses for
## species populations (as they're provided per unit biomass) 

#ef.level is set to prey as efficiencies provided in these datasets depends on prey identitites 

species.fluxes

#Stability of the food web of fluxes is returned by the stability.value function 

##requires growth.rates vector
growth.rates = rep(NA, dim(species$mat)[1])
growth.rates[colSums(species$mat) == 0] = 0.5
growth.rates

#Compute species per unit biomass metabolic rates using the
## metabolic theory general allometric equation 
losses = 0.71 * species$bodymasses^(-0.25)

#Stability of food web 
#More Negative = More Stable

## Maximum eigenvalue of the Jacobian matrix of a Lotka voltera like system of equations 
## Can get all eigenvalues and eigenvectors if you want 
##The more negative, the more stable
stability.value(val.mat = species.fluxes, species$biomasses, losses, species$efficiencies, growth.rates, bioms.prefs = TRUE, bioms.losses = TRUE, ef.level = "prey")


#Use a scalar value multiplied against the losses term to see ## which value gives a stability of 0 (move from stable to unstable)
make.stability(val.mat = species.fluxes, species$biomasses, losses, species$efficiencies, growth.rates, bioms.prefs = TRUE, bioms.losses = TRUE, ef.level = "prey")

#function doesn't apply in this case 

#### Sensitivity Analysis ####
rm(list=ls())
attach(species.level)
load("species.level.RData")
set.seed(12)
losses = 0.71*bodymasses^-0.25

#create vectors to store the standard deviation of c.v. for each 
## uncertainty level 
sd.cvs.eff = c()
sd.cvs.los = c() 
sd.cvs.mat = c() 

for (var in seq(0, 0.12, 0.01)){ 
  cat( 'var: ', var, '\n')
  # for efficienceis 
  res = sensitivity(fluxing, "efficiencies", var, 50, 
                    mat = mat, 
                    biomasses = biomasses, 
                    losses = losses, 
                    efficiencies = efficiencies)
  sd.cvs.eff = c(sd.cvs.eff, mean(res[[2]], na.rm = T))
  
  # for losses
  res = sensitivity(fluxing, "losses", var, 50, 
                    mat = mat, 
                    biomasses = biomasses, 
                    losses = losses, 
                    efficiencies = efficiencies)
  sd.cvs.los = c(sd.cvs.los, mean(res[[2]], na.rm = T))
  
  # for prefrences 
  res = sensitivity(fluxing, "mat", var, 50, 
                    mat = mat, 
                    biomasses = biomasses, 
                    losses = losses, 
                    efficiencies = efficiencies)
  sd.cvs.mat = c(sd.cvs.mat, mean(res[[2]], na.rm = T))
}

# Plot 
windows(height=6, width=7)
par(mai=c(1,1.25,0.25,0.25))

seq1 <- seq(from=0, to=0.10, by=0.01)
plot(sd.cvs.eff ~ seq1, xlim = c(0, 0.12), 
     xlab = 'variation in parameters', 
     ylab = 'observed departure to original results', 
     pch = 16)
points(sd.cvs.los ~ seq1, col = 'red', pch = 16) 
points(sd.cvs.mat ~ seq1, col = 'green', pch = 16)
abline(a = 0, b = 1, lty = 2) 

legend('topleft', 
       legend = c('efficiency', 'metabolic losses', 
                  'species preferences'), 
       col = c('black', 'red', 'green'), 
       pt.cex=1.5, bty='n', 
       pt.bg = c('black', 'red', 'green'), pch = 21 )
