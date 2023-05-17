# Pigment analysis
## J.C. Montes-Herrera

library(lme4)
library(rsq)
library(lmerTest)

#setwd("C:/Users/jcmontes/OneDrive - University of Tasmania/01_Projects_Drive/Imaging_spectroscopy/Phenotyping_macroalgae/results/Hyperspectral_Imaging/res_tables")

path <- "C:\\Users\\jcmontes\\OneDrive - University of Tasmania\\01_Projects_Drive\\Imaging_spectroscopy\\Phenotyping_macroalgae\\results\\Hyperspectral_Imaging\\res_tables\\"

df<- read.csv(paste0(path, "PE-pigments-indices.csv"))
#df<- read.csv(paste0(path, "CHL-pigments-indices.csv"))
#df<- read.csv(paste0(path, "PHAEO-pigments-indices.csv"))

# looping through data
# Get names of spectral indices (index)
index <- colnames(df)[2:length(colnames(df))]

# Result matrix
res <- matrix(NA, nrow = 4, ncol = length(index),
              dimnames = list("stats" = c("slope","intercept","R2","pvalue"), "index" = index))

for(iCol in colnames(df)[4:length(colnames(df))])
{
  
  my_lme<- lmer(PE_area ~ df[,iCol] + (1 | Site), data=df) # mixed model, with site as random effect
  #my_lm<- lm(PE_area ~ df[,iCol], data=df)
  #op4 <- lmer(PE_area ~ dd_PE569 + (1|Site), data=df, REML=T)
  
  res[1, iCol] <- summary(my_lme)$coefficients[2] # m value, slope
  res[2, iCol] <- summary(my_lme)$coefficients[1] # b value, intercept
  res[3, iCol] <- rsq(my_lme, adj=TRUE)[[1]] # R^2
  res[4, iCol] <- anova(my_lme)[[6]] #(Pr(>F))
  #res[4, iCol] <- summary(my_lme)$coefficients[9]
  #res[4,iCol] <- anova(my_lme)[["Pr(>F)"]][[1]] #p value
}

res

write.table(res,file=(paste0(path, "PE-LMM.csv")))
#write.matrix(res,file="CHL-LMM.csv")
