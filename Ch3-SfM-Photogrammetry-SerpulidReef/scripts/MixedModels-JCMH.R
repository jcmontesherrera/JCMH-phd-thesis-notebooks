# R code used in:
# Towards reproducible analysis of benthos structural complexity: 
# A case study on Antarctic polychaete reefs using action cameras and remotely operated vehicles

# Authors: J.C. Montes-Herrera, G. Johnstone, J. Stark, N. Hill, V. Cummings, V. Lucieer

# Contact: juancarlos.montesherrera@utas.edu.au

###### Linear mixed models ######

library(ggplot2)
library(lme4)
library(rsq)
library(lmerTest)

path<- "C:\\Users\\jcmontes\\Documents\\GitHub\\ant_biogenic_structures\\data\\"

# Load file:
df<- read.csv(paste0(path, "for_LMM_merged_metrics.csv"))


# ----

# Fractal dimension (D32) linear models:
d32_pc_lm<- lm(fd32 ~ Polychaete.Colonies, data=df)
d32_mud_lm<- lm(fd32 ~ Mud, data=df)
d32_bpt_lm<- lm(fd32 ~ Broken.Polychaete.Tubes, data=df)

# Profile Curvature linear models:
profc_pc_lm<- lm(mean_profile_curvature ~ Polychaete.Colonies, data=df)
profc_mud_lm<- lm(mean_profile_curvature ~ Mud, data=df)
profc_bpt_lm<- lm(mean_profile_curvature ~ Broken.Polychaete.Tubes, data=df)

# Richness (rs)
rs_pc_lm <- lm(richness ~ Polychaete.Colonies, data=df)
rs_mud_lm<- lm(richness ~ Mud, data=df)
rs_bpt_lm<- lm(richness ~ Broken.Polychaete.Tubes, data=df)

# Shannon Index
sw_pc_lm <- lm(shannon.idx ~ Polychaete.Colonies, data=df)
sw_mud_lm<- lm(shannon.idx ~ Mud, data=df)
sw_bpt_lm<- lm(shannon.idx ~ Broken.Polychaete.Tubes, data=df)

# ----

# D32 linear mixed model with Transect as a random effect:
d32_pc_lme <- lmer(fd32 ~ Polychaete.Colonies + (1 | Transect),
                   data = df)
d32_mud_lme <- lmer(fd32 ~ Mud + (1 | Transect),
                    data = df)
d32_bpt_lme <- lmer(fd32 ~ Broken.Polychaete.Tubes + (1 | Transect),
                    data = df)


# Profile curvature linear mixed model with Transect as a random effect:
profc_pc_lme <- lmer(mean_profile_curvature ~ Polychaete.Colonies + (1 | Transect),
                     data = df)
profc_mud_lme <- lmer(mean_profile_curvature ~ Mud + (1 | Transect),
                      data = df)
profc_bpt_lme <- lmer(mean_profile_curvature ~ Broken.Polychaete.Tubes + (1 | Transect),
                      data = df)
profc_rs_lme <- lmer(mean_profile_curvature ~ richness + (1 | Transect),
                      data = df)

# Richness linear mixed model with Transect as a random effect:
rs_pc_lme <- lmer(richness ~ Polychaete.Colonies + (1 | Transect),
                     data = df)
rs_mud_lme <- lmer(richness ~ Mud + (1 | Transect),
                      data = df)
rs_bpt_lme <- lmer(richness ~ Broken.Polychaete.Tubes + (1 | Transect),
                      data = df)

# Shannon-Weiner Index linear mixed model with Transect as a random effect:
sw_pc_lme <- lmer(shannon.idx ~ Polychaete.Colonies + (1 | Transect),
                  data = df)
sw_mud_lme <- lmer(shannon.idx ~ Mud + (1 | Transect),
                   data = df)
sw_bpt_lme <- lmer(shannon.idx ~ Broken.Polychaete.Tubes + (1 | Transect),
                   data = df)


## ---- D32 linear mixed model with Random intercept and slope

d32_pc_rslme <- lmer(fd32 ~ Polychaete.Colonies + (1 + Polychaete.Colonies|Transect), 
                      data = df)

d32_mud_rslme <- lmer(fd32 ~ Mud + (1 + Mud|Transect), 
                      data = df)

d32_bpt_rslme <- lmer(fd32 ~ Broken.Polychaete.Tubes + (1 + Broken.Polychaete.Tubes|Transect), 
                       data = df)


## ---- Profile Curvature linear mixed model with Random intercept and slope

profc_pc_rslme <- lmer(mean_profile_curvature ~ Polychaete.Colonies + 
                        (1 + Polychaete.Colonies |Transect), 
                      data = df)

profc_mud_rslme <- lmer(mean_profile_curvature ~ Mud + 
                        (1 + Mud |Transect), 
                      data = df)

profc_bpt_rslme <- lmer(mean_profile_curvature ~ Broken.Polychaete.Tubes + 
                        (1 + Broken.Polychaete.Tubes|Transect), 
                      data = df)

## ---- Richness linear mixed model with Random intercept and slope

rs_pc_rslme <- lmer(richness ~ Polychaete.Colonies + 
                         (1 + Polychaete.Colonies |Transect), 
                       data = df)

rs_mud_rslme <- lmer(richness ~ Mud + 
                          (1 + Mud |Transect), 
                        data = df)

rs_bpt_rslme <- lmer(richness ~ Broken.Polychaete.Tubes + 
                          (1 + Broken.Polychaete.Tubes|Transect), 
                        data = df)

## ---- Shannon Index linear mixed model with Random intercept and slope

sw_pc_rslme <- lmer(shannon.idx ~ Polychaete.Colonies + 
                      (1 + Polychaete.Colonies |Transect), 
                    data = df)

sw_mud_rslme <- lmer(shannon.idx ~ Mud + 
                       (1 + Mud |Transect), 
                     data = df)

sw_bpt_rslme <- lmer(shannon.idx ~ Broken.Polychaete.Tubes + 
                       (1 + Broken.Polychaete.Tubes|Transect), 
                     data = df)



# ----
# Model decision based on AIC

anova(sw_pc_rslme, sw_pc_lme, sw_pc_lm)

## Performance and goodness-of-fit
performance::r2(sw_pc_lme)

# Plots
#Polychaete.Colonies
#mean_profile_curvature
#Broken.Polychaete.Tubes
#shannon.idx

library(ggeffects)
pred_mm <- ggpredict(sw_pc_lme, terms = c("Polychaete.Colonies")) # extract prediction dataframe

p <- (ggplot(pred_mm) + 
    geom_line(aes(x = x, y = predicted), color="black", size=1) +          # slope
    geom_ribbon(aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error), 
                fill = "darkgrey", alpha = 0.5) +  # error band
    geom_point(data = df,                      # adding the raw data (scaled values)
               aes(x = Polychaete.Colonies, y = shannon.idx, color = Transect), size=2) + 
    labs(x = "Polychaete.Colonies", y = "ShannonIDX", 
         title = "Polychaete.Colonies - Model 2") + 
    theme(legend.box.background = element_blank()) +
    scale_x_continuous(limits = c(0, 100), breaks = seq(0, 90, by = 30))
)

p

#Save figure
tiff(file="C:\\Users\\jcmontes\\Desktop\\sw-pc_M2.png", 
     units="in", width=5, height=4, res=300)
p
dev.off()


