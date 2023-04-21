## MFA of immunology and parasitology datasets 

library("FactoMineR")
library("factoextra")
library("ggthemes")
library("vegan")

# Original data already pre-processed from previous analyses
load("IgGSerumBiplotData.RData") 
load("IgASerumBiplotData.RData")
load("IgGMucusBiplotData.RData")
load("IgAMucusBiplotData.RData")

names(dat.supp.IgGSerum) <- c("APY","MEP","CF","MIF","SAA","TGH","ASP","ES20","cFEC.red","WB.red")
names(dat.supp.IgASerum) <- c("APY","MEP","CF","MIF","SAA","TGH","ASP","ES20","cFEC.red","WB.red")
names(dat.supp.IgGMucus) <- c("APY","MEP","CF","MIF","SAA","TGH","ASP","ES20","cFEC.red","WB.red")
names(dat.supp.IgAMucus) <- c("APY","MEP","CF","MIF","SAA","TGH","ASP","ES20","cFEC.red","WB.red")

trial <- read.csv("Levels_T_IgG_Serum.csv")[,1]
tmp <- which(trial%in%c("1","2a","2b","3","5"))
trial <- as.factor(trial[tmp])
trial <- droplevels(trial)

group <- read.csv("Levels_T_IgG_Serum.csv")[,3]
group <- as.factor(group[tmp])

# Extract parasitology variables
para <- dat.supp.IgGSerum[,9:10]

dat1 <- dat.supp.IgGSerum[,-c(9,10)]
dat2 <- dat.supp.IgASerum[,-c(9,10)]
dat3 <- dat.supp.IgGMucus[,-c(9,10)]
dat4 <- dat.supp.IgAMucus[,-c(9,10)]

names(dat1) <- paste0("GS.",names(dat1))
names(dat2) <- paste0("AS.",names(dat2))
names(dat3) <- paste0("GM.",names(dat3))
names(dat4) <- paste0("AM.",names(dat4))

dat <- cbind(dat1,dat2,dat3,dat4,para)
grn <- c(ncol(dat1),ncol(dat2),ncol(dat3),ncol(dat4),ncol(para))

dat.mfa <- MFA(dat,group=grn,type=c("s","s","s","s","s"),ncp=2,
               name.group=c("IgG.Serum","IgA.Serum","IgG.Mucus","IgA.Mucus","Parasitology"),graph = FALSE)

fviz_mfa_var(dat.mfa,
             geom.var=c("arrow","text"),
             addEllipses = FALSE,
             repel=TRUE,
             title="") +
  #labs(x="MFA1 (48.8%)",y="MFA2 (17.4%)") +
  scale_colour_colorblind() +
  theme_bw() +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=14),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12)) +
  theme(legend.position="bottom", legend.direction="vertical",
        legend.box="vertical", legend.margin=margin()) +
  theme(plot.margin = unit(c(-8,1.5,0,0.1), "cm")) +
  guides(col=guide_legend("Dataset",nrow=3,byrow = TRUE))

fviz_mfa_ind(dat.mfa,
             habillage = group,
             geom="",
             addEllipses = FALSE,
             repel=TRUE,
             title="") +
  #labs(x="MFA1 (48.8%)",y="MFA2 (17.4%)") +
  scale_fill_colorblind() +
  scale_colour_colorblind() +
  scale_shape_manual(values=c(16,15,17,3,7)) +
  geom_point(aes(shape=trial,col=group),size=2) +
  theme_bw() +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=14),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12)) +
  theme(legend.position="bottom", legend.direction="vertical",
        legend.box="vertical", legend.margin=margin(),
        legend.box.just="left") +
  theme(plot.margin = unit(c(0.75,0.3,8.8,-1), "cm")) +
  guides(col=guide_legend("Treatment group",nrow=1),fill=guide_legend("Treatment group",nrow=1),
         shape=guide_legend("Trial",nrow=1,byrow = TRUE))


## RDA 1: Is total IgG response in serum influencing IgG response in abomasal mucus?

# Control for Group and Trial effects
dat1plus <- data.frame(dat1,Group=group,Trial=trial)

# Partial RDA controlling for trial and treatment group effects
rda.fit <- rda(dat3 ~ GS.APY + GS.MEP + GS.CF + GS.MIF + GS.SAA + GS.TGH + GS.ASP + GS.ES20 + Condition(Group + Trial),
               data=dat1plus, scale = TRUE)

rda.fit # Variance partition
# Global test of significance of the RDA result
set.seed(1234)
anova.cca(rda.fit)
# Test of all RDA  axes
set.seed(1234)
anova.cca(rda.fit, by = "axis")
# Goodness of fit
RsquareAdj(rda.fit)
rda.fit.sum <- summary(rda.fit)
# Variability explained by RDA axis
rda.fit.sum$concont$importance

## RDA 2: Is total IgA response in abomasal mocus influencing IgA response in serum?

# Control for Group and Trial effects
dat4plus <- data.frame(dat4,Group=group,Trial=trial)

# Partial RDA controlling for trial and treatment group effects
rda.fit <- rda(dat2 ~ AM.APY + AM.MEP + AM.CF + AM.MIF + AM.SAA + AM.TGH + AM.ASP + AM.ES20 + Condition(Group + Trial),
               data=dat4plus, scale = TRUE)

rda.fit # Variance partition
# Global test of significance of the RDA result
set.seed(1234)
anova.cca(rda.fit)
# Test of all RDA  axes
set.seed(1234)
anova.cca(rda.fit, by = "axis")
# Goodness of fit
RsquareAdj(rda.fit)
rda.fit.sum <- summary(rda.fit)
# Variability explained by RDA axis
rda.fit.sum$concont$importance
