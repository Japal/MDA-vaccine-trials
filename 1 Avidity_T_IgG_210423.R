# Total IgG Avidity index data analysis

dat <- read.csv("Avidity_T_IgG.csv", as.is = FALSE)

names(dat) <- c("Trial","ID","Group","Age","Sex","Wgain",
                "cFEC","cFEC.red","WB","WB.red",
                "APY","MEP","CF","MIF","SAA","TGH",
                "ASP","ES20")

avi <- dat[,11:18]

# Impute a few NAs
library("mice")
avi <- complete(mice(avi, m=1, maxit = 50, method = 'pmm', seed = 500))
log.avi <- log(avi)

## Multivariate linear model testing for differences in (logged) avidity measures by sex and weight gain controlling for trial effect

# (Age is the same for all animal in a trial, so confounded with trial)

library(car)
avi.mlm <- lm(cbind(APY,MEP,CF,MIF,SAA,TGH,ASP,ES20) ~
                dat$Sex + dat$Wgain + dat$Trial, data=log.avi)
Anova(avi.mlm)
# Residuals
library("heplots")
cqplot(avi.mlm)

## PCA of avidity

library("FactoMineR")
library("factoextra")

# Avidities in log scale (large outlier in trial 4b, row 70, removed)
dat.supp <- data.frame(log.avi[-70,],dat[-70,c(8,10)]) 
dat.out <- dat[-70,]

res.pca <- PCA(dat.supp, quanti.sup = 9:10, graph = FALSE)

# Tailored biplot with different symbols by trial
fviz_pca_biplot(res.pca, col.quanti.sup="darkgreen",
                col.var = "blue3",geom.var=c("arrow","text"),geom.sup="point",
                col.ind = "black",
                geom.ind="",
                repel=FALSE,
                title="") +
  scale_shape_manual(values=c(16,15,17,3,4,5,7)) +
  geom_point(aes(shape=dat.out$Trial),color="black",size=2) +
  theme_bw() +
  # labs(x="PC1 (41.1%)",y="PC2 (16.8%)") +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=14),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12)) +
  theme(legend.position="bottom", legend.direction="horizontal") +
  guides(shape=guide_legend("Trial",ncol=7))

# Overall association of reduction in parasitology variables with avidity measures
avi.mlm <- lm(cbind(APY,MEP,CF,MIF,SAA,TGH,ASP,ES20) ~
                cFEC.red + WB.red, data=dat.supp)
Anova(avi.mlm)

# Association of reduction in parasitology variables with biplot axes by regression
log.avi.pca <- prcomp(dat.supp[,1:8])
avi.red.lm <- lm(dat.supp[,9] ~ log.avi.pca$x[,1:2])
summary(avi.red.lm) # R2 = 0.01817, 1.81% variability in cFEC.red explained by PC1 and PC2
avi.red.lm <- lm(dat.supp[,10] ~ log.avi.pca$x[,1:2])
summary(avi.red.lm) # R2 = 0.01573, 1.57% variability in WB.red explained by PC1 and PC2

cor(dat.supp) # Correlations between any reduction variable and indexes are < 0.2
