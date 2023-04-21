# Levels of total IgA in serum data analysis 

dat <- read.csv("Levels_T_IgA_Serum.csv", as.is = FALSE)

names(dat) <- c("Trial","ID","Group",
                "cFEC","cFEC.red","WB","WB.red",
                "APY","MEP","CF","MIF","SAA","TGH",
                "ASP","ES20")

# Remove trial 4 data
dat <- dat[!(dat$Trial%in%c("4a","4b")),]
dat$Trial <- droplevels(dat$Trial)

lev <- dat[,8:15]
log.lev <- log(lev) # Logged data

# Examine overall difference in covariance matrices between groups
library("heplots")
res <- boxM(cbind(APY,MEP,CF,MIF,SAA,TGH,ASP,ES20) ~ dat$Group, data = log.lev)
res
plot(res)

## Parametric MANOVA testing of differences in T IgG levels between treatment groups controlling for trial
library("car")
lev.aov <- lm(cbind(APY,MEP,CF,MIF,SAA,TGH,ASP,ES20) ~ dat$Group*dat$Trial, data = log.lev)
Anova(lev.aov)

# PERMANOVA
library("vegan")
adonis(cbind(APY,MEP,CF,MIF,SAA,TGH,ASP,ES20) ~ dat$Group*dat$Trial, method = "euclidean", data = log.lev)

# MANOVA for heterogeneous covariance matrices (although in our case, as balanced, PERMANOVA should work fine)

# James' test (function used only allows one factor at a time)
library("Compositional")
tmp <- data.frame(log.lev, Group = dat$Group, Trial = dat$Trial)
maovjames(as.matrix(tmp[,1:8]), tmp[,9]) # By group

# Modified ANOVA-type statistic (MATS)
library("MANOVA.RM")
mats <- MANOVA.wide(cbind(APY,MEP,CF,MIF,SAA,TGH,ASP,ES20) ~ Group*Trial, CPU = 1, seed = 234, data = tmp, iter=1000)
mats

## PCA of T IgA levels

library("FactoMineR")
library("factoextra")
library("ggthemes")

dat.supp <- data.frame(log.lev,dat[,c(5,7)])

res.pca <- PCA(dat.supp, quanti.sup = 9:10, scale.unit = TRUE, graph = FALSE)

fviz_pca_biplot(res.pca, col.quanti.sup="darkgreen",
                habillage = dat$Group,
                col.var = "blue3",geom.var=c("arrow","text"),
                geom.ind="",
                addEllipses = TRUE,
                repel=TRUE,
                title="")  +
  # labs(x="PC1 (86.5%)",y="PC2 (4.7%)") +
  scale_fill_colorblind() +
  scale_colour_colorblind() +
  scale_shape_manual(values=c(16,15,17,3,7)) +
  geom_point(aes(shape=dat$Trial,col=dat$Group),size=2) +
  theme_bw() +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=14),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12)) +
  theme(legend.position="bottom", legend.direction="horizontal",
        legend.box="vertical", legend.margin=margin()) +
  guides(col=guide_legend("Treatment group"),fill=guide_legend("Treatment group"),
         shape=guide_legend("Trial"))


# Overall association of reduction in parasitology variables with responses
avi.mlm <- lm(cbind(APY,MEP,CF,MIF,SAA,TGH,ASP,ES20) ~
                cFEC.red + WB.red, data=dat.supp)
Anova(avi.mlm)

# Association of reduction in parasitology variables with biplot axes by regression
log.avi.pca <- prcomp(dat.supp[,1:8])
avi.red.lm <- lm(dat.supp[,9] ~ log.avi.pca$x[,1:2])
summary(avi.red.lm) # R2 = 0.0009176, 0.09% variability in cFEC.red explained by PC1 and PC2
avi.red.lm <- lm(dat.supp[,10] ~ log.avi.pca$x[,1:2])
summary(avi.red.lm) # R2 = 0.006564, 0.66% variability in WB.red explained by PC1 and PC2

cor(dat.supp) # Correlations between any reduction variable and indexes are < 0.15

## Quadratic Discriminant analysis

library(caret)

set.seed(1244)
fitControl <- trainControl(## 10-fold CV
  method = "repeatedcv",
  number = 10,
  ## repeated ten times
  repeats = 5)

training <- data.frame(dat.supp[,1:8], Group = dat$Group)

lev.qda <- train(Group ~ ., data = training, method = "qda", trControl = fitControl,verbose = FALSE)
lev.qda

varImp(lev.qda)
