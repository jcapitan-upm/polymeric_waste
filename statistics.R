rm(list = ls())

####
# Shore D surface hardness analysis
####

data <- read.csv("Dureza.csv", sep = ";")
data$MORTERO <- as.factor(data$MORTERO)
data$EDAD <- as.factor(data$EDAD)

# Plot figure with dosage and age as factors

library("ggpubr")
ggline(data, x = "MORTERO", y = "DUREZA", color = "EDAD", add=c("mean_se","jitter")) +
  labs(y = "Shore D surface hardness") + labs(x = "Composite material") + 
  scale_color_discrete(name = "Age", labels = c("28 days", "90 days", "180 days")) + 
  scale_x_discrete(labels=c('REF', 'CW-25%', 'CW-50%', 'CW-75%', 'CW-100%'))

# Two-way ANOVA with interaction: effect of dosage

res.aov <- aov(DUREZA ~ MORTERO*EDAD, data = data)
summary(res.aov)

# Test of the normality of the residuals

aov_residuals <- residuals(object = res.aov)
ks.test(unique(aov_residuals), pnorm, 0, sd(aov_residuals))

# Test for homoscedasticity (with a logarithmic transformation)

library("car")
leveneTest(DUREZA ~ MORTERO*EDAD, data = data)
leveneTest(log(DUREZA) ~ MORTERO*EDAD, data = data)

# Bootstrap to account for heteroscedasticity

library(twowaytests)
gpTwoWay(DUREZA ~ MORTERO*EDAD, data = data, method = "gPB")

# Pairwise comparisons: for fixed ages

pairwise.t.test(data$DUREZA[which(data$EDAD=="028D")], 
                data$MORTERO[which(data$EDAD=="028D")], p.adjust.method = "bonferroni")
pairwise.t.test(data$DUREZA[which(data$EDAD=="090D")], 
                data$MORTERO[which(data$EDAD=="090D")], p.adjust.method = "bonferroni")
pairwise.t.test(data$DUREZA[which(data$EDAD=="180D")], 
                data$MORTERO[which(data$EDAD=="180D")], p.adjust.method = "bonferroni")

# Pairwise comparisons: for fixed dosages

pairwise.t.test(data$DUREZA[which(data$MORTERO=="000CW")], 
                data$EDAD[which(data$MORTERO=="000CW")], p.adjust.method = "bonferroni")
pairwise.t.test(data$DUREZA[which(data$MORTERO=="025CW")], 
                data$EDAD[which(data$MORTERO=="025CW")], p.adjust.method = "bonferroni")
pairwise.t.test(data$DUREZA[which(data$MORTERO=="050CW")], 
                data$EDAD[which(data$MORTERO=="050CW")], p.adjust.method = "bonferroni")
pairwise.t.test(data$DUREZA[which(data$MORTERO=="075CW")], 
                data$EDAD[which(data$MORTERO=="075CW")], p.adjust.method = "bonferroni")
pairwise.t.test(data$DUREZA[which(data$MORTERO=="100CW")], 
                data$EDAD[which(data$MORTERO=="100CW")], p.adjust.method = "bonferroni")


####
# Compression test analysis
####

data_comp <- read.csv("Compresion.csv", sep = ";")
data_comp$MORTERO <- as.factor(data_comp$MORTERO)
data_comp$EDAD <- as.factor(data_comp$EDAD)

# Plot figure with dosage and age as factors

ggline(data_comp, x = "MORTERO", y = "COMPRESION", color = "EDAD", add=c("mean_se","jitter")) +
  labs(y = "Compressive strength (N/mm2)") + labs(x = "") + 
  scale_color_discrete(name = "Age", labels = c("28 days", "90 days", "180 days")) + 
  scale_x_discrete(labels=c('REF', 'CW-25%', 'CW-50%', 'CW-75%', 'CW-100%'))

# Two-way ANOVA with interaction

res.aov.c <- aov(COMPRESION ~ MORTERO*EDAD, data = data_comp)
summary(res.aov.c)

# Test of the normality of the residuals

aov_residuals.c <- residuals(object = res.aov.c)
ks.test(unique(aov_residuals.c), pnorm, 0, sd(aov_residuals.c))

# Test for homoscedasticity (with a logarithmic transformation)

leveneTest(COMPRESION ~ MORTERO*EDAD, data = data_comp)
leveneTest(log(COMPRESION) ~ MORTERO*EDAD, data = data_comp)

# Bootstrap to account for the absence of normality, with log transformation

data_comp_log <- data_comp
data_comp_log$COMPRESION <- log(data_comp_log$COMPRESION)

gpTwoWay(COMPRESION ~ MORTERO*EDAD, data = data_comp_log, method = "gPB")

# Pairwise comparisons: for fixed ages

pairwise.t.test(data_comp_log$COMPRESION[which(data_comp_log$EDAD=="028D")], 
                data_comp_log$MORTERO[which(data_comp_log$EDAD=="028D")], p.adjust.method = "bonferroni")
pairwise.t.test(data_comp_log$COMPRESION[which(data_comp_log$EDAD=="090D")], 
                data_comp_log$MORTERO[which(data_comp_log$EDAD=="090D")], p.adjust.method = "bonferroni")
pairwise.t.test(data_comp_log$COMPRESION[which(data_comp_log$EDAD=="180D")], 
                data_comp_log$MORTERO[which(data_comp_log$EDAD=="180D")], p.adjust.method = "bonferroni")

# Pairwise comparisons: for fixed dosages

pairwise.t.test(data_comp_log$COMPRESION[which(data_comp_log$MORTERO=="000CW")], 
                data_comp_log$EDAD[which(data_comp_log$MORTERO=="000CW")], p.adjust.method = "bonferroni")
pairwise.t.test(data_comp_log$COMPRESION[which(data_comp_log$MORTERO=="025CW")], 
                data_comp_log$EDAD[which(data_comp_log$MORTERO=="025CW")], p.adjust.method = "bonferroni")
pairwise.t.test(data_comp_log$COMPRESION[which(data_comp_log$MORTERO=="050CW")], 
                data_comp_log$EDAD[which(data_comp_log$MORTERO=="050CW")], p.adjust.method = "bonferroni")
pairwise.t.test(data_comp_log$COMPRESION[which(data_comp_log$MORTERO=="075CW")], 
                data_comp_log$EDAD[which(data_comp_log$MORTERO=="075CW")], p.adjust.method = "bonferroni")
pairwise.t.test(data_comp_log$COMPRESION[which(data_comp_log$MORTERO=="100CW")], 
                data_comp_log$EDAD[which(data_comp_log$MORTERO=="100CW")], p.adjust.method = "bonferroni")


####
# Flexural test analysis
####

data_flex <- read.csv("Flexion.csv", sep = ";")
data_flex$MORTERO <- as.factor(data_flex$MORTERO)
data_flex$EDAD <- as.factor(data_flex$EDAD)

# Plot figure with dosage and age as factors

ggline(data_flex, x = "MORTERO", y = "FLEXION", color = "EDAD", add=c("mean_se","jitter")) +
  labs(y = "Flexural strength (N/mm2)") + labs(x = "") + 
  scale_color_discrete(name = "Age", labels = c("28 days", "90 days", "180 days")) + 
  scale_x_discrete(labels=c('REF', 'CW-25%', 'CW-50%', 'CW-75%', 'CW-100%'))

# Two-way ANOVA with interaction

res.aov.f <- aov(FLEXION ~ MORTERO*EDAD, data = data_flex)
summary(res.aov.f)

# Test of the normality of the residuals

aov_residuals.f <- residuals(object = res.aov.f)
ks.test(unique(aov_residuals.f), pnorm, 0, sd(aov_residuals.f))

# Test for homoscedasticity (with a logarithmic transformation)

leveneTest(FLEXION ~ MORTERO*EDAD, data = data_flex)

# Pairwise comparisons: for fixed ages

pairwise.t.test(data_flex$FLEXION[which(data_flex$EDAD=="028D")], 
                data_flex$MORTERO[which(data_flex$EDAD=="028D")], p.adjust.method = "bonferroni")
pairwise.t.test(data_flex$FLEXION[which(data_flex$EDAD=="090D")], 
                data_flex$MORTERO[which(data_flex$EDAD=="090D")], p.adjust.method = "bonferroni")
pairwise.t.test(data_flex$FLEXION[which(data_flex$EDAD=="180D")], 
                data_flex$MORTERO[which(data_flex$EDAD=="180D")], p.adjust.method = "bonferroni")

# Pairwise comparisons: for fixed dosages

pairwise.t.test(data_flex$FLEXION[which(data_flex$MORTERO=="000CW")], 
                data_flex$EDAD[which(data_flex$MORTERO=="000CW")], p.adjust.method = "bonferroni")
pairwise.t.test(data_flex$FLEXION[which(data_flex$MORTERO=="025CW")], 
                data_flex$EDAD[which(data_flex$MORTERO=="025CW")], p.adjust.method = "bonferroni")
pairwise.t.test(data_flex$FLEXION[which(data_flex$MORTERO=="050CW")], 
                data_flex$EDAD[which(data_flex$MORTERO=="050CW")], p.adjust.method = "bonferroni")
pairwise.t.test(data_flex$FLEXION[which(data_flex$MORTERO=="075CW")], 
                data_flex$EDAD[which(data_flex$MORTERO=="075CW")], p.adjust.method = "bonferroni")
pairwise.t.test(data_flex$FLEXION[which(data_flex$MORTERO=="100CW")], 
                data_flex$EDAD[which(data_flex$MORTERO=="100CW")], p.adjust.method = "bonferroni")
