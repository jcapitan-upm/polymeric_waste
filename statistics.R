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

# Bootstrap to account for heteroscedasticity and the lack of normality

library(twowaytests)
out <- gpTwoWay(DUREZA ~ MORTERO*EDAD, data = data, method = "gPB")

# Bootstrapping considering only separate measurements for different test tubes
# Loop for 100 samples of the measurement taken in each test tube

subsets <- c()

for (dosage in c("000CW","025CW","050CW","075CW","100CW")){
  for (age in c("028D","090D","180D")){
    
    subsets <- rbind(subsets, data$DUREZA[which(data$MORTERO==dosage & data$EDAD==age)])
  }
}

count <- rep(0,3)

for (i in 1:100){
  
  print(i)
  data_sampled <- c()

  index <- 1
  
  for (dosage in c("000CW","025CW","050CW","075CW","100CW")){
    for (age in c("028D","090D","180D")){
      
      sdval <- 0
      
      while (sdval == 0) { # To avoid three equal values at the same level
        
        s1 <- dqsample(subsets[index, 1:10], size = 1)
        data_sampled <- rbind(data_sampled, c(dosage, age, s1))
        s2 <- dqsample(subsets[index, 11:20], size = 1)
        data_sampled <- rbind(data_sampled, c(dosage, age, s2))
        s3 <- dqsample(subsets[index, 21:30], size = 1)
        data_sampled <- rbind(data_sampled, c(dosage, age, s3))
        L <- length(data_sampled[,3])
        sdval <- sd(data_sampled[(L-2):L,3])
        if (sdval > 0){
          break
        }
      }
      
      index <- index + 1
    }
  }
  
  data_sampled1 <- as.data.frame(data_sampled)
  names(data_sampled1) <- c("MORTERO","EDAD","DUREZA")
  data_sampled1$MORTERO <- as.factor(data_sampled1$MORTERO)
  data_sampled1$EDAD <- as.factor(data_sampled1$EDAD)
  data_sampled1$DUREZA <- as.numeric(data_sampled1$DUREZA)
  
  out <- gpTwoWay(DUREZA ~ MORTERO*EDAD, data = data_sampled1, method = "gPB")
  
  # Count cases in which effects are not rejected
  nrej <- which(out$output$P.value>0.05)
  
  if (length(nrej) > 0) {
    count[nrej] <- count[nrej] + 1
  }
}

# Number of cases (out of 100) that one of the three terms in the model is non-significant
count

# Print the last sample

ggline(data_sampled1, x = "MORTERO", y = "DUREZA", color = "EDAD", add=c("mean_se","jitter")) +
  labs(y = "Shore D surface hardness") + labs(x = "") + 
  scale_color_discrete(name = "Age", labels = c("28 days", "90 days", "180 days")) + 
  scale_x_discrete(labels=c('REF', 'CW-25%', 'CW-50%', 'CW-75%', 'CW-100%'))

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
shapiro.test(aov_residuals.c)

# Test for homoscedasticity (with a logarithmic transformation)

leveneTest(COMPRESION ~ MORTERO*EDAD, data = data_comp)
leveneTest(log(COMPRESION) ~ MORTERO*EDAD, data = data_comp)

# Bootstrap to account for the absence of normality, with log transformation

data_comp_log <- data_comp
data_comp_log$COMPRESION <- log(data_comp_log$COMPRESION)

gpTwoWay(COMPRESION ~ MORTERO*EDAD, data = data_comp_log, method = "gPB")

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
shapiro.test(aov_residuals.f)

# Test for homoscedasticity (with a logarithmic transformation)

leveneTest(FLEXION ~ MORTERO*EDAD, data = data_flex)
