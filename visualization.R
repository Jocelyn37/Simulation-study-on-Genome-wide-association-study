###### Visualization and analysis ############
library(ggplot2)
#library(ggbump)
library(dplyr)

################# Type1error ###################################################
#### sample 350
setwd("/Users/jocelyn/Desktop/BiolabIntern/results+graph/type1error/sample350/Simulation_sig_SPA_1billion_5e-8")
type1error_SPA <- read.csv("Type1_error.csv", header = TRUE)
setwd("/Users/jocelyn/Desktop/BiolabIntern/results+graph/type1error/sample350/Simulation_sig_withoutSPA")
type1error_4tests <- read.csv("Type1_error.csv", header = TRUE)

type1error <- cbind(type1error_4tests, type1error_SPA$SPACox)
colnames(type1error)
names(type1error)[names(type1error) == "type1error_SPA$SPACox"] <- "SPACox"
type1error$MAF <- factor(type1error$MAF)
type1error$X <- NULL
type1error <- type1error[,c(1,2,3,4,5,6,7,9,8)]

type1error_mod <- type1error %>%
  mutate(Shape_para = recode(Shape_para, "1.0" = "Shape_par = 1", "1.5" = "Shape_par = 1.5"))

type1error_mod <- type1error_mod %>%
  mutate(Cens_Prob = recode(Cens_Prob, "0.25" = "CR = 0.25", "0.50" = "CR = 0.50"))

t1e_plot <- ggplot(data = type1error_mod, aes(x = MAF)) +
  geom_line(aes(y = Score, group = 1, color = "Score"), size = 1, linetype = "dashed") +
  geom_point(aes(y = Score, color = "Score"), size = 3) +
  geom_line(aes(y = Wald, group = 1, color = "Wald"), size = 1, linetype = "dashed") +
  geom_point(aes(y = Wald, color = "Wald"), size = 3) +
  geom_line(aes(y = LikeRatio, group = 1, color = "LikeRatio"), size = 1, linetype = "dashed") +
  geom_point(aes(y = LikeRatio, color = "LikeRatio"), size = 3) +
  geom_line(aes(y = Firth, group = 1, color = "Firth"), size = 1, linetype = "dashed") +
  geom_point(aes(y = Firth, color = "Firth"), size = 3) +
  geom_line(aes(y = SPACox, group = 1, color = "SPACox"), size = 1, linetype = "dashed") + 
  geom_point(aes(y = SPACox, color = "SPACox"), size = 3) + scale_y_log10(limits = c(1e-9, 1e-4)) +
  geom_hline(yintercept = 5e-8, linetype = "dashed") +
  scale_color_manual(name = "Methods", values = c("Score" = "steelblue", "Wald" = "hotpink", "LikeRatio" = "seagreen", 
                                                  "Firth" = "orange", "SPACox" = "lightcoral"))

t1e_plot + facet_grid(Shape_para ~ Cens_Prob) + 
  ylab("Type 1 Error Rates") + ggtitle("Sample Size = 350") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5))

# Performance
#1
perfor <- type1error
perfor$X <- NULL
perfor[5:9] <- abs(log(perfor[5:9]/5e-8))

#2
perfor2<- type1error
perfor2$X <- NULL
perfor2[5:9] <- log(perfor2[5:9]/5e-8)
perfor2[5:9][perfor2[5:9]<0] <- 0

#3 valid/non-valid
perfor3 <- type1error

for(i in 1:20){
  for(j in 5:9){
    prop_result <- prop.test(perfor3[i,j]*1e9, 1e9, conf.level = 0.95)
    left1 <- prop_result$conf.int[1]
    right1 <- prop_result$conf.int[2]
    perfor3[i,j] <- as.numeric(5.5e-8 >left1)
  }
}

####### Read power data and data precossing(below) and get powerCox_1, powerCox_2(shape1,1.5)
#sample size 350
MAF <- c(0.5, 0.2, 0.1, 0.05, 0.02)
SHAPE <- c(1, 1.5)
CENS_PROB <- c(0.25, 0.5)

perfor3_mean <- perfor3
rank_power <- perfor3

#shape=1.0
powerCox_1$X <- NULL
for(cens in CENS_PROB){
  for(maf in MAF){
    power_per <- powerCox_1[(powerCox_1$MAF==maf)&(powerCox_1$Cens_Prob==cens),]
    power_per <- power_per[power_per$coef != 0.0,]
    power_mean <- colMeans(power_per[5:9])
    perfor3_mean[(perfor3_mean$MAF==maf)&(perfor3_mean$Cens_Prob==cens)&(perfor3_mean$Shape_para==1.0),5:9] <- power_mean
    
    rank_power[(rank_power$MAF==maf)&(rank_power$Cens_Prob==cens)&(rank_power$Shape_para==1.0),5:9] <- rank(-power_mean, ties.method = "min")
  }
}

#shape=1.5
powerCox_2$X <- NULL
for(cens in CENS_PROB){
  for(maf in MAF){
    power_per <- powerCox_2[(powerCox_1$MAF==maf)&(powerCox_1$Cens_Prob==cens),]
    power_per <- power_per[power_per$coef != 0.0,]
    power_mean <- colMeans(power_per[5:9])
    perfor3_mean[(perfor3_mean$MAF==maf)&(perfor3_mean$Cens_Prob==cens)&(perfor3_mean$Shape_para==1.5),5:9] <- power_mean
    
    rank_power[(rank_power$MAF==maf)&(rank_power$Cens_Prob==cens)&(rank_power$Shape_para==1.5),5:9] <- rank(-power_mean, ties.method = "min")
  }
}


write.csv(perfor3, "/Users/jocelyn/Desktop/BiolabIntern/bump_chart/350valid.csv", row.names = TRUE)
write.csv(rank_power, "/Users/jocelyn/Desktop/BiolabIntern/bump_chart/350rank.csv", row.names = TRUE)

# rank_power <- perfor3_mean
# scenarios<- rank_power[(rank_power$Shape_para==1.0)&(rank_power$Cens_Prob==0.25),]
# rank_power[5:9] <- t(apply(-scenarios[5:9], 1, rank, ties.method='min'))


#### sample 800
setwd("/Users/jocelyn/Desktop/BiolabIntern/results+graph/type1error/sample800/Simulation_sig_SPA_1billion_2")
type1error_SPA <- read.csv("Type1_error.csv", header = TRUE)
setwd("/Users/jocelyn/Desktop/BiolabIntern/results+graph/type1error/sample800/Simulation_sig_withoutSPA_2")
type1error_4tests <- read.csv("Type1_error.csv", header = TRUE)

type1error <- cbind(type1error_4tests, type1error_SPA$SPACox)
colnames(type1error)
names(type1error)[names(type1error) == "type1error_SPA$SPACox"] <- "SPACox"
type1error$MAF <- factor(type1error$MAF)
type1error$X <- NULL
type1error <- type1error[,c(1,2,3,4,5,6,7,9,8)]

type1error_mod <- type1error %>%
  mutate(Shape_para = recode(Shape_para, "1.0" = "Shape_par = 1", "1.5" = "Shape_par = 1.5"))

type1error_mod <- type1error_mod %>%
  mutate(Cens_Prob = recode(Cens_Prob, "0.25" = "CR = 0.25", "0.50" = "CR = 0.50"))


t1e_plot + facet_grid(Shape_para ~ Cens_Prob) + 
  ylab("Type 1 Error Rates") + ggtitle("Sample Size = 800") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5))



#3 valid/non-valid
perfor3 <- type1error

for(i in 1:20){
  for(j in 5:9){
    prop_result <- prop.test(perfor3[i,j]*1e9, 1e9, conf.level = 0.95)
    left1 <- prop_result$conf.int[1]
    right1 <- prop_result$conf.int[2]
    perfor3[i,j] <- as.numeric(5.5e-8 >left1)
  }
}

####### Read power data and data precossing(below) and get powerCox_1, powerCox_2(shape1,1.5)
#sample size 800
MAF <- c(0.5, 0.2, 0.1, 0.05, 0.02)
SHAPE <- c(1, 1.5)
CENS_PROB <- c(0.25, 0.5)

perfor3_mean <- perfor3
rank_power <- perfor3

#shape=1.0
powerCox_1$X <- NULL
for(cens in CENS_PROB){
  for(maf in MAF){
    power_per <- powerCox_1[(powerCox_1$MAF==maf)&(powerCox_1$Cens_Prob==cens),]
    power_per <- power_per[power_per$coef != 0.0,]
    power_mean <- colMeans(power_per[5:9])
    perfor3_mean[(perfor3_mean$MAF==maf)&(perfor3_mean$Cens_Prob==cens)&(perfor3_mean$Shape_para==1.0),5:9] <- power_mean
    
    rank_power[(rank_power$MAF==maf)&(rank_power$Cens_Prob==cens)&(rank_power$Shape_para==1.0),5:9] <- rank(-power_mean, ties.method = "min")
  }
}

#shape=1.5
powerCox_2$X <- NULL
for(cens in CENS_PROB){
  for(maf in MAF){
    power_per <- powerCox_2[(powerCox_1$MAF==maf)&(powerCox_1$Cens_Prob==cens),]
    power_per <- power_per[power_per$coef != 0.0,]
    power_mean <- colMeans(power_per[5:9])
    perfor3_mean[(perfor3_mean$MAF==maf)&(perfor3_mean$Cens_Prob==cens)&(perfor3_mean$Shape_para==1.5),5:9] <- power_mean
    
    rank_power[(rank_power$MAF==maf)&(rank_power$Cens_Prob==cens)&(rank_power$Shape_para==1.5),5:9] <- rank(-power_mean, ties.method = "min")
  }
}


write.csv(perfor3, "/Users/jocelyn/Desktop/BiolabIntern/bump_chart/800valid.csv", row.names = TRUE)
write.csv(rank_power, "/Users/jocelyn/Desktop/BiolabIntern/bump_chart/800rank.csv", row.names = TRUE)




#### sample 1600
setwd("/Users/jocelyn/Desktop/BiolabIntern/results+graph/type1error/sample1600")
type1error <- read.csv("Type1_error.csv", header = TRUE)
type1error$MAF <- factor(type1error$MAF)

type1error_mod <- type1error %>%
  mutate(Shape_para = recode(Shape_para, "1.0" = "Shape_par = 1", "1.5" = "Shape_par = 1.5"))

type1error_mod <- type1error_mod %>%
  mutate(Cens_Prob = recode(Cens_Prob, "0.25" = "CR = 0.25", "0.50" = "CR = 0.50"))


t1e_plot + facet_grid(Shape_para ~ Cens_Prob) + 
  ylab("Type 1 Error Rates") + ggtitle("Sample Size = 1600") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5))



#3 valid/non-valid
perfor3 <- type1error
perfor3$X <- NULL

for(i in 1:20){
  for(j in 5:9){
    prop_result <- prop.test(perfor3[i,j]*1e9, 1e9, conf.level = 0.95)
    left1 <- prop_result$conf.int[1]
    right1 <- prop_result$conf.int[2]
    perfor3[i,j] <- as.numeric(5.5e-8 >left1)
  }
}

####### Read power data and data precossing(below) and get powerCox_1, powerCox_2(shape1,1.5)
#sample size 1600
MAF <- c(0.5, 0.2, 0.1, 0.05, 0.02)
SHAPE <- c(1, 1.5)
CENS_PROB <- c(0.25, 0.5)

perfor3_mean <- perfor3
rank_power <- perfor3

#shape=1.0
powerCox_1$X <- NULL
for(cens in CENS_PROB){
  for(maf in MAF){
    power_per <- powerCox_1[(powerCox_1$MAF==maf)&(powerCox_1$Cens_Prob==cens),]
    power_per <- power_per[power_per$coef != 0.0,]
    power_mean <- colMeans(power_per[5:9])
    perfor3_mean[(perfor3_mean$MAF==maf)&(perfor3_mean$Cens_Prob==cens)&(perfor3_mean$Shape_para==1.0),5:9] <- power_mean
    
    rank_power[(rank_power$MAF==maf)&(rank_power$Cens_Prob==cens)&(rank_power$Shape_para==1.0),5:9] <- rank(-power_mean, ties.method = "min")
  }
}

#shape=1.5
powerCox_2$X <- NULL
for(cens in CENS_PROB){
  for(maf in MAF){
    power_per <- powerCox_2[(powerCox_1$MAF==maf)&(powerCox_1$Cens_Prob==cens),]
    power_per <- power_per[power_per$coef != 0.0,]
    power_mean <- colMeans(power_per[5:9])
    perfor3_mean[(perfor3_mean$MAF==maf)&(perfor3_mean$Cens_Prob==cens)&(perfor3_mean$Shape_para==1.5),5:9] <- power_mean
    
    rank_power[(rank_power$MAF==maf)&(rank_power$Cens_Prob==cens)&(rank_power$Shape_para==1.5),5:9] <- rank(-power_mean, ties.method = "min")
  }
}


write.csv(perfor3, "/Users/jocelyn/Desktop/BiolabIntern/bump_chart/1600valid.csv", row.names = TRUE)
write.csv(rank_power, "/Users/jocelyn/Desktop/BiolabIntern/bump_chart/1600rank.csv", row.names = TRUE)


###################### Power ###################################################
###########nogroup
######350
COEF <- seq(-1, 1, by = 0.1)
for(coef in COEF){
  filepath <- paste('/Users/jocelyn/Desktop/BiolabIntern/Project/Simulation code/power_results/result_nogroup/coef_', coef, sep = "")
  filename <- paste(filepath, '/powerCox.csv', sep = "")
  powerCox <- read.csv(filename, header = TRUE)
  powerCox$coef <- rep(coef, 20)
  
  if(coef == -1){
    powerCox_1 <- powerCox[powerCox$Shape_para == 1.0, ]
    powerCox_2 <- powerCox[powerCox$Shape_para == 1.5, ]
  }
  else{
    powerCox_1_new <- powerCox[powerCox$Shape_para == 1.0, ]
    powerCox_2_new <- powerCox[powerCox$Shape_para == 1.5, ]
    powerCox_1 <- rbind(powerCox_1, powerCox_1_new)
    powerCox_2 <- rbind(powerCox_2, powerCox_1_new)
  }
}

### Shape para = 1.0
colnames(powerCox_1)
powerCox_1$MAF <- factor(powerCox_1$MAF)

# powerCox_1_mod <- powerCox_1 %>%
#   mutate(Shape_para = recode(Shape_para, "1.0" = "Shape_par = 1", "1.5" = "Shape_par = 1.5"))

powerCox_1_mod <- powerCox_1 %>%
  mutate(Cens_Prob = recode(Cens_Prob, "0.25" = "CR = 0.25", "0.50" = "CR = 0.50"))
powerCox_1_mod <- powerCox_1_mod %>%
  mutate(MAF = recode(MAF, "0.02" = "MAF = 0.02", "0.05" = "MAF = 0.05", "0.1" = "MAF = 0.1","0.2" = "MAF = 0.2","0.5" = "MAF = 0.5"))

t1e_plot <- ggplot(data = powerCox_1_mod, aes(x = coef)) +
  geom_line(aes(y = Score, group = 1, color = "Score"), size = 0.7, linetype = "dashed") +
  geom_point(aes(y = Score, color = "Score"), size = 2) +
  geom_line(aes(y = Wald, group = 1, color = "Wald"), size = 0.7, linetype = "dashed") +
  geom_point(aes(y = Wald, color = "Wald"), size = 2) +
  geom_line(aes(y = LikeRatio, group = 1, color = "LikeRatio"), size = 0.7, linetype = "dashed") +
  geom_point(aes(y = LikeRatio, color = "LikeRatio"), size = 2) +
  geom_line(aes(y = Firth, group = 1, color = "Firth"), size = 0.7, linetype = "dashed") +
  geom_point(aes(y = Firth, color = "Firth"), size = 2) +
  geom_line(aes(y = SPACox, group = 1, color = "SPACox"), size = 0.7, linetype = "dashed") + 
  geom_point(aes(y = SPACox, color = "SPACox"), size = 2) + 
  ylim(0,1) +
  scale_color_manual(name = "Methods", values = c("Score" = "steelblue", "Wald" = "hotpink", "LikeRatio" = "seagreen", 
                                                  "Firth" = "orange", "SPACox" = "lightcoral"))

t1e_plot + facet_grid(Cens_Prob ~ MAF) + 
  ylab("Power") + ggtitle("Shape para = 1.0") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5))


### Shape para = 1.5
colnames(powerCox_2)
powerCox_2$MAF <- factor(powerCox_2$MAF)

# powerCox_1_mod <- powerCox_1 %>%
#   mutate(Shape_para = recode(Shape_para, "1.0" = "Shape_par = 1", "1.5" = "Shape_par = 1.5"))

powerCox_2_mod <- powerCox_2 %>%
  mutate(Cens_Prob = recode(Cens_Prob, "0.25" = "CR = 0.25", "0.50" = "CR = 0.50"))
powerCox_2_mod <- powerCox_2_mod %>%
  mutate(MAF = recode(MAF, "0.02" = "MAF = 0.02", "0.05" = "MAF = 0.05", "0.1" = "MAF = 0.1","0.2" = "MAF = 0.2","0.5" = "MAF = 0.5"))

t1e_plot <- ggplot(data = powerCox_2_mod, aes(x = coef)) +
  geom_line(aes(y = Score, group = 1, color = "Score"), size = 0.7, linetype = "dashed") +
  geom_point(aes(y = Score, color = "Score"), size = 2) +
  geom_line(aes(y = Wald, group = 1, color = "Wald"), size = 0.7, linetype = "dashed") +
  geom_point(aes(y = Wald, color = "Wald"), size = 2) +
  geom_line(aes(y = LikeRatio, group = 1, color = "LikeRatio"), size = 0.7, linetype = "dashed") +
  geom_point(aes(y = LikeRatio, color = "LikeRatio"), size = 2) +
  geom_line(aes(y = Firth, group = 1, color = "Firth"), size = 0.7, linetype = "dashed") +
  geom_point(aes(y = Firth, color = "Firth"), size = 2) +
  geom_line(aes(y = SPACox, group = 1, color = "SPACox"), size = 0.7, linetype = "dashed") + 
  geom_point(aes(y = SPACox, color = "SPACox"), size = 2) + 
  ylim(0,1) +
  scale_color_manual(name = "Methods", values = c("Score" = "steelblue", "Wald" = "hotpink", "LikeRatio" = "seagreen", 
                                                  "Firth" = "orange", "SPACox" = "lightcoral"))

t1e_plot + facet_grid(Cens_Prob ~ MAF) + 
  ylab("Power") + ggtitle("Shape para = 1.5") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5))

######800
COEF <- seq(-1, 1, by = 0.1)
for(coef in COEF){
  filepath <- paste('/Users/jocelyn/Desktop/BiolabIntern/results+graph/power/result_800/coef_', coef, sep = "")
  filename <- paste(filepath, '/powerCox.csv', sep = "")
  powerCox <- read.csv(filename, header = TRUE)
  powerCox$coef <- rep(coef, 20)
  
  if(coef == -1){
    powerCox_1 <- powerCox[powerCox$Shape_para == 1.0, ]
    powerCox_2 <- powerCox[powerCox$Shape_para == 1.5, ]
  }
  else{
    powerCox_1_new <- powerCox[powerCox$Shape_para == 1.0, ]
    powerCox_2_new <- powerCox[powerCox$Shape_para == 1.5, ]
    powerCox_1 <- rbind(powerCox_1, powerCox_1_new)
    powerCox_2 <- rbind(powerCox_2, powerCox_1_new)
  }
}

### Shape para = 1.0
colnames(powerCox_1)
# powerCox_1 <- powerCox_1[powerCox_1$MAF==0.05,]
powerCox_1$MAF <- factor(powerCox_1$MAF)

# powerCox_1_mod <- powerCox_1 %>%
#   mutate(Shape_para = recode(Shape_para, "1.0" = "Shape_par = 1", "1.5" = "Shape_par = 1.5"))

powerCox_1_mod <- powerCox_1 %>%
  mutate(Cens_Prob = recode(Cens_Prob, "0.25" = "CR = 0.25", "0.50" = "CR = 0.50"))
powerCox_1_mod <- powerCox_1_mod %>%
  mutate(MAF = recode(MAF, "0.02" = "MAF = 0.02", "0.05" = "MAF = 0.05", "0.1" = "MAF = 0.1","0.2" = "MAF = 0.2","0.5" = "MAF = 0.5"))

# powerCox_1_mod <- powerCox_1_mod %>%
#   mutate(MAF = recode(MAF,  "0.05" = "MAF = 0.05"))

t1e_plot <- ggplot(data = powerCox_1_mod, aes(x = coef)) +
  geom_line(aes(y = Score, group = 1, color = "Score"), size = 0.7, linetype = "dashed") +
  geom_point(aes(y = Score, color = "Score"), size = 2) +
  geom_line(aes(y = Wald, group = 1, color = "Wald"), size = 0.7, linetype = "dashed") +
  geom_point(aes(y = Wald, color = "Wald"), size = 2) +
  geom_line(aes(y = LikeRatio, group = 1, color = "LikeRatio"), size = 0.7, linetype = "dashed") +
  geom_point(aes(y = LikeRatio, color = "LikeRatio"), size = 2) +
  geom_line(aes(y = Firth, group = 1, color = "Firth"), size = 0.7, linetype = "dashed") +
  geom_point(aes(y = Firth, color = "Firth"), size = 2) +
  geom_line(aes(y = SPACox, group = 1, color = "SPACox"), size = 0.7, linetype = "dashed") + 
  geom_point(aes(y = SPACox, color = "SPACox"), size = 2) + 
  ylim(0,1) +
  scale_color_manual(name = "Methods", values = c("Score" = "steelblue", "Wald" = "hotpink", "LikeRatio" = "seagreen", 
                                                  "Firth" = "orange", "SPACox" = "lightcoral"))

t1e_plot + facet_grid(Cens_Prob ~ MAF) + 
  ylab("Power") + ggtitle("Shape para = 1.0") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5))

### Shape para = 1.5
colnames(powerCox_2)
powerCox_2$MAF <- factor(powerCox_2$MAF)

# powerCox_1_mod <- powerCox_1 %>%
#   mutate(Shape_para = recode(Shape_para, "1.0" = "Shape_par = 1", "1.5" = "Shape_par = 1.5"))

powerCox_2_mod <- powerCox_2 %>%
  mutate(Cens_Prob = recode(Cens_Prob, "0.25" = "CR = 0.25", "0.50" = "CR = 0.50"))
powerCox_2_mod <- powerCox_2_mod %>%
  mutate(MAF = recode(MAF, "0.02" = "MAF = 0.02", "0.05" = "MAF = 0.05", "0.1" = "MAF = 0.1","0.2" = "MAF = 0.2","0.5" = "MAF = 0.5"))

t1e_plot <- ggplot(data = powerCox_2_mod, aes(x = coef)) +
  geom_line(aes(y = Score, group = 1, color = "Score"), size = 0.7, linetype = "dashed") +
  geom_point(aes(y = Score, color = "Score"), size = 2) +
  geom_line(aes(y = Wald, group = 1, color = "Wald"), size = 0.7, linetype = "dashed") +
  geom_point(aes(y = Wald, color = "Wald"), size = 2) +
  geom_line(aes(y = LikeRatio, group = 1, color = "LikeRatio"), size = 0.7, linetype = "dashed") +
  geom_point(aes(y = LikeRatio, color = "LikeRatio"), size = 2) +
  geom_line(aes(y = Firth, group = 1, color = "Firth"), size = 0.7, linetype = "dashed") +
  geom_point(aes(y = Firth, color = "Firth"), size = 2) +
  geom_line(aes(y = SPACox, group = 1, color = "SPACox"), size = 0.7, linetype = "dashed") + 
  geom_point(aes(y = SPACox, color = "SPACox"), size = 2) + 
  ylim(0,1) +
  scale_color_manual(name = "Methods", values = c("Score" = "steelblue", "Wald" = "hotpink", "LikeRatio" = "seagreen", 
                                                  "Firth" = "orange", "SPACox" = "lightcoral"))

t1e_plot + facet_grid(Cens_Prob ~ MAF) + 
  ylab("Power") + ggtitle("Shape para = 1.5") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5))


######1600
COEF <- seq(-1, 1, by = 0.1)
for(coef in COEF){
  filepath <- paste('/Users/jocelyn/Desktop/BiolabIntern/results+graph/power/result_1600/coef_', coef, sep = "")
  filename <- paste(filepath, '/powerCox.csv', sep = "")
  powerCox <- read.csv(filename, header = TRUE)
  powerCox$coef <- rep(coef, 20)
  
  if(coef == -1){
    powerCox_1 <- powerCox[powerCox$Shape_para == 1.0, ]
    powerCox_2 <- powerCox[powerCox$Shape_para == 1.5, ]
  }
  else{
    powerCox_1_new <- powerCox[powerCox$Shape_para == 1.0, ]
    powerCox_2_new <- powerCox[powerCox$Shape_para == 1.5, ]
    powerCox_1 <- rbind(powerCox_1, powerCox_1_new)
    powerCox_2 <- rbind(powerCox_2, powerCox_1_new)
  }
}

### Shape para = 1.0
colnames(powerCox_1)
powerCox_1$MAF <- factor(powerCox_1$MAF)

# powerCox_1_mod <- powerCox_1 %>%
#   mutate(Shape_para = recode(Shape_para, "1.0" = "Shape_par = 1", "1.5" = "Shape_par = 1.5"))

powerCox_1_mod <- powerCox_1 %>%
  mutate(Cens_Prob = recode(Cens_Prob, "0.25" = "CR = 0.25", "0.50" = "CR = 0.50"))
powerCox_1_mod <- powerCox_1_mod %>%
  mutate(MAF = recode(MAF, "0.02" = "MAF = 0.02", "0.05" = "MAF = 0.05", "0.1" = "MAF = 0.1","0.2" = "MAF = 0.2","0.5" = "MAF = 0.5"))

t1e_plot <- ggplot(data = powerCox_1_mod, aes(x = coef)) +
  geom_line(aes(y = Score, group = 1, color = "Score"), size = 0.7, linetype = "dashed") +
  geom_point(aes(y = Score, color = "Score"), size = 2) +
  geom_line(aes(y = Wald, group = 1, color = "Wald"), size = 0.7, linetype = "dashed") +
  geom_point(aes(y = Wald, color = "Wald"), size = 2) +
  geom_line(aes(y = LikeRatio, group = 1, color = "LikeRatio"), size = 0.7, linetype = "dashed") +
  geom_point(aes(y = LikeRatio, color = "LikeRatio"), size = 2) +
  geom_line(aes(y = Firth, group = 1, color = "Firth"), size = 0.7, linetype = "dashed") +
  geom_point(aes(y = Firth, color = "Firth"), size = 2) +
  geom_line(aes(y = SPACox, group = 1, color = "SPACox"), size = 0.7, linetype = "dashed") + 
  geom_point(aes(y = SPACox, color = "SPACox"), size = 2) + 
  ylim(0,1) +
  scale_color_manual(name = "Methods", values = c("Score" = "steelblue", "Wald" = "hotpink", "LikeRatio" = "seagreen", 
                                                  "Firth" = "orange", "SPACox" = "lightcoral"))

t1e_plot + facet_grid(Cens_Prob ~ MAF) + 
  ylab("Power") + ggtitle("Shape para = 1.0") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5))

### Shape para = 1.5
colnames(powerCox_2)
powerCox_2$MAF <- factor(powerCox_2$MAF)

# powerCox_1_mod <- powerCox_1 %>%
#   mutate(Shape_para = recode(Shape_para, "1.0" = "Shape_par = 1", "1.5" = "Shape_par = 1.5"))

powerCox_2_mod <- powerCox_2 %>%
  mutate(Cens_Prob = recode(Cens_Prob, "0.25" = "CR = 0.25", "0.50" = "CR = 0.50"))
powerCox_2_mod <- powerCox_2_mod %>%
  mutate(MAF = recode(MAF, "0.02" = "MAF = 0.02", "0.05" = "MAF = 0.05", "0.1" = "MAF = 0.1","0.2" = "MAF = 0.2","0.5" = "MAF = 0.5"))

t1e_plot <- ggplot(data = powerCox_2_mod, aes(x = coef)) +
  geom_line(aes(y = Score, group = 1, color = "Score"), size = 0.7, linetype = "dashed") +
  geom_point(aes(y = Score, color = "Score"), size = 2) +
  geom_line(aes(y = Wald, group = 1, color = "Wald"), size = 0.7, linetype = "dashed") +
  geom_point(aes(y = Wald, color = "Wald"), size = 2) +
  geom_line(aes(y = LikeRatio, group = 1, color = "LikeRatio"), size = 0.7, linetype = "dashed") +
  geom_point(aes(y = LikeRatio, color = "LikeRatio"), size = 2) +
  geom_line(aes(y = Firth, group = 1, color = "Firth"), size = 0.7, linetype = "dashed") +
  geom_point(aes(y = Firth, color = "Firth"), size = 2) +
  geom_line(aes(y = SPACox, group = 1, color = "SPACox"), size = 0.7, linetype = "dashed") + 
  geom_point(aes(y = SPACox, color = "SPACox"), size = 2) + 
  ylim(0,1) +
  scale_color_manual(name = "Methods", values = c("Score" = "steelblue", "Wald" = "hotpink", "LikeRatio" = "seagreen", 
                                                  "Firth" = "orange", "SPACox" = "lightcoral"))

t1e_plot + facet_grid(Cens_Prob ~ MAF) + 
  ylab("Power") + ggtitle("Shape para = 1.5") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5))







################################################################################
###########group
COEF <- seq(-1, 1, by = 0.1)
for(coef in COEF){
  filepath <- paste('/Users/jocelyn/Desktop/Biolab/Project/Simulation code/power_results/results_group/coef_', coef, sep = "")
  filename <- paste(filepath, '/ powerCox.csv', sep = "")
  powerCox <- read.csv(filename, header = TRUE)
  powerCox$coef <- rep(coef, 20)
  
  if(coef == -1){
    powerCox_1 <- powerCox[powerCox$Shape_para == 1.0, ]
    powerCox_2 <- powerCox[powerCox$Shape_para == 1.5, ]
  }
  else{
    powerCox_1_new <- powerCox[powerCox$Shape_para == 1.0, ]
    powerCox_2_new <- powerCox[powerCox$Shape_para == 1.5, ]
    powerCox_1 <- rbind(powerCox_1, powerCox_1_new)
    powerCox_2 <- rbind(powerCox_2, powerCox_1_new)
  }
}

### Shape para = 1.0
colnames(powerCox_1)
powerCox_1$MAF <- factor(powerCox_1$MAF)

# powerCox_1_mod <- powerCox_1 %>%
#   mutate(Shape_para = recode(Shape_para, "1.0" = "Shape_par = 1", "1.5" = "Shape_par = 1.5"))

powerCox_1_mod <- powerCox_1 %>%
  mutate(Cens_Prob = recode(Cens_Prob, "0.25" = "CR = 0.25", "0.50" = "CR = 0.50"))
powerCox_1_mod <- powerCox_1_mod %>%
  mutate(MAF = recode(MAF, "0.02" = "MAF = 0.25", "0.05" = "MAF = 0.50", "0.1" = "MAF = 0.1","0.2" = "MAF = 0.2","0.5" = "MAF = 0.5"))

t1e_plot <- ggplot(data = powerCox_1_mod, aes(x = coef)) +
  geom_line(aes(y = Score, group = 1, color = "Score"), size = 0.7, linetype = "dashed") +
  geom_point(aes(y = Score, color = "Score"), size = 2) +
  geom_line(aes(y = Wald, group = 1, color = "Wald"), size = 0.7, linetype = "dashed") +
  geom_point(aes(y = Wald, color = "Wald"), size = 2) +
  geom_line(aes(y = LikeRatio, group = 1, color = "LikeRatio"), size = 0.7, linetype = "dashed") +
  geom_point(aes(y = LikeRatio, color = "LikeRatio"), size = 2) +
  geom_line(aes(y = Firth, group = 1, color = "Firth"), size = 0.7, linetype = "dashed") +
  geom_point(aes(y = Firth, color = "Firth"), size = 2) +
  geom_line(aes(y = SPACox, group = 1, color = "SPACox"), size = 0.7, linetype = "dashed") + 
  geom_point(aes(y = SPACox, color = "SPACox"), size = 2) + 
  #scale_y_log10(limits = c(1e-9, 1e-2)) +
  ylim(0,1e-2) +
  scale_color_manual(name = "Methods", values = c("Score" = "steelblue", "Wald" = "hotpink", "LikeRatio" = "seagreen", 
                                                  "Firth" = "orange", "SPACox" = "lightcoral"))

t1e_plot + facet_grid(Cens_Prob ~ MAF) + 
  ylab("Power") + ggtitle("Shape para = 1.0") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5))

### Shape para = 1.5
colnames(powerCox_2)
powerCox_2$MAF <- factor(powerCox_2$MAF)

# powerCox_1_mod <- powerCox_1 %>%
#   mutate(Shape_para = recode(Shape_para, "1.0" = "Shape_par = 1", "1.5" = "Shape_par = 1.5"))

powerCox_2_mod <- powerCox_2 %>%
  mutate(Cens_Prob = recode(Cens_Prob, "0.25" = "CR = 0.25", "0.50" = "CR = 0.50"))
powerCox_2_mod <- powerCox_2_mod %>%
  mutate(MAF = recode(MAF, "0.02" = "MAF = 0.25", "0.05" = "MAF = 0.50", "0.1" = "MAF = 0.1","0.2" = "MAF = 0.2","0.5" = "MAF = 0.5"))

t1e_plot <- ggplot(data = powerCox_2_mod, aes(x = coef)) +
  geom_line(aes(y = Score, group = 1, color = "Score"), size = 0.7, linetype = "dashed") +
  geom_point(aes(y = Score, color = "Score"), size = 2) +
  geom_line(aes(y = Wald, group = 1, color = "Wald"), size = 0.7, linetype = "dashed") +
  geom_point(aes(y = Wald, color = "Wald"), size = 2) +
  geom_line(aes(y = LikeRatio, group = 1, color = "LikeRatio"), size = 0.7, linetype = "dashed") +
  geom_point(aes(y = LikeRatio, color = "LikeRatio"), size = 2) +
  geom_line(aes(y = Firth, group = 1, color = "Firth"), size = 0.7, linetype = "dashed") +
  geom_point(aes(y = Firth, color = "Firth"), size = 2) +
  geom_line(aes(y = SPACox, group = 1, color = "SPACox"), size = 0.7, linetype = "dashed") + 
  geom_point(aes(y = SPACox, color = "SPACox"), size = 2) + 
  ylim(0,1) +
  scale_color_manual(name = "Methods", values = c("Score" = "steelblue", "Wald" = "hotpink", "LikeRatio" = "seagreen", 
                                                  "Firth" = "orange", "SPACox" = "lightcoral"))

t1e_plot + facet_grid(Cens_Prob ~ MAF) + 
  ylab("Power") + ggtitle("Shape para = 1.5") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5))


####


