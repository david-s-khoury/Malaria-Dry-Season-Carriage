# Analysis of the simulated data of individuals 

# This file contains the code for the following data analysis and plots:
# 1. Comparison of carriers and non-carriers previous exposure
# 2. Sankey diagram of the change of infection status after 1 year of surveillance
# 3. Survival curves for comparison of homogeneous and heterogeneous infection risk
# 4. Time to first infection for homogeneous exposure with 10 different FOIs


# load packages and data:
library(survival)
library(survminer)
library(R.matlab)
library(ggplot2)
library(riverplot)


##########################################################################################################
###################### 1. Comparison of carriers and non-carriers previous exposure ######################
##########################################################################################################

# load data: 
data <- readMat("~/Simulated data/Data/Data-2-car-vs-non-car.mat")
data <- as.data.frame(data[[1]])
names(data) <- c("FOI","age","carrier","num_bites","last_bite_time")
data$carrier <- as.factor(data$carrier)
levels(data$carrier) <- c("non-carrier","carrier")
data <- dplyr::mutate(data,event=!is.na(last_bite_time))

### Number of bites in the previous season for different FOIs: (Fig. S5)
ggplot(data,aes(x=carrier,y=num_bites)) +
  geom_boxplot() +
  facet_wrap(~FOI) + 
  labs(x="", y="Number of bites", title="Number of infectious bites in the previous transmission season") +
  ylim(0,max(data$num_bites,na.rm = TRUE)*1.07) +
  theme_bw()

# boxplot for number of bites in the previous season for a FOI of 0.024 bites per day: (Fig. 4a)
ggplot(data[data$FOI==0.024,],aes(x=carrier,y=num_bites)) +
  geom_boxplot() +
  labs(x="", y="Number of bites", title="Number of infectious bites in the previous transmission season") +
  ylim(0,max(data$num_bites,na.rm = TRUE)*1.07) +
  theme_bw()

# p-values of comparison with a Wilcoxon rank-sum test:
w_p_fois <- matrix(0,10,1)
for(i in 1:10){
  if(any(data$carrier=="carrier" & data$FOI==unique(data$FOI)[i])){
    w_p_fois[i] <- wilcox.test(data$num_bites[data$carrier=="non-carrier" & data$FOI==unique(data$FOI)[i]],
                               data$num_bites[data$carrier=="carrier" & data$FOI==unique(data$FOI)[i]],alternative = "less")$p.value
  }else{
    w_p_fois[i] <- NaN
  }
}


### Time since the last bite survival curves for different FOIs: (Fig. 4b and Fig. S6)
for(i in 1:10){
  fit <- survfit(Surv(last_bite_time,event) ~ carrier, data = data[data$FOI==unique(data$FOI)[i],])
  # summary(fit)
  if(i==1){
    ggsurvplot(fit, data = data[data$FOI==unique(data$FOI)[i],], title=paste("Time since last infectious bite (FOI = ",unique(data$FOI)[i],")",sep=""), pval = TRUE, 
               legend.title = "", legend.labs = c("Non-carriers"), legend="bottom",break.time.by=100,
               # conf.int = TRUE,risk.table = TRUE,tables.height = 0.2,tables.theme = theme_cleantable(),break.time.by=50,xlim=c(0,230),
               ylab="Fraction without a bite",xlab="Time [days]",xlim=c(0,365))
  }else{
    ggsurvplot(fit, data = data[data$FOI==unique(data$FOI)[i],], title=paste("Time since last infectious bite (FOI = ",unique(data$FOI)[i],")",sep=""), pval = TRUE, 
               legend.title = "", legend.labs = c("Non-carriers","Carriers"), legend="bottom",break.time.by=100,
               # conf.int = TRUE,risk.table = TRUE,tables.height = 0.2,tables.theme = theme_cleantable(),break.time.by=50,xlim=c(0,230),
               ylab="Fraction without a bite",xlab="Time [days]",xlim=c(0,365))
  }
  ggsave(paste("~/OneDrive - UNSW/R/DrySeason/Homogeneous simulation carriers vs non-carriers/Survival-curve-data-2-",i,".pdf",sep=""), width = 6, height = 4)
}


##############################################################################################################################
###################### 2. Sankey diagram of the change of infection status after 1 year of surveillance ######################
##############################################################################################################################

# Define sample data frame (using the numbers from the MATLAB simulation):
nodes <- c("A","B","C","D","E","G")
edges <- list(A=list(D=15,E=47,G=12),B=list(D=26,E=87,G=71),C=list(D=6,E=84,G=652))
# Plot Sankey diagram: (Fig. S10)
r <- makeRiver( nodes, edges, node_xpos= c(1,1,1,2,2,2),
                node_labels= c( A="RDT+ (treated)", B="RDT-PCR+", C="PCR-", D="RDT+ (treated)",E="RDT-PCR+", G="PCR-"),
                node_styles= list( A= list(col="red"),B=list(col="yellow"),C=list(col= "blue"),
                                   D= list(col="red"),E=list(col="yellow"),G=list(col="blue")));
plot(r)


###############################################################################################################################
###################### 3. Survival curves for comparison of homogeneous and heterogeneous infection risk ######################
###############################################################################################################################

# Survival curves for the time from the first follow-up to PCR-detectable infection for carriers and non-carriers


### Homogeneous infection risk (Fig. 5g)
data <- readMat("~/Simulated data/Data/Survival curve data/survival-curve-data-3-6.mat")

time <- t(data[[1]])
event <- t(data[[4]]) # 0=right censored , 1=event , 2=left censored
rdt_pos <- data[[6]]
pcr_pos <- data[[2]]
group <- ifelse(pcr_pos==0,0,ifelse(rdt_pos==0,1,2))

surv_data <- as.data.frame(cbind(time,event,group))
names(surv_data) <- c("time","event","group")
surv_data2 <- dplyr::filter(surv_data,event%in%c(0,1)) # only events and right-censored data
surv_data3 <- dplyr::filter(surv_data2,group%in%c(0,2)) # only carriers (RDT+) and non-carriers (PCR-)

# Plot survival curves for RDT pos. and PCR neg only
fit <- survfit(Surv(time,event) ~ group, data = surv_data3)
summary(fit)
ggsurvplot(fit, data = surv_data3, title="Time from first follow-up to PCR positivity", pval = TRUE, 
           legend.title = "", legend.labs = c("Non-carriers","Carriers (treated)"), legend="bottom",
           # conf.int = TRUE,risk.table = TRUE,tables.height = 0.2,tables.theme = theme_cleantable(),break.time.by=50,xlim=c(0,230),
           ylab="Probability of being malaria free",xlab="Time [days]")
# ggsave("~/Survival-curve-data-3-6.pdf", width = 6, height = 4)


### Heterogeneous infection risk (Fig. 5h and Fig. S8a)
# load data:
data <- readMat("~/Simulated data/Data/Survival curve data/survival-curve-data-1.mat")

time <- t(data[[1]])
event <- t(data[[4]]) # 0=right censored , 1=event , 2=left censored
rdt_pos <- data[[6]]
pcr_pos <- data[[2]]
group <- ifelse(pcr_pos==0,0,ifelse(rdt_pos==0,1,2))

surv_data <- as.data.frame(cbind(time,event,group))
names(surv_data) <- c("time","event","group")
surv_data2 <- dplyr::filter(surv_data,event%in%c(0,1)) # only events and right-censored data
surv_data3 <- dplyr::filter(surv_data2,group%in%c(0,2)) # only carriers (RDT+) and non-carriers (PCR-)

# Plot survival curves for RDT pos. and PCR neg only
fit <- survfit(Surv(time,event) ~ group, data = surv_data3)
summary(fit)
ggsurvplot(fit, data = surv_data3, title="Time from first follow-up to PCR positivity", pval = TRUE, 
           legend.title = "", legend.labs = c("Non-carriers","Carriers (treated)"), legend="bottom",
           # conf.int = TRUE,risk.table = TRUE,tables.height = 0.2,tables.theme = theme_cleantable(),break.time.by=50,xlim=c(0,230),
           ylab="Probability of being malaria free",xlab="Time [days]")
# ggsave("~/Survival-curve-data-1.pdf", width = 6, height = 4)


########################################################################################################################
###################### 4. Time to first infection for homogeneous exposure with 10 different FOIs ######################
########################################################################################################################

# Time from first follow-up to positive PCR result for a population of 1,000 individuals with the same FOI (one of 10 
# different FOIs) but different age (Fig. S7 and p-values in Table S2).

all_fois <- c(0.004,0.008,0.012,0.016,0.02,0.024,0.028,0.032,0.036,0.04)

for(i in 1:10){
  filename <- paste("~/Simulated data/Data/Survival curve data/survival-curve-data-3-",i,".mat",sep="")
  data <- readMat(filename)
  
  foi_tmp <- all_fois[i]
  time <- t(data[[1]])
  event <- t(data[[4]]) # 0=right censored , 1=event , 2=left censored
  rdt_pos <- data[[6]]
  pcr_pos <- data[[2]]
  group <- ifelse(pcr_pos==0,0,ifelse(rdt_pos==0,1,2))
  
  surv_data <- as.data.frame(cbind(time,event,group))
  names(surv_data) <- c("time","event","group")
  surv_data2 <- dplyr::filter(surv_data,event%in%c(0,1)) # only events and right-censored data
  surv_data3 <- dplyr::filter(surv_data2,group%in%c(0,2)) # only carriers (RDT+) and non-carriers (PCR-)
  
  # Plot survival curves for RDT pos. and PCR neg only
  if(i==2){
    fit <- survfit(Surv(time,event) ~ group, data = surv_data3)
    summary(fit)
    ggsurvplot(fit, data = surv_data3, title=paste("FOI ",foi_tmp," bites per day",sep=""), pval = TRUE, 
               legend.title = "", legend.labs = c("Non-carriers"), legend="bottom",
               conf.int = FALSE,ylab="Fraction without infection",xlab="Time [days]")
    ggsave(paste("~/Survival-curve-data-3-",i,".pdf",sep=""), width = 6, height = 4)
  }else{
    fit <- survfit(Surv(time,event) ~ group, data = surv_data3)
    summary(fit)
    ggsurvplot(fit, data = surv_data3, title=paste("FOI ",foi_tmp," bites per day",sep=""), pval = TRUE, 
               legend.title = "", legend.labs = c("Non-carriers","Carriers (treated)"), legend="bottom",
               # conf.int = TRUE,risk.table = TRUE,tables.height = 0.2,tables.theme = theme_cleantable(),break.time.by=50,xlim=c(0,230),
               ylab="Fraction without infection",xlab="Time [days]")
    ggsave(paste("~/Survival-curve-data-3-",i,".pdf",sep=""), width = 6, height = 4)
  }
}






