# Analysis of the data from the longitudinal study in Mali

# This file contains the code for the following data analysis and plots:
# 1. Data overview: number of individuals by status of infection
# 2. Age and status of infection
# 3. Time from 1st follow-up to first positive PCR
# 4. Time from 1st follow-up to clinical malaria
# 5. Time from enrollment to clinical malaria
# 6. Time from enrollment to first positive PCR


# load packages and data:
library(survival)
library(survminer)
library(dplyr)
library(readxl)
library(graphics)

Workbook <- read_excel("Mali-2012-data.xlsx")
Workbook_full <- read_excel("Mali-2012-data (full).xlsx",range = cell_cols(c(NA,128)))


############################################################################################################
###################### 1. Data overview: number of individuals by status of infection ######################
############################################################################################################


# date of enrollment:
sort(unique(Workbook$`VisitDate of enrolm`))

data <- Workbook
data <- unique(data) # remove duplicates

# Number of individuals by infection status at the end of the dry season (for Fig. 6)

# number of individuals with RDT+PCR+
data_pp <- unique(data[data$Interaction=="pospos" | data$RdtPcr=="pospos",])
length(unique(data_pp$subj_id))-as.numeric(any(is.na(data_pp$subj_id)))

# number of individuals with RDT-PCR+:
data_pn <- unique(data[data$Interaction=="posneg",])
length(unique(data_pn$subj_id))-as.numeric(any(is.na(data_pn$subj_id)))

# number of individuals with RDT-PCR-:
data_nn <- unique(data[data$Interaction=="negneg" | data$RdtPcr=="negneg",])
length(unique(data_nn$subj_id))


############################################################################################################
###################################### 2. Age and status of infection ######################################
############################################################################################################

age_par <- data.frame(id = Workbook$`subj_id of Time2PCR+`,result=Workbook$Interaction,age=Workbook$age)
age_par$result <- ifelse(age_par$result == "pospos", "RDT+PCR+", ifelse(age_par$result== "posneg","RDT-PCR+","RDT-PCR-"))
age_par <- dplyr::mutate(age_par,PCR=ifelse(result=="RDT-PCR-","neg","pos"))
age_par <- dplyr::mutate(age_par,age_group=ifelse(age_par$age <=1, 1,ifelse(age_par$age <=2,2,ifelse(age_par$age <=3,3,
           ifelse(age_par$age <=4,4,ifelse(age_par$age <=5,5,ifelse(age_par$age <=6,6,ifelse(age_par$age <=7,7,
           ifelse(age_par$age <=8,8,ifelse(age_par$age <=9,9,ifelse(age_par$age <=10,10,12)))))))))))

# Remove duplicates:
age_par <- unique(age_par)

# dplyr::count(age_par,result) # How many individuals are there for each status of infection?
# dplyr::count(age_par,age_group) # Age distribution

# Fraction of individuals with parasites that have a parasite density below the level of detection for RDT:
sum(age_par$result=="RDT-PCR+",na.rm = TRUE)/(sum(age_par$result=="RDT+PCR+",na.rm = TRUE)+sum(age_par$result=="RDT-PCR+",na.rm = TRUE))

# Fraction of RDT and PCR positive by age
age_res <- table(age=age_par$age_group,age_par$result)
age_res <- rbind(RDTpos=age_res[,3],PCRpos=age_res[,2]+age_res[,3],all=age_res[,1]+age_res[,2]+age_res[,3])
age_res_frac <- rbind(RDTpos=age_res[1,]/age_res[3,],PCRpos=age_res[2,]/age_res[3,])

# Barplot for RDT/PCR results by age and age distribution boxplots for carriers and non-carriers (Fig. S12)
par(mfrow=c(1,2))
barplot(age_res_frac,main=expression(paste("Fraction of ",RDT^+phantom(1),"and ",PCR^+phantom(1),"by age")),
        xlab="Age [years]",ylab="Fraction",ylim=c(0,0.7),
        names.arg = c(expression(phantom(.)<=1),"2","3","4","5","6","7","8","9","10",">10"),
        legend.text = c(expression(RDT^+phantom(1)),expression(PCR^+phantom(1))), beside = TRUE)
boxplot(age_par$age[which(age_par$result%in%c("RDT+PCR+","RDT-PCR-"))]~age_par$result[which(age_par$result%in%c("RDT+PCR+","RDT-PCR-"))],
        main=expression(paste("Age distribution of ",PCR^+phantom(.),"and ",PCR^-phantom(.))),xlab="Status of infection (May 2012)",
        ylab="Age [years]",ylim=c(0,12.2),names=c(expression(PCR^-phantom(.)),expression(RDT^+phantom(.))))
par(mfrow=c(1,1))

# Comparison of the median age with the Wilcoxon rank-sum test:
age_par <- age_par[which(age_par$result%in%c("RDT+PCR+","RDT-PCR-")),]
wilcox.test(age_par$age~age_par$PCR, alternative="less", data=age_par)


##############################################################################################################################
###################################### 3. Time from 1st follow-up to first positive PCR ######################################
##############################################################################################################################


# PCR results:
PCRres <- dplyr::select(Workbook,'PCR+ t1','PCR+ t2','PCR* t3','PCR+ t4','PCR+ t5','PCR+ t6','PCR+ t7',
                        'PCR+ t8','PCR+ t9','PCR+ t10','PCR+ t11','PCR+ t12')
# Visit dates:
VisitDate <- dplyr::select(Workbook,'VisitDate of ENROLMENT','VisitDate of untitled 2','T02VisitDate',
                           'T01VisitDate','T04VisitDate','T05VisitDate','T06VisitDate','T07VisitDate',
                           'T08VisitDate','T09VisitDate','T10VisitDate','T11VisitDate','T12VisitDate')
# Time from enrollment to the different visits:
time2visit <- data.frame(id = Workbook$id,
                         dplyr::transmute(VisitDate,V1 = as.numeric(VisitDate$`VisitDate of untitled 2`-VisitDate$`VisitDate of ENROLMENT`)),
                         dplyr::transmute(VisitDate,V2 = as.numeric(VisitDate$T02VisitDate-VisitDate$`VisitDate of ENROLMENT`)),
                         dplyr::transmute(VisitDate,V3 = as.numeric(VisitDate$T01VisitDate-VisitDate$`VisitDate of ENROLMENT`)),
                         dplyr::transmute(VisitDate,V4 = as.numeric(VisitDate$T04VisitDate-VisitDate$`VisitDate of ENROLMENT`)),
                         dplyr::transmute(VisitDate,V5 = as.numeric(VisitDate$T05VisitDate-VisitDate$`VisitDate of ENROLMENT`)),
                         dplyr::transmute(VisitDate,V6 = as.numeric(VisitDate$T06VisitDate-VisitDate$`VisitDate of ENROLMENT`)),
                         dplyr::transmute(VisitDate,V7 = as.numeric(VisitDate$T07VisitDate-VisitDate$`VisitDate of ENROLMENT`)),
                         dplyr::transmute(VisitDate,V8 = as.numeric(VisitDate$T08VisitDate-VisitDate$`VisitDate of ENROLMENT`)),
                         dplyr::transmute(VisitDate,V9 = as.numeric(VisitDate$T09VisitDate-VisitDate$`VisitDate of ENROLMENT`)),
                         dplyr::transmute(VisitDate,V10 = as.numeric(VisitDate$T10VisitDate-VisitDate$`VisitDate of ENROLMENT`)),
                         dplyr::transmute(VisitDate,V11 = as.numeric(VisitDate$T11VisitDate-VisitDate$`VisitDate of ENROLMENT`)),
                         dplyr::transmute(VisitDate,V12 = as.numeric(VisitDate$T12VisitDate-VisitDate$`VisitDate of ENROLMENT`)))

# Individuals with inconsistencies in the follow-up dates:
ind <- rep(0,length(time2visit$V1))
for(i in 1:length(time2visit$V1)){
  tmp <- c(FALSE,!is.na(time2visit[i,2:length(time2visit[i,])]))
  if(any(diff(as.numeric(time2visit[i,tmp]))<=0)){
    ind[i] <- i
  }
}
ind1 <- ind[ind>0] 

# Maximal number of days allowed between two PCR results: (otherwise consider the data as censored)
tmax <- 20 

VisitDate2 <- VisitDate[,2:13]
PCRres2 <- PCRres[,2:12]

ind <- rep(0,length(PCRres2$`PCR+ t12`)) # index for considered individuals
event <- rep(0,length(PCRres2$`PCR+ t12`)) # 1 if event, 0 if censored
times <- rep(0,length(PCRres2$`PCR+ t12`)) # time to first PCR+ or censored

for (i in 1:length(Workbook$id)){
  tmp <- time2visit[i,3:length(time2visit[i,])]-time2visit[i,2:(length(time2visit[i,])-1)] # time between consecutive follow-up visits
  if(!is.na(time2visit[i,2])){ # If the first visit date is missing, we do not know the time from first follow-up to first PCR+, so we exclude these cases.
    if(!(i%in%ind1)){ 
      k <- length(PCRres2[i,])
    }else{
      k <- which(tmp<=0)[1]-1 # if there are inconsistencies in the visit dates, we only consider the data before the first inconsistency
    }
    if(any(PCRres2[i,1:k]=='+',na.rm = TRUE)){ # if there is a positive result
      j <- which(PCRres2[i,1:k]=='+')[1] # first positive PCR result
      if(!any(is.na(PCRres2[i,1:(j-1)])) && max(tmp[1:j],na.rm = TRUE)<tmax){ # no missing PCR result before first positive and max. 20 days between visits
        ind[i] <- ind[i]+1
        times[i] <- as.numeric(VisitDate2[i,j+1]-VisitDate2[i,1])
        event[i] <- 1
      }else if (any(is.na(PCRres2[i,1:(j-1)]))){ # if there are NAs before the first positive PCR result
        tmp2 <- tmp
        h <- 0 
        for (l in 1:(sum(is.na(PCRres2[i,1:(j-1)])))){
          tmp2[which(is.na(PCRres2[i,1:k]))[l]+1] <- tmp2[which(is.na(PCRres2[i,1:k]))[l]+1]+tmp2[which(is.na(PCRres2[i,1:k]))[l]] # time between visits with known results
          tmp2[which(is.na(PCRres2[i,1:k]))[l]] <- 0
          if(is.na(tmp2[which(is.na(PCRres2[i,1:k]))[l]+1])){
            h <- 1 # there is a visit date missing
          }
        }
        if(h==0 && max(tmp2[1:j],na.rm = TRUE)<tmax){ # no missing visit date & time between visits with known results is maximally 20 days
          ind[i] <- ind[i]+1
          times[i] <- as.numeric(VisitDate2[i,j+1]-VisitDate2[i,1])
          event[i] <- 1
        }
      }
    }
    if(!any(is.na(PCRres2[i,1:k])) && all(PCRres2[i,1:k]=='-') && max(tmp[1:(j-1)],na.rm = TRUE)<tmax){ # if there are only negative PCR results
      ind[i] <- ind[i]+1
      times[i] <- as.numeric(VisitDate2[i,length(VisitDate2[i,])]-VisitDate2[i,1])
    }
    if(!is.na(PCRres2[i,1]) && PCRres2[i,1]=='-'){ # if there are only negative results before missing results
      if(any(is.na(PCRres2[i,1:k]))){ 
        j <- which(is.na(PCRres2[i,1:k]))[1]-1
        if(all(PCRres2[i,1:j]=='-')){
          ind[i] <- ind[i]+1
          times[i] <- as.numeric(VisitDate2[i,j+1]-VisitDate2[i,1])
          if(is.na(times[i]) && j>1){
            j <- j-1
            times[i] <- as.numeric(VisitDate2[i,j]-VisitDate2[i,1])
          }
        }
      }
    }
  }
}

new_data <- data.frame(id=Workbook$id[as.logical(ind)],time=times[as.logical(ind)],
                       event=event[as.logical(ind)],PCR=Workbook$`PCR+may`[as.logical(ind)],
                       age=Workbook$age[as.logical(ind)],gender=Workbook$gender[as.logical(ind)])

new_data <- new_data[-which(Workbook$`PCR+ t1`[which(Workbook$id%in%new_data$id)]=='+'),] # exclude individuals with positive PCR result at the first visit
new_data <- dplyr::mutate(new_data,age_group=ifelse(new_data$age <=3, "(0,3]",ifelse(new_data$age <=6,"(3,6]",ifelse(new_data$age <=9,"(6,9]",">9"))))

# Remove duplicates:
new_data <- unique(new_data)

# Kaplan-Meier curve for RDT-PCR- and RDT+PCR+ (treated) (Fig. S16a)
fit1 <- survfit(Surv(time,event) ~ PCR, data = new_data)
summary(fit1)
ggsurvplot(fit1, data = new_data, title=expression(paste("Time from first follow-up to ",PCR^+phantom(.))), pval = TRUE, 
           legend.title = "Parasites at enrollment in May 2012:", legend.labs = c("PCR-","RDT+ (treated)"), legend="bottom",
           ylab=expression(paste("Probability of being ",PCR^-phantom(.))),xlab="Time [days]")

# Cox PH model for RDT-PCR- and RDT+PCR+ (treated) including age (Table S4)
fit.coxph1 <- coxph(Surv(time,event) ~ PCR+age, data = new_data)
summary(fit.coxph1)

# Kaplan-Meier curve for different ages: (Fig. S17)
fit2 <- survfit(Surv(time,event) ~ age_group, data = new_data)
summary(fit2)
ggsurvplot(fit2, data = new_data, title=expression(paste("Time from first follow-up to ",PCR^+phantom(.))), pval = TRUE, 
           legend.title = "Age [years]:",legend.labs = c("(0,3]","(3,6]","(6,9]",">9"),legend="bottom",
           ylab=expression(paste("Probability of being ",PCR^-phantom(.))),xlab="Time [days]")

# Consider only individuals with age > 3
# Kaplan-Meier curve for different ages: (Fig. S18)
new_data2 <- dplyr::filter(new_data,age>3)
fit3 <- survfit(Surv(time,event) ~ age_group, data = new_data2)
summary(fit3)
ggsurvplot(fit3, data = new_data2, title=expression(paste("Time from first follow-up to ",PCR^+phantom(.))), pval = TRUE, 
           legend.title = "Age [years]:",legend.labs = c("(3,6]","(6,9]",">9"),legend="bottom",
           ylab=expression(paste("Probability of being ",PCR^-phantom(.))),xlab="Time [days]")

# Cox PH model: (Table S5)
fit3.coxph <- coxph(Surv(time,event) ~ PCR+age, data = new_data2)
summary(fit3.coxph)


#######################################################################################################################
################################### 4. Time from 1st follow-up to clinical malaria ####################################
#######################################################################################################################

time2mal <- data.frame(id = Workbook$`subj_id of Time2PCR+`,result=Workbook$Interaction,enr=Workbook$`VisitDate of ENROLMENT`,
                       mal=Workbook$`VisitDate of 1stMAL`,lastvisit=Workbook$T12VisitDate,age=Workbook$age)

time2mal$result <- ifelse(time2mal$result == "pospos", "RDT+PCR+", ifelse(time2mal$result== "posneg","RDT-PCR+","RDT-PCR-"))
time2mal <- dplyr::mutate(time2mal,PCR=ifelse(result=="RDT-PCR-","neg","pos"))
time2mal <- dplyr::mutate(time2mal,age_group = ifelse(time2mal$age <=3, "(0,3]",ifelse(time2mal$age <=6,"(3,6]",ifelse(time2mal$age <=9,"(6,9]",">9"))))
time2mal <- dplyr::mutate(time2mal,age_sq = time2mal$age^2)
time2mal <- dplyr::mutate(time2mal,firstvisit=Workbook$`VisitDate of untitled 2`)

# add missing lastvisit and firstvisit information, if possible:
time2mal$lastvisit[is.na(time2mal$lastvisit)] <- apply(Workbook[is.na(time2mal$lastvisit),c(1+2*(3:14))], 1, 
                                                       function(x) ifelse( !all(is.na(x)), max(x, na.rm=T), NA))  
time2mal$firstvisit[is.na(time2mal$firstvisit)] <- apply(Workbook[is.na(time2mal$firstvisit),c(1+2*(3:14))], 1, 
                                                         function(x) ifelse( !all(is.na(x)), min(x, na.rm=T), NA))  

# remove rows that contain only NA:
time2mal <- time2mal[rowSums(is.na(time2mal[ ,c(1:5,10)]))<6, ]

# Remove duplicates:
time2mal <- unique(time2mal) # remove duplicates
if(length(time2mal$id) != length(unique(time2mal$id))){ # remove rows with the same id but other different values ### TODO: all these rows differ only in the date for clinical malaria => use the earlier date
  time2mal <- time2mal[!(time2mal$id%in%time2mal$id[duplicated(time2mal$id)]),]
}

# Visit Date 1st Mal and Visit Date 1st Mal 2 comparison:
any(unique(Workbook$`subj_id of Time2PCR+`[which(Workbook$`VisitDate of 1stMAL`-Workbook$`VisitDate of 1stMAL 2`!=0)])%in%time2mal$id)

# Event or censored:
time2mal <- dplyr::mutate(time2mal,event=ifelse(is.na(mal),0,1))

# time from first follow-up to malaria in days:
time <- floor(as.numeric(difftime(time2mal$mal, time2mal$firstvisit, units = "days")))
time2mal <- dplyr::mutate(time2mal,time=time)
time2mal[time2mal$event==0,] <- dplyr::mutate(time2mal[time2mal$event==0,],time=lastvisit-firstvisit)
time2mal$time <- floor(as.numeric(time2mal$time))

# clinical malaria before the first follow-up -> exclude these individuals:
time2mal <- dplyr::filter(time2mal,time>0)

# Kaplan-Meier curve for time to malaria, only RDT-PCR- and RDT+PCR+ (treated) (Fig. S16b)
tmp <- time2mal[time2mal$result%in%c("RDT-PCR-","RDT+PCR+"),]
fit3 <- survfit(Surv(time,event) ~ result, data = tmp)
summary(fit3)
ggsurvplot(fit3, data = tmp, title="Time to clinical malaria", pval = TRUE, 
           legend.title = "Parasites at enrollment in May 2012:", legend.labs = c("PCR-","RDT+ (treated)"), 
           legend="bottom", ylab="Probability of being malaria free",xlab="Time [days]")


########################################################################################################################
##################################### 5. Time from enrollment to clinical malaria ###################################### 
########################################################################################################################


time2mal <- data.frame(id = Workbook$`subj_id of Time2PCR+`,result=Workbook$Interaction,enr=Workbook$`VisitDate of ENROLMENT`,
                       mal=Workbook$`VisitDate of 1stMAL`,lastvisit=Workbook$T12VisitDate,age=Workbook$age)

time2mal$result <- ifelse(time2mal$result == "pospos", "RDT+PCR+", ifelse(time2mal$result== "posneg","RDT-PCR+","RDT-PCR-"))
time2mal <- dplyr::mutate(time2mal,PCR=ifelse(result=="RDT-PCR-","neg","pos"))
time2mal <- dplyr::mutate(time2mal,age_group = ifelse(time2mal$age <=3, "(0,3]",ifelse(time2mal$age <=6,"(3,6]",ifelse(time2mal$age <=9,"(6,9]",">9"))))

# Add missing lastvisit information if possible:
time2mal$lastvisit[is.na(time2mal$lastvisit)] <- apply(Workbook[is.na(time2mal$lastvisit),c(2*(4:15))], 1, 
                                                       function(x) ifelse( !all(is.na(x)), max(x, na.rm=T), NA))  

# Remove duplicates and rows with the same id but otherwise different values:
time2mal <- unique(time2mal) # remove duplicates
if(length(time2mal$id) != length(unique(time2mal$id))){
  time2mal <- time2mal[!(time2mal$id%in%time2mal$id[duplicated(time2mal$id)]),]
}

# Event or censored:
time2mal <- dplyr::mutate(time2mal,event=ifelse(is.na(mal),0,1))

# time from enrollment to clinical malaria in days:
time2mal <- dplyr::mutate(time2mal,time=mal-enr)
# If no clinical malaria episode was observed, then we consider the time from enrollment to the last visit as malaria-free and have a right-censored event:
time2mal[time2mal$event==0,] <- dplyr::mutate(time2mal[time2mal$event==0,],time=lastvisit-enr) 
time2mal$time <- floor(time2mal$time) # round to full days

# Kaplan-Meier curve for RDT-PCR- and RDT+PCR+ (treated) only (Fig. 7B)
tmp <- time2mal[time2mal$result%in%c("RDT-PCR-","RDT+PCR+"),]
fit2 <- survfit(Surv(time,event) ~ result, data = tmp)
summary(fit2)
ggsurvplot(fit2, data = tmp, title="Time to clinical malaria", pval = TRUE, 
           legend.title = "Parasites at enrollment in May 2012:", legend.labs = c("PCR-","RDT+ (treated)"), 
           legend="bottom", ylab="Probability of being malaria free",xlab="Time [days]")


##########################################################################################################################
##################################### 6. Time from enrollment to first positive PCR ######################################
##########################################################################################################################


# PCR results:
PCRres <- dplyr::select(Workbook,'PCR+ t1','PCR+ t2','PCR* t3','PCR+ t4','PCR+ t5','PCR+ t6','PCR+ t7',
                        'PCR+ t8','PCR+ t9','PCR+ t10','PCR+ t11','PCR+ t12')
# Visit dates:
VisitDate <- dplyr::select(Workbook,'VisitDate of ENROLMENT','VisitDate of untitled 2','T02VisitDate',
                           'T01VisitDate','T04VisitDate','T05VisitDate','T06VisitDate','T07VisitDate',
                           'T08VisitDate','T09VisitDate','T10VisitDate','T11VisitDate','T12VisitDate')
# Time from enrollment to the different visits:
time2visit <- data.frame(id = Workbook$id,
                         dplyr::transmute(VisitDate,V1 = as.numeric(VisitDate$`VisitDate of untitled 2`-VisitDate$`VisitDate of ENROLMENT`)),
                         dplyr::transmute(VisitDate,V2 = as.numeric(VisitDate$T02VisitDate-VisitDate$`VisitDate of ENROLMENT`)),
                         dplyr::transmute(VisitDate,V3 = as.numeric(VisitDate$T01VisitDate-VisitDate$`VisitDate of ENROLMENT`)),
                         dplyr::transmute(VisitDate,V4 = as.numeric(VisitDate$T04VisitDate-VisitDate$`VisitDate of ENROLMENT`)),
                         dplyr::transmute(VisitDate,V5 = as.numeric(VisitDate$T05VisitDate-VisitDate$`VisitDate of ENROLMENT`)),
                         dplyr::transmute(VisitDate,V6 = as.numeric(VisitDate$T06VisitDate-VisitDate$`VisitDate of ENROLMENT`)),
                         dplyr::transmute(VisitDate,V7 = as.numeric(VisitDate$T07VisitDate-VisitDate$`VisitDate of ENROLMENT`)),
                         dplyr::transmute(VisitDate,V8 = as.numeric(VisitDate$T08VisitDate-VisitDate$`VisitDate of ENROLMENT`)),
                         dplyr::transmute(VisitDate,V9 = as.numeric(VisitDate$T09VisitDate-VisitDate$`VisitDate of ENROLMENT`)),
                         dplyr::transmute(VisitDate,V10 = as.numeric(VisitDate$T10VisitDate-VisitDate$`VisitDate of ENROLMENT`)),
                         dplyr::transmute(VisitDate,V11 = as.numeric(VisitDate$T11VisitDate-VisitDate$`VisitDate of ENROLMENT`)),
                         dplyr::transmute(VisitDate,V12 = as.numeric(VisitDate$T12VisitDate-VisitDate$`VisitDate of ENROLMENT`)))

ind <- rep(0,length(PCRres$`PCR+ t1`)) # index for considered individuals
event <- rep(0,length(PCRres$`PCR+ t1`)) # 1 if event, 0 if censored
times <- rep(0,length(PCRres$`PCR+ t1`)) # time to first PCR+ or censored

for (i in 1:length(Workbook$id)){
  if(any(PCRres[i,]=='+',na.rm = TRUE)){ # if there is a positive PCR result
    j <- which(PCRres[i,]=='+')[1] # first positive PCR result
    if(!any(is.na(PCRres[i,1:j]))){ # if there are no missing PCR results
      ind[i] <- ind[i]+1
      times[i] <- as.numeric(VisitDate[i,j+1]-VisitDate[i,1])
      event[i] <- 1
    }
  }
  if(!any(is.na(PCRres[i,])) && all(PCRres[i,]=='-')){ # if all PCR results are negative
    ind[i] <- ind[i]+1
    times[i] <- as.numeric(VisitDate[i,13]-VisitDate[i,1])
  }
  if(!is.na(PCRres[i,1]) && PCRres[i,1]=='-'){ # at least one negative PCR result: censored at the last negative PCR result
    if(any(is.na(PCRres[i,]))){
      j <- which(is.na(PCRres[i,]))[1]-1
      if(all(PCRres[i,1:j]=='-')){
        ind[i] <- ind[i]+1
        times[i] <- as.numeric(VisitDate[i,j+1]-VisitDate[i,1])
        if(is.na(times[i]) && j>1){ # if there is a missing date, then take the last known date
          j <- j-1
          times[i] <- as.numeric(VisitDate[i,j]-VisitDate[i,1])
        }
      }
    }
  }
}

new_data <- data.frame(id=Workbook$id[as.logical(ind)],time=times[as.logical(ind)],
                       event=event[as.logical(ind)],PCR=Workbook$`PCR+may`[as.logical(ind)],age=Workbook$age[as.logical(ind)])
new_data <- dplyr::mutate(new_data,age_group=ifelse(new_data$age <=3, "(0,3]",ifelse(new_data$age <=6,"(3,6]",ifelse(new_data$age <=9,"(6,9]",">9"))))

# Remove duplicates and rows with the same id but otherwise different values:
new_data <- unique(new_data)

# Kaplan-Meier curve for RDT-PCR- and RDT+PCR+ (treated) (Fig. 7A)
fit1 <- survfit(Surv(time,event) ~ PCR, data = new_data)
summary(fit1)
ggsurvplot(fit1, data = new_data, title=expression(paste("Time from enrollment to ",PCR^+phantom(.))), pval = TRUE, 
           legend.title = "Parasites at enrollment in May 2012:", legend.labs = c("PCR-","RDT+ (treated)"), legend="bottom",
           ylab=expression(paste("Probability of being ",PCR^-phantom(.))),xlab="Time [days]")

# Cox PH model for RDT-PCR- and RDT+PCR+ (treated) including age (Table S3)
fit.coxph <- coxph(Surv(time,event) ~ PCR+age, data = new_data)
summary(fit.coxph)

# Kaplan-Meier curve for different age groups: (Fig. S15)
fit2 <- survfit(Surv(time,event) ~ age_group, data = new_data)
summary(fit2)
ggsurvplot(fit2, data = new_data, title=expression(paste("Time from enrollment to ",PCR^+phantom(.))), pval = TRUE, 
           legend.title = "Age [years]:",legend.labs = c("(0,3]","(3,6]","(6,9]",">9"),legend="bottom",
           ylab=expression(paste("Probability of being ",PCR^-phantom(.))),xlab="Time [days]")


### Fit an exponential distribution with left-censoring:
t_cen <- 70 # time limit for left-censoring (in days)
new_data4 <- new_data
new_data4$event[new_data$time<t_cen] <- -1 # event=-1 indicates left-censoring
fun_neg <- Vectorize(function(x){-log(x)*sum(new_data4$event[!is.na(new_data4$PCR) & new_data4$PCR=="neg"]==1)+
    x*sum(new_data4$time[!is.na(new_data4$PCR) & new_data4$PCR=="neg" & new_data4$event%in%c(0,1)])-
    sum(log(1-exp(-x*new_data4$time[!is.na(new_data4$PCR) & new_data4$PCR=="neg" & new_data4$event==-1])))})
fun_pos <- Vectorize(function(x){-log(x)*sum(new_data4$event[!is.na(new_data4$PCR) & new_data4$PCR=="pos"]==1)+
    x*sum(new_data4$time[!is.na(new_data4$PCR) & new_data4$PCR=="pos" & new_data4$event%in%c(0,1)])-
    sum(log(1-exp(-x*new_data4$time[!is.na(new_data4$PCR) & new_data4$PCR=="pos" & new_data4$event==-1])))})
seq(1.57e-2,1.59e-2,by=1e-7)[which.min(fun_pos(seq(1.57e-2,1.59e-2,by=1e-7)))] # MLE for 70 days censoring limit and RDT+PCR+ (treated) patients: 0.0158633
seq(3.5e-3,4.5e-3,by=1e-7)[which.min(fun_neg(seq(3.5e-3,4.5e-3,by=1e-7)))] # MLE for 70 days censoring limit and RDT-PCR- patients: 0.0042027

# visualize fit: (Fig. S13)
fit1 <- survfit(Surv(time,event) ~ PCR, data = new_data)
lam_pos <- 0.0158633
lam_neg <- 0.0042027
plot(c(0,fit1$time[1:fit1$strata[1]]),c(1,fit1$surv[1:fit1$strata[1]]),type="s",lwd=2,col="red",ylim=c(0,1),las=1,xlab="Time [days]", 
     ylab = "Fraction without infection", main="Time from enrollment to first infection")
lines(c(0,fit1$time[(fit1$strata[1]+1):length(fit1$surv)]),c(1,fit1$surv[(fit1$strata[1]+1):length(fit1$surv)]),lwd=2,col="blue",type="s")
lines(seq(0,300,by=1),exp(-lam_neg*seq(0,300,by=1)),lty=2,col="red",lwd=2)
lines(seq(0,300,by=1),exp(-lam_pos*seq(0,300,by=1)),lty=2,col="blue",lwd=2)
legend(130,1,legend=c("Carriers","Non-carriers"),col=c("blue","red"),lwd=2,lty=1)

# biting rate per year:
lam_pos*365 # carriers
lam_neg*365 # non-carriers
lam_pos/lam_neg # biting rate of carriers compared to non-carriers

frac_car <- 100/579 # fraction of carriers (also fraction of bites for carriers in the case of homogeneous distribution of infectious bites)
frac_non_car <- 424/579
# fraction of bites for carriers in the heterogeneous case:
lam_pos*frac_car/(lam_pos*frac_car+lam_neg*(1-frac_car)) # assuming that RDT-PCR+ have the same infection risk as non-carriers
# commparison of bites for carriers in the heterogeneous and the homogeneous case:
(lam_pos*frac_car/(lam_pos*frac_car+lam_neg*(1-frac_car)))/frac_car
# time to x% of the mosquito population being infected in the heterogeneous case compared to the homogeneous case:
frac_car/(lam_pos*frac_car/(lam_pos*frac_car+lam_neg*(1-frac_car)))

