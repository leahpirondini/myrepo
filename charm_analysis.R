### CHARM

# install packages

# load packages
library(readstata13)
library(gtsummary)
library(tidyverse)
library(survival)

# set working directory
setwd("C:/Users/leahp/OneDrive/Documents/JOB_LSHTM/Covariate_Adjustment/Data/CHARM")

# load data
covars <- read.dta13("mscbase.dta")
hospx <- read.dta13("hosp_x.dta")
endpoints0 <- read.dta13("mscendpoint0.dta")
endpoints1 <- read.dta13("mscendpoint1.dta")
endpoints2 <- read.dta13("mscendpoint2.dta")
endpoints3 <- read.dta13("mscendpoint3.dta")
endpoints4 <- read.dta13("mscendpoint4.dta")
endpoints5 <- read.dta13("mscendpoint5.dta")
endpoints6 <- read.dta13("mscendpoint6.dta")
endpoints7 <- read.dta13("mscendpoint7.dta")
endpoints8 <- read.dta13("mscendpoint8.dta")
endpoints9 <- read.dta13("mscendpoint9.dta")
endpoints10 <- read.dta13("mscendpoint10.dta")
endpoints11 <- read.dta13("mscendpoint11.dta")
endpoints12 <- read.dta13("mscendpoint12.dta")
endpoints13 <- read.dta13("mscendpoint13.dta")

# data checks
any(duplicated(covars$patid))
any(duplicated(endpoints0$patid))
table(endpoints0$study)


any(duplicated(hospx$patid))
length(unique(hospx[["patid"]]))
unique(hospx$study)
hospx.unique <- unique(hospx[c("patid","study","center","country")])

# add treatment group and study id to covars data
treatdta <- endpoints0 %>%
  select(treat, patid, study)
covars <- left_join(covars, treatdta, by = "patid")

# label variables
covars$sex <- factor(covars2$sex,
                     levels = c(1,2),
                     labels = c("Male", "Female"))
covars$smoke <- factor(covars2$smoke,
                       levels = c(1,2,3),
                       labels = c("Non-smoker", "Previous smoker", "Current smoker"))
covars$nyha <- factor(covars2$nyha,
                      levels = c(2,3,4),
                      labels = c("II", "III", "IV"))
covars$agecat <- factor(covars2$agecat,
                        levels = c(1,2,3,4,5,6,7,8),
                        labels = c("<50", "50-54","55-59","60-64","65-69","70-74",
                                   "75-79",">=80"))
covars$region6 <- factor(covars2$region6,
                         levels = c(1,2,3,4,5,6),
                         labels = c("US", "Canada", "Scandinavia", "Eastern Europe",
                                    "Rest of europe", "Rest of world"))
covars$race <- factor(covars2$race,
                      levels = c(1,2,3,4,5,6,7),
                      labels = c("European origin", "Black", "South Asian", "Arab/Middle East", 
                                 "Oriental", "Malay", "Other"))

# summarise baseline data in table 1 format
names(covars)

baseline.vars <- subset(covars, study=="SH-AHS-0007",
                        select=c("treat", "sex", "ef", "smoke", "nyha", "hrate", "sbp", "dbp", "weightkg",
                                 "heightcm", "chftime", "chftimey", "bmi", "age", "agecat", "region6",
                                 "race", "smokstop", "smokno", "smokyear", "oedema"))

baseline.vars %>% tbl_summary(by = treat,
                              missing = "always",
                              missing_text="(Missing)",
                              statistic = list(all_continuous() ~ "{mean} ({sd})"),
                              digits = all_continuous() ~ 1) %>%
  add_overall() %>%
  bold_labels()

baseline.meds <- subset(covars, study=="SH-AHS-0007",
                        select=c("treat", "medh1","medh2","medh3","medh4","medh5","medh6","medh7","medh8",
                                 "medh9","medh10","medh11","medh12","medh13","medh14","medh15","medh16", "medh17",
                                 "hfdurq1", "hfdurq2", "hfdurq3", "hfdurq4", "hfdurq5", "hfdurq6",
                                 "hfdurq7", "hfdurq8", "hfdurq9", "hfdurq10","hfdurq11",
                                 "hfcriq0", "hfcriq1", "hfcriq2", "hfcriq3", "hfcriq4",
                                 "hfcriq5", "hfcriq6", "hfcriq7", "hfcriq8", "hfcriq9", "hfcriq10",
                                 "ecgq1", "ecgq2", "ecgq3", "ecgq4", "ecgq5", "ecgq6", "ecgq9",
                                 "cmed1", "cmed2", "cmed3", "cmed4", "cmed5", "cmed6", "cmed7", 
                                 "cmed8", "cmed9", "cmed10", "cmed11", "cmed12", "cmed13", "cmed14"))

baseline.meds %>% tbl_summary(by = treat,
                              missing = "always",
                              missing_text="(Missing)",
                              statistic = list(all_continuous() ~ "{mean} ({sd})"),
                              digits = all_continuous() ~ 1) %>%
  add_overall() %>%
  bold_labels()

# subset of preserve trial data
preserve <- subset(endpoints0, study=="SH-AHS-0007")
preserve <- left_join(preserve, covars, by = "patid")

# subset of added trial data
added <- subset(endpoints0, study=="SH-AHS-0006")
added <- left_join(added, covars, by = "patid")

# subset of alternative trial data
alternative <- subset(endpoints0, study=="SH-AHS-0003")
alternative <- left_join(alternative, covars, by = "patid")

# unadjusted trial stats
table(preserve$event)
with(preserve, table(event, tr))
table(added$event)
with(added, table(event, tr))
table(alternative$event)
with(alternative, table(event, tr))

unadj.cox.preserve = coxph(Surv(time,event)~tr,data=preserve)
summary(unadj.cox.preserve)

unadj.cox.added = coxph(Surv(time,event)~tr,data=added)
summary(unadj.cox.added)

unadj.cox.alternative = coxph(Surv(time,event)~tr,data=alternative)
summary(unadj.cox.alternative)

# adjusted analyses
univar.cox = coxph(Surv(time,event)~age,data=preserve)
univar.cox = coxph(Surv(time,event)~sex,data=preserve)
univar.cox = coxph(Surv(time,event)~smoke,data=preserve)
univar.cox = coxph(Surv(time,event)~nyha,data=preserve)
univar.cox = coxph(Surv(time,event)~sbp,data=preserve)
univar.cox = coxph(Surv(time,event)~bmi,data=preserve)
univar.cox = coxph(Surv(time,event)~region6,data=preserve)
univar.cox = coxph(Surv(time,event)~race,data=preserve)

summary(univar.cox)

adj.cox.preserve = coxph(Surv(time,event)~tr+age+region6,data=preserve)
summary(adj.cox.preserve)

table(preserve$nyha,preserve$agecat)

adj.cox.preserve = coxph(Surv(time,event)~tr+female+race+age+nyha
                         +ef+hrate+sbp+dbp+bmi+smoke
                         +medh1+medh2+medh3+medh4+medh5
                         +medh6+medh7+medh8+medh9+medh10+medh11+medh12
                         +hfdurq1+hfdurq2+hfdurq3
                         +cmed13+cmed2+cmed3+cmed14+cmed1+cmed5+cmed6+cmed7+cmed10
                         +hfcriq8
                         ,data=preserve)
summary(adj.cox.preserve)




