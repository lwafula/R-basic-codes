# TASK: 
# 1.	fit an MNL (later: mixed logit) on the information in designs1s3 
# (all choice sets of size 3) (call this dataset design1) , RID is identification
# updates:
## 1. do not delete those that did not finish the experiment: we want to assess 
### preference differences between those that finished and those that did not finish the
#### experiment

## 2. fit a discrete time hazard. Compare estimates with a self-written code. The self-
### written code, if it works fine, will be included in a program for estimating a DCE +
### dropout code

rm(list = ls())

setwd('C:/KU_Leuven/DCE/People/VMartina/DeborahStreet')

# data 4 design 1 ####
d1 <- data.table::data.table(openxlsx::read.xlsx('data/m_ds_forshaad_design1s3.xlsx',
                                                 sheet = 1, startRow = 1, colNames = TRUE))
# unique entries
length(unique(d1[['RID']])) # 156

# number of choice sets answered by RID
d1[, `:=`(fnd = .N), by= RID]
table(d1[, fnd]) ; length(unique(d1[fnd==18, RID])) # 131 completed 18 choice sets, 

# sample completers ####
# for completers, delete 100 to remain with a comparable number as for non-completers
set.seed(07022020)
idfin <- sample(unique(d1[fnd==18, RID]), 31, replace = FALSE)
sort(idfin)

# keep either those that did not finish or those sampled above
d1 <- d1[(fnd<18|RID%in%idfin),]
length(unique(d1[['RID']])) # 56 i.e 25 + 31

# indicator of finished or not finished the study
d1[, `:=`(FinI=ifelse(fnd==18,1,0))]
apply(table(unique(d1[,.(RID, FinI)])),2,sum)

# EDA
table(d1[['SEQ']]-d1[['SCENARIO']])
table(d1[['DESIGN_ROW']])

# 5 attributes, 3 alternatives

d1<- d1[,!(c(4,21))]

d1<- d1[order(RID, SCENARIO)]

# change column names for easy reshaping to long format
dfnames<- colnames(d1)[-(c(1:4,20))]
nnames <- c()

for(i in 1:length(dfnames)){
        nnames[i] <- 
                paste0(strsplit(dfnames[i],"_")[[1]][c(2)],"_",strsplit(dfnames[i],"_")[[1]][c(1)])
}

colnames(d1) <- c(colnames(d1)[1:4],nnames, colnames(d1)[20])

# convert to long format
colA = paste("x1_a", 1:3, sep = ""); colB = paste("x2_a", 1:3, sep = "")
colC = paste("x3_a", 1:3, sep = ""); colD = paste("x4_a", 1:3, sep = "")
colE = paste("x5_a", 1:3, sep = "")
d1long = data.table::melt(d1, measure = list(colA, colB, colC, colD, colE), 
                          value.name = c("x1", "x2", "x3", "x4", "x5"))
d1long <- d1long[order(RID, SCENARIO)]
d1long[, `:=`(choice = 1.*(pref1==variable))]

# create dummies for the atributes
d1long[, `:=` (x1_d2 = 1.*(x1==2), x1_d3 = 1.*(x1==3), 
               x2_d2 = 1.*(x2==2), x2_d3 = 1.*(x2==3),
               x3_d2 = 1.*(x3==2), x3_d3 = 1.*(x3==3),
               x4_d2 = 1.*(x4==2), x4_d3 = 1.*(x4==3),
               x5_d2 = 1.*(x5==2), x5_d3 = 1.*(x5==3))]

# d1long <- d1long[,-(6:10)]

# choice data by RID
M <- max(as.numeric(unique(d1long[,variable])))
S <- max(as.numeric(unique(d1long[,SCENARIO])))

# create new IDs
apply(table(unique(d1long[, RID, by= .(FinI)])),1, sum)
ids = unique(d1long[, RID])
d1long <- d1long[, PID:= match(RID, table = ids)]

I <- max(as.numeric(unique(d1long[,PID])))
table(unique(d1long[, PID, by= .(FinI)])[, FinI])

## dropout indicator ####
d1long[, `:=`(stop_xp = ifelse(is.na(d1long[,choice]),1,0))]
d1long.event <- unique(d1long[,.(PID, FinI, SCENARIO, stop_xp)])[order(PID)]

m.event <- glm(stop_xp ~ SCENARIO-1, family = 'binomial', data = d1long.event)
summary(m.event)

# BUILD: LOGIT MODEL ####
# events
evntSI <- array(rep(1, I*S), c(S,I))
for (i in 1:I) {
        for(s in 1:S){
                evntSI[s,i] <- ifelse(s %in% d1long.event[PID==i, SCENARIO],
                                      d1long.event[PID==i & SCENARIO==s, stop_xp], 
                                      evntSI[s,i])
        }
        
}

# variable
int.SI <- array(rep(1, I*S), c(S,I))
setID.SI <- array(rep(0, I*S), c(S,I))
for (i in 1:I) {
        for(s in 1:S){
                setID.SI[s,i] <- ifelse(s %in% d1long.event[PID==i, SCENARIO],
                                        d1long.event[PID==i & SCENARIO==s, SCENARIO], 
                                        setID.SI[s,i])
        }
        
}

P <- 1
logist<-function(x){1/(1+exp(-x))} # sigmoid function
fbin <- function(param){
        
        # prior
        # prior <- sum(log(dnorm(param)))
        
        beta <- param[1:P]
        # beta.SI <- ((beta%o%rep(1,S))%o%rep(1,I))
        
        # 
        pin.SI <- logist(setID.SI*beta)
        
        #probabilities
        temp <- exp(apply(evntSI*log(pin.SI)+(1-evntSI)*log(1-pin.SI),2, sum))
        # f <- -2*(prior + sum(log(temp)))
        f <- -sum(log(temp))
}

# tbeta<- c(-2, 0.2)
tbeta<- c(-2)
m<-optim(tbeta,fbin, hessian=TRUE, method = "BFGS")
m$par
c('Est' = m$par, 'Std. Error' = sqrt(1/m$hessian), 
      'z-value' = m$par/sqrt(1/m$hessian), 'Pr(>|z|)' = 2*pnorm(m$par/sqrt(1/m$hessian)))
summary(m.event)$coefficients

# TWO VARIABLES
# ********************************
m.event <- glm(stop_xp ~ SCENARIO, family = 'binomial', data = d1long.event)
summary(m.event)

d1long.event[, `:=`(I=1)]
P <- 2
X.PSI <- array(rep(0, P*S*I), c(P,S,I))

for (i in 1:I) {
        for(s in 1:S){
                X.PSI[1,s,i] <- ifelse(s %in% d1long.event[PID==i, SCENARIO],
                                       as.numeric(d1long.event[PID==i & SCENARIO==s, 5]), X.PSI[1,s,i]);
                X.PSI[2,s,i] <- ifelse(s %in% d1long.event[PID==i, SCENARIO],
                                       as.numeric(d1long.event[PID==i & SCENARIO==s, 3]), 
                                       X.PSI[2,s,i])
        }
}

f2.1 <- function(param){
        
        # prior
        # prior <- sum(log(dnorm(param)))
        
        beta <- param[1:P]
        beta.PSI <- ((beta%o%rep(1,S))%o%rep(1,I))
        
        # 
        pin.SI <- logist(apply(X.PSI*beta.PSI, c(2,3), sum))
        
        #probabilities
        temp <- exp(apply(evntSI*log(pin.SI)+(1-evntSI)*log(1-pin.SI),2, sum))
        # f <- -2*(prior + sum(log(temp)))
        f <- -sum(log(temp))
}

tbeta<- c(-2,-0.2)

m2.1<-optim(tbeta,f2.1, hessian=TRUE, method = "BFGS")
m2.1$par
cbind('Est' = m2.1$par, 'Std. Error' = sqrt(diag(solve(m2.1$hessian))), 
  'z-value' = m2.1$par/ sqrt(diag(solve(m2.1$hessian))), 
  'Pr(>|z|)' = 2*pnorm(m2.1$par/ sqrt(diag(solve(m2.1$hessian)))))
summary(m.event)$coefficients
