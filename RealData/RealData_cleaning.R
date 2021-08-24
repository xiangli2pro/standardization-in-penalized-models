##---------------------------------------------------------------

## Data cleaning of NAMCS data

## Author: Xiang Li
## Date: 05-13-2021

## Input: namcs_2016.csv (2016 NAMCS data requested from CDC https://www.cdc.gov/nchs/ahcd/index.htm)

## Output: 
## X_Y.rda (cleaned predictors and outcome)
## NAMCS16_Clean.RData (cleaned predictors and outcome with weights)

##---------------------------------------------------------------




## import library
library(data.table)
library(tidyverse)
library(dmm)

## import data
coded16_full<-fread("namcs_2016.csv")

## select predictors
dt16sub<-coded16_full %>% 
  dplyr::select(CPSUM,CSTRATM,PATWT,PHYSWT,PATCODE,PHYCODE,opioid_use,
                #PATIENT GENERAL
                VMONTH,VDAYR,AGE,SEX,RACERETH,PAYTYPER,USETOBAC,EVERTOBAC,
                INJURY72,INJPOISAD,PRIMCARE,REFER,SENBEFOR,PASTVIS,MAJOR, 
                #PATIENT CONDITION
                ETOHAB,ALZHD,ARTHRTIS,ASTHMA,ASTH_SEV,ASTH_CON,ADD,AUTISM,
                CEBVD,CKD,COPD,CHF,CAD,DEPRN,diabetes,ESRD,HEPB,HEPC,HPE,HIV,
                HYPLIPID,HTN,OBESITY,OSA,OSTPRSIS,SUBSTAB,NOCHRON,TOTCHRON, 
                #SERVICES PROVIDED
                SERVICES,ETOH,BREAST,DEPRESS,DVS,FOOT,NEURO,PELVIC,RECTAL,
                RETINAL,SKIN,SUBST,BMP,CBC,CHLAMYD,CMP,CREAT,BLDCX,TRTCX,
                URNCX,OTHCX,GLUCOSE,GCT,HGBA,HEPTEST,HIVTEST,HPVDNA,CHOLEST,
                HEPATIC,PAP,PREGTEST,PSA,STREP,THYROID,URINE,VITD,ANYIMAGE,
                BONEDENS,CATSCAN,ECHOCARD,OTHULTRA,MAMMO,MRI,XRAY,OTHIMAGE,
                AUDIO,BIOPSY,BIOPROV,CARDIAC,COLON,CRYO,EKG,EEG,EMG,EXCISION,
                FETAL,PEAK,SIGMOID,SPIRO,TONO,TBTEST,EGD,SIGCOLON,CSW,CAM,DME,
                HOMEHLTH,MENTAL,OCCUPY,PT,PSYCHOTH,RADTHER,WOUND,ETOHED,ASTHMAED,
                ASTHMAP,DIAEDUC,DIETNUTR,EXERCISE,FAMPLAN,GRWTHDEV,INJPREV,STDPREV,
                STRESMGT,SUBSTED,TOBACED,WTREDUC,SERVCNT,ALLSERV,MED,#NUMMED#,
                numother,
                NUMNEW,NUMCONT,NOPROVID,PHYS,PHYSASST,NPNMW,RNLPN,MHP,OTHPROV,
                PROVNONE,TIMEMD,NODISP,RETREFPHY,REFOTHMD,RETAPPT1,RETAPPT2,RETAPPT3,
                RETUNSP,RETNEED,ERADMHOS,OTHDISP,
                #SERVICES PROVIDED
                SPECR,SPECCAT,MDDO,RETYPOFFR,SOLO,EMPSTAT,OWNSR,PATEVEN,
                NHVISR,HOMVISR,HOSVISR,TELCONR,ECONR,ECONTR,ACEPTNEW,PHYSCOMP,
                COMPPROD,COMPSAT,COMPQUAL,COMPDROF,COMPFIN,COMPUNK,COMPREF,
                SASDAPPT,SDAPPTPCT,APPTTIME,REGIONOFF,MSA)

## opioid_use as outcome
dt16sub$opioid_use<-factor(dt16sub$opioid_use)

##
dt16sub <-  dt16sub %>% mutate(
  
  AGE = as.numeric(case_when(
    AGE == "92 years or older" ~ "92",
    TRUE ~ AGE
  )),
  
  ## PASTVIS (Past Visits in last 12 months)
  PASTVIS = as.numeric(case_when(
    PASTVIS == "35 or more (except psychiatry)" ~ "35",
    PASTVIS == "Not applicable" ~ "999",
    TRUE ~ PASTVIS
  )),
  
  ## TOTCHRON (Total number of chronic conditions)
  TOTCHRON = as.numeric(case_when(
    TOTCHRON == "Entire item blank" ~ "999",
    TRUE ~ TOTCHRON
  )),
  
  ## SDAPPTPCT (What percent of daily visits are same day appointments)
  SDAPPTPCT = as.numeric(case_when(
    SDAPPTPCT == "Refused to answer" ~ "666",
    SDAPPTPCT == "Blank" ~ "999",
    TRUE ~ SDAPPTPCT
  )),
  
  diabetes = as.factor(case_when(
    diabetes == 1 ~ "Yes",
    diabetes == 0 ~ "No"
  ))
) %>% mutate(
  
  PASTVIS = factor(case_when(
    PASTVIS == 999 ~ "999",
    (PASTVIS >= 0) & (PASTVIS < 2) ~ "0-1",
    (PASTVIS >= 2) & (PASTVIS < 6) ~ "2-5",
    (PASTVIS >= 6) & (PASTVIS < 10) ~ "6-9",
    (PASTVIS >= 10) & (PASTVIS < 15) ~ "10-14",
    (PASTVIS >= 15) & (PASTVIS <= 35) ~ "15-35",
    (PASTVIS > 35) ~ ">=36"
  ),levels=c("0-1","2-5","6-9","10-14","15-35",">=36","999")),
  
  TOTCHRON = factor(case_when(
    TOTCHRON == 999 ~ "999",
    TOTCHRON == 0 ~ "0",
    TOTCHRON == 1 ~ "1",
    TOTCHRON == 2 ~ "2",
    TOTCHRON == 3 ~ "3",
    (TOTCHRON >= 4) & (TOTCHRON <= 6) ~ "4-6",
    TOTCHRON >= 7 ~ ">=7"
  ),levels = c("0","1","2","3","4-6",">=7","999")),
  
  ## SERVCNT (Total number of services ordered or provided including vital signs)
  SERVCNT = factor(case_when(
    (SERVCNT >=0) & (SERVCNT <=1) ~ "0-1",
    (SERVCNT >=2) & (SERVCNT <=3) ~ "2-3",
    (SERVCNT >=4) & (SERVCNT <=5) ~ "4-5",
    (SERVCNT >=6) & (SERVCNT <=10) ~ "6-10",
    (SERVCNT >=11) ~ ">=11"
  ), levels = c("0-1", "2-3", "4-5","6-10",">=11")),
  
  SDAPPTPCT = factor(case_when(
    SDAPPTPCT == 999 ~ "999",
    SDAPPTPCT == 666 ~ "Refused to answer",
    SDAPPTPCT < 25 ~ "<25",
    SDAPPTPCT >= 25 ~ ">=25"
  ),levels = c("<25", ">=25", "999", "Refused to answer"))
)

## TIMEMD (time spent with physicians in mins)
## change 1 hour 20 mins to 80 mins
timemd <- as.character(dt16sub$TIMEMD)
timemd <- strsplit(timemd, ' ')

getTime <- function(x) {
  result <- 0
  if ("hour" %in% x) {
    result <- as.numeric(x[which(x == "hour") - 1])*60
  }
  if("minutes" %in% x) {
    result<-result+as.numeric(x[which(x=="minutes")-1])
  }
  result
}

dt16sub$TIMEMD <- sapply(timemd,getTime) 


#####################################
##Set reference level        ########
#####################################

dt16sub <- dt16sub %>% mutate(
  # from start model
  RACERETH = relevel(as.factor(RACERETH),'Non-Hispanic White'),
  USETOBAC=relevel(as.factor(USETOBAC),'Not current'),
  EVERTOBAC=relevel(as.factor(EVERTOBAC),'Never'),
  INJPOISAD=relevel(as.factor(INJPOISAD),"No, visit is not related to injury/trauma, overdose/poisoning, or adverse effect of medical/surgical treatment"),
  MAJOR=relevel(as.factor(MAJOR),'Preventive care'),
  ARTHRTIS=relevel(as.factor(ARTHRTIS),'No'),
  SPECR=relevel(as.factor(SPECR),'General/family practice'),
  diabetes=relevel(as.factor(diabetes),'No'),
  CAD=relevel(as.factor(CAD),'No'),
  SUBSTAB=relevel(as.factor(SUBSTAB),'No'),
  TOTCHRON=relevel(as.factor(TOTCHRON),'0'),
  HTN=relevel(as.factor(HTN),'No'),
  DEPRN=relevel(as.factor(DEPRN),'No'),
  VMONTH=relevel(as.factor(VMONTH),'November'),
  PHYSCOMP=relevel(as.factor(PHYSCOMP),'Mix of salary and share of billings or other measures of performance'),
  HYPLIPID=relevel(as.factor(HYPLIPID),'No'),
  
  # patient general
  PAYTYPER=relevel(as.factor(PAYTYPER),"Private insurance"),
  INJURY72=relevel(as.factor(INJURY72),"No"),
  PRIMCARE=relevel(as.factor(PRIMCARE),"No"),
  REFER=relevel(as.factor(REFER),"No"),
  
  # patients conditions
  ASTH_SEV=relevel(as.factor(ASTH_SEV),"None recorded"),
  ASTH_CON=relevel(as.factor(ASTH_CON),"None recorded"),
  NOCHRON=relevel(as.factor(NOCHRON),"No"),
  
  # services provided
  SERVICES=relevel(as.factor(SERVICES),"No services were reported" ),
  BIOPROV=relevel(as.factor(BIOPROV),"No"),
  ALLSERV=relevel(as.factor(ALLSERV),"No services were reported" ),
  MED=relevel(as.factor(MED),"No medications were reported"),
  SERVCNT=relevel(as.factor(SERVCNT),"4-5"),
  
  # physican relatedx
  SPECCAT=relevel(as.factor(SPECCAT),"Primary care specialty"),
  MDDO=relevel(as.factor(MDDO),"D.O. - Doctor of Osteopathy"),
  EMPSTAT=relevel(as.factor(EMPSTAT),"Full-owner"),
  OWNSR=relevel(as.factor(OWNSR),"Physician or physician group"),
  RETYPOFFR=relevel(as.factor(RETYPOFFR),"Private solo or group practice"),
  HOMVISR=relevel(as.factor(HOMVISR),"No"),
  TELCONR=relevel(as.factor(TELCONR),"No"),
  ACEPTNEW=relevel(as.factor(ACEPTNEW),"No"),
  COMPPROD=relevel(as.factor(COMPPROD),"Box is not marked"),
  COMPSAT=relevel(as.factor(COMPSAT),"Box is not marked"),
  COMPQUAL=relevel(as.factor(COMPQUAL),"Box is not marked"),
  COMPDROF=relevel(as.factor(COMPDROF),"Box is not marked"),
  COMPFIN=relevel(as.factor(COMPFIN),"Box is not marked"),
  COMPUNK=relevel(as.factor(COMPUNK),"Box is not marked"),
  COMPREF=relevel(as.factor(COMPREF),"Box is not marked"),
  SASDAPPT=relevel(as.factor(SASDAPPT),"No"),
  APPTTIME=relevel(as.factor(APPTTIME),"Within 1 week"),
  REGIONOFF=relevel(as.factor(REGIONOFF),"South"),
  MSA=relevel(as.factor(MSA),"MSA (Metropolitan Statistical Area)"),
  
  # newly added 0130
  SIGCOLON=relevel(as.factor(SIGCOLON),"No"),
  SEX=relevel(as.factor(SEX),"Female"),
  SENBEFOR=relevel(as.factor(SENBEFOR),"No, new patient"),
  NOPROVID=relevel(as.factor(NOPROVID),"No categories marked"),
  NODISP=relevel(as.factor(NODISP),"No categories marked"),
  
  VDAYR=relevel(as.factor(VDAYR),"Friday"),
  SOLO=relevel(as.factor(SOLO),"Non-solo"),
  PATEVEN=relevel(as.factor(PATEVEN),"No"),
  NHVISR=relevel(as.factor(NHVISR),"No"),
  HOSVISR=relevel(as.factor(HOSVISR),"No"),
  ECONR=relevel(as.factor(ECONR),"No"),
  ECONTR=relevel(as.factor(ECONTR),"No")
)


str(dt16sub, list.len=ncol(dt16sub))


##################
## set Y and X
##################

## check the class of all the variables
## some are factors, some are numeric and characters
# dt_str <- summary.default(dt16sub) %>% as.data.frame() %>% 
#   dplyr::group_by(Var1) %>%  tidyr::spread(key = Var2, value = Freq)


# set response Y
Y <- dt16sub$opioid_use 
Y <- ifelse(Y==0,-1,1)
# table(Y)

# set predictors X
# removed survey structure related variable and varaible that are highly correlated
X <- dt16sub %>%
  dplyr::select(-c(CPSUM,CSTRATM,PATWT,PHYSWT,PATCODE,PHYCODE,opioid_use,NUMCONT,NUMNEW,MED,SPECCAT,COLON))

# recode sdapptpct; split 999 to 0 and missing(2483), remove sdadappt.
X <- X %>% mutate(
  SDAPPTPCT = factor(case_when(
    (SDAPPTPCT == "999") & (SASDAPPT == "No") ~ "0",
    (SDAPPTPCT == "999") & (SASDAPPT == "999") ~ "999",
    (SDAPPTPCT == "<25") ~ "<25",
    (SDAPPTPCT == ">=25") ~ ">=25",
    (SDAPPTPCT == "Refused to answer") ~ "Refused to answer"
  ),levels = c("<25","0","999",">=25","Refused to answer"))
) %>% dplyr::select(-SASDAPPT)


## check the class and reference level of the variables
X_str <- summary.default(X) %>% as.data.frame() %>% 
  dplyr::group_by(Var1) %>%  tidyr::spread(key = Var2, value = Freq)

X_str$Unique <- apply(X,2,function(x) length(unique(x)))
X_str$Yes_No <- apply(X,2,function(x) sum(c("No","Yes") %in% unique(x)))

## convert 2 level character to 2 level factor with "No" as reference
cha_index <- which(X_str$Mode=="character" & X_str$Unique==2 & X_str$Yes_No==2)
X[,cha_index] <- lapply(X[,cha_index], function(x) factor(x,levels = c("No","Yes")))

## recheck the structure of X
X_str_new <- summary.default(X) %>% as.data.frame() %>% 
  dplyr::group_by(Var1) %>%  tidyr::spread(key = Var2, value = Freq) 

X_str_new$Num_levels <- apply(X,2,function(x) length(unique(x)))
X_str_new$Yes_No <- apply(X,2,function(x) sum(c("No","Yes") %in% unique(x)))
X_str_new$Ref_level <- lapply(X,function(x) ifelse(is.factor(x),levels(x)[1],"numeric"))


# save(X_str_new,file = "X_str_new.rda")

save(Y,X,file = "X_Y.rda")

save(dt16sub, file = "NAMCS16_Clean.RData")

