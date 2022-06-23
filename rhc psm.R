rhc <- read.csv("rhc.csv")
# Change the Age variable into categories  below 50, [50,60), [60,70), [70,80), above 80
# categorizing a continuous variable is not recommended.
rhc$age <- cut(rhc$age,breaks=c(-Inf, 50, 60, 70, 80, Inf),right=FALSE)
# Re-order the levels of race to white, black and other
rhc$race <- factor(rhc$race, levels=c("white","black","other"))
# merging disease categories
rhc$cat1 <- as.character(rhc$cat1)
rhc$cat1[rhc$cat1 == "Lung Cancer"] <- "Other"
rhc$cat1[rhc$cat1 == "COPD"] <- "Other"
rhc$cat1[rhc$cat1 == "Coma"] <- "Other"
rhc$cat1[rhc$cat1 == "Cirrhosis"] <- "Other"
rhc$cat1[rhc$cat1 == "Colon Cancer"] <- "Other"
rhc$cat1[rhc$cat1 == "MOSF w/Malignancy"] <- "MOSF"
rhc$cat1[rhc$cat1 == "MOSF w/Sepsis"] <- "MOSF"
rhc$cat1 <- as.factor(rhc$cat1)
# Change the baseline for gender to Male
rhc$sex <- as.factor(rhc$sex)
rhc$sex <- relevel(rhc$sex, ref = "Male")
# Regroup the levels for disease categories to "ARF","CHF","MOSF","Other".
levels(rhc$ca) <- c("Metastatic","None","Localized (Yes)")
# Rename the levels of "ca" (Cancer) to "Metastatic","None" and "Localized (Yes)" 
rhc$ca <- factor(rhc$ca, levels=c("None","Localized (Yes)","Metastatic"))
# re-order the levels to "None","Localized (Yes)" and "Metastatic"
rhc$ca <- factor(rhc$ca, levels=c("None","Localized (Yes)","Metastatic"))
# create a new variable called "numcom" to count number of comorbidities illness for each person  (12 categories)
rhc$numcom <- rhc$cardiohx + rhc$chfhx + rhc$dementhx + rhc$psychhx + 
  rhc$chrpulhx + rhc$renalhx + rhc$liverhx + rhc$gibledhx + rhc$malighx + 
  rhc$immunhx + rhc$transhx +rhc$amihx
rhc2 <- rhc[c("age","sex", "race","cat1", "ca", "dnr1", "aps1",
              "surv2md1","numcom","adld3p","das2d3pc","temp1",
              "hrt1","meanbp1","resp1","wblc1","pafi1","paco21",
              "ph1","crea1","alb1","scoma1","swang1", "death")]
names(rhc2) <- c("age","sex", "race","Disease.category", "Cancer", 
                 "DNR.status", "APACHE.III.score", "Pr.2mo.survival",
                 "No.of.comorbidity","ADLs.2wk.prior","DASI.2wk.prior",
                 "Temperature","Heart.rate","Blood.pressure",
                 "Respiratory.rate","WBC.count","PaO2.by.FIO2","PaCO2",
                 "pH","Creatinine","Albumin","GComa.Score","RHC", "Death")
dim(rhc2)
rhc2$age <- factor(rhc2$age, levels = c("[-Inf,50)","[50,60)","[60,70)",
                                        "[70,80)","[80, Inf)"), 
                   ordered = TRUE)
levels(rhc2$age)
# Assess missing values
require(DataExplorer)
plot_missing(rhc2) 
# simplifying
rhc2$ADLs.2wk.prior <- NULL
rhc2$Cancer <- NULL
analytic.data0 <- rhc2 
rm(rhc2)
dim(analytic.data0)
table(analytic.data0$RHC)
table(analytic.data0$Death)
# inducing some bias in the study!!
analytic.data0$ID <- 1:nrow(analytic.data0)
# Younger age and no treated and did not survive
id1 <- analytic.data0$ID[analytic.data0$RHC!="RHC" & analytic.data0$age =="[-Inf,50)" & analytic.data0$Death=="Yes"]
# Female and not treated and did not survive
id2 <- analytic.data0$ID[analytic.data0$RHC!="RHC" & analytic.data0$sex !="Male" & analytic.data0$Death=="Yes"]
# Other race (other than white and black) and not treated and did not survive
id3 <- analytic.data0$ID[analytic.data0$RHC!="RHC" & analytic.data0$race =="other" & analytic.data0$Death=="Yes"]
# Abnormal heart rate and not treated and did not survive
id4 <- analytic.data0$ID[analytic.data0$RHC!="RHC" & analytic.data0$Heart.rate < 70 & analytic.data0$Heart.rate > 110 & analytic.data0$Death=="Yes"]
idx <- unique(c(id1,id2,id3,id4))
length(idx)
set.seed(123)
# take a random sample of the above group
exclude.id <- sample(idx, 
                     size = round(length(idx)*3/4), 
                     replace = FALSE)
head(sort(exclude.id))
# exclude the selected sample from the analytic data
analytic.data <- analytic.data0[ !analytic.data0$ID %in% exclude.id, ]
head(sort(analytic.data$ID))
table(analytic.data$RHC)
table(analytic.data$Death)
dim(analytic.data)
##########################################################################

## Analysis strategy: matching RHC patients with non-RHC patients
# I Outcome Death (Y )
# I Death at any time up to 180 Days
# I Treatment swang1 (A: Swan-Ganz catheter)
# I Whether or not a patient received a RHC
# I Covariate list: L (age, sex, race , . . .)
names(analytic.data)
library(tableone)
# 2 x 2 table
tab0 <- CreateTableOne(vars = "RHC",
                       data = analytic.data,
                       strata = "Death")
print(tab0, showAllLevels = TRUE)
table(analytic.data$RHC)
table(analytic.data$Death)

baselinevars <- c("age","sex", "race")
# Table 1
tab1 <- CreateTableOne(vars = baselinevars,
                       data = analytic.data,
                       strata = "Death", includeNA = TRUE,
                       test = TRUE, smd = FALSE)
print(tab1, showAllLevels = FALSE, smd = FALSE)

## crude regression

# adjust the exposure variable (primary interest)
fit0 <- glm(I(Death=="Yes")~RHC,
            family=binomial, data = analytic.data)
require(Publish)
publish(fit0)

## adjusted regression
# adjust the exposure variable + demographics
fit1 <- glm(I(Death=="Yes")~RHC + age + sex + race,
            family=binomial, data = analytic.data)
publish(fit1)

# adjust the exposure variable + adjustment variables
baselinevars <- c("age","sex", "race","Disease.category",
                  "DNR.status", "APACHE.III.score",
                  "Pr.2mo.survival","No.of.comorbidity",
                  "DASI.2wk.prior","Temperature",
                  "Heart.rate", "Blood.pressure",
                  "Respiratory.rate", "WBC.count",
                  "PaO2.by.FIO2","PaCO2","pH",
                  "Creatinine","Albumin","GComa.Score")
out.formula <- as.formula(paste("I(Death=='Yes')", "~ RHC +",
                                paste(baselinevars,
                                      collapse = "+")))
out.formula
fit2 <- glm(out.formula,
            family=binomial, data = analytic.data)
publish(fit2)

## Assumptions check 

## plot(fit2, which =1) for linear regression

plot(fit2, which =3)
plot(fit2, which =4)
# 35:00
# 
# How sure are you about the model-specification?
#   I Interaction?
#   I Polynomial?
#   I Potential solution?
#   I Exact Matching



### start to show that exact matching will hurt you

## Exact Matching: 2 variables
var.comb <- do.call('paste0',
                    analytic.data[, c('race', 'sex')])
length(table(var.comb))
table(var.comb)
table(analytic.data$RHC,var.comb)

require(MatchIt)
# exact match by sex and race (?non parametric solution)
m.out = matchit (RHC=="RHC" ~ sex + race,
                 data = analytic.data,
                 method = "exact")
m.out
summary(m.out)


var.comb <- do.call('paste0',
                    analytic.data[, c('race', 'sex', 'age')])
length(table(var.comb))
table(analytic.data$RHC,var.comb=="otherMale[80, Inf)")
table(analytic.data$RHC,var.comb=="otherMale[80, Inf)")
table(analytic.data$RHC,var.comb=="otherFemale[80, Inf)")

# exact match by age, sex and race
m.out = matchit (RHC=="RHC" ~ age + sex + race,
                 data = analytic.data,
                 method = "exact")
m.out
print(m.out)
summary(m.out)
matched.data <- match.data(m.out)
dim(matched.data)
nrow(analytic.data)-nrow(matched.data) # subjects deleted
# Not taking into account of matched sets
fit1m <- glm(I(Death=="Yes")~RHC,
             family=binomial, data = matched.data)
publish(fit1m)
m.out = matchit (RHC=="RHC" ~ age + sex + race +
                   Disease.category + DNR.status,
                 data = analytic.data,
                 method = "exact")
m.out
matched.data <- match.data(m.out)
dim(matched.data)

fit2m <- glm(I(Death=="Yes")~RHC,
             family=binomial, data = matched.data)
publish(fit2m)

m.out = matchit (RHC=="RHC" ~ age + sex + race +
                   Disease.category + DNR.status+
                   Heart.rate, # continuous
                 data = analytic.data,
                 method = "exact")
m.out

m.out = matchit (RHC=="RHC" ~ age + sex + race +
                   Disease.category + DNR.status+
                   Heart.rate + Blood.pressure +
                   Temperature,
                 data = analytic.data,
                 method = "exact")
m.out
matched.data <- match.data(m.out)
dim(matched.data)

nrow(analytic.data)-nrow(matched.data) # subjects deleted

fit3m <- glm(I(Death=="Yes")~RHC,
             family=binomial, data = matched.data)
publish(fit3m)

## if you match your dataset you 39:00

matched.data[1:2,]
colnames(matched.data)

mytest = matched.data[1:2,c(1:5,10:12,21,22)]

##################
# SOLUTION IS PSM
##################

baselinevars
