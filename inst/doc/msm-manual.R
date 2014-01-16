### R code from vignette source 'msm-manual.Rnw'

###################################################
### code chunk number 1: msm-manual.Rnw:23-25
###################################################
version <- gsub("Version: +", "",
               packageDescription("msm", lib.loc="../../..")$Version)


###################################################
### code chunk number 2: msm-manual.Rnw:30-31
###################################################
cat(version)


###################################################
### code chunk number 3: msm-manual.Rnw:34-35
###################################################
cat(format(Sys.time(), "%d %B, %Y"))


###################################################
### code chunk number 4: msm-manual.Rnw:843-844
###################################################
options(width = 60)


###################################################
### code chunk number 5: msm-manual.Rnw:879-880
###################################################
library(msm)


###################################################
### code chunk number 6: msm-manual.Rnw:939-940
###################################################
cav[1:21,]


###################################################
### code chunk number 7: msm-manual.Rnw:950-951
###################################################
statetable.msm(state, PTNUM, data=cav)


###################################################
### code chunk number 8: msm-manual.Rnw:1002-1006
###################################################
twoway4.q  <-  rbind ( c(0, 0.25, 0, 0.25),
                       c(0.166, 0, 0.166, 0.166),
                       c(0, 0.25, 0, 0.25),
                       c(0, 0, 0, 0) )


###################################################
### code chunk number 9: msm-manual.Rnw:1040-1042
###################################################
twoway4.crude.q  <- crudeinits.msm(state ~ years, PTNUM, data=cav, 
                                   qmatrix=twoway4.q)


###################################################
### code chunk number 10: msm-manual.Rnw:1066-1068
###################################################
cav.msm <- msm( state ~ years, subject=PTNUM, data = cav,
                    qmatrix = twoway4.q, death = 4 )


###################################################
### code chunk number 11: msm-manual.Rnw:1098-1099 (eval = FALSE)
###################################################
## help(optim)


###################################################
### code chunk number 12: msm-manual.Rnw:1121-1122
###################################################
cav.msm


###################################################
### code chunk number 13: msm-manual.Rnw:1161-1163
###################################################
cavsex.msm <- msm( state ~ years, subject=PTNUM, data = cav,
                     qmatrix = twoway4.q, death = 4, covariates = ~ sex)


###################################################
### code chunk number 14: msm-manual.Rnw:1169-1170
###################################################
cavsex.msm


###################################################
### code chunk number 15: msm-manual.Rnw:1185-1187
###################################################
qmatrix.msm(cavsex.msm, covariates=list(sex=0)) # Male
qmatrix.msm(cavsex.msm, covariates=list(sex=1)) # Female


###################################################
### code chunk number 16: msm-manual.Rnw:1199-1202 (eval = FALSE)
###################################################
## cavsex.msm <- msm( state ~ years, subject=PTNUM, data = cav,
##                   qmatrix = twoway4.q, death = 4,
##                   covariates = list("1-2" = ~ sex, "3-4" = ~sex) )


###################################################
### code chunk number 17: msm-manual.Rnw:1213-1217 (eval = FALSE)
###################################################
## cav3.msm <- msm( state ~ years, subject=PTNUM, data = cav,
##                 qmatrix = twoway4.q, death = 4,
##                 covariates = ~ sex, 
##                 constraint = list(sex=c(1,2,3,1,2,3,2)) )


###################################################
### code chunk number 18: msm-manual.Rnw:1253-1257 (eval = FALSE)
###################################################
## cav4.msm <- msm( state ~ years, subject=PTNUM, data = cav,
##                 qmatrix = twoway4.q, death = 4,
##                 control = list(trace=2, REPORT=1),
##                 fixedpars = c(6, 7) )


###################################################
### code chunk number 19: msm-manual.Rnw:1296-1297
###################################################
pmatrix.msm(cav.msm, t=10)


###################################################
### code chunk number 20: msm-manual.Rnw:1326-1327
###################################################
sojourn.msm(cav.msm)


###################################################
### code chunk number 21: msm-manual.Rnw:1339-1340
###################################################
pnext.msm(cav.msm)


###################################################
### code chunk number 22: msm-manual.Rnw:1365-1366
###################################################
totlos.msm(cav.msm)


###################################################
### code chunk number 23: msm-manual.Rnw:1383-1384
###################################################
qratio.msm(cav.msm, ind1=c(2,1), ind2=c(1,2))


###################################################
### code chunk number 24: msm-manual.Rnw:1395-1396
###################################################
hazard.msm(cavsex.msm)


###################################################
### code chunk number 25: msm-manual.Rnw:1405-1406 (eval = FALSE)
###################################################
## qmatrix.msm(cav.msm)


###################################################
### code chunk number 26: msm-manual.Rnw:1415-1416 (eval = FALSE)
###################################################
## qmatrix.msm(cavsex.msm, covariates = 0)


###################################################
### code chunk number 27: msm-manual.Rnw:1421-1422 (eval = FALSE)
###################################################
## qmatrix.msm(cavsex.msm, covariates = list(sex = 1))


###################################################
### code chunk number 28: msm-manual.Rnw:1448-1449
###################################################
plot(cav.msm, legend.pos=c(8, 1))


###################################################
### code chunk number 29: msm-manual.Rnw:1685-1687
###################################################
options(digits=3)
prevalence.msm(cav.msm, times=seq(0,20,2))


###################################################
### code chunk number 30: msm-manual.Rnw:1689-1690
###################################################
plot.prevalence.msm(cav.msm, mintime=0, maxtime=20)


###################################################
### code chunk number 31: msm-manual.Rnw:1823-1826
###################################################
options(digits=2)
pearson.msm(cav.msm, timegroups=2, 
            transitions=c(1,2,3,4,5,6,7,8,9,9,9,10))


###################################################
### code chunk number 32: msm-manual.Rnw:1951-1963
###################################################
oneway4.q <- rbind(c(0, 0.148, 0, 0.0171),
                     c(0, 0, 0.202, 0.081),
                     c(0, 0, 0, 0.126),
                     c(0, 0, 0, 0))
ematrix <- rbind(c(0, 0.1, 0, 0),
                 c(0.1, 0, 0.1, 0),
                 c(0, 0.1, 0, 0),
                 c(0, 0, 0, 0))
cavmisc.msm <- msm(state ~ years, subject = PTNUM, data = cav,
                   qmatrix = oneway4.q, ematrix = ematrix, death = 4,
                   obstrue = firstobs, method="BFGS")
cavmisc.msm


###################################################
### code chunk number 33: msm-manual.Rnw:1991-1995
###################################################
cavmiscsex.msm <- msm(state ~ years, subject = PTNUM, data = cav,
                      qmatrix = oneway4.q, ematrix = ematrix, 
                      death = 4, misccovariates = ~sex, 
                      obstrue=firstobs, method="BFGS")


###################################################
### code chunk number 34: msm-manual.Rnw:1997-1998
###################################################
cavmiscsex.msm


###################################################
### code chunk number 35: msm-manual.Rnw:2018-2020
###################################################
ematrix.msm(cavmiscsex.msm, covariates=list(sex=0))
ematrix.msm(cavmiscsex.msm, covariates=list(sex=1))


###################################################
### code chunk number 36: msm-manual.Rnw:2067-2069
###################################################
pearson.msm(cavmisc.msm, timegroups=2, 
            transitions=c(1,2,3,4,5,6,7,8,9,9,9,10))


###################################################
### code chunk number 37: msm-manual.Rnw:2115-2117
###################################################
vit <- viterbi.msm(cavmisc.msm)
vit[vit$subject==100103,]


###################################################
### code chunk number 38: msm-manual.Rnw:2315-2316
###################################################
three.q <- rbind(c(0, exp(-6), exp(-9)), c(0, 0, exp(-6)), c(0, 0, 0))


###################################################
### code chunk number 39: msm-manual.Rnw:2334-2346
###################################################
hmodel1 <- list(hmmNorm(mean=100, sd=16), hmmNorm(mean=54, sd=18),
                hmmIdent(999))

fev1.msm <- msm(fev ~ days, subject=ptnum, data=fev, qmatrix=three.q,
                death=3, hmodel=hmodel1,
                hcovariates=list(~acute, ~acute, NULL),
                hcovinits = list(-8, -8, NULL),
                hconstraint = list(acute = c(1,1)), method="BFGS")

fev1.msm

sojourn.msm(fev1.msm)


###################################################
### code chunk number 40: msm-manual.Rnw:2591-2592 (eval = FALSE)
###################################################
## help(msm)


###################################################
### code chunk number 41: msm-manual.Rnw:2600-2601 (eval = FALSE)
###################################################
## help.start()


