### R code from vignette source 'msm-manual.Rnw'

###################################################
### code chunk number 1: msm-manual.Rnw:20-22
###################################################
version <- gsub("Version: +", "",
                packageDescription("msm", lib.loc="../../..")$Version)


###################################################
### code chunk number 2: msm-manual.Rnw:27-28
###################################################
cat(version)


###################################################
### code chunk number 3: msm-manual.Rnw:31-32
###################################################
cat(format(Sys.time(), "%d %B, %Y"))


###################################################
### code chunk number 4: msm-manual.Rnw:842-843
###################################################
options(width = 60)


###################################################
### code chunk number 5: msm-manual.Rnw:878-879
###################################################
library(msm)


###################################################
### code chunk number 6: msm-manual.Rnw:921-922
###################################################
data(cav)


###################################################
### code chunk number 7: msm-manual.Rnw:943-944
###################################################
cav[1:21,]


###################################################
### code chunk number 8: msm-manual.Rnw:954-955
###################################################
statetable.msm(state, PTNUM, data=cav)


###################################################
### code chunk number 9: msm-manual.Rnw:1006-1010
###################################################
twoway4.q  <-  rbind ( c(0, 0.25, 0, 0.25),
                       c(0.166, 0, 0.166, 0.166),
                       c(0, 0.25, 0, 0.25),
                       c(0, 0, 0, 0) )


###################################################
### code chunk number 10: msm-manual.Rnw:1039-1040
###################################################
twoway4.crude.q  <- crudeinits.msm(state ~ years, PTNUM, data=cav, qmatrix=twoway4.q)


###################################################
### code chunk number 11: msm-manual.Rnw:1065-1067
###################################################
cav.msm <- msm( state ~ years, subject=PTNUM, data = cav,
                    qmatrix = twoway4.q, death = 4 )


###################################################
### code chunk number 12: msm-manual.Rnw:1095-1096 (eval = FALSE)
###################################################
## help(optim)


###################################################
### code chunk number 13: msm-manual.Rnw:1118-1119
###################################################
cav.msm


###################################################
### code chunk number 14: msm-manual.Rnw:1157-1159
###################################################
cavsex.msm <- msm( state ~ years, subject=PTNUM, data = cav,
                     qmatrix = twoway4.q, death = 4, covariates = ~ sex)


###################################################
### code chunk number 15: msm-manual.Rnw:1165-1166
###################################################
cavsex.msm


###################################################
### code chunk number 16: msm-manual.Rnw:1181-1183
###################################################
qmatrix.msm(cavsex.msm, covariates=list(sex=0)) # Male
qmatrix.msm(cavsex.msm, covariates=list(sex=1)) # Female


###################################################
### code chunk number 17: msm-manual.Rnw:1194-1197 (eval = FALSE)
###################################################
## cav3.msm <- msm( state ~ years, subject=PTNUM, data = cav,
##                    qmatrix = twoway4.q, death = 4,
##                    covariates = ~ sex, constraint = list(sex=c(1,2,3,1,2,3,2)) )


###################################################
### code chunk number 18: msm-manual.Rnw:1233-1237 (eval = FALSE)
###################################################
## cav4.msm <- msm( state ~ years, subject=PTNUM, data = cav,
##                 qmatrix = twoway4.q, death = 4,
##                 control = list(trace=2, REPORT=1),
##                 fixedpars = c(6, 7) )


###################################################
### code chunk number 19: msm-manual.Rnw:1276-1277
###################################################
pmatrix.msm(cav.msm, t=10)


###################################################
### code chunk number 20: msm-manual.Rnw:1306-1307
###################################################
sojourn.msm(cav.msm)


###################################################
### code chunk number 21: msm-manual.Rnw:1319-1320
###################################################
pnext.msm(cav.msm)


###################################################
### code chunk number 22: msm-manual.Rnw:1345-1346
###################################################
totlos.msm(cav.msm)


###################################################
### code chunk number 23: msm-manual.Rnw:1363-1364
###################################################
qratio.msm(cav.msm, ind1=c(2,1), ind2=c(1,2))


###################################################
### code chunk number 24: msm-manual.Rnw:1375-1376
###################################################
hazard.msm(cavsex.msm)


###################################################
### code chunk number 25: msm-manual.Rnw:1385-1386 (eval = FALSE)
###################################################
## qmatrix.msm(cav.msm)


###################################################
### code chunk number 26: msm-manual.Rnw:1395-1396 (eval = FALSE)
###################################################
## qmatrix.msm(cavsex.msm, covariates = 0)


###################################################
### code chunk number 27: msm-manual.Rnw:1401-1402 (eval = FALSE)
###################################################
## qmatrix.msm(cavsex.msm, covariates = list(sex = 1))


###################################################
### code chunk number 28: msm-manual.Rnw:1428-1429
###################################################
plot(cav.msm, legend.pos=c(8, 1))


###################################################
### code chunk number 29: msm-manual.Rnw:1648-1650
###################################################
options(digits=3)
prevalence.msm(cav.msm, times=seq(0,20,2))


###################################################
### code chunk number 30: msm-manual.Rnw:1652-1653
###################################################
plot.prevalence.msm(cav.msm, mintime=0, maxtime=20)


###################################################
### code chunk number 31: msm-manual.Rnw:1786-1788
###################################################
options(digits=2)
pearson.msm(cav.msm, timegroups=2, transitions=c(1,2,3,4,5,6,7,8,9,9,9,10))


###################################################
### code chunk number 32: msm-manual.Rnw:1911-1923
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
### code chunk number 33: msm-manual.Rnw:1951-1954
###################################################
cavmiscsex.msm <- msm(state ~ years, subject = PTNUM, data = cav,
                          qmatrix = oneway4.q, ematrix = ematrix, death = 4,
                          misccovariates = ~sex, obstrue=firstobs, method="BFGS")


###################################################
### code chunk number 34: msm-manual.Rnw:1956-1957
###################################################
cavmiscsex.msm


###################################################
### code chunk number 35: msm-manual.Rnw:1977-1979
###################################################
ematrix.msm(cavmiscsex.msm, covariates=list(sex=0))
ematrix.msm(cavmiscsex.msm, covariates=list(sex=1))


###################################################
### code chunk number 36: msm-manual.Rnw:1988-1989
###################################################
odds.msm(cavmiscsex.msm)


###################################################
### code chunk number 37: msm-manual.Rnw:2026-2027
###################################################
pearson.msm(cavmisc.msm, timegroups=2, transitions=c(1,2,3,4,5,6,7,8,9,9,9,10))


###################################################
### code chunk number 38: msm-manual.Rnw:2073-2075
###################################################
vit <- viterbi.msm(cavmisc.msm)
vit[vit$subject==100103,]


###################################################
### code chunk number 39: msm-manual.Rnw:2261-2262
###################################################
data(fev)


###################################################
### code chunk number 40: msm-manual.Rnw:2267-2268
###################################################
three.q <- rbind(c(0, exp(-6), exp(-9)), c(0, 0, exp(-6)), c(0, 0, 0))


###################################################
### code chunk number 41: msm-manual.Rnw:2286-2295
###################################################
hmodel1 <- list(hmmNorm(mean=100, sd=16), hmmNorm(mean=54, sd=18), hmmIdent(999))

fev1.msm <- msm(fev ~ days, subject=ptnum, data=fev, qmatrix=three.q, death=3, hmodel=hmodel1,
                hcovariates=list(~acute, ~acute, NULL), hcovinits = list(-8, -8, NULL),
                hconstraint = list(acute = c(1,1)), method="BFGS")

fev1.msm

sojourn.msm(fev1.msm)


###################################################
### code chunk number 42: msm-manual.Rnw:2328-2339
###################################################
hmodel2 <- list(hmmMETNorm(mean=100, sd=16, sderr=8, lower=80, upper=Inf, meanerr=0),
                hmmMETNorm(mean=54, sd=18, sderr=8, lower=0, upper=80, meanerr=0),
                hmmIdent(999))

fev2.msm <- msm(fev ~ days, subject=ptnum, data=fev, qmatrix=three.q, death=3, hmodel=hmodel2,
                hcovariates=list(~acute, ~acute, NULL), hcovinits = list(-8, -8, NULL),
                hconstraint = list(sderr = c(1,1), acute = c(1,1)), method="BFGS")

fev2.msm

sojourn.msm(fev2.msm)


###################################################
### code chunk number 43: msm-manual.Rnw:2357-2366
###################################################
keep <- fev$ptnum==1 & fev$fev<999
plot(fev$days[keep], fev$fev[keep], type="l",
ylab=expression(paste("% baseline ", FEV[1])), xlab="Days after transplant")
vit <- viterbi.msm(fev2.msm)[keep,]
(max1 <- max(vit$time[vit$fitted==1]))
(min2 <- min(vit$time[vit$fitted==2]))
abline(v = mean(max1,min2), lty=2)
text(max1 - 500, 50, "STATE 1")
text(min2 + 500, 50, "STATE 2")


###################################################
### code chunk number 44: msm-manual.Rnw:2397-2405
###################################################
oneway4.q <- rbind(c(0, 0.148, 0, 0.0171),
                   c(0, 0, 0.202, 0.081),
                   c(0, 0, 0, 0.126),
                   c(0, 0, 0, 0))
cavmisc.msm <- msm(state ~ years, subject = PTNUM, data = cav,
                     hmodel = list (hmmCat(c(0.9, 0.1, 0, 0)), hmmCat(c(0.1, 0.8, 0.1, 0)), hmmCat(c(0, 0.1, 0.9, 0)), hmmIdent(4)),
                     qmatrix = oneway4.q, obstrue=firstobs, death = 4, method="BFGS")
cavmisc.msm


###################################################
### code chunk number 45: msm-manual.Rnw:2528-2529 (eval = FALSE)
###################################################
## help(msm)


###################################################
### code chunk number 46: msm-manual.Rnw:2537-2538 (eval = FALSE)
###################################################
## help.start()


