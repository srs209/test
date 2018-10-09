library(bimixt)
library(car)

library(pls)
library(ranger)

##functions for boxcox
##functions for boxcox transformation
for.trans <- function(x){
  lam <- powerTransform(x[x>1e-30], family="bcPower")$lambda[[1]]
  test1 <- boxcox(x+1, lambda = lam )
  return(list("var" = test1, "lam"=lam))
}

back.trans <- function(x,lam){   #x is a transformed vector, while is an untransformed vector
  x1 <- boxcox.inv(x, lambda = lam) - 1
  return(x1)
}



##log transformation
for.log <- function(x){
  trans.tmp <- log10(x+1)
  return(trans.tmp)
}

back.log <- function(x){
  back.tmp <- (10^x) - 1 
  return(back.tmp)
}


##NOT RUN ****** 
#do transformation and save it
##OC
load("C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/full_data_trun6000/sub.calib.oc.RData")
load("C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/full_data_trun6000/sub.valid.oc.RData")
sub.calib.oc$oc <- sub.calib.oc$OC
sub.valid.oc$oc <- sub.valid.oc$OC
sub.calib.oc$bc.oc <- for.trans(sub.calib.oc$oc)$var
sub.calib.oc$log.oc <- for.log(sub.calib.oc$oc)
lam <- for.trans(sub.calib.oc$OC)$lam
lam
save(sub.calib.oc, file = "C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/newResults/trans/data/sub.calib.oc.RData")
save(sub.valid.oc, file = "C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/newResults/trans/data/sub.valid.oc.RData")


##pH
load("C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/full_data_trun6000/sub.calib.ph.RData")
load("C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/full_data_trun6000/sub.valid.ph.RData")
sub.calib.ph$bc.ph <- for.trans(sub.calib.ph$ph)$var
sub.calib.ph$log.ph <- for.log(sub.calib.ph$ph)
lam <- for.trans(sub.calib.ph$ph)$lam
lam
save(sub.calib.ph, file = "C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/newResults/trans/data/sub.calib.ph.RData")
save(sub.valid.ph, file = "C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/newResults/trans/data/sub.valid.ph.RData")
##END NOT RUN **********

##get squareroot transformation results from other folder
## we only run no transformation and bc transformation for oc and ph


##start building models by loading data
load("C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/newResults/trans/data/sub.calib.oc.RData")
load("C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/newResults/trans/data/sub.valid.oc.RData")
load("C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/newResults/trans/data/sub.calib.ph.RData")
load("C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/newResults/trans/data/sub.valid.ph.RData")
##1. no transformation model --OC
## total 4 models -- sbl, rf, cubist, plsr

#1a plsr OC Untrans
fit.plsr.untrans.oc <- plsr(oc~spc, ncomp=20, data = sub.calib.oc, valid="CV")
save(fit.plsr.untrans.oc, file = "C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/newResults/trans/plsr/fit.plsr.untrans.oc.RData")
ncomp.onesigma <- selectNcomp(fit.plsr.untrans.oc, method = "onesigma", plot = TRUE,ylim = c(0, .2))
ncomp.onesigma
val.pred <- predict(fit.plsr.untrans.oc, newdata = sub.valid.oc$spc, ncomp=ncomp.onesigma)
pls.valid.pred.oc <- data.frame(sub.valid.oc$oc, val.pred)
names(pls.valid.pred.oc) <- c("obs", "pred")
write.csv(pls.valid.pred.oc, file ="C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/newResults/trans/plsr/pls.untrans.valid.oc.csv")

#1b rf OC untrans
X1 <- data.frame(sub.calib.oc$spc)
Y1 <- sub.calib.oc$oc
X1 <- X1[!is.na(Y1), ]
Y1 <- Y1[!is.na(Y1)]

tmp.calib <- cbind(Y1,X1)

fit.rf.untrans.oc <- ranger(Y1~., data = tmp.calib, quantreg=TRUE, keep.inbag=TRUE, num.trees=150)
save(fit.rf.untrans.oc, file = "C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/newResults/trans/rf/fit.rf.untrans.oc.RData")

pred.new <- predict(fit.rf.untrans.oc, data = sub.valid.oc$spc, type = "se")
valid.pred.oc <- data.frame(sub.valid.oc$oc,pred.new$predictions,pred.new$se)
names(valid.pred.oc) <- c("obs", "pred", "se")
write.csv(valid.pred.oc, file ="C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/newResults/trans/rf/rf.untrans.valid.oc.csv")


#1c sbl OC untrans
load("/home/sdangal/test/localreg/sub.calib.oc.RData")
load("/home/sdangal/test/localreg/sub.valid.oc.RData")

Xu <- sub.valid.oc$spc
Yu <- sub.valid.oc$oc 
Yr <- sub.calib.oc$oc
Xr <- sub.calib.oc$spc
Xu <- Xu[!is.na(Yu),]
Yu <- Yu[!is.na(Yu)]
Xr <- Xr[!is.na(Yr),]
Yr <- Yr[!is.na(Yr)]

ctrl <- mblControl(sm = 'pc', pcSelection = list('opc', 50),
                   valMethod = 'loc_crossval',center=TRUE,scale=FALSE,allowParallel=FALSE)

sbl.untrans.oc <- mbl(Yr = Yr, Xr = Xr, Yu = Yu, Xu = Xu,
                      mblCtrl = ctrl,
                      dissUsage = 'none',
                      k = seq(40, 100, by = 20),
                      method = 'pls', pls.c = 6)
save(sbl.untrans.oc, file = "/home/sdangal/test/localreg/sbl.untrans.oc.RData")
obs <- sbl.untrans.oc$results$Nearest_neighbours_40$yu.obs
pred <- sbl.untrans.oc$results$Nearest_neighbours_40$pred
valid.pred.oc <- data.frame(obs, pred)
write.csv(valid.pred.oc, file = "/home/sdangal/test/localreg/sbl.untrans.valid.oc.csv")
rm(sub.calib.oc, sub.valid.oc, sbl.untrans.oc)






#1d cubist OC untrans
load("/home/sdangal/test/localreg/sub.calib.oc.RData")
load("/home/sdangal/test/localreg/sub.valid.oc.RData")
resp <- sub.calib.oc$oc
sub.calib.oc <- sub.calib.oc[!is.na(resp),]
resp <- resp[!is.na(resp)]
cub.untrans.oc <- cubist(x=sub.calib.oc$spc, y = resp)
save(cub.untrans.oc, file = "/home/sdangal/test/localreg/cub.untrans.oc.RData")
valid.pred <- predict(cub.untrans.oc, sub.valid.oc$spc)
valid.pred.oc <- data.frame(sub.valid.oc$oc, valid.pred)
write.csv(valid.pred.oc, file = "/home/sdangal/test/localreg/cub.untrans.valid.oc.csv")


##2. no transformation model --ph
## total 4 models -- sbl, rf, cubist, plsr

#2a plsr ph Untrans
fit.plsr.untrans.ph <- plsr(ph~spc, ncomp=20, data = sub.calib.ph, valid="CV")
save(fit.plsr.untrans.ph, file = "C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/newResults/trans/plsr/fit.plsr.untrans.ph.RData")
ncomp.onesigma <- selectNcomp(fit.plsr.untrans.ph, method = "onesigma", plot = TRUE,ylim = c(0, .2))
ncomp.onesigma
val.pred <- predict(fit.plsr.untrans.ph, newdata = sub.valid.ph$spc, ncomp=ncomp.onesigma)
pls.valid.pred.ph <- data.frame(sub.valid.ph$ph, val.pred)
names(pls.valid.pred.ph) <- c("obs", "pred")
write.csv(pls.valid.pred.ph, file ="C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/newResults/trans/plsr/pls.untrans.valid.ph.csv")

#2b rf ph untrans
X1 <- data.frame(sub.calib.ph$spc)
Y1 <- sub.calib.ph$ph
X1 <- X1[!is.na(Y1), ]
Y1 <- Y1[!is.na(Y1)]

tmp.calib <- cbind(Y1,X1)

fit.rf.untrans.ph <- ranger(Y1~., data = tmp.calib, quantreg=TRUE, keep.inbag=TRUE, num.trees=150)
save(fit.rf.untrans.ph, file = "C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/newResults/trans/rf/fit.rf.untrans.ph.RData")

pred.new <- predict(fit.rf.untrans.ph, data = sub.valid.ph$spc, type = "se")
valid.pred.ph <- data.frame(sub.valid.ph$ph,pred.new$predictions,pred.new$se)
names(valid.pred.ph) <- c("obs", "pred", "se")
write.csv(valid.pred.ph, file ="C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/newResults/trans/rf/rf.untrans.valid.ph.csv")


#2c sbl ph untrans
load("/home/sdangal/test/localreg/sub.calib.ph.RData")
load("/home/sdangal/test/localreg/sub.valid.ph.RData")

Xu <- sub.valid.ph$spc
Yu <- sub.valid.ph$ph 
Yr <- sub.calib.ph$ph
Xr <- sub.calib.ph$spc
Xu <- Xu[!is.na(Yu),]
Yu <- Yu[!is.na(Yu)]
Xr <- Xr[!is.na(Yr),]
Yr <- Yr[!is.na(Yr)]

ctrl <- mblControl(sm = 'pc', pcSelection = list('opc', 50),
                   valMethod = 'loc_crossval',center=TRUE,scale=FALSE,allowParallel=FALSE)

sbl.untrans.ph <- mbl(Yr = Yr, Xr = Xr, Yu = Yu, Xu = Xu,
                      mblCtrl = ctrl,
                      dissUsage = 'none',
                      k = seq(40, 100, by = 20),
                      method = 'pls', pls.c = 6)
save(sbl.untrans.ph, file = "/home/sdangal/test/localreg/sbl.untrans.ph.RData")
obs <- sbl.untrans.ph$results$Nearest_neighbours_40$yu.obs
pred <- sbl.untrans.ph$results$Nearest_neighbours_40$pred
valid.pred.ph <- data.frame(obs, pred)
write.csv(valid.pred.ph, file = "/home/sdangal/test/localreg/sbl.untrans.valid.ph.csv")
rm(sub.calib.ph, sub.valid.ph, sbl.untrans.ph)





#2d cubist ph untrans
load("/home/sdangal/test/localreg/sub.calib.ph.RData")
load("/home/sdangal/test/localreg/sub.valid.ph.RData")
resp <- sub.calib.ph$ph
sub.calib.ph <- sub.calib.ph[!is.na(resp),]
resp <- resp[!is.na(resp)]
cub.untrans.ph <- cubist(x=sub.calib.ph$spc, y = resp)
save(cub.untrans.ph, file = "/home/sdangal/test/localreg/cub.untrans.ph.RData")
valid.pred <- predict(cub.untrans.ph, sub.valid.ph$spc)
valid.pred.ph <- data.frame(sub.valid.ph$ph, valid.pred)
write.csv(valid.pred.ph, file = "/home/sdangal/test/localreg/cub.untrans.valid.ph.csv")


##3. box cox transformation model --OC
## total 4 models -- sbl, rf, cubist, plsr

#3a plsr OC bc
fit.plsr.bc.oc <- plsr(bc.oc~spc, ncomp=20, data = sub.calib.oc, valid="CV")
save(fit.plsr.bc.oc, file = "C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/newResults/trans/plsr/fit.plsr.bc.oc.RData")
ncomp.onesigma <- selectNcomp(fit.plsr.bc.oc, method = "onesigma", plot = TRUE,ylim = c(0, .2))
ncomp.onesigma
val.pred <- predict(fit.plsr.bc.oc, newdata = sub.valid.oc$spc, ncomp=ncomp.onesigma)
##backtrans
val.pred <- back.trans(val.pred, -0.02363157)
pls.valid.pred.oc <- data.frame(sub.valid.oc$oc, val.pred)
names(pls.valid.pred.oc) <- c("obs", "pred")
write.csv(pls.valid.pred.oc, file ="C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/newResults/trans/plsr/pls.bc.valid.oc.csv")

#3b random forest model for bc trans
X1 <- data.frame(sub.calib.oc$spc)
Y1 <- sub.calib.oc$bc.oc
X1 <- X1[!is.na(Y1), ]
Y1 <- Y1[!is.na(Y1)]
tmp.calib <- cbind(Y1,X1)
fit.rf.bc.oc <- ranger(Y1~., data = tmp.calib, quantreg=TRUE, keep.inbag=TRUE, num.trees=150)
save(fit.rf.bc.oc, file = "C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/newResults/trans/rf/fit.rf.bc.oc.RData")

pred.new <- predict(fit.rf.bc.oc, data = sub.valid.oc$spc, type = "se")
pred.new <- back.trans(pred.new$predictions, -0.02363157)
valid.pred.oc <- data.frame(sub.valid.oc$oc,pred.new)
names(valid.pred.oc) <- c("obs", "pred")
write.csv(valid.pred.oc, file ="C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/newResults/trans/rf/rf.bc.valid.oc.csv")


#3c sbl model for bc trans oc
load("/home/sdangal/test/localreg/sub.calib.oc.RData")
load("/home/sdangal/test/localreg/sub.valid.oc.RData")

Xu <- sub.valid.oc$spc
Yu <- sub.valid.oc$bc.oc 
Yr <- sub.calib.oc$bc.oc
Xr <- sub.calib.oc$spc
Xu <- Xu[!is.na(Yu),]
Yu <- Yu[!is.na(Yu)]
Xr <- Xr[!is.na(Yr),]
Yr <- Yr[!is.na(Yr)]

ctrl <- mblControl(sm = 'pc', pcSelection = list('opc', 50),
                   valMethod = 'loc_crossval',center=TRUE,scale=FALSE,allowParallel=FALSE)

sbl.bc.oc <- mbl(Yr = Yr, Xr = Xr, Yu = Yu, Xu = Xu,
                 mblCtrl = ctrl,
                 dissUsage = 'none',
                 k = seq(40, 100, by = 20),
                 method = 'pls', pls.c = 6)
save(sbl.bc.oc, file = "/home/sdangal/test/localreg/sbl.bc.oc.RData")
obs <- sbl.bc.oc$results$Nearest_neighbours_40$yu.obs
pred <- sbl.bc.oc$results$Nearest_neighbours_40$pred
obs <- back.trans(obs, -0.02363157)
pred <- back.trans(pred, -0.02363157)
valid.pred.oc <- data.frame(obs, pred)
write.csv(valid.pred.oc, file = "/home/sdangal/test/localreg/sbl.bc.valid.oc.csv")
rm(sub.calib.oc, sub.valid.oc, sbl.bc.oc)




#3d cubist model for bc trans oc
load("/home/sdangal/test/localreg/sub.calib.oc.RData")
load("/home/sdangal/test/localreg/sub.valid.oc.RData")
resp <- sub.calib.oc$bc.oc
sub.calib.oc <- sub.calib.oc[!is.na(resp),]
resp <- resp[!is.na(resp)]
cub.bc.oc <- cubist(x=sub.calib.oc$spc, y = resp)
save(cub.bc.oc, file = "/home/sdangal/test/localreg/cub.bc.oc.RData")
valid.pred <- predict(cub.bc.oc, sub.valid.oc$spc)
valid.new <- back.trans(valid.new, -0.02363157)
valid.pred.oc <- data.frame(sub.valid.oc$oc, valid.pred)
write.csv(valid.pred.oc, file = "/home/sdangal/test/localreg/cub.bc.valid.oc.csv")



##4. box cox transformation model --ph
## total 4 models -- sbl, rf, cubist, plsr

#4a plsr
fit.plsr.bc.ph <- plsr(bc.ph~spc, ncomp=20, data = sub.calib.ph, valid="CV")
save(fit.plsr.bc.ph, file = "C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/newResults/trans/plsr/fit.plsr.bc.ph.RData")
ncomp.onesigma <- selectNcomp(fit.plsr.bc.ph, method = "onesigma", plot = TRUE,ylim = c(0, .2))
ncomp.onesigma
val.pred <- predict(fit.plsr.bc.ph, newdata = sub.valid.ph$spc, ncomp=ncomp.onesigma)
##backtrans
val.pred <- back.trans(val.pred, 0.6186485)
pls.valid.pred.ph <- data.frame(sub.valid.ph$ph, val.pred)
names(pls.valid.pred.ph) <- c("obs", "pred")
write.csv(pls.valid.pred.ph, file ="C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/newResults/trans/plsr/pls.bc.valid.ph.csv")

#4b random forest model
X1 <- data.frame(sub.calib.ph$spc)
Y1 <- sub.calib.ph$bc.ph
X1 <- X1[!is.na(Y1), ]
Y1 <- Y1[!is.na(Y1)]

tmp.calib <- cbind(Y1,X1)

fit.rf.bc.ph <- ranger(Y1~., data = tmp.calib, quantreg=TRUE, keep.inbag=TRUE, num.trees=150)
save(fit.rf.bc.ph, file = "C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/newResults/trans/rf/fit.rf.bc.ph.RData")

pred.new <- predict(fit.rf.bc.ph, data = sub.valid.ph$spc, type = "se")
pred.new <- back.trans(pred.new$predictions, 0.6186485)
valid.pred.ph <- data.frame(sub.valid.ph$ph,pred.new)
names(valid.pred.ph) <- c("obs", "pred")
write.csv(valid.pred.ph, file ="C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/newResults/trans/rf/rf.bc.valid.ph.csv")



#4c sbl ph bc

load("/home/sdangal/test/localreg/sub.calib.ph.RData")
load("/home/sdangal/test/localreg/sub.valid.ph.RData")

Xu <- sub.valid.ph$spc
Yu <- sub.valid.ph$bc.ph 
Yr <- sub.calib.ph$bc.ph
Xr <- sub.calib.ph$spc
Xu <- Xu[!is.na(Yu),]
Yu <- Yu[!is.na(Yu)]
Xr <- Xr[!is.na(Yr),]
Yr <- Yr[!is.na(Yr)]

ctrl <- mblControl(sm = 'pc', pcSelection = list('opc', 50),
                   valMethod = 'loc_crossval',center=TRUE,scale=FALSE,allowParallel=FALSE)

sbl.bc.ph <- mbl(Yr = Yr, Xr = Xr, Yu = Yu, Xu = Xu,
                 mblCtrl = ctrl,
                 dissUsage = 'none',
                 k = seq(40, 100, by = 20),
                 method = 'pls', pls.c = 6)
save(sbl.bc.ph, file = "/home/sdangal/test/localreg/sbl.bc.ph.RData")
obs <- sbl.bc.ph$results$Nearest_neighbours_40$yu.obs
pred <- sbl.bc.ph$results$Nearest_neighbours_40$pred
obs <- back.trans(obs, 0.6186485)
pred <- back.trans(pred, 0.6186485)
valid.pred.ph <- data.frame(obs, pred)
write.csv(valid.pred.ph, file = "/home/sdangal/test/localreg/sbl.bc.valid.ph.csv")
rm(sub.calib.ph, sub.valid.ph, sbl.bc.ph)




#4d cubist ph bc
load("/home/sdangal/test/localreg/sub.calib.ph.RData")
load("/home/sdangal/test/localreg/sub.valid.ph.RData")
resp <- sub.calib.ph$bc.ph
sub.calib.ph <- sub.calib.ph[!is.na(resp),]
resp <- resp[!is.na(resp)]
cub.bc.ph <- cubist(x=sub.calib.ph$spc, y = resp)
save(cub.bc.ph, file = "/home/sdangal/test/localreg/cub.bc.ph.RData")
valid.pred <- predict(cub.bc.ph, sub.valid.ph$spc)
valid.new <- back.trans(valid.new, 0.6186485)
valid.pred.ph <- data.frame(sub.valid.ph$ph, valid.pred)
write.csv(valid.pred.ph, file = "/home/sdangal/test/localreg/cub.bc.valid.ph.csv")


##5. log transformation model --OC
## total 4 models -- sbl, rf, cubist, plsr

#5a plsr OC Untrans
fit.plsr.log.oc <- plsr(log.oc~spc, ncomp=20, data = sub.calib.oc, valid="CV")
save(fit.plsr.log.oc, file = "C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/newResults/trans/plsr/fit.plsr.log.oc.RData")
ncomp.onesigma <- selectNcomp(fit.plsr.log.oc, method = "onesigma", plot = TRUE,ylim = c(0, .2))
ncomp.onesigma
val.pred <- predict(fit.plsr.log.oc, newdata = sub.valid.oc$spc, ncomp=ncomp.onesigma)
##backtrans
val.pred <- back.log(val.pred)
pls.valid.pred.oc <- data.frame(sub.valid.oc$oc, val.pred)
names(pls.valid.pred.oc) <- c("obs", "pred")
write.csv(pls.valid.pred.oc, file ="C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/newResults/trans/plsr/pls.log.valid.oc.csv")

#5b random forest model for log trans
X1 <- data.frame(sub.calib.oc$spc)
Y1 <- sub.calib.oc$log.oc
X1 <- X1[!is.na(Y1), ]
Y1 <- Y1[!is.na(Y1)]
tmp.calib <- cbind(Y1,X1)
fit.rf.log.oc <- ranger(Y1~., data = tmp.calib, quantreg=TRUE, keep.inbag=TRUE, num.trees=150)
save(fit.rf.log.oc, file = "C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/newResults/trans/rf/fit.rf.log.oc.RData")

pred.new <- predict(fit.rf.log.oc, data = sub.valid.oc$spc, type = "se")
pred.new <- back.log(pred.new$predictions)
valid.pred.oc <- data.frame(sub.valid.oc$oc,pred.new)
names(valid.pred.oc) <- c("obs", "pred")
write.csv(valid.pred.oc, file ="C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/newResults/trans/rf/rf.log.valid.oc.csv")


#5c sbl model for log trans oc

load("/home/sdangal/test/localreg/sub.calib.oc.RData")
load("/home/sdangal/test/localreg/sub.valid.oc.RData")

Xu <- sub.valid.oc$spc
Yu <- sub.valid.oc$log.oc 
Yr <- sub.calib.oc$log.oc
Xr <- sub.calib.oc$spc
Xu <- Xu[!is.na(Yu),]
Yu <- Yu[!is.na(Yu)]
Xr <- Xr[!is.na(Yr),]
Yr <- Yr[!is.na(Yr)]

ctrl <- mblControl(sm = 'pc', pcSelection = list('opc', 50),
                   valMethod = 'loc_crossval',center=TRUE,scale=FALSE,allowParallel=FALSE)

sbl.log.oc <- mbl(Yr = Yr, Xr = Xr, Yu = Yu, Xu = Xu,
                  mblCtrl = ctrl,
                  dissUsage = 'none',
                  k = seq(40, 100, by = 20),
                  method = 'pls', pls.c = 6)
save(sbl.log.oc, file = "/home/sdangal/test/localreg/sbl.log.oc.RData")
obs <- sbl.log.oc$results$Nearest_neighbours_40$yu.obs
pred <- sbl.log.oc$results$Nearest_neighbours_40$pred
obs <- back.log(obs)
pred <- back.log(pred)
valid.pred.oc <- data.frame(obs, pred)
write.csv(valid.pred.oc, file = "/home/sdangal/test/localreg/sbl.log.valid.oc.csv")
rm(sub.calib.oc, sub.valid.oc, sbl.log.oc)





#5d cubist model for log trans oc
load("/home/sdangal/test/localreg/sub.calib.oc.RData")
load("/home/sdangal/test/localreg/sub.valid.oc.RData")
resp <- sub.calib.oc$log.oc
sub.calib.oc <- sub.calib.oc[!is.na(resp),]
resp <- resp[!is.na(resp)]
cub.log.oc <- cubist(x=sub.calib.oc$spc, y = resp)
save(cub.log.oc, file = "/home/sdangal/test/localreg/cub.log.oc.RData")
valid.pred <- predict(cub.log.oc, sub.valid.oc$spc)
valid.new <- back.log(valid.new)
valid.pred.oc <- data.frame(sub.valid.oc$oc, valid.pred)
write.csv(valid.pred.oc, file = "/home/sdangal/test/localreg/cub.log.valid.oc.csv")



##5. log transformation model --ph
## total 4 models -- sbl, rf, cubist, plsr

#6a plsr
fit.plsr.log.ph <- plsr(log.ph~spc, ncomp=20, data = sub.calib.ph, valid="CV")
save(fit.plsr.log.ph, file = "C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/newResults/trans/plsr/fit.plsr.log.ph.RData")
ncomp.onesigma <- selectNcomp(fit.plsr.log.ph, method = "onesigma", plot = TRUE,ylim = c(0, .2))
ncomp.onesigma
val.pred <- predict(fit.plsr.log.ph, newdata = sub.valid.ph$spc, ncomp=ncomp.onesigma)
##backtrans
val.pred <- back.log(val.pred)
pls.valid.pred.ph <- data.frame(sub.valid.ph$ph, val.pred)
names(pls.valid.pred.ph) <- c("obs", "pred")
write.csv(pls.valid.pred.ph, file ="C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/newResults/trans/plsr/pls.log.valid.ph.csv")

#6b random forest model
X1 <- data.frame(sub.calib.ph$spc)
Y1 <- sub.calib.ph$log.ph
X1 <- X1[!is.na(Y1), ]
Y1 <- Y1[!is.na(Y1)]

tmp.calib <- cbind(Y1,X1)

fit.rf.log.ph <- ranger(Y1~., data = tmp.calib, quantreg=TRUE, keep.inbag=TRUE, num.trees=150)
save(fit.rf.log.ph, file = "C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/newResults/trans/rf/fit.log.bc.ph.RData")

pred.new <- predict(fit.rf.log.ph, data = sub.valid.ph$spc, type = "se")
pred.new <- back.log(pred.new$predictions)
valid.pred.ph <- data.frame(sub.valid.ph$ph,pred.new)
names(valid.pred.ph) <- c("obs", "pred")
write.csv(valid.pred.ph, file ="C:/Users/sdangal/Documents/ModelChoice_17SEPT2018/newResults/trans/rf/rf.log.valid.ph.csv")



#6c sbl ph log
load("/home/sdangal/test/localreg/sub.calib.ph.RData")
load("/home/sdangal/test/localreg/sub.valid.ph.RData")

Xu <- sub.valid.ph$spc
Yu <- sub.valid.ph$log.ph 
Yr <- sub.calib.ph$log.ph
Xr <- sub.calib.ph$spc
Xu <- Xu[!is.na(Yu),]
Yu <- Yu[!is.na(Yu)]
Xr <- Xr[!is.na(Yr),]
Yr <- Yr[!is.na(Yr)]

ctrl <- mblControl(sm = 'pc', pcSelection = list('opc', 50),
                   valMethod = 'loc_crossval',center=TRUE,scale=FALSE,allowParallel=FALSE)

sbl.log.ph <- mbl(Yr = Yr, Xr = Xr, Yu = Yu, Xu = Xu,
                  mblCtrl = ctrl,
                  dissUsage = 'none',
                  k = seq(40, 100, by = 20),
                  method = 'pls', pls.c = 6)
save(sbl.log.ph, file = "/home/sdangal/test/localreg/sbl.log.ph.RData")
obs <- sbl.log.ph$results$Nearest_neighbours_40$yu.obs
pred <- sbl.log.ph$results$Nearest_neighbours_40$pred
obs <- back.log(obs)
pred <- back.log(pred)
valid.pred.ph <- data.frame(obs, pred)
write.csv(valid.pred.ph, file = "/home/sdangal/test/localreg/sbl.log.valid.ph.csv")
rm(sub.calib.ph, sub.valid.ph, sbl.log.ph)




#6d cubist ph log
load("/home/sdangal/test/localreg/sub.calib.ph.RData")
load("/home/sdangal/test/localreg/sub.valid.ph.RData")
resp <- sub.calib.ph$log.ph
sub.calib.ph <- sub.calib.ph[!is.na(resp),]
resp <- resp[!is.na(resp)]
cub.log.ph <- cubist(x=sub.calib.ph$spc, y = resp)
save(cub.log.ph, file = "/home/sdangal/test/localreg/cub.log.ph.RData")
valid.pred <- predict(cub.log.ph, sub.valid.ph$spc)
valid.pred <- back.log(valid.pred)
valid.pred.ph <- data.frame(sub.valid.ph$ph, valid.pred)
write.csv(valid.pred.ph, file = "/home/sdangal/test/localreg/cub.log.valid.ph.csv")








