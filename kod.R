library(foreign)
library(RJDemetra)
library(tidyverse)
library(tidyquant)
library(forecast)
library(tseries)
library(lmtest)
library(xts)
library(urca)
library(aTSA)


niesezonowe <- read.csv('niesezonowe.csv',
                        header=TRUE)
sezonowe <- read.csv('sezonowe.csv',
                      header=TRUE)

niesezonowe$Date <- as.Date(niesezonowe$Date)
sezonowe$Date <- as.Date(sezonowe$Date)

ts_niesezonowe <-  ts(data=niesezonowe$Value, start=c(2000,1), end=c(2023, 2), frequency=12)
ts_niesezonowe.xts = as.xts(ts_niesezonowe)
ts_sezonowe <- ts(data=sezonowe$Value, start = c(2000, 1), end=c(2022, 12), frequency=12)
ts_sezonowe.xts = as.xts(ts_sezonowe)

# wykresy
plot(ts_niesezonowe, 
     type='l', 
     col='black',
     xlab='Data',
     ylab='Cena kawy (U.S. cents/pound)',
     lwd=1)


plot(ts_sezonowe, 
     type='l', 
     col='black',
     xlab='Data',
     ylab='Liczba pasażerów kolei w USA (w milionach)',
     lwd=1)

##########################
### SZEREG NIESEZONOWY ###
##########################



# dekompozycja szeregu
niesezonowy_da = decompose(ts_niesezonowe, 'additive')
plot(niesezonowy_da)

niesezonowy_dm = decompose(ts_niesezonowe, 'multiplicative')
plot(niesezonowy_dm)

forecast::ggsubseriesplot(ts_niesezonowe, ylab = "")

# ADF

source("funs/TESTDF.R")
x <- ur.df(na.omit(ts_niesezonowe), type = c('drift'), lags=0)
summary(x)


I1_ts_niesezonowe <- diff.xts(ts_niesezonowe)
d.nonseas.is <- window(I1_ts_niesezonowe, end=c(2022,11))
nonseas.is <- window(ts_niesezonowe, end=c(2022,11))

testdf(variable = nonseas.is, ADF_type="c", ADF_max_order = 3,BG_max_order = 4)

testdf(variable = d.nonseas.is, ADF_type="nc", ADF_max_order = 2,BG_max_order = 4)

plot(d.nonseas.is, 
     type='l', 
     col='black',
     xlab='Data',
     ylab='Cena kawy (U.S. cents/pound)',
     lwd=1)

kpss.test <- ur.kpss(I1_ts_niesezonowe, type=c('mu'))
summary(kpss.test)
# ARIMA

Box.test(I1_ts_niesezonowe, type = "Ljung-Box", lag = 24)
Box.test(I1_ts_niesezonowe, type = "Box-Pierce", lag = 24)


ACF <- acf(d.nonseas.is,
    lag.max = 24,
    na.action = na.pass,
    plot=FALSE)

ACF$lag <- ACF$lag * 12

plot(ACF,
     xlim=c(1,24),
     ylim=c(-0.5,0.5))

PACF <- pacf(d.nonseas.is,
    na.action=na.pass,
    plot=FALSE)

PACF$lag <- PACF$lag * 12

plot(PACF,
     xlim=c(1,24),
     ylim=c(-0.5,0.5),
     ylab='PACF')

# ARIMA

d.nonseas.is <- window(I1_ts_niesezonowe, end=c(2022,11))
nonseas.is <- window(ts_niesezonowe, end=c(2022,11))

nobs <- length(nonseas.is)

#######################################
arima111 <- arima(nonseas.is,
                  order = c(1, 1, 1),
                  xreg=1:nobs)

arima111
coeftest(arima111) # MA(1) i AR(1) nieistotne

# badanie reszt
ACF <- acf(resid(arima111),
           lag.max = 24,
           na.action = na.pass,
           plot=FALSE)

ACF$lag <- ACF$lag * 12

plot(ACF,
     xlim=c(1,24),
     ylim=c(-0.5,0.5))

PACF <- pacf(resid(arima111),
             na.action=na.pass,
             plot=FALSE)

PACF$lag <- PACF$lag * 12

plot(PACF,
     xlim=c(1,24),
     ylim=c(-0.5,0.5),
     ylab='PACF')

Box.test(resid(arima111), type = "Ljung-Box", lag = 24)
Box.test(resid(arima111), type = "Box-Pierce", lag = 24)
#######################################
# ARIMA(1,1,0)
arima110 <- arima(nonseas.is,
                  order = c(1, 1, 0))
arima110
coeftest(arima110)

# test LR
# H0 phi1 = 0
teststat <- 2*(as.numeric(logLik(arima111))- as.numeric(logLik(arima110)))
teststat
pchisq(teststat, df=1, lower.tail = FALSE) # brak podstaw do odrzucenia H0

# badanie reszt
ACF <- acf(resid(arima110),
           lag.max = 24,
           na.action = na.pass,
           plot=FALSE)

ACF$lag <- ACF$lag * 12

plot(ACF,
     xlim=c(1,24),
     ylim=c(-0.5,0.5))

PACF <- pacf(resid(arima110),
             na.action=na.pass,
             plot=FALSE)

PACF$lag <- PACF$lag * 12

plot(PACF,
     xlim=c(1,24),
     ylim=c(-0.5,0.5),
     ylab='PACF')

Box.test(resid(arima110), type = "Ljung-Box", lag = 24)
Box.test(resid(arima110), type = "Box-Pierce", lag = 24)

#######################################
# ARIMA(0,1,1)
arima011 <- arima(nonseas.is,
                  order = c(0, 1, 1))
arima011
coeftest(arima011)

# test LR
# H0 alfa1 = 0
teststat <- 2*(as.numeric(logLik(arima111))- as.numeric(logLik(arima011)))
teststat
pchisq(teststat, df=1, lower.tail = FALSE) # brak podstaw do odrzucenia H0

# badanie reszt
ACF <- acf(resid(arima011),
           lag.max = 24,
           na.action = na.pass,
           plot=FALSE)

ACF$lag <- ACF$lag * 12

plot(ACF,
     xlim=c(1,24),
     ylim=c(-0.5,0.5))

PACF <- pacf(resid(arima011),
             na.action=na.pass,
             plot=FALSE)

PACF$lag <- PACF$lag * 12

plot(PACF,
     xlim=c(1,24),
     ylim=c(-0.5,0.5),
     ylab='PACF')

Box.test(resid(arima011), type = "Ljung-Box", lag = 24)
Box.test(resid(arima011), type = "Box-Pierce", lag = 24)

#######################################
# ARIMA(0,1,0)
arima010 <- arima(nonseas.is,
                  order = c(0, 1, 0))

arima010
#coeftest(arima010)

# test LR
# H0 alfa1 = phi1 = 0
teststat <- 2*(as.numeric(logLik(arima111))- as.numeric(logLik(arima010)))
teststat
pchisq(teststat, df=1, lower.tail = FALSE) # odrzucamy H0

# badanie reszt
ACF <- acf(resid(arima010),
           lag.max = 24,
           na.action = na.pass,
           plot=FALSE)

ACF$lag <- ACF$lag * 12

plot(ACF,
     xlim=c(1,24),
     ylim=c(-0.5,0.5))

PACF <- pacf(resid(arima010),
             na.action=na.pass,
             plot=FALSE)

PACF$lag <- PACF$lag * 12

plot(PACF,
     xlim=c(1,24),
     ylim=c(-0.5,0.5),
     ylab='PACF')

Box.test(resid(arima010), type = "Ljung-Box", lag = 24)
Box.test(resid(arima010), type = "Box-Pierce", lag = 24)

#######################################

AIC(arima111, arima011, arima110)
BIC(arima111, arima011, arima110)

Box.test(resid(arima111), type = "Ljung-Box", lag = 24)
Box.test(resid(arima011), type = "Ljung-Box", lag = 24)
Box.test(resid(arima110), type = "Ljung-Box", lag = 24)

#######################################


arima.best.AIC <- auto.arima(nonseas.is,
                             d = 1,             # parameter d w modelu ARIMA
                             max.p = 4,         # p = maksymalna wartosc
                             max.q = 4,         # q = maksymalna wartosc
                             max.order = 8,     # suma p+q
                             start.p = 0,       # Wartosc startowa dla p
                             start.q = 0,       # Wartosc startowa dla q
                             ic = "aic",        # Wybor modelu na podstawie kryterium informcyjne
                             stepwise = FALSE,  # jezeli FALSE rozwaza wszystkie modeli
                             allowdrift = TRUE, # model zawiera stalą
                             trace = TRUE)      # wyswietlenie rozwazonych modeli

coeftest(arima.best.AIC)

arima110 <- arima(as.numeric(nonseas.is), # zmienna zależna
                    order = c(1, 1, 0),  # rzędy (p,d,q)
)

forecast110  =  predict(arima110,n.ahead = 3) 

nonseas.prog=ts_niesezonowe[276:278]

plot(as.numeric(ts_niesezonowe), xlim=c(250,279),xlab='Numer obserwacji',
     ylab='Cena kawy (U.S. cents/pound)', type='l')
abline(v = 276, lty = 2, col = "gray")
lines(forecast110$pred, col = "red", lwd = 2)
lines(forecast110$pred + 2 * forecast110$se, col = "red", lty = 3)
lines(forecast110$pred - 2 * forecast110$se, col = "red", lty = 3)

source("funs/functions.R")

errors(1,1,1,TRUE,nonseas.prog, ts_niesezonowe, 3)

### Modele ekstrapolacyjne

holt <- HoltWinters(nonseas.is,gamma=FALSE)
holt.pred <- predict(holt,n.ahead = 3,prediction.interval = TRUE)

plot(holt,xlab='Data',
     ylab='Cena kawy (U.S. cents/pound)')
holt

arima011 <- arima(nonseas.is, # zmienna zależna
                  order = c(0, 1, 1),  # rzędy (p,d,q)
)

forecast011  =  predict(arima011,n.ahead = 3) 




plot(nonseas.is, xlim=c(2020.5, 2023.4),xlab='Data',
     ylab='Cena kawy (U.S. cents/pound)', type='l')
abline(v = 2022.9167, lty = 2, col = "gray")
#lines(forecast011$pred, col = "red", lwd = 2)
#lines(forecast011$pred + 2 * forecast110$se, col = "red", lty = 3)
#lines(forecast011$pred - 2 * forecast110$se, col = "red", lty = 3)
lines(holt.pred[, 1], col = "blue")
lines(holt.pred[, 2], col = "blue", lty = 2)
lines(holt.pred[, 3], col = "blue", lty = 2)
legend(2020.5, 120, legend=c("ARIMA(0,1,1)", "Holt"),
       col=c("red", "blue"), lty=1:1, cex=0.4, pt.cex=1)

holt.pred
errors(1,1,1,FALSE,holt.pred[,1], ts_niesezonowe, 3)

cpi.holt<-holt(ts(nonseas.is), h = 3)

cpi.holt
# najlepszy model - ARIMA(0,1,1)
arima011 <- arima(as.numeric(nonseas.is), # zmienna zależna
                  order = c(0, 1, 1),  # rzędy (p,d,q)
)

forecast011 <- predict(arima011, n.ahead= 3)
forecast110 <- predict(arima110, n.ahead= 3)

summary(cpi.holt)

errors(1,1,1,FALSE,holt.pred[,1], ts_niesezonowe, 3)

cpi.holt<-holt(as.numeric(nonseas.is), h = 3)

plot(ts(ts_niesezonowe.xts), xlim=c(250,279),xlab='Obserwacja',
     ylab='Cena kawy (U.S. cents/pound)')
abline(v = 276, lty = 2, col = "gray")
lines(forecast011$pred, col = "red", lwd = 2)
lines(cpi.holt$mean, col = "blue", lwd = 2)
legend(250, 120, legend=c("ARIMA(0,1,1)", "Holt"),
       col=c("red", "blue"), lty=1:1, cex=0.5)

#######################
### SZEREG SEZONOWY ###
#######################

sezonowe$Value <- sezonowe$Value/1000000

ts_sezonowe <- ts(data=sezonowe$Value, start = c(2000, 1), end=c(2022, 12), frequency=12)
ts_sezonowe.xts = as.xts(ts_sezonowe)

sezonowy_da = decompose(ts_sezonowe, 'additive')
plot(sezonowy_da)

sezonowy_dm = decompose(ts_sezonowe, 'multiplicative')
plot(sezonowy_dm)

forecast::ggsubseriesplot(ts_sezonowe, ylab = "")


ts_sezonowe_full = ts_sezonowe
ts_sezonowe_oos <- ts_sezonowe[265:276]
ts_sezonowe <- ts_sezonowe[1:264]

ts_sezonowe
ts_sezonowe_oos

dsez = diff.xts(ts_sezonowe, lag=1)
plot(dsez, type='l')

d12dsez = diff.xts(dsez, lag=12)


ACF <- acf(ts_sezonowe,
           lag.max = 36,
           xlim=c(2,36),
           na.action = na.pass,
           plot=TRUE)

# ACF$lag <- ACF$lag * 12

#plot(ACF,
#     xlim=c(2,36),
#     ylim=c(-0.9,0.9))

PACF <- pacf(d12dsez,
             na.action=na.pass,
             xlim=c(2,36),
             lag.max = 36,
             plot=TRUE)

#PACF$lag <- PACF$lag * 12

#plot(PACF,
#     xlim=c(1,36),
#     ylim=c(-0.9,0.9),
#     ylab='PACF')

# Test Dickey-Haszy-Fullera

d12sez = diff.xts(ts_sezonowe, lag=12)
lag12sez = lag.xts(ts_sezonowe, k=12) 

model1 = lm(d12sez~0+lag12sez)
summary(model1)
bg1 <- bgtest(model1, order=1)
bg1

# test ADHF
lagd12sez = lag.xts(d12sez, k=1)

model2=lm(d12sez~0+lag12sez+lagd12sez)
summary(model2)

bg1 <- bgtest(model2, order=1)
bg1

lag2d12sez = lag.xts(d12sez, k=2)

model3 <- lm(d12sez~0+lag12sez+lagd12sez+lag2d12sez)
summary(model3)

bg1 <- bgtest(model3, order=1)
bg1
bg2 <- bgtest(model3, order=2)
bg2

lag3d12sez = lag.xts(d12sez, k=3)

model4 <- lm(d12sez~0+lag12sez+lagd12sez+lag2d12sez+lag3d12sez)
summary(model4)

bg1 <- bgtest(model4, order=1)
bg1
bg2 <- bgtest(model4, order=2)
bg2
bg3 <- bgtest(model4, order=3)
bg3
bg4 <- bgtest(model4, order=4)
bg4
bg5 <- bgtest(model4, order=5)
bg5

# stat. test: -0.840
# stat krytyczna: -5.83
# brak podstaw do odrzucenia H0
d12sez = diff.xts(ts_sezonowe, lag=12)
plot(d12dsez, type='l')
testdf(d12dsez, 'nc', 3, 4)
kpss.test <- ur.kpss(d12dsez, type=c('mu'))
summary(kpss.test)

nobs = length(ts_sezonowe)

#### SARIMA estymacja czesci sezonowej

arima010311 <- arima(ts_sezonowe,
                     # rzędy (p,d,q)
                     order = c(0, 1, 0),
                     # rzędy sezonowe (P,D,Q)
                     seasonal = list(order = c(3, 1, 1),
                                     # częstotliwość danych (12 dla danych miesięcznych)
                                     period = 12)       # dodatkowe regresory - stala
)
arima010211 <- arima(ts_sezonowe,
                     # rzędy (p,d,q)
                     order = c(0, 1, 0),
                     # rzędy sezonowe (P,D,Q)
                     seasonal = list(order = c(2, 1, 1),
                                     # częstotliwość danych (12 dla danych miesięcznych)
                                     period = 12)       # dodatkowe regresory - stala
)
teststat<- 2*(as.numeric(logLik(arima010311))-as.numeric(logLik(arima010211)))
teststat

pchisq(teststat, df=1, lower.tail = FALSE )
# brak podstaw do odrzucenia H0

arima010111 <- arima(ts_sezonowe,
                     # rzędy (p,d,q)
                     order = c(0, 1, 0),
                     # rzędy sezonowe (P,D,Q)
                     seasonal = list(order = c(1, 1, 1),
                                     # częstotliwość danych (12 dla danych miesięcznych)
                                     period = 12)       # dodatkowe regresory - stala
)
teststat<- 2*(as.numeric(logLik(arima010311))-as.numeric(logLik(arima010111)))
teststat

pchisq(teststat, df=1, lower.tail = FALSE )
# brak podstaw do odrzucenia H0

arima010011 <- arima(ts_sezonowe,
                     # rzędy (p,d,q)
                     order = c(0, 1, 0),
                     # rzędy sezonowe (P,D,Q)
                     seasonal = list(order = c(0, 1, 1),
                                     # częstotliwość danych (12 dla danych miesięcznych)
                                     period = 12)       # dodatkowe regresory - stala
)
teststat<- 2*(as.numeric(logLik(arima010311))-as.numeric(logLik(arima010011)))
teststat

pchisq(teststat, df=1, lower.tail = FALSE )
# brak podstaw do odrzucenia H0


arima010010 <- arima(ts_sezonowe,
                     # rzędy (p,d,q)
                     order = c(0, 1, 0),
                     # rzędy sezonowe (P,D,Q)
                     seasonal = list(order = c(0, 1, 0),
                                     # częstotliwość danych (12 dla danych miesięcznych)
                                     period = 12)       # dodatkowe regresory - stala
)
teststat<- 2*(as.numeric(logLik(arima010311))-as.numeric(logLik(arima010010)))
teststat

pchisq(teststat, df=1, lower.tail = FALSE )
# odrzucamy H0

sarima_test
coeftest(sarima_test)


AIC(arima010311, arima010211, arima010111, arima010011, arima010010)
#czesc sezonowa SARIMA(0,1,0)(0,1,1) - AIC
# wartosci BIC
BIC(arima010311, arima010211, arima010111, arima010011, arima010010) 


# SARIMA  (3,1,1) - SAR(3) nieist.
# SARIMA  (2,1,1) - SAR(2) nieist
# SARIMA  (1,1,1) - SAR(1) nieist

#### SARIMA estymacja czesci regularnej

# 2,1,2
arima212011 <- arima(ts_sezonowe,
                     # rzędy (p,d,q)
                     order = c(2, 1, 2),
                     # rzędy sezonowe (P,D,Q)
                     seasonal = list(order = c(0, 1, 1),
                                     # częstotliwość danych (12 dla danych miesięcznych)
                                     period = 12)
)

arima212011
coeftest(arima212011)


par(mfrow = c(2, 1))
acf(resid(arima212011), lag.max = 36,
    ylim = c(-0.4, 0.4), xlim=c(2,36), lwd = 4, col = "red")
pacf(resid(arima212011), lag.max = 36,
     lwd = 4, col = "red")
par(mfrow=c(1, 1))

Box.test(resid(arima212011), type = "Ljung-Box", lag = 36)
Box.test(resid(arima212011), type = "Box-Pierce", lag = 36)

# 1,1,2 OK
arima112011 <- arima(ts_sezonowe,
                     # rzędy (p,d,q)
                     order = c(1, 1, 2),
                     # rzędy sezonowe (P,D,Q)
                     seasonal = list(order = c(0, 1, 1),
                                     # częstotliwość danych (12 dla danych miesięcznych)
                                     period = 12)
)

arima112011 
coeftest(arima112011)


par(mfrow = c(2, 1))
acf(resid(arima112011), lag.max = 36,
    ylim = c(-0.4, 0.4), xlim=c(2,36), lwd = 4, col = "red")
pacf(resid(arima112011), lag.max = 36,
     lwd = 4, col = "red")
par(mfrow=c(1, 1))

Box.test(resid(arima112011), type = "Ljung-Box", lag = 36)
Box.test(resid(arima112011), type = "Box-Pierce", lag = 36)

teststat<- 2*(as.numeric(logLik(arima212011))-as.numeric(logLik(arima112011)))
teststat

pchisq(teststat, df=1, lower.tail = FALSE )
# brak podstaw do odrzucenia H0

# 0,1,2 OK
arima012011 <- arima(ts_sezonowe,
                     # rzędy (p,d,q)
                     order = c(0, 1, 2),
                     # rzędy sezonowe (P,D,Q)
                     seasonal = list(order = c(0, 1, 1),
                                     # częstotliwość danych (12 dla danych miesięcznych)
                                     period = 12)
)

arima012011 
coeftest(arima012011)


par(mfrow = c(2, 1))
acf(resid(arima012011), lag.max = 36,
    ylim = c(-0.4, 0.4), xlim=c(2,36), lwd = 4, col = "red")
pacf(resid(arima012011), lag.max = 36,
     lwd = 4, col = "red")
par(mfrow=c(1, 1))

Box.test(resid(arima012011), type = "Ljung-Box", lag = 36)
Box.test(resid(arima012011), type = "Box-Pierce", lag = 36)

teststat<- 2*(as.numeric(logLik(arima212011))-as.numeric(logLik(arima012011)))
teststat

pchisq(teststat, df=1, lower.tail = FALSE )
# brak podstaw do odrzucenia H0

# 0,1,1 MEH OK
arima011011 <- arima(ts_sezonowe,
                     # rzędy (p,d,q)
                     order = c(0, 1, 1),
                     # rzędy sezonowe (P,D,Q)
                     seasonal = list(order = c(0, 1, 1),
                                     # częstotliwość danych (12 dla danych miesięcznych)
                                     period = 12)
)

arima011011 
coeftest(arima011011)


par(mfrow = c(2, 1))
acf(resid(arima011011), lag.max = 36,
    ylim = c(-0.4, 0.4), xlim=c(2,36), lwd = 4, col = "red")
pacf(resid(arima011011), lag.max = 36,
     lwd = 4, col = "red")
par(mfrow=c(1, 1))

Box.test(resid(arima011011), type = "Ljung-Box", lag = 36)
Box.test(resid(arima011011), type = "Box-Pierce", lag = 36)

teststat<- 2*(as.numeric(logLik(arima212011))-as.numeric(logLik(arima011011)))
teststat

pchisq(teststat, df=1, lower.tail = FALSE )
# brak podstaw do odrzucenia H0

# 0,1,0 NAH
arima010011 <- arima(ts_sezonowe,
                     # rzędy (p,d,q)
                     order = c(0, 1, 0),
                     # rzędy sezonowe (P,D,Q)
                     seasonal = list(order = c(0, 1, 1),
                                     # częstotliwość danych (12 dla danych miesięcznych)
                                     period = 12)
)

arima010011 
coeftest(arima010011)


par(mfrow = c(2, 1))
acf(resid(arima010011), lag.max = 36,
    ylim = c(-0.4, 0.4), xlim=c(2,36), lwd = 4, col = "red")
pacf(resid(arima010011), lag.max = 36,
     lwd = 4, col = "red")
par(mfrow=c(1, 1))

Box.test(resid(arima010011), type = "Ljung-Box", lag = 36)
Box.test(resid(arima010011), type = "Box-Pierce", lag = 36)

teststat<- 2*(as.numeric(logLik(arima212011))-as.numeric(logLik(arima010011)))
teststat

pchisq(teststat, df=1, lower.tail = FALSE )
# odrzucamy H0

# 2,1,1 MEH OK
arima211011 <- arima(ts_sezonowe,
                     # rzędy (p,d,q)
                     order = c(2, 1, 1),
                     # rzędy sezonowe (P,D,Q)
                     seasonal = list(order = c(0, 1, 1),
                                     # częstotliwość danych (12 dla danych miesięcznych)
                                     period = 12)
)

arima211011 
coeftest(arima211011)


par(mfrow = c(2, 1))
acf(resid(arima211011), lag.max = 36,
    ylim = c(-0.4, 0.4), xlim=c(2,36), lwd = 4, col = "red")
pacf(resid(arima211011), lag.max = 36,
     lwd = 4, col = "red")
par(mfrow=c(1, 1))

Box.test(resid(arima211011), type = "Ljung-Box", lag = 36)
Box.test(resid(arima211011), type = "Box-Pierce", lag = 36)

teststat<- 2*(as.numeric(logLik(arima212011))-as.numeric(logLik(arima211011)))
teststat

pchisq(teststat, df=1, lower.tail = FALSE )
# brak podstaw do odrzucenia H0

# 2,1,0 MEH OK
arima210011 <- arima(ts_sezonowe,
                     # rzędy (p,d,q)
                     order = c(2, 1, 0),
                     # rzędy sezonowe (P,D,Q)
                     seasonal = list(order = c(0, 1, 1),
                                     # częstotliwość danych (12 dla danych miesięcznych)
                                     period = 12)
)

arima210011 
coeftest(arima210011)


par(mfrow = c(2, 1))
acf(resid(arima210011), lag.max = 36,
    ylim = c(-0.4, 0.4), xlim=c(2,36), lwd = 4, col = "red")
pacf(resid(arima210011), lag.max = 36,
     lwd = 4, col = "red")
par(mfrow=c(1, 1))

Box.test(resid(arima210011), type = "Ljung-Box", lag = 36)
Box.test(resid(arima210011), type = "Box-Pierce", lag = 36)

teststat<- 2*(as.numeric(logLik(arima212011))-as.numeric(logLik(arima210011)))
teststat

pchisq(teststat, df=1, lower.tail = FALSE )
# brak podstaw do odrzucenia H0

# 1,1,1 MEH OK nieist
arima111011 <- arima(ts_sezonowe,
                     # rzędy (p,d,q)
                     order = c(1, 1, 1),
                     # rzędy sezonowe (P,D,Q)
                     seasonal = list(order = c(0, 1, 1),
                                     # częstotliwość danych (12 dla danych miesięcznych)
                                     period = 12)
)

arima111011 
coeftest(arima111011)


par(mfrow = c(2, 1))
acf(resid(arima111011), lag.max = 36,
    ylim = c(-0.4, 0.4), xlim=c(2,36), lwd = 4, col = "red")
pacf(resid(arima111011), lag.max = 36,
     lwd = 4, col = "red")
par(mfrow=c(1, 1))

Box.test(resid(arima111011), type = "Ljung-Box", lag = 36)
Box.test(resid(arima111011), type = "Box-Pierce", lag = 36)

teststat<- 2*(as.numeric(logLik(arima212011))-as.numeric(logLik(arima111011)))
teststat

pchisq(teststat, df=1, lower.tail = FALSE )
# brak podstaw do odrzucenia H0

# 1,1,0 MEH OK 
arima110011 <- arima(ts_sezonowe,
                     # rzędy (p,d,q)
                     order = c(1, 1, 0),
                     # rzędy sezonowe (P,D,Q)
                     seasonal = list(order = c(0, 1, 1),
                                     # częstotliwość danych (12 dla danych miesięcznych)
                                     period = 12)
)

arima110011 
coeftest(arima110011)


par(mfrow = c(2, 1))
acf(resid(arima110011), lag.max = 36,
    ylim = c(-0.4, 0.4), xlim=c(2,36), lwd = 4, col = "red")
pacf(resid(arima110011), lag.max = 36,
     lwd = 4, col = "red")
par(mfrow=c(1, 1))

Box.test(resid(arima110011), type = "Ljung-Box", lag = 36)
Box.test(resid(arima110011), type = "Box-Pierce", lag = 36)

teststat<- 2*(as.numeric(logLik(arima212011))-as.numeric(logLik(arima110011)))
teststat

pchisq(teststat, df=1, lower.tail = FALSE )
# brak podstaw do odrzucenia H0


# 2,1,2
# 1,1,2 OK
# dw 0,1,2 OK - MA(2) nieist
# 0,1,1 MEH OK p-value < 10%
# dw 0,1,0 NAH
# dw 2,1,1 MEH OK AR(2,1) MA(1) nieist
# dw 2,1,0 MEH OK AR(2) nieist
# dw 1,1,1 MEH OK nieist AR(1) MA(1) nieist
# 1,1,0 MEH OK

AIC(arima212011, arima112011, arima011011, arima110011)
BIC(arima212011, arima112011, arima011011, arima110011) 
Box.test(resid(arima212011), type = "Ljung-Box", lag = 36)
Box.test(resid(arima112011), type = "Ljung-Box", lag = 36)
Box.test(resid(arima011011), type = "Ljung-Box", lag = 36)
Box.test(resid(arima110011), type = "Ljung-Box", lag = 36)

### PROGNOZOWANIE

# 212
forecast110 <- predict(arima110011, n.ahead = 12)

forecast110
str(forecast110)

ts.plot(as.numeric(ts_sezonowe_full),
        xlim = c(200,280),
        ylab = "Liczba pasażerów kolei w USA (w milionach)",
        xlab='Numer obserwacji')
# pocztek okresu prognozy
abline(v = 132, lty = 2, col = "gray")
lines(forecast110$pred, col = "red", lwd = 2)
abline(v = 265, lty = 2, col = "gray")
lines(forecast110$pred + 2 * forecast110$se, col = "red", lty = 3)
lines(forecast110$pred - 2 * forecast110$se, col = "red", lty = 3)


source("funs/functions.R")
list = list(arima212011,arima112011,arima011011,arima110011)
errors_sarima(list, ts_sezonowe_full, 12)

HWadd <- HoltWinters(ts(ts_sezonowe, frequency = 12, start = c(2000, 1)),seasonal = "additive")
HWadd.forecast <- predict(HWadd,n.ahead = 12,prediction.interval = TRUE)
plot(HWadd, lwd=2)
HWadd
ts(ts_sezonowe, frequency = 12, start = c(2000, 1))

plot(ts_sezonowe_full,
     xlim=c(2018, 2023.1),
     ylab = "Liczba pasażerów kolei w USA (w milionach)",
     xlab='Data')
lines(HWadd.forecast[, 1], col = "blue") # prognoza
lines(HWadd.forecast[, 2], col = "blue", lty = 2) # dolna granica przedziału ufności
lines(HWadd.forecast[, 3], col = "blue", lty = 2) # górna granica przedziału ufności
abline(v = c(2022, 01), lty=2)

errors(0,0,0, FALSE, HWadd.forecast[, 1], ts_sezonowe_full, 12)

arima212011graph <- arima(ts(ts_sezonowe, frequency = 12, start = c(2000, 1)),
                                         # rzędy (p,d,q)
                                         order = c(2, 1, 2),
                                         # rzędy sezonowe (P,D,Q)
                                         seasonal = list(order = c(0, 1, 1),
                                                         # częstotliwość danych (12 dla danych miesięcznych)
                                                         period = 12)
)
forecast212 <- predict(arima212011graph, n.ahead = 12)

plot(ts_sezonowe_full, xlim=c(2018, 2023.1),
     ylab = "Liczba pasażerów kolei w USA (w milionach)",
     xlab='Data')
abline(v = c(2022,01), lty = 2, col = "gray")
lines(forecast212$pred, col = "red", lwd = 2)
lines(forecast212$pred + 2 * forecast212$se, col = "red", lty = 3)
lines(forecast212$pred - 2 * forecast212$se, col = "red", lty = 3)
lines(HWadd.forecast[, 1], col = "blue")
lines(HWadd.forecast[, 2], col = "blue", lty = 2)
lines(HWadd.forecast[, 3], col = "blue", lty = 2)
legend(2018, 680, legend=c("SARIMA(2,1,2)(0,1,1)", "Holt-Winters"),
       col=c("red", "blue"), lty=1:1, cex=0.4, pt.cex=1)

