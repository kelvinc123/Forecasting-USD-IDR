rm(list=ls())

library(quantmod)
library(ggplot2)
library(forecast)
library(vars)
library(lmtest)
library(timeSeries)
library(rugarch)
library(tseries)

#Get necessary dataset from Yahoo Finance
usd_idr <- getSymbols("IDR=X",auto.assign = FALSE, from = as.Date("01/01/2002", format = "%m/%d/%Y"))[,6]
usd_idr <- to.monthly(usd_idr)[,4]

sp500 <- getSymbols("^GSPC",auto.assign = FALSE, from = as.Date("01/01/2002", format = "%m/%d/%Y"))[,6]
sp500 <- to.monthly(sp500)[,4]

#Convert Data to time series and calculate returns
usd_idr_ts <- ts(usd_idr, start = 2002, freq = 12)
sp500_ts <- ts(sp500, start = 2002, freq = 12)

usd_idr_ret=100*(diff(usd_idr_ts)/lag(usd_idr_ts, -1))
sp500_ret=100*(diff(sp500_ts)/lag(sp500_ts, -1))

#Plot Data
plot((sp500_ret-mean(sp500_ret))/sd(sp500_ret), main = "S&P500", ylab="Returns")
lines((usd_idr_ret-mean(usd_idr_ret))/sd(usd_idr_ret), col = "red")

#ACF and PACF plot
tsdisplay(usd_idr_ret, main = "USD/IDR Returns")
acf(usd_idr_ret, main="USD/IDR Returns")
pacf (usd_idr_ret, main="USD/IDR Returns")

tsdisplay(sp500_ret, main="S&P500 Returns")
acf(sp500_ret, main="S&P500 Returns")
pacf(sp500_ret, main="S&P500 Returns")



# Model -------------------------------------------------------------------

# specify the index
t = seq(2002 + (1/12), 2019 + (4/12), length = length(usd_idr_ret))

# periodic
mod_period_usdidr <- lm(usd_idr_ret ~ I(sin(2*pi*t)) + I(cos(2*pi*t)))
summary(mod_period)
plot(usd_idr_ret)
lines(mod_period$fit ~ t, col = "red")

# Seasonal
mod_seasonal_usdidr <- tslm(usd_idr_ret ~ season)
summary(mod_seasonal)
plot(usd_idr_ret)
lines(mod_seasonal$fit, col = "red")

# arma
mod_arma_22_usdidr <- arima(usd_idr_ret, order = c(2,0,2)) # based on the acf pacf
summary(mod_arma_22_usdidr)

plot(usd_idr_ret)
lines(mod_arma_22_usdidr$fit, col = "red")

summary(mod_arma_22_usdidr)

arima_summary <- function(model){
  coef_arma <- model$coef
  se_arma <- sqrt(diag(model$var.coef))
  t_arma <- coef_arma / se_arma
  p_arma <- pnorm(abs(t_arma), lower.tail = FALSE)
  cbind("coef" = coef_arma, "se" = se_arma, "t" = t_arma, "p" = round(p_arma, 5))
}

arima_summary(mod_arma_22_usdidr)


m = 4
AIC_mat <- matrix(NA, ncol = m, nrow = m)

for(i in 1:m){
  for(j in 1:m){
    AIC_mat[i,j] <- summary(arma(usd_idr_ret, order = c(i,j)))$aic
  }
}
AIC_mat

mod_arma_32_usdidr <- arima(usd_idr_ret, order = c(3,0,2)) # based on the for loop
arima_summary(mod_arma_32_usdidr) # but ar 3 is not stats sig under alpha = 0.05, so we use arma 22 instead

#auto arima
mod_autoarima_usdidr = auto.arima(usd_idr_ret)
arima_summary(mod_autoarima_usdidr)

#AIC BIC
AIC(mod_period_usdidr, mod_seasonal_usdidr, mod_arma_22_usdidr, mod_arma_32_usdidr, mod_autoarima_usdidr)
BIC(mod_period_usdidr, mod_seasonal_usdidr, mod_arma_22_usdidr, mod_arma_32_usdidr, mod_autoarima_usdidr)

#Choose mod_arma22
tsdisplay(mod_arma_22_usdidr$res)

#Analyze GARCH
tsdisplay(mod_arma_22_usdidr$res*mod_arma_22_usdidr$res)

#use for loop to find best GARCH(p,q)
m <- 8
garch_mse_mat <- matrix(NA, ncol = m, nrow = m)
garch_mae_mat <- matrix(NA, ncol = m, nrow = m)
for (i in 1:m){
  for(j in 1:m){
    model <- ugarchspec(
      variance.model = list(model = "sGARCH", garchOrder = c(i, j)),
      mean.model = list(armaOrder = c(2, 2), include.mean = TRUE),
      distribution.model = "sstd")
    
    modelfit <- ugarchfit(spec=model,data=usd_idr_ret)
    
    a <- attributes(modelfit)
    
    garch_mse_mat[i,j] <- MSE(usd_idr_ret, a$fit$fitted.values, digits = 5)
    garch_mae_mat[i,j] <- MAE(usd_idr_ret, a$fit$fitted.values, digits = 5)
    
  }
}

garch_mae_mat
garch_mse_mat

which(garch_mae_mat == min(garch_mae_mat), arr.ind = TRUE)
which(garch_mse_mat == min(garch_mse_mat), arr.ind = TRUE)

#GARCH(4,1) smalest error
garch_mod_11=ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(4, 1)),
  mean.model = list(armaOrder = c(2, 2), include.mean = TRUE),
  distribution.model = "sstd")

garchfit11 = ugarchfit(spec=garch_mod_11, data=usd_idr_ret)
garchfit11
plot(garchfit11)

att_garch_usdidr = attributes(garchfit11)
tsdisplay(att_garch_usdidr$fit$res)
att <- attributes(mod_arma_22_usdidr)

plot(usd_idr_ret)
lines(att_garch_usdidr$fit$fitted.values ~ t, col = "green")
lines(fitted(mod_arma_22_usdidr), col = "red")
plot(forecast(mod_arma_22_usdidr))
modelfor = ugarchforecast(garchfit11, data = NULL, n.ahead = 12, n.roll = 0, out.sample = 0)
plot(modelfor)

#Respective Residuals
tsplot(usdidr_mod$res ~ usdidr_mod$fit, main = "USD/IDR Model Residual",pch = 20, xlab = "Fitted values", ylab = "Residuals")
tsplot(sp500_mod$res ~ sp500_mod$fit, main = "S&P500 Model Residual",pch = 20, xlab = "Fitted values", ylab = "Residuals")

plot(usd_idr_ret, xlim = c(2002, 2021))
lines(ugarch_usdidr_mod, col = "red")
lines(predict(ugarch_usdidr_mod, n.ahead = 12)$pred, col = "red")

plot(usd_idr_ret, xlim = c(2002, 2021))
lines(usdidr_mod$fit, col = "red")
lines(predict(usdidr_mod, n.ahead = 12)$pred, col = "red")

plot(sp500_ret, xlim = c(2002, 2021))
lines(sp500_mod$fit, col = "green")
lines(predict(sp500_mod, n.ahead = 12)$pred, col = "green")

#ACF PACF Residuals
acf(usdidr_mod$residuals)
pacf(usdidr_mod$residuals)

acf(sp500_mod$residuals)
pacf(sp500_mod$residuals)

#CUSUM
plot(efp(usdidr_mod$res ~ 1, type = "Rec-CUSUM"))
plot(efp(sp500_mod$res ~ 1, type = "Rec-CUSUM"))

# recursive residual
rec_res_usdidr <- recresid(usdidr_mod$res ~ 1)
rec_res_sp500 <- recresid(sp500_mod$res ~ 1)

plot(rec_res_usdidr, main = "Recursive Residuals for USD/IDR", ylab = "Res", xlab = "Time")
plot(rec_res_sp500, main = "Recursive Residuals for S&P500", ylab = "Res", xlab = "Time")

#VAR
varbind = cbind(sp500_ret, usd_idr_ret)

#Choose Best model by AIC BIC
row = 20
aicbic = matrix(NA, ncol=3, nrow=row)
colnames(aicbic) = c("p", "AIC", "BIC")

for(i in 1:row){
  aicbic[i, 1] = i
  aicbic[i, 2] = AIC(VAR(varbind, p=i))
  aicbic[i, 3] = BIC(VAR(varbind, p=i))
}

aicbic
plot(aicbic[,1], aicbic[,2], xlab = "p", ylab="AIC", main="AIC plot")
plot(aicbic[,1], aicbic[,3], xlab = "p", ylab="BIC", main="BIC plot")

#Choose VAR(1)
varmod = VAR(varbind, p=1)
summary(varmod)

#Granger test 1 (H0 = USD/IDR does not explain SP500)
grangertest(sp500_ret ~ usd_idr_ret, order=1)

#Granger test 2 (H0 = SP500 does not explain USD/IDR)
grangertest(usd_idr_ret ~ sp500_ret, order=1)


