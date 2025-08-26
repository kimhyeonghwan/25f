
source("C:/Users/haoxing/Desktop/Work/circuit_breaker/Replication package/Empirical Analysis/Initialize.R")

setwd("C:/Users/haoxing/Desktop/Work/circuit_breaker/Replication package/Empirical Analysis/data")

# Load Data
Y <- 2013
for (m in seq(5,12)){
  Name <- paste("DataReg", Y,m, "New.rds", sep =  "")
  if (m == 5) {
    Data <- read_rds(Name)
    print(dim(Data)[1])
  } else {
    Data_Temp <-  read_rds(Name)
    print(dim(Data_Temp)[1])
    Data <- rbind(Data, Data_Temp)
  }
  print(paste("Year: ", Y))
  print(paste("Month: ", m))
}

for (Y  in (2014:2020)) {
  for (m in seq(1,12)){
    Name <- paste("DataReg", Y,m, "New.rds", sep =  "")
    Data_Temp <-  read_rds(Name)
    print(dim(Data_Temp)[1])
    Data <- rbind(Data, Data_Temp)
    print(paste("Year: ", Y))
    print(paste("Month: ", m))
  }
}

rm(Data_Temp)

### Prepare variables

NormFactor <- 60 * 6.5

# Make daily
Data <- Data %>%
  mutate(TV = sqrt(NormFactor) * TV,
         AdjTV_5S = sqrt(NormFactor) * AdjTV_5S,
         AdjTV_10S = sqrt(NormFactor) * AdjTV_10S,
         Ret = NormFactor * Ret,
  )




# Remove contract change day and the previous day (uncomment for volume data only)

# Data <- Data %>% filter(!(Day %in% c(14,15) & Month == 3 & Year == 2013), !(Day %in% c(20,21) & Month == 6 & Year == 2013), !(Day %in% c(19,20) & Month %in% c(9,12) & Year == 2013),
#                        !(Day %in% c(20,21) & Month == 3 & Year == 2014), !(Day %in% c(19,20) & Month == 6 & Year == 2014), !(Day %in% c(18,19) & Month %in% c(9,12) & Year == 2014),
#                        !(Day %in% c(19,20) & Month == 3 & Year == 2015), !(Day %in% c(18,19) & Month == 6 & Year == 2015), !(Day %in% c(17,18) & Month %in% c(9,12) & Year == 2015),
#                        !(Day %in% c(17,18) & Month == 3 & Year == 2016), !(Day %in% c(16,17) & Month == 6 & Year == 2016), !(Day %in% c(15,16) & Month %in% c(9,12) & Year == 2016),
#                        !(Day %in% c(16,17) & Month == 3 & Year == 2017), !(Day %in% c(15,16) & Month == 6 & Year == 2017), !(Day %in% c(14,15) & Month %in% c(9,12) & Year == 2017),
#                        !(Day %in% c(15,16) & Month == 3 & Year == 2018), !(Day %in% c(14,15) & Month == 6 & Year == 2018), !(Day %in% c(20,21) & Month %in% c(9,12) & Year == 2018),
#                        !(Day %in% c(14,15) & Month == 3 & Year == 2019), !(Day %in% c(20,21) & Month == 6 & Year == 2019), !(Day %in% c(19,20) & Month %in% c(9,12) & Year == 2019),
#                        !(Day %in% c(19,20) & Month == 3 & Year == 2020), !(Day %in% c(18,19) & Month == 6 & Year == 2020), !(Day %in% c(17,18) & Month %in% c(9,12) & Year == 2020), )

# uncomment to remove 4 circuit breaker triggerred days

# Data <- Data %>% filter(!(Day %in% c(9, 12, 16, 18) & Month == 3 & Year == 2020)) 


# filter out data with Emini price limit is triggered

Data <- Data %>%
  filter(PDiff >0)


# Make date FE
Data <- Data %>%
  mutate(YY = factor(paste(Year,Month,Day)))

# Create rolling averages of prices and volumes
Data <- Data %>% 
  arrange(Year, Month, Day, TimeMin) %>%
  mutate(Price1Day = roll_mean(PriceM, 60 * 6.5),
         Price1Month = roll_mean(PriceM, 60 * 6.5 * 21),
         Price1Hour = roll_mean(PriceM, 60),
         Vol1Day = roll_mean(VolumeM, 60 * 6.5),
         Vol1Month = roll_mean(VolumeM, 60 * 6.5 * 21)
  )


# Create monthly average, excluding today, for denom of PLev
DataDen <- Data %>%
  select(Year, Month, Day, TimeMin, Price1Month, Price1Day, Vol1Month, Vol1Day) %>%
  arrange(Year, Month, Day, TimeMin) %>%
  group_by(Year, Month, Day) %>%
  slice(n()) %>%
  ungroup() %>%
  mutate(PriceDen = dplyr::lag(Price1Month),
         VolDen = dplyr::lag(Vol1Month),
         VolDenDay = dplyr::lag(Vol1Day)) %>%
  select(Year, Month, Day, PriceDen, VolDen,  VolDenDay)

Data <- Data %>%
  left_join(DataDen)


Data <- Data %>%
  mutate(PLevMonthly = Price1Hour/PriceDen) %>%
  mutate(PLev_Square = (PLevMonthly-1)^2)


# Create lags, leads, and normalized values
for (Lag in seq(1, 45, 1)) {
  Name <- paste("AdjTV_5S_Lag", Lag, sep = "")
  Data <- Data %>%
    mutate(!!sym(Name) := dplyr::lag(AdjTV_5S,Lag))
} 

for (Lag in seq(1, 19, 1)) {
  Name <- paste("TV_Lag", Lag, sep = "")
  Data <- Data %>%
    mutate(!!sym(Name) := dplyr::lag(TV,Lag))
}

for (Lag in seq(1, 45, 1)) {
  Name <- paste("AdjTV_10S_Lag", Lag, sep = "")
  Data <- Data %>%
    mutate(!!sym(Name) := dplyr::lag(AdjTV_10S,Lag))
}

for (Lag in seq(1, 10, 1)) {
  Name <- paste("Skew_Lag", Lag, sep = "")
  Data <- Data %>%
    mutate(!!sym(Name) := dplyr::lag(Skew_ACJV,Lag))
}

for (Lag in seq(1, 10, 1)) {
  Name <- paste("Ret_Lag", Lag, sep = "")
  Data <- Data %>%
    mutate(!!sym(Name) := dplyr::lag(Ret,Lag))
}

Data <- Data %>%
  mutate(Ret_Lead = dplyr::lead(Ret,1))       


Data <- Data %>%
  mutate(VolNormalizedDay = VolumeM / VolDenDay - 1) 

for (Lag in seq(1, 40, 1)) {
  Name <- paste("VolNormalizedDay_Lag", Lag, sep = "")
  Data <- Data %>%
    group_by(Year, Month, Day) %>%
    mutate(!!sym(Name) := dplyr::lag(VolNormalizedDay,Lag)) %>%
    ungroup()
}


# Create indicators
for (Dum in seq(0.01, 0.07, 0.01)) {
  Name <- paste("Pre", Dum, sep = "")
  Data <- Data %>%
    mutate(!!sym(Name) := (PDiff < Dum))
}

# filter out data after 3:25pm

Data <- Data %>%
  filter(Pre325 == 1)


### Summary statistics
Data_summary <- Data %>%
  select(YY, TimeMin, Year, Month, Day, TV, AdjTV_5S, AdjTV_10S, PDiff, Skew_ACJV, Ret, VolNormalizedDay, PLevMonthly, PLev_Square)

summary(Data_summary)

print(paste("TV std: ", sd(Data_summary$TV, na.rm= T)))
print(paste("AdjTV_5S std: ", sd(Data_summary$AdjTV_5S, na.rm= T)))
print(paste("AdjTV_10S std: ", sd(Data_summary$AdjTV_10S, na.rm= T)))
print(paste("PDiff std: ", sd(Data_summary$PDiff, na.rm= T)))
print(paste("Skew std: ", sd(Data_summary$Skew_ACJV, na.rm= T)))
print(paste("Return std: ", sd(Data_summary$Ret, na.rm= T)))
print(paste("VolNormalizedDay std: ", sd(Data_summary$VolNormalizedDay, na.rm= T)))
print(paste("PLevMonthly std: ", sd(Data_summary$PLevMonthly, na.rm= T)))

quantile(Data_summary$TV, probs = c(0.01, 0.99), na.rm=T)
quantile(Data_summary$AdjTV_5S, probs = c(0.01, 0.99), na.rm=T)
quantile(Data_summary$AdjTV_10S, probs = c(0.01, 0.99), na.rm=T)
quantile(Data_summary$PDiff, probs = c(0.01, 0.99), na.rm=T)
quantile(Data_summary$Skew_ACJV, probs = c(0.01, 0.99), na.rm=T)
quantile(Data_summary$Ret, probs = c(0.01, 0.99), na.rm=T)
quantile(Data_summary$VolNormalizedDay, probs = c(0.01, 0.99), na.rm=T)
quantile(Data_summary$PLevMonthly, probs = c(0.01, 0.99), na.rm=T)

write.csv(summary(Data_summary),"summary_stat.csv")

rm(Data_summary)

Data_cov <- Data %>%
  select(AdjTV_5S, Skew_ACJV, Ret_Lead, VolNormalizedDay, PDiff, PLevMonthly, PLev_Square) %>%
  drop_na()

Cor_var <- cor(Data_cov)

tidy_Cor_var <- tidy(Cor_var)
write.csv(tidy_Cor_var, "Cor_var.csv")

Data_cov <- Data %>%
  select(VolNormalizedDay, PDiff, PLevMonthly, PLev_Square) %>%
  drop_na()

Cor_var <- cor(Data_cov)

tidy_Cor_var <- tidy(Cor_var)
write.csv(tidy_Cor_var, "Cor_vol_var.csv")



### Left hand side variable analysis (Can jump to regressions directly)
Data_ARIMA <- Data %>%
  filter(!is.na(AdjTV_5S), !is.na(VolNormalizedDay), !is.na(Skew), !is.na(Ret), !is.na(Ret_Lead),
        !is.na(PLevMonthly), !is.na(TV), !is.na(AdjTV_10S))


AdjTV_res <- felm(AdjTV_5S ~ 1|YY + TimeFac|0|0, Data_ARIMA)
TV_res <- felm(TV ~ 1|YY + TimeFac|0|0, Data_ARIMA)
AdjTV_10S_res <- felm(AdjTV_10S ~ 1|YY + TimeFac|0|0, Data_ARIMA)
VolNormalizedDay_res <- felm(VolNormalizedDay ~ 1|YY + TimeFac|0|0, Data_ARIMA)
Skew_res <- felm(Skew_ACJV ~ 1|YY + TimeFac|0|0, Data_ARIMA)
Ret_Lead_res <- felm(Ret_Lead ~ 1|YY + TimeFac|0|0, Data_ARIMA)
PDiff_res <- felm(PDiff ~ 1|YY + TimeFac|0|0, Data_ARIMA)
Lev_res <- felm(PLevMonthly ~ 1|YY + TimeFac|0|0, Data_ARIMA)

Data_ARIMA <- Data_ARIMA %>%
  mutate(AdjTV_5S_res = AdjTV_res$residuals,
         AdjTV_10S_res = AdjTV_10S_res$residuals,
         TV_res = TV_res$residuals,
         Volume_res = VolNormalizedDay_res$residuals,
         Skew_res = Skew_res$residuals,
         Ret_Lead_res = Ret_Lead_res$residuals,
         PLevMonthly_res = Lev_res$residuals
  ) %>%
  filter(!is.na(AdjTV_5S_res), !is.na(Volume_res), !is.na(Skew_res), !is.na(Ret_Lead_res), !is.na(PLevMonthly_res))


# Volatility

AdjTV_ARIMA <- auto.arima(Data_ARIMA$AdjTV_5S_res, max.p=50, max.d=0, max.q=0)

summary(AdjTV_ARIMA)

TV_ARIMA <- auto.arima(Data_ARIMA$TV_res, max.p=50, max.d=0, max.q=0)

summary(TV_ARIMA)

AdjTV_10S_ARIMA <- auto.arima(Data_ARIMA$AdjTV_10S_res, max.p=50, max.d=0, max.q=0)

summary(AdjTV_10S_ARIMA)

# Volume

Volume_ARIMA <- auto.arima(Data_ARIMA$Volume_res, max.p=50, max.d=0, max.q=0)

summary(Volume_ARIMA)

# Skewness

Skew_ARIMA <- auto.arima(Data_ARIMA$Skew_res, max.p=50, max.d=0, max.q=0)

summary(Skew_ARIMA)

# Return

Ret_ARIMA <- auto.arima(Data_ARIMA$Ret_Lead_res, max.p=50, max.d=0, max.q=0)

summary(Ret_ARIMA)


### Bin-Scatter Plots

D2 <- Data %>%
  filter(PDiff >=0, !is.na(AdjTV_5S)) 


for (Bin in seq(0.03, 0.14, 0.01)) {
  Name <- paste("Bin", Bin, sep = "")
  D2 <- D2 %>%
    mutate(!!sym(Name) := ((PDiff < Bin) & (PDiff >= Bin-0.01)))
}

D2 <- D2 %>%
 mutate(Bin0.0075 := ((PDiff < 0.0075) & (PDiff >=0))) %>%
 mutate(Bin0.0125 := ((PDiff < 0.0125) & (PDiff >=0.0075))) %>%
 mutate(Bin0.02 := ((PDiff < 0.02) & (PDiff >=0.0125))) #%>%


Reg <- felm(AdjTV_5S ~ PLevMonthly+PLev_Square+Bin0.0075+Bin0.0125
            +Bin0.02+Bin0.03+Bin0.04+Bin0.05+Bin0.06+Bin0.07+Bin0.08+Bin0.09+Bin0.1+
              Bin0.11+Bin0.12+Bin0.13+Bin0.14
            |TimeFac+YY|0|0, D2)

summary(Reg)


R <- NeweyWest(Reg, 29)
R_diag <- diag(R)
R_bin <- R_diag[3:length(R_diag)]
Betas <- Reg$coefficients
XX <- eye(length(Betas))
Fit_all <- XX %*% Betas
Fit <- Fit_all[3:length(Betas)]
SDs <- sqrt(R_bin)
Upper <- Fit + 2 * SDs
Lower <- Fit - 2 * SDs


ggplot(tibble(Fit=Fit,Upper=Upper,Lower=Lower), aes(x = c(0.0075,0.0125,seq(0.02,0.14,0.01)), y = Fit)) + 
  geom_point(col='red') + 
  geom_ribbon(aes(ymin =Upper, ymax = Lower), alpha = 0.1) +
  xlab("Emini_PDiff") +
  ylab("TSTV")



### Regressions

# Volatility

Vol_Data <- Data %>%
  select(TimeMin,  YY, TimeFac,
         AdjTV_5S, PLevMonthly, PLev_Square, PDiff, Pre0.07, Pre0.06, Pre0.05, Pre0.04, Pre0.03, Pre0.02, Pre0.01, 
         AdjTV_5S_Lag1, AdjTV_5S_Lag2, AdjTV_5S_Lag3, AdjTV_5S_Lag4, AdjTV_5S_Lag5,
         AdjTV_5S_Lag6, AdjTV_5S_Lag7, AdjTV_5S_Lag8, AdjTV_5S_Lag9, AdjTV_5S_Lag10,
         AdjTV_5S_Lag11, AdjTV_5S_Lag12, 
         TV, TV_Lag1, TV_Lag2, TV_Lag3, TV_Lag4, TV_Lag5, TV_Lag6, TV_Lag7, TV_Lag8, TV_Lag9,
         TV_Lag10, TV_Lag11, TV_Lag12, TV_Lag13, TV_Lag14, TV_Lag15, TV_Lag16, TV_Lag17, TV_Lag18, TV_Lag19,
         AdjTV_10S, AdjTV_10S_Lag1, AdjTV_10S_Lag2, AdjTV_10S_Lag3, AdjTV_10S_Lag4, AdjTV_10S_Lag5,
         AdjTV_10S_Lag6, AdjTV_10S_Lag7, AdjTV_10S_Lag8, AdjTV_10S_Lag9, AdjTV_10S_Lag10,
         AdjTV_10S_Lag11, AdjTV_10S_Lag12, AdjTV_10S_Lag13, AdjTV_10S_Lag14, AdjTV_10S_Lag15,
         AdjTV_10S_Lag16, AdjTV_10S_Lag17, AdjTV_10S_Lag18, AdjTV_10S_Lag19, AdjTV_10S_Lag20,
         AdjTV_10S_Lag21, AdjTV_10S_Lag22, AdjTV_10S_Lag23, AdjTV_10S_Lag24, AdjTV_10S_Lag25, AdjTV_10S_Lag26,
         AdjTV_10S_Lag27, AdjTV_10S_Lag28, AdjTV_10S_Lag29, AdjTV_10S_Lag30, AdjTV_10S_Lag31, AdjTV_10S_Lag32,
         AdjTV_10S_Lag33, AdjTV_10S_Lag34, AdjTV_10S_Lag35, AdjTV_10S_Lag36, AdjTV_10S_Lag37, AdjTV_10S_Lag38, 
         AdjTV_10S_Lag39, AdjTV_10S_Lag40
  ) %>%
  drop_na() %>%
  mutate(PDiff_below7 = pmin(PDiff-0.07,0)) %>%
  mutate(PDiff_below6 = pmin(PDiff-0.06,0)) %>%
  mutate(PDiff_below5 = pmin(PDiff-0.05,0)) %>%
  mutate(PDiff_below4 = pmin(PDiff-0.04,0)) %>%
  mutate(PDiff_below3 = pmin(PDiff-0.03,0)) %>%
  mutate(PDiff_below2 = pmin(PDiff-0.02,0)) %>%
  mutate(PDiff_below1 = pmin(PDiff-0.01,0)) %>%
  mutate(PDiff_above7 = pmax(PDiff-0.07,0))


Reg1_Adj5S <- felm(AdjTV_5S ~ PLevMonthly + PLev_Square 
                   + PDiff
                   + PDiff_below2
                   + AdjTV_5S_Lag1 + AdjTV_5S_Lag2 + AdjTV_5S_Lag3
                   + AdjTV_5S_Lag4 + AdjTV_5S_Lag5 + AdjTV_5S_Lag6
                   + AdjTV_5S_Lag7 + AdjTV_5S_Lag8 + AdjTV_5S_Lag9
                   + AdjTV_5S_Lag10 + AdjTV_5S_Lag11 
                   |TimeFac+YY|0|0, Vol_Data)

summary(Reg1_Adj5S)

acf(Reg1_Adj5S$residuals)

auto.arima(Reg1_Adj5S$residuals, max.p=5, max.d=1, max.q=2)

acf(Reg1_Adj5S$residuals*Vol_Data$PDiff)

R1<-sqrt(diag(NeweyWest(Reg1_Adj5S, 29)))

R1

tidy_Reg1_Adj5S <- tidy(Reg1_Adj5S)
write.csv(tidy_Reg1_Adj5S, "Reg1_Adj5S.csv")

# Reg TV
Reg1_TV <- felm(TV ~ PLevMonthly + PLev_Square 
                + PDiff
                + PDiff_below2
                + TV_Lag1 + TV_Lag2 + TV_Lag3
                + TV_Lag4 + TV_Lag5 + TV_Lag6
                + TV_Lag7 + TV_Lag8 + TV_Lag9
                |YY + TimeFac|0|0, Vol_Data)

summary(Reg1_TV)

acf(Reg1_TV$residuals)

auto.arima(Reg1_TV$residuals, max.p=5, max.d=1, max.q=2)

acf(Reg1_TV$residuals*Vol_Data$SPX_PDiff)

R1_TV<-sqrt(diag(NeweyWest(Reg1_TV, 29)))

R1_TV

tidy_Reg1_TV <- tidy(Reg1_TV)
write.csv(tidy_Reg1_TV, "Reg1_TV.csv")

# Reg_AdjTV_10S

Reg1_Adj10S <- felm(AdjTV_10S ~ PLevMonthly + PLev_Square 
                    + PDiff
                    + PDiff_below2
                    + AdjTV_10S_Lag1 + AdjTV_10S_Lag2 + AdjTV_10S_Lag3
                    + AdjTV_10S_Lag4 + AdjTV_10S_Lag5 + AdjTV_10S_Lag6
                    + AdjTV_10S_Lag7 + AdjTV_10S_Lag8 + AdjTV_10S_Lag9
                    + AdjTV_10S_Lag10 + AdjTV_10S_Lag11 + AdjTV_10S_Lag12
                    + AdjTV_10S_Lag13 + AdjTV_10S_Lag14 + AdjTV_10S_Lag15
                    |TimeFac+YY|0|0, Vol_Data)

summary(Reg1_Adj10S)

acf(Reg1_Adj10S$residuals)

auto.arima(Reg1_Adj10S$residuals, max.p=5, max.d=1, max.q=2)

acf(Reg1_Adj10S$residuals*Vol_Data$SPX_PDiff)

R1_10S<-sqrt(diag(NeweyWest(Reg1_Adj10S, 33)))

R1_10S

tidy_Reg1_Adj10S <- tidy(Reg1_Adj10S)
write.csv(tidy_Reg1_Adj10S, "Reg1_Adj10S.csv")

# Volume

Volume_Data <- Data %>%
  select(TimeMin, TimeFac, YY,
         VolNormalizedDay, PDiff, PLevMonthly, PLev_Square, Pre0.07, Pre0.06, Pre0.05, Pre0.04, Pre0.03, Pre0.02, Pre0.01,
         VolNormalizedDay_Lag1, VolNormalizedDay_Lag2, VolNormalizedDay_Lag3,
         VolNormalizedDay_Lag4, VolNormalizedDay_Lag5, VolNormalizedDay_Lag6,
         VolNormalizedDay_Lag7, VolNormalizedDay_Lag8, VolNormalizedDay_Lag9,
         VolNormalizedDay_Lag10, VolNormalizedDay_Lag11, VolNormalizedDay_Lag12, 
         VolNormalizedDay_Lag13, VolNormalizedDay_Lag14, VolNormalizedDay_Lag15,
         VolNormalizedDay_Lag16, VolNormalizedDay_Lag17, VolNormalizedDay_Lag18,
         VolNormalizedDay_Lag19, VolNormalizedDay_Lag20, VolNormalizedDay_Lag21,
         VolNormalizedDay_Lag22, VolNormalizedDay_Lag23, VolNormalizedDay_Lag24,
         VolNormalizedDay_Lag25, VolNormalizedDay_Lag26, VolNormalizedDay_Lag27,
         VolNormalizedDay_Lag28, VolNormalizedDay_Lag29, VolNormalizedDay_Lag30,
         VolNormalizedDay_Lag31, VolNormalizedDay_Lag32, 
         VolNormalizedDay_Lag33
  ) %>%
  drop_na() %>%
  mutate(PDiff_above10 = pmax(PDiff-0.1,0)) %>%
  mutate(PDiff_below4 = pmin(PDiff-0.04,0)) %>%
  mutate(PDiff_below3 = pmin(PDiff-0.03,0)) %>%
  mutate(PDiff_below2 = pmin(PDiff-0.02,0)) %>%
  mutate(PDiff_below1 = pmin(PDiff-0.01,0))

#Volume_Data <- Volume_Data %>%
#  filter(PDiff <=0.07)

Reg2Day <- felm(VolNormalizedDay ~ PLevMonthly + PLev_Square
                + PDiff
                #+ PDiff_above10
                + PDiff_below2
                + VolNormalizedDay_Lag1 + VolNormalizedDay_Lag2 + VolNormalizedDay_Lag3
                + VolNormalizedDay_Lag4 + VolNormalizedDay_Lag5 + VolNormalizedDay_Lag6
                + VolNormalizedDay_Lag7 + VolNormalizedDay_Lag8 + VolNormalizedDay_Lag9
                + VolNormalizedDay_Lag10 + VolNormalizedDay_Lag11 + VolNormalizedDay_Lag12 
                + VolNormalizedDay_Lag13 + VolNormalizedDay_Lag14 + VolNormalizedDay_Lag15
                + VolNormalizedDay_Lag16 + VolNormalizedDay_Lag17 + VolNormalizedDay_Lag18
                + VolNormalizedDay_Lag19 + VolNormalizedDay_Lag20 
                + VolNormalizedDay_Lag21
                + VolNormalizedDay_Lag22 + VolNormalizedDay_Lag23 + VolNormalizedDay_Lag24
                + VolNormalizedDay_Lag25 + VolNormalizedDay_Lag26 + VolNormalizedDay_Lag27
                + VolNormalizedDay_Lag28 + VolNormalizedDay_Lag29 + VolNormalizedDay_Lag30
                + VolNormalizedDay_Lag31 + VolNormalizedDay_Lag32 
                + VolNormalizedDay_Lag33
                |YY+TimeFac|0|0, Volume_Data)

summary(Reg2Day)

acf(Reg2Day$residuals)

auto.arima(Reg2Day$residuals, max.p=5, max.d=1, max.q=2)

acf(Reg2Day$residuals*Volume_Data$PDiff)


R2 <-sqrt(diag(NeweyWest(Reg2Day, 29)))

R2

tidy_Reg2Day <- tidy(Reg2Day)
write.csv(tidy_Reg2Day, "Reg2Day.csv")


# Skewness

Skew_Data <- Data %>%
  select(TimeMin, YY, TimeFac,
         Skew_ACJV, PLevMonthly, PLev_Square, PDiff, Pre0.07, Pre0.06, Pre0.05, Pre0.04, Pre0.03, Pre0.02, Pre0.01,
         Skew_Lag1, Skew_Lag2, Skew_Lag3, Skew_Lag4, Skew_Lag5,
         Skew_Lag6, Skew_Lag7, Skew_Lag8, Skew_Lag9) %>%
  drop_na() %>%
  mutate(PDiff_below7=pmin(PDiff-0.07,0)) %>%
  mutate(PDiff_below5=pmin(PDiff-0.05,0)) %>%
  mutate(PDiff_below4=pmin(PDiff-0.04,0)) %>%
  mutate(PDiff_below3=pmin(PDiff-0.03,0)) %>%
  mutate(PDiff_below2=pmin(PDiff-0.02,0)) 


Reg3 <- felm(Skew_ACJV ~  PLevMonthly + PLev_Square
             + PDiff
             + PDiff_below2
             +Skew_Lag1 + Skew_Lag2 + Skew_Lag3 + Skew_Lag4 + Skew_Lag5
             +Skew_Lag6 + Skew_Lag7 
             |TimeFac+YY|0|0, Skew_Data)

summary(Reg3)

acf(Reg3$residuals)

auto.arima(Reg3$residuals, max.p=5, max.d=1, max.q=2)

acf(Reg3$residuals*Skew_Data$PDiff)

R3 <-sqrt(diag(NeweyWest(Reg3, 29)))

R3

tidy_Reg3 <- tidy(Reg3)
write.csv(tidy_Reg3, "Reg3.csv")

# Return

Return_Data <- Data %>%
  select(TimeMin,
         Ret_Lead, PLevMonthly, PLev_Square, PDiff, Pre0.07, Pre0.06, Pre0.05, Pre0.04, Pre0.03, Pre0.02, Pre0.01, YY, TimeFac,
         Ret, Ret_Lag1, Ret_Lag2, Ret_Lag3, Ret_Lag4, Ret_Lag5,
         Ret_Lag6 , Ret_Lag7, Ret_Lag8, Ret_Lag9
  ) %>%
  drop_na() %>%
  mutate(PDiff_above10=pmax(PDiff-0.1,0)) %>%
  mutate(PDiff_below7=pmin(PDiff-0.07,0)) %>%
  mutate(PDiff_below5=pmin(PDiff-0.05,0)) %>%
  mutate(PDiff_below4=pmin(PDiff-0.04,0)) %>%
  mutate(PDiff_below3=pmin(PDiff-0.03,0)) %>%
  mutate(PDiff_below2=pmin(PDiff-0.02,0)) %>%
  mutate(PDiff_below1=pmin(PDiff-0.01,0))


Reg4 <- felm(Ret_Lead ~ PLevMonthly + PLev_Square
             + PDiff 
               + PDiff_below2
               + Ret 
             + Ret_Lag1 + Ret_Lag2 + Ret_Lag3 
             + Ret_Lag4 + Ret_Lag5
             + Ret_Lag6 + Ret_Lag7 
             + Ret_Lag8 + Ret_Lag9
             |YY+TimeFac|0|0, Return_Data)

summary(Reg4)

pacf(Reg4$residuals)

auto.arima(Reg4$residuals, max.p=5, max.d=1, max.q=3)

acf(Reg4$residuals*Return_Data$PDiff)

R4 <-sqrt(diag(NeweyWest(Reg4, 29)))

R4

tidy_Reg4 <- tidy(Reg4)
write.csv(tidy_Reg4, "Reg4.csv")
