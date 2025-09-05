# running Global setting
setwd("/Users/hwan/KAIST/25fall/thesis")
source("1.setting.R")

data_min <- read_csv(gzfile("data/data_min2.csv.gz"))
data_close <- read_csv(gzfile("data/data_close2.csv.gz"))

# Thresholds base price
data_close <- data_close %>% 
  arrange(ISU_CD, TRD_DD) %>% 
  group_by(ISU_CD) %>% 
  mutate(P_close_lag = dplyr::lag(P_close)) %>% 
  drop_na(P_close_lag) %>% 
  ungroup() %>% 
  select(TRD_DD, ISU_CD, P_close_lag) %>% 
  filter()

# Convert to daily-scale & make thresholds
data_min <- data_min %>% 
  left_join(., data_close, by=c("TRD_DD", "ISU_CD")) %>% 
  drop_na(P_close_lag) %>% 
  select(-AvgRV_5s, -AvgRV_10s, -Skew) %>% 
  mutate(RV = 390 * RV,
         TSRV_5s = 390 * TSRV_5s,
         TSRV_10s = 390 * TSRV_10s,
         Ret_min = 390 * Ret_min) %>% 
  mutate(DTCB=if_else(T_min>=1450,
                      (P_low-0.8*P_close_lag)/P_close_lag,
                      if_else(CB_grp==0,
                              (P_low-0.92*P_close_lag)/P_close_lag,
                              (P_low-0.15*P_close_lag)/P_close_lag)))
  
# data_min %>% filter(DTCB<0) %>% view()

# Test
error_test <- data_min %>% group_by(TRD_DD,PROD_ID, ISU_CD) %>% summarise(cnt=1) %>% group_by(TRD_DD, PROD_ID) %>% summarise(cnt=sum(cnt))
sum(error_test$cnt>1) # zero!
# for 1day, 1 isu -> lead month base & second month ONLY expiry day
data_min %>% filter((P_low-P_close_lag)/P_close_lag>0.15) # NO thresholds level 2 (15%)

year_first_day <- data_close %>% 
  mutate(TRD_YY=substr(TRD_DD,1,4)) %>% 
  group_by(TRD_YY, TRD_DD) %>% 
  summarise() %>% 
  slice(1) %>% 
  filter(substr(TRD_DD,5,6)=='01')

# Create 1h & 21day Price, 1day(6.5h) Volume 
# 더미를 만들어 붙이고 롤섬해야함
Price_1h <- data_min %>% 
  group_by(PROD_ID) %>% 
  mutate(T_min_1h_lag=if_else(TRD_DD %in% year_first_day$TRD_DD,
                              if_else(T_min>=1100,T_min-100,T_min+430),
                              if_else(T_min>=1000,T_min-100,T_min+530))) %>% 
  mutate(TRD_DD_1h_lag=if_else(T_min_1h_lag>T_min,dplyr::lag(TRD_DD),TRD_DD))
  mutate(Avg_VWAP_1h_lag=roll_meanr(VWAP_min, 60))

Volume_1day <- data_min %>% 
  group_by(TRD_DD, PROD_ID) %>% 
  summarise(Avg_Volume_1day=sum(Volume)/390, .groups = "drop") %>% 
  arrange(PROD_ID, TRD_DD) %>% 
  group_by(PROD_ID) %>% 
  mutate(Avg_Volume_1day_lag=dplyr::lag(Avg_Volume_1day))

data_min %>% mutate(T_hour=T_min %/% 100) %>% group_by(TRD_DD, PROD_ID, T_hour) %>% 
  summarise(cnt=n()) %>% filter(cnt>31, cnt<60) %>% view()
