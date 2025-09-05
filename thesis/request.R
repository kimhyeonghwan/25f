rm(list=ls())

# Set local environment
setwd("/Users/hwan/KAIST/25fall/thesis/data")
o <- odbcConnect(dsn = "oracle", uid = "USFO_INF", pwd = "vkfksgksmf03!",
                 believeNRows = FALSE, DBMSencoding = "")

# Import library
library(tidyverse)
library(data.table)
library(e1071)
library(RODBC)

# Import csv
start_date <- 20150501
end_date <- 20181231

target_isu <- read_csv("target_isu.csv.gz") %>% 
  filter(TRD_DD>=start_date, TRD_DD<=end_date)
target_isu_wide <- target_isu %>%
  select(TRD_DD,PROD_ID,EXPMM_NO,ISU_CD) %>% 
  pivot_wider(names_from = c("PROD_ID","EXPMM_NO"), values_from = "ISU_CD")

year_first_day <- target_isu_wide %>% 
  mutate(TRD_YY=substr(TRD_DD,1,4)) %>% 
  group_by(TRD_YY) %>% 
  slice(1) %>% 
  ungroup() %>% 
  select(TRD_DD)

# CB history : remove sec-return after CB triggered
cb_history <- tibble(
  TRD_DD=c(20200313, 20200313, 20200319, 20200319, 20240805, 20240805, 20160212),
  PROD_ID=c("KRDRVFUK2I", "KRDRVFUKQI", "KRDRVFUK2I", "KRDRVFUKQI", "KRDRVFUK2I", "KRDRVFUKQI", "KRDRVFUKQI"),
  CB_event=c(104353, 90403, 120553, 120533, 141432, 135602, 115504)
)

# Price-limit history : remove sec_return during 5min+1sec after triggered
limit_history <- tibble(
  TRD_DD=c(20200323, 20240805, 20200318, 20190805, 20170612, 20231107),
  PROD_ID=c("KRDRVFUK2I", "KRDRVFUK2I", "KRDRVFUKQI", "KRDRVFUKQI", "KRDRVFUKQI", "KRDRVFUKQI"),
  Limit_start=c(090001, 133656, 152451, 151722, 090022, 121333),
  Limit_end=c(090501, 134156, 152951, 152222, 090522, 121833)
)

# Make min-level data & variables using transaction-level data

result_min <- list()
result_close <- list()
loop_i <- 1

for(i in 1:NROW(target_isu_wide)){
  day <- target_isu_wide$TRD_DD[i]
  isu <- c(target_isu_wide$KRDRVFUK2I_1[i],target_isu_wide$KRDRVFUKQI_1[i],
           target_isu_wide$KRDRVFUK2I_2[i],target_isu_wide$KRDRVFUKQI_2[i])
  
  # transaction level data to second level data
  tmp_second <- sqlQuery(o,paste0("SELECT TRD_DD, TRD_TM, ISU_CD, TRD_PRC, TRDVOL, TRD_NO",
                                  " FROM VWSV_DRV_TRD3",
                                  " WHERE trd_dd='",day,"'",
                                  " AND BRD_ID='G1' AND SESS_ID='40'")) %>%
    tibble() %>% 
    mutate(across(c(TRD_DD, TRD_TM, TRD_NO), as.integer),
           across(c(TRD_PRC, TRDVOL), as.double),
           ISU_CD=as.character(ISU_CD)) %>% 
    filter(ISU_CD %in% isu) %>%
    mutate(T_min=TRD_TM %/% 100000,
           T_sec=TRD_TM %/% 1000) %>%
    arrange(TRD_DD,ISU_CD,TRD_NO) %>%
    group_by(TRD_DD,ISU_CD,SESS_ID,T_min,T_sec) %>%
    summarise(VWAP=weighted.mean(TRD_PRC,TRDVOL),
              TRDVOL=sum(TRDVOL),
              .groups = "drop")
  
  # first day every year & 수능
  if(target_isu_wide$TRD_DD[i] %in% year_first_day){
    tmp_second <- tmp_second %>% mutate(day_grp=1)
  }else if(max(tmp_second$T_min)>1600){
    tmp_second <- tmp_second %>% mutate(day_grp=2,
                                        T_min=T_min-100, T_sec=T_sec-10000)
  }else{tmp_second <- tmp_second %>% mutate(day_grp=0)}
  
  tmp_second <- tmp_second %>% 
    filter(T_min>=900, T_min<=1530) %>% 
    group_by(TRD_DD, ISU_CD) %>% 
    arrange(TRD_DD, ISU_CD, T_sec) %>% 
    mutate(Ret_sec=log(VWAP)-log(dplyr::lag(VWAP))) %>%
    ungroup() %>% 
    # Treat CB & Price-limit history in Real market
    mutate(PROD_ID=if_else(substr(ISU_CD,5,6)=="01","KRDRVFUK2I","KRDRVFUKQI")) %>%
    left_join(.,cb_history, by=c("TRD_DD", "PROD_ID")) %>%
    left_join(.,limit_history, by=c("TRD_DD", "PROD_ID")) %>%
    mutate(CB_event=replace_na(CB_event, 999999),
           Limit_start=replace_na(Limit_start,999999),
           Limit_end=replace_na(Limit_end,999999)) %>% 
    mutate(CB_grp=if_else(T_sec>=CB_event,1,0),
           Limit_grp=if_else(T_sec>=Limit_start&T_sec<=Limit_end,1,0)) %>% 
    mutate(Ret_sec=if_else(CB_grp==1&dplyr::lag(CB_grp)==0,as.numeric(NA),Ret_sec)) %>% 
    mutate(Ret_sec=if_else(Limit_grp==1,as.numeric(NA),Ret_sec))
  
  tmp_min <- tmp_second %>% 
    group_by(TRD_DD, PROD_ID, ISU_CD, T_min) %>% 
    summarise(Ret_min=sum(Ret_sec, na.rm = TRUE),
              TSRV_5s=0,
              RDSkew=mean(Ret_sec^3, na.rm = TRUE) / (mean(Ret_sec^2, na.rm = TRUE)^(3/2)),
              P_low=min(VWAP),
              Volume=sum(TRDVOL),
              RV=sum(Ret_sec^2, na.rm = TRUE),
              TSRV_10s=0,
              AvgRV_5s=0,
              AvgRV_10s=0,
              Skew=skewness(Ret_sec, na.rm = TRUE),
              VWAP_min=weighted.mean(VWAP, TRDVOL),
              Count=n(),
              CB_grp=dplyr::last(CB_grp),
              day_grp=dplyr::first(day_grp),
              .groups = "drop")
  
  # Commute Autocorrelation of return's each min-level group
  tmp_min_auto <- tmp_second %>% 
    select(TRD_DD, ISU_CD, T_min, Ret_sec) %>% 
    group_by(TRD_DD, ISU_CD, T_min) %>% 
    mutate(Ret_sec_lag=dplyr::lag(Ret_sec)) %>% 
    summarise(Avg_Ret_lag=mean(Ret_sec_lag, na.rm = TRUE),
              Avg_Ret=mean(Ret_sec, na.rm = TRUE),
              Ret_cross=mean(Ret_sec * Ret_sec_lag, na.rm = TRUE),
              SD_Ret_lag=sd(Ret_sec_lag, na.rm = TRUE),
              SD_Ret=sd(Ret_sec, na.rm = TRUE),
              .groups = "drop") %>% 
    mutate(AutoCov = (Ret_cross - Avg_Ret_lag*Avg_Ret) / (SD_Ret_lag * SD_Ret)) %>% 
    select(TRD_DD, ISU_CD, T_min, AutoCov)
  
  tmp_min <- tmp_min %>% 
    left_join(., tmp_min_auto, by=c("TRD_DD", "ISU_CD", "T_min"))
  
  # calculate 5-second TSRV(Two Scale Realized Volatility)
  for(k in 0:4){
    ttmp_subsample <- tmp_second %>% 
      mutate(SubGrp=(T_sec%%100-k) %/% 5) %>% 
      filter(SubGrp>=0, SubGrp<11) %>% 
      group_by(TRD_DD, ISU_CD, T_min, SubGrp) %>% 
      summarise(Ret_sub=sum(Ret_sec, na.rm = TRUE), .groups = "drop") %>% 
      group_by(TRD_DD, ISU_CD, T_min) %>% 
      summarise(RV_sub=sum(Ret_sub^2, na.rm = TRUE), .groups = "drop")
    
    tmp_min <- tmp_min %>% 
      left_join(ttmp_subsample, by=c("TRD_DD","ISU_CD","T_min")) %>% 
      mutate(AvgRV_5s=AvgRV_5s+RV_sub/5) %>% # 5-times subsampling 
      select(-RV_sub)
    
    rm(ttmp_subsample)
  }
  
  tmp_min <- tmp_min %>% 
    mutate(TSRV_5s=(AvgRV_5s-(11/60)*RV)/(1-11/60)) # remove noise & scaling to minute level
  
  # calculate 10-second TSRV(Two Scale Realized Volatility)
  for(k in 0:9){
    ttmp_subsample <- tmp_second %>% 
      mutate(SubGrp=(T_sec%%100-k) %/% 10) %>% 
      filter(SubGrp>=0, SubGrp<5) %>% 
      group_by(TRD_DD, ISU_CD, T_min, SubGrp) %>% 
      summarise(Ret_sub=sum(Ret_sec, na.rm = TRUE), .groups = "drop") %>% 
      group_by(TRD_DD, ISU_CD, T_min) %>% 
      summarise(RV_sub=sum(Ret_sub^2, na.rm = TRUE), .groups = "drop")
    
    tmp_min <- tmp_min %>% 
      left_join(ttmp_subsample, by=c("TRD_DD","ISU_CD","T_min")) %>% 
      mutate(AvgRV_10s=AvgRV_10s+RV_sub/10) %>% # 10-times subsampling 
      select(-RV_sub)
    
    rm(ttmp_subsample)
  }
  
  tmp_min <- tmp_min %>% 
    mutate(TSRV_10s=(AvgRV_10s-(5/60)*RV)/(1-5/60)) # remove noise & scaling to original
  
  tmp_daily <- tmp_second %>% 
    filter(T_sec>=152930) %>% 
    group_by(TRD_DD, ISU_CD) %>% 
    summarise(P_close=weighted.mean(VWAP,TRDVOL),
              .groups = "drop")
  
  result_min[[loop_i]] <- tmp_min
  result_close[[loop_i]] <- tmp_daily
  
  print(paste0("Complete ",loop_i,"th : ",day))
  rm(tmp_second, tmp_min, tmp_daily, tmp_min_auto)
  
  loop_i <- loop_i+1
}

data_min <- bind_rows(result_min)
data_close <- bind_rows(result_close)

fwrite(data_min %>% data.table(),"request_min.csv.gz", compress = "gzip")
fwrite(data_close %>% data.table(),"request_close.csv.gz", compress = "gzip")
