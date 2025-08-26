# running Global setting
source("1.setting.R")

# set parameter
master_strt <- 2018
master_end <- 2023
target_strt <- 20190101
target_end <- 20231231

# daily data
master <- tibble()
for(i in master_strt:master_end){
  tmp <- read_csv(gzfile(paste0("data/DRV_Master/",i,"_Futures_Master.csv.gz"))) %>%
    tibble() %>% 
    filter(PROD_ID %in% c("KRDRVFUK2I","KRDRVFUKQI"),
           EXPMM_NO<=2, substr(ISU_CD,4,4)==1)
  master <- master %>% bind_rows(tmp)
}

# expiry day & (T-1) day
lsttrd_dd <- master %>% 
  mutate(TRD_MM=TRD_DD %/% 100) %>% 
  filter(PROD_ID=="KRDRVFUK2I", EXPMM_NO==1, as.integer(substr(TRD_MM,5,6)) %% 3==0) %>% 
  select(TRD_MM, TRD_DD,remain_dys) %>% 
  arrange(TRD_MM, remain_dys) %>% 
  group_by(TRD_MM) %>% 
  slice(1:2) %>% 
  ungroup() %>% 
  select(TRD_DD)

master_lead <- master %>% filter(EXPMM_NO==1)
master_second <- master %>% filter(EXPMM_NO==2)

lsttrd_second <- lsttrd_dd %>% left_join(master_second)

# target is lead month + second month(ONLY expiry day)
target_isu <- master_lead %>%
  union_all(lsttrd_second) %>% 
  filter(TRD_DD>=target_strt,
         TRD_DD<=target_end) %>% 
  select(TRD_DD,PROD_ID,ISU_CD,TDD_CLSPRC,ACC_TRDVOL,ACC_TRDVAL,EXPMM_NO,EXPMM,remain_dys) %>% 
  arrange(TRD_DD,PROD_ID,EXPMM_NO)

target_isu_wide <- target_isu %>%
  select(TRD_DD,PROD_ID,EXPMM_NO,ISU_CD) %>% 
  pivot_wider(names_from = c("PROD_ID","EXPMM_NO"), values_from = "ISU_CD")

# using transaction level data -> second level data
for(k in 2019:2019){
  tmp_target <- target_isu_wide %>% filter(substr(TRD_DD,1,4)==as.character(k))
  tmp_second <- tibble()
  for(i in 1:NROW(tmp_target)){
    tmp <- read_csv(gzfile(paste0("data/DRV_TRD/",tmp_target$TRD_DD[i],"_DRV_TRD.csv.gz"))) %>%
      tibble() %>% 
      filter(ISU_CD %in% c(tmp_target$KRDRVFUK2I_1[i],tmp_target$KRDRVFUKQI_1[i],tmp_target$KRDRVFUK2I_2[i],tmp_target$KRDRVFUKQI_2[i])) %>% 
      filter(BRD_ID=="G1") %>% 
      mutate(T_min=TRD_TM %/% 100000,
             T_sec=TRD_TM %/% 1000) %>% 
      arrange(TRD_DD,ISU_CD,TRD_NO) %>% 
      group_by(TRD_DD,ISU_CD,SESS_ID,T_min,T_sec) %>% 
      summarise(VWAP=weighted.mean(TRD_PRC,TRDVOL),
                TRDVOL=sum(TRDVOL),
                AVG_PRC=mean(TRD_PRC),
                MID_PRC=dplyr::last((ASK_STEP1_BSTORD_PRC+BID_STEP1_BSTORD_PRC)/2),
                MICRO_PRC=dplyr::last((ASK_STEP1_BSTORD_RQTY*BID_STEP1_BSTORD_PRC+BID_STEP1_BSTORD_RQTY*ASK_STEP1_BSTORD_PRC)/(ASK_STEP1_BSTORD_RQTY+BID_STEP1_BSTORD_RQTY)),
                OPN=dplyr::first(TRD_PRC),
                HG=max(TRD_PRC),
                LW=min(TRD_PRC),
                CLS=dplyr::last(TRD_PRC)) %>% 
      ungroup()
    tmp_second <- tmp_second %>% bind_rows(tmp)
  }
  
  fwrite(tmp_second %>% data.table(), paste0(k,"_target_second.csv.gz"), compress = "gzip")
}

# 신년 개장일-개장만 1시간늦게, 수능일-1시간순연 따로 처리 필요..
# 특이한 날을 따로 발라내서 9~1530으로 시간을 맞추자
k200_sample <- target_second %>% 
  filter(substr(ISU_CD,5,6)=="01",SESS_ID==40, T_min<1530) %>% 
  select(TRD_DD,ISU_CD,T_min,T_sec,VWAP,TRDVOL)

# last 30second VWAP, (P_t close)
# 샘플늘려봐야혀
daily_threshold <- k200_sample %>% 
  filter(T_sec>=152930, T_sec<=152959) %>% 
  group_by(TRD_DD,ISU_CD) %>% 
  summarise(P_close=weighted.mean(VWAP,TRDVOL)) %>% 
  ungroup() %>% 
  mutate(lag_P_close=dplyr::lag(P_close))


