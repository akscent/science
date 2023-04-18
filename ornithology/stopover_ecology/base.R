library(dplyr)
library(tidyr)
library(tidyverse)
library(lubridate)
library(readxl)
library(openxlsx)
base<- base %>%
  separate("Date_capture", c("Y_dc" , "M_dc", "D_dc"), sep = "-", remove = FALSE ) %>%
  separate("Date_banding", c("Y_db" , "M_db", "D_db"), sep = "-", remove = FALSE ) %>% 
  mutate(Jday_c = yday(as.Date(base$Date_capture)), .after = "D_dc") %>% 
  mutate(Season_c = ifelse(as.numeric(M_dc)<8, "S", "A"), .after = "Jday_c") %>% 
  mutate(Jday_b = yday(as.Date(base$Date_banding)), .after = "D_db") %>% 
  mutate(Season_b = ifelse(as.numeric(M_db)<8, "S", "A"), .after = "Jday_b")

names<-c('Key', 'Gene', 'G', 'count_rec', 'ring', 'Date_capture', "Y_dc" , "M_dc", "D_dc", 'Jday_c', 'Season_c', 
         'Date_recapture', "Y_dr" , "M_dr", "D_dr", 'Jday_r', 'Season_r', 'Recaptures', 'Date_recapture', 'Date_death', 'Genus', 'Sp', 'Ru_name', 'Count', 'Sex', 'Euring', 'Age', 'Pneumatization', 'Fat', 'Muscle', 'Weight', "beakfw", 'beak', 'wing_min', 'wing_max', 'Tarsus', 'Tarsus2', 'Tail', 'Tail_rost', 'Head', 'Net', 'Net2', 'Biotop', 'count_f_bio', 'Nets_number', 'Nets_nimber_big', 'Q_nets', 'Number', 'time', 'weather', 'other', 'collection', 'Place_capture', 'Place_recapture', 'Taper', 'Latitude', 'Longitude',
         'Body_length', 'Mouth_beak', 'Height_beak', 'Width_beak', 'Width_beak_nostril', 'Height_beak_nostril', 'Tip_wing', 'Al>BVKPM', 
         'Al<BVKPM',	'1Pm>BVKPM', '1Pm<BVKPM' ,'1Pm=BVKPM',	'Alula', '1Pm', '2Pm'	,'3Пм',	'4Pm',	'5Pm',	'6Pm',	'7Pm',	'8Pm',	'9Pm',	'10Pm',	'11Pm',	'2Pm>6Pm',	'2Pm<6Pm',	'2Pm=6Pm',	'2Pm>7Pm',	'2Pm<7Pm',	'2ПМ=7Pm',	'5Pm>6Pm',	'5Pm<6Pm',	'Length_white_brow',	'Length_black_brow',	'Code_tars')
names(base)<-names

write.xlsx(base, "base_01.xlsx", na = "NA")
base %>%
  mutate(Day_time_c = sum(Date_capture,time, na.rm = FALSE), .after = "Date_capture")
