# set load
setwd("E:\\AT\\knockoff相关资料\\ukbcode\\data")
library(data.table)
library(dplyr)

#### read data ####
#ICD10 ICD10_date birth_year birth_month death_age lost_date
data_TTE <- as.data.table(fread("E:\\AT\\knockoff相关资料\\ukbcode\\data\\ukbb_phenotype0710.csv", 
                                select = c("eid", 
                                           paste0("41270-0.", 0:0),   #ICD10
                                           paste0("41280-0.", 0:242),   #ICD10_date
                                           "34-0.0",   #birth_year
                                           "52-0.0",   #birth_month
                                           "40007-0.0", "40007-1.0",   #death_age
                                           "191-0.0",   #lost_date
                                           "53-0.0", "53-1.0", "53-2.0", "53-3.0",   #attend_year
                                           "55-0.0", "55-1.0", "55-2.0", "55-3.0"   #attend_month
                                )))

# Convert date format(xxxxyear-xxmonth-xxday)
# 1.birth_date
#data_TTE$birth_date <- paste("X34-0.0", sprintf("%02d", "52-0.0"), sep = "-")
#data_TTE$birth_date <- as.Date(paste(data_TTE$birth_date, "01", sep = "-"), format = "%Y-%m-%d")
month_map <- c("January"=1, "February"=2, "March"=3, "April"=4, "May"=5, "June"=6,
               "July"=7, "August"=8, "September"=9, "October"=10, "November"=11, "December"=12)

data_TTE[, birth_date := as.Date(paste(
  `34-0.0`, 
  month_map[`52-0.0`], 
  "01", 
  sep="-"
))] 

head(data_TTE)

# 2.lost_date
data_TTE$lost_date <- NA
data_TTE$lost_date <- format(as.Date(data_TTE[["191-0.0"]]), "%Y-%m-%d")
table(data_TTE$lost_date)

# 3.latest date
#data_TTE$latest_date <- do.call(pmax, c(data_TTE[c(paste0("41280-0.", 0:242), "53-0.0", "53-1.0", "53-2.0", "53-3.0",
                                                   #"lost_date")], na.rm = TRUE))
col_names <- c(paste0("41280-0.", 0:242), "53-0.0", "53-1.0", "53-2.0", "53-3.0", "lost_date")  
if (all(col_names %in% names(data_TTE))) {  
  date_cols <- data_TTE[, ..col_names, with = FALSE]  
  date_cols <- lapply(date_cols, as.Date, format = "%Y-%m-%d")  
  data_TTE[, latest_date := do.call(pmax, c(as.list(date_cols), na.rm = TRUE))]  
} else {  
  stop("Not all column names are present in data_TTE")  
}

# 4.survival_months(death age*12)
data_TTE$survival_months <- ifelse(!is.na(data_TTE[["40007-0.0"]]),
                                   round(as.numeric(data_TTE[["40007-0.0"]]*12)),
                                   round(as.numeric(data_TTE[["40007-1.0"]]*12)))

#change values of IDC10
# 保留原始描述信息
data_TTE[,`41270-0.0_description` := `41270-0.0`]
library(data.table)
library(stringr)


# 拆分合并的ICD字符串
split_icd <- function(x) {
  if(is.na(x)) return(NA)
  codes <- unlist(strsplit(x, "\\|"))
  codes <- sub("^([A-Z][0-9]{2}\\.[0-9]).*", "\\1", codes) # 提取类似"I10"或"E11.9"的代码
  codes
}

# 应用到所有行
data_TTE[, icd_codes := lapply(`41270-0.0`, split_icd)]

# 展开到各列（最多242个诊断）
max_codes <- max(sapply(data_TTE$icd_codes, length))
for(i in 0:(max_codes-1)){
  data_TTE[, paste0("41270-0.",i) := sapply(icd_codes, function(x) ifelse(length(x)>i, x[i+1], NA))]
}

for (i in 0:242) {
  column_name <- paste0("41270-0.", i)
  pattern <- "^(.{3})(\\d+)$"
  data_TTE[[column_name]] <- gsub(pattern, "\\1.\\2", data_TTE[[column_name]])
}

#### match_code_I10 ####
## definite Phecode and survival time
## phecode
data_TTE$Phecode_I10 <- 0

for (i in 0:242) {
  icd_col <- paste0("41270-0.", i)
  time_col <- paste0("41280-0.", i)
  data_TTE$Phecode_I10 <- ifelse(!is.na(data_TTE[[icd_col]]) & !is.na(data_TTE[[time_col]]) & grepl("I10", data_TTE[[icd_col]]), 1, data_TTE$Phecode_I10)
}
#table(data_TTE$Phecode_I10)

## calculate survival time
# Phecode==0 (censoring)
data_TTE$survival_time_I10 <- ifelse(data_TTE$Phecode_I10 == 0,
                                     ifelse(!is.na(data_TTE$survival_months),
                                            pmax(data_TTE$survival_months, 
                                                 round(as.numeric(difftime(as.Date(data_TTE$latest_date, format="%Y-%m-%d"), 
                                                                           as.Date(data_TTE$birth_date, format="%Y-%m-%d"), 
                                                                           units = "days") / 30))),
                                            round(as.numeric(difftime(as.Date(data_TTE$latest_date, format="%Y-%m-%d"), 
                                                                      as.Date(data_TTE$birth_date, format="%Y-%m-%d"), 
                                                                      units = "days") / 30))),
                                     NA)

# Phecode==1
for (i in 0:242) {
  icd_col <- paste0("41270-0.", i)
  time_col <- paste0("41280-0.", i)
  
  indexes <- grep("I10", data_TTE[[icd_col]])
  
  if (length(indexes) > 0) {
    data_TTE$survival_time_I10[indexes] <- round(as.numeric(difftime(as.Date(data_TTE[[time_col]][indexes]), 
                                                                 as.Date(data_TTE$birth_date[indexes], "%Y-%m-%d"), 
                                                                 units = "days") / 30))
  }
}


#### match_code_Type 2 Diabetes(T2D) (E11.2- E11.9) ####
## definite Phecode and survival time
## phecode
# definite Phecode and survival time
data_TTE$Phecode_T2D <- 0

for (i in 0:242) {
  icd_col <- paste0("41270-0.", i)
  time_col <- paste0("41280-0.", i)
  
  data_TTE$Phecode_T2D <- ifelse(!is.na(data_TTE[[icd_col]]) & !is.na(data_TTE[[time_col]]) & 
                                           grepl("E11.2|E11.3|E11.4|E11.5|E11.6|E11.7|E11.8|E11.9", data_TTE[[icd_col]]), 1, data_TTE$Phecode_T2D)
}
#table(data_TTE$Phecode_T2D)

# Phecode==0 (censoring)
data_TTE$survival_time_T2D <- ifelse(data_TTE$Phecode_T2D == 0,
                                     ifelse(!is.na(data_TTE$survival_months),
                                            pmax(data_TTE$survival_months, 
                                                 round(as.numeric(difftime(as.Date(data_TTE$latest_date, format="%Y-%m-%d"), 
                                                                           as.Date(data_TTE$birth_date, format="%Y-%m-%d"), 
                                                                           units = "days") / 30))),
                                            round(as.numeric(difftime(as.Date(data_TTE$latest_date, format="%Y-%m-%d"), 
                                                                      as.Date(data_TTE$birth_date, format="%Y-%m-%d"), 
                                                                      units = "days") / 30))),
                                     NA)
# Phecode==1(failure)
for (i in 0:242) {
  icd_col <- paste0("41270-0.", i)
  time_col <- paste0("41280-0.", i)
  
  indexes <- grep("E11.2|E11.3|E11.4|E11.5|E11.6|E11.7|E11.8|E11.9", data_TTE[[icd_col]])
  
  if (length(indexes) > 0) {
    data_TTE$survival_time_T2D[indexes] <- round(as.numeric(difftime(min(as.Date(data_TTE[[time_col]][indexes])), 
                                                                             as.Date(data_TTE$birth_date[indexes], "%Y-%m-%d"), 
                                                                             units = "days") / 30))
  }
}

# save
data_TTE_df <- data_TTE[, .SD, .SDcols = !"icd_codes"]
write.table(data_TTE_df,file="ukbb_phenotype0714.csv",row.names=FALSE,col.names=TRUE,sep=",")




