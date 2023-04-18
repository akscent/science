#Достать данные из базы данных "bd"
#Связать базу данных с R
library(RODBC)
library(dplyr)
db<-"C:/Users/user/OneDrive/Документы/ds/bd.accdb"
con2 <- odbcConnectAccess2007(db)
#Узнать имена всех входящих таблиц, в т.ч. запросов
sqlTables(con2)$TABLE_NAME
#Подготовка датафрэйма. Запрос списка видов
listsp<-"SELECT DISTINCT Speciese FROM Morphology_sp"
#Вывести список видов как вектор
#Либо самостоятельно записать список видов как вектор в переменную
splist<- sqlQuery(con2, listsp)[,1] %>% as.vector()
#Для каждого вида сделать запрос и создать соответствующую таблицу
for (i in splist) {
  assign(paste0(i), sqlQuery(con2, paste("SELECT * FROM Morphology_sp WHERE Speciese = '",i,"';", sep ="")))
  
}

#Либо Записать запросы отдельно для нужного вида, если предусмотрено их разделение в разные таблицы
# Например 
# qry1<-"SELECT * FROM Morphology_sp WHERE Speciese = 'Ocyris tristrami';"
#Выбрать данные из подключенной базы данных в отдельные датафреймы по видам
# O_tristrami<-sqlQuery(con2, qry1)
