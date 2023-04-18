#Подготовка данных
#Замена запятых точками и преобразлвание числовых переменных в класс numeric а факторных в класс factor
#speciese<-list(`t.Ocyris aureolus`, `t.Ocyris chrysophrys`, `t.Ocyris personatus`, `t.Ocyris pusillus`, 
#               `t.Ocyris rusticus`, `t.Ocyris rutilus`, `t.Ocyris spodocephalus`, `t.Ocyris tristrami`, 
#              `t.Ocyris variabilis`, `t.Ocyris yessoensis`)
#создаём функцию для замены точками запятых в определенных столбцах
funt_df_asnum <- function(df) {
  #здесь suppressWarnings заставляет принимать ошибки по возникновению NA, а sapply возвращает вектор данных, а не список, как lapply
  df[,13:25]<-sapply(df[,13:25],function(x) suppressWarnings(as.numeric(sub(",", ".", x, fixed = TRUE))))
  #и вернём датафрейм
  return(df)
}
#тоже создать список доступных датафреймов в оккружении, название которых наичнается с Ocyris
#nms <- ls(pattern = "^Ocyris*", envir = e) #
#добавить переменную, которая будет означать глобальное окружение
e <- .GlobalEnv
#запустить цикл (loop) для функции 
for(i in splist) e[[i]] <- funt_df_asnum(e[[i]])
#сделать тоже самое дальше
#for (i in splist) {
#  output[[paste("'",i,"'", sep ="")[,13:25]]]<- sapply(paste("'",i,"'", sep ="")[,13:25], function(x) as.numeric(sub(",", ".", x, fixed = TRUE)))
#}

`Ocyris tristrami`[,13:25]<-sapply(`Ocyris tristrami`[,13:25], function(x) as.numeric(sub(",", ".", x, fixed = TRUE)))
`Ocyris tristrami`[,33:41]<-sapply(`Ocyris tristrami`[,33:41], function(x) as.numeric(sub(",", ".", x, fixed = TRUE)))
`Ocyris tristrami`[,11:12]<-sapply(`Ocyris tristrami`[,11:12], function(x) as.factor(x))
data<-`Ocyris tristrami`
data<-data[, c(1,2,3,4,6,7,9,15,13,10,26,28,30,31,32,27,5,8,14,29,11,12,16,17,18,19,20,21,23,25,33,34,35,36,37,38,39,40,41)]

#_______________________________________________________________________________________________________
library(ggplot2)
library(ggcorrplot)
library(effectsize)
library(stargazer)
library(sjPlot)
library(gtsummary)
library(ggstatsplot)
library(ggsci)
library("ggdensity")
library(plotly)
library(ggfortify)
library(missMDA)
library(tidyr)
library(dplyr)
library(hrbrthemes)
library(ggridges)
library(viridis)
library(hrbrthemes)
library(flextable)
library(MVTests)
library(PCAtools)
library(FactoMineR)
library(DESeq2)
library(factoextra)
library(Factoshiny)
library(gt)
library(ggfortify)
library(psych)
library(MASS)
library(tidyverse)
library(modelr)
library(janitor)
library(skimr)
library(broom)
library(corrplot)
library(class)
library(klaR)
library(ROCR)
library(yarrr)
library(multcomp)
library(devtools)
library(stats)
library(nsprcomp)
library(tibble)
library(ggplotgui)
library(reshape2)
library(broom)

#_____________________________
#Удаление выбросов в данных

# create detect outlier function
detect_outlier<-function(x) {
  # calculate first quantile
  Quantile1<- quantile(x, probs=.25)
  # calculate third quantile
  Quantile3<- quantile(x, probs=.75)
  # calculate inter quartile range
  IQR=Quantile3-Quantile1
  # return true or false
  x>Quantile3 + (IQR*1.5) | x<Quantile1 - (IQR*1.5)
}
# create remove outlier function
remove_outlier<- function(dataframe,
                            columns=names(dataframe)) {
  # for loop to traverse in columns vector
  for (col in columns) {
    # remove observation if it satisfies outlier function
    dataframe<- dataframe[!detect_outlier(dataframe[[col]]), ]
  }
  # return dataframe
  print("Remove outliers")
  print(dataframe)
}
#Удаление NA ( 7-8 ПМ)
df<-data[,-40]
df<-df[,-40]
#
df_cca<-na.omit(df)
df_pm<-df_pm[,-3]
df_pm<-df_pm[,-3]
df_pm<-df_pm[,-3]
df_pm<-df_pm[,-3]
df_pm<-df_pm[,-3]
df_pm<-df_pm[,-3]
df_pm<-df_pm[,-3]
df_pm<-df_pm[,-3]
#удалить ненужные ст
df_pm<-na.omit(data)
#Удаление выбросов и новый фрейм данных df_pm, получилось всего 367 наблюдений
df_pm<-remove_outlier(df_pm, c("Alula"  ,  "1Pm",
                        "2Pm"  ,    "3Pm"  ,    "4Pm"  ,    "5Pm"   ,   "6Pm" ))
df_pm$Age[df_pm$Age=='juv']<-'HY'
df_pm$Age[df_pm$Age=='sad']<-'HY'
df_pm$Age[df_pm$Age=='ad']<-'AHY'
df_pm<-mutate(df_pm, Sex_Age=paste(df_pm$Sex,df_pm$Age))
#Корреляционная таблица
M <- cor(df_pm[,4:18], method = "pearson")
ggcorrplot(M, hc.order = TRUE, type = "lower",
           colors = c("white","yellow","purple" ), lab = TRUE)
#____________________________________________________________________________________________________________________
#Сделаем тоже самое для промеров без ПМ, чтобы оценить больше наблюдений
df<-na.omit(data[,20:30])
df<-remove_outlier(df, c("Weight",   "beakfw",   "beak" ,    "wing_min", "wing_max", "Tarsus"  , "Tail"   ,
                               "Head"))
M1 <- cor(df[,4:11], method = "pearson")
ggcorrplot(M1, hc.order = TRUE, type = "lower",
           colors = c("white","yellow","purple" ), lab = TRUE)
#Удление неопределенных пола и возраста
df<- df %>% filter(Sex!='0')
df<- df %>% filter(Age!='0')
#преобразование в длинную переменную
df_long <- gather(df[,2:11], key = "measure", value = "measurement", Weight, beakfw, beak, wing_min, wing_max, Tarsus, Tail, Head) 
df_long<-mutate(df_long, Sex_Age=paste(df_long$Sex,df_long$Age))

#Провекра на многомерную нормальность
library(mvnormtest)
mshapiro.test(t(df[,4:11]))
#Если нет многомерной нормальности, делаем попарный АНОВА, группируемый по измерениям
df_long %>% 
  group_by(measure) %>%
  group_modify(~broom::tidy(TukeyHSD(aov(measurement~Sex_Age, data = .x)))) %>% as.data.frame() %>% 
  filter(adj.p.value < .05)
#Если нет значимых различий с juv, заменяем в данных juv & sad на HY, ad на AHY
df_long$Age[df_long$Age=='juv']<-'HY'
df_long$Age[df_long$Age=='sad']<-'HY'
df_long$Age[df_long$Age=='ad']<-'AHY'
df_long<-mutate(df_long[,-5], Sex_Age=paste(df_long$Sex,df_long$Age))
#
df$Age[df$Age=='juv']<-'HY'
df$Age[df$Age=='sad']<-'HY'
df$Age[df$Age=='ad']<-'AHY'
df<-mutate(df[,-12], Sex_Age=paste(df$Sex,df$Age))
#И снова повторяем АНОВА с фильтрованием по значимым сравнениям
table.aov<-df_long %>% 
  group_by(measure) %>%
  group_modify(~broom::tidy(TukeyHSD(aov(measurement~Sex_Age, data = .x)))) %>% as.data.frame() %>% 
  filter(adj.p.value < .05)
#Удобная интерактивная таблица
library(reactable)
reactable(table.aov, groupBy = 'measure')
#Визуализация различий
ggplot(df_long, aes(x = measure, y = measurement, color = Sex, fill=Age)) +
  geom_boxplot() +
  facet_wrap(~measure, scale="free")
grouped_ggbetweenstats(
  plot.type = "boxviolin",
  type = 'p',
  data  = df_long,
  x     = Sex_Age,
  xlab = 'Половозрастная группа',
  ylab = 'Измерение',
  y     = measurement,
  outlier.shape = NA,
  pairwise.display = "s",
  p.adjust.method  = "hochberg",
  palette          = "default_jama",
  package          = "ggsci",
  grouping.var = measure,
  plotgrid.args = list(nrow = 3, ncol = 3)
)
# И только по полу
grouped_ggbetweenstats(
  plot.type = "boxviolin",
  type = 'p',
  data  = df_long,
  x     = Sex,
  xlab = 'Пол',
  ylab = 'Измерение',
  y     = measurement,
  outlier.shape = NA,
  pairwise.display = "s",
  p.adjust.method  = "hochberg",
  palette          = "default_jama",
  package          = "ggsci",
  grouping.var = measure,
  plotgrid.args = list(nrow = 3, ncol = 3)
)
#Исходя из размера эффекта выбираем самые подходящие для классификационного моделирования признаки: wing, tail, weight, head
#необходимо оценить cca признаки по длинам ПМ, df_pm будет использоваться для моделирования
#проверить ALULA на различия и другие промеры, доступные в общей бд, если они есть
df_pm<-df_pm %>% mutate(C1 = 4.705 * 10 - 4*`1Pm`^0.209*`2Pm`^0.200*`3Pm`^0.198*`4Pm`^0.195*`5Pm`^0.192*`6Pm`^0.192, .after = `6Pm`)
df_pm<-df_pm %>% mutate(C2 = 3.332 * `1Pm`^-3.49*`2Pm`^-1.816*`3Pm`^-0.893*`4Pm`^-0.003*`5Pm`^0.829*`6Pm`^1.351, .after = C1)
df_pm<-df_pm %>% mutate(C3 = 0.0879*`1Pm`^-6.231*`2Pm`^1.683*`3Pm`^4.033*`4Pm`^4.721*`5Pm`^3.955*`6Pm`^1.349, .after = C2)
#здесь сезон S редставлен всего 2 записями, их можно удалить
df_pm<-filter(df_pm, Season_c!='S')
#И снова повторяем АНОВА с фильтрованием по значимым сравнениям
df_pm_l<-gather(df_pm, key = "measure", value = "measurement", Alula, `1Pm`, `2Pm`, `3Pm`, `4Pm`, `5Pm`, `6Pm`, C1, C2, C3) 
table.aov.pm<-df_pm_l %>% 
  group_by(measure) %>%
  group_modify(~broom::tidy(TukeyHSD(aov(measurement~Sex, data = .x)))) %>% as.data.frame() %>% 
  filter(adj.p.value < .05)
reactable(table.aov.pm, groupBy = 'measure')
grouped_ggbetweenstats(
  plot.type = "boxviolin",
  type = 'p',
  data  = df_pm_l,
  x     = Sex,
  xlab = 'Пол',
  ylab = 'Измерение',
  y     = measurement,
  outlier.shape = NA,
  pairwise.display = "s",
  p.adjust.method  = "hochberg",
  palette          = "default_jama",
  package          = "ggsci",
  grouping.var = measure,
  plotgrid.args = list(nrow = 3, ncol = 4)
)
#
#Данные для моделирования
df_mod<-data[20:37]
df_mod<- df_mod %>% filter(Sex!='0')
df_mod<-df_mod %>% mutate(C1 = 4.705 * 10 - 4*`1Pm`^0.209*`2Pm`^0.200*`3Pm`^0.198*`4Pm`^0.195*`5Pm`^0.192*`6Pm`^0.192, .after = `6Pm`)
df_mod<-df_mod %>% mutate(C2 = 3.332 * `1Pm`^-3.49*`2Pm`^-1.816*`3Pm`^-0.893*`4Pm`^-0.003*`5Pm`^0.829*`6Pm`^1.351, .after = C1)
df_mod<-df_mod %>% mutate(C3 = 0.0879*`1Pm`^-6.231*`2Pm`^1.683*`3Pm`^4.033*`4Pm`^4.721*`5Pm`^3.955*`6Pm`^1.349, .after = C2)
df_mod<- df_mod %>% filter(Sex!='0')
df_mod<- df_mod %>% filter(Age!='0')
df_mod$Age[df_mod$Age=='juv']<-'HY'
df_mod$Age[df_mod$Age=='sad']<-'HY'
df_mod$Age[df_mod$Age=='ad']<-'AHY'
#

#________________________________________________________________________________________________
#моделирование
#EFA Факторный анализ
fa.parallel(M,n.obs=112,fa="both",n.iter=100, main = "Screeplots with parrallel analysis")
fa<-fa(M,nfactors=2,rotate="none",fm="pa")
fa
fa.promax<-fa(M,nfactors=2,rotate="promax",fm="pa")
fa.promax
fa.diagram(fa.promax,simple=FALSE, main = "Факторные нагрузки")

#Создание таблиц сравнения между мужчинами и женщинами

add_ES_sad <- filter(df_mod, Age == 'HY') %>% select(Sex, Head, Tail, Tarsus, Weight, wing_max, beak, beakfw, Alula, `1Pm`, `2Pm`, `3Pm`, `4Pm`, `5Pm`, `6Pm`, C1, C2, C3)%>% 
  tbl_summary(by = Sex, 
              label = list(Head ~ "Голова, см", Tail ~ "Хвост, см", Tarsus ~ "Цевка, см", Weight ~ "Вес, г", wing_max ~ "Крыло, см", beak ~"Клюв от ноздри, см",  `beakfw` ~ "Клюв от лба, см", `Alula` ~ "Алула, см", `1Pm` ~ "1 ПМ, см" , `2Pm` ~ "2 ПМ, см" , `3Pm` ~ "3 ПМ, см" , `4Pm` ~ "4 ПМ, см" , `5Pm` ~ "5 ПМ, см" , `6Pm` ~ "6 ПМ, см", C1 ~ 'Изометричный размер', C2 ~ 'Заострённость', C3 ~ 'Выпуклость'),
              statistic = list(all_continuous2() ~ "{mean} ({sd})"),  missing = "no")  %>% add_difference(test = list(all_continuous()~"t.test", all_continuous()~"cohens_d")) %>% add_p(test = everything() ~ t.test) %>% bold_p() %>% 
              modify_header (label = "**Характеристика**") 
add_ES_sad
add_ES_ad <- filter(df_mod, Age == 'AHY') %>% select(Sex, Head, Tail, Tarsus, Weight, wing_max, beak, beakfw, Alula, `1Pm`, `2Pm`, `3Pm`, `4Pm`, `5Pm`, `6Pm`, C1, C2, C3)%>% 
  tbl_summary(by = Sex, 
              label = list(Head ~ "Голова, см", Tail ~ "Хвост, см", Tarsus ~ "Цевка, см", Weight ~ "Вес, г", wing_max ~ "Крыло, см", beak ~"Клюв от ноздри, см",  `beakfw` ~ "Клюв от лба, см", `Alula` ~ "Алула, см", `1Pm` ~ "1 ПМ, см" , `2Pm` ~ "2 ПМ, см" , `3Pm` ~ "3 ПМ, см" , `4Pm` ~ "4 ПМ, см" , `5Pm` ~ "5 ПМ, см" , `6Pm` ~ "6 ПМ, см", C1 ~ 'Изометричный размер', C2 ~ 'Заострённость', C3 ~ 'Выпуклость'),
              statistic = list(all_continuous2() ~ "{mean} ({sd})"),  missing = "no")  %>% add_difference(test = list(all_continuous()~"t.test", all_continuous()~"cohens_d")) %>% add_p(test = everything() ~ t.test) %>% bold_p() %>% 
  modify_header (label = "**Характеристика**") 
add_ES_ad

#В заострённости и выпуклости не обнражуено различий у взрослых, только в размерах крыла, значит пм можно не учитывать в модели
#Однако обнаружено у молодых особей
tbl_Es<- as_gt(tbl_merge(tbls = list(add_ES_sad, add_ES_ad), tab_spanner = c("**Первый год жизни**", "**Второй и более год жизни**")))
tbl_Es
Mpm <- cor(na.omit(df_mod[,13:21]), method = "pearson")
ggcorrplot(Mpm, hc.order = TRUE, type = "lower",
           colors = c("white","yellow","purple" ), lab = TRUE)
#В общем ПМ можно исключать
#______________________________________________________________________________________________________

#Соедигнить таблицы потом по видам
#______________________________________________________________________________________________________

box<-boxM(cbind(data[,13], data[,15], df_mod[,21]), df_mod$Sex) 
#Chi-Sq (approx.) = 17.961, df = 6, p-value = 0.006331
box
plot(box)
#результаты указывают на то, что ковариационные матрицы не отличаются, т.е. можножно применять LDA
#нарушение допущения однородности дисперсионно-ковариационных матриц

sample <- sample(c(TRUE, FALSE), nrow(data), replace=TRUE, prob=c(0.7,0.3))
train <- data[sample, ]
train<-clean_names(train)
test <- data[!sample, ]
test<-clean_names(test)
#QDA
modelall<-qda(sex~tail+wing_max+x2pm, data = train)

partimat(as.factor(sex)~tail + wing_max + x2pm, data = train, method = 'qda', nplots.hor = 3, main = "Квадратичный дисперсионный анализ")

p2_train_all <- predict(modelall, train)$class
tab2_train_all <- table(Predicted = p2_train_all, Actual = train$sex)
tab2_train_all
qda_train_accuracy_all <- sum(diag(tab2_train_all))/sum(tab2_train_all) * 100
qda_train_accuracy_all
#79.88615%
#LDA
ldamodelall<-lda(sex~tail+wing_max+x2pm, data = train)

partimat(as.factor(sex)~tail + wing_max + x2pm, data = train, method = 'lda', nplots.hor = 3, main = "Квадратичный дисперсионный анализ")

plot(ldamodelall)
coef(ldamodelall)
p2_train_all_lda <- predict(ldamodelall, train)$class
tab2_train_all_lda <- table(Predicted = p2_train_all_lda, Actual = train$sex)
tab2_train_all_lda
qda_train_accuracy_all_lda <- sum(diag(tab2_train_all_lda))/sum(tab2_train_all_lda) * 100
qda_train_accuracy_all_lda
#80.0759
#RDA
library(mlbench)
library(caret)
library(glmnet)
library(klaR)
cv_5_grid = trainControl(method = "cv", number = 5)

trrdamodelall<-train(sex~tail+wing_max+x2pm, data = train, trControl = cv_5_grid)
partimat(as.factor(sex)~tail + wing_max + x2pm, data = train, method = 'rda', nplots.hor = 3, main = "Квадратичный дисперсионный анализ")

preproc.param <- train %>%
  preProcess(method = c("center", "scale"))
train.transformed <- preproc.param %>% predict(train)
test.transformed <- preproc.param %>% predict(test)

rdamodelall<-vegan::rda(sex~tail+wing_max+x2pm, data = train)
rdamodelall

p2_train_all_rda <- predict(rdamodelall, train)$class
tab2_train_all_rda <- table(Predicted = p2_train_all_rda, Actual = train$sex)
tab2_train_all_rda
qda_train_accuracy_all_rda <- sum(diag(tab2_train_all_rda))/sum(tab2_train_all_rda) * 100
qda_train_accuracy_all_rda
#80.83491

#
library(ggvegan)

#MDA
library(mda)
# Fit the model
modelallmda <- mda(sex~tail+wing_max+x2pm, data = train)
modelallmda
# Make predictions
predicted.classes <- modelallmda %>% predict(test)
# Model accuracy
mean(predicted.classes == test$sex)
#0.7783251%
plot(modelallmda, data = train, group = "pred")

#FDA
# Fit the model
modelallfda <- fda(sex~tail+wing_max+x2pm, data = train)
modelallfda
# Make predictions
predicted.classes <- modelallfda %>% predict(test)
# Model accuracy
mean(predicted.classes == test$sex)
#0.7832512%

plot(modelallfda, data = train, group = "pred")
#Хвост+Крыло
#QDA
modelt_w<-qda(sex~tail+wing_max, data = train)

#LDA
ldamodelt_w<-lda(sex~tail+wing_max, data = data)
p2_train_all_lda <- predict(ldamodelt, data)$class
tab2_train_all_lda <- table(Predicted = p2_train_all_lda, Actual = data$sex)
tab2_train_all_lda
train_accuracy_all_lda <- sum(diag(tab2_train_all_lda))/sum(tab2_train_all_lda) * 100
train_accuracy_all_lda

#Крыло + 2ПМ
#QDA
modelw_pm<-qda(sex~wing_max+x2pm, data = train)
#LDA
ldamodelw_pm<-lda(sex~wing_max+x2pm, data = train)
#Хвост + 2ПМ
#QDA
modelt_pm<-qda(sex~tail+x2pm, data = train)
#LDA
ldamodelt_pm<-lda(sex~tail+x2pm, data = train)
#Хвост
#QDA
modelt<-qda(sex~tail, data = train)
#LDA
ldamodelt<-lda(sex~tail, data = train)
#Крыло
#QDA
modelw<-qda(sex~wing_max, data = train)
p2_train_all_lda <- predict(modelw, train)$class
tab2_train_all_lda <- table(Predicted = p2_train_all_lda, Actual = train$sex)
tab2_train_all_lda
train_accuracy_all_lda <- sum(diag(tab2_train_all_lda))/sum(tab2_train_all_lda) * 100
train_accuracy_all_lda
#81.97343
coef(modelw)

#LDA
ldamodelw<-lda(sex~wing_max, data = train)
#2ПМ
#QDA
modelpm<-qda(sex~x2pm, data = train)
#LDA
ldamodelpm<-lda(sex~x2pm, data = train)

#AUC for all
allqda <- predict(modelall, newdata = test)
alllda <- predict(ldamodelall, newdata = test)
t_wqda<- predict(modelt_w, newdata = test)
t_wlda <- predict(ldamodelt_w, newdata = test)
w_pmqda<- predict(modelw_pm, newdata = test)
w_pmlda <- predict(ldamodelw_pm, newdata = test)
t_pmqda<- predict(modelt_pm, newdata = test)
t_pmlda <- predict(ldamodelt_pm, newdata = test)
t_qda<- predict(modelt, newdata = test)
t_lda <- predict(ldamodelt, newdata = test)
w_qda<- predict(modelw, newdata = test)
w_lda <- predict(ldamodelw, newdata = test)
pm_qda<- predict(modelpm, newdata = test)
pm_lda <- predict(ldamodelpm, newdata = test)

#
library(ROCR)
library(pROC)
rallqda<-roc(test$sex, allqda$posterior[,2])
ralllda <- roc(test$sex, alllda$posterior[,2])
rt_wqda<- roc(test$sex, t_wqda$posterior[,2])
rt_wlda <- roc(test$sex, t_wlda$posterior[,2])
rw_pmqda<- roc(test$sex, w_pmqda$posterior[,2])
rw_pmlda <- roc(test$sex, w_pmlda$posterior[,2])
rt_pmqda<- roc(test$sex, t_pmqda$posterior[,2])
rt_pmlda <- roc(test$sex, t_pmlda$posterior[,2])
rt_qda<- roc(test$sex, t_qda$posterior[,2])
rt_lda <- roc(test$sex, t_lda$posterior[,2])
rw_qda<- roc(test$sex, w_qda$posterior[,2])
rw_lda <- roc(test$sex, w_lda$posterior[,2])
rpm_qda<- roc(test$sex, pm_qda$posterior[,2])
rpm_lda <- roc(test$sex, pm_lda$posterior[,2])
list(rallqda$auc,
     ralllda$auc,
     rt_wqda$auc,
     rt_wlda$auc,
     rw_pmqda$auc,
     rw_pmlda$auc,
     rt_pmqda$auc,
     rt_pmlda$auc,
     rt_qda$auc,
     rt_lda$auc,
     rw_qda$auc,
     rw_lda$auc,
     rpm_qda$auc,
     rpm_lda$auc)

#par(mfrow=c(1, 2))
prediction(allqda$posterior[,2], test$sex) %>%
  performance(measure = "tpr", x.measure = "fpr") %>% 
  plot(main = "ROC", xlab = "Доля ложных положительных оценок", ylab = "Доля правильных положительных оценок", col = 1)

prediction(alllda$posterior[,2], test$sex) %>%
  performance(measure = "tpr", x.measure = "fpr") %>%
  plot(col = 2, add = TRUE)
prediction(t_wqda$posterior[,2], test$sex) %>%
  performance(measure = "tpr", x.measure = "fpr") %>%
  plot(col = 3, add = TRUE)
prediction(t_wlda$posterior[,2], test$sex) %>%
  performance(measure = "tpr", x.measure = "fpr") %>%
  plot(col = 4, add = TRUE)
prediction(w_pmqda$posterior[,2], test$sex) %>%
  performance(measure = "tpr", x.measure = "fpr") %>%
  plot(col = 5, add = TRUE)
prediction(w_pmlda$posterior[,2], test$sex) %>%
  performance(measure = "tpr", x.measure = "fpr") %>%
  plot(col = 6, add = TRUE)
prediction(t_pmqda$posterior[,2], test$sex) %>%
  performance(measure = "tpr", x.measure = "fpr") %>%
  plot(col = 7, add = TRUE)
prediction(t_pmlda$posterior[,2], test$sex) %>%
  performance(measure = "tpr", x.measure = "fpr") %>%
  plot(col = 8, add = TRUE)
prediction(t_qda$posterior[,2], test$sex) %>%
  performance(measure = "tpr", x.measure = "fpr") %>%
  plot(col = 9, add = TRUE)
prediction(t_lda$posterior[,2], test$sex) %>%
  performance(measure = "tpr", x.measure = "fpr") %>%
  plot(col = 10, add = TRUE)
prediction(w_qda$posterior[,2], test$sex) %>%
  performance(measure = "tpr", x.measure = "fpr") %>%
  plot(col = 11, add = TRUE)
prediction(w_lda$posterior[,2], test$sex) %>%
  performance(measure = "tpr", x.measure = "fpr") %>%
  plot(col = 12, add = TRUE)
prediction(pm_qda$posterior[,2], test$sex) %>%
  performance(measure = "tpr", x.measure = "fpr") %>%
  plot(col = 13, add = TRUE)
prediction(pm_lda$posterior[,2], test$sex) %>%
  performance(measure = "tpr", x.measure = "fpr") %>%
  plot(col = 14, add = TRUE)


legend(0.85, 0.8, c(rallqda$auc,
                     ralllda$auc,
                     rt_wqda$auc,
                     rt_wlda$auc,
                     rw_pmqda$auc,
                     rw_pmlda$auc,
                     rt_pmqda$auc,
                     rt_pmlda$auc,
                     rt_qda$auc,
                     rt_lda$auc,
                     rw_qda$auc,
                     rw_lda$auc,
                     rpm_qda$auc,
                     rpm_lda$auc),1:14)


#итак лучшая модель - это lda крыло+2пм+хвост
plot(ldamodelall)

coef(ldamodelall)
#tail     0.03897862
#wing_max 0.32550768
#x2pm     0.19319842

p2_train_all_lda <- predict(ldamodelall, train)$class
tab2_train_all_lda <- table(Predicted = p2_train_all_lda, Actual = train$sex)
tab2_train_all_lda
train_accuracy_all_lda <- sum(diag(tab2_train_all_lda))/sum(tab2_train_all_lda) * 100
train_accuracy_all_lda
#80.0759

#лучшая модель по одному признаку - это lda крыло
#Протестируем ее на полных данных df_na

df_na<-df[!(is.na(df$wing)), ]
df_na<-df_na[!(is.na(df_na$sex)), ]
summary(df_na)

df_na<-filter(df_na, wing > 60)
df_na<-filter(df_na, wing < 75)

samplena <- sample(c(TRUE, FALSE), nrow(df_na), replace=TRUE, prob=c(0.7,0.3))
traina <- df_na[samplena, ]
testna <- df_na[!samplena, ] 

#test predicted
names(df_na)[names(df_na) == 'wing'] <- 'wing_max'

test.predicted.w <- predict(ldamodelt_w, newdata = df_na)
prediction(test.predicted.w$posterior[,2], df_na$sex) %>%
  performance(measure = "auc") %>%
  .@y.values
## auc
## 0.8250813
coef(ldamodelt_w)
plot(ldamodelt_w)
#wing_max 0.4920731
p2_train_all_lda <- predict(ldamodelt_w, data)$class
tab2_train_all_lda <- table(Predicted = p2_train_all_lda, Actual = data$sex)
tab2_train_all_lda
train_accuracy_all_lda <- sum(diag(tab2_train_all_lda))/sum(tab2_train_all_lda) * 100
train_accuracy_all_lda

par(mfrow=c(1,2))
ggplot(test.predicted.w$x, aes(test.predicted.w$x[,1], test.predicted.w$class, color=df_na$sex)) + 
  geom_point() +  theme_minimal() + scale_y_discrete(expand= c(0, 0.3), name = "Предсказанный пол") +
  scale_x_continuous(name = "X - переменная") +
  scale_colour_brewer("Пол", type = "qual", palette = 6) 
prediction(test.predicted.w$posterior[,2], df_na$sex) %>%
  performance(measure = "tpr", x.measure = "fpr") %>%
  plot(main = "Пол=Крыло, AUC = 0.835069", xlab = "Доля ложных положительных оценок", ylab = "Доля правильных положительных оценок")
#
#
#_____________________________________________________________________________________________________



##Дифференциальная миграция

df_na_sex<-df[!(is.na(df$age_sex)), ]
df_na_sex<-df_na_sex[!(is.na(df_na_sex$sex)), ]
df_na_sex<-df_na_sex[!(is.na(df_na_sex$age)), ]
ggplot(df_na_sex, aes(x = `day`, y = `year`, color = age, fill = age)) + geom_density_ridges(size = 0.25) +
  scale_y_discrete(expand= c(0, 0)) +
  scale_x_continuous(expand = c(0, 0), breaks = waiver(), minor_breaks = waiver(), n.breaks = 14, name = "День от начала года") +
  scale_fill_manual(values = c("#D55E0050", "#0072B250"), labels = c("Молодые", "Взрослые"), name = "Возраст") +
  scale_color_manual(values = c("#D55E00", "#0072B2"), guide = "none") +
  
  coord_cartesian(clip = "off") +
  guides(fill = guide_legend(
    override.aes = list(
      fill = c("#D55E00A0", "#0072B2A0"),
      color = NA, point_color = NA)
  )
  ) +
  ggtitle("Диффиренциальная миграция") + labs(y = "Год") +
  theme_ridges(center = TRUE)
#

ggbetweenstats( plot.type = "box", 
                data  = df_na_sex,
                x     = age_sex,
                xlab = "Половозрастная группа",
                ylab = "День от начала года",
                y     = day,
                outlier.shape = NA,
                pairwise.display = "s",
                p.adjust.method  = "hochberg",
                palette          = "default_jama",
                package          = "ggsci",
                plotgrid.args    = list(nrow = 1),
                title = "Половозрастные различия в днях отлова"
)
ggbetweenstats( plot.type = "box", 
                data  = df_na,
                x     = age,
                xlab = "Возрастная группа",
                ylab = "День от начала года",
                y     = day,
                outlier.shape = NA,
                pairwise.display = "s",
                p.adjust.method  = "hochberg",
                palette          = "default_jama",
                package          = "ggsci",
                plotgrid.args    = list(nrow = 1),
                title = "Возрастные различия в днях отлова"
)
ggscatterstats( data  = df1d,
                x     = day,
                xlab = "День",
                ylab = "С3",
                y     = c3,
                type = "parametric",
                conf.level = 0.95,
                marginal = TRUE,
                xfill = "#009E73",
                yfill = "#D55E00",
                point.args = list(size = 3, alpha = 0.4, stroke = 0, na.rm = TRUE),
                point.width.jitter = 0,
                point.height.jitter = 0,
                point.label.args = list(size = 3, max.overlaps = 1e+06),
                smooth.line.args = list(size = 1.5, color = "blue", method = "lm", formula = y ~ x,
                                        na.rm = TRUE), 
)
grouped_ggbetweenstats(
  plot.type = "boxviolin", 
  data  = df_na,
  x     = age,
  xlab = "Возрастная группа",
  ylab = "День от начала года",
  y     = day,
  outlier.shape = NA,
  pairwise.display = "s",
  p.adjust.method  = "hochberg",
  palette          = "default_jama",
  package          = "ggsci",
  grouping.var = year,
  plotgrid.args = list(nrow = 5, ncol = 5)
)
grouped_ggscatterstats(
  data = df1d,
  type = "parametric",
  x = day,
  y = c3,
  xlab = "День от начала года",
  ylab = "С3",
  conf.level = 0.95,
  grouping.var = sex,
  plotgrid.args = list(nrow = 1, ncol = 2),
  xfill = "#009E73",
  yfill = "#D55E00",
  point.args = list(size = 3, alpha = 0.4, stroke = 0, na.rm = TRUE),
  point.width.jitter = 0,
  point.height.jitter = 0,
  point.label.args = list(size = 3, max.overlaps = 1e+06),
  smooth.line.args = list(size = 1.5, color = "blue", method = "lm", formula = y ~ x,
                          na.rm = TRUE)
)
grouped_ggscatterstats(
  data = df_na,
  type = "parametric",
  x = day,
  y = wing,
  xlab = "День от начала года",
  ylab = "Крыло, см",
  conf.level = 0.95,
  grouping.var = sex,
  plotgrid.args = list(nrow = 2, ncol = 1),
  xfill = "#009E73",
  yfill = "#D55E00",
  point.args = list(size = 3, alpha = 0.4, stroke = 0, na.rm = TRUE),
  point.width.jitter = 1,
  point.height.jitter = 1,
  point.label.args = list(size = 3, max.overlaps = 1e+06),
  smooth.line.args = list(size = 1.5, color = "blue", method = "lm", formula = y ~ x,
                          na.rm = TRUE)
)

cleaner.aov <- aov(formula = c3 ~ fat555,
                   data = df1d)
summary(cleaner.aov)
#              Df Sum Sq Mean Sq F value Pr(>F)    
# age_sex        3  19319    6440   46.82 <2e-16 ***
#  Residuals   6720 924232     138                   
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
TukeyHSD(cleaner.aov)
results_test <- TukeyHSD(cleaner.aov, conf.level=.95)
results_matrix <- as.matrix (results_test) 
df_res <- as.data.frame(results_matrix[1])
plot(results_matrix, col= ifelse(df_res[,4]<0.05, 'red', 'black'))



cleaner.aov <- aov(formula = day ~ age_sex,
                   data = df_na_sex)
summary(cleaner.aov)

TukeyHSD(cleaner.aov)
results_test <- TukeyHSD(cleaner.aov, conf.level=.95)
results_matrix <- as.matrix (results_test) 
df_res <- as.data.frame(results_matrix[1])
plot(results_matrix, col= ifelse(df_res[,4]<0.05, 'red', 'black'))

cleaner.aov1 <- aov(formula = day ~ year,
                   data = gt)
summary(cleaner.aov1)

TukeyHSD(cleaner.aov1)
results_test <- TukeyHSD(cleaner.aov1, conf.level=.95)
results_matrix <- as.matrix (results_test) 
df_res <- as.data.frame(results_matrix[1])
plot(results_matrix, col= ifelse(df_res[,4]<0.05, 'red', 'black'))
