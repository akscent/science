#Графический анализ
library(readxl)
морфология <- read_excel("C:/Users/user/Рабочий стол/Moio govno/FNC/Database/Banding-station/Litovka/морфология.xlsx", 
                         sheet = "morph", col_types = c("text", 
                                                        "text", "numeric", "numeric", "text", 
                                                        "text", "text", "text", "numeric", 
                                                        "numeric", "numeric", "numeric", 
                                                        "numeric", "numeric", "numeric", 
                                                        "numeric", "text", "text", "numeric", 
                                                        "numeric", "numeric", "numeric", 
                                                        "numeric", "numeric", "numeric"))

library(ggplot2)
library(ggcorrplot)
library(effectsize)
library(stargazer)
library(sjPlot)
library(MorphoTools2)
library(gtsummary)
library(ggstatsplot)
library(ggsci)
library(effectsize)
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
library(sjPlot)
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
library(plotly)
library(stats)
library(nsprcomp)


data<-морфология
M <- cor(df1d[,2:23], method = "kendall")
corrplot(M)

ggcorrplot(M, hc.order = TRUE, type = "lower",
           colors = c("white","yellow","purple" ), lab = TRUE)

interpret_r(M, rules = "funder2019")
#Далее нужно выделить коррелирующие переменные в отдельный список или пару списков, PCA
pca_rut <- prcomp(cbind(data[,9:16], data[,19:25]), center = TRUE, scale. = TRUE, )

p <- autoplot(pca_rut, data = data, colour = 'sex')
ggbiplot(pca_rut)
ggplotly(p)
biplot(pca_rut)
print(pca_rut)
pairs.panels(pca_rut$x,
             gap=0,
             bg = c("red", "blue")[data$sex],
             pch=21)
autoplot(pca_rut, data = data, colour = 'sex',
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3)


#PCA графики


X <- cbind(data[,9:16], data[,19:25])
prin_comp <- prcomp(X)
explained_variance_ratio <- summary(prin_comp)[["importance"]]['Proportion of Variance',]
explained_variance_ratio <- 100 * explained_variance_ratio
components <- prin_comp[["x"]]
components <- data.frame(components)
components <- cbind(components, data$sex)
components$PC3 <- -components$PC3
components$PC2 <- -components$PC2

axis = list(showline=FALSE,
            zeroline=FALSE,
            gridcolor='#ffff',
            ticklen=4,
            titlefont=list(size=13))

fig <- components %>%
  plot_ly()  %>%
  add_trace(
    type = 'splom',
    dimensions = list(
      list(label=paste('PC 1 (',toString(round(explained_variance_ratio[1],1)),'%)',sep = ''), values=~PC1),
      list(label=paste('PC 2 (',toString(round(explained_variance_ratio[2],1)),'%)',sep = ''), values=~PC2),
      list(label=paste('PC 3 (',toString(round(explained_variance_ratio[3],1)),'%)',sep = ''), values=~PC3),
      list(label=paste('PC 4 (',toString(round(explained_variance_ratio[4],1)),'%)',sep = ''), values=~PC4)
    ),
    color = ~data$sex, colors = c('#636EFA50','#EF553B50')
  ) %>%
  layout(
    legend=list(title=list(text='color')),
    hovermode='closest',
    dragmode= 'select',
    plot_bgcolor='rgba(240,240,240, 0.95)',
    xaxis=list(domain=NULL, showline=F, zeroline=F, gridcolor='#ffff', ticklen=4),
    yaxis=list(domain=NULL, showline=F, zeroline=F, gridcolor='#ffff', ticklen=4),
    xaxis2=axis,
    xaxis3=axis,
    xaxis4=axis,
    yaxis2=axis,
    yaxis3=axis,
    yaxis4=axis
  )

fig

#EFA
fa.parallel(M,n.obs=112,fa="both",n.iter=100, main = "Screeplots with parrallel analysis")
fa<-fa(M,nfactors=2,rotate="none",fm="pa")
fa
fa.promax<-fa(M,nfactors=2,rotate="promax",fm="pa")
fa.promax
fa.diagram(fa.promax,simple=FALSE, main = "Факторные нагрузки")

# Возможно это стоит убрать
wing<- df %>%
  group_by(age_sex) %>% 
  summarize(
    q_05 = quantile(wing, 0.05, na.rm = TRUE),
    median = median(wing, na.rm = TRUE),
    q_95 = quantile(wing, 0.95, na.rm = TRUE), mean = mean(wing, na.rm = TRUE),
    SD = sd(wing, na.rm = TRUE), Cv = sd(wing, na.rm = TRUE)/mean(wing, na.rm = TRUE))




ggbetweenstats(
  data  = data,
  x     = H_sex,
  xlab = "Половозрастная группа",
  ylab = "Длина, см",
  y     = Alula,
  outlier.shape = NA,
  pairwise.display = "s",
  p.adjust.method  = "hochberg",
  palette          = "default_jama",
  package          = "ggsci",
  plotgrid.args    = list(nrow = 1),
  title = "Половозрастные различия в ALULA, см"
)
ggbetweenstats(
  data  = data,
  x     = H_sex,
  xlab = "Половозрастная группа",
  ylab = "Длина, см",
  y     = "1pm",
  outlier.shape = NA,
  pairwise.display = "s",
  p.adjust.method  = "hochberg",
  palette          = "default_jama",
  package          = "ggsci",
  plotgrid.args    = list(nrow = 1),
  title = "Половозрастные различия в длине 1 ПМ, см"
)
ggbetweenstats(
  data  = data,
  x     = H_sex,
  xlab = "Половозрастная группа",
  ylab = "Длина, см",
  y     = '2pm',
  outlier.shape = NA,
  pairwise.display = "s",
  p.adjust.method  = "hochberg",
  palette          = "default_jama",
  package          = "ggsci",
  plotgrid.args    = list(nrow = 1),
  title = "Половозрастные различия в длине 2 ПМ, см"
)
ggbetweenstats(
  data  = data,
  x     = H_sex,
  xlab = "Половозрастная группа",
  ylab = "Длина, см",
  y     = '3pm',
  outlier.shape = NA,
  pairwise.display = "s",
  p.adjust.method  = "hochberg",
  palette          = "default_jama",
  package          = "ggsci",
  plotgrid.args    = list(nrow = 1),
  title = "Половозрастные различия в длине 3 ПМ, см"
)
ggbetweenstats(
  data  = data,
  x     = H_sex,
  xlab = "Половозрастная группа",
  ylab = "Длина, см",
  y     = '4pm',
  outlier.shape = NA,
  pairwise.display = "s",
  p.adjust.method  = "hochberg",
  palette          = "default_jama",
  package          = "ggsci",
  plotgrid.args    = list(nrow = 1),
  title = "Половозрастные различия в длине 4 ПМ, см"
)
ggbetweenstats(
  data  = data,
  x     = H_sex,
  xlab = "Половозрастная группа",
  ylab = "Длина, см",
  y     = '5pm',
  outlier.shape = NA,
  pairwise.display = "s",
  p.adjust.method  = "hochberg",
  palette          = "default_jama",
  package          = "ggsci",
  plotgrid.args    = list(nrow = 1),
  title = "Половозрастные различия в длине 5 ПМ, см"
)
ggbetweenstats(
  data  = data,
  x     = H_sex,
  xlab = "Половозрастная группа",
  ylab = "Длина, см",
  y     = '6pm',
  outlier.shape = NA,
  pairwise.display = "s",
  p.adjust.method  = "hochberg",
  palette          = "default_jama",
  package          = "ggsci",
  plotgrid.args    = list(nrow = 1),
  title = "Половозрастные различия в длине 6 ПМ, см"
)
ggbetweenstats(
  data  = data,
  x     = H_sex,
  xlab = "Половозрастная группа",
  ylab = "Длина, см",
  y     = wing_min,
  outlier.shape = NA,
  pairwise.display = "s",
  p.adjust.method  = "hochberg",
  palette          = "default_jama",
  package          = "ggsci",
  plotgrid.args    = list(nrow = 1),
  title = "Половозрастные различия в мин. длине крыла, см"
)
ggbetweenstats(
  data  = gt,
  x     = year,
  xlab = "Год",
  ylab = "День от начала года",
  y     = day,
  outlier.shape = NA,
  pairwise.display = "s",
  p.adjust.method  = "hochberg",
  palette          = "default_jama",
  package          = "ggsci",
  plotgrid.args    = list(nrow = 1),
  title = "Распределение дней отловов от года"
) + geom_point()


#Создание таблиц по возрасту между мужчинами и женщинами
data1<-`морфология`
df_ahy<-filter(data1, H == "AHY")
df_hy<-filter(data1, H == "HY")

add_ES_hy <- df_data %>% select(sex, head, tail, tarsus, weight, wing_min, beak_min, beak_max, Al, `al<primary_coverts`, `1pm`, `2pm`, `3pm`, `4pm`, `5pm`, `6pm`)%>% 
  tbl_summary(by = sex, 
              label = list(head ~ "Голова, см", tail ~ "Хвост, см", tarsus ~ "Цевка, см", weight ~ "Вес, г", wing_min ~ "Крыло, см", beak_min ~"Клюв от ноздри, см",  `beak_max` ~ "Клюв от лба, см", `Al` ~ "Алула, см", `al<primary_coverts` ~ "Алула<БВКПМ, см" `1pm` ~ "1 ПМ, см" , `2pm` ~ "2 ПМ, см" , `3pm` ~ "3 ПМ, см" , `4pm` ~ "4 ПМ, см" , `5pm` ~ "5 ПМ, см" , `6pm` ~ "6 ПМ, см" ),
              statistic = list(all_continuous() ~ "{mean} ({sd})"),  missing = "no")  %>% add_difference(test = list(all_continuous()~"t.test", all_continuous()~"cohens_d")) %>% add_p(test = everything() ~ t.test) %>% bold_p() %>% 
              modify_header (label = "**Характеристика**") 
add_ES_hy
add_ES_ahy <- df_ahy %>% select(sex, head_length, tail, tarsus, weight, wing_min, wing_max, beak_min, beak_max, Alula, `1pm`, `2pm`, `3pm`, `4pm`, `5pm`, `6pm`)%>% 
  tbl_summary(by = sex, 
              label = list(head_length ~ "Голова, см", tail ~ "Хвост, см", tarsus ~ "Цевка, см", weight ~ "Вес, г", wing_min ~ "Расправленное крыло, см", wing_min ~ "Крыло, см", beak_min ~"Клюв от ноздри, см",  beak_max ~ "Клюв от лба, см", Alula ~ "АЛУЛА, см", `1pm` ~ "1 ПМ, см" , `2pm` ~ "2 ПМ, см" , `3pm` ~ "3 ПМ, см" , `4pm` ~ "4 ПМ, см" , `5pm` ~ "5 ПМ, см" , `6pm` ~ "6 ПМ, см" ),
              statistic = list(all_continuous() ~ "{mean} ({sd})"),  missing = "no")  %>% add_difference(test = list(all_continuous()~"t.test", all_continuous()~"cohens_d")) %>% add_p(test = everything() ~ t.test) %>% bold_p() %>% 
  modify_header (label = "**Характеристика**") 
tbl1<-tbl_merge(tbls = list(add_ES_hy, add_ES_ahy), tab_spanner = c("**Первый год жизни**", "**Второй и более год жизни**"))

tbl1 %>%    # build gtsummary table
  as_gt() %>%             # convert to gt table
  gt::gtsave(             # save table as image
    filename = "my_table_image.png"
  )

#Наибольшие различия в крыле, 2 пм и хвосте. При этом крыло между ними высокая корреляция

box<-boxM(cbind(data[,13], data[,15], data[,21]), data$sex) 
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

cleaner.aov1 <- aov(formula = day ~ year,
                   data = gt)
summary(cleaner.aov1)
#              Df Sum Sq Mean Sq F value Pr(>F)    
# age_sex        3  19319    6440   46.82 <2e-16 ***
#  Residuals   6720 924232     138                   
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
TukeyHSD(cleaner.aov1)
results_test <- TukeyHSD(cleaner.aov1, conf.level=.95)
results_matrix <- as.matrix (results_test) 
df_res <- as.data.frame(results_matrix[1])
plot(results_matrix, col= ifelse(df_res[,4]<0.05, 'red', 'black'))
