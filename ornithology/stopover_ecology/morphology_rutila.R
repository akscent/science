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

df<-морфология
df_na<-df[!(is.na(df$head)), ]
df_na<-df_na[!(is.na(df_na$tail)), ]
df_na<-df_na[!(is.na(df_na$tarsus)), ]
df_na<-df_na[!(is.na(df_na$height)), ]
df_na<-df_na[!(is.na(df_na$wing)), ]
df_na<-df_na[!(is.na(df_na$beak)), ]
df_na<-df_na[!(is.na(df_na$age_sex)), ]

M <- cor(df_na[, 7:12])
ggcorrplot(M, hc.order = TRUE, type = "lower",
           colors = c("white","yellow","purple" ), lab = TRUE)

interpret_r(M, rules = "funder2019")
wing<- df %>%
  group_by(age_sex) %>% 
  summarize(
    q_05 = quantile(wing, 0.05, na.rm = TRUE),
    median = median(wing, na.rm = TRUE),
    q_95 = quantile(wing, 0.95, na.rm = TRUE), mean = mean(wing, na.rm = TRUE),
    SD = sd(wing, na.rm = TRUE), Cv = sd(wing, na.rm = TRUE)/mean(wing, na.rm = TRUE))
manovamdl<-manova(cbind(head, tail, tarsus, height, wing, beak) ~ age_sex, data = df_na)
sum_manova<-summary(manovamdl)
summary(manovamdl) #Есть статистически значимая разница
sum_aov<-summary.aov(manovamdl)
summary.aov(manovamdl) #Есть статистически значимое влияние на все, кроме клюва

ggbetweenstats(
  data  = df,
  x     = age_sex,
  xlab = "Половозрастная группа",
  ylab = "Длина, см",
  y     = wing,
  outlier.shape = NA,
  pairwise.display = "s",
  p.adjust.method  = "hochberg",
  palette          = "default_jama",
  package          = "ggsci",
  plotgrid.args    = list(nrow = 1),
  title = "Половозрастные различия в длине крыла, см"
)
ggbetweenstats(
  data  = df,
  x     = age_sex,
  xlab = "Половозрастная группа",
  ylab = "Длина, см",
  y     = head,
  outlier.shape = NA,
  pairwise.display = "s",
  p.adjust.method  = "hochberg",
  palette          = "default_jama",
  package          = "ggsci",
  plotgrid.args    = list(nrow = 1),
  title = "Половозрастные различия в длине головы, см"
)
ggbetweenstats(
  data  = df,
  x     = age_sex,
  xlab = "Половозрастная группа",
  ylab = "Длина, см",
  y     = tail,
  outlier.shape = NA,
  pairwise.display = "s",
  p.adjust.method  = "hochberg",
  palette          = "default_jama",
  package          = "ggsci",
  plotgrid.args    = list(nrow = 1),
  title = "Половозрастные различия в длине хвоста, см"
)
ggbetweenstats(
  data  = df,
  x     = age_sex,
  xlab = "Половозрастная группа",
  ylab = "Длина, см",
  y     = tarsus,
  outlier.shape = NA,
  pairwise.display = "s",
  p.adjust.method  = "hochberg",
  palette          = "default_jama",
  package          = "ggsci",
  plotgrid.args    = list(nrow = 1),
  title = "Половозрастные различия в длине цевки, см"
)
ggbetweenstats(
  data  = df,
  x     = age_sex,
  xlab = "Половозрастная группа",
  ylab = "Длина, см",
  y     = beak,
  outlier.shape = NA,
  pairwise.display = "s",
  p.adjust.method  = "hochberg",
  palette          = "default_jama",
  package          = "ggsci",
  plotgrid.args    = list(nrow = 1),
  title = "Половозрастные различия в длине клюва, см"
)
gghistostats(
  data       = df,
  x          = wing,
  title      = "Длина крыла, см",
  test.value = 12,
  binwidth   = 1
)
gghistostats(
  data       = df,
  x          = tail,
  title      = "Длина хвоста, см",
  test.value = 12,
  binwidth   = 1
)
gghistostats(
  data       = df,
  x          = tarsus,
  title      = "Длина цевки, см",
  test.value = 12,
  binwidth   = 1
)
gghistostats(
  data       = df,
  x          = head,
  title      = "Длина головы, см",
  test.value = 12,
  binwidth   = 1
)
gghistostats(
  data       = df,
  x          = beak,
  title      = "Длина клюва, см",
  test.value = 12,
  binwidth   = 1
)
t.test(wing ~ sex, data = df) #тест стьюдента на различия

cohens_d(wing ~ sex, data = df, pooled_sd = FALSE) #Размер эффекта Коэна wing
interpret_cohens_d(d = 0.96, rules = "cohen1988")
#@Article{, title = {{e}ffectsize: Estimation of Effect Size Indices and Standardized Parameters}, author = {Mattan S. Ben-Shachar and Daniel Lüdecke and Dominique Makowski},
  #year = {2020},
  #journal = {Journal of Open Source Software},
  #volume = {5},
  #number = {56},
  #pages = {2815},
  #publisher = {The Open Journal},
  #doi = {10.21105/joss.02815},
  #url = {https://doi.org/10.21105/joss.02815}


pca_rut <- prcomp(df_na[7:12], scale. = TRUE)

p <- autoplot(pca_rut, data = df_na, colour = 'age')

ggplotly(p)


pc_na<-mutate(df_na, as.data.frame(pca_rut$x))

ggplot(pc_na, aes(PC1, PC2, fill = age_sex)) + coord_equal() +
  geom_hdr()
ggplot(pc_na, aes(PC1, PC2, fill = age_sex)) + coord_equal() +
  geom_hdr() +  facet_wrap(vars(age_sex))

hy<-filter(pc_na, age == "HY")
ahy<-filter(pc_na, age == "AHY")

ggplotly(p)


pc_na<-mutate(df_na, as.data.frame(pca_rut$x))

ggplot(pc_na, aes(PC1, PC2, fill = age_sex)) + coord_equal() +
  geom_hdr()
ggplot(pc_na, aes(PC1, PC2, fill = age_sex)) + coord_equal() +
  geom_hdr() +  facet_wrap(vars(age_sex))


d <- df_na %>%
  ggplot( aes(x=day, fill=age)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("#69b3a2", "#404080")) +
  theme_ipsum() +
  labs(x = "День в году", y = "Плотность", fill="Возраст")



ggplot(df_na) +
 aes(x = day, fill = age_sex) +
 geom_histogram(bins = 80L) +
 scale_fill_hue(direction = 1) +
 labs(x = "День в году", y = "Количество первичных отловов, шт", 
 title = "Рыжая овсянка", subtitle = "Ocyris rutilus") +
 theme_minimal() +
 theme(legend.position = "bottom") +
 facet_wrap(vars(year))
ggplot(df_na) + aes(day, fill = sex) + geom_density(position = "stack", color="#e9ecef", alpha=0.6) +  scale_fill_manual(values=c("#69b3a2", "#404080")) + theme_ipsum() + labs(x = "День в году", y = "Плотность", fill="")
ggplot(df_na) + aes(day, fill = age) + geom_density(position = "stack", color="#e9ecef", alpha=0.6) +  scale_fill_manual(values=c("#69b3a2", "#404080")) + theme_ipsum() + labs(x = "День в году", y = "Плотность", fill="")
g <- df_na %>%
  ggplot( aes(x=day, fill=sex)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("#69b3a2", "#404080")) +
  theme_ipsum() +
  labs(x = "День в году", y = "Плотность", fill="Пол")

#Создание таблиц по возрасту между мужчинами и женщинами
df_ahy<-filter(df_hy, age == "AHY")
df_hy<-filter(df_hy, age == "HY")
#Без NA
df_ahy_na<-filter(df_na, age == "AHY")
df_hy_na<-filter(df_na, age == "HY")

my_EStest <- function(data, variable, by, ...) {
  # Cohen's D
  d <- effsize::cohen.d(data[[variable]] ~ as.factor(data[[by]]), 
                        conf.level=.90, pooled=TRUE, paired=FALSE, 
                        hedges.correction=TRUE) 
  # Formatting statistic with CI
  est <- style_sigfig(d$estimate)
  ci <- style_sigfig(d$conf.int) %>% paste(collapse = ", ")
  
  # returning estimate with CI together
  str_glue("{est} ({ci})")
}
add_ES_hy <- df_hy %>% select(sex, head, tail, tarsus, height, wing, beak)%>% 
  tbl_summary(by = sex, 
              label = list(head ~ "Голова, см", tail ~ "Хвост, см", tarsus ~ "Цевка, см", height ~ "Вес, г", wing ~ "Крыло, см", beak ~ "Клюв, см"),
              statistic = list(all_continuous() ~ "{mean} ({sd})"),  missing = "no")  %>% add_n()  %>% add_p(test = everything() ~ t.test) %>% bold_p() %>%
  add_stat(
    fns = everything() ~ my_EStest,
    fmt_fun = NULL) %>% 
    modify_header (list(
      add_stat_1 ~ "**ES (90% CI)**", 
      all_stat_cols() ~ "**{level}**"
    ), label = "**Характеристика**", p.value = "**P**"
  ) %>%
  modify_footnote( add_stat_1 ~ "Cohen's D (90% CI)")
add_ES_hy
add_ES_ahy <- df_ahy %>% select(sex, head, tail, tarsus, height, wing, beak) %>% 
  tbl_summary(by = sex, 
              label = list(head ~ "Голова, см", tail ~ "Хвост, см", tarsus ~ "Цевка, см", height ~ "Вес, г", wing ~ "Крыло, см", beak ~ "Клюв, см"),
              statistic = list(all_continuous() ~ "{mean} ({sd})"),  missing = "no")  %>% add_n()  %>% add_p(test = everything() ~ t.test) %>% bold_p() %>%
  add_stat(
    fns = everything() ~ my_EStest,
    fmt_fun = NULL) %>% 
  modify_header (list(
    add_stat_1 ~ "**ES (90% CI)**",
    all_stat_cols() ~ "**{level}**"
  ), label = "**Характеристика**", p.value = "**P**"
  ) %>%
  modify_footnote(add_stat_1 ~ "Cohen's D (90% CI)")
add_ES_ahy
tbl1<-tbl_merge(tbls = list(add_ES_hy, add_ES_ahy), tab_spanner = c("**Первый год жизни**", "**Второй и более год жизни**"))

tbl1

BoxM(df_na[7:12], df_na$sex) #результаты указывают на то, что ковариационные матрицы отличаются, т.е. нужно применять QDA
#нарушение допущения однородности дисперсионно-ковариационных матриц
#PCA анализ

pca_res<-prcomp(df_na[7:12])
autoplot(pca_res, data = df_na, colour = 'age_sex')
#ничего не дал
#QDA
#ahy
sample <- sample(c(TRUE, FALSE), nrow(df_ahy_na), replace=TRUE, prob=c(0.7,0.3))
trainahy <- df_ahy_na[sample, ]
testahy <- df_ahy_na[!sample, ] 
model1<-qda(sex~tail+wing, data = trainahy)
model2<-lda(sex~tail+wing, data = trainahy)
plot(model2)
partimat(sex~tail+wing, data = trainahy, method = 'lda')
partimat(sex~tail+wing, data = trainahy, method = 'qda')
p2_train_all <- predict(model1, trainahy)$class
tab2_train_all <- table(Predicted = p2_train_all, Actual = trainahy$sex)
tab2_train_all
# Actual
#Predicted  F  M
         #F 21 13
         #M 11 26
qda_train_accuracy_all <- sum(diag(tab2_train_all))/sum(tab2_train_all) * 100
qda_train_accuracy_all
#66,1972%

#только крыло
model1w<-qda(sex~wing, data = trainahy)
model2w<-lda(sex~wing, data = trainahy)
plot(model2w)
p2_train_allw <- predict(model1w, trainahy)$class
tab2_train_allw <- table(Predicted = p2_train_allw, Actual = trainahy$sex)
tab2_train_allw
qda_train_accuracy_allw <- sum(diag(tab2_train_allw))/sum(tab2_train_allw) * 100
qda_train_accuracy_allw
#70.4225%
#все переменные
model1all<-qda(sex~tail+wing+tarsus+head+beak+height, data = trainahy)
model2all<-lda(sex~tail+wing+tarsus+head+beak+height, data = trainahy)
plot(model2all)
partimat(sex~tail+wing+tarsus+head+beak+height, data = trainahy, method = 'qda')
p2_train_all <- predict(model1all, trainahy)$class
tab2_train_all <- table(Predicted = p2_train_all, Actual = trainahy$sex)
tab2_train_all
qda_train_accuracy_all <- sum(diag(tab2_train_all))/sum(tab2_train_all) * 100
qda_train_accuracy_all
#74.6479
#test predicted
test.predicted.qda_ahy1 <- predict(model1, newdata = testahy)
test.predicted.qda_ahy2 <- predict(model1w, newdata = testahy)
test.predicted.qda_ahy3 <- predict(model1all, newdata = testahy)
par(mfrow=c(1, 3))
prediction(test.predicted.qda_ahy1$posterior[,2], testahy$sex) %>%
  performance(measure = "tpr", x.measure = "fpr") %>%
  plot(main = "Пол=Крыло+Хвост", xlab = "Доля ложных положительных оценок", ylab = "Доля правильных положительных оценок")

prediction(test.predicted.qda_ahy2$posterior[,2], testahy$sex) %>%
  performance(measure = "tpr", x.measure = "fpr") %>%
  plot(main = "Пол=Крыло", xlab = "Доля ложных положительных оценок", ylab = "Доля правильных положительных оценок")
prediction(test.predicted.qda_ahy3$posterior[,2], testahy$sex) %>%
  performance(measure = "tpr", x.measure = "fpr") %>%
  plot(main = "Пол=Крыло+Хвост+Цевка+Голова+Клюв+Вес", xlab = "Доля ложных положительных оценок", ylab = "Доля правильных положительных оценок")
# model 1 w+t
prediction(test.predicted.qda_ahy1$posterior[,2], testahy$sex) %>%
  performance(measure = "auc") %>%
  .@y.values
## [[1]]
## 0.786458

# model 2 w
prediction(test.predicted.qda_ahy2$posterior[,2], testahy$sex) %>%
  performance(measure = "auc") %>%
  .@y.values
## [[1]]
## 0.835069

#model 3 all
prediction(test.predicted.qda_ahy3$posterior[,2], testahy$sex) %>%
  performance(measure = "auc") %>%
  .@y.values
## 0.552083


#hy
samplehy <- sample(c(TRUE, FALSE), nrow(df_hy_na), replace=TRUE, prob=c(0.7,0.3))
trainhy <- df_hy_na[samplehy, ]
testhy <- df_hy_na[!samplehy, ] 
model1hy<-qda(sex~tail+wing, data = trainhy)
model2hy<-lda(sex~tail+wing, data = trainhy)
plot(model2hy)
partimat(sex~tail+wing, data = trainhy, method = 'lda')
partimat(sex~tail+wing, data = trainhy, method = 'qda')
p2_train_all <- predict(model1, trainhy)$class
tab2_train_all <- table(Predicted = p2_train_all, Actual = trainhy$sex)
tab2_train_all
qda_train_accuracy_all <- sum(diag(tab2_train_all))/sum(tab2_train_all) * 100
qda_train_accuracy_all
#75,9459%

#только крыло
model1hyw<-qda(sex~wing, data = trainhy)
model2hyw<-lda(sex~wing, data = trainhy)
plot(model2hyw)
p2_train_allw <- predict(model1hyw, trainhy)$class
tab2_train_allw <- table(Predicted = p2_train_allw, Actual = trainhy$sex)
tab2_train_allw
qda_train_accuracy_allw <- sum(diag(tab2_train_allw))/sum(tab2_train_allw) * 100
qda_train_accuracy_allw
#83.2432%
#все переменные
model1hyall<-qda(sex~tail+wing+tarsus+head+beak+height, data = trainhy)
model2hyall<-lda(sex~tail+wing+tarsus+head+beak+height, data = trainhy)
plot(model2hyall)
partimat(sex~tail+wing+tarsus+head+beak+height, data = trainhy, method = 'qda')
p2_train_all <- predict(model1all, trainhy)$class
tab2_train_all <- table(Predicted = p2_train_all, Actual = trainhy$sex)
tab2_train_all
qda_train_accuracy_all <- sum(diag(tab2_train_all))/sum(tab2_train_all) * 100
qda_train_accuracy_all
#71.0811
#AUC for hy
test.predicted.qda_ahy01 <- predict(model1hy, newdata = testhy)
test.predicted.qda_ahy02 <- predict(model1hyw, newdata = testhy)
test.predicted.qda_ahy03 <- predict(model1hyall, newdata = testhy)
par(mfrow=c(1, 3))
prediction(test.predicted.qda_ahy01$posterior[,2], testhy$sex) %>%
  performance(measure = "tpr", x.measure = "fpr") %>%
  plot(main = "Пол = Крыло + Хвост", xlab = "Доля ложных положительных оценок", ylab = "Доля правильных положительных оценок")

prediction(test.predicted.qda_ahy02$posterior[,2], testhy$sex) %>%
  performance(measure = "tpr", x.measure = "fpr") %>%
  plot(main = "Пол = Хвост", xlab = "Доля ложных положительных оценок", ylab = "Доля правильных положительных оценок")
prediction(test.predicted.qda_ahy03$posterior[,2], testhy$sex) %>%
  performance(measure = "tpr", x.measure = "fpr") %>%
  plot(main = "Пол=Крыло+Хвост+Цевка+Голова+Клюв+Вес", xlab = "Доля ложных положительных оценок", ylab = "Доля правильных положительных оценок")
# model 1 w+t
prediction(test.predicted.qda_ahy01$posterior[,2], testhy$sex) %>%
  performance(measure = "auc") %>%
  .@y.values
## [[1]]
## 0.904653

# model 2 w
prediction(test.predicted.qda_ahy02$posterior[,2], testhy$sex) %>%
  performance(measure = "auc") %>%
  .@y.values
## [[1]]
## 0.901725

#model 3 all
prediction(test.predicted.qda_ahy03$posterior[,2], testhy$sex) %>%
  performance(measure = "auc") %>%
  .@y.values
## 0.909163

##Результирующая модель для всей выборки по полу
sample <- sample(c(TRUE, FALSE), nrow(df_na), replace=TRUE, prob=c(0.7,0.3))
train <- df_na[sample, ]
test <- df_na[!sample, ] 
model<-qda(sex~wing, data = train)
plot(model)
p2_train_all <- predict(model, train)$class
tab2_train_all <- table(Predicted = p2_train_all, Actual = train$sex)
tab2_train_all
qda_train_accuracy_all <- sum(diag(tab2_train_all))/sum(tab2_train_all) * 100
qda_train_accuracy_all
# 82.1918%
test.predicted.qda <- predict(model, newdata = test)
test.predicted.qda
qda.cm <- table(test$sex, test.predicted.qda$class)
list(QDA_model = qda.cm %>% prop.table() %>% round(3))
qda.pred.adj <- ifelse(test.predicted.qda$posterior[, 2] > .20, "Yes", "No")
list(QDA_model = table(test$sex, qda.pred.adj))
#   No Yes
# F 53  56
# M  1  87

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
