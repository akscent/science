library(nscancor)
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
data<-морфология

library(tidyverse)
library(dplyr)
library(janitor)
library(Seurat)

require("plyr")
require("abind")
library(nsprcomp)
trans.spca <- spca(X, ncomp = 10, center = TRUE, scale = TRUE)

data1<-filter(data1, x1pm>14)
data2 <- data1 %>%
  as_tibble() %>%
  mutate(log(data1[, 17:23]))


C <- daply(cbind(data2[, 5], data2[, 17:23]), "sex", function(x) cov(x[, -ncol(x)]))

C <- aperm(C, c(2, 3, 1)) # put the 1st dimension to the end
dim(C)
dimnames(C)

mod <- cpc(C)
str(mod)


scpca_sim <- scPCA(target = data1[, 17:23],
                   background = data2[, 17:23],
                   n_centers = 4, cv = 3, 
                   penalties = exp(seq(log(0.01), log(0.5), length.out = 10)),
                   alg = "var_proj")
scpca_d <- scpca_sim$x %>%
  as_tibble()
colnames(scpca_d) <- c("scPC1", "scPC2", "label")
p_scpca <- ggplot(scpca_d, aes(x = scPC1, y = scPC2, colour = data1$h_sex)) +
  geom_point(alpha = 0.5) +
  ggtitle("scPCA of Simulated Data") +
  theme_minimal()
p_scpca

RunSPCA(data2[, 17:23],
        assay = NULL,
        npcs = 3,
        reduction.key = "SPC_",
        graph = NULL,
        verbose = FALSE,
        seed.use = 42)
