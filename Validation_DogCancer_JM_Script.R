
##Name: Joao Marcos G. Barbosa         
##ORCID address: https://orcid.org/0000-0003-3560-5838



                        ###############################################################################

####Variable Selection - Cancer discriminaton####
##Data visualization using the 18VOMs (previously selected by GA-PLS)

 
df <- read.csv("GA_DataSet_2023.csv", sep=",",header = TRUE)
Num <- df[,8:25]
Num.mat <- as.matrix(Num)                        
Cancer <- as.factor(df$Cancer)
Sex <- as.factor(df$Sex)
Age <- as.factor(df$Age)
Therapy <- as.factor(df$Therapy)
Cancer_Type <- as.factor(df$Cancer_Type)


### HCA Ward Hamming
library(e1071)
library(ape)
dist.ham<-as.dist(hamming.distance(Num.mat))
groupCodes <- as.character(df$Cancer)
rownames(Num.mat) <- make.unique(groupCodes)
colorCodes <- c(Y="Red", N="Blue")
hc <- hclust(dist.ham, method = "ward.D")
d <- as.dendrogram(hc)
require(dendextend)
labels_colors(d) <- colorCodes[groupCodes][order.dendrogram(d)]
par(cex=1.0)
plot(d)


d %>% set("branches_k_color", value = c("red", "blue"), k = 2) %>% plot()  

colors <- c("blue", "red")
clus2 <- cutree(d,2)
plot(as.phylo(d), type ="fan", tip.color = colors [clus2],
     label.offset = 1, cex = 0.7)




library(e1071) #hamming distance
dist.ham<-as.dist(hamming.distance(Num.mat))
require(ade4) # multivariate analysis pco
mds <- dudi.pco(dist.ham, scannf=F)
require(ggplot2) # fancy plotting
ggplot(data.frame(mds$li, df),
       aes(x=A1,y=A2,col=Cancer, shape = Therapy))+
  geom_point(size=4, alpha=1) +
  scale_color_manual(values=c(Y="Red", N="Blue"))+
  geom_hline(yintercept=0, col="darkgrey")+ 
  geom_vline(xintercept=0, col="darkgrey")+
  scale_y_continuous(name='PCO2')+
  scale_x_continuous(name='PCO1')+
  coord_fixed() +
  theme_bw() +
  theme(text=element_text(size = 16, family="serif"))

ggsave(filename = "PCO_18.png",width = 10,height = 5, dpi = 300, limitsize = F)




###Validation of the 18VOMs selected by GA-PLS (Chromosome 1)
library(caret)
library(pls)
df <- read.csv("Chromosome1.csv", sep=",",header = TRUE) 
Num <- df[,7:25]
Num$Cancer <- factor(Num$Cancer)
Num.mat <- as.matrix(Num)


ind <- sample(2, nrow(Num), replace = T, prob = c(0.8, 0.2))
train <- Num[ind == 1,]
test <- Num[ind == 2,]
str(test)

set.seed(1)
ctrl <- trainControl(method = "repeatedcv",
                     number = 10,
                     repeats = 10,
                     classProbs = TRUE,
                     summaryFunction = twoClassSummary)

set.seed(2)
plsFit <- train(
  Cancer ~ .,
  data = train,
  method = "pls",
  trControl = ctrl,
  metric = "ROC",
  tuneLength = 5,
  tuneGrid = expand.grid(ncomp = 1:5))

plsFit
varImp(plsFit)

ggplot(plsFit)+ theme_bw()+theme(text=element_text(size = 16, family="serif"))

ggsave("Chromosome_1.tiff")

p1 <- predict(plsFit, train)
confusionMatrix(p1, train$Cancer)

p2 <- predict(plsFit, test)
confusionMatrix(p2, test$Cancer)



### GA-PLS Chromosome 2
#Use the Binary dataset file containing the 128VOMs detected in the study 
df <- read.csv("DataSet_dogs_50CA.csv", sep=",",header = TRUE) 
Cancer <- as.factor(df$Cancer)
Num <- df[,8:135]
Num.mat <- as.matrix(Num)
Mx <- df[,7:135]

library("plsVarSel")
library('caret')
set.seed(4)
X <- Num.mat
y <- Cancer
xga <- ga_pls(y, X, GA.threshold = 10, iters = 10, popSize = 100)
Xga <- Mx[,xga$ga.selection]  
Xga$Cancer <- Mx$Cancer

Xga$Cancer <- factor(Xga$Cancer)


ind <- sample(2, nrow(Xga), replace = T, prob = c(0.8, 0.2))
train <- Xga[ind == 1,]
test <- Xga[ind == 2,]
str(test)


set.seed(3)
ctrl <- trainControl(method = "repeatedcv",
                     number = 10,
                     repeats = 10,
                     classProbs = TRUE,
                     summaryFunction = twoClassSummary)

set.seed(5)
plsFit <- train(
  Cancer ~ .,
  data = train,
  method = "pls",
  trControl = ctrl,
  metric = "ROC",
  tuneLength = 5,
  tuneGrid = expand.grid(ncomp = 1:5))

plsFit
varImp(plsFit)

ggplot(plsFit)+ theme_bw()+theme(text=element_text(size = 16, family="serif"))

ggsave("Chromosome_2.tiff")


p1 <- predict(plsFit, train)
confusionMatrix(p1, train$Cancer)

p2 <- predict(plsFit, test)
confusionMatrix(p2, test$Cancer)

#####GA-PLS Chromosome 3
set.seed(6)
X <- Num.mat
y <- Cancer
xga <- ga_pls(y, X, GA.threshold = 10, iters = 100, popSize = 100)
Xga <- Mx[,xga$ga.selection]  
Xga$Cancer <- Mx$Cancer
Xga$Cancer <- factor(Xga$Cancer)


ind <- sample(2, nrow(Xga), replace = T, prob = c(0.8, 0.2))
train <- Xga[ind == 1,]
test <- Xga[ind == 2,]
str(test)


set.seed(7)
plsFit <- train(
  Cancer ~ .,
  data = train,
  method = "pls",
  trControl = ctrl,
  metric = "ROC",
  tuneLength = 5,
  tuneGrid = expand.grid(ncomp = 1:5))

plsFit
varImp(plsFit)

ggplot(plsFit)+ theme_bw()+theme(text=element_text(size = 16, family="serif"))

ggsave("Chromosome_3.png", dpi = 300, limitsize = F)


p1 <- predict(plsFit, train)
confusionMatrix(p1, train$Cancer)

p2 <- predict(plsFit, test)
confusionMatrix(p2, test$Cancer)

#####GA-PLS Chromosome 4
set.seed(8)
X <- Num.mat
y <- Cancer
xga <- ga_pls(y, X, GA.threshold = 10, iters = 100, popSize = 100)
Xga <- Mx[,xga$ga.selection]  
Xga$Cancer <- Mx$Cancer
Xga$Cancer <- factor(Xga$Cancer)


ind <- sample(2, nrow(Xga), replace = T, prob = c(0.8, 0.2))
train <- Xga[ind == 1,]
test <- Xga[ind == 2,]
str(test)


set.seed(9)
plsFit <- train(
  Cancer ~ .,
  data = train,
  method = "pls",
  trControl = ctrl,
  metric = "ROC",
  tuneLength = 5,
  tuneGrid = expand.grid(ncomp = 1:5))

plsFit
varImp(plsFit)

ggplot(plsFit)+ theme_bw()+theme(text=element_text(size = 16, family="serif"))

ggsave("Chromosome_4.png", dpi = 300, limitsize = F)


p1 <- predict(plsFit, train)
confusionMatrix(p1, train$Cancer)

p2 <- predict(plsFit, test)
confusionMatrix(p2, test$Cancer)


####GA-PLS Chromosome 5
set.seed(10)
X <- Num.mat
y <- Cancer
xga <- ga_pls(y, X, GA.threshold = 10, iters = 100, popSize = 100)
Xga <- Mx[,xga$ga.selection]  
Xga$Cancer <- Mx$Cancer
Xga$Cancer <- factor(Xga$Cancer)


ind <- sample(2, nrow(Xga), replace = T, prob = c(0.8, 0.2))
train <- Xga[ind == 1,]
test <- Xga[ind == 2,]
str(test)


set.seed(11)
plsFit <- train(
  Cancer ~ .,
  data = train,
  method = "pls",
  trControl = ctrl,
  metric = "ROC",
  tuneLength = 5,
  tuneGrid = expand.grid(ncomp = 1:5))

plsFit
varImp(plsFit)

ggplot(plsFit)+ theme_bw()+theme(text=element_text(size = 16, family="serif"))

ggsave("Chromosome_5.png", dpi = 300, limitsize = F)


p1 <- predict(plsFit, train)
confusionMatrix(p1, train$Cancer)

p2 <- predict(plsFit, test)
confusionMatrix(p2, test$Cancer)

####GA-PLS Chromosome 6
set.seed(12)
X <- Num.mat
y <- Cancer
xga <- ga_pls(y, X, GA.threshold = 10, iters = 100, popSize = 100)
Xga <- Mx[,xga$ga.selection]  
Xga$Cancer <- Mx$Cancer
Xga$Cancer <- factor(Xga$Cancer)


ind <- sample(2, nrow(Xga), replace = T, prob = c(0.8, 0.2))
train <- Xga[ind == 1,]
test <- Xga[ind == 2,]
str(test)


set.seed(13)
plsFit <- train(
  Cancer ~ .,
  data = train,
  method = "pls",
  trControl = ctrl,
  metric = "ROC",
  tuneLength = 5,
  tuneGrid = expand.grid(ncomp = 1:5))

plsFit
varImp(plsFit)

ggplot(plsFit)+ theme_bw()+theme(text=element_text(size = 16, family="serif"))

ggsave("Chromosome_6.png", dpi = 300, limitsize = F)


p1 <- predict(plsFit, train)
confusionMatrix(p1, train$Cancer)

p2 <- predict(plsFit, test)
confusionMatrix(p2, test$Cancer)


####GA-PLS Chromosome 7
set.seed(14)
X <- Num.mat
y <- Cancer
xga <- ga_pls(y, X, GA.threshold = 10, iters = 100, popSize = 100)
Xga <- Mx[,xga$ga.selection]  
Xga$Cancer <- Mx$Cancer
Xga$Cancer <- factor(Xga$Cancer)


ind <- sample(2, nrow(Xga), replace = T, prob = c(0.8, 0.2))
train <- Xga[ind == 1,]
test <- Xga[ind == 2,]
str(test)


set.seed(15)
plsFit <- train(
  Cancer ~ .,
  data = train,
  method = "pls",
  trControl = ctrl,
  metric = "ROC",
  tuneLength = 5,
  tuneGrid = expand.grid(ncomp = 1:5))

plsFit
varImp(plsFit)

ggplot(plsFit)+ theme_bw()+theme(text=element_text(size = 16, family="serif"))

ggsave("Chromosome_7.png", dpi = 300, limitsize = F)


p1 <- predict(plsFit, train)
confusionMatrix(p1, train$Cancer)

p2 <- predict(plsFit, test)
confusionMatrix(p2, test$Cancer)

####GA-PLS Chromosome 8
set.seed(19)
X <- Num.mat
y <- Cancer
xga <- ga_pls(y, X, GA.threshold = 10, iters = 100, popSize = 100)
Xga <- Mx[,xga$ga.selection]  
Xga$Cancer <- Mx$Cancer
Xga$Cancer <- factor(Xga$Cancer)


ind <- sample(2, nrow(Xga), replace = T, prob = c(0.8, 0.2))
train <- Xga[ind == 1,]
test <- Xga[ind == 2,]
str(test)


set.seed(20)
plsFit <- train(
  Cancer ~ .,
  data = train,
  method = "pls",
  trControl = ctrl,
  metric = "ROC",
  tuneLength = 5,
  tuneGrid = expand.grid(ncomp = 1:5))

plsFit
varImp(plsFit)

ggplot(plsFit)+ theme_bw()+theme(text=element_text(size = 16, family="serif"))

ggsave("Chromosome_8.png", dpi = 300, limitsize = F)


p1 <- predict(plsFit, train)
confusionMatrix(p1, train$Cancer)

p2 <- predict(plsFit, test)
confusionMatrix(p2, test$Cancer)


####GA-PLS Chromosome 9
set.seed(21)
X <- Num.mat
y <- Cancer
xga <- ga_pls(y, X, GA.threshold = 10, iters = 100, popSize = 100)
Xga <- Mx[,xga$ga.selection]  
Xga$Cancer <- Mx$Cancer
Xga$Cancer <- factor(Xga$Cancer)


ind <- sample(2, nrow(Xga), replace = T, prob = c(0.8, 0.2))
train <- Xga[ind == 1,]
test <- Xga[ind == 2,]
str(test)


set.seed(22)
plsFit <- train(
  Cancer ~ .,
  data = train,
  method = "pls",
  trControl = ctrl,
  metric = "ROC",
  tuneLength = 5,
  tuneGrid = expand.grid(ncomp = 1:5))

plsFit
varImp(plsFit)

ggplot(plsFit)+ theme_bw()+theme(text=element_text(size = 16, family="serif"))

ggsave("Chromosome_9.png", dpi = 300, limitsize = F)


p1 <- predict(plsFit, train)
confusionMatrix(p1, train$Cancer)

p2 <- predict(plsFit, test)
confusionMatrix(p2, test$Cancer)


####GA-PLS Chromosome 10
set.seed(24)
X <- Num.mat
y <- Cancer
xga <- ga_pls(y, X, GA.threshold = 10, iters = 100, popSize = 100)
Xga <- Mx[,xga$ga.selection]  
Xga$Cancer <- Mx$Cancer
Xga$Cancer <- factor(Xga$Cancer)


ind <- sample(2, nrow(Xga), replace = T, prob = c(0.8, 0.2))
train <- Xga[ind == 1,]
test <- Xga[ind == 2,]
str(test)


set.seed(25)
plsFit <- train(
  Cancer ~ .,
  data = train,
  method = "pls",
  trControl = ctrl,
  metric = "ROC",
  tuneLength = 5,
  tuneGrid = expand.grid(ncomp = 1:5))

plsFit
varImp(plsFit)

ggplot(plsFit)+ theme_bw()+theme(text=element_text(size = 16, family="serif"))

ggsave("Chromosome_10.png", dpi = 300, limitsize = F)


p1 <- predict(plsFit, train)
confusionMatrix(p1, train$Cancer)

p2 <- predict(plsFit, test)
confusionMatrix(p2, test$Cancer)


###Validation of the 16 VOMs selected 5/10 by GA-PLS (Chromosome 11)
library(caret)
library(pls)
df <- read.csv("ChromosomeXI.csv", sep=",",header = TRUE) 
Num <- df[,7:23]
Num$Cancer <- factor(Num$Cancer)
Num.mat <- as.matrix(Num)


ind <- sample(2, nrow(Num), replace = T, prob = c(0.8, 0.2))
train <- Num[ind == 1,]
test <- Num[ind == 2,]
str(test)

set.seed(26)
ctrl <- trainControl(method = "repeatedcv",
                     number = 10,
                     repeats = 10,
                     classProbs = TRUE,
                     summaryFunction = twoClassSummary)

set.seed(27)
plsFit <- train(
  Cancer ~ .,
  data = train,
  method = "pls",
  trControl = ctrl,
  metric = "ROC",
  tuneLength = 5,
  tuneGrid = expand.grid(ncomp = 1:5))

plsFit
varImp(plsFit)

ggplot(plsFit)+ theme_bw()+theme(text=element_text(size = 16, family="serif"))

ggsave("Chromosome_11.png", dpi = 300, limitsize = F)

p1 <- predict(plsFit, train)
confusionMatrix(p1, train$Cancer)

p2 <- predict(plsFit, test)
confusionMatrix(p2, test$Cancer)


###Validation of the 7 VOMs selected 8/10 by GA-PLS (Chromosome 12)
library(caret)
library(pls)
df <- read.csv("ChromosomeXII.csv", sep=",",header = TRUE) 
Num <- df[,7:14]
Num$Cancer <- factor(Num$Cancer)
Num.mat <- as.matrix(Num)


ind <- sample(2, nrow(Num), replace = T, prob = c(0.8, 0.2))
train <- Num[ind == 1,]
test <- Num[ind == 2,]
str(test)

set.seed(28)
ctrl <- trainControl(method = "repeatedcv",
                     number = 10,
                     repeats = 10,
                     classProbs = TRUE,
                     summaryFunction = twoClassSummary)

set.seed(29)
plsFit <- train(
  Cancer ~ .,
  data = train,
  method = "pls",
  trControl = ctrl,
  metric = "ROC",
  tuneLength = 5,
  tuneGrid = expand.grid(ncomp = 1:5))

plsFit
varImp(plsFit)

ggplot(plsFit)+ theme_bw()+theme(text=element_text(size = 16, family="serif"))

ggsave("Chromosome_12.png", dpi = 300, limitsize = F)

p1 <- predict(plsFit, train)
confusionMatrix(p1, train$Cancer)

p2 <- predict(plsFit, test)
confusionMatrix(p2, test$Cancer)
