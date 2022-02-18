##R-code Husband & McIntyre 2022

library(multcomp)
library(multcompView)
library(dplyr)
library(tidyr)
library(ggplot2)
library(PerformanceAnalytics)
library(corrplot)
library(vegan)
library(ggpubr)
library(BiodiversityR)
library(labdsv)
library(lattice)
library(MASS)
library(MVA)
library(optpart)
library(picante)
library(stats)
library(ggrepel)
library(factoextra)
library(ggforce)
library(rcartocolor)

setwd("C://Users/whale/OneDrive/Desktop/GitHub_R")
abio2020<-read.csv("2020_Abiotic.csv", header = TRUE, row.names = 1) #raw data
abio2021<-read.csv("2021_Abiotic.csv", header = TRUE, row.names = 1) #raw data

abio2020_n<-read.csv("2020_AbioNorm.csv", header = TRUE, row.names = 1) 
#normalized data
abio2021_n<-read.csv("2021_AbioNorm.csv", header = TRUE, row.names = 1) 
#normalized data

A_2020_n<-abio2020_n[,2:14] #select abiotic fields 
A_2021_n<-abio2021_n[,2:14] #select abiotic fields 

#Correlation Matrix to remove highly correlated variables 
M1 = cor(A_2020_n)
corrplot(M1, method = 'number') #remove Cond, TDS, NH3

M2 = cor(A_2021_n)
corrplot(M2, method = 'number') #remove NH3, TDS

abio2020_rem<-read.csv("2020_AbioNorm_rem.csv", header = TRUE, row.names = 1) 
#removed correlated variables
abio2021_rem<-read.csv("2021_AbioNorm_rem.csv", header = TRUE, row.names = 1) 
#removed correlated variables 

abio2020_rem_anosim<-read.csv("2020_AbioNorm_rem_anosim.csv", header = TRUE, 
                              row.names = 1) #collapsed multi test sites
abio2021_rem_anosim<-read.csv("2021_AbioNorm_rem_anosim.csv", header = TRUE, 
                              row.names = 1) #collapsed multi test sites

#analysis of similarity (Anosim)
comm2020<-read.csv("comm2020.csv", header = TRUE, row.names = 1) 
#binary site x species matrix
comm2021<-read.csv("comm2021.csv", header = TRUE, row.names = 1) 
#binary site x species matrix

ano2020 = anosim(comm2020, abio2020_rem_anosim$Type, distance= "bray", 
                 permutations = 9999)
ano2020

ano2021 = anosim(comm2021, abio2021_rem_anosim$Type, distance= "bray", 
                 permutations = 9999)
ano2021

#Manova with abiotic variables 
manova2020<-manova(cbind(NO3, SO2, P, NO2, NTU, Basin_Fill,
                                 Veg, Na, H2O_Temp, pH)
                           ~abio2020_rem$Type, data=abio2020_rem)
summary(manova2020, intercept=TRUE)
summary.aov(manova2020)

manova2021<-manova(cbind(NO3, SO2, P, NO2, Basin_Fill,
                         Veg, Na, H2O_Temp, pH, Cond)
                   ~abio2021_rem$Type, data=abio2021_rem)
summary(manova2021, intercept=TRUE)
summary.aov(manova2021)

#posthoc ANOVA and Tukey's test for 2021 data
class(abio2021_rem$Type)
abio2021_rem$Type<-as.factor(abio2021_rem$Type)
class(abio2021_rem$Type)

TypeA_norm_rem21<-abio2021_rem[,1] #wetland type 

#Sulfate 
anova<-aov(SO2 ~ TypeA_norm_rem21, data = abio2021_rem)
tukey<-TukeyHSD(anova)
print(tukey)

cld<-multcompLetters4(anova, tukey)
Tk<-group_by(abio2021_rem, abio2021_rem$Type) %>%
  summarise(mean=mean(SO2), quant = quantile(SO2, probs = 0.75)) %>%
  arrange(desc(mean))
as.factor(abio2021_rem$Type)

cld<-as.data.frame.list(cld$Type)
Tk$cld <-cld$Letters

S<-ggplot(abio2021_rem, aes(TypeA_norm_rem21, abio2021_rem$SO2))+
  geom_boxplot()+
  labs(x= "Wetland Type", y = "Sulfate")+
  geom_text(Tk, mapping = aes(x = Tk$`abio2021_rem$Type`, y = Tk$quant, label = cld), 
            vjust = 1.1, hjust = 1, size =5)
S


#P 
anova2<-aov(P ~ TypeA_norm_rem21, data = abio2021_rem)
tukey2<-TukeyHSD(anova2)
print(tukey2)

cld2<-multcompLetters4(anova2, tukey2)
Tk2<-group_by(abio2021_rem, abio2021_rem$Type) %>%
  summarise(mean=mean(P), quant = quantile(P, probs = 0.75)) %>%
  arrange(desc(mean))
as.factor(abio2021_rem$Type)

cld2<-as.data.frame.list(cld2$Type)
Tk2$cld2 <-cld2$Letters

P<- ggplot(abio2021_rem, aes(TypeA_norm_rem21, P))+
  geom_boxplot()+
  labs(x= "Wetland Type", y = "Phosphorous")+
  geom_text(Tk2, mapping = aes(x = Tk2$`abio2021_rem$Type`, y = Tk2$quant, label = cld2), 
            vjust = -1, hjust = -1, size =5)  
P

#NO2
anova3<-aov(NO2 ~ TypeA_norm_rem21, data = abio2021_rem)
tukey3<-TukeyHSD(anova3)
print(tukey3)

cld3<-multcompLetters4(anova3, tukey3)
Tk3<-group_by(abio2021_rem, abio2021_rem$Type) %>%
  summarise(mean=mean(NO2), quant = quantile(NO2, probs = 0.75)) %>%
  arrange(desc(mean))
as.factor(abio2021_rem$Type)
cld3<-as.data.frame.list(cld3$Type)
Tk3$cld3 <-cld3$Letters

N<-ggplot(abio2021_rem, aes(TypeA_norm_rem21, NO2))+
  geom_boxplot()+
  labs(x= "Wetland Type", y = "Nitrate")+
  geom_text(Tk3, mapping = aes(x = Tk3$`abio2021_rem$Type`, y = Tk3$quant, label = cld3), 
            vjust = -1, hjust = -1, size =5)  
N

#% Edge Vegetation 
anova6<-aov(Veg ~ TypeA_norm_rem21, data = abio2021_rem)
tukey6<-TukeyHSD(anova6)
print(tukey6)

cld6<-multcompLetters4(anova6, tukey6)
Tk6<-group_by(abio2021_rem, abio2021_rem$Type) %>%
  summarise(mean=mean(Veg), quant = quantile(Veg, probs = 0.75)) %>%
  arrange(desc(mean))
as.factor(abio2021_rem$Type)
cld6<-as.data.frame.list(cld6$Type)
Tk6$cld6 <-cld6$Letters

V<-ggplot(abio2021_rem, aes(TypeA_norm_rem21, abio2021_rem$Veg))+
  geom_boxplot()+
  labs(x= "Wetland Type", y = "Prop. Edge Vegetation")+
  geom_text(Tk6, mapping = aes(x = Tk6$`abio2021_rem$Type`, y = Tk6$quant, label = cld6), 
            vjust = -1, hjust = -1, size =5)
V

#Salinity
anova7<-aov(Na ~ TypeA_norm_rem21, data = abio2021_rem)
tukey7<-TukeyHSD(anova7)
print(tukey7)

cld7<-multcompLetters4(anova7, tukey7)
Tk7<-group_by(abio2021_rem, abio2021_rem$Type) %>%
  summarise(mean=mean(Na), quant = quantile(Na, probs = 0.75)) %>%
  arrange(desc(mean))
as.factor(abio2021_rem$Type)
cld7<-as.data.frame.list(cld7$Type)
Tk7$cld7 <-cld7$Letters

Na<-ggplot(abio2021_rem, aes(TypeA_norm_rem21, Na))+
  geom_boxplot()+
  labs(x= "Wetland Type", y = "Salinity")+
  geom_text(Tk7, mapping = aes(x = Tk7$`abio2021_rem$Type`, y = Tk7$quant, label = cld7), 
            vjust = -1, hjust = -1, size =5)
Na

#NMDS
community<-read.csv("SiteXspp.csv", header = TRUE, row.names = 1)
set.seed(15)
MDS1 <- metaMDS(community, distance="bray", binary=TRUE, k=2, try = 999, zerodist="add")
plot(MDS1)
MDS1

#extract site scores
data.scores <- as.data.frame(scores(MDS1))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$site <- rownames(data.scores)  # create a column of site names, from the rownames of data.scores
data.scores$Type = abio2021_rem_anosim$Type
head(data.scores)  #look at the data

#extract species scores
species.scores <- as.data.frame(scores(MDS1, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
head(species.scores)  #look at the data

theme.custom<- theme_bw()+ 
  theme(axis.title = element_text(size=15), #text larger
        panel.border = element_rect(size =2)) #border thicker

ggplot() + 
  geom_point(data=data.scores, aes(x=NMDS1, y=NMDS2, colour=abio2021_rem_anosim$Type), alpha= 0.7, size=3, shape= 17)+
  theme.custom+
  stat_ellipse(data=data.scores, aes(x=NMDS1, y=NMDS2,color=abio2021_rem_anosim$Type),type = "t", lwd=2)+
  scale_color_manual(values = c("#882255", "#999933", "#44AA99", "#6699CC"))

dev.off() # if 'invalid graphics state' error appears 
citation(package = "rcartocolor")

#Diversity Metrics etc. 
specpool(community)
all.sp<-specpool(community, pool = abio2021_rem_anosim$Type)
all.sp

gamma<- specpool(community)
gamma

#Species Accumulation Curve
Accum.1<-accumcomp(community, y=abio2021_rem_anosim, factor = 'Type', method = 'exact', conditioned = FALSE, plotit = FALSE)

accum.long1 <- accumcomp.long(Accum.1, ci=NA, label.freq=5)

theme.custom<- theme_bw()+ 
  theme(axis.title = element_text(size=15), #text larger
        panel.border = element_rect(size =2)) #border thicker

plotgg1 <- ggplot(data=accum.long1, aes(x = Sites, y = Richness, ymax = UPR, ymin = LWR)) +
  geom_line(aes(colour=Grouping), size=2)+
  scale_color_manual(values = c("#882255", "#999933", "#44AA99", "#6699CC"))+
  geom_ribbon(aes(colour=Grouping, fill=after_scale(alpha(colour, 0.3))), 
              show.legend=FALSE)+ theme.custom
plotgg1


#beta diversity

all2 <- rowsum(community, group = abio2021_rem_anosim$Type)
betadiver(all2, method="sor")

all3 <- rowsum(community, group = abio2021_rem_anosim$Type)
betadiver(all3, method="j")

gb.b <- betadiver(community, method = "w")
gb.ano <- anosim(gb.b, abio2021_rem_anosim$Type)
gb.ano 

gb.b2 <- betadiver(community, method = "sor")
gb.ano2 <- anosim(gb.b2, abio2021_rem_anosim$Type)
gb.ano2

#raw data boxplots for significant manova variables 


theme(axis.title = element_text(size = 20))   

abio2020
A<- ggplot(abio2020, aes(Type,Na))+
  geom_boxplot()+
  labs(y = "Salinity (ppm)")+
  theme.custom
A + theme(axis.title=element_text(size=10))

B<-ggplot(abio2020, aes(Type,SO2))+
  geom_boxplot()+
  labs(y = "Sulfate (g/mL)")+
  theme.custom
B + theme(axis.title=element_text(size=10))

C<-ggplot(abio2020, aes(Type,Veg))+
  geom_boxplot()+
  labs(x="Wetland Type", y = "% Edge Veg.")+
  theme.custom
C + theme(axis.title=element_text(size=10))

D<-ggplot(abio2020, aes(Type,Basin_Fill))+
  geom_boxplot()+
  labs(x="Wetland Type", y = "% Basin Fill")+
  theme.custom
D + theme(axis.title=element_text(size=10))


abio2021

E<- ggplot(abio2021, aes(Type,Na))+
  geom_boxplot()+
  labs(y = "Salinity (ppm)")+
  theme.custom
E + theme(axis.title=element_text(size=10))

F <-ggplot(abio2021, aes(Type,SO2))+
  geom_boxplot()+
  labs(y = "Sulfate (g/mL)")+
  theme.custom
F + theme(axis.title=element_text(size=10))

G<-ggplot(abio2021, aes(Type,Veg))+
  geom_boxplot()+
  labs(y = "% Edge Veg.")+
  theme.custom
G + theme(axis.title=element_text(size=10))

I<-ggplot(abio2021, aes(Type,P))+
  geom_boxplot()+
  labs(x="Wetland Type", y = "Phosphorous (g/mL)")+
  theme.custom
I + theme(axis.title=element_text(size=10))

H<- ggplot(abio2021, aes(Type,NO2))+
  geom_boxplot()+
  labs(x="Wetland Type", y = "Nitrite (g/mL)")+
  theme.custom
H + theme(axis.title=element_text(size=10))

ggarrange(A, B, C, D, E, F, G, H, I, 
          labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"))


