
library(rstatix)
library(ggplot2)
library(ggpubr)
library(emmeans) #library(emmeans) 
basedir <- '/Users/WIWIFH/Dropbox/Projects/EXCOOL/data'; # set base directory 
setwd(basedir) # set working directory
(list.files()) # see the files
#sub_i <- 2;
#run_i <- 2; # run number 
#csvfile <- paste(basedir,"/","EXC_task_Sub-EX00", (sub_i), "_run0",(run_i),".csv", sep = "")
#data<- read.csv(csvfile, header=TRUE) # Load CSV files

# LOAD CSV FILES per participants
subjname <- "EX002" # subject number
file_list <- list.files(path = basedir, pattern=subjname, full.names = TRUE,ignore.case = TRUE)## file list 

final<-NULL ## 빈 변수 (NULL) 선언 
for(i in 1:length(file_list)){ ## file_list의 1번째방부터 6번째(길이)까지    
     file<-read.csv(file_list[i]) ## file_list의 i번째 방에 있는 csv 파일 read 
     final<-rbind(final,file)  
     ## i가 1일 때 : 빈 값과 1번째 csv rbind -> 1번째 csv 
     ## i가 2일 때 : 1번째 csv와 2번째 csv rbind -> 1번째+2번째 csv 
     ## i가 3일 때 : (1번째+2번째 csv)와 세번째 csv rbind -> 1번째+2번째+3번째 csv 
     cat("\n",i) 
}
data <- final
# example 
# <- "A$B"는 자료A에 포함된 변수B(또는 행렬A의 열B)를 의미
barplot(data$ratings,horiz=TRUE)
barplot(data$ratings,horiz=FALSE, main = "Aggrement ratings")
plot(data$Wvalence, data$ratings, main = "Aggrement ratings")

# Load the vioplot library
#install.packages("viridis")
# Libraries

colnames(data) 
# Draw the violin plot
# Grouped
data<-data.frame(data)
data2 <- transform(data, 
                   TGroup = cut(Tcondition, breaks = c(0,1,2,3), include.lowest = FALSE, 
                                   right = TRUE, 
                                   labels = c("SELF", "FRIEND", "CELEB")),
                   WGroup = cut(Wvalence, breaks = c(-2,-1,1), include.lowest = FALSE, 
                                 right = TRUE, 
                                 labels = c("NEG", "POS")),
                   RatingsBIN = cut(ratings, breaks = c(-1,0,1), include.lowest = FALSE, 
                                   right = TRUE, 
                                   labels = c("NEG", "POS"))
                  )
# Summary statistics
data2 %>%
     group_by(TGroup, WGroup) %>%
     get_summary_stats(ratings, type = "mean_sd")

#ggplot: visualization
res_grap<-ggplot(data2,                         
       aes(x = TGroup,
           y = ratings,
           fill = WGroup)) +
     geom_violin(trim = FALSE, position="dodge", alpha=0.5)+
     geom_dotplot(binaxis='y', stackdir='center',position=position_dodge(1)) +
     scale_fill_manual(values=c("#67a9cf", "#ef8a62", "#56B4E9")) + 
     labs(title="Ratings by Valence by Target",x="Target", y = "Ratings (-0.5~+0.5)")+
     theme_classic()

sepres<-data2 %>%
     group_by(RunN, TGroup, WGroup) %>%
     get_summary_stats(ratings, type = "mean_sd")
sepres <-data.frame(sepres)

#geom_point(shape=21, color="black", fill="#69b3a2", size=6) +
res_grap3<-ggplot(sepres, aes(x=RunN, y=mean,shape=WGroup)) +
     theme_classic() +
     ggtitle("assd") + 
     labs(title="Ratings by Run",x="RUN No.", y = "Ratings (-0.5~+0.5)")
res_grap3 + geom_point(aes(colour = factor(TGroup), shape = factor(WGroup)),size=4,alpha = 0.8) + 
     geom_line( aes(colour = factor(TGroup)),linetype = 1,alpha = 0.8,size=0.8)


# statistical test      
res.aov <- data2 %>% anova_test(ratings ~ WGroup + TGroup + TGroup*WGroup)
res.aov


pwc <- data2 %>% 
     group_by(WGroup) %>%
     emmeans_test(ratings ~ TGroup, p.adjust.method = "bonferroni") 
pwc




data2_pos <- subset(data2,data2$WGroup=="POS") # When word is positive 
data2_neg <- subset(data2,data2$WGroup=="NEG") # When word is negative

(m_Diff_pos <- tapply(data2_pos$ratings, data2_pos$TGroup, mean))
(m_Diff_neg <- tapply(data2_neg$ratings, data2_neg$TGroup, mean))
barplot(m_Diff_pos) # 평균의 막대그래프
barplot(m_Diff) # 평균의 막대그래프

plot(data2$RatingsBIN)
table(data2$RatingsBIN)
table(data2$RatingsBIN)/length(data2$RatingsBIN) # 상대 빈도

# ANOVA
two.way <- aov(ratings ~ WGroup + TGroup + TGroup*WGroup, data = data2)
summary(two.way)

#
# LOAD CSV FILES ALL participants
file_list <- list.files(path = basedir, pattern="*.csv", full.names = TRUE,ignore.case = TRUE)## file list 

final<-NULL ## 빈 변수 (NULL) 선언 
for(i in 1:length(file_list)){ ## file_list의 1번째방부터 6번째(길이)까지    
     file<-read.csv(file_list[i]) ## file_list의 i번째 방에 있는 csv 파일 read 
     temp_subnames<-regmatches(file_list[i],gregexpr("Sub-.....", file_list[i])[1])
     temp_df<-data.frame(rep(temp_subnames[[1]],times=length(file$RunN)))
     colnames(temp_df)<-"subjID"
     file <- cbind(temp_df,file)
     final<-rbind(final,file)  
     ## i가 1일 때 : 빈 값과 1번째 csv rbind -> 1번째 csv 
     ## i가 2일 때 : 1번째 csv와 2번째 csv rbind -> 1번째+2번째 csv 
     ## i가 3일 때 : (1번째+2번째 csv)와 세번째 csv rbind -> 1번째+2번째+3번째 csv 
     cat("\n",i) 
}
#temp_subnames<-regmatches(file_list[2],gregexpr("Sub-.....", file_list[2])[1])
data <- final
# Grouped
data<-data.frame(data)
data2 <- transform(data, 
                   TGroup = cut(Tcondition, breaks = c(0,1,2,3), include.lowest = FALSE, 
                                right = TRUE, 
                                labels = c("SELF", "FRIEND", "CELEB")),
                   WGroup = cut(Wvalence, breaks = c(-2,-1,1), include.lowest = FALSE, 
                                right = TRUE, 
                                labels = c("NEG", "POS")),
                   RatingsBIN = cut(ratings, breaks = c(-1,0,1), include.lowest = FALSE, 
                                    right = TRUE, 
                                    labels = c("NEG", "POS"))
)
# Summary statistics
data2 %>%
     group_by(TGroup, WGroup,subjID) %>%
     get_summary_stats(ratings, type = "mean_sd")

#ggplot: visualization
res_grap<-ggplot(data2,                         
                 aes(x = subjID,
                     y = ratings,
                     fill = WGroup)) +
     geom_violin(trim = FALSE, position="dodge", alpha=0.5)+
     geom_dotplot(binaxis='y', stackdir='center',position=position_dodge(1),dotsize=0.4) +
     scale_fill_manual(values=c("#67a9cf", "#ef8a62", "#56B4E9")) + 
     labs(title="Ratings by Word Valence",x="Target", y = "Ratings (-0.5~+0.5)")+
     theme_classic()
res_grap

res_grap2<-ggplot(data2,                         
                 aes(x = subjID,
                     y = ratings,
                     fill = TGroup)) +
     geom_violin(trim = FALSE, position="dodge", alpha=0.5)+
     geom_dotplot(binaxis='y', stackdir='center',position=position_dodge(1),dotsize=0.4) +
     scale_fill_manual(values=c("#4daf4a","#984ea3", "#ff7f00", "#56B4E9")) + 
     labs(title="Ratings by Target",x="Target", y = "Ratings (-0.5~+0.5)")+
     theme_classic()
res_grap2


# statistical test      
library(lme4) # for the analysis
library(RColorBrewer) # needed for some extra colours in one of the graphs
library(lmerTest)# to get p-value estimations that are not part of the standard lme4 packages
library(effects) 

m1 <- lmer(ratings ~  1 + TGroup + WGroup + TGroup:WGroup + ( 1 + TGroup + WGroup + WGroup:TGroup |subjID),
                           data    = data2,na.action = na.omit) #to run the model

m2 <- lmer(ratings ~  1 + Tcondition + Wvalence + Tcondition:Wvalence + ( 1 |subjID),
           data    = data2,na.action = na.omit) #to run the model

model <- summary(m1)
model2 <- summary(m2)
#Plot model data
plot(allEffects(m2))
plot(predictorEffects(m2, ~ Tcondition,xlevels=2))
plot(predictorEffects(m2, ~ Wvalence,xlevels=3))


