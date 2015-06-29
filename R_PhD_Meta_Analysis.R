#libraries
library(metafor)
library(ggplot2)

theme_set(theme_bw(base_size=18))

#extract bioclim variables
#setwd("S:\\users\\kambach\\PhD\\Project_I_Meta_analysis_BD_and_stability\\Data\\R")
#data = read.csv("Coding_Table_Shannon_Diversity.csv",header=TRUE, sep=";",dec=",",na.strings="NA")
#clim = getData('worldclim', var = 'bio', res = 2.5, lon = -100, lat = 45)
#climate = extract(clim$bio1,cbind(data$GPS_DEZ_X,data$GPS_DEZ_Y),method="bilinear")
#data = cbind(data,climate)
#names(data)[35] = "bio1"
#data$bio1 = data$bio / 10
#write.table(data,"Coding_Table_Shannon_Diversity_with_clim.txt",sep="\t",dec=",")

#--------------start------------------

setwd("C:\\Users\\kambach\\Desktop\\aktuelle Arbeiten\\Meta-analysis Herbivory\\Code\\meta-analysis-herbivory")
data = read.table("Coding_Table_Shannon_Diversity_with_clim.txt",header=TRUE, sep="\t",dec=",",na.strings="NA")
data.frame(names(data))
str(data)


#functions
as.numeric.factor = function(x) {as.numeric(levels(x))[x]}
fisher.to.r = function(x){round(((exp(1)^(2*x))-1)/((exp(1)^(2*x))+1),6)}


#calculate z-score
data$R_VALUE = (data$BETA_ESTIMATE*sqrt(data$R_SQUARE))/ sqrt(data$BETA_ESTIMATE^2)
data$Z_SCORE = 0.5 * log((1 + data$R_VALUE)/(1-data$R_VALUE))
data$VARIANCE_ESTIMATE = 1/(data$PLOT_NUMBER-3)



#order and clean the table and construct a unique study identifier
data = data[order(data$Z_SCORE),]
data$UNIQUE_ID <- paste(data$AUTHOR,data$YEAR,data$FOREST_NAME,data$SAMPLING_ON, sep=" ")
data$UNIQUE_ID <- paste(data$UNIQUE_ID,data$HERBIVORE_SPECIES, sep=": ")
data$UNIQUE_ID <- paste(data$UNIQUE_ID,data$FOCAL_SPECIES, sep=" on ")

#exclude non-shannon treee diversities
dat.filtered = subset(data, TREE_DIV_METRIC %in% "shannon")
#write.table(dat.filtered,"for_the_map.txt",sep=";",dec=".",quote=FALSE,row.names=FALSE)

# convert GPS_DEZ_Y to just positive values
#data$GPS_DEZ_Y <- sqrt(data$GPS_DEZ_Y^2)



#heterogeneity between herbivore damage and abundance--------
#####################################################
dat.damage.abundance <- subset(dat.filtered, !(VARIANCE_ESTIMATE == Inf)  & VARIANCE_ESTIMATE >= 0 & !(Z_SCORE %in% NaN) & !(Z_SCORE %in% Inf)
                     & !(Z_SCORE %in% NA) 
                     & HERBIVORE_GUILD %in% "insects"
                     & !(DAMAGE_OR_HERBIVORES %in% "richness"))

#write.table(dat.damage.abundance,"data_frame_damage_abundance_insects.csv", sep=";",dec=",",row.names=F)

#without difference between damage, abundance and richness
rma.damage.abundance.empty <- rma.mv(Z_SCORE, VARIANCE_ESTIMATE,method="REML",
                            random=~ list(1|factor(CASE_ID))+list(1|factor(STUDY_ID)),
                            data=dat.damage.abundance)
rma.damage.abundance.empty
forest(rma.damage.abundance.empty)
funnel(rma.damage.abundance.empty)
regtest(rma.damage.abundance.empty)


#plot with ggplot2
plot.data.empty = rbind(data.frame("CASE_ID" = as.factor("grand mean"),
                                    "ES" =rma.damage.abundance.empty$b,
                                    "SE"=rma.damage.abundance.empty$se),
                         data.frame("CASE_ID" = as.factor("spaceholder"), #space holder
                                    "ES" = 0,
                                    "SE"= 0),
                         data.frame("CASE_ID" =factor(dat.damage.abundance$CASE_ID,levels=dat.damage.abundance$CASE_ID),
                                    "ES" =rma.damage.abundance.empty$yi,
                                    "SE" =sqrt(rma.damage.abundance.empty$vi)/sqrt(dat.damage.abundance$PLOT_NUMBER)))

plot.data.empty = merge(plot.data.empty,dat.damage.abundance[,c(1,8)],by.x="CASE_ID",by.y="CASE_ID",all.x=TRUE)
plot.data.empty$DAMAGE_OR_HERBIVORES = c("placeholder","placeholder",as.character(plot.data.empty$DAMAGE_OR_HERBIVORES[3:59]))

#plot.data.damage$TYPE1 = c("black",rep(NA,nrow(plot.data.empty) -1))
#plot.data.damage$TYPE1 = c("white",rep(NA,nrow(plot.data.empty) -1))
#plot.data.damage$TYPE2 = c(NA,NA,rep("grey",nrow(plot.data.empty) -2))

levels(plot.data.empty$CASE_ID)
theme_set(theme_bw(base_size=18))

svg(filename="herbivore_damage_abundance.svg",height=12,width=6)  
plot.empty = ggplot(data=(plot.data.empty)) 
plot.empty = plot.empty + geom_blank(aes(x=factor(CASE_ID),y=ES)) + coord_flip() + 
             geom_hline(aes(x=0), lty=2,size=0.5,colour="black") 
plot.empty = plot.empty + geom_point(aes(x=factor(CASE_ID),y=ES,size=plot.data.empty$SE,colour=(DAMAGE_OR_HERBIVORES)))
plot.empty = plot.empty + scale_colour_manual(values = c("grey","black","white"))
plot.empty = plot.empty + scale_size(range=c(1,8))
plot.full = plot.empty + geom_pointrange(aes(x=1,y=rma.damage.abundance.empty$b[1,1],ymin=rma.damage.abundance.empty$b[1,1] - 1.96*rma.damage.abundance.empty$se,ymax= rma.damage.abundance.empty$b + 1.96*rma.damage.abundance.empty$se),
            shape=3, size=1.2,linetype=1)
plot.full +  ylim(c(-2,2)) + theme(legend.position="none",axis.ticks=element_blank()) +
labs(y="Fisher's z",x="Study case ID")
graphics.off()


#for herbivore species richness
dat.richness <- subset(dat.filtered, !(VARIANCE_ESTIMATE == Inf)  & VARIANCE_ESTIMATE >= 0 & !(Z_SCORE %in% NaN) & !(Z_SCORE %in% Inf)
                               & !(Z_SCORE %in% NA) 
                               & HERBIVORE_GUILD %in% "insects"
                               & (DAMAGE_OR_HERBIVORES %in% "richness"))

rma.richness <- rma.mv(Z_SCORE, VARIANCE_ESTIMATE,method="REML",
                       random=~ list(1|factor(CASE_ID)),
                       data=dat.richness)
rma.richness

plot.data.empty = rbind(data.frame("CASE_ID" = as.factor("grand mean"),
                                   "ES" =rma.richness$b,
                                   "SE"=rma.richness$se),
                        data.frame("CASE_ID" = as.factor("spaceholder"), #space holder
                                   "ES" = 0,
                                   "SE"= 0),
                        data.frame("CASE_ID" =factor(dat.richness$CASE_ID,levels=dat.richness$CASE_ID),
                                   "ES" =rma.richness$yi,
                                   "SE" =sqrt(rma.richness$vi)/sqrt(dat.richness$PLOT_NUMBER)))

plot.data.empty = merge(plot.data.empty,dat.richness[,c(1,8)],by.x="CASE_ID",by.y="CASE_ID",all.x=TRUE)
plot.data.empty$DAMAGE_OR_HERBIVORES = c("placeholder","placeholder",as.character(plot.data.empty$DAMAGE_OR_HERBIVORES[3:8]))

#plot.data.damage$TYPE1 = c("black",rep(NA,nrow(plot.data.empty) -1))
#plot.data.damage$TYPE1 = c("white",rep(NA,nrow(plot.data.empty) -1))
#plot.data.damage$TYPE2 = c(NA,NA,rep("grey",nrow(plot.data.empty) -2))

levels(plot.data.empty$CASE_ID)
theme_set(theme_bw(base_size=18))

svg(filename="herbivore_richness2.svg",height=3,width=6)  
plot.empty = ggplot(data=(plot.data.empty)) 
plot.empty = plot.empty + geom_blank(aes(x=factor(CASE_ID),y=ES)) + coord_flip() + 
  geom_hline(aes(x=0), lty=2,size=0.5,colour="black") 
plot.empty = plot.empty + geom_point(aes(x=factor(CASE_ID),y=ES,size=plot.data.empty$SE,colour=(DAMAGE_OR_HERBIVORES)))
plot.empty = plot.empty + scale_colour_manual(values = c("white","black","white"))
plot.empty = plot.empty + scale_size(range=c(2,8))
plot.full = plot.empty + geom_pointrange(aes(x=1,y=rma.richness$b[1,1],ymin=rma.richness$b[1,1] - 1.96*rma.richness$se,ymax= rma.richness$b + 1.96*rma.richness$se),
                                         shape=3, size=1.2,linetype=1)
plot.full +  ylim(c(-2,2)) + theme(legend.position="none",axis.ticks=element_blank()) +
  labs(y="Fisher's z",x="Study case ID")
graphics.off()




#with difference between damage and abundance
rma.damage.abundance.full <- rma.mv(Z_SCORE, VARIANCE_ESTIMATE,method="REML",
                               mods=~ factor(DAMAGE_OR_HERBIVORES)*bio1,
                               random=~ list(1|factor(CASE_ID))+list(1|factor(STUDY_ID)),
                               data=dat.damage.abundance)
rma.damage.abundance.full
forest(rma.damage.abundance.full)

rma.damage.abundance.bio1 <- rma.mv(Z_SCORE, VARIANCE_ESTIMATE,method="REML",
                                    mods=~ bio1,
                                    random=~ list(1|factor(CASE_ID))+list(1|factor(STUDY_ID)),
                                    data=dat.damage.abundance)

rma.damage.abundance.bio1
funnel(rma.damage.abundance.bio1)
regtest(rma.damage.abundance.bio1)

#Comparing the AICc between all three models based on ML
rma.damage.abundance.empty$fit.stats
rma.damage.abundance.full$fit.stats
rma.damage.abundance.bio1$fit.stats


#R2
Z_fitted = fitted(rma.damage.abundance.bio1)
model_var = sum((dat.damage.abundance$Z_SCORE - Z_fitted)^2)
dat_var = sum((dat.damage.abundance$Z_SCORE - mean(dat.damage.abundance$Z_SCORE))^2)
Z_r2 = 1- model_var/dat_var
Z_r2

#fitstats
rma.empty <- rma.mv(Z_SCORE, VARIANCE_ESTIMATE,method="ML",
                    random=~ list(1|factor(CASE_ID))+list(1|factor(STUDY_ID)),
                    data=dat.damage.abundance)
rma.bio1 <- rma.mv(Z_SCORE, VARIANCE_ESTIMATE,method="ML",
                   mods =~ bio1, 
                   random=~ list(1|factor(CASE_ID))+list(1|factor(STUDY_ID)),
                   control=list(optmethod="BFGS"),
                   data=dat.damage.abundance)
fitstats(rma.empty)
fitstats(rma.bio1)




##trim and fill 
trimfill(rma.damage.abundance.empty,side="right",verbose=TRUE )
#funnel plot
funnel(trim.abundance)
#fail safe number
fsn(Z_SCORE, VARIANCE_ESTIMATE,data=dat.damage.abundance,
    type="Rosenthal",alpha=0.05,)  


rma.damage.abundance.bio1


dat.damage.abundance$Z_SCORE_fitted = rma.damage.abundance.bio1$b[1]+(rma.damage.abundance.bio1$b[2]*dat.damage.abundance$bio1)

svg(filename="herbivore_damage_abundance_bio1.svg",height =6, width=6)
plot.empty = ggplot(data=dat.damage.abundance)
plot.empty = plot.empty + geom_point(aes(x=bio1, y=Z_SCORE,size=VARIANCE_ESTIMATE),colour="#7e7e7e")
plot.empty = plot.empty + geom_line(aes(x=bio1,y=Z_SCORE_fitted)) + scale_size(range=(c(3,6)))
plot.full = plot.empty + geom_line(aes(x=bio1,y=0),linetype="dashed")
plot.full + theme(legend.position="none",axis.ticks=element_blank()) + 
  xlim(c(min(dat.damage.abundance$bio1),max(dat.damage.abundance$bio1))) + ylim(c(-2,2)) +
  labs(y="Fisher's z",x="Mean annual temperature") +
  scale_x_continuous(breaks = c(0,5,10,15,20,25)) 
graphics.off()









summary(test1)
rma.damage.abundance.bio1
forest(rma.damage.abundance.bio1)
profile(rma.damage.abundance.bio1,tau2=1)
profile(rma.damage.abundance.bio1,rho=1)
funnel(rma.damage.abundance.bio1)
qqnorm(rma.damage.abundance.bio1)


fitstats(rma.damage.abundance.empty)
fitstats(rma.damage.abundance.bio1)
fitstats(rma.damage.abundance.full)







rma.damage.abundance.full <- rma.mv(Z_SCORE, VARIANCE_ESTIMATE,method="REML",
                                    mods=~ factor(DAMAGE_OR_HERBIVORES) + bio1,
                                    random=~ CASE_ID|STUDY_ID,
                                    data=dat.damage.abundance)
rma.damage.abundance.full
forest(rma.damage.abundance.full)



#-> no difference between damage and abundance



#split the dataset
dat.damage = subset(dat.damage.abundance, DAMAGE_OR_HERBIVORES %in% "damage")
dat.abundance = subset(dat.damage.abundance, DAMAGE_OR_HERBIVORES %in% "abundance")
dat.richness = subset(dat.damage.abundance, DAMAGE_OR_HERBIVORES %in% "richness")

#grand mean for all three aspects
#--------------------------------
rma.damage <- rma(Z_SCORE, VARIANCE_ESTIMATE,method="REML",
                  #random = ~ CASE_ID|STUDY_ID,
                  data=dat.damage)
rma.damage
forest(rma.damage)
#-> no grand mean of zero


rma.abundance <- rma(Z_SCORE, VARIANCE_ESTIMATE,method="REML",
                      # random=~ FOREST_NAME|STUDY_ID,
                        data=dat.abundance)
rma.abundance
forest(rma.abundance)
#-> negative grand mean
funnel(rma.abundance)
#trim and fill 
trimfill(rma.abundance,side="right",verbose=TRUE )
#funnel plot
funnel(trim.abundance)
#fail safe number
fsn(Z_SCORE, VARIANCE_ESTIMATE,data=dat.abundance,
    type="Rosenthal",alpha=0.05,)  

rma.richness <- rma(Z_SCORE, VARIANCE_ESTIMATE,method="REML",
                    #   random=~  factor(CASE_ID)|factor(FOREST_NAME),
                        data=dat.richness)
rma.richness
forest(rma.richness)
#-> Grand mean of zero


#influence of MAT
#----------------
rma.damage.bio1 <- rma(Z_SCORE, VARIANCE_ESTIMATE,method="REML",
                      #    random=~ 1|CASE_ID,
                          mods =~ bio1,
                          data=dat.damage)
rma.damage.bio1
#-> positive influence of MAT

rma.abundance.bio1 <- rma.mv(Z_SCORE, VARIANCE_ESTIMATE,method="REML",
                             #random=~ 1|CASE_ID,
                             mods =~ bio1,
                             data=dat.abundance)
rma.abundance.bio1
#-> positive influence of MAT


#----plotting area-----------------

#damage forest plot
plot.data.damage = rbind(data.frame("CASE_ID" = as.factor("grand mean"),
                                    "ES" =rma.damage$b,
                                    "SE"=rma.damage$se),
                         data.frame("CASE_ID" = as.factor("spaceholder"), #space holder
                                    "ES" = 0,
                                    "SE"= 0),
                         data.frame("CASE_ID" =factor(dat.damage$CASE_ID,levels=dat.damage$CASE_ID),
                                    "ES" =rma.damage$yi,
                                    "SE" =sqrt(rma.damage$vi)/sqrt(dat.damage$PLOT_NUMBER)))

plot.data.damage$TYPE1 = c("black",rep(NA,nrow(plot.data.damage) -1))
plot.data.damage$TYPE2 = c(NA,NA,rep("grey",nrow(plot.data.damage) -2))

levels(plot.data.damage$CASE_ID)
theme_set(theme_bw(base_size=18))

svg(filename="herbivore_damage.svg",height=11,width=6)
plot.empty = ggplot(data=(plot.data.damage)) 
plot.empty = plot.empty + geom_blank(aes(x=factor(CASE_ID),y=ES)) + coord_flip() + geom_hline(aes(x=0), lty=2,size=0.5,colour="black") 
plot.empty = plot.empty + geom_point(aes(x=factor(CASE_ID),y=ES,size=5-(plot.data.damage$SE*5),colour=(TYPE2)))
plot.empty = plot.empty + scale_colour_manual(values = c("#7e7e7e"))
plot.full = plot.empty + geom_pointrange(aes(x=1,y=rma.damage$b[1,1],ymin=rma.damage$b[1,1] - 1.96*rma.damage$se,ymax= rma.damage$b + 1.96*rma.damage$se),
                                         shape=3, size=1.2,linetype=1)
plot.full +  ylim(c(-2,2)) + theme(legend.position="none",axis.ticks=element_blank(),axis.text.y = element_text(colour="grey")) +
  labs(y="Fisher's z",x="Study case ID", title="a) Herbivore damage")
graphics.off()


#abundance forest plot
plot.data.abundance = rbind(data.frame("CASE_ID" = as.factor("grand mean"),
                                    "ES" =rma.abundance$b,
                                    "SE"=rma.abundance$se),
                         data.frame("CASE_ID" = as.factor("spaceholder"), #space holder
                                    "ES" = 0,
                                    "SE"= 0),
                         data.frame("CASE_ID" =factor(dat.abundance$CASE_ID,levels=dat.abundance$CASE_ID),
                                    "ES" =rma.abundance$yi,
                                    "SE" =sqrt(rma.abundance$vi)/sqrt(dat.abundance$PLOT_NUMBER)))

plot.data.abundance$TYPE1 = c("black",rep(NA,nrow(plot.data.abundance) -1))
plot.data.abundance$TYPE2 = c(NA,NA,rep("grey",nrow(plot.data.abundance) -2))

levels(plot.data.abundance$CASE_ID)
theme_set(theme_bw(base_size=18))

svg(filename="herbivore_abundance.svg",height=8,width=6)
plot.empty = ggplot(data=(plot.data.abundance)) 
plot.empty = plot.empty + geom_blank(aes(x=factor(CASE_ID),y=ES)) + coord_flip() + geom_hline(aes(x=0), lty=2,size=0.5,colour="black") 
plot.empty = plot.empty + geom_point(aes(x=factor(CASE_ID),y=ES,size=5-(plot.data.abundance$SE*5),colour=(TYPE2)))
plot.empty = plot.empty + scale_colour_manual(values = c("#7e7e7e"))
plot.full = plot.empty + geom_pointrange(aes(x=1,y=rma.abundance$b[1,1],ymin=rma.abundance$b[1,1] - 1.96*rma.abundance$se,ymax= rma.abundance$b + 1.96*rma.abundance$se),
                                         shape=3, size=1.2,linetype=1)
plot.full +  ylim(c(-2,2)) + theme(legend.position="none",axis.ticks=element_blank(),axis.text.y = element_text(colour="grey")) +
  labs(y="Fisher's z",x="Study case ID", title="b) Herbivore abundance")
graphics.off()


#damage ~ bio1
dat.damage$Z_SCORE_fitted = rma.damage.bio1$b[1]+(rma.damage.bio1$b[2]*dat.damage$bio1)

svg(filename="herbivore_damage_bio1.svg",height =6, width=6)
plot.empty = ggplot(data=dat.damage)
plot.empty = plot.empty + geom_point(aes(x=bio1, y=Z_SCORE,size=VARIANCE_ESTIMATE),colour="#7e7e7e")
plot.empty = plot.empty + geom_line(aes(x=bio1,y=Z_SCORE_fitted)) + scale_size(range=(c(3,6)))
plot.full = plot.empty + geom_line(aes(x=bio1,y=0),linetype="dashed")
plot.full + theme(legend.position="none",axis.ticks=element_blank()) + 
  xlim(c(min(dat.damage.abundance$bio1),max(dat.damage.abundance$bio1))) + ylim(c(-2,2)) +
  labs(y="Fisher's z",x="Mean annual temperature", title="a) Herbivore damage")
graphics.off()


#abundance ~ bio1
dat.abundance$Z_SCORE_fitted = rma.abundance.bio1$b[1]+(rma.abundance.bio1$b[2]*dat.abundance$bio1)

svg(filename="herbivore_abundance_bio1.svg",height=6,width=6)
plot.empty = ggplot(data=dat.abundance)
plot.empty = plot.empty + geom_point(aes(x=bio1, y=Z_SCORE,size=VARIANCE_ESTIMATE),colour="#7e7e7e")
plot.empty = plot.empty + geom_line(aes(x=bio1,y=Z_SCORE_fitted)) + scale_size(range=(c(3,6)))
plot.full = plot.empty + geom_line(aes(x=bio1,y=0),linetype="dashed")
plot.full + theme(legend.position="none",axis.ticks=element_blank()) + 
  xlim(c(min(dat.damage.abundance$bio1),max(dat.damage.abundance$bio1))) + ylim(c(-2,2)) +
  labs(y="Fisher's z",x="Mean annual temperature", title="b) Herbivore abundance")
graphics.off()


#diversity analysis ------------------------------------
#####################################################
dat.richness <- subset(dat.filtered, DAMAGE_OR_HERBIVORES %in% "richness" & !(VARIANCE_ESTIMATE == Inf)  & VARIANCE_ESTIMATE >= 0 & !(Z_SCORE %in% NaN)& !(Z_SCORE %in% Inf))


#random effect model without moderators
rma.richness <- rma.mv(Z_SCORE, VARIANCE_ESTIMATE,
                       random=~ STUDY_ID|CASE_ID,
                       method="REML", data=dat.richness)

forest(rma.richness)
summary(rma.richness)

#try ploting with ggplot2
plot.data.richness = rbind(data.frame("CASE_ID" = as.factor("grand mean"),
                                      "ES" =rma.richness$b,
                                      "SE"=rma.richness$se),
                           data.frame("CASE_ID" = as.factor("spaceholder"), #space holder
                                      "ES" = 0,
                                      "SE"= 0),
                           data.frame("CASE_ID" =factor(dat.richness$CASE_ID,levels=dat.richness$CASE_ID),
                                      "ES" =rma.richness$yi,
                                      "SE" =sqrt(rma.richness$vi)/sqrt(dat.richness$PLOT_NUMBER)))

plot.data.richness$TYPE1 = c("black",rep(NA,nrow(plot.data.richness) -1))
plot.data.richness$TYPE2 = c(NA,NA,rep("grey",nrow(plot.data.richness) -2))



#plot.data.richness$CASE_ID = factor(plot.data.richness$CASE_ID,levels=(c(plot.data.richness$CASE_ID)))
#plot.data.richness$TYPE = factor(plot.data.richness$TYPE,levels=rev(c(plot.data.richness$TYPE)))
levels(plot.data.richness$CASE_ID)

theme_set(theme_bw(base_size=18))

svg(filename="herbivore_richness.svg",height=3,width=6)
plot.empty = ggplot(data=(plot.data.richness)) 
plot.empty = plot.empty + geom_blank(aes(x=factor(CASE_ID),y=ES)) + coord_flip() + geom_hline(aes(x=0), lty=2,size=0.5,colour="black") 
plot.empty = plot.empty + geom_point(aes(x=factor(CASE_ID),y=ES,size=5-(plot.data.richness$SE*5),colour=(TYPE2)))
plot.empty = plot.empty + scale_colour_manual(values = c("#7e7e7e"))
plot.full = plot.empty + geom_pointrange(aes(x=1,y=rma.richness$b[1,1],ymin=rma.richness$b[1,1] - 1.96*rma.richness$se,ymax= rma.richness$b + 1.96*rma.richness$se),
                                         shape=3, size=1.2,linetype=1)
plot.full +  ylim(c(-1.6,1.6)) + theme(legend.position="none",axis.ticks=element_blank(),axis.text.y = element_text(colour="grey")) +
  labs(y="Fisher's z",x="Study case ID", title="c) Herbivore species richness")
graphics.off()




















#looking for omnibus test of moderators
rma.damage.abundance <- rma(Z_SCORE, VARIANCE_ESTIMATE,method="REML",
                            mods= ~ factor(DAMAGE_OR_HERBIVORES),
                            data=dat.damage.abundance)

rma.damage.abundance
forest(rma.damage.abundance)
permutest(rma.damage.abundance, exact=FALSE, iter=1000, progbar=TRUE,
          retpermdist=FALSE)


#damage analysis ------------------------------------
#####################################################
dat.damage.clear <- subset(dat.filtered, DAMAGE_OR_HERBIVORES %in% "damage" & !(VARIANCE_ESTIMATE == Inf)  & VARIANCE_ESTIMATE >= 0 & !(Z_SCORE %in% NaN) & !(Z_SCORE %in% Inf)
                     & !(Z_SCORE %in% NA))
#dat.damage.clear = dat.damage[row.names(na.omit(dat.damage[,c(23,24,25,28)])),]
#dat.damage.clear$NATURAL_OR_PLANTATIONS = factor(dat.damage.clear$NATURAL_OR_PLANTATIONS)
#dat.damage.clear$SURROUNDING = factor(dat.damage.clear$SURROUNDING)
#dat.damage.clear$HERBIVORE_GUILD = factor(dat.damage.clear$HERBIVORE_GUILD)


#random effect model without moderators
rma.damage.all <- rma.mv(Z_SCORE, VARIANCE_ESTIMATE,
                      mods= ~ bio1,
                      random = ~ 1|STUDY_ID,
                      data=dat.damage.clear,method="REML")
forest(rma.damage.all)
summary(rma.damage.all)


#try ploting with ggplot2
plot.data.damage = rbind(data.frame("CASE_ID" = as.factor("grand mean"),
                                    "ES" =rma.damage.all$b,
                                    "SE"=rma.damage.all$se),
                         data.frame("CASE_ID" = as.factor("spaceholder"), #space holder
                                    "ES" = 0,
                                    "SE"= 0),
                         data.frame("CASE_ID" =factor(dat.damage.clear$CASE_ID,levels=dat.damage.clear$CASE_ID),
                                    "ES" =rma.damage.all$yi,
                                    "SE" =sqrt(rma.damage.all$vi)/sqrt(dat.damage.clear$PLOT_NUMBER)))
                                   
plot.data.damage$TYPE1 = c("black",rep(NA,nrow(plot.data.damage) -1))
plot.data.damage$TYPE2 = c(NA,NA,rep("grey",nrow(plot.data.damage) -2))



#plot.data.damage$CASE_ID = factor(plot.data.damage$CASE_ID,levels=(c(plot.data.damage$CASE_ID)))
#plot.data.damage$TYPE = factor(plot.data.damage$TYPE,levels=rev(c(plot.data.damage$TYPE)))
levels(plot.data.damage$CASE_ID)

theme_set(theme_bw(base_size=18))

svg(filename="herbivore_damage.svg",height=10,width=6)
plot.empty = ggplot(data=(plot.data.damage)) 
plot.empty = plot.empty + geom_blank(aes(x=factor(CASE_ID),y=ES)) + coord_flip() + geom_hline(aes(x=0), lty=2,size=0.5,colour="black") 
plot.empty = plot.empty + geom_point(aes(x=factor(CASE_ID),y=ES,size=5-(plot.data.damage$SE*5),colour=(TYPE2)))
plot.empty = plot.empty + scale_colour_manual(values = c("grey"))
plot.full = plot.empty + geom_pointrange(aes(x=1,y=rma.damage.all$b[1,1],ymin=rma.damage.all$b[1,1] - 1.96*rma.damage.all$se,ymax= rma.damage.all$b + 1.96*rma.damage.all$se),
                             shape=3, size=1.2,linetype=1)
plot.full +  ylim(c(-1.3,1.4)) + theme(legend.position="none",axis.ticks=element_blank(),axis.text.y = element_text(colour="grey")) +
            labs(y="Fisher's z",x="Study", title="a) Herbivore damage")
graphics.off()



#publication bias
#################
#funnel plot
funnel(rma.damage.all)
#trim and fill 
trimfill(rma.damage.all,side="right",verbose=TRUE )
#funnel plot
funnel(trim.damage)
#fail safe number
fsn(Z_SCORE, VARIANCE_ESTIMATE,data=dat.damage,
    type="Rosenthal",alpha=0.05,)  


rma.resid.damage = rstandard.rma.mv(rma.damage.all)
plot(rma.resid.damage$resid ~ dat.damage$Z_SCORE)
abline(0,1)
plot(density(rma.resid.damage$resid,breaks=50))


# the same for moderator bio1
#run 1000 models with unique study_ids picked by chance
dat.damage = dat.damage.clear
set.seed(42)
sample.damage = data.frame("GRAND_MEAN"=0,"GRAND_P_VAL"=0,"GRAND_LB"=0,"GRAND_UB"=0,
                           "QM"=0, "R2"=0, "COMBI"=0,
                           "BIO1_MEAN"=0,"BIO1_P_VAL"=0,"BIO1_LB"=0,"BIO1_UB"=0,
                           "QM"=0, "R2"=0, "COMBI"=0)
sample.damage = sample.damage[-1,]
dat.damage.clear$UNIQUE_IDENTIFIER = paste(dat.damage.clear$FOREST_NAME,dat.damage.clear$FOREST_NR)
rma.without.mods = rma.mv(Z_SCORE, VARIANCE_ESTIMATE,
                         mods= ~ 1,
                         random = ~ 1|STUDY_ID,
                         data=dat.damage.unique.by.chance,method="ML")
  
i=1
for(i in 1:10000){
  dat.damage.sort.by.chance = dat.damage.clear[sample(rownames(dat.damage.clear),size = nrow(dat.damage.clear), replace = FALSE),]
  dat.damage.unique.by.chance = dat.damage.sort.by.chance[-c(which(duplicated(dat.damage.sort.by.chance$UNIQUE_IDENTIFIER))),]
  #plot(dat.damage.clear$Z_SCORE ~ dat.damage.clear$bio1)
  #points(dat.damage.unique.by.chance$Z_SCORE ~ dat.damage.unique.by.chance$bio1,col="red")
  combi = paste(as.character(sort(rownames(dat.damage.unique.by.chance))),collapse=" ")
  
  rma.temp <- rma.uni(Z_SCORE, VARIANCE_ESTIMATE,
                  mods= ~ bio1,
                  data=dat.damage.unique.by.chance,method="ML",control=list(maxiter=1000))  
  
  sample.damage = data.frame(rbind(sample.damage,data.frame(
    "GRAND_MEAN"=rma.temp$b[1],
    "GRAND_P_VAL"=rma.temp$pval[1],
    "GRAND_LB"=rma.temp$ci.lb[1],
    "GRAND_UB"=rma.temp$ci.ub[1],
    "BIO1_MEAN"=rma.temp$b[2],
    "BIO1_P_VAL"=rma.temp$pval[2],
    "BIO1_LB"=rma.temp$ci.lb[2],
    "BIO1_UB"=rma.temp$ci.ub[2],
    "QM"=rma.temp$QMp, 
    "R2" = rma.temp$R2,
    "COMBI"=paste(as.character(sort(rownames(dat.damage.unique.by.chance))),collapse=" "))))
  
  print(i)
}
sample.damage = sample.damage[order(sample.damage$GRAND_MEAN),]



svg(filename="damage_grand_mean_sampled.svg",height=3,width=8)
p = ggplot(data=sample.damage) + theme(axis.title.x= element_text(size=20),axis.title.y= element_text(size=20), axis.text.x =element_text(size=20), axis.text.y =element_text(size=20) )
p = p + geom_blank(aes(x=c(1:dim(sample.damage)[1]),y=GRAND_MEAN)) 
p = p + geom_pointrange(data=sample.damage,aes(x=c(1:dim(sample.damage)[1]),y=GRAND_MEAN,ymin= GRAND_LB, ymax = GRAND_UB),colour="dark grey") 
p = p + geom_smooth(data=sample.damage,aes(x=c(1:dim(sample.damage)[1]),y=GRAND_MEAN),colour="#666666",size=2) 
p = p + geom_hline(aes(yintercept=0),colour="black",linetype="dashed",cex=1)
p.sample.damage = p  + ylab("Grand mean effect size") + xlab("Permutations") + labs(title = "b) Herbivore damage")
p.sample.damage
graphics.off()


svg(filename="damage_bio1_sampled.svg",height=3,width=8)
p = ggplot(data=sample.damage) + theme(axis.title.x= element_text(size=20),axis.title.y= element_text(size=20), axis.text.x =element_text(size=20), axis.text.y =element_text(size=20) )
p = p + geom_blank(aes(x=c(1:dim(sample.damage)[1]),y=BIO1_MEAN)) 
p = p + geom_pointrange(data=sample.damage,aes(x=c(1:dim(sample.damage)[1]),y=BIO1_MEAN,ymin= BIO1_LB, ymax = BIO1_UB),colour="dark grey") 
p = p + geom_smooth(data=sample.damage,aes(x=c(1:dim(sample.damage)[1]),y=BIO1_MEAN),colour="#666666",size=2) 
p = p + geom_hline(aes(yintercept=0),colour="black",linetype="dashed",cex=1)
p.sample.damage.bio1 = p  + ylab("Grand mean effect size") + xlab("Permutations") + labs(title = "b) Herbivore damage")
p.sample.damage.bio1
graphics.off()



# meta-regression: WITH ALL POSSIBLE MODEL COMBINATIONS
#construct all possible formulas for meta-regression models
vars <- c("bio1","HERBIVORE_GUILD","NATURAL_OR_PLANTATIONS","SURROUNDING","bio1 * HERBIVORE_GUILD")
dat.damage.clear$NATURAL_OR_PLANTATIONS = factor(dat.damage.clear$NATURAL_OR_PLANTATIONS)
dat.damage.clear$HERBIVORE_GUILD = factor(dat.damage.clear$HERBIVORE_GUILD)
dat.damage.clear$SURROUNDING = factor(dat.damage.clear$SURROUNDING)
dat.damage.clear$FORMER_LAND_USE = factor(dat.damage.clear$FORMER_LAND_USE)

indexes<-unique(apply(combinations(length(vars), length(vars), repeats=T), 1, unique))
gen.form<-function(x) as.formula(paste('~ ',paste( vars[x],collapse='+')))
formulas<-lapply(indexes, gen.form)
formulas[[length(formulas)+1]] = as.formula( ~ 1)

all.models.temp = data.frame("NR"=0,"N"=0,"FORMULA"=0,"AIC"=0,"AICC"=0,"BIC"=0)
all.models.temp$bio1 = 0
all.models.temp$NATURAL_OR_PLANTATIONSplantation = 0
all.models.temp$NATURAL_OR_PLANTATIONSforest = 0
all.models.temp$HERBIVORE_GUILDmammals = 0
all.models.temp$HERBIVORE_GUILDinsects = 0
all.models.temp$SURROUNDINGgrassland = 0
all.models.temp$SURROUNDINGforest = 0
all.models.temp$SURROUNDINGmosaic = 0
all.models.temp$SURROUNDINGplantation = 0
all.models.temp$bio1HERBIVORES = 0
all.models.temp$SE_bio1 = 0
all.models.temp$SE_NATURAL_OR_PLANTATIONSplantation = 0
all.models.temp$SE_NATURAL_OR_PLANTATIONSforest = 0
all.models.temp$SE_HERBIVORE_GUILDmammals = 0
all.models.temp$SE_HERBIVORE_GUILDinsects = 0
all.models.temp$SE_SURROUNDINGgrassland = 0
all.models.temp$SE_SURROUNDINGforest = 0
all.models.temp$SE_SURROUNDINGmosaic = 0
all.models.temp$SE_SURROUNDINGplantation = 0
all.models.temp$SE_bio1HERBIVORES = 0



#normalize data
#dat.damage.clear$bio1 = dat.damage.clear$bio1 / max(dat.damage.clear$bio1)
                                                    
                                                    
i=1
for(i in 1:length(formulas)){
  rma.fit <- rma.mv(Z_SCORE, VARIANCE_ESTIMATE,
                    mods= formulas[[i]],
                    random = ~ 1|STUDY_ID,
                    data=dat.damage.clear,method="ML")
  
  all.models.temp = rbind(all.models.temp,data.frame("NR"=i,
                                                     "N"=rma.fit$k,
                                                     "FORMULA"=paste("~ ",as.character(formulas[[i]][2]),sep=""),
                                                     "AIC"=fitstats(rma.fit)[3],
                                                     "AICC"=fitstats(rma.fit)[5],
                                                     "BIC"=fitstats(rma.fit)[4],
                                                     "bio1"=  if(length(which(names(coef(rma.fit)) %in% "bio1")) >=1){coef(rma.fit)[which(names(coef(rma.fit)) %in% "bio1")]}else(0),
                                                     "NATURAL_OR_PLANTATIONSplantation" = if(length(which(names(coef(rma.fit)) %in% "NATURAL_OR_PLANTATIONSplantation")) >=1){coef(rma.fit)[which(names(coef(rma.fit)) %in% "NATURAL_OR_PLANTATIONSplantation")]}else(0),
                                                     "NATURAL_OR_PLANTATIONSforest" = if(length(which(names(coef(rma.fit)) %in% "NATURAL_OR_PLANTATIONSforest")) >=1){coef(rma.fit)[which(names(coef(rma.fit)) %in% "NATURAL_OR_PLANTATIONSforest")]}else(0),
                                                     "HERBIVORE_GUILDmammals" =if(length(which(names(coef(rma.fit)) %in% "HERBIVORE_GUILDmammals")) >=1){coef(rma.fit)[which(names(coef(rma.fit)) %in% "HERBIVORE_GUILDmammals")]}else(0),
                                                     "HERBIVORE_GUILDinsects" = if(length(which(names(coef(rma.fit)) %in% "HERBIVORE_GUILDinsects")) >=1){coef(rma.fit)[which(names(coef(rma.fit)) %in% "HERBIVORE_GUILDinsects")]}else(0),
                                                     "SURROUNDINGgrassland" = if(length(which(names(coef(rma.fit)) %in% "SURROUNDINGgrassland")) >=1){coef(rma.fit)[which(names(coef(rma.fit)) %in% "SURROUNDINGgrassland")]}else(0),
                                                     "SURROUNDINGforest" =  if(length(which(names(coef(rma.fit)) %in% "SURROUNDINGforest")) >=1){coef(rma.fit)[which(names(coef(rma.fit)) %in% "SURROUNDINGforest")]}else(0),
                                                     "SURROUNDINGmosaic" = if(length(which(names(coef(rma.fit)) %in% "SURROUNDINGmosaic")) >=1){coef(rma.fit)[which(names(coef(rma.fit)) %in% "SURROUNDINGmosaic")]}else(0),
                                                     "SURROUNDINGplantation" = if(length(which(names(coef(rma.fit)) %in% "SURROUNDINGplantation")) >=1){coef(rma.fit)[which(names(coef(rma.fit)) %in% "SURROUNDINGplantation")]}else(0),
                                                     "bio1HERBIVORES" = if(length(which(names(coef(rma.fit)) %in% "bio1:HERBIVORE_GUILDmammals")) >=1){coef(rma.fit)[which(names(coef(rma.fit)) %in% "bio1:HERBIVORE_GUILDmammals")]}else(0),
                                                     
                                                        "SE_bio1"=  if(length(which(names(coef(rma.fit)) %in% "bio1")) >=1){rma.fit$se[which(names(coef(rma.fit)) %in% "bio1")]}else(0),
                                                        "SE_NATURAL_OR_PLANTATIONSplantation" = if(length(which(names(coef(rma.fit)) %in% "NATURAL_OR_PLANTATIONSplantation")) >=1){rma.fit$se[which(names(coef(rma.fit)) %in% "NATURAL_OR_PLANTATIONSplantation")]}else(0),
                                                        "SE_NATURAL_OR_PLANTATIONSforest" = if(length(which(names(coef(rma.fit)) %in% "NATURAL_OR_PLANTATIONSforest")) >=1){rma.fit$se[which(names(coef(rma.fit)) %in% "NATURAL_OR_PLANTATIONSforest")]}else(0),
                                                        "SE_HERBIVORE_GUILDmammals" =if(length(which(names(coef(rma.fit)) %in% "HERBIVORE_GUILDmammals")) >=1){rma.fit$se[which(names(coef(rma.fit)) %in% "HERBIVORE_GUILDmammals")]}else(0),
                                                        "SE_HERBIVORE_GUILDinsects" = if(length(which(names(coef(rma.fit)) %in% "HERBIVORE_GUILDinsects")) >=1){rma.fit$se[which(names(coef(rma.fit)) %in% "HERBIVORE_GUILDinsects")]}else(0),
                                                        "SE_SURROUNDINGgrassland" = if(length(which(names(coef(rma.fit)) %in% "SURROUNDINGgrassland")) >=1){rma.fit$se[which(names(coef(rma.fit)) %in% "SURROUNDINGgrassland")]}else(0),
                                                        "SE_SURROUNDINGforest" =  if(length(which(names(coef(rma.fit)) %in% "SURROUNDINGforest")) >=1){rma.fit$se[which(names(coef(rma.fit)) %in% "SURROUNDINGforest")]}else(0),
                                                        "SE_SURROUNDINGmosaic" = if(length(which(names(coef(rma.fit)) %in% "SURROUNDINGmosaic")) >=1){rma.fit$se[which(names(coef(rma.fit)) %in% "SURROUNDINGmosaic")]}else(0),
                                                        "SE_SURROUNDINGplantation" = if(length(which(names(coef(rma.fit)) %in% "SURROUNDINGplantation")) >=1){rma.fit$se[which(names(coef(rma.fit)) %in% "SURROUNDINGplantation")]}else(0),
                                                        "SE_bio1HERBIVORES" = if(length(which(names(coef(rma.fit)) %in% "bio1:HERBIVORE_GUILDmammals")) >=1){rma.fit$se[which(names(coef(rma.fit)) %in% "bio1:HERBIVORE_GUILDmammals")]}else(0)))

                                                        
}

#remove first line
all.models.damage = all.models.temp[-1,]
#remove all factors in the model to which the remaining fixed factors refer
all.models.damage= all.models.damage[,-c(9,11,13,19,21,23)]
all.models.damage = all.models.damage[order(all.models.damage$AICC),]
all.models.damage

#remove duplicates from interaction term
dupli = c()
i=1
for( i in 1:dim(all.models.damage)[1]){
    if(length(which(textcnt(all.models.damage[i,3],n=1L,method="string") == 2) >=1) >=1){dupli = c(dupli,i)}
}
all.models.damage = all.models.damage[-c(dupli),]

#calculate the relative evidence weight for each model
all.models.damage$LIKELIHOOD = exp(-(all.models.damage$AICC - all.models.damage$AICC[1]) / 2)
all.models.damage$AKAIKE_WEIGHT = all.models.damage$LIKELIHOOD / sum(all.models.damage$LIKELIHOOD)


#calculate weighted estimates and variances for all factors
data.frame(names(all.models.damage))
all.models.damage[,c(14:20)] = (all.models.damage[,c(14:20)] * sqrt(all.models.damage$N))^2


importances.damage = data.frame("bio1" = sum(all.models.damage$AKAIKE_WEIGHT[which(!(all.models.damage$bio1 == 0))]),
                               "NATURAL_OR_PLANTATIONS" = sum(all.models.damage$AKAIKE_WEIGHT[which(!(all.models.damage$NATURAL_OR_PLANTATIONSplantation == 0))]),
                               "HERBIVORE_GUILD" = sum(all.models.damage$AKAIKE_WEIGHT[which(!(all.models.damage$HERBIVORE_GUILDmammals == 0))]),
                               "SURROUNDING" = sum(all.models.damage$AKAIKE_WEIGHT[which(!(all.models.damage$SURROUNDINGgrassland == 0))]),
                               "bio1HERBIVORES" = sum(all.models.damage$AKAIKE_WEIGHT[which(!(all.models.damage$bio1HERBIVORES == 0))]))
                               

estimates.damage = data.frame("bio1" = sum(all.models.damage$AKAIKE_WEIGHT * all.models.damage$bio1),
                              "NATURAL_OR_PLANTATIONSplantation" = sum(all.models.damage$AKAIKE_WEIGHT * all.models.damage$NATURAL_OR_PLANTATIONSplantation),
                              "HERBIVORE_GUILDmammals" = sum(all.models.damage$AKAIKE_WEIGHT * all.models.damage$HERBIVORE_GUILDmammals),
                              "SURROUNDINGgrassland" = sum(all.models.damage$AKAIKE_WEIGHT * all.models.damage$SURROUNDINGgrassland),
                              "SURROUNDINGmosaic" = sum(all.models.damage$AKAIKE_WEIGHT * all.models.damage$SURROUNDINGmosaic),
                              "SURROUNDINGplantation" = sum(all.models.damage$AKAIKE_WEIGHT * all.models.damage$SURROUNDINGplantation),
                              "bio1HERBIVORES" = sum(all.models.damage$AKAIKE_WEIGHT * all.models.damage$bio1HERBIVORES))


#compute the variances
variances.damage = data.frame("bio1" = sum(all.models.damage$AKAIKE_WEIGHT*(all.models.damage$SE_bio1+((all.models.damage$bio1 - mean(all.models.damage$bio1))^2))),
                              "NATURAL_OR_PLANTATIONSplantation" = sum(all.models.damage$AKAIKE_WEIGHT*(all.models.damage$SE_NATURAL_OR_PLANTATIONSplantation+((all.models.damage$NATURAL_OR_PLANTATIONSplantation - mean(all.models.damage$NATURAL_OR_PLANTATIONSplantation))^2))),
                              "HERBIVORE_GUILDmammals" = sum(all.models.damage$AKAIKE_WEIGHT*(all.models.damage$SE_HERBIVORE_GUILDmammals+((all.models.damage$HERBIVORE_GUILDmammals - mean(all.models.damage$HERBIVORE_GUILDmammals))^2))),
                              "SURROUNDINGgrassland" = sum(all.models.damage$AKAIKE_WEIGHT*(all.models.damage$SE_SURROUNDINGgrassland+((all.models.damage$SURROUNDINGgrassland - mean(all.models.damage$SURROUNDINGgrassland)^2)))),
                              "SURROUNDINGmosaic" = sum(all.models.damage$AKAIKE_WEIGHT*(all.models.damage$SE_SURROUNDINGmosaic+((all.models.damage$SURROUNDINGmosaic - mean(all.models.damage$SURROUNDINGmosaic))^2))),
                              "SURROUNDINGplantation" = sum(all.models.damage$AKAIKE_WEIGHT*(all.models.damage$SE_SURROUNDINGplantation+((all.models.damage$SURROUNDINGplantation - mean(all.models.damage$SURROUNDINGplantation))^2))),
                              "bio1HERBIVORES" = sum(all.models.damage$AKAIKE_WEIGHT*(all.models.damage$SE_bio1HERBIVORES+((all.models.damage$bio1HERBIVORES - mean(all.models.damage$bio1HERBIVORES))^2))))


head(all.models.damage[,c(3,5)])

round(importances.damage,3)
fisher.to.r(estimates.damage)
fisher.to.r(variances.damage)


#model fit for the best model
best.model.damage = rma.mv(Z_SCORE, VARIANCE_ESTIMATE,
                    mods= ~ bio1,
                    random = ~ 1|STUDY_ID,
                    data=dat.damage.clear,method="ML")

best.model.damage.no.mods = rma.mv(Z_SCORE, VARIANCE_ESTIMATE,
                           #mods= ~ factor(ECOREGION),
                           random = ~ 1|STUDY_ID,
                           data=dat.damage.clear,method="ML")
best.model.damage

((best.model.damage.no.mods$sigma2 - best.model.damage$sigma2)/best.model.damage.no.mods$sigma2)^2

theme_set(theme_bw(base_size=18))

svg(filename="damage_bio_correlation.svg",height=5,width=6)
p = ggplot() + theme(legend.position="none",axis.ticks = element_blank())
p = p + geom_point(aes(x=dat.damage.clear$bio1,y=fisher.to.r(dat.damage.clear$Z_SCORE),size=(1/dat.damage.clear$VARIANCE_ESTIMATE/max(dat.damage.clear$VARIANCE_ESTIMATE))),colour="#666666")
p= p + scale_size(range=c(3,8))  + ylim(c(-1,1)) + xlim(c(0,30))
p = p + geom_line(aes(x=dat.damage.clear$bio1,y=fitted(best.model.damage)))
p = p + xlab("Mean annual temperature") + ylab("Correlation coefficient") + ggtitle("a) Herbivore damage")
p
graphics.off()


plot(Z_SCORE ~ bio1, data=dat.damage.clear)
abline(lm(fitted(best.model.damage) ~ dat.damage.clear$bio1))
plot(residuals(best.model.damage) ~ dat.damage.clear$bio1)
shapiro.test(residuals(best.model.damage))


#meta-regression with random effect model, method is ML because REML is not comparable between different models
rma.fit.damage1 <- rma.mv(Z_SCORE, VARIANCE_ESTIMATE,
                          mods= ~ bio1 + NATURAL_OR_PLANTATIONS + AGE_FOREST + factor(HERBIVORE_GUILD) + factor(SURROUNDING),
                          random= ~ 1 | FOREST_NAME,
                          data=dat.damage.clear,method="REML")
summary(rma.fit.damage1)

#threre is no call stored in the rma.fit object, so I have to reduce complexity without the update() function
#exclude AGE_FOREST
rma.fit.damage2 <- rma.mv(Z_SCORE, VARIANCE_ESTIMATE,
                          mods= ~ bio1 + NATURAL_OR_PLANTATIONS  + factor(HERBIVORE_GUILD) + factor(SURROUNDING),
                          random= ~ 1 | FOREST_NAME,
                          data=dat.damage.clear,method="REML")
summary(rma.fit.damage2)

#exclude bio1
rma.fit.damage3 <- rma.mv(Z_SCORE, VARIANCE_ESTIMATE,
                          mods= ~ bio1 + NATURAL_OR_PLANTATIONS + factor(HERBIVORE_GUILD),
                          random= ~ 1 | FOREST_NAME,
                          data=dat.damage.clear,method="REML")
summary(rma.fit.damage3)

#exclude NATURAL_OR_PLANTATIONS
rma.fit.damage3 <- rma.mv(Z_SCORE, VARIANCE_ESTIMATE,
                          mods= ~ bio1 + factor(HERBIVORE_GUILD),
                          random= ~ 1 | FOREST_NAME,
                          data=dat.damage.clear,method="REML")
summary(rma.fit.damage3)

#exclude HERBIVORE_GUILD
rma.fit.damage3 <- rma.mv(Z_SCORE, VARIANCE_ESTIMATE,
                          mods= ~  bio1,
                          random= ~ 1 | FOREST_NAME,
                          data=dat.damage.clear,method="REML")
summary(rma.fit.damage3)




rma.fit.damage3[[1]] = fisher.to.r(rma.fit.damage3[[1]])
rma.fit.damage3[[2]] = fisher.to.r(rma.fit.damage3[[1]])
rma.fit.damage3[[5]] = fisher.to.r(rma.fit.damage3[[5]])
rma.fit.damage3[[6]] = fisher.to.r(rma.fit.damage3[[6]])

summary(rma.fit.damage3)




#try ploting with ggplot2
theme_set(theme_bw(base_size=18))


plot.data.damage = cbind(dat.damage,ES=rma.damage.all$yi,SE=sqrt(rma.damage.all$vi)/sqrt(dat.damage$PLOT_NUMBER))
plot.data.damage.forest = plot.data.damage[plot.data.damage$NATURAL_OR_PLANTATIONS %in% "forest",]
plot.data.damage.plantation = plot.data.damage[plot.data.damage$NATURAL_OR_PLANTATIONS %in% "plantation",]
plot.data.damage.forest = plot.data.damage.forest[order(plot.data.damage.forest$Z_SCORE),]
plot.data.damage.plantation = plot.data.damage.plantation[order(plot.data.damage.plantation$Z_SCORE),]
plot.data.damage.all = plot.data.damage[order(plot.data.damage$Z_SCORE),]



#over-all effect  line
rma.temp.damage.all <- rma(yi=plot.data.damage.all$Z_SCORE, vi=plot.data.damage.all$VARIANCE_ESTIMATE,method="EB")
summary.line = data.frame(CASE_ID = "grand mean", STUDY_ID =NA,FOREST_NR=NA,FOREST_NAME="summary",AUTHOR="summary",
                          YEAR=NA,PEER_OR_GREY="summary",UNIQUE_ID="summary",DAMAGE_OR_HERBIVORES="summary",
                          UNITY="summary",HERBIVORE_GUILD="summary",HERBIVORE_SPECIES="summary",TREE_DIV_METRIC="summary",
                          SAMPLING_ON="summary",BETA_ESTIMATE=NA,R_SQUARE=NA,
                          R_SQUARE_ADJUSTED=NA,P_VALUE=NA,PLOT_NUMBER=NA,
                          MIN_RICH=NA,MAX_RICH=NA,MAX_DAMAGE_FOR_STANDARD=NA,
                          NATURAL_OR_PLANTATIONS="summary",SURROUNDING="summary",AGE_FOREST=NA,FORMER_LAND_USE="summary",
                          GPS_DEZ_X=NA,GPS_DEZ_Y=NA,ALTITUDE=NA,ECOREGION="summary",
                          FOCAL_SPECIES="summary",JUST_FOCAL="summary",PLOT_SIZE=NA,X=NA,
                          R_VALUE = NA, Z_SCORE=NA,VARIANCE_ESTIMATE=NA,
                          ES=rma.temp.damage.all$b,SE=rma.temp.damage.all$se)              

#forest effect size
#ma.temp.damage.forest <- rma(yi=plot.data.damage.forest$Z_SCORE, vi=plot.data.damage.forest$VARIANCE_ESTIMATE,method="EB")
#plot.data.damage.forest = rbind(plot.data.damage.forest,data.frame(CASE_ID = "forest mean", STUDY_ID =NA,FOREST_NR=NA,FOREST_NAME="summary",AUTHOR="summary",
#                                                             YEAR=NA,PEER_OR_GREY="summary",UNIQUE_ID="summary",DAMAGE_OR_HERBIVORES="summary",
#                                                             UNITY="summary",HERBIVORE_GUILD="summary",HERBIVORE_SPECIES="summary",TREE_DIV_METRIC="summary",
#                                                             SAMPLING_ON="summary",BETA_ESTIMATE=NA,R_SQUARE=NA,
#                                                             R_SQUARE_ADJUSTED=NA,P_VALUE=NA,PLOT_NUMBER=NA,
#                                                             MIN_RICH=NA,MAX_RICH=NA,MAX_DAMAGE_FOR_STANDARD=NA,
#                                                             NATURAL_OR_PLANTATIONS="summary",SURROUNDING="summary",AGE_FOREST=NA,FORMER_LAND_USE="summary",
#                                                             GPS_DEZ_X=NA,GPS_DEZ_Y=NA,ALTITUDE=NA,ECOREGION="summary",
#                                                             FOCAL_SPECIES="summary",JUST_FOCAL="summary",PLOT_SIZE=NA,X=NA,
#                                                             R_VALUE = NA, Z_SCORE=NA,VARIANCE_ESTIMATE=NA,
#                                                             ES=rma.temp.damage.forest$b,SE=rma.temp.damage.forest$se))              
#dummy-line
#dummy.line = data.frame(CASE_ID = " ", STUDY_ID =NA,FOREST_NR=NA,FOREST_NAME="dummy",AUTHOR="dummy",
#                        YEAR=NA,PEER_OR_GREY="dummy",UNIQUE_ID="dummy",DAMAGE_OR_HERBIVORES="dummy",
#                        UNITY="dummy",HERBIVORE_GUILD="dummy",HERBIVORE_SPECIES="dummy",TREE_DIV_METRIC="dummy",
#                        SAMPLING_ON="dummy",BETA_ESTIMATE=NA,R_SQUARE=NA,
#                        R_SQUARE_ADJUSTED=NA,P_VALUE=NA,PLOT_NUMBER=NA,
#                        MIN_RICH=NA,MAX_RICH=NA,MAX_DAMAGE_FOR_STANDARD=NA,
#                        NATURAL_OR_PLANTATIONS="dummy",SURROUNDING="dummy",AGE_FOREST=NA,FORMER_LAND_USE="dummy",
#                        GPS_DEZ_X=NA,GPS_DEZ_Y=NA,ALTITUDE=NA,ECOREGION="dummy",
#                        FOCAL_SPECIES="dummy",JUST_FOCAL="dummy",PLOT_SIZE=NA,X=NA,
#                        R_VALUE = NA, Z_SCORE=NA,VARIANCE_ESTIMATE=NA,
#                        ES=0,SE=0)              

#plantations effect size
#rma.temp.damage.plantation <- rma(yi=plot.data.damage.plantation$Z_SCORE, vi=plot.data.damage.plantation$VARIANCE_ESTIMATE,method="EB")
#plot.data.damage.plantation = rbind(plot.data.damage.plantation,data.frame(CASE_ID = "plantation mean", STUDY_ID =NA,FOREST_NR=NA,FOREST_NAME="summary",AUTHOR="summary",
#                                                                   YEAR=NA,PEER_OR_GREY="summary",UNIQUE_ID="summary",DAMAGE_OR_HERBIVORES="summary",
#                                                                   UNITY="summary",HERBIVORE_GUILD="summary",HERBIVORE_SPECIES="summary",TREE_DIV_METRIC="summary",
#                                                                   SAMPLING_ON="summary",BETA_ESTIMATE=NA,R_SQUARE=NA,
#                                                                   R_SQUARE_ADJUSTED=NA,P_VALUE=NA,PLOT_NUMBER=NA,
#                                                                   MIN_RICH=NA,MAX_RICH=NA,MAX_DAMAGE_FOR_STANDARD=NA,
#                                                                   NATURAL_OR_PLANTATIONS="summary",SURROUNDING="summary",AGE_FOREST=NA,FORMER_LAND_USE="summary",
#                                                                   GPS_DEZ_X=NA,GPS_DEZ_Y=NA,ALTITUDE=NA,ECOREGION="summary",
#                                                                   FOCAL_SPECIES="summary",JUST_FOCAL="summary",PLOT_SIZE=NA,X=NA,
#                                                                   R_VALUE = NA, Z_SCORE=NA,VARIANCE_ESTIMATE=NA,
#                                                                   ES=rma.temp.damage.plantation$b,SE=rma.temp.damage.plantation$se))              

#creat CASE_ID with unique names and CASE_ID2 to code 
#plot.data.damage.separated = rbind(plot.data.damage.forest,dummy.line,plot.data.damage.plantation,dummy.line,summary.line)
#plot.data.damage.separated$CASE_ID2 = c(rep("Forest",10),NA,NA,rep("Plantation",33),NA,NA,NA)
#plot.data.damage.separated$CASE_ID3 = c(rep(NA,10),"Study",NA,rep(NA,33),"Study",NA,"Summary")

#plot.data.damage.separated$CASE_ID[which(plot.data.damage.separated$CASE_ID %in% " ")][2] = "  "
#plot.data.damage.separated$CASE_ID = factor(plot.data.damage.separated$CASE_ID, levels=rev(c(as.character(plot.data.damage.separated$CASE_ID))))


#build from scratch
#svg(filename="herbivore_damage.svg",height=9,width=6)
#plot.empty = ggplot(data=plot.data.damage.separated) 
#plot.empty = plot.empty + geom_blank(aes(x=CASE_ID,y=ES)) + coord_flip() + geom_hline(aes(x=0), lty=2,size=1,colour="grey") 
#plot.full = plot.empty + geom_point(data=plot.data.damage.separated,aes(x=factor(CASE_ID),y=ES,size=5-(plot.data.damage.separated$SE*5), colour=plot.data.damage.separated$CASE_ID2)) 
#plot.full = plot.full + geom_pointrange(data=plot.data.damage.separated,aes(x=factor(CASE_ID),y=ES,ymin=ES-1.96*SE,ymax=ES+1.96*SE, colour=plot.data.damage.separated$CASE_ID3),size=1.2,pch=3) + scale_colour_manual(values=c("grey","black","black","black"))
#plot.full + ylim(c(-1.5, 1.5)) + theme(legend.position="none") + labs(y="Fisher's z-transformed correlation coefficient",x="Study",title ="(a) Herbivore Damage")
#variance estimate lines are out of plot borders
#graphics.off()
#dev.copy2eps(file="herbivore_damage.eps", height=9, width=6)
#dev.copy2pdf(file="herbivore_damage.pdf", height=9, width=6)
#svg(filename="herbivore_damage.svg")












#abundance analysis -----------------------------------
#####################################################
dat.abundance.clear <- subset(dat.filtered, DAMAGE_OR_HERBIVORES %in% "abundance" & !(VARIANCE_ESTIMATE == Inf)  & VARIANCE_ESTIMATE >= 0 & !(Z_SCORE %in% NaN)
                        & !(Z_SCORE %in% Inf))
#dat.abundance.clear = dat.abundance[row.names(na.omit(dat.abundance[,c(23,24,25,28)])),]
#dat.abundance.clear$NATURAL_OR_PLANTATIONS = factor(dat.abundance.clear$NATURAL_OR_PLANTATIONS)
#dat.abundance.clear$SURROUNDING = factor(dat.abundance.clear$SURROUNDING)
#dat.abundance.clear$HERBIVORE_GUILD = factor(dat.abundance.clear$HERBIVORE_GUILD)


#random effect model without moderators
rma.abundance.all <- rma.mv(Z_SCORE, VARIANCE_ESTIMATE,
                         mods = ~ bio1,
                         random = ~ 1|STUDY_ID,
                         method="REML",data=dat.abundance.clear)

forest(rma.abundance.all)
summary(rma.abundance.all)



#try ploting with ggplot2

theme_set(theme_bw(base_size=18))

plot.data.abundance = rbind(data.frame("CASE_ID" = as.factor("grand mean"),
                                    "ES" =rma.abundance.all$b,
                                    "SE"=rma.abundance.all$se),
                         data.frame("CASE_ID" = as.factor("spaceholder"), #space holder
                                    "ES" = 0,
                                    "SE"= 0),
                         data.frame("CASE_ID" =factor(dat.abundance.clear$CASE_ID,levels=dat.abundance.clear$CASE_ID),
                                    "ES" =rma.abundance.all$yi,
                                    "SE" =sqrt(rma.abundance.all$vi)/sqrt(dat.abundance.clear$PLOT_NUMBER)))

plot.data.abundance$TYPE1 = c("black",rep(NA,nrow(plot.data.abundance) -1))
plot.data.abundance$TYPE2 = c(NA,NA,rep("grey",nrow(plot.data.abundance) -2))



#plot.data.damage$CASE_ID = factor(plot.data.damage$CASE_ID,levels=(c(plot.data.damage$CASE_ID)))
#plot.data.damage$TYPE = factor(plot.data.damage$TYPE,levels=rev(c(plot.data.damage$TYPE)))
levels(plot.data.abundancedamage$CASE_ID)

theme_set(theme_bw(base_size=18))

svg(filename="herbivore_abundance.svg",height=7,width=6)
plot.empty = ggplot(data=(plot.data.abundance)) 
plot.empty = plot.empty + geom_blank(aes(x=factor(CASE_ID),y=ES)) + coord_flip() + geom_hline(aes(x=0), lty=2,size=0.5,colour="black") 
plot.empty = plot.empty + geom_point(aes(x=factor(CASE_ID),y=ES,size=5-(plot.data.abundance$SE*5),colour=(TYPE2)))
plot.empty = plot.empty + scale_colour_manual(values = c("grey"))
plot.full = plot.empty + geom_pointrange(aes(x=1,y=rma.abundance.all$b[1,1],ymin=rma.abundance.all$b[1,1] - 1.96*rma.abundance.all$se,ymax= rma.abundance.all$b + 1.96*rma.abundance.all$se),
                                         shape=3, size=1.2,linetype=1)
plot.full +  ylim(c(-2,2)) + theme(legend.position="none",axis.ticks=element_blank(),axis.text.y = element_text(colour="grey")) +
  labs(y="Fisher's z",x="Study", title="b) Herbivore abundance")
graphics.off()



#publication bias
#################
#funnel plot
funnel(rma.abundance.all)
#trim and fill 
trim.abundance = trimfill(rma.abundance.all,verbose=TRUE )
trim.abundance
#funnel plot
funnel(trim.abundance)
#fail safe number
fsn(Z_SCORE, VARIANCE_ESTIMATE,data=dat.abundance,
    type="Rosenthal",alpha=0.05,)



rma.resid.abundance = rstandard.rma.mv(rma.abundance.all)
plot(rma.resid.abundance$resid ~ dat.abundance$Z_SCORE)
abline(0,1)
plot(density(rma.resid.abundance$resid,breaks=50))


#new approach:
#sample one case study per forest_id  (1000x) and display the distribution of the grand mean effect size

# create data set
#run 1000 models with unique study_ids picked by chance
dat.abundance = dat.abundance.clear
set.seed(42)
sample.abundance = data.frame("GRAND_MEAN"=0,"GRAND_P_VAL"=0,"GRAND_LB"=0,"GRAND_UB"=0,
                           "QM"=0, "R2"=0, "COMBI"=0,
                           "BIO1_MEAN"=0,"BIO1_P_VAL"=0,"BIO1_LB"=0,"BIO1_UB"=0,
                           "QM"=0, "R2"=0, "COMBI"=0)
sample.abundance = sample.abundance[-1,]
dat.abundance.clear$UNIQUE_IDENTIFIER = paste(dat.abundance.clear$FOREST_NAME,dat.abundance.clear$FOREST_NR)
rma.without.mods = rma.mv(Z_SCORE, VARIANCE_ESTIMATE,
                          mods= ~ 1,
                          random = ~ 1|STUDY_ID,
                          data=dat.abundance.unique.by.chance,method="ML")

i=1
for(i in 1:10000){
  dat.abundance.sort.by.chance = dat.abundance.clear[sample(rownames(dat.abundance.clear),size = nrow(dat.abundance.clear), replace = FALSE),]
  dat.abundance.unique.by.chance = dat.abundance.sort.by.chance[-c(which(duplicated(dat.abundance.sort.by.chance$UNIQUE_IDENTIFIER))),]
  #plot(dat.abundance.clear$Z_SCORE ~ dat.abundance.clear$bio1)
  #points(dat.abundance.unique.by.chance$Z_SCORE ~ dat.abundance.unique.by.chance$bio1,col="red")
  combi = paste(as.character(sort(rownames(dat.abundance.unique.by.chance))),collapse=" ")
  
  rma.temp <- rma.uni(Z_SCORE, VARIANCE_ESTIMATE,
                      mods= ~ bio1,
                      data=dat.abundance.unique.by.chance,method="DL",control=list(maxiter=1000))  
  
  sample.abundance = data.frame(rbind(sample.abundance,data.frame(
    "GRAND_MEAN"=rma.temp$b[1],
    "GRAND_P_VAL"=rma.temp$pval[1],
    "GRAND_LB"=rma.temp$ci.lb[1],
    "GRAND_UB"=rma.temp$ci.ub[1],
    "BIO1_MEAN"=rma.temp$b[2],
    "BIO1_P_VAL"=rma.temp$pval[2],
    "BIO1_LB"=rma.temp$ci.lb[2],
    "BIO1_UB"=rma.temp$ci.ub[2],
    "QM"=rma.temp$QMp, 
    "R2" = rma.temp$R2,
    "COMBI"=paste(as.character(sort(rownames(dat.abundance.unique.by.chance))),collapse=" "))))
  
  print(i)
}
sample.abundance = sample.abundance[order(sample.abundance$GRAND_MEAN),]



svg(filename="abundance_grand_mean_sampled.svg",height=3,width=8)
p = ggplot(data=sample.abundance) + theme(axis.title.x= element_text(size=20),axis.title.y= element_text(size=20), axis.text.x =element_text(size=20), axis.text.y =element_text(size=20) )
p = p + geom_blank(aes(x=c(1:dim(sample.abundance)[1]),y=GRAND_MEAN)) 
p = p + geom_pointrange(data=sample.abundance,aes(x=c(1:dim(sample.abundance)[1]),y=GRAND_MEAN,ymin= GRAND_LB, ymax = GRAND_UB),colour="dark grey") 
p = p + geom_smooth(data=sample.abundance,aes(x=c(1:dim(sample.abundance)[1]),y=GRAND_MEAN),colour="#666666",size=2) 
p = p + geom_hline(aes(yintercept=0),colour="black",linetype="dashed",cex=1)
p.sample.abundance = p  + ylab("Grand mean effect size") + xlab("Permutations") + labs(title = "b) Herbivore abundance")
p.sample.abundance
graphics.off()


svg(filename="abundance_bio1_sampled.svg",height=3,width=8)
p = ggplot(data=sample.abundance) + theme(axis.title.x= element_text(size=20),axis.title.y= element_text(size=20), axis.text.x =element_text(size=20), axis.text.y =element_text(size=20) )
p = p + geom_blank(aes(x=c(1:dim(sample.abundance)[1]),y=BIO1_MEAN)) 
p = p + geom_pointrange(data=sample.abundance,aes(x=c(1:dim(sample.abundance)[1]),y=BIO1_MEAN,ymin= BIO1_LB, ymax = BIO1_UB),colour="dark grey") 
p = p + geom_smooth(data=sample.abundance,aes(x=c(1:dim(sample.abundance)[1]),y=BIO1_MEAN),colour="#666666",size=2) 
p = p + geom_hline(aes(yintercept=0),colour="black",linetype="dashed",cex=1)
p.sample.abundance.bio1 = p  + ylab("bio1 estimate") + xlab("Permutations") + labs(title = "b) Herbivore abundance")
p.sample.damage.bio1
graphics.off()





#remove the UNIQUE_IDENTIFIER
dat.abundance = dat.abundance[,-(length(names(dat.abundance)))]


#influence of climate
#----------------------
rma.abundance.bio1 <- rma(Z_SCORE, VARIANCE_ESTIMATE,
                       mods = ~ bio1,
                       data=dat.abundance.clear,method="REML")

rma.abundance.bio1


# meta-regression: WITH ALL POSSIBLE MODEL COMBINATIONS
#construct all possible formulas for meta-regression models
vars <- c("bio1","NATURAL_OR_PLANTATIONS","SURROUNDING")
dat.abundance.clear$NATURAL_OR_PLANTATIONS = factor(dat.abundance.clear$NATURAL_OR_PLANTATIONS)
dat.abundance.clear$HERBIVORE_GUILD = factor(dat.abundance.clear$HERBIVORE_GUILD)
dat.abundance.clear$FORMER_LAND_USE = factor(dat.abundance.clear$FORMER_LAND_USE)
dat.abundance.clear$SURROUNDING = factor(dat.abundance.clear$SURROUNDING)
indexes<-unique(apply(combinations(length(vars), length(vars), repeats=T), 1, unique))
gen.form<-function(x) as.formula(paste('~',paste( vars[x],collapse='+')))
formulas<-lapply(indexes, gen.form)
formulas[[length(formulas)+1]] = as.formula( ~ 1)


all.models.temp = data.frame("NR"=0,"N"=0,"FORMULA"=0,"AIC"=0,"AICC"=0,"BIC"=0)
all.models.temp$bio1 = 0
all.models.temp$NATURAL_OR_PLANTATIONSplantation = 0
all.models.temp$NATURAL_OR_PLANTATIONSforest = 0
all.models.temp$SURROUNDINGgrassland = 0
all.models.temp$SURROUNDINGforest = 0
all.models.temp$SURROUNDINGmosaic = 0
all.models.temp$SURROUNDINGplantation = 0
all.models.temp$SE_bio1 = 0
all.models.temp$SE_NATURAL_OR_PLANTATIONSplantation = 0
all.models.temp$SE_NATURAL_OR_PLANTATIONSforest = 0
all.models.temp$SE_SURROUNDINGgrassland = 0
all.models.temp$SE_SURROUNDINGforest = 0
all.models.temp$SE_SURROUNDINGmosaic = 0
all.models.temp$SE_SURROUNDINGplantation = 0


#normalize data
#dat.abundance.clear$bio1 = dat.abundance.clear$bio1 / max(dat.abundance.clear$bio1)




i=1
for(i in 1:length(formulas)){
  rma.fit <- rma.mv(Z_SCORE, VARIANCE_ESTIMATE,
                    mods= formulas[[i]],
                    random= ~ 1|STUDY_ID,
                    data=dat.abundance.clear,method="ML")
  all.models.temp = rbind(all.models.temp,data.frame("NR"=i,
                                                     "N"=rma.fit$k,
                                                     "FORMULA"=paste("~",as.character(formulas[[i]][2]),sep=""),
                                                     "AIC"=fitstats(rma.fit)[3],
                                                     "AICC"=fitstats(rma.fit)[5],
                                                     "BIC"=fitstats(rma.fit)[4],
                                                     "bio1"=  if(length(which(names(coef(rma.fit)) %in% "bio1")) >=1){coef(rma.fit)[which(names(coef(rma.fit)) %in% "bio1")]}else(0),
                                                     "NATURAL_OR_PLANTATIONSplantation" = if(length(which(names(coef(rma.fit)) %in% "NATURAL_OR_PLANTATIONSplantation")) >=1){coef(rma.fit)[which(names(coef(rma.fit)) %in% "NATURAL_OR_PLANTATIONSplantation")]}else(0),
                                                     "NATURAL_OR_PLANTATIONSforest" = if(length(which(names(coef(rma.fit)) %in% "NATURAL_OR_PLANTATIONSforest")) >=1){coef(rma.fit)[which(names(coef(rma.fit)) %in% "NATURAL_OR_PLANTATIONSforest")]}else(0),
                                                     "SURROUNDINGgrassland" = if(length(which(names(coef(rma.fit)) %in% "SURROUNDINGgrassland")) >=1){coef(rma.fit)[which(names(coef(rma.fit)) %in% "SURROUNDINGgrassland")]}else(0),
                                                     "SURROUNDINGforest" =  if(length(which(names(coef(rma.fit)) %in% "SURROUNDINGforest")) >=1){coef(rma.fit)[which(names(coef(rma.fit)) %in% "SURROUNDINGforest")]}else(0),
                                                     "SURROUNDINGmosaic" = if(length(which(names(coef(rma.fit)) %in% "SURROUNDINGmosaic")) >=1){coef(rma.fit)[which(names(coef(rma.fit)) %in% "SURROUNDINGmosaic")]}else(0),
                                                     "SURROUNDINGplantation" = if(length(which(names(coef(rma.fit)) %in% "SURROUNDINGplantation")) >=1){coef(rma.fit)[which(names(coef(rma.fit)) %in% "SURROUNDINGplantation")]}else(0),
                                                     
                                                     
                                                     "SE_bio1"=  if(length(which(names(coef(rma.fit)) %in% "bio1")) >=1){rma.fit$se[which(names(coef(rma.fit)) %in% "bio1")]}else(0),
                                                     "SE_NATURAL_OR_PLANTATIONSplantation" = if(length(which(names(coef(rma.fit)) %in% "NATURAL_OR_PLANTATIONSplantation")) >=1){rma.fit$se[which(names(coef(rma.fit)) %in% "NATURAL_OR_PLANTATIONSplantation")]}else(0),
                                                     "SE_NATURAL_OR_PLANTATIONSforest" = if(length(which(names(coef(rma.fit)) %in% "NATURAL_OR_PLANTATIONSforest")) >=1){rma.fit$se[which(names(coef(rma.fit)) %in% "NATURAL_OR_PLANTATIONSforest")]}else(0),
                                                     "SE_SURROUNDINGgrassland" = if(length(which(names(coef(rma.fit)) %in% "SURROUNDINGgrassland")) >=1){rma.fit$se[which(names(coef(rma.fit)) %in% "SURROUNDINGgrassland")]}else(0),
                                                     "SE_SURROUNDINGforest" =  if(length(which(names(coef(rma.fit)) %in% "SURROUNDINGforest")) >=1){rma.fit$se[which(names(coef(rma.fit)) %in% "SURROUNDINGforest")]}else(0),
                                                     "SE_SURROUNDINGmosaic" = if(length(which(names(coef(rma.fit)) %in% "SURROUNDINGmosaic")) >=1){rma.fit$se[which(names(coef(rma.fit)) %in% "SURROUNDINGmosaic")]}else(0),
                                                     "SE_SURROUNDINGplantation" = if(length(which(names(coef(rma.fit)) %in% "SURROUNDINGplantation")) >=1){rma.fit$se[which(names(coef(rma.fit)) %in% "SURROUNDINGplantation")]}else(0)))
}

#remove first line
all.models.abundance = all.models.temp[-1,]
#remove all factors in the model to which the remaining fixed factors refer
all.models.abundance= all.models.abundance[,-c(9,10,11,13,16,17,18,20)]
all.models.abundance = all.models.abundance[order(all.models.abundance$AICC),]
all.models.abundance


#calculate the relative evidence weight for each model
all.models.abundance$LIKELIHOOD = exp(1)^(-(all.models.abundance$AICC - all.models.abundance$AICC[1]) / 2)
all.models.abundance$AKAIKE_WEIGHT = all.models.abundance$LIKELIHOOD / sum(all.models.abundance$LIKELIHOOD)


#calculate weighted estimates and variances for all factors
all.models.abundance[,c(10:12)] = (all.models.abundance[,c(10:12)] * sqrt(all.models.abundance$N))^2

data.frame(names(all.models.abundance))

importances.abundance = data.frame("bio1" = sum(all.models.abundance$AKAIKE_WEIGHT[which(!(all.models.abundance$bio1 == 0))]),
                                "NATURAL_OR_PLANTATIONS" = sum(all.models.abundance$AKAIKE_WEIGHT[which(!(all.models.abundance$NATURAL_OR_PLANTATIONSplantation == 0))]),
                                "SURROUNDING" = sum(all.models.abundance$AKAIKE_WEIGHT[which(!(all.models.abundance$SURROUNDINGmosaic == 0))]))
                                

estimates.abundance = data.frame("bio1" = sum(all.models.abundance$AKAIKE_WEIGHT * all.models.abundance$bio1),
                              "NATURAL_OR_PLANTATIONSplantation" = sum(all.models.abundance$AKAIKE_WEIGHT * all.models.abundance$NATURAL_OR_PLANTATIONSplantation),
                              "SURROUNDINGmosaic" = sum(all.models.abundance$AKAIKE_WEIGHT * all.models.abundance$SURROUNDINGmosaic))
                              

variances.abundance = data.frame("bio1" = sum(all.models.abundance$AKAIKE_WEIGHT*(all.models.abundance$SE_bio1+((all.models.abundance$bio1 - mean(all.models.abundance$bio1))^2))),
                              "NATURAL_OR_PLANTATIONSplantation" = sum(all.models.abundance$AKAIKE_WEIGHT*(all.models.abundance$SE_NATURAL_OR_PLANTATIONSplantation+((all.models.abundance$NATURAL_OR_PLANTATIONSplantation - mean(all.models.abundance$NATURAL_OR_PLANTATIONSplantation))^2))),
                              "SURROUNDINGmosaic" = sum(all.models.abundance$AKAIKE_WEIGHT*(all.models.abundance$SE_SURROUNDINGmosaic+((all.models.abundance$SURROUNDINGmosaic - mean(all.models.abundance$SURROUNDINGmosaic))^2))))

round(importances.abundance,3)
fisher.to.r(estimates.abundance)
fisher.to.r(variances.abundance)

#model fit for the best model
best.model.abundance = rma.mv(Z_SCORE, VARIANCE_ESTIMATE,
                           mods= ~ bio1,
                           random = ~ 1|STUDY_ID,
                           data=dat.abundance.clear,method="ML")

best.model.abundance.no.mods = rma.mv(Z_SCORE, VARIANCE_ESTIMATE,
                                   #mods= ~ factor(ECOREGION),
                                   random = ~ 1|STUDY_ID,
                                   data=dat.abundance.clear,method="ML")

best.model.abundance
best.model.abundance.no.mods
#pseudo R2
((best.model.abundance.no.mods$sigma2 - best.model.abundance$sigma2)/best.model.abundance.no.mods$sigma2)^2


theme_set(theme_bw(base_size=18))

svg(filename="abundance_bio_correlation.svg",height=5,width=6)
p = ggplot() + theme(legend.position="none",axis.ticks = element_blank())
p = p + geom_point(aes(x=dat.abundance.clear$bio1,y=fisher.to.r(dat.abundance.clear$Z_SCORE),size=(1/dat.abundance.clear$VARIANCE_ESTIMATE/max(dat.abundance.clear$VARIANCE_ESTIMATE))),colour="#666666")
p= p + scale_size(range=c(3,8)) + ylim(c(-1,1)) + xlim(c(0,30))
p = p + geom_line(aes(x=dat.abundance.clear$bio1,y=fitted(best.model.abundance)))
p = p + xlab("Mean annual temperature") + ylab("Correlation coefficient") + ggtitle("b) Herbivore abundance")
p
graphics.off()







plot(Z_SCORE ~ bio1, data=dat.abundance.clear)
abline(lm(fitted(best.model.abundance) ~ dat.abundance.clear$bio1))
plot(residuals(best.model.abundance) ~ dat.abundance.clear$bio1)
shapiro.test(residuals(best.model.abundance))





#model fit for the best model
best.model.abundance = rma(Z_SCORE, VARIANCE_ESTIMATE,
                        mods= ~ bio1,
                        data=dat.abundance.clear,method="ML")

best.model.abundance

#get a pseudo R2
reduced.abundance = rma(Z_SCORE, VARIANCE_ESTIMATE,
                     data=dat.abundance.clear,method="ML") 
reduced.abundance
#-> cannot estimate R2 because Heterogeneity (tau2 is even larger in the case with moderators)






#meta-regression with random effect model, method is ML because REML is not comparable between different models
#no HERBIVORW_GUILD and no FORMER_LAND_USE because of too small sample size
rma.fit.abundance1 <- rma.mv(Z_SCORE, VARIANCE_ESTIMATE,
                          mods= ~ bio1 + NATURAL_OR_PLANTATIONS + AGE_FOREST + SURROUNDING,
                          random= ~ 1 | FOREST_NAME,
                          data=dat.abundance.clear,method="ML")
summary(rma.fit.abundance1)

#threre is no call stored in the rma.fit object, so I have to reduce complexity without the update() function
#exclude NATURAL_OR_PLANTATIONS
rma.fit.abundance2 <- rma.mv(Z_SCORE, VARIANCE_ESTIMATE,
                          mods= ~ bio1 + AGE_FOREST + SURROUNDING ,
                          random= ~ 1 | FOREST_NAME,
                          data=dat.abundance.clear,method="ML")
summary(rma.fit.abundance2)

#exclude SURROUNDING
rma.fit.abundance3 <- rma.mv(Z_SCORE, VARIANCE_ESTIMATE,
                          mods= ~ bio1,
                          random= ~ 1 | FOREST_NAME,
                          data=dat.abundance.clear,method="ML")
summary(rma.fit.abundance3)



rma.fit.abundance3[[1]] = fisher.to.r(rma.fit.abundance3[[1]])
rma.fit.abundance3[[2]] = fisher.to.r(rma.fit.abundance3[[1]])
rma.fit.abundance3[[5]] = fisher.to.r(rma.fit.abundance3[[5]])
rma.fit.abundance3[[6]] = fisher.to.r(rma.fit.abundance3[[6]])

summary(rma.fit.abundance3)

#try ploting with ggplot2
theme_set(theme_bw(base_size=18))

plot.data.abundance = cbind(dat.abundance,ES=rma.abundance.all$yi,SE=sqrt(rma.abundance.all$vi)/sqrt(dat.abundance$PLOT_NUMBER))

# data.frame(Plantation = dat.abundance$NATURAL_OR_PLANTATIONS, Ecoregion= dat.abundance$ECOREGION,ES=rma.abundance$yi,SE=sqrt(rma.abundance$vi),Type="Study",Study=dat.abundance$UNIQUE_ID)
plot.data.abundance = plot.data.abundance[order(plot.data.abundance$ES),]
plot.data.abundance = plot.data.abundance[,-c(34:53)]

plot.data.abundance = rbind(plot.data.abundance,data.frame(CASE_ID="grand mean",STUDY_ID =NA,FOREST_NR=NA,FOREST_NAME="summary",AUTHOR="summary",
                                                     YEAR=NA,PEER_OR_GREY="summary",DAMAGE_OR_HERBIVORES="summary",HERBIVORE_SPECIFICITY=NA,
                                                     UNITY="summary",HERBIVORE_GUILD="summary",HERBIVORE_SPECIES="summary",TREE_DIV_METRIC="summary",
                                                     SAMPLING_ON="summary",BETA_ESTIMATE=NA,R_SQUARE=NA,
                                                     R_SQUARE_ADJUSTED=NA,P_VALUE=NA,PLOT_NUMBER=NA,
                                                     MIN_RICH=NA,MAX_RICH=NA,MAX_DAMAGE_FOR_STANDARD=NA,
                                                     NATURAL_OR_PLANTATIONS=NA,SURROUNDING="summary",AGE_FOREST=NA,FORMER_LAND_USE="summary",
                                                     GPS_DEZ_X=NA,GPS_DEZ_Y=NA,ALTITUDE=NA,ECOREGION="summary",
                                                     FOCAL_SPECIES="summary",JUST_FOCAL="summary",PLOT_SIZE=NA,R_VALUE=NA,
                                                     Z_SCORE=NA,VARIANCE_ESTIMATE=NA,ES=rma.abundance.all$b,SE=rma.abundance.all$se))                


plot.data.abundance$CASE_ID = factor(plot.data.abundance$CASE_ID,levels=rev(c(plot.data.abundance$CASE_ID)))
plot.data.abundance$VEC_COL1 = c(rep("Study",27),NA)
plot.data.abundance$VEC_COL2 = c(rep(NA,27),"Grand Mean")

svg(filename="herbivore_abundance.svg",height=6,width=6)
plot.empty = ggplot(data=plot.data.abundance) 
plot.empty = plot.empty + geom_blank(aes(x=factor(CASE_ID),y=ES)) + coord_flip() + geom_hline(aes(x=0), lty=2,size=1,colour="grey") 
plot.full = plot.empty + geom_point(data=plot.data.abundance,aes(x=factor(CASE_ID),y=ES,size=5-(plot.data.abundance$SE*5), colour=plot.data.abundance$NATURAL_OR_PLANTATIONS)) 
plot.full = plot.full + geom_pointrange(data=plot.data.abundance,aes(x=factor(CASE_ID),y=ES,ymin=ES-1.96*SE,ymax=ES+1.96*SE, colour=plot.data.abundance$VEC_COL2),size=1.2,pch=3) + scale_colour_manual(values=c("grey","black","black"))
plot.full.abundance= plot.full +  ylim(c(-3,3)) + theme(legend.position="none") +labs(y="Fisher's z-transformed correlation coefficient",x="Study", title="b) Herbivore Abundance")
plot.full.abundance
graphics.off()






#diversity analysis ------------------------------------
#####################################################
dat.richness.clear <- subset(dat.filtered, DAMAGE_OR_HERBIVORES %in% "richness" & !(VARIANCE_ESTIMATE == Inf)  & VARIANCE_ESTIMATE >= 0 & !(Z_SCORE %in% NaN)& !(Z_SCORE %in% Inf))


#random effect model without moderators
rma.richness.all <- rma.uni(Z_SCORE, VARIANCE_ESTIMATE,
                           method="REML", data=dat.richness.clear)

forest(rma.richness.all)
summary(rma.richness.all)

#try ploting with ggplot2
plot.data.richness = rbind(data.frame("CASE_ID" = as.factor("grand mean"),
                                    "ES" =rma.richness.all$b,
                                    "SE"=rma.richness.all$se),
                         data.frame("CASE_ID" = as.factor("spaceholder"), #space holder
                                    "ES" = 0,
                                    "SE"= 0),
                         data.frame("CASE_ID" =factor(dat.richness.clear$CASE_ID,levels=dat.richness.clear$CASE_ID),
                                    "ES" =rma.richness.all$yi,
                                    "SE" =sqrt(rma.richness.all$vi)/sqrt(dat.richness.clear$PLOT_NUMBER)))

plot.data.richness$TYPE1 = c("black",rep(NA,nrow(plot.data.richness) -1))
plot.data.richness$TYPE2 = c(NA,NA,rep("grey",nrow(plot.data.richness) -2))



#plot.data.richness$CASE_ID = factor(plot.data.richness$CASE_ID,levels=(c(plot.data.richness$CASE_ID)))
#plot.data.richness$TYPE = factor(plot.data.richness$TYPE,levels=rev(c(plot.data.richness$TYPE)))
levels(plot.data.richness$CASE_ID)

theme_set(theme_bw(base_size=18))

svg(filename="herbivore_richness.svg",height=3,width=6)
plot.empty = ggplot(data=(plot.data.richness)) 
plot.empty = plot.empty + geom_blank(aes(x=factor(CASE_ID),y=ES)) + coord_flip() + geom_hline(aes(x=0), lty=2,size=0.5,colour="black") 
plot.empty = plot.empty + geom_point(aes(x=factor(CASE_ID),y=ES,size=5-(plot.data.richness$SE*5),colour=(TYPE2)))
plot.empty = plot.empty + scale_colour_manual(values = c("grey"))
plot.full = plot.empty + geom_pointrange(aes(x=1,y=rma.richness.all$b[1,1],ymin=rma.richness.all$b[1,1] - 1.96*rma.richness.all$se,ymax= rma.richness.all$b + 1.96*rma.richness.all$se),
                                         shape=3, size=1.2,linetype=1)
plot.full +  ylim(c(-1.6,1.6)) + theme(legend.position="none",axis.ticks=element_blank(),axis.text.y = element_text(colour="grey")) +
  labs(y="Fisher's z",x="Study", title="c) Herbivore richness")
graphics.off()







#publication bias
#################
#funnel plot
funnel(rma.richness.all)
#trim and fill 
trimfill(rma.richness.all,verbose=TRUE )
#funnel plot
funnel(trim.richness)
#fail safe number
fsn(Z_SCORE, VARIANCE_ESTIMATE,data=dat.richness.clear,
    type="Rosenthal",alpha=0.05,)
















theme_set(theme_bw(base_size=18))

plot.data.richness = cbind(dat.richness,ES=rma.richness.all$yi,SE=sqrt(rma.richness.all$vi)/sqrt(dat.richness$PLOT_NUMBER))

# data.frame(Plantation = dat.richness$NATURAL_OR_PLANTATIONS, Ecoregion= dat.richness$ECOREGION,ES=rma.richness$yi,SE=sqrt(rma.richness$vi),Type="Study",Study=dat.richness$UNIQUE_ID)
plot.data.richness = plot.data.richness[order(plot.data.richness$ES),]
plot.data.richness = plot.data.richness[,-c(34:53)]

plot.data.richness = rbind(plot.data.richness,data.frame(CASE_ID="grand mean",STUDY_ID =NA,FOREST_NR=NA,FOREST_NAME="summary",AUTHOR="summary",
                                                           YEAR=NA,PEER_OR_GREY="summary",DAMAGE_OR_HERBIVORES="summary",
                                                           UNITY="summary",HERBIVORE_GUILD="summary",HERBIVORE_SPECIES="summary", HERBIVORE_SPECIFICITY =NA, TREE_DIV_METRIC="summary",
                                                           SAMPLING_ON="summary",BETA_ESTIMATE=NA,R_SQUARE=NA,
                                                           R_SQUARE_ADJUSTED=NA,P_VALUE=NA,PLOT_NUMBER=NA,
                                                           MIN_RICH=NA,MAX_RICH=NA,MAX_DAMAGE_FOR_STANDARD=NA,
                                                           NATURAL_OR_PLANTATIONS=NA,SURROUNDING="summary",AGE_FOREST=NA,FORMER_LAND_USE="summary",
                                                           GPS_DEZ_X=NA,GPS_DEZ_Y=NA,ALTITUDE=NA,ECOREGION="summary",
                                                           FOCAL_SPECIES="summary",JUST_FOCAL="summary",PLOT_SIZE=NA,R_VALUE=NA,
                                                           Z_SCORE=NA,VARIANCE_ESTIMATE=NA,UNIQUE_ID = "summary",ES=rma.richness.all$b,SE=rma.richness.all$se))                


plot.data.richness$CASE_ID = factor(plot.data.richness$CASE_ID,levels=rev(c(plot.data.richness$CASE_ID)))
plot.data.richness$VEC_COL1 = c(rep("Study",6),NA)
plot.data.richness$VEC_COL2 = c(rep(NA,6),"Grand Mean")

svg(filename="herbivore_richness.svg",height=3,width=6)
plot.empty = ggplot(data=plot.data.richness) 
plot.empty = plot.empty + geom_blank(aes(x=factor(CASE_ID),y=ES)) + coord_flip() + geom_hline(aes(x=0), lty=2,size=1,colour="grey") 
plot.full = plot.empty + geom_point(data=plot.data.richness,aes(x=factor(CASE_ID),y=ES,size=5-(plot.data.richness$SE*5), colour=plot.data.richness$NATURAL_OR_PLANTATIONS)) 

plot.full = plot.full + geom_pointrange(data=plot.data.richness,aes(x=factor(CASE_ID),y=ES,ymin=ES-1.96*SE,ymax=ES+1.96*SE, colour=plot.data.richness$VEC_COL2),size=1.2,pch=3) + scale_colour_manual(values=c("grey","black","black"))
plot.full.diversity = plot.full +  ylim(c(-2,2)) + theme(legend.position="none") +labs(y="Fisher's z-transformed correlation coefficient",x="Study", title="c) Herbivore Species Richness")
plot.full.diversity
graphics.off()











#------------------------draw bar-charts to represent the weightened importance of moderating factors
###################################################################################################

importances.damage
importances.abundance

estimates.damage
estimates.abundance

variances.damage
variances.abundance

importances = data.frame("Driving_factors" =rep(c("MAT","Forest vs plantation","Age of tree stand","Herbivore taxon","Surrounding landscape","Interaction \n MAT-Herbivore taxon"),2),
                         "Importance" = c(0-importances.damage$bio1, 0-importances.damage$NATURAL_OR_PLANTATIONS,0-importances.damage$AGE_FOREST,0-importances.damage$HERBIVORE_GUILD,0-importances.damage$SURROUNDING,0-importances.damage$bio1HERBIVORES,importances.abundance$bio1,importances.abundance$NATURAL_OR_PLANTATIONS,importances.abundance$AGE_FOREST,NA,importances.abundance$SURROUNDING,NA),
                         "Herbivory_aspect" = c(rep("Herbivore damage",6),rep("Herbivore abundance",6)))

importances.a = importances[c(1:6),]
importances.b = importances[c(6:12),]

importances.a$Driving_factors = factor(importances.a$Driving_factors,rev(as.character(importances.a$Driving_factors)))
importances.b$Driving_factors = factor(importances.b$Driving_factors,rev(as.character(importances.b$Driving_factors)))


svg(filename="Weighted importances.svg",height=4,width=10)
p = ggplot() + theme_bw(base_size=25)
p = p + geom_bar(data=importances.b,aes(y=Importance,x=factor(Driving_factors)),stat="identity",fill="dark grey",colour="black")
p = p + geom_bar(data= importances.a,aes(y=Importance,x=order(Driving_factors)),stat="identity",fill="white",colour="black")
p = p + annotate("text", x=1, y=-0.8, label ="Herbivore \n damage",size=7.0) + annotate("text", x=1, y=0.8, label ="Herbivore \n abundance",size=7.0)
p = p + annotate("text", x= 3, y=0.2,label="NA",size=7.0) + annotate("text", x= 1, y=0.2,label="NA",size=7.0)
p + scale_y_continuous(limits=c(-1,1)) + coord_flip() + xlab("")
graphics.off()


#calculate the models with the most important factors
#####################################################

dat.damage$bio1 = dat.damage$bio1 / max(dat.damage$bio1)
rma.damage.mods = rma.uni(Z_SCORE, VARIANCE_ESTIMATE,method="REML",
                          mods = ~ bio1 + factor(HERBIVORE_GUILD),
                          data=dat.damage)

rma.damage.mods


dat.abundance$bio1 = dat.abundance$bio1 / max(dat.abundance$bio1)
rma.abundance.mods = rma.uni(Z_SCORE, VARIANCE_ESTIMATE,method="REML",
                          mods = ~ bio1,
                          data=dat.abundance)
rma.abundance.mods





#multiplot for the sampling of damage and abundance
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


svg(filename="damage_and abundance_sampled.svg",height=10,width=10)
multiplot(p.sample.damage,p.sample.abundance)
graphics.off()


svg(filename="all_forest_plots.svg",height=10,width=10)
multiplot(plot.full.damage, plot.full.abundance,plot.full.diversity)
graphics.off()











#damage
p = ggplot(data=dat.damage)
p = ggplot(data=dat.abundance)


p + geom_point(aes(x=bio1,y=Z_SCORE,colour=factor(HERBIVORE_GUILD)),size=5) + geom_smooth(aes(x=bio1, y=Z_SCORE))

p + geom_point(aes(x=bio1,y=Z_SCORE),subset=.(bio1 >= 200 )) + geom_smooth(aes(x=bio1, y=Z_SCORE),subset=.(bio1 >= 200 ))
p + geom_point(aes(x=bio1,y=Z_SCORE,colour=factor(HERBIVORE_GUILD)),subset=.(NATURAL_OR_PLANTATIONS %in% "plantation"),size=5) + geom_smooth(aes(x=bio1, y=Z_SCORE),subset=.(NATURAL_OR_PLANTATIONS %in% "plantation"))
p + geom_point(aes(x=bio1,y=Z_SCORE),subset=.(HERBIVORE_GUILD %in% "insects")) + geom_smooth(aes(x=bio1, y=Z_SCORE),subset=.(HERBIVORE_GUILD %in% "insects"))

p + geom_point(aes(x=factor(HERBIVORE_GUILD),y=Z_SCORE))





abundance.overall <- rma(Z_SCORE, VARIANCE_ESTIMATE,method="REML",
                            data=dat.abundance,level=95,subset=((JUST_FOCAL %in% "pooled")))
rma.damage.insects <- rma(Z_SCORE, VARIANCE_ESTIMATE,method="REML",
                          mods = ~ factor(HERBIVORE_GUILD)*bio1,
                          data=dat.damage)



dat.test = subset(dat.filtered, !(DAMAGE_OR_HERBIVORES %in% "richness") & !(VARIANCE_ESTIMATE ==  Inf))


test1 = rma.uni(Z_SCORE, VARIANCE_ESTIMATE,method="REML",
                  data=dat.test)


test2 <- rma(Z_SCORE, VARIANCE_ESTIMATE,method="REML",
                    mods = ~ factor(DAMAGE_OR_HERBIVORES),
                    data=dat.test)

aov(test1,test2)

test.damage2 <- rma(Z_SCORE, VARIANCE_ESTIMATE,method="ML",
                     mods = ~ DAMAGE_OR_HERBIVORES,
                     data=dat.damage)





rma.abundance.mammals <- rma(Z_SCORE, VARIANCE_ESTIMATE,method="ML",
                            data=dat.abundance,level=85,subset=(HERBIVORE_GUILD %in% "mammals"))
rma.abundance.plantation <- rma(Z_SCORE, VARIANCE_ESTIMATE,method="ML",
                                data=dat.abundance,level=85,subset=(NATURAL_OR_PLANTATIONS =="plantation"))


rma.abundance.low <- rma(Z_SCORE, VARIANCE_ESTIMATE,method="REML",
                         data=dat.abundance,level=85,subset=(bio1 <= median(bio1)))
rma.abundance.high <- rma(Z_SCORE, VARIANCE_ESTIMATE,method="REML",
                          data=dat.abundance,level=85,subset=(bio1 >= median(bio1)))
rma.damage.low <- rma(Z_SCORE, VARIANCE_ESTIMATE,method="REML",
                      data=dat.damage,level=85,subset=(bio1 <= median((bio1)))
rma.damage.high <- rma(Z_SCORE, VARIANCE_ESTIMATE,method="REML",
                       data=dat.damage,level=85,subset=(bio1 >= median(bio1)))




#test area
#----------
