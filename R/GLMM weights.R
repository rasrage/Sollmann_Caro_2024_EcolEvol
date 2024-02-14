################################################################################
####### analyze coconut crab weights ###########################################

## Caro and Sollmann, Spatio-temporal metapopulation trends: the coconut crabs 
## of Zanzibar. Ecology and Evolution.

## Script implements analyses described under Methods - Weights
## and recreates figure 3 (weight panel)

rm(list=ls())

library(lme4)
library(MuMIn)
library(lmerTest)
library(DHARMa)
library(glmmTMB)
library(readxl)
library(writexl)
library(ggplot2)
library(cowplot)

#read in data
dat.all<-readRDS('Interim data products/Weight_data.rds')

##some data conversion incl sqrt for weights
##set Ngao East to 99 (to separate from Ngao West)
dat.all$Area[dat.all$Name == 'Ngao East']<-99

##Ngao East and Manta Cave should be removed as they only have 1 survey with
##captures
rem<-unique(dat.all$Area[dat.all$Name == 'Ngao East' |dat.all$Name == 'Manta Cave'])
dat<-dat.all[!dat.all$Area %in% rem,]

dat$Area<-as.factor(dat$Area)
dat$Interval.mo<-dat$Interval/365
dat$ID<-as.factor(1:nrow(dat))
dat$y<-sqrt(dat$Weight)


## basic data summaries
nrow(dat) #total sample size
nrow(dat[dat$Year<2023,]) #data points before 2023 when individual ID was determined
sum(dat$Sex ==1)/nrow(dat) #percent females
mean(dat$Weight[dat$Sex ==1]); sd(dat$Weight[dat$Sex ==1])
mean(dat$Weight[dat$Sex ==0]); sd(dat$Weight[dat$Sex ==0])



## base model: site random effect, sex
## sex=1 is females (smaller)
basemod<-glmmTMB(y~(1|Area)+Sex+Interval.mo, data=dat,
              dispformula =~Area)

summary(basemod)

## covariates from previous analyses: protection, hotel, exploitation
step1covs<-c("Protected","Hotel","Exploit")

mod1<-glmmTMB(y~(1|Area)+Sex+Interval.mo+Protected, data=dat,
              dispformula =~Area)
summary(mod1)

mod2<-glmmTMB(y~(1|Area)+Sex+Interval.mo+Hotel, data=dat,
              dispformula =~Area)
summary(mod2)

mod3<-glmmTMB(y~(1|Area)+Sex+Interval.mo+Exploit, data=dat,
              dispformula =~Area)
summary(mod3)

model.sel(basemod, mod1, mod2, mod3)

## both protection and exploitation are similar
## so add model with both

mod4<-glmmTMB(y~(1|Area)+Sex+Interval.mo + Exploit+Protected, data=dat,
              dispformula =~Area)
summary(mod4)
aic.out<-AIC(basemod, mod1, mod2, mod3, mod4)$AIC
betas<-sapply(lapply(list(basemod, mod1, mod2, mod3, mod4), fixef), function(x)x[[1]][4])
names(betas)<-NULL
se.beta<-sapply(list(basemod, mod1, mod2, mod3, mod4),function(x)
  summary(x)$coef$cond[, 2, drop = FALSE][4]
)
##note: manuscript table requires manual addition of second beta from mod4
##      as this code only extracts one beta per model
summary(mod4)$coef$cond[5,1:2]

out<-model.sel(basemod, mod1, mod2, mod3, mod4, rank=AIC)

tab3<-data.frame(Variable=c('NULL', step1covs, 'Protected+Exploited')[order(aic.out)],
                 dAIC=round(out$delta, dig=2),
                 beta=round(betas[order(aic.out)], dig=3),
                 SE=round(se.beta[order(aic.out)], dig=3))

write_xlsx(tab3, 'Manuscript/MainEffectsWeights_AIC.xlsx')


##set model with both important main effect as base trend model
modtrend<-mod4


### add 3 predictors as interactions with trend

modtrendb<-glmmTMB(y~(1|Area)+Sex+Exploit+Protected+Interval.mo+
                     Interval.mo:Exploit, data=dat,
                   dispformula =~Area)
summary(modtrendb)

modtrendc<-glmmTMB(y~(1|Area)+Sex+Exploit+Protected+Interval.mo+
                  Interval.mo:Protected, data=dat,
                  dispformula =~Area)
summary(modtrendc)

modtrendd<-glmmTMB(y~(1|Area)+Sex+Exploit+Protected+Interval.mo+
                  Interval.mo:Hotel, data=dat,
                  dispformula =~Area)

summary(modtrendd)


out2<-model.sel(modtrend,modtrendb, modtrendc,modtrendd, rank=AIC)
aic.out2<-AIC(modtrend,modtrendb, modtrendc,modtrendd)$AIC
betas2<-sapply(lapply(list(modtrend,modtrendb, modtrendc,modtrendd), fixef), function(x)x[[1]][6])
names(betas2)<-NULL
se.beta2<-sapply(list(modtrend,modtrendb, modtrendc,modtrendd),function(x)
  summary(x)$coef$cond[, 2, drop = FALSE][6]
)

tab2<-data.frame(Variable=c('NULL', step1covs)[order(aic.out2)],
                 dAIC=round(out2$delta, dig=2),
                 beta=round(betas2[order(aic.out2)], dig=3),
                 SE=round(se.beta2[order(aic.out2)], dig=3))

write_xlsx(tab2, 'Manuscript/InteractionEffectsWeights_AIC.xlsx')
summary(modtrend)$coef$cond[3,]

##no evidence for any interactions with trend, and trend is almost 0
## therefore, no separate analysis for data-rich sites

## check model fit for top model
simulationOutput <- simulateResiduals(fittedModel = mod4)
plot(simulationOutput)

jpeg('Manuscript/ResidualPlotWeightsFull.jpg', width=18, height = 10, units = 'cm',
     res=600)
plot(simulationOutput)
dev.off()
## looks good except for minor patterns in residuals



################################################################################
##### prevalence of large/small individuals ####################################

### get 25th, 75th percentiles from all weights

perc<-quantile(dat$Weight, p=c(0.25, 0.75))

dat$small<-dat$large<-0

dat$small[dat$Weight<=perc[1]]<-1
dat$large[dat$Weight>=perc[2]]<-1

#### small individuals
## base model: site random effect, sex, trend
basemod<-glmmTMB(small~(1|Area)+Sex+Interval.mo, data=dat,
                 family='binomial')
summary(basemod)

## from previous analyses: protection, hotel, exploitation
mod1<-glmmTMB(small~(1|Area)+Sex+Interval.mo+Protected, data=dat,
              family='binomial')
summary(mod1)

mod2<-glmmTMB(small~(1|Area)+Sex+Interval.mo+Hotel, data=dat,
              family='binomial')
summary(mod2)

mod3<-glmmTMB(small~(1|Area)+Sex+Interval.mo+Exploit, data=dat,
              family='binomial')
summary(mod3)

out<-model.sel(basemod, mod1, mod2, mod3, rank=AIC)

###protection is the only thing that matters
aic.out<-AIC(basemod, mod1, mod2, mod3)$AIC
betas<-sapply(lapply(list(basemod, mod1, mod2, mod3), fixef), function(x)x[[1]][4])
names(betas)<-NULL
se.beta<-sapply(list(basemod, mod1, mod2, mod3),function(x)
  summary(x)$coef$cond[, 2, drop = FALSE][4]
)

tab<-data.frame(Variable=c('NULL', step1covs)[order(aic.out)],
                 dAIC=round(out$delta, dig=2),
                 beta=round(betas[order(aic.out)], dig=3),
                 SE=round(se.beta[order(aic.out)], dig=3))

write_xlsx(tab, 'Manuscript/MainEffectsSMALL_AIC.xlsx')

### Interactions
mod1a<-glmmTMB(small~(1|Area)+Sex+Interval.mo+Protected +Interval.mo:Protected, 
               data=dat,
              family='binomial')
mod1b<-glmmTMB(small~(1|Area)+Sex+Interval.mo+Protected +Interval.mo:Hotel, 
               data=dat,
               family='binomial')
mod1c<-glmmTMB(small~(1|Area)+Sex+Interval.mo+Protected +Interval.mo:Exploit, 
               data=dat,
               family='binomial')

out2<-model.sel(mod1, mod1a, mod1b, mod1c, rank=AIC)
aic.out2<-AIC(mod1, mod1a, mod1b, mod1c)$AIC
betas2<-sapply(lapply(list(mod1, mod1a, mod1b, mod1c), fixef), function(x)x[[1]][5])
names(betas2)<-NULL
se.beta2<-sapply(list(mod1, mod1a, mod1b, mod1c),function(x)
  summary(x)$coef$cond[, 2, drop = FALSE][5]
)

tab2<-data.frame(Variable=c('NULL', step1covs)[order(aic.out2)],
                dAIC=round(out2$delta, dig=2),
                beta=round(betas2[order(aic.out2)], dig=3),
                SE=round(se.beta2[order(aic.out2)], dig=3))

write_xlsx(tab2, 'Manuscript/InteractionEffectsSMALL_AIC.xlsx')

##no evidence for any interaction effects, so look at top main effects model
summary(mod1)$coef$cond[3,] #trend, close to 0
summary(mod1)$coef$cond[2,] #effect of being female

## check model fit
simulationOutput <- simulateResiduals(fittedModel = mod1)
plot(simulationOutput)
##everything looks good

jpeg('Manuscript/ResidualPlotSmallFull.jpg', width=18, height = 10, units = 'cm',
     res=600)
plot(simulationOutput)
dev.off()

################################################################################
#### same procedure for large individuals

## base model: site random effect, sex, trend
basemod<-glmmTMB(large~(1|Area)+Sex+Interval.mo, data=dat,
                 family='binomial')
summary(basemod)

## from previous analyses: protection, hotel, exploitation

mod1<-glmmTMB(large~(1|Area)+Sex+Interval.mo+Protected, data=dat,
              family='binomial')
summary(mod1)

mod2<-glmmTMB(large~(1|Area)+Sex+Interval.mo+Hotel, data=dat,
              family='binomial')
summary(mod2)

# Data too sparse to include exploitation as a predictor
# mod3<-glmmTMB(large~(1|Area)+Sex+Interval.mo+Exploit, data=dat,
#               family='binomial')
# summary(mod3)

out<-model.sel(basemod, mod1, mod2,rank=AIC)

###both hotel and protection matter, but protection is similar to null model
aic.out<-AIC(basemod, mod1, mod2)$AIC
betas<-sapply(lapply(list(basemod, mod1, mod2), fixef), function(x)x[[1]][4])
names(betas)<-NULL
se.beta<-sapply(list(basemod, mod1, mod2),function(x)
  summary(x)$coef$cond[, 2, drop = FALSE][4]
)

tab<-data.frame(Variable=c('NULL', step1covs[1:2])[order(aic.out)],
                dAIC=round(out$delta, dig=2),
                beta=round(betas[order(aic.out)], dig=3),
                SE=round(se.beta[order(aic.out)], dig=3))

write_xlsx(tab, 'Manuscript/MainEffectsLARGE_AIC.xlsx')

##hotel and protected have similar support as main effects, but 'protected' is 
##not 2 dAIC better than null model, so move forward with hotel only

### Interactions
mod1a<-glmmTMB(large~(1|Area)+Sex+Interval.mo+Hotel +Interval.mo:Protected, 
               data=dat,
               family='binomial')
mod1b<-glmmTMB(large~(1|Area)+Sex+Interval.mo+Hotel +Interval.mo:Hotel, 
               data=dat,
               family='binomial')

out2<-model.sel(mod2, mod1a, mod1b,rank=AIC)
aic.out2<-AIC(mod2, mod1a, mod1b)$AIC
betas2<-sapply(lapply(list(mod2, mod1a, mod1b), fixef), function(x)x[[1]][5])
names(betas2)<-NULL
se.beta2<-sapply(list(mod2, mod1a, mod1b),function(x)
  summary(x)$coef$cond[, 2, drop = FALSE][5]
)

tab2<-data.frame(Variable=c('NULL', step1covs[1:2])[order(aic.out2)],
                 dAIC=round(out2$delta, dig=2),
                 beta=round(betas2[order(aic.out2)], dig=3),
                 SE=round(se.beta2[order(aic.out2)], dig=3))

write_xlsx(tab2, 'Manuscript/InteractionEffectsLARGE_AIC.xlsx')
###no evidence for an interaction effect, so look at trend of top main effects model
summary(mod2)$coef$cond[3,1] #trend, near 0
summary(mod2)$coef$cond[2,1] #effect fo being female

## check model fit
simulationOutput <- simulateResiduals(fittedModel = mod2)
plot(simulationOutput)
## some evidence of lack of fit, but likely only significant because of large
## sample size
jpeg('Manuscript/ResidualPlotLargeFull.jpg', width=18, height = 10, units = 'cm',
     res=600)
plot(simulationOutput)
dev.off()


################################################################################
########### recreate figure 3 (weight panel) ###################################

## excel sheet is created in CPUE analyis
trends<-read_xlsx('Trend estimates good data sites.xlsx')
nam.order<-c(trends$Site, unique(dat$Name)[c(8,11,13,12)])

##add limits for large and small
lg<-perc[2]
sm<-perc[1]
hline<-data.frame(y=c(lg, sm))

g1<-list()

for (i in 1:length(nam.order)){
  
  dat.sub<-dat[dat$Name == nam.order[i],c("Weight", "Interval.mo")]
  
  ##for leftmost panels, make left margin wider
  ##for rightmost panels make right margin wider
  marg.vec<-c(1,.3,.5,.3)
  if(i %in% c(1,5,9,13 ))
    marg.vec<-c(1,.3,.5,1)
  if(i %in% c(4,8,12 ))
    marg.vec<-c(1,1,.5,.3)
  
  #set max high enough that you can display limit for large individuals
  upper<-ifelse(max(dat.sub$Weight)<1, 1, max(dat.sub$Weight))
  
  #calculate median weight
  hline2<-data.frame(y=mean(dat.sub$Weight))
  
  g1[[i]]<-ggplot(dat.sub, aes(x = Interval.mo, y = Weight)) +
    ##data points
    geom_point(color='lightcoral', size=0.8) +
    labs(x = ifelse(i %in% c(10:13 ), 'Time', ''), 
         y = ifelse(i %in% c(1,5,9,13 ), "Weight", '')) +
    coord_cartesian(xlim=c(0, 7), ylim=c(0, upper)) +
    scale_x_continuous(breaks=seq(0,7, 1)) +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
    geom_hline(data=hline,aes(yintercept = y),
               linetype='dashed', color = 'darkgrey', 
               linewidth=0.5) +
    geom_hline(data=hline2,aes(yintercept = y),
               color = 'black', 
               linewidth=0.5) +
    theme_bw()+
    theme(plot.margin = unit(marg.vec, "points"),
          plot.title = element_text(size = 10,
                                    margin=margin(b = 2, unit = "pt")),
          axis.title = element_text(size = 10),
          axis.text.x = element_text(size = 8))+
    ggtitle(label=nam.order[i])
  
}

##combine

p<-plot_grid(plotlist=g1, align='hv',nrow=4, ncol=4, byrow=TRUE)

jpeg('Weight test figure.jpg', res=600, width=18, height=19, units = 'cm')
p
dev.off()



################################################################################
#### recruitment: presence of very small over time #############################

## recreates Appendix Figure A12

vsmall<-quantile(dat$Weight, p=0.05)
dat$Vsmall<-0
dat$Vsmall[dat$Weight<=vsmall]<-1

##extract calendar year/season 
xx<-sapply(strsplit(dat$unique, '_'), function(x)x[2])
xx<-gsub(' ', '', xx)
dat$Season<-as.factor(xx)

##subset to good data sites only
keep<-c(1,2,3,4,5,6, 7,9, 20)
dat.good<-dat[dat$Area %in% keep,]
nvs<-table(dat.good$Name,dat.good$Season, dat.good$Vsmall)

##calculate proportion very small in each season
pvs<-round(nvs[,,2]/apply(nvs, 1:2, sum), dig=3)

##chronologically, wet comes before dry, so change order of columns
nam<-dimnames(pvs)[[2]][c(1:3, 5, 4, 6, 8,7, 10, 9)]

pvs<-pvs[,nam]

##make abbreviated names as labels
lab<-c('CH', 'FU', 'KI', 'KO', 'JO', 'MA', 'MI', 'SW', 'VE')

p.df<-data.frame(Name=rep(rownames(pvs), ncol(pvs)),
                 Label=rep(lab, ncol(pvs)),
                 Season=rep(colnames(pvs), each=nrow(pvs)),
                 Proportion= c(pvs))
p.df$Season2=as.factor(substr(p.df$Season, 5, 8))
p.df$Year=as.factor(substr(p.df$Season, 1, 4))
p.df$x = rep(1:9, ncol(pvs))
#p.df$new.label<-rep(letters[1:9], ncol(pvs))

#remove non-sampled surveys
p.df<-p.df[!is.na(p.df$Proportion),]

##try as grid by season
p2<-ggplot(p.df, aes(x = x, y = Proportion, shape=Name)) +
  geom_point(size=3)+
  theme_bw()+
  facet_wrap(.~Season, ncol=4) +
  scale_shape_manual(values = 0:8)+ 
  theme(axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        axis.title.x = element_blank(),
        legend.position="bottom") 

jpeg('Recruitment figure.jpg', res=600, width=18, height=16, units = 'cm')
p2
dev.off()

