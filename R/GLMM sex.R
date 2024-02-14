################################################################################
#### analyze sex to see if there are trends over time ######################

## Caro and Sollmann, Spatio-temporal metapopulation trends: the coconut crabs 
## of Zanzibar. Ecology and Evolution.

## Script implements analyses described under Methods - Proportion of females
## and recreates figure 4 (p(female) panel)

rm(list=ls())

library(lme4)
library(merTools)
library(MuMIn)
library(lmerTest)
library(DHARMa)
library(glmmTMB)
library(readxl)
library(writexl)
library(msm)
library(ggplot2)
library(cowplot)

#read in data
dat.all<-readRDS('Interim data products/Weight_data.rds')

##some data conversion 
##set Ngao East to 99 to distinguish from Ngao West
dat.all$Area[dat.all$Name == 'Ngao East']<-99

##Ngao East and Manta Cave should be removed as they only have 1 survey with
##captures
rem<-unique(dat.all$Area[dat.all$Name == 'Ngao East' |dat.all$Name == 'Manta Cave'])
dat<-dat.all[!dat.all$Area %in% rem,]

dat$Area<-as.factor(dat$Area)
dat$Interval.mo<-dat$Interval/365
dat$ID<-as.factor(1:nrow(dat))
dat$y<-dat$Sex

## base model: site random effect, sex
## sex=1 is females (smaller)
basemod<-glmmTMB(y~(1|Area)+Interval.mo, data=dat,
              family='binomial')

summary(basemod)

## important covariates from analyses of weight, CPUE
step1covs<-c("Protected","Agriculture","Exploit")

mod1<-glmmTMB(y~(1|Area)+Interval.mo+Protected, data=dat,
              family='binomial')
summary(mod1)

mod2<-glmmTMB(y~(1|Area)+Interval.mo+Agriculture, data=dat,
              family='binomial')
summary(mod2)

mod3<-glmmTMB(y~(1|Area)+Interval.mo+Exploit, data=dat,
              family='binomial')
summary(mod3)

model.sel(basemod, mod1, mod2, mod3)

aic.out<-AIC(basemod, mod1, mod2, mod3)$AIC
betas<-sapply(lapply(list(basemod, mod1, mod2, mod3), fixef), function(x)x[[1]][3])
names(betas)<-NULL
se.beta<-sapply(list(basemod, mod1, mod2, mod3),function(x)
  summary(x)$coef$cond[, 2, drop = FALSE][3]
)


out<-model.sel(basemod, mod1, mod2, mod3,rank=AIC)

tab3<-data.frame(Variable=c('NULL', step1covs)[order(aic.out)],
                 dAIC=round(out$delta, dig=2),
                 beta=round(betas[order(aic.out)], dig=3),
                 SE=round(se.beta[order(aic.out)], dig=3))

write_xlsx(tab3, 'Manuscript/MainEffectsSex_AIC.xlsx')
## base model is in similar to top model (exploitation) and effects of top model 
## are uncertain and not sign.

## use base model as best main effects model
modtrend<-basemod

### add 3 predictors as interaction with trend

modtrendb<-glmmTMB(y~(1|Area)+Interval.mo+Interval.mo:Protected, data=dat,
                   family='binomial')
summary(modtrendb)

modtrendc<-glmmTMB(y~(1|Area)+Interval.mo+Interval.mo:Agriculture, data=dat,
                   family='binomial')
summary(modtrendc)


modtrendd<-glmmTMB(y~(1|Area)+Interval.mo+Interval.mo:Exploit, data=dat,
                   family='binomial')

summary(modtrendd)


out2<-model.sel(modtrend,modtrendb, modtrendc,modtrendd, rank=AIC)
aic.out2<-AIC(modtrend,modtrendb, modtrendc,modtrendd)$AIC
betas2<-sapply(lapply(list(modtrend,modtrendb, modtrendc,modtrendd), fixef), function(x)x[[1]][3])
names(betas2)<-NULL
se.beta2<-sapply(list(modtrend,modtrendb, modtrendc,modtrendd),function(x)
  summary(x)$coef$cond[, 2, drop = FALSE][3]
)

tab2<-data.frame(Variable=c('NULL', step1covs)[order(aic.out2)],
                 dAIC=round(out2$delta, dig=2),
                 beta=round(betas2[order(aic.out2)], dig=3),
                 SE=round(se.beta2[order(aic.out2)], dig=3))

write_xlsx(tab2, 'Manuscript/InteractionEffectsSex_AIC.xlsx')
## important interaction of exploitation with trend

## look at trends exploited and unexploited
summary(modtrendd)$coef$cond[2,1:2] #unexploited
sum(summary(modtrendd)$coef$cond[2:3,1]) #trend exploited
##SE for trend in exploited: the variance of the sum is the sum of the variances
sqrt(0.03229846^2 + 0.06945286^2 )


##note: only 2 exploited sites, so could be other factors as well

## check model fit
simulationOutput <- simulateResiduals(fittedModel = modtrendd)
plot(simulationOutput)
## somepatterns in residuals
jpeg('Manuscript/ResidualPlotSexFull.jpg', width=18, height = 10, units = 'cm',
     res=600)
plot(simulationOutput)
dev.off()



################################################################################
############# analysis of data-rich sites ######################################

## subset to good data sites
incl.area<-c(1,2,3,4,5,6, 7,9, 20)
dat.good<-dat[dat$Area %in% incl.area,]
incl.name<-unique(dat.good$Name)[pmatch(incl.area, unique(dat.good$Area))]

## read in trajectory categories
traj<-read_xlsx('Trend estimates good data sites.xlsx')

mod.f<-glmmTMB(y~1+Interval.mo*Area, data=dat.good,
               family='binomial')

##check for model fit - generally fits, some patterns in residuals
simulationOutput <- simulateResiduals(fittedModel = mod.f)
jpeg('Manuscript/ResPlotSexSubset.jpg', width=18, height = 10, units = 'cm',
     res=600)
plot(simulationOutput)
dev.off()


###### get trends in p(females) on logit scale; also baseline proportions

trend.mat<-prop.mat<-matrix(NA, length(incl.area), 4)
vcv<-vcov(mod.f)[[1]]
ses<-sqrt(diag(vcv))

#adjust indicator depending on model structure
i1<-which(names(fixef(mod.f)[[1]]) == 'Interval.mo')
i2<-which(names(fixef(mod.f)[[1]]) == 'Interval.mo:Area2')-2
b1<-which(names(fixef(mod.f)[[1]]) == '(Intercept)')
b2<-which(names(fixef(mod.f)[[1]]) == 'Area2')-2
  
trend.mat[1,1]<-fixef(mod.f)[[1]][i1]
trend.mat[1,2]<-ses[i1]
trend.mat[1,3]<-fixef(mod.f)[[1]][i1] - 1.96*ses[i1]
trend.mat[1,4]<-fixef(mod.f)[[1]][i1] + 1.96*ses[i1]

prop.mat[1,1]<-plogis(fixef(mod.f)[[1]][b1])
prop.mat[1,2]<-deltamethod(~exp(x1)/(exp(x1)+1), mean=fixef(mod.f)[[1]][b1], 
                           cov = vcv[b1,b1])
prop.mat[1,3]<-plogis(fixef(mod.f)[[1]][b1] - 1.96*ses[b1])
prop.mat[1,4]<-plogis(fixef(mod.f)[[1]][b1] + 1.96*ses[b1])

for (j in 2:9){
  trend.mat[j,1]<-fixef(mod.f)[[1]][i1] + fixef(mod.f)[[1]][(j+i2)]
  trend.mat[j,2]<-deltamethod(~x1 + x2, mean=fixef(mod.f)[[1]][c(i1, (j+i2))], 
                              cov = vcv[c(i1, (j+i2)),c(i1, (j+i2))])
  trend.mat[j,3]<-(fixef(mod.f)[[1]][i1] + fixef(mod.f)[[1]][(j+i2)]) - 1.96*trend.mat[j,2]
  trend.mat[j,4]<-(fixef(mod.f)[[1]][i1] + fixef(mod.f)[[1]][(j+i2)]) + 1.96*trend.mat[j,2]
  
  prop.mat[j,1]<-plogis(fixef(mod.f)[[1]][b1] + fixef(mod.f)[[1]][(j+b2)])
  prop.mat[j,2]<-deltamethod(~exp(x1 + x2)/(1+ exp(x1 + x2)), 
                             mean=fixef(mod.f)[[1]][c(b1, (j+b2))], 
                              cov = vcv[c(b1, (j+b2)),c(b1, (j+b2))])
  se.link<-deltamethod(~x1 + x2, mean=fixef(mod.f)[[1]][c(b1, (j+b2))], 
                             cov = vcv[c(b1, (j+b2)),c(b1, (j+b2))])
  prop.mat[j,3]<-plogis((fixef(mod.f)[[1]][b1] + fixef(mod.f)[[1]][(j+b2)]) - 1.96*se.link)
  prop.mat[j,4]<-plogis((fixef(mod.f)[[1]][b1] + fixef(mod.f)[[1]][(j+b2)]) + 1.96*se.link)
  
}

trend.mat<-round(trend.mat, dig=2)
prop.mat<-round(prop.mat, dig=2)

##write out trend estimates

trends<-data.frame(Site=incl.name, Code=incl.area, Trend=trend.mat[,1],
                   SE=trend.mat[,2], CI.lower=trend.mat[,3],
                   CI.upper=trend.mat[,4])
##say: above 1.1 - growing; below 0.9 - declining; in between - stable
trends$Category<-'stable'
trends$Category[trends$Trend<=-0.1]<-'declining'
trends$Category[trends$Trend>=0.1]<-'increasing'

trends[,c('Proportion', 'SE.Prop.', 'LCI.Prop', 'UCI.Prop')]<-prop.mat

##order: incr, stable, decl
trends$Category<-factor(trends$Category, levels=c('increasing', 'stable', 'declining'))
trends<-trends[order(trends$Category),]
write_xlsx(trends, 'Manuscript/Trend p(females) good data sites.xlsx')


################################################################################
### make panel of trends in p(female) (Figure 4)################################

## re-fit top models using glmer (to use merTools)
##fit top model for all sites 
top.all<-glmer(y~(1|Area)+Interval.mo+Interval.mo:Exploit, data=dat,
               family='binomial',
               control=glmerControl(optimizer="bobyqa",
                                    optCtrl=list(maxfun=1e6)))
##fit top model for data-rich sites
top.sub<-glm(y~1+Interval.mo*Area, data=dat.good,
             family='binomial')

##make predictions intervals and combine in a single data frame

##for data-poor sites
##create ranges of values for covariates in model
all.area<-as.numeric(as.character(sort(unique(dat$Area))))
excl.area<-all.area[!all.area %in% incl.area]

new.int<-seq(0,7,0.1)
new.ar<-as.factor(rep(excl.area, each=length(new.int)))
new.expl<-rep(dat$Exploit[pmatch(excl.area, dat$Area)], each=length(new.int)) 

newdat=data.frame(Area=new.ar,
                  Interval.mo=rep(new.int, length(excl.area)),
                  Exploit = new.expl)

pi<-predictInterval(top.all, newdata=newdat,
                    level=0.95, n.sims=10000,
                    stat = "median", type="linear.prediction",
                    include.resid.var = FALSE)

##for data-rich sites
new.int<-seq(0,7,0.1)
new.ar<-as.factor(rep(incl.area, each=length(new.int)))

newdat=data.frame(Interval.mo=rep(new.int, length(incl.area)),
                  Area=new.ar)

pi.rich.pre<-predict(top.sub, newdata=newdat, se.fit = TRUE, type='link')
pi.rich<-data.frame(fit=pi.rich.pre$fit, 
                    upr=pi.rich.pre$fit + 1.96*pi.rich.pre$se.fit, 
                    lwr=pi.rich.pre$fit - 1.96*pi.rich.pre$se.fit)

##combine and add info for plotting
pi.all<-rbind(pi, pi.rich)
pi.all$Area<-c(rep(excl.area, each=length(new.int)), 
               rep(incl.area, each=length(new.int)))
pi.all$Name<-dat$Name[pmatch(pi.all$Area, dat$Area, duplicates.ok = TRUE)]
pi.all$Time<-c(rep(new.int, length(excl.area)), 
               rep(new.int, length(incl.area)))

##use same order as for CPUE
##sort: data rich sites first, increasing-stable-decreasing, then data-poor sites
## excel sheet is created in CPUE analysis
trends<-read_xlsx('Trend estimates good data sites.xlsx')

ar.order<-pmatch(c(excl.area, incl.area), trends$Code)
ar.order[is.na(ar.order)]<-c(10,11,13,12)
pi.all$Order<-rep(ar.order, each = length(new.int))

##finally, sort by Order first, then Time
pi.all<-pi.all[with(pi.all, order(Order, Time)),]

##plogis!
pi.all[,c('fit','upr', 'lwr')]<-apply(pi.all[,c('fit','upr', 'lwr')], 2, plogis)
nam.order<-unique(pi.all$Name)

## plot separately per site, then combine

g1<-list()

for (i in 1:length(nam.order)){
  pi.sub<-pi.all[pi.all$Order == i, ]
  dat.sub<-dat[dat$Name == unique(pi.sub$Name),]
  
  ##summarize by survey
  
  nsv<-unique(dat.sub$Interval.mo)
  sdf<-data.frame(prop = rep(NA, length(nsv)), Interval.mo = nsv)
  for (gg in 1:length(nsv)){
    subsub<-dat.sub[dat.sub$Interval.mo == nsv[gg],]
    sdf$prop[gg]<-sum(subsub$Sex)/nrow(subsub)
  }
  
  ##for leftmost panels, make left margin wider
  ##for rightmost panels make right margin wider
  marg.vec<-c(1,.3,.5,.3)
  if(i %in% c(1,5,9,13 ))
    marg.vec<-c(1,.3,.5,1)
  if(i %in% c(4,8,12 ))
    marg.vec<-c(1,1,.5,.3)
  
  upper<-ifelse(max(c(pi.sub$upr, max(sdf$prop)))<0.1, 0.1, 
                max(c(pi.sub$upr, max(sdf$prop))))
  
  g1[[i]]<-ggplot(pi.sub, aes(x = Time, y = fit)) +
    ##data points
    geom_point(data = sdf, aes(x = Interval.mo, y = prop), col='lightcoral') +
    ##estimated trends
    geom_line(color = "black",  aes(group=Name)) + # alpha = .3,
    labs(x = ifelse(i %in% c(10:13 ), 'Time', ''), 
         y = ifelse(i %in% c(1,5,9,13 ), "p(female)", '')) +
    coord_cartesian(ylim=c(0, upper)) +
    geom_line(aes(x = Time, y = upr), color='darkgrey') + #lty=2
    geom_line(aes(x = Time, y = lwr), color='darkgrey') + 
    scale_x_continuous(breaks=seq(0,7, 1)) +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
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

jpeg('p(female) test figure.jpg', res=600, width=18, height=19, units = 'cm')
p
dev.off()
