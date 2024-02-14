################################################################################
#### Analysis of coconut crab CPUE data for all sites ##########################

## Caro and Sollmann, Spatio-temporal metapopulation trends: the coconut crabs 
## of Zanzibar. Ecology and Evolution.

## Script implements analyses described under Methods - Population trends
## and recreates figure 2 (CPUE panel)

rm(list=ls())

library(DHARMa)
library(MuMIn)
library(msm)
library(writexl)#
library(glmmTMB)
library(merTools)
library(readxl)
library(ggplot2)
library(cowplot)

##read in prepped data
dat<-readRDS('Interim data products/CPUE.trend.data.rds') 
##note: dat$y holds counts of crabs captured for each nightly visit

##do some covariate processing
##convert time since first survey from days to years
dat$Interval.mo<-dat$Interval/365
##convert site/subpopulation to factor
dat$Area.f<-as.factor(dat$New.Area)
#scale shehia level rain data
dat$RainDay.sc<-scale(as.numeric(dat$RainDayS))
dat$RainWeek.sc<-scale(as.numeric(dat$RainWeekS))
#binary rain/no rain on sampling day
dat$Rain.bin<-0
dat$Rain.bin[dat$RainDayS>0]<-1
##percent moon
dat$Moon.<-dat$Moon./100
#binary moon phase
dat$Moon.bin<-0
dat$Moon.bin[which(dat$Moon %in% c(12:16))]<-1

#fix missing data entry in one effort measure
dat$ActualSPmin[dat$ActualSPmin == '?']<-NA
dat$ActualSPmin<-as.numeric(dat$ActualSPmin)

##extract calendar year/season to check for overall environmental effects
xx<-sapply(strsplit(dat$unique, '_'), function(x)x[2])
xx<-gsub(' ', '', xx)
dat$Season<-as.factor(xx)

#add individual data point ID for random effect (lack of fit)
dat$ID<-as.factor(1:nrow(dat))

## some data summaries for manuscript: number of nights and people in field
##set RealTeam to 2 for Lisi's data
dat$RealTeam[is.na(dat$RealTeam)]<-2
nrow(dat)
mean(dat$RealTeam, na.rm=TRUE)
sd(dat$RealTeam, na.rm=TRUE)
range(dat$RealTeam, na.rm=TRUE)

################################################################################
###### Step 1 - test 3 alternative offset effort measures ######################

##subset to data for which all 3 measures are available
dat.sub<-dat[which(!is.na(dat$ActualSPmin)),]

##Infield: time in field, 
##SearchEst: time in field minus time measuring, 
##ActualSPmin: people searching times time searching

##basic model structure accounts for repeated surveys at sites, (1|Area.f),
##and repeated visits nested in surveys,  (1|Area.f:SurveyV3)

m1.1<-glmmTMB(y~(1|Area.f)+(1|Area.f:SurveyV3)+offset(log(Infield)), 
            family='poisson', data=dat.sub)
m1.2<-glmmTMB(y~(1|Area.f)+(1|Area.f:SurveyV3)+offset(log(SearchEst)), 
            family='poisson', data=dat.sub)
m1.3<-glmmTMB(y~(1|Area.f)+(1|Area.f:SurveyV3)+offset(log(ActualSPmin)), 
            family='poisson', data=dat.sub)


eff.tab<-AIC(m1.1, m1.2, m1.3)
eff.tab[order(eff.tab$AIC),]
out<-model.sel(m1.1, m1.2, m1.3, rank=AIC)

##make and write out model selection table for Appendix
tab<-data.frame(Effort=c('Infield','SearchEst','ActualSPmin')[order(eff.tab$AIC)],
                dAIC=round(out$delta, dig=2))
write_xlsx(tab, 'Manuscript/Effort_AIC.xlsx')

## move forward with Infield as offset  

################################################################################
#### Step 1a: do we need a season random effect? ###############################
m2<-glmmTMB(y~(1|Area.f)+(1|Area.f:SurveyV3)+offset(log(Infield)), 
            family='poisson', data=dat)

m2s<-glmmTMB(y~(1|Area.f)+(1|Area.f:SurveyV3)+(1|Season)+offset(log(Infield)), 
            family='poisson', data=dat)
AIC(m2, m2s)
##nope, model is worse

################################################################################
#### Step 2: check for visit level variables affecting capture rates ###########

## test daily/weekly and binary rain, moon percent and binary moon phase
## use subset of data for which all variables are available (%moon missing for
## two visits)
dat.sub<-dat[which(!is.na(dat$Moon.)),]

##'null' model
m2<-glmmTMB(y~(1|Area.f)+(1|Area.f:SurveyV3)+offset(log(Infield)), 
          family='poisson', data=dat.sub)

m2a<-glmmTMB(y~(1|Area.f)+(1|Area.f:SurveyV3)+offset(log(Infield)) + RainDay.sc, 
             family='poisson', data=dat.sub)
#summary(m2a)

m2b<-glmmTMB(y~(1|Area.f)+(1|Area.f:SurveyV3)+offset(log(Infield)) + RainWeek.sc, 
             family='poisson', data=dat.sub)
#summary(m2b)

m2c<-glmmTMB(y~(1|Area.f)+(1|Area.f:SurveyV3)+offset(log(Infield))+ Moon., 
           family='poisson', data=dat.sub)
#summary(m2c)

m2d<-glmmTMB(y~(1|Area.f)+(1|Area.f:SurveyV3)+offset(log(Infield))+ Moon.bin, 
           family='poisson', data=dat.sub)
#summary(m2d)

m2e<-glmmTMB(y~(1|Area.f)+(1|Area.f:SurveyV3)+offset(log(Infield))+ Rain.bin, 
           family='poisson', data=dat.sub)
#summary(m2e)

##compare via AIC
varr<-c('RainDay.sc', 'RainWeek.sc', 'Moon.', 'Moon.bin',
        'Rain.bin','Null')
#get beta estimates
betas<-sapply(lapply(list(m2a, m2b, m2c, m2d, m2e, m2), fixef), function(x)x$cond[2])
names(betas)<-NULL
se.beta<-sapply(list(m2a, m2b, m2c, m2d, m2e, m2),function(x)
  summary(x)$coef$cond[, 2, drop = FALSE][2]
)
var.tab<-AIC(m2a, m2b, m2c, m2d, m2e, m2)

out2<-model.sel(m2a, m2b, m2c, m2d, m2e, m2, rank = AIC)

tab2<-data.frame(Variable=varr[order(var.tab$AIC)],
             dAIC=round(out2$delta, dig=2),
             beta=round(betas[order(var.tab$AIC)], dig=3),
             SE=round(se.beta[order(var.tab$AIC)], dig=3))

write_xlsx(tab2, 'Manuscript/VisitVariables_AIC.xlsx')

##patterns hold: only binary moon phase is important

################################################################################
##### Step 3: add overall trend to create base model ###########################

## add overall trend, fit to full data
m2base<-glmmTMB(y~(1|Area.f)+(1|Area.f:SurveyV3)+offset(log(Infield))+ 
              Moon.bin + Interval.mo, 
            family='poisson', data=dat)
summary(m2base)

################################################################################
#### Step 4: test for main effects affecting baseline CPUE #####################

## use only top 3 covariates from previous CPUE models:
## Protected, Village and Agriculture
## Use only binary versions of these covariates

covs<-c("Agriculture", "Protected", "Village")

m2ha<-list()
for (jj in 1:length(covs)){
  
  m2ha[[jj]]<-glmmTMB(formula(paste0('y~(1|Area.f)+(1|Area.f:SurveyV3)+offset(log(Infield))+ 
              Moon.bin + Interval.mo +', covs[jj])), 
                    family='poisson', data=dat)
}
##add model without any main effects for comparison
m2ha[[4]]<-m2base
aic.out<-sapply(m2ha, AIC)
#get beta estimates
betas2<-sapply(lapply(m2ha, fixef), function(x)x$cond[4])
names(betas2)<-NULL
se.beta2<-sapply(m2ha,function(x)
  summary(x)$coef$cond[, 2, drop = FALSE][4]
)

out3<-model.sel(m2ha, rank=AIC)

tab3<-data.frame(Variable=c(covs, 'NULL')[order(aic.out)],
                 dAIC=round(out3$delta, dig=2),
                 beta=round(betas2[order(aic.out)], dig=3),
                 SE=round(se.beta2[order(aic.out)], dig=3))

write_xlsx(tab3, 'Manuscript/MainEffects_AIC.xlsx')
#summary(m2ha[[1]])

##top model is 'Agriculture', so keep that and then only add interaction, not main
## effect of all covariates

options(warn=0)
m2hi<-list()
for (jj in 1:length(covs)){
  
  m2hi[[jj]]<-glmmTMB(formula(paste0('y~(1|Area.f)+(1|Area.f:SurveyV3)+offset(log(Infield))+ 
              Moon.bin + Agriculture +Interval.mo+Interval.mo:', covs[jj])), 
                    family='poisson', data=dat)
}
##add top model from main effects for comparison
m2hi[[4]]<-m2ha[[1]]
aic.out2<-sapply(m2hi, AIC)
out4<-model.sel(m2hi, rank=AIC)
#get beta estimates
betas3<-sapply(lapply(m2hi, fixef), function(x)x$cond[5])
names(betas3)<-NULL
se.beta3<-sapply(m2hi,function(x)
  summary(x)$coef$cond[, 2, drop = FALSE][5]
)

tab4<-data.frame(Variable=c(covs, 'NULL')[order(aic.out2)],
                 dAIC=round(out4$delta, dig=2),
                 beta=round(betas3[order(aic.out2)], dig=3),
                 SE=round(se.beta3[order(aic.out2)], dig=3))

write_xlsx(tab4, 'Manuscript/InteractionEffects_AIC.xlsx')


##only interaction trend x protected is better than no-interaction model
## use delta method to calculate trend in protected and unprotected sites
fixef(m2hi[[2]])$cond[4]
vcv<-vcov(m2hi[[2]])[[1]]
sum(fixef(m2hi[[2]])$cond[4:5])
deltamethod(~x1+x2, mean = fixef(m2hi[[2]])$cond[4:5],
            cov = vcv[4:5,4:5])

################################################################################
#### Step 5: check best model for lack of fit  #################################

##check best model for lack of fit with dharma
topmodel<-glmmTMB(y~(1|Area.f)+(1|Area.f:SurveyV3)+offset(log(Infield))+ 
                  Moon.bin + Agriculture+Interval.mo+Interval.mo:Protected, 
                family='poisson', data=dat)
simulationOutput <- simulateResiduals(fittedModel = m2hi[[2]])
jpeg('Manuscript/ResPlotOriginal.jpg', width=18, height = 10, units = 'cm',
     res=600)
plot(simulationOutput)
dev.off()
##there is some evidence for lack of fit and patterns in residuals

##test for zero-inflation
jpeg('Manuscript/ZeroInflPlotOriginal.jpg', width=18, height = 10, units = 'cm',
     res=600)
testZeroInflation(simulationOutput)
dev.off()
## no evidence for zero inflation


## try running with negative binomial 
topmodel.nb<-glmmTMB(y~(1|Area.f)+(1|Area.f:SurveyV3)+offset(log(Infield))+ 
                  Moon.bin + Agriculture+Interval.mo+Interval.mo:Protected, 
                data=dat,family='nbinom2')
summary(topmodel.nb)
##very large dispersion parameter, so no evidence of overdispersion
simulationOutput <- simulateResiduals(fittedModel = topmodel.nb)
jpeg('Manuscript/ResPlotNegBin.jpg', width=18, height = 10, units = 'cm',
     res=600)
plot(simulationOutput)
dev.off()
## does not address lack of fit


## try adding visit level random effect
topmodel2<-glmmTMB(y~(1|Area.f)+(1|Area.f:SurveyV3)+(1|ID)+offset(log(Infield))+ 
                  Moon.bin + Agriculture+Interval.mo+Interval.mo:Protected, 
                family='poisson', data=dat)
simulationOutput <- simulateResiduals(fittedModel = topmodel2)
jpeg('Manuscript/ResPlotIndRaneff.jpg', width=18, height = 10, units = 'cm',
     res=600)
plot(simulationOutput)
dev.off()
## does not address lack of fit


################################################################################
#### Step 6: Analysis of 9 data-rich sites  ####################################

## Subset to the following sites:
## SD, Kigomasha, Misali, Fundo, Verani, Chumbe (good data)
## Jombe, matumbini, Kokota (fair data)

## code for aras to keep
incl.area<-c(1,2,3,4,5,6, 7,9, 20)

dat.good<-dat[dat$Area %in% incl.area,]
incl.name<-unique(dat.good$Name)[pmatch(incl.area, unique(dat.good$Area))]

## fit model with fixed site effect (main and interaction with trend)

## just out of curiosity, compare to a model using site as random effect
mr<-glmmTMB(y~(1|Area.f:SurveyV3)+offset(log(Infield))+ 
              Moon.bin + Interval.mo + (1+Interval.mo|Area.f), 
            family='poisson', data=dat.good)

mf<-glmmTMB(y~(1|Area.f:SurveyV3)+offset(log(Infield))+ 
              Moon.bin + Interval.mo*Area.f, 
            family='poisson', data=dat.good)

AIC(mr, mf)
## fixed effect model much better

##check for model fit
simulationOutput <- simulateResiduals(fittedModel = mf)
jpeg('Manuscript/ResPlotSubset.jpg', width=18, height = 10, units = 'cm',
     res=600)
plot(simulationOutput)
dev.off()
##this model fits well except for some pattern in residuals at high predicted 
##values



##calculate trend with SE and Wald 95%CI for all 6 sites
trend.mat<-matrix(NA, length(incl.area), 4)
vcv<-vcov(mf)[[1]]
ses<-sqrt(diag(vcv))

#adjust indicator depending on model structure
i1<-which(names(fixef(mf)[[1]]) == 'Interval.mo')
i2<-which(names(fixef(mf)[[1]]) == 'Interval.mo:Area.f2')-2

trend.mat[1,1]<-exp(fixef(mf)[[1]][i1])
trend.mat[1,2]<-deltamethod(~exp(x1), mean=fixef(mf)[[1]][i1], cov = vcv[i1,i1])
trend.mat[1,3]<-exp(fixef(mf)[[1]][i1] - 1.96*ses[i1])
trend.mat[1,4]<-exp(fixef(mf)[[1]][i1] + 1.96*ses[i1])

for (j in 2:9){
  trend.mat[j,1]<-exp(fixef(mf)[[1]][i1] + fixef(mf)[[1]][(j+i2)])
  trend.mat[j,2]<-deltamethod(~exp(x1 + x2), mean=fixef(mf)[[1]][c(i1, (j+i2))], 
                              cov = vcv[c(i1, (j+i2)),c(i1, (j+i2))])
  se.link<-deltamethod(~x1+x2, mean=fixef(mf)[[1]][c(i1, (j+i2))],
                       cov = vcv[c(i1, (j+i2)),c(i1, (j+i2))])
  trend.mat[j,3]<-exp((fixef(mf)[[1]][i1] + fixef(mf)[[1]][(j+i2)]) - 1.96*se.link)
  trend.mat[j,4]<-exp((fixef(mf)[[1]][i1] + fixef(mf)[[1]][(j+i2)]) + 1.96*se.link)
}

trend.mat<-round(trend.mat, dig=2)


##write out trend estimates

trends<-data.frame(Site=incl.name, Code=incl.area, Trend=trend.mat[,1],
                   SE=trend.mat[,2], CI.lower=trend.mat[,3],
                   CI.upper=trend.mat[,4])
##say: above 1.1 - growing; below 0.9 - declining; in between - stable
trends$Category<-'stable'
trends$Category[trends$Trend<=0.9]<-'declining'
trends$Category[trends$Trend>=1.1]<-'increasing'

##order: incr, stable, decl
trends$Category<-factor(trends$Category, levels=c('increasing', 'stable', 'declining'))
trends<-trends[order(trends$Category),]
write_xlsx(trends, 'Trend estimates good data sites.xlsx')



################################################################################
#### Step 7: Effect of education campaign  #####################################

#### for sites with 2020 campaign, check CPUE before and after campaign

camp<-c(1,2,3,4,5,20)
##I believe wet comes before dry season, chronologically
before<-c(paste0(2016:2020, 'dry'), paste0(2016:2020, 'wet'))
after<-c(paste0(2021:2023, 'dry'), paste0(2021:2023, 'wet'))

dat.camp<-dat.good[dat.good$Code %in% camp,]
dat.camp$after<-0
dat.camp$after[dat.camp$Season %in% after]<-1

m0<-glmmTMB(y~(1|Area.f:SurveyV3)+offset(log(Infield))+ 
              Moon.bin + Interval.mo*Area.f, 
            family='poisson', data=dat.camp)
mc<-glmmTMB(y~(1|Area.f:SurveyV3)+offset(log(Infield))+ 
              Moon.bin + Area.f+after, 
            family='poisson', data=dat.camp)
AIC(m0, mc)
summary(mc)
##model with campaign effect it is worse, campaign effect is near 0


################################################################################
#### Step 8: Make CPUE panel (Figure 2)  #######################################

##re-fit models with glmer in order to use predictInterval function from 
##merTools
##Estimates are identical to glmmtmb estimates

##fit top model for all sites 
top.all<-glmer(y~(1|Area.f)+(1|Area.f:SurveyV3)+offset(log(Infield))+ 
                 Moon.bin + Agriculture+Interval.mo+Interval.mo:Protected, 
               family='poisson', data=dat,
               control=glmerControl(optimizer="bobyqa",
                                    optCtrl=list(maxfun=1e6)))
##fit top model for data-rich sites
top.sub<-glmer(y~(1|Area.f:SurveyV3)+offset(log(Infield))+ 
                 Moon.bin + Interval.mo*Area.f, 
               family='poisson', data=dat.good,
               control=glmerControl(optimizer="bobyqa",
                                    optCtrl=list(maxfun=1e6)))

##make predictions intervals and combine in a single data frame

##for data-poor sites (from overall model)
all.area<-sort(unique(dat$New.Area))
excl.area<-all.area[!all.area %in% incl.area]
##create range of values for covariates in top model for full data
new.int<-seq(0,7,0.1)
new.ar<-as.factor(rep(excl.area, each=length(new.int)))
new.moon<-rep(0, length(new.ar))
new.inf<-rep(1, length(new.ar))
new.ag<-rep(dat$Agriculture[pmatch(excl.area, dat$New.Area)], each=length(new.int))
new.prot<-rep(dat$Protected[pmatch(excl.area, dat$New.Area)], each=length(new.int)) 

newdat=data.frame(Moon.bin=new.moon, 
                  Interval.mo=rep(new.int, length(excl.area)),
                  Infield=new.inf,
                  Area.f=new.ar,
                  Agriculture = new.ag,
                  Protected = new.prot,
                  SurveyV3=as.factor(rep(0,length(new.ar) )))

pi<-predictInterval(top.all, newdata=newdat,
                    level=0.95, n.sims=10000,
                    stat = "median", type="linear.prediction",
                    include.resid.var = FALSE)

##for data-rich sites (from data-rich model)
new.int<-seq(0,7,0.1)
new.ar<-as.factor(rep(incl.area, each=length(new.int)))
new.moon<-rep(0, length(new.ar))
new.inf<-rep(1, length(new.ar))

newdat=data.frame(Moon.bin=new.moon, 
                  Interval.mo=rep(new.int, length(incl.area)),
                  Infield=new.inf,
                  Area.f=new.ar,
                  SurveyV3=as.factor(rep(0,length(new.ar) )))

pi.rich<-predictInterval(top.sub, newdata=newdat,
                         level=0.95, n.sims=10000,
                         stat = "median", type="linear.prediction",
                         include.resid.var = FALSE)

##combine and add info for plotting
pi.all<-rbind(pi, pi.rich)
pi.all$Area<-c(rep(excl.area, each=length(new.int)), 
               rep(incl.area, each=length(new.int)))
pi.all$Name<-dat$Name[pmatch(pi.all$Area, dat$New.Area, duplicates.ok = TRUE)]
pi.all$Time<-c(rep(new.int, length(excl.area)), 
               rep(new.int, length(incl.area)))

##sort: data rich sites first, increasing-stable-decreasing, then data-poor sites
##xlsx read in here is created earlier in this script

trends<-read_xlsx('Trend estimates good data sites.xlsx')

ar.order<-pmatch(c(excl.area, incl.area), trends$Code)
##odd order to make plotting w names work
ar.order[is.na(ar.order)]<-c(10,11,15,13,14,12)
pi.all$Order<-rep(ar.order, each = length(new.int))

##finally, sort by Order first, then Time
pi.all<-pi.all[with(pi.all, order(Order, Time)),]

##exponentiate!
pi.all[,c('fit','upr', 'lwr')]<-exp(pi.all[,c('fit','upr', 'lwr')])

##plot per area, then combine (necessary due to different axes)
##meaning, only 1,5,9,13 have y axis label (all have axes, diff scales)
##only 13-15 have x axis, only 15 has x axis label 
##if upper limit is <0.1, use 0.1 as upper y axis limit
##campaign, removal if applicable
campaign<-read_xlsx('data-raw/Conservation Campaign dates.xlsx', sheet='PerSite')
campaign<-as.data.frame(campaign)
nam.order<-unique(pi.all$Name)
camp.order<-campaign[pmatch(nam.order, campaign$...1),]

g1<-list()

for (i in 1:15){
  pi.sub<-pi.all[pi.all$Order == i, ]
  dat.sub<-dat[dat$Name == unique(pi.sub$Name),]
  dat.sub$cpue<-dat.sub$y/dat.sub$Infield
  ##if campaign, calculate where in graph it should be
  if(!is.na(camp.order[i,'C.Year'])){
    mind<-min(dat.sub$midpoint)
    cdate<-as.Date(paste0(camp.order[i, 'C.Year'], '-', camp.order[i, 'C.Month'], '-', '01'))
    campx<-as.numeric(cdate-mind)/365
    vline<-data.frame(y=campx)} else {
      vline<-data.frame(y=NA)
    }
  
  if(campx<0) stop('HELP')
  
  ##for leftmost panels, make left margin wider
  ##for rightmost panels make right margin wider
  marg.vec<-c(1,.3,.5,.3)
  if(i %in% c(1,5,9,13 ))
    marg.vec<-c(1,.3,.5,1)
  if(i %in% c(4,8,12 ))
    marg.vec<-c(1,1,.5,.3)
  
  upper<-ifelse(max(c(pi.sub$upr, max(dat.sub$cpue)))<0.1, 0.1, 
                max(c(pi.sub$upr, max(dat.sub$cpue))))
  yr.dat<-data.frame(x=0.7, y=upper, 
                     Year=as.character(min(dat.sub$Year)),
                     Year2=substr(as.character(min(dat.sub$Year)), 3,4))
  g1[[i]]<-ggplot(pi.sub, aes(x = Time, y = fit)) +
    ##data points
    geom_point(data = dat.sub, aes(x = Interval.mo, y = cpue), col='lightcoral') +
    ##estimated trends
    geom_line(color = "black",  aes(group=Name)) + # alpha = .3,
    labs(x = ifelse(i %in% c(12:15), 'Time', ''), 
         y = ifelse(i %in% c(1,5,9,13 ), "CPUE", '')) +
    coord_cartesian(ylim=c(0, upper)) +
    geom_line(aes(x = Time, y = upr), color='darkgrey') + #lty=2
    geom_line(aes(x = Time, y = lwr), color='darkgrey') + 
    scale_x_continuous(breaks=seq(0,7, 1)) +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
    theme_bw()+
    #ggtitle(label=paste0(nam.order[i],  " ('", yr.dat$Year2,')'))
    ggtitle(label=nam.order[i])
  
  ##line for campaign
  if(!is.na(vline$y)){
    g1[[i]]<-g1[[i]] + 
      geom_vline(data=vline,aes(xintercept = y),
                 linetype='dashed', color = 'darkgrey', 
                 linewidth=0.5)}
  
  ## for Jombe, arrows for removals
  
  if(unique(dat.sub$Name) == 'KP Jombe'){
    
    mind<-min(dat.sub$midpoint)
    rdate<-as.Date(paste0(2020:2022, '-01-01'))
    remx<-as.numeric(rdate-mind)/365
    
    g1[[i]]<-g1[[i]] + 
      geom_segment(aes(x=remx[1]-0.5, y=0.01, xend=remx[1], yend=0.02),
                   arrow = arrow(length=unit(.2, 'cm')), lwd=1, col='darkgrey')+ 
      geom_segment(aes(x=remx[2]-0.5, y=0.01, xend=remx[2], yend=0.02),
                   arrow = arrow(length=unit(.2, 'cm')), lwd=1, col='darkgrey')+ 
      geom_segment(aes(x=remx[3]-0.5, y=0.01, xend=remx[3], yend=0.02),
                   arrow = arrow(length=unit(.2, 'cm')), lwd=1, col='darkgrey')
  }
  
  ###for all, add year of first survey
  g1[[i]]<-g1[[i]] + 
    geom_text(data=yr.dat,aes(x=x, y=y, label=Year), size=3)+
    theme(plot.margin = unit(marg.vec, "points"),
          plot.title = element_text(size = 10,
                                    margin=margin(b = 2, unit = "pt")),
          axis.title = element_text(size = 10),
          axis.text.x = element_text(size = 8))
}


# geom_hline(data=v_line2,aes(yintercept = yintercept), 
#            linetype='dashed', color = 'orangered', size=0.5)+
# theme_bw(base_size = 15)+
# theme(panel.spacing.x	= unit(0.1, "lines"))+
# theme(plot.margin = unit(c(5.5,.2,5.5,5.5), "points"))+
# theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 

##combine

p<-plot_grid(plotlist=g1, align='hv',nrow=4, ncol=4, byrow=TRUE)

jpeg('CPUE test figure.jpg', res=600, width=18, height=19, units = 'cm')
p
dev.off()


### make addl plots for Chumbe and Misali (Appendix figure A11)
sitc<-read_xlsx('data-raw/Report data.xlsx', sheet = 'Chumbe')
sitm<-read_xlsx('data-raw/Report data.xlsx', sheet = 'Misali')

##combine both sites, exclude data from latest reports (part of formal analysis)
sit<-rbind(sitc[1:4, c('Year', 'New caught', 'Effort', 'Reference')],
           sitm[1:6, c('Year', 'New caught', 'Effort', 'Reference')])
sit$Area<-c(rep('Chumbe', 4), rep('Misali', 6))
sit<-as.data.frame(sit)
sit$Year<-as.numeric(sit$Year)
sit$CPUE<-sit$`New caught`/sit$Effort

g1<-ggplot(sit, aes(x = Year, y = CPUE, color = Area)) +
  ##data points
  geom_point() +
  #labs(x = ifelse(i %in% c(10:13 ), 'Time', ''), 
  #     y = ifelse(i %in% c(1,5,9,13 ), "Weight", '')) +
  #coord_cartesian(xlim=c(0, 7), ylim=c(0, upper)) +
  scale_x_continuous(breaks=seq(1999,2018, 2)) +
  # scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  # geom_hline(data=hline,aes(yintercept = y),
  #            linetype='dashed', color = 'darkgrey', 
  #            linewidth=0.5) +
  theme_bw()+ 
  theme(legend.position="bottom")


jpeg('SIT trends.jpg', width=15, height=10, unit='cm', res=600)
g1
dev.off()
