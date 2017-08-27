# create suitability layers for vegetable leaf miner from temperature-based mortality and development rate

library(tidyr)
library(raster)
library(sp)
library(ggplot2)
library(plyr)
library(maptools)  ## For wrld_simpl
library(rgdal)
library(rasterVis)
library(gridExtra) # also loads grid for lattice multi-plot
# global
source('./develop.fun.R')

setwd("C:\\Users\\james\\Documents\\git\\leafminer")

# development rate, 1/d, as a function of temperature, C
# https://cipotato.org/riskatlasforafrica/7-3-10-liriomyza-sativae-mujica-et-al/

# egg development
eggDev =function(temp){
  p  = 0.411
  To = 295.9
  Tl = 179.804
  Th = 313.431
  Ha = 16312.606
  Hl =-1180.651
  Hh = 526984.333
  x=temp+273 # convert to Kelvin
 dev.rate.out = p*(x/(To))*exp((Ha/1.987)*((1/To)-(1/x)))/(1+exp((Hl/1.987)*((1/Tl)-(1/x)))+exp((Hh/1.987)*((1/Th)-(1/x))))
  return(dev.rate.out)
}
# temp.C	dev.1.d^-1
eggdf<-data.frame(matrix(
  c(15.0877	,	0.133765	,
    20.0877	,	0.412082	,
    25.0877	,	0.343042	,
    30.0439	,	0.550162	,
    35.0439	,	1.08306	)
  ,ncol=2,byrow = TRUE))
names(eggdf)<-c('temp.C', "dev.rate.1/d")
temps = seq(0,55, length=1000)
png('plots/egg_dev.png',height = 4, width = 6, units = 'in', res = 300)
  plot(temps,eggDev(temps),'l', xlab='temp, (C)',ylab="egg dev. rate (1/d)", main = 'egg')
  points(eggdf)
dev.off()

# larval development 
larvaDev=function(temp){
  Dmin = 2.743
  Topt = 33.519
  K    = 0.118
  x = temp
  dev.rate.out = 2/(Dmin*(exp(K*(x-Topt))+exp(-K*(x-Topt)))) 
  return(dev.rate.out)
}
# temp.C	dev.1.d^-1
larvadf<-data.frame(matrix(
  c(15.0753	,	0.0798295	,
  20.0448	,	0.153306	,
  25.0138	,	0.219979	,
  30.0304	,	0.344861	,
  35.04	,	0.357107	)
,ncol=2,byrow = TRUE))
names(larvadf)<-c('temp.C', "dev.rate.1/d")
temps = seq(0,55, length=1000)
png('plots/larva_dev.png',height = 4, width = 6, units = 'in', res = 300)
plot(temps,larvaDev(temps),'l', xlab='temp, (C)',ylab="larval dev. rate (1/d)", main = 'larva')
points(larvadf)
dev.off()

# pupal devlopment
pupaDev=function(temp){
 p  = 0.049
 To = 288.52
 Tl = 284.349
 Th = 308.51
 Ha = 16661.255
 Hl =-192502.313
 Hh = 25785.133
 x=temp+273 # convert to Kelvin
 dev.rate.out = p*(x/(To))*exp((Ha/1.987)*((1/To)-(1/x)))/(1+exp((Hl/1.987)*((1/Tl)-(1/x)))+exp((Hh/1.987)*((1/Th)-(1/x))))
 return(dev.rate.out)
}
# temp.C	dev.1.d^-1
pupadf<-data.frame(matrix(
  c(15.0842	,	0.0534868	,
    20.0974	,	0.0726102	,
    25.0652	,	0.100795	,
    30.0306	,	0.145808	)
  ,ncol=2,byrow = TRUE))
names(pupadf)<-c('temp.C', "dev.rate.1/d")
temps = seq(0,55, length=1000)
png('plots/pupal_dev.png',height = 4, width = 6, units = 'in', res = 300)
plot(temps,pupaDev(temps),'l', xlab='temp, (C)',ylab="pupal dev. rate (1/d)", main='pupa')
points(pupadf)
dev.off()

#save dev fun plot accross all stages
eggdf$stage  <-'egg'
larvadf$stage<-'larva'
pupadf$stage <-'pupa'
obsdf<-rbind(eggdf,larvadf,pupadf)
names(obsdf)<-c('temp','dev.rate', 'stage')
preddfw<-data.frame(temp=temps,
                   egg = eggDev(temps), 
                   larva = larvaDev(temps),
                   pupa = pupaDev(temps)
                   )
preddf<-gather(preddfw,stage,dev.rate, -temp)
ggplot()+geom_line(data=preddf, aes(temp, dev.rate, colour = stage))+
  geom_point(data=obsdf, aes(temp, dev.rate, colour = stage))+
  ylab('development rate (1/d)')+
  xlab('temperature (C)')+
  theme_classic()
ggsave('plots/dev.rate.vs.temp.png')

# mortality functions
# egg mortality
eggMort=function(temp){
  a=  0.00122
  b= -0.071 
  c=  1.19
  x = temp
  mort.rate.out = a*x^2+b*x + c
  mort.rate.out[mort.rate.out>1]=1
  mort.rate.out[mort.rate.out<0]=0
  return(mort.rate.out)
}
# temp.C, mortality (%)
eggmortdf<-data.frame(matrix(
  c(15.0883	,	43.1914	,
    20.0973	,	20.19	,
    24.9955	,	15.1761	,
    30.0929	,	24.3397	,
    35.0759	,	15.5152	)
  ,ncol=2,byrow = TRUE))
names(eggmortdf)<-c('temp.C', "mort")
temps = seq(0,55, length=1000)
plot(temps,100*eggMort(temps),'l', xlab='temp, (C)',ylab="egg mort. (%)")
points(eggmortdf)

# larval mortality
larvaMort=function(temp){
  a= 0.002 
  b=-0.083 
  c= 0.977
  x = temp
  mort.rate.out = a*x^2+b*x + c
  mort.rate.out[mort.rate.out>1]=1
  mort.rate.out[mort.rate.out<0]=0
  return(mort.rate.out)
}
# temp (C), mortality (%)
larvamortdf<-data.frame(matrix(
  c(15.0962	,	14.2063	,
    20.0996	,	15.1365	,
    25.0151	,	4.1942	,
    30.1222	,	21.4109	,
    35.2386	,	41.2151	)
  ,ncol=2,byrow = TRUE))
names(larvamortdf)<-c('temp.C', "mort")
temps = seq(0,55, length=1000)
plot(temps,100*larvaMort(temps),'l', xlab='temp, (C)',ylab="larval mort. (%)")
points(df)

#pupal mortality
pupalMort=function(temp){
  Topt = 24.724 
  B    = 2.945
  H    = 0.077
  x = temp
  # typo in online website, 
  mort.rate.out =  1-1/exp(H*(1+exp(-(x-Topt)/B))*(1+exp(-(Topt-x)/B)))
  # mort.rate.out[mort.rate.out>1]=1
  # mort.rate.out[mort.rate.out<0]=0
  return(mort.rate.out)
}
# temp (C), mortality (%)
pupamortdf<-data.frame(matrix(
  c(15.0563	,	91.8697	,
    20.0215	,	36.67	,
    25.0578	,	39.7017	,
    30.0921	,	38.3127	,
    35.065	,	100.034	)
  ,ncol=2,byrow = TRUE))
names(pupamortdf)<-c('temp.C', "mort")
temps = seq(0,55, length=1000)
plot(temps,100*pupalMort(temps),'l', xlab='temp, (C)',ylab="pupal mort. (%)")
points(pupamortdf)


#save mortality plot accross all stages
eggmortdf$stage  <-'egg'
larvamortdf$stage<-'larva'
pupamortdf$stage <-'pupa'
mortobsdf<-rbind(eggmortdf,larvamortdf,pupamortdf)
names(mortobsdf)<-c('temp','mortality', 'stage')
mortpreddfw<-data.frame(temp=temps,
                    egg = eggMort(temps), 
                    larva = larvaMort(temps),
                    pupa = pupalMort(temps)
)
mortpreddf<-gather(mortpreddfw,stage,mortality, -temp)
ggplot()+geom_line(data=mortpreddf, aes(temp, 100*mortality, colour = stage))+
  geom_point(data=mortobsdf, aes(temp, mortality, colour = stage))+
  ylab('mortality (%)')+
  xlab('temperature (C)')+
  theme_classic()
ggsave('plots/mort.vs.temp.png')

# male senescence rate
adultSenM=function(temp){
  trid=384432.141 
  Tmax=38.567 
  Tmin=21.1 
  D=159910855.531 
  Dt=0.017 
  Smin=0.233
  x=temp
  mort.rate.out = trid*(((x-Tmin)^2)/((x-Tmin)^2+D)-exp(-(Tmax-(x-Tmin))/Dt))+Smin
  return(mort.rate.out)
}
# temp, (C), sen. rate (1/d)
senMdf<-data.frame(matrix(
  c(15.1261	,	0.289135	,
    20.0452	,	0.228492	,
    25.0736	,	0.233406	,
    30.0254	,	0.392278	)
  ,ncol=2,byrow = TRUE))
names(senMdf)<-c('temp, (C)', "sen rate, 1/d")

# female senescence rate
adultSenF=function(temp){
  b1=0.028 
  b2=0.057
  x=temp
  mort.rate.out = b1*exp(b2*x)
  return(mort.rate.out)
}
# temp, (C), sen. rate (1/d)
senFdf<-data.frame(matrix(
  c(15.0132	,	0.0881459	,
    20.0121	,	0.0805471	,
    25.0641	,	0.0942249	,
    30.0258	,	0.171733)
  ,ncol=2,byrow = TRUE))
names(senFdf)<-c('temp, (C)', "mort. (%)")
temps = seq(0,55, length=1000)
png('plots/senescence_temp.png',height = 4, width = 6, units = 'in', res = 300)
plot(temps,adultSenF(temps),'l', col = 'red',
     xlab='temp, (C)',ylab="sen. rate (1/d)")
points(senFdf, col = 'red')
lines(temps,adultSenM(temps),'l', col='blue')
points(senMdf, col = 'blue')
legend( 40,.2,c("male","female"), col = c('blue','red'), lty = rep(1,2))
dev.off()

# female fecundity
fecFun=function(temp){
  b1=-64.662 
  b2=-1.334 
  b3=32.137
  x=temp
  mort.rate.out = exp(b1 + b2*x + b3*log(x))
  return(mort.rate.out)
}
# temp C,, total eggs laid
df<-data.frame(matrix(
  c(15.0447	,	-0.244035	,
    20.1251	,	141.279	,
    25.2081	,	220.029	,
    30.1191	,	110.095	,
    34.9846	,	-0.567472	)
  ,ncol=2,byrow = TRUE))
names(df)<-c('temp, (C)', "mort. (%)")
temps = seq(0,55, length=1000)
png('plots/female_fecundity_vs_temp.png',height = 4, width = 6, units = 'in', res = 300)
plot(temps,fecFun(temps),'l', xlab='temp, C',ylab="total eggs laid, eggs/ind")
points(df)
dev.off()

# intrinsic rate of increase
r.d<-as.data.frame(read.csv('C:/Users/james/Dropbox (cesar)/1303CR2 - Leafminer/Literature/ecophysiology/Mujica et al/intrinsic_rate_of_increase.csv'))
r.d$temp2<-r.d$temp^2
lm1<-lm(r ~ temp + temp2, data= r.d[r.d$source == 'Mujica2016',])
preds<-data.frame(temp<-seq(0,55, length=1000),temp2<-seq(0,55, length=1000)^2)
preds$preds<-predict(lm1,preds)
ggplot()+geom_point(data=r.d, aes(temp,r,colour=source))+geom_line(data=preds,aes(temp, preds))
r.fun=function(temp){
  a = -1.089141111       
  b =  0.093601269         
  c = -0.001751445  
  x=temp
  r.out = a + b*temp + c*temp^2 
  return(r.out)
}

########################################

# load temperature data, raster stack of 365 days of the year
TMIN<-brick('E:/AWAP_daily_downloaded/mu_Tmin_for_DOY_ag10.tif')
TMAX<-brick('E:/AWAP_daily_downloaded/mu_Tmax_for_DOY_ag10.tif')
# TMIN<-brick('E:/AWAP_daily_downloaded/mu_Tmin_for_DOY.tif')
# TMAX<-brick('E:/AWAP_daily_downloaded/mu_Tmax_for_DOY.tif')

# creat QLD spatial polygon
if(!exists('mapdata')){
  # load map data each id (state) has a certain number of pieces so select what you need to speed computation
  aus <- readOGR(dsn = "C:/Users/james/Dropbox (Personal)/Programming/R/Maps/Australia by state ABS", layer = "STE11aAust")
  aus_smpl<-rgeos::gSimplify(spgeom = aus,tol = 0.005,topologyPreserve = FALSE)
  plot(aus_smpl)
  mapdata <- fortify(aus_smpl)
}

SPDF = SpatialPolygonsDataFrame(aus_smpl, data.frame(state = 0:8, row.names = c("0","1", "2", "3", "4", "5", "6", "7", "8")))
png('plots/Australia.png',height = 4, width = 6, units = 'in', res = 300)
plot(SPDF)
dev.off()
QLDmask <- subset(SPDF, state==2)
QLDmask<-spTransform(QLDmask, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
# crop temperature stacks, or load if high res
Tmin = TMIN
Tmax = TMAX

# if(TMIN@ncols==886){
#   Tmin <- raster('E:/AWAP_daily_downloaded/mu_Tmin_for_DOY_QLD')
#   Tmax <- raster('E:/AWAP_daily_downloaded/mu_Tmax_for_DOY_QLD')
# }else{
#   Tmin<-crop(TMIN,QLDmask)
#   Tmax<-crop(TMAX,QLDmask)
# }

# plot Seisia conditions
seisia<- dismo::geocode('Seisia, QLD')
seisiaTmin<-raster::extract(TMIN, cbind(seisia$longitude,seisia$latitude))
seisiaTmax<-raster::extract(TMAX, cbind(seisia$longitude,seisia$latitude))
seisiaDfw<-data.frame(date=as.Date('2017-01-01')+0:364,
               Tmin = as.numeric(seisiaTmin),
               Tmax = as.numeric(seisiaTmax))
seisiaDf<-gather(seisiaDfw, key, temp,-date)

# plot Cairns conditions
cairns<- dismo::geocode('Cairns, QLD')
cairnsTmin<-raster::extract(TMIN, cbind(cairns$longitude,cairns$latitude))
cairnsTmax<-raster::extract(TMAX, cbind(cairns$longitude,cairns$latitude))
cairnsDfw<-data.frame(date=as.Date('2017-01-01')+0:364,
                      Tmin = as.numeric(cairnsTmin),
                      Tmax = as.numeric(cairnsTmax))
cairnsDf<-gather(cairnsDfw, key, temp,-date)


# plot Melbourne conditions
melb<- dismo::geocode('Melbourne, VIC')
melbTmin<-raster::extract(TMIN, cbind(melb$longitude,melb$latitude))
melbTmax<-raster::extract(TMAX, cbind(melb$longitude,melb$latitude))
melbDfw<-data.frame(date=as.Date('2017-01-01')+0:364,
                      Tmin = as.numeric(melbTmin),
                      Tmax = as.numeric(melbTmax))
melbDf<-gather(melbDfw, key, temp,-date)

# get lower and upper temp limits for each life stage
# assume <5% cohort survival to be critical limit
x1<-rep(min(seisiaDf$date),3)
x2<-rep(max(seisiaDf$date),3)
crit.lim = 0.9
y1<-c(uniroot(function(x)   eggMort(x) - crit.lim, interval = c(0,30))$root,
      uniroot(function(x) larvaMort(x) - crit.lim, interval = c(0,30))$root,
      uniroot(function(x)  pupalMort(x) - crit.lim, interval = c(0,30))$root)

y2<-c(uniroot(function(x)   eggMort(x) - crit.lim, interval = c(30,60))$root,
      uniroot(function(x) larvaMort(x) - crit.lim, interval = c(30,60))$root,
      uniroot(function(x)  pupalMort(x) - crit.lim, interval = c(25,60))$root)

# plot Seisia
boxes<-data.frame(x1=x1,x2=x2, y1 = y1, y2 =  y2, stage = c('egg','larva','pupa'))
ggplot()+
  geom_rect(data=boxes, aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=stage),alpha=0.2)+
  geom_line(data=seisiaDf,aes(x=date, y=temp, colour=key), size = 2)+ggtitle("Seisia, QLD")+
  scale_x_date(date_breaks = "1 month", date_labels = "%b")+theme_classic()
ggsave('plots/temp_Seisia.png')
ggplot()+
  geom_line(data=seisiaDf,aes(x=date, y=temp, colour=key), size = 2)+ggtitle("Seisia, QLD")+
  scale_x_date(date_breaks = "1 month", date_labels = "%b")+theme_classic()
ggsave('plots/temp_Seisia_no_temprange.png')

# plot Cairns
ggplot()+
  geom_rect(data=boxes, aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=stage),alpha=0.2)+
  geom_line(data=cairnsDf,aes(x=date, y=temp, colour=key), size = 2)+ggtitle("Cairns, QLD")+
  scale_x_date(date_breaks = "1 month", date_labels = "%b")+theme_classic()
ggsave('plots/temp_Cairns.png')

# plot Melbourne
ggplot()+
  geom_rect(data=boxes, aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=stage),alpha=0.2)+
  geom_line(data=melbDf,aes(x=date, y=temp, colour=key), size = 2)+ggtitle("Melbourne, VIC")+
  scale_x_date(date_breaks = "1 month", date_labels = "%b")+theme_classic()
ggsave('plots/temp_Melbourne.png')

#plot box 
ggplot()+
  geom_rect(data=boxes, aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=stage),alpha=0.2)+
  scale_x_date(date_breaks = "1 month", date_labels = "%b")+
  theme_classic()+theme( panel.background = element_rect(fill = "transparent",colour = NA), # or theme_blank()
                         panel.grid.minor = element_blank(), 
                         panel.grid.major = element_blank(),
                         plot.background = element_rect(fill = "transparent",colour = NA))
ggsave('plots/temp_pref.png')

# simulate development
startDate<-as.Date('2017-01-01')
startDay<-as.numeric(format(startDate,'%j'))
startStage<-rep(1,Tmin@ncols*Tmin@nrows)
# TMAX =Tmax +0.00001
# TMIN =Tmin +0.00001
dev.fun = list(egg = eggDev, larva=larvaDev, pupa=pupaDev)
mort.fun = list(egg = eggMort, L1 = larvaMort, pupa = pupalMort)

# run simulation
data<-develop(Tmax,Tmin, startDay, startStage, dev.fun, mort.fun, r.fun)
# save(data,file='hires_temp_simulations.Rdata')
# blank raster
r<- Tmax[[1]]
projection(r)<-'+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
projection(aus_smpl)<-'+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
## crop and mask
QLDplot<-function(data, title, hide.legend = FALSE, legRange=NULL, plotAus = FALSE){
  r[]<-data
  if(plotAus){
    r1 <- mask(r, aus_smpl)
  }else{
    r1 <- crop(r, QLDmask)
    r1 <- mask(r1, QLDmask)
  }
  colr <- colorRampPalette(brewer.pal(11, 'RdYlGn'))
  if(is.null(legRange)){
    minLeg<-floor(min(r1[],na.rm=TRUE))
    maxLeg<-ceiling(max(r1[],na.rm=TRUE))
  }else{
    if(length(legRange)!=2)stop('legRange min max not specified')
    minLeg<-legRange[1]
    maxLeg<-legRange[2]
  }
  plotout<-levelplot(r1, 
            margin=FALSE,                       # suppress marginal graphics
            colorkey=ifelse(hide.legend,FALSE,list(
              space='bottom',                   # plot legend at bottom
              labels=list(at=seq(minLeg,maxLeg,length = 3), font=4)      # legend ticks and labels 
            )),
            main = title,
            par.settings=list(
              axis.line=list(col='transparent') # suppress axes and legend outline
            ),
            scales=list(draw=FALSE),            # suppress axis labels
            col.regions=colr,                   # colour ramp,
            xlab = NULL, ylab=NULL,
            at=seq(minLeg, maxLeg, len=101))           # colour ramp breaks
    
  if(!plotAus)plotout<-plotout+layer(sp.polygons(QLDmask, lwd=1, col = 'black'))           # add oregon SPDF with latticeExtra::layer
  if( plotAus)plotout<-plotout+layer(sp.polygons(aus_smpl, lwd=1, col = 'black'))           # add oregon SPDF with latticeExtra::layer
  plotout
}


# save plot of generations per year
png('plots/generations.png',height = 4, width = 4, units = 'in', res = 300)
print(QLDplot(data$genCount, title = 'generations per year'))
dev.off()

# save plot of survival index
png('plots/survival_index.png',height = 4, width = 12, units = 'in', res = 300)
plot1<-(QLDplot(data$surv[,1], title = 'egg survival index', hide.legend = TRUE, plotAus = TRUE))
plot2<-(QLDplot(data$surv[,2], title = 'larval survival index',hide.legend = TRUE))
plot3<-(QLDplot(data$surv[,3], title = 'pupal survival index',hide.legend = TRUE))
grid.arrange(plot1,plot2, plot3,ncol=3)
dev.off()

# plot monthly plot of dev potential
rs<-stack(r,r,r,r,r,r,r,r,r,r,r,r)
month.names<-c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')
for(month in 1:12){
  rs[[month]][]<-data$popSizeMonthly[,month]  
  png(paste0('plots/pop_growth/pop_growth_',stringr::str_pad(month,pad = '0',width=2),month.names[month],'.png'),
      height = 4, width = 4, units = 'in', res = 300)
  print(QLDplot(data$popSizeMonthly[,month], title = paste0('pop. growth potential (',month.names[month],')'), legRange=c(0,70), plotAus = TRUE))
  dev.off()
}
monthdata<-data.frame(
  month = 1:12,
  month.names = month.names,
  cairns = as.numeric(raster::extract(rs,cbind(cairns$longitude,cairns$latitude))[1,]),
  seisia = as.numeric(raster::extract(rs,cbind(seisia$longitude,seisia$latitude))[1,])
)
monthdf<-gather(monthdata,location, pop.growth, -month, -month.names)
monthdf$month.names<-factor(monthdf$month.names, levels = unique(monthdf$month.names[order(monthdf$month)]))
ggplot()+geom_line(data=monthdf,aes(factor(month.names),pop.growth,group=location, colour=location))+
  ylab('monthly pop. growth (factor increase)')+
  xlab('')+theme_classic()
ggsave('plots/pop.growth.seisia.cairns.png')
# save plot of popsize 
print(QLDplot(data$crit.mort[,1], title = 'egg critical mortality index'))
# print(QLDplot(data$crit.mort[,2], title = 'larval critical mortality index'))
print(QLDplot(data$crit.mort[,3], title = 'pupal critical mortality index'))







