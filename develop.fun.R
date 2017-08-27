develop<-function(Tmax, Tmin, startDay, startStage, dev.funs, mort.funs, r.funs){
  # arguments:
  # Tmax: a raster stack or list of daily max temperatures for 365 days [365][D]
  # Tmin: a raster stack or list of daily min temperatures for 365 days [365][D]
  # startStage: a vector specifying the stage of development at each deme [D]
  # startDay: a scalar specifying the day developmnetal stages were observed [1]
  # 'name': the name of the bug
  # 'dev.funs': temperature dependent temperature function for each stage

  # output
  # genCount: count of generations 
  # some checks
  if(length(startDay)!=1)stop('startDay must have length 1')
  if(length(Tmax[[1]])!=length(startStage))stop('number of locations must match length of startStage')
  if(length(names(Tmax))!=length(names(Tmin)))stop('size of climate data Tmin and Tmax must match')
  if(length(dev.funs)!=length(mort.funs))stop('dev.funs and mort.funs must match')
  
  # some useful variables
  D <- length(Tmax[[1]])
  S <- length(dev.funs)
  
  # hourly temp as a function of Tmin, Tmax, time (of day),d
  hrtemp<-function(TMIN,TMAX,hr){
    # Tmax: daily maximum temp, C
    # Tmin: daily minimum temp, C
    # hr: hour of day, h
    (TMAX-TMIN)/2*sin(2*pi*(hr - 8)/24 )  + (TMAX +TMIN)/2
  }
  
  hr.dev <- function(yearHour, TMIN, TMAX, curStage, dev.funs){
    # calculate hourly development accross all demes
    hr<-yearHour%%24
    if(length(hr)>1)stop('more than one hour supplied')
    if(all(c(length(TMIN), length(TMAX)) != length(curStage)) )stop('size of TMAX, TMIN, curStage must match')
    S <- length(dev.funs) # stages
    D <- length(TMIN)
    
    # get temp for given hr accross all demes
    temp<-hrtemp(TMIN,TMAX,hr)
    # intialise dev matrix 
    dev<-rep(0,D)
    # iterate through stages
    for(s in 1:S){
      # stage specific function for temp dependence of dev. rate
      dev[curStage==s]<-dev.funs[[s]](temp[curStage==s])/24 # increment dev during hr for stage
    }
    return(dev) # return hourly development
  }
  
  # initial some state vectors
  curStage<-startStage # current stage of insect
  curDev  <-rep(0, D) # proportion completed of current stage
  genCount<-rep(0, D) # count of annual generations
  popSizeMonthly <- array(0,c(D,12)) # store pop size for each month
  popSize<-rep(1,D) # starting popsize
  surv<-array(0,c(D,S)) # sum of all hourly surv index
  crit.mort<-array(0,c(D,S)) # sum of all hourly mort index > 0.95
  
  # forward simulation
  h = ((startDay-1)*24)
  day = startDay # intialise day
  while(day<=365){
    date=as.Date('2017-01-01')+day-1
    TMIN<-Tmin[[day]][] # index cannot be 0
    TMAX<-Tmax[[day]][]
    for (dayhour in 0:23){
      h<-(day-1)*24+dayhour
      curDev<-curDev+hr.dev(h, TMIN, TMAX, curStage, dev.funs) # increment dev
      curDev[is.na(curDev)]<-0 # na's to 0
      curStage[curDev>1]<-curStage[curDev>1] + 1 # curStage
      D_st <- which(curDev>1) # demes with stage transitions
      S_st <- curStage[curDev>1] # new stages 
      curDev[curDev>1]<-0 # if new stage, set dev to zero
      # if current stage is higher than final stage, increment generation count
      genCount[curStage==(S+1)]<- genCount[curStage==(S+1)] + 1
      curStage[curStage==(S+1)]<-1 # if current stage is higher than final stage, reset to 1
      
      # get temp for given hr accross all demes
      temp<-hrtemp(TMIN,TMAX,dayhour)
      # iterate through stages
      for(s in 1:S){
        s_surv<-1-mort.funs[[s]](temp)
        surv[,s]<-surv[,s]+s_surv
        # count the number of hours where survival index drops below 5%
        crit.mort[s_surv<0.05,s]<-crit.mort[s_surv<0.05,s]+1
      }
      # if first day of month, set popSize to one
      if((format(date,'%d')=='01')){
        popSize<-rep(1,D)  
      }
      popSize<-popSize+popSize*r.funs(temp)/24
    }
    # if 28th day of month, store pop size for month
    if(format(date,'%d')=='28'){
      month<-as.numeric(format(date,'%m'))
      popSizeMonthly[,month]<-popSize  
    }
    day = day+1 # increment day
  }
  
  return(list(genCount=genCount,surv=surv,crit.mort=crit.mort, popSizeMonthly=popSizeMonthly))
}