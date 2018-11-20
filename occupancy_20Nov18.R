## ensure that the working directory has list of India's birds with scientific names 
## (just a safety mechanism for the function to work for small subsets, needs to be enabled if required)
## only need to input data, the species of interest and the complete list of India's bird species
## also groupspecs if required (a dataframe with all relevant list level info), it is defaulted to data

expandbyspecies = function(data, species)
{
  require(tidyverse)
  
  data = data %>%
    mutate(timegroups = as.character(year)) %>%
    mutate(timegroups = ifelse(year < 1990, "before 1990", timegroups)) %>%
    mutate(timegroups = ifelse(year >= 1990 & year <= 1999, "1990-1999", timegroups)) %>%
    mutate(timegroups = ifelse(year > 1999 & year <= 2005, "2000-2005", timegroups)) %>%
    mutate(timegroups = ifelse(year > 2005 & year <= 2010, "2006-2010", timegroups)) %>%
    mutate(timegroups = ifelse(year > 2010 & year <= 2013, "2011-2013", timegroups)) %>%
    mutate(timegroups = ifelse(year == 2014, "2014", timegroups)) %>%
    mutate(timegroups = ifelse(year == 2015, "2015", timegroups)) %>%
    mutate(timegroups = ifelse(year == 2016, "2016", timegroups)) %>%
    mutate(timegroups = ifelse(year == 2017, "2017", timegroups)) %>%
    mutate(timegroups = ifelse(year == 2018, "2018", timegroups))
  
  data$timegroups = as.factor(data$timegroups)
  
  ## considers only complete lists
  
  checklistinfo = data %>%
    distinct(gridg1,gridg2,gridg3,gridg4,gridg5,DISTRICT,ST_NM,
             LOCALITY.ID,LOCALITY.TYPE,LATITUDE,LONGITUDE,OBSERVATION.DATE,TIME.OBSERVATIONS.STARTED,
             OBSERVER.ID,PROTOCOL.TYPE,DURATION.MINUTES,EFFORT.DISTANCE.KM,NUMBER.OBSERVERS,ALL.SPECIES.REPORTED,
             group.id,month,year,day,week,fort,LOCALITY.HOTSPOT,no.sp,timegroups)
  
  checklistinfo = checklistinfo %>%
    filter(ALL.SPECIES.REPORTED == 1) %>%
    group_by(group.id) %>% slice(1) %>% ungroup
  
  ## expand data frame to include all bird species in every list
  
  expanded = checklistinfo
  expanded$COMMON.NAME = species
  
  ## join the two, deal with NAs next
  
  expanded = left_join(expanded,data)
  
  ## deal with NAs
  
  expanded = expanded %>% mutate(OBSERVATION.COUNT = replace(OBSERVATION.COUNT, is.na(OBSERVATION.COUNT), "0"))
  
  
  expanded = expanded %>%
    mutate(OBSERVATION.NUMBER=OBSERVATION.COUNT,
           OBSERVATION.COUNT=replace(OBSERVATION.COUNT, OBSERVATION.COUNT != "0", "1"))
  
  
  
  expanded$OBSERVATION.COUNT = as.numeric(expanded$OBSERVATION.COUNT)
  expanded$OBSERVATION.NUMBER = as.numeric(expanded$OBSERVATION.NUMBER)
  
  return(expanded)
}

############################################################

## occupancy analyses for bird abundance/range
###Requires tidyverse, reshape2, data.table and unmarked####

occufreq = function(data, species, resolution)
{
  require(tidyverse)
  require(reshape2)
  require(data.table)
  require(unmarked)
  
  # create dataframe to store occupancy and detection proabability estimates across species and spatial resolutions
  est = array(data=NA,dim=c(length(species),3,length(resolution)),
              dimnames=list(species,c("detprob","occ","occ.se"),c("20km","40km","60km","80km")))
  
  for(s in 1:length(species))
  {selexp = expandbyspecies(data,species[s])
  
  eff = quantile(selexp$EFFORT.DISTANCE.KM, 0.95, na.rm=TRUE)
  selexp = selexp %>%
  filter(EFFORT.DISTANCE.KM < eff, year > 2013)
  selexp = selexp[sample(1:nrow(selexp)),]
  
  selexp$month[selexp$month %in% c(11,12,1,2)] = "Win"
  selexp$month[selexp$month %in% c(3,4,5,6)] = "Sum"
  selexp$month[selexp$month %in% c(7,8,9,10)] = "Mon"
  
  for(r in 1:length(resolution))
  {
    if(r == 1)
    {lpg = selexp %>%
     group_by(gridg1) %>% summarize(lpg = n())
     listcutoff = quantile(lpg$lpg, 0.95, na.rm=TRUE)
      selexp = selexp %>% 
      arrange(gridg1) %>%
      mutate(gridg = as.numeric(gridg1)) %>%
      group_by(gridg) %>% mutate(group.id = 1:n())}
    
    if(r == 2)
    {lpg = selexp %>%
      group_by(gridg2) %>% summarize(lpg = n())
      listcutoff = quantile(lpg$lpg, 0.95, na.rm=TRUE)
      selexp = selexp %>% 
      arrange(gridg2) %>%
      mutate(gridg = as.numeric(gridg2)) %>%
      group_by(gridg) %>% mutate(group.id = 1:n())}
    
    if(r == 3)
    {lpg = selexp %>%
      group_by(gridg3) %>% summarize(lpg = n())
      listcutoff = quantile(lpg$lpg, 0.95, na.rm=TRUE)
      selexp = selexp %>% 
      arrange(gridg3) %>%
      mutate(gridg = as.numeric(gridg3)) %>%
      group_by(gridg) %>% mutate(group.id = 1:n())}
    
    if(r == 4)
    {lpg = selexp %>%
      group_by(gridg4) %>% summarize(lpg = n())
      listcutoff = quantile(lpg$lpg, 0.95, na.rm=TRUE)
      selexp = selexp %>% 
      arrange(gridg4) %>%
      mutate(gridg = as.numeric(gridg4)) %>%
      group_by(gridg) %>% mutate(group.id = 1:n())}
    
    setDT(selexp)
    
    det = dcast(selexp, gridg ~ group.id, value.var = "OBSERVATION.COUNT")
    cov.month = dcast(selexp, gridg ~ group.id, value.var = "month")
    cov.nosp = dcast(selexp, gridg ~ group.id, value.var = "no.sp")
    
    det = setDF(det)
    cov.month = setDF(cov.month)
    cov.nosp = setDF(cov.nosp)
    
    det = det[,1:listcutoff]
    cov.month = cov.month[,1:listcutoff]
    cov.nosp = cov.nosp[,1:listcutoff]
    
    umf = unmarkedFrameOccu(y=det[,-1], siteCovs =NULL, obsCovs = list(cov1 = cov.nosp[,-1], cov2 = cov.month[,-1]))
    
    occ_det = occu(~cov1+cov2 ~1, data=umf)
    g = backTransform(occ_det, type="state")
    
    newdat = data.frame(cov1=20, cov2=factor(c("Mon","Win","Sum")))
    f = predict(occ_det, newdata = newdat, type = "det")
    f = mean(f$Predicted)
     
  est[s,"detprob",r] =  f
  est[s,"occ",r] = g@estimate
  }
}  
  return(est)
}

species = c("Indian Peafowl","Ashy Prinia")
resolution = c("g1","g2","g3", "g4")

occufreq(data, species, resolution)