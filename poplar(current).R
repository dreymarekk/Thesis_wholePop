######################################################################################################################
#
# Poplar Island Mark-Recapture Analysis
# Survival Rate Estimation
# 5/23/2018
# David Jenkins
#
######################################################################################################################


library(RMark)
library(reshape)
library(dplyr)
library(ggplot2)

makeHistory<-function(x){                     
  k<-ncol(x)                                  
  n<-nrow(x)                                    
  out<-array(dim=n)                             
  
  for (i in 1:n){                               
    y<-(x[i,]>0)*1                              
    out[i]<-paste(y, collapse = "")}            
  
  return(out)                                   
  
}                                


#setwd("")


raw<-read.csv("wild3.csv")
raw<-raw[-6003:-6012,]

#df of first capture per animal
initcaps<-distinct(raw,pit, .keep_all = TRUE)

#convert dates to date obj
lct <- Sys.getlocale("LC_TIME"); Sys.setlocale("LC_TIME", "C")
raw$do_cap<-as.Date(raw$do_cap, "%m/%d/%Y")
Sys.setlocale("LC_TIME", lct)
rm(lct)

#extract yoc
raw<-mutate(raw, year = format(do_cap, "%Y"))

#group all ages for which there is no yob
raw$age[is.na(raw$age)] = 25

#shape enc hist object
junk<-melt(mutate(raw, detect = 1),id.var=c("pit","year"),measure.var="detect")
y=cast(junk, pit~year)


#replace any integers >1 with 1
for (i in 1:nrow(y)){
  for (j in 2:15){
    if(y[i,j] > 1){
      y[i,j]<-1
    }
  }
}

init_sort<-initcaps[order(initcaps$pit),]
y_sort<-y[order(y$pit),]

#attach covariates and temporary variables

y<-cbind(y_sort, init_sort$type, init_sort$yob, init_sort$sex, init_sort$age)
names(y)[16:19] <- c("type", "yob", "sex", "age")

#remove any MH where yob is unknown
y<-y[-which(y$type == "MH" & is.na(y$yob)),]

#give all MH initial encounters in yob
for (i in 1:nrow(y)){
  if (y[i, 16] == "MH"){
    
    if(y[i,17] == 2004){
      y[i,2]<-1
    }
    else if(y[i,17] == 2005){
      y[i,3]<-1
    }
    else if(y[i,17] == 2006){
      y[i,4]<-1
    }
    else if(y[i,17] == 2007){
      y[i,5]<-1
    }
    else if(y[i,17] == 2008){
      y[i,6]<-1
    }
    else if(y[i,17] == 2009){
      y[i,7]<-1
    }
    else if(y[i,17] == 2010){
      y[i,8]<-1
    }
    else if(y[i,17] == 2011){
      y[i,9]<-1
    }
    else if(y[i,17] == 2012){
      y[i,10]<-1
    }
    else if(y[i,17] == 2013){
      y[i,11]<-1
    }
    else if(y[i,17] == 2014){
      y[i,12]<-1
    }
    else if(y[i,17] == 2015){
      y[i,13]<-1
    }
    else if(y[i,17] == 2016){
      y[i,14]<-1
    }
  }
}

#give MH initial age 0
for (i in 1:nrow(y)){
  if (y[i,16] == "MH"){
    y[i,19]<-"0"
  }
}
#################################################################################################################
#At this point analysis from HS.R must have been run up until construction of hs encounter history data frame
#If it has not, run HS.R up until indicated point and then save object as z.hs
#
#Note to self: this is very sloppy programming and should be changed at a later date
#Response to self: has been changed, z.hs now .csv in project directory. Ignore above warning
#################################################################################################################
read.csv("z.hs.csv")
z.hs<-mutate(z.hs, age = 0)
library(tibble)
z.hs<-add_column(z.hs, 0, .before = "2005")
names(z.hs)[2]<-"2004"

x<-rbind(y,z.hs)
x$age<-as.integer(capt.hist$age)

capt.hist<-data.frame(ch = makeHistory(x[,2:15]))

capt.hist<-cbind(capt.hist, x[,18:19])

# fix sex errors
for(i in 1:nrow(capt.hist)){
  
  if (capt.hist[i,2] == "j" | is.na(capt.hist[i,2]))
    capt.hist[i,2] = "J"
}

capt.hist$sex<-factor(capt.hist$sex)

capt.hist.processed<-process.data(capt.hist,model="CJS",begin.time=2004, groups=c("sex", "age"), age.var=2, 
                                       initial.age=c(0,1,2,3,4,5,6,7,8,9,10,11,12,14,25))

ddl<-make.design.data(capt.hist.processed, 
                           parameters = list(Phi=list(age.bins = c(0,4,8,40)),p=list(age.bins = c(0,4,8,40))), 
                           right = FALSE)
#Model definitions
Phi.dot = list(formula=~1)
Phi.time = list(formula=~time)
Phi.sex = list(formula=~sex)
Phi.age = list(formula=~age)
Phi.sex.age =list(formula=~sex*age)
Phi.sex.time = list(formula=~sex*time)
Phi.sex.time.age = list(formula=~sex*time*age)
Phi.age.time = list(formula=~age*time)

p.dot = list(formula=~1)
p.time = list(formula=~time)
p.sex = list(formula=~sex)
p.age = list(formula=~age)
p.sex.age =list(formula=~sex*age)
p.sex.time = list(formula=~sex*time)
p.sex.time.age = list(formula=~sex*time*age)
p.age.time = list(formula=~age*time)

model.list=create.model.list(model="CJS")
results=mark.wrapper(model.list, data=capt.hist.processed,ddl=ddl)
results
PIMS(results$Phi.sex.time.age.p.sex.time.age, parameter = "Phi")


#############################################################
#Extract and plot results


beta<-mutate(results$Phi.sex.time.age.p.sex.time.age$results$real, 
             index = 1:nrow(results$Phi.sex.time.age.p.sex.time.age$results$real))
phi<-filter(beta, index<= 117)

pim.translation<-results$Phi.sex.time.age.p.sex.time.age$simplify$pim.translation

design<-ddl$Phi
for(i in 1:3185){
  design[i,1] <- pim.translation[i]
}

ref<-select(design, par.index, age, time, sex)
ref<-group_by(ref, par.index)
ref<-distinct(ref, .keep_all = TRUE)

ggbest.survival<-bind_cols(select(phi, estimate, se, lcl, ucl), select(ref, age,time,sex))


ggplot(ggbest.survival, aes (x=age, y=estimate, color=sex))+
  geom_point(position = position_dodge(width=0.5))+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=0, position = position_dodge(width=0.5))+
  facet_wrap(~time)+
  theme_bw()

############################################################
#Without time effect

beta2<-mutate(results$Phi.sex.age.p.sex.time.age$results$real, 
             index = 1:nrow(results$Phi.sex.age.p.sex.time.age$results$real))
phi2<-filter(beta2, index<= 9)

pim.translation2<-results$Phi.sex.age.p.sex.time.age$simplify$pim.translation

design2<-ddl$Phi
for(i in 1:3185){
  design[i,1] <- pim.translation2[i]
}

ref<-select(design, par.index, age, sex)
ref<-group_by(ref, par.index)
ref<-distinct(ref, .keep_all = TRUE)

ggbest.survival2<-bind_cols(select(phi2, estimate, se, lcl, ucl), select(ref, age,sex))

ggbest.survival2$sex <-factor(ggbest.survival2$sex, levels = c("J", "M", "F"))


Sex<- list(
  'J' = "Sexually Immature",
  'M' = "Male",
  'F' = "Female"
)

sex_labeller<- function(variable, value){
  return(Sex[value])
}


ggplot(ggbest.survival2, aes(x=age, y=estimate))+facet_wrap(~sex, labeller = sex_labeller)+
  geom_point()+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=.1)+
  theme_bw()+scale_y_continuous(limits = c(0,1))+
  scale_shape_manual(values = c(1, 17, 16))+
  ylab(expression(Apparent~Survival~(phi)))+xlab("Age Class")
#+scale_color_manual(values = c("black", "blue", "red"))
