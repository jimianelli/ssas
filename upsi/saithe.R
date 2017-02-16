source("../R/tools.R")

#----------------Run some single runs w/ different selectivity options plot----------------------------
fit <- Do_Run(Nretro=0,yrs_sel_change=seq(3,55,1))
fit <- Do_Run(Nretro=0,yrs_sel_change=seq(3,55,25))
names(fit[[1]])
#compute RMSE of surveys
(mean((log(fit[[1]]$U1obs/fit[[1]]$U1hat))^2))^.5
(mean((log(fit[[1]]$U2obs/fit[[1]]$U2hat))^2))^.5
#Make results in dataframe
df <- data.frame(Year=1961:2015,retro=rep(0,55),SSB=fit[[1]]$SSB, SSB.sd=fit[[1]]$SSB.sd, Fbar=fit[[1]]$Fbar,Fbar.sd=fit[[1]]$Fbar.sd, rec=fit[[1]]$rec,rec.sd=fit[[1]]$rec.sd )

#----------------Profile over different plus-group selectivity options---
nsel <- Do_Run(Nretro=0,rn="nSelAges_",nselages=5,yrs_sel_change=seq(3,55,5))
nsel
length(nsel[[1]]$SSB)
df <- data.frame(Year=1961:2015,nsel=rep(5,55),SSB=nsel[[1]]$SSB, SSB.sd=nsel[[1]]$SSB.sd, Fbar=nsel[[1]]$Fbar,Fbar.sd=nsel[[1]]$Fbar.sd, rec=nsel[[1]]$rec,rec.sd=nsel[[1]]$rec.sd )
for (i in 1:6){
	nselages <- i + 5
  nsel[i] <- Do_Run(Nretro=0,rn="nSelAges_",nselages=nselages,yrs_sel_change=seq(3,55,5))
  df        <- rbind(df,data.frame(Year=1961:2015,nsel=rep(nselages,55),SSB=nsel[[i]]$SSB, SSB.sd=nsel[[i]]$SSB.sd, Fbar=nsel[[i]]$Fbar,Fbar.sd=nsel[[i]]$Fbar.sd, rec=nsel[[i]]$rec,rec.sd=nsel[[i]]$rec.sd ))
}
df$nselages <- as.factor(df$nselages)
ggplot(df,aes(x=Year,y=SSB,color=nselages)) + geom_line(size=2) + mytheme + 

#----------------Run a retrospective and plot----------------------------
retro <- Do_Run(Nretro=10,yrs_sel_change=seq(3,55,5))
df <- data.frame(Year=1961:2015,retro=rep(0,55),SSB=retro[[1]]$SSB, SSB.sd=retro[[1]]$SSB.sd, Fbar=retro[[1]]$Fbar,Fbar.sd=retro[[1]]$Fbar.sd, rec=retro[[1]]$rec,rec.sd=retro[[1]]$rec.sd )
for (i in 0:10){
df <- rbind(df,data.frame(Year=1961:(2015-i),retro=rep(i,55-i),SSB=retro[[i+1]]$SSB,
	SSB.sd=retro[[i+1]]$SSB.sd, Fbar=retro[[i+1]]$Fbar,Fbar.sd=retro[[i+1]]$Fbar.sd, rec=retro[[i+1]]$rec,rec.sd=retro[[i+1]]$rec.sd )
 )
}
tail(df)
df$retro <- as.factor(df$retro)
ggplot(df,aes(x=Year,y=SSB,color=retro)) + geom_line(size=2) + mytheme + xlim(c(1970,2016)) + ylim(c(0,220))

