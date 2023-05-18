library(devtools)
#set working directory to workbook location
path<-dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)
#load the large matrix A' from c-code compiled package 
install_github("schoolmasterd/sigma-matrix")
library(Sigma)

#create distribution stuff#
#create pdf
fq<-function(mu,sig,y)-(exp(-((mu-log(1/(1-y))-log(y))^2/(2*sig^2)))*(-1/(1-y)-1/y))/(sqrt(2*pi)*sig)
#create cdf
Fq<-function(mu,sig,y)sapply(y,function(y)pnorm(log(1/(1-y))+log(y),mu,sig))
#create rng
randFq<-function(x){
  rnd<-runif(1)
  est<-integrate(function(y)y*fq(x[1],x[2],y),lower = 0,upper = 1)
  ans=uniroot(function(y)Fq(x[1],x[2],y)-rnd,interval = c(0,1))$root
  return(ans)
}

#get fitted parmaters from both species
svPars <- read.csv("Data/SilverParams.csv",header = F)
bhPars <-read.csv("Data/BigHeadParams.csv",header = F)

#find means
Ebh=mapply(function(mu,sd)integrate(function(y)y*fq(mu,sd,y),lower = 0,upper = 1)$value,bhPars[,1],bhPars[,2])
Esv=mapply(function(mu,sd)integrate(function(y)y*fq(mu,sd,y),lower = 0,upper = 1)$value,svPars[,1],svPars[,2])


#order of arguments expected in A'
arg<-"q12,q23,q34,q45,q56,q21,q32,q43,q54,q65,b0"

#try it with BH data and b0=0
do.call('Sigma',as.list(c(Ebh,0)))

#use the randFq function to grab random values for each movement prob
pars<-unlist(apply(bhPars, 1, list), recursive = FALSE)

#now use this construction to get random values
sapply(pars,randFq)

#setup weight vector to modify barriers
len<-100
wts<-rep(1,10)
names(wts)<-unlist(strsplit(arg,","))[1:10]

#vector of recruitment rates to use
b0_vec<-seq(0,1,.01)

barrier<-function(nm,amt=0.5,len=10,pars){
  #create vector to modify for behavioral barrier
  wts<-rep(1,10)
  names(wts)<-c("q12","q23","q34","q45","q56","q21","q32","q43","q54","q65")
  #modify the wts of the movement probs in nm by the proportion amt
  wts[nm]<-amt
  
  #list for results
  ans<-list()
  
  #loop through recruitment rates
  for(i in 1:length(b0_vec)){
   temp<-rep(0,len)
   temp_null<-rep(0,len)
   tempeig<-rep(0,len)
   #loop through resamples
   for(j in 1:len){
     
    #sample distribution of movement parameters
    par<-sapply(pars,randFq)
    
    #create a post-barrier deterrent matrix (A_twiddle) and pre-barrier displacment matrix (A_null)
    A_twiddle<-do.call('Sigma',as.list(c(par*wts,b0_vec[i])))
    A_null<-do.call('Sigma',as.list(c(par,b0_vec[i])))
    
    #get eigensystem for A_twiddle
    eigs_twiddle<-eigen(A_twiddle)
    
    #get eigensystem for A_null
    eigs_null<-eigen(A_null)
    #get stable stage distributions
    ssd_twiddle<-eigs_twiddle$vectors[,1]/sum(eigs_twiddle$vectors[,1])
    ssd_null<-eigs_null$vectors[,1]/sum(eigs_null$vectors[,1])
    #calculate pct. change in metacommunity growth rate
    tempeig[j]<-eigs_twiddle$values[1]/eigs_null$values[1]-1
    #get the elements to calc pct change in Dresden Island per cap growth rate
    temp[j]<-A_twiddle[6,]%*%ssd_twiddle
    temp_null[j]<-A_null[6,]%*%ssd_null
   }
  ans[[i]]<-list(local=temp/temp_null-1,meta=tempeig)
  }
 return(ans)
}

#modify each upstream movement prob. resampling 500 times for Big Head Carp (bh)
q21_bh_ans<-barrier("q21",len = 500,pars=unlist(apply(bhPars, 1, list), recursive = FALSE))
q32_bh_ans<-barrier("q32",len = 500,pars=unlist(apply(bhPars, 1, list), recursive = FALSE))
q43_bh_ans<-barrier("q43",len = 500,pars=unlist(apply(bhPars, 1, list), recursive = FALSE))
q54_bh_ans<-barrier("q54",len = 500,pars=unlist(apply(bhPars, 1, list), recursive = FALSE))
q65_bh_ans<-barrier("q65",len = 500,pars=unlist(apply(bhPars, 1, list), recursive = FALSE))

#modify each upstream movement prob. resampling 500 times for Silver Carp (s)
q21_s_ans<-barrier("q21",len = 500,pars=unlist(apply(svPars, 1, list), recursive = FALSE))
q32_s_ans<-barrier("q32",len = 500,pars=unlist(apply(svPars, 1, list), recursive = FALSE))
q43_s_ans<-barrier("q43",len = 500,pars=unlist(apply(svPars, 1, list), recursive = FALSE))
q54_s_ans<-barrier("q54",len = 500,pars=unlist(apply(svPars, 1, list), recursive = FALSE))
q65_s_ans<-barrier("q65",len = 500,pars=unlist(apply(svPars, 1, list), recursive = FALSE))

#set color for plotting uncertainty polygons
cols<-rgb(t(col2rgb("grey44")),alpha = 128,maxColorValue = 255)


#plotter function for Dresden Island result
overlay_pltr<-function(data_name,data_name2){
    dat<-data.frame(sapply(get(data_name),'[',1))
    dat2<-data.frame(sapply(get(data_name2),'[',1))
    ylims<-c(-.1+min(apply(dat,2,mean)-1.96*apply(dat,2,sd)),.1+max(apply(dat,2,mean)+1.96*apply(dat,2,sd)))
    ylims2<-c(-.1+min(apply(dat2,2,mean)-1.96*apply(dat2,2,sd)),.1+max(apply(dat2,2,mean)+1.96*apply(dat2,2,sd)))
    ylims<-c(min(c(ylims[1],ylims2[1])),max(c(ylims[2],ylims2[2])))
    ylims[2]<-0
    std<-predict(smooth.spline(x = b0_vec, y=apply(dat,2,sd),spar = .6))$y
    mn<-predict(smooth.spline(x = b0_vec, y=apply(dat,2,mean),spar = .6))$y
    std2<-predict(smooth.spline(x = b0_vec, y=apply(dat2,2,sd),spar = .6))$y
    mn2<-predict(smooth.spline(x = b0_vec, y=apply(dat2,2,mean),spar = .6))$y
    plot(b0_vec,mn,ylim=ylims,type = "l",lwd=2,xlab="",ylab = "",bty="n",cex.axis=1.25)
    lines(b0_vec,mn2,lwd=2,lty=2,type = "l")
    polygon(c(b0_vec,rev(b0_vec)),
            c(mn-1.96*std,
              rev(mn+1.96*std)),col=cols,border = NA)
    polygon(c(b0_vec,rev(b0_vec)),
            c(mn2-1.96*std2,
              rev(mn2+1.96*std2)),col=cols,border = NA)
  }
#set up vectors to feed to plotter function
dat_nms<-c("q21_bh_ans","q32_bh_ans","q43_bh_ans","q54_bh_ans","q65_bh_ans")
dat_nms2<-c("q21_s_ans","q32_s_ans","q43_s_ans","q54_s_ans","q65_s_ans")

####create Figure 2####
pdf("Figure2.pdf")
par(mfrow=c(3,2),mar=c(1,2,2,1),oma=c(5,5,1,1))
for(i in 1:5){
  overlay_pltr(dat_nms[i],dat_nms2[i])
  mtext(bquote(phi[.(i+1)*','*.(i)]),3,padj = 1,cex = 1.25)
  mtext(paste0("(",letters[i],")"),adj=-.2,padj=.6)
}
plot(0:1,0:1,yaxt="n",xaxt='n',ylab='',xlab='',bty='n',type='n')
legend("center",legend = c("Bighead","Silver"),lty=c(1,2),lwd=2,bty = 'n',cex=1.5)
mtext(bquote('Recruitment Rate ('*b[0]*')'),side = 1,outer = TRUE,cex=1.5,padj =2 )
mtext(bquote('Pct. Change Growth Rate (%'*Delta*e^tilde(r[1])*')'),side = 2,outer = TRUE,cex=1.5,padj=-0.8 )
dev.off()

#create plotter function for metapopulation result
overlay_pltr_meta<-function(data_name,data_name2){
  dat<-data.frame(sapply(get(data_name),'[',2))
  dat2<-data.frame(sapply(get(data_name2),'[',2))
  ylims<-c(-.01+min(apply(dat,2,mean)-1.96*apply(dat,2,sd)),.01+max(apply(dat,2,mean)+1.96*apply(dat,2,sd)))
  ylims2<-c(-.01+min(apply(dat2,2,mean)-1.96*apply(dat2,2,sd)),.01+max(apply(dat2,2,mean)+1.96*apply(dat2,2,sd)))
  ylims<-c(min(c(ylims[1],ylims2[1])),max(c(ylims[2],ylims2[2])))
  ylims<-c(-.01,.015)
  std<-predict(smooth.spline(x = b0_vec, y=apply(dat,2,sd),spar = .6))$y
  mn<-predict(smooth.spline(x = b0_vec, y=apply(dat,2,mean),spar = .6))$y
  std2<-predict(smooth.spline(x = b0_vec, y=apply(dat2,2,sd),spar = .6))$y
  mn2<-predict(smooth.spline(x = b0_vec, y=apply(dat2,2,mean),spar = .6))$y
  plot(b0_vec,mn,ylim=ylims,type = "l",lwd=2,xlab="",ylab = "",bty="n",cex.axis=1.25)
  lines(b0_vec,mn2,lwd=2,lty=2,type = "l")
  polygon(c(b0_vec,rev(b0_vec)),
          c(mn-1.96*std,
            rev(mn+1.96*std)),col=cols,border = NA)
  polygon(c(b0_vec,rev(b0_vec)),
          c(mn2-1.96*std2,
            rev(mn2+1.96*std2)),col=cols,border = NA)
}

####create Figure 3####
pdf("Figure3.pdf")
par(mfrow=c(3,2),mar=c(1,3,2,1),oma=c(5,5,1,1))
#par(mfrow=c(3,2),pty = "s",mar=c(4,2,1,1))
for(i in 1:5){
  overlay_pltr_meta(dat_nms[i],dat_nms2[i])
  mtext(bquote(phi[.(i+1)*','*.(i)]),3,padj = 1,cex = 1.25)
  mtext(paste0("(",letters[i],")"),adj=-.2,padj=.6)
}
plot(0:1,0:1,yaxt="n",xaxt='n',ylab='',xlab='',bty='n',type='n')
legend("center",legend = c("Bighead","Silver"),lty=c(1,2),lwd=2,bty = 'n',cex=1.5)
mtext(bquote('Recruitment Rate ('*b[0]*')'),side = 1,outer = TRUE,cex=1.5,padj =2 )
mtext(bquote('Pct. Change Metapopulation Growth Rate (%'*Delta*tilde(lambda)[1]*')'),side = 2,outer = TRUE,cex=1.5,padj=-1 )
dev.off()

#### TRANSIENT DYNAMICS ####

#function to multiply matrix by itself many times 
mat_pwr<-function(x,power){
  fun<-paste(rep("x",power),collapse="%*%")
  return(eval(str2expression(fun)))
}

#create function similar in structure to one above but calculates elements of 
#transient dynamics
barrier_t_dyn<-function(nm,len=10,amt=0.5,pars){
  
  wts<-rep(1,10)
  names(wts)<-c("q12","q23","q34","q45","q56","q21","q32","q43","q54","q65")
  ans<-list()
  wts[nm]<-amt
  for(i in 1:length(b0_vec)){
    temp_s1<-rep(0,len)
    temp_life<-rep(0,len)
    
    for(j in 1:len){
      par<-sapply(pars,randFq)
      A_twiddle<-do.call('Sigma',as.list(c(par*wts,b0_vec[i])))
      A_null<-do.call('Sigma',as.list(c(par,b0_vec[i])))
      eigs_twiddle<-eigen(A_twiddle)
      eigs_null<-eigen(A_null)
      #get pre-deterrent barrier stable spacial distribution
      n0_hat<-eigs_null$vectors[,1]/sum(eigs_null$vectors[,1])
      
      #get maximum of first 30 steps of Dresden Island population after post-barrier
      max_t_dres<-max(sapply(1:30,function(t)mat_pwr(A_twiddle,t)%*%n0_hat/norm(mat_pwr(A_twiddle,t)%*%n0_hat))[6,])/n0_hat[6]-1
      
      #get estimate of timescale of trans. dynamics
      avg_life<-1/log(eigs_twiddle$values[1]/eigs_twiddle$values[2])
      
      temp_life[j]<-avg_life
      temp_s1[j]<-max_t_dres
     
    }
    ans[[i]]<-list(avg_life=temp_life,s1=temp_s1)
  }
  return(ans)
}

#modify each upstream movement prob. resampling 500 times for Big Head Carp (bh) 
q21_bh_ans_t<-barrier_t_dyn("q21",len = 500,pars=unlist(apply(bhPars, 1, list), recursive = FALSE))
q32_bh_ans_t<-barrier_t_dyn("q32",len = 500,pars=unlist(apply(bhPars, 1, list), recursive = FALSE))
q43_bh_ans_t<-barrier_t_dyn("q43",len = 500,pars=unlist(apply(bhPars, 1, list), recursive = FALSE))
q54_bh_ans_t<-barrier_t_dyn("q54",len = 500,pars=unlist(apply(bhPars, 1, list), recursive = FALSE))
q65_bh_ans_t<-barrier_t_dyn("q65",len = 500,pars=unlist(apply(bhPars, 1, list), recursive = FALSE))

#modify each upstream movement prob. resampling 500 times for Silver Carp (s) 
q21_s_ans_t<-barrier_t_dyn("q21",len = 500,pars=unlist(apply(svPars, 1, list), recursive = FALSE))
q32_s_ans_t<-barrier_t_dyn("q32",len = 500,pars=unlist(apply(svPars, 1, list), recursive = FALSE))
q43_s_ans_t<-barrier_t_dyn("q43",len = 500,pars=unlist(apply(svPars, 1, list), recursive = FALSE))
q54_s_ans_t<-barrier_t_dyn("q54",len = 500,pars=unlist(apply(svPars, 1, list), recursive = FALSE))
q65_s_ans_t<-barrier_t_dyn("q65",len = 500,pars=unlist(apply(svPars, 1, list), recursive = FALSE))


cols<-rgb(t(col2rgb("grey44")),alpha = 128,maxColorValue = 255)
par_names<-c("q12","q23","q34","q45","q56","q21","q32","q43","q54","q65")
dat<-data.frame(sapply(get(dat_nms[1]),'[',1))
dat2<-data.frame(sapply(get(dat_nms2[1]),'[',1))

#create plotter function for transient timescale
pltr_trans_meta<-function(data_name,data_name2){
  dat<-data.frame(sapply(get(data_name),'[',1))
  dat2<-data.frame(sapply(get(data_name2),'[',1))
  ylims<-c(-.1+min(apply(dat,2,mean)-1.96*apply(dat,2,sd)),.1+max(apply(dat,2,mean)+1.96*apply(dat,2,sd)))
  ylims2<-c(-.1+min(apply(dat2,2,mean)-1.96*apply(dat2,2,sd)),.1+max(apply(dat2,2,mean)+1.96*apply(dat2,2,sd)))
  ylims<-c(min(c(ylims[1],ylims2[1])),max(c(ylims[2],ylims2[2])))
  ylims[2]<-80
  std<-predict(smooth.spline(x = b0_vec, y=apply(dat,2,sd),spar = .6))$y
  mn<-predict(smooth.spline(x = b0_vec, y=apply(dat,2,mean),spar = .6))$y
  std2<-predict(smooth.spline(x = b0_vec, y=apply(dat2,2,sd),spar = .6))$y
  mn2<-predict(smooth.spline(x = b0_vec, y=apply(dat2,2,mean),spar = .6))$y
  plot(b0_vec,mn,ylim=ylims,type = "l",lwd=2,xlab="",ylab = "",bty="n",cex.axis=1.25)
  lines(b0_vec,mn2,lwd=2,lty=2,type = "l")
  polygon(c(b0_vec,rev(b0_vec)),
          c(mn-1.96*std,
            rev(mn+1.96*std)),col=cols,border = NA)
  polygon(c(b0_vec,rev(b0_vec)),
          c(mn2-1.96*std2,
            rev(mn2+1.96*std2)),col=cols,border = NA)
}

#set up data to feed to plotter functions
dat_nms<-c("q21_bh_ans_t","q32_bh_ans_t","q43_bh_ans_t","q54_bh_ans_t","q65_bh_ans_t")
dat_nms2<-c("q21_s_ans_t","q32_s_ans_t","q43_s_ans_t","q54_s_ans_t","q65_s_ans_t")

####Figure 4####
pdf("Figure4.pdf")
par(mfrow=c(3,2),mar=c(1,2,2,1),oma=c(5,5,1,1))
#par(mfrow=c(3,2),pty = "s",mar=c(4,2,1,1))
for(i in 1:5){
  pltr_trans_meta(dat_nms[i],dat_nms2[i])
  mtext(bquote(phi[.(i+1)*','*.(i)]),3,padj = 2,cex = 1.25)
  mtext(paste0("(",letters[i],")"),adj=-.2,padj=.6)
  
}
plot(0:1,0:1,yaxt="n",xaxt='n',ylab='',xlab='',bty='n',type='n')
legend("center",legend = c("Bighead","Silver"),lty=c(1,2),lwd=2,bty = 'n',cex=1.5)
mtext(bquote('Recruitment Rate ('*b[0]*')'),side = 1,outer = TRUE,cex=1.5,padj =2 )
mtext("Transient Timescale (Years)",side = 2,outer = TRUE,cex=1.5,padj=-2 )
dev.off()

#create plotter function for Dresden Island transient dynamics
pltr_trans_dres<-function(data_name,data_name2){
  dat<-data.frame(sapply(get(data_name),'[',2))
  dat2<-data.frame(sapply(get(data_name2),'[',2))
  ylims<-c(-.01+min(apply(dat,2,mean)-1.96*apply(dat,2,sd)),.01+max(apply(dat,2,mean)+1.96*apply(dat,2,sd)))
  ylims2<-c(-.01+min(apply(dat2,2,mean)-1.96*apply(dat2,2,sd)),.01+max(apply(dat2,2,mean)+1.96*apply(dat2,2,sd)))
  ylims<-c(min(c(ylims[1],ylims2[1])),max(c(ylims[2],ylims2[2])))
  ylims<-c(-.2,0)
  std<-predict(smooth.spline(x = b0_vec, y=apply(dat,2,sd),spar = .6))$y
  mn<-predict(smooth.spline(x = b0_vec, y=apply(dat,2,mean),spar = .6))$y
  std2<-predict(smooth.spline(x = b0_vec, y=apply(dat2,2,sd),spar = .6))$y
  mn2<-predict(smooth.spline(x = b0_vec, y=apply(dat2,2,mean),spar = .6))$y
  plot(b0_vec,mn,ylim=ylims,type = "l",lwd=2,xlab="",ylab = "",bty="n",cex.axis=1.25)
  lines(b0_vec,mn2,lwd=2,lty=2,type = "l")
  polygon(c(b0_vec,rev(b0_vec)),
          c(mn-1.96*std,
            rev(mn+1.96*std)),col=cols,border = NA)
  polygon(c(b0_vec,rev(b0_vec)),
          c(mn2-1.96*std2,
            rev(mn2+1.96*std2)),col=cols,border = NA)
  abline(h=1,lty=3)
}

####Figure 5####
pdf("Figure5.pdf")
par(mfrow=c(3,2),mar=c(1,2,3,1),oma=c(5,5,1,1))
for(i in 1:5){
  pltr_trans_dres(dat_nms[i],dat_nms2[i])
  mtext(bquote(phi[.(i+1)*','*.(i)]),3,padj = .5,cex = 1.25)
  mtext(paste0("(",letters[i],")"),adj=-.2,padj=.6)
  
}
plot(0:1,0:1,yaxt="n",xaxt='n',ylab='',xlab='',bty='n',type='n')
legend("center",legend = c("Bighead","Silver"),lty=c(1,2),lwd=2,bty = 'n',cex=1.5)
mtext(bquote('Recruitment Rate ('*b[0]*')'),side = 1,outer = TRUE,cex=1.5,padj =2 )
mtext(bquote('Pct. Change Max. Pop. Dresden Island (%'*Delta*'max('*hat(n)[1]*')'[t]*")"),side = 2,outer = TRUE,cex=1.5,padj=-1 )
dev.off()

####create Figure 6####
#calculate averages percent change in per-capita growth rate of Dresden Island population
dat_nms<-c("q21_bh_ans","q32_bh_ans","q43_bh_ans","q54_bh_ans","q65_bh_ans")
dat_nms2<-c("q21_s_ans","q32_s_ans","q43_s_ans","q54_s_ans","q65_s_ans")

bh_avg<-NULL
sv_avg<-NULL

for(i in 1:5){
  dat<-data.frame(sapply(get(dat_nms[i]),'[',1))
  dat2<-data.frame(sapply(get(dat_nms2[i]),'[',1))
  bh_avg[i]<-mean(apply(dat,2,mean))
  sv_avg[i]<-mean(apply(dat2,2,mean))
}


pdf("Figure6.pdf")
par(oma=c(0,1,0,0))
plot(c(Ebh[6:10],Esv[6:10]),c(bh_avg,sv_avg),pch=21,bg=rep(c("black","grey45"),each=5),bty="n",
     xlab="Movement Probability",ylab="",cex.lab=1.25,cex.axis=1.25,cex=1.25)
mtext(bquote('Pct. Change Growth Rate (%'*Delta*e^tilde(r)[1]*')'),2,padj=-2,cex=1.25)
for(i in 1:5)text(Ebh[6:10][i],bh_avg[i],bquote(phi[.(i+1)*','*.(i)]),pos=1,col = "black")
for(i in 1:5)text(Esv[6:10][i],sv_avg[i],bquote(phi[.(i+1)*','*.(i)]),pos = 4,col="grey45")
legend("topleft",legend = c("Bighead","Silver"),pch=21,pt.bg =c("black","grey45"),bty = "n")
dev.off()


