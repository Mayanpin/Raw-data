
if (T) {
  dir.create("00_rawdatas/TARGET")
  dir.create("00_rawdatas/GEO",recursive = T)
  dir.create("00_pre_datas/TARGET")
  dir.create("00_pre_datas/GEO",recursive = T)
  dir.create("reference")
  dir.create("files")
  dir.create("Figures")
  
}
library(openxlsx)
library(tidyselect, lib.loc = "/home/R36library251mount/library")
library("hgu133plus2.db")
library(WGCNA, lib.loc = "/home/R36library251mount/library")
library(glmnet)
library(reshape2)
library(ggpubr)
library(ggsci)
library(maftools)
library(tidyr)
library(pheatmap)
library(clusterProfiler)
library(ComplexHeatmap)
library(GSVA)
library(GSEABase)
library(corrplot)
library(survival)
library(survminer)
library(colorspace)
options(stringsAsFactors = F)

custom_theme <- function() {
  theme_survminer() %+replace%
    theme(text = element_text(family = 'Times'),panel.grid = element_blank())
}
my_violin=function(dat,group,group_cols=ggsci::pal_aaas()(10),#test_method='kruskal.test',
                   fill= "Group",label=c("p.format",'p.signif')[1],
                   xlab='',ylab='',title='',x.size=10,y.size=10,legend.position='top'){
  
  # dat = tcga.b.cell$GeneSet,
  # group = tcga.subtype$Cluster
  data=data.frame(Group=group,value=dat)
  data=crbind2DataFrame(data)
  data=melt(data)
  data=data[which(!is.na(data[,1])),]
  if(length(names(table(group)))>2){
    test_method='kruskal.test'
  }else{
    test_method='wilcox.test'
  }
  p=ggplot(data,aes(x=Group, y=value,fill=Group)) +
    geom_violin()+  
    scale_fill_manual(values = group_cols)+
    geom_boxplot(width=0.2,position=position_dodge(0.8),outlier.colour = NA,fill="white")+
    theme_classic(base_size = 20)+labs(x=xlab,y=ylab,title = title)+
    ggpubr::stat_compare_means(aes(group=Group), label = label, method =test_method)+
    theme(legend.position = legend.position,axis.text = element_text(color = 'black'),
          title = element_text(size = 12),
          axis.text.x = element_text(size = 12),axis.title.x = element_text(size = 15),
          axis.text.y = element_text(size = 12),axis.title.y = element_text(size = 12))
  return(p)
}

my_mutiboxplot=function(dat,group,group_cols=ggsci::pal_aaas()(10),
                        #test_method=c('t.test','wilcox.test','anova','kruskal.test')[4],
                        bw=T,xlab='',ylab='score',title='',size=10,angle = 45, hjust = 1,
                        legend.position='top',fill='group',notch=F){
  # dat=tcga.est[tcga.subtype.cli$Samples,]
  # group=tcga.subtype.cli$Cluster
  dat.bind=cbind(dat,Cluster=group)
  dat.bind=crbind2DataFrame(dat.bind)
  dat.melt=melt(dat.bind)
  #data=data[which(!is.na(data[,1])),]
  colnames(dat.melt)=c('Group','type','value')
  if(length(names(table(group)))>2){
    test_method='kruskal.test'
  }else{
    test_method=c('wilcox.test','t.test')[1]
  }
  p=dat.melt %>%
    ggplot(aes(x=type, y=value,fill=Group)) +
    geom_boxplot(notch = notch) +  
    scale_fill_manual(values =group_cols)+   #
    ggpubr::stat_compare_means(aes(group=Group), label = "p.signif", method = test_method)+
    labs(x="", y = ylab, fill =fill,title =title) +
    #theme_light()+
    theme_bw()+
    #theme_classic()
    theme(legend.position = legend.position,                 #
          plot.title = element_text(hjust = 0.5),
          text = element_text(family = 'Times',size = size),
          axis.text.x = element_text(size = size,angle = angle, hjust = hjust)) # 
  return(p)
}

sig_boxplot_t<-function(dat,leg,ylab,xlab='',palette=ggsci::pal_lancet()(10)[3:4]){
  dat=na.omit(dat)
  colnames(dat)=c('group','gene')
  dat=dat[order(dat$group),]
  all.combn=combn(as.character(unique(dat$group)),2)
  my_comparisons=lapply(seq_len(ncol(all.combn)), function(i) all.combn[,i])
  pp=ggboxplot(dat, 
               x='group', y='gene', color = 'group',
               palette =  palette,
               short.panel.labs = T,outlier.shape = NA)+
    stat_compare_means(comparisons=my_comparisons,method="t.test",label = "p.signif")+
    ylab(ylab)+xlab(xlab)+labs(color=leg)
  return(pp)
}
sig_boxplot_w<-function(dat,leg,ylab,xlab='',palette=ggsci::pal_lancet()(10)[3:4]){
  dat=na.omit(dat)
  colnames(dat)=c('group','gene')
  dat=dat[order(dat$group),]
  all.combn=combn(as.character(unique(dat$group)),2)
  my_comparisons=lapply(seq_len(ncol(all.combn)), function(i) all.combn[,i])
  pp=ggboxplot(dat, 
               x='group', y='gene', color = 'group',
               palette =  palette,
               short.panel.labs = T,outlier.shape = NA)+
    stat_compare_means(comparisons=my_comparisons,method="wilcox.test",label = "p.signif")+
    ylab(ylab)+xlab(xlab)+labs(color=leg)
  return(pp)
}
mg_violin_1=function(data,xangle=0,ylab='value',xlab='',leg.title='Group',test_method='anova',cmp_test_method='t.test',
                     legend.pos='r',melt=F,jitter=T,ylim=NULL,show_compare=NULL,point_size=NULL,group_col){
  library(ggplot2)
  if(is.null(ylim)){
    
  }
  if(melt){
    data_m=data
    colnames(data_m)=c('Group','value')
  }else{
    data_m=reshape2::melt(data)
    colnames(data_m)=c('Group','value')
  }
  if(!is.null(ylim)){
    data_m$value[data_m$value>ylim[2]]<-NA
  }
  data_m=data_m[which(!is.na(data_m[,1])),]
  if(xangle==0){
    tx=element_text(colour="black",family="Times")
  }else{
    tx=element_text(angle=xangle,hjust = 1,colour="black",family="Times")
  }
  
  pos='right'
  if(is.null(legend.pos)){
    pos='none'
  }else if(legend.pos=='tr'){
    pos=c(1,1)
  }else if(legend.pos=='br'){
    pos=c(1,0)
  }else if(legend.pos=='tl'){
    pos=c(0,1)
  }else if(legend.pos=='bl'){
    pos=c(0,0)
  }else if(legend.pos=='t'){
    pos='top'
  }else if(legend.pos=='r'){
    pos='right'
  }else if(legend.pos=='b'){
    pos='bottom'
  }
  uni.group=unique(data_m[,1])
  ct=length(uni.group)
  
  p1<-ggplot(data_m,aes(x=Group,y=value))+geom_violin(alpha=0.7)
  # if(ct<=4){
  #   p1=p1+ggsci::scale_fill_lancet()
  # }else if(ct<=10){
  #   p1=p1+ggsci::scale_fill_npg(name=leg.title)
  # }else if(ct<=20){
  #   p1=p1+ggsci::scale_fill_d3(palette = "category20",name=leg.title)
  # }else if(ct<=30){
  #   cbPalette=c(ggsci::pal_npg("nrc", alpha = 0.6)(10),ggsci::pal_d3("category20", alpha = 0.6)(20))
  #   p1=p1+scale_fill_manual(values=cbPalette[1:ct])
  # }else if(ct<=38){
  #   cbPalette=c(ggsci::pal_lancet()(10)
  #               ,ggsci::pal_npg("nrc", alpha = 0.6)(10)
  #               ,ggsci::pal_d3("category20", alpha = 0.6)(20)
  #               ,ggsci::pal_nejm("default", alpha = 0.6)(8))
  #   p1=p1+scale_fill_manual(values=cbPalette[1:ct])
  # }
  p1=p1+scale_fill_manual(values=group_col)
  if(jitter){
    if(is.null(point_size)){
      p1<-p1+geom_jitter(alpha=0.3,col='black',show.legend=FALSE,width = 0.2)
    }else{
      p1<-p1+geom_jitter(alpha=0.3,col='black',show.legend=FALSE,width = 0.2,size=point_size)
    }
  }
  
  p1=p1+theme_bw()+geom_boxplot(width=0.2,aes(fill=Group),outlier.shape = NA)
  p1=p1+theme(axis.text.x=tx, 
              axis.text.y=element_text(family="Times",face="plain"), 
              axis.title.y=element_text(family="Times",face="plain"), 
              legend.text=element_text(face="plain", family="Times", colour="black" 
              ),
              legend.title=element_text(face="plain", family="Times", colour="black"
              ),
              legend.justification=pos, legend.position=pos
              ,legend.background = element_rect(fill = NA, colour = NA)
  )+ylab(ylab)+xlab(xlab)
  til=''
  if(test_method=='anova'){
    if(length(unique(data_m[,1]))<3){
      x1=data_m[,2][which(data_m[,1]==unique(data_m[,1])[1])]
      x2=data_m[,2][which(data_m[,1]==unique(data_m[,1])[2])]
      pv=t.test(x1,x2)$p.value 
      til=paste0('t-tests p=',signif(pv,2))
    }else{
      fit <- aov(value~Group, data = data_m)
      pv=summary(fit)[[1]][5][[1]]
      fv=summary(fit)[[1]][4][[1]]
      til=paste0('ANOVA tests p=',signif(pv,2))
    }
  }else{
    if(length(unique(data_m[,1]))<3){
      x1=data_m[,2][which(data_m[,1]==unique(data_m[,1])[1])]
      x2=data_m[,2][which(data_m[,1]==unique(data_m[,1])[2])]
      pv=wilcox.test(x1,x2)$p.value 
      til=paste0('wilcox.tests p=',signif(pv,2))
    }else{
      fit=kruskal.test(value~Group, data = data_m)
      pv=fit$p.value
      til=paste0('Kruskal-Wallis test p=',signif(pv,2))
    }
  }
  p1=p1+ggtitle(til) 
  if(!is.null(ylim)){
    p1=p1+ylim(ylim)
  }
  if(is.null(show_compare)){
    if(length(uni.group)>5){
      show_compare=F
    }else{
      show_compare=T
    }
  }
  if(show_compare){
    comps=list()
    for(i in 1:(length(uni.group)-1)){
      for(j in (i+1):length(uni.group)){
        comps=c(comps,list(c(uni.group[i],uni.group[j])))
      }
    }
    p1=p1+ggpubr::stat_compare_means(comparisons = comps,method = cmp_test_method,label= "p.signif")
  }
  return(p1)
}
coxFun <- function(dat){
  library(survival)
  colnames(dat)=c('time','status','gene')
  fmla=as.formula("Surv(time,status)~gene")
  cox=coxph(fmla,data=dat)
  p=summary(cox)[[7]][5]
  result=c(p,summary(cox)[[8]][1],summary(cox)[[8]][3],summary(cox)[[8]][4])
  return(result)
}
get_riskscore.lasso<-function(dat,os,os.time,labels=c('A','B')){
  library(glmnet)
  set.seed(2021)
  fit1=glmnet(as.matrix(dat)
              ,cbind(time=os.time,
                     status=os)
              ,family="cox"
              ,nlambda=100
              , alpha=1) 
  
  cv.fit<-cv.glmnet(as.matrix(dat)
                    ,cbind(time=os.time,
                           status=os)
                    ,family="cox"
                    ,nfolds = 10
                    ,nlambda=100
                    , alpha=1)
  sig.coef <- coefficients(cv.fit,s=cv.fit$lambda.min)[which(coefficients(cv.fit,s=cv.fit$lambda.min)[,1]!=0),1]
  #print(cv.fit$lambda.min)
  #length(names(sig.coef))
  #10
  mg_plot_lasso <- function(fit,cv_fit,lambda=NULL,show_text=T,figLabels=c('A','B')){
    if(is.null(lambda)){
      lmda=cv_fit$lambda.min
    }else{
      lmda=lambda
    }
    fit.coef=fit$beta[(apply(fit$beta,1,function(x){
      return(sum(x!=0))
    })>0),]
    
    fit.coef=as.matrix(fit.coef)
    colnames(fit.coef)=fit$lambda
    #fit$lambda==cv_fit$lambda
    library(ggplot2)
    dat=data.table::melt(t(as.matrix(fit.coef)))
    dat_z=dat[which(dat$value==0),]
    dat=dat[which(dat$value!=0),]
    dat.sv=rbind()
    for (u in unique(dat_z[,2])) {
      t.z=dat_z[which(dat_z[,2]==u),1]
      t.zx=max(t.z)
      dat.sv=rbind(dat.sv,c(t.zx,u,0))
      t.zn=min(t.z)
      if(t.zx!=t.zn){
        dat.sv=rbind(dat.sv,c(t.zn,u,0))
      }
    }
    colnames(dat.sv)=colnames(dat_z)
    #dat_z=dat_z[dat_z[,2]%in%names(which(fit.coef[,which(fit$lambda==lmda)]!=0)),]
    dat=crbind2DataFrame(rbind(dat,dat.sv))
    mn=min(-log(dat$Var1))
    mx=max(-log(dat$Var1))
    if(show_text){
      mx=(mx-mn)*0.1+mx
    }
    p=ggplot(dat, aes(x=-log(Var1), y=value,colour=Var2))+geom_line()+theme_bw()+theme(legend.position = "none")
    p=p+coord_cartesian(xlim=c(mn, mx))+xlab('-ln(lambda)')+ylab('Coefficients')
    if(show_text){
      fl=fit.coef[which(fit.coef[,which(fit$lambda==lmda)]!=0),ncol(fit.coef)]
      for_label=data.frame(Var1=rep(min(dat$Var1),length(fl)),Var2=names(fl),value=fl)
      p=p+ggrepel::geom_label_repel(
        aes(label = Var2,color=Var2),
        data = for_label,hjust = 0
      )
    }
    p=p+geom_vline(aes(xintercept=-log(lmda)), colour="#BB0000", linetype="dashed")
    p=p+annotate('text',x=-log(lmda),y=min(dat[,3]),label=paste0('lambda=',round(lmda,4)))
    tgc=data.frame(lambda=cv_fit$lambda,cvm=cv_fit$cvm,cvup=cv_fit$cvup,cvlo=cv_fit$cvlo,cvsd=cv_fit$cvsd
                   ,col=ifelse(cv_fit$lambda>=cv_fit$lambda.min&cv_fit$lambda<=cv_fit$lambda.1se,ifelse(cv_fit$lambda==lmda,'A','C'),'B'))
    p1=ggplot(tgc, aes(x=log(lambda), y=cvm)) + xlab('ln(lambda)')+ ylab('Parial Likelihood Deviance')+
      geom_errorbar(aes(ymin=cvm-cvsd, ymax=cvm+cvsd)) +
      geom_point(aes(colour=col))
    p1=p1+theme_bw()+theme(legend.position = "none")
    gal=ggpubr::ggarrange(p,p1, ncol = 2, nrow = 1
                          #,align = "hv"
                          ,labels = figLabels)
    return(gal)
  }
  lasso.pdf <- mg_plot_lasso(fit1,
                             cv.fit,
                             show_text=T,
                             figLabels=labels)
  return(list(lasso.gene=names(sig.coef),lambda.min=cv.fit$lambda.min,plot=lasso.pdf))
}
get_riskscore<-function(dat,os,os.time,step=T,direction=c("both", "backward", "forward")[1]){
  tcga_dat1 <- cbind(time=os.time,
                     status=os,
                     dat)
  tcga_dat1=crbind2DataFrame(tcga_dat1)
  colnames(tcga_dat1)=gsub('-','__',colnames(tcga_dat1))
  gene111=gsub('-','__',colnames(dat))
  fmla <- as.formula(paste0("Surv(time, status) ~"
                            ,paste0(gene111,collapse = '+')))
  cox <- coxph(fmla, data =as.data.frame(tcga_dat1))
  if(step==T){
    cox1 <- step(cox,direction =direction)
  }else{
    cox1=cox
  }
  lan <- coef(cox1)
  #round(lan, 3)
  genes <- names(cox1$coefficients)
  mult_results=paste0(round(lan, 3), '*', names(lan),collapse = '+')
  risk.tcga <- as.numeric(lan%*%as.matrix(t(tcga_dat1[,genes])))
  
  data_gene_score_final<-tcga_dat1
  data_gene_score_final$Samples<-rownames(data_gene_score_final)
  data_gene_score_final$riskscore=risk.tcga
  data_gene_score_final$riskscorez=mosaic::zscore(risk.tcga)
  #
  optimalCutoff <- survminer::surv_cutpoint(data.frame(time=data_gene_score_final$time/365,
                                                       event=data_gene_score_final$status,
                                                       risk=data_gene_score_final$riskscore), 
                                            time = "time", event = "event",variables = c("risk"))
  optimalCutoff=optimalCutoff$cutpoint$cutpoint[1]
  #print(optimalCutoff)
  #optimalCutoff=median(data_gene_score_final$riskscore)
  #optimalCutoff=0
  data_gene_score_final$Risk=ifelse(data_gene_score_final$riskscore>optimalCutoff,'High','Low')
  table(data_gene_score_final$Risk)
  data_gene_score_final$cutoff=optimalCutoff
  return(list(result=data_gene_score_final,module.gene=genes,model=mult_results))
}
bioForest=function(rt=null,col){
  #
  gene <- rownames(rt)
  hr <- sprintf("%.3f",rt[,2])
  hrLow  <- sprintf("%.3f",rt[,3])
  hrHigh <- sprintf("%.3f",rt[,4])
  Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
  pVal <- ifelse(rt[,1]<0.001, "<0.001", sprintf("%.3f", rt[,1]))
  
  #
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(3,2.5))
  
  #
  xlim = c(0,3)
  par(mar=c(4,2.5,2,1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
  text(3.1,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3.1,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1)
  
  #
  par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
  xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
  abline(v=1,col="black",lty=2,lwd=2)
  boxcolor = ifelse(as.numeric(hr) > 1, col[2], col[1])
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.5)
  axis(1)
}
#
pie_compare_plot<-function(dat,gname,group_cols){
  library(dplyr)
  library(ggplot2)
  g_n<-as.character(names(table(dat[,gname])))
  vname <- setdiff(colnames(dat), gname)
  pie.gname=data.frame()
  fisher.p <- c()
  for (i in vname) {
    tmp <- table(dat[,gname], dat[,i])
    p <- signif(chisq.test(tmp)$p.value,digits=4)
    names(p) <- i
    fisher.p <- c(fisher.p, p)
    pie.dat <- 
      tmp %>% as.data.frame() %>% group_by(Var1) %>% mutate(Pct = Freq/sum(Freq)) %>% as.data.frame()
    pie.gname=rbind.data.frame(pie.gname,pie.dat)
  }
  #plot--
  vname_col<-list()
  for(i in 1:length(vname)){
    col_i<-length(as.character(unique(na.omit(dat[,vname[i]]))))
    vname_col[[i]] <- alpha(group_cols[i], 
                            sample(x = 1:2,size = col_i,replace = F )/2)
  }
  
  #names(vname_col)=vname
  col_num<-ncol(dat)
  row_num<-length(g_n)
  row_nums=1+2*row_num+2
  #c(1:col_num,2*row_num*col_num+1)
  #第一行
  nums1=c()
  for (i in 1:col_num){
    nums1=c(nums1,rep(i,3))
  }
  #
  nums2=c()
  for (j in 1:row_num){
    nums21=c()
    for (i in (col_num*j+1):((1+j)*col_num)){
      nums21=c(nums21,rep(i,3))
    }
    nums2=c(nums2,nums21,nums21)
  }
  
  #
  nums3=c()
  #(1+row_num)*col_num+1,(2+row_num)*col_num
  for (i in (((1+row_num)*col_num+1):((2+row_num)*col_num))){
    nums3=c(nums3,rep(i,3))
  }
  #
  nums4=c(rep(((2+row_num)*col_num)+1,col_num*3))
  nums=c(nums1,nums2,nums3,nums4)
  showLayout <- F 
  
  layout(matrix(nums,byrow = T,nrow = row_nums))
  if(showLayout) {
    layout.show(n = ((2+row_num)*col_num)+1) #
  }
  par(bty="n", mgp = c(0,0,0), mar = c(0,0,0,0), lwd = 2) 
  plot(1,1,
       xlab = "",xaxt = "n", # 
       ylab = "",yaxt = "n") # 
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4]
       ,col = "white"
  ) # 
  text((par("usr")[1]+par("usr")[2])/2, # 
       (par("usr")[3]+par("usr")[4])/2,
       gname,cex = 2, col = "black") #
  #
  for (i in 1:length(vname)){
    plot(1,1,
         xlab = "",xaxt = "n", 
         ylab = "",yaxt = "n") 
    rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4]
         ,col = "white"
    ) 
    text((par("usr")[1]+par("usr")[2])/2, 
         (par("usr")[3]+par("usr")[4])/2,
         vname[i],cex = 2, col = "black") 
  }
  #
  for (i in 1:length(g_n)){
    plot(1,1,
         xlab = "",xaxt = "n", # 
         ylab = "",yaxt = "n") # 
    rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4]
         ,col = "white"
    ) # 
    text((par("usr")[1]+par("usr")[2])/2,
         (par("usr")[3]+par("usr")[4])/2,
         paste0(g_n[i],"\n(n = ",as.numeric(table(dat[,gname])[i]),")")
         ,cex = 2
         , col = "black") 
    for (j in 1:length(vname)){
      aa=as.character(unique(dat[,vname[j]]))
      pie.gname1=pie.gname[pie.gname$Var1==g_n[i],]
      pie.gname1=pie.gname[pie.gname$Var1==g_n[i] & pie.gname$Var2 %in% aa,]
      pie(pie.gname1$Pct, 
          col = vname_col[[j]], 
          border = "white",  
          radius = 1, 
          labels = NA,
          init.angle = 90)
      symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)
      if(j==length(vname)){
        abline(v = par("usr")[2], col = "black")
      }
    }
    
  }
  plot(1,1,
       xlab = "",xaxt = "n",
       ylab = "",yaxt = "n")
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4]
       ,col = "white"
  ) 
  for(i in vname){
    plot(1,1,col = "white",
         xlab = "",xaxt = "n", # 
         ylab = "",yaxt = "n") # 
    text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")
         (par("usr")[3]+par("usr")[4])/2,
         paste0("p = ",fisher.p[i]),cex = 1.5, col = "black") # 
    abline(h = par("usr")[3], col = "black")
  }
  abline(v = par("usr")[2], col = "black")
  plot(0,0,col = "white",
       xlab = "",xaxt = "n", # 
       ylab = "",yaxt = "n") # 
  leg_nam=c()
  for (ii in 1:length(vname)){
    aa=as.character(unique(dat[,vname[ii]]))
    pie.gname1=pie.gname[pie.gname$Var2 %in% aa,]
    leg_nam=c(leg_nam,as.character(unique(pie.gname1$Var2)))
  }
  vname_col_name=unlist(vname_col)
  
  legend(x='center',
         legend = leg_nam,
         fill = vname_col_name,
         border = NA, 
         bty = "n", 
         cex = 1.2,
         x.intersp = 0.05,
         y.intersp = 1,
         text.width = 0.075, 
         horiz = T)
}
getGeneFC=function(gene.exp,group,ulab=ulab,dlab=dlab){
  degs_C1_C3=mg_limma_DEG(gene.exp, 
                          group,
                          ulab=ulab,
                          dlab=dlab)
  
  ## 
  degs_C1_C3_sig<-degs_C1_C3$DEG[which(degs_C1_C3$DEG$adj.P.Val <= 1),]
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(stringr)
  degs_C1_C3_sig_gene<-rownames(degs_C1_C3_sig)
  
  degs_gene_entz=bitr(degs_C1_C3_sig_gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
  degs_gene_entz <- dplyr::distinct(degs_gene_entz,SYMBOL,.keep_all=TRUE)
  
  gene_df <- data.frame(logFC=degs_C1_C3_sig$logFC,
                        SYMBOL = rownames(degs_C1_C3_sig))
  gene_df <- merge(gene_df,degs_gene_entz,by="SYMBOL")
  head(gene_df)
  
  geneList<-gene_df$logFC
  names(geneList)=gene_df$ENTREZID 
  head(geneList)
  
  geneList=sort(geneList,decreasing = T)
  head(geneList) ## 
  return(geneList)
}


####00.#######################

#TARGET####
os_data<-read.delim('00_rawdatas/TARGET/Merge_RNA_seq_FPKM_TARGET-OS.txt',sep='\t',header = T,row.names = 1,check.names = F)
head(rownames(os_data))
range(os_data)
###fpkm tpm
# os_data=mg_FPKM2TPMs(os_data)
os_data=log2(os_data+1)
range(os_data)
#
target_exp <- exp_ensg2symbol(os_data)
#
genecode=read.delim('data/GeneTag.genecode.v32.txt',sep='\t',header = T)
table(genecode$TYPE)
mrna_genecode=genecode[which(genecode$TYPE=='protein_coding'),]
head(mrna_genecode)
target_exp=target_exp[rownames(target_exp)%in%mrna_genecode$SYMBOL,]
dim(target_exp)
#18448    88
#
target_cli=read.delim('00_rawdatas/TARGET/TARGET_OS_ClinicalData_Discovery_20210520.txt',sep='\t',header = T,stringsAsFactors = F)
colnames(target_cli)
cli.chose=c("TARGET.USI","Gender","Age.at.Diagnosis.in.Days","Disease.at.diagnosis",
            "Overall.Survival.Time.in.Days", "Vital.Status")
target_cli=target_cli[,cli.chose]
colnames(target_cli)=c('Samples','Gender','Age','Metaststic','OS.time','OS')
target_cli$Samples=paste0(target_cli$Samples,'-01')
rownames(target_cli)=target_cli$Samples
range(target_cli$Age)
target_cli$Age=target_cli$Age/365
fivenum(target_cli$Age)
target_cli$Age1=ifelse(target_cli$Age>15,'>15','<=15')
table(target_cli$Metaststic)
target_cli$Metaststic[target_cli$Metaststic=='Metastatic (confirmed)'|target_cli$Metaststic=='Metastatic']<-'YES'
target_cli$Metaststic[target_cli$Metaststic=='Non-metastatic (confirmed)' | target_cli$Metaststic=='Non-metastatic (Confirmed)']<-'NO'
target_cli$OS.time/365
target_cli <- target_cli %>% drop_na(OS.time)
#& target_cli$OS.time<10*365
dim(target_cli)
table(is.na(target_cli$OS.time))
###
target_cli=target_cli[target_cli$OS.time>0,]
table(target_cli$OS)
target_cli$OS=ifelse(target_cli$OS=='Alive',0,1)
com_sam=intersect(rownames(target_cli),colnames(target_exp))
target_cli=target_cli[com_sam,]
target_exp=target_exp[,com_sam]
dim(target_exp)
#18998    84
saveRDS(target_exp,file='00_pre_datas/TARGET/target_os_fpkm.RDS')
saveRDS(target_cli,file='00_pre_datas/TARGET/target_cli_os.RDS')
write.table(target_cli,'00_pre_datas/TARGET/target_cli.txt',quote = F,row.names = F,sep='\t')



####GSE21257####
GSE21257_exp2 <- as.data.frame(read.delim('00_rawdatas/GEO/GSE21257_OS_Gene_GPL10295.txt', 
                                         header = T, check.names = F,stringsAsFactors = F, row.names = 1))
GSE21257_cli <- read.delim('00_rawdatas/GEO/GSE21257_OS_SanmleInfo_GPL10295.txt',header = T, check.names = F,stringsAsFactors = F,)
colnames(GSE21257_cli)
GSE21257_cli <- GSE21257_cli[, c("rownames(geo.samp)", "times", "STATUS", 
                                 "1:AGE", "1:GENDER", "1:TUMOR LOCATION",
                                 "1:GROUP")]
colnames(GSE21257_cli) <- c('Samples', 'OS.time', 'OS', 'Age','Gender','Site', 'Metastases')
which(is.na(GSE21257_cli)==TRUE)
range(GSE21257_cli$OS.time)
GSE21257_cli$OS.time <- GSE21257_cli$OS.time * 30
table(GSE21257_cli$OS)
GSE21257_cli$OS <- ifelse(GSE21257_cli$OS == 'Alive', 0, 1)
GSE21257_cli$Age <- gsub(' months', '', GSE21257_cli$Age)
range(GSE21257_cli$Age)
GSE21257_cli$Age <- as.numeric(GSE21257_cli$Age) / 12
table(GSE21257_cli$Metastases)
GSE21257_cli$Metastases <- ifelse(GSE21257_cli$Metastases == 'No metastases', 'NO', 'YES')
rownames(GSE21257_cli)=GSE21257_cli$Samples
#GSE21257_cli=GSE21257_cli[GSE21257_cli$OS.time/365<10,]
GSE21257_exp=GSE21257_exp[,GSE21257_cli$Samples]
dim(GSE21257_exp)
#24998    53
range(GSE21257_exp)
saveRDS(GSE21257_cli,file = "00_pre_datas/GEO/GSE21257_cli.RDS")
saveRDS(GSE21257_exp,file = "00_pre_datas/GEO/GSE21257_exp.RDS")
write.table(GSE21257_cli,'00_pre_datas/GSE21257_cli.txt',quote = F,row.names = F,sep='\t')


#01.######
dir.create('01_GRGs')
GRGs.signature=read.xlsx('data/GRGs_GENES_PMID35808804.xlsx')
GRGs.signature=unique(GRGs.signature$GT)
length(GRGs.signature)
#185


###########
target.GRGs.score=t(ssGSEAScore_by_genes(gene.exp = target_exp,genes = GRGs.signature))
GSE39058.GRGs.score=t(ssGSEAScore_by_genes(gene.exp = GSE39058_exp,genes = GRGs.signature))

target.cli.merge=data.frame(target_cli,GRGs=target.GRGs.score[target_cli$Samples,1])
target.cli.merge=crbind2DataFrame(target.cli.merge)
head(target.cli.merge)
library(survival)
library(survminer)
library(RColorBrewer)
GRGs_score.col<- brewer.pal(12, "Set3")[6:5]
cutoff<-surv_cutpoint(target.cli.merge,time="OS.time",event="OS",variables='GRGs')
target.cli.merge$group <- ifelse(target.cli.merge$GRGs > cutoff$cutpoint$cutpoint, 'High', 'Low')
###################GRGs###########
target.GRGs_score.cli.use=target.cli.merge[order(target.cli.merge$GRGs),]
colnames(target.GRGs_score.cli.use)[9]='GRGs score group'
cli.colors=list()
for(i in c(2,4,7,9)){
  var=target.GRGs_score.cli.use[i]
  var.clas=names(table(target.GRGs_score.cli.use[,i]))
  var.color=ggsci::pal_d3()(10)[c(7,10,5,6)][1:length(var.clas)]
  names(var.color)=var.clas
  cli.colors=c(cli.colors,list(var.color))
}
names(cli.colors)=colnames(target.GRGs_score.cli.use)[c(2,4,7,9)]



library(ComplexHeatmap)

target.riskscore.barplot=columnAnnotation(GRGs_score = anno_barplot(target.GRGs_score.cli.use$GRGs,
                                                                    baseline =median(target.GRGs.score),
                                                                    bar_width=1,
                                                                    gp=gpar(fill=ifelse(target.GRGs_score.cli.use$GRGs>median(target.GRGs.score),
                                                                                        GRGs_score.col[1],GRGs_score.col[2])))
                                          ,annotation_name_side ='left'
                                          ,height =unit(4,'inches'))
draw(target.riskscore.barplot)

target.cli.heatmap=columnAnnotation(df = target.GRGs_score.cli.use[,c(2,4,7,9)]
                                    ,annotation_name_side='left'
                                    ,annotation_height =unit(4,'inches')
                                    ,col = cli.colors)

ht_list=target.riskscore.barplot %v% target.cli.heatmap
pdf('01_GRGs/fig1a.pdf',height = 8,width = 13,onefile = F)
ht_list
dev.off()





######GRGs####
km.target.os.GRGs=ggsurvplot(fit=survfit( Surv(OS.time/365, OS) ~ group,data = target.cli.merge),
           data=target.cli.merge,
           conf.int = F,pval = T,fun = "pct",risk.table =T, size = 0.7,
           surv.median.line = 'hv',
           linetype = c("solid", "dashed","strata")[1],
           palette = GRGs_score.col,
           #legend = c('top', 'bottom', 'left', 'right', 'none')[5],
           legend = c(0.8,0.85), # 
           legend.title = "")
km.target.os.GRGs=mg_merge_plot(km.target.os.GRGs$plot,km.target.os.GRGs$table,nrow=2,heights = c(3,1),align = 'v')
km.target.os.GRGs
ggsave("01_GRGs/target.km.OS_GRGs.pdf",height = 6,width = 6)

################################
target.GRGs.geneList=getGeneFC(gene.exp=target_exp[,target.cli.merge$Samples],group=target.cli.merge$group
                               ,ulab='High',dlab = 'Low')
h.all.gmt<-read.gmt("data/h.all.v2023.2.Hs.entrez.gmt")
target.GRGs.hallmark.gsea<-GSEA(target.GRGs.geneList,TERM2GENE = h.all.gmt,seed=T)
library(enrichplot)
library(ggplot2)
source('data/dotplotGsea.R')
gsea.dotplot=dotplotGsea(data = target.GRGs.hallmark.gsea,topn = 10,
                         order.by = 'NES')
pdf('01_GRGs/GSEA.pdf',height = 10,width = 9)
gsea.dotplot$plot+theme(text = element_text(family = 'Times'))+
  scale_y_discrete(labels=function(x) str_remove(x,"HALLMARK_"))
dev.off()

gsea.dotplot$df$Description


####02.#################
dir.create('02_Cluster')
#########
pre.genes=GRGs.signature
# pre.genes=intersect(tcga.cluster.degs ,rownames(tcga.degs))
length(pre.genes)#175
target_model_data=t(target_exp[pre.genes,target_cli$Samples])
colnames(target_model_data)=gsub('-','_',colnames(target_model_data))
target_model_data=merge(data.frame(Samples=target_cli$Samples,OS=target_cli$OS,OS.time=target_cli$OS.time),
                      data.frame(Samples=rownames(target_model_data),target_model_data),by='Samples')
rownames(target_model_data)=target_model_data$Samples
target_model_data=target_model_data[,-1]
target_model_data=crbind2DataFrame(target_model_data)
dim(target_model_data)
#  84 177
####
gene_sig_cox=cox_batch(dat = target_exp[pre.genes,target_cli$Samples],
                   time = target_cli$OS.time,event = target_cli$OS)
table(gene_sig_cox$p.value<.01)
# FALSE  TRUE 
#  168     7 
gene_sig_cox_res=gene_sig_cox[gene_sig_cox$p.value < 0.01, ]
gene_sig_cox_res$Type=ifelse(gene_sig_cox_res$HR>1,'Risk','Portect')
sig_fit_gene=rownames(gene_sig_cox_res)
length(sig_fit_gene)#7
writeMatrix(gene_sig_cox,outpath = '02_Cluster/sig_cox_res.txt')
writeMatrix(gene_sig_cox,outpath = 'files/sig_cox_res.txt')


########
pdf('02_Cluster/sig_cox_bioForest.pdf',height = 8,width = 8)
bioForest(rt = gene_sig_cox_res,col=c('blue','red'))
dev.off()



##############
library(ConsensusClusterPlus)
clusterAlg_name=c('hc','pam','km','kmdist')[1]
distance_name=c('pearson','spearman','euclidean','binary','maximum','canberra','minkowski')[1]
consen_gene=sig_fit_gene
#####TARGET#######
target_consen_data=as.matrix(target_exp[intersect(consen_gene,rownames(target_exp)),])
target_consen_data=t(scale(t(target_consen_data),scale = T))
#target_consen_data=sweep(target_consen_data,1,apply(target_consen_data, 1, mean))
#target_consen_data=sweep(target_consen_data,1,apply(target_consen_data, 1, median))
#target_consen_data=as.dist(1-cor(target_consen_data,method = 'pearson'))
target_consen_data=as.matrix(target_consen_data)
dim(target_consen_data)
target_clust_subtype <- ConsensusClusterPlus(target_consen_data
                                             , maxK = 10, reps = 500, pItem = 0.8
                                             , pFeature = 1
                                             , title = "Target_subtype"
                                             , clusterAlg = clusterAlg_name
                                             , distance = distance_name
                                             , plot = "pdf"
                                             , writeTable = T
                                             , seed = 123456)
k=2
target.subtype <- data.frame(Samples = names(target_clust_subtype[[k]]$consensusClass),
                             Cluster=target_clust_subtype[[k]]$consensusClass)
target.subtype$Cluster=paste0('C',target.subtype$Cluster)
writeMatrix(target.subtype,'02_Cluster/target.subtype.txt')
writeMatrix(target.subtype,'files/target.subtype.txt')

table(target.subtype$Cluster)

subcluster.color= c('gold1','deepskyblue3')
fig2a=ggplotKMCox(data.frame(time = target_cli[rownames(target.subtype),]$OS.time/365
                             , event = target_cli[rownames(target.subtype),]$OS
                             , target.subtype$Cluster)
                  ,add_text = '',title = 'Target'
                  ,palette = c('gold1','deepskyblue3'),labs = c('C1','C2'))
fig2a
ggsave("02_Cluster/fig2b.pdf",height = 6,width = 6)




###03.##########################
dir.create('03_DEGs')
############target
library(openxlsx)
p_fit=0.05
fc_fit=log2(1.5)
####c1 vs c2
target.degs.clucter=mg_limma_DEG(exp = target_exp[rownames(target_exp),target.subtype$Samples],group = as.character(target.subtype$Cluster),ulab = 'C1',dlab = 'C2')
target.degs.clucter$Summary
target.degs.clucter.up=rownames(target.degs.clucter$DEG[which(target.degs.clucter$DEG$adj.P.Val<p_fit & target.degs.clucter$DEG$logFC>fc_fit),])
target.degs.clucter.down=rownames(target.degs.clucter$DEG[which(target.degs.clucter$DEG$adj.P.Val<p_fit & target.degs.clucter$DEG$logFC<(-fc_fit)),])
length(target.degs.clucter.up)#108
length(target.degs.clucter.down)#173
target.cluster.degs=unique(c(target.degs.clucter.up,target.degs.clucter.down))
length(target.cluster.degs)
#281



write.xlsx(list(C1_vs_C2=cbind(gene=target.cluster.degs,
                               target.degs.clucter$DEG[which(target.degs.clucter$DEG$adj.P.Val<p_fit & abs(target.degs.clucter$DEG$logFC)>fc_fit),])
),
'03_DEGs/cluster_degs.xlsx',overwrite = T)

c1.degs.enrichment=mg_clusterProfiler(genes = target.degs.clucter.up)


write.xlsx(list(GO_BP=c1.degs.enrichment$GO_BP@result,
                GO_CC=c1.degs.enrichment$GO_CC@result,
                GO_MF=c1.degs.enrichment$GO_MF@result,
                KEGG=c1.degs.enrichment$KEGG@result),
           '03_DEGs/C1_DEGs_enrichment.xlsx',overwrite = T)


c2.degs.enrichment=mg_clusterProfiler(genes = target.degs.clucter.down)
fig3=list()
fig3[[1]]=barplot(c1.degs.enrichment$GO_BP)+scale_y_discrete(labels=function(y)str_wrap(y,width = 25))
fig3[[2]]=barplot(c2.degs.enrichment$KEGG)+scale_y_discrete(labels=function(y)str_wrap(y,width = 25))
fig3[[3]]=dotplot(c2.degs.enrichment$GO_BP)+scale_y_discrete(labels=function(y)str_wrap(y,width = 25))
fig3[[4]]=barplot(c2.degs.enrichment$GO_CC)+scale_y_discrete(labels=function(y)str_wrap(y,width = 25))

write.xlsx(list(GO_BP=c2.degs.enrichment$GO_BP@result,
                GO_CC=c2.degs.enrichment$GO_CC@result,
                GO_MF=c2.degs.enrichment$GO_MF@result,
                KEGG=c2.degs.enrichment$KEGG@result),
           '03_DEGs/C2_DEGs_enrichment.xlsx',overwrite = T)



fig3=mg_merge_plot(fig3,labels = LETTERS[1:4],common.legend = T,nrow=2,ncol=2,heights = c(1,1),hjust =0.5)
savePDF('03_DEGs/Fig3.pdf',fig3,height = 15,width = 20)


###########
dir.create('04_model')
target_cli=crbind2DataFrame(target_cli)
which(is.na(target_cli)==TRUE)
gene.cox=cox_batch(dat = target_exp[target.cluster.degs,target_cli$Samples],
                   time = target_cli$OS.time,event = target_cli$OS)
table(gene.cox$p.value<.01)
# FALSE  TRUE 
# 202    79 
gene.cox.fit=gene.cox[gene.cox$p.value<.01,]
pre.genes=rownames(gene.cox.fit)


target_model_data <- cbind(target_cli[, c("OS.time", "OS")],
                           t(target_exp[pre.genes, target_cli$Samples]))
colnames(target_model_data) <- gsub('-', '_', colnames(target_model_data))
dim(target_model_data)
# 84 81

#####LASSO####
library(glmnet)
set.seed(2024)
fit1=glmnet(as.matrix(target_model_data[,-c(1,2)])
            ,cbind(time=target_model_data$OS.time,
                   status=target_model_data$OS)
            ,family="cox"
            ,nlambda=100
            , alpha=1) 

cv.fit<-cv.glmnet(as.matrix( target_model_data[,-c(1,2)])
                  ,cbind(time=target_model_data$OS.time,
                         status=target_model_data$OS)
                  ,family="cox"
                  ,nfolds = 10
                  ,nlambda=100
                  , alpha=1)

sig.coef <- coefficients(cv.fit,s=cv.fit$lambda.min)[which(coefficients(cv.fit,s=cv.fit$lambda.min)[,1]!=0),1]
print(cv.fit$lambda.min)
length(names(sig.coef))
#16
pdf('04_model/LASSO.pdf',height = 5,width = 10,onefile = F)
par(mfrow=c(1,2))
plot(fit1)
plot(cv.fit)
dev.off()
#########
fmla <- as.formula(paste0("Surv(OS.time, OS) ~"
                          ,paste0(names(sig.coef),collapse = '+')))
cox <- coxph(fmla, data =as.data.frame(target_model_data))
cox=step(cox)

lan <- coef(cox)
lan
paste0(round(lan, 3), '*', names(lan),collapse = '+')
# "0.728*CORT+-0.999*LPAR5+-0.674*CEBPA+-0.386*MYH10+-1.465*MAGEA11"
write.csv(data.frame(gene=names(lan),coef=as.numeric(lan)),'04_model/gene_coef.csv',row.names = F)
module.coxforest=ggforest(cox, data = target_model_data, 
                          main = "Hazardratio", fontsize =1.0, 
                          noDigits = 2)
module.coxforest
ggsave('04_model/coxforest.pdf',module.coxforest,height = 5,width = 9)

risktype.col=c('lightpink2','turquoise3')



#####target#####
risk.target=as.numeric(lan%*%as.matrix(t(target_model_data[target_cli.merge$Samples,names(lan)])))
target.risktype.cli=data.frame(target_cli.merge,Riskscore=risk.target)
#######
target.data.point <- surv_cutpoint(target.risktype.cli, time = "OS.time", event = "OS",
                                   variables = 'Riskscore')
target.cutoff <- as.numeric(summary(target.data.point)[1])
target.cutoff
target.risktype.cli$Risktype=ifelse(target.risktype.cli$Riskscore>target.cutoff,'High','Low')
#target.risktype.cli$Risktype=ifelse(target.risktype.cli$Riskscore>median(target.risktype.cli$Riskscore),'High','Low')
target.km.OS=ggsurvplot(fit=survfit(Surv(OS.time/365, OS) ~ Risktype,
                                    data = target.risktype.cli),
                        data=target.risktype.cli,
                        conf.int = F,pval = T,risk.table = T, 
                        fun = "pct",size = 1,surv.median.line = 'hv',
                        title='target',legend.title='Risktype',
                        legend.labs = c('High','Low'),
                        linetype = c("solid", "dashed","strata")[1],
                        palette = risktype.col,ylab='Overall Survival(OS)',
                        legend=c(0.85,0.8),#
                        ggtheme = theme_classic())
target.km.OS=mg_merge_plot(target.km.OS$plot,target.km.OS$table,nrow=2,heights = c(3,1),align = 'v')
target.km.OS
ggsave("04_model/target.km.OS.pdf",height = 6,width = 6)
target.roc=ggplotTimeROC(target.risktype.cli$OS.time,
                         target.risktype.cli$OS,
                         target.risktype.cli$Riskscore,mks = c(1,3,5))
target.roc
ggsave("04_model/target.roc.pdf",height = 6,width = 6)
target.pca <- prcomp(t(target_exp[names(lan),target.risktype.cli$Samples]), scale=T)


# install.packages("ggbiplot")
library(ggbiplot)
ggbiplot(target.pca, scale=1, groups = target.risktype.cli$Risktype,
         ellipse = T,ellipse.prob=0.5, circle = F,var.axes=F) +
  scale_color_manual(values = risktype.col) + 
  # xlim(-5, 5) + ylim(-5,5) +
  theme_light() +theme(text = element_text(family = 'Times'))+
  theme(legend.direction = 'horizontal', legend.position = 'top') +
  xlab('PCA1') + ylab('PCA2')
ggsave('04_model/risktype_PCA.pdf',height = 5,width = 5)


my_mutiboxplot(t(target_exp[names(lan),target.risktype.cli$Samples]),
               group = target.risktype.cli$Risktype,legend.pos = 'top',group_cols = risktype.col,
               ylab = 'Gene Expression',angle = 0,hjust = .5)
ggsave('04_model/modelgene_expr.pdf',height = 5,width = 8)

  
##05.################# 
dir.create('05_External_data_validation')

######GSE21257#####
dir.create('03_External_data_validation')
dim(GSE21257_exp)
GSE21257_gene <- c("CORT", "GPR92","CEBPA","MYH10","MAGEA11")
GSE21257_model_data <- cbind(GSE21257_cli[, c("OS.time", "OS")],
                             t(GSE21257_exp[rownames(GSE21257_exp) %in% GSE21257_gene , GSE21257_cli$Samples]))
colnames(GSE21257_model_data) <- gsub('-', '_', colnames(GSE21257_model_data))

risk.GSE21257=as.numeric(lan%*%as.matrix(t(GSE21257_model_data[GSE21257_cli$Samples,GSE21257_gene ])))
GSE21257.risktype.cli=data.frame(GSE21257_cli,Riskscore=risk.GSE21257)
#######
GSE21257.data.point <- surv_cutpoint(GSE21257.risktype.cli, time = "OS.time", event = "OS",
                                     variables = 'Riskscore')
GSE21257.cutoff <- as.numeric(summary(GSE21257.data.point)[1])
GSE21257.cutoff
GSE21257.risktype.cli$Risktype=ifelse(GSE21257.risktype.cli$Riskscore>GSE21257.cutoff,'High','Low')
#GSE21257.risktype.cli$Risktype=ifelse(GSE21257.risktype.cli$Riskscore>median(GSE21257.risktype.cli$Riskscore),'High','Low')
GSE21257.km.OS=ggsurvplot(fit=survfit(Surv(OS.time/365, OS) ~ Risktype,
                                      data = GSE21257.risktype.cli),
                          data=GSE21257.risktype.cli,
                          conf.int = F,pval = T,risk.table = T, 
                          fun = "pct",size = 1,surv.median.line = 'hv',
                          title='GSE21257',legend.title='Risktype',
                          legend.labs = c('High','Low'),
                          linetype = c("solid", "dashed","strata")[1],
                          palette = risktype.col,ylab='Overall Survival(OS)',
                          legend=c(0.85,0.8),#
                          ggtheme = theme_classic())
GSE21257.km.OS=mg_merge_plot(GSE21257.km.OS$plot,GSE21257.km.OS$table,nrow=2,heights = c(3,1),align = 'v')
GSE21257.km.OS
ggsave("05_External_data_validation/GSE21257.km.OS.pdf",height =6,width = 6)
GSE21257.roc=ggplotTimeROC(GSE21257.risktype.cli$OS.time,
                           GSE21257.risktype.cli$OS,
                           GSE21257.risktype.cli$Riskscore,mks = c(1,3,5))
GSE21257.roc
ggsave("05_External_data_validation/GSE21257.roc.pdf",height =6,width = 6)


###06.####################
dir.create('06_risktype.immu')

#########
target.est=immu_estimate(target_exp)
save(target.est,file ='06_risktype.immu/target.est.RData')
load('06_risktype.immu/target.est.RData')
fig6a=mg_PlotMutiBoxplot(target.est[target.risktype.cli$Samples,]
                         , group = target.risktype.cli$Risktype
                         , legend.pos = 'top'
                         #, test_method = 'anova'
                         , test_method = 't.test'
                         , add = 'boxplot'
                         , ylab = 'Score'
                         , group_cols = risktype.col)
#########
target.immu.mcp=immu_MCPcounter(target_exp,isTCGA = F)
fig6b=mg_PlotMutiBoxplot(target.immu.mcp[target.risktype.cli$Samples,]
                   , group = target.risktype.cli$Risktype
                   , legend.pos = 'top'
                   #, test_method = 'anova'
                   , test_method = 't.test'
                   , add = 'boxplot'
                   , ylab = 'Score'
                   , group_cols = risktype.col)
#####ssGSEA####
# target.immu.ssgsea=immu_ssgsea(target_exp,isTCGA = F)
# save(target.immu.ssgsea,file ='05_risktype.immu/target.immu.ssgsea.RData')
load('06_risktype.immu/target.immu.ssgsea.RData')
fig6c=mg_PlotMutiBoxplot(target.immu.ssgsea[target.risktype.cli$Samples,]
                         , group = target.risktype.cli$Risktype
                         #, legend.pos = 'top'
                         #, test_method = 'anova'
                         , test_method = 't.test'
                         , add = 'boxplot'
                         , ylab = 'Score'
                         , group_cols = risktype.col)
# #####CIBERSORT#####
# #target.ciber=immu_CIBERSORT(exp_data = target_exp)
# #save(target.ciber,file='06_risktype.immu/target.ciber.RData')
# target.ciber <- readRDS("06_risktype.immu/cibersort_target_os.RDS")
# target.ciber <- target.ciber[intersect(target.risktype.cli$Samples,rownames(target.ciber)),]
# dim(target.ciber)
# # #
# fig6d=mg_PlotMutiBoxplot(target.ciber
#                          , group = target.risktype.cli[rownames(target.ciber),'Risktype']
#                          , legend.pos = 'top'
#                          #, test_method = 'anova'
#                          , test_method = 't.test'
#                          , add = 'boxplot'
#                          , ylab = 'Score'
#                          , group_cols = risktype.col)




fig6ab=mg_merge_plot(fig5a,fig5b,labels = c('A','B'),align = 'h',ncol = 2,nrow = 1,common.legend = T,widths = c(1,2))
fig6abc=mg_merge_plot(mg_merge_plot(fig5a,fig5b,widths = c(1,2),align = 'h',common.legend = T,labels = c('A','B')),
              fig5c,ncol = 1,nrow=2,labels = c('','C'))
savePDF('06_risktype.immu/fig6abc.pdf',fig6abc,height = 15,width = 15)


####07.###########
dir.create('07_module.sig')
##############
fig7a=list()
fig7a[[1]]=my_violin(dat = target.risktype.cli$Riskscore,group = target.cli.merge$Age1,
                     group_cols = pal_d3()(10)[c(7,10,5,6)],ylab = 'Riskscore')
fig7a[[2]]=my_violin(dat = target.risktype.cli$Riskscore,group = target.cli.merge$Gender,
                     group_cols = pal_d3()(10)[c(7,10,5,6)],ylab = 'Riskscore')
fig7a[[3]]=my_violin(dat = target.risktype.cli$Riskscore,group = target.cli.merge$Metaststic,
                     group_cols = pal_d3()(10)[c(7,10,5,6)],ylab = 'Riskscore')


fig7a=mg_merge_plot(fig7a,labels = LETTERS[1:3],common.legend = F,nrow=1,ncol=3,hjust =0.5)
savePDF('07_module.sig/fig7a.pdf',fig7a,height = 5,width = 12)



#07-2.####
target_cox_datas=target.risktype.cli
colnames(target_cox_datas)
target_cox_datas$Risktype <- factor(target_cox_datas$Risktype,levels = c('Low','High'))
table(target_cox_datas$Age1)
table(target_cox_datas$Metaststic)
#Age
Age_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Age,
                             data=target_cox_datas))
Age_sig_cox_dat <- data.frame(Names=rownames(Age_sig_cox[[8]]),
                              HR = round(Age_sig_cox[[7]][,2],3),
                              lower.95 = round(Age_sig_cox[[8]][,3],3),
                              upper.95 = round(Age_sig_cox[[8]][,4],3),
                              p.value=round(Age_sig_cox[[7]][,5],3))
Age_sig_cox_dat

#Gender
Gender_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Gender,
                                data=target_cox_datas))
Gender_sig_cox_dat <- data.frame(Names=rownames(Gender_sig_cox[[8]]),
                                 HR = round(Gender_sig_cox[[7]][,2],3),
                                 lower.95 = round(Gender_sig_cox[[8]][,3],3),
                                 upper.95 = round(Gender_sig_cox[[8]][,4],3),
                                 p.value=round(Gender_sig_cox[[7]][,5],3))
Gender_sig_cox_dat

#Metaststic
Metaststic_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Metaststic,
                                    data=target_cox_datas))
Metaststic_sig_cox_dat <- data.frame(Names=rownames(Metaststic_sig_cox[[8]]),
                                     HR = round(Metaststic_sig_cox[[7]][,2],3),
                                     lower.95 = round(Metaststic_sig_cox[[8]][,3],3),
                                     upper.95 = round(Metaststic_sig_cox[[8]][,4],3),
                                     p.value=round(Metaststic_sig_cox[[7]][,5],3))
Metaststic_sig_cox_dat

#risktype
risktype_sig_cox<-summary(coxph(formula=Surv(OS.time, OS)~Risktype,
                                data=target_cox_datas))
risktype_sig_cox_dat <- data.frame(Names=rownames(risktype_sig_cox[[8]]),
                                   HR = round(risktype_sig_cox[[7]][,2],3),
                                   lower.95 = round(risktype_sig_cox[[8]][,3],3),
                                   upper.95 = round(risktype_sig_cox[[8]][,4],3),
                                   p.value=round(risktype_sig_cox[[7]][,5],3))
risktype_sig_cox_dat

#riskscore
riskscore_sig_cox<-summary(coxph(formula=Surv(OS.time, OS)~Riskscore,
                                 data=target_cox_datas))
riskscore_sig_cox_dat <- data.frame(Names=rownames(riskscore_sig_cox[[8]]),
                                    HR = round(riskscore_sig_cox[[7]][,2],3),
                                    lower.95 = round(riskscore_sig_cox[[8]][,3],3),
                                    upper.95 = round(riskscore_sig_cox[[8]][,4],3),
                                    p.value=round(riskscore_sig_cox[[7]][,5],3))
riskscore_sig_cox_dat

sig_cox_dat <- rbind(Age_sig_cox_dat,
                     Gender_sig_cox_dat,
                     Metaststic_sig_cox_dat,
                     riskscore_sig_cox_dat)
data.sig <- data.frame(Names=sig_cox_dat$Names,
                       p.value=sig_cox_dat$p.value,
                       sig_cox_dat$HR,
                       sig_cox_dat$lower.95,
                       sig_cox_dat$upper.95)
data.sig <- crbind2DataFrame(data.sig)
data.sig$Names
rownames(data.sig) <- c("Age",
                        "Gender",
                        "Metaststic",
                        "Riskscore")
data.sig$Names <- rownames(data.sig)
data.sig
data.sig$p.value=ifelse(data.sig$p.value==0,'<0.001',data.sig$p.value)
pdf('07_module.sig/Fig07_sig.pdf', width = 8, height = 5,onefile = F)
mg_forestplot_v2(data.sig,xlog = T,colgap = 8,lineheight = 10)
dev.off()

#
muti_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~ Age  + Gender + Metaststic  + Riskscore, data=target_cox_datas))
muti_cox_dat <- data.frame(Names=rownames(muti_sig_cox[[8]]),
                           HR = round(muti_sig_cox[[7]][,2],3),
                           lower.95 = round(muti_sig_cox[[8]][,3],3),
                           upper.95 = round(muti_sig_cox[[8]][,4],3),
                           p.value=round(muti_sig_cox[[7]][,5],3))
data.muti <- data.frame(Names=muti_cox_dat$Names,
                        p.value=muti_cox_dat$p.value,
                        muti_cox_dat$HR,
                        muti_cox_dat$lower.95,
                        muti_cox_dat$upper.95)
data.muti <- crbind2DataFrame(data.muti)
rownames(data.muti) <- c("Age",
                         "Gender",
                         "Metaststic",
                         "Riskscore")
data.muti$Names <- rownames(data.muti)
data.muti
data.muti$p.value=ifelse(data.muti$p.value==0,'<0.001',data.muti$p.value)
pdf('07_module.sig/Fig07_muti.pdf', width = 8, height = 5,onefile = F)
mg_forestplot_v2(data.muti,xlog = T,colgap = 8,lineheight = 10)
dev.off()


#07-3.####

####
saveRDS(target_cox_datas,file ='07_module.sig/target_cox_datas.RDS' )
# mg_nomogram=function(clinical_riskscore,os,status,title='Nomogram',
#                      quick=T,mks = c(1,3,5)){
#   #clinical_riskscore=dat1[,3:5]
#   #os=dat1[,1]
#   #status=dat1[,2]
#   #sum(is.na(norm.stat.al[,3]))
#   norm.stat.al=data.frame(clinical_riskscore,time=os,status=status)
#   norm.stat.al=as.data.frame(norm.stat.al)
#   library(rms)
#   env <- globalenv()
#   env$MG_Grobal_DDSet <- rms::datadist(norm.stat.al) 
#   options(datadist='MG_Grobal_DDSet')
#   fmla <- as.formula(paste0("Surv(time, status) ~",paste0(colnames(clinical_riskscore),collapse = '+')))
#   cox2 <- cph(fmla, data = norm.stat.al,surv = T,x = T,y = T)
#   #summary(cox2)
#   #surv=Survival(cox2)
#   
#   fp <- predict(cox2)
#   cindex=getC_index(fp,norm.stat.al$time,norm.stat.al$status)
#   #cindex.orig=1-rcorr.cens(fp,Surv(norm.stat.al$time,norm.stat.al$status))
#   cut.time=c()
#   if(quantile(os[!is.na(os)])['75%']<12){
#     cut.time=mks
#   }else if(quantile(os[!is.na(os)])['75%']<365){
#     cut.time=c(12*mks[1],12*mks[2],12*mks[3])
#   }else{
#     cut.time=c(365*mks[1],365*mks[2],365*mks[3])
#   }
#   cut.time=cut.time[which(cut.time<quantile(os,seq(0,1,0.01))['100%'])]
#   print(cut.time)
#   #regplot(cox2)
#   #  print(regplot(cox3#
#   #,observation=pbc[2,] #
#   #
#   #              ,title=title
#   #              ,failtime = cut.time
#   #              ,prfail = TRUE #
#   #              ,showP = T #
#   #              ,droplines = F#
#   #,colors = mg_colors[1:3] #
#   #,rank="decreasing") #
#   #,interval="confidence"
#   #,rank="decreasing"
#   #,clickable=T
#   #              ,points=TRUE)) #
#   
#   #  plot(nom)
#   surv=Survival(cox2)
#   survs=list()
#   cal_all=list()
#   for(i in 1:length(cut.time)){
#     f1<-cph(formula = fmla,data=norm.stat.al,x=T,y=T,surv = T,na.action=na.delete,time.inc = cut.time[i]) 
#     cal1<-calibrate(f1, cmethod="KM", method="boot",u=cut.time[i],m=floor(sum(f1$n)/3)) 
#     cal_all=c(cal_all,list(cal1))
#     #    surv0 <- function(x)surv(cut.time[i],lp=x) 
#     #    survs=c(survs,list(surv0))
#   }
#   #f2<-cph(formula = fmla,data=norm.stat.al,x=T,y=T,surv = T,na.action=na.delete,time.inc = cut.time[2]) 
#   #cal3<-calibrate(f2, cmethod="KM", method="boot",u=cut.time[2],m=100,B=200) 
#   #f3<-cph(formula = fmla,data=norm.stat.al,x=T,y=T,surv = T,na.action=na.delete,time.inc = cut.time[3]) 
#   #cal5<-calibrate(f3, cmethod="KM", method="boot",u=cut.time[3],m=100,B=200) 
#   if(length(cut.time)==1){
#     surv1 <- function(x)surv(cut.time[1],lp=x) 
#     survs=list(surv1)
#   }else if(length(cut.time)==2){
#     surv1 <- function(x)surv(cut.time[1],lp=x) 
#     surv2 <- function(x)surv(cut.time[2],lp=x) 
#     survs=list(surv1,surv2)
#   }else if(length(cut.time)==3){
#     surv1 <- function(x)surv(cut.time[1],lp=x) 
#     surv2 <- function(x)surv(cut.time[2],lp=x) 
#     surv3 <- function(x)surv(cut.time[3],lp=x) 
#     survs=list(surv1,surv2,surv3)
#   }
#   nom=nomogram(cox2,fun=survs,lp= F
#                ,funlabel=c(paste0(mks[1], '-Year Survival'),
#                            paste0(mks[2], '-Year Survival'),
#                            paste0(mks[3], '-Year Survival'))[1:length(cut.time)]
#                ,maxscale=100
#                ,fun.at=seq(0,1,0.2)
#   )
#   
#   if(!quick){
#     cal_all=list()
#     for(i in 1:length(cut.time)){
#       cal1=get_best_calibrate(cox2,cut.time[i])
#       cal_all=c(cal_all,list(cal1))
#     }
#     #cal3=get_best_calibrate(cox2,cut.time[2])
#     #cal5=get_best_calibrate(cox2,cut.time[3])
#   }
#   lay2 <- customLayout::lay_new(matrix(1:2))
#   lay1 <- customLayout::lay_new(matrix(1:1))
#   cl <- customLayout::lay_bind_col(lay2, lay1, widths = c(1, 1.5)) 
#   #customLayout::lay_show(cl)
#   customLayout::lay_set(cl) 
#   
#   plot(cal_all[[1]],lwd = 2,lty = 1,errbar.col = mg_colors[1], bty = "l",xlim = c(0,1),ylim= c(0,1),xlab = "Nomogram-prediced OS (%)"
#        ,ylab = "Observed OS (%)",col = mg_colors[1],cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6,subtitles = F)
#   #lines(cal_all[[1]][,c('mean.predicted',"KM")],type = 'b', lwd = 1, col = mg_colors[1], pch = 16)
#   mtext("")
#   if(length(cal_all)>1){
#     for(i in 2:length(cal_all)){
#       plot(cal_all[[i]],lwd = 2,lty = 1,errbar.col = mg_colors[i],xlim = c(0,1),ylim= c(0,1),col = mg_colors[i],add = T)
#       #lines(cal_all[[i]][,c('mean.predicted',"KM")],type = 'b', lwd = 1, col = mg_colors[i], pch = 16)
#     }
#     #plot(cal3,lwd = 2,lty = 0,errbar.col = mg_colors[3],xlim = c(0,1),ylim= c(0,1),col = mg_colors[3],add = T)
#     #lines(cal3[,c('mean.predicted',"KM")],type = 'b', lwd = 1, col = mg_colors[3], pch = 16)
#   }
#   abline(0,1, lwd = 2, lty = 3, col = 'black')
#   legend("topleft", legend = c(paste0(mks[1], '-Year'),
#                                paste0(mks[2], '-Year'),
#                                paste0(mks[3], '-Year'))[1:length(cut.time)],col =mg_colors[1:length(cut.time)],lwd = 2,cex = 1.2,bty = "n")
#   
#   fp <- predict(cox2)
#   cindex=getC_index(fp,norm.stat.al$time,norm.stat.al$status)
#   dca_dat=data.frame(Nomogram=fp,status=norm.stat.al$status)
#   fp.al=cbind(fp)
#   for(i in 1:ncol(clinical_riskscore)){
#     fmla1 <- as.formula(paste0("Surv(time, status) ~",colnames(clinical_riskscore)[i]))
#     cox1 <- cph(fmla1, data = norm.stat.al,surv = T,x = T,y = T)
#     fp1 <- predict(cox1)
#     fp.al=cbind(fp.al,fp1)
#   }
#   colnames(fp.al)=c('Nomogram',colnames(clinical_riskscore))
#   fp.al=as.data.frame(fp.al)
#   fp.al$status=norm.stat.al$status
#   mg_plotDCA(fp.al$status
#              ,c('Nomogram',colnames(clinical_riskscore))
#              ,c('Nomogram',colnames(clinical_riskscore)),fp.al)
#   #plot(cal1,xlim=c(0,1),ylim=c(0,1))
#   #plot(cal3,xlim=c(0,1),ylim=c(0,1))
#   #plot(cal5,xlim=c(0,1),ylim=c(0,1))
#   plot(nom)
#   options(datadist=NULL)
#   return(list(Mod=cox2,Cindex=cindex,CutTime=cut.time))
# }
# library(devtools)
# install_github("zzawadz/customLayout")
# library(customLayout)
# pdf('07_module.sig/Fig07_Nomogram.pdf', width = 12, height = 10)
nom.plot=mg_nomogram(clinical_riskscore = data.frame(RiskScore=target_cox_datas$Riskscore,
                                                     Metaststic=target_cox_datas$Metaststic),
                     os = target_cox_datas$OS.time,
                     status = target_cox_datas$OS,mks = c(1,3,5))
dev.off()
# 
# 
# mg_nomogram_buti(nom.plot$Mod,cut.time = c(1*365,3*365,5*365))
# 
#08#####
dir.create('08_Immune_treatment')

##IMvigor210######
library("IMvigor210CoreBiologies")
data(cds)
pheno<-pData(cds)
head(pheno)
exper_tpm=mg_get_immu_pd1_treament_exp()
exper_id=exper_tpm$fpkm
exper_id$symbol=rownames(exper_id)
rownames(exper_id)<-exper_id$symbol
exper_id$symbol<-NULL
range(exper_id)
exper_id_use<-log2(exper_id+1)
dim(exper_id_use)
# 31085   348
range(exper_id_use)
#rownames(exper_id_use)=gsub('-','__',rownames(exper_id_use))
exper_id_use[1:5,1:5]
exper_id_use=exper_id_use[,rownames(pheno[which(pheno$binaryResponse!='NA'),])]
dim(exper_id_use)
# 31085   298
pheno=pheno[which(pheno$binaryResponse!='NA'),]

IMvigor210_model_data=data.frame(OS=pheno$censOS,OS.time = pheno$os,
                                 t(exper_id_use[intersect(pre.genes,rownames(exper_id_use)),rownames(pheno)]))
head(IMvigor210_model_data)

 fmla.IMvigor210 <- as.formula(paste0("Surv(OS.time, OS) ~"
                                      ,paste0(names(lan),collapse = '+')))
cox.IMvigor210 <- coxph(fmla.IMvigor210, data =as.data.frame(IMvigor210_model_data))
IMvigor210_lan <- coef(cox.IMvigor210)

#IMvigor210_lan <- lan
risk.imv210=as.numeric(IMvigor210_lan%*%as.matrix(t(IMvigor210_model_data[rownames(pheno),names(IMvigor210_lan)])))

imv210.risktype.cli=cbind.data.frame(pheno,Riskscore=risk.imv210)
#######
imv210.data.point <- surv_cutpoint(imv210.risktype.cli, time = "os", event = "censOS",
                                   variables = 'Riskscore')
imv210.cutoff <- as.numeric(summary(imv210.data.point)[1])
imv210.cutoff
imv210.risktype.cli$Risktype=ifelse(risk.imv210>imv210.cutoff,'High','Low')

imv210.roc=ggplotTimeROC(imv210.risktype.cli$os,
                         imv210.risktype.cli$censOS,
                         imv210.risktype.cli$Riskscore,mks = c(1,2,3,4,5))
imv210.roc


imv210.km=ggsurvplot(fit=survfit(Surv(os, censOS) ~ Risktype,
                                 data = imv210.risktype.cli),
                     data=imv210.risktype.cli,
                     conf.int = F,pval = T,risk.table = T, 
                     fun = "pct",size = 1,surv.median.line = 'hv',
                     title='IMvigor210',legend.title='Risktype',
                     legend.labs = c('High','Low'),
                     linetype = c("solid", "dashed","strata")[1],
                     palette = risktype.col,ylab='Overall Survival(OS)',
                     legend=c(0.85,0.8),#
                     ggtheme = custom_theme())
imv210.km=mg_merge_plot(imv210.km$plot,imv210.km$table,nrow=2,heights = c(3,1),align = 'v')
imv210.km



head(imv210.risktype.cli)
table(imv210.risktype.cli$binaryResponse)
imv210.bar=plotMutiBar(table(imv210.risktype.cli$binaryResponse,imv210.risktype.cli$Risktype))
imv210.bar

imv210.boxplot=ggplot(imv210.risktype.cli,aes(x=`Best Confirmed Overall Response`,y=Riskscore,fill=`Best Confirmed Overall Response`))+
  geom_boxplot()+stat_compare_means(aes(group=`Best Confirmed Overall Response`), label = 'p.format', method = 'kruskal.test')+
  theme(legend.position = 'none',text = element_text(family = 'Times',size = 12))
imv210.boxplot


##GSE78220###############
library(stringr)
GSE78220_cli <- getGEOSampleData('GSE78220')
GSE78220_cli1 <- GSE78220_cli[, c("Acc", "Title", "anti-pd-1 response", 
                                  "overall survival (days)", "vital status")]
colnames(GSE78220_cli1) <- c(c("Samples", "Title", "Rresponse", 
                               "OS.time", "OS"))
GSE78220_cli1 <- na.omit(GSE78220_cli1)
table(GSE78220_cli1$OS)
GSE78220_cli1$OS <- ifelse(GSE78220_cli1$OS == 'Alive', 0, 1)
rownames(GSE78220_cli1) <- GSE78220_cli1$Title

GSE78220_exp <- openxlsx::read.xlsx('00_rawdatas/GSE78220/GSE78220_PatientFPKM.xlsx',
                                    sheet = 1)
rownames(GSE78220_exp) <- GSE78220_exp$Gene
GSE78220_exp <- GSE78220_exp[, -1]
colnames(GSE78220_exp) <- str_split_fixed(colnames(GSE78220_exp), '\\.', 3)[, 1]
boxplot(GSE78220_exp[, 1:5])

GSE78220_exp <- log2(GSE78220_exp + 1)
boxplot(GSE78220_exp[, 1:5])

GSE78220_model_data <- cbind(GSE78220_cli1[, c("OS.time", 'OS')],
                             t(GSE78220_exp[, rownames(GSE78220_cli1)]))
GSE78220.genes <- intersect(names(lan), colnames(GSE78220_model_data))
GSE78220.genes

GSE78220_model_data=GSE78220_model_data[,c('OS.time','OS',GSE78220.genes)]

# GSE78220.genes=gsub('-','__',GSE78220.genes)
# colnames(GSE78220_model_data)=gsub('-','__',colnames(GSE78220_model_data))

fmla.GSE78220 <- as.formula(paste0("Surv(OS.time, OS) ~"
                                   ,paste0(GSE78220.genes,collapse = '+')))
cox.GSE78220 <- coxph(fmla.GSE78220, data =as.data.frame(GSE78220_model_data))
GSE78220_lan <- coef(cox.GSE78220)

risk.GSE78220=as.numeric(GSE78220_lan%*%as.matrix(t(GSE78220_model_data[GSE78220_cli1$Title,names(GSE78220_lan)])))

GSE78220_model_data$RS <- risk.GSE78220
GSE78220.data.point <- surv_cutpoint(GSE78220_model_data, time = "OS.time", event = "OS",
                                     variables = 'RS')
GSE78220.cutoff <- as.numeric(summary(GSE78220.data.point)[1])
GSE78220.cutoff

GSE78220.roc <- ggplotTimeROC(GSE78220_model_data$OS.time / 365,
                              GSE78220_model_data$OS,
                              risk.GSE78220,
                              mks = c(1,2,2.5))


GSE78220.risktype.cli=data.frame(GSE78220_cli1,
                                 Riskscore=risk.GSE78220,
                                 Risktype=ifelse(risk.GSE78220>=GSE78220.cutoff,'High','Low'))
GSE78220.km <-ggsurvplot(fit=survfit(Surv(OS.time/365, OS) ~ Risktype,
                                     data = GSE78220.risktype.cli),
                         data=GSE78220.risktype.cli,
                         conf.int = F,pval = T,risk.table = T, 
                         fun = "pct",size = 1,surv.median.line = 'hv',
                         title='GSE78220',legend.title='Risktype',
                         legend.labs = c('High','Low'),
                         linetype = c("solid", "dashed","strata")[1],
                         palette = risktype.col,ylab='Overall Survival(OS)',
                         legend=c(0.85,0.85),#
                         ggtheme = custom_theme()) 

GSE78220.km=mg_merge_plot(GSE78220.km$plot,GSE78220.km$table,nrow=2,heights = c(3,1),align = 'v')
GSE78220.km

GSE78220.bar=plotMutiBar(table(GSE78220.risktype.cli$Rresponse,GSE78220.risktype.cli$Risktype))
GSE78220.bar


GSE78220.boxplot=ggplot(GSE78220.risktype.cli,aes(x=Rresponse,y=Riskscore,fill=Rresponse))+
  geom_boxplot()+stat_compare_means(aes(group=Rresponse), label = 'p.format', method = 'kruskal.test')+
  theme(legend.position = 'none',text = element_text(family = 'Times',size = 12))+xlab('Response')
GSE78220.boxplot


p=mg_merge_plot(imv210.km,imv210.bar,imv210.boxplot,GSE78220.km,GSE78220.bar,GSE78220.boxplot,
                ncol=3,nrow=2,labels = LETTERS[1:6])
ggsave('08_Immune_treatment/Fig8.pdf',p,height = 10,width = 15)

save.image(file = 'OS_GRGs.RData')
