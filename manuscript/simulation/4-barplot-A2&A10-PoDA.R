library(ggplot2)
library(foreach)
library(ggpubr)
library(ggnewscale)
library(gridExtra)
library(reshape2)
# A2
pseu <- readRDS("~/manuscript/simulation/Results/zero/A2_pseudo/A2pseudo5.rds")
imp <- readRDS("~/manuscript/simulation/Results/zero/A2_pseudo/A2imputation.rds")
LN<-readRDS("~/manuscript/simulation/Results/zero/A2_pseudo/A2multLN.rds")
Repl <- readRDS("~/manuscript/simulation/Results/zero/A2_pseudo/A2multRepl.rds")
#

# A10
pseu <- readRDS("~/manuscript/simulation/Results/zero/A10_pseudo/A10pseudo5.rds")
imp <- readRDS("~/manuscript/simulation/Results/zero/A10_pseudo/A10imputation.rds")
LN<-readRDS("~/manuscript/simulation/Results/zero/A10_pseudo/A10multLN.rds")
Repl <- readRDS("~/manuscript/simulation/Results/zero/A10_pseudo/A10multRepl.rds")
#



pseu_poda_power=cbind(pseu[,1],imp[,1],LN[,1],Repl[,1])
pseu_poda_fdr=cbind(pseu[,6],imp[,6],LN[,6],Repl[,6])
pseu_poda=cbind(pseu_poda_power,pseu_poda_fdr)


colind=c(2,1,3,4)
colind=c(colind,colind+4)
simlist=t(pseu_poda[1:3,colind])
simlist=simlist*100
rownames(simlist)=c("power.imputation","power.pseudo","power.LN","power.Repl",
                    "fdr.imputation","fdr.pseudo","fdr.LN","fdr.Repl")



##########bar plot#########
f=function(i){
  index=i
  data.frame(class=rep(c("Power(%)","Empirical FDR(%)"),each=4),
             method=rep(c('adaptive  ','naive ','multLN ','multRepl '),2),
             value=simlist[,index])
}

power.fdr=foreach(i=1:3, .combine=rbind) %do% f(i)



out=data.frame(sig=rep(c(50,100,200),each=8),
               q=rep(c(0.5,5),each=4),
               power.fdr)


colnames(out)=c("sig.prob","q","class","method","value")
out$sig.prob=factor(out$sig.prob, levels =c("50","100","200"))
out$class=factor(out$class, levels =c("Power(%)","Empirical FDR(%)"))
out$method=factor(out$method,levels=c('adaptive  ','naive ','multLN ','multRepl '))
out$q=factor(out$q, levels =c("0.5","5"))
line=data.frame(class=factor(c("Power(%)","Empirical FDR(%)"), levels =c("Power(%)","Empirical FDR(%)")),y=c(NA,5))





p1 = ggplot(out[out$q=="5",], aes(x=sig.prob, y=value, fill=method))+ theme_bw()
ptest1 = p1+geom_bar(stat="identity",position=position_dodge(0.7),width = 0.6)+
  scale_fill_manual(values=c('red','#8074C8','#00A08799','#F0C284'))+
  facet_grid(vars(class), labeller = labeller(.cols = label_both))+ylim(c(0,100))+
  geom_hline(data= line, aes(yintercept=y),linetype = "dashed")+
  theme(legend.position = "top", legend.title = element_blank(),legend.text = element_text(size = 6),
        legend.key.size = unit(0.3, "cm"),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6) ,
        panel.spacing = unit(0.5, "lines"))+
  guides(fill = guide_legend(ncol = 6))+
  ggnewscale::new_scale_fill()+
  geom_bar(data=out[out$q=="0.5",], aes(x=sig.prob, y=value, fill=method),stat="identity",position=position_dodge(0.7),
           width = 0.6,
           show.legend=FALSE)+
  scale_fill_manual(values=c('red','#8074C8','#00A08799','#F0C284'))+
  facet_grid(vars(class), labeller = labeller(.cols = label_both))+ylim(c(0,100))+
  labs( y = ' ',x=expression(n))
ptest1


colind=c(2,1,3,4)
colind=c(colind,colind+4)
simlist=t(pseu_poda[4:6,colind])
simlist=simlist*100
rownames(simlist)=c("power.imputation","power.pseudo","power.LN","power.Repl",
                    "fdr.imputation","fdr.pseudo","fdr.LN","fdr.Repl")



##########bar plot#########
f=function(i){
  index=i
  data.frame(class=rep(c("Power(%)","Empirical FDR(%)"),each=4),
             method=rep(c('adaptive  ','naive ','multLN ','multRepl '),2),
             value=simlist[,index])
}

power.fdr=foreach(i=1:3, .combine=rbind) %do% f(i)



out=data.frame(sig=rep(c(50,100,200),each=8),
               q=rep(c(0.5,5),each=4),
               power.fdr)


colnames(out)=c("sig.prob","q","class","method","value")
out$sig.prob=factor(out$sig.prob, levels =c("50","100","200"))
out$class=factor(out$class, levels =c("Power(%)","Empirical FDR(%)"))
out$method=factor(out$method,levels=c('adaptive  ','naive ','multLN ','multRepl '))
out$q=factor(out$q, levels =c("0.5","5"))
line=data.frame(class=factor(c("Power(%)","Empirical FDR(%)"), levels =c("Power(%)","Empirical FDR(%)")),y=c(NA,5))






p1 = ggplot(out[out$q=="5",], aes(x=sig.prob, y=value, fill=method))+ theme_bw()
ptest2 = p1+geom_bar(stat="identity",position=position_dodge(0.7),width = 0.6)+
  scale_fill_manual(values=c('red','#8074C8','#00A08799','#F0C284'))+
  facet_grid(vars(class), labeller = labeller(.cols = label_both))+ylim(c(0,100))+
  geom_hline(data= line, aes(yintercept=y),linetype = "dashed")+
  theme(legend.position = "top", legend.title = element_blank(),legend.text = element_text(size = 6),
        legend.key.size = unit(0.3, "cm"),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6) ,
        panel.spacing = unit(0.5, "lines"))+
  guides(fill = guide_legend(ncol = 6))+
  ggnewscale::new_scale_fill()+
  geom_bar(data=out[out$q=="0.5",], aes(x=sig.prob, y=value, fill=method),stat="identity",show.legend=FALSE,
           position=position_dodge(0.7),width = 0.6)+
  scale_fill_manual(values=c('red','#8074C8','#00A08799','#F0C284'))+
  facet_grid(vars(class), labeller = labeller(.cols = label_both))+ylim(c(0,100))+
  labs( y = ' ',x=expression(n))


colind=c(2,1,3,4)
colind=c(colind,colind+4)
simlist=t(pseu_poda[7:9,colind])
simlist=simlist*100
rownames(simlist)=c("power.imputation","power.pseudo","power.LN","power.Repl",
                    "fdr.imputation","fdr.pseudo","fdr.LN","fdr.Repl")



##########bar plot#########
f=function(i){
  index=i
  data.frame(class=rep(c("Power(%)","Empirical FDR(%)"),each=4),
             method=rep(c('adaptive  ','naive ','multLN ','multRepl '),2),
             value=simlist[,index])
}

power.fdr=foreach(i=1:3, .combine=rbind) %do% f(i)



out=data.frame(sig=rep(c(50,100,200),each=8),
               q=rep(c(0.5,5),each=4),
               power.fdr)


colnames(out)=c("sig.prob","q","class","method","value")
out$sig.prob=factor(out$sig.prob, levels =c("50","100","200"))
out$class=factor(out$class, levels =c("Power(%)","Empirical FDR(%)"))
out$method=factor(out$method,levels=c('adaptive  ','naive ','multLN ','multRepl '))
out$q=factor(out$q, levels =c("0.5","5"))
line=data.frame(class=factor(c("Power(%)","Empirical FDR(%)"), levels =c("Power(%)","Empirical FDR(%)")),y=c(NA,5))




p1 = ggplot(out[out$q=="5",], aes(x=sig.prob, y=value, fill=method))+ theme_bw()
ptest3 = p1+geom_bar(stat="identity",position=position_dodge(0.7),width = 0.6)+
  scale_fill_manual(values=c('red','#8074C8','#00A08799','#F0C284'))+
  facet_grid(vars(class), labeller = labeller(.cols = label_both))+ylim(c(0,100))+
  geom_hline(data= line, aes(yintercept=y),linetype = "dashed")+
  theme(legend.position = "top", legend.title = element_blank(),legend.text = element_text(size = 6),
        legend.key.size = unit(0.3, "cm"),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6) ,
        panel.spacing = unit(0.5, "lines"),strip.text = element_text(size = 12))+
  guides(fill = guide_legend(ncol = 6))+
  ggnewscale::new_scale_fill()+
  geom_bar(data=out[out$q=="0.5",], aes(x=sig.prob, y=value, fill=method),stat="identity",show.legend=FALSE,
           position=position_dodge(0.7),width = 0.6)+
  scale_fill_manual(values=c('red','#8074C8','#00A08799','#F0C284'))+
  facet_grid(vars(class), labeller = labeller(.cols = label_both))+ylim(c(0,100))+
  labs( y = ' ',x=expression(n))


ptest1 <- ptest1 + theme(plot.margin = margin(1, 1, 1, 1))
ptest2 <- ptest2 + theme(plot.margin = margin(1, 1, 1, 1))
ptest3 <- ptest3 + theme(plot.margin = margin(1, 1, 1, 1))




ptest1 <- ptest1 +
  theme(
    axis.title.x = element_blank(),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.5, "cm"),
    axis.title.y = element_text(size = 10),
    plot.margin = margin(1, 1, 1, 1),
    strip.background = element_blank(),
    strip.text = element_blank()
  )

ptest2 <- ptest2 +
  theme(
    axis.title.x = element_blank(),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.5, "cm"),
    plot.margin = margin(1, 1, 1, 1),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank()
  )

ptest3 <- ptest3 +
  theme(
    axis.title.x = element_blank(),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.5, "cm"),
    plot.margin = margin(1, 1, 1, 1),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank()
  )


p.comb <- ggarrange(
  ptest1, ptest2, ptest3,
  ncol = 3,
  common.legend = TRUE,
  legend = "top"
) +
  theme(
    axis.title.y = element_text(size = 12),
    legend.text = element_text(size = 20),
    legend.key.size = unit(10, "cm")
  )

p.comb

p.comb <- annotate_figure(
  p.comb,
  bottom = text_grob(expression(n), size = 12)
)

# ggsave("~/manuscript/simulation/Plots/zero/A2new.png", p.comb, width = 4*3.5, height =2.5*3.5, limitsize = FALSE ,bg = "white")
