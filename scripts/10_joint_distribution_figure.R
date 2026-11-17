# Figure 2 contract: show the full joint distribution, including unknown history.
# One quantitative panel; each point uses all valid HAV serology as denominator.
# Period-specific survey-logit 95% CIs; no pooled national estimate.
library(ggplot2)
z<-read.csv('outputs/joint/joint_distribution.csv')
ps<-unique(z$period); cats<-unique(z$category)
z$category<-factor(z$category,levels=rev(cats));z$period<-factor(z$period,levels=ps)
labels<-c('Anti-HAV negative: No doses','Anti-HAV negative: <2 doses, excluding zero','Anti-HAV negative: At least 2 doses','Anti-HAV negative: Unknown history','Anti-HAV positive: No doses','Anti-HAV positive: <2 doses, excluding zero','Anti-HAV positive: At least 2 doses','Anti-HAV positive: Unknown history')
p<-ggplot(z,aes(x=estimate*100,y=category,colour=period,shape=period))+geom_vline(xintercept=seq(0,60,10),colour='#E5E7EB',linewidth=.3)+geom_errorbar(aes(xmin=ci_low*100,xmax=ci_high*100),orientation='y',position=position_dodge(width=.45),width=.16,linewidth=.5)+geom_point(position=position_dodge(width=.45),size=2)+scale_y_discrete(labels=setNames(labels,cats))+scale_x_continuous(limits=c(0,60),breaks=seq(0,60,10),expand=expansion(mult=c(.01,.01)))+scale_colour_manual(values=c('#286A8D','#B46A32'))+scale_shape_manual(values=c(16,17))+labs(x='Survey-weighted proportion of all valid HAV results (%)',y=NULL,colour=NULL,shape=NULL)+theme_classic(base_size=10,base_family='Arial')+theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),axis.text=element_text(colour='black',size=9),axis.title.x=element_text(size=10,margin=margin(t=10)),legend.position='bottom',legend.text=element_text(size=9),plot.margin=margin(12,12,6,6),legend.margin=margin(t=8))
# Single plot area: multi-panel alignment is not applicable.
ggsave('outputs/figures/Figure_2.png',p,width=183,height=125,units='mm',dpi=300,bg='white')
ggsave('outputs/figures/Figure_2.tiff',p,width=183,height=125,units='mm',dpi=600,compression='lzw',bg='white')
svglite::svglite('outputs/figures/Figure_2.svg',width=183/25.4,height=125/25.4)
print(p)
dev.off()
