# Plot the survey-weighted joint distribution by period.
library(ggplot2)
d <- read.csv('outputs/joint/joint_distribution.csv')
parts <- do.call(rbind,strsplit(d$category,' / ',fixed=TRUE))
d$result <- factor(parts[,1], levels=c('Positive','Negative'))
d$dose <- factor(parts[,2],levels=c('No doses','Fewer than 2 excluding zero','At least 2','Unknown'),labels=c('Zero','<2, excl.\nzero','At least 2','Unknown'))
d$panel <- factor(d$period,levels=unique(d$period),labels=c('2017\u2013March 2020 (N = 716)','August 2021\u2013August 2023 (N = 376)'))
stopifnot(all(abs(tapply(d$estimate,d$panel,sum)-1)<1e-8))
g <- ggplot(d,aes(dose,100*estimate,fill=result))+geom_col(width=.68)+facet_wrap(~panel,nrow=1)+scale_fill_manual(values=c(Positive='#C47B39',Negative='#286A8D'),breaks=c('Negative','Positive'))+scale_y_continuous(breaks=seq(0,70,10),limits=c(0,72),expand=c(0,0))+labs(x='Self-reported hepatitis A vaccine doses',y='All adults with valid serology (%)',fill='Total anti-HAV')+theme_minimal(base_size=11,base_family='Arial')+theme(panel.grid.major.x=element_blank(),panel.grid.minor=element_blank(),strip.text=element_text(face='bold'),legend.position='bottom',axis.text=element_text(colour='black'))
ggsave('outputs/figures/Figure_2.tiff',g,width=183,height=110,units='mm',dpi=600,compression='lzw',bg='white')

ggsave('outputs/figures/Figure_2.png',g,width=183,height=110,units='mm',dpi=300,bg='white')
ggsave('outputs/figures/Figure_2.svg',g,width=183,height=110,units='mm',device=svglite::svglite,bg='white')
