#Cox-Knockoff 
common_TTE <- read.table("E:/AT/knockoff相关资料/Cox konckoff/Cox_Knockoff/result_knockoff_Common_TTE.txt",header=F)
fdr_level=common_TTE$V4
power <- common_TTE$V5
fdr <- common_TTE$V6
Results_Cox_knock=cbind(fdr_level,power,fdr)

#SPACox-BH 
#common_TTE_SPACox_BH <- read.table("./Cox_Knockoff/result.BH_common_TTE.txt",header=F)
common_TTE_SPACox_BH <- read.table("E:/AT/knockoff相关资料/Cox konckoff/Cox_Knockoff/result.BH_common_TTE.txt",header=F)
power_SPACox_BH <- common_TTE_SPACox_BH$V5
fdr_SPACox_BH <- common_TTE_SPACox_BH$V6
Results_SPACox_BH=cbind(fdr_level,power_SPACox_BH,fdr_SPACox_BH)

library(ggplot2)
library(RColorBrewer)

cols=brewer.pal(9, "Paired")[c(2,9)]
spline.power <- as.data.frame(spline(fdr_level, power))
spline.power_SPACox_BH <- as.data.frame(spline(fdr_level, power_SPACox_BH))
Power=ggplot(Results_Cox_knock, aes(fdr_level, power)) +
  geom_smooth(data = spline.power, aes(x=x, y=y,color="Cox-Knockoff"), size=2, se = FALSE)+
  geom_smooth(data = spline.power_SPACox_BH, aes(x=x, y=y,color="SPACox-BH"), size=2, se = FALSE)+
  scale_y_continuous(limits = c(-0.05,1))+scale_x_continuous(limits =  c(0, 0.2))+
  labs(x = "Target FDR", y = "Power")+
  theme(text = element_text(size = 35),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Remove panel background
        panel.background = element_rect(fill = "white", linetype = "solid",color = "black"),
        # Add axis line
        axis.line = element_line(color = "black"),
        legend.position = 'top',
        legend.text = element_text(size = 35),plot.margin=unit(c(0.5,0.5,-0.1,0.5), "cm"))+
  guides(color = guide_legend(nrow = 1,byrow=TRUE)) +
  scale_color_manual(name = "", 
                     breaks = c("Cox-Knockoff","SPACox-BH"),
                     values = c(`Cox-Knockoff` =cols[1],`SPACox-BH` =cols[2]))
Power

spline.fdr <- as.data.frame(spline(fdr_level, fdr))
spline.fdr_SPACox_BH <- as.data.frame(spline(fdr_level, fdr_SPACox_BH))
FDR=ggplot(Results_Cox_knock, aes(fdr_level, fdr)) +
  geom_smooth(data = spline.fdr, aes(x=x, y=y), col=cols[1],se = FALSE, size=2)+
  geom_smooth(data = spline.fdr_SPACox_BH, aes(x=x, y=y), col=cols[2],se = FALSE, size=2)+
  geom_abline(aes(slope=1, intercept=0, colour="Expected FDR"), linetype="dotdash",size=2)+
  scale_y_continuous(limits = c(-0.05,0.85))+scale_x_continuous(limits =  c(0, 0.2))+
  labs(x = "Target FDR", y = "Observed FDR")+ 
  theme(text = element_text(size = 35),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Remove panel background
        panel.background = element_rect(fill = "white", linetype = "solid",colour = "black"),
        # Add axis line
        axis.line = element_line(colour = "black"),
        legend.position =c(0.2, 0.36),
        legend.text = element_text(size = 25),plot.margin=unit(c(-0.1,0.5,0.5,0.5), "cm"))+
  guides(color = guide_legend(nrow =1)) +
  scale_color_manual(name = "", values = c(`Expected FDR` = "gray"))
FDR

library(ggpubr)
arrange <- ggarrange(Power,FDR, align = "v",ncol = 1, nrow = 2,heights = c(0.65,0.35))
arrange

Cox_Knockoff<-annotate_figure(arrange,
                                top = text_grob("Power and FDR of Cox-knockoff", size = 40)
)
Cox_Knockoff

png(file=paste0('Simulations_Cox_Knockoff.png'),width=1200,height=800,units='px',pointsize =40)
Cox_Knockoff
dev.off()

