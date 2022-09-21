require(fitur)
require("ggstance")
require(ggplot2)
require(Hmisc)
require(grid)



#reading dataset
bias <- read.table("dataset_figure2", 
                    header = TRUE)




###############################################################
##################Figure 2####################################
##############################################################


bias$rho <- ordered(bias$rho)


#adding labels
bias$rho <- factor(bias$rho, labels = c("rho == -0.8", "rho == 0.2", "rho == 0.5","rho ==0.8"))

bias$par <- factor(bias$par, labels = c("hat(beta)[11]", "hat(beta)[12]", "hat(beta)[21]","hat(beta)[22]","hat(phi)[11]","hat(phi)[12]","hat(phi)[21]","hat(phi)[22]","hat(rho)"))




p1 <- ggplot(data=bias[bias$par == "hat(beta)[11]",],
       aes(x =as.factor(n),y = bias, shape=as.factor(n)))+geom_hline(yintercept = 0               ,linetype = "dotted")+
   geom_pointrange(aes(ymin = lwr, ymax = upr,shape=as.factor(n                        )),color="grey40")+
     ylab("")+ylim(-0.2,0.2)+xlab(expression(hat(beta)[11]))+
   facet_grid(~rho, 
              labeller = label_parsed) +
      theme(legend.title = element_text(size = 18),
            legend.text=element_text(size=18),
            axis.text = element_text(color="gray50",size = 10),
            plot.margin = unit(c(0,0.5,-0.2,0), "cm"),
            legend.position="top",
            plot.title=element_text(size=18,face="bold"),
         panel.background = element_rect(fill = "white", colour = "grey40"),
         axis.text.y=element_blank(),
         axis.ticks.y=element_blank(),
         axis.title.x=element_text("",face="bold", size=18),
          axis.title.y=element_text(angle=360,vjust = 0.5,face="bold", size=18),
         axis.text.x=element_text(face="bold"),
         axis.title=element_text(size=12,face="bold"),
         strip.text.y = element_text(size=12,angle=180,face="bold"),strip.text.x = element_text(size=18,face="bold"))+
   scale_shape_manual("Sample size",values=c(21, 22,23,24))+
   coord_flip() 


p2 <- ggplot(data=bias[bias$par == "hat(beta)[12]",],
             aes(x = par:as.factor(n),y = bias, ymin = lwr, ymax = upr,shape=as.factor(n)))+
   geom_pointrange(aes(shape=as.factor(n)),size=0.5,color="grey40")+geom_hline(yintercept = 0,linetype = "dotted")+
   ylab("")+xlab(expression(hat(beta)[12]))+
   facet_grid(~rho, switch = "y",
              labeller = label_parsed,scales= 'free') +ylim(-0.2,0.2)+
   theme(axis.text = element_text(color="gray50",size = 10),plot.margin = unit(c(-0.2,0.5,-0.2,0), "cm"),legend.position="",legend.direction = 'horizontal',plot.title=element_text(size=16,face="bold"),
         panel.background = element_rect(fill = "white", colour = "grey40"),
         axis.text.y=element_blank(),
         axis.ticks.y=element_blank(),
         axis.title.x=element_text("",face="bold", size=18),
         axis.title.y=element_text(angle=360,vjust = 0.5,face="bold", size=18),
         axis.text.x=element_text(face="bold"),
         axis.title=element_text(size=12,face="bold"),
         strip.background = element_blank(),
         strip.text.x = element_blank())+scale_shape_manual("",values=c(21, 22,23,24))+
   coord_flip()


p3 <- ggplot(data=bias[bias$par == "hat(beta)[21]",],
             aes(x = par:as.factor(n),y = bias, ymin = lwr, ymax = upr,shape=as.factor(n)))+
   geom_pointrange(aes(shape=as.factor(n)),size=0.5,color="grey40")+geom_hline(yintercept = 0,linetype = "dotted")+
   ylab("")+xlab(expression(hat(beta)[21]))+
   facet_grid(~rho, switch = "y",
              labeller = label_parsed,scales= 'free') +ylim(-2,2)+
   theme(axis.text = element_text(color="gray50",size = 10),plot.margin = unit(c(-0.2,0.5,-0.2,0), "cm"),legend.position="",legend.direction = 'horizontal',plot.title=element_text(size=16,face="bold"),
         panel.background = element_rect(fill = "white", colour = "grey40"),
         axis.text.y=element_blank(),
         axis.ticks.y=element_blank(),
         axis.title.x=element_text("",face="bold", size=18),
         axis.title.y=element_text(angle=360,vjust = 0.5,face="bold", size=18),
         axis.text.x=element_text(face="bold"),
         axis.title=element_text(size=12,face="bold"),
         strip.background = element_blank(),
         strip.text.x = element_blank())+scale_shape_manual("",values=c(21, 22,23,24))+
   coord_flip()


p4 <- ggplot(data=bias[bias$par == "hat(beta)[22]",],
             aes(x = par:as.factor(n),y = bias, ymin = lwr, ymax = upr,shape=as.factor(n)))+
   geom_pointrange(aes(shape=as.factor(n)),size=0.5,color="grey40")+geom_hline(yintercept = 0,linetype = "dotted")+
   ylab("")+xlab(expression(hat(beta)[22]))+
   facet_grid(~rho, switch = "y",
              labeller = label_parsed,scales= 'free') +ylim(-25,25)+
   theme(axis.text = element_text(color="gray50",size = 10),plot.margin = unit(c(-0.2,0.5,-0.2,0), "cm"),legend.position="",legend.direction = 'horizontal',plot.title=element_text(size=16,face="bold"),
         panel.background = element_rect(fill = "white", colour = "grey40"),
         axis.text.y=element_blank(),
         axis.ticks.y=element_blank(),
         axis.title.x=element_text("",face="bold", size=18),
         axis.title.y=element_text(angle=360,vjust = 0.5,face="bold", size=18),
         axis.text.x=element_text(face="bold"),
         axis.title=element_text(size=12,face="bold"),
         strip.background = element_blank(),
         strip.text.x = element_blank())+scale_shape_manual("",values=c(21, 22,23,24))+
   coord_flip()


p5 <- ggplot(data=bias[bias$par == "hat(phi)[11]",],
             aes(x = par:as.factor(n),y = bias, ymin = lwr, ymax = upr,shape=as.factor(n)))+
   geom_pointrange(aes(shape=as.factor(n)),size=0.5,color="grey40")+geom_hline(yintercept = 0,linetype = "dotted")+
   ylab("")+xlab(expression(hat(Phi)[11]))+
   facet_grid(~rho, switch = "y",
              labeller = label_parsed,scales= 'free') +ylim(-0.5,0.5)+
   theme(axis.text = element_text(color="gray50",size =10),plot.margin = unit(c(-0.2,0.5,-0.2,0), "cm"),legend.position="",legend.direction = 'horizontal',plot.title=element_text(size=16,face="bold"),
         panel.background = element_rect(fill = "white", colour = "grey40"),
         axis.text.y=element_blank(),
         axis.ticks.y=element_blank(),
         axis.title.x=element_text("",face="bold", size=18),
         axis.title.y=element_text(angle=360,vjust = 0.5,face="bold", size=18),
         axis.text.x=element_text(face="bold"),
         axis.title=element_text(size=12,face="bold"),
         strip.background = element_blank(),
         strip.text.x = element_blank())+scale_shape_manual("",values=c(21, 22,23,24))+
   coord_flip()


p6 <- ggplot(data=bias[bias$par == "hat(phi)[12]" ,],
             aes(x = par:as.factor(n),y = bias, ymin = lwr, ymax = upr,shape=as.factor(n)))+
   geom_pointrange(aes(shape=as.factor(n)),size=0.5,color="grey40")+geom_hline(yintercept = 0,linetype = "dotted")+
   ylab("")+xlab(expression(hat(Phi)[12]))+
   facet_grid(~rho, switch = "y",
              labeller = label_parsed,scales= 'free') +ylim(-2,2)+
   theme(axis.text = element_text(color="gray50",size = 10),plot.margin = unit(c(-0.2,0.5,-0.2,0), "cm"),legend.position="",legend.direction = 'horizontal',plot.title=element_text(size=16,face="bold"),
         panel.background = element_rect(fill = "white", colour = "grey40"),
         axis.text.y=element_blank(),
         axis.ticks.y=element_blank(),
         axis.title.x=element_text("",face="bold", size=18),
         axis.title.y=element_text(angle=360,vjust = 0.5,face="bold", size=18),
         axis.text.x=element_text(face="bold"),
         axis.title=element_text(size=12,face="bold"),
         strip.background = element_blank(),
         strip.text.x = element_blank())+scale_shape_manual("",values=c(21, 22,23,24))+
   coord_flip()


p7 <- ggplot(data=bias[bias$par == "hat(phi)[21]",],
             aes(x = par:as.factor(n),y = bias, ymin = lwr, ymax = upr,shape=as.factor(n)))+
   geom_pointrange(aes(shape=as.factor(n)),size=0.5,color="grey40")+geom_hline(yintercept = 0,linetype = "dotted")+
   ylab("")+xlab(expression(hat(Phi)[21]))+
   facet_grid(~rho, switch = "y",
              labeller = label_parsed,scales= 'free') +ylim(-1,1)+
   theme(axis.text = element_text(color="gray50",size = 10),plot.margin = unit(c(-0.2,0.5,-0.2,0), "cm"),legend.position="",legend.direction = 'horizontal',plot.title=element_text(size=16,face="bold"),
         panel.background = element_rect(fill = "white", colour = "grey40"),
         axis.text.y=element_blank(),
         axis.ticks.y=element_blank(),
         axis.title.x=element_text("",face="bold", size=18),
         axis.title.y=element_text(angle=360,vjust = 0.5,face="bold", size=18),
         axis.text.x=element_text(face="bold"),
         axis.title=element_text(size=12,face="bold"),
         strip.background = element_blank(),
         strip.text.x = element_blank())+scale_shape_manual("",values=c(21, 22,23,24))+
   coord_flip()


p8 <- ggplot(data=bias[bias$par == "hat(phi)[22]",],
             aes(x = par:as.factor(n),y = bias, ymin = lwr, ymax = upr,shape=as.factor(n)))+
   geom_pointrange(aes(shape=as.factor(n)),size=0.5,color="grey40")+geom_hline(yintercept = 0,linetype = "dotted")+
   ylab("")+xlab(expression(hat(Phi)[22]))+
   facet_grid(~rho, switch = "y",
              labeller = label_parsed,scales= 'free') +ylim(-8,8)+
   theme(axis.text = element_text(color="gray50",size = 10),plot.margin = unit(c(-0.2,0.5,-0.2,0), "cm"),legend.position="",legend.direction = 'horizontal',plot.title=element_text(size=16,face="bold"),
         panel.background = element_rect(fill = "white", colour = "grey40"),
         axis.text.y=element_blank(),
         axis.ticks.y=element_blank(),
         axis.title.x=element_text("",face="bold", size=18),
         axis.title.y=element_text(angle=360,vjust = 0.5,face="bold", size=18),
         axis.text.x=element_text(face="bold"),
         axis.title=element_text(size=12,face="bold"),
         strip.background = element_blank(),
         strip.text.x = element_blank())+scale_shape_manual("",values=c(21,22,23,24))+
   coord_flip()


p9 <- ggplot(data=bias[bias$par == "hat(rho)" ,],
             aes(x = par:as.factor(n),y = bias, ymin = lwr, ymax = upr,shape=as.factor(n)))+
   geom_pointrange(aes(shape=as.factor(n)),size=0.5,color="grey40")+geom_hline(yintercept = 0,linetype = "dotted")+
   ylab("Bias")+xlab(expression(hat(rho)))+
   facet_grid(~rho, switch = "y",
              labeller = label_parsed,scales= 'free') +ylim(-0.3,0.3)+
   theme(axis.text = element_text(color="gray50",size = 10),plot.margin = unit(c(-0.2,0.5,0,0.4), "cm"),legend.position="",legend.direction = 'horizontal',plot.title=element_text(size=16,face="bold"),
         panel.background = element_rect(fill = "white", colour = "grey40"),
         axis.text.y=element_blank(),
         axis.ticks.y=element_blank(),
         axis.title.x=element_text("",face="bold", size=18),
         axis.title.y=element_text(angle=360,vjust = 0.5,face="bold", size=18),
         axis.text.x=element_text(face="bold"),
         axis.title=element_text(size=12,face="bold"),
         strip.background = element_blank(),
         strip.text.x = element_blank())+scale_shape_manual("",values=c(21, 22,23,24))+
   coord_flip()

   
require(gridExtra)
g <- grid.arrange(arrangeGrob(p1, ncol=1, nrow=1),
             arrangeGrob(p2, ncol=1, nrow=1),
             arrangeGrob(p3, ncol=1, nrow=1),
             arrangeGrob(p4, ncol=1, nrow=1),
             arrangeGrob(p5, ncol=1, nrow=1),
             arrangeGrob(p6, ncol=1, nrow=1),
             arrangeGrob(p7, ncol=1, nrow=1),
             arrangeGrob(p8, ncol=1, nrow=1),
             arrangeGrob(p9, ncol=1, nrow=1),
                     heights=c(1.6,1,1,1,1,1,1,1,1)) 



ggsave(file="Figure_2.eps", g, width = 30, height = 40, units = "cm") 



###############################################################
##################Coverage Rate - Figure 3####################
##############################################################
coverage <- read.table("dataset_figure3", 
                   header = TRUE)

coverage$par <- factor(coverage$par, labels = c("hat(beta)[11]", "hat(beta)[12]", "hat(beta)[21]","hat(beta)[22]","hat(phi)[11]","hat(phi)[12]","hat(phi)[21]","hat(phi)[22]","hat(rho)"))

b <- ggplot(data = coverage, mapping = aes(y = coverage, x = n ,group = as.factor(rho),colour=as.factor(rho))) +
   geom_line(size=0.8) +xlab("Sample size")+ylab("Coverage Rate \n")+
   geom_point(size=0.8) +
   facet_wrap(~par, labeller = label_parsed,ncol = 4)+ scale_colour_grey(name = "", labels = c(expression(rho==-0.8),expression(rho==0.2),expression(rho==0.5),expression(rho==0.8)),start = 0.05, end = 0.85)+theme(axis.text = element_text(color="gray50",size = 12,face="bold"),legend.position="top",legend.text=element_text(size=18),axis.title=element_text(size=15,face="bold"),strip.text.x = element_text(size = 18), panel.background = element_rect(fill = "white", colour = "grey40"))+geom_hline(yintercept = 0.95,linetype = "dotted")+ylim(0.7,1)



ggsave(file="Figure_3.png", b, width = 25, height = 30, units = "cm") 

