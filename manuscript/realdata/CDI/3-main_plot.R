library(ggstar)
library(ggpubr)
library(patchwork)
library(ggplot2)


df <- expand.grid(
  Pseudocount = c(as.character(seq(0.1, 1, 0.1)),"adaImpute","multLN","multRepl"),
  method = c(" Clostridium_XI", " Clostridium_XlVa"),
  dataset = c("PoDA", "LinDA")
)

df$status <- c(# poda
               rep("detected",13), # 15 Clostridium_XI
               "detected",rep("undetected",12), # 17 Clostridium_XlVa
               #linda
               rep("detected",13),# 15 Clostridium_XI
               rep("detected",5),rep("undetected",5),"detected","undetected","detected") # 17 Clostridium_XlVa



df$method <- factor(df$method, levels = c(" Clostridium_XI", " Clostridium_XlVa"))
df1=df[1:26,]
df2=df[27:52,]



df3 <- expand.grid(
  Pseudocount = c(as.character(seq(0.1, 1, 0.1)),"adaImpute","multLN","multRepl"),
  method = c("Alistipes","Anaerofustis","Clostridium_XI","Collinsella","Eubacterium","Lactococcus","Oscillibacter",
             "Propionibacterium","Ralstonia","Sphingomonas"),
  dataset = c("ANCOM-BC")
)

# imp 1 2 3 4 5 6 7 8 9 10 LN Repl
df3$status <- c(rep("detected",3),rep("undetected",9), "undetected", # 3 Alistipes
                rep("detected",2),rep("undetected",10),"undetected", # 5 Anaerofustis
                rep("detected",13), # 15 Clostridium_XI
                rep("detected",3),rep("undetected",9),"undetected", # 20 Collinsella
                rep("detected",3),rep("undetected",9), "undetected",# 26 Eubacterium
                rep("detected",4),rep("undetected",8),"undetected", # 32 Lactococcus
                rep("detected",3),rep("undetected",9),"undetected", # 34 Oscillibacter
                rep("detected",3),rep("undetected",9),"undetected", # 40 Propionibacterium
                rep("detected",3),rep("undetected",9), "undetected",# 41 Ralstonia
                rep("detected",3),rep("undetected",9),"undetected") # 45 Sphingomonas







p6 = ggplot() +
  geom_star(data = df1, aes(x = Pseudocount, y = method, fill = status), size = 3,
            starshape = 15) +
  theme_bw() +
  theme(panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y.right = element_text(angle = 0,
                                         hjust = 0,
                                         vjust = 0.5),
        plot.caption = element_text(hjust = 0.5),
        legend.title = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10)) +
  theme(panel.grid = element_blank()) +
  labs(y = expression("PoDA")) + labs(x = "Zero-replacement strategy") +
  scale_fill_manual(values = c("#2a83a2", "#ffe8d560")) +
  scale_x_discrete(
    limits = c("0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7",
               "0.8", "0.9", "1", "adaImpute","multLN", "multRepl")
  )
p6



p7 <- ggplot() +
  geom_star(data = df2, aes(x = Pseudocount, y = method, fill = status),
            size = 3, starshape = 15) +
  theme_bw() +
  theme(
    panel.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 10),
    legend.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 10),
    axis.text.y.right = element_text(
      angle = 0, hjust = 0, vjust = 0.5, size = 12
    ),
    plot.caption = element_text(hjust = 0.5),
    panel.grid = element_blank()
  ) +
  labs(y = expression("LinDA")) +
  scale_x_discrete(
    limits = c("0.1", "0.2", "0.3", "0.4", "0.5", "0.6",
               "0.7", "0.8", "0.9", "1", "adaImpute","multLN", "multRepl")
  ) +
  scale_fill_manual(values = c("#2a83a2", "#ffe8d560"))

p7


p71 <- ggplot() +
  geom_star(data = df3, aes(x = Pseudocount, y = method, fill = status),
            size = 3, starshape = 15) +
  theme_bw() +
  theme(
    panel.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(hjust = 0.5,size = 10),
    axis.ticks = element_blank(),
    legend.title = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 10),
    axis.text.y.right = element_text(angle = 0, hjust = 0, vjust = 0.5, size = 8),
    panel.grid = element_blank()
  ) +
  labs(y = "ANCOM-BC") +
  scale_x_discrete(
    limits = c("0.1", "0.2", "0.3", "0.4", "0.5", "0.6",
               "0.7", "0.8", "0.9", "1", "adaImpute","multLN", "multRepl")
  ) +
  scale_fill_manual(values = c("#2a83a2", "#ffe8d560"))

p71


ggarrange(p7,p6,ncol=1,heights  = c(1,1),common.legend = TRUE,legend="top")

p8=ggarrange(p71,p7,p6,ncol=1,heights  = c(5, 1,1.45),common.legend = TRUE,legend="top")
p8
p8_wrapped <- wrap_elements(p8)

p_fin11 <- p_fin + coord_cartesian(xlim = c(-17, 19.6), ylim = c(-10, 5))
p_fin_wrapped <- wrap_elements( p_fin11)

p_fin_wrapped + p8_wrapped +
  plot_layout(ncol = 1) +plot_annotation(tag_levels = "a")



final_plot <-p_fin_wrapped + p8_wrapped +
  plot_layout(ncol = 1,heights = c(1, 1.5)) +plot_annotation(tag_levels = "a")

final_plot

# ggsave("~/manuscript/realdata/Plots/cdi_main.png", final_plot, width = 9.5, height = 9.5,bg="white")
