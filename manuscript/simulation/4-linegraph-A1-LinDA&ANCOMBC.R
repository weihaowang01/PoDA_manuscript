library(ggplot2)
library(foreach)
library(ggpubr)
library(ggnewscale)
library(gridExtra)
library(reshape2)
aimp <- readRDS("~/manuscript/simulation/Results/zero/A1_pseudo/A1imputation.ds")
a1 <- readRDS("~/manuscript/simulation/Results/zero/A1_pseudo/A1pseudo1.rds")
a2 <- readRDS("~/manuscript/simulation/Results/zero/A1_pseudo/A1pseudo2.rds")
a3 <- readRDS("~/manuscript/simulation/Results/zero/A1_pseudo/A1pseudo3.rds")
a4 <- readRDS("~/manuscript/simulation/Results/zero/A1_pseudo/A1pseudo4.rds")
a5 <- readRDS("~/manuscript/simulation/Results/zero/A1_pseudo/A1pseudo5.rds")
a6 <- readRDS("~/manuscript/simulation/Results/zero/A1_pseudo/A1pseudo6.rds")
a7 <- readRDS("~/manuscript/simulation/Results/zero/A1_pseudo/A1pseudo7.rds")
a8 <- readRDS("~/manuscript/simulation/Results/zero/A1_pseudo/A1pseudo8.rds")
a9 <- readRDS("~/manuscript/simulation/Results/zero/A1_pseudo/A1pseudo9.rds")
a10 <- readRDS("~/manuscript/simulation/Results/zero/A1_pseudo/A1pseudo10.rds")
aLN<-readRDS("~/manuscript/simulation/Results/zero/A1_pseudo/A1multLN.rds")
aRepl<-readRDS("~/manuscript/simulation/Results/zero/A1_pseudo/A1multirepl.rds")
# #LinDA
# pseu_poda_power=cbind(a1[,2],a2[,2],a3[,2],a4[,2],a5[,2],aimp[,2],
#                       a6[,2],a7[,2],a8[,2],a9[,2],a10[,2],aLN[,2],aRepl[,2])
# pseu_poda_fdr=cbind(a1[,7],a2[,7],a3[,7],a4[,7],a5[,7],aimp[,7],
#                     a6[,7],a7[,7],a8[,7],a9[,7],a10[,7],aLN[,7],aRepl[,7])

# #ACBC
pseu_poda_power=cbind(a1[,4],a2[,4],a3[,4],a4[,4],a5[,4],
                      a6[,4],a7[,4],a8[,4],a9[,4],a10[,4],aimp[,4],aLN[,4],aRepl[,4])
pseu_poda_fdr=cbind(a1[,9],a2[,9],a3[,9],a4[,9],a5[,9],
                    a6[,9],a7[,9],a8[,9],a9[,9],a10[,9],aimp[,9],aLN[,9],aRepl[,9])

pseu_poda=cbind(pseu_poda_power,pseu_poda_fdr)

A=pseu_poda
simlist=t(A[1:3,])
simlist=simlist*100
rownames(simlist)=c("power.p1","power.p2","power.p3","power.p4","power.p5","power.p6","power.p7","power.p8","power.p9","power.p10",
                    "power.imp","power.LN","power.Repl",
                    "fdr.p1","fdr.p2","fdr.p3","fdr.p4","fdr.p5","fdr.p6","fdr.p7","fdr.p8","fdr.p9","fdr.p10",
                    "fdr.imp","fdr.LN","fdr.Repl")



A1=A[1:3,]*100
rownames(A1) <- c("n=50", "n=100", "n=200")
colnames(A1) <-rep(c(as.character(seq(0.1, 1, by = 0.1)),"adaImpute","multLN","multRepl"), 2)

power <- A1[1:3, 1:13]
fdr   <- A1[1:3, 14:26]


df_power <- melt(power)
df_fdr   <- melt(fdr)


df_power$class <- "Power(%)"
df_fdr$class   <- "Empirical FDR(%)"


df_all <- rbind(df_power, df_fdr)
colnames(df_all) <- c("method", "setting", "value", "class")


df_all$method <- factor(df_all$method, levels = unique(df_all$method))
df_all$setting <- factor(df_all$setting, levels = c(as.character(seq(0.1, 1, by = 0.1)),"adaImpute","multLN","multRepl"))
df_all$class <- factor(df_all$class, levels = c("Power(%)", "Empirical FDR(%)"))


line_df <- data.frame(class = "Empirical FDR(%)", yintercept = 5)
line_df$class <- factor(line_df$class, levels = levels(df_all$class))


p1=ggplot(df_all, aes(x = setting, y = value, color = method, group = method)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  facet_grid(rows = vars(class), labeller = labeller(.cols = label_both)) +
  ylim(0, NA) +
  geom_hline(data = line_df, aes(yintercept = yintercept),
             linetype = "dashed", color = "red", inherit.aes = FALSE) +
  theme_bw() +
  labs(x = "Setting", y = NULL, title = NULL) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    strip.text = element_text(size = 10),
    legend.title = element_blank()
  )+
  theme(
    legend.position = "top",
    axis.title.x = element_blank(),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.5, "cm"),
    axis.text = element_text(size = 6),
    plot.margin = margin(1, 1, 1, 1),
    strip.background = element_blank(),
    strip.text = element_blank()
  )


p1




A1=A[4:6,]*100
rownames(A1) <- c("n=50", "n=100", "n=200")


colnames(A1) <-rep(c(as.character(seq(0.1, 1, by = 0.1)),"adaImpute","multLN","multRepl"), 2)

power <- A1[1:3, 1:13]
fdr   <- A1[1:3, 14:26]

df_power <- melt(power)
df_fdr   <- melt(fdr)


df_power$class <- "Power(%)"
df_fdr$class   <- "Empirical FDR(%)"



df_all <- rbind(df_power, df_fdr)
colnames(df_all) <- c("method", "setting", "value", "class")


df_all$method <- factor(df_all$method, levels = unique(df_all$method))
df_all$setting <- factor(df_all$setting, levels = c(as.character(seq(0.1, 1, by = 0.1)),"adaImpute","multLN","multRepl"))
df_all$class <- factor(df_all$class, levels = c("Power(%)", "Empirical FDR(%)"))


line_df <- data.frame(class = "Empirical FDR(%)", yintercept = 5)
line_df$class <- factor(line_df$class, levels = levels(df_all$class))


p2=ggplot(df_all, aes(x = setting, y = value, color = method, group = method)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  facet_grid(rows = vars(class), labeller = labeller(.cols = label_both)) +
  ylim(0, NA) +
  geom_hline(data = line_df, aes(yintercept = yintercept),
             linetype = "dashed", color = "red", inherit.aes = FALSE) +
  theme_bw() +
  labs(x = "Setting", y = NULL, title = NULL) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    strip.text = element_text(size = 10),
    legend.title = element_blank()
  )+
  theme(
    axis.text = element_text(size = 6),
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

p2


A1=A[7:9,]*100
rownames(A1) <- c("n=50", "n=100", "n=200")


colnames(A1) <-rep(c(as.character(seq(0.1, 1, by = 0.1)),"adaImpute ","multLN","multRepl"), 2)

power <- A1[1:3, 1:13]
fdr   <- A1[1:3, 14:26]

df_power <- melt(power)
df_fdr   <- melt(fdr)


df_power$class <- "Power(%)"
df_fdr$class   <- "Empirical FDR(%)"


df_all <- rbind(df_power, df_fdr)
colnames(df_all) <- c("method", "setting", "value", "class")


df_all$method <- factor(df_all$method, levels = unique(df_all$method))
df_all$setting <- factor(df_all$setting, levels = c(as.character(seq(0.1, 1, by = 0.1)),"adaImpute ","multLN","multRepl"))
df_all$class <- factor(df_all$class, levels = c("Power(%)", "Empirical FDR(%)"))

line_df <- data.frame(class = "Empirical FDR(%)", yintercept = 5)
line_df$class <- factor(line_df$class, levels = levels(df_all$class))


p3=ggplot(df_all, aes(x = setting, y = value, color = method, group = method)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  facet_grid(rows = vars(class), labeller = labeller(.cols = label_both)) +
  ylim(0, NA) +
  geom_hline(data = line_df, aes(yintercept = yintercept),
             linetype = "dashed", color = "red", inherit.aes = FALSE) +
  theme_bw() +
  labs(x = "Setting", y = NULL, title = NULL) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    strip.text = element_text(size = 13),
    legend.title = element_blank()
  )+
  theme(
    axis.title.x = element_blank(),
    axis.text = element_text(size = 6) ,
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.5, "cm"),
    plot.margin = margin(1, 1, 1, 1),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank()
  )



p3



ptest1 <- p1 + theme(plot.margin = margin(1, 1, 1, 1))
ptest2 <- p2 + theme(plot.margin = margin(1, 1, 1, 1))
ptest3 <- p3 + theme(plot.margin = margin(1, 1, 1, 1))

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
  bottom = text_grob("Zero-replacement strategy", size = 12)
)



# ggsave("~/manuscript/simulation/Plots/zero/A1new_acbc.png", p.comb, width = 4*4, height =2.5*4, limitsize = FALSE ,bg = "white")
