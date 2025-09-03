library(foreach)
library(ggplot2)
library(ggpubr)
library(ggnewscale)


A <- readRDS("~/weihaowang/PoDAsuppsim/simulation/SE/unbalanceplnSE.rds")
simlist <- t(A[1:3,])
colind <- c(1,2,5,4,3)
colind <- c(colind, colind+5)
error=simlist[10+colind,]*100
simlist <- simlist[colind,]
simlist <- simlist*100
rownames(simlist) <- c("power.IPOD","power.LinDA","power.locom","power.ANCOMBC","power.ANCOMBC2",
                       "fdr.IPOD","fdr.LinDA","fdr.locom","fdr.ANCOMBC","fdr.ANCOMBC2")

f <- function(i) {
  index <- i
  data.frame(
    class = rep(c("Power(%)","Empirical FDR(%)"), each=5),
    method = rep(c('PoDA','LinDA','LOCOM','ANCOM-BC','ANCOM-BC2'), 2),
    value = simlist[, index],
    se = error[, index]  # 添加误差列
  )
}


power.fdr <- foreach(i=1:3, .combine=rbind) %do% f(i)


out <- data.frame(
  sig = rep(c(50,100,200), each=10),
  q = rep(c(0.5,5), each=5),
  power.fdr
)

colnames(out) <- c("sig.prob","q","class","method","value","se")  # 包含误差列
out$sig.prob <- factor(out$sig.prob, levels = c("50","100","200"))
out$class <- factor(out$class, levels = c("Power(%)","Empirical FDR(%)"))
out$method <- factor(out$method, levels = c('PoDA','LinDA','LOCOM','ANCOM-BC','ANCOM-BC2'))
out$q <- factor(out$q, levels = c("0.5","5"))
line <- data.frame(
  class = factor(c("Power(%)","Empirical FDR(%)"), levels = c("Power(%)","Empirical FDR(%)")),
  y = c(NA,5)
)
out$mu <- factor(rep("nu=0.05", nrow(out)))


ptest1 <- ggplot(out, aes(x=sig.prob, y=value, fill=method)) +  # 主数据集使用全部数据
  theme_bw() +
  geom_bar(
    data = out[out$q=="5",],
    stat="identity",
    position=position_dodge(0.7),
    width=0.6
  ) +
  geom_errorbar(
    data = out[out$q=="5",],
    aes(ymin = value - se, ymax = value + se, color = after_scale(fill)),
    position = position_dodge(0.7),
    width = 0.2
  ) +
  geom_bar(
    data = out[out$q=="0.5",],
    stat="identity",
    position=position_dodge(0.7),
    width=0.6
  ) +
  geom_errorbar(
    data = out[out$q=="0.5",],
    aes(ymin = value - se, ymax = value + se, color = after_scale(fill)),
    position = position_dodge(0.7),
    width = 0.2
  ) +
  scale_fill_manual(values = c('red','#8074C8','#00A08799','#F0C284','#E64B3599')) +
  scale_color_identity() +  # 确保颜色与填充色一致
  facet_grid(vars(class), labeller = labeller(.cols = label_both)) +
  ylim(c(0,100)) +
  geom_hline(data = line, aes(yintercept = y), linetype = "dashed") +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 6),
    legend.key.size = unit(0.3, "cm"),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 6),
    panel.spacing = unit(0.5, "lines")
  ) +
  guides(fill = guide_legend(ncol = 6)) +
  labs(y = ' ', x = expression(n))

ptest1


simlist <- t(A[4:6,])
colind <- c(1,2,5,4,3)
colind <- c(colind, colind+5)
error=simlist[10+colind,]*100
simlist <- simlist[colind,]
simlist <- simlist *100
rownames(simlist) <- c("power.IPOD","power.LinDA","power.locom","power.ANCOMBC","power.ANCOMBC2",
                       "fdr.IPOD","fdr.LinDA","fdr.locom","fdr.ANCOMBC","fdr.ANCOMBC2")







f <- function(i) {
  index <- i
  data.frame(
    class = rep(c("Power(%)","Empirical FDR(%)"), each=5),
    method = rep(c('PoDA','LinDA','LOCOM','ANCOM-BC','ANCOM-BC2'), 2),
    value = simlist[, index],
    se = error[, index]  # 添加误差列
  )
}

# 生成包含误差的数据
power.fdr <- foreach(i=1:3, .combine=rbind) %do% f(i)

# 完善数据框
out <- data.frame(
  sig = rep(c(50,100,200), each=10),
  q = rep(c(0.5,5), each=5),
  power.fdr
)

colnames(out) <- c("sig.prob","q","class","method","value","se")  # 包含误差列
out$sig.prob <- factor(out$sig.prob, levels = c("50","100","200"))
out$class <- factor(out$class, levels = c("Power(%)","Empirical FDR(%)"))
out$method <- factor(out$method, levels = c('PoDA','LinDA','LOCOM','ANCOM-BC','ANCOM-BC2'))
out$q <- factor(out$q, levels = c("0.5","5"))
line <- data.frame(
  class = factor(c("Power(%)","Empirical FDR(%)"), levels = c("Power(%)","Empirical FDR(%)")),
  y = c(NA,5)
)
out$mu <- factor(rep("nu=0.05", nrow(out)))


ptest2 <- ggplot(out, aes(x=sig.prob, y=value, fill=method)) +
  theme_bw() +
  geom_bar(
    data = out[out$q=="5",],
    stat="identity",
    position=position_dodge(0.7),
    width=0.6
  ) +
  geom_errorbar(
    data = out[out$q=="5",],
    aes(ymin = value - se, ymax = value + se, color = after_scale(fill)),
    position = position_dodge(0.7),
    width = 0.2
  ) +
  geom_bar(
    data = out[out$q=="0.5",],
    stat="identity",
    position=position_dodge(0.7),
    width=0.6
  ) +
  geom_errorbar(
    data = out[out$q=="0.5",],
    aes(ymin = value - se, ymax = value + se, color = after_scale(fill)),
    position = position_dodge(0.7),
    width = 0.2
  ) +
  scale_fill_manual(values = c('red','#8074C8','#00A08799','#F0C284','#E64B3599')) +
  scale_color_identity() +  # 确保颜色与填充色一致
  facet_grid(vars(class), labeller = labeller(.cols = label_both)) +
  ylim(c(0,100)) +
  geom_hline(data = line, aes(yintercept = y), linetype = "dashed") +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 6),
    legend.key.size = unit(0.3, "cm"),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 6),
    panel.spacing = unit(0.5, "lines")
  ) +
  guides(fill = guide_legend(ncol = 6)) +
  labs(y = ' ', x = expression(n))




simlist <- t(A[7:9,])
colind <- c(1,2,5,4,3)
colind <- c(colind, colind+5)
error=simlist[10+colind,]*100
simlist <- simlist[colind,]
simlist <- simlist *100
rownames(simlist) <- c("power.IPOD","power.LinDA","power.locom","power.ANCOMBC","power.ANCOMBC2",
                       "fdr.IPOD","fdr.LinDA","fdr.locom","fdr.ANCOMBC","fdr.ANCOMBC2")


f <- function(i) {
  index <- i
  data.frame(
    class = rep(c("Power(%)","Empirical FDR(%)"), each=5),
    method = rep(c('PoDA','LinDA','LOCOM','ANCOM-BC','ANCOM-BC2'), 2),
    value = simlist[, index],
    se = error[, index]
  )
}


power.fdr <- foreach(i=1:3, .combine=rbind) %do% f(i)


out <- data.frame(
  sig = rep(c(50,100,200), each=10),
  q = rep(c(0.5,5), each=5),
  power.fdr
)

colnames(out) <- c("sig.prob","q","class","method","value","se")
out$sig.prob <- factor(out$sig.prob, levels = c("50","100","200"))
out$class <- factor(out$class, levels = c("Power(%)","Empirical FDR(%)"))
out$method <- factor(out$method, levels = c('PoDA','LinDA','LOCOM','ANCOM-BC','ANCOM-BC2'))
out$q <- factor(out$q, levels = c("0.5","5"))
line <- data.frame(
  class = factor(c("Power(%)","Empirical FDR(%)"), levels = c("Power(%)","Empirical FDR(%)")),
  y = c(NA,5)
)
out$mu <- factor(rep("nu=0.05", nrow(out)))


ptest3 <- ggplot(out, aes(x=sig.prob, y=value, fill=method)) +  # 主数据集使用全部数据
  theme_bw() +
  geom_bar(
    data = out[out$q=="5",],
    stat="identity",
    position=position_dodge(0.7),
    width=0.6
  ) +
  geom_errorbar(
    data = out[out$q=="5",],
    aes(ymin = value - se, ymax = value + se, color = after_scale(fill)),
    position = position_dodge(0.7),
    width = 0.2
  ) +
  geom_bar(
    data = out[out$q=="0.5",],
    stat="identity",
    position=position_dodge(0.7),
    width=0.6
  ) +
  geom_errorbar(
    data = out[out$q=="0.5",],
    aes(ymin = value - se, ymax = value + se, color = after_scale(fill)),
    position = position_dodge(0.7),
    width = 0.2
  ) +
  scale_fill_manual(values = c('red','#8074C8','#00A08799','#F0C284','#E64B3599')) +
  scale_color_identity() +
  facet_grid(vars(class), labeller = labeller(.cols = label_both)) +
  ylim(c(0,100)) +
  geom_hline(data = line, aes(yintercept = y), linetype = "dashed") +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 6),
    legend.key.size = unit(0.3, "cm"),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 6),
    panel.spacing = unit(0.5, "lines")
  ) +
  guides(fill = guide_legend(ncol = 6)) +
  labs(y = ' ', x = expression(n))

ptest1 <- ptest1 + theme(plot.margin = margin(1, 1, 1, 1))
ptest2 <- ptest2 + theme(plot.margin = margin(1, 1, 1, 1))
ptest3 <- ptest3 + theme(plot.margin = margin(1, 1, 1, 1))



# 调整每个子图的样式
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
    strip.text = element_text(size = 12),
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
  bottom = text_grob(expression(n), size = 12) # 设置横坐标标签
)

ggsave("~/weihaowang/PoDAsuppsim/errorbar画图/A12unbalancepln.png", p.comb, width = 4*3.5, height =2.5*3.5, limitsize = FALSE ,bg = "white")

