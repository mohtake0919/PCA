
library(tidyverse)
library(readxl)
library(FactoMineR)
library(factoextra)


circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(tibble(x = xx, y = yy))
}


d3 = tibble(sheet = c("Raw data_Growth period",
                 "Raw data_Maturation period",
                 "Raw data_full year")) %>%
  mutate(data = map(sheet, read_xlsx, path = "主成分分析_190924.xlsx"))

d3 = d3 %>% separate(sheet, letters[1:3]) %>%
  select(group = c, data) %>%
  unnest()

library(vegan)

X = d3 %>% select(Vmax, Km, `Vmax/Km`, `P uptake`, `P demand_max`, `P demand_in situ`)

manova(cbind(Vmax, Km, `Vmax/Km`, `P uptake`) ~
         生長速度 + 全リン + PO4濃度 + NO3濃度 + PAR + 水温 + WW + group, data = d3) %>%
  summary()


tmp1 = get_pca_var(d3out) %>%  pluck("contrib") %>% as_tibble(rownames = "Variable") %>%
  gather("Dim", "Contribution", -Variable)
tmp2 = get_pca_var(d3out) %>%  pluck("cor") %>% as_tibble(rownames = "Variable") %>%
  gather("Dim", "Correlation", -Variable)

N = length(unique(tmp1$Variable))

d3_x = full_join(tmp1, tmp2)
d3_x = d3_x %>% mutate(C = ifelse(Contribution > 100 / N,
                           "Above average contribution",
                           "Below average contribution"))

p1 = d3_x %>% filter(str_detect(Dim, "Dim.1")) %>%
  ggplot(aes(x = reorder(Variable, -Correlation), y = Correlation)) +
  geom_line(aes(group = 1), linetype = "dotted") +
  geom_point(aes(shape = C, fill = C), size = 3) +
  geom_hline(yintercept= 0, linetype = "dashed") +
  scale_y_continuous("Dimension 1", limits = c(-1,1)) +
  scale_shape_manual(values = c(19, 21)) +
  scale_fill_manual(values = c("black", "white")) +
  theme(legend.background=element_blank(),
        legend.position=c(1,1),
        legend.justification=c(1,1),
        legend.title=element_blank())

p2 = d3_x %>% filter(str_detect(Dim, "Dim.2")) %>%
  ggplot(aes(x = reorder(Variable, -Correlation), y = Correlation)) +
  geom_line(aes(group = 1), linetype = "dotted") +
  geom_point(aes(shape = C, fill = C), size = 3) +
  geom_hline(yintercept= 0, linetype = "dashed") +
  scale_y_continuous("Dimension 1", limits = c(-1,1)) +
  scale_shape_manual(values = c(19, 21)) +
  scale_fill_manual(values = c("black", "white")) +
  theme(legend.background=element_blank(),
        legend.position=c(1,1),
        legend.justification=c(1,1),
        legend.title=element_blank())

p3 = d3_x %>% filter(str_detect(Dim, "Dim.3")) %>%
  ggplot(aes(x = reorder(Variable, -Correlation), y = Correlation)) +
  geom_line(aes(group = 1), linetype = "dotted") +
  geom_point(aes(shape = C, fill = C), size = 3) +
  geom_hline(yintercept= 0, linetype = "dashed") +
  scale_y_continuous("Dimension 1", limits = c(-1,1)) +
  scale_shape_manual(values = c(19, 21)) +
  scale_fill_manual(values = c("black", "white")) +
  theme(legend.background=element_blank(),
        legend.position=c(1,1),
        legend.justification=c(1,1),
        legend.title=element_blank())
p4 = d3_x %>% filter(str_detect(Dim, "Dim.4")) %>%
  ggplot(aes(x = reorder(Variable, -Correlation), y = Correlation)) +
  geom_line(aes(group = 1), linetype = "dotted") +
  geom_point(aes(shape = C, fill = C), size = 3) +
  geom_hline(yintercept= 0, linetype = "dashed") +
  scale_y_continuous("Dimension 1", limits = c(-1,1)) +
  scale_shape_manual(values = c(19, 21)) +
  scale_fill_manual(values = c("black", "white")) +
  theme(legend.background=element_blank(),
        legend.position=c(1,1),
        legend.justification=c(1,1),
        legend.title=element_blank())

cowplot::plot_grid(p1,p2,p3, p4)



pcacoord = get_pca_var(d3out) %>%  pluck("coord") %>% as_tibble(rownames = "Variable")
pcacoord = pcacoord %>%
  mutate(d1 = 1.1*sqrt(Dim.1^2 + Dim.2^2) * sin(atan2(Dim.1, Dim.2)),
         d2 = 1.1*sqrt(Dim.1^2 + Dim.2^2) * cos(atan2(Dim.1, Dim.2)),
         d3 = 1.1*sqrt(Dim.3^2 + Dim.4^2) * sin(atan2(Dim.3, Dim.4)),
         d4 = 1.1*sqrt(Dim.3^2 + Dim.4^2) * cos(atan2(Dim.3, Dim.4)),
         a1 = atan2(Dim.2, Dim.1),
         a2 = atan2(Dim.4, Dim.3)) %>%
  mutate(adj1 = ifelse(abs(a1) > pi/2, a1 - pi, a1),
         adj2 = ifelse(abs(a2) > pi/2, a2 - pi, a2))

pcacoord %>% select(Variable, adj)

eval = d3out %>% get_eig()

xlab = str_glue("Dim. 1 ({round(eval[1,2], 1)}%)")
xlab2 = str_glue("Dim. 3 ({round(eval[3,2], 1)}%)")
ylab = str_glue("Dim. 2 ({round(eval[2,2], 1)}%)")
ylab2 = str_glue("Dim. 4 ({round(eval[4,2], 1)}%)")


p1 = ggplot(pcacoord) +
  geom_vline(xintercept=0, linetype = "dashed") +
  geom_hline(yintercept=0, linetype = "dashed") +
  geom_segment(aes(x = 0, y = 0,
                   xend = Dim.1, yend = Dim.2),
               arrow=arrow(angle = 20, length = unit(5, "mm"), type="closed")) +
  geom_text(aes(x = d1, y = d2, label = Variable,
                angle = 360*adj1/(2*pi)),
            hjust = 1,
            data = pcacoord %>% filter(abs(pcacoord$a1) > pi/2)) +
  geom_text(aes(x = d1, y = d2, label = Variable,
                angle = 360*adj1/(2*pi)),
            hjust = 0,
            data = pcacoord %>% filter(abs(pcacoord$a1) < pi/2)) +
  geom_path(aes(x = x, y = y),
            data = circleFun(diameter = 2)) +
  scale_x_continuous(xlab, limits = c(-1.2, 1.2)) +
  scale_y_continuous(ylab, limits = c(-1.2, 1.2)) +
  coord_equal()

p2 = ggplot(pcacoord) +
  geom_vline(xintercept=0, linetype = "dashed") +
  geom_hline(yintercept=0, linetype = "dashed") +
  geom_segment(aes(x = 0, y = 0,
                   xend = Dim.3, yend = Dim.4),
               arrow=arrow(angle = 20, length = unit(5, "mm"), type="closed")) +
  geom_text(aes(x = d3, y = d4, label = Variable,
                angle = 360*adj2/(2*pi)),
            hjust = 1,
            data = pcacoord %>% filter(abs(pcacoord$a2) > pi/2)) +
  geom_text(aes(x = d3, y = d4, label = Variable,
                angle = 360*adj2/(2*pi)),
            hjust = 0,
            data = pcacoord %>% filter(abs(pcacoord$a2) < pi/2)) +
  geom_path(aes(x = x, y = y),
            data = circleFun(diameter = 2)) +
  scale_x_continuous(xlab2, limits = c(-1.2, 1.2)) +
  scale_y_continuous(ylab2, limits = c(-1.2, 1.2)) +
  coord_equal()
cowplot::plot_grid(p1, p2)


