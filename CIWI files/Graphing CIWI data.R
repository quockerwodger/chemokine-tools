library(tidyverse)
library(xml2)

# Read in data

variable_table <- as_tibble(read.csv("CKBP CIWI for thesis/Calculated means per variable.csv"))
position_table <- as_tibble(read.csv("CKBP CIWI for thesis/Calculated means per consensus position.csv"))
raw_scores <- as_tibble(read.csv("CKBP CIWI for thesis/All raw scores.csv"))

# Reformat data

variable_table <- variable_table %>%
  rename("Chemokine" = "X") %>%
  mutate(
    Chemokine = factor(Chemokine),
    Binder = factor(Binder)
  )
position_table <- position_table %>%
  rename("Chemokine" = "X") %>%
  mutate(
    Chemokine = factor(Chemokine),
    Binder = factor(Binder)
  )
raw_scores <- raw_scores %>%
  rename("Variable" = "X") %>%
  mutate(
    Variable = factor(Variable),
    Binder = factor(Binder)
  )

# Convert to graphable tables
position_longer <- pivot_longer(position_table, 2:42, names_to = "Position", values_to = "Value")
variable_longer <- pivot_longer(variable_table, 2:4, names_to = "Variable", values_to = "Value")
raw_longer <- pivot_longer(raw_scores, 2:42, names_to = "Position", values_to = "Value")
raw_longer <- raw_longer %>%
  mutate(
    Position = factor(Position)
  )

# Plotting positional means
position_boxplot <- ggplot(position_longer, aes(y = Value, color = Binder))+
  geom_boxplot()+
  facet_wrap(vars(Position), scales = "free_y")
position_boxplot
position_density <- ggplot(position_longer, aes(x = Value, color = Binder))+
  geom_density()+
  facet_wrap(vars(Position), scales = "free_y")
position_density
position_histogram <- ggplot(position_longer, aes(x = Value, y = after_stat((density/10)*100), fill = Binder))+
  geom_histogram(binwidth = 0.1, alpha = 0.5, position = "identity")+
  facet_wrap(vars(Position), scales = "free_y")
position_histogram
position_bar <- ggplot(position_longer, aes(x = Value, fill = Binder))+
  geom_bar(position = "dodge")+
  facet_wrap(vars(Position))
position_bar
position_col <- ggplot(position_longer, aes(x = Position, y = Value, fill = Binder))+
  geom_col(position = "dodge")
position_col

# Plotting variable means
variable_density <- ggplot(variable_longer, aes(x = Value, color = Binder))+
  geom_density()+
  facet_wrap(vars(Variable), scales = "free")
variable_density
variable_boxplot <- ggplot(variable_longer, aes(y = Value, color = Binder))+
  geom_boxplot()+
  facet_wrap(vars(Variable), scales = "free")
variable_boxplot
variable_histogram <- ggplot(variable_longer, aes(x = Value, y = after_stat((density/10)*100), fill = Binder))+
  geom_histogram(binwidth = 0.1, alpha = 0.5, position = "identity")+
  facet_wrap(vars(Variable), scales = "free_y")
variable_histogram
variable_bar <- ggplot(variable_longer, aes(x = Value, fill = Binder))+
  geom_bar(position = "dodge")+
  facet_wrap(vars(Variable))
variable_bar
variable_col <- ggplot(variable_longer, aes(x = Variable, y = Value, fill = Binder))+
  geom_col(position = "dodge")
variable_col

# Creating raw score summary table
raw_summary <- raw_longer %>%
  group_by(Variable, Position, Binder) %>%
  summarise(mean_value = mean(Value),
            value_se = sd(Value)/sqrt(n()),
            n_chemokines = n())%>%
  ungroup()
raw_summary$Position <- factor(raw_summary$Position, levels = c("X0", "X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10", "X11", "X12", "X13", "X14", "X15", "X16", "X17", "X18", "X19", "X20", "X21", "X22", "X23", "X24", "X25", "X26", "X27", "X28", "X29", "X30", "X31", "X32", "X33", "X34", "X35", "X36", "X37", "X38", "X39", "X40"))
raw_summary <- raw_summary %>%
  mutate(
    Position = as.numeric(Position)
  )

# Plotting raw scores
raw_density <- ggplot(raw_longer, aes(x = Value, color = Binder))+
  geom_density()+
  facet_grid(cols = vars(Position), rows = vars(Variable), scales = "free")
raw_density
raw_boxplot <- ggplot(raw_longer, aes(x = Value, color = Binder))+
  geom_boxplot()+
  facet_grid(cols = vars(Position), rows = vars(Variable), scales = "free")
raw_boxplot
raw_histogram <- ggplot(raw_longer, aes(x = Value, y = after_stat((count/10)*100), fill = Binder))+
  geom_histogram(binwidth = 1.0, alpha = 0.5, position = "identity")+
  facet_grid(cols = vars(Position), rows = vars(Variable), scales = "free")
raw_histogram
raw_bar <- ggplot(raw_longer, aes (x = Value, fill = Binder))+
  geom_bar(position = "dodge")+
  facet_grid(rows = vars(Variable), scales = "free")
raw_bar
raw_col <- ggplot(raw_summary, aes(x = Position, y = mean_value, fill = Binder))+
  geom_col(position=position_dodge(0.8), show.legend = TRUE, linewidth = 0.3)+
  labs(x = "Residue", y = "Mean match")+
  #geom_errorbar(aes(x=res, ymax=Involved_mean+Involved_se, ymin = Involved_mean-Involved_se, colour = "black"), position=position_dodge(0.8))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), axis.text = element_text(size = 15), axis.title = element_text(size = 20))+
  scale_fill_manual(values = c("#0dbbbf", "#dd3224"))+
  scale_colour_manual(values = c("black"))+
  scale_y_continuous(expand = c(0,0), limits = c(0, 1.05))+
  scale_x_continuous(expand = c(0,0), limits = c(0, 42))+
  theme(panel.background = element_blank(),
        panel.grid.major.y = element_line(colour = "black", linewidth = 0.5, linetype = "dashed"),
        panel.border = element_rect(colour = "black", linewidth =0.75, fill = NA),
        text = element_text(size = 11, colour = "black"),
        axis.text = element_text(size = 11, colour = "black"),
        axis.ticks.x= element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank()
  )+
  geom_vline(xintercept = c(0.5, 1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5,14.5,15.5,16.5,17.5,18.5,19.5,20.5,21.5,22.5,23.5,24.5,25.5,26.5,27.5,28.5,29.5,30.5,31.5,32.5,33.5,34.5,35.5,36.5,37.5,38.5,39.5,40.5, 41.5), linetype = "dashed", colour = "black", linewidth = 0.5)+

  geom_errorbar(aes(ymax = mean_value+value_se, ymin = mean_value-value_se), position = "dodge")+
  #facet_wrap(vars(Variable), scales = "free_x")
  facet_wrap(nrow = 3, vars(Variable), scales = "free_x")
raw_col

#KD plots
raw_longer$Position <- factor(raw_longer$Position, levels = c("X0", "X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10", "X11", "X12", "X13", "X14", "X15", "X16", "X17", "X18", "X19", "X20", "X21", "X22", "X23", "X24", "X25", "X26", "X27", "X28", "X29", "X30", "X31"))
raw_longer <- raw_longer %>%
  mutate(
    Position = as.numeric(Position),
    Value = factor(Value)
  )
KD_density <- ggplot(raw_longer, aes(x = KD, fill = Value))+
  geom_density(alpha = 0.5)+
  facet_grid(cols = vars(Position), rows = vars(Variable))
KD_density
KD_boxplot <- ggplot(raw_longer, aes (y = KD, colour = Value))+
  geom_boxplot()+
  #geom_point(aes(x = Value), position = "jitter")+
  facet_grid(cols = vars(Position), rows = vars(Variable))
KD_boxplot
KD_scatter <- ggplot(raw_longer, aes(y = KD, x = Position, colour = Value))+
  geom_point(position = "jitter")+
  facet_grid(rows = vars(Variable), scales = 'free')
KD_scatter
KD_summary <- raw_longer %>%
  group_by(Variable, Position, Value)%>%
  summarise(mean_KD = mean(KD),
            KD_se = sd(KD)/sqrt(n()),
            n_chemokines = n())%>%
  ungroup()
KD_column <- ggplot(KD_summary, aes(y = mean_KD, x = Position, fill = Value))+
  geom_col(position = "dodge")+
  geom_errorbar(aes(ymax = mean_KD+KD_se, ymin = mean_KD-KD_se), position = "dodge")+
  facet_grid(rows = vars(Variable), scales = 'free')
KD_column
KD_line <- ggplot(KD_summary, aes(y = mean_KD, x = Position, colour = Value))+
  geom_line(linewidth = 1)+
  facet_grid(rows = vars(Variable), scales = 'free')
KD_line
KD_wrap <- ggplot(KD_summary, aes(y = mean_KD, x = Position, fill = Value))+
  geom_col(position = "dodge")+
  geom_errorbar(aes(ymax = mean_KD+KD_se, ymin = mean_KD-KD_se), position = "dodge")+
  facet_wrap(vars(Variable), scales = 'free')
ggsave("CKBP CIWI for thesis/KD density.png", KD_density, height = 30, width = 40, units = "cm", dpi = 500)
ggsave("CKBP CIWI for thesis/KD boxplot.png", KD_boxplot, height = 30, width = 40, units = "cm", dpi = 500)
ggsave("CKBP CIWI for thesis/KD scatter.png", KD_scatter, height = 30, width = 40, units = "cm", dpi = 500)
ggsave("CKBP CIWI for thesis/KD column.png", KD_column, height = 30, width = 40, units = "cm", dpi = 500)
ggsave("CKBP CIWI for thesis/KD line.png", KD_line, height = 30, width = 40, units = "cm", dpi = 500)
ggsave("CKBP CIWI for thesis/KD column wrapped.png", KD_wrap, height = 30, width = 40, units = "cm", dpi = 500)

# Save plots for means

ggsave("CKBP CIWI for thesis/Position boxplot.png", position_boxplot, height = 20, width = 20, units = "cm", dpi = 500)
ggsave("CKBP CIWI for thesis/Position density.png", position_density, height = 20, width = 20, units = "cm", dpi = 500)
ggsave("CKBP CIWI for thesis/Position histogram.png", position_histogram, height = 20, width = 20, units = "cm", dpi = 500)
ggsave("CKBP CIWI for thesis/Position bar chart.png", position_bar, height = 20, width = 20, units = "cm", dpi = 500)
ggsave("CKBP CIWI for thesis/Position columns.png", position_col, height = 20, width = 20, units = "cm", dpi = 500)
ggsave("CKBP CIWI for thesis/Variable boxplot.png", variable_boxplot, height = 20, width = 20, units = "cm", dpi = 500)
ggsave("CKBP CIWI for thesis/Variable density.png", variable_density, height = 20, width = 20, units = "cm", dpi = 500)
ggsave("CKBP CIWI for thesis/Variable histogram.png", variable_histogram, height = 20, width = 20, units = "cm", dpi = 500)
ggsave("CKBP CIWI for thesis/Variable bar chart.png", variable_bar, height = 20, width = 20, units = "cm", dpi = 500)
ggsave("CKBP CIWI for thesis/Variable columns.png", variable_col, height = 20, width = 20, units = "cm", dpi = 500)
ggsave("CKBP CIWI for thesis/Raw boxplot.png", raw_boxplot, height = 20, width = 20, units = "cm", dpi = 500)
ggsave("CKBP CIWI for thesis/Raw density.png", raw_density, height = 20, width = 20, units = "cm", dpi = 500)
ggsave("CKBP CIWI for thesis/Raw histogram.png", raw_histogram, height = 20, width = 20, units = "cm", dpi = 500)
ggsave("CKBP CIWI for thesis/Raw bar chart.png", raw_bar, height = 20, width = 20, units = "cm", dpi = 500)
ggsave("CKBP CIWI for thesis/Raw columns.png", raw_col, height = 25, width = 25, units = "cm", dpi = 500)

