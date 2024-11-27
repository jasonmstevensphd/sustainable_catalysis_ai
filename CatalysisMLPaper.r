library(tidyverse)
library(readxl)

df <- read_xlsx('2024-08-02_CatalysisML.xlsx')

plot <- ggplot(df, aes(x = reorder(Ligand, `Predicted Yield`, mean), y = `Predicted Yield`))+
  geom_point(aes(color = Substrate, shape = Solvent), size = 2, alpha = 0.8)+
  #scale_color_viridis_d()+
  ggtitle(paste("Predicted Reaction Yield for Substrates 1a-c"))+
  labs(x = paste("Ligand"), y = paste("Predicted Yield"))+
  theme(plot.background = element_rect(fill="white"))+
  theme(panel.background = element_rect(fill="white", colour="grey50"))+
  theme(plot.title = element_text(face = "bold", 
                                  size = 18,
                                  color="#BE2BBB"))+
  theme(axis.title = element_text(face = "bold", size = 16))+
  theme(axis.text.x = element_text(vjust = 1,
                                   hjust = 1,
                                   size=12,
                                   angle = 60))+
  theme(axis.text.y = element_text(size=12))

plot #

avg_yield_sub <- df %>% group_by(Substrate) %>% summarise(avg_yield = mean(`Predicted Yield`))
avg_yield_solv <- df %>% group_by(Solvent) %>% summarise(avg_yield = mean(`Predicted Yield`))
avg_yield_Ligand <- df %>% group_by(Ligand) %>% summarise(avg_yield = mean(`Predicted Yield`)) %>% arrange(desc(avg_yield))



df_1 <- read_xlsx('NiBorylationData.xlsx')
df_1 <- df_1 %>% filter(!is.na(AP_PDT))

plot_1 <- ggplot(df_1, aes(x = AP_PDT))+
  geom_histogram(aes(fill = Ligand), color = "black", bins = 10)+
  #stat_bin(binwidth = 10, geom = "text", aes(label = ..count..), vjust = -0.5)+
  #scale_fill_viridis_d()+
  ggtitle(paste("Ligand performance distribution for the HTE study"))+
  labs(x = paste("Area Percent Range"), y = paste("Counts"))+
  #scale_x_continuous(breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100), limits = c(-0, 100))+
  theme(plot.background = element_rect(fill="white"))+
  theme(panel.background = element_rect(fill="white", colour="grey50"))+
  theme(plot.title = element_text(face = "bold", 
                                  size = 18,
                                  color="#BE2BBB"))+
  theme(axis.title = element_text(face = "bold", size = 16))+
  theme(axis.text.x = element_text(vjust = 1,
                                   hjust = 1,
                                   size=12,
                                   angle = 60))+
  theme(axis.text.y = element_text(size=12))

plot_1 #

print(max(df_1$AP_PDT))
print(min(df_1$AP_PDT))

plot_2 <- ggplot(df_1, aes(x = reorder(Ligand, AP_PDT, mean), y = AP_PDT))+
  geom_boxplot()+
  geom_jitter(position = position_jitter(width = 0.1, height = 0.05),
              aes(color = Solvent), alpha = 0.8)+
  ggtitle(paste("High Throughput Experimentation"))+
  labs(subtitle = 'Ligand and Solvent Screening',
       x = paste("Ligand"), y = paste("Product (UPLC Area Percent)"))+
  theme(plot.background = element_rect(fill="white"))+
  theme(panel.background = element_rect(fill="white", colour="grey50"))+
  theme(plot.title = element_text(face = "bold", 
                                  size = 18,
                                  color="#BE2BBB"))+
  theme(axis.title = element_text(face = "bold", size = 16))+
  theme(axis.text.x = element_text(vjust = 1,
                                   hjust = 1,
                                   size=12,
                                   angle = 60))+
  theme(axis.text.y = element_text(size=12))
plot_2

df_1_ligand <- df_1 %>% group_by(Ligand) %>% summarize(Avg = mean(AP_PDT), Max = max(AP_PDT)) %>%
  arrange(desc(Avg))


df_2 <- df_1 %>% filter(!is.na(Mordred_Prediction))

df_b <- df %>% filter(Substrate == "1b")

df_2 <- df_2 %>% left_join(df_b, by = c("Ligand", "Solvent"))

df_2$error <- abs(df_2$AP_PDT-df_2$`Predicted Yield`)
df_2$sq_error <- df_2$error^2

df_2_ligand <- df_2 %>% group_by(Ligand) %>% summarize(avg_error = mean(error))

rmse <- sqrt(sum(df_2$sq_error)/nrow(df_2))

model <- lm(df_2$`Predicted Yield` ~ 0 + df_2$AP_PDT)
summary(model)
r2 <- round(summary(model)$r.squared, 2)

plot_3 <- ggplot(df_2, aes(x = AP_PDT, y = `Predicted Yield`))+
  geom_point(aes(color = Ligand, shape = Solvent), size = 4, alpha = 0.8)+
  geom_abline(intercept = 0, slope = 1)+
  geom_smooth(method=lm, se=FALSE, formula=y~x-1, linetype = "dashed")+
  xlim(0, 100)+
  ylim(0, 100)+
  ggtitle(paste("Machine Learning Model Performance"))+
  labs(x = paste("Actual Area Percent"), y = paste("Predicted Yield"))+
  annotate("text", x=85, y=65, label= paste0("r2 = ", r2)) +
  theme(plot.background = element_rect(fill="white"))+
  theme(panel.background = element_rect(fill="white", colour="grey50"))+
  theme(plot.title = element_text(face = "bold", 
                                  size = 18,
                                  color="#BE2BBB"))+
  theme(axis.title = element_text(face = "bold", size = 16))+
  theme(axis.text.x = element_text(vjust = 1,
                                   hjust = 1,
                                   size=12))+
  theme(axis.text.y = element_text(size=12))
plot_3

df_3 <- read.csv("Borylation_Updated (6).csv") %>% filter(Yield != "PENDING")

df_3$Yield <- as.numeric(df_3$Yield)
df_3$Cost <- as.numeric(df_3$Cost)
df_3$Round <- as.factor(df_3$Round)

df_3 <- df_3 %>% rename(`EDBO+ Selected` = Jay_Selected)

df_3$`EDBO+ Selected` <- ifelse(df_3$`EDBO+ Selected` == 1, "No", "Yes")

plot_4 <- ggplot(df_3, aes(x = Round, y = Yield))+
  geom_jitter(position = position_jitter(width = 0.1, height = 0.05),
              aes(size = Cost, color = `EDBO+ Selected`), alpha = 0.8)+
  #scale_color_viridis_d()+
  ggtitle(paste("EBDO Performance by Round"))+
  labs(subtitle = 'Ni Catalyzed Borylation with Ni(PPh3)2Cl2',
       x = paste("Experimentation Round"), y = paste("Solution Yield"))+
  theme(plot.background = element_rect(fill="white"))+
  theme(panel.background = element_rect(fill="white", colour="grey50"))+
  theme(plot.title = element_text(face = "bold", 
                                  size = 18,
                                  color="#BE2BBB"))+
  theme(axis.title = element_text(face = "bold", size = 16))+
  theme(axis.text.x = element_text(vjust = 1,
                                   hjust = 1,
                                   size=12))+
  theme(axis.text.y = element_text(size=12))
plot_4


df_3 <- df_3 %>% arrange(Round, Yield)


plot_4 <- ggplot(df_1, aes(x = Column, y = Row))+
  geom_point(aes(color = AP_PDT), alpha = 0.8, shape = 15, size = 20)+
  scale_color_viridis_c()+
  theme(plot.background = element_rect(fill="white"))+
  theme(panel.background = element_rect(fill="white", colour="grey50"))+
  theme(plot.title = element_blank())+
  theme(axis.title = element_blank())+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_blank())
plot_4


# Prediction Analysis -----------------------------------------------------

pred_ind <- c(1028, 4, 3, 2, 1, 0, 876, 540, 518, 488, 477, 274, 1280, 776, 574,
              479, 386, 111, 688, 667, 658, 257, 144, 31)

pred_1 <- read_csv('pred_Borylation_Updated.csv') %>% filter(Index %in% pred_ind) %>%
  select(Index, Yield_predicted_mean, Yield_predicted_variance,
         Cost_predicted_mean, Cost_predicted_variance) %>%
  rename(post_Rd_1_pred_yield = Yield_predicted_mean) %>%
  rename(post_Rd_1_var_yield = Yield_predicted_variance) %>%
  rename(post_Rd_1_pred_cost = Cost_predicted_mean) %>%
  rename(post_Rd_1_var_cost = Cost_predicted_variance)
pred_2 <- read_csv('pred_Borylation_Updated (1).csv') %>% filter(Index %in% pred_ind) %>%
  select(Index, Yield_predicted_mean, Yield_predicted_variance,
         Cost_predicted_mean, Cost_predicted_variance) %>%
  rename(post_Rd_2_pred_yield = Yield_predicted_mean) %>%
  rename(post_Rd_2_var_yield = Yield_predicted_variance)%>%
  rename(post_Rd_2_pred_cost = Cost_predicted_mean) %>%
  rename(post_Rd_2_var_cost = Cost_predicted_variance)
pred_3 <- read_csv('pred_Borylation_Updated (2).csv') %>% filter(Index %in% pred_ind) %>%
  rename(post_Rd_3_pred_yield = Yield_predicted_mean) %>%
  rename(post_Rd_3_var_yield = Yield_predicted_variance)%>%
  rename(post_Rd_3_pred_cost = Cost_predicted_mean) %>%
  rename(post_Rd_3_var_cost = Cost_predicted_variance)

pred <- pred_3 %>% left_join(pred_1, by = 'Index')
pred <- pred %>% left_join(pred_2, by = 'Index')

pred_yield <- pred %>% select(Index, Round, Yield, post_Rd_1_pred_yield,
                                  post_Rd_2_pred_yield, post_Rd_3_pred_yield)


pred_yield$Yield <- as.numeric(pred_yield$Yield)

pred_yield$Entry <- 1:nrow(pred_yield)

plot_5 <- ggplot(pred_yield)+
  geom_point(aes(x = Entry, y = Yield), size = 5, color = '#BE2BBB', shape = 19)+
  geom_point(aes(x = Entry, y = post_Rd_1_pred_yield), size = 3, color = '#595454')+
  geom_point(aes(x = Entry, y = post_Rd_2_pred_yield), size = 3, color = '#595454', shape = 17)+
  geom_point(aes(x = Entry, y = post_Rd_3_pred_yield), size = 3, color = '#595454', shape = 18)+
  scale_shape_manual(name = "Legend",
                     labels = c("Experimental Yield","Post-Rd 1 Prediction",
                                "Post-Rd 2 Prediction", "Post-Rd 3 Prediction"),
                     values = c(19, 15, 17, 18))

plot_5

