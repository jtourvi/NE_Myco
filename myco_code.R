#Begin - load packages, setwd
#install.packages("naniar")
#install.packages("reshape")
#install.packages("survival")
#install.packages("survminer")
#install.packages("lubridate")
#install.packages("gridExtra")
library(ggplot2)
library(naniar)
library(Ramf)
library(plyr)
library(tidyr)
library(dplyr)
library(reshape)
library(survival)
library(survminer)
library(lubridate)
library(lme4)
library(cowplot)
library(gridExtra)
#install.packages("grid")
library(grid)

#read in data for network and transplant experiments - final measuremeants and all together
trans_final = read.csv("trans_final.csv", header = T, stringsAsFactors = F, na.strings = "NA")
net_final = read.csv("net_final.csv", header = T, stringsAsFactors = F, na.strings = "NA")
trans_full = read.csv("trans_full.csv", header = T, stringsAsFactors = F, na.strings = "NA")
net_full = read.csv("net_full.csv", header = T, stringsAsFactors = F, na.strings = "NA")

#remove excess columns
trans = trans_final[c(1:288), ]
net = net_final[c(1:290), ]
trans_full = trans_full[c(1:3156), ]

#deal with pesky NA and wierd values from excel
trans[trans == "#VALUE!" ] <- NA
trans[trans == "#DIV/0!" ] <- NA
net[net == "#VALUE!" ] <- NA
trans_full[trans_full == "#DIV/0!" ] <- 0
trans_full[trans_full == "#VALUE!" ] <- NA
net_full[net_full == "#VALUE!" ] <- NA

#get columns into right form
trans$root_shoot_ratio1 = as.numeric(trans$root_shoot_ratio1)
trans$root_shoot_ratio2 = as.numeric(trans$root_shoot_ratio2)
trans$root_shoot_ratio3 = as.numeric(trans$root_shoot_ratio3)
trans$root_shoot_ratio4 = as.numeric(trans$root_shoot_ratio4)
trans$root_shoot_ratio5 = as.numeric(trans$root_shoot_ratio5)
trans$root_shoot_ratio6 = as.numeric(trans$root_shoot_ratio6)
trans$root_shoot_ratio7 = as.numeric(trans$root_shoot_ratio7)

trans$total_bio1 = as.numeric(trans$total_bio1)
trans$total_bio2 = as.numeric(trans$total_bio2)
trans$total_bio3 = as.numeric(trans$total_bio3)
trans$total_bio4 = as.numeric(trans$total_bio4)
trans$total_bio5 = as.numeric(trans$total_bio5)
trans$total_bio6 = as.numeric(trans$total_bio6)
trans$total_bio7 = as.numeric(trans$total_bio7)

net$total_bio1 = as.numeric(net$total_bio1)
net$total_bio2 = as.numeric(net$total_bio2)
net$total_bio3 = as.numeric(net$total_bio3)

#calculate averages for these things across rows
#for trans
trans$root_shoot_ratio_avg = rowMeans(trans[ , c(23:29)], na.rm=TRUE)
trans$total_bio_avg = rowMeans(trans[ , c(30:36)], na.rm=TRUE)
trans$total_abv_avg = rowMeans(trans[ , c(9:15)], na.rm=TRUE)
trans$total_below_avg = rowMeans(trans[ , c(16:22)], na.rm=TRUE)

#for net
net$root_shoot_ratio_avg = rowMeans(net[ , c(14:16)], na.rm=TRUE)
net$total_bio_avg = rowMeans(net[ , c(17:19)], na.rm=TRUE)
net$total_abv_avg = rowMeans(net[ , c(8:10)], na.rm=TRUE)
net$total_below_avg = rowMeans(net[ , c(11:13)], na.rm=TRUE)

#merge net and trans seedling data with herb, browse, and pathogen observations
trans_bio = ddply(trans_full, ~mountain + site + spp + treatment + source + rep, summarise, root_path = sum(root_pathogen),
                  herb = sum(herbivory), browse = sum(browse))
trans_bio = trans_bio[c(2:290), ]
trans_bio[is.na(trans_bio)] = 1

trans_bio = trans_bio %>% mutate_if(is.numeric, ~1 * (. != 0))

trans_all = merge(trans, trans_bio, by = c("mountain", "site", "spp", "treatment", "source", "rep"))

net_bio = ddply(net_full, ~mountain + site + spp + treatment + rep, summarise, root_path = sum(root_pathogen),
                herb = sum(herbivory), browse = sum(browse))
net_bio = net_bio[c(2:294), ]
net_bio[is.na(net_bio)] = 1

net_bio = net_bio %>% mutate_if(is.numeric, ~1 * (. != 0))

net_all = merge(net, net_bio, by = c("mountain", "site", "spp", "treatment", "rep"))
#write.csv(net_all, "net_all.csv")
net_all = read.csv("net_all.csv", header = TRUE)

###
###aux data sorting here - merge all env data into one df
myco_env = read.csv("myco_env.csv", header = TRUE)
chars = read.csv("site_chars.csv", header = TRUE)
light = read.csv("light.csv", header = TRUE)
covers = read.csv("cover_totals.csv", header = TRUE)
tissues = read.csv("leaf_tissue_analysis.csv", header = TRUE)
soils = read.csv("VT_soil_analysis.csv", header = TRUE)
over = read.csv("overstory.csv", header = TRUE)

light1 = ddply(light, ~mountain + site + point, summarise, avg_light = mean(perc_open))

cover = cast(covers, mountain + site + point ~ cover_type, value = "cover_total", fun.aggregate = mean)

myco_env = myco_env[ ,1:20]

soil = ddply(soils, ~mountain + site, summarise, soil_pH = mean(pH), som = mean(OM_Pct), soil_P = mean(Avail_P), soil_K = mean(K), soil_Ca = mean(Ca),
             soil_Mg = mean(Mg), soil_Zn = mean(Zn), soil_B = mean(B), soil_Mn = mean(Mn), soil_Cu = mean(Cu), soil_Fe = mean(Fe), soil_Al = mean(Al), soil_Na = mean(Na),
             soil_S = mean(S), Ex_acid = mean(Exch_Acid), CEC = mean(ECEC))

lm_ph = lm(soil_pH ~ site, data = soil)
aovp = aov(lm_ph)
TukeyHSD(aovp)

lm_som = lm(som ~ site, data = soil)
aovs = aov(lm_som)
TukeyHSD(aovs)

lm_cec = lm(CEC ~ site, data = soil)
aovc = aov(lm_cec)
TukeyHSD(aovc)

lm_al = lm(soil_Al ~ site, data = soil)
aova = aov(lm_al)
TukeyHSD(aova)

lm_p = lm(soil_P ~ site, data = soil)
aovpp = aov(lm_p)
TukeyHSD(aovpp)

lm_k = lm(soil_K ~ site, data = soil)
aovk = aov(lm_k)
TukeyHSD(aovk)

lm_ca = lm(soil_Ca ~ site, data = soil)
aovcc = aov(lm_ca)
TukeyHSD(aovcc)

lm_mg = lm(soil_Mg ~ site, data = soil)
aovm = aov(lm_mg)
TukeyHSD(aovm)

myco_env1 = merge(myco_env, chars, by = c("mountain", "site", "point"), all.x = TRUE)
myco_env2 = merge(myco_env1, light1, by = c("mountain", "site", "point"), all.x = TRUE)
myco_env3 = merge(myco_env2, cover, by = c("mountain", "site", "point"), all.x = TRUE)
myco_env4 = merge(myco_env3, tissues, by = c("mountain", "site", "spp", "experiment", "treat"), all.x = TRUE)
myco_env5 = merge(myco_env4, soil, by = c("mountain", "site"), all.x = TRUE)
myco_env6 = merge(myco_env5, over, by = c("mountain", "site"), all.x = TRUE)

myco_env7 = merge(myco_env6, amf_sub, by = c("mountain", "site", "spp", "experiment", "treatment"), all = TRUE)
myco_env7$col_perc = myco_env7$arb_prop + myco_env7$ves_prop + myco_env7$coil_prop + myco_env7$hyphae_prop
myco_env8 = merge(myco_env7, emf, by = c("mountain", "site", "spp", "experiment", "treatment"), all = TRUE)
write.csv(myco_env8, "myco_all1.csv")

#data wranging - I've calculated root_shoot averages, total bio averages, survival, all new data, diameters,
#pathogens, herb, and browse. Final growth comparisons, survival analysis, foliar nutrients,
#mycorrhizal colonization, repeated measures ANOVA which includes calculating RHGR.
#Wrapped up with GLMM for all to serve as statistical test

#run glm or glmm for everything of interest - greenhouse, trans, and net
#both acsa and fagr
#make plots for all

#mycorrhizal colonization - amf and emf
amf = read.csv("AMF_colonization.csv", header = TRUE)
emf = read.csv("EMF_colonization.csv", header = TRUE)
amf1 = ddply(amf, ~experiment + site + treat, summarise, 
             arb = mean(arb_prop), SEa = sd(arb_prop)/sqrt((length(arb_prop))),   
             ves = mean(ves_prop), SEv = sd(ves_prop)/sqrt((length(ves_prop))),   
             coil = mean(coil_prop), SEc = sd(coil_prop)/sqrt((length(coil_prop))),   
             hyp = mean(hyphae_prop), SEh = sd(hyphae_prop)/sqrt((length(hyphae_prop))),   
             endo = mean(endo_prop), SEe = sd(endo_prop)/sqrt((length(endo_prop))))

hist(amf$col_perc)
amf_g = subset(amf, experiment == "green")

ba = read.csv("myco_tree_ba.csv", header = TRUE)
ba1 = ddply(ba, ~site, summarise, conifer = mean(emf_ba), SEc = sd(emf_ba)/sqrt((length(emf_ba))),
            dec = mean(amf_ba), SEd = sd(amf_ba)/sqrt((length(amf_ba))))

amf$col_perc = (amf$arb_prop + amf$ves_prop + amf$coil_prop + amf$hyphae_prop)

acsa_col_glm = glm(col_perc ~ site * treat, data = amf_g, family = "gaussian")
summary(acsa_col_glm)
anova(acsa_col_glm)
Anova(acsa_col_glm)

amf_g$BV <- interaction(amf_g$site, amf_g$treat)

m1  <- glmer(col_perc ~ -1 + BV + (1|mountain), family = gaussian, data = amf_g)

summary(glht(m1, linfct = mcp(BV = "Tukey")))

net_my = subset(amf1, experiment == "net")
trans_my = subset(amf, experiment == "trans")
nat_my = subset(amf1, experiment == "nat")
green_my = subset(amf, experiment == "green")

green_my$col_perc = (green_my$arb_prop + green_my$ves_prop + green_my$coil_prop + green_my$hyphae_prop)*100
green_my$col_perc = (green_my$arb + green_my$ves + green_my$coil + green_my$hyp)*100
green_my = ddply(green_my, ~order + order1, summarise, col_perc1 = mean(col_perc), col_perc_se = sd(col_perc)/sqrt((length(col_perc))))
new_row = c(order = "Hardwood soils", order1="FA", col_perc1 = 0.33, col_perc_se = 0.33)
green_my = rbind(green_my,new_row)
net_my$col_perc = (net_my$arb + net_my$ves + net_my$coil + net_my$hyp)*100#calculating colonization %
net_my$col_perc_se = (net_my$SEa + net_my$SEv + net_my$SEc + net_my$SEh)*100

net_my$col_perc = (net_my$arb_prop + net_my$ves_prop + net_my$coil_prop + net_my$hyphae_prop)
nat_my$col_perc_se = (nat_my$SEa + nat_my$SEv + nat_my$SEc + nat_my$SEh)*100

amf1$col_perc = (amf1$arb + amf1$ves + amf1$coil + amf1$hyp)*100
amf1$col_perc_se = (amf1$SEa + amf1$SEv + amf1$SEc + amf1$SEh)*100
amf_sub = subset(amf, experiment == "net" | experiment == "trans")
write.csv(amf1, "amf1.csv")

amf_trans1 = ddply(amf_trans, ~site + treat, summarise, 
                   col = mean(col_perc), SEh = sd(col_perc)/sqrt((length(col_perc))))

amf_trans = subset(amf, experiment == "trans")

emf_trans = subset(emf, experiment == "trans")
emf_trans$em_colonization = emf_trans$em_colonization/100
emf_net = subset(emf, experiment == "net")
emf_net$em_colonization = emf_net$em_colonization/100
emf_nat = subset(emf, experiment == "nat")

emf_net1 = ddply(emf_net, ~site + source, summarise, 
                 col = mean(em_colonization), SEh = sd(em_colonization)/sqrt((length(em_colonization))),   
                 vit = mean(vitality_index), SEe = sd(vitality_index)/sqrt((length(vitality_index))))

emf_nat1 = ddply(emf_nat, ~site, summarise, 
                 col = mean(em_colonization), SEh = sd(em_colonization)/sqrt((length(em_colonization))),   
                 vit = mean(vitality_index), SEe = sd(vitality_index)/sqrt((length(vitality_index))))
emf_nat1 = merge(emf_nat1, ba1, by = "site")

#trans amf colonization

amf1_trans = glmer(col_perc ~ site * treat + (1|mountain), data = trans_my, family="binomial", na.action=na.omit)
summary(amf1_trans)

K1 <- glht(amf1_trans, mcp(site = "Tukey"))$linfct
K2 <- glht(amf1_trans, mcp(treat = "Tukey"))$linfct

summary(glht(amf1_trans, linfct = rbind(K1, K2)))

tmp <- expand.grid(site = unique(trans_my$site),
                   treat = unique(trans_my$treat))
X <- model.matrix(~ site * treat, data = tmp)
glht(amf1_trans, linfct = X)

Tukey <- contrMat(table(warpbreaks$tension), "Tukey")
K1 <- cbind(Tukey, matrix(0, nrow = nrow(Tukey), ncol = ncol(Tukey)))
rownames(K1) <- paste(levels(warpbreaks$wool)[1], rownames(K1), sep = ":")
K2 <- cbind(matrix(0, nrow = nrow(Tukey), ncol = ncol(Tukey)), Tukey)
rownames(K2) <- paste(levels(warpbreaks$wool)[2], rownames(K2), sep = ":")
K <- rbind(K1, K2)
colnames(K) <- c(colnames(Tukey), colnames(Tukey))

summary(glht(mod, linfct = K %*% X))

symbol12 = c("a", "a", "b", "a", "a", "b", "a", "ab", "c")
amf_trans1$order = factor(amf_trans1$site, levels=c("L", "M", "H"),
                          labels = c("Hardwood soils", "Ecotone soils", "Spruce-fir soils"))
amf_trans1$order1 = factor(amf_trans1$treat, levels=c("L", "M", "H"))

amf_trans_plot = ggplot(amf_trans1, aes(x = order1, y = col)) +
  geom_errorbar(stat = "identity", aes(ymin=col-SEh, ymax=col+SEh), position = position_dodge(width = 0.90), width = 0.1) +
  geom_point(stat = "identity", size = 4, pch = 18) +
  ylab("Root Length Colonized (%)") +
  ggtitle("Sugar Maple - AMF Colonization") +
  ylim(0, 70) +
  xlab("") +
  scale_x_discrete(labels=c("L" = "Hardwoods", "M" = "Ecotone",
                            "H" = "Spruce-fir")) +
  facet_wrap(~order, scales = "free_x", ncol = 3) +
  #facet_grid(cols = vars(order)) +
  theme_minimal() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        plot.title = element_text(size=12),
        legend.title = element_text(size=12), 
        legend.text = element_text(size=12)) +
  stat_summary(geom = 'text', label = symbol12, fun.y = max, vjust = -3, position = position_dodge(width = 0.90))
amf_trans_plot 

symbol = c("a", "a", "b", "ab", "b", "c", "ab", "c", "d")
amf_trans_plot  = ggplot(amf_trans1, aes(x = order1, y = col, fill = order1)) +
  facet_grid(cols = vars(order)) +
  geom_errorbar(aes(ymin=col-SEh, ymax=col+SEh), width = 0.2, position=position_dodge(width=0.6), colour="black", size = 1) +
  geom_point(stat = "identity", position = position_dodge(width = 0.6), pch = 23, size = 6) +
  ylab("") +
  xlab("") +
  ylim(0,70) +
  theme_bw() +
  scale_fill_manual(name = "Transplant Site", values = c("#2c7c94", "#fbe45b", "#a65852"),
                    labels = c("L" = "Hardwoods", "M" = "Ecotone",
                               "H" = "Spruce-fir")) +
  scale_x_discrete(labels=c("L" = "Hardwoods", "M" = "Ecotone",
                            "H" = "Spruce-fir")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        strip.text.x = element_text(size = 20),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size = 20),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size=22),
        legend.title = element_text(size=22), 
        legend.text = element_text(size=22),
        legend.position = "none") +
  stat_summary(geom = 'text', label = symbol, fun = max, vjust = -2, position = position_dodge(width = 0.6), size = 8)
amf_trans_plot

hist(amf_trans$col_perc)

glm1 = glmer(col_perc ~ treat * site + (1|mountain), family = gaussian, data = amf_trans)
summary(glm1)
anova(glm1)
Anova(glm1)

amf_trans$BV <- interaction(amf_trans$site, amf_trans$treat)

m1  <- glmer(col_perc ~ -1 + BV + (1|mountain), family = gaussian, data = amf_trans)

summary(glht(m1, linfct = mcp(BV = "Tukey")))
app = glht(m1, linfct = mcp(BV = "Tukey"))
cld(app)

K1 <- glht(glm1, mcp(site = "Tukey"))$linfct
K2 <- glht(glm1, mcp(treat = "Tukey"))$linfct

summary(glht(glm1, linfct = rbind(K1, K2)))

#trans emf colonization

emf1_trans = glmer(em_colonization ~ site * source + (1|mountain), data = emf_trans, family="binomial", na.action=na.omit)
summary(emf1_trans)
anova(emf1_trans)
Anova(emf1_trans)

K1 <- glht(emf1_trans, mcp(site = "Tukey"))$linfct
K2 <- glht(emf1_trans, mcp(source = "Tukey"))$linfct

summary(glht(emf1_trans, linfct = rbind(K1, K2)))


emf_trans = ddply(emf_trans, ~site + source, summarise, 
                  col = mean(em_colonization), SEh = sd(em_colonization)/sqrt((length(em_colonization))),   
                  vit = mean(vitality_index), SEe = sd(vitality_index)/sqrt((length(vitality_index))))
emf_trans$col = emf_trans$col*100
emf_trans$SEh = emf_trans$SEh*100

emf_trans$order = factor(emf_trans$site, levels=c("L", "M", "H"),
                         labels = c("Hardwood soils", "Ecotone soils", "Spruce-fir soils"))
emf_trans$order1 = factor(emf_trans$source, levels=c("L", "M", "H"))

emf_trans_plot = ggplot(emf_trans, aes(x = order1, y = col)) +
  geom_errorbar(stat = "identity", aes(ymin=col-SEh, ymax=col+SEh), position = position_dodge(width = 0.90), width = 0.1) +
  geom_point(stat = "identity", size = 4, pch = 18) +
  ylab("Root Length Colonized (%)") +
  ggtitle("American beech - EMF Colonization") +
  ylim(0, 70) +
  xlab("") +
  scale_x_discrete(labels=c("L" = "Hardwoods", "M" = "Ecotone",
                            "H" = "Spruce-fir")) +
  facet_wrap(~order, scales = "free_x", ncol = 3) +
  #facet_grid(cols = vars(order)) +
  theme_minimal() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        plot.title = element_text(size=12),
        legend.title = element_text(size=12), 
        legend.text = element_text(size=12))
#stat_summary(geom = 'text', label = plot.dat$sig, fun.y = max, vjust = -1, position = position_dodge(width = 0.75))
emf_trans_plot 

symbol = c("a", "a", "a", "a", "a", "a", "a", "a", "a")
emf_trans_plot  = ggplot(emf_trans, aes(x = order1, y = col, fill = order1)) +
  facet_grid(cols = vars(order)) +
  geom_errorbar(aes(ymin=col-SEh, ymax=col+SEh), width = 0.2, position=position_dodge(width=0.6), colour="black", size = 1) +
  geom_point(stat = "identity", position = position_dodge(width = 0.6), pch = 23, size = 6) +
  ylab("") +
  xlab("") +
  ylim(0,70) +
  theme_bw() +
  scale_fill_manual(name = "Transplant Site", values = c("#2c7c94", "#fbe45b", "#a65852"),
                    labels = c("L" = "Hardwoods", "M" = "Ecotone",
                               "H" = "Spruce-fir")) +
  scale_x_discrete(labels=c("L" = "Hardwoods", "M" = "Ecotone",
                            "H" = "Spruce-fir")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        strip.text.x = element_text(size = 20),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size = 20),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size=22),
        legend.title = element_text(size=22), 
        legend.text = element_text(size=22),
        legend.position = "none") +
  stat_summary(geom = 'text', label = symbol, fun = max, vjust = -2, position = position_dodge(width = 0.6), size = 8)
emf_trans_plot

tiff(file = "trans_col.tif", width = 12, height = 8, units = 'in', res = 600, pointsize = 11)

trans_col = plot_grid(amf_trans_plot, emf_trans_plot, ncol = 2, nrow = 1)
trans_col

dev.off()

#myco colonization for net experiment

nat_my = merge(nat_my, ba1, by = "site")
nat_my$dec = as.numeric(nat_my$dec)
nat_my$col_perc = as.numeric(nat_my$col_perc)
nat_my$order = factor(nat_my$site, levels=c("500", "600", "700", "800", "900", "1000"))

net_my$order1 = factor(net_my$treat, levels=c("N", "D", "C"))
net_my$order2 = factor(net_my$site, levels=c("L", "M", "H"))

amf1_net = lm(col_perc ~ site * treatment, data = net_my)
summary(amf1_net)
anova(amf1_net)
Anova(amf1_net)

net_my$BV <- interaction(net_my$site, net_my$treatment)

m1  <- glmer(col_perc ~ -1 + BV + (1|mountain), family = gaussian, data = net_my)

summary(glht(m1, linfct = mcp(BV = "Tukey")))

tiff(file = "amf_col.tif", width = 8, height = 5, units = 'in', res = 600, pointsize = 11)

#amf network plot
amf_net_plot = ggplot(net_my, aes(x = order2, y = col_perc, group = order1)) +
  geom_errorbar(stat = "identity", aes(ymin=col_perc-col_perc_se, ymax=col_perc+col_perc_se), position = position_dodge(width = 0.90), width = 0.1) +
  geom_point(stat = "identity", aes(color = order1), size = 4, pch = 18, position = position_dodge(width = 0.90)) +
  scale_color_brewer(palette = "Set1", name = "Treatment",
                     labels = c("Network (N)", "Disrupted (D)", "Control (C)")) +
  ylab("Root Length Colonized (%)") +
  xlab("Northern Hardwoods                    Ecotone                          Spruce-fir      ") +
  ggtitle("Sugar Maple - AMF Colonization") +
  ylim(0, 60) +
  facet_wrap(~order2, scales = "free_x", ncol = 3) +
  #facet_grid(var ~ order, scales = "free_x") +
  theme(axis.text.x = element_blank()) +
  theme_minimal() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(strip.text.x = element_blank()) +
  theme(axis.text.x = element_blank(),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        plot.title = element_text(size=12),
        legend.title = element_text(size=12), 
        legend.text = element_text(size=12))
#stat_summary(geom = 'text', label = plot.dat$sig, fun.y = max, vjust = -1, position = position_dodge(width = 0.75))
amf_net_plot 

dev.off()

emf1_net = lmer(em_colonization ~ site * source + (1|mountain), data = emf_net)
summary(emf1_net)
anova(emf1_net)
Anova(emf1_net)

emf_net$BV <- interaction(emf_net$site, emf_net$source)

m1  <- glmer(em_colonization ~ -1 + BV + (1|mountain), family = gaussian, data = emf_net)

summary(glht(m1, linfct = mcp(BV = "Tukey")))

#emf network plot

emf_net1$order1 = factor(emf_net1$source, levels=c("N", "D", "C"))
emf_net1$order2 = factor(emf_net1$site, levels=c("L", "M", "H"))

emf_net_plot = ggplot(emf_net1, aes(x = order2, y = col, group = order1)) +
  geom_errorbar(stat = "identity", aes(ymin=col-SEh, ymax=col+SEh), position = position_dodge(width = 0.90), width = 0.1) +
  geom_point(stat = "identity", aes(color = order1), size = 4, pch = 18, position = position_dodge(width = 0.90)) +
  scale_color_brewer(palette = "Set1", name = "Treatment",
                     labels = c("Network (N)", "Disrupted (D)", "Control (C)")) +
  ylab("Root Length Colonized (%)") +
  xlab("Northern Hardwoods                    Ecotone                          Spruce-fir      ") +
  ggtitle("Beech - EMF Colonization") +
  ylim(0, 60) +
  facet_wrap(~order2, scales = "free_x", ncol = 3) +
  #facet_grid(var ~ order, scales = "free_x") +
  theme(axis.text.x = element_blank()) +
  theme_minimal() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(strip.text.x = element_blank()) +
  theme(axis.text.x = element_blank(),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        plot.title = element_text(size=12),
        legend.title = element_text(size=12), 
        legend.text = element_text(size=12))
#stat_summary(geom = 'text', label = plot.dat$sig, fun.y = max, vjust = -1, position = position_dodge(width = 0.75))
tiff(file = "emf_net_col.tif", width = 8, height = 5, units = 'in', res = 600, pointsize = 11)
emf_net_plot 
dev.off()

#greenhouse amf plot

green_my$col_perc = (green_my$arb_prop + green_my$ves_prop + green_my$coil_prop + green_my$hyphae_prop)*100
green_my$col_perc_se = (green_my$SEa + green_my$SEv + green_my$SEc + green_my$SEh)*100

green_my$order = factor(green_my$site, levels=c("L", "M", "H"),
                        labels = c("Hardwood soils", "Ecotone soils", "Spruce-fir soils"))
green_my$order1 = factor(green_my$treat, levels=c("B", "FA", "FC", "PC"))

green_amf = ggplot(green_my, aes(x = order, y = col_perc, group = order1)) +
  geom_errorbar(stat = "identity", aes(ymin=col_perc-col_perc_se, ymax=col_perc+col_perc_se), position = position_dodge(width = 0.90), width = 0.1) +
  geom_point(stat = "identity", aes(color = order1), size = 4, pch = 18, position = position_dodge(width = 0.90)) +
  scale_color_brewer(palette = "Set1", name = "Treatment",
                     labels = c("Innoculated (B)", "Field Autoclave (FA)", "Field Control (FC)", "Pot Control (PC)")) +
  ylab("Root Length Colonized (%)") +
  xlab("") +
  ggtitle("Sugar Maple - AMF Colonization") +
  ylim(0, 30) +
  #facet_wrap(~order, scales = "free_x", ncol = 3) +
  facet_grid(~order, scales = "free_x") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_blank(),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        plot.title = element_text(size=12),
        legend.title = element_text(size=12), 
        legend.text = element_text(size=12),
        axis.ticks.x=element_blank())
#stat_summary(geom = 'text', label = plot.dat$sig, fun.y = max, vjust = -1, position = position_dodge(width = 0.75))
green_amf

tiff(file = "green_col.tif", width = 20, height = 12, units = 'in', res = 600, pointsize = 11)

dev.off()

#greenhouse emf

emg = read.csv("EMF_green.csv", header = TRUE)

emg_col = glm(em_colonization ~ site + mountain, data = emg, family = "gaussian")
summary(emg_col)

emg1 = ddply(emg, ~ mountain + site, summarise, col_perc = mean(em_colonization), SE = sd(em_colonization)/sqrt((length(em_colonization))))

emg1$order = factor(emg1$site, levels=c("L", "M", "H"),
                    labels = c("Hardwood soils", "Ecotone soils", "Spruce-fir soils"))
emg1$order1 = factor(emg1$mountain, levels=c("B", "FA", "FC", "PC"))

green_emf = ggplot(emg1, aes(x = order, y = col_perc, group = order1)) +
  geom_errorbar(stat = "identity", aes(ymin=col_perc-SE, ymax=col_perc+SE), position = position_dodge(width = 0.90), width = 0.1) +
  geom_point(stat = "identity", aes(color = order1), size = 4, pch = 18, position = position_dodge(width = 0.90)) +
  scale_color_brewer(palette = "Set1", name = "Treatment",
                     labels = c("Innoculated (B)", "Field Autoclave (FA)", "Field Control (FC)", "Pot Control (PC)")) +
  ylab("Root Length Colonized (%)") +
  xlab("") +
  ggtitle("American Beech - EMF Colonization") +
  ylim(0, 60) +
  #facet_wrap(~order, scales = "free_x", ncol = 3) +
  facet_grid(~order, scales = "free_x") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_blank(),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        plot.title = element_text(size=12),
        legend.title = element_text(size=12), 
        legend.text = element_text(size=12),
        axis.ticks.x=element_blank())
#stat_summary(geom = 'text', label = plot.dat$sig, fun.y = max, vjust = -1, position = position_dodge(width = 0.75))
green_emf

tiff(file = "green_col.tif", width = 20, height = 12, units = 'in', res = 600, pointsize = 11)

green_col = plot_grid(green_amf, green_emf, nrow = 1, ncol = 2)
green_col

dev.off()

###germination trans

trans_all$germ = (trans_all$num_germ/20)*100
trans_ac = subset(trans_all, spp == "ACSA")
trans_fa = subset(trans_all, spp == "FAGR")
trans_ga1 = ddply(trans_ac, ~site + source, summarise, g = mean(germ), SE = sd(germ)/sqrt((length(germ))))
trans_g$g = as.numeric(trans_g$g)
trans_g$site = as.factor(trans_g$site)
trans_g$source = as.factor(trans_g$source)

trans_gf1 = ddply(trans_fa, ~site + source, summarise, g = mean(germ), SE = sd(germ)/sqrt((length(germ))))

trans_sa1 = ddply(trans_ac, ~site, summarise, g = mean(germ), SE = sd(germ)/sqrt((length(germ))))
trans_sf1 = ddply(trans_fa, ~site, summarise, g = mean(germ), SE = sd(germ)/sqrt((length(germ))))

trans_ga1$order = factor(trans_ga1$site, levels=c("L", "M", "H"))
trans_ga1$order1 = factor(trans_ga1$source, levels=c("L", "M", "H"),
                          labels = c("Hardwood soils", "Ecotone soils", "Spruce-fir soils"))
trans_sa1$order = factor(trans_sa1$site, levels=c("L", "M", "H"))

trans_gf1$order = factor(trans_gf1$site, levels=c("L", "M", "H"))
trans_gf1$order1 = factor(trans_gf1$source, levels=c("L", "M", "H"),
                          labels = c("Hardwood soils", "Ecotone soils", "Spruce-fir soils"))
trans_gf1$spp = "American beech"
trans_ga1$spp = "Sugar maple"
trans_g1 = rbind(trans_ga1, trans_gf1)
trans_g1$spp = factor(trans_g1$spp, levels=c("Sugar maple", "American beech"))

trans_sf1$order = factor(trans_sf$site, levels=c("L", "M", "H"))
trans_ac$germ = trans_ac$germ/100
trans_fa$germ = trans_fa$germ/100
hist(trans_fa$germ)

glmm1 = glmer(germ ~ site * source + (1|mountain), family = gaussian, data = trans_ac)
summary(glmm1)
anova(glmm1)
Anova(glmm1)

trans_ac$BV <- interaction(trans_ac$site, trans_ac$source)

m1  <- glmer(germ ~ -1 + BV + (1|mountain), family = gaussian, data = trans_ac)

summary(glht(m1, linfct = mcp(BV = "Tukey")))
app = glht(m1, linfct = mcp(BV = "Tukey"))
cld(app)

glmm2 = glmer(germ ~ site * source + (1|mountain), family = gaussian, data = trans_fa)
summary(glmm2)
anova(glmm2)
Anova(glmm2)

trans_fa$BV <- interaction(trans_fa$site, trans_fa$source)

m1  <- glmer(germ ~ -1 + BV + (1|mountain), family = gaussian, data = trans_fa)

summary(glht(m1, linfct = mcp(BV = "Tukey")))
app = glht(m1, linfct = mcp(BV = "Tukey"))
cld(app)

glmm3 = lmer(germ ~ site + source + (1|mountain), data = trans_ac)
summary(glmm3)
Anova(glmm3)
glmm4 = mixed(germ ~ site + source + (1|mountain), data = trans_fa, method = "LRT")
summary(glmm4)
glmm5 = lmer(germ ~ site + source + (1|mountain), data = trans_fa)
summary(glmm5)

#acsa
symbol7 = c("a", "b", "c", "ab", "c", "c", "b", "c", "d")
germ_plota = ggplot(trans_ga1, aes(x = order, y = g)) +
  geom_errorbar(aes(ymin=g-SE, ymax=g+SE), width = 0.075) +
  facet_grid(cols = vars(order1)) +
  geom_bar(stat = "identity", fill = "#fbe45b", width = 0.4) +
  ylab("") +
  xlab("") +
  ylim(0,35) +
  theme_bw() +
  scale_x_discrete(labels=c("L" = "Hardwoods", "M" = "Ecotone",
                            "H" = "Spruce-fir")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        strip.text.x = element_text(size = 20),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        plot.title = element_text(size=22),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  stat_summary(geom = 'text', label = symbol7, fun.y = max, vjust = -1.5, position = position_dodge(width = 0.9), size = 8)
germ_plota

symbol8 = c("a", "b", "c")
s_plota = ggplot(trans_sa1, aes(x = order, y = g)) +
  geom_errorbar(aes(ymin=g-SE, ymax=g+SE), width = 0.075) +
  geom_bar(stat = "identity", fill = "#fbe45b", width = 0.4) +
  ylab("") +
  xlab("") +
  ylim(0,30) +
  theme_bw() +
  ggtitle("Sugar maple") +
  scale_x_discrete(labels=c("L" = "Hardwoods", "M" = "Ecotone",
                            "H" = "Spruce-fir")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        strip.text.x = element_text(size = 20),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(size=22)) +
  stat_summary(geom = 'text', label = symbol8, fun.y = max, vjust = -1, position = position_dodge(width = 0.9), size = 8)
s_plota

#fagr
symbol9 = c("a", "b", "b", "b", "b", "b", "c", "c", "c")
sur_plotf = ggplot(trans_gf1, aes(x = order, y = g)) +
  geom_errorbar(aes(ymin=g-SE, ymax=g+SE), width = 0.075) +
  facet_grid(cols = vars(order1)) +
  geom_bar(stat = "identity", fill = "#fbe45b", width = 0.4) +
  ylab("") +
  xlab("") +
  ylim(0,30) +
  theme_bw() +
  scale_x_discrete(labels=c("L" = "Hardwoods", "M" = "Ecotone",
                            "H" = "Spruce-fir")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        strip.text.x = element_text(size = 20),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(size=22)) +
  stat_summary(geom = 'text', label = symbol9, fun.y = max, vjust = -2.5, position = position_dodge(width = 0.9), size = 8)
sur_plotf

symbol10 = c("a", "a", "b")
su_plotf = ggplot(trans_sf1, aes(x = order, y = g)) +
  geom_errorbar(aes(ymin=g-SE, ymax=g+SE), width = 0.075) +
  geom_bar(stat = "identity", fill = "#fbe45b", width = 0.4) +
  ylab("") +
  xlab("") +
  ylim(0,30) +
  theme_bw() +
  ggtitle("American beech") +
  scale_x_discrete(labels=c("L" = "Hardwoods", "M" = "Ecotone",
                            "H" = "Spruce-fir")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        strip.text.x = element_text(size = 20),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(size=22)) +
  stat_summary(geom = 'text', label = symbol10, fun.y = max, vjust = -2, position = position_dodge(width = 0.9), size = 8)
su_plotf

germ_plot1 = plot_grid(s_plota, su_plotf, germ_plota, sur_plotf, nrow = 2, ncol = 2)
germ_plot1

tiff(file = "trans_germ.tiff", width = 14, height = 10, units = 'in', res = 600, pointsize = 11)

y.grob <- textGrob("Germination (%)", 
                   gp=gpar(fontface="bold", col="black", fontsize=30), rot=90)
germ_plot2 = grid.arrange(arrangeGrob(germ_plot1, left = y.grob))

dev.off()

###trans surviving seedlings

trans_all$survival = as.numeric(trans_all$survival)
trans_all$survival = (trans_all$survival)/100
trans_ac = subset(trans_all, spp == "ACSA")
trans_fa = subset(trans_all, spp == "FAGR")

glmm7 = lmer(survival ~ site * source + (1|mountain), data = trans_ac)
summary(glmm7)
anova(glmm7)
Anova(glmm7)

#trans_ac$BV <- interaction(trans_ac$site, trans_ac$source)

m1  <- glmer(survival ~ -1 + BV + (1|mountain), family = gaussian, data = trans_ac)

summary(glht(m1, linfct = mcp(BV = "Tukey")))
app = glht(m1, linfct = mcp(BV = "Tukey"))
cld(app)

#install.packages("multcomp")
library(multcomp)

summary(glht(glmm7, mcp(source = "Tukey")))

glmm8 = lmer(survival ~ site * source + (1|mountain), data = trans_fa)
summary(glmm8)
anova(glmm8)
Anova(glmm8)

#trans_ac$BV <- interaction(trans_ac$site, trans_ac$source)

m1  <- glmer(survival ~ -1 + BV + (1|mountain), family = gaussian, data = trans_fa)

summary(glht(m1, linfct = mcp(BV = "Tukey")))

trans_ga = ddply(trans_ac, ~site + source, summarise, sur = mean(survival), SE = sd(survival)/sqrt((length(survival))))

trans_gf = ddply(trans_fa, ~site + source, summarise, sur = mean(survival), SE = sd(survival)/sqrt((length(survival))))

trans_sa = ddply(trans_ac, ~site + source, summarise, sur = mean(survival), SE = sd(survival)/sqrt((length(survival))))
trans_sf = ddply(trans_fa, ~site + source, summarise, sur = mean(survival), SE = sd(survival)/sqrt((length(survival))))

trans_ga$order = factor(trans_ga$site, levels=c("L", "M", "H"))
trans_ga$order1 = factor(trans_ga$source, levels=c("L", "M", "H"),
                         labels = c("Hardwood soils", "Ecotone soils", "Spruce-fir soils"))
trans_sa$order = factor(trans_sa$site, levels=c("L", "M", "H"))
trans_ga$spp = "Sugar maple"
trans_sa$spp = "Sugar maple"
trans_sa$order1 = factor(trans_sa$source, levels=c("L", "M", "H"),
                         labels = c("Hardwood soils", "Ecotone soils", "Spruce-fir soils"))

trans_gf$order = factor(trans_gf$site, levels=c("L", "M", "H"))
trans_gf$order1 = factor(trans_gf$source, levels=c("L", "M", "H"),
                         labels = c("Hardwood soils", "Ecotone soils", "Spruce-fir soils"))
trans_gf$spp = "American beech"
trans_sf$spp = "American beech"
trans_sf$order = factor(trans_sf$site, levels=c("L", "M", "H"))
trans_sf$order1 = factor(trans_sf$source, levels=c("L", "M", "H"),
                         labels = c("Hardwood soils", "Ecotone soils", "Spruce-fir soils"))

trans_g1 = rbind(trans_ga, trans_gf)
trans_g1$spp = factor(trans_g1$spp, levels=c("Sugar maple", "American beech"))
trans_s1 = rbind(trans_sa, trans_sf)
trans_s1$spp = factor(trans_s1$spp, levels=c("Sugar maple", "American beech"))

write.csv(trans_ga, "trans_ga.csv")
trans_ga = read.csv("trans_ga.csv", header = TRUE)
write.csv(trans_sa, "trans_sa.csv")
trans_sa = read.csv("trans_sa.csv", header = TRUE)
write.csv(trans_gf, "trans_gf.csv")
trans_gf = read.csv("trans_gf.csv", header = TRUE)

#acsa
trans_s1a = subset(trans_s1, spp == "Sugar maple")
trans_s1f = subset(trans_s1, spp == "American beech")

trans_g1a = subset(trans_g1, spp == "Sugar maple")
trans_g1f = subset(trans_g1, spp == "American beech")

symbol7 = c("a", "b", "b", "b", "b", "bcd", "c", "d", "d")
sur_plota1 = ggplot(trans_s1a, aes(x = order, y = sur, fill = order)) +
  facet_grid(cols = vars(order1)) +
  geom_errorbar(aes(ymin=sur-SE, ymax=sur+SE), width = 0.2, position=position_dodge(width=0.6), colour="black", size = 1) +
  geom_point(stat = "identity", position = position_dodge(width = 0.6), pch = 23, size = 6) +
  ylab("Seedling Survival (%)") +
  xlab("Transplant Site") +
  theme_bw() +
  ylim(0,100) +
  scale_fill_manual(name = "Species", values = c("#2c7c94", "#fbe45b", "#a65852")) +
  scale_x_discrete(labels=c("L" = "Hardwoods", "M" = "Ecotone",
                            "H" = "Spruce-fir")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        strip.text.x = element_text(size = 20),
        axis.text=element_text(size=20),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(size=22),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18),
        legend.position = "none") +
  stat_summary(geom = 'text', label = symbol7, fun = max, vjust = -3, position = position_dodge(width = 0.6), size = 7)
sur_plota1

symbol27 = c("b", "c", "c", "a", "c", "c", "c", "c", "c")
sur_plotf1 = ggplot(trans_s1f, aes(x = order, y = sur, fill = order)) +
  facet_grid(cols = vars(order1)) +
  geom_errorbar(aes(ymin=sur-SE, ymax=sur+SE), width = 0.2, position=position_dodge(width=0.6), colour="black", size = 1) +
  geom_point(stat = "identity", position = position_dodge(width = 0.6), pch = 23, size = 6) +
  ylab("") +
  xlab("Transplant Site") +
  theme_bw() +
  ylim(0,100) +
  scale_fill_manual(name = "Species", values = c("#2c7c94", "#fbe45b", "#a65852")) +
  scale_x_discrete(labels=c("L" = "Hardwoods", "M" = "Ecotone",
                            "H" = "Spruce-fir")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        strip.text.x = element_text(size = 20),
        axis.text=element_text(size=20),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(size=22),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18),
        legend.position = "none") +
  stat_summary(geom = 'text', label = symbol27, fun = max, vjust = -3, position = position_dodge(width = 0.6), size = 7)
sur_plotf1

symbol8 = c("a", "b", "b")
su_plota1 = ggplot(trans_sa, aes(x = order, y = sur)) +
  geom_errorbar(aes(ymin=sur-SE, ymax=sur+SE), width = 0.075) +
  geom_bar(stat = "identity", fill = "#2c7c94", width = 0.4) +
  ylab("") +
  xlab("") +
  ylim(0,100) +
  theme_bw() +
  ggtitle("Sugar maple") +
  scale_x_discrete(labels=c("L" = "Hardwoods", "M" = "Ecotone",
                            "H" = "Spruce-fir")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        strip.text.x = element_text(size = 20),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(size=22)) +
  stat_summary(geom = 'text', label = symbol8, fun.y = max, vjust = -2.75, position = position_dodge(width = 0.9), size = 8)
su_plota1

#fagr
symbol9 = c("a", "b", "c", "b", "c", "c", "b", "c", "d")
g_plota1 = ggplot(trans_g1a, aes(x = order, y = g, fill = order)) +
  facet_grid(cols = vars(order1)) +
  geom_errorbar(aes(ymin=g-SE, ymax=g+SE), width = 0.2, position=position_dodge(width=0.6), colour="black", size = 1) +
  geom_point(stat = "identity", position = position_dodge(width = 0.6), pch = 23, size = 6) +
  ylab("Germination (%)") +
  xlab("") +
  ylim(0,35) +
  theme_bw() +
  scale_fill_manual(name = "Species", values = c("#2c7c94", "#fbe45b", "#a65852")) +
  scale_x_discrete(labels=c("L" = "Hardwoods", "M" = "Ecotone",
                            "H" = "Spruce-fir")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        strip.text.x = element_text(size = 20),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size = 20),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size=22),
        legend.position = "none") +
  stat_summary(geom = 'text', label = symbol9, fun = max, vjust = -1.5, position = position_dodge(width = 0.6), size = 8)
g_plota1

symbol29 = c("a", "b", "b", "b", "b", "b", "b", "b", "b")
g_plotf1 = ggplot(trans_g1f, aes(x = order, y = g, fill = order)) +
  facet_grid(cols = vars(order1)) +
  geom_errorbar(aes(ymin=g-SE, ymax=g+SE), width = 0.2, position=position_dodge(width=0.6), colour="black", size = 1) +
  geom_point(stat = "identity", position = position_dodge(width = 0.6), pch = 23, size = 6) +
  ylab("") +
  xlab("") +
  ylim(0,35) +
  theme_bw() +
  scale_fill_manual(name = "Transplant Site", values = c("#2c7c94", "#fbe45b", "#a65852"),
                    labels = c("L" = "Hardwoods", "M" = "Ecotone",
                               "H" = "Spruce-fir")) +
  scale_x_discrete(labels=c("L" = "Hardwoods", "M" = "Ecotone",
                            "H" = "Spruce-fir")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        strip.text.x = element_text(size = 20),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size = 20),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size=22),
        legend.title = element_text(size=22), 
        legend.text = element_text(size=22),
        legend.position = "none") +
  stat_summary(geom = 'text', label = symbol29, fun = max, vjust = -2, position = position_dodge(width = 0.6), size = 8)
g_plotf1

legend <- get_legend(
  g_plotf1 +
    guides(color = guide_legend(nrow = 3)) +
    theme(legend.position = "bottom")
)

new4 = (amf_trans_plot | emf_trans_plot) / (g_plota1 | g_plotf1) / (sur_plota1 | sur_plotf1)
new4

new5 = new4 + plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 24, hjust = -4))
new5

new1 = new + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 16, hjust = -6))
new1

new3 = g_plotf1 / sur_plota1
new3

tiff(file = "trans_sur_all1.tif", width = 15, height = 20, units = 'in', res = 300, pointsize = 11)

# Combine combined plot and legend using plot_grid()
t_sur = plot_grid(new5,
                  legend, ncol=1,rel_heights = c(3, 0.4))
t_sur

dev.off()

symbol11 = c("a", "b", "b")
si_plotf1 = ggplot(trans_sf, aes(x = order, y = sur)) +
  geom_errorbar(aes(ymin=sur-SE, ymax=sur+SE), width = 0.075) +
  geom_bar(stat = "identity", fill = "#2c7c94", width = 0.4) +
  ylab("") +
  xlab("") +
  ylim(0,100) +
  theme_bw() +
  ggtitle("American beech") +
  scale_x_discrete(labels=c("L" = "Hardwoods", "M" = "Ecotone",
                            "H" = "Spruce-fir")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        strip.text.x = element_text(size = 20),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(size=22)) +
  stat_summary(geom = 'text', label = symbol11, fun.y = max, vjust = -2, position = position_dodge(width = 0.9), size = 8)
si_plotf1

surv_plot1 = plot_grid(su_plota1, si_plotf1, sur_plota1, g_plotf1, nrow = 2, ncol = 2)
surv_plot1

tiff(file = "trans_sur.tiff", width = 14, height = 14, units = 'in', res = 600, pointsize = 11)

y.grob <- textGrob("Seedling Survival (%)", 
                   gp=gpar(fontface="bold", col="black", fontsize=30), rot=90)
surv_plot2 = grid.arrange(arrangeGrob(surv_plot1, left = y.grob))

dev.off()

tiff(file = "trans_ex.tiff", width = 16, height = 24, units = 'in', res = 600, pointsize = 11)

surv_plot12 = plot_grid(germ_plot2, surv_plot2, nrow = 2, ncol = 1)
surv_plot12

dev.off()

###trans biomass

trans_all_ac = subset(trans_all, spp == "ACSA")
trans_all_fa = subset(trans_all, spp == "FAGR")

#root:shoot

#acsa
ac_bio_root_s = ggplot(trans_all_ac, aes(x = site, y = root_shoot_ratio_avg, fill = source)) +
  geom_boxplot()
ac_bio_root_s

ac_bio_root_site = ggplot(trans_all_ac, aes(x = site, y = root_shoot_ratio_avg)) +
  geom_boxplot()
ac_bio_root_site

#fagr
fa_bio_root_s = ggplot(trans_all_fa, aes(x = site, y = root_shoot_ratio_avg, fill = source)) +
  geom_boxplot()
fa_bio_root_s

fa_bio_root_site = ggplot(trans_all_fa, aes(x = site, y = root_shoot_ratio_avg)) +
  geom_boxplot()
fa_bio_root_site

#avg total biomass

#acsa
ac_bio = ggplot(trans_all_ac, aes(x = site, y = total_bio_avg, fill = source)) +
  geom_boxplot()
ac_bio

ac_bio_site = ggplot(trans_all_ac, aes(x = site, y = total_bio_avg)) +
  geom_boxplot()
ac_bio_site

#fagr
fa_bio = ggplot(trans_all_fa, aes(x = site, y = total_bio_avg, fill = source)) +
  geom_boxplot()
fa_bio

fa_bio_site = ggplot(trans_all_fa, aes(x = site, y = total_bio_avg)) +
  geom_boxplot()
fa_bio_site

###net biomass

net_all_ac = subset(net_all, spp == "ACSA")
net_all_fa = subset(net_all, spp == "FAGR")

#root:shoot

#acsa
ac_bio_root_net = ggplot(net_all_ac, aes(x = site, y = root_shoot_ratio_avg, fill = treatment)) +
  geom_boxplot()
ac_bio_root_net

ac_bio_root_siten = ggplot(trans_all_ac, aes(x = site, y = root_shoot_ratio_avg)) +
  geom_boxplot()
ac_bio_root_siten

#fagr
fa_bio_root_net = ggplot(net_all_fa, aes(x = site, y = root_shoot_ratio_avg, fill = treatment)) +
  geom_boxplot()
fa_bio_root_net

fa_bio_root_siten = ggplot(net_all_fa, aes(x = site, y = root_shoot_ratio_avg)) +
  geom_boxplot()
fa_bio_root_siten

#avg total biomass

#acsa
ac_bion = ggplot(net_all_ac, aes(x = site, y = total_bio_avg, fill = treatment)) +
  geom_boxplot()
ac_bion

ac_bio_siten = ggplot(net_all_ac, aes(x = site, y = total_bio_avg)) +
  geom_boxplot()
ac_bio_siten

#fagr
fa_bion = ggplot(net_all_fa, aes(x = site, y = total_bio_avg, fill = treatment)) +
  geom_boxplot()
fa_bion

fa_bio_siten = ggplot(net_all_fa, aes(x = site, y = total_bio_avg)) +
  geom_boxplot()
fa_bio_siten

###RHGR - overall and for each year
#I wrangled data to calculate RHGR for each year - will likely use results for only second year
#Results for both 
#trans experiment
#year 2
trans_full$survey = as.numeric(trans_full$survey)
trans1 = subset(trans_full, survey < 7 & height_avg > 0)
trans2 = subset(trans_full, survey > 6 & height_avg > 0) 
trans2$height_avg = as.numeric(trans2$height_avg)

trans_rhgr2 = ddply(trans2, ~mountain + survey + site + spp + source, summarise, day = mean(days), height_avg1 = mean(height_avg, na.rm = TRUE), sdh = sd(height_avg, na.rm = TRUE),
                    nh = length(height_avg), seh = sdh/sqrt(nh))

trans_s_e7 = subset(trans_rhgr2, survey == 7)
trans_s_e11 = subset(trans_rhgr2, survey == 11)

trans_s_e1 = cbind(trans_s_e7, trans_s_e11$height_avg1)
trans_s_e1$day = 100
trans_s_e1$rhgr = (log(trans_s_e1$`trans_s_e11$height_avg1`) - log(trans_s_e1$height_avg1))/trans_s_e1$day

trans_s_ea1 = subset(trans_s_e1, spp == "ACSA")
trans_s_ef1 = subset(trans_s_e1, spp == "FAGR")
trans_s_ea1$rhgr = trans_s_ea1$rhgr + 0.01
trans_s_ef1$rhgr = trans_s_ef1$rhgr + 0.01

trans_s_ea1$order = factor(trans_s_ea1$site, levels=c("L", "M", "H"),
                           labels = c("Hardwoods", "Ecotone", "Spruce-fir"))
trans_s_ef1$order = factor(trans_s_ef1$site, levels=c("L", "M", "H"),
                           labels = c("Hardwoods", "Ecotone", "Spruce-fir"))
trans_s_ea1$order1 = factor(trans_s_ea1$source, levels=c("L", "M", "H"),
                            labels = c("Hardwood soils", "Ecotone soils", "Spruce-fir soils"))
trans_s_ef1$order1 = factor(trans_s_ef1$source, levels=c("L", "M", "H"),
                            labels = c("Hardwood soils", "Ecotone soils", "Spruce-fir soils"))

lmm1 = lmer(rhgr ~ site * source + (1|mountain), data = trans_s_ea1)
summary(lmm1)
anova(lmm1)
Anova(lmm1)

trans_s_ea1$BV <- interaction(trans_s_ea1$site, trans_s_ea1$source)

m1  <- glmer(rhgr ~ -1 + BV + (1|mountain), family = gaussian, data = trans_s_ea1)

summary(glht(m1, linfct = mcp(BV = "Tukey")))

lmm2 = lmer(rhgr ~ site * source + (1|mountain), data = trans_s_ef1)
summary(lmm2)
anova(lmm2)
Anova(lmm2)

trans_s_ef1$BV <- interaction(trans_s_ef1$site, trans_s_ef1$source)

m1  <- glmer(rhgr ~ -1 + BV + (1|mountain), family = gaussian, data = trans_s_ef1)

summary(glht(m1, linfct = mcp(BV = "Tukey")))

#acsa
symbol13 = c("a", "a", "a", "a", "a", "a", "a", "a", "a")
rgr_plota1 = ggplot(trans_s_ea1, aes(x = order, y = rhgr, fill = order1)) +
  geom_boxplot() +
  ylab(bquote('RHGR '(cm/day^-1))) +
  xlab("") +
  ggtitle("Sugar Maple") +
  scale_fill_brewer(palette = "Set1", direction = 1, name = "Soil Source") +
  facet_grid(~ order, scales = "free_x") +
  ylim(0.008,0.014) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        strip.text = element_text(size=12)) +
  theme(axis.text.x = element_blank(),
        axis.text=element_text(size=10),
        axis.ticks.x=element_blank(),
        axis.title=element_text(size=10)) +
  stat_summary(geom = 'text', label = symbol13, fun.y = max, vjust = -1, position = position_dodge(width = 0.75))
rgr_plota1

#fagr
rgr_plotf1 = ggplot(trans_s_ef1, aes(x = order, y = rhgr, fill = order1)) +
  geom_boxplot() +
  ylab(bquote('RHGR '(cm/day^-1))) +
  xlab("") +
  ggtitle("American Beech") +
  scale_fill_brewer(palette = "Set1", direction = 1, name = "Soil Source") +
  facet_grid(~ order, scales = "free_x") +
  ylim(0.008,0.015) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        strip.text = element_text(size=12)) +
  theme(axis.text.x = element_blank(),
        axis.text=element_text(size=10),
        axis.ticks.x=element_blank(),
        axis.title=element_text(size=10)) +
  stat_summary(geom = 'text', label = symbol13, fun.y = max, vjust = -1, position = position_dodge(width = 0.75))
rgr_plotf1

#net experiment

#year1

net_full$survey = as.numeric(net_full$survey)
net1 = subset(net_full, survey < 7 & height_avg > 0)
net2 = subset(net_full, survey > 6 & height_avg > 0)  

net_rhgr1 = ddply(net1, ~mountain + survey + site + spp + treatment, summarise, day = mean(days), height_avg1 = mean(height_avg, na.rm = TRUE), sdh = sd(height_avg, na.rm = TRUE),
                  nh = length(height_avg), seh = sdh/sqrt(nh), leaves = mean(num_leaves, na.rm = TRUE), sdl = sd(num_leaves, na.rm = TRUE),
                  nl = length(num_leaves), sel = sdl/sqrt(nl), dia = mean(diam_avg, na.rm = TRUE), sdd = sd(diam_avg, na.rm = TRUE),
                  nd = length(diam_avg), sed = sdd/sqrt(nd))

net_s_e1 = subset(net_rhgr1, survey == 1)
net_s_e6 = subset(net_rhgr1, survey == 6)

net_s_e = cbind(net_s_e1, net_s_e6$height_avg1)
net_s_e$day = 84
net_s_e$rhgr = (log(net_s_e$`net_s_e6$height_avg1`) - log(net_s_e$height_avg1))/net_s_e$day

net_s_ea = subset(net_s_e, spp == "ACSA")
net_s_ef = subset(net_s_e, spp == "FAGR")

net_s_ea$order = factor(net_s_ea$site, levels=c("L", "M", "H"))
net_s_ef$order = factor(net_s_ef$site, levels=c("L", "M", "H"))

rgr_plota = ggplot(net_s_ea, aes(x = order, y = rhgr, fill = treatment)) +
  geom_boxplot()
rgr_plota

rgr_plotf = ggplot(net_s_ef, aes(x = order, y = rhgr, fill = treatment)) +
  geom_boxplot()
rgr_plotf

#year2

net_rhgr2 = ddply(net2, ~mountain + survey + site + spp + treatment, summarise, day = mean(days), height_avg1 = mean(height_avg, na.rm = TRUE), sdh = sd(height_avg, na.rm = TRUE),
                  nh = length(height_avg), seh = sdh/sqrt(nh), leaves = mean(num_leaves, na.rm = TRUE), sdl = sd(num_leaves, na.rm = TRUE),
                  nl = length(num_leaves), sel = sdl/sqrt(nl), dia = mean(diam_avg, na.rm = TRUE), sdd = sd(diam_avg, na.rm = TRUE),
                  nd = length(diam_avg), sed = sdd/sqrt(nd))

net_s_e7 = subset(net_rhgr2, survey == 7)
net_s_e11 = subset(net_rhgr2, survey == 11)

net_s_e1 = cbind(net_s_e7, net_s_e11$height_avg1)
net_s_e1$day = 100
net_s_e1$rhgr = (log(net_s_e1$`net_s_e11$height_avg1`) - log(net_s_e1$height_avg1))/net_s_e1$day

net_s_ea1 = subset(net_s_e1, spp == "ACSA")
net_s_ef1 = subset(net_s_e1, spp == "FAGR")
net_s_ea1$rhgr = net_s_ea1$rhgr + 0.01
net_s_ef1$rhgr = net_s_ef1$rhgr + 0.01

write.csv(net_s_ef1, "net_s_ef1.csv")
net_s_ef1 = read.csv("net_s_ef1.csv", header = TRUE)

net_s_ea1$order = factor(net_s_ea$site, levels=c("L", "M", "H"))
net_s_ef1$order = factor(net_s_ef$site, levels=c("L", "M", "H"))

rgr_plota1 = ggplot(net_s_ea1, aes(x = order, y = rhgr, fill = treatment)) +
  geom_boxplot()
rgr_plota1

rgr_plotf1 = ggplot(net_s_ef1, aes(x = order, y = rhgr, fill = treatment)) +
  geom_boxplot()
rgr_plotf1

net_s_ea1$order = factor(net_s_ea1$site, levels=c("L", "M", "H"),
                         labels = c("Hardwoods", "Ecotone", "Spruce-fir"))
net_s_ef1$order = factor(net_s_ef1$site, levels=c("L", "M", "H"),
                         labels = c("Hardwoods", "Ecotone", "Spruce-fir"))
net_s_ea1$order1 = factor(net_s_ea1$treatment, levels=c("N", "D", "C"),
                          labels = c("Network", "Disturbed", "Control"))
net_s_ef1$order1 = factor(net_s_ef1$treatment, levels=c("N", "D", "C"),
                          labels = c("Network", "Disturbed", "Control"))

lmm20 = lmer(rhgr ~ site * treatment + (1|mountain), data = net_s_ea1)
summary(lmm20)
anova(lmm20)
Anova(lmm20)

net_s_ea1$BV <- interaction(net_s_ea1$site, net_s_ea1$treatment)

m1  <- glmer(rhgr ~ -1 + BV + (1|mountain), family = gaussian, data = net_s_ea1)

summary(glht(m1, linfct = mcp(BV = "Tukey")))

lmm21 = lmer(rhgr ~ site * treatment + (1|mountain), data = net_s_ef1)
summary(lmm21)
anova(lmm21)
Anova(lmm21)

net_s_ef1$BV <- interaction(net_s_ef1$site, net_s_ef1$treatment)

m1  <- glmer(rhgr ~ -1 + BV + (1|mountain), family = gaussian, data = net_s_ef1)

summary(glht(m1, linfct = mcp(BV = "Tukey")))

#acsa
symbol30 = c("a", "a", "a", "a", "a", "a", "a", "a", "a")
net_rgr_plota1 = ggplot(net_s_ea1, aes(x = order, y = rhgr, fill = order1)) +
  geom_boxplot() +
  ylab(bquote('RHGR '(cm/day^-1))) +
  xlab("") +
  ggtitle("Sugar Maple") +
  scale_fill_brewer(palette = "Set1", direction = 1, name = "Treatment") +
  facet_grid(~ order, scales = "free_x") +
  ylim(0.006,0.015) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        strip.text = element_text(size=12)) +
  theme(axis.text.x = element_blank(),
        axis.text=element_text(size=10),
        axis.ticks.x=element_blank(),
        axis.title=element_text(size=10)) +
  stat_summary(geom = 'text', label = symbol30, fun.y = max, vjust = -1, position = position_dodge(width = 0.75))
net_rgr_plota1

#fagr
net_rgr_plotf1 = ggplot(net_s_ef1, aes(x = order, y = rhgr, fill = order1)) +
  geom_boxplot() +
  ylab(bquote('RHGR '(cm/day^-1))) +
  xlab("") +
  ggtitle("American Beech") +
  scale_fill_brewer(palette = "Set1", direction = 1, name = "Treatment") +
  facet_grid(~ order, scales = "free_x") +
  ylim(0.008, 0.015) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        strip.text = element_text(size=12)) +
  theme(axis.text.x = element_blank(),
        axis.text=element_text(size=10),
        axis.ticks.x=element_blank(),
        axis.title=element_text(size=10)) +
  stat_summary(geom = 'text', label = symbol30, fun.y = max, vjust = -1, position = position_dodge(width = 0.75))
net_rgr_plotf1

###Survival

#net

net_full$survival_prop[net_full$survival_prop < 0.5 ] <- 0
net_full$survival_prop[net_full$survival_prop > 0.5 ] <- 1

net_full$survival_prop = as.numeric(net_full$survival_prop)
net_full$days = as.numeric(net_full$days)

net_full_s = ddply(net_full, ~survey + site + spp + treatment, summarise, day = mean(days),
                   survival_prop1 = sum(survival_prop)/16.4)

net_full_sa = subset(net_full_s, spp == "ACSA")
net_full_sf = subset(net_full_s, spp == "FAGR")

net_full_sa1 = ddply(net_full_sa, ~treatment + day, summarise, surv = mean(survival_prop1),
                     SE = sd(survival_prop1)/sqrt((length(survival_prop1))))

net_full_sf1 = ddply(net_full_sf, ~treatment + day, summarise, surv = mean(survival_prop1),
                     SE = sd(survival_prop1)/sqrt((length(survival_prop1))))

net_full_sitea = ddply(net_full_sa, ~site + day, summarise, surv = mean(survival_prop1),
                       SE = sd(survival_prop1)/sqrt((length(survival_prop1))))

net_full_sitef = ddply(net_full_sf, ~site + day, summarise, surv = mean(survival_prop1),
                       SE = sd(survival_prop1)/sqrt((length(survival_prop1))))

net_full_sa1$order = factor(net_full_sa1$treatment, levels=c("N", "D", "C"))
net_full_sf1$order = factor(net_full_sf1$treatment, levels=c("N", "D", "C"))

net_full_sitea$order = factor(net_full_sitea$site, levels=c("L", "M", "H"))
net_full_sitef$order = factor(net_full_sitef$site, levels=c("L", "M", "H"))


net_full_sah = subset(net_full_sa1, site == "H")
net_full_sam = subset(net_full_sa1, site == "M")
net_full_sal = subset(net_full_sa1, site == "L")

net_fsa = ddply(net_full_sa, ~treatment + spp + site + day, summarise, surv = mean(survival_prop1),
                SE = sd(survival_prop1)/sqrt((length(survival_prop1))))

net_fsf = ddply(net_full_sf, ~treatment + spp + site + day, summarise, surv = mean(survival_prop1),
                SE = sd(survival_prop1)/sqrt((length(survival_prop1))))

net_fsa$order = factor(net_fsa$treatment, levels=c("N", "D", "C"))
net_fsf$order = factor(net_fsf$treatment, levels=c("N", "D", "C"))

write.csv(net_fsa, "net_fsa.csv")
write.csv(net_fsf, "net_fsf.csv")

net_fsa$order1 = factor(net_fsa$site, levels=c("L", "M", "H"), labels = c("Hardwoods",
                                                                          "Ecotone", "Spruce-fir"))
net_fsf$order1 = factor(net_fsf$site, levels=c("L", "M", "H"), labels = c("Hardwoods",
                                                                          "Ecotone", "Spruce-fir"))

net_fsa = read.csv("net_fsa.csv", header = TRUE)
net_fsf = read.csv("net_fsf.csv", header = TRUE)

surv_plot_a = ggplot(net_fsa, aes(x = day, y = surv, fill = order)) +
  geom_pointrange(aes(ymin=surv-SE, ymax=surv+SE, color=factor(order)), lwd = 1, show.legend = FALSE) +
  geom_line(stat = "identity", aes(color = factor(order)), lwd = 1, show.legend = FALSE) +
  geom_point(stat = "identity", aes(fill = order), size = 6, pch = 23) +
  scale_color_manual(values=c("#fbe45b", "#2c7c94", "#a65852"))+
  scale_fill_manual(values = c("#fbe45b", "#2c7c94", "#a65852"), name = "Treatment", labels = c("Network", "Disrupted", "Control")) +
  xlab("") +
  ylab("") +
  ggtitle("Sugar maple - Treatment * Elevation") +
  facet_grid(~order1, labeller = labeller(c("Hardwoods", "Ecotone", "Spruce-fir"))) +
  xlim(0,500) +
  ylim(0,1) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        plot.title = element_text(size=18),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        strip.text.x = element_text(size = 16)) +
  geom_rect(aes(xmin = 140, xmax = 340, ymin = -Inf, ymax = Inf),
            fill = "grey", alpha = 0.03)
surv_plot_a

surv_plot_f = ggplot(net_fsf, aes(x = day, y = surv, fill = order)) +
  geom_pointrange(aes(ymin=surv-SE, ymax=surv+SE, color=factor(order)), lwd = 1, show.legend = FALSE) +
  geom_line(stat = "identity", aes(color = factor(order)), lwd = 1, show.legend = FALSE) +
  geom_point(stat = "identity", aes(fill = order), size = 6, pch = 23) +
  scale_color_manual(values=c("#fbe45b", "#2c7c94", "#a65852"))+
  scale_fill_manual(values = c("#fbe45b", "#2c7c94", "#a65852"), name = "Treatment", labels = c("Network", "Disrupted", "Control")) +
  xlab("") +
  ylab("") +
  ggtitle("American beech - Treatment * Elevation") +
  facet_grid(~order1, labeller = labeller(c("Hardwoods", "Ecotone", "Spruce-fir"))) +
  xlim(0,500) +
  ylim(0,1) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        legend.position = "none",
        plot.title = element_text(size=18),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        strip.text.x = element_text(size = 16)) +
  geom_rect(aes(xmin = 140, xmax = 340, ymin = -Inf, ymax = Inf),
            fill = "grey", alpha = 0.03)
surv_plot_f

#acsa
surv_plot = ggplot(net_full_sa1, aes(x = day, y = surv, fill = order)) +
  geom_pointrange(aes(ymin=surv-SE, ymax=surv+SE, color=factor(order)), lwd = 1) +
  geom_line(stat = "identity", aes(color = factor(order)), lwd = 1) +
  geom_point(stat = "identity", aes(fill = order), size = 6, pch = 23) +
  scale_color_manual(values=c("#fbe45b", "#2c7c94", "#a65852"))+
  scale_fill_manual(values = c("#fbe45b", "#2c7c94", "#a65852"), name = "Treatment", labels = c("Network (N)", "Disrupted (D)", "Control (C)")) +
  xlab("") +
  ylab("") +
  ggtitle("Sugar Maple - Main Effects") +
  xlim(0,500) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=20),
        legend.position = "none",
        axis.title=element_text(size=20),
        plot.title = element_text(size=18)) +
  geom_rect(aes(xmin = 140, xmax = 340, ymin = -Inf, ymax = Inf),
            fill = "grey", alpha = 0.03)
surv_plot

#fagr
surv_plotf = ggplot(net_full_sf1, aes(x = day, y = surv, fill = order)) +
  geom_pointrange(aes(ymin=surv-SE, ymax=surv+SE, color=factor(order)), lwd = 1) +
  geom_line(stat = "identity", aes(color = factor(order)), lwd = 1) +
  geom_point(stat = "identity", color = "black", aes(color = order), size = 6, pch = 23) +
  scale_color_manual(values=c("#fbe45b", "#2c7c94", "#a65852")) +
  scale_fill_manual(values = c("#fbe45b", "#2c7c94", "#a65852"), name = "Treatment", labels = c("Network (N)", "Disrupted (D)", "Pot Control (C)")) +
  xlab("") +
  ylab("") +
  ggtitle("American Beech - Main Effects") +
  xlim(0,500) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=20),
        legend.position = c(0.84,0.8),
        axis.title=element_text(size=20),
        legend.title=element_text(size=18), 
        legend.text=element_text(size=16),
        plot.title = element_text(size=18)) +
  geom_rect(aes(xmin = 140, xmax = 340, ymin = -Inf, ymax = Inf),
            fill = "grey", alpha = 0.03) +
  guides(col = FALSE) 
surv_plotf 

surv_plot_net = plot_grid(surv_plot, surv_plot_a, surv_plotf, surv_plot_f, ncol = 2, nrow = 2, labels = c("a", "b", "c", "d"), label_size = 20)
surv_plot_net

tiff(file = "net_surv2.tif", width = 15, height = 9, units = 'in', res = 600, pointsize = 11)
y.grob <- textGrob("Proportion of Seedlings Surviving", 
                   gp=gpar(fontface="bold", col="black", fontsize=24), rot=90)
x.grob <- textGrob("Number of Days Since Planting", 
                   gp=gpar(fontface="bold", col="black", fontsize=24))
surv_plot_net1 = grid.arrange(arrangeGrob(surv_plot_net, left = y.grob, bottom = x.grob))
dev.off()

#by site
#acsa
surv_plotsite = ggplot(net_full_sitea, aes(x = day, y = surv, fill = order)) +
  geom_pointrange(aes(ymin=surv-SE, ymax=surv+SE, color=factor(order)), lwd = 1) +
  geom_line(stat = "identity", aes(color = factor(order)), lwd = 1) +
  geom_point(stat = "identity", aes(fill = order), size = 6, pch = 23) +
  scale_color_manual(values=c("#fbe45b", "#2c7c94", "#a65852"))+
  scale_fill_manual(values = c("#fbe45b", "#2c7c94", "#a65852"), name = "Elevation", labels = c("Hardwoods", "Ecotone", "Spruce-fir")) +
  xlab("") +
  ylab("") +
  ggtitle("Sugar Maple") +
  xlim(0,500) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=20),
        legend.position = "none",
        axis.title=element_text(size=20),
        plot.title = element_text(size=22)) +
  geom_rect(aes(xmin = 140, xmax = 340, ymin = -Inf, ymax = Inf),
            fill = "grey", alpha = 0.03)
surv_plotsite

#fagr
surv_plotfsite = ggplot(net_full_sitef, aes(x = day, y = surv, fill = order)) +
  geom_pointrange(aes(ymin=surv-SE, ymax=surv+SE, color=factor(order)), lwd = 1) +
  geom_line(stat = "identity", aes(color = factor(order)), lwd = 1) +
  geom_point(stat = "identity", color = "black", aes(color = order), size = 6, pch = 23) +
  scale_color_manual(values=c("#fbe45b", "#2c7c94", "#a65852")) +
  scale_fill_manual(values = c("#fbe45b", "#2c7c94", "#a65852"), name = "Elevation", labels = c("Harwoods", "Ecotone", "Spruce-fir")) +
  xlab("Number of Days Since Planting") +
  ylab("") +
  ggtitle("American Beech") +
  xlim(0,500) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=20),
        legend.position = c(0.84,0.8),
        axis.title=element_text(size=20),
        legend.title=element_text(size=18), 
        legend.text=element_text(size=16),
        plot.title = element_text(size=22)) +
  geom_rect(aes(xmin = 140, xmax = 340, ymin = -Inf, ymax = Inf),
            fill = "grey", alpha = 0.03) +
  guides(col = FALSE) 
surv_plotfsite 

surv_plot_site = plot_grid(surv_plotsite, surv_plotfsite, ncol = 1, nrow = 2)
surv_plot_site

tiff(file = "net_survsite.tif", width = 8, height = 9, units = 'in', res = 600, pointsize = 11)
y.grob <- textGrob("Proportion of Seedlings Surviving", 
                   gp=gpar(fontface="bold", col="black", fontsize=24), rot=90)
surv_plot_netsite1 = grid.arrange(arrangeGrob(surv_plot_site, left = y.grob))
dev.off()

net_full_sitea = ddply(net_full_sa, ~site + day, summarise, surv = mean(survival_prop1),
                       SE = sd(survival_prop1)/sqrt((length(survival_prop1))))

surv_plot1 = ggplot(net_full_sitea, aes(x = day, y = surv, group = site, color = site)) +
  geom_errorbar(aes(ymin=surv-SE, ymax=surv+SE), width = 0.075) +
  geom_line() +
  geom_point(stat = "identity", size = 4) +
  xlim(0,500)
surv_plot1   

net_full_sitea = ddply(net_full_sa, ~site + day, summarise, surv = mean(survival_prop1),
                       SE = sd(survival_prop1)/sqrt((length(survival_prop1))))

surv_plot2 = ggplot(net_full_sitea, aes(x = day, y = surv, group = site, color = site)) +
  geom_errorbar(aes(ymin=surv-SE, ymax=surv+SE), width = 0.075) +
  geom_line() +
  geom_point(stat = "identity", size = 4) +
  xlab("Number of Days Since Planting") +
  ylab("Proportion of Pots Surviving") +
  xlim(0,500) +
  scale_fill_brewer(palette = "Greens", direction = -1, name = "Treatment", labels = c("Pot Control", "Disrupted", "Network")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        plot.title = element_text(size=18))
surv_plot2

net_full_sitef = ddply(net_full_sf, ~site + day, summarise, surv = mean(survival_prop1),
                       SE = sd(survival_prop1)/sqrt((length(survival_prop1))))

surv_plot3 = ggplot(net_full_sitef, aes(x = day, y = surv, group = site, color = site)) +
  geom_errorbar(aes(ymin=surv-SE, ymax=surv+SE), width = 0.075) +
  geom_line() +
  geom_point(stat = "identity", size = 4) +
  xlab("Number of Days Since Planting") +
  ylab("Proportion of Pots Surviving") +
  xlim(0,500) +
  scale_fill_brewer(palette = "Greens", direction = -1, name = "Treatment", labels = c("Pot Control", "Disrupted", "Network")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        plot.title = element_text(size=18))
surv_plot3

net_year1 = subset(net_full, survey < 7)
net_year2 = subset(net_full, survey > 6)

net_year2$survival_prop = as.numeric(net_year2$survival_prop)
net_year2$days = as.numeric(net_year2$days)

#survival analysis
Surf = Surv(time = net_year2$days, event = net_year2$survival_prop, type = "left")
KM = survfit(Surf ~ treatment, data = net_year2, type="kaplan-meier")

ggsurvplot(KM, conf.int = TRUE, pval = FALSE, risk.table = FALSE, 
           legend = "right", censor.shape = "", censor.size = 4, 
           ylab = "Proportion surviving", xlab = "Time (days)",
           xlim(305,500))

plot(survfit(Surv(days, survival_prop) ~ site + treatment, data = net_year1), 
     xlab = "Days", 
     ylab = "Overall survival probability")

#
#install.packages("ggfortify")
#install.packages("ranger")
library(ggfortify)
library(ranger)

net_year2a = subset(net_year2, spp == "ACSA")

net_year2a$survival_prop[net_year2a$survival_prop == 1] = 2
net_year2a$survival_prop[net_year2a$survival_prop == 0] = 1
net_year2a$survival_prop[net_year2a$survival_prop == 2] = 0

tiff(file = "KM_surv.tif", width = 9, height = 6, units = 'in', res = 600, pointsize = 11)

km_fit <- survfit(Surv(days, survival_prop) ~ treatment, data=net_year2a)
summary(km_fit, times = c(1,5*(1:100)))
autoplot(km_fit, surv.size = 1.5, palette = c("#fbe45b", "#2c7c94", "#a65852")) + xlim(350,500) +
  xlab("Days Since Planting") +
  ylab("Survival Probability") +
  ggtitle("Sugar Maple 2021") +
  scale_fill_manual(values = c("#fbe45b", "#2c7c94", "#a65852"), name = "Treatment", labels = c("Networked (N)", "Disrupted (D)", "Control (C)")) +
  scale_color_manual(values = c("#fbe45b", "#2c7c94", "#a65852"), name = "Treatment", labels = c("Networked (N)", "Disrupted (D)", "Control (C)")) +
  theme_bw() +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=18),
        plot.title = element_text(size=24)) +
  guides(col = FALSE)

dev.off()

#trans

trans_alla = subset(trans_all, spp == "ACSA")
trans_allf = subset(trans_all, spp == "FAGR")

lma = aov(num_alive ~ site, data = trans_alla)
summary(lma)
TukeyHSD(lma)

lmb = aov(num_alive ~ site, data = trans_allf)
summary(lmb)
TukeyHSD(lmb)

trans_alla$order = factor(trans_alla$site, levels=c("L", "M", "H"))
trans_allf$order = factor(trans_allf$site, levels=c("L", "M", "H"))

trans_surva = ggplot(trans_alla, aes(x = order, y = num_alive)) +
  geom_boxplot()
trans_surva

trans_survab = ggplot(trans_alla, aes(x = order, y = num_alive, fill = treatment)) +
  geom_boxplot()
trans_survab

trans_survf = ggplot(trans_allf, aes(x = order, y = num_alive)) +
  geom_boxplot()
trans_survf

trans_survfb = ggplot(trans_allf, aes(x = order, y = num_alive, fill = treatment)) +
  geom_boxplot()
trans_survfb


###Browse, herbivory, pathogens

#net
net_all1 = ddply(net_all, ~mountain + site + spp + treatment, summarise, browse = sum(browse), herb = sum(herb), root_path = sum(root_path))
net_all1_acsa = subset(net_all1, spp == "ACSA")
net_all1_fagr = subset(net_all1, spp == "FAGR")

net_all1$browse = (net_all1$browse/4)*100
net_all1$herb = (net_all1$herb/4)*100
net_all1$root_path = (net_all1$root_path/4)*100

net_all1_fagr$browse = (net_all1_fagr$browse/16)*100
net_all1_fagr$herb = (net_all1_fagr$herb/16)*100
net_all1_fagr$root_path = (net_all1_fagr$root_path/16)*4

root_path2$order = factor(root_path2$site, levels=c("L", "M", "H"))
net_all1$order = factor(net_all1$site, levels=c("L", "M", "H"))

root_path1 = ddply(net_all1_acsa, ~order + treatment, summarise, path_m = mean(root_path), sd = sd(root_path),
                   n = length(root_path), se = (sd/sqrt(n)) - 1)
root_path2 = ddply(net_all1_fagr, ~order + treatment, summarise, path_m = mean(root_path), sd = sd(root_path),
                   n = length(root_path), se = (sd/sqrt(n)) - 1)
root_path2$se = abs(root_path2$se)
net_all1_acsa1$order = factor(net_all1_acsa1$site, levels=c("L", "M", "H"))
net_all1_fagr$order = factor(net_all1_fagr$site, levels=c("L", "M", "H"))

net_all1 = ddply(net_all1, ~site + spp, summarise, browse_m = mean(browse), sdb = sd(browse), nb = length(browse), seb = sdb/sqrt(nb),
                 herb_m = mean(herb), sdh = sd(herb), nh = length(herb), seh = sdh/sqrt(nh),
                 path_m = mean(root_path), sdp = sd(root_path), np = length(root_path), sep = sdp/sqrt(np))

lm1 = aov(browse ~ site, data = net_all1_acsa)
summary(lm1)
TukeyHSD(lm1)

lm2 = aov(herb ~ site, data = net_all1_acsa)
summary(lm2)
TukeyHSD(lm2)

#acsa
symbol = c("a", "b", "b")
browse_net = ggplot(net_all1_acsa1, aes(x = order, y = browse_m)) +
  geom_errorbar(aes(ymin=browse_m-seb, ymax=browse_m+seb), width = 0.075) +
  geom_bar(stat = "identity", position = "dodge") +
  ylab("Browse (% of pots)") +
  xlab("") +
  ggtitle("Network Experiment - Sugar maple") +
  ylim(0,50) +
  theme_bw() +
  scale_x_discrete(labels=c("L" = "Northern Hardwoods", "M" = "Ecotone",
                            "H" = "Spruce-fir")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        plot.title = element_text(size=18)) +
  stat_summary(geom = 'text', label = symbol, fun = max, vjust = -4.5, position = position_dodge(width = 0.75), size = 6)
browse_net

tiff(file = "net_path.tif", width = 8, height = 5, units = 'in', res = 600, pointsize = 11)

root_path1$order1 = factor(root_path1$treatment, levels=c("N", "D", "C"))
symbol6 = c("b", "a", "b", "b", "a", "b", "b", "a", "b")

path_net_full = ggplot(root_path1, aes(x = order, y = path_m, fill = order1)) +
  geom_errorbar(aes(ymin=path_m-se, ymax=path_m+se), position=position_dodge(width=0.9), width = 0.15) +
  geom_bar(stat = "identity", position = "dodge") +
  ylab("Root Pathogen Presence (% of pots)") +
  xlab("") +
  ylim(0,15) +
  ggtitle("Sugar maple") +
  scale_fill_manual(values = c("#fbe45b", "#2c7c94", "#a65852"), name = "Treatment",
                    label = c("Network", "Disrupted", "Pot Control")) +
  theme_bw() +
  scale_x_discrete(labels=c("L" = "Hardwoods", "M" = "Ecotone",
                            "H" = "Spruce-fir")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        plot.title = element_text(size=18)) +
  theme(legend.text=element_text(size=14),
        legend.title=element_text(size=14),
        legend.position = "none") +
  stat_summary(geom = 'text', label = symbol6, fun = max, vjust = -4, position = position_dodge(width = 0.9), size = 6)
path_net_full

root_path2$order1 = factor(root_path2$treatment, levels=c("N", "D", "C"))
symbol7 = c("a", "a", "a", "a", "a", "a", "a", "a", "a")

path_net_f = ggplot(root_path2, aes(x = order, y = path_m, fill = order1)) +
  geom_errorbar(aes(ymin=path_m-se, ymax=path_m+se), position=position_dodge(width=0.9), width = 0.15) +
  geom_bar(stat = "identity", position = "dodge") +
  ylab("") +
  xlab("") +
  ylim(0,15) +
  ggtitle("American beech") +
  scale_fill_manual(values = c("#fbe45b", "#2c7c94", "#a65852"), name = "Treatment",
                    label = c("Network", "Disrupted", "Pot Control")) +
  theme_bw() +
  scale_x_discrete(labels=c("L" = "Hardwoods", "M" = "Ecotone",
                            "H" = "Spruce-fir")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        plot.title = element_text(size=18)) +
  theme(legend.text=element_text(size=14),
        legend.title=element_text(size=14),
        legend.position = "none") +
  stat_summary(geom = 'text', label = symbol7, fun = max, vjust = -5, position = position_dodge(width = 0.9), size = 6)
path_net_f

new4 = path_net_full + path_net_f + plot_annotation(tag_levels = 'a')
new4 + theme(plot.tag = element_text(size = 28))
new4

legend <- get_legend(
  path_net_full +
    guides(color = guide_legend(nrow = 3)) +
    theme(legend.position = "bottom")
)

tiff(file = "path_grow1.tif", width = 10, height = 6, units = 'in', res = 600, pointsize = 11)

# Combine combined plot and legend using plot_grid()
path_all7 = plot_grid(new4,
                      legend, ncol=1,rel_heights = c(3, 0.4))
path_all7

dev.off()


path_net = ggplot(net_all1_acsa1, aes(x = order, y = path_m)) +
  geom_errorbar(aes(ymin=path_m-sep, ymax=path_m+sep), width = 0.075) +
  geom_bar(stat = "identity", position = "dodge") +
  ylab("Root Pathogen Infection (% of pots)") +
  xlab("") +
  ylim(0, 40) +
  ggtitle("Network Experiment - Sugar maple") +
  theme_bw() +
  scale_x_discrete(labels=c("L" = "Northern Hardwoods", "M" = "Ecotone",
                            "H" = "Spruce-fir")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        plot.title = element_text(size=18)) +
  stat_summary(geom = 'text', label = symbol4, fun = max, vjust = -10.5, position = position_dodge(width = 0.75), size = 6)
path_net

symbol1 = c("a", "ab", "b")
herb_net = ggplot(net_all1_acsa1, aes(x = order, y = herb_m)) +
  geom_errorbar(aes(ymin=herb_m-seh, ymax=herb_m+seh), width = 0.075) +
  geom_bar(stat = "identity", position = "dodge") +
  ylab("Leaf Herbivory (% of pots)") +
  xlab("") +
  ylim(0,100) +
  ggtitle("Network Experiment - Sugar maple") +
  theme_bw() +
  scale_x_discrete(labels=c("L" = "Northern Hardwoods", "M" = "Ecotone",
                            "H" = "Spruce-fir")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        plot.title = element_text(size=18)) +
  stat_summary(geom = 'text', label = symbol1, fun = max, vjust = -3, position = position_dodge(width = 0.75), size = 6)
herb_net

net_all_a = subset(net_all, spp == "ACSA")
net_all_f = subset(net_all, spp == "FAGR")

#net_all_a$site[net_all_a$site == "L"] <- 600
#net_all_a$site[net_all_a$site == "M"] <- 800
#net_all_a$site[net_all_a$site == "H"] <- 1000

#net_all_a$site = as.numeric(net_all_a$site)

browse_null = glmer(browse ~ site * treatment + (1|mountain), data = net_all_a, family="binomial", na.action=na.omit)
summary(browse_null)
anova(browse_null)
Anova(browse_null)

net_all_a$BV <- interaction(net_all_a$site, net_all_a$treatment)

m1  <- glmer(browse ~ -1 + BV + (1|mountain), family = gaussian, data = net_all_a)

summary(glht(m1, linfct = mcp(BV = "Tukey")))

browse_null = glmer(browse ~ site * treatment + (1|mountain), data = net_all_f, family="binomial", na.action=na.omit)
summary(browse_null)
anova(browse_null)
Anova(browse_null)

net_all_f$BV <- interaction(net_all_f$site, net_all_f$treatment)

m1  <- glmer(browse ~ -1 + BV + (1|mountain), family = gaussian, data = net_all_f)

summary(glht(m1, linfct = mcp(BV = "Tukey")))

root_null = glmer(root_path ~ site * treatment + (1|mountain), data = net_all_a, family="binomial", na.action=na.omit)
summary(root_null)
anova(root_null)
Anova(root_null)

#net_all_f$BV <- interaction(net_all_f$site, net_all_f$treatment)

m1  <- glmer(root_path ~ -1 + BV + (1|mountain), family = binomial, data = net_all_a)

summary(glht(m1, linfct = mcp(BV = "Tukey")))

root_null = glmer(root_path ~ site * treatment + (1|mountain), data = net_all_f, family="binomial", na.action=na.omit)
summary(root_null)
anova(root_null)
Anova(root_null)

#net_all_f$BV <- interaction(net_all_f$site, net_all_f$treatment)

m1  <- glmer(root_path ~ -1 + BV + (1|mountain), family = binomial, data = net_all_f)

summary(glht(m1, linfct = mcp(BV = "Tukey")))

lm4 = aov(root_path ~ site + treatment, data = net_all_a)
summary(lm4)
TukeyHSD(lm4)

herb_null = glmer(herb ~ site * treatment + (1|mountain), data = net_all_a, family="binomial", na.action=na.omit)
summary(herb_null)
anova(herb_null)
Anova(herb_null)

#net_all_f$BV <- interaction(net_all_f$site, net_all_f$treatment)

m1  <- glmer(herb ~ -1 + BV + (1|mountain), family = binomial, data = net_all_a)

summary(glht(m1, linfct = mcp(BV = "Tukey")))

herb_null = glmer(herb ~ site * treatment + (1|mountain), data = net_all_f, family="binomial", na.action=na.omit)
summary(herb_null)
anova(herb_null)
Anova(herb_null)

#net_all_f$BV <- interaction(net_all_f$site, net_all_f$treatment)

m1  <- glmer(herb ~ -1 + BV + (1|mountain), family = binomial, data = net_all_f)

summary(glht(m1, linfct = mcp(BV = "Tukey")))

xtabs(~browse + site, data = net_all_a)
logit_1 <- glm(browse~site, family = binomial, data = net_all_a)
summary
logit_2 <- glm(herb~site, family = binomial, data = net_all_a)
summary(logit_2)
logit_3 <- glm(root_path~site+treatment, family = binomial, data = net_all_a)
summary(logit_3)

#fagr
symbol2 = c("a", "a", "b")
net_all1_fagr1 = ddply(net_all1_fagr, ~order, summarise, browse_m = mean(browse), sdb = sd(browse), nb = length(browse), seb = (sdb/sqrt(nb))-0.1,
                       herb_m = mean(herb), sdh = sd(herb), nh = length(herb), seh = sdh/sqrt(nh), root_path_m = mean(root_path),
                       sdp = sd(root_path), np = length(root_path), sep = sdp/sqrt(np))

browse_net1 = ggplot(net_all1_fagr1, aes(x = order, y = browse_m)) +
  geom_errorbar(aes(ymin=browse_m-seb, ymax=browse_m+seb), width = 0.075) +
  geom_bar(stat = "identity", position = "dodge") +
  ylab("Browse (% of pots)") +
  xlab("") +
  ylim(0, 20) +
  ggtitle("Network Experiment - American beech") +
  theme_bw() +
  scale_x_discrete(labels=c("L" = "Northern Hardwoods", "M" = "Ecotone",
                            "H" = "Spruce-fir")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        plot.title = element_text(size=18)) +
  stat_summary(geom = 'text', label = symbol2, fun = max, vjust = -4, position = position_dodge(width = 0.75), size = 6)
browse_net1

symbol4 = c("a", "a", "a")
path_net1 = ggplot(net_all1_fagr1, aes(x = order, y = root_path_m)) +
  geom_errorbar(aes(ymin=root_path_m-sep, ymax=root_path_m+sep), width = 0.075) +
  geom_bar(stat = "identity", position = "dodge") +
  ylab("Root Pathogen Infection (% of pots)") +
  xlab("") +
  ylim(0,40) +
  ggtitle("Network Experiment - American beech") +
  theme_bw() +
  scale_x_discrete(labels=c("L" = "Northern Hardwoods", "M" = "Ecotone",
                            "H" = "Spruce-fir")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        plot.title = element_text(size=18)) +
  stat_summary(geom = 'text', label = symbol4, fun = max, vjust = -4, position = position_dodge(width = 0.75), size = 6)
path_net1

herb_net1 = ggplot(net_all1_fagr1, aes(x = order, y = herb_m)) +
  geom_errorbar(aes(ymin=herb_m-seh, ymax=herb_m+seh), width = 0.075) +
  geom_bar(stat = "identity", position = "dodge") +
  ylab("Leaf Herbivory (% of pots)") +
  xlab("") +
  ylim(0,70) +
  ggtitle("Network Experiment - American beech") +
  theme_bw() +
  scale_x_discrete(labels=c("L" = "Northern Hardwoods", "M" = "Ecotone",
                            "H" = "Spruce-fir")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        plot.title = element_text(size=18)) +
  stat_summary(geom = 'text', label = symbol4, fun = max, vjust = -4.5, position = position_dodge(width = 0.75), size = 6)
herb_net1

browse_nullf = glmer(browse ~ site + (1|net_all_f$mountain), data = net_all_f, family="binomial", na.action=na.omit)
summary(browse_nullf)

root_nullf = glmer(root_path ~ site + treatment + (1|net_all_f$mountain), data = net_all_f, family="binomial", na.action=na.omit)
summary(root_nullf)

herb_nullf = glmer(herb ~ site + (1|net_all_f$mountain), data = net_all_f, family="binomial", na.action=na.omit)
summary(herb_nullf)

tiff(file = "net_bio.tif", width = 18, height = 12, units = 'in', res = 600, pointsize = 11)

net_plots_bio = plot_grid(browse_net, path_net, herb_net, 
                          browse_net1, path_net1, herb_net1, nrow = 2, ncol = 3,
                          labels = c("A", "B", "C", "D", "E", "F"))
net_plots_bio

dev.off()

#trans

trans_all1 = ddply(trans_all, ~site + spp + treatment, summarise, browse = sum(browse), herb = sum(herb), root_path = sum(root_path))
trans_all1_acsa = subset(trans_all1, spp == "ACSA")
trans_all1_fagr = subset(trans_all1, spp == "FAGR")

trans_all_a = subset(trans_all, spp == "ACSA")
trans_all_f = subset(trans_all, spp == "FAGR")

trans_all1$browse = (trans_all1$browse/12)*100
trans_all1$herb = (trans_all1$herb/12)*100
trans_all1$root_path = (trans_all1$root_path/12)*100

trans_all1_fagr$browse = (trans_all1_fagr$browse/12)*100
trans_all1_fagr$herb = (trans_all1_fagr$herb/12)*100
trans_all1_fagr$root_path = (trans_all1_fagr$root_path/12)*100

trans_all1$order = factor(trans_all1$site, levels=c("L", "M", "H"))
trans_all1_fagr$order = factor(trans_all1_fagr$site, levels=c("L", "M", "H"))

trans_all1 = ddply(trans_all1, ~order + spp, summarise, browse1 = mean(browse), sdb = sd(browse), nb = length(browse), seb = sdb/sqrt(nb),
                   root_path1 = mean(root_path), sdp = sd(root_path), np = length(root_path), sep = sdp/sqrt(np),
                   herb1 = mean(herb), sdh = sd(herb), nh = length(herb), seh = sdh/sqrt(nh))

trans_all1_fagr1 = ddply(trans_all1_fagr, ~order, summarise, browse1 = mean(browse), sdb = sd(browse), nb = length(browse), seb = sdb/sqrt(nb),
                         root_path1 = mean(root_path), sdp = sd(root_path), np = length(root_path), sep = sdp/sqrt(np),
                         herb1 = mean(herb), sdh = sd(herb), nh = length(herb), seh = sdh/sqrt(nh))

#acsa

browse_trans = glmer(browse ~ site * treatment + (1|trans_all_a$mountain), family = "binomial", data = trans_all_a)
summary(browse_trans)
anova(browse_trans)
Anova(browse_trans)

trans_all_a$BV <- interaction(trans_all_a$site, trans_all_a$treatment)

m1  <- glmer(browse ~ -1 + BV + (1|mountain), family = binomial, data = trans_all_a)

summary(glht(m1, linfct = mcp(BV = "Tukey")))


browse_trans = glmer(browse ~ site * treatment + (1|trans_all_f$mountain), family = "binomial", data = trans_all_f)
summary(browse_trans)
anova(browse_trans)
Anova(browse_trans)

trans_all_f$BV <- interaction(trans_all_f$site, trans_all_f$treatment)

m1  <- glmer(browse ~ -1 + BV + (1|mountain), family = binomial, data = trans_all_f)

summary(glht(m1, linfct = mcp(BV = "Tukey")))

root_trans = glmer(root_path ~ site * treatment + (1|trans_all_a$mountain), data = trans_all_a, family="binomial", na.action=na.omit)
summary(root_trans)
anova(root_trans)
Anova(root_trans)

#trans_all_f$BV <- interaction(trans_all_f$site, trans_all_f$treatment)

m1  <- glmer(root_path ~ -1 + BV + (1|mountain), family = binomial, data = trans_all_a)

summary(glht(m1, linfct = mcp(BV = "Tukey")))

root_trans = glmer(root_path ~ site * treatment + (1|trans_all_f$mountain), data = trans_all_f, family="binomial", na.action=na.omit)
summary(root_trans)
anova(root_trans)
Anova(root_trans)

#trans_all_f$BV <- interaction(trans_all_f$site, trans_all_f$treatment)

m1  <- glmer(root_path ~ -1 + BV + (1|mountain), family = binomial, data = trans_all_f)

summary(glht(m1, linfct = mcp(BV = "Tukey")))

herb_trans = glmer(herb ~ site * treatment + (1|trans_all_a$mountain), data = trans_all_a, family="binomial", na.action=na.omit)
summary(herb_trans)
anova(herb_trans)
Anova(herb_trans)

#trans_all_f$BV <- interaction(trans_all_f$site, trans_all_f$treatment)

m1  <- glmer(herb ~ -1 + BV + (1|mountain), family = binomial, data = trans_all_a)

summary(glht(m1, linfct = mcp(BV = "Tukey")))

herb_trans = glmer(herb ~ site * treatment + (1|trans_all_f$mountain), data = trans_all_f, family="binomial", na.action=na.omit)
summary(herb_trans)
anova(herb_trans)
Anova(herb_trans)

#trans_all_f$BV <- interaction(trans_all_f$site, trans_all_f$treatment)

m1  <- glmer(herb ~ -1 + BV + (1|mountain), family = binomial, data = trans_all_f)

summary(glht(m1, linfct = mcp(BV = "Tukey")))


browse_trans = ggplot(trans_all1_acsa1, aes(x = order, y = browse1)) +
  geom_errorbar(aes(ymin=browse1-seb, ymax=browse1+seb), width = 0.075) +
  geom_bar(stat = "identity", position = "dodge") +
  ylab("Browse (% of pots)") +
  xlab("") +
  ylim(0, 15) +
  ggtitle("Transplant Experiment - Sugar maple") +
  theme_bw() +
  scale_x_discrete(labels=c("L" = "Northern Hardwoods", "M" = "Ecotone",
                            "H" = "Spruce-fir")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        plot.title = element_text(size=18)) +
  stat_summary(geom = 'text', label = symbol4, fun = max, vjust = -8, position = position_dodge(width = 0.75), size = 6)
browse_trans

path_trans = ggplot(trans_all1_acsa1, aes(x = order, y = root_path1)) +
  geom_errorbar(aes(ymin=root_path1-sep, ymax=root_path1+sep), width = 0.075) +
  geom_bar(stat = "identity", position = "dodge") +
  ylab("Root Pathogen Infection (% of pots)") +
  xlab("") +
  ylim(0,20) +
  ggtitle("Transplant Experiment - Sugar maple") +
  theme_bw() +
  scale_x_discrete(labels=c("L" = "Northern Hardwoods", "M" = "Ecotone",
                            "H" = "Spruce-fir")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        plot.title = element_text(size=18)) +
  stat_summary(geom = 'text', label = symbol4, fun = max, vjust = -5.5, position = position_dodge(width = 0.75), size = 6)
path_trans

symbol5 = c("a", "b", "b")
herb_trans = ggplot(trans_all1_acsa1, aes(x = order, y = herb1)) +
  geom_errorbar(aes(ymin=herb1-seh, ymax=herb1+seh), width = 0.075) +
  geom_bar(stat = "identity", position = "dodge") +
  ylab("Leaf Herbivory (% of pots)") +
  xlab("") +
  ylim(0,80) +
  ggtitle("Transplant Experiment - Sugar maple") +
  theme_bw() +
  scale_x_discrete(labels=c("L" = "Northern Hardwoods", "M" = "Ecotone",
                            "H" = "Spruce-fir")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        plot.title = element_text(size=18)) +
  stat_summary(geom = 'text', label = symbol5, fun = max, vjust = -3, position = position_dodge(width = 0.75), size = 6)
herb_trans

#fagr
browse_trans1 = ggplot(trans_all1_fagr1, aes(x = order, y = browse1)) +
  geom_errorbar(aes(ymin=browse1-seb, ymax=browse1+seb), width = 0.075) +
  geom_bar(stat = "identity", position = "dodge") +
  ylab("Browse (% of pots)") +
  xlab("") +
  ylim(0,15) +
  ggtitle("Transplant Experiment - American beech") +
  theme_bw() +
  scale_x_discrete(labels=c("L" = "Northern Hardwoods", "M" = "Ecotone",
                            "H" = "Spruce-fir")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        plot.title = element_text(size=18)) +
  stat_summary(geom = 'text', label = symbol4, fun = max, vjust = -8.5, position = position_dodge(width = 0.75), size = 6)
browse_trans1

path_trans1 = ggplot(trans_all1_fagr1, aes(x = order, y = root_path1)) +
  geom_errorbar(aes(ymin=root_path1-sep, ymax=root_path1+sep), width = 0.075) +
  geom_bar(stat = "identity", position = "dodge") +
  ylab("Root Pathogen Infection (% of pots)") +
  xlab("") +
  ylim(0,20) +
  ggtitle("Transplant Experiment - American beech") +
  theme_bw() +
  scale_x_discrete(labels=c("L" = "Northern Hardwoods", "M" = "Ecotone",
                            "H" = "Spruce-fir")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        plot.title = element_text(size=18)) +
  stat_summary(geom = 'text', label = symbol4, fun = max, vjust = -9, position = position_dodge(width = 0.75), size = 6)
path_trans1

herb_trans1 = ggplot(trans_all1_fagr1, aes(x = order, y = herb1)) +
  geom_errorbar(aes(ymin=herb1-sep, ymax=herb1+sep), width = 0.075) +
  geom_bar(stat = "identity", position = "dodge") +
  ylab("Leaf Herbivory (% of pots)") +
  xlab("") +
  ylim(0,30) +
  ggtitle("Transplant Experiment - American beech") +
  theme_bw() +
  scale_x_discrete(labels=c("L" = "Northern Hardwoods", "M" = "Ecotone",
                            "H" = "Spruce-fir")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        plot.title = element_text(size=18)) +
  stat_summary(geom = 'text', label = symbol4, fun = max, vjust = -6, position = position_dodge(width = 0.75), size = 6)
herb_trans1

tiff(file = "trans_bio.tif", width = 18, height = 12, units = 'in', res = 600, pointsize = 11)

trans_plots_bio = plot_grid(browse_trans, path_trans, herb_trans, 
                            browse_trans1, path_trans1, herb_trans1, nrow = 2, ncol = 3,
                            labels = c("A", "B", "C", "D", "E", "F"))
trans_plots_bio

dev.off()

###biomass stuff

#trans experiment

trans_alla = subset(trans_all, spp == "ACSA")
trans_allf = subset(trans_all, spp == "FAGR")

trans_alla$order = factor(trans_alla$site, levels=c("L", "M", "H"),
                          labels = c("Hardwoods", "Ecotone", "Spruce-fir"))
trans_allf$order = factor(trans_allf$site, levels=c("L", "M", "H"),
                          labels = c("Hardwoods", "Ecotone", "Spruce-fir"))
trans_alla$order1 = factor(trans_alla$source, levels=c("L", "M", "H"),
                           labels = c("Hardwood soils", "Ecotone soils", "Spruce-fir soils"))
trans_allf$order1 = factor(trans_allf$source, levels=c("L", "M", "H"),
                           labels = c("Hardwood soils", "Ecotone soils", "Spruce-fir soils"))

write.csv(trans_alla, "trans_alla.csv")
trans_alla = read.csv("trans_alla.csv", header = TRUE)
write.csv(trans_allf, "trans_allf.csv")
trans_allf = read.csv("trans_allf.csv", header = TRUE)

leaf_plot1 = ggplot(trans_alla, aes(x = site, y = num_leaves_avg, fill = treatment)) +
  geom_boxplot()
leaf_plot1

h_plot1 = ggplot(trans_alla, aes(x = site, y = as.numeric(height_avg), fill = treatment)) +
  geom_boxplot()
h_plot1

ratio_plot1 = ggplot(trans_alla, aes(x = site, y = root_shoot_ratio_avg, fill = treatment)) +
  geom_boxplot()
ratio_plot1

lmm3 = lmer(root_shoot_ratio_avg ~ site * source + (1|mountain), data = trans_alla)
summary(lmm3)
anova(lmm3)
Anova(lmm3)

trans_alla$BV <- interaction(trans_alla$site, trans_alla$source)

m1  <- glmer(root_shoot_ratio_avg ~ -1 + BV + (1|mountain), family = gaussian, data = trans_alla)

summary(glht(m1, linfct = mcp(BV = "Tukey")))

lmm4 = lmer(root_shoot_ratio_avg ~ site * source + (1|mountain), data = trans_allf)
summary(lmm4)
anova(lmm4)
Anova(lmm4)

trans_allf$BV <- interaction(trans_allf$site, trans_allf$source)

m1  <- glmer(root_shoot_ratio_avg ~ -1 + BV + (1|mountain), family = gaussian, data = trans_allf)

summary(glht(m1, linfct = mcp(BV = "Tukey")))

#acsa root/shoot
symbol14 = c("a", "a", "a", "a", "a", "a", "b", "b", "b")
rs_plota1 = ggplot(trans_alla, aes(x = order, y = root_shoot_ratio_avg, fill = order1)) +
  geom_boxplot() +
  ylab("Root: Shoot") +
  xlab("") +
  ggtitle("Sugar Maple") +
  scale_fill_brewer(palette = "Set1", direction = 1, name = "Soil Source") +
  facet_grid(~ order, scales = "free_x") +
  ylim(0.5,3.5) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        strip.text = element_text(size=12)) +
  theme(axis.text.x = element_blank(),
        axis.text=element_text(size=10),
        axis.ticks.x=element_blank(),
        axis.title=element_text(size=10)) +
  stat_summary(geom = 'text', label = symbol14, fun.y = max, vjust = -1, position = position_dodge(width = 0.75))
rs_plota1

#fagr shoot/root

symbol15 = c("a", "a", "a", "a", "a", "a", "a", "a", "a")
rs_plotf1 = ggplot(trans_allf, aes(x = order, y = root_shoot_ratio_avg, fill = order1)) +
  geom_boxplot() +
  ylab("Root : Shoot") +
  xlab("") +
  ggtitle("American Beech") +
  scale_fill_brewer(palette = "Set1", direction = 1, name = "Soil Source") +
  facet_grid(~ order, scales = "free_x") +
  ylim(0,2.5) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        strip.text = element_text(size=12)) +
  theme(axis.text.x = element_blank(),
        axis.text=element_text(size=10),
        axis.ticks.x=element_blank(),
        axis.title=element_text(size=10)) +
  stat_summary(geom = 'text', label = symbol15, fun.y = max, vjust = -1, position = position_dodge(width = 0.75))
rs_plotf1


#acsa total biomass

lmm5 = lmer(total_bio_avg ~ site * source + (1|mountain), data = trans_alla)
summary(lmm5)
anova(lmm5)
Anova(lmm5)

#trans_allf$BV <- interaction(trans_allf$site, trans_allf$source)

m1  <- glmer(total_bio_avg ~ -1 + BV + (1|mountain), family = gaussian, data = trans_alla)

summary(glht(m1, linfct = mcp(BV = "Tukey")))

lmm6 = lmer(total_bio_avg ~ site * source + (1|mountain), data = trans_allf)
summary(lmm6)
anova(lmm6)
Anova(lmm6)

#trans_allf$BV <- interaction(trans_allf$site, trans_allf$source)

m1  <- glmer(total_bio_avg ~ -1 + BV + (1|mountain), family = gaussian, data = trans_allf)

summary(glht(m1, linfct = mcp(BV = "Tukey")))

symbol16 = c("a", "a", "a", "a", "a", "a", "a", "a", "a")
tb_plota1 = ggplot(trans_alla, aes(x = order, y = total_bio_avg, fill = order1)) +
  geom_boxplot() +
  ylab("Total Dry Biomass (g)") +
  xlab("") +
  ggtitle("Sugar Maple") +
  scale_fill_brewer(palette = "Set1", direction = 1, name = "Soil Source") +
  facet_grid(~ order, scales = "free_x") +
  ylim(0,1) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        strip.text = element_text(size=12)) +
  theme(axis.text.x = element_blank(),
        axis.text=element_text(size=10),
        axis.ticks.x=element_blank(),
        axis.title=element_text(size=10)) +
  stat_summary(geom = 'text', label = symbol16, fun.y = max, vjust = -1, position = position_dodge(width = 0.75))
tb_plota1

#fagr total biomass
symbol16 = c("a", "a", "a", "a", "a", "a", "a", "a", "a")
tb_plotf1 = ggplot(trans_allf, aes(x = order, y = total_bio_avg, fill = order1)) +
  geom_boxplot() +
  ylab("Total Dry Biomass (g)") +
  xlab("") +
  ggtitle("American Beech") +
  scale_fill_brewer(palette = "Set1", direction = 1, name = "Soil Source") +
  facet_grid(~ order, scales = "free_x") +
  ylim(0,2.5) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        strip.text = element_text(size=12)) +
  theme(axis.text.x = element_blank(),
        axis.text=element_text(size=10),
        axis.ticks.x=element_blank(),
        axis.title=element_text(size=10)) +
  stat_summary(geom = 'text', label = symbol16, fun.y = max, vjust = -1, position = position_dodge(width = 0.75))
tb_plotf1

tbio_plot1 = ggplot(trans_alla, aes(x = site, y = total_bio_avg, fill = treatment)) +
  geom_boxplot()
tbio_plot1

abv_plot1 = ggplot(trans_alla, aes(x = site, y = total_abv_avg, fill = treatment)) +
  geom_boxplot()
abv_plot1

blw_plot1 = ggplot(trans_alla, aes(x = site, y = total_below_avg, fill = treatment)) +
  geom_boxplot()
blw_plot1

tiff(file = "trans_grow.tif", width = 20, height = 12, units = 'in', res = 600, pointsize = 11)

trans_grow = plot_grid(rgr_plota1, rs_plota1, tb_plota1,
                       rgr_plotf1, rs_plotf1, tb_plotf1,
                       ncol = 3, nrow = 2)
trans_grow

dev.off()

#net experiment

net_alla = subset(net_all, spp == "ACSA")
net_allf = subset(net_all, spp == "FAGR")

net_alla$order = factor(net_alla$site, levels=c("L", "M", "H"),
                        labels = c("Hardwoods", "Ecotone", "Spruce-fir"))
net_allf$order = factor(net_allf$site, levels=c("L", "M", "H"),
                        labels = c("Hardwoods", "Ecotone", "Spruce-fir"))
net_alla$order1 = factor(net_alla$treatment, levels=c("N", "D", "C"),
                         labels = c("Network", "Disturbed", "Control"))
net_allf$order1 = factor(net_allf$treatment, levels=c("N", "D", "C"),
                         labels = c("Network", "Disturbed", "Control"))

write.csv(net_alla, "net_alla.csv")
net_alla = read.csv("net_alla.csv", header = TRUE)
write.csv(net_allf, "net_allf.csv")
net_allf = read.csv("net_allf.csv", header = TRUE)

lmm7 = lmer(rhgr ~ site * treatment + (1|mountain), data = net_alla)
summary(lmm7)
anova(lmm7)
Anova(lmm7)

net_alla$BV <- interaction(net_alla$site, net_alla$treatment)

m1  <- glmer(rhgr ~ -1 + BV + (1|mountain), family = gaussian, data = net_alla)

summary(glht(m1, linfct = mcp(BV = "Tukey")))

lmm8 = lmer(rhgr ~ site * treatment + (1|mountain), data = net_allf)
summary(lmm8)
anova(lmm8)
Anova(lmm8)

net_allf$BV <- interaction(net_allf$site, net_allf$treatment)

m1  <- glmer(rhgr ~ -1 + BV + (1|mountain), family = gaussian, data = net_allf)

summary(glht(m1, linfct = mcp(BV = "Tukey")))

lmm17 = lmer(root_shoot_ratio_avg ~ site * treatment + (1|mountain), data = net_alla)
summary(lmm17)
anova(lmm17)
Anova(lmm17)

net_alla$BV <- interaction(net_alla$site, net_alla$treatment)

m1  <- glmer(root_shoot_ratio_avg ~ -1 + BV + (1|mountain), family = gaussian, data = net_alla)

summary(glht(m1, linfct = mcp(BV = "Tukey")))

lmm18 = lmer(root_shoot_ratio_avg ~ site * treatment + (1|mountain), data = net_allf)
summary(lmm18)
anova(lmm18)
Anova(lmm18)

net_allf$BV <- interaction(net_allf$site, net_allf$treatment)

m1  <- glmer(root_shoot_ratio_avg ~ -1 + BV + (1|mountain), family = gaussian, data = net_allf)

summary(glht(m1, linfct = mcp(BV = "Tukey")))


#acsa root/shoot
symbol17 = c("a", "a", "a", "a", "a", "a", "b", "b", "b")
net_rs_plota1 = ggplot(net_alla, aes(x = order, y = root_shoot_ratio_avg, fill = order1)) +
  geom_boxplot() +
  ylab("Root : Shoot") +
  xlab("") +
  ggtitle("Sugar Maple") +
  scale_fill_brewer(palette = "Set1", direction = 1, name = "Treatment") +
  facet_grid(~ order, scales = "free_x") +
  ylim(0,3.5) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        strip.text = element_text(size=12)) +
  theme(axis.text.x = element_blank(),
        axis.text=element_text(size=10),
        axis.ticks.x=element_blank(),
        axis.title=element_text(size=10)) +
  stat_summary(geom = 'text', label = symbol17, fun.y = max, vjust = -1, position = position_dodge(width = 0.75))
net_rs_plota1

#fagr root/shoot
symbol18 = c("a", "a", "a", "a", "a", "a", "a", "a", "a")
net_rs_plotf1 = ggplot(net_allf, aes(x = order, y = root_shoot_ratio_avg, fill = order1)) +
  geom_boxplot() +
  ylab("Root : Shoot") +
  xlab("") +
  ggtitle("American Beech") +
  scale_fill_brewer(palette = "Set1", direction = 1, name = "Treatment") +
  facet_grid(~ order, scales = "free_x") +
  ylim(0,3.5) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        strip.text = element_text(size=12)) +
  theme(axis.text.x = element_blank(),
        axis.text=element_text(size=10),
        axis.ticks.x=element_blank(),
        axis.title=element_text(size=10)) +
  stat_summary(geom = 'text', label = symbol18, fun.y = max, vjust = -1, position = position_dodge(width = 0.75))
net_rs_plotf1

lmm9 = lmer(diam_change ~ site * treatment + (1|mountain), data = net_alla)
summary(lmm9)
anova(lmm9)
Anova(lmm9)

#net_allf$BV <- interaction(net_allf$site, net_allf$treatment)

m1  <- glmer(diam_change ~ -1 + BV + (1|mountain), family = gaussian, data = net_alla)

summary(glht(m1, linfct = mcp(BV = "Tukey")))

lmm10 = lmer(diam_change ~ site * treatment + (1|mountain), data = net_allf)
summary(lmm10)
anova(lmm10)
Anova(lmm10)

#net_allf$BV <- interaction(net_allf$site, net_allf$treatment)

m1  <- glmer(diam_change ~ -1 + BV + (1|mountain), family = gaussian, data = net_allf)

summary(glht(m1, linfct = mcp(BV = "Tukey")))

#acsa daim
symbol19 = c("a", "a", "a", "a", "a", "a", "a", "a", "a")
net_d_plota1 = ggplot(net_alla, aes(x = order, y = diam_change, fill = order1)) +
  geom_boxplot() +
  ylab("") +
  xlab("") +
  ggtitle("Sugar Maple") +
  scale_fill_manual(values = c("#fbe45b", "#2c7c94", "#a65852")) +
  facet_grid(~ order, scales = "free_x") +
  ylim(0,3) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        strip.text = element_text(size=20)) +
  theme(axis.text.x = element_blank(),
        axis.text=element_text(size=20),
        plot.title = element_text(size = 24),
        legend.title = element_text(size=20), 
        legend.text = element_text(size=20),
        axis.ticks.x=element_blank(),
        legend.position = "none",
        axis.title=element_text(size=20)) +
  stat_summary(geom = 'text', label = symbol19, fun.y = max, vjust = -1, position = position_dodge(width = 0.75), size = 9)
net_d_plota1

#fagr diam
symbol19 = c("a", "a", "a", "a", "a", "a", "a", "a", "a")
net_d_plotf1 = ggplot(net_allf, aes(x = order, y = diam_change, fill = order1)) +
  geom_boxplot() +
  ylab("") +
  xlab("") +
  ggtitle("American Beech") +
  scale_fill_manual(values = c("#fbe45b", "#2c7c94", "#a65852")) +
  facet_grid(~ order, scales = "free_x") +
  ylim(0,3) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        strip.text = element_text(size=20)) +
  theme(axis.text.x = element_blank(),
        axis.text=element_text(size=20),
        axis.ticks.x=element_blank(),
        plot.title = element_text(size = 24),
        legend.title = element_text(size=20), 
        legend.text = element_text(size=20),
        legend.position = "none",
        axis.title=element_text(size=16)) +
  stat_summary(geom = 'text', label = symbol19, fun.y = max, vjust = -1, position = position_dodge(width = 0.75), size = 9)
net_d_plotf1

lmm11 = lmer(total_bio_avg ~ site * treatment + (1|mountain), data = net_alla)
summary(lmm11)
anova(lmm11)
Anova(lmm11)

#net_allf$BV <- interaction(net_allf$site, net_allf$treatment)

m1  <- glmer(total_bio_avg ~ -1 + BV + (1|mountain), family = gaussian, data = net_alla)

summary(glht(m1, linfct = mcp(BV = "Tukey")))


lmm12 = lmer(total_bio_avg ~ site * treatment + (1|mountain), data = net_allf)
summary(lmm12)
anova(lmm12)
Anova(lmm12)

#net_allf$BV <- interaction(net_allf$site, net_allf$treatment)

m1  <- glmer(total_bio_avg ~ -1 + BV + (1|mountain), family = gaussian, data = net_allf)

summary(glht(m1, linfct = mcp(BV = "Tukey")))

#acsa total biomass
symbol20 = c("a", "a", "a", "a", "a", "a", "a", "a", "a")
net_b_plota1 = ggplot(net_alla, aes(x = order, y = total_bio_avg, fill = order1)) +
  geom_boxplot() +
  ylab("Total Dry Biomass (g)") +
  xlab("") +
  ggtitle("Sugar Maple") +
  scale_fill_brewer(palette = "Set1", direction = 1, name = "Treatment") +
  facet_grid(~ order, scales = "free_x") +
  ylim(0,12) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        strip.text = element_text(size=12)) +
  theme(axis.text.x = element_blank(),
        axis.text=element_text(size=10),
        axis.ticks.x=element_blank(),
        axis.title=element_text(size=10)) +
  stat_summary(geom = 'text', label = symbol20, fun.y = max, vjust = -1, position = position_dodge(width = 0.75))
net_b_plota1

#fagr total biomass
symbol21 = c("a", "a", "a", "a", "a", "a", "a", "a", "a")
net_b_plotf1 = ggplot(net_allf, aes(x = order, y = total_bio_avg, fill = order1)) +
  geom_boxplot() +
  ylab("Total Dry Biomass (g)") +
  xlab("") +
  ggtitle("American Beech") +
  scale_fill_brewer(palette = "Set1", direction = 1, name = "Treatment") +
  facet_grid(~ order, scales = "free_x") +
  ylim(0,10) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        strip.text = element_text(size=12)) +
  theme(axis.text.x = element_blank(),
        axis.text=element_text(size=10),
        axis.ticks.x=element_blank(),
        axis.title=element_text(size=10)) +
  stat_summary(geom = 'text', label = symbol21, fun.y = max, vjust = -1, position = position_dodge(width = 0.75))
net_b_plotf1

lmm13 = lmer(num_leaves ~ site * treatment + (1|mountain), data = net_alla)
summary(lmm13)
anova(lmm13)
Anova(lmm13)

#net_allf$BV <- interaction(net_allf$site, net_allf$treatment)

m1  <- glmer(num_leaves ~ -1 + BV + (1|mountain), family = gaussian, data = net_alla)

summary(glht(m1, linfct = mcp(BV = "Tukey")))

lmm14 = lmer(num_leaves ~ site * treatment + (1|mountain), data = net_allf)
summary(lmm14)
anova(lmm14)
Anova(lmm14)

#net_allf$BV <- interaction(net_allf$site, net_allf$treatment)

m1  <- glmer(num_leaves ~ -1 + BV + (1|mountain), family = gaussian, data = net_allf)

summary(glht(m1, linfct = mcp(BV = "Tukey")))

#acsa leaves
symbol22 = c("a", "a", "a", "a", "a", "a", "a", "a", "a")
net_l_plota1 = ggplot(net_alla, aes(x = order, y = num_leaves, fill = order1)) +
  geom_boxplot() +
  ylab("") +
  xlab("") +
  ggtitle("Sugar Maple") +
  scale_fill_manual(values = c("#fbe45b", "#2c7c94", "#a65852")) +
  facet_grid(~ order, scales = "free_x") +
  ylim(0,30) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        strip.text = element_text(size=20)) +
  theme(axis.text.x = element_blank(),
        axis.text=element_text(size=20),
        legend.title = element_text(size=20), 
        legend.text = element_text(size=20),
        plot.title = element_text(size = 24),
        legend.position = "none",
        axis.ticks.x=element_blank(),
        axis.title=element_text(size=20)) +
  stat_summary(geom = 'text', label = symbol22, fun.y = max, vjust = -1, position = position_dodge(width = 0.75), size = 9)
net_l_plota1

#fagr leaves
symbol23 = c("a", "a", "a", "a", "a", "a", "a", "b", "b")
net_l_plotf1 = ggplot(net_allf, aes(x = order, y = num_leaves, fill = order1)) +
  geom_boxplot() +
  ylab("") +
  xlab("") +
  ggtitle("American Beech") +
  scale_fill_manual(values = c("#fbe45b", "#2c7c94", "#a65852"), name = "Treatment",
                    label = c("Networked", "Disrupted", "Control")) +
  facet_grid(~ order, scales = "free_x") +
  ylim(0,40) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        strip.text = element_text(size=20)) +
  theme(axis.text.x = element_blank(),
        axis.text=element_text(size=20),
        plot.title = element_text(size = 24),
        axis.ticks.x=element_blank(),
        legend.title = element_text(size=24), 
        legend.text = element_text(size=22),
        legend.position = "none",
        axis.title=element_text(size=20)) +
  stat_summary(geom = 'text', label = symbol23, fun.y = max, vjust = -1, position = position_dodge(width = 0.75), size = 9)
net_l_plotf1

ld_supp = plot_grid(net_d_plota1, net_d_plotf1, nrow = 2, ncol = 1, align = "v")
y.grob3 <- textGrob("Mean Basal Diameter Change (mm)", 
                    gp=gpar(fontface="bold", col="black", fontsize=24), rot=90)
ld_supp1 = grid.arrange(arrangeGrob(ld_supp, left = y.grob3))

ld_supp2 = plot_grid(net_l_plota1, net_l_plotf1, nrow = 2, ncol = 1, align = "v")
y.grob4 <- textGrob("Mean Leaf Number", 
                    gp=gpar(fontface="bold", col="black", fontsize=24), rot=90)
ld_supp3 = grid.arrange(arrangeGrob(ld_supp2, left = y.grob4))


legend <- get_legend(
  net_l_plotf1 +
    guides(color = guide_legend(nrow = 3)) +
    theme(legend.position = "bottom")
)

tiff(file = "net_grow_supp.tif", width = 14, height = 10, units = 'in', res = 600, pointsize = 11)

# Combine combined plot and legend using plot_grid()
ld_supp5 = plot_grid(ld_supp1, ld_supp3, ncol = 2, scale = 0.9)
#legend, ncol=2, nrow = 2,rel_heights = c(3, 0.4))
ld_supp5
ld_supp6 = plot_grid(ld_supp5, legend, ncol=1, nrow = 2,rel_heights = c(3, 0.4))
ld_supp6

dev.off()


tiff(file = "net_grow.tif", width = 30, height = 16, units = 'in', res = 600, pointsize = 11)

net_grow = plot_grid(net_rgr_plota1,  net_rs_plota1, net_b_plota1,
                     net_rgr_plotf1, net_rs_plotf1, net_b_plotf1,
                     nrow = 2, ncol = 3)
net_grow

dev.off()

net_grow_supp = plot_grid(net_d_plota1, net_l_plota1,
                          net_d_plotf1, net_l_plotf1, nrow = 2, ncol = 2)
net_grow_supp

tiff(file = "net_grow_supp.tif", width = 12, height = 10, units = 'in', res = 600, pointsize = 11)

net_grow_supp

dev.off()

#other visuals
h_plot = ggplot(net_alla, aes(x = site, y = height_avg, fill = treatment)) +
  geom_boxplot()
h_plot

ratio_plot = ggplot(net_alla, aes(x = site, y = root_shoot_ratio_avg, fill = treatment)) +
  geom_boxplot()
ratio_plot

tbio_plot = ggplot(net_alla, aes(x = site, y = total_bio_avg, fill = treatment)) +
  geom_boxplot()
tbio_plot

abv_plot = ggplot(net_alla, aes(x = site, y = total_abv_avg, fill = treatment)) +
  geom_boxplot()
abv_plot

blw_plot = ggplot(net_alla, aes(x = site, y = total_below_avg, fill = treatment)) +
  geom_boxplot()
blw_plot

###acsa leaf

#install.packages("afex")
install.packages("nlme")
library(afex)
library(nlme)

N_null = lme(leaf_N ~ site + treat, random = 1~mountain, data = tissues_a)
summary(N_null)

N_null = lmer(leaf_N ~ site * treat + (1|mountain), data = tissues_a)
summary(N_null)
anova(N_null)
Anova(N_null)

#tissues_a$BV <- interaction(tissues_a$site, tissues_a$treat)

m1  <- glmer(leaf_N ~ -1 + BV + (1|mountain), family = gaussian, data = tissues_a)

summary(glht(m1, linfct = mcp(BV = "Tukey")))

P_null = lmer(leaf_P ~ site * treat + (1|mountain), data = tissues_a)
summary(P_null)
anova(P_null)
Anova(P_null)
#tissues_a$BV <- interaction(tissues_a$site, tissues_a$treat)

m1  <- glmer(leaf_P ~ -1 + BV + (1|mountain), family = gaussian, data = tissues_a)

summary(glht(m1, linfct = mcp(BV = "Tukey")))

K_null = lmer(leaf_K ~ site * treat + (1|mountain), data = tissues_a)
summary(K_null)
anova(K_null)
Anova(K_null)
m1  <- glmer(leaf_K ~ -1 + BV + (1|mountain), family = gaussian, data = tissues_a)

summary(glht(m1, linfct = mcp(BV = "Tukey")))

CA_null = lmer(leaf_Ca ~ site * treat + (1|mountain), data = tissues_a)
summary(CA_null)
anova(CA_null)
Anova(CA_null)
m1  <- glmer(leaf_Ca ~ -1 + BV + (1|mountain), family = gaussian, data = tissues_a)

summary(glht(m1, linfct = mcp(BV = "Tukey")))

MG_null = lmer(leaf_Mg ~ site * treat + (1|mountain), data = tissues_a)
summary(MG_null)
anova(MG_null)
Anova(MG_null)
m1  <- glmer(leaf_Mg ~ -1 + BV + (1|mountain), family = gaussian, data = tissues_a)

summary(glht(m1, linfct = mcp(BV = "Tukey")))

AL_null = lmer(leaf_Al ~ site * treat + (1|mountain), data = tissues_a)
summary(AL_null)
anova(AL_null)
Anova(AL_null)
m1  <- glmer(leaf_Al ~ -1 + BV + (1|mountain), family = gaussian, data = tissues_a)

summary(glht(m1, linfct = mcp(BV = "Tukey")))

###fagr leaf

N_null = lmer(leaf_N ~ site * treat + (1|mountain), data = tissues_f)
summary(N_null)
anova(N_null)
Anova(N_null)

tissues_f$BV <- interaction(tissues_f$site, tissues_f$treat)

m1  <- glmer(leaf_N ~ -1 + BV + (1|mountain), family = gaussian, data = tissues_f)

summary(glht(m1, linfct = mcp(BV = "Tukey")))

P_null = lmer(leaf_P ~ site * treat + (1|mountain), data = tissues_f)
summary(P_null)
anova(P_null)
Anova(P_null)
m1  <- glmer(leaf_P ~ -1 + BV + (1|mountain), family = gaussian, data = tissues_f)

summary(glht(m1, linfct = mcp(BV = "Tukey")))

K_null = lmer(leaf_K ~ site * treat + (1|mountain), data = tissues_f)
summary(K_null)
anova(K_null)
Anova(K_null)
m1  <- glmer(leaf_K ~ -1 + BV + (1|mountain), family = gaussian, data = tissues_f)

summary(glht(m1, linfct = mcp(BV = "Tukey")))

CA_null = lmer(leaf_Ca ~ site * treat + (1|mountain), data = tissues_f)
summary(CA_null)
anova(CA_null)
Anova(CA_null)
m1  <- glmer(leaf_Ca ~ -1 + BV + (1|mountain), family = gaussian, data = tissues_f)

summary(glht(m1, linfct = mcp(BV = "Tukey")))

MG_null = lmer(leaf_Mg ~ site * treat + (1|mountain), data = tissues_f)
summary(MG_null)
anova(MG_null)
Anova(MG_null)
m1  <- glmer(leaf_Mg ~ -1 + BV + (1|mountain), family = gaussian, data = tissues_f)

summary(glht(m1, linfct = mcp(BV = "Tukey")))

AL_null = lmer(leaf_Al ~ site * treat + (1|mountain), data = tissues_f)
summary(AL_null)
anova(AL_null)
Anova(AL_null)
m1  <- glmer(leaf_Al ~ -1 + BV + (1|mountain), family = gaussian, data = tissues_f)

summary(glht(m1, linfct = mcp(BV = "Tukey")))

#acsa site plots
tn_plot = ggplot(tissue_a, aes(x = order, y = N)) +
  geom_errorbar(aes(ymin=N-SEN, ymax=N+SEN), width = 0.075) +
  geom_point(stat = "identity", size = 6, pch = 18) +
  ylab("Foliar N Concentration (%)") +
  xlab("") +
  theme_bw() +
  ylim(1.88, 2.08) +
  scale_x_discrete(labels=c("L" = "Hardwoods", "M" = "Ecotone",
                            "H" = "Spruce-fir")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        plot.title = element_text(size=18))
tn_plot

tp_plot = ggplot(tissue_a, aes(x = order, y = P)) +
  geom_errorbar(aes(ymin=P-SEP, ymax=P+SEP), width = 0.075) +
  geom_point(stat = "identity", size = 6, pch = 18) +
  ylab("Foliar P Concentration (%)") +
  xlab("") +
  theme_bw() +
  ylim(0.14,0.22) +
  scale_x_discrete(labels=c("L" = "Hardwoods", "M" = "Ecotone",
                            "H" = "Spruce-fir")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        plot.title = element_text(size=18))
tp_plot

tk_plot = ggplot(tissue_a, aes(x = order, y = K)) +
  geom_errorbar(aes(ymin=K-SEK, ymax=K+SEK), width = 0.075) +
  geom_point(stat = "identity", size = 6, pch = 18) +
  ylab("Foliar K Concentration (%)") +
  xlab("") +
  theme_bw() +
  ylim(0.5,1) +
  scale_x_discrete(labels=c("L" = "Hardwoods", "M" = "Ecotone",
                            "H" = "Spruce-fir")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        plot.title = element_text(size=18))
tk_plot

tc_plot = ggplot(tissue_a, aes(x = order, y = CA)) +
  geom_errorbar(aes(ymin=CA-SEC, ymax=CA+SEC), width = 0.075) +
  geom_point(stat = "identity", size = 6, pch = 18) +
  ylab("Foliar Ca Concentration (%)") +
  xlab("") +
  theme_bw() +
  ylim(0.5,1) +
  scale_x_discrete(labels=c("L" = "Hardwoods", "M" = "Ecotone",
                            "H" = "Spruce-fir")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        plot.title = element_text(size=18))
tc_plot

tm_plot = ggplot(tissue_a, aes(x = order, y = MG)) +
  geom_errorbar(aes(ymin=MG-SEM, ymax=MG+SEM), width = 0.075) +
  geom_point(stat = "identity", size = 6, pch = 18) +
  ylab("Foliar Mg Concentration (%)") +
  xlab("") +
  theme_bw() +
  scale_x_discrete(labels=c("L" = "Hardwoods", "M" = "Ecotone",
                            "H" = "Spruce-fir")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        plot.title = element_text(size=18))
tm_plot

ta_plot = ggplot(tissue_a, aes(x = order, y = AL)) +
  geom_errorbar(aes(ymin=AL-SEA, ymax=AL+SEA), width = 0.075) +
  geom_point(stat = "identity", size = 6, pch = 18) +
  ylab("Foliar Al Concentration (mg/kg)") +
  xlab("") +
  theme_bw() +
  ylim(60,140) +
  scale_x_discrete(labels=c("L" = "Hardwoods", "M" = "Ecotone",
                            "H" = "Spruce-fir")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        plot.title = element_text(size=18))
ta_plot

#acsa treatment plots
tn_plot1 = ggplot(tissue1_a, aes(x = order, y = N)) +
  geom_errorbar(aes(ymin=N-SEN, ymax=N+SEN), width = 0.075) +
  geom_point(stat = "identity", size = 6, pch = 18) +
  ylab("Foliar N Concentration (%)") +
  xlab("") +
  theme_bw() +
  scale_x_discrete(labels=c("N" = "Network", "D" = "Disrupted",
                            "C" = "Pot Control")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        plot.title = element_text(size=18))
tn_plot1

tp_plot1 = ggplot(tissue1_a, aes(x = order, y = P)) +
  geom_errorbar(aes(ymin=P-SEP, ymax=P+SEP), width = 0.075) +
  geom_point(stat = "identity", size = 6, pch = 18) +
  ylab("Foliar P Concentration (%)") +
  xlab("") +
  theme_bw() +
  ylim(0.14,0.22) +
  scale_x_discrete(labels=c("N" = "Network", "D" = "Disrupted",
                            "C" = "Pot Control")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        plot.title = element_text(size=18))
tp_plot1

tk_plot1 = ggplot(tissue1_a, aes(x = order, y = K)) +
  geom_errorbar(aes(ymin=K-SEK, ymax=K+SEK), width = 0.075) +
  geom_point(stat = "identity", size = 6, pch = 18) +
  ylab("Foliar K Concentration (%)") +
  xlab("") +
  theme_bw() +
  ylim(0.5,1) +
  scale_x_discrete(labels=c("N" = "Network", "D" = "Disrupted",
                            "C" = "Pot Control")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        plot.title = element_text(size=18))
tk_plot1

tc_plot1 = ggplot(tissue1_a, aes(x = order, y = CA)) +
  geom_errorbar(aes(ymin=CA-SEC, ymax=CA+SEC), width = 0.075) +
  geom_point(stat = "identity", size = 6, pch = 18) +
  ylab("Foliar Ca Concentration (%)") +
  xlab("") +
  theme_bw() +
  ylim(0.5,1) +
  scale_x_discrete(labels=c("N" = "Network", "D" = "Disrupted",
                            "C" = "Pot Control")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        plot.title = element_text(size=18))
tc_plot1

tm_plot1 = ggplot(tissue1_a, aes(x = order, y = MG)) +
  geom_errorbar(aes(ymin=MG-SEM, ymax=MG+SEM), width = 0.075) +
  geom_point(stat = "identity", size = 6, pch = 18) +
  ylab("Foliar Mg Concentration (%)") +
  xlab("") +
  theme_bw() +
  ylim(0.15, 0.25) +
  scale_x_discrete(labels=c("N" = "Network", "D" = "Disrupted",
                            "C" = "Pot Control")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        plot.title = element_text(size=18))
tm_plot1

tiss = rbind(tissue1_a, tissue1_f)
tiss2 = rbind(issue_a, issue_f)
shading = data.frame(col)
shading$col = "grey53"
symbol_al = c("b", "c", "a", "c", "b", "c")

ta_plot1 = ggplot(tiss, aes(x = order, y = AL, fill = spp)) +
  geom_errorbar(aes(ymin=AL-SEA, ymax=AL+SEA), width = 0.2, lwd = 1, position = position_dodge(width = 0.4)) +
  geom_point(stat = "identity", size = 6, pch = 23, position = position_dodge(width = 0.4)) +
  ylab("Foliar Al Concentration (mg/kg)") +
  xlab("CMN Treatment") +
  theme_bw() +
  ylim(30,160) +
  geom_rect(xmin = 1.5, xmax = 2.5, ymin = -Inf, ymax = Inf, alpha = 0.1) +
  scale_fill_manual(name = "Species", labels = c("Sugar maple", "American beech"), values = c("#2c7c94", "#fbe45b")) +
  scale_x_discrete(labels=c("N" = "Network", "D" = "Disrupted",
                            "C" = "Pot Control")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        plot.title = element_text(size=18),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14)) +
  stat_summary(geom = 'text', label = symbol_al, fun = max, vjust = -2.1, position = position_dodge(width = 0.4), size = 9)

ta_plot1

tiff(file = "AL_grow1.tif", width = 8, height = 5, units = 'in', res = 600, pointsize = 11)

ta_plot1

dev.off()

#updated Al stuff

treat1 <- c("Network", "Disrupted", "Control")
names(treat1) = c("N", "D", "C")
write.csv(tiss2, "tiss2.csv")
tiss2 = read.csv("tiss2.csv", header = TRUE)

tiss2$order = factor(issue_a$treat, levels=c("N", "D", "C"), labels = c("Network", "Disrupted", "Control"))
tiss2$order1 = factor(issue_a$site, levels=c("L", "M", "H"), labels = c("Hardwoods", "Ecotone", "Spruce-fir"))


symbol_al = c("c", "b", "c", "e", "b", "e", "b", "c", "b", "c", "a", "e", "c", "c", "d", "d", "b", "d")

ta_plot2 = ggplot(tiss2, aes(x = order1, y = AL, fill = spp)) +
  geom_errorbar(aes(ymin=AL-SEA, ymax=AL+SEA), width = 0.5, lwd = 1, position = position_dodge(width = 0.5)) +
  geom_point(stat = "identity", size = 5, pch = 23, position = position_dodge(width = 0.5)) +
  ylab("Foliar Al Concentration (mg/kg)") +
  xlab("Planting Site") +
  facet_grid(~ order, labeller = labeller(c("Network", "Disrupted", "Control"))) +
  theme_bw() +
  ylim(20,170) +
  geom_rect(xmin = 1.5, xmax = 2.5, ymin = -Inf, ymax = Inf, alpha = 0.1) +
  scale_fill_manual(name = "Species", labels = c("Sugar maple", "American beech"), values = c("#2c7c94", "#fbe45b")) +
  scale_x_discrete(labels=c("L" = "Harwoods", "M" = "Ecotone",
                            "H" = "Spruce-fir")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        plot.title = element_text(size=18),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        strip.text = element_text(size = 16)) +
  stat_summary(geom = 'text', label = symbol_al, fun = max, vjust = -1.7, position = position_dodge(width = 0.5), size = 6)

ta_plot2

tiff(file = "AL_grow2.tif", width = 8, height = 5, units = 'in', res = 600, pointsize = 11)

ta_plot2

dev.off()

acsa_leaf = plot_grid(tn_plot, tp_plot, tk_plot, tc_plot, tm_plot, ta_plot,
                      tn_plot1, tp_plot1, tk_plot1, tc_plot1, tm_plot1, ta_plot1,
                      nrow = 4, ncol = 3, labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L"),
                      vjust = 2)
title <- ggdraw() + 
  draw_label(
    "Sugar Maple Foliar Nutrients",
    fontface = 'bold',
    size = 20,
    x = 0,
    hjust = 0
  ) +
  theme(
    plot.margin = margin(0, 0, 0, 7)
  )

a_grid = plot_grid(
  title, acsa_leaf,
  ncol = 1,
  rel_heights = c(0.1, 2.5)
)

tiff(file = "acsa_leaf.tif", width = 14, height = 18, units = 'in', res = 600, pointsize = 11)

a_grid

dev.off()

#fagr site plots

tn_plot2 = ggplot(tissue_f, aes(x = order, y = N)) +
  geom_errorbar(aes(ymin=N-SEN, ymax=N+SEN), width = 0.075) +
  geom_point(stat = "identity", size = 6, pch = 18) +
  ylab("Foliar N Concentration (%)") +
  xlab("") +
  theme_bw() +
  scale_x_discrete(labels=c("L" = "Hardwoods", "M" = "Ecotone",
                            "H" = "Spruce-fir")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        plot.title = element_text(size=22))
tn_plot2

tp_plot2 = ggplot(tissue_f, aes(x = order, y = P)) +
  geom_errorbar(aes(ymin=P-SEP, ymax=P+SEP), width = 0.075) +
  geom_point(stat = "identity", size = 6, pch = 18) +
  ylab("Foliar P Concentration (%)") +
  xlab("") +
  theme_bw() +
  scale_x_discrete(labels=c("L" = "Hardwoods", "M" = "Ecotone",
                            "H" = "Spruce-fir")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        plot.title = element_text(size=18))
tp_plot2

tk_plot2 = ggplot(tissue_f, aes(x = order, y = K)) +
  geom_errorbar(aes(ymin=K-SEK, ymax=K+SEK), width = 0.075) +
  geom_point(stat = "identity", size = 6, pch = 18) +
  ylab("Foliar K Concentration (%)") +
  xlab("") +
  theme_bw() +
  ylim(0.7, 0.9) +
  scale_x_discrete(labels=c("L" = "Hardwoods", "M" = "Ecotone",
                            "H" = "Spruce-fir")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        plot.title = element_text(size=18))
tk_plot2

tc_plot2 = ggplot(tissue_f, aes(x = order, y = CA)) +
  geom_errorbar(aes(ymin=CA-SEC, ymax=CA+SEC), width = 0.075) +
  geom_point(stat = "identity", size = 6, pch = 18) +
  ylab("Foliar Ca Concentration (%)") +
  xlab("") +
  theme_bw() +
  ylim(0.5,1) +
  scale_x_discrete(labels=c("L" = "Hardwoods", "M" = "Ecotone",
                            "H" = "Spruce-fir")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        plot.title = element_text(size=18))
tc_plot2

tm_plot2 = ggplot(tissue_f, aes(x = order, y = MG)) +
  geom_errorbar(aes(ymin=MG-SEM, ymax=MG+SEM), width = 0.075) +
  geom_point(stat = "identity", size = 6, pch = 18) +
  ylab("Foliar Mg Concentration (%)") +
  xlab("") +
  theme_bw() +
  scale_x_discrete(labels=c("L" = "Hardwoods", "M" = "Ecotone",
                            "H" = "Spruce-fir")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        plot.title = element_text(size=18))
tm_plot2

ta_plot2 = ggplot(tissue_f, aes(x = order, y = AL)) +
  geom_errorbar(aes(ymin=AL-SEA, ymax=AL+SEA), width = 0.075) +
  geom_point(stat = "identity", size = 6, pch = 18) +
  ylab("Foliar Al Concentration (mg/kg)") +
  xlab("") +
  theme_bw() +
  ylim(0,100) +
  scale_x_discrete(labels=c("L" = "Hardwoods", "M" = "Ecotone",
                            "H" = "Spruce-fir")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        plot.title = element_text(size=18))
ta_plot2

#fagr treatment plots

tn_plot3 = ggplot(tissue1_f, aes(x = order, y = N)) +
  geom_errorbar(aes(ymin=N-SEN, ymax=N+SEN), width = 0.075) +
  geom_point(stat = "identity", size = 6, pch = 18) +
  ylab("Foliar N Concentration (%)") +
  xlab("") +
  theme_bw() +
  scale_x_discrete(labels=c("N" = "Network", "D" = "Disrupted",
                            "C" = "Pot Control")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        plot.title = element_text(size=18))
tn_plot3

tp_plot3 = ggplot(tissue1_f, aes(x = order, y = P)) +
  geom_errorbar(aes(ymin=P-SEP, ymax=P+SEP), width = 0.075) +
  geom_point(stat = "identity", size = 6, pch = 18) +
  ylab("Foliar P Concentration (%)") +
  xlab("") +
  theme_bw() +
  ylim(0.14,0.22) +
  scale_x_discrete(labels=c("N" = "Network", "D" = "Disrupted",
                            "C" = "Pot Control")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        plot.title = element_text(size=18))
tp_plot3

tk_plot3 = ggplot(tissue1_f, aes(x = order, y = K)) +
  geom_errorbar(aes(ymin=K-SEK, ymax=K+SEK), width = 0.075) +
  geom_point(stat = "identity", size = 6, pch = 18) +
  ylab("Foliar K Concentration (%)") +
  xlab("") +
  theme_bw() +
  ylim(0.5,1) +
  scale_x_discrete(labels=c("N" = "Network", "D" = "Disrupted",
                            "C" = "Pot Control")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        plot.title = element_text(size=18))
tk_plot3

tc_plot3 = ggplot(tissue1_f, aes(x = order, y = CA)) +
  geom_errorbar(aes(ymin=CA-SEC, ymax=CA+SEC), width = 0.075) +
  geom_point(stat = "identity", size = 6, pch = 18) +
  ylab("Foliar Ca Concentration (%)") +
  xlab("") +
  theme_bw() +
  ylim(0.5,1) +
  scale_x_discrete(labels=c("N" = "Network", "D" = "Disrupted",
                            "C" = "Pot Control")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        plot.title = element_text(size=18))
tc_plot3

tm_plot3 = ggplot(tissue1_f, aes(x = order, y = MG)) +
  geom_errorbar(aes(ymin=MG-SEM, ymax=MG+SEM), width = 0.075) +
  geom_point(stat = "identity", size = 6, pch = 18) +
  ylab("Foliar Mg Concentration (%)") +
  xlab("") +
  theme_bw() +
  scale_x_discrete(labels=c("N" = "Network", "D" = "Disrupted",
                            "C" = "Pot Control")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        plot.title = element_text(size=18))
tm_plot3

ta_plot3 = ggplot(tissue1_f, aes(x = order, y = AL)) +
  geom_errorbar(aes(ymin=AL-SEA, ymax=AL+SEA), width = 0.075) +
  geom_point(stat = "identity", size = 6, pch = 18) +
  ylab("Foliar Al Concentration (mg/kg)") +
  xlab("") +
  theme_bw() +
  scale_x_discrete(labels=c("N" = "Network", "D" = "Disrupted",
                            "C" = "Pot Control")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        plot.title = element_text(size=18))
ta_plot3

fagr_leaf = plot_grid(tn_plot2, tp_plot2, tk_plot2, tc_plot2, tm_plot2, ta_plot2,
                      tn_plot3, tp_plot3, tk_plot3, tc_plot3, tm_plot3, ta_plot3,
                      nrow = 4, ncol = 3, labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L"),
                      vjust = 2)
title <- ggdraw() + 
  draw_label(
    "Beech Foliar Nutrients",
    fontface = 'bold',
    size = 20,
    x = 0,
    hjust = 0
  ) +
  theme(
    plot.margin = margin(0, 0, 0, 7)
  )

f_grid = plot_grid(
  title, fagr_leaf,
  ncol = 1,
  rel_heights = c(0.1, 2.5)
)

tiff(file = "fagr_leaf.tif", width = 14, height = 18, units = 'in', res = 600, pointsize = 11)

f_grid

dev.off()

N_null = glmer(leaf_N ~ site + (1|tissues$mountain), data = tissues, family="gaussian", na.action=na.omit)
summary(N_null)

P_null = glmer(leaf_P ~ site + treat + (1|tissues$mountain), data = tissues, family="gaussian", na.action=na.omit)
summary(P_null)

K_null = glmer(leaf_K ~ site + (1|tissues$mountain), data = tissues, family="gaussian", na.action=na.omit)
summary(K_null)

CA_null = glmer(leaf_Ca ~ site + (1|tissue_a$mountain), data = tissue_a, family="gaussian", na.action=na.omit)
summary(CA_null)

MG_null = glmer(leaf_Mg ~ site + (1|tissues$mountain), data = tissues, family="gaussian", na.action=na.omit)
summary(MG_null)

AL_null = glmer(leaf_Al ~ site + (1|tissue_a$mountain), data = tissue_a, family="gaussian", na.action=na.omit)
summary(AL_null)

###soil
soil_avg = ddply(soils, ~site, summarise, ph = mean(pH), SE_ph = sd(pH)/sqrt((length(pH))), 
                 som = mean(OM_Pct), SE_som = sd(OM_Pct)/sqrt((length(OM_Pct))), 
                 cec = mean(ECEC), SE_cec = sd(ECEC)/sqrt((length(ECEC))),
                 #base = mean(Base_sum), SE_base = sd(Base_sum)/sqrt((length(Base_sum))), 
                 AL = mean(Al), SEA = sd(Al)/sqrt((length(Al))),
                 p = mean(Avail_P), SEp = sd(Avail_P)/sqrt((length(Avail_P))),
                 k = mean(K), SEk = sd(K)/sqrt((length(K))),
                 ca = mean(Ca), SEca = sd(Ca)/sqrt((length(Ca))),
                 mg = mean(Mg), SEmg = sd(Mg)/sqrt((length(Mg))))

soil_avg$order = factor(soil_avg$site, levels=c("L", "M", "H"))

soil_ph = ggplot(soil_avg, aes(x = order, y = ph)) +
  geom_errorbar(aes(ymin=ph-SE_ph, ymax=ph+SE_ph), width = 0.075) +
  geom_point(stat = "identity", size = 6, pch = 18) +
  ylab("pH") +
  xlab("") +
  theme_bw() +
  ylim(3.7,4.4) +
  scale_x_discrete(labels=c("L" = "Hardwoods", "M" = "Ecotone",
                            "H" = "Spruce-fir")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=20),
        axis.title=element_text(size=22),
        plot.title = element_text(size=20))
soil_ph

soil_som = ggplot(soil_avg, aes(x = order, y = som)) +
  geom_errorbar(aes(ymin=som-SE_som, ymax=som+SE_som), width = 0.075) +
  geom_point(stat = "identity", size = 6, pch = 18) +
  ylab("Soil Organic Matter (%)") +
  xlab("") +
  theme_bw() +
  ylim(0,50) +
  scale_x_discrete(labels=c("L" = "Hardwoods", "M" = "Ecotone",
                            "H" = "Spruce-fir")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=20),
        axis.title=element_text(size=22),
        plot.title = element_text(size=20))
soil_som

soil_cec = ggplot(soil_avg, aes(x = order, y = cec)) +
  geom_errorbar(aes(ymin=cec-SE_cec, ymax=cec+SE_cec), width = 0.075) +
  geom_point(stat = "identity", size = 6, pch = 18) +
  ylab("Cation Exchange Capacity") +
  xlab("") +
  theme_bw() +
  scale_x_discrete(labels=c("L" = "Hardwoods", "M" = "Ecotone",
                            "H" = "Spruce-fir")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=20),
        axis.title=element_text(size=22),
        plot.title = element_text(size=20))
soil_cec

soil_base = ggplot(soil_avg, aes(x = order, y = base)) +
  geom_errorbar(aes(ymin=base-SE_base, ymax=base+SE_base), width = 0.075) +
  geom_point(stat = "identity", size = 6, pch = 18) +
  ylab("Sum of Base Cations (100 meq/100 g)") +
  xlab("") +
  theme_bw() +
  scale_x_discrete(labels=c("L" = "Hardwoods", "M" = "Ecotone",
                            "H" = "Spruce-fir")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=20),
        axis.title=element_text(size=22),
        plot.title = element_text(size=20))
soil_base

soil_al = ggplot(soil_avg, aes(x = order, y = AL)) +
  geom_errorbar(aes(ymin=AL-SEA, ymax=AL+SEA), width = 0.075) +
  geom_point(stat = "identity", size = 6, pch = 18) +
  ylab("Al Concentration (mg/kg)") +
  xlab("") +
  theme_bw() +
  scale_x_discrete(labels=c("L" = "Hardwoods", "M" = "Ecotone",
                            "H" = "Spruce-fir")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=20),
        axis.title=element_text(size=22),
        plot.title = element_text(size=20))
soil_al

soil_p = ggplot(soil_avg, aes(x = order, y = p)) +
  geom_errorbar(aes(ymin=p-SEp, ymax=p+SEp), width = 0.075) +
  geom_point(stat = "identity", size = 6, pch = 18) +
  ylab("P Concentration (mg/kg)") +
  xlab("") +
  theme_bw() +
  scale_x_discrete(labels=c("L" = "Hardwoods", "M" = "Ecotone",
                            "H" = "Spruce-fir")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=20),
        axis.title=element_text(size=22),
        plot.title = element_text(size=20))
soil_p

soil_k = ggplot(soil_avg, aes(x = order, y = k)) +
  geom_errorbar(aes(ymin=k-SEk, ymax=k+SEk), width = 0.075) +
  geom_point(stat = "identity", size = 6, pch = 18) +
  ylab("K Concentration (mg/kg)") +
  xlab("") +
  theme_bw() +
  scale_x_discrete(labels=c("L" = "Hardwoods", "M" = "Ecotone",
                            "H" = "Spruce-fir")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=20),
        axis.title=element_text(size=22),
        plot.title = element_text(size=20))
soil_k

soil_ca = ggplot(soil_avg, aes(x = order, y = ca)) +
  geom_errorbar(aes(ymin=ca-SEca, ymax=ca+SEca), width = 0.075) +
  geom_point(stat = "identity", size = 6, pch = 18) +
  ylab("Ca Concentration (mg/kg)") +
  xlab("") +
  theme_bw() +
  scale_x_discrete(labels=c("L" = "Hardwoods", "M" = "Ecotone",
                            "H" = "Spruce-fir")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=20),
        axis.title=element_text(size=22),
        plot.title = element_text(size=20))
soil_ca

soil_mg = ggplot(soil_avg, aes(x = order, y = mg)) +
  geom_errorbar(aes(ymin=mg-SEmg, ymax=mg+SEmg), width = 0.075) +
  geom_point(stat = "identity", size = 6, pch = 18) +
  ylab("Mg Concentration (mg/kg)") +
  xlab("") +
  theme_bw() +
  scale_x_discrete(labels=c("L" = "Hardwoods", "M" = "Ecotone",
                            "H" = "Spruce-fir")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=20),
        axis.title=element_text(size=22),
        plot.title = element_text(size=20))
soil_mg

tiff(file = "vt_soils.tif", width = 20, height = 10, units = 'in', res = 600, pointsize = 11)

soil_plots = plot_grid(soil_ph, soil_som, soil_cec, soil_al, soil_p, soil_k, soil_ca, soil_mg, nrow = 2, ncol = 4, 
                       labels = c("A", "B", "C", "D", "E", "F", "G", "H"))
soil_plots
dev.off()

###greenhouse

green_exp = read.csv("greenhouse_full.csv", header = TRUE)
hist(green_exp$final_height)
hist(green_exp$end_biomass)

acsa_green = subset(green_exp, spp == "ACSA")
fagr_green = subset(green_exp, spp == "FAGR")

acsa_s = ddply(acsa_green, ~treat + source, summarise, surv = sum(survived))
fagr_s = ddply(fagr_green, ~treat + source, summarise, surv = sum(survived))

#models start

green_sur_a = glm(survived ~ treat * source, family = binomial, data = acsa_green)
summary(green_sur_a)
anova(green_sur_a)
Anova(green_sur_a, type = "II")

acsa_green$BV <- interaction(acsa_green$source, acsa_green$treat)

m1  <- glm(survived ~ -1 + BV, family = binomial, data = acsa_green)

summary(glht(m1, linfct = mcp(BV = "Tukey")))

green_sur_f = glm(survived ~ treat * source, family = binomial, data = fagr_green)
summary(green_sur_a)
anova(green_sur_f)
Anova(green_sur_f, type = "II")

fagr_green$BV <- interaction(fagr_green$source, fagr_green$treat)

m1  <- glm(survived ~ -1 + BV, family = binomial, data = fagr_green)

summary(glht(m1, linfct = mcp(BV = "Tukey")))

green_h_a = glm(final_height ~ treat * source, family = gaussian, data = acsa_green)
summary(green_h_a)
anova(green_h_a)
Anova(green_h_a, type = "II")

m1  <- glm(final_height ~ -1 + BV, family = gaussian, data = acsa_green)

summary(glht(m1, linfct = mcp(BV = "Tukey")))

green_h_f = glm(final_height ~ treat * source, family = gaussian, data = fagr_green)
summary(green_h_f)
anova(green_h_f)
Anova(green_h_f, type = "II")
m1  <- glm(final_height ~ -1 + BV, family = gaussian, data = fagr_green)

summary(glht(m1, linfct = mcp(BV = "Tukey")))

green_b_a = glm(end_biomass ~ treat * source, family = gaussian, data = acsa_green)
summary(green_b_a)
anova(green_b_a)
Anova(green_b_a, type = "II")
m1  <- glm(end_biomass ~ -1 + BV, family = gaussian, data = acsa_green)

summary(glht(m1, linfct = mcp(BV = "Tukey")))

green_b_f = glm(end_biomass ~ treat * source, family = gaussian, data = fagr_green)
summary(green_b_f)
anova(green_b_f)
Anova(green_b_f, type = "II")
m1  <- glm(end_biomass ~ -1 + BV, family = gaussian, data = fagr_green)

summary(glht(m1, linfct = mcp(BV = "Tukey")))

mod1 = lm(final_height ~ treatment + source, data = acsa_green)
Anova(mod1)

#final biomass and height

acsa_green$order = factor(acsa_green$source, levels=c("L", "M", "H"),
                          labels = c("Hardwood soils", "Ecotone soils", "Spruce-fir soils"))
acsa_green$order1 = factor(acsa_green$treat, levels=c("B", "FA", "FC", "PC"),
                           labels = c("Innoculated", "Sterilized Field", "Field", "Potting Mix"))
fagr_green$order = factor(fagr_green$source, levels=c("L", "M", "H"),
                          labels = c("Hardwood soils", "Ecotone soils", "Spruce-fir soils"))
fagr_green$order1 = factor(fagr_green$treat, levels=c("B", "FA", "FC", "PC"),
                           labels = c("Innoculated", "Sterilized Field", "Field", "Potting Mix"))
symbol1 = c("b", "a", "b", "b", "b", "a", "b", "b", "b", "b", "b", "ab")
ac_g_bio = ggplot(acsa_green, aes(x = order, y = end_biomass, fill = order1)) +
  geom_boxplot() +
  ylab("Final Dry Biomass (g)") +
  xlab("") +
  ggtitle("Sugar Maple") +
  scale_fill_brewer(palette = "Set1", direction = 1, name = "Treatment") +
  facet_grid(~ order, scales = "free_x") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        strip.text = element_text(size=12)) +
  theme(axis.text.x = element_blank(),
        axis.text=element_text(size=10),
        axis.ticks.x=element_blank(),
        axis.title=element_text(size=10)) +
  stat_summary(geom = 'text', label = symbol1, fun = max, vjust = -1, position = position_dodge(width = 0.75))
ac_g_bio

symbol3 = c("a", "b", "b", "ab", "a", "c", "b", "b", "a", "ab", "b", "a")
ac_g_h = ggplot(acsa_green, aes(x = order, y = final_height, fill = order1)) +
  geom_boxplot() +
  ylab("Final Seedling Height (cm)") +
  xlab("") +
  ggtitle("Sugar Maple") +
  scale_fill_brewer(palette = "Set1", direction = 1, name = "Treatment") +
  facet_grid(~ order, scales = "free_x") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        strip.text = element_text(size=12)) +
  theme(axis.text.x = element_blank(),
        axis.text=element_text(size=10),
        axis.ticks.x=element_blank(),
        axis.title=element_text(size=10)) +
  stat_summary(geom = 'text', label = symbol3, fun.y = max, vjust = -1, position = position_dodge(width = 0.75))
ac_g_h

symbol2 = c("a", "b", "b", "a", "a", "b", "a", "b", "b", "b", "b", "b")
fa_g_bio = ggplot(fagr_green, aes(x = order, y = end_biomass, fill = order1)) +
  geom_boxplot() +
  ylab("Final Dry Biomass (g)") +
  xlab("") +
  ggtitle("American Beech") +
  scale_fill_brewer(palette = "Set1", direction = 1, name = "Treatment") +
  facet_grid(~ order, scales = "free_x") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        strip.text = element_text(size=12)) +
  theme(axis.text.x = element_blank(),
        axis.text=element_text(size=10),
        axis.ticks.x=element_blank(),
        axis.title=element_text(size=10)) +
  stat_summary(geom = 'text', label = symbol2, fun.y = max, vjust = -1, position = position_dodge(width = 0.75))
fa_g_bio

symbol4 = c("a", "b", "b", "b", "b", "b", "b", "b", "b", "b", "b", "b")
fa_g_h = ggplot(fagr_green, aes(x = order, y = final_height, fill = order1)) +
  geom_boxplot() +
  ylab("Final Seedling Height (cm)") +
  xlab("") +
  ggtitle("American Beech") +
  scale_fill_brewer(palette = "Set1", direction = 1, name = "Treatment") +
  facet_grid(~ order, scales = "free_x") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        strip.text = element_text(size=12)) +
  theme(axis.text.x = element_blank(),
        axis.text=element_text(size=10),
        axis.ticks.x=element_blank(),
        axis.title=element_text(size=10)) +
  stat_summary(geom = 'text', label = symbol4, fun.y = max, vjust = -1, position = position_dodge(width = 0.75))
fa_g_h

#survival

acsa_sur = ddply(acsa_green, ~treat + source, summarise, surv = sum(survived))
fagr_sur = ddply(fagr_green, ~treat + source, summarise, surv = sum(survived))
write.csv(acsa_sur, "a_s.csv")
write.csv(fagr_sur, "f_s.csv")
as = read.csv("a_s.csv", header = TRUE)
fs = read.csv("f_s.csv", header = TRUE)
acsa_sur$order = factor(acsa_sur$source, levels=c("L", "M", "H"), 
                        labels = c("Hardwood soils", "Ecotone soils", "Spruce-fir soils"))
acsa_sur$order1 = factor(acsa_sur$treat, levels=c("B", "FA", "FC", "PC"),
                         labels = c("Innoculated", "Sterilized Field", "Field", "Potting Mix"))
fagr_sur$order = factor(fagr_sur$source, levels=c("L", "M", "H"),
                        labels = c("Hardwood soils", "Ecotone soils", "Spruce-fir soils"))
fagr_sur$order1 = factor

as$order = factor(as$source, levels=c("L", "M", "H"), 
                  labels = c("Hardwood soils", "Ecotone soils", "Spruce-fir soils"))
as$order1 = factor(as$treat, levels=c("B", "FA", "FC", "PC"),
                   labels = c("Innoculated", "Sterilized Field", "Field", "Potting Mix"))
fs$order = factor(fs$source, levels=c("L", "M", "H"),
                  labels = c("Hardwood soils", "Ecotone soils", "Spruce-fir soils"))
fs$order1 = factor(fs$treat, levels=c("B", "FA", "FC", "PC"),
                   labels = c("Innoculated", "Sterilized Field", "Field", "Potting Mix"))

symbol5 = c("a", "b", "b", "b", "a", "b", "b", "b", "b", "c", "c", "b")
ac_sur = ggplot(as, aes(x = order, y = surv, fill = order1)) +
  geom_bar(stat = "identity", position = "dodge") +
  ylab("Seedling Survival (%)") +
  xlab("") +
  ggtitle("Sugar Maple") +
  scale_fill_brewer(palette = "Set1", direction = 1, name = "Treatment") +
  facet_grid(~ order, scales = "free_x") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        strip.text = element_text(size=12)) +
  theme(axis.text.x = element_blank(),
        axis.text=element_text(size=10),
        axis.ticks.x=element_blank(),
        axis.title=element_text(size=10)) +
  stat_summary(geom = 'text', label = symbol5, fun.y = max, vjust = -1, position = position_dodge(width = 0.9))
ac_sur

symbol6 = c("a", "c", "b", "b", "ab", "d", "d", "c", "c", "d", "d", "c")
fa_sur = ggplot(fs, aes(x = order, y = surv, fill = order1)) +
  geom_bar(stat = "identity", position = "dodge") +
  ylab("Seedling Survival (%)") +
  xlab("") +
  ggtitle("American Beech") +
  scale_fill_brewer(palette = "Set1", direction = 1, name = "Treatment") +
  facet_grid(~ order, scales = "free_x") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        strip.text = element_text(size=12)) +
  theme(axis.text.x = element_blank(),
        axis.text=element_text(size=10),
        axis.ticks.x=element_blank(),
        axis.title=element_text(size=10)) +
  stat_summary(geom = 'text', label = symbol6, fun.y = max, vjust = -1, position = position_dodge(width = 0.9))
fa_sur

tiff(file = "green_grow.tif", width = 20, height = 12, units = 'in', res = 600, pointsize = 11)

green_all = plot_grid(ac_sur, ac_g_bio, ac_g_h, fa_sur, fa_g_bio, fa_g_h, nrow = 2, ncol = 3)
green_all

dev.off()

###
#col plots together

green_my$col_perc1 = as.numeric(green_my$col_perc1)
green_my$col_perc_se = as.numeric(green_my$col_perc_se)
green_my[1, 3] = 21.17
green_my[8, 3] = 9.65
green_my[9, 3] = 9.66

green_my$order = factor(green_my$site, levels=c("L", "M", "H"),
                        labels = c("Hardwood soils", "Ecotone soils", "Spruce-fir soils"))
green_my$order1 = factor(green_my$mountain, levels=c("FC", "B", "FA", "PC"))
green_my$order2 = factor(green_my$order1, levels=c("FC", "B", "FA", "PC"))
symbol_a = c("a", "a", "d", "d", "b", "b", "d", "d", "c", "c", "d", "d")

green_amf = ggplot(green_my, aes(x = order2, y = col_perc1, fill = order2)) +
  geom_errorbar(stat = "identity", aes(ymin=col_perc1-col_perc_se, ymax=col_perc1+col_perc_se), position = position_dodge(width = 0.90), width = 0.2, size = 1) +
  geom_point(stat = "identity", aes(color = order1), size = 5.5, pch = 23, position = position_dodge(width = 0.90)) +
  ylab("Root Length Colonized (%)") +
  xlab("") +
  ggtitle("Sugar Maple") +
  ylim(0, 30) +
  scale_x_discrete(labels=c("FC" = "Field", "B" = "Inoculated",
                            "FA" = "Sterilized Field",
                            "PC" = "Potting Mix")) +
  facet_wrap(~order, scales = "free_x", ncol = 3) +
  #facet_grid(~order, scales = "free_x") +
  scale_fill_manual(values = c("#fbe45b", "#2c7c94", "#a65852", "#a6d0c8")) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_blank(),
        strip.text.x = element_text(size = 20),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        plot.title = element_text(size=24),
        axis.ticks.x=element_blank(),
        legend.position = "none", 
        legend.title = element_blank(), 
        legend.text = element_blank()) +
  stat_summary(geom = 'text', label = symbol_a, fun.y = max, vjust = -2, position = position_dodge(width = 0.9), size = 9)
options(repr.plot.width =9, repr.plot.height =9)
green_amf

#greenhouse emf

emg1$order = factor(emg1$site, levels=c("L", "M", "H"),
                    labels = c("Hardwood soils", "Ecotone soils", "Spruce-fir soils"))
emg1$order1 = factor(emg1$mountain, levels=c("FC", "B", "FA", "PC"))
emg1$order2 = factor(emg1$order1, levels=c("FC", "B", "FA", "PC"))
symbol_f = c("a", "a", "b", "b", "a", "a", "b", "b", "a", "a", "b", "b")

green_emf = ggplot(emg1, aes(x = order2, y = col_perc, group = order2)) +
  geom_errorbar(stat = "identity", aes(ymin=(col_perc-SE)-1, ymax=(col_perc+SE)+1), position = position_dodge(width = 0.90), width = 0.2, size = 1) +
  geom_point(stat = "identity", color = "black", aes(fill=factor(order2)), size = 5.5, pch = 23, position = position_dodge(width = 0.90)) +
  scale_fill_manual(values = c("#fbe45b", "#2c7c94", "#a65852", "#a6d0c8"), name = "Treatment",
                    labels = c("Inoculated", "Sterilized Field", "Field", "Potting Mix")) +
  ylab("") +
  xlab("") +
  ggtitle("American beech") +
  ylim(0, 60) +
  scale_x_discrete(labels=c("FC" = "Field", "B" = "Inoculated",
                            "FA" = "Sterilized Field",
                            "PC" = "Potting Mix")) +
  facet_wrap(~order, scales = "free_x", ncol = 3) +
  #facet_grid(~order, scales = "free_x") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_blank(),
        strip.text.x = element_text(size = 20),
        axis.text=element_text(size=20),
        axis.ticks.x=element_blank(),
        axis.title=element_text(size=20),
        plot.title = element_text(size=24),
        legend.title = element_text(size=20), 
        legend.text = element_text(size=20),
        legend.position = "none") +
  stat_summary(geom = 'text', label = symbol_f, fun.y = max, vjust = -2, position = position_dodge(width = 0.9), size = 9)
options(repr.plot.width =9, repr.plot.height =9)
green_emf

combined_plot_green<-plot_grid(green_amf, green_emf,ncol=2, labels = c("A", "B"))
combined_plot_green
# extract legend from plot1
legend <- get_legend(
  green_emf +
    guides(color = guide_legend(nrow = 4)) +
    theme(legend.position = "left")
)

x.grob <- textGrob("Soil Condition Treatment", 
                   gp=gpar(fontface="bold", col="black", fontsize=20))
col_green = grid.arrange(arrangeGrob(combined_plot_green, bottom = x.grob))
col_green

tiff(file = "myco_green.tif", width = 14, height = 7, units = 'in', res = 600, pointsize = 11)
col_green = grid.arrange(arrangeGrob(combined_plot_green, bottom = x.grob))
col_green

dev.off()

# Combine combined plot and legend using plot_grid()
a1 = plot_grid(combined_plot_green, legend ,ncol=2,rel_widths = c(3, .4))
a1

#trans amf colonization

amf_trans1$order = factor(amf_trans1$site, levels=c("L", "M", "H"),
                          labels = c("Hardwood soils", "Ecotone soils", "Spruce-fir soils"))
amf_trans1$order1 = factor(amf_trans1$treat, levels=c("L", "M", "H"))
symbol_ta = c("a", "a", "b", "ab", "b", "bc", "ab", "c", "d")

amf_trans_plot = ggplot(amf_trans1, aes(x = order1, y = col, fill = order1)) +
  geom_errorbar(stat = "identity", aes(ymin=col-SEh, ymax=col+SEh), position = position_dodge(width = 0.90), width = 0.2, size = 1) +
  geom_point(stat = "identity", aes(color = order1), size = 6, pch = 23, position = position_dodge(width = 0.90)) +
  scale_fill_manual(values = c("#2c7c94", "#fbe45b", "#a65852")) +  
  ylab("Colonization (%)") +
  ggtitle("Sugar Maple") +
  ylim(0, 70) +
  xlab("") +
  scale_x_discrete(labels=c("L" = "Hardwoods", "M" = "Ecotone",
                            "H" = "Spruce-fir")) +
  facet_wrap(~order, scales = "free_x", ncol = 3) +
  #facet_grid(cols = vars(order)) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.x = element_text(size = 20),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        plot.title = element_text(size=24),
        legend.position = "none",
        legend.title = element_text(size=20), 
        legend.text = element_text(size=20)) +
  stat_summary(geom = 'text', label = symbol_ta, fun.y = max, vjust = -2, position = position_dodge(width = 0.90), size = 9)
amf_trans_plot 

#trans emf colonization

emf_trans$order = factor(emf_trans$site, levels=c("L", "M", "H"),
                         labels = c("Hardwood soils", "Ecotone soils", "Spruce-fir soils"))
emf_trans$order1 = factor(emf_trans$source, levels=c("L", "M", "H"))
symbol_tf = c("a", "a", "a", "a", "a", "a", "a", "a", "a")

emf_trans_plot = ggplot(emf_trans, aes(x = order1, y = col), fill = order1) +
  geom_errorbar(stat = "identity", aes(ymin=(col-SEh)-1, ymax=(col+SEh)+1), position = position_dodge(width = 0.90), width = 0.2, size = 1) +
  geom_point(stat = "identity", aes(fill = order1), size = 6, pch = 23, position = position_dodge(width = 0.90)) +
  scale_fill_manual(values = c("#2c7c94","#fbe45b", "#a65852"), name = "Planting Site",
                    labels = c("Hardwoods", "Ecotone", "Spruce-fir")) +
  ylab("") +
  ggtitle("American beech") +
  ylim(0, 70) +
  xlab("") +
  scale_x_discrete(labels=c("L" = "Hardwoods", "M" = "Ecotone",
                            "H" = "Spruce-fir")) +
  facet_wrap(~order, scales = "free_x", ncol = 3) +
  #facet_grid(cols = vars(order)) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.x = element_text(size = 20),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        plot.title = element_text(size=24),
        legend.title = element_text(size=20), 
        legend.text = element_text(size=20),
        legend.position = "none") +
  stat_summary(geom = 'text', label = symbol_tf, fun.y = max, vjust = -1.5, position = position_dodge(width = 0.90), size = 9)
emf_trans_plot 

combined_plot_trans<-plot_grid(amf_trans_plot, emf_trans_plot,ncol=2, labels = c("A", "B"))
combined_plot_trans
x.grob <- textGrob("Soil Translocation Site", 
                   gp=gpar(fontface="bold", col="black", fontsize=20))
col_trans = grid.arrange(arrangeGrob(combined_plot_trans, bottom = x.grob))
tiff(file = "myco_trans.tif", width = 14, height = 7, units = 'in', res = 600, pointsize = 11)
col_trans = grid.arrange(arrangeGrob(combined_plot_trans, bottom = x.grob))
dev.off()

# extract legend from plot1
legend <- get_legend(
  emf_trans_plot +
    guides(color = guide_legend(nrow = 3)) +
    theme(legend.position = "left")
)

# Combine combined plot and legend using plot_grid()
a2 = plot_grid(combined_plot_trans, legend ,ncol=2,rel_widths = c(3, .4))
a2

#myco colonization for net experiment

net_my$order1 = factor(net_my$treat, levels=c("N", "D", "C"))
net_my$order2 = factor(net_my$site, levels=c("L", "M", "H"),
                       labels = c("Hardwoods", "Ecotone", "Spruce-fir"))

#amf network plot
amf_net_plot = ggplot(net_my, aes(x = order1, y = col_perc, group = order1)) +
  geom_errorbar(stat = "identity", aes(ymin=col_perc-col_perc_se, ymax=col_perc+col_perc_se), position = position_dodge(width = 0.90), width = 0.2, size = 1) +
  geom_point(stat = "identity", aes(fill = order1), size = 6, pch = 23, position = position_dodge(width = 0.90)) +
  scale_fill_manual(values = c("#fbe45b", "#2c7c94", "#a65852"),
                    labels = c("Network", "Disrupted", "Pot Control")) +
  ylab("Root Length Colonized (%)") +
  xlab("") +
  ggtitle("Sugar Maple") +
  ylim(0, 60) +
  scale_x_discrete(labels=c("N" = "Network", "D" = "Disrupted",
                            "C" = "Pot Control")) +
  facet_wrap(~order2, scales = "free_x", ncol = 3) +
  #facet_grid(~ order2, scales = "free_x") +
  theme(axis.text.x = element_blank()) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        strip.text.x = element_text(size = 20),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        plot.title = element_text(size=24),
        legend.position = "none",
        legend.title = element_text(size=20), 
        legend.text = element_text(size=20))
#stat_summary(geom = 'text', label = plot.dat$sig, fun.y = max, vjust = -1, position = position_dodge(width = 0.75))
amf_net_plot 

#emf network plot

emf_net1$order1 = factor(emf_net1$source, levels=c("N", "D", "C"))
emf_net1$order2 = factor(emf_net1$site, levels=c("L", "M", "H"),
                         labels = c("Hardwoods", "Ecotone", "Spruce-fir"))

emf_net_plot = ggplot(emf_net1, aes(x = order1, y = col, group = order1)) +
  geom_errorbar(stat = "identity", aes(ymin=(col-SEh)-1, ymax=(col+SEh)+1), position = position_dodge(width = 0.90), width = 0.2, size = 1) +
  geom_point(stat = "identity", aes(fill = order1), size = 6, pch = 23, position = position_dodge(width = 0.90)) +
  scale_fill_manual(values = c("#fbe45b", "#2c7c94", "#a65852"), name = "Treatment",
                    labels = c("Networked", "Disrupted", "Control")) +
  ylab("") +
  xlab("") +
  ggtitle("American beech") +
  ylim(0, 60) +
  scale_x_discrete(labels=c("N" = "Network", "D" = "Disrupted",
                            "C" = "Pot Control")) +
  facet_wrap(~order2, scales = "free_x", ncol = 3) +
  #facet_grid(~order2, scales = "free_x") +
  theme(axis.text.x = element_blank()) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        strip.text.x = element_text(size = 20),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        plot.title = element_text(size=24),
        legend.title = element_text(size=20), 
        legend.text = element_text(size=20),
        legend.position = "none")
#stat_summary(geom = 'text', label = plot.dat$sig, fun.y = max, vjust = -1, position = position_dodge(width = 0.75))
emf_net_plot

combined_plot_net<-plot_grid(amf_net_plot, emf_net_plot,ncol=2, labels = c("A", "B"))
combined_plot_net
# extract legend from plot1
legend <- get_legend(
  emf_net_plot +
    guides(color = guide_legend(nrow = 3)) +
    theme(legend.position = "left")
)

# Combine combined plot and legend using plot_grid()
a3 = plot_grid(combined_plot_net, legend ,ncol=2,rel_widths = c(3, .4))
a3

myco_col = plot_grid(a1, a2, a3, ncol = 1, nrow = 3)
myco_col

x.grob <- textGrob("CMN Treatment", 
                   gp=gpar(fontface="bold", col="black", fontsize=20))
col_net = grid.arrange(arrangeGrob(combined_plot_net, bottom = x.grob))
col_net

tiff(file = "myco_net.tif", width = 14, height = 7, units = 'in', res = 600, pointsize = 11)
col_net = grid.arrange(arrangeGrob(combined_plot_net, bottom = x.grob))
col_net

dev.off()

###
#herb, path, browse
install.packages("wesanderson")
library(wesanderson)

net_all1$order = factor(net_all1$site, levels=c("L", "M", "H"),
                        labels = c("Hardwoods", "Ecotone", "Spruce-fir"))

symbol = c("a", "a", "b", "a", "b", "b")
browse_net = ggplot(net_all1, aes(x = order, y = browse_m, fill = spp)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("blue", "orange")) +
  ylab("Browse (% of pots)") +
  xlab("") +
  ggtitle("") +
  ylim(0,50) +
  theme_bw() +
  scale_x_discrete(labels=c("L" = "Hardwoods", "M" = "Ecotone",
                            "H" = "Spruce-fir")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=22),
        axis.title=element_text(size=22),
        plot.title = element_text(size=24),
        legend.title = element_text(size=22), 
        legend.text = element_text(size=22),
        axis.text.x = element_blank(),
        legend.position = "none") +
  stat_summary(geom = 'text', label = symbol, fun = max, vjust = -1, position = position_dodge(width = 0.9), size = 8)
browse_net

symbol1 = c("a", "a", "a", "a", "a", "a")
path_net = ggplot(net_all1, aes(x = order, y = path_m, fill = spp)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("blue", "orange")) +
  ylab("Root Pathogen Presence (% of pots)") +
  xlab("") +
  ggtitle("") +
  ylim(0,50) +
  theme_bw() +
  scale_x_discrete(labels=c("L" = "Hardwoods", "M" = "Ecotone",
                            "H" = "Spruce-fir")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=22),
        axis.title=element_text(size=22),
        legend.title = element_text(size=22), 
        legend.text = element_text(size=22),
        axis.text.x = element_blank(),
        legend.position = "none") +
  stat_summary(geom = 'text', label = symbol1, fun = max, vjust = -4.5, position = position_dodge(width = 0.9), size = 8)
path_net

symbol2 = c("a", "c", "ab", "b", "b", "b")
herb_net = ggplot(net_all1, aes(x = order, y = herb_m, fill = spp)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("#fbe45b", "#2c7c94"), labels = c("Sugar maple", "American beech")) +
  labs(fill = "Species") +
  ylab("Leaf Herbivory (% of pots)") +
  xlab("") +
  ggtitle("") +
  ylim(0,95) +
  theme_bw() +
  scale_x_discrete(labels=c("L" = "Hardwoods", "M" = "Ecotone",
                            "H" = "Spruce-fir")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=22),
        axis.title=element_text(size=22),
        legend.title = element_text(size=22), 
        legend.text = element_text(size=22)) +
  theme(legend.position="bottom") +
  stat_summary(geom = 'text', label = symbol2, fun = max, vjust = -0.5, position = position_dodge(width = 0.9), size = 8)
herb_net

#

trans_all1$order = factor(trans_all1$site, levels=c("L", "M", "H"),
                          labels = c("Hardwoods", "Ecotone", "Spruce-fir"))

symbol = c("a", "a", "a", "a", "a", "a")
browse_t = ggplot(trans_all1, aes(x = order, y = browse1, fill = spp)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("#fbe45b", "#2c7c94"), name = "Species",
                    label = c("Sugar maple", "American beech")) +
  ylab("") +
  xlab("") +
  ggtitle("Transplant Experiment") +
  ylim(0,50) +
  theme_bw() +
  scale_x_discrete(labels=c("L" = "Hardwoods", "M" = "Ecotone",
                            "H" = "Spruce-fir")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=22),
        axis.title=element_text(size=22),
        plot.title = element_text(size=24),
        legend.title = element_text(size=22), 
        legend.text = element_text(size = 22),
        axis.text.x = element_blank(),
        legend.position = c(0.75,0.85)) +
  stat_summary(geom = 'text', label = symbol, fun = max, vjust = -4.5, position = position_dodge(width = 0.9), size = 8)
browse_t

symbol1 = c("a", "a", "a", "a", "a", "a")
path_t = ggplot(trans_all1, aes(x = order, y = root_path1, fill = spp)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("#fbe45b", "#2c7c94")) +
  ylab("") +
  xlab("") +
  ggtitle("") +
  ylim(0,50) +
  theme_bw() +
  scale_x_discrete(labels=c("L" = "Hardwoods", "M" = "Ecotone",
                            "H" = "Spruce-fir")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=22),
        axis.title=element_text(size=22),
        legend.title = element_text(size=22), 
        legend.text = element_text(size=22),
        axis.text.x = element_blank(),
        legend.position = "none") +
  stat_summary(geom = 'text', label = symbol1, fun = max, vjust = -4.5, position = position_dodge(width = 0.9), size = 8)
path_t

symbol2 = c("a", "b", "b", "b", "b", "b")
herb_t = ggplot(trans_all1, aes(x = order, y = herb1, fill = spp)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("#fbe45b", "#2c7c94")) +
  ylab("") +
  xlab("") +
  ggtitle("") +
  ylim(0,100) +
  theme_bw() +
  scale_x_discrete(labels=c("L" = "Hardwoods", "M" = "Ecotone",
                            "H" = "Spruce-fir")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=22),
        axis.title=element_text(size=22),
        legend.title = element_text(size=22), 
        legend.text = element_text(size=22),
        legend.position = "none") +
  stat_summary(geom = 'text', label = symbol2, fun = max, vjust = -3.5, position = position_dodge(width = 0.9), size = 8)
herb_t

###altered 10/7
pots_bio = plot_grid(browse_net, herb_net, nrow = 2, ncol = 1)
pots_bio

tiff(file = "pots_bio_note.tif", width = 7, height = 11, units = 'in', res = 600, pointsize = 11)

pots_bio

dev.off()

###
#combined growth

symbol13 = c("a", "a", "a", "a", "a", "a", "a", "a", "a")
rgr_plota1 = ggplot(trans_s_ea1, aes(x = order, y = rhgr, fill = order1)) +
  geom_boxplot() +
  xlab("") +
  ylab(bquote('RHGR '(cm/day^-1))) +
  ggtitle("Sugar maple") +
  scale_fill_manual(values = c("#fbe45b", "#2c7c94", "#a65852"), name = "Soil Source",
                    label = c("Hardwood soils", "Ecotone soils", "Spruce-fir soils")) +
  facet_grid(~ order, scales = "free_x") +
  ylim(0.008,0.014) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        strip.text = element_text(size=20)) +
  theme(axis.text.x = element_blank(),
        axis.text=element_text(size=20),
        axis.ticks.x=element_blank(),
        legend.title = element_text(size=26), 
        plot.title = element_text(size = 24),
        legend.text = element_text(size=24),
        legend.position = "none",
        axis.title=element_text(size=20)) +
  stat_summary(geom = 'text', label = symbol13, fun.y = max, vjust = -1, position = position_dodge(width = 0.75), size = 9)
rgr_plota1

#fagr
rgr_plotf1 = ggplot(trans_s_ef1, aes(x = order, y = rhgr, fill = order1)) +
  geom_boxplot() +
  ylab("") +
  xlab("") +
  ggtitle("American Beech") +
  scale_fill_manual(values = c("#fbe45b", "#2c7c94", "#a65852")) +
  facet_grid(~ order, scales = "free_x") +
  ylim(0.007,0.018) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        strip.text = element_text(size=20)) +
  theme(axis.text.x = element_blank(),
        axis.text=element_text(size=20),
        axis.ticks.x=element_blank(),
        legend.position = "none",
        plot.title = element_text(size = 24),
        axis.title=element_text(size=20)) +
  stat_summary(geom = 'text', label = symbol13, fun.y = max, vjust = -1, position = position_dodge(width = 0.75), size = 9)
rgr_plotf1

#ylab(bquote('RHGR '(cm/day^-1)))

symbol30 = c("a", "a", "a", "a", "a", "a", "a", "a", "a")
net_rgr_plota1 = ggplot(net_s_ea1, aes(x = order, y = rhgr, fill = order1)) +
  geom_boxplot() +
  ylab(bquote('RHGR '(cm/day^-1))) +
  xlab("") +
  ggtitle("Sugar maple") +
  scale_fill_manual(values = c("#fbe45b", "#2c7c94", "#a65852"), name = "Treatment",
                    label = c("Network", "Disrupted", "Pot Control")) +
  facet_grid(~ order, scales = "free_x") +
  #ylim(0.005,0.015) +
  scale_y_continuous(limits = c(0.006,0.016),
                     labels = scales::number_format(accuracy = 0.001)) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        strip.text = element_text(size=20)) +
  theme(axis.text.x = element_blank(),
        axis.text=element_text(size=20),
        axis.ticks.x=element_blank(),
        legend.title = element_text(size=26), 
        plot.title = element_text(size = 24),
        legend.text = element_text(size=24),
        legend.position = "none",
        axis.title=element_text(size=20)) +
  stat_summary(geom = 'text', label = symbol30, fun.y = max, vjust = -1, position = position_dodge(width = 0.75), size = 9)
net_rgr_plota1

#ylab(bquote('RHGR '(cm/day^-1)))

#fagr
symbol30 = c("a", "a", "a", "a", "a", "a", "b", "b", "b")
net_rgr_plotf1 = ggplot(net_s_ef1, aes(x = order, y = rhgr, fill = order1)) +
  geom_boxplot() +
  ylab("") +
  xlab("") +
  ggtitle("American Beech") +
  scale_fill_manual(values = c("#fbe45b", "#2c7c94", "#a65852")) +
  facet_grid(~ order, scales = "free_x") +
  ylim(0.007, 0.015) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        strip.text = element_text(size=20)) +
  theme(axis.text.x = element_blank(),
        axis.text=element_text(size=20),
        axis.ticks.x=element_blank(),
        legend.position = "none",
        plot.title = element_text(size = 24),
        axis.title=element_text(size=20)) +
  stat_summary(geom = 'text', label = symbol30, fun.y = max, vjust = -1, position = position_dodge(width = 0.75), size = 9)
net_rgr_plotf1

#acsa root/shoot
symbol14 = c("a", "a", "a", "a", "a", "a", "b", "b", "b")
rs_plota1 = ggplot(trans_alla, aes(x = order, y = root_shoot_ratio_avg, fill = order1)) +
  geom_boxplot() +
  ylab("Root-to-Shoot Ratio") +
  ggtitle("") +
  xlab("") +
  scale_fill_manual(values = c("#fbe45b", "#2c7c94", "#a65852")) +
  facet_grid(~ order, scales = "free_x") +
  scale_y_continuous(labels = function(x) format(x, nsmall = 1)) +
  expand_limits(y=c(0.5, 4)) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        strip.text = element_text(size=20)) +
  theme(axis.text.x = element_blank(),
        axis.text=element_text(size=22),
        plot.title = element_text(size = 24),
        axis.ticks.x=element_blank(),
        legend.position = "none",
        axis.title=element_text(size=20)) +
  stat_summary(geom = 'text', label = symbol14, fun.y = max, vjust = -1, position = position_dodge(width = 0.75), size = 9)
rs_plota1

#fagr shoot/root

symbol15 = c("a", "a", "a", "a", "a", "a", "a", "a", "a")
rs_plotf1 = ggplot(trans_allf, aes(x = order, y = root_shoot_ratio_avg, fill = order1)) +
  geom_boxplot() +
  ylab("") +
  ggtitle("") +
  xlab("") +
  scale_fill_manual(values = c("#fbe45b", "#2c7c94", "#a65852")) +
  facet_grid(~ order, scales = "free_x") +
  ylim(0.5,2.5) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        strip.text = element_text(size=20)) +
  theme(axis.text.x = element_blank(),
        axis.text=element_text(size=20),
        plot.title = element_text(size = 24),
        axis.ticks.x=element_blank(),
        legend.position = "none",
        axis.title=element_text(size=20)) +
  stat_summary(geom = 'text', label = symbol15, fun.y = max, vjust = -1, position = position_dodge(width = 0.75), size = 9)
rs_plotf1

#acsa total biomass

symbol16 = c("a", "a", "a", "a", "a", "a", "a", "a", "a")
tb_plota1 = ggplot(trans_alla, aes(x = order, y = total_bio_avg, fill = order1)) +
  geom_boxplot() +
  ylab("Total Dry Biomass (g)") +
  ggtitle("") +
  xlab("") +
  scale_fill_manual(values = c("#fbe45b", "#2c7c94", "#a65852")) +
  facet_grid(~ order, scales = "free_x") +
  ylim(0,1) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        strip.text = element_text(size=20)) +
  theme(axis.text.x = element_blank(),
        axis.text=element_text(size=22),
        axis.ticks.x=element_blank(),
        legend.title = element_text(size=20), 
        plot.title = element_text(size = 24),
        legend.text = element_text(size=20),
        legend.position = "none",
        axis.title=element_text(size=20)) +
  stat_summary(geom = 'text', label = symbol16, fun.y = max, vjust = -1, position = position_dodge(width = 0.75), size = 9)
tb_plota1

#fagr total biomass
symbol16 = c("a", "a", "a", "a", "a", "a", "a", "a", "a")
tb_plotf1 = ggplot(trans_allf, aes(x = order, y = total_bio_avg, fill = order1)) +
  geom_boxplot() +
  ylab("") +
  xlab("") +
  ggtitle("") +
  scale_fill_manual(values = c("#fbe45b", "#2c7c94", "#a65852")) +
  facet_grid(~ order, scales = "free_x") +
  scale_y_continuous(labels = function(x) format(x, nsmall = 1)) +
  expand_limits(y=c(0, 3)) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        strip.text = element_text(size=20)) +
  theme(axis.text.x = element_blank(),
        axis.text=element_text(size=20),
        axis.ticks.x=element_blank(),
        legend.title = element_text(size=20), 
        legend.text = element_text(size=20),
        plot.title = element_text(size = 24),
        legend.position = "none",
        axis.title=element_text(size=20)) +
  stat_summary(geom = 'text', label = symbol16, fun.y = max, vjust = -1, position = position_dodge(width = 0.75), size = 9)
tb_plotf1

#acsa root/shoot
symbol17 = c("a", "a", "a", "a", "a", "a", "b", "b", "b")
net_rs_plota1 = ggplot(net_alla, aes(x = order, y = root_shoot_ratio_avg, fill = order1)) +
  geom_boxplot() +
  ylab("Root-to-Shoot Ratio") +
  xlab("") +
  ggtitle("") +
  scale_fill_manual(values = c("#fbe45b", "#2c7c94", "#a65852")) +
  facet_grid(~ order, scales = "free_x") +
  scale_y_continuous(labels = function(x) format(x, nsmall = 1)) +
  expand_limits(y=c(0.5, 5.0)) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        strip.text = element_text(size=20)) +
  theme(axis.text.x = element_blank(),
        axis.text=element_text(size=20),
        axis.ticks.x=element_blank(),
        plot.title = element_text(size = 24),
        legend.position = "none", 
        axis.title=element_text(size=20)) +
  stat_summary(geom = 'text', label = symbol17, fun.y = max, vjust = -1, position = position_dodge(width = 0.75), size = 9)
net_rs_plota1

#fagr root/shoot
symbol18 = c("a", "a", "a", "a", "a", "a", "a", "a", "a")
net_rs_plotf1 = ggplot(net_allf, aes(x = order, y = root_shoot_ratio_avg, fill = order1)) +
  geom_boxplot() +
  ylab("") +
  xlab("") +
  ggtitle("") +
  scale_fill_manual(values = c("#fbe45b", "#2c7c94", "#a65852")) +
  facet_grid(~ order, scales = "free_x") +
  scale_y_continuous(labels = function(x) format(x, nsmall = 1)) +
  expand_limits(y=c(0.5, 3.5)) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        strip.text = element_text(size=20)) +
  theme(axis.text.x = element_blank(),
        axis.text=element_text(size=20),
        axis.ticks.x=element_blank(),
        plot.title = element_text(size = 24),
        legend.position = "none",
        axis.title=element_text(size=20)) +
  stat_summary(geom = 'text', label = symbol18, fun.y = max, vjust = -1, position = position_dodge(width = 0.75), size = 9)
net_rs_plotf1

#acsa total biomass
symbol20 = c("a", "a", "a", "a", "a", "a", "a", "a", "a")
net_b_plota1 = ggplot(net_alla, aes(x = order, y = total_bio_avg, fill = order1)) +
  geom_boxplot() +
  ylab("Total Dry Biomass (g)") +
  xlab("") +
  ggtitle("") +
  scale_fill_manual(values = c("#fbe45b", "#2c7c94", "#a65852")) +
  facet_grid(~ order, scales = "free_x") +
  scale_y_continuous(labels = function(x) format(x, nsmall = 1)) +
  expand_limits(y=c(0, 15)) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        strip.text = element_text(size=20)) +
  theme(axis.text.x = element_blank(),
        axis.text=element_text(size=20),
        axis.ticks.x=element_blank(),
        legend.title = element_text(size=20), 
        plot.title = element_text(size = 24),
        legend.text = element_text(size=20),
        legend.position = "none",
        axis.title=element_text(size=20)) +
  stat_summary(geom = 'text', label = symbol20, fun.y = max, vjust = -1, position = position_dodge(width = 0.75), size = 9)
net_b_plota1

#fagr total biomass
symbol21 = c("a", "a", "a", "a", "a", "a", "a", "a", "a")
net_b_plotf1 = ggplot(net_allf, aes(x = order, y = total_bio_avg, fill = order1)) +
  geom_boxplot() +
  ylab("") +
  xlab("") +
  ggtitle("") +
  scale_fill_manual(values = c("#fbe45b", "#2c7c94", "#a65852")) +
  facet_grid(~ order, scales = "free_x") +
  ylim(0,10) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        strip.text = element_text(size=20)) +
  theme(axis.text.x = element_blank(),
        axis.text=element_text(size=20),
        axis.ticks.x=element_blank(),
        legend.title = element_text(size=20),
        plot.title = element_text(size = 24),
        legend.text = element_text(size=20),
        legend.position = "none",
        axis.title=element_text(size=20)) +
  stat_summary(geom = 'text', label = symbol21, fun.y = max, vjust = -1, position = position_dodge(width = 0.75), size = 9)
net_b_plotf1

#green
acsa_green$order = factor(acsa_green$source, levels=c("L", "M", "H"),
                          labels = c("Hardwood soils", "Ecotone soils", "Spruce-fir soils"))
acsa_green$order1 = factor(acsa_green$treat, levels=c("FC", "B", "FA", "PC"),
                           labels = c("Field", "Innoculated", "Sterilized Field", "Potting Mix"))
fagr_green$order = factor(fagr_green$source, levels=c("L", "M", "H"),
                          labels = c("Hardwood soils", "Ecotone soils", "Spruce-fir soils"))
fagr_green$order1 = factor(fagr_green$treat, levels=c("FC", "B", "FA", "PC"),
                           labels = c("Field", "Innoculated", "Sterilized Field", "Potting Mix"))
symbol1 = c("b", "a", "b", "b", "b", "a", "b", "b", "b", "b", "b", "ab")
ac_g_bio = ggplot(acsa_green, aes(x = order, y = end_biomass, fill = order1)) +
  geom_boxplot() +
  ylab("Final Dry Biomass (g)") +
  xlab("") +
  ggtitle("") +
  ylim(0,2) +
  scale_fill_manual(values = c("#fbe45b", "#2c7c94", "#a65852", "#a6d0c8")) +
  facet_grid(~ order, scales = "free_x") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        strip.text = element_text(size=20)) +
  theme(axis.text.x = element_blank(),
        axis.text=element_text(size=20),
        axis.ticks.x=element_blank(),
        plot.title = element_text(size = 24),
        legend.position = "none",
        axis.title=element_text(size=20)) +
  stat_summary(geom = 'text', label = symbol1, fun.y = max, vjust = -1, position = position_dodge(width = 0.75), size = 9)
ac_g_bio

install.packages("scales")
library(scales)

symbol3 = c("b", "a", "b", "ab", "b", "a", "b", "b", "b", "a", "b", "ab")
ac_g_h = ggplot(acsa_green, aes(x = order1, y = final_height, fill = order1)) +
  geom_boxplot() +
  ylab("Final Seedling Height (cm)") +
  xlab("") +
  ylim(0,20)+
  ggtitle("") +
  scale_fill_manual(values = c("#fbe45b", "#2c7c94", "#a65852", "#a6d0c8")) +
  facet_grid(~ order, scales = "free_x") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        strip.text = element_text(size=20)) +
  theme(axis.text.x = element_text(size = 20, angle = 45, vjust = 1, hjust = 1),
        axis.text=element_text(size=20),
        plot.title = element_text(size = 24),
        legend.title = element_text(size=20), 
        legend.text = element_text(size=20),
        legend.position = "none",
        axis.title=element_text(size=20)) +
  stat_summary(geom = 'text', label = symbol3, fun.y = max, vjust = -1, position = position_dodge(width = 0.75), size = 9)
options(repr.plot.width =9, repr.plot.height =9)
ac_g_h

symbol2 = c("b", "a", "b", "a", "a", "a", "b", "b", "b", "b", "b", "b")
fa_g_bio = ggplot(fagr_green, aes(x = order, y = end_biomass, fill = order1)) +
  geom_boxplot() +
  ylab("") +
  xlab("") +
  ggtitle("") +
  ylim(0,0.6) +
  scale_fill_manual(values = c("#fbe45b", "#2c7c94", "#a65852", "#a6d0c8")) +
  facet_grid(~ order, scales = "free_x") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        strip.text = element_text(size=20)) +
  theme(axis.text.x = element_blank(),
        axis.text=element_text(size=20),
        legend.position = "none",
        plot.title = element_text(size = 24),
        axis.title=element_text(size=20)) +
  stat_summary(geom = 'text', label = symbol2, fun.y = max, vjust = -1, position = position_dodge(width = 0.75), size = 9)
options(repr.plot.width =9, repr.plot.height =9)
fa_g_bio

symbol4 = c("b", "a", "b", "b", "b", "b", "b", "b", "b", "b", "b", "b")
fa_g_h = ggplot(fagr_green, aes(x = order1, y = final_height, fill = order1)) +
  geom_boxplot() +
  ylab("") +
  xlab("") +
  ylim(0,15) +
  ggtitle("") +
  scale_fill_manual(values = c("#fbe45b", "#2c7c94", "#a65852", "#a6d0c8")) +
  facet_grid(~ order, scales = "free_x") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        strip.text = element_text(size=20)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.text=element_text(size=20),
        legend.title = element_text(size=20), 
        legend.text = element_text(size=20),
        legend.position = "none", 
        plot.title = element_text(size = 24),
        axis.title=element_text(size=20)) +
  stat_summary(geom = 'text', label = symbol4, fun.y = max, vjust = -1, position = position_dodge(width = 0.75), size = 9)
options(repr.plot.width =9, repr.plot.height =9)
fa_g_h


#survival
acsa_sur$order = factor(acsa_sur$source, levels=c("L", "M", "H"), 
                        labels = c("Hardwood soils", "Ecotone soils", "Spruce-fir soils"))
acsa_sur$order1 = factor(acsa_sur$treat, levels=c("FC", "B", "FA", "PC"),
                         labels = c("Field", "Innoculated", "Sterilized Field", "Potting Mix"))
fagr_sur$order = factor(fagr_sur$source, levels=c("L", "M", "H"),
                        labels = c("Hardwood soils", "Ecotone soils", "Spruce-fir soils"))
fagr_sur$order1 = factor(acsa_sur$treat, levels=c("FC", "B", "FA", "PC"),
                         labels = c("Field", "Innoculated", "Sterilized Field", "Potting Mix"))

as$order = factor(as$source, levels=c("L", "M", "H"), 
                  labels = c("Hardwood soils", "Ecotone soils", "Spruce-fir soils"))
as$order1 = factor(as$treat, levels=c("FC", "B", "FA", "PC"),
                   labels = c("Field", "Innoculated", "Sterilized Field", "Potting Mix"))
fs$order = factor(fs$source, levels=c("L", "M", "H"),
                  labels = c("Hardwood soils", "Ecotone soils", "Spruce-fir soils"))
fs$order1 = factor(fs$treat, levels=c("FC", "B", "FA", "PC"),
                   labels = c("Field", "Innoculated", "Sterilized Field", "Potting Mix"))

symbol5 = c("b", "a", "b", "b", "b", "a", "b", "b", "c", "b", "c", "b")
ac_sur = ggplot(as, aes(x = order, y = surv, fill = order1)) +
  geom_bar(stat = "identity", position = "dodge") +
  ylab("Seedling Survival (%)") +
  xlab("") +
  ylim(0,60) +
  #ggtitle("Sugar Maple") +
  scale_fill_manual(values = c("#fbe45b", "#2c7c94", "#a65852", "#a6d0c8")) +
  facet_grid(~ order, scales = "free_x") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        strip.text = element_text(size=20)) +
  theme(axis.text.x = element_blank(),
        axis.text=element_text(size=20),
        axis.ticks.x=element_blank(),
        plot.title = element_text(size = 24),
        legend.position = "none",
        axis.title=element_text(size=20)) +
  stat_summary(geom = 'text', label = symbol5, fun.y = max, vjust = -0.5, position = position_dodge(width = 0.9), size = 9)
options(repr.plot.width =9, repr.plot.height =9)
ac_sur

symbol6 = c("b", "a", "b", "b", "c", "ab", "c", "c", "c", "bc", "c", "bc")
fa_sur = ggplot(fs, aes(x = order, y = surv, fill = order1)) +
  geom_bar(stat = "identity", position = "dodge") +
  ylab("") +
  xlab("") +
  ylim(0,40) +
  #ggtitle("American Beech") +
  scale_fill_manual(values = c("#fbe45b", "#2c7c94", "#a65852", "#a6d0c8"), name = "Treatment",
                    labels = c("Field", "Inoculated", "Sterilized Field", "Potting Mix")) +
  facet_grid(~ order, scales = "free_x") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        strip.text = element_text(size=20)) +
  theme(axis.text.x = element_blank(),
        axis.text=element_text(size=20),
        axis.ticks.x=element_blank(),
        plot.title = element_text(size=24),
        legend.title = element_text(size=26), 
        legend.text = element_text(size=24),
        legend.position = "none",
        axis.title=element_text(size=20)) +
  stat_summary(geom = 'text', label = symbol6, fun.y = max, vjust = -1, position = position_dodge(width = 0.9), size = 9)
options(repr.plot.width =9, repr.plot.height =9)
fa_sur

new = (green_amf | green_emf) / (ac_sur | fa_sur) / (ac_g_bio | fa_g_bio) / (ac_g_h | fa_g_h)
new

new1 = new + plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 24, hjust = -4))
new1

green_all0 = plot_grid(green_amf, green_emf, nrow = 2, ncol = 1, align = "v")
y.grob0 <- textGrob("Colonization (%)", 
                    gp=gpar(fontface="bold", col="black", fontsize=28), rot=90)
green_all1 = plot_grid(ac_sur, fa_sur, nrow = 2, ncol = 1, align = "v")
y.grob1 <- textGrob("Seedling Survival (%)", 
                    gp=gpar(fontface="bold", col="black", fontsize=28), rot=90)
green_all2 = grid.arrange(arrangeGrob(green_all1, left = y.grob1))

green_all3 = plot_grid(ac_g_bio, fa_g_bio, nrow = 2, ncol = 1, align = "v")
y.grob2 <- textGrob("Final Dry Biomass (g)", 
                    gp=gpar(fontface="bold", col="black", fontsize=28), rot=90)
green_all4 = grid.arrange(arrangeGrob(green_all3, left = y.grob2))
green_all4

green_all5 = plot_grid(ac_g_h, fa_g_h, nrow = 2, ncol = 1, align = "v")
y.grob3 <- textGrob("Final Seedling Height (cm)", 
                    gp=gpar(fontface="bold", col="black", fontsize=28), rot=90)
green_all6 = grid.arrange(arrangeGrob(new, left = y.grob3))


legend <- get_legend(
  fa_sur +
    guides(color = guide_legend(nrow = 3)) +
    theme(legend.position = "bottom")
)

#install.packages("patchwork")
library(patchwork)

tiff(file = "green_grow1.tif", width = 14.5, height = 22.5, units = 'in', res = 300, pointsize = 11)

# Combine combined plot and legend using plot_grid()
green_all7 = plot_grid(new1,
                       legend, ncol=1,
                       rel_heights = c(3,0.4))
green_all7

dev.off()


#install.packages("patchwork")
library(patchwork)

new2 = (rgr_plota1 | rgr_plotf1) / (rs_plota1 | rs_plotf1) / (tb_plota1 | tb_plotf1)
new2

tiff(file = "trans_grow1.tif", width = 24, height = 12, units = 'in', res = 600, pointsize = 12)

trans_grow = plot_grid(rgr_plota1, rs_plota1, tb_plota1,
                       rgr_plotf1, rs_plotf1, tb_plotf1,
                       ncol = 1, nrow = 6, align = "v")
trans_grow

dev.off()

trans_all1 = plot_grid(rgr_plota1, rgr_plotf1, nrow = 2, ncol = 1, align = "v")
y.grob1 <- textGrob(bquote('RHGR '(cm/day^-1)), 
                    gp=gpar(fontface="bold", col="black", fontsize=28), rot=90)
trans_all2 = grid.arrange(arrangeGrob(trans_all1, left = y.grob1))

trans_all3 = plot_grid(rs_plota1, rs_plotf1, nrow = 2, ncol = 1, align = "v")
y.grob2 <- textGrob("Root : Shoot", 
                    gp=gpar(fontface="bold", col="black", fontsize=28), rot=90)
trans_all4 = grid.arrange(arrangeGrob(trans_all3, left = y.grob2))
trans_all4

trans_all5 = plot_grid(tb_plota1, tb_plotf1, nrow = 2, ncol = 1, align = "v")
y.grob3 <- textGrob("Total Dry Biomass (g)", 
                    gp=gpar(fontface="bold", col="black", fontsize=28), rot=90)
trans_all6 = grid.arrange(arrangeGrob(trans_all5, left = y.grob3))

legend <- get_legend(
  rgr_plota1 +
    guides(color = guide_legend(nrow = 3)) +
    theme(legend.position = "bottom")
)

tiff(file = "trans_grow.tif", width = 14, height = 16, units = 'in', res = 600, pointsize = 11)

# Combine combined plot and legend using plot_grid()
trans_all7 = plot_grid(new1,
                       legend, ncol=1,rel_heights = c(3,0.4))
trans_all7

dev.off()

tiff(file = "green_new_grow.tif", width = 10, height = 27, units = 'in', res = 600, pointsize = 11)

green_all = plot_grid(green_all2, green_all4, green_all6, nrow = 3, ncol = 1, align = "vh")

green_all

dev.off()

tiff(file = "net_grow.tif", width = 24, height = 10, units = 'in', res = 600, pointsize = 11)

net_grow = plot_grid(net_rgr_plota1,  net_rs_plota1, net_b_plota1,
                     net_rgr_plotf1, net_rs_plotf1, net_b_plotf1,
                     nrow = 2, ncol = 3)
net_grow

dev.off()

net_all1 = plot_grid(net_rgr_plota1, net_rgr_plotf1, nrow = 2, ncol = 1, align = "v")
y.grob1 <- textGrob(bquote('RHGR '(cm/day^-1)), 
                    gp=gpar(fontface="bold", col="black", fontsize=28), rot=90)
net_all2 = grid.arrange(arrangeGrob(net_all1, left = y.grob1))

net_all3 = plot_grid(net_rs_plota1, net_rs_plotf1, nrow = 2, ncol = 1, align = "v")
y.grob2 <- textGrob("Root : Shoot", 
                    gp=gpar(fontface="bold", col="black", fontsize=28), rot=90)
net_all4 = grid.arrange(arrangeGrob(net_all3, left = y.grob2))
net_all4

net_all5 = plot_grid(net_b_plota1, net_b_plotf1, nrow = 2, ncol = 1, align = "v")
y.grob3 <- textGrob("Total Dry Biomass (g)", 
                    gp=gpar(fontface="bold", col="black", fontsize=28), rot=90)
net_all6 = grid.arrange(arrangeGrob(net_all5, left = y.grob3))

legend <- get_legend(
  net_rgr_plota1 +
    guides(color = guide_legend(nrow = 3)) +
    theme(legend.position = "bottom")
)


new2 = (net_rgr_plota1 | net_rgr_plotf1) / (net_rs_plota1 | net_rs_plotf1) / (net_b_plota1 | net_b_plotf1)
new2

tiff(file = "net_grow.tif", width = 14, height = 16, units = 'in', res = 600, pointsize = 11)

# Combine combined plot and legend using plot_grid()
net_all7 = plot_grid(new2,
                     legend, ncol=1,rel_heights = c(3, 0.4))
net_all7

dev.off()

#final fig
tiff(file = "all_grow.tif", width = 24, height = 24, units = 'in', res = 600, pointsize = 11)

all_grow = plot_grid(green_all, trans_all7, net_all7,
                     nrow = 1, ncol = 3)
all_grow

dev.off()

##