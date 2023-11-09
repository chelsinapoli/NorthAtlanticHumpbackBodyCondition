# Purpose -----------------------------------------------------------------

# purpose: to investigate the differences between regional body condition
# of the north Atlantic humpback whales

# CN

# created: nov 26 2022
# last modified: october 30 2023

# Library Setup -------------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(car)
library(smatr)
library(dplyr)
library(npreg)
library(stringr)
library(ggformula)
library(multcomp)
library(readxl)
library(ggsignif)
library(plyr)
library(coin)
library(FSA)
library(bootstrap)
library(sm)
library(lsr)
library(lme4)
# library(ggsmatr)

set.seed(03282020)


# Loading Data ------------------------------------------------------------

# setting the working directory to major data folder
all_data <- read.csv("./Napoli_et_al_RawData.csv")

# Add Year ----------------------------------------------------------------

# add year 

all_data$Date = format(as.Date(all_data$Date, "%m/%d/%y"), "20%y-%m-%d")
all_data$Year <- format(as.Date(all_data$Date),"%Y")

# Adding Juvenile v Adult -------------------------------------------------

# add categorization of juvenile vs adult 
all_data <- all_data %>%
  mutate(age_class = ifelse(Total_Length > 11.47, "Adult", "Juvenile")) # cut off of adult being 11.47 meters as used in Bierlich. (2022)


# Adding in Average Mass --------------------------------------------------

all_data <- all_data %>%
  mutate(mass = Body_Volume*1035.3)


# Regional BCI & BCI Plots  ---------------------------------------------------------------------

# setting group colors
# red for stony brook, blue for iceland, orange-ish for gulf of maine for unb, grey for historical, green for greenland
site_colors <- c("New York Bight" = "firebrick3", "Iceland" = "lightskyblue", "Gulf of Maine" = "goldenrod2", "Greenland" = "forestgreen")

# creating the line used as a guide for the plots 
vol_by_length <- nls(Body_Volume ~ a*(Total_Length^b),
                     data = all_data, 
                     start = list(a = 0.5, b = 0.2))

bci_equation <- lm(log(all_data$Body_Volume) ~ log(all_data$Total_Length)) # equation to find expected body volume from Christiansen et al (2020)

all_data$logfit <- fitted(bci_equation) # logfit 
all_data$fit <- exp(all_data$logfit) # expected value 
all_data$BCI <- (all_data$Body_Volume-all_data$fit)/all_data$fit # equation for BCI 

# Violin Plot for BCI
all_data$`Foraging Region` <- all_data$Foraging.Region

all_data$`Foraging Region` <- factor(all_data$`Foraging Region`, levels = c("New York Bight", "Gulf of Maine", "Greenland", "Iceland"))

BCI_violin <- ggplot(all_data, aes(x = `Foraging Region`, y = BCI, fill = `Foraging Region`, alpha = 0.75, color = `Foraging Region`)) + # all regional data, with foraging region as x groups and BCI on y axis. Colors determined by foraging group
  geom_violin(fill = "transparent", draw_quantiles = 0.5, color = "black") + # formatting for violin plot
  geom_violin() + # telling ggplot to make a violin plot
  geom_dotplot(binaxis = 'y', 
               stackratio = 1, 
               stackgroups = TRUE, 
               stackdir = 'centerwhole', 
               dotsize = 0.3, 
               alpha = 2, 
               binpositions = "bygroup", 
               stroke = NA,
               method = "histodot") +
  labs(x = "Foraging Region", y = "Body Condition Index") + # labels
  ggtitle("Body Condition Index by Region") + # title
  scale_fill_manual(values = site_colors) + # legend 
  scale_color_manual(values = site_colors) +  # legend 
  scale_alpha(guide = "none") +  # legend
  guides(fill="none") +
  guides(color = "none") +
  theme(legend.position = "none") + # title formatting
  geom_signif(comparisons = list(c("Gulf of Maine", "New York Bight")), map_signif_level=TRUE, textsize = 5, fontface = "bold", y_position =  0.62, aes(color = "black"), alpha = 0.8) + # making significance stars between NYB and GOM
  geom_signif(comparisons = list(c("Gulf of Maine", "Iceland")), map_signif_level=TRUE, textsize = 5, fontface = "bold", y_position =  0.56, aes(color = "black"), alpha = 0.8) + # making significance stars between GOM and Iceland
  geom_signif(comparisons = list(c("Gulf of Maine", "Greenland")), map_signif_level=TRUE, textsize = 5, fontface = "bold", y_position =  0.49, aes(color = "black"), alpha = 0.8) + # making significance stars between NYB and Greenland
  theme_bw() + # making theme black and white, removing the default grey background 
  theme(axis.text=element_text(size=12), title=element_text(size=12))
BCI_violin # plot the plot 

# ANCOVA using DOY, Year ---------------------------------------------------------------------

# ANOVA to control for age class
all_data$Foraging_Region <- as.factor(all_data$`Foraging Region`)
agemodel <- aov(BCI ~ age_class + Error(Foraging_Region), data = all_data)

# ANCOVA using DOY, Year ---------------------------------------------------------------------

# ANCOVA to control for DOY and assess any differences in our Body Volume Metric
all_data$Foraging_Region <- as.factor(all_data$`Foraging Region`)
ancova_resid <- aov(BCI ~ Foraging_Region + DOY + Year, data=all_data)
Anova(ancova_resid, type = "III")

# tukey post hoc 
postHocs_resid <- glht(ancova_resid, linfct = mcp(Foraging_Region = "Tukey"))
summary(postHocs_resid)

# effect size 
etaSquared(ancova_resid)

# Regional SMA ---------------------------------------------------------------------

sma_sites <- sma(mass ~ Total_Length +Foraging_Region, log="xy", data = all_data, multcomp = TRUE, multcompmethod = "adjusted") # sma with body volume and total length, and foraging region as a covariate
summary(sma_sites) # summary

# using base R to plot it, which is terrible :) 

# plotting SMA with legend
# note that base R assigns colors in alphabetical order, so they MUST be assigned Greenland, GOM, Iceland, NYB! Otherwise the colors will not match with the foraging region! 
#plot(sma_sites, 
#main = "Standardized Major Axes Regression by Foraging Region", 
#xlab = "Total Length (m)", ylab = "Body Volume (cubic meters)", 
#col = c("Greenland" = "forestgreen", "GOM" = "goldenrod2", "Iceland" = "skyblue2", "NYB" = "firebrick3"),
#cex.main=1.25, 
#cex.lab=1.3)
#legend("topleft", legend=c("NYB", "Iceland", "Gulf of Maine", "Greenland"),
#col=c("firebrick3", "skyblue2","goldenrod2", "forestgreen"), 
#lty=1, cex=0.8)

# using ggplot to plot it :) Adapted from the ggSMATR package! 

# function taken from ggsmatr package but that package is a little messed up and doesn't like the log = "xy" function part of SMAtr

# this function just makes everything a factor 
as_factorALL<- function(dataframe){
  
  df<- dataframe
  df[sapply(df, is.character)] <- lapply(df[sapply(df, is.character)],
                                         as.factor)
  return(df)
}

data<- as_factorALL(sma_sites$data) # this inputs data as the same data as the smatr package uses it 

# define colors
grps.df = unique(data[["Foraging_Region"]])
station_cols = scales::hue_pal()(length(grps.df))
names(station_cols) <- unique(data[["Foraging_Region"]])

fitsall <- sma_sites$groupsummary |> dplyr::select(group , r2, pval, Slope, Int)
fitsall <- fitsall |>  dplyr::mutate_if(is.numeric, round, digits = 3)
fitsall

# extraer coeficientes de SMA, vars:  dataframe  groups   xvar # direct comment from ggsmatr 
bb <- data.frame(stats::coef(sma_sites))

bb <- bb |>
  tibble::rownames_to_column(var = "group")

# calcular minimo y maximo del eje x para cada sitio # direct comment from ggsmatr
bb2 <- data |> 
  dplyr::select(Foraging_Region, Total_Length) |>  
  dplyr::group_by(Foraging_Region) |>
  dplyr::summarise(dplyr::across(.cols = Total_Length,
                                 .fns = list(min = min, max = max),
                                 .names = "{fn}_x"))
names(bb2)[1] <- 'group'


#unir dataframe
bb3 <- base::merge(bb,bb2, by= "group")

#calcular min y max de y con intercepto y elevacion
bb4 <- bb3 |>
  dplyr::mutate(min_y = (slope*min_x) + elevation) |>
  dplyr::mutate(max_y = (slope*max_x) + elevation)
data<- stats::na.omit(data)

df.pl <- cbind(bb4, fitsall[,-1])

print(df.pl[,c(1,8,9)])

df.pl.log <- df.pl

gg_sma_log_plot <- ggplot(data = sma_sites$data, ggplot2::aes(x = Total_Length, y = mass, color= Foraging_Region, fill= Foraging_Region)) +
  ggplot2::geom_point(data = sma_sites$data, aes(x = Total_Length, y = mass, color= Foraging_Region, fill= Foraging_Region), size=2, alpha = 0.5) +
  ggplot2::geom_segment(data = df.pl.log ,ggplot2::aes(x= min_x, xend= max_x, y=min_y, yend=max_y, colour= group),inherit.aes = FALSE, linewidth = 1) +
  scale_color_manual(values = site_colors) +
  labs(title = "(b)", 
       x = "Log (Total Length) (m)", 
       y = bquote(' Log (Mass) '(kg)),
       color = "Foraging Region", 
       fill = "Foraging Region") +
  theme_bw()

gg_sma_log_plot


# Violin Plot for BCI Historical Data

# Plotting Total Length v Body Volume 
tl_v_bv_plot <- ggplot(all_data, aes(x = Total_Length, y = Body_Volume)) + # all historical data with total length as x and body volume as y 
  geom_point(data = all_data %>% filter(Foraging.Region == "Iceland"), aes(x = Total_Length, y = Body_Volume, colour = "Iceland"), size = 2, alpha = 0.75) + # icelandic data combined and loaded in with the assigned color 
  geom_point(data = all_data %>% filter(Foraging.Region == "New York Bight"), aes(x = Total_Length, y = Body_Volume, colour = "New York Bight"), size = 2, alpha = 0.75) + # NYB data loaded in with assigned color 
  geom_point(data = all_data %>% filter(Foraging.Region == "Gulf of Maine"), aes(x = Total_Length, y = Body_Volume, colour = "Gulf of Maine"), size = 2, alpha = 0.75) + # AW and UNB data loaded in with assigned coor 
  geom_point(data = all_data %>% filter(Foraging.Region == "Greenland"), aes(x = Total_Length, y = Body_Volume, colour = "Greenland"), size = 2, alpha = 0.75) + # greenland data loaded in with the assigned color 
  geom_line(aes(y = exp(fitted(bci_equation)))) + # adding in the expected value line
  labs(x = "Total Length (m)", y = bquote('Body Volume '(m^3)), color = "Legend") + # title and labels 
  ggtitle("(a)") + # title 
  scale_fill_manual(values = site_colors, name = "Foraging Region", breaks = c("New York Bight", "Gulf of Maine", "Greenland", "Iceland", "Historical")) + # legend 
  scale_color_manual(values = site_colors, name = "Foraging Region", breaks = c("New York Bight", "Gulf of Maine", "Greenland", "Iceland", "Historical")) + # legend 
  theme_bw() # + # theme changed to black and white to get rid of the default ggplot 
#theme(axis.text=element_text(size=15), title=element_text(size=15), legend.text=element_text(size=10), legend.title=element_text(size=11)) # labels formatting 
tl_v_bv_plot # plot the plot

# Fixing Up The Plots -----------------------------------------------------


gg_sma_log_plot + theme(legend.position = "bottom", legend.text=element_text(size=19), 
                        legend.title = element_text(size=19),
                        #legend.title=element_text(size=11, face = "bold"),
                        plot.subtitle=element_text(size=22),
                        axis.text = element_text(size=19),
                        axis.title = element_text(size=22),
                        plot.title=element_text(face = "bold", size = 22)) + guides(color = guide_legend(nrow = 2))

tl_v_bv_plot + theme(legend.position = "none", 
                          legend.text=element_text(size=19), 
                          legend.title = element_text(size=19),
                          #legend.title=element_text(size=11, face = "bold"),
                          plot.subtitle=element_text(size=22),
                          axis.title = element_text(size=22),
                          axis.text = element_text(size=19),
                          plot.title=element_text(face = "bold", size = 22)) + guides(color = guide_legend(nrow = 3))


# 12 M Whale Analysis -----------------------------------------------------

whale_avglength <- 12 # this is the whale length that we're going with
sma_length <- log10(whale_avglength) # whale length in the same log unit as in the SMA analysis 
nyb_elev <- bb$elevation[1] # elevation of nyb from the SMA 
gom_elev <- bb$elevation[2] # elevation of GOM from the SMA
greenland_elev <- bb$elevation[3] # elevation of greenland from the SMA 
iceland_elev <- bb$elevation[4] # elevation of iceland from the SMA 
slope <- bb$slope[1] # common slope from the SMA 

greenland_avlogmass <- sma_length*slope + greenland_elev # finding average log mass for greenland for a 12 m whale
gom_avlogmass <- sma_length*slope + gom_elev # finding average log mass for GOM for a 12 m whale 
nyb_avlogmass <- sma_length*slope + nyb_elev # finding average log mass for nyb for a 12 m whale
iceland_avlogmass <- sma_length*slope + iceland_elev # finding average log mass for iceland for a 12 m whale 

greenland_mass <- 10^greenland_avlogmass
gom_mass <- 10^gom_avlogmass
nyb_mass <- 10^nyb_avlogmass
iceland_mass <-10^iceland_avlogmass

otherregions_mass <- mean(c(greenland_mass, nyb_mass, iceland_mass))
mass_diff <- gom_mass-otherregions_mass

