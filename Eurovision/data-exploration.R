# Code to do that data analysis on the Eurovision data

library(ggplot2)
library(ggrepel)
library(foreach)
library(here)

library(rnaturalearth)
library(rnaturalearthdata)
library(sf)

library(dplyr)

library(network)
library(ggraph)

library(GGally)

plotdir <- here("Eurovision", "plots")

# Read in all the csv files with response and covariates
Y <- read.csv(here("Eurovision", "2015data", "votes-ranks.csv"), row.names=1)
Xc.gender <- read.csv(here("Eurovision", "2015data", "covariate-gender.csv"), row.names=1, stringsAsFactors=FALSE)
Xd.lang <- read.csv(here("Eurovision", "2015data", "covariate-language.csv"), row.names=1)
Xd.contig <- read.csv(here("Eurovision", "2015data", "covariate-contig.csv"), row.names=1)
Xc.odds <- read.csv(here("Eurovision", "2015data", "covariate-odds.csv"), row.names=1)
Xc.pop <- read.csv(here("Eurovision", "2015data", "covariate-population.csv"), row.names=1)
Xc.gdp <- read.csv(here("Eurovision", "2015data", "covariate-gdp.csv"), row.names=1)

# Which countries are we restricting the analysis to?
countries.in.final <- row.names(Xc.gender)[Xc.gender[,"Female"] != "#N/A"]

Xd <- array(
  c(#as.matrix(Xd.lang[countries.in.final, countries.in.final]),
    as.matrix(Xd.contig[countries.in.final, countries.in.final])),
  dim=c(length(countries.in.final), length(countries.in.final), 1)
)
dimnames(Xd)[[3]] <- list("Contig")
#Xc <- data.matrix(cbind(
#  Xc.gender[countries.in.final,],
#  Xc.odds[countries.in.final,"Median16"]
#))
#colnames(Xc) <- c(colnames(Xc.gender), "MedianOdds")
Xc <- data.matrix(data.frame(LogMedianOdds=log(Xc.odds[countries.in.final,"Median16"]),
                             LogPopulation=log(Xc.pop[countries.in.final,"Pop2015"]),
                             LogGDP=log(Xc.gdp[countries.in.final,"GDPpc2015"])))


################################################################################
# Make a map of the participating countries

world <- ne_countries(returnclass="sf", scale="medium")

highlighted.countries <- countries.in.final |>
  replace(countries.in.final=="UnitedKingdom", "United Kingdom") |>
  replace(countries.in.final=="Serbia", "Republic of Serbia")

mapplot <- world |>
  mutate(eurovisionfinalist = ifelse(admin %in% highlighted.countries, "Y", "N")) |>
  ggplot() +
  geom_sf(aes(fill = eurovisionfinalist)) +
  scale_fill_manual(values=c(N="#ffffff", Y="#7755ff")) +
  theme_void() +
  theme(legend.position = "none")

################################################################################
# Plot the vote network

ncountries <- length(countries.in.final)

votenet <- network(as.matrix(Y[countries.in.final,countries.in.final]), directed=TRUE)
votenet %v% "nodename" <- countries.in.final
votenet %e% "rank" <- as.matrix(Y[countries.in.final,countries.in.final])
points <- as.matrix(Y[countries.in.final,countries.in.final])
points[points == 10] <- 12
points[points == 9] <- 10
votenet %e% "points" <- points
sender <- matrix(rep(countries.in.final, ncountries), ncol=ncountries)
votenet %e% "voter" <- sender
votenet %e% "song" <- t(sender)

arrowsize <- 8
radius <- 2
labelpad <- 2

voteplot <- ggraph(votenet, layout="manual", 
       x = radius * cos(2*pi*seq_len(ncountries)/ncountries),
       y = radius * sin(2*pi*seq_len(ncountries)/ncountries)) +
  geom_edge_fan(
    aes(edge_width=rank, edge_alpha=rank,
        start_cap=label_rect(voter, padding = margin(labelpad, labelpad, labelpad, labelpad, "mm")), 
        end_cap=label_rect(song, padding = margin(labelpad, labelpad, labelpad, labelpad, "mm"))),
    arrow = arrow(type = "closed", length=unit(arrowsize, "pt"))) +
  geom_node_label(aes(label=nodename)) +
  theme_graph(base_family="serif") +
  scale_edge_width_continuous(range=c(0.1,1.2), trans=scales::exp_trans(1.4)) +
  scale_edge_alpha_continuous(range=c(0.1,1), trans=scales::exp_trans(1.1)) +
  theme(legend.position = "none")


################################################################################
# Make a network of the country contiguity covariate

ncountries <- length(countries.in.final)
countryborders <- as.matrix(Xd.contig[countries.in.final, countries.in.final])
diag(countryborders) <- 0
notislands <- names(which(rowSums(countryborders) > 0))
islands <- setdiff(countries.in.final, notislands)
countryborders <- countryborders[notislands, notislands]
bordernet <- network(countryborders, directed=FALSE)
bordernet %v% "nodename" <- rownames(countryborders)

labelpad <- 2

borderplot <- ggraph(bordernet, layout="igraph", algorithm="dh") +
  geom_edge_link() +
  geom_node_label(aes(label=nodename))+
  theme_graph(base_family="serif")

################################################################################
# Plot covariates

covariateplot <- ggpairs(as.data.frame(Xc), diag=list(continuous="densityDiag")) +
  theme_bw(base_family = "serif")

################################################################################
# Save figures

fig.width <- 9
fig.height <- 6
fig.res <- 150

png(file.path(plotdir, "eurovision-votes-network.png"), width=9, height=6, units="in", res=fig.res)
print(voteplot)
dev.off()

png(file.path(plotdir, "eurovision-country-map.png"), width=10, height=6, units="in", res=fig.res)
print(mapplot)
dev.off()

png(file.path(plotdir, "eurovision-country-borders.png"), width=10, height=5, units="in", res=fig.res)
print(borderplot)
dev.off()

png(file.path(plotdir, "eurovision-covariates.png"), width=6, height=6, units="in", res=fig.res)
print(covariateplot)
dev.off()