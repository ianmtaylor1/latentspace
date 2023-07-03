library(here)
outdir <- here("manuscript", "plots")


# Figure 1 - illustration of network and spatial dependence
library(usmap)
library(network)
library(ggplot2)
library(ggraph)
library(latex2exp)
library(gridExtra)

# Make a map of the states
statesmapplot <- plot_usmap("states", include=c(.mountain, "CA", "WA", "OR"), labels=TRUE)

# Make a network of state border relationships
statesnames <- c("WA", "OR", "CA", "ID", "NV", "MT", "WY", "UT", "CO", "AZ", "NM")
statesmat <- matrix(rep(0, length(statesnames)^2), nrow=length(statesnames))
colnames(statesmat) <- rownames(statesmat) <- statesnames
statesmat["WA", c("OR", "ID")] <- 1
statesmat["OR", c("WA", "ID", "NV", "CA")] <- 1
statesmat["CA", c("OR", "NV", "AZ")] <- 1
statesmat["ID", c("WA", "OR", "NV", "UT", "WY", "MT")] <- 1
statesmat["NV", c("CA", "AZ", "UT", "ID", "OR")] <- 1
statesmat["MT", c("ID", "WY")] <- 1
statesmat["WY", c("MT", "ID", "UT", "CO")] <- 1
statesmat["UT", c("ID", "WY", "CO", "AZ", "NV")] <- 1
statesmat["CO", c("WY", "UT", "NM")] <- 1
statesmat["AZ", c("CA", "NV", "UT", "NM")] <- 1
statesmat["NM", c("CO", "AZ")] <- 1

statesnet <- network(statesmat, directed = FALSE)
#statesnet %v% "stateobsname" <- sapply(statesnames, function(x) { TeX(paste0("y_{", x, "}")) } )
statesnet %v% "stateobsname" <- sapply(statesnames, function(x) paste0("y[", x, "]"))

statesnetplot <- ggraph(statesnet) +
  geom_edge_link() +
  geom_node_label(aes(label=stateobsname), parse=TRUE) +
  theme_graph()
statesnetplot

# Make a simple directed complete graph on 4 nodes
examplemat <- matrix(1, nrow=4, ncol=4)
colnames(examplemat) <- rownames(examplemat) <- c(1, 2, 3, 4)
examplenet <- network(examplemat, directed=TRUE, loops = TRUE)

nodenames <- c(1,2,3,4)
examplenet %v% "nodename" <- nodenames
for (i in nodenames) {
  for (j in nodenames) {
    examplenet[i,j, names.eval="obsname"] <- paste0("y[", i, "][", j, "]")
  }
}

arrowsize <- 6
edgegap <- 12

examplenetplot <- ggraph(examplenet, layout="manual", x=c(1, -1, -1, 1), y=c(1, 1, -1, -1)) +
  geom_edge_fan(aes(label=obsname), 
                label_parse = TRUE, angle_calc="along", label_dodge=unit(-8, "pt"),
                arrow = arrow(type = "closed", length=unit(arrowsize, "pt")), 
                end_cap=circle(edgegap, "pt"), start_cap=circle(edgegap, "pt")) +
  geom_edge_loop(aes(label=obsname, direction=90 * from - 45),
                 label_parse = TRUE, angle_calc="along", label_dodge=unit(-8, "pt"),
                 arrow = arrow(type = "closed", length=unit(arrowsize, "pt")), 
                 end_cap=circle(edgegap, "pt"), start_cap=circle(edgegap, "pt")) +
  geom_node_label(aes(label=nodename)) +
  theme_graph()
examplenetplot

# Depict the dependence relationships between edges of that graph
pairsep=","

edgenames <- c(outer(nodenames, nodenames, function(x,y) {paste(x, y, sep=pairsep)}))
edgerelmat <- matrix(0, nrow=length(edgenames), ncol=length(edgenames))
rownames(edgerelmat) <- colnames(edgerelmat) <- edgenames
edgerelnet <- network(edgerelmat, directed = FALSE, loops = FALSE, multiple = FALSE)

# Populate sender and receiver relationships
for (i in nodenames) {
  for (j in nodenames) {
    for (k in nodenames) {
      edge1 <- paste(i,j, sep=pairsep)
      # Same sender
      if (k != j) {
        edge2 <- paste(i, k, sep=pairsep)
        edgerelnet[edge1, edge2, names.eval="common", add.edges=TRUE] <- "sender"
      }
      # Same receiver
      if (k != i) {
        edge2 <- paste(k, j, sep=pairsep)
        edgerelnet[edge1, edge2, names.eval="common", add.edges=TRUE] <- "receiver"
      }
    }
  }
}
# Apply labels
edgerelnet %v% "obsname" <- sapply(edgerelnet %v% "vertex.names", function(x) paste0("y[", strsplit(x, pairsep)[[1]][1], "][", strsplit(x, pairsep)[[1]][2], "]"))
edgerelnet %v% "sender" <- sapply(edgerelnet %v% "vertex.names", function(x) as.numeric(strsplit(x, pairsep)[[1]][1]))
edgerelnet %v% "receiver" <- sapply(edgerelnet %v% "vertex.names", function(x) as.numeric(strsplit(x, pairsep)[[1]][2]))

senderoffset <- 0
receiveroffset <- pi/12
sendermag <- 1
receivermag <- 3

edgerelplot <- ggraph(edgerelnet, layout="manual", 
                      x = sendermag * cos(senderoffset + sender * pi / 2) + receivermag * cos(receiveroffset + receiver * pi / 2), 
                      y = sendermag * sin(senderoffset + sender * pi / 2) + receivermag * sin(receiveroffset + receiver * pi / 2)) +
  geom_edge_link(aes(color=common)) +
  geom_node_label(aes(label=obsname), parse=TRUE) +
  theme_graph() +
  scale_color_brewer(palette="Dark2") +
  theme(legend.position="none") 
edgerelplot

png(file.path(outdir, "fig1-example-dependene.png"), width=10, height=8, units="in", res = 150)
grid.arrange(statesmapplot, statesnetplot, examplenetplot, edgerelplot,
             nrow=2)
dev.off()