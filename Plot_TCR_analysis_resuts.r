library(ggplot2)
library(viridisLite)
library(dplyr)
library(tidyverse)
library(RColorBrewer)

# Load the confidential TCR sequence information from an external R script
source("C:/Users/Desktop/TCR_analysis/TCR_parameters.r") 

#Directory path
directory_path <- "C:/Users/Desktop/TCR_analysis"

all_matrix_beta <- read_csv(file.path(directory_path, "all_matrix_beta_with_markers.csv"))

# Classify TCR beta sequences into two groups: "Studied_TCR" and "Other_TCR" based on a specific peptide sequence
all_matrix_beta$group <- ifelse(all_matrix_beta$cdr3_b_aa_search == studied_beta_TCR, "Studied_TCR", "Other_TCR") 

# Remove rows where the distance is -1 (when TCRs are compared against themselves, the distance is -1)
all_matrix_beta <- subset(all_matrix_beta, dist != -1)

all_matrix_alpha <- read_csv(file.path(directory_path, "all_matrix_alpha_with_markers.csv"))

# Classify TCR alpha sequences into two groups: "Studied_TCR" and "Other_TCR" based on a specific peptide sequence
all_matrix_alpha$group <- ifelse(all_matrix_alpha$cdr3_a_aa_search == studied_alpha_TCR, "Studied_TCR", "Other_TCR") 

# Remove rows where the distance is -1 (when TCRs are compared against themselves, the distance is -1)
all_matrix_alpha <- subset(all_matrix_alpha, dist != -1)

# Plot for TCR beta
plot_beta <- ggplot(all_matrix_beta,aes(x=dist,fill=group))+
            geom_histogram(aes(y=asin(sqrt(after_stat(density)))),
                        alpha=0.5,position='identity',binwidth=3)+
            scale_fill_manual(values = c(rgb(0,0,1,0.5) , rgb(1, 0, 0, 0.5)), 
                            labels = c(studied_beta_TCR, "All CDR3s"), 
                            name = "") +
            ylab("Transformed Relative Frequency of CDR3 Distance Scores") +
            xlab("Distance Score") +
            ggtitle(expression("Distribution of TCR"~beta~"Distance Scores"))+
            theme(plot.title = element_text(size = 25, hjust = 0.5, face = "bold"),
            legend.text = element_text(size = 20, face= "bold"),
            axis.text = element_text(size = 18, face= "bold"),
            axis.title = element_text(size = 20, face= "bold"))

# Save the TCR beta plot as PNG and EPS files
ggsave(file.path(directory_path, "TCR_beta_arcsin_transf.png"), plot = plot_beta, width = 16, height = 12, units = "in")
ggsave(file.path(directory_path, "TCR_beta_arcsin_transf.eps"), plot = plot_beta, width = 16, height = 12, units = "in", device=cairo_ps)

# Plot for TCR alpha
plot_alpha <- ggplot(all_matrix_alpha,aes(x=dist,fill=group))+
            geom_histogram(aes(y=asin(sqrt(..density..))),
                        alpha=0.5,position='identity',binwidth=3)+
            scale_fill_manual(values = c(rgb(0,0,1,0.5) , rgb(1, 0, 0, 0.5)), 
                            labels = c(studied_alpha_TCR, "All CDR3s"),
                            name = "") +
            ylab("Transformed Relative Frequency of CDR3 Distance Scores") +
            xlab("Distance Score") +
            ggtitle(expression("Distribution of TCR"~alpha~"Distance Scores"))+
            theme(plot.title = element_text(size = 25, hjust = 0.5, face = "bold"),
            legend.text = element_text(size = 20, face= "bold"),
            axis.text = element_text(size = 18, face= "bold"),
            axis.title = element_text(size = 20, face= "bold"))

# Save the TCR beta plot as PNG and EPS files
ggsave(file.path(directory_path, "TCR_alpha_arcsin_transf.png"), plot = plot_alpha, width = 16, height = 12, units = "in")
ggsave(file.path(directory_path, "TCR_alpha_arcsin_transf.eps"), plot = plot_alpha, width = 16, height = 12, units = "in", device=cairo_ps)




