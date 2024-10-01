# Install required packages if not already installed
# install.packages("pheatmap")

# Load the pheatmap package
library(pheatmap)
library(reshape2)
library(grid)  

# Create the dataframe for AD/PD heatmap
data <- data.frame(
  Disease = c("AD_sexcombined", "AD_sexcombined", "AD_sexcombined", "AD_sexcombined", "AD_sexcombined", "AD_sexcombined",
              "AD_male", "AD_male", "AD_male", "AD_male", "AD_male", "AD_male",
              "AD_female", "AD_female", "AD_female", "AD_female", "AD_female", "AD_female",
              "PD_sexcombined", "PD_sexcombined", "PD_sexcombined", "PD_sexcombined", "PD_sexcombined", "PD_sexcombined",
              "PD_male", "PD_male", "PD_male", "PD_male", "PD_male", "PD_male",
              "PD_female", "PD_female", "PD_female", "PD_female", "PD_female", "PD_female"),
  Traits = c("Hematuria", "Microalbumin", "Creatinine", "Potassium", "Sodium", "eGFR",
             "Hematuria", "Microalbumin", "Creatinine", "Potassium", "Sodium", "eGFR",
             "Hematuria", "Microalbumin", "Creatinine", "Potassium", "Sodium", "eGFR",
             "Hematuria", "Microalbumin", "Creatinine", "Potassium", "Sodium", "eGFR",
             "Hematuria", "Microalbumin", "Creatinine", "Potassium", "Sodium", "eGFR",
             "Hematuria", "Microalbumin", "Creatinine", "Potassium", "Sodium", "eGFR"),
  Rg = c(0.1702, -0.0401, 0.0651, 0.1156, 0.0245, 0.0222,
         -0.1443, -0.2499, -0.0995, -0.0705, 0.0408, 0.0158,
         0.0624, -0.2769, 0.0034, -0.0705, -0.1937, 0.0059,
         -0.1218, 0.0873, 0.1124, 0.0915, 0.1238, 0.0236,
         -0.1959, -0.1173, -0.0601, 0.0725, -0.0878, 0.012,
         -0.3245, -0.1497, -0.2061, -0.1495, -0.1869, 0.0668),
  P_value = c(0.1444, 0.7217, 0.2199, 0.0763, 0.6343, 0.6867,
              0.6121, 0.4076, 0.4931, 0.5305, 0.7852, 0.9005,
              0.7705, 0.3427, 0.9714, 0.53305, 0.0757, 0.94,
              0.0991, 0.2253, 0.0014, 0.0208, 0.0009, 0.4872,
              0.1071, 0.2959, 0.2312, 0.2891, 0.0909, 0.8002,
              0.0056, 0.4406, 0.0003, 0.0077, 0.0004, 0.213)
)

# Remove rows with missing Rg values
data <- data[!is.na(data$Rg), ]

# Add significance asterisks based on P_value
data$Significance <- ifelse(data$P_value < 0.001, "***", 
                            ifelse(data$P_value < 0.01, "**", 
                                   ifelse(data$P_value < 0.05, "*", "")))

# Reshape the data into a matrix for pheatmap
heatmap_data <- reshape2::acast(data, Traits ~ Disease, value.var = "Rg")
significance_data <- acast(data, Traits ~ Disease, value.var = "Significance")

# Reorder the columns (AD/PD sex-combined first, then AD/PD males, then AD/PD females)
ordered_cols <- c("AD_sexcombined", "PD_sexcombined", 
                  "AD_male", "PD_male", 
                  "AD_female", "PD_female")

# Ensure columns are in the correct order
heatmap_data <- heatmap_data[, ordered_cols]
significance_data <- significance_data[, ordered_cols]

# Define a color palette with 0 as white, negative values as blue, and positive as red
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)

# Create the heatmap
max_abs_value <- max(abs(heatmap_data), na.rm = TRUE)
pheatmap(heatmap_data, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         display_numbers = custom_display(heatmap_data, significance_data), 
         color = color_palette, 
         breaks = seq(-max_abs_value, max_abs_value, length.out = 101), 
         main = "Genetic Correlations (Rg) between Kidney Traits and Neurodegenerative Diseases",
         fontsize_number = 10)

# Create the data frame for PD and No Proxy PD
data <- data.frame(
  Disease = c("PD_sexcombined", "PD_sexcombined", "PD_sexcombined", "PD_sexcombined", "PD_sexcombined", "PD_sexcombined",
              "PD_male", "PD_male", "PD_male", "PD_male", "PD_male", "PD_male",
              "PD_female", "PD_female", "PD_female", "PD_female", "PD_female", "PD_female",
              "No_proxy_PD_male", "No_proxy_PD_male", "No_proxy_PD_male", "No_proxy_PD_male", "No_proxy_PD_male", "No_proxy_PD_male",
              "No_proxy_PD_female", "No_proxy_PD_female", "No_proxy_PD_female", "No_proxy_PD_female", "No_proxy_PD_female", "No_proxy_PD_female"),
  Traits = c("Hematuria", "Microalbumin", "Creatinine", "Potassium", "Sodium", "eGFR",
             "Hematuria", "Microalbumin", "Creatinine", "Potassium", "Sodium", "eGFR",
             "Hematuria", "Microalbumin", "Creatinine", "Potassium", "Sodium", "eGFR",
             "Hematuria", "Microalbumin", "Creatinine", "Potassium", "Sodium", "eGFR",
             "Hematuria", "Microalbumin", "Creatinine", "Potassium", "Sodium", "eGFR"),
  Rg = c(-0.1218, 0.0873, 0.1124, 0.0915, 0.1238, 0.0236,
         -0.1959, -0.1173, -0.0601, 0.0725, -0.0878, 0.012,
         -0.3245, -0.1497, -0.2061, -0.1495, -0.1869, 0.0668,
         -0.1681, -0.0573, -0.0466, 0.0557, -0.1197, -0.0003,
         -0.3343, -0.2015, -0.2125, -0.1214, -0.2013, 0.0955),
  P_value = c(0.0991, 0.2253, 0.0014, 0.0208, 0.0009, 0.4872,
              0.1071, 0.2959, 0.2312, 0.2891, 0.0909, 0.8002,
              0.0056, 0.4406, 0.0003, 0.0077, 0.0004, 0.213,
              0.1742, 0.6188, 0.3881, 0.4405, -2.1369, 0.995,
              0.0092, 0.3709, 0.0004, 0.0446, 9.8933e-05, 0.1009)
)

# Remove rows with missing Rg values
data <- data[!is.na(data$Rg), ]

# Add significance asterisks based on P_value
data$Significance <- ifelse(data$P_value < 0.001, "***", 
                            ifelse(data$P_value < 0.01, "**", 
                                   ifelse(data$P_value < 0.05, "*", "")))

# Reshape the data into a matrix for pheatmap
heatmap_data <- acast(data, Traits ~ Disease, value.var = "Rg")
significance_data <- acast(data, Traits ~ Disease, value.var = "Significance")

# Reorder the columns (PD first, No Proxie PD next)
ordered_cols <- c("PD_sexcombined", "PD_male", "PD_female", 
                  "No_proxy_PD_male", "No_proxy_PD_female")

# Ensure columns are in the correct order
heatmap_data <- heatmap_data[, ordered_cols]
significance_data <- significance_data[, ordered_cols]

# Define a color palette with 0 as white, negative values as blue, and positive as red
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)

# Create heatmap
max_abs_value <- max(abs(heatmap_data), na.rm = TRUE)
pheatmap(heatmap_data, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         display_numbers = custom_display(heatmap_data, significance_data), 
         color = color_palette, 
         breaks = seq(-max_abs_value, max_abs_value, length.out = 101), 
         main = "Genetic Correlations (Rg) between Kidney Traits and PD/No Proxy PD",
         fontsize_number = 10)
