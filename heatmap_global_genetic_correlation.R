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
         -0.1443, -0.2499, -0.0995, -0.216, 0.0408, 0.0158,
         0.0624, -0.2769, 0.0034, -0.0705, -0.1937, 0.0059,
         -0.1218, 0.0873, 0.1124, 0.0915, 0.1238, 0.0236,
         -0.1959, -0.1173, -0.0601, 0.0725, -0.0878, 0.012,
         -0.3245, -0.1497, -0.2061, -0.1495, -0.1869, 0.0668),
  P_value = c(0.1444, 0.7217, 0.2199, 0.0763, 0.6343, 0.6867,
              0.6121, 0.4076, 0.4931, 0.2463, 0.7852, 0.9005,
              0.7705, 0.3427, 0.9714, 0.53305, 0.0757, 0.94,
              0.0991, 0.2253, 0.0014, 0.0208, 0.0009, 0.4872,
              0.1071, 0.2959, 0.2312, 0.2891, 0.0909, 0.8002,
              0.0056, 0.4406, 0.0003, 0.0077, 0.0004, 0.213)
)

# Remove rows with missing Rg values
data <- data[!is.na(data$Rg), ]

# Add significance asterisks based on P_value
data$Significance <- ifelse(data$P_value < 0.008, "***", 
                            ifelse(data$P_value < 0.01, "**", 
                                   ifelse(data$P_value < 0.05, "*", "")))

# Reshape the data into a matrix for pheatmap
heatmap_data <- reshape2::acast(data, Traits ~ Disease, value.var = "Rg")
significance_data <- acast(data, Traits ~ Disease, value.var = "Significance")

# Reorder the traits and diseases
ordered_traits <- c("eGFR", "Hematuria", "Creatinine", "Microalbumin", "Potassium", "Sodium")
ordered_diseases <- c("AD_sexcombined", "AD_male", "AD_female", "PD_sexcombined", "PD_male", "PD_female")

# Ensure columns are in the correct order
heatmap_data <- heatmap_data[ordered_traits, ordered_diseases]
significance_data <- significance_data[ordered_traits, ordered_diseases]

# Define a color palette with 0 as white, negative values as blue, and positive as red
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)

# Function for displaying numbers with significance
custom_display <- function(values, significance) {
  matrix(sprintf("%.3f%s", values, significance), nrow = nrow(values), dimnames = dimnames(values))
}

# Set a fixed range for the legend
fixed_range <- 0.5

# Create the heatmap
pheatmap(heatmap_data, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         display_numbers = custom_display(heatmap_data, significance_data), 
         color = color_palette, 
         breaks = seq(-fixed_range, fixed_range, length.out = 101), 
         number_color = "black",
         gaps_col = c(3),
         labels_col = c("Sex-combined", "Male", "Female", "Sex-combined", "Male", "Female"),  #Change the labels for different data used
         main = " ",
         fontsize_number = 10)

### Heatmap for PD with and No proxy

# Create the data frame for PD and No Proxy PD
data <- data.frame(
  Disease = c("PD_sexcombined", "PD_sexcombined", "PD_sexcombined", "PD_sexcombined", "PD_sexcombined", "PD_sexcombined",
              "NP_PD_sexcombined", "NP_PD_sexcombined", "NP_PD_sexcombined", "NP_PD_sexcombined", "NP_PD_sexcombined", "NP_PD_sexcombined",
              "PD_male", "PD_male", "PD_male", "PD_male", "PD_male", "PD_male",
              "NP_PD_male", "NP_PD_male", "NP_PD_male", "NP_PD_male", "NP_PD_male", "NP_PD_male",
              "PD_female", "PD_female", "PD_female", "PD_female", "PD_female", "PD_female",
              "NP_PD_female", "NP_PD_female", "NP_PD_female", "NP_PD_female", "NP_PD_female", "NP_PD_female"),
  Traits = c(rep(c("Creatinine", "eGFR", "Hematuria", "Microalbumin", "Potassium", "Sodium"), 6)),
  Rg = c(0.1124, 0.0236, -0.1218, 0.0873, 0.0915, 0.1238,  # PD_sexcombined
         0.0871, 0.0262, -0.1160, 0.1027, 0.1062, 0.0726,  # No_proxy_PD_sexcombined
         -0.0601, 0.012, -0.1959, -0.1173, 0.0725, -0.0878,  # PD_male
         -0.0466, -0.0003, -0.1681, -0.0573, 0.0557, -0.1197,  # No_proxy_PD_male
         -0.2061, 0.0668, -0.3245, -0.1497, -0.1495, -0.1869,  # PD_female
         -0.2125, 0.0955, -0.3343, -0.2015, -0.1214, -0.2013),  # No_proxy_PD_female
  P_value = c(0.0014, 0.4872, 0.0991, 0.2253, 0.0208, 0.0009,  # PD_sexcombined
              0.0115, 0.419, 0.1501, 0.1376, 0.0066, 0.0545,  # No_proxy_PD_sexcombined
              0.2312, 0.8002, 0.1071, 0.2959, 0.2891, 0.0909,  # PD_male
              0.3881, 0.995, 0.1742, 0.6188, 0.4405, 0.0326,  # No_proxy_PD_male
              0.0003, 0.213, 0.0056, 0.4406, 0.0077, 0.0004,  # PD_female
              0.0004, 0.1009, 0.0092, 0.3709, 0.0446, 9.8933e-05)  # No_proxy_PD_female
)

# Remove rows with missing Rg values
data <- data[!is.na(data$Rg), ]

# Add significance asterisks based on P_value
data$Significance <- ifelse(data$P_value < 0.008, "***", 
                            ifelse(data$P_value < 0.01, "**", 
                                   ifelse(data$P_value < 0.05, "*", "")))

# Reshape the data into a matrix for pheatmap
heatmap_data <- acast(data, Traits ~ Disease, value.var = "Rg")
significance_data <- acast(data, Traits ~ Disease, value.var = "Significance")

# Reorder the diseases and traits in the desired order
ordered_traits <- c("eGFR", "Hematuria", "Creatinine", "Microalbumin", "Potassium", "Sodium")
ordered_diseases <- c("PD_sexcombined", "NP_PD_sexcombined", "PD_male", "NP_PD_male", "PD_female", "NP_PD_female")

# Ensure columns are in the correct order
heatmap_data <- heatmap_data[ordered_traits, ordered_diseases]
significance_data <- significance_data[ordered_traits, ordered_diseases]

# Define a color palette with 0 as white, negative values as blue, and positive as red
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)

fixed_range <- 0.50

custom_display <- function(values, significance) {
  # Combine the values and significance annotations
  formatted_data <- matrix(sprintf("%.3f%s", values, significance), 
                           nrow = nrow(values), 
                           dimnames = dimnames(values))
  return(formatted_data)
}

# Create heatmap
pheatmap(heatmap_data, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         display_numbers = custom_display(heatmap_data, significance_data), 
         color = color_palette, 
         breaks = seq(-fixed_range, fixed_range, length.out = 101), 
         fontsize_number = 10,
         number_color = "black",  # Set number color to black
         labels_col = c("Proxy", "No Proxy", "Proxy", "No Proxy", "Proxy", "No Proxy"), # Add the custom labels for columns
         angle_col = 0,
         main = " "
)
