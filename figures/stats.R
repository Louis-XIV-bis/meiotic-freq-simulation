#!/usr/bin/env Rscript

## Author : Louis OLLIVIER ~ louis.ollivier@universite-paris-saclay.fr
## Université Paris-Saclay
## Lab : LISN ~ UMR9015 ~ BIOINFO team 

###### Package initialization  ----------------------------------------

if (!require('readr', quietly = T)) install.packages('readr');
if (!require('tibble', quietly = T)) install.packages('tibble');
if (!require('dplyr', quietly = T)) install.packages('dplyr');
if (!require('dunn.test', quietly = T)) install.packages('dunn.test');

library(readr)
library(tibble)
library(dplyr)
library(dunn.test)


########################################################################
# Function that upload and format the data

rm(list=ls())

# Load full dataset
full_dataset = read_csv('../data/t_pichr2_full_merged.csv') %>%
  mutate(alpha = as.factor(1 / GR)) %>% 
  mutate(t = t - 2000)
full_dataset

########################################################################
# Different functions used 

# Function to test if the distribution of a variable is Gaussian for each combination of `group_var` and `alpha`
check_normality <- function(data, group_var, test_var) {
  # Perform Shapiro-Wilk test for normality for each combination of `group_var` and `alpha`
  normality_results <- data %>%
    group_by(!!sym(group_var), alpha) %>%
    summarize(
      shapiro_p_value = shapiro.test(!!sym(test_var))$p.value,  # Shapiro-Wilk test for the specified variable
      is_gaussian = ifelse(shapiro.test(!!sym(test_var))$p.value > 0.05, "Yes", "No")  # Normality check
    )
  
  # Return the results
  return(normality_results)
}

perform_tests <- function(data, group_var, value_var, group_col) {
  # Check if the grouping column exists
  if (!(group_var %in% names(data))) stop("Group variable not found in data.")
  if (!(value_var %in% names(data))) stop("Value variable not found in data.")
  
  # Get unique values of the grouping variable
  unique_values <- unique(data[[group_col]])
  
  # Loop through each unique value and perform Kruskal-Wallis test
  for (value in unique_values) {
    # Filter data for the current grouping value
    filtered_data <- data %>% filter(!!sym(group_col) == value)
    
    # Perform Kruskal-Wallis test
    test_result <- kruskal.test(as.formula(paste(value_var, "~", group_var)), data = filtered_data)
    
    # Print the results
    cat("Results for", group_col, "=", value, ":\n")
    print(test_result)
    
    # Check if p-value is less than 0.05
    if (test_result$p.value < 0.05) {
      cat("There's a significant difference between some groups\n")
      
      # Perform Dunn's test for pairwise comparisons
      cat("Dunn's test results:\n")
      dunn.test(filtered_data[[value_var]], filtered_data[[group_var]])
      
    } else {
      cat("No significant difference between groups\n")
    }
    
    cat("\n")
  }
}
########################################################################
# Statistical tests 

###### Fig 1 ###########################################################

data_neutral = read_csv('../data/pi_neutral.csv') %>%
  mutate(alpha = as.factor(1 / GR)) 

pi_neutral_alpha_1 <- data_neutral %>% filter(alpha == 1) %>% pull(mean)
shapiro.test(pi_neutral_alpha_1) #  p-value = 0.759 ==> gaussian

pi_neutral_alpha_0.01 <- data_neutral %>% filter(alpha == 0.01) %>% pull(mean)
shapiro.test(pi_neutral_alpha_0.01) #  p-value = 0.768 ==> gaussian

t.test(pi_neutral_alpha_1, pi_neutral_alpha_0.01, alternative = "less")
# p-value = 0.0001139 ==> pi alpha 1 < pi alpha 0.01

4.007228e-5 / 4.115257e-5
###### Fig 2 ###########################################################


### pi chr 2 alpha = 0.01 vs alpha == 1 ###
data_fig2 = full_dataset %>% 
  filter(s == 0.05 & rho == '5e-08' & rho_scaled != 'rho_fixe', h == 0.5 & (alpha == 1 | alpha == 0.01)) 
data_fig2
check_normality(data_fig2, "window", "pi")
# window alpha shapiro_p_value is_gaussian
# <chr>  <fct>           <dbl> <chr>      
# 1 chr2   0.01           0.298  Yes        
# 2 chr2   1              0.0265 No       

# non parametric test to be exact because only 1 gaussian 
pi_alpha_1 <- data_fig2 %>% filter(alpha == 1) %>% pull(pi)
mean(pi_alpha_1)
pi_alpha_0.01 <- data_fig2 %>% filter(alpha == 0.01) %>% pull(pi)
mean(pi_alpha_0.01)
wilcox.test(pi_alpha_1, pi_alpha_0.01, alternative = "less")
# p-value < 2.2e-16 ==> pi_alpha=1 != pi_alpha=0.01

### T chr 2 alpha = 0.01 vs alpha == 1 ###

check_normality(data_fig2, "window", "t")
# window alpha shapiro_p_value is_gaussian
# <chr>  <fct>           <dbl> <chr>      
# 1 chr2   0.01         3.93e-11 No         
# 2 chr2   1            1.28e- 7 No         

# non parametric test to be exact 
t_alpha_1 <- data_fig2 %>% filter(alpha == 1) %>% pull(t)
t_alpha_0.01 <- data_fig2 %>% filter(alpha == 0.01) %>% pull(t)
wilcox.test(t_alpha_1, t_alpha_0.01)
# p-value = 0.02631 > t_alpha=1 != t_alpha=0.01 ==> difference in mean = 20 generation

rm(pi_alpha_1, pi_alpha_0.01, t_alpha_1, t_alpha_0.01, data_fig2)

###### Fig 3 ######################

### Compare pi / pi0 for a given s for each alpha and the opposite ###
data_s = full_dataset %>%
  filter(rho == '5e-08' & rho_scaled != 'rho_fixe', h == 0.5) %>%
  group_by(s) %>%
  mutate(pi_normalized = pi / mean(pi[alpha == 1], na.rm = TRUE)) %>%
  ungroup()
data_s

check_normality(data_s, "s", "pi_normalized")

# s alpha shapiro_p_value is_gaussian
# <dbl> <fct>           <dbl> <chr>      
# 1  0.02 0.01          0.205   Yes        
# 2  0.02 0.02          0.219   Yes        
# 3  0.02 0.1           0.00164 No         
# 4  0.02 1             0.00964 No         
# 5  0.05 0.01          0.298   Yes        
# 6  0.05 0.02          0.946   Yes        
# 7  0.05 0.1           0.372   Yes        
# 8  0.05 1             0.0265  No         
# 9  0.1  0.01          0.524   Yes        
# 10  0.1  0.02          0.0632  Yes        
# 11  0.1  0.1           0.202   Yes        
# 12  0.1  1             0.455   Yes        

# Will do kruskal walis and dunn.test because not all of them are gaussian 

# Compare pi for each alpha for a given s  
perform_tests(data_s, group_var = "alpha", value_var = "pi_normalized", group_col = "s") 

# Compare pi for each s for a given alpha 
perform_tests(data_s, group_var = "s", value_var = "pi_normalized", group_col = "alpha")

### Compare t for a given s for each alpha and the opposite ###
check_normality(data_s, "s", "t")

# s alpha shapiro_p_value is_gaussian
# <dbl> <fct>           <dbl> <chr>      
# 1  0.02 0.01         6.34e- 8 No         
# 2  0.02 0.02         1.57e-13 No         
# 3  0.02 0.1          3.59e- 9 No         
# 4  0.02 1            8.66e-12 No         
# 5  0.05 0.01         3.93e-11 No         
# 6  0.05 0.02         3.42e- 8 No         
# 7  0.05 0.1          1.03e- 8 No         
# 8  0.05 1            1.28e- 7 No         
# 9  0.1  0.01         6.89e- 3 No         
# 10  0.1  0.02         1.36e-10 No         
# 11  0.1  0.1          6.52e- 7 No         
# 12  0.1  1            7.16e-10 No    

# Will do kruskal walis and dunn.test because not all of them are gaussian 
# Compare t for each alpha for a given s  
perform_tests(data_s, group_var = "alpha", value_var = "t", group_col = "s") 

# Compare t for each s for a given alpha 
perform_tests(data_s, group_var = "s", value_var = "t", group_col = "alpha")

###### Fig 4 ######################

### Compare pi for a given rho for each alpha and the opposite ###
data_rho = full_dataset %>%
  filter(s ==  0.05 & rho_scaled != 'rho_fixe', h == 0.5) %>%
  group_by(rho) %>%
  mutate(pi_normalized = pi / mean(pi[alpha == 1], na.rm = TRUE)) %>%
  ungroup()

check_normality(data_rho, "rho", "pi_normalized")

# rho alpha shapiro_p_value is_gaussian
# <dbl> <fct>           <dbl> <chr>      
# 1 0.00000001 0.01         0.106    Yes        
# 2 0.00000001 0.02         0.812    Yes        
# 3 0.00000001 0.1          0.644    Yes        
# 4 0.00000001 1            0.196    Yes        
# 5 0.00000005 0.01         0.298    Yes        
# 6 0.00000005 0.02         0.946    Yes        
# 7 0.00000005 0.1          0.372    Yes        
# 8 0.00000005 1            0.0265   No         
# 9 0.0000001  0.01         0.000159 No         
# 10 0.0000001  0.02         0.843    Yes        
# 11 0.0000001  0.1          0.189    Yes        
# 12 0.0000001  1            0.0711   Yes   

# Will do kruskal walis and dunn.test because not all of them are gaussian 

# Compare pi for each alpha for a given rho  
perform_tests(data_rho, group_var = "alpha", value_var = "pi_normalized", group_col = "rho") 

# Compare pi for each rho for a given alpha 
perform_tests(data_rho, group_var = "rho", value_var = "pi_normalized", group_col = "alpha")

### Compare t for a given rho for each alpha and the opposite ###
check_normality(data_rho, "rho", "t")

# rho alpha shapiro_p_value is_gaussian
# <dbl> <fct>           <dbl> <chr>      
#   1 0.00000001 0.01         1.68e- 7 No         
# 2 0.00000001 0.02         9.76e- 4 No         
# 3 0.00000001 0.1          4.26e-12 No         
# 4 0.00000001 1            1.63e-10 No         
# 5 0.00000005 0.01         3.93e-11 No         
# 6 0.00000005 0.02         3.42e- 8 No         
# 7 0.00000005 0.1          1.03e- 8 No         
# 8 0.00000005 1            1.28e- 7 No         
# 9 0.0000001  0.01         5.12e- 8 No         
# 10 0.0000001  0.02         5.27e- 9 No         
# 11 0.0000001  0.1          3.58e- 9 No         
# 12 0.0000001  1            4.57e- 6 No 

# Will do kruskal walis and dunn.test because not all of them are gaussian 
# Compare t for each alpha for a given rho  
perform_tests(data_rho, group_var = "alpha", value_var = "t", group_col = "rho") 

# Compare t for each rho for a given alpha 
perform_tests(data_rho, group_var = "rho", value_var = "t", group_col = "alpha")

###### Fig 5 ######################
### Compare pi for a given rho_var for each alpha and the opposite ###
data_rho_m = full_dataset %>%
  filter((h == 0.5 & s == 0.05) & (rho == '5e-08' | rho_scaled == 'rho_fixe')) %>% 
  mutate(rho_scaled = case_when(rho_scaled == 'rho_fixe' ~ 'rho_fixe', TRUE ~ 'rho_var')) %>%
  group_by(rho_scaled) %>%
  mutate(pi_normalized = pi / mean(pi[alpha == 1], na.rm = TRUE)) %>%
  ungroup()
data_rho_m

check_normality(data_rho_m, "rho_scaled", "pi_normalized")

# rho_scaled alpha shapiro_p_value is_gaussian
# <chr>      <fct>           <dbl> <chr>      
# 1 rho_fixe   0.01       0.00000103 No         
# 2 rho_fixe   0.02       0.00152    No         
# 3 rho_fixe   0.1        0.574      Yes        
# 4 rho_fixe   1          0.213      Yes        
# 5 rho_var    0.01       0.298      Yes        
# 6 rho_var    0.02       0.946      Yes        
# 7 rho_var    0.1        0.372      Yes        
# 8 rho_var    1          0.0265     No         

# Will do kruskal walis and dunn.test because not all of them are gaussian 

# Compare pi for each alpha for a given rho_scaled  
perform_tests(data_rho_m, group_var = "alpha", value_var = "pi_normalized", group_col = "rho_scaled") 

# Compare pi for each rho_scaled for a given alpha 
perform_tests(data_rho_m, group_var = "rho_scaled", value_var = "pi_normalized", group_col = "alpha")

### Compare t for a given rho for each alpha and the opposite ###
check_normality(data_rho_m, "rho_scaled", "t")

# rho_scaled alpha shapiro_p_value is_gaussian
# <chr>      <fct>           <dbl> <chr>      
# 1 rho_fixe   0.01         3.31e- 8 No         
# 2 rho_fixe   0.02         1.74e- 9 No         
# 3 rho_fixe   0.1          1.17e- 7 No         
# 4 rho_fixe   1            1.83e- 7 No         
# 5 rho_var    0.01         3.93e-11 No         
# 6 rho_var    0.02         3.42e- 8 No         
# 7 rho_var    0.1          1.03e- 8 No         
# 8 rho_var    1            1.28e- 7 No 

# Will do kruskal walis and dunn.test because not all of them are gaussian 
# Compare t for each alpha for a given rho_scaled  
perform_tests(data_rho_m, group_var = "alpha", value_var = "t", group_col = "rho_scaled") 

# Compare t for each rho_scaled for a given alpha 
perform_tests(data_rho_m, group_var = "rho_scaled", value_var = "t", group_col = "alpha")


###### Fig 6 ######################

### Compare pi for a given rho for each alpha and the opposite ###
data_h = full_dataset %>%
  filter(s ==  0.05 & rho_scaled != 'rho_fixe' & rho == '5e-08') %>%
  group_by(h) %>%
  mutate(pi_normalized = pi / mean(pi[alpha == 1], na.rm = TRUE)) %>%
  ungroup()
data_h
check_normality(data_h, "h", "pi_normalized")

# h alpha shapiro_p_value is_gaussian
# <dbl> <fct>           <dbl> <chr>      
# 1   0.2 0.01      0.0152      No         
# 2   0.2 0.02      0.000000187 No         
# 3   0.2 0.1       0.241       Yes        
# 4   0.2 1         0.899       Yes        
# 5   0.5 0.01      0.298       Yes        
# 6   0.5 0.02      0.946       Yes        
# 7   0.5 0.1       0.372       Yes        
# 8   0.5 1         0.0265      No         
# 9   0.8 0.01      0.421       Yes        
# 10   0.8 0.02      0.978       Yes        
# 11   0.8 0.1       0.0215      No         
# 12   0.8 1         0.772       Yes   

# Will do kruskal walis and dunn.test because not all of them are gaussian 

# Compare pi for each alpha for a given h  
perform_tests(data_h, group_var = "alpha", value_var = "pi_normalized", group_col = "h") 

# Compare pi for each h for a given alpha 
perform_tests(data_h, group_var = "h", value_var = "pi_normalized", group_col = "alpha")

### Compare t for a given h for each alpha and the opposite ###
check_normality(data_h, "h", "t")

# h alpha shapiro_p_value is_gaussian
# <dbl> <fct>           <dbl> <chr>      
#   1   0.2 0.01         1.23e-14 No         
# 2   0.2 0.02         2.68e-11 No         
# 3   0.2 0.1          2.65e-13 No         
# 4   0.2 1            2.33e-12 No         
# 5   0.5 0.01         3.93e-11 No         
# 6   0.5 0.02         3.42e- 8 No         
# 7   0.5 0.1          1.03e- 8 No         
# 8   0.5 1            1.28e- 7 No         
# 9   0.8 0.01         7.08e-16 No         
# 10   0.8 0.02         3.94e-14 No         
# 11   0.8 0.1          3.60e-14 No         
# 12   0.8 1            6.38e-13 No     

# Will do kruskal walis and dunn.test because not all of them are gaussian 
# Compare t for each alpha for a given h  
perform_tests(data_h, group_var = "alpha", value_var = "t", group_col = "h") 

# Compare t for each h for a given alpha 
perform_tests(data_h, group_var = "h", value_var = "t", group_col = "alpha")


