#Set working directory
setwd("/Users/zhaolei/Desktop/psmc_restart_1_1")

#Load required packages
library(ggplot2)
library(dplyr)
library(mgcv)
library(patchwork)
library(scales)

samples <- "sample"
#
estimate_data <- do.call(rbind, lapply(samples, function(sample) {
  file_path <- paste0("plot/", sample, "_plot.0.txt")
  data <- read.table(file_path, header = FALSE, sep = "\t")[, 1:2]
  colnames(data) <- c("t", "Ne")
  data$Ne <- data$Ne
  data$Sample <- sample
  data$Type <- "Estimate"
  return(data)
}))

bootstrap_files <- paste0("bootstrap/bootstrap_plot.", 0:99, ".txt")
bootstrap_data <- do.call(rbind, lapply(bootstrap_files, function(file_path) {
  
  #Extract sample ID from file names
  i <- gsub(".*bootstrap_plot\\.(\\d+)\\.txt", "\\1", file_path)
  
  data <- read.table(file_path, header = FALSE, sep = "\t")[, 1:2]
  colnames(data) <- c("t", "Ne")
  data$Ne <- data$Ne
  data$Sample <- paste0("Bootstrap_", i) 
  data$Type <- "Bootstrap"  
  return(data)
}))

combined_data <- rbind(estimate_data, bootstrap_data)

#===================== >= 1000 Retain boundary points =====================
filtered_data <- combined_data %>%
  filter(t >= 10000, is.finite(t), is.finite(Ne))
#===================================================================

#Sort by Sample and time (t), and use lead(t) to define step transitions
filtered_data <- filtered_data %>%
  arrange(Sample, t) %>%
  group_by(Sample) %>%
  mutate(
    StepX = lead(t)
  ) %>%
  # ===================== Assign StepX to the last row to avoid data loss =====================
mutate(StepX = ifelse(is.na(StepX), t * 2, StepX)) %>%
  # ===========================================================================
ungroup()

# Define x-axis limits
x_min <- min(filtered_data$t)
x_max <- max(filtered_data$StepX)

#  Define major tick marks and labels
major_breaks_x <- 10^(floor(log10(x_min)):ceiling(log10(x_max)))
labels_x <- parse(text = paste0("10^", floor(log10(x_min)):ceiling(log10(x_max))))

# Define minor tick marks
minor_breaks_x <- unlist(lapply(floor(log10(x_min)):ceiling(log10(x_max)), function(i) seq(1, 9, 1) * 10^i))
minor_breaks_x <- minor_breaks_x[minor_breaks_x >= x_min & minor_breaks_x <= x_max]

peak_left <- 11161
peak_right <- 17953
peak_ne <- 1.46
# Draw main plot
plot1 <- ggplot(filtered_data, aes(color = Type)) +
  #geom_vline(xintercept = c(peak_left, peak_right), 
  #           linetype = "dashed", color = "blue", size = 0.5) +
  # Add rectangle highlighting peak region
  #annotate("rect", 
           #xmin = peak_left, xmax = peak_right,
           #ymin = peak_ne * 0.9, ymax = peak_ne * 1.1,
           #fill = "blue", alpha = 0.3) + 
  # Add rectangle annotation
  #annotate("rect", fill = "#377EB8", xmin = 19000, xmax = 26500, ymin = 1, ymax = Inf, alpha = 0.2) +
  # annotate("rect", fill = "grey",    xmin = 31998.73, xmax = 75682.07, ymin = 1, ymax = Inf, alpha = 0.2) +
  # Add LOESS smoothing curve (using raw t values)
  #geom_smooth(aes(x = t, y = Ne, group = NULL), method = "loess", se = TRUE,
              #color = "black", size = 1, na.rm = TRUE) +
  # Plot step function using StepX (Ne changes only when reaching the next t)
  geom_step(aes(x = t, y = Ne, group = Sample), size = 1, alpha = 1, na.rm = TRUE) +
  theme_minimal() +
  labs(y = expression("Effective size ("*10^4*")")) +
  scale_x_log10(
    breaks       = major_breaks_x,
    minor_breaks = minor_breaks_x,
    labels       = labels_x
  ) +
  scale_y_log10(
    breaks       = c(1, 10),
    labels       = c(1, 10),
    minor_breaks = unlist(lapply(0:2, function(i) seq(1, 9, 1) * 10^i)),
    limits       = c(0.8, 15),
    expand       = c(0, 0)
  ) +
  annotation_logticks(
    sides  = "bl",
    scaled = TRUE,
    short  = unit(0.2, "cm"),
    mid    = unit(0.4, "cm"),
    long   = unit(0.6, "cm"),
    color  = "black"
  ) +
  theme(
    plot.title        = element_text(size = 24, face = "bold", hjust = 0.5),
    axis.text.x       = element_text(size = 20, face = "bold", colour = "black"),
    axis.title.x      = element_text(size = 22, face = "bold", colour = "black"),
    #axis.ticks.x      = element_blank(),
    axis.text.y       = element_text(size = 20, face = "bold", colour = "black"),
    axis.title.y      = element_text(size = 22, face = "bold", colour = "black"),
    panel.border      = element_rect(colour = "black", fill = NA, size = 1),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    axis.ticks.length = unit(0.3, "cm"),
    axis.ticks        = element_line(size = 1, colour = "black"),
    legend.position   = c(0, 1),
    legend.justification = c("left", "top"),
    legend.text       = element_text(size = 18),
    legend.key.size   = unit(2, "lines"),
    legend.title      = element_blank()
  ) +
  guides(color = guide_legend(ncol = 3)) +
  scale_color_manual(
    values = c(
      "Estimate"  = "#E83A29", "Bootstrap" = "#F095B3"
    ),
    breaks = c(
      "Estimate", "Bootstrap"
    )
  ) +
  # ===================== Truncate the left side of the x-axis to 10^3  =====================
coord_cartesian(xlim = c(10000, max(filtered_data$StepX)))
print(plot1)

diff_ne <- diff(filtered_data[filtered_data$Type == "Estimate",]$Ne)

# Provide detailed information on continuous segments
zero_runs <- rle(abs(diff_ne) <= 1e-6)
zero_positions <- which(zero_runs$values == TRUE)
for (pos in zero_positions) {
  start <- ifelse(pos == 1, 1, sum(zero_runs$lengths[1:(pos-1)]) + 1)
  end <- start + zero_runs$lengths[pos] - 1
  cat(sprintf("Plateau interval: diff index [%d, %d], raw data index [%d, %d]\n", 
              start, end, start, end+1))
} 
# Extract Estimate data
estimate_data <- filtered_data[filtered_data$Type == "Estimate", ]

# Define variables
t_values <- estimate_data$t          
stepx_values <- estimate_data$StepX  
ne_values <- estimate_data$Ne       
diff_ne <- diff(ne_values)          

#  Function: Identify left and right boundaries of peaks/valleys(use the minimum left boundary and the maximum right boundary)
get_full_endpoints <- function(peak_idx, ne_values, t_values, stepx_values, is_peak = TRUE) {
  n <- length(ne_values)
  
  left_idx <- peak_idx
  if (is_peak) {
    
    k <- peak_idx - 1
    while (k >= 1 && ne_values[k] <= ne_values[k+1]) {
      left_idx <- k
      k <- k - 1
    }
  } else {
  
    k <- peak_idx - 1
    while (k >= 1 && ne_values[k] >= ne_values[k+1]) {
      left_idx <- k
      k <- k - 1
    }
  }
  
  right_idx <- peak_idx
  if (is_peak) {
   
    k <- peak_idx + 1
    while (k <= n && ne_values[k] <= ne_values[k-1]) {
      right_idx <- k
      k <- k + 1
    }
  } else {
    
    k <- peak_idx + 1
    while (k <= n && ne_values[k] >= ne_values[k-1]) {
      right_idx <- k
      k <- k + 1
    }
  }
  
  return(list(
    left_idx = left_idx,
    right_idx = right_idx,
    left_time = t_values[left_idx],
    right_time = stepx_values[right_idx],
    ne_value = ne_values[peak_idx]
  ))
}

# Store all peak and valley results
results <- data.frame()

zero_runs <- rle(abs(diff_ne) <= 1e-6)
zero_positions <- which(zero_runs$values == TRUE)

processed_indices <- c()

if (length(zero_positions) > 0) {
  for (pos in zero_positions) {
    
    start <- ifelse(pos == 1, 1, sum(zero_runs$lengths[1:(pos-1)]) + 1)
    end <- start + zero_runs$lengths[pos] - 1
    
    left_idx <- start
    right_idx <- end + 1
     
    processed_indices <- c(processed_indices, left_idx:right_idx)
    
    # Compute left and right differences
    left_diff <- if(start > 1) diff_ne[start - 1] else NA
    right_diff <- if(end < length(diff_ne)) diff_ne[end + 1] else NA
    
    # Identify peaks and valleys
    is_peak <- FALSE
    is_valley <- FALSE
    
    if (!is.na(left_diff) && !is.na(right_diff)) {
      is_peak <- (left_diff > 0 && right_diff < 0)
      is_valley <- (left_diff < 0 && right_diff > 0)
    } else if (!is.na(left_diff)) {
      is_peak <- (left_diff > 0)
      is_valley <- (left_diff < 0)
    } else if (!is.na(right_diff)) {
      is_peak <- (right_diff < 0)
      is_valley <- (right_diff > 0)
    }
    
    if (is_peak || is_valley) {
     
      # Determine complete boundaries (use the minimum left boundary and maximum right boundary)
      endpoints <- get_full_endpoints(left_idx, ne_values, t_values, stepx_values, is_peak)
      
      results <- rbind(results, data.frame(
        type = ifelse(is_peak, "peak", "valley"),
        is_plateau = TRUE,
        peak_idx_start = left_idx,
        peak_idx_end = right_idx,
        peak_t_start = t_values[left_idx],
        peak_t_end = stepx_values[right_idx],
        ne_value = ne_values[left_idx],
        plateau_length = right_idx - left_idx + 1,
        left_endpoint_idx = endpoints$left_idx,
        right_endpoint_idx = endpoints$right_idx,
        left_endpoint_time = endpoints$left_time,
        right_endpoint_time = endpoints$right_time
      ))
    }
  }
}

for (i in 2:(length(ne_values)-1)) {

  if (i %in% processed_indices) next
  
  # Determine whether the segment is a peak or a valley
  is_peak <- (ne_values[i-1] < ne_values[i] && ne_values[i] > ne_values[i+1])
  is_valley <- (ne_values[i-1] > ne_values[i] && ne_values[i] < ne_values[i+1])
  
  if (is_peak || is_valley) {
    # Obtain complete boundaries (use the minimum left endpoint and maximum right endpoint)
    endpoints <- get_full_endpoints(i, ne_values, t_values, stepx_values, is_peak)
    
    results <- rbind(results, data.frame(
      type = ifelse(is_peak, "peak", "valley"),
      is_plateau = FALSE,
      peak_idx_start = i,
      peak_idx_end = i,
      peak_t_start = t_values[i],
      peak_t_end = stepx_values[i],
      ne_value = ne_values[i],
      plateau_length = 1,
      left_endpoint_idx = endpoints$left_idx,
      right_endpoint_idx = endpoints$right_idx,
      left_endpoint_time = endpoints$left_time,
      right_endpoint_time = endpoints$right_time
    ))
  }
}

# Output results
cat("\n========== All detected peaks and valleys with their full left and right boundaries ==========\n\n")
for (i in 1:nrow(results)) {
  r <- results[i, ]
  type_name <- ifelse(r$type == "peak", 
                      ifelse(r$is_plateau, "peak plateau", "single peak"), 
                      ifelse(r$is_plateau, "valley plateau", "single valley"))
  
  cat(sprintf("%s %d:\n", type_name, i))
  cat(sprintf("  # Peak and valley positions: t=[%.2f, %.2f], Ne=%.4f\n", 
              r$peak_t_start, r$peak_t_end, r$ne_value))
  cat("\n")
}
