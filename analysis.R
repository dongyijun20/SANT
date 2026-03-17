# ==============================================================================
# Part 1: Environment Setup and Data Loading
# ==============================================================================

if(!require(table1)) install.packages("table1")
if(!require(ggplot2)) install.packages("ggplot2")
if(!require(gridExtra)) install.packages("gridExtra")
if(!require(ggsignif)) install.packages("ggsignif")
if(!require(ggpubr)) install.packages("ggpubr")
if(!require(pROC)) install.packages("pROC")
if(!require(tidyr)) install.packages("tidyr")
if(!require(dplyr)) install.packages("dplyr")

library(table1)
library(ggplot2)
library(gridExtra)
library(ggsignif)
library(ggpubr)
library(pROC)
library(tidyr)
library(dplyr)

if(!dir.exists("figures")) dir.create("figures", recursive = TRUE)
if(!dir.exists("figures/main")) dir.create("figures/main", recursive = TRUE)
if(!dir.exists("figures/supplementary")) dir.create("figures/supplementary", recursive = TRUE)

# Publication style helpers (Arial, title=10, body=8)
plot_font_family <- "Arial"
can_use_arial <- tryCatch({
  tf <- tempfile(fileext = ".pdf")
  grDevices::pdf(tf, family = "Arial")
  plot.new()
  dev.off()
  TRUE
}, error = function(e) FALSE)
if(!can_use_arial) plot_font_family <- "sans"

theme_pub_a4 <- function() {
  theme_classic(base_family = plot_font_family, base_size = 8) +
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 8, hjust = 0.5),
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 8),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 8)
    )
}

a4_width <- 8.27
# Use A4 width with figure-specific heights (avoid over-stretching panels)
fig2_height <- 7.2
fig3_height <- 9.2
fig4_height <- 10.0
fig5_height <- 8.0
figS_height <- 10.0

arrange_with_panel_labels <- function(plot_list, ncol = 1, nrow = NULL) {
  n_plot <- length(plot_list)
  if(is.null(nrow)) nrow <- ceiling(n_plot / ncol)
  ggpubr::ggarrange(
    plotlist = plot_list,
    ncol = ncol,
    nrow = nrow,
    labels = LETTERS[seq_len(n_plot)],
    font.label = list(size = 18, face = "bold", family = plot_font_family),
    hjust = -0.15,
    vjust = 1.2,
    align = "hv"
  )
}

arrange_two_plus_one <- function(pA, pB, pC) {
  p_blank <- ggplot() + theme_void()
  ggpubr::ggarrange(
    plotlist = list(pA, pB, pC, p_blank),
    ncol = 2,
    nrow = 2,
    labels = c("A", "B", "C", ""),
    font.label = list(size = 18, face = "bold", family = plot_font_family),
    hjust = -0.15,
    vjust = 1.2,
    align = "hv"
  )
}

arrange_two_plus_one_compact <- function(pA, pB, pC, bottom_height = 0.82) {
  p_blank <- ggplot() + theme_void()
  top_row <- ggpubr::ggarrange(
    pA, pB,
    ncol = 2,
    labels = c("A", "B"),
    font.label = list(size = 18, face = "bold", family = plot_font_family),
    hjust = -0.15,
    vjust = 1.2,
    align = "hv"
  )
  bottom_row <- ggpubr::ggarrange(
    pC, p_blank,
    ncol = 2,
    labels = c("C", ""),
    font.label = list(size = 18, face = "bold", family = plot_font_family),
    hjust = -0.15,
    vjust = 1.2,
    align = "hv"
  )
  ggpubr::ggarrange(top_row, bottom_row, ncol = 1, heights = c(1, bottom_height))
}

# Load data
df_sub <- read.csv("standard_data.csv", stringsAsFactors = FALSE)

# ==============================================================================
# Part 2: Data Cleaning (strictly exclude Unknown)
# ==============================================================================

# --- A. Basic grouping ---
df_sub$Group <- factor(df_sub$Group, levels=c("SANT", "Splenic Lymphoma"))
df_sub$Gender <- factor(df_sub$Gender, levels=c("M", "F"), labels=c("Male", "Female"))

# --- B. Immunohistochemistry (IHC) ---
# Rule: values other than Positive/Negative (including empty and Unknown) are set to NA
process_ihc_na <- function(x) {
  # Normalize text first
  x <- as.character(x)
  # Convert numeric encodings to text labels (compatibility handling)
  if(any(grepl("^[0-9.]+$", x[!is.na(x)]))) {
    x_num <- as.numeric(x)
    x <- ifelse(x_num > 0, "Positive", "Negative")
  }
  
  # Core step: convert Unknown and empty strings to NA
  x[x == "Unknown" | x == "" | is.na(x)] <- NA
  
  # Set factor levels (excluding Unknown)
  return(factor(x, levels=c("Positive", "Negative")))
}

ihc_vars <- c("CD31", "CD34", "CD8", "SMA", "CD68", "CD20", "CD3", "EBV_EBER")
for(var in ihc_vars) {
  if(var %in% names(df_sub)) {
    df_sub[[paste0(var, "_Status")]] <- process_ihc_na(df_sub[[var]])
  }
}

df_sub$Ki67_Index <- as.numeric(df_sub$Ki67_Index)

# --- C. CT imaging features ---

# 1. CT value and density
df_sub$CT_Value_HU <- as.numeric(df_sub$CT_Value_HU)

df_sub$CT_Density_Type[df_sub$CT_Density_Type == "Unknown" | df_sub$CT_Density_Type == ""] <- NA
df_sub$CT_Density_Type <- factor(df_sub$CT_Density_Type, 
                                 levels=c("Hypodense", "Isodense", "Heterogeneous"))

# 2. Enhancement features
# Prefer the detailed enhancement column
if("CT_Enhancement_Detailed" %in% names(df_sub)){
  enh_col <- "CT_Enhancement_Detailed"
} else {
  enh_col <- "CT_Enhancement_Type" # Backward compatibility with legacy column name
}

# Key point: do not include Unknown as a factor level
df_sub$CT_Enhancement_Display <- df_sub[[enh_col]]
df_sub$CT_Enhancement_Display[df_sub$CT_Enhancement_Display == "Unknown" | df_sub$CT_Enhancement_Display == ""] <- NA

# Set factor order (retain only meaningful categories)
df_sub$CT_Enhancement_Display <- factor(df_sub$CT_Enhancement_Display, 
                                        levels=c("Spoke-wheel", "Progressive", "Delayed", 
                                                 "Iso-enhancing", "Hypo-enhancing", 
                                                 "Heterogeneous", "None"))

# --- D. Ultrasound imaging features ---

# 1. Echogenicity
if("US_Echo_Pattern" %in% names(df_sub)) {
  df_sub$US_Echo_Pattern[df_sub$US_Echo_Pattern == "Unknown" | df_sub$US_Echo_Pattern == ""] <- NA
  # Main categories are assumed to be Hypoechoic and Heterogeneous; add others (e.g., Isoechoic) if needed
  df_sub$US_Echo_Pattern <- factor(df_sub$US_Echo_Pattern, levels=c("Hypoechoic", "Heterogeneous"))
}

# 2. Blood flow
if("US_Blood_Flow_Category" %in% names(df_sub)) {
  # Clean values
  df_sub$US_Blood_Flow_Category[df_sub$US_Blood_Flow_Category == "Sparse"] <- "None/Sparse"
  df_sub$US_Blood_Flow_Category[df_sub$US_Blood_Flow_Category == "Unknown" | df_sub$US_Blood_Flow_Category == ""] <- NA
  
  df_sub$US_Blood_Flow_Category <- factor(df_sub$US_Blood_Flow_Category,
                                          levels=c("Peripheral", "None/Sparse"))
}

# 3. Morphology and boundary
df_sub$US_Shape[df_sub$US_Shape == "Unknown" | df_sub$US_Shape == ""] <- NA
df_sub$US_Shape <- factor(df_sub$US_Shape, levels=c("Regular", "Irregular"))

df_sub$US_Boundary_Clear[df_sub$US_Boundary_Clear == "Unknown" | df_sub$US_Boundary_Clear == ""] <- NA
df_sub$US_Boundary_Clear <- factor(df_sub$US_Boundary_Clear, levels=c("Clear", "Unclear"))


# ==============================================================================
# Part 3: Generate Table 1 (Exclude Unknown from P-value)
# ==============================================================================

# Set labels
label(df_sub$Age)                    <- "Age (years)"
label(df_sub$Gender)                 <- "Gender"
label(df_sub$Gross_Spleen_Size_cm)   <- "Spleen Size (Max Diameter, cm) *"
label(df_sub$Gross_Spleen_Weight_g)  <- "Spleen Weight (g)"
label(df_sub$Gross_Lesion_Size_cm)   <- "Lesion Size (Gross, cm)"
label(df_sub$Ki67_Index)             <- "Ki-67 Index (%) *"

# IHC Labels
label(df_sub$CD31_Status) <- "CD31 (Endothelial)"
label(df_sub$CD34_Status) <- "CD34 (Endothelial)"
label(df_sub$CD8_Status)  <- "CD8 (Littoral/T-cell)"
label(df_sub$SMA_Status)  <- "SMA (Stromal)"
label(df_sub$CD68_Status) <- "CD68 (Histiocytic)"
label(df_sub$CD20_Status) <- "CD20 (B-cell)"
label(df_sub$CD3_Status)  <- "CD3 (T-cell)"
if("EBV_EBER_Status" %in% names(df_sub)) label(df_sub$EBV_EBER_Status) <- "EBV (EBER)"

# Imaging Labels
label(df_sub$CT_Value_HU)            <- "CT Value (Plain, HU)"
label(df_sub$CT_Density_Type)        <- "CT Density Type"
label(df_sub$CT_Enhancement_Display) <- "CT Enhancement Pattern"
if("US_Echo_Pattern" %in% names(df_sub)) label(df_sub$US_Echo_Pattern) <- "US Echogenicity"
label(df_sub$US_Blood_Flow_Category) <- "US Blood Flow Pattern"
label(df_sub$US_Shape)               <- "US Morphology"
label(df_sub$US_Boundary_Clear)      <- "US Boundary"

# Define p-value function (automatically ignores NA)
pvalue_fmt <- function(x, ...) {
  y <- unlist(x)
  g <- factor(rep(1:length(x), times=sapply(x, length)))
  
  if(is.numeric(y)) {
    # Continuous variable: Wilcoxon test (non-normal) or t-test
    tryCatch({
      # Ignore NA
      if(length(na.omit(y)) < 2) return(NA) 
      p <- wilcox.test(y ~ g)$p.value
    }, error=function(e) { NA })
  } else {
    # Categorical variable: Fisher's exact test
    tryCatch({
      # table() ignores NA by default, so this computes p-values on non-missing data
      p <- fisher.test(table(y, g), simulate.p.value=TRUE)$p.value
    }, error=function(e) { NA })
  }
  c("", sub("<", "&lt;", format.pval(p, digits=3, eps=0.001)))
}

# Generate table
tab1 <- table1(~ Age + Gender + 
                 # Gross Pathology
                 Gross_Spleen_Size_cm + Gross_Spleen_Weight_g + Gross_Lesion_Size_cm + 
                 # IHC (Pathology)
                 Ki67_Index + CD31_Status + CD34_Status + CD8_Status + SMA_Status + 
                 CD68_Status + CD20_Status + CD3_Status + 
                 # CT Features
                 CT_Value_HU + CT_Density_Type + CT_Enhancement_Display + 
                 # US Features
                 US_Echo_Pattern + US_Blood_Flow_Category + US_Shape + US_Boundary_Clear | Group, 
               data=df_sub, 
               extra.col=list(`P-value`=pvalue_fmt),
               # render.missing=NULL still displays Missing rows (table1 default behavior)
               # This shows Unknown/missing counts without affecting p-values
               overall=FALSE)

# Output
print(tab1)
# Save CSV
write.csv(as.data.frame(tab1), "Table1_Output_Corrected.csv", row.names = FALSE)


# ==============================================================================
# Part 4: Plots (kept as-is for presentation)
# ==============================================================================
# (Because the main update was Table 1, the plotting structure is reused)
# Reload data for plotting to avoid NA side effects
plot_data <- read.csv("standard_data.csv", stringsAsFactors = FALSE)
plot_data$Group <- factor(plot_data$Group, levels=c("SANT", "Splenic Lymphoma"))

# A. Boxplots
create_sig_boxplot <- function(data, y_var, y_label, title) {
  sub_data <- data[!is.na(data[[y_var]]), ]
  test <- wilcox.test(sub_data[[y_var]] ~ sub_data$Group)
  p_val <- test$p.value
  
  sig_label <- ifelse(p_val < 0.001, "***", ifelse(p_val < 0.01, "**", ifelse(p_val < 0.05, "*", paste0("p=", round(p_val, 3)))))
  y_max <- max(sub_data[[y_var]])
  
  ggplot(sub_data, aes(x = Group, y = .data[[y_var]], fill = Group)) +
    geom_boxplot(width = 0.5, alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
    annotate("text", x = 1.5, y = y_max * 1.1, label = sig_label, size = 5, fontface = "bold") +
    annotate("segment", x = 1, xend = 2, y = y_max * 1.05, yend = y_max * 1.05) +
    scale_fill_manual(values = c("#3498DB", "#E74C3C")) +
    theme_pub_a4() +
    labs(title = title, y = y_label, x = "") + 
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold"))
}

# Generate subplot panels
p1 <- create_sig_boxplot(plot_data, "Age", "Age (years)", "A. Age")
p2 <- create_sig_boxplot(plot_data, "Gross_Spleen_Size_cm", "Size (cm)", "B. Spleen Size")
p3 <- create_sig_boxplot(plot_data, "Gross_Spleen_Weight_g", "Weight (g)", "C. Spleen Weight")
p4 <- create_sig_boxplot(plot_data, "Ki67_Index", "Ki-67 Index (%)", "D. Ki-67 Index")

# Preserve key quantitative comparison panels for main Figure 2
p_sig_age <- p1
p_sig_spleen <- p2
p_sig_weight <- p3
p_sig_ki67 <- p4

# Save boxplots and main Figure 2 (A4 portrait)
combined_box <- grid.arrange(p1, p2, p3, p4, nrow = 2)
fig2_main <- arrange_with_panel_labels(list(p_sig_age, p_sig_spleen, p_sig_weight, p_sig_ki67), ncol = 2, nrow = 2)
ggsave("Figure_Boxplots.png", combined_box, width = 10, height = 10, dpi = 300)
ggsave("figures/Figure_Significant_Comparisons.pdf", fig2_main, width = a4_width, height = fig2_height)
ggsave("figures/main/Figure2_Significant_Comparisons.pdf", fig2_main, width = a4_width, height = fig2_height)

# B. ROC curve (16 cm cutoff)
df_roc <- plot_data %>% filter(!is.na(Gross_Spleen_Size_cm))
df_roc$outcome <- ifelse(df_roc$Group == "Splenic Lymphoma", 1, 0)
roc_obj <- roc(df_roc$outcome, df_roc$Gross_Spleen_Size_cm, levels=c(0, 1), direction="<", quiet=TRUE)
target_cutoff <- 16.0
specific_coords <- coords(roc_obj, x = target_cutoff, input = "threshold", ret = c("threshold", "specificity", "sensitivity"), transpose = FALSE)

roc_plot <- ggroc(roc_obj, legacy.axes = TRUE, color="#E74C3C", linewidth = 1.1) +
  theme_pub_a4() +
  annotate("segment", x = 0, xend = 1, y = 0, yend = 1, color="grey60", linetype="dashed") +
  annotate("point", x = 1 - specific_coords$specificity, y = specific_coords$sensitivity, color = "#2C3E50", size = 4, shape=18) +
  annotate("text", x = 0.65, y = 0.25, label = paste0("Cutoff > ", round(specific_coords$threshold, 1), " cm\nSens = ", round(specific_coords$sensitivity*100,1), "%\nSpec = ", round(specific_coords$specificity*100,1), "%"), size = 4, fontface = "bold") +
  labs(title = "ROC: Spleen Size", x = "1 - Specificity", y = "Sensitivity")

ggsave("Figure_ROC.png", roc_plot, width = 6, height = 6, dpi = 300)
ggsave("figures/Figure_ROC_Spleen_Size.pdf", roc_plot, width = 6.5, height = 6)

# C. Diagnostic scatter plot
p_scatter <- ggplot(plot_data, aes(x = Gross_Lesion_Size_cm, y = Gross_Spleen_Size_cm)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 16.0, ymax = Inf, fill = "#E74C3C", alpha = 0.1) + 
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 16.0, fill = "#3498DB", alpha = 0.1) + 
  geom_point(aes(color = Group, shape = Group), size = 4, alpha = 0.8) +
  geom_hline(yintercept = 16.0, linetype = "solid", color = "black", size = 0.8) +
  scale_color_manual(values = c("#3498DB", "#E74C3C")) +
  theme_pub_a4() +
  labs(title = "Spleen vs Lesion Size", subtitle = "Cutoff > 16.0 cm separates groups", x = "Lesion Size (cm)", y = "Spleen Size (cm)") +
  theme(legend.position = "top", plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("Figure_Scatter.png", p_scatter, width = 8, height = 6, dpi = 300)

# ==============================================================================
# Part 5: Correlation analysis
# ==============================================================================

# Ensure required package is loaded
if(!require(ggpubr)) install.packages("ggpubr")
library(ggpubr)

# Reload data (latest version)
plot_data <- read.csv("standard_data.csv", stringsAsFactors = FALSE)
plot_data$Group <- factor(plot_data$Group, levels=c("SANT", "Splenic Lymphoma"))

# --- Plot A: Measurement accuracy check (imaging vs pathology) ---
# Goal: show preoperative CT/US lesion size agrees with gross pathology measurement.
p_corr1 <- ggplot(plot_data, aes(x = Imaging_Lesion_Size_cm, y = Gross_Lesion_Size_cm)) +
  geom_point(aes(color = Group, shape = Group), size = 3, alpha = 0.8) +
  # Add y=x reference line (dashed)
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
  # Add linear fit with light gray confidence region
  geom_smooth(method = "lm", color = "black", fill = "lightgrey", alpha = 0.2, size = 0.8) +
  # Add correlation coefficient (R and p-value)
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top", size = 5) +
  scale_color_manual(values = c("#3498DB", "#E74C3C")) +
  theme_pub_a4() +
  labs(title = "A. Measurement Accuracy",
       subtitle = "Strong correlation validates preoperative imaging",
       x = "Imaging Lesion Size (cm)",
       y = "Gross Pathology Lesion Size (cm)") +
  theme(legend.position = "top",
        legend.background = element_rect(fill = "white", color = "black"))

# --- Plot B: Tumor burden effect (lesion size vs spleen size) ---
# Goal: illustrate possible mechanism differences between groups.
# SANT (blue): expected positive correlation (larger lesion -> larger spleen; mass effect).
# Lymphoma (red): expected weak/no correlation (small lesion can still coexist with splenomegaly; diffuse infiltration).
p_corr2 <- ggplot(plot_data, aes(x = Gross_Lesion_Size_cm, y = Gross_Spleen_Size_cm)) +
  geom_point(aes(color = Group, shape = Group), size = 3, alpha = 0.8) +
  # Fit trend lines by group
  geom_smooth(aes(color = Group, fill = Group), method = "lm", alpha = 0.15) +
  scale_color_manual(values = c("#3498DB", "#E74C3C")) +
  scale_fill_manual(values = c("#3498DB", "#E74C3C")) +
  # Add group-wise correlation annotations
  stat_cor(aes(color = Group), method = "pearson", label.x.npc = "center", size = 4) +
  theme_pub_a4() +
  labs(title = "B. Tumor Burden Effect",
       subtitle = "SANT shows mass effect; Lymphoma shows diffuse enlargement",
       x = "Lesion Size (cm)",
       y = "Spleen Size (cm)") +
  theme(legend.position = "none")

# --- Combine and save ---
pdf("Figure_Correlations.pdf", width = 12, height = 6)
grid.arrange(p_corr1, p_corr2, nrow = 1)
dev.off()
pdf("figures/Correlation_Plots.pdf", width = 12, height = 6)
grid.arrange(p_corr1, p_corr2, nrow = 1)
dev.off()

# ==============================================================================
# Visualization: Proportions with Missing Values
# ==============================================================================

# --- 1. Load packages ---

library(ggplot2)
library(gridExtra)
library(dplyr)
library(tidyr)
library(scales)

# --- 2. Load data ---
df_plot <- read.csv("standard_data.csv", stringsAsFactors = FALSE)

# --- 3. Data preprocessing (explicitly mark Missing) ---
# Key step: convert all NA/empty/Unknown values to "Missing"

# Set group
df_plot$Group <- factor(df_plot$Group, levels=c("SANT", "Splenic Lymphoma"))

# Define helper for missing-value filling
fill_missing <- function(x) {
  x <- as.character(x)
  x_trim <- trimws(x)
  missing_tokens <- c("", "unknown", "na", "n/a", "nan", "null", "missing")
  x[is.na(x) | tolower(x_trim) %in% missing_tokens] <- "Missing"
  return(x)
}

# Apply cleaning in batch
# Select columns for proportion plots
cols_to_plot <- c("Gender", 
                  "CD31", "CD34", "CD8", "SMA", "CD68", "CD20", "CD3",
                  "CT_Enhancement_Detailed", 
                  "US_Blood_Flow_Category", "US_Echo_Pattern") # Include if echogenicity column exists

# Check columns exist to avoid errors
valid_cols <- cols_to_plot[cols_to_plot %in% names(df_plot)]

for(col in valid_cols) {
  df_plot[[col]] <- fill_missing(df_plot[[col]])
}

# --- 4. Optimize category ordering for display ---
# Set factor levels per variable so Missing is consistently positioned (gray)

# IHC: Positive, Negative, Missing
ihc_cols <- c("CD31", "CD34", "CD8", "SMA", "CD68", "CD20", "CD3")
for(col in ihc_cols) {
  if(col %in% names(df_plot)) {
    # Simple text conversion: 1 -> Positive, 0 -> Negative (for numeric-coded source data)
    # After fill_missing values are characters; handle cases like "1"/"0"
    df_plot[[col]][df_plot[[col]] == "1" | df_plot[[col]] == "1.0"] <- "Positive"
    df_plot[[col]][df_plot[[col]] == "0" | df_plot[[col]] == "0.0"] <- "Negative"
    
    df_plot[[col]] <- factor(df_plot[[col]], levels=c("Positive", "Negative", "Missing"))
  }
}

# CT Enhancement
if("CT_Enhancement_Detailed" %in% names(df_plot)) {
  # Keep long labels organized for cleaner plotting
  df_plot$CT_Enhancement_Detailed <- factor(df_plot$CT_Enhancement_Detailed,
                                            levels=c("Spoke-wheel", "Progressive", "Delayed", "Iso-enhancing", "Hypo-enhancing", "Heterogeneous", "None", "Missing"))
}

# US Blood Flow
if("US_Blood_Flow_Category" %in% names(df_plot)) {
  df_plot$US_Blood_Flow_Category[df_plot$US_Blood_Flow_Category == "Sparse"] <- "None/Sparse"
  df_plot$US_Blood_Flow_Category <- factor(df_plot$US_Blood_Flow_Category,
                                           levels=c("Peripheral", "Central", "None/Sparse", "Missing"))
}

# --- 5. Define general plotting function ---
plot_proportion <- function(data, var, title) {
  # Safety normalization to ensure NA/Missing are unified even after downstream recoding
  data[[var]] <- fill_missing(data[[var]])
  
  # Force Missing to be the top segment in stacked bars (last level)
  if(is.factor(data[[var]])) {
    lv <- levels(data[[var]])
    lv_non_missing <- lv[lv != "Missing"]
  } else {
    lv_non_missing <- sort(unique(data[[var]][data[[var]] != "Missing"]))
  }
  data[[var]] <- factor(data[[var]], levels = c(lv_non_missing, "Missing"))
  
  # Compute within-group proportions
  plot_data <- data %>%
    group_by(Group, .data[[var]]) %>%
    summarise(Count = n(), .groups = 'drop') %>%
    group_by(Group) %>%
    mutate(Prop = Count / sum(Count))
  
  ggplot(plot_data, aes(x = Group, y = Prop, fill = .data[[var]])) +
    geom_bar(stat = "identity", position = "fill", width = 0.7) +
    # Add percentage labels (show only >5% to reduce clutter)
    geom_text(aes(label = ifelse(Prop > 0.05, scales::percent(Prop, accuracy = 1), "")),
              position = position_fill(vjust = 0.5), size = 4, color = "white", fontface = "bold") +
    scale_y_continuous(labels = scales::percent) +
    # Color scheme: Missing in gray; observed categories in color
    scale_fill_manual(values = c(
      "Positive" = "#E74C3C", "Negative" = "#3498DB", 
      "Spoke-wheel" = "#E74C3C", "Progressive" = "#E67E22", "Delayed" = "#F1C40F",
      "Hypo-enhancing" = "#3498DB", "None" = "#2980B9", "Iso-enhancing" = "#95A5A6", "Heterogeneous" = "#8E44AD",
      "Peripheral" = "#E74C3C", "Central" = "#3498DB", "None/Sparse" = "#95A5A6",
      "Missing" = "grey80"
    )) +
    theme_pub_a4() +
    labs(title = title, x = "", y = "Proportion", fill = "") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
          legend.position = "right",
          axis.text.x = element_text(face = "bold"))
}

# --- 6. Generate plots ---

# A. Key IHC markers
p1 <- plot_proportion(df_plot, "CD31", "A. CD31 (Vascular)")
p2 <- plot_proportion(df_plot, "CD34", "B. CD34 (Vascular)")
p3 <- plot_proportion(df_plot, "CD20", "C. CD20 (B-cell)")
p4 <- plot_proportion(df_plot, "CD8",  "D. CD8")

# B. Imaging features
# Plot only if the column exists
plot_list_img <- list()
if("CT_Enhancement_Detailed" %in% names(df_plot)) {
  p_ct <- plot_proportion(df_plot, "CT_Enhancement_Detailed", "E. CT Enhancement")
  plot_list_img[[1]] <- p_ct
}
if("US_Blood_Flow_Category" %in% names(df_plot)) {
  p_us <- plot_proportion(df_plot, "US_Blood_Flow_Category", "F. US Blood Flow")
  plot_list_img[[2]] <- p_us
}

# --- 7. Combine and save ---
# Figure 3 should include all proportion panels (no supplementary proportion files)
all_prop_list <- list(p1, p2, p3, p4)
if(exists("p_ct")) all_prop_list[[length(all_prop_list) + 1]] <- p_ct
if(exists("p_us")) all_prop_list[[length(all_prop_list) + 1]] <- p_us

if(length(all_prop_list) >= 2) {
  ncol_prop <- 2
  nrow_prop <- ceiling(length(all_prop_list) / ncol_prop)
  fig3_main <- arrange_with_panel_labels(all_prop_list, ncol = ncol_prop, nrow = nrow_prop)
  ggsave("figures/main/Figure3_Key_Proportions.pdf", fig3_main, width = a4_width, height = fig3_height)
}

# Avoid errors when no graphics device is open
if(!is.null(dev.list())) dev.off()

# ==============================================================================
# Part 6: Additional analysis figures
# ==============================================================================

if(!dir.exists("figures")) dir.create("figures", recursive = TRUE)

# Reload data to keep this section self-contained and consistent
df_extra <- read.csv("standard_data.csv", stringsAsFactors = FALSE)
df_extra$Group <- factor(df_extra$Group, levels = c("SANT", "Splenic Lymphoma"))
df_extra$outcome <- ifelse(df_extra$Group == "Splenic Lymphoma", 1, 0)

# --- A. ROC comparison across continuous variables ---
roc_vars <- c("Age", "Gross_Spleen_Size_cm", "Gross_Spleen_Weight_g", "Ki67_Index", "CT_Value_HU", "Imaging_Lesion_Size_cm")
roc_vars <- roc_vars[roc_vars %in% names(df_extra)]

safe_build_roc <- function(data, var_name) {
  sub <- data[!is.na(data[[var_name]]) & !is.na(data$outcome), c("outcome", var_name)]
  if(nrow(sub) < 8) return(NULL)
  if(length(unique(sub$outcome)) < 2) return(NULL)
  if(length(unique(sub[[var_name]])) < 3) return(NULL)
  
  roc_obj <- tryCatch(
    pROC::roc(sub$outcome, sub[[var_name]], levels = c(0, 1), direction = "<", quiet = TRUE),
    error = function(e) NULL
  )
  if(is.null(roc_obj)) return(NULL)
  
  auc_val <- as.numeric(pROC::auc(roc_obj))
  ci_auc <- tryCatch(as.numeric(pROC::ci.auc(roc_obj)), error = function(e) c(NA, NA, NA))
  
  auc_label <- paste0(
    var_name, " (AUC=", round(auc_val, 3),
    if(all(!is.na(ci_auc))) paste0(", 95%CI ", round(ci_auc[1], 3), "-", round(ci_auc[3], 3)) else "",
    ")"
  )
  
  list(
    roc_obj = roc_obj,
    summary = data.frame(
      Variable = var_name,
      N = nrow(sub),
      AUC = auc_val,
      AUC_CI_L = ci_auc[1],
      AUC_CI_U = ci_auc[3],
      stringsAsFactors = FALSE
    ),
    legend_label = auc_label
  )
}

roc_res <- lapply(roc_vars, function(v) safe_build_roc(df_extra, v))
roc_res <- roc_res[!sapply(roc_res, is.null)]

if(length(roc_res) > 0) {
  roc_summary <- bind_rows(lapply(roc_res, function(x) x$summary))
  legend_labels <- setNames(
    sapply(roc_res, function(x) x$legend_label),
    sapply(roc_res, function(x) x$summary$Variable)
  )
  roc_list <- setNames(
    lapply(roc_res, function(x) x$roc_obj),
    sapply(roc_res, function(x) x$summary$Variable)
  )
  
  # Use ggroc(list) so each ROC is drawn in the correct threshold order
  p_roc_multi <- ggroc(roc_list, legacy.axes = TRUE, linewidth = 1.1) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey60") +
    scale_color_brewer(palette = "Dark2", labels = legend_labels) +
    coord_equal() +
    theme_pub_a4() +
    labs(
      title = "Supplementary Figure A. ROC Comparison Across Continuous Predictors",
      x = "1 - Specificity",
      y = "Sensitivity",
      color = ""
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "bottom"
    )
  
  ggsave("figures/Figure_ROC_Comparison.pdf", p_roc_multi, width = 9, height = 7)
  ggsave("figures/supplementary/Figure_ROC_Comparison.pdf", p_roc_multi, width = 9, height = 7)
  write.csv(roc_summary[order(-roc_summary$AUC), ], "figures/ROC_Comparison_Summary.csv", row.names = FALSE)
}

# --- B. Categorical-feature OR forest (lymphoma direction) ---
feature_binary <- list(
  Gender_Male = c("Gender", "M"),
  CD31_Positive = c("CD31", "Positive"),
  CD34_Positive = c("CD34", "Positive"),
  CD8_Positive = c("CD8", "Positive"),
  CD20_Positive = c("CD20", "Positive"),
  CD3_Positive = c("CD3", "Positive"),
  CT_Density_Heterogeneous = c("CT_Density_Type", "Heterogeneous"),
  CT_Enhancement_Progressive = c("CT_Enhancement_Detailed", "Progressive"),
  US_BloodFlow_Peripheral = c("US_Blood_Flow_Category", "Peripheral"),
  US_Shape_Irregular = c("US_Shape", "Irregular"),
  US_Boundary_Unclear = c("US_Boundary_Clear", "Unclear")
)

safe_fisher_or <- function(data, var_col, target_level, feature_name) {
  if(!(var_col %in% names(data))) return(NULL)
  x <- as.character(data[[var_col]])
  x[x == "" | x == "Unknown"] <- NA
  y <- data$outcome
  keep <- !is.na(x) & !is.na(y)
  if(sum(keep) < 6) return(NULL)
  
  x_bin <- ifelse(x[keep] == target_level, 1, 0)
  y_bin <- y[keep]
  
  tb <- table(x_bin, y_bin)
  if(!all(dim(tb) == c(2, 2))) return(NULL)
  
  ft <- tryCatch(fisher.test(tb), error = function(e) NULL)
  if(is.null(ft)) return(NULL)
  
  data.frame(
    Feature = feature_name,
    N = sum(keep),
    OR = unname(ft$estimate),
    CI_L = ft$conf.int[1],
    CI_U = ft$conf.int[2],
    P = ft$p.value,
    stringsAsFactors = FALSE
  )
}

or_res <- lapply(names(feature_binary), function(nm) {
  item <- feature_binary[[nm]]
  safe_fisher_or(df_extra, item[1], item[2], nm)
})
or_res <- bind_rows(or_res)

if(nrow(or_res) > 0) {
  # Clip OR=0/Inf for visualization to avoid log-scale plotting errors
  lower_clip <- 0.05
  upper_clip <- 20
  
  or_res <- or_res %>%
    mutate(
      log2OR = log2(OR),
      OR_plot = pmin(pmax(OR, lower_clip), upper_clip),
      CI_L_plot = pmin(pmax(CI_L, lower_clip), upper_clip),
      CI_U_plot = pmin(pmax(CI_U, lower_clip), upper_clip),
      Feature = factor(Feature, levels = Feature[order(OR_plot)])
    )
  
  p_or <- ggplot(or_res, aes(x = Feature, y = OR_plot, ymin = CI_L_plot, ymax = CI_U_plot)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
    geom_pointrange(aes(color = P < 0.05), size = 0.6) +
    scale_y_log10() +
    scale_color_manual(values = c("TRUE" = "#E74C3C", "FALSE" = "#2C3E50"), guide = "none") +
    coord_flip() +
    theme_pub_a4() +
    labs(
      title = "Supplementary Figure B. Odds Ratios for Splenic Lymphoma",
      subtitle = "Univariable Fisher exact test (reference OR=1, clipped for display)",
      x = "",
      y = "Odds Ratio (log scale)"
    ) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  ggsave("figures/Figure_OR_Forest.pdf", p_or, width = 8, height = 6)
  write.csv(or_res %>% arrange(P), "figures/OR_Forest_Summary.csv", row.names = FALSE)
}

# ==============================================================================
# Part 7: Combined Model ROC + Bootstrap Stability
# ==============================================================================

candidate_predictors <- c("Ki67_Index", "Gross_Spleen_Size_cm", "Age", "Imaging_Lesion_Size_cm", "Gross_Spleen_Weight_g")
candidate_predictors <- candidate_predictors[candidate_predictors %in% names(df_extra)]

build_model_dataset <- function(data, predictors) {
  keep_cols <- c("outcome", predictors)
  sub <- data[complete.cases(data[, keep_cols]), keep_cols, drop = FALSE]
  if(nrow(sub) < 10) return(NULL)
  if(length(unique(sub$outcome)) < 2) return(NULL)
  sub
}

# Prefer 3-variable combinations, then 2-variable combinations; pick the one with largest usable N
choose_best_combo <- function(data, candidates, k) {
  if(length(candidates) < k) return(NULL)
  combos <- combn(candidates, k, simplify = FALSE)
  
  best <- NULL
  best_n <- -1
  for(v in combos) {
    sub <- build_model_dataset(data, v)
    if(is.null(sub)) next
    if(nrow(sub) > best_n) {
      best <- list(vars = v, data = sub)
      best_n <- nrow(sub)
    }
  }
  best
}

best3 <- choose_best_combo(df_extra, candidate_predictors, 3)
best2 <- choose_best_combo(df_extra, candidate_predictors, 2)

model_choice <- if(!is.null(best3)) best3 else best2
if(!is.null(best3) && !is.null(best2) && nrow(best2$data) > (nrow(best3$data) + 3)) {
  # If 2-variable model has much better data coverage, prefer robustness over complexity
  model_choice <- best2
}

model_vars <- if(!is.null(model_choice)) model_choice$vars else NULL
model_df <- if(!is.null(model_choice)) model_choice$data else NULL

if(!is.null(model_df)) {
  model_formula <- as.formula(paste("outcome ~", paste(model_vars, collapse = " + ")))
  fit_model <- glm(model_formula, data = model_df, family = binomial())
  
  pred_prob <- predict(fit_model, type = "response")
  roc_model <- pROC::roc(model_df$outcome, pred_prob, levels = c(0, 1), direction = "<", quiet = TRUE)
  auc_model <- as.numeric(pROC::auc(roc_model))
  ci_model <- tryCatch(as.numeric(pROC::ci.auc(roc_model)), error = function(e) c(NA, NA, NA))
  model_coords <- tryCatch(
    pROC::coords(
      roc_model,
      x = "best",
      best.method = "youden",
      ret = c("threshold", "specificity", "sensitivity"),
      transpose = FALSE
    ),
    error = function(e) NULL
  )
  
  p_model_roc <- ggroc(roc_model, legacy.axes = TRUE, color = "#8E44AD", linewidth = 1.2) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey60") +
    coord_equal() +
    theme_pub_a4() +
    labs(
      title = "Supplementary Figure C. Combined Model ROC",
      subtitle = paste0(
        "Predictors: ", paste(model_vars, collapse = " + "), 
        " | AUC=", round(auc_model, 3),
        if(all(!is.na(ci_model))) paste0(" (95%CI ", round(ci_model[1], 3), "-", round(ci_model[3], 3), ")") else ""
      ),
      x = "1 - Specificity",
      y = "Sensitivity"
    ) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  if(!is.null(model_coords)) {
    p_model_roc <- p_model_roc +
      annotate(
        "point",
        x = 1 - model_coords$specificity,
        y = model_coords$sensitivity,
        color = "#2C3E50",
        size = 3.5,
        shape = 18
      ) +
      annotate(
        "text",
        x = 0.62,
        y = 0.22,
        label = paste0(
          "Cutoff > ", round(model_coords$threshold, 3),
          "\nSens = ", round(model_coords$sensitivity * 100, 1), "%",
          "\nSpec = ", round(model_coords$specificity * 100, 1), "%"
        ),
        size = 2.8,
        fontface = "bold"
      )
  }
  
  ggsave("figures/Figure_ROC_Combined_Model.pdf", p_model_roc, width = 7, height = 6)
  
  # Bootstrap robustness: refit on bootstrap sample and evaluate AUC on original modeling set
  set.seed(20260317)
  n_boot <- 1000
  auc_boot <- rep(NA_real_, n_boot)
  n_model <- nrow(model_df)
  
  for(i in seq_len(n_boot)) {
    idx <- sample(seq_len(n_model), size = n_model, replace = TRUE)
    boot_df <- model_df[idx, , drop = FALSE]
    
    if(length(unique(boot_df$outcome)) < 2) next
    
    fit_boot <- tryCatch(glm(model_formula, data = boot_df, family = binomial()), error = function(e) NULL)
    if(is.null(fit_boot)) next
    
    pred_boot_on_orig <- tryCatch(predict(fit_boot, newdata = model_df, type = "response"), error = function(e) NULL)
    if(is.null(pred_boot_on_orig)) next
    
    roc_boot <- tryCatch(pROC::roc(model_df$outcome, pred_boot_on_orig, levels = c(0, 1), direction = "<", quiet = TRUE), error = function(e) NULL)
    if(is.null(roc_boot)) next
    
    auc_boot[i] <- as.numeric(pROC::auc(roc_boot))
  }
  
  auc_boot <- auc_boot[!is.na(auc_boot)]
  if(length(auc_boot) > 50) {
    boot_summary <- data.frame(
      Metric = c("Apparent_AUC", "Bootstrap_Median", "Bootstrap_CI_L", "Bootstrap_CI_U", "Bootstrap_Valid_N"),
      Value = c(
        auc_model,
        median(auc_boot),
        quantile(auc_boot, 0.025),
        quantile(auc_boot, 0.975),
        length(auc_boot)
      ),
      stringsAsFactors = FALSE
    )
    
    boot_df_plot <- data.frame(AUC = auc_boot)
    p_boot_plot <- ggplot(boot_df_plot, aes(x = AUC)) +
      geom_histogram(bins = 30, fill = "#3498DB", color = "white", alpha = 0.85) +
      geom_vline(xintercept = auc_model, color = "#E74C3C", linetype = "dashed", linewidth = 1) +
      geom_vline(xintercept = quantile(auc_boot, c(0.025, 0.975)), color = "#2C3E50", linetype = "dotted", linewidth = 0.9) +
      theme_pub_a4() +
      labs(
        title = "Supplementary Figure D. Bootstrap AUC Stability",
        subtitle = paste0(
          "n=", n_model, ", bootstrap=", length(auc_boot),
          " | median AUC=", round(median(auc_boot), 3),
          ", 95% percentile CI=", round(quantile(auc_boot, 0.025), 3), "-", round(quantile(auc_boot, 0.975), 3)
        ),
        x = "AUC (bootstrap-refit model on original dataset)",
        y = "Count"
      ) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    
    ggsave("figures/Figure_Bootstrap_AUC_Stability.pdf", p_boot_plot, width = 8, height = 6)
    ggsave("figures/supplementary/Figure_Bootstrap_AUC_Stability.pdf", p_boot_plot, width = 8, height = 6)
    write.csv(boot_summary, "figures/Bootstrap_AUC_Summary.csv", row.names = FALSE)
  }
}

# ==============================================================================
# Part 8: Final layout exports (Main vs Supplementary)
# ==============================================================================

# Figure 4 (main): diagnostic framework in portrait A4 (no combined model panel)
if(exists("p_corr1") && exists("roc_plot") && exists("p_scatter")) {
  p4a_main <- p_corr1 +
    theme_pub_a4() +
    theme(legend.position = "top") +
    labs(
      title = "Measurement Correlation",
      subtitle = "Imaging lesion size vs gross pathology size"
    )
  
  p4b_main <- roc_plot +
    theme_pub_a4() +
    labs(
      title = "Single-Marker ROC",
      subtitle = "Spleen size with 16 cm clinical cutoff"
    )
  
  p4c_main <- p_scatter +
    theme_pub_a4() +
    labs(
      title = "Lesion vs Spleen Size",
      subtitle = "16 cm threshold highlights splenic-burden separation"
    )
  
  fig4_main <- arrange_two_plus_one(p4a_main, p4b_main, p4c_main)
  ggsave("figures/main/Figure4_Diagnostic_Framework.pdf", fig4_main, width = a4_width, height = fig4_height)
}

# ==============================================================================
# Part 9: Exploratory model utilities for supplementary figures
# ==============================================================================
if(exists("fit_model") && !is.null(model_df) && !is.null(model_vars)) {
  # Panel A: points-based nomogram representation
  coef_vec <- coef(fit_model)[model_vars]
  coef_abs <- abs(coef_vec)
  if(max(coef_abs, na.rm = TRUE) == 0) coef_abs <- rep(1, length(coef_abs))
  
  nom_df <- do.call(rbind, lapply(model_vars, function(v) {
    rng <- range(model_df[[v]], na.rm = TRUE)
    beta <- unname(coef_vec[v])
    max_pts <- abs(beta) * abs(diff(rng))
    data.frame(
      Variable = v,
      MinValue = rng[1],
      MaxValue = rng[2],
      MaxPoints = max_pts,
      stringsAsFactors = FALSE
    )
  }))
  
  if(max(nom_df$MaxPoints, na.rm = TRUE) > 0) {
    nom_df$MaxPoints <- 100 * nom_df$MaxPoints / max(nom_df$MaxPoints, na.rm = TRUE)
  }
  nom_df$Variable <- factor(nom_df$Variable, levels = rev(nom_df$Variable))
  
  p5a_nomogram <- ggplot(nom_df, aes(y = Variable)) +
    geom_segment(aes(x = 0, xend = 1, yend = Variable), linewidth = 0.8, color = "#2C3E50") +
    geom_text(aes(x = 0, label = sprintf("min %.2f", MinValue)), hjust = 1.05, vjust = 1.2, size = 2.6, family = plot_font_family) +
    geom_text(aes(x = 1, label = sprintf("max %.2f", MaxValue)), hjust = -0.05, vjust = 1.2, size = 2.6, family = plot_font_family) +
    geom_point(aes(x = 0.5, size = MaxPoints), color = "#8E44AD", alpha = 0.85) +
    geom_text(aes(x = 0.5, label = paste0(round(MaxPoints, 1), " pts")), vjust = -1.2, size = 2.5, family = plot_font_family) +
    scale_size_continuous(range = c(2.5, 7), guide = "none") +
    coord_cartesian(xlim = c(-0.2, 1.2)) +
    theme_pub_a4() +
    labs(
      title = "Nomogram Representation",
      subtitle = "Predictor contribution scaled to 0-100 points",
      x = "Predictor range",
      y = ""
    )
  
  # Panel B: bootstrap calibration
  pred_app <- predict(fit_model, type = "response")
  obs_app <- model_df$outcome
  grid_prob <- seq(0.01, 0.99, by = 0.01)
  
  lo_app <- loess(obs_app ~ pred_app, span = 0.75, degree = 1)
  cal_app <- pmin(pmax(predict(lo_app, newdata = data.frame(pred_app = grid_prob)), 0), 1)
  
  n_boot_cal <- 300
  optimism_mat <- matrix(NA_real_, nrow = n_boot_cal, ncol = length(grid_prob))
  for(b in seq_len(n_boot_cal)) {
    idx <- sample(seq_len(nrow(model_df)), size = nrow(model_df), replace = TRUE)
    boot_df <- model_df[idx, , drop = FALSE]
    if(length(unique(boot_df$outcome)) < 2) next
    
    fit_b <- tryCatch(glm(model_formula, data = boot_df, family = binomial()), error = function(e) NULL)
    if(is.null(fit_b)) next
    
    p_boot <- tryCatch(predict(fit_b, newdata = boot_df, type = "response"), error = function(e) NULL)
    p_test <- tryCatch(predict(fit_b, newdata = model_df, type = "response"), error = function(e) NULL)
    if(is.null(p_boot) || is.null(p_test)) next
    
    lo_boot <- tryCatch(loess(boot_df$outcome ~ p_boot, span = 0.75, degree = 1), error = function(e) NULL)
    lo_test <- tryCatch(loess(model_df$outcome ~ p_test, span = 0.75, degree = 1), error = function(e) NULL)
    if(is.null(lo_boot) || is.null(lo_test)) next
    
    curve_boot <- pmin(pmax(predict(lo_boot, newdata = data.frame(p_boot = grid_prob)), 0), 1)
    curve_test <- pmin(pmax(predict(lo_test, newdata = data.frame(p_test = grid_prob)), 0), 1)
    optimism_mat[b, ] <- curve_boot - curve_test
  }
  
  optimism_mean <- colMeans(optimism_mat, na.rm = TRUE)
  optimism_mean[is.na(optimism_mean)] <- 0
  cal_corrected <- pmin(pmax(cal_app - optimism_mean, 0), 1)
  
  cal_df <- data.frame(
    Predicted = rep(grid_prob, 3),
    Observed = c(grid_prob, cal_app, cal_corrected),
    Curve = rep(c("Ideal", "Apparent", "Bias-corrected"), each = length(grid_prob)),
    stringsAsFactors = FALSE
  )
  
  p5b_cal <- ggplot(cal_df, aes(x = Predicted, y = Observed, color = Curve, linetype = Curve)) +
    geom_line(linewidth = 0.9) +
    scale_color_manual(values = c("Ideal" = "grey45", "Apparent" = "#E74C3C", "Bias-corrected" = "#2C3E50")) +
    theme_pub_a4() +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(
      title = "Calibration",
      subtitle = "Bootstrap-corrected calibration curve",
      x = "Predicted probability",
      y = "Observed probability",
      color = "",
      linetype = ""
    )
  
  # Panel C: decision curve analysis
  thresh <- seq(0.05, 0.95, by = 0.01)
  n_obs <- nrow(model_df)
  prev <- mean(model_df$outcome)
  pred_model <- pred_app
  
  nb_model <- sapply(thresh, function(t) {
    pred_pos <- pred_model >= t
    tp <- sum(pred_pos & model_df$outcome == 1)
    fp <- sum(pred_pos & model_df$outcome == 0)
    (tp / n_obs) - (fp / n_obs) * (t / (1 - t))
  })
  nb_all <- prev - (1 - prev) * (thresh / (1 - thresh))
  nb_none <- rep(0, length(thresh))
  
  dca_df <- data.frame(
    Threshold = rep(thresh, 3),
    NetBenefit = c(nb_model, nb_all, nb_none),
    Strategy = rep(c("Model", "Treat-all", "Treat-none"), each = length(thresh)),
    stringsAsFactors = FALSE
  )
  
  p5c_dca <- ggplot(dca_df, aes(x = Threshold, y = NetBenefit, color = Strategy, linetype = Strategy)) +
    geom_line(linewidth = 0.9) +
    scale_color_manual(values = c("Model" = "#8E44AD", "Treat-all" = "#E67E22", "Treat-none" = "grey45")) +
    theme_pub_a4() +
    labs(
      title = "Decision Curve Analysis",
      subtitle = "Net benefit across threshold probabilities",
      x = "Threshold probability",
      y = "Net benefit",
      color = "",
      linetype = ""
    )
  
  fig5_supp <- arrange_two_plus_one_compact(p5a_nomogram, p5b_cal, p5c_dca, bottom_height = 0.78)
  ggsave(
    "figures/supplementary/Supplementary_FigureS4_Clinical_Utility_Exploratory.pdf",
    fig5_supp,
    width = a4_width,
    height = fig5_height
  )
  
  # Remove outdated main-figure export to avoid manuscript confusion
  old_fig5_main <- "figures/main/Figure5_Clinical_Utility.pdf"
  if(file.exists(old_fig5_main)) file.remove(old_fig5_main)
}

# Writing guidance file to align legend style and result flow
layout_notes <- c(
  "Main Figure Sequence Recommendation",
  "Figure 1: Representative imaging cases (clinical phenotype).",
  "Figure 2: Significant comparisons (age, spleen size, spleen weight, Ki-67).",
  "Figure 3: Unified feature proportions (all proportion panels; Missing/NA harmonized).",
  "Figure 4: Diagnostic framework = measurement correlation + spleen-size ROC + lesion-vs-spleen threshold scatter.",
  "",
  "Supplementary Figures",
  "S1: ROC comparison across continuous predictors (single panel, readability-prioritized).",
  "S2: Odds-ratio forest + bootstrap AUC stability (adaptive stacked layout).",
  "S3: Correlation analysis extension (measurement + burden effects).",
  "S4: Exploratory model utility panels (nomogram + calibration + decision curve analysis).",
  "",
  "Results paragraph order",
  "1) Qualitative imaging distinction -> 2) quantitative group differences -> 3) diagnostic framework without multivariable model dependence -> 4) exploratory model utilities in supplementary."
)
writeLines(layout_notes, "figures/main/Figure_Layout_and_Legend_Notes.txt")

# Clean legacy supplementary export to prevent stale-file confusion
legacy_supp <- "figures/supplementary/Supplementary_Figure1_Diagnostic_Extensions.pdf"
if(file.exists(legacy_supp)) file.remove(legacy_supp)

# ==============================================================================
# Part 10: Supplementary figure layout (adaptive, readability-first)
# ==============================================================================
if(exists("p_roc_multi")) {
  s1_roc <- p_roc_multi +
    theme_pub_a4() +
    labs(
      title = "Supplementary Figure S1. ROC Comparison Across Continuous Predictors",
      subtitle = "Single-panel layout to preserve curve and legend readability"
    )
  ggsave(
    "figures/supplementary/Supplementary_FigureS1_ROC_Comparison.pdf",
    s1_roc,
    width = a4_width,
    height = 6.0
  )
}

if(exists("p_or") && exists("p_boot_plot")) {
  s2_or <- p_or +
    theme_pub_a4() +
    labs(
      title = "Supplementary Figure S2A. Odds-Ratio Forest",
      subtitle = "Univariable categorical-feature effects"
    )
  s2_boot <- p_boot_plot +
    theme_pub_a4() +
    labs(
      title = "Supplementary Figure S2B. Bootstrap AUC Stability",
      subtitle = "Internal robustness of combined model"
    )
  s2_combo <- ggpubr::ggarrange(
    s2_or, s2_boot,
    ncol = 1,
    labels = c("A", "B"),
    font.label = list(size = 18, face = "bold", family = plot_font_family),
    hjust = -0.15,
    vjust = 1.2,
    heights = c(1.15, 1.0),
    align = "v"
  )
  ggsave(
    "figures/supplementary/Supplementary_FigureS2_OR_Bootstrap.pdf",
    s2_combo,
    width = a4_width,
    height = 8.2
  )
}

if(exists("p_corr2") && exists("p_scatter")) {
  s3_corr <- ggpubr::ggarrange(
    p_corr2 + theme_pub_a4() + labs(title = "Tumor Burden Effect by Group"),
    p_scatter + theme_pub_a4() + labs(title = "Lesion vs Spleen Size with 16 cm Threshold"),
    ncol = 1,
    labels = c("A", "B"),
    font.label = list(size = 18, face = "bold", family = plot_font_family),
    hjust = -0.15,
    vjust = 1.2,
    heights = c(1, 1)
  )
  ggsave(
    "figures/supplementary/Supplementary_FigureS3_Correlations.pdf",
    s3_corr,
    width = a4_width,
    height = 8.2
  )
}