
# Venn/Euler Diagram -----------------------------------------------------------
library(eulerr)

# Simplified example with 3 variables
fit <- euler(c(
  "physchem_eucl" = 0.1595,
  "coor_haver_log" = 0.0850,
  "plant_bray" = 0.0067,
  "physchem_eucl&coor_haver_log" = 0.1415,
  "physchem_eucl&plant_bray" = 0.0320,
  "coor_haver_log&plant_bray" = 0,  # small negative; can be set to 0
  "physchem_eucl&coor_haver_log&plant_bray" = 0.1511
))

plot(fit, quantities = TRUE, fills = list(fill = c("#440154FF", "#31688EFF", "#35B779FF")))


# Stacked Barplot of Contributions ---------------------------------------------
# Create a data frame of main contributors (you can adjust as needed)
df <- data.frame(
  Component = c("Unique: physchem_eucl", "Unique: coor_haver_log", "Unique: plant_bray", 
                "Shared: physchem_eucl & coor_haver_log", "Shared: physchem_eucl & plant_bray",
                "Shared: physchem_eucl & coor_haver_log & plant_bray"),
  R2 = c(0.1595, 0.0850, 0.0067, 0.1415, 0.0320, 0.1511)
)

# Plot
ggplot(df, aes(x = "", y = R2, fill = Component)) +
  geom_bar(stat = "identity", width = 0.5) +
  coord_polar("y") +
  theme_void() +
  labs(title = "Partitioned R² Contributions", fill = "Component") +
  scale_fill_brewer(palette = "Set3")


# Heatmap or Matrix Visualization ----------------------------------------------
# Copy-paste your lognormal commonality results here
r2_data <- tribble(
  ~Component, ~R2,
  "Unique: physchem_eucl", 0.1595,
  "Unique: coor_haver_log", 0.0850,
  "Unique: plant_bray", 0.0067,
  "Unique: p_trait_eucl", -0.0004,
  "Unique: physchem_eucl:plant_bray", 0.0113,
  "Shared: physchem_eucl & coor_haver_log", 0.1415,
  "Shared: physchem_eucl & plant_bray", 0.0320,
  "Shared: coor_haver_log & plant_bray", -0.0058,
  "Shared: physchem_eucl & p_trait_eucl", 0.0002,
  "Shared: coor_haver_log & p_trait_eucl", -0.0005,
  "Shared: plant_bray & p_trait_eucl", 0.0000,
  "Shared: physchem_eucl & physchem_eucl:plant_bray", 0.0164,
  "Shared: coor_haver_log & physchem_eucl:plant_bray", -0.0054,
  "Shared: plant_bray & physchem_eucl:plant_bray", -0.0052,
  "Shared: p_trait_eucl & physchem_eucl:plant_bray", -0.0001,
  "Shared: physchem_eucl & coor_haver_log & plant_bray", 0.1511,
  "Shared: physchem_eucl & coor_haver_log & p_trait_eucl", 0.0006,
  "Shared: physchem_eucl & plant_bray & p_trait_eucl", 0.0006,
  "Shared: coor_haver_log & plant_bray & p_trait_eucl", 0.0001,
  "Shared: physchem_eucl & coor_haver_log & physchem_eucl:plant_bray", 0.0914,
  "Shared: physchem_eucl & plant_bray & physchem_eucl:plant_bray", -0.0226,
  "Shared: coor_haver_log & plant_bray & physchem_eucl:plant_bray", 0.0092,
  "Shared: physchem_eucl & p_trait_eucl & physchem_eucl:plant_bray", 0.0002,
  "Shared: coor_haver_log & p_trait_eucl & physchem_eucl:plant_bray", 0.0000,
  "Shared: plant_bray & p_trait_eucl & physchem_eucl:plant_bray", 0.0000,
  "Shared: physchem_eucl & coor_haver_log & plant_bray & p_trait_eucl", 0.0025,
  "Shared: physchem_eucl & coor_haver_log & plant_bray & physchem_eucl:plant_bray", -0.0950,
  "Shared: physchem_eucl & coor_haver_log & p_trait_eucl & physchem_eucl:plant_bray", 0.0000,
  "Shared: physchem_eucl & plant_bray & p_trait_eucl & physchem_eucl:plant_bray", -0.0001,
  "Shared: coor_haver_log & plant_bray & p_trait_eucl & physchem_eucl:plant_bray", 0.0002,
  "Shared: All", -0.0002
)

# Clean up names and split predictors for plotting
r2_long <- r2_data %>%
  separate_rows(Component, sep = "[:&]+") %>%
  mutate(Component = str_trim(str_remove(Component, "Unique:|Shared:"))) %>%
  group_by(across()) %>%
  mutate(Combination = paste(sort(Component), collapse = " + ")) %>%
  ungroup() %>%
  dplyr::select(Combination, R2) %>%
  distinct() %>%
  arrange(desc(R2)) %>%
  mutate(R2_category = case_when(
    R2 > 0.05 ~ "High",
    R2 > 0.01 ~ "Moderate",
    R2 > 0    ~ "Low",
    R2 <= 0   ~ "Negative"
  ))


# Plot heatmap
ggplot(r2_long, aes(x = reorder(Combination, R2), y = 1, fill = R2)) +
  geom_tile(color = "white") +
  scale_fill_viridis(option = "C", direction = -1, name = "R²",
                     limits = c(min(r2_long$R2), max(r2_long$R2))) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title = "Commonality Partitioning (R² Contributions)",
       subtitle = "GLMM Hierarchical Partitioning: Lognormal model",
       x = "Predictor Combination", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 10),
        panel.grid = element_blank())




# Matrix ----------
library(tidyverse)

# Base data
r2_data <- tribble(
  ~Component, ~R2,
  "Unique: physchem_eucl", 0.1595,
  "Unique: coor_haver_log", 0.0850,
  "Unique: plant_bray", 0.0067,
  "Unique: p_trait_eucl", -0.0004,
  "Unique: physchem_eucl:plant_bray", 0.0113,
  "Shared: physchem_eucl & coor_haver_log", 0.1415,
  "Shared: physchem_eucl & plant_bray", 0.0320,
  "Shared: coor_haver_log & plant_bray", -0.0058,
  "Shared: physchem_eucl & p_trait_eucl", 0.0002,
  "Shared: coor_haver_log & p_trait_eucl", -0.0005,
  "Shared: plant_bray & p_trait_eucl", 0.0000,
  "Shared: physchem_eucl & physchem_eucl:plant_bray", 0.0164,
  "Shared: coor_haver_log & physchem_eucl:plant_bray", -0.0054,
  "Shared: plant_bray & physchem_eucl:plant_bray", -0.0052,
  "Shared: p_trait_eucl & physchem_eucl:plant_bray", -0.0001,
  "Shared: physchem_eucl & coor_haver_log & plant_bray", 0.1511,
  "Shared: physchem_eucl & coor_haver_log & p_trait_eucl", 0.0006,
  "Shared: physchem_eucl & plant_bray & p_trait_eucl", 0.0006,
  "Shared: coor_haver_log & plant_bray & p_trait_eucl", 0.0001,
  "Shared: physchem_eucl & coor_haver_log & physchem_eucl:plant_bray", 0.0914,
  "Shared: physchem_eucl & plant_bray & physchem_eucl:plant_bray", -0.0226,
  "Shared: coor_haver_log & plant_bray & physchem_eucl:plant_bray", 0.0092,
  "Shared: physchem_eucl & p_trait_eucl & physchem_eucl:plant_bray", 0.0002,
  "Shared: coor_haver_log & p_trait_eucl & physchem_eucl:plant_bray", 0.0000,
  "Shared: plant_bray & p_trait_eucl & physchem_eucl:plant_bray", 0.0000,
  "Shared: physchem_eucl & coor_haver_log & plant_bray & p_trait_eucl", 0.0025,
  "Shared: physchem_eucl & coor_haver_log & plant_bray & physchem_eucl:plant_bray", -0.0950,
  "Shared: physchem_eucl & coor_haver_log & p_trait_eucl & physchem_eucl:plant_bray", 0.0000,
  "Shared: physchem_eucl & plant_bray & p_trait_eucl & physchem_eucl:plant_bray", -0.0001,
  "Shared: coor_haver_log & plant_bray & p_trait_eucl & physchem_eucl:plant_bray", 0.0002,
  "Shared: All", -0.0002
)

# List of variables to track
variables <- c("physchem_eucl", "coor_haver_log", "plant_bray", "p_trait_eucl", "physchem_eucl:plant_bray")

# Expand the data for heatmap matrix
r2_matrix <- r2_data %>%
  mutate(Predictors = str_remove(Component, "Unique: |Shared: ")) %>%
  mutate(Predictors = str_split(Predictors, " & ")) %>%
  unnest(Predictors) %>%
  mutate(Predictors = str_trim(Predictors)) %>%
  mutate(included = 1) %>%
  pivot_wider(names_from = Predictors, values_from = included, values_fill = 0) %>%
  mutate(R2 = r2_data$R2,
         Component = r2_data$Component) %>%
  dplyr::select(Component, R2, all_of(variables))

# Convert to long format for plotting
heatmap_data <- r2_matrix %>%
  pivot_longer(cols = all_of(variables), names_to = "Variable", values_to = "Included")

# Plot the matrix-style heatmap
ggplot(heatmap_data, aes(x = Variable, y = fct_reorder(Component, R2))) +
  geom_tile(aes(fill = Included), color = "white") +
  geom_text(aes(label = ifelse(Included == 1, "✓", "")), size = 3.5, color = "black") +
  scale_fill_gradient(low = "white", high = "#219EBC") +
  theme_minimal(base_size = 11) +
  labs(
    title = "Matrix Heatmap of Predictor Inclusion in R² Components",
    x = NULL,
    y = "R² Component",
    fill = "Included"
  )



# PEDAL 4 Venn -----------------------------------------------------------------
library(ggvenn)

# Manually create sets for symbolic Venn layout
venn_sets <- list(
  "physchem_eucl" = 1:10,
  "coor_haver_log" = 5:15,
  "plant_bray" = 10:20,
  "p_trait_eucl" = 15:25,
  "physchem_eucl:plant_bray" = 20:30
)

# Plot the 5-set Venn diagram (not proportional)
ggvenn(
  venn_sets,
  fill_color = c("#8ECAE6", "#219EBC", "#FFB703", "#FB8500", "#D4A373"),
  stroke_size = 0.5,
  set_name_size = 5
) +
  theme_void() +
  annotate("text", x = 2.5, y = 2.5, label = "0.1595", size = 5) +  # physchem_eucl
  annotate("text", x = -2.5, y = 2.5, label = "0.0850", size = 5) + # coor_haver_log
  annotate("text", x = -2.5, y = -2.5, label = "0.0067", size = 5) +# plant_bray
  annotate("text", x = 2.5, y = -2.5, label = "-0.0004", size = 5) +# p_trait_eucl
  annotate("text", x = 0, y = 3.8, label = "0.0113", size = 5)      # physchem_eucl:plant_bray



# -------------
library(ggplot2)
library(ggforce)
library(tibble)
library(dplyr)

# Define 5 unique variables and their R² values
r2_unique <- tribble(
  ~label, ~r2,
  "physchem_eucl",       "0.1595",
  "coor_haver_log",      "0.0850",
  "plant_bray",          "0.0067",
  "p_trait_eucl",        "-0.0004",
  "physchem_eucl:plant_bray", "0.0113"
)

# Circle (petal) layout
n <- nrow(r2_unique)
radius <- 2.2             # radius of each circle
inner_radius <- 1.7       # how close circles are to center (smaller = more overlap)

# Calculate positions around a circle
angles <- seq(0, 2*pi, length.out = n+1)[- (n+1)]
r2_unique <- r2_unique %>%
  mutate(
    x = inner_radius * cos(angles),
    y = inner_radius * sin(angles),
    r = radius
  )

# Plot overlapping circular petals
ggplot() +
  geom_circle(data = r2_unique,
              aes(x0 = x, y0 = y, r = r, fill = label),
              alpha = 0.45, color = "black") +
  geom_text(data = r2_unique,
            aes(x = x, y = y + 1.2, label = label),
            size = 4, fontface = "bold") +
  geom_text(data = r2_unique,
            aes(x = x, y = y - 1.2, label = r2),
            size = 4) +
  scale_fill_manual(values = c(
    "physchem_eucl" = "#8ECAE6",
    "coor_haver_log" = "#219EBC",
    "plant_bray" = "#FFB703",
    "p_trait_eucl" = "#FB8500",
    "physchem_eucl:plant_bray" = "#D4A373"
  )) +
  coord_fixed() +
  theme_void() +
  theme(legend.position = "none")


# ------------------------------------
library(ggplot2)
library(ggforce)
library(tibble)
library(dplyr)

# Base petal data
petals <- tribble(
  ~label,                     ~r2,       ~x,   ~y,   ~r,  ~label_x, ~label_y, ~r2_x, ~r2_y,
  "physchem_eucl",            " 0.1595",  0,    1.7, 2.2,  0,        3.2,      0,     0.1,
  "coor_haver_log",           " 0.0850",  1.7,  0.5, 2.2,  3.1,      0.3,      1.7,  -0.8,
  "plant_bray",               " 0.0067",  1,   -1.5, 2.2,  2.0,     -3.2,      1,    -2.9,
  "p_trait_eucl",             "-0.0004", -1,   -1.5, 2.2, -2.6,     -2.2,     -1.3,  -2.9,
  "physchem_eucl:plant_bray", " 0.0113", -1.7,  0.5, 2.2, -3.4,      0.2,     -1.7,  -0.8
)

# Plot
ggplot() +
  geom_circle(data = petals, aes(x0 = x, y0 = y, r = r, fill = label),
              alpha = 0.45, color = "black") +
  geom_text(data = petals, aes(x = label_x, y = label_y, label = label),
            size = 4.5, fontface = "bold") +
  geom_text(data = petals, aes(x = r2_x, y = r2_y, label = r2),
            size = 5) +
  scale_fill_manual(values = c(
    "physchem_eucl" =            "#8ECAE6",
    "coor_haver_log" =           "#219EBC",
    "plant_bray" =               "#FFB703",
    "p_trait_eucl" =             "#FB8500",
    "physchem_eucl:plant_bray" = "#D4A373"
  )) +
  coord_fixed() +
  theme_void() +
  theme(legend.position = "none")


# pedal 5 shape ---------
library(ggplot2)
library(ggforce)
library(tibble)
library(dplyr)

# Define original radius positions
petals <- tribble(
  ~label, ~r2, ~x, ~y,
  "physchem_eucl",       "0.1595",  0,    3,
  "coor_haver_log",      "0.0850",  2.85, 0.93,
  "plant_bray",          "0.0067",  1.76, -2.43,
  "p_trait_eucl",        "-0.0004", -1.76, -2.43,
  "physchem_eucl:plant_bray", "0.0113", -2.85, 0.93
)

# Shrink factor for tight overlap
shrink <- 0.6

petals <- petals %>%
  mutate(
    x = x * shrink,
    y = y * shrink
  ) %>%
  rowwise() %>%
  mutate(
    angle = atan2(-y, -x),  # face toward center
    a = 2.5,                # major axis (length)
    b = 1.2,                # minor axis (width)
    label_x = x * 1.9,      # position for label text
    label_y = y * 1.9,
    r2_x = x * 1.2,         # position for R² value
    r2_y = y * 1.2
  )

# Plot
ggplot() +
  geom_ellipse(data = petals, 
               aes(x0 = x, y0 = y, a = a, b = b, angle = angle, fill = label),
               alpha = 0.45, color = "black") +
  geom_text(data = petals, aes(x = label_x, y = label_y, label = label),
            size = 4.5, fontface = "bold") +
  geom_text(data = petals, aes(x = r2_x, y = r2_y, label = r2),
            size = 5) +
  scale_fill_manual(values = c(
    "physchem_eucl" = "#8ECAE6",
    "coor_haver_log" = "#219EBC",
    "plant_bray" = "#FFB703",
    "p_trait_eucl" = "#FB8500",
    "physchem_eucl:plant_bray" = "#D4A373"
  )) +
  coord_fixed() +
  theme_void() +
  theme(legend.position = "none")



# Pedal on the other side --------
library(ggplot2)
library(ggforce)
library(tibble)
library(dplyr)

# Define petal data
petals <- tribble(
  ~label, ~r2, ~x, ~y,
  "physchem_eucl",       "0.1595",  0,    3,
  "coor_haver_log",      "0.0850",  2.85, 0.93,
  "plant_bray",          "0.0067",  1.76, -2.43,
  "p_trait_eucl",        "-0.0004", -1.76, -2.43,
  "physchem_eucl:plant_bray", "0.0113", -2.85, 0.93
)

# Deeper overlap (more inward)
shrink <- 0.35  # Shrink more for deeper central overlap

petals <- petals %>%
  mutate(
    x = x * shrink,
    y = y * shrink
  ) %>%
  rowwise() %>%
  mutate(
    angle = atan2(-y, -x),  # face toward center
    a = 2.5,                # major axis (length)
    b = 1.2,                # minor axis (width)
    label_x = x * 3.5,      # move labels further out
    label_y = y * 3.5,
    r2_x = x * 1.5,         # and r2 values closer in
    r2_y = y * 1.5
  )

# Plot
ggplot() +
  geom_ellipse(data = petals, 
               aes(x0 = x, y0 = y, a = a, b = b, angle = angle, fill = label),
               alpha = 0.45, color = "black") +
  geom_text(data = petals, aes(x = label_x, y = label_y, label = label),
            size = 4.5, fontface = "bold") +
  geom_text(data = petals, aes(x = r2_x, y = r2_y, label = r2),
            size = 5) +
  scale_fill_manual(values = c(
    "physchem_eucl" = "#8ECAE6",
    "coor_haver_log" = "#219EBC",
    "plant_bray" = "#FFB703",
    "p_trait_eucl" = "#FB8500",
    "physchem_eucl:plant_bray" = "#D4A373"
  )) +
  coord_fixed() +
  theme_void() +
  theme(legend.position = "none")


# Shift----

library(tidyverse)
library(ggforce)

# Your original petals data
petals <- tribble(
  ~label, ~r2, ~x, ~y,
  "physchem_eucl",       "0.1595",  0,    3,
  "coor_haver_log",      "0.0850",  2.85, 0.93,
  "plant_bray",          "0.0067",  1.76, -2.43,
  "p_trait_eucl",        "-0.0004", -1.76, -2.43,
  "physchem_eucl:plant_bray", "0.0113", -2.85, 0.93
)

# Deeper overlap (more inward)
shrink <- 0.35  # Shrink to bring petals inward

petals <- petals %>%
  mutate(
    x = x * shrink,
    y = y * shrink
  ) %>%
  rowwise() %>%
  mutate(
    angle = atan2(-y, -x),  # face petals toward center
    a = 2.5,                # petal major axis length
    b = 1.2,                # petal minor axis width
    # Compute shift left relative to petal orientation:
    # left direction vector is perpendicular to angle: (cos(angle + pi/2), sin(angle + pi/2))
    shift_dist = 0.8,  # how far to shift left — tweak this
    x_shifted = x + shift_dist * cos(angle + pi/2),
    y_shifted = y + shift_dist * sin(angle + pi/2),
    # Label and r2 positions relative to new petal center
    label_x = x_shifted * 3.5,
    label_y = y_shifted * 3.5,
    r2_x = x_shifted * 1.5,
    r2_y = y_shifted * 1.5
  ) %>%
  ungroup()

# Plot with shifted petals
ggplot() +
  geom_ellipse(data = petals, 
               aes(x0 = x_shifted, y0 = y_shifted, a = a, b = b, angle = angle, fill = label),
               alpha = 0.45, color = "black") +
  geom_text(data = petals, aes(x = label_x, y = label_y, label = label),
            size = 4.5, fontface = "bold") +
  geom_text(data = petals, aes(x = r2_x, y = r2_y, label = r2),
            size = 5) +
  scale_fill_manual(values = c(
    "physchem_eucl" = "#8ECAE6",
    "coor_haver_log" = "#219EBC",
    "plant_bray" = "#FFB703",
    "p_trait_eucl" = "#FB8500",
    "physchem_eucl:plant_bray" = "#D4A373"
  )) +
  coord_fixed() +
  theme_void() +
  theme(legend.position = "none")




# TEARS ----------
7