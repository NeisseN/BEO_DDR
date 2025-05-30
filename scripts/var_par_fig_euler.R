library(eulerr)
library(grid)
library(viridis)

# Create the Euler diagram object
p <- plot(mod_vp_eu,
          quantities = FALSE,  # we'll add quantities manually
          fills = list(fill = viridis(4), alpha = 0.6),
          labels = FALSE,
          main = NULL)

# Start a new plotting page and scale it
grid.newpage()
vp <- viewport(width = 0.8, height = 0.8)
pushViewport(vp)
grid.draw(p)

# Add title and main category labels
grid.text("Variation Partitioning", x = 0.5,  y = 1.05, gp = gpar(fontsize = 16, fontface = "bold"))
grid.text("Dispersal",              x = 0.05, y = 0.05, gp = gpar(fontface = "bold"))
grid.text("Abiotic NP",             x = 0.95, y = 0.05, gp = gpar(fontface = "bold"))
grid.text("Active Biotic NP",       x = 0.05, y = 0.95, gp = gpar(fontface = "bold"))
grid.text("Passive Biotic NP",      x = 0.95, y = 0.95, gp = gpar(fontface = "bold"))
grid.text("Residuals = 0.58",       x = 0.7,  y = -0.05)

# === Add the manual numbers from your data (adjust positions as needed) ===

# Unique areas
grid.text("0.02 **", x = 0.21, y = 0.12)  # Distance only
grid.text("0.01", x = 0.84, y = 0.28)  # Abiotic NP only
grid.text("0.07 ***", x = 0.22, y = 0.72)  # Active Biotic NP only
grid.text("0.01", x = 0.78, y = 0.72)  # Passive Biotic NP only

# 2-way overlaps
grid.text("0.03", x = 0.50, y = 0.25)  # Distance & Abiotic
grid.text("0.13", x = 0.35, y = 0.5)   # Distance & Active
grid.text("0.02", x = 0.65, y = 0.5)   # Distance & Passive
grid.text("0.01", x = 0.65, y = 0.6)   # Abiotic & Passive
grid.text("0.06", x = 0.5, y = 0.75)   # Active & Passive

# 3-way overlaps
grid.text("0.02", x = 0.45, y = 0.35)  # Distance + Abiotic + Active
grid.text("0.00", x = 0.55, y = 0.35)  # Distance + Abiotic + Passive
grid.text("0.01", x = 0.4, y = 0.6)    # Distance + Active + Passive
grid.text("0.04", x = 0.55, y = 0.65)  # Abiotic + Active + Passive

# 4-way overlap
grid.text("0.01", x = 0.5, y = 0.5, gp = gpar(fontface = "bold"))