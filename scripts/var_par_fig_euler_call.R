library(eulerr)
library(grid)
library(viridis)

# Create a function to generate the full custom Euler diagram
plot_custom_euler <- function(mod_vp_eu) {
  # Create the Euler diagram object
  p <- plot(mod_vp_eu,
            quantities = FALSE,
            fills = list(fill = viridis(4), alpha = 0.6),
            labels = FALSE,
            main = NULL)
  
  # Start a new page and draw scaled plot
  grid.newpage()
  vp <- viewport(width = 0.8, height = 0.8)
  pushViewport(vp)
  grid.draw(p)
  
  # Custom annotations
  grid.text("Variation Partitioning", x = 0.5,  y = 1.05, gp = gpar(fontsize = 16, fontface = "bold"))
  grid.text("Dispersal",              x = 0.05, y = 0.05, gp = gpar(fontface = "bold"))
  grid.text("Abiotic NP",             x = 0.95, y = 0.05, gp = gpar(fontface = "bold"))
  grid.text("Active Biotic NP",       x = 0.05, y = 0.95, gp = gpar(fontface = "bold"))
  grid.text("Passive Biotic NP",      x = 0.95, y = 0.95, gp = gpar(fontface = "bold"))
  grid.text("Residuals = 0.58",       x = 0.7,  y = -0.05)
  
  # Unique fractions
  grid.text("0.02 **",  x = 0.05, y = 0.15)
  grid.text("0.01",     x = 0.95, y = 0.15)
  grid.text("0.07 ***", x = 0.05, y = 0.85)
  grid.text("0.01",     x = 0.95, y = 0.85)
  
  # 2-way overlaps
  grid.text("0.13", x = 0.32, y = 0.37)
  grid.text("0.01", x = 0.93, y = 0.41)
  grid.text("0.06", x = 0.68, y = 0.72)
  
  # 3-way overlaps
  grid.text("0.02", x = 0.57, y = 0.26)
  grid.text("0.00", x = 0.7, y = 0.22)
  grid.text("0.01", x = 0.51, y = 0.57)
  grid.text("0.04", x = 0.73, y = 0.44)
  
  # 4-way overlap
  grid.text("0.01", x = 0.58, y = 0.42)
}
