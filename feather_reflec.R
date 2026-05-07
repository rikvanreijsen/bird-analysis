# Reset Commands
closeAllConnections()
rm(list = ls())

# Programs
require(colorspace); require(ggplot2); require(pavo)

# Data Preparation
setwd("~/Desktop/R/Data")
reflec <- read.csv2("~/Desktop/R/Data/reflecdata.csv")
reflec <- as.rspec(reflec, whichwl = 1, inter = TRUE, lim = c(300, 700))
reflec.sm <- procspec(reflec, opt = "smooth")
reflec.mean <- aggspec(reflec.sm, by = 3)
reflec.mn.sm <- prospec(reflec.mean, opt = "smooth")
summary(reflec.mn.sm)

# Plots
plots(reflec.mn.sm, col = spec2rgb(reflec.mn.sm))

##### BIRD VIEW
# Construct a Model Using a Bird Viewer
reflec.mn.sm.visbybird <- vismodel(reflec.mn.sm, visual = 'avg.uv', illum = 'D65', vonkries = TRUE, 
                                   relative = FALSE, achromatic = 'none')

# Model Reflectance Spectra in a Color Space
reflec.mn.sm.vis.tsc <- colspace(reflec.mn.sm.visbybird, space = 'tcs')
colorplaette <- spec2rgb(reflec.mn.sm)
colorpalette
plot(reflec.mn.sm.vis.tsc, pch = 21, bg = colorpalette)

##### HUMAN VIEW #####
# Construct a Model Using Human Viewer
reflec.mn.sm.visbyhuman <- vismodel(reflec.mn.sm, visual = 'cie10', illum = 'D65', vonkries = TRUE, 
                                    relative = FALSE, achromatic = 'none')

# Model Reflectance Spectra in a Color Space
reflec.mn.sm.vis.ciexyz <- colspace(reflec.mn.sm.visbyhuman, space = 'ciexyz')
plot(reflec.mn.sm.vis.ciexyz, pch = 21, bg = colorpalette)
