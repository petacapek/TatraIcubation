pkgs <- c("ggplot2", 
          'openxlsx', 
          'dplyr', 
          'brms', 
          'arm',
          'multcomp',
          'DescTools',
          'reshape2',
          'ggsignif',
          'gridExtra',
          'metafor',
          'ggnewscale',
          'nlme',
          'nlraa',
          'emmeans',
          'multcompView',
          'tidyverse',
          "lubridate")

# install any that are missing, then load all
missing <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(missing)) install.packages(missing)
invisible(lapply(pkgs, library, character.only = TRUE))