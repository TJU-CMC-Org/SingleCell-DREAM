# Distmap installation
dependencies <- c("devtools", "DistMap", "data.table", "doMC", "glmnet", "ggplot2", "ggfortify", "R.utils")
new.packages <- dependencies[!(dependencies %in% installed.packages()[,"Package"])]
if(length(new.packages)) {
    install.packages(new.packages)
}

