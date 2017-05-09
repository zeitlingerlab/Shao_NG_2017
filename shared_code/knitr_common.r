library(pander)
panderOptions("table.split.table", Inf)

options(knitr.project_name = "Transcription initiation and Pol II pausing")

figure_path <- function(filename="") {
  file.path(getOption("knitr.figure_dir"), filename)
}

# Format number with commas
pn <- function(i, ...) {
  prettyNum(i, big.mark=",", ...)
}

# Wrap output and code
options(width=80)

# Force knitr to stop evaluation when an error is encountered
knitr::opts_chunk$set(error=FALSE)

# Don't show code blocks by default
knitr::opts_chunk$set(echo=FALSE)

# Don't show warnings by default
knitr::opts_chunk$set(warning=FALSE)

# Don't show messages by default
knitr::opts_chunk$set(message=FALSE)

# Don't reformat R code
knitr::opts_chunk$set(tidy=FALSE)

# Set up figure defaults 
knitr::opts_chunk$set(fig.width=7, fig.height=5, fig.path=figure_path(), dev=c("png", "pdf"))

# Create output directory if it doesn't exist
if(!file.exists(getOption("knitr.figure_dir"))) dir.create(getOption("knitr.figure_dir"))

source("shared_code/caching.r")

