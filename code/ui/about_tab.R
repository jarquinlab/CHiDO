# Created by: Francisco Gonzalez
# Last updated by: Francisco Gonzalez
# Last updated: 06/09/2024

source("code/functions/utils.R")

### UI instructions ----
### Supporting objects ----

about_md <- markdown::renderMarkdown(file="README.Rmd")

### Main tab ----

about_tab <- tabItem(
  tabName = "about",
  h1("CHiDO (Characterization & Integration of Driven Omics)", class="tab-header"),
  HTML(about_md)
)