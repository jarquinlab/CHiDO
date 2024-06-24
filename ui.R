# Created by: Francisco Gonzalez
# Last updated by: Francisco Gonzalez
# Last updated: 06/09/2024

source("ui/data_tab.R")
source("ui/model_tab.R")
source("ui/cv_tab.R")
source("ui/results_tab.R")
source("ui/about_tab.R")

body <- dashboardBody(
  # Enable JavaScript
  useShinyjs(),
  # Use CSS and JS files in www/ directory
  tags$head(
    tags$link(rel="stylesheet",type="text/css",href="custom.css"),
    tags$script(src="custom.js")
  ),
  # Tabs for dashboard
  tabItems(
    data_tab,
    model_tab,
    cv_tab,
    results_tab,
    about_tab
  )
)

sidebar <- dashboardSidebar(
 sidebarMenu(
   id = "tabs",
   menuItem("Upload Omics", tabName="data"),
   menuItem("Create Model(s)", tabName="model"),
   menuItem("Train/Validate Model(s)", tabName="validate"),
   menuItem("View CV results", tabName="results"),
   menuItem("About", tabName="about")
 ) 
)

header <- dashboardHeader(
  title = div(tags$img(src="logo_jarquin_uf.png", class="logo-image"), "CHiDO")
)

ui <- dashboardPage(
  title = "CHiDO",
  skin = "black",
  header,
  sidebar,
  body
)