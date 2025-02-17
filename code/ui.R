# CHiDO is a no-code platform to integrate multi-omics data to build, train and test
# linear mixed models for identifying candidates for desired GxE interactions.
#
# Copyright (C) 2024 Francisco Gonzalez, Diego Jarquin, and Julian Garcia
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# Created by: Francisco Gonzalez
# Last updated by: Francisco Gonzalez
# Last updated: 06/09/2024

source("code/ui/mo-prediction/data_tab.R")
source("code/ui/mo-prediction/model_tab.R")
source("code/ui/mo-prediction/cv_tab.R")
source("code/ui/mo-prediction/results_tab.R")
source("code/ui/gei-analysis/gei_data_tab.R")
source("code/ui/gei-analysis/runana_tab.R")
source("code/ui/gei-analysis/gei_results_tab.R")
source("code/ui/about_tab.R")

body <- dashboardBody(
  # Enable JavaScript
  useShinyjs(),
  # Use CSS and JS files in www/ directory
  tags$head(
    tags$link(rel="stylesheet",type="text/css",href="custom.css"),
    tags$script(src="custom2.js")
  ),
  # Tabs for dashboard
  tabItems(
    data_tab,
    model_tab,
    cv_tab,
    results_tab,
    gei_data_tab,
    runana_tab,
    gei_results_tab,
    about_tab
  )
)

sidebar <- dashboardSidebar(
  sidebarMenu(
    id = "tabs",
    menuItem("Multi-omic Prediction",tabName="prediction", startExpanded = TRUE,
             menuSubItem("Upload Data", tabName="data", icon = icon("arrow-right")),
             menuSubItem("Create Model(s)", tabName="model", icon = icon("arrow-right")),
             menuSubItem("Train/Validate Model(s)", tabName="validate", icon=icon("arrow-right")),
             menuSubItem("View CV results", tabName="results", icon = icon("arrow-right"))),
   
    menuItem("GEI Analysis", tabName="gei",startExpanded = FALSE,
      menuSubItem("Upload Data", tabName="gei_data", icon = icon("arrow-right")),
      menuSubItem("Run model", tabName="run_gei_tab", icon = icon("arrow-right")),
      menuSubItem("View results", tabName="gei_results", icon = icon("arrow-right"))),
    
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