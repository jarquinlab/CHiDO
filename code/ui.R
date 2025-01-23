# CHiDO is a no-code platform to integrate multi-omics data to build, train and test
# linear mixed models for identifying candidates for desired GxE interactions.
#
# Copyright (C) 2025 Francisco Gonzalez, Diego Jarquin, and Julian Garcia, Vitor Sagae
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

source("code/ui/data_tab.R")
source("code/ui/model_tab.R")
source("code/ui/cv_tab.R")
source("code/ui/results_tab.R")
source("code/ui/about_tab.R")

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

ui <- fluidPage(
  tags$head(tags$link(rel = "icon", type = "image/png", href = "favicon-32x32.png")),
  dashboardPage(
    title = "CHiDO",
    skin = "black",
    header,
    sidebar,
    body
  )
)


