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

source("code/functions/utils.R")

### UI instructions ----
### Supporting objects ----

metrics_choices_gei<-list(
  "bammi"=c("Trace mu","Trace tau","Trace Lambda PC1","Trace Lambda PC2","AMMI Plot"),
  "bstab"=c("Stabdist Plot"),
  "bwaas"=c("Waas Plot"),
  "bwaasy"=c("Waasy Plot")
)

gei_results_view_panel <- tabsetPanel(
  tabPanel(
    "Tables",
    div(style="margin-top: 15px;", uiOutput("gei_results_menu")),
    div(style="overflow-y: scroll;", DTOutput("gei_results_table"))
    ),
  tabPanel(
    style="height:800px;",
    "Visualizations",
    div(style="padding: 10px;", selectInput("gei_metrics", "Select a visual:", 
                choices = NULL)
        ),
   div(style="align-items:center;", uiOutput("plot_or_message_gei"))
    
  )
)

# Download results

### Main tab ----

gei_results_tab <- tabItem(
  tabName = "gei_results",
  h1("View and Download Results", class="tab-header"),
  fluidRow(
  box(
    class = "upload-box",
    width = 12,
    div(downloadButton("download_results_gei", "Download Results"),
        style="margin-top: 10px; margin-bottom: 10px;"),
    gei_results_view_panel
    )
  )
)