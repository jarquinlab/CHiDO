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

# Created by: Vitor Sagae
# Last Updated By: Vitor Sagae
# Last Updated: 08/01/2025

library(shiny)
library(shinydashboard)
library(DT)
library(shinyjqui)
library(shinyjs)

source("code/functions/utils.R")

### UI instructions ----

add_gei_instructions <- "<strong>Instructions:</strong><br>Select the file
  from your computer to upload and identify the necessary columns. For example, for the <strong>phenotypic response file</strong> 
  , you must specify at least the following columns: (1) genotypes/lines 
  IDs, (2) environment/plot IDs, and (3) the target trait. Optionally, you can
  specify a sample/individual-level ID if there was data collected at such 
  level. All fields marked with an asterisk (*) are required."

data_view_instructions <- "To view the data you have uploaded for this session,
  scroll down and ensure the file name, label and other inputted values are
  correct. You can click on a column and press <strong>View</strong> to view its
  contents in a separate window, or press <strong>Delete Row(s)</strong> to
  delete one or more omics selected.<br><br>"

### Supporting objects ----


# Upload phenotypic response file box
upload_y_box_gei <- box(
  width=8,
  class = "gei-upload-box",
  fileInput("gei_y_file", "Upload phenotype response (Y) file:", multiple=FALSE),
  tags$hr(class="separator"),
  uiOutput("gei_y_panel"),
  # Upload button
  div(class="action-btn", div(class="btn-contents",
                              actionButton("gei_y_load", "Upload"))
      )
)

# Data view panel
data_view_panel_gei <- tabsetPanel(
  tabPanel(
    "Uploaded Data",
    div(style="padding: 5px; float: right;", conditionalPanel(
      condition = "input.table_rows_selected.length > 0",
      actionButton("gei_delete_rows", "Delete Row(s)"))
    ),
    div(style="margin-top: 20px;", DTOutput("gei_table"))
  ),
  tabPanel(
    "View Data",
    selectInput("gei_labels", "Select column to view:", choices=NULL),
    div(style="margin-top: 20px;", DTOutput("gei_table_to_view"))
  )
)

### Main tab ----

gei_data_tab <- tabItem(
  tabName = "gei_data",
  h1("Upload Data", class="tab-header"),
  text_to_html(add_gei_instructions),
  fluidRow(
    upload_y_box_gei
  ),
  tags$hr(class="separator"),
  text_to_html(data_view_instructions),
  data_view_panel_gei
)


