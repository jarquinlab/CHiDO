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
# Last Updated By: Francisco Gonzalez
# Last Updated: 01/10/2025

library(shiny)
library(shinydashboard)
library(DT)
library(shinyjqui)
library(shinyjs)

source("code/functions/utils.R")

### UI instructions ----

add_omics_instructions <- "<strong>Instructions:</strong><br>Select the file(s)
  from your computer to upload and identify the necessary columns to link the file
  to others. For example, for the <strong>phenotypic response file</strong> 
  (left), you must specify at least the following columns: (1) genotypes/lines 
  IDs, (2) environment/plot IDs, and (3) the target trait. Optionally, you can
  specify a sample/inidividual-level ID if there was data collected at such 
  level. For omics uploads (right), you must specify the ID column and ID type
  to link the data to your phenotypic response data. You can modify the 
  <em>Reference Label</em> at your discretion. Make sure to mark the correct 
  <em>Linkage Type</em> before uploading the data. All fields marked with an 
  asterisk (*) are required."

data_view_instructions <- "To view the data you have uploaded for this session,
  scroll down and ensure the file name, label and other inputted values are
  correct. You can click on a column and press <strong>View</strong> to view its
  contents in a separate window, or press <strong>Delete Row(s)</strong> to
  delete one or more omics selected.<br><br>"

### Supporting objects ----

# Dropdown items for data type selection
omic_types_gen <- c("genomic markers", "phenomic markers", "environmental markers",
                "pedigree data", "high-throughput data", "other")

omic_types_cga <- c("genomic markers", "phenomic markers", "environmental markers",
                "pedigree data", "high-throughput data", "other")

# Upload phenotypic response file box
upload_y_box <- box(
  width = 6,
  class = "upload-box",
  fileInput("y_file", "Upload phenotype response (Y) file:", multiple=FALSE),
  tags$div(style="position: relative; top: -20px;",
  radioButtons(inputId="modtype",label="Select model type:",choices=c("Genotype level", "Combining ability"),selected="Genotype level",inline=T)),
  uiOutput("y_panel"),
  tags$hr(class="separator"),
  # Upload button
  div(class="action-btn", div(class="btn-contents",
                              actionButton("y_load", "Upload"))
      )
)

# Omics upload box
upload_omics_box <- box(
  width = 6, 
  class = "upload-box",
  selectInput("data_type", "Select data type:", omic_types_gen, selected="genomic markers"),
  fileInput("file","Choose file to upload:", multiple=FALSE),
  tags$hr(style="margin-top: -20px;"),
  uiOutput("data_panel"),
  checkboxInput("check_data","Check data consistency",value = F),
  tags$hr(class="separator"),
  div(class="action-btn", div(class="btn-contents",
                              actionButton("upload", "Upload"))
      )
)

# Data view panel
data_view_panel <- tabsetPanel(
  tabPanel(
    "Uploaded Data",
    div(style="padding: 5px; float: right;", conditionalPanel(
      condition = "input.table_rows_selected.length > 0",
      actionButton("delete_rows", "Delete Row(s)"))
    ),
    div(style="margin-top: 20px;", DTOutput("omics_table"))
  ),
  tabPanel(
    "View Data",
    selectInput("omics_labels", "Select omic to view:", choices=NULL),
    div(style="margin-top: 20px;", DTOutput("table_to_view"))
  )
)

### Main tab ----

data_tab <- tabItem(
  tabName = "data",
  h1("Upload Data", class="tab-header"),
  text_to_html(add_omics_instructions),
  fluidRow(
    upload_y_box,
    upload_omics_box
  ),
  tags$hr(class="separator"),
  text_to_html(data_view_instructions),
  data_view_panel
)


