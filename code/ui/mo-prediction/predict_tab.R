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
# Last Updated By: Francisco Gonzalez
# Last Updated: 06/10/2024

library(shiny)
library(shinydashboard)
library(DT)
library(shinyjqui)
library(shinyjs)

source("code/functions/utils.R")

### UI instructions ----

add_omics_instructions2 <- "<strong>Instructions:</strong><br>Instructions Placeholder"

data_view_instructions2 <- "Instructions Placeholder<br><br>"

### Supporting objects ----

# Dropdown items for data type selection
omic_types_gen <- c("genomic markers", "phenomic markers", "environmental markers",
                "pedigree data", "high-throughput data", "other")

omic_types_cga <- c("genomic markers","genomic relationship matrix", "phenomic markers", "environmental markers",
                "pedigree data", "high-throughput data", "other")

omic_types_hp<-c("genomic markers","phenomic markers","pedigree data","high-throughput data","other")

# long select model
select_model2 <- box(
  width = 15, 
  class = "upload-box",
  selectInput("predict_model", "Select a model:", choices = NULL, width = "100%"),
)

upload_y_box2 <- box(
  width = 6,
  class = "upload-box",
  panel<-div(fileInput("pred_y_file", "Upload the phenotypes you wish to predict:", multiple=FALSE),
             create_input_box("EID column:", "eid_col2", 1),
             create_input_box("GID column:", "gid_col2", 2),
             create_input_box("UID column:", "uid_col2", 3),
             div(class="action-btn", div(class="btn-contents",
                                         actionButton("pred_y_load", "Upload")))
  )
)

# Omics upload box
upload_omics_box2 <- box(
  width = 6, 
  class = "upload-box",
  selectInput("omics_new_features", "Select data type:", choices = NULL, selected=NULL),
  fileInput("feature_file","Choose file to upload:", multiple=FALSE),
  tags$hr(style="margin-top: -20px;"),
  uiOutput("data_panel2"),
  checkboxInput("check_data2","Check data consistency",value = F),
  tags$hr(class="separator"),
  div(class="action-btn", div(class="btn-contents",
                              actionButton("feature_upload_btn", "Upload"))
  ),
  div(style = "text-align: center;", # will only enable if all the omics are provided for predictions. 
      actionButton("run_predictions", "Run Predictions", class = "btn-secondary", disabled = TRUE))
)

data_view_panel2 <- tabsetPanel(
  tabPanel(
    "Uploaded Data",
    div(style="padding: 5px; float: right;", conditionalPanel(
      condition = "input.table_rows_selected.length > 0",
      actionButton("delete_rows2", "Delete Row(s)"))
    ),
    div(style="margin-top: 20px;", DTOutput("predict_prompt"))
  ),
  tabPanel(
    "View Data",
    selectInput("omics_labels2", "Select omic to view:", choices=NULL),
    div(style="margin-top: 20px;", DTOutput("table_to_view2"))
  )
)
### Main tab ----

predict_tab <- tabItem(
  tabName = "predict",
  h1("Predict New Phenotypes", class="tab-header"),
  text_to_html(add_omics_instructions2),
  select_model2,
  fluidRow(
    upload_y_box2,
    upload_omics_box2
  ),
  tags$hr(class="separator"),
  text_to_html(data_view_instructions2),
  data_view_panel2
)

