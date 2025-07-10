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
pred_instructions <- tagList(
    p("CHiDO builds and trains ", strong("linear models"), " using a ", strong("supervised"), 
      " machine learning framework. Understanding the two main phases will help you get the most out of this tool:"),
    
    tags$ul(
      tags$li(strong("Training Phase (TRN):"), " During this stage, models learn from both your phenotypic data and omics information. Cross-validation is used here to help us evaluate how well the models perform."),
      tags$li(strong("Prediction Phase (PRD):"), " Models apply what they've learned during the training phase to new sets of omics data, allowing us to predict unknown phenotypes.")
    ),
    
    p("This tab is specifically designed to guide you through the ", strong("Prediction (PRD) phase"),"."),
  )

interpretation_guide_ui <- tagList(
  h2("Results"),
  
  p("The ", strong("Results tab"), " provides access to two CSV files for each of the models:"),
  tags$ul(
    tags$li(code("prediction_prompt.csv"), ": contains the predictions (yHat) for the provided prompt. The column ", strong("Predictable"),
            "indicates whether the phenotype could be predicted using the provided omics."),
    tags$li(code("prediction_main.csv"), " contains all the relevant information in the model, including both TRN and PRD sets and is used for ", strong("visualizations"))
  ),
  tags$hr(),
  
  p("The ", strong("Visualization tab"), " is designed to help you understand the predicted phenotypes:
    "),

  h3("Top Genotypes in Known Environments (CV1-like)"),
  p("This plot highlights genotypes that fall into the extreme ", strong("p% of predicted values"), 
    ", allowing the detection of promising candidates that are part of the prediction (PRD) set and therefore were never tested in the TRN environments."),
  
  h3("Top Genotypes in PRD Environment Across TRN Environments (CV0-like)"),
  p("This plot displays the performance of the ", strong("top k lines"), " (identified as best in a specific ", strong("new environment"), ") across ", 
    strong("all TRN environments"), "providing a measure of similarity between environments in the same TPE."),
  tags$hr()
)

### Supporting objects ----

# Dropdown items for data type selection
omic_types_gen <- c("genomic markers", "phenomic markers", "environmental markers",
                "pedigree data", "high-throughput data", "other")

omic_types_cga <- c("genomic markers","genomic relationship matrix", "phenomic markers", "environmental markers",
                "pedigree data", "high-throughput data", "other")

omic_types_hp<-c("genomic markers","phenomic markers","pedigree data","high-throughput data","other")

select_model2 <- box(
  width = 15,
  class = "upload-box",
  title = "Select Prediction Model", 
  
  p("Choose the", strong("pre-trained model"), "you wish to use for your predictions from the list below. This selection will determine the required input features."), 
  
  selectInput("predict_model", "Select a model:", choices = NULL, width = "100%"),
)

# --- UI definition for upload_y_box2 ---
upload_y_box2 <- box(
  width = 6,
  class = "upload-box", 
  title = "Prepare New Phenotype Predictions", 
  
  p(strong("Upload a CSV file"), " containing the identifiers for the phenotypes you wish to predict. This file should include columns for ", strong("Environment (EID), Genotype / Line (GID),"), " and a ", strong("Compound / Unit (UID).")),
  tags$hr(), 
  
    h4("1. Select Your Prediction Prompt File"), 
    fileInput(
      "pred_y_file",
      label = "Choose file to upload:", 
      multiple = FALSE
    ),
    h4("2. Map ID Columns in Your File"), 
    p("Specify the exact column indices in your uploaded CSV file that correspond to the required identifiers:"),
    create_input_box("Environment ID (EID):", "eid_col2", 1), 
    create_input_box("Genotype / Line ID (GID):", "gid_col2", 2),
    create_input_box("Compound / Unit ID (UID):", "uid_col2", 3),
  tags$hr(), 
    div(class="action-btn", 
        div(class="btn-contents",
            actionButton(
              "pred_y_load",
              "Upload"#,
              #icon = icon("upload") 
            )
        )
    )
  )

# Omics upload box with improved UX and styling
upload_omics_box2 <- box(
  width = 6,
  class = "upload-box",
  title = "Upload Omics Features for Prediction", 
  
  p(strong("Upload a CSV file"), 
    " containing the new features that will be used to predict phenotypes. Each file should contain the necessary data for a particular omics type."),
  tags$hr(), # Separator after introduction
  
  h4("1. Select Omics Data Type"), 
  selectInput("omics_new_features", "Select data type:", choices = NULL, selected=NULL),
  
  h4("2. Upload Your Feature File"), 
  fileInput("feature_file", "Choose file to upload:", multiple=FALSE),
  tags$hr(),
  uiOutput("data_panel2"),
  
  # Upload button with consistent styling
  div(class="action-btn",
      div(class="btn-contents",
          actionButton("feature_upload_btn", "Upload")
      )
  ),
  
  # Run Predictions button (disabled by default)
  div(style = "text-align: center; margin-top: 20px;", 
      actionButton("run_predictions", "Run Predictions", class = "btn-secondary", disabled = TRUE, 
                   icon("bolt")))
)

data_view_4_viz <- box(
  width = 6, 
  class = "upload-box",
  div(style = "padding: 10px;",
  uiOutput("predict_viz_file"),
  tags$hr(style="margin-top: -20px;"),
  uiOutput("predict_cv_options")
  )
)

data_view_panel2 <- tabsetPanel(
  # tabPanel(
  #   "Data Sources",
  #   div(style="padding: 5px; float: right;", conditionalPanel(
  #     condition = "input.table_rows_selected.length > 0",
  #     actionButton("delete_rows2", "Delete Row(s)"))
  #   ),
  #   div(style="margin-top: 20px;", DTOutput("predict_prompt"))
  # ),
  # tabPanel(
  #   "View Data",
  #   selectInput("omics_labels2", "Select omic to view:", choices=NULL),
  #   div(style="margin-top: 20px;", DTOutput("table_to_view2"))
  # ),
  tabPanel(
    "Results",
    div(style="margin-top: 20px;", uiOutput("predict_results_menu")),
    div(style="overflow-y: scroll;", DTOutput("predict_results_table"))
  ),
  tabPanel(
    "Visualization", 
    fluidRow(
      data_view_4_viz,
      box(
        width = 6,
        height = "620px",
        title = '',
        div(plotOutput("predict_plot", height = "550px")))
    )
  )
)

### Main tab ----
predict_tab <- tabItem(
  tabName = "predict",
  h1("Predict New Phenotypes", class="tab-header"),
  pred_instructions,
  select_model2,
  fluidRow(
    upload_y_box2,
    upload_omics_box2
  ),
  tags$hr(class="separator"),
  interpretation_guide_ui,
  data_view_panel2
)

