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
# Last Updated: 06/09/2024

source("code/functions/utils.R")

### UI instructions ----

cv_instructions <- "With the model(s) you saved, you can begin performing
  cross-validations to train/test them. Customize the hyperparameters and data
  pre-processing steps as necessary. The settings do not reset once a cross-validation
  is completed, so please re-adjust if you want to use different configurations."

### Supporting objects ----

params_box <- box(
  width = 6,
  class = "upload-box",
  h3("Hyperparameter Tuning", class="tab-header"),
  div(style="display:flex; align-items:center; margin-top:-20px; padding:5px;",
      tags$label("Number of iterations:", style="width:200px;"),
      div(style="margin-left:auto; min-width:60px; max-width: 80px;",
          numericInput("n_iter_input", "",value=5000, max=50000))
      ),
  sliderInput("n_iter_slider", label=NULL, min=0, max=50000, value=5000, step=1000),
  div(style="display:flex; align-items:center; margin-top: -10px; padding:5px;",
      tags$label("Burn-in rate:", style="width:200px;"),
      div(style="margin-left:auto; min-width:60px; max-width:80px;",
          numericInput("burn_in_input","",value=500,max=5000))
      ),
  sliderInput("burn_in_slider", label=NULL, min=0, max=5000, value=500, step=100),
  tags$hr(),
  div(style="margin-top:20px;", checkboxInput("user_set_seed", label="Set seed?", value=FALSE),
      conditionalPanel(
        condition="input.user_set_seed==true",
        div(style="display:flex; align-items:center; margin-top:-30px; padding:3px;",
            tags$label("Seed value:", style="width:200px;"),
            div(style="margin-left:auto; min-width:60px; max-width:80px;",
                numericInput("set_seed","", value=1)))))
)

preprocess_box <- box(
  width = 6,
  class = "upload-box",
  h3("Data Pre-processing", class="tab-header"),
  div(style="margin-top:20px;", helpText("Specify number of folds for cross-validation")),
  div(style="display:flex; align-items:center; margin-top:-20px; padding:5px;",
      tags$label("Folds:", style="width:200px;"),
      div(style="margin-left:auto; min-width:60px; max-width:80px;",
          numericInput("folds","", value=5, min=1, max=99))),
  helpText("Specify how much missing data can be tolerated in your omics"),
  div(style="display:flex; align-items:center; margin-top:-20px; padding:5px;",
      tags$label("NaN threshold (%):", style="width:200px;"),
      div(style="margin-left:auto; min-width:60px; max-width: 80px;",
          numericInput("nan_freq","",value=0, min=0, max=99))),
  helpText("Select any data transformations needed prior to model training"),
  column(6, 
         checkboxInput("std","Standardizing", value=FALSE),
         checkboxInput("ctr","Centering", value=FALSE)
  )
)

cv_box <- fluidRow(
  tags$head(tags$script(src = "custom.js")),
  box(
    width = 12,
    class = "upload-box",
    h3("Cross-validation", class="tab-header"),
    column(
      12,
      selectInput("saved_models", "Select a model:", choices=NULL, width="100%"),
      div(style="padding: 5px;", fluidRow(
        column(4,
               tags$b("Select desired CV scheme(s):"),
               checkboxInput("cv1", "CV1", value=FALSE),
               checkboxInput("cv2", "CV2", value=FALSE),
               checkboxInput("cv0", "CV0", value=FALSE),
               checkboxInput("cv00", "CV00", value=FALSE),
               tags$br(),
               actionButton("run_cv", "Run")
               ),
        column(8, verbatimTextOutput("cv_logs", placeholder=TRUE))
      )
      )
  )
  )
)

### Main tab ----

cv_tab <- tabItem(
  tabName = "validate",
  h1("Train & Validate", class="tab-header"),
  text_to_html(cv_instructions),
  # Tuning
  fluidRow(
    width = 12,
    params_box,
    preprocess_box
  ),
  # Cross validations
  cv_box
)