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

model_instructions <- "<strong>Available Omics:</strong> These are the uploaded
  data sources you can use to create linear mixed models. To add an omic as a 
  model term, drag the label to the <strong>Terms to Add</strong> box. You can add
  them as individual terms or place multiple in the box to create interaction terms.
  <br><br>
  <strong>Term(s) to Add:</strong> Add the terms to the model equation by pressing
  'Add' or discard them by pressing 'Clear'.<br><br>
  <strong>Model Preview:</strong> This is the current model formula. Press 'Save'
   to save the model using or press 'Clear' to remove it entirely. The 'Previous' 
  button removes only the last term added to the model equation."


### Supporting objects ----

all_omics_box <- column(
  6,
  h4("Available Omics", style="margin-top:0px; margin-bottom:0px;"),
  div(style="margin-bottom: 46px;", htmlOutput("omics_grid"))
)

terms_to_add_box <- column(
  6,
  h4("Term(s) to Add", style="margin-bottom:0px; margin-top:0px;"),
  orderInput("added_terms", "", items=NULL, placeholder="Drag labels here"),
  div(style="margin-top:10px;", actionButton("clear_terms","Clear"),
      actionButton("add_terms","Add"))
)

model_prev_box <- box(
  width = 12,
  class = "upload-box",
  h4("Model Preview", class="tab-header"),
  fluidRow(
    column(8, textInput("typed_formula", "Already know the formula you want? Type it here (e.g., E+G+E*G)", value="")),
    column(4, div(style="margin-top:25px;", actionButton("add_typed_formula", "Add Formula")))
  ),
  tags$hr(),
  div(style="padding: 10px; 1px solid #ddd; text-aling:center; margin-bottom:58px;
      margin-top: 62px;", uiOutput("model_formula")),
  actionButton("clear_model", "Clear"),
  actionButton("prev_model","Previous"),
  actionButton("save_model","Save")
)

### Main tab ----

model_tab <- tabItem(
  tabName = "model",
  h1("Model Assembly", class="tab-header"),
  text_to_html(model_instructions),
  fluidRow(
    box(
      width = 12,
      class = "upload-box",
      all_omics_box,
      terms_to_add_box
    )
  ),
  tags$br(),
  fluidRow(
    model_prev_box
  )
)