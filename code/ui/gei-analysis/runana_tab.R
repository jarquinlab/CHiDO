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

# Created by: Vitor Sagae
# Last Updated By: Vitor Sagae
# Last updated: 02/27/2025

source("code/functions/utils.R")

### UI instructions ----

ruana_instructions <- "With the data you have uploaded, you can begin performing
  the following analyses: Bayesian AMMI (Crossa et al., 2011) and Bayesian weighted average of absolute scores including or no the response variable (BWAAS and BWAASY) (Olivoto et al., 2019; Nascimento et al., 2025). 
  Customize the hyperparameters as necessary. The settings do not reset once the analyses
  are completed, so please re-adjust if you want to use different configurations."

### Supporting objects ----

gei_params_box <- box(
  width = 6,
  class = "params-box",
  h3("Hyperparameter Tuning", class="tab-header"),
  div(style="display:flex; align-items:center; margin-top:-20px; padding:5px;",
      tags$label("Number of iterations:", style="width:200px;"),
      div(style="margin-left:auto; min-width:60px; max-width: 80px;",
          numericInput("gei_n_iter_input", "",value=1000, max=100000))
      ),
  sliderInput("gei_n_iter_slider", label=NULL, min=0, max=100000, value=1000, step=1000),
  div(style="display:flex; align-items:center; margin-top: -10px; padding:5px;",
      tags$label("Burn-in rate:", style="width:200px;"),
      div(style="margin-left:auto; min-width:60px; max-width:80px;",
          numericInput("gei_burn_in_input","",value=100,max=10000))
      ),
  sliderInput("gei_burn_in_slider", label=NULL, min=0, max=10000, value=100, step=100),
  div(style="display:flex; align-items:center; margin-top: -10px; padding:5px;",
      tags$label("Thin rate:", style="width:200px;"),
      div(style="margin-left:auto; min-width:60px; max-width:80px;",
          numericInput("gei_thin_input","",value=5,max=100))
  ),
  sliderInput("gei_thin_slider", label=NULL, min=0, max=100, value=5, step=5),
  tags$hr(),
  fluidRow(
    div(
      style = "display: flex; justify-content: space-between; align-items: center",
      
      div(
        style = "flex: 1; padding-right: 30px; margin-left: 30px;",
        checkboxInput("gei_user_set_seed", label = "Set seed?", value = FALSE),
        conditionalPanel(
          condition = "input.gei_user_set_seed==true",
          div(
            style = "display: flex; align-items: center; margin-top: -30px; padding: 3px;",
            tags$label("Seed value:", style = "margin-right: 10px;"),
            div(
              style = "min-width:60px; max-width:70px;",
              numericInput("gei_set_seed", "", value = 1)
            )))),
      
      div(
        style = "flex: 1; padding-left: 10px;",
        div(
          style = "display: flex; align-items: center; margin-top: 10px; padding: 3px;",
          tags$label("Number of chains:", style = "margin-right: 10px;"),
          div(
            style = "min-width:60px; max-width:80px;",
            numericInput("gei_chains", "", value = 2, min = 2, max = 50),
            tags$script("$('#gei_chains').prop('disabled',true);")
          ))))
  )
)

analysis_box <- box(
  width = 6,
  class = "analysis-box",
  h3("Analyses", class="tab-header"),
  helpText("Select analyses to be performed under bayesian approach"),
  column(6, 
         checkboxInput("bammi","AMMI model", value=FALSE),
         checkboxInput("bstab","STABDIST index", value=FALSE),
         checkboxInput("bwaas","WAAS index", value=FALSE),
         checkboxInput("bwaasy","WAASY index", value=FALSE),
         uiOutput("weigh_panel"),
         actionButton("run_gei", "Run")
  )
)

gei_box <- fluidRow(
  tags$head(tags$script(src = "custom2.js")),
  box(
    width = 12,
    class = "log-box",
    h3("Processing", class="tab-header"),
    column(
      12,
      div(style="padding: 5px;", fluidRow(
        
        column(12, verbatimTextOutput("gei_logs", placeholder=TRUE))
      )
      )
  )
  )
)

### Main tab ----

runana_tab <- tabItem(
  tabName = "run_gei_tab",
  h1("Run analysis", class="tab-header"),
  text_to_html(ruana_instructions),
  # Tuning
  fluidRow(
    width = 12,
    gei_params_box,
    analysis_box
  ),
  # Processing analysis
  gei_box
)