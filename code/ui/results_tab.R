# Created by: Francisco Gonzalez
# Last updated by: Francisco Gonzalez
# Last updated: 06/09/2024

source("code/functions/utils.R")

### UI instructions ----
### Supporting objects ----

metrics_choices <- c("Accuracy (by environment)", "Accuracy (overall)",
                     "RMSE (by environment)", "RMSE (overall)",
                     "Variance components")

results_view_panel <- tabsetPanel(
  tabPanel(
    "Tables",
    div(style="margin-top: 15px;", uiOutput("results_menu")),
    div(style="overflow-y: scroll;", DTOutput("results_table"))
    ),
  tabPanel(
    style="height:800px;",
    "Visualizations",
    div(style="padding: 10px;", selectInput("metrics", "Select a visual:", 
                choices = metrics_choices, selected = metrics_choices[[1]])
        ),
   div(style="align-items:center;", uiOutput("plot_or_message"))
    
  )
)

# Download results

### Main tab ----

results_tab <- tabItem(
  tabName = "results",
  h1("View and Download Results", class="tab-header"),
  fluidRow(
  box(
    class = "upload-box",
    width = 12,
    div(downloadButton("download_results", "Download Results"),
        style="margin-top: 10px; margin-bottom: 10px;"),
    results_view_panel
    )
  )
)