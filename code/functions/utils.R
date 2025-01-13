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

library(shiny)

# Convert text into a Shiny compliant format
text_to_html <- function(text) {
  return(HTML(paste(text,"<br><br>")))
}

# Numeric input boxes to define ID columns in Add Omics tab
create_input_box <- function(label_text, input_id, default_value = NULL) {
  div(style = "display: flex; align-items: center; margin-top: -30px; padding: 3px;", 
      tags$label(label_text, style = "width: 200px;"),
      div(style = "margin-left: auto; min-width: 60px; max-width: 100px;", 
          numericInput(input_id, "", value = default_value, min = 1, max = 999))
  )
}

create_panel <- function(data_type, label) {
  id_col <- tolower(label)
  
  panel <- conditionalPanel(
    condition = paste0("input.data_type == '", data_type, "'"),
    # Label
    div(style="display: flex; align-items: center; margin-top: -20px; padding: 2px;", 
        tags$label("Reference Label:", style="width:200px;"),
        div(style="margin-left: auto; min-width: 60px; max-width: 70px;",
            textInput("label", "", label))
    ),
    # ID
    div(style="display: flex; align-items: center; margin-top: -20px; padding: 2px;", 
        tags$label("Linkage column (relation to Y):", style="width: 200px;"),
        div(style="margin-left: auto; min-width: 60px; max-width: 100px;", 
            numericInput("id_col", "", value = 1, min = 1, max = 999))
    ),
    selectInput("linkage_type", "Linkage type", 
                c("Environment ID", "Genotype / Line ID", "Compound ID / UID")
                )
  )
  return(panel)
}

get_num_fields_from_list <- function(config) {
  return(names(config)[sapply(config, function(x) is.numeric(x) && 
                                !is.null(x) && 
                                !is.na(x))]
  )
}

logging <- function(session, message="", t = 0.5) {
  if (message == "") { 
    message == "clearLogs"
  } else {
    message = paste0(message,"\n")
  }
  
  session$sendCustomMessage("updateLog", message)
  Sys.sleep(t)
}