library(shiny)
library(shinydashboard)
library(shinyjqui)
library(shinyjs)
library(DT)
library(ggplot2)
library(gridExtra)
library(BGLR)
library(dplyr)

source("code/ui.R")
source("code/server.R")

shinyApp(ui,server)