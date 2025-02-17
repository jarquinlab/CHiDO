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
# Last updated: 06/10/2024

### Libraries ----

library(shiny)
library(shinyjqui)
library(dplyr)
library(ggplot2)
library(BGLR)
library(gridExtra)

### Loading source files ----

# Import functions for server
source("code/functions/utils.R")
source("code/functions/data.R")
source("code/functions/model.R")
source("code/functions/metrics.R")
source("code/functions/gei_model.R")
source("code/functions/gei_metrics.R")

create_temp_dir <- function() {
  tmpdir <- tempdir()
  if(!dir.exists(tmpdir)) { dir.create(tmpdir, recursive=TRUE) }
  return(tmpdir)
}

server <- function(input, output, session) {
  
  ### Global ----
  
  # Set maximm file size to 50MB # I Changed to 100 mb to test
  options(shiny.maxRequestSize = 100 * 1024 ^ 2)
  
  # Create temporary directory for session files
  tmpdir <- create_temp_dir()
  
  # Reactive variables to maintain sesssion
  phen_data <- reactiveVal(list())
  data_sources <- reactiveVal(list())
  omics_labels <- reactiveVal(list())
  data_types <- reactiveVal(list())
  file_counter <- reactiveVal(0)
  model_terms <- reactiveVal(character(0))
  saved_models <- reactiveVal(list())
  logs <- reactiveVal("")
  table_data <- reactiveVal(data.frame(
    Label = character(),
    File = character(),
    Type = character(),
    ID = numeric(),
    Linkage = character(),
    stringsAsFactors = FALSE
  ))
  
  ### Data tab ----
  
  # Generate a dynamic upload panel based on model type selected
  output$y_panel<-renderUI({
    if (input$modtype%in%"Genotype level"){
      panel<-div(create_input_box('*Genotype ID column:', "gid_col", 1),
                 create_input_box("*Environment ID column:", "eid_col", 2),
                 create_input_box("Individual ID (UID) column:", "uid_col"),
                 create_input_box("*Target trait column:", "trait_col", 3))
    }else{
      
      panel<-div(create_input_box('*Genotype ID column:', "gid_col", 1),
                 create_input_box('*Parent Group 1 ID column:', "g1id_col", 2),
                 create_input_box('*Parent Group 2 ID column:', "g2id_col", 3),
                 create_input_box("*Environment ID column:", "eid_col", 4),
                 create_input_box("Individual ID (UID) column:", "uid_col"),
                 create_input_box("*Target trait column:", "trait_col", 5))
    }
  })
  
  # Generate a dynamic upload panel based on omic type selected
  output$data_panel <- renderUI({
    if (input$modtype=="Genotype level"){
      if (input$data_type %in% c("high-throughput data", "other")) {
        panel <- create_panel(input$data_type, input$label)
      } 
      else {
        
        label_id <- ifelse(input$data_type == "environmental markers",
                           "*Environment ID",
                           "*Genotype / Line ID")
        
        panel <- div(
          # Label
          div(
            style = "display: flex; align-items: center; margin-top: -27px; padding: 4px;",
            tags$label('*Reference label:', style = "width:200px;"),
            div(style = "margin-left: auto; min-width: 60px; max-width: 70px;",
                textInput("label", "", ""))
          ),
          # Genome or Env ID
          div(
            style = "display: flex; align-items: center; margin-top: -27px; padding: 4px;",
            tags$label(label_id, style = "width: 200px;"),
            div(
              style = "margin-left: auto; min-width: 60px; max-width: 100px;",
              numericInput("id_col", "", value = 1, min = 1, max = 999)
            )
          )
        )
      }
    }
    else{
      #if (input$data_type %in% c("high-throughput data", "other")) {
      id_col <- tolower(input$label)
      
      panel <- conditionalPanel(
        condition = paste0("input.data_type == '", input$data_type, "'"),
        # Label
        div(style="display: flex; align-items: center; margin-top: -20px; padding: 2px;", 
            tags$label("Reference Label:", style="width:200px;"),
            div(style="margin-left: auto; min-width: 60px; max-width: 70px;",
                textInput("label", "", input$label))
        ),
        # ID
        div(style="display: flex; align-items: center; margin-top: -20px; padding: 2px;", 
            tags$label("Linkage column (relation to Y):", style="width: 200px;"),
            div(style="margin-left: auto; min-width: 60px; max-width: 100px;", 
                numericInput("id_col", "", value = 1, min = 1, max = 999))
        ),
        selectInput("linkage_type", "Linkage type", 
                    c("Environment ID", "Parent Group 1 ID", "Parent Group 2 ID","Parent Groups ID", "Compound ID / UID")
        )
      )
      #} 
    }
  })
  
  # Switch omic options according to model type
  observeEvent(input$modtype,{
    omic<-if (input$modtype=="Genotype level"){
      omic_types_gen
    }else{
      omic_types_cga
    }
    updateSelectInput(session,"data_type",choices=omic)
  })
  
  # Switch panel defaults based on data type
  observe({
    switch(
      input$data_type,
      "genomic markers" = {
        updateTextInput(session, "label", value="G")
      },
      "phenomic markers" = {
        updateTextInput(session, "label", value="P")
      },
      "pedigree data" = {
        updateTextInput(session, "label", value="A")
      },
      "environmental markers" = {
        updateTextInput(session, "label", value="W")
      },
      "high-throughput data" = {
        updateTextInput(session, "label", value="F")
      },
      "other" = {
        updateTextInput(session, "label", value="O")
      }
    )
  })
  
  # Loading Y file after upload button is pressed
  observeEvent(input$y_load, {
    # Create omic metadata object
    if(input$modtype=="Genotype level"){
      config <- create_omic_config(input$y_file, "Y", "phenotypic",input$modtype, input$trait_col, 
                                   gid_col = input$gid_col, eid_col = input$eid_col,
                                   uid_col = input$uid_col)
    }else{
      config <- create_omic_config(input$y_file, "Y", "phenotypic",input$modtype, input$trait_col, 
                                   gid_col = input$gid_col,g1id_col = input$g1id_col, g2id_col = input$g2id_col, eid_col = input$eid_col,
                                   uid_col = input$uid_col)
    }
    
    if(is.null(config)) {
      showNotification("No file was uploaded, try again.", type="error")
      return()
    }
    
    verified <- verify_omic_object(config, phen = TRUE,modtype=input$modtype)
    
    if(verified) {
      
      # Add data to metadata config 
      config$data <- load_data(input$y_file$datapath)
      
      # Grab values from reactive objects
      phen_data(config)
      curr_data <- data_sources()
      
      if (input$modtype=="Genotype level"){
        for (omic in c("E","L")) {
          
          if(omic == "E") {
            data_type = "Environment IDs"
            link = "Environment ID"
            id_col = input$eid_col
            uid_col=input$uid_col
          } else {
            data_type = "Line IDs"
            link = "Genotype / Line ID"
            id_col = input$gid_col
            uid_col=input$uid_col
          }
          
          # Get config for ID column
          curr_omic <- create_omic_config(input$y_file, omic, data_type,input$modtype,id_col,
                                          link = link, uid_col=uid_col)
          # Load data as data.frame for that respective column
          curr_omic$data <- config$data[, id_col, drop=FALSE]
          
          if(!is.na(input$uid_col)){
          curr_omic$uid<- config$data[, uid_col,drop=FALSE]
          }
          
          curr_data[[curr_omic$label]] <- curr_omic
          
          # Add omic to preview table
          curr_table <- table_data()
          updated_table <- update_table(curr_table, curr_data[[curr_omic$label]])
          table_data(updated_table)
        }}else{
          for (omic in c("E","L","L1","L2")) {
            
            if(omic == "E") {
              data_type = "Environment IDs"
              link = "Environment ID"
              id_col = input$eid_col
              uid_col=input$uid_col
            } 
            if(omic == "L") {
              data_type = "Line IDs"
              link = "Genotype / Line ID"
              id_col = input$gid_col
              uid_col=input$uid_col
            } 
            if (omic == "L1"){
              data_type = "Line Group 1 IDs"
              link = "Parent Group 1 ID"
              id_col = input$g1id_col
              uid_col=input$uid_col
            } 
            if (omic == "L2") {
              data_type = "Line Group 2 IDs"
              link = "Parent Group 2 ID"
              id_col = input$g2id_col
              uid_col=input$uid_col
            }
            
            
            # Get config for ID column
            curr_omic <- create_omic_config(input$y_file, omic, data_type,input$modtype,id_col,
                                            link = link, uid_col=uid_col)
            
            # Load data as data.frame for that respective column
            curr_omic$data <- config$data[, id_col, drop=FALSE]
            
            if(!is.na(uid_col)){
              curr_omic$uid<- config$data[, uid_col,drop=FALSE]
            }
            
            curr_data[[curr_omic$label]] <- curr_omic
            
            # Add omic to preview table
            curr_table <- table_data()
            updated_table <- update_table(curr_table, curr_data[[curr_omic$label]])
            table_data(updated_table)
          }}
      
      # Save data added for E and L
      data_sources(curr_data)
      # Update avalable omics in model assembly page
      updateSelectInput(session, "omics_labels", choices = names(data_sources()))
      showNotification("Data was successfully uploaded!", type="message")
      
    } else {
      showNotification("Error: check the data entered in the upload box.", type="error")
      return()
    }
  })
  
  
  # Add omics data to session
  observeEvent(input$upload, {
    if (input$modtype=="Genotype level"){
      
      # Create omic metadata object
      config <- create_omic_config(input$file, input$label, input$data_type,input$modtype, 
                                   input$id_col, link = input$linkage_type)
      
      if(is.null(config)) {
        showNotification("No file was uploaded, try again.", type="error")
        return()
      }
      
      if(input$check_data & is.null(input$y_file)){
        #showNotification("To check data consistency you need to upload phenotype response file (Y).", type="error")
        return()
      }
      
      verified <- verify_omic_object(config)
      
      if(verified) {
        
        # Add data to metadata config 
        config$data <- load_data(input$file$datapath)
        
        # Grab values from reactive objects
        curr_data <- data_sources()
        curr_types <- data_types()
        curr_table <- table_data()
        
        # New omic type, register and increment file counter
        if (!config$type %in% curr_types) {
          file_counter(file_counter() + 1)
          data_types(c(curr_types, config$type))
        }
        
        # Check if the label already exists and remove the corresponding row
        existing_label_idx <- which(curr_table$Label == config$label)
        if (length(existing_label_idx) > 0) {
          curr_table <- curr_table[-existing_label_idx, ]
        }
        
        # Add omic with data to available data sources
        curr_data[[config$label]] <- config
        data_sources(curr_data)
        
        # Add omic to preview table
        updated_table <- update_table(curr_table, curr_data[[config$label]])
        table_data(updated_table)
        
        # Update avalable omics in model assembly page
        updateSelectInput(session, "omics_labels", choices = names(data_sources()))
        showNotification("Data was successfully uploaded!", type="message")
        
      } else {
        showNotification("Error: check the data entered in the upload box.", type="error")
        return()
      }
    }else{
      if (!input$linkage_type=="Parent Groups ID"){
        # Create omic metadata object
        config <- create_omic_config(input$file, input$label, input$data_type,input$modtype, 
                                     input$id_col, link = input$linkage_type)
        if(is.null(config)) {
          showNotification("No file was uploaded, try again.", type="error")
          return()
        }
        
        if(input$check_data & is.null(input$y_file)){
          #showNotification("To check data consistency you need to upload phenotype response file (Y).", type="error")
          return()
        }
        
        verified <- verify_omic_object(config)
        
        if(verified) {
          
          # Add data to metadata config 
          config$data <- load_data(input$file$datapath)
          
          # Grab values from reactive objects
          curr_data <- data_sources()
          curr_types <- data_types()
          curr_table <- table_data()
          
          # New omic type, register and increment file counter
          if (!config$type %in% curr_types) {
            file_counter(file_counter() + 1)
            data_types(c(curr_types, config$type))
          }
          
          # Check if the label already exists and remove the corresponding row
          existing_label_idx <- which(curr_table$Label == config$label)
          if (length(existing_label_idx) > 0) {
            curr_table <- curr_table[-existing_label_idx, ]
          }
          
          # Add omic with data to available data sources
          curr_data[[config$label]] <- config
          data_sources(curr_data)
          
          # Add omic to preview table
          updated_table <- update_table(curr_table, curr_data[[config$label]])
          table_data(updated_table)
          
          # Update avalable omics in model assembly page
          updateSelectInput(session, "omics_labels", choices = names(data_sources()))
          showNotification("Data was successfully uploaded!", type="message")
          
        } else {
          showNotification("Error: check the data entered in the upload box.", type="error")
          return()
        }
      }
    }
  })
  
  # Duplicate omic in case of combining ability model with one genomic marker information for both groups
  observeEvent(input$upload, {
    if(!input$modtype=="Genotype level"){
      if (input$linkage_type=="Parent Groups ID"){
        for(omic in c(paste(input$label,"1",sep=""),paste(input$label,"2",sep=""))){
          if (omic == paste(input$label,"1",sep="")) {
            data_type = "genomic markers"
            link = "Parent Group 1 ID"
            id_col = input$id_col
          }
          if (omic == paste(input$label,"2",sep="")) {
            data_type = "genomic markers"
            link = "Parent Group 2 ID"
            id_col = input$id_col
          }
          # Create omic metadata object
          
          config <- create_omic_config(input$file, omic, data_type,input$modtype, 
                                       id_col, link=link)
          
          if(is.null(config)) {
            showNotification("No file was uploaded, try again.", type="error")
            return()
          }
          
          if(input$check_data & is.null(input$y_file)){
            #showNotification("To check data consistency you need to upload phenotype response file (Y).", type="error")
            return()
          }
          
          verified <- verify_omic_object(config)
          if(verified) {
            
            # Add data to metadata config 
            config$data <- load_data(input$file$datapath)
            
            # Grab values from reactive objects
            curr_data <- data_sources()
            curr_types <- data_types()
            curr_table <- table_data()
            
            # New omic type, register and increment file counter
            if (!config$type %in% curr_types) {
              file_counter(file_counter() + 1)
              data_types(c(curr_types, config$type))
            }
            
            # Check if the label already exists and remove the corresponding row
            existing_label_idx <- which(curr_table$Label == config$label)
            if (length(existing_label_idx) > 0) {
              curr_table <- curr_table[-existing_label_idx, ]
            }
            
            # Add omic with data to available data sources
            curr_data[[config$label]] <- config
            data_sources(curr_data)
            
            # Add omic to preview table
            updated_table <- update_table(curr_table, curr_data[[config$label]])
            table_data(updated_table)
            
          }else {
            showNotification("Error: check the data entered in the upload box.", type="error")
            return()
          }
        }
        
        # Update avalable omics in model assembly page
        updateSelectInput(session, "omics_labels", choices = names(data_sources()))
        showNotification("Data was successfully uploaded!", type="message")
      }
    }
  }
  )
  
  # set function to check data consistency
  observeEvent(input$upload, {
    
    if(input$check_data){ # Only works when check box is TRUE
      if(is.null(input$y_file)){
        showNotification("To check data consistency you first need to upload phenotype response file (Y).", type="error")
        return()
      }
      
      if(input$data_type=="environmental markers"){ # Function to check consistency between environmental markers and environments
        ids<-data_sources()[["E"]]$data[,1]
        update_data_sources<-data_sources()
        
        update_data_sources[[input$label]]$data <- update_data_sources[[input$label]]$data[update_data_sources[[input$label]]$data[[input$id_col]] %in% ids, ]
        
        data_sources(update_data_sources)
        if (is.na(update_data_sources[[input$label]]$data[1,1]) ){
          showNotification("Please check omic ID column and/or linkage type.", type="error")
        }else{
          data_sources(update_data_sources)
          showNotification("Data was successfully checked!", type="message")
        }
      }else{
      if(input$modtype=="Genotype level"){ # Check at Genotype level
        
        # Function to check genotypes
        update_data_sources<-data_sources()
        
        if(data_sources()[["L"]]$linkage_type=="Compound ID / UID"){
        ids<-data_sources()[["L"]]$uid[,1]
        }else{
        ids<-data_sources()[["L"]]$data[,1]
        }
        
        update_data_sources[[input$label]]$data <- update_data_sources[[input$label]]$data[update_data_sources[[input$label]]$data[[input$id_col]] %in% ids, ]
        
        if (is.na(update_data_sources[[input$label]]$data[1,1]) ){
          showNotification("Please check omic ID column and/or linkage type.", type="error")
        }else{
          data_sources(update_data_sources)
          showNotification("Data was successfully checked!", type="message")
        }
      }else{ # Check at parent level
        
        update_data_sources<-data_sources()
        if(input$linkage_type=="Parent Group 1 ID"){
          ids<-data_sources()[["L1"]]$data[,1]
          update_data_sources[[input$label]]$data <- update_data_sources[[input$label]]$data[update_data_sources[[input$label]]$data[[input$id_col]] %in% ids, ]
        }
        else if(input$linkage_type=="Parent Group 2 ID"){
          ids<-data_sources()[["L2"]]$data[,1]
          update_data_sources[[input$label]]$data <- update_data_sources[[input$label]]$data[update_data_sources[[input$label]]$data[[input$id_col]] %in% ids, ]
        }
        
        else if(input$linkage_type=="Parent Groups ID"){
          ids<-c(data_sources()[["L1"]]$data[,1],data_sources()[["L2"]]$data[,1])
          update_data_sources[[paste(input$label,"1",sep="")]]$data <- update_data_sources[[paste(input$label,"1",sep="")]]$data[update_data_sources[[paste(input$label,"1",sep="")]]$data[[input$id_col]] %in% ids, ]
          update_data_sources[[paste(input$label,"2",sep="")]]$data <- update_data_sources[[paste(input$label,"2",sep="")]]$data[update_data_sources[[paste(input$label,"2",sep="")]]$data[[input$id_col]] %in% ids, ]
          
        }
        else if(input$linkage_type=="Compound ID / UID"){
          ids<-c(data_sources()[["L"]]$uid[,1],data_sources()[["L"]]$uid[,1])
          update_data_sources[[input$label]]$data <- update_data_sources[[input$label]]$data[update_data_sources[[input$label]]$data[[input$id_col]] %in% ids, ]
          
        }
        if (nrow(update_data_sources[[length(update_data_sources)]]$data)<=0 ){
          showNotification("Please check omic ID column and/or linkage type.", type="error")
        }else{
          data_sources(update_data_sources)
          showNotification("Data was successfully checked!", type="message")
        }
      }
      }
    }
  })
  
  # Display data in the available omics table
  output$omics_table <- renderDataTable({
    datatable(table_data(),
              editable = FALSE,
              options = list(
                # Can only select 1 row at a time
                select = list(style="single", target="row"),
                lengthChange = FALSE,
                searching = FALSE
              ))
  })
  
  # Delete a row from the Available Omics table
  observeEvent(input$delete_rows, {
    select_rows <- input$omics_table_rows_selected
    
    if(length(select_rows) > 0) {
      # Create new table object without these rows
      curr_table <- table_data()[-select_rows,]
      # Delete the associated omics
      del_omics <- table_data()[select_rows, 1]
      curr_data <- data_sources()
      curr_data <- curr_data[!(names(curr_data) %in% del_omics)]
      # Update table
      table_data(curr_table)
      data_sources(curr_data)
    }
  })
  
  # Consistently update omics to view their data 
  output$omics_labels <- renderUI({
    selectInput("omics_labels","Select omic to view:",choices=names(data_sources()))
  })
  
  # Display data for a specific omic
  output$table_to_view <- renderDT({
    omic <- input$omics_labels
    DT::datatable(data_sources()[[omic]]$data, options = list(
      searching = TRUE,
      scroller = TRUE,
      scrollX = TRUE
    ))
  })
  
  ### Model tab ----
  
  # Display all available omics as draggable icons
  output$omics_grid <- renderUI({
    # Get list of labels
    curr_data <- data_sources()
    labels <- names(curr_data)
    # Make labels draggable icons
    orderInput("terms","",items=labels,as_source=TRUE, connect="added_terms", style="width:75%;")
  })
  
  # Add the typed in formula 
  observeEvent(input$add_typed_formula, {
    formula <- strsplit(gsub(" ", "", input$typed_formula), "\\+")[[1]]
    
    # Validate formula
    curr_data <- data_sources()
    labels <- names(curr_data)
    
    valid <- validate_formula(formula, labels)
    
    if(valid) {
      # Update model terms with formula
      model_terms(formula)
      # Clear text entry box
      updateTextInput(session, "typed_formula", value="")
    } else {
      showNotification("The formula includes unidentified terms, check the available omics",
                       type = "error")
      updateTextInput(session, "typed_formula", value="")
      return()
    }
    
  })
  
  # Show the model preview
  output$model_formula <- renderUI({
    terms <- model_terms()
    
    if(length(terms) > 0) {
      paste(terms, collapse = " + ")
    } else {
      "No terms added"
    }
  })
  
  # Add the selected terms into the model equation
  observeEvent(input$add_terms, {
    new_terms <- unique(input$added_terms)
    curr_terms <- model_terms()
    latest_terms <- curr_terms
    
    # Ensure new_terms is not NULL or empty
    if(!is.null(new_terms) && length(new_terms) > 0) {
      if(length(new_terms) == 1) {
        # If only one new term and it's not in curr_terms, add it
        if(!(new_terms %in% curr_terms)) {
          latest_terms <- c(curr_terms, new_terms)
        }
      } else if(length(new_terms) > 1) {
        # For multiple new terms, create an interaction term
        int_term <- paste0(sort(new_terms), collapse = "*")
        # Ensure interaction term is a single string and not already in curr_terms
        if(length(int_term) == 1 && !int_term %in% curr_terms) {
          latest_terms <- c(curr_terms, int_term)
        }
      }
    }
    
    # Update reactive variable
    model_terms(latest_terms)
    
    # Reset source orderInput
    curr_data <- data_sources()
    labels <- names(curr_data)
    updateOrderInput(session, "terms", items=labels)
    # Clear items in added_terms box
    updateOrderInput(session, "added_terms", items=list())
  })
  
  # Delete terms to NOT add them into the model equation
  observeEvent(input$clear_terms, {
    # Reset the source orderInput
    ds <- data_sources()
    labels <- names(ds)
    updateOrderInput(session, inputId = "terms", items = labels)
    
    # Clear the items in the "added_terms" orderInput
    updateOrderInput(session, inputId = "added_terms", items = list())
  })
  
  # Clear terms from model formula
  observeEvent(input$clear_model, {
    model_terms(character(0))
  })
  
  # Go back one step when creating model
  observeEvent(input$prev_model, {
    curr_terms <- model_terms()
    if (length(curr_terms) > 0) {
      model_terms(curr_terms[-length(curr_terms)])
    }
  })
  
  # Save the model formula
  observeEvent(input$save_model, {
    curr_terms <- model_terms()
    default_name <- paste0(sapply(curr_terms, function(x) gsub("[*]", "", x)), collapse = "+")
    
    if (length(curr_terms) > 0) {
      showModal(modalDialog(
        title = "Save Model",
        textInput("saved_model", "Name for model", default_name),
        actionButton("save_name_btn", "Save"),
        easyClose = TRUE,
        footer = NULL
      ))
    }
  })
  
  # Give saved model a name to recall it; opens up as a pop up
  observeEvent(input$save_name_btn, {
    if (!is.null(input$saved_model) && input$saved_model != "") {
      curr_saved <- saved_models()
      curr_saved[[input$saved_model]] <- model_terms()
      saved_models(curr_saved)
      removeModal()
      
      # Update models available for validation
      updateSelectInput(session, "saved_models", choices = names(saved_models()))
      
      showNotification(paste0(input$saved_model, " was successfully saved!"),
                       type = "message")
    }
  })
  
  ### Train / validate tab ----
  
  # Update available models
  output$saved_models <- renderUI({
    selectInput("saved_models",
                "Select a saved model",
                choices = names(saved_models()))
  })
  
  # Update and sync the slider and input fields for hyperparameters
  observeEvent(input$n_iter_input, {
    updateSliderInput(session, "n_iter_slider", value = input$n_iter_input)
  })
  observeEvent(input$n_iter_slider, {
    updateNumericInput(session, "n_iter_input", value = input$n_iter_slider)
  })
  observeEvent(input$burn_in_input, {
    updateSliderInput(session, "burn_in_slider", value = input$burn_in_input)
  })
  observeEvent(input$burn_in_slider, {
    updateNumericInput(session, "burn_in_input", value = input$burn_in_slider)
  })
  
  # Run the selected cross-validation(s)
  observeEvent(input$run_cv, {
    
    # Restart logs
    logging(session)
    
    ### 1/6 -- Gather all parameters ###
    logging(session, "[Step 1/6] Checking parameters.")
    
    # Get the model selected
    model_selected <- input$saved_models
    model_equation <- lapply(saved_models()[[model_selected]], reformat_model)
    
    logging(session, paste0("Model selected: ", model_selected))
    logging(session, paste0(
      "Model equation: ", 
      paste0(saved_models()[[model_selected]], collapse="+"))
    )
    
    # Tuning parameters
    nIter <- input$n_iter_input
    burnIn <- input$burn_in_input
    std <- input$std
    ctr <- input$ctr
    wht <- FALSE #input$wht
    esc <- FALSE #input$esc
    nan_freq <- input$nan_freq
    folds <- input$folds
    prop_maf_j <- input$prop_maf_j
    seed <- input$set_seed
    
    # CV methods
    cv1 <- input$cv1
    cv2 <- input$cv2
    cv0 <- input$cv0
    cv00 <- input$cv00
    
    # Data for training/testing
    ds <- data_sources()
    # Separate Y data from the rest
    ydata <- phen_data()
    # Get unique list of labels
    labels <- unique(unlist(model_equation))
    
    
    ### 2/6 -- Creating matrices ###
    logging(session, "Step [2/6] Creating G matrices.")
    
    # Loop through unique matrices to create G matrices
    for (label in labels) {
      # Get the ID column shared by omic and Y data
      curr_term <- ds[[label]]
      curr_term$y_col <- get_join_id(curr_term, ydata[!names(ydata) %in% "data"])
      if (is.null(curr_term$y_col)) {
        showNotification("There is a linkage error between your data, please check your ID columns",
                         type = "error")
        return()
      }
      # Create incidence matrix using share ID by Y and term
      logging(session, paste0("1 ==> Preparing Z matrix for: ", label))
      
      # Perform necessary data processing
      logging(session, paste0("2 ==> Pre-procesing data for: ", label))
      
      # Create G matrices 
      logging(session, paste0("3 ==> Creating G matrix for: ", label))
      
      result <- create_g_matrix(curr_term, ydata, wht, ctr, std, nan_freq, prop_maf_j)
      G <- result[["G"]]
      EVD <- result[["EVD"]]
      print(dim(G))
      print(dim(EVD$vectors))
      # Save the results as RDA files
      logging(session, paste0("4 ==> Saving results for: ", label))
      
      outdir <- file.path(tmpdir, label)
      if (!dir.exists(outdir)) { dir.create(outdir, recursive = TRUE)}
      
      save(G, file = file.path(outdir, "G.rda"))
      save(EVD, file = file.path(outdir, "EVD.rda"))
      
      # Plot EVD
      #logging(session, paste0("5 ==> Generating EVD plots for: ", label))
      #get_evd_plots(EVD, file.path("tmp",label))
      
    }
    
    ### 3/6 -- Creating interaction matrices, if applicable ###
    logging(session, "Step [3/6] Creating interaction matrices, if applicable.")
    
    for (i in seq_along(model_equation)) {
      # Generate interaction matrix for explicitly defined interactions
      if (length(model_equation[[i]]) > 1 & is.vector(model_equation[[i]])) {
        term <- paste(model_equation[[i]], collapse = "_")
        result <- create_interaction_matrix(model_equation[[i]], tmpdir)
        
        # Save the results as RDA files
        logging(session, paste0("1 ==> Saving results for: ", term))
        
        outdir <- file.path(tmpdir, term)
        if (!dir.exists(outdir)) { dir.create(outdir, recursive = TRUE)}
        
        G <- result[["G"]]
        EVD <- result[["EVD"]]
        
        save(G, file=file.path(outdir,"G.rda"))
        save(EVD, file=file.path(outdir,"EVD.rda"))
        
        #logging(session, paste0("2 ==> Generating EVD plots for: ", term))
        #get_evd_plots(EVD, outdir)
        
      }
    }
    
    ### 4/6 -- Preparing data for cross-validation ###
    logging(session, "Step [4/6] Assigning folds for cross validation")
    cv_ids <- prep_data_for_cv(ydata, folds, cv1, cv2, cv0, cv00)
    #write.csv(cv_ids,file="cv_ids.csv",row.names=F)
    # Get terms to add as column names
    all_terms <- sapply(model_equation, function(x) paste(x, collapse="_"))
    
    ### 5/6 -- Performing cross-validations ###
    logging(session, "Step [5/6] Performing cross-validations")
    
    #####################
    
    set.seed(seed)
    
    outdir <- file.path(tmpdir, model_selected)
    if (!dir.exists(outdir)) { dir.create(outdir, recursive = TRUE) }
    
    if (cv00) { 
      logging(session, "CV00")
      cv1_preds <- get_predictions(cv_ids, tmpdir, model_selected, all_terms, "cv00", 
                                   folds, nIter=nIter, burnIn=burnIn)
    }
    if (cv1) { 
      logging(session, "CV1")
      cv1_preds <- get_predictions(cv_ids, tmpdir, model_selected, all_terms, "cv1", 
                                   folds, nIter=nIter, burnIn=burnIn)
    }
    if (cv2) { 
      logging(session, "CV2")
      cv1_preds <- get_predictions(cv_ids, tmpdir, model_selected, all_terms, "cv2", 
                                   folds, nIter=nIter, burnIn=burnIn)
    }
    if (cv0) { 
      logging(session, "CV0")
      cv1_preds <- get_predictions(cv_ids, tmpdir, model_selected, all_terms, "cv0", 
                                   folds, nIter=nIter, burnIn=burnIn)
    }
    
    logging(session, "Cross-validation(s) were completed successfully!")
    
    ### Get full data variance components ###
    logging(session, "Step [6/6] Calculating variance components.")
    get_predictions(cv_ids, tmpdir, model_selected, all_terms, folds=folds, nIter=nIter,
                    burnIn=burnIn, cv=NULL)
    
    ### Create visualizations ###
    logging(session, "Calculating evaluation metrics: ")
    get_model_accuracy(3, file.path(tmpdir,"output"), FALSE)
    get_model_accuracy(3, file.path(tmpdir,"output"), TRUE)
    logging(session, "==> Accuracy")
    get_model_rmse(3, file.path(tmpdir,"output"), FALSE)
    get_model_rmse(3, file.path(tmpdir,"output"), TRUE)
    logging(session, "==> RMSE")
    get_model_varcomps(file.path(tmpdir,"output"))
    logging(session, "==> Variance components")
    
    ### Complete, notify user ###
    logging(session, "END: Training and validation complete!")
    
    Sys.sleep(1)
    showModal(
      modalDialog(title = "Notification",
                  "Cross-validation(s) completed! Please go to the 'Results' tab to view outcome",
                  easyClose = TRUE,
                  footer = tagList(actionButton("go_to_results", "View Results")))
    )
  })
  
  observeEvent(input$go_to_results, {
    removeModal()
    updateTabItems(session, "tabs", "results")
  })
  
  # Render log messages reactively
  output$cv_logs <- renderText({
    invalidateLater(1000, session)  # Check for updates every second
    logs()
  })
  
  ### Results tab ----
  
  # Reactive poll to check for new files every 5 seconds
  files_reactive <- reactivePoll(
    5000, 
    session,
    checkFunc = function() {
      list.files(file.path(tmpdir, "output"), pattern = "\\.csv$", recursive = TRUE)
    },
    valueFunc = function() {
      list.files(file.path(tmpdir, "output"), pattern = "\\.csv$", recursive = TRUE)
    })
  
  # Generate dynamic UI for file menu
  output$results_menu <- renderUI({
    files <- files_reactive()
    if (length(files) > 0) {
      selectInput("selected_file", "Choose a file:", choices = files)
    } else {
      p("No CSV files found.")
    }
  })
  
  # Render the data table based on the selected file
  output$results_table <- renderDT({
    req(input$selected_file)
    file_path <- file.path(tmpdir, "output", input$selected_file)
    
    if (file.exists(file_path)) {
      DT::datatable(read.csv(file_path))
    } else {
      DT::datatable(data.frame())
    }
  })
  
  # Download CV data
  output$download_results <- downloadHandler(
    filename = function() {
      paste0("results_", Sys.Date(), ".zip")
    },
    content = function(file) {
      zip::zipr(file, file.path(tmpdir,"output"))
    }
  )
  
  # Data visualizations
  output$plot <- renderImage({
    req(input$metrics)
    
    image_list <- list(
      "Accuracy (by environment)" = "accuracy_by_env.png",
      "Accuracy (overall)" = "accuracy_overall.png",
      "RMSE (by environment)" = "rmse_by_env.png",
      "RMSE (overall)" = "rmse_overall.png",
      "Variance components" = "varcomps.png"
    )
    
    filename <- image_list[[input$metrics]]
    img_path <- file.path(tmpdir, "output", "visualizations", filename)
    
    # Check if the selected image file exists
    if (file.exists(img_path)) {
      list(src = img_path, contentType = "image/png", alt = "Plot", style = "min-width: 70%; max-width: 100%; min-height: 500px; max-height: 600px;")
    } else {
      return(NULL)
    }
    
  }, deleteFile = FALSE)
  
  output$plot_or_message <- renderUI({
    image_list <- list(
      "Accuracy (by environment)" = "accuracy_by_env.png",
      "Accuracy (overall)" = "accuracy_overall.png",
      "RMSE (by environment)" = "rmse_by_env.png",
      "RMSE (overall)" = "rmse_overall.png",
      "Variance components" = "varcomps.png"
    )
    
    filename <- image_list[[input$metrics]]
    img_path <- file.path(tmpdir, "output", "visualizations", filename)
    
    if (file.exists(img_path)) {
      imageOutput("plot")
    } else {
      renderText("No plots available yet")
    }
  })
  
  # AMMI and Bayesian indexes (GEI analyses)
  
  ### GEI Data tab ----
  # Generate a dynamic upload panel to select columns
  output$gei_y_panel<-renderUI({
      panel<-div(create_input_box('*Genotype ID column:', "gid_col", 1),
                 create_input_box("*Environment ID column:", "eid_col", 2),
                 create_input_box("*Individual ID (UID) column:", "uid_col",3),
                 create_input_box("*Target trait column:", "trait_col", 4))
  })
  
  # Loading Y file after upload button is pressed
  observeEvent(input$gei_y_load, {
    config <- create_omic_config(input$gei_y_file, "Y", "phenotypic",modtype="Genotype level", input$trait_col, 
                                 gid_col = input$gid_col, eid_col = input$eid_col,
                                 uid_col = input$uid_col)
    config$AMMI<-TRUE
    # Check if data was uploaded.
    if(is.null(config)) {
      showNotification("No file was uploaded, try again.", type="error")
      return()
    }
    
    # Verify uploaded data
    verified <- verify_omic_object(config, phen = TRUE,modtype="Genotype level")

    if(verified) {
      
      # Add data to metadata config 
      config$data <- load_data(input$gei_y_file$datapath)

      # Grab values from reactive objects
      phen_data(config)
      curr_data <- data_sources()
      
      if (is.na(input$uid_col)){
        for (omic in c("E","L")) {
          
          if(omic == "E") {
            data_type = "Environment IDs"
            link = "Environment ID"
            id_col = input$eid_col
          } else {
            data_type = "Line IDs"
            link = "Genotype / Line ID"
            id_col = input$gid_col
          }
 
          # Get config for ID column
          curr_omic <- create_omic_config(input$gei_y_file, omic, data_type,modtype="Genotype level",id_col,
                                          link = link)
          # Load data as data.frame for that respective column
          curr_omic$data <- config$data[, id_col, drop=FALSE]

          curr_data[[curr_omic$label]] <- curr_omic
          
          # Add omic to preview table
          curr_table <- table_data()
          updated_table <- update_table(curr_table, curr_data[[curr_omic$label]])
          table_data(updated_table)
        }}else{
          for (omic in c("E","L","R")) {
            
            if(omic == "E") {
              data_type = "Environment IDs"
              link = "Environment ID"
              id_col = input$eid_col
            } else if (omic == "L") {
              data_type = "Line IDs"
              link = "Genotype / Line ID"
              id_col = input$gid_col
            } else{
              data_type = "Line Replicates"
              link = "Genotype / Line ID"
              id_col = input$uid_col
            }
            
            # Get config for ID column
            curr_omic <- create_omic_config(input$gei_y_file, omic, data_type,modtype="Genotype level",id_col,
                                            link = link)
            # Load data as data.frame for that respective column
            curr_omic$data <- config$data[, id_col, drop=FALSE]
            
            curr_data[[curr_omic$label]] <- curr_omic
            
            # Add omic to preview table
            curr_table <- table_data()
            updated_table <- update_table(curr_table, curr_data[[curr_omic$label]])
            table_data(updated_table)
          }
        }
      
      # Save data added for E and L
      data_sources(curr_data)
      # Update avalable omics in model assembly page
      updateSelectInput(session, "gei_labels", choices = names(data_sources()))
      showNotification("Data was successfully uploaded!", type="message")
      
    } else {
      showNotification("Error: check the data entered in the upload box.", type="error")
      return()
    }
    
    }
  )
  
  # Display data in the available table
  output$gei_table <- renderDataTable({
    datatable(table_data()[,1:4],
              editable = FALSE,
              options = list(
                # Can only select 1 row at a time
                select = list(style="single", target="row"),
                lengthChange = FALSE,
                searching = FALSE
              ))
  })
  
  # Delete a row from the Available table
  observeEvent(input$gei_delete_rows, {
    select_rows <- input$gei_table_rows_selected
    
    if(length(select_rows) > 0) {
      # Create new table object without these rows
      curr_table <- table_data()[-select_rows,]
      # Delete the associated omics
      del_gei <- table_data()[select_rows, 1]
      curr_data <- data_sources()
      curr_data <- curr_data[!(names(curr_data) %in% del_gei)]
      # Update table
      table_data(curr_table)
      data_sources(curr_data)
    }
  })
  
  # Consistently update omics to view their data 
  output$gei_labels <- renderUI({
    selectInput("gei_labels","Select omic to view:",choices=names(data_sources()))
  })
  
  # Display data for a specific omic
  output$gei_table_to_view <- renderDT({
    gei <- input$gei_labels
    DT::datatable(data_sources()[[gei]]$data, options = list(
      searching = TRUE,
      scroller = TRUE,
      scrollX = TRUE
    ))
  })
  
  ### GEI Run analysis tab ----
  
  # Update and sync the slider and input fields for hyperparameters
  observeEvent(input$gei_n_iter_input, {
    updateSliderInput(session, "gei_n_iter_slider", value = input$gei_n_iter_input)
  })
  observeEvent(input$gei_n_iter_slider, {
    updateNumericInput(session, "gei_n_iter_input", value = input$gei_n_iter_slider)
  })
  observeEvent(input$gei_burn_in_input, {
    updateSliderInput(session, "gei_burn_in_slider", value = input$gei_burn_in_input)
  })
  observeEvent(input$gei_burn_in_slider, {
    updateNumericInput(session, "gei_burn_in_input", value = input$gei_burn_in_slider)
  })
  
  # Generate a dynamic weigh panel if waasy is selected
  output$weigh_panel <- renderUI({
    if(input$bwaasy){
      
    # weights
      panel <- div(
        style="width: 200%;",
        tags$head(
          tags$style(HTML("
      .irs-min, 
      .irs-max { 
        display: none !important; 
      }
    ")) 
        ),
        tags$label('*Weights between response variable and absolute scores:', 
                   style = "display: block; margin-bottom: 5px;"),
        div(
          style = "margin-top: 20px; width: 100%; position: relative;", # Define largura e posicionamento relativo
          sliderInput(
            inputId = "weights",
            label = NULL, # Remove label padrão
            min = 0,
            max = 100,
            value = 50,
            step = 1
          ),
          # Labels personalizadas no lugar dos valores padrão
          div(
            style = "display: flex; justify-content: space-between; position: absolute; width: 100%; top: -12px;",
            tags$span("Response variable", style = "font-size: 12px; color: black;"), # Label personalizada no lugar do 0
            tags$span("Absolute scores", style = "font-size: 12px; color: black;")    # Label personalizada no lugar do 100
          )
        )
      )

    }
    }
  )
  
  
  # Run the selected analyses
  observeEvent(input$run_gei, {
    
    # Restart logs
    logging(session,target="gei_logs")
    
    ### 1/6 -- Gather all parameters ###
    logging(session, "[Step 1/6] Checking tuning parameters.",target="gei_logs")
    
    # Tuning parameters
    nIter <- input$gei_n_iter_input
    burnIn <- input$gei_burn_in_input
    thin <- input$gei_thin_input
    seed <- input$gei_set_seed
    trait_col<-input$trait_col
    gid_col<-input$gid_col
    eid_col<-input$eid_col
    uid_col<-input$uid_col
    chains<-input$gei_chains
    
    # GEI Methods
    ammi <- input$bammi
    waas <- input$bwaas
    waasy <- input$bwaasy
    stabdist<- input$bstab
    
    if(waasy){ammi=TRUE}
    if(waas){ammi=TRUE}
    if(stabdist){ammi=TRUE}
    
    # Data for run analysis
    ds <- data_sources()
    
    # Building DF to pre-analisys
    ydata <- phen_data()
    
    # Get list of labels
    labels <- names(df)
    
    ### 2/6 -- Obtaining initial parameters ###
    logging(session, "Step [2/6] Obtaining initial parameters.",target="gei_logs")
    
    if (!ammi) {
      showNotification("Please check at least one procedure to be executed.",
                       type = "error")
      return()
    }
    
    # Fit initial linear model to obtain a residual value
    logging(session, paste("1 ==> Fitting an initial bilinear model to obtain a residual estimative."),target="gei_logs")
      
    ### 3/6 -- Generating Gibbs sampler chains ###
    logging(session, "Step [3/6] Generating Gibbs sampler chains.",target="gei_logs")
    
    set.seed(seed)
    outdir <- file.path(tmpdir,"output/AMMI")
    if (!dir.exists(outdir)) { dir.create(outdir, recursive = TRUE)}
    
    chain_list<-list()
    for(chain in 1:chains){
    bilinear<-bilinear_model(ydata$data, trait_col, eid_col, gid_col, uid_col)
      
    exit1<-NULL
    nombre1 <- file.path(outdir,paste0("chain_",chain,"_output.csv",sep=""))
    write.table(exit1,file=nombre1,append = FALSE, quote = TRUE, sep = " , ",row.names = FALSE,col.names = FALSE)
    
    logging(session, paste0(chain," ==> Generating chain number: ", chain),target="gei_logs")
    
      # Values for hyper-parameters 
      sigma2_u<<-1*10^15
      alfa_sig_ap<<-0
      beta_sig_ap<<-0
      mu_tao_ap<<-rep(0,r)
      sig_tao_ap<<-(1*10^15)
      mu_del_ap<<-rep(0,cc)
      sig_del_ap<<-(1*10^15)
      mu_lamb_ap<<-0
      sig_lamb_ap<<-1*10^15
      
      # Initial values for bilinear terms 
      alphas0 <- bilinear[[12]]
      gammas0 <- bilinear[[13]]
      lambdas0 <- bilinear[[11]]
      
    times<-1
    gibbs_chain(exit1,nombre1,nIter,times,alphas0,gammas0,lambdas0)
    
    chain_list[[chain]]<-read.csv(nombre1,header=FALSE)
    }
    
    logging(session, "Chains were sucessfully generated!",target="gei_logs")

    
    ### 4/6 -- Conducting analyses ###
    logging(session, "Step [4/6] Conducting analyses.",target="gei_logs")
    logging(session, paste0("1 ==> Removing burn-in iterations. "),target="gei_logs")
    chain_list<-burnin_filter(chain_list,burnIn,r,cc,k)

    # TODO dependence and stationary tests
    
    logging(session, paste0("2 ==> Applying thin. "),target="gei_logs")
    chain_list<-thin_filter(chain_list,nIter,burnIn,thin)

    logging(session, paste0("3 ==> Merging chains. "),target="gei_logs")
    chains<-do.call(rbind,chain_list)
    # Save chains
    
    write.table(chains, file=file.path(outdir, 'chains.csv'),quote=TRUE,sep=",",row.names=F,col.names=TRUE)

    mu    <<-  chains[,1]
    tau   <<-  chains[,2]
    tao <<-  chains[,3:(2+r)]
    beta  <-  chains[,(2+r+1):(2+r+cc)]
    
    alpha <- chains[,(2+r+cc+1):(2+r+cc+r*k)] #U
    gamma <- chains[,(2+r+cc+r*k+1):(2+r+cc+r*k+cc*k)] #VT
    lambda <- chains[,(2+r+cc+r*k+cc*k+1):(2+r+cc+r*k+cc*k+k)] #lambda
    
    saving_parameters(outdir,mu,tau,tao,beta,lambda,alpha,gamma)
    
    if(waas){
    logging(session, paste0("4 ==> WAAS index. "),target="gei_logs")
      
    out <- file.path(tmpdir,"output/WAAS")
    if (!dir.exists(out)) { dir.create(out, recursive = TRUE)}
        
    waas<-waas_index(outdir,bilinear[[4]])
    #write.csv(waas,file="waas.csv",row.names=F)
    write.csv(waas,file=file.path(out, 'waas.csv'),row.names=F)
    waas_plot(waas,tmpdir,bilinear[[4]])}
    if(waasy){
    logging(session, paste0("4 ==> WAASY index. "),target="gei_logs")
    
    out <- file.path(tmpdir,"output/WAASY")
    if (!dir.exists(out)) { dir.create(out, recursive = TRUE)}
    
    wts<-input$weights
    waasy<-waasy_index(outdir,bilinear[[4]],wts)
    #write.csv(waasy,file="waasy.csv",row.names=F)
    write.csv(waasy,file=file.path(out, 'waasy.csv'),row.names=F)
    
    waasy_plot(waasy,tmpdir,bilinear[[4]])
    }
    if(stabdist){stabdist<-stabdist_index(outdir,bilinear[[4]])
    logging(session, paste0("4 ==> STABDIST index. "),target="gei_logs")
    
    out <- file.path(tmpdir,"output/STABDIST")
    if (!dir.exists(out)) { dir.create(out, recursive = TRUE)}
    
    means<-mean_gen(outdir,bilinear[[4]])
    stabdist_merged<-data.frame(means,stabdist[,-1])
    stabdist_merged$Genotype<-tratamientos
    write.csv(stabdist_merged,file=file.path(out, 'stabdist.csv'),row.names=F)
    
    stabdist_plot(stabdist_merged,tmpdir,bilinear[[4]])}
    
    #hpd<-hpd_intervals(chains,outdir)
    #write.table(hpd,file="hpd.csv",row.names=F,col.names=T,sep=",")
    logging(session, paste0("5 ==> Calculating high posterior density intervals. "),target="gei_logs")
    hpd_intervals(chains,outdir)
    logging(session, paste0("6 ==> Obtaining estimated parameters. "),target="gei_logs")
 
    scores<-average_scores(outdir,k=2,r,cc)
    
    logging(session, paste0("7 ==> Analyses were sucessfully done!. "),target="gei_logs")
    

    ### 5/6 -- Building plots. ###
    logging(session, "Step [5/5] Building plots.",target="gei_logs")
    
    logging(session, paste0("1 ==> Trace plots. "),target="gei_logs")
    
    trace_plots(tmpdir,tau,mu,lambda)
    
    if(input$bammi){
      logging(session, paste0("2 ==> Biplot. "),target="gei_logs")
      bayesian_biplot(tmpdir,scores,lambda,r)}
    
    ### Complete, notify user ###
    logging(session, "END: The analyses are completed!")
    
        Sys.sleep(1)
        showModal(
          modalDialog(title = "Notification",
                      "Analyse(s) completed! Please go to the 'Results' tab to view outcome",
                      easyClose = TRUE,
                      footer = tagList(actionButton("go_to_results_gei", "View Results")))
        )
  })

  observeEvent(input$go_to_results_gei, {
    removeModal()
    updateTabItems(session, "tabs", "gei_results")
  })
  
  # Render log messages reactively
  output$gei_logs <- renderText({
    invalidateLater(1000, session)  # Check for updates every second
    logs()
  })
  
  ### Results tab ----

  # Reactive poll to check for new files every 5 seconds
  #gei_files_reactive <- reactivePoll(
    #5000, 
    #session,
    #checkFunc = function() {
    #  list.files(file.path(tmpdir, "AMMI"), pattern = "\\.csv$", recursive = TRUE)
    #},
    #valueFunc = function() {
    #  list.files(file.path(tmpdir, "AMMI"), pattern = "\\.csv$", recursive = TRUE)
    #})

  # Generate dynamic UI for file menu
  output$gei_results_menu <- renderUI({
    files_gei <- files_reactive()
    if (length(files_gei) > 0) {
      selectInput("gei_selected_file", "Choose a file:", choices = files_gei)
    } else {
      p("No CSV files found.")
    }
  })
  
  # Render the data table based on the selected file
  output$gei_results_table <- renderDT({
    req(input$gei_selected_file)
    file_path <- file.path(tmpdir, "output",input$gei_selected_file)
    
    if (file.exists(file_path)) {
      DT::datatable(read.csv(file_path))
    } else {
      DT::datatable(data.frame())
    }
  })
  
  # Download analyses results
  output$download_results_gei <- downloadHandler(
    filename = function() {
      paste0("results_", Sys.Date(), ".zip")
    },
    content = function(file) {
      zip::zipr(file, file.path(tmpdir,"output"))
    }
  )
  
  # Data visualizations
  
  observeEvent(input$run_gei, {
    selected_categories <- c()
    
    if (input$bammi) selected_categories <- c(selected_categories, "bammi")
    if (input$bwaas) selected_categories <- c(selected_categories, "bwaas")
    if (input$bwaasy) selected_categories <- c(selected_categories, "bwaasy")
    if (input$bstab) selected_categories <- c(selected_categories, "bstab")
    
    combined_options <- c(unlist(metrics_choices_gei[selected_categories],use.names=F))

    updateSelectInput(
      inputId = "gei_metrics",
      choices = combined_options
    )

  })

  
  output$plot_gei <- renderImage({
    req(input$gei_metrics)
    
    image_list <- list(
      "Trace mu" = "trace_mu.png",
      "Trace tau" = "trace_tau.png",
      "Trace Lambda PC1" = "trace_lambda1.png",
      "Trace Lambda PC2" = "trace_lambda2.png",
      "AMMI Plot" = "bayesian_ammi.png",
      "Waas Plot" = "waas_plot.png",
      "Waasy Plot"= "waasy_plot.png",
      "Stabdist Plot"= "stabdist_plot.png"
    )
    
    filename <- image_list[[input$gei_metrics]]

    img_path <- file.path(tmpdir, "output", "visualizations", filename)
    # Check if the selected image file exists
    if (file.exists(img_path)) {
      list(src = img_path, contentType = "image/png", alt = "Plot", style = "min-width: 70%; max-width: 100%; min-height: 500px; max-height: 600px;")
    } else {
      return(NULL)
    }
    
  }, deleteFile = FALSE)
  
  output$plot_or_message_gei <- renderUI({
    image_list <- list(
      "Trace mu" = "trace_mu.png",
      "Trace tau" = "trace_tau.png",
      "Trace Lambda PC1" = "trace_lambda1.png",
      "Trace Lambda PC2" = "trace_lambda2.png",
      "AMMI Plot" = "bayesian_ammi.png",
      "Waas Plot" = "waas_plot.png",
      "Waasy Plot"= "waasy_plot.png",
      "Stabdist Plot"= "stabdist_plot.png"
    )
    
    filename <- image_list[[input$gei_metrics]]
    img_path <- file.path(tmpdir, "output", "visualizations", filename)

    if (file.exists(img_path)) {
      imageOutput("plot_gei")
    } else {
      renderText("No plots available yet")
    }
  })
  
  # Clean up temporary directory when the session ends
  session$onSessionEnded(function() {
    unlink(tmpdir, recursive = TRUE)
  })
  
}