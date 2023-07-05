library(shiny)
library(xml2)
library(shinyWidgets)
library(magrittr)
library(shinyBS)
library(XML)

xml_file <- "/Users/fs8/Desktop/Project/MAAT/src/system_models/leaf/leaf_default.xml"
xml_data <- read_xml(xml_file)

leaf_static <- read_xml("/Users/fs8/Desktop/Project/MAAT/src/system_models/leaf/init_files/leaf_user_static.xml")

leaf_pars <- xml_data %>% xml_find_first(".//pars/leaf")
leaf_pars_name <- xml_children(leaf_pars) %>% xml_name()
leaf_pars_value <- xml_children(leaf_pars) %>% xml_text()

leaf_pars_ggc <- xml_find_all(leaf_pars, ".//*[count(*) > 1]") %>% xml_name()

ui <- fluidPage(
  titlePanel(""),
  fluidRow(
    column(2, 
           wellPanel(
             prettyRadioButtons(
               "parts",
               "Please Select",
               choices = "pars",
               status = "danger",
               fill = TRUE
             )
           )
    ), 
    column(10, 
           bsCollapse(id = "collapse_ex", open = "Parameters", multiple = TRUE,
                      bsCollapsePanel("Parameters w/Sub-Script", 
                                      selectInput("gc_pars", "Select Parameter", choices = leaf_pars_ggc),
                                      uiOutput("grandchild_textboxes"),
                                      style = "success", 
                                      fluidRow(
                                        column(12, align = "right",
                                               actionButton("pss_update_button", "Update")
                                        )
                                      )
                      ),
                      bsCollapsePanel("Parameters", 
                                      uiOutput("ui"), 
                                      style = "info", 
                                      fluidRow(
                                        column(12, align = "right",
                                               actionButton("p_update_button", "Update")
                                        )
                                      )
                      )
           )
    )
  )
)

server <- function(input, output, session) {
  pars_gc_names <- reactiveValues(names = NULL)
  
  observeEvent(input$parts, {
    leaf_pars <- leaf_pars_name
    
    output$ui <- renderUI({
      fluidRow(
        lapply(leaf_pars, function(node_name) {
          if (node_name %in% leaf_pars_ggc) {
            return(NULL)
          }
          if (node_name %in% leaf_pars_name) {
            node_index <- match(node_name, leaf_pars_name)
            default_value <- leaf_pars_value[node_index]
            
            column(6, textInput(node_name, label = node_name, value = default_value))
          } 
        }) 
      )
    })
    
    output$grandchild_textboxes <- renderUI({
      pars_gc_options <- input$gc_pars # grandchild_options gets the selected parameter from the dropdown menu
      
      pars_grandchilds <- xml_data %>% xml_find_first(paste0(".//pars/leaf/", pars_gc_options))
      pars_gc_names$names <- xml_children(pars_grandchilds) %>% xml_name() # Update the reactive value
      pars_gc_values <- xml_children(pars_grandchilds) %>% xml_text()
      
      lapply(pars_gc_names$names, function(pars_gc_nodes) {
        pars_node_index <- match(pars_gc_nodes, pars_gc_names$names)
        pars_def_value <- pars_gc_values[pars_node_index]
        
        column(6, textInput(pars_gc_nodes, label = pars_gc_nodes, value = pars_def_value))
      })
    })
  })
  
  observeEvent(input$pss_update_button, {
    updated_xml_data <- xml_data
    pars_gc_options <- input$gc_pars
    
    # Retrieve the input data from the "Parameters w/Sub-Script" panel
    input_data <- lapply(pars_gc_names$names, function(pars_gc_nodes) {
      if (pars_gc_nodes %in% pars_gc_names$names) {
        input_value <- input[[pars_gc_nodes]]
        return(data.frame(node_name = pars_gc_nodes, value = input_value))
      }
    })
    
    # Remove NULL elements from the list
    input_data <- input_data[!sapply(input_data, is.null)]
    
    # Combine the list of data frames into a single data frame
    input_data_df <- do.call(rbind, input_data)
    
    for (i in seq_len(nrow(input_data_df))) {
      # Access individual rows of the data frame
      input_item <- input_data_df[i, ]
      
      # Access values from the input_item row
      node_name <- input_item$node_name
      value <- input_item$value
      
      # Find the corresponding XML node and update its text
      xpath <- paste0(".//pars/leaf/", pars_gc_options, "/", node_name)
      # print(xpath)
      node <- xml_find_first(updated_xml_data, xpath)
      
      if (!is.null(node)) {
        xml_text(node) <- value
      }
    }
    
    new_xml_file <- "/Users/fs8/Desktop/Project/test/updated_leaf.xml"
    write_xml(updated_xml_data, new_xml_file)
    
    showModal(modalDialog(
      title = "Update Successful",
      "The XML file has been updated.",
      easyClose = TRUE
    ))
  })
  
  observeEvent(input$p_update_button, {
    # Create a copy of the XML data
    updated_xml_data <- xml_data
    
    # Update the values in the XML data based on the user input
    for (node_name in leaf_pars_name) {
      if (node_name %in% leaf_pars_ggc) {
        next
      }
      
      new_value <- input[[node_name]]
      
      # Find the corresponding XML node and update its text
      node <- xml_find_first(updated_xml_data, paste0(".//pars/leaf/", node_name))
      
      if (!is.null(node)) {
        xml_text(node) <- new_value
      }
    }
    
    # Save the updated XML data to a new file
    new_xml_file <- "/Users/fs8/Desktop/Project/test/updated_leaf.xml"
    write_xml(updated_xml_data, new_xml_file)

    # Show a success message or perform any additional actions
    showModal(modalDialog(
      title = "Update Successful",
      "The XML file has been updated.",
      easyClose = TRUE
    ))
  })
}

shinyApp(ui = ui, server = server)
