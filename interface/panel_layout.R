library(shiny)
library(xml2)
library(shinyWidgets)
library(magrittr)
library(shinyBS)
library(XML)
library(glue)

# Read the XML file and extract the options
xml_data <- read_xml("/Users/fs8/Desktop/Project/MAAT/src/system_models/leaf/leaf_default.xml")

# Find parent nodes
parent_nodes <- xml_data %>% xml_children() %>% xml_name() 

# Find grandchild node names for each parent node
fnames <- xml_data %>% xml_find_first(".//fnames/leaf")
fnames_names <- xml_children(fnames) %>% xml_name()
fnames_values <- xml_children(fnames)  %>% xml_text()
fnames_gc <- xml_find_all(fnames, ".//*[count(*) > 1]") %>% xml_name()

pars <- xml_data %>% xml_find_first(".//pars/leaf")
pars_names <- xml_children(pars) %>% xml_name()
pars_values <- xml_children(pars) %>% xml_text()
pars_gc <- xml_find_all(pars, ".//*[count(*) > 1]") %>% xml_name()

env <- xml_data %>% xml_find_first(".//env/leaf")
env_names <- xml_children(env) %>% xml_name()
env_values <- xml_children(env) %>% xml_text()
env_gc <- xml_find_all(env, ".//*[count(*) > 1]") %>% xml_name()

# Define the UI
ui <- fluidPage(
  titlePanel(""),
  fluidRow(
    column(2, 
           wellPanel(
             prettyRadioButtons(
               "parts",
               "Please Select",
               choices = parent_nodes,
               status = "danger",
               fill = TRUE
               )
             )
           ), 
    column(10, 
           bsCollapse(id = "collapse_ex", open = "Parameters", multiple = TRUE,
                      bsCollapsePanel("Parameters w/Sub-Script", uiOutput("gc_params_ui"), uiOutput("gc_textboxes"), style = "success", 
                                      fluidRow(column(12, align = "right", actionButton("fir_update_button", "Update")))), 
                      bsCollapsePanel("Parameters", uiOutput("ui"), style = "info", 
                                      fluidRow(column(12, align = "right", actionButton("sec_update_button", "Update"))))
           )
    )
  )
)

server <- function(input, output, session) {
  pars_gc_names <- reactiveValues(names = NULL)
  fnames_gc_names <- reactiveValues(names = NULL)

  observeEvent(input$parts, {
    parent_node <- input$parts
    all_grandchild_nodes <- switch(parent_node, "fnames" = fnames_names, "pars" = pars_names, "env" = env_names, NULL)
    
    output$ui <- renderUI({
      fluidRow(
        lapply(all_grandchild_nodes, function(node_name) {
          if (parent_node == "env") {
            output$gc_params_ui <- renderUI({})
            output$gc_textboxes <- renderUI({print("Nothing here")})
            if (node_name %in% env_names) {
              node_index <- match(node_name, env_names)
              default_value <- env_values[node_index]
              
              column(6, textInput(node_name, label = node_name, value = default_value))
            }
          } else if (parent_node == "pars") {
            output$gc_params_ui <- renderUI({
              selectInput("gc_params", "Select a Parameter", choices = pars_gc)
            })
            output$gc_textboxes <- renderUI({
              selected_option <- input$gc_params
              
              pars_gc_nodes <- xml_data %>% xml_find_first(paste0(".//pars/leaf/", selected_option))
              pars_gc_names$names <- xml_children(pars_gc_nodes) %>% xml_name()
              pars_gc_values <- xml_children(pars_gc_nodes) %>% xml_text()
              
              lapply(pars_gc_names$names, function(pars_nodes) {
                pars_node_index <- match(pars_nodes, pars_gc_names$names)
                pars_node_value <- pars_gc_values[pars_node_index]
                
                column(6, textInput(pars_nodes, label = pars_nodes, value = pars_node_value))
              }) 
            })
            if (node_name %in% pars_gc) {
              return(NULL)
            }
            if (node_name %in% pars_names) {
              node_index <- match(node_name, pars_names)
              default_value <- pars_values[node_index]
              
              column(6, textInput(node_name, label = node_name, value = default_value))
            }
          } else if (parent_node == "fnames") {
            output$gc_params_ui <- renderUI({
              selectInput("gc_params", "Select a Parameter", choices = fnames_gc)
            })
            output$gc_textboxes <- renderUI({
              selected_option <- input$gc_params
              
              if (!is.null(selected_option)) {
                fnames_gc_nodes <- xml_data %>% xml_find_first(paste0(".//fnames/leaf/", selected_option))
                fnames_gc_names$names <- xml_children(fnames_gc_nodes) %>% xml_name()
                fnames_gc_values <- xml_children(fnames_gc_nodes) %>% xml_text()
                
                lapply(fnames_gc_names$names, function(fnames_nodes) {
                  fnames_node_index <- match(fnames_nodes, fnames_gc_names$names)
                  fnames_node_value <- fnames_gc_values[fnames_node_index]
                  
                  column(6, textInput(fnames_nodes, label = fnames_nodes, value = fnames_node_value))
                })
              } 
            })
            if (node_name %in% fnames_gc) {
              return(NULL)
            }
            if (node_name %in% fnames_names) {
              node_index <- match(node_name, fnames_names)
              default_value <- fnames_values[node_index]
              
              column(6, textInput(node_name, label = node_name, value = default_value))
            }
          }
        })
      )
    })
  })

  observeEvent(input$fir_update_button, {
    parent_node <- input$parts
    updated_xml_data <- xml_data
    selected_option <- input$gc_params
    
    if (parent_node == "fnames") {
      for (node_name in fnames_gc_names$names) {
        new_value <- input[[node_name]]
        node <- xml_find_first(updated_xml_data, paste0(".//fnames/leaf/", selected_option, "/", node_name))
        
        if (!is.null(node)) {
          xml_text(node) <- new_value
        }
      }
    }
    
    if (parent_node == "pars") {
      for (node_name in pars_gc_names$names) {
        new_value <- input[[node_name]]
      
        node <- xml_find_first(updated_xml_data, paste0(".//pars/leaf/", selected_option, "/", node_name))

        if (!is.null(node)) {
          xml_text(node) <- new_value
        }
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

  observeEvent(input$sec_update_button, {
    parent_node <- input$parts
    updated_xml_data <- xml_data

    if (parent_node == "pars") {
      for (node_name in pars_names) {
        if (node_name %in% pars_gc) {next}
        new_value <- input[[node_name]]

        node <- xml_find_first(updated_xml_data, paste0(".//pars/leaf/", node_name))

        if (!is.null(node)) {
          xml_text(node) <- new_value
        }
      }
    }
    if (parent_node == "fnames") {
      for (node_name in fnames_names) {
        if (node_name %in% fnames_gc) {next}
        new_value <- input[[node_name]]

        node <- xml_find_first(updated_xml_data, paste0(".//fnames/leaf/", node_name))

        if (!is.null(node)) {
          xml_text(node) <- new_value
        }
      }
    }
    if (parent_node == "env") {
      for (node_name in env_names) {
        if (node_name %in% env_gc) {next}
        new_value <- input[[node_name]]

        node <- xml_find_first(updated_xml_data, paste0(".//env/leaf/", node_name))

        if (!is.null(node)) {
          xml_text(node) <- new_value
        }
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
}



shinyApp(ui = ui, server = server)
