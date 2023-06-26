library(shiny)
library(xml2)
library(shinyWidgets)
library(magrittr)
library(shinyBS)
library(XML)

# Read the XML file and extract the options
xml_file <- "/Users/fs8/Desktop/Project/MAAT/src/system_models/leaf/leaf_default.xml"
xml_data <- read_xml(xml_file)

# Read the leaf_options XML file and extract the options 
leaf_options <- "/Users/fs8/Desktop/Project/MAAT/src/system_models/leaf/leaf_options.xml"
leaf_xml_options <- read_xml(leaf_options)
# xml_structure(leaf_xml_options)

# Find parent nodes
top_parent <- xml_data %>% xml_find_first(".")
parent_nodes <- xml_data %>% xml_children() %>% xml_name() 

# Find grandchild node names for each parent node
fnames_leaf_node <- xml_data %>% xml_find_first(".//fnames/leaf")
fnames_grandchild_node_names <- xml_children(fnames_leaf_node) %>% xml_name()
fnames_grandchild_node_values <- xml_children(fnames_leaf_node)  %>% xml_text()

pars_leaf_node <- xml_data %>% xml_find_first(".//pars/leaf")
pars_grandchild_node_names <- xml_children(pars_leaf_node) %>% xml_name()
pars_grandchild_node_values <- xml_children(pars_leaf_node) %>% xml_text()

pars_ggc_names <- xml_find_all(pars_leaf_node, ".//*[count(*) > 1]") %>% xml_name()

env_leaf_node <- xml_data %>% xml_find_first(".//env/leaf")
env_grandchild_node_names <- xml_children(env_leaf_node) %>% xml_name()
env_grandchild_node_values <- xml_children(env_leaf_node) %>% xml_text()


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
                                      selectInput("grandchild_params", "Select Parameter", choices = pars_ggc_names),
                                      uiOutput("grandchild_textboxes"),
                                      style = "success"),
                      bsCollapsePanel("Parameters", uiOutput("ui"), style = "info")
           )
    )
  )
)

server <- function(input, output, session) {
  observeEvent(input$parts, {
    pars_leaf_node <- pars_grandchild_node_names
    
    output$ui <- renderUI({
      fluidRow(
        lapply(pars_leaf_node, function(node_name) {
          if (node_name %in% pars_ggc_names) {
            return(NULL)
          }
          if (node_name %in% pars_grandchild_node_names) {
            node_index <- match(node_name, pars_grandchild_node_names)
            default_value <- pars_grandchild_node_values[node_index]
            
            column(6, textInput(node_name, label = node_name, value = default_value))
          }  
        }) 
      )
    })
    
    output$grandchild_textboxes <- renderUI({
      grandchild_param <- input$grandchild_params # grandchild_param gets the selected parameter from the dropdown menu
      
      reftemp <- xml_data %>% xml_find_first(paste0(".//pars/leaf/", grandchild_param))
      reftemp_names <- xml_children(reftemp) %>% xml_name()
      reftemp_values <- xml_children(reftemp) %>% xml_text()
      
      # grandchild_node_values <- get_grandchild_node_values(grandchild_param)  # Function to get the values based on the selected parameter
      
      pars_ggc <- reftemp_names
      
      lapply(pars_ggc, function(pars_node_name) {
        if (pars_node_name %in% reftemp_names) {
          pars_node_index <- match(pars_node_name, reftemp_names)
          pars_def_value <- reftemp_values[pars_node_index]
          
          column(6, textInput(pars_node_name, label = pars_node_name, value = pars_def_value))
        }
      })
    })
  })
  
  get_grandchild_node_values <- function(param) {
    switch(param,
           "reftemp" = list(rd = 25, vcmax = 25, jmax = 25, tpu = 25, k_pepc = 25, Kc = 25, Ko = 25, gstar = 25, tau = 25),
           "atref" = list(rd = 2, vcmax = 50, jmax = 100, tpu = 5, k_pepc = 7e+05, Kc = 40.49, Ko = 27.84, gstar = 4.325, tau = 2.6, vomax = 0),
           "Ha" = list(rd = 69830, vcmax = 69830, jmax = 100280, tpu = 69830, k_pepc = 69830, Kc = 79430, Ko = 36380, gstar = 37830, tau = -41572, vomax = 60110),
           "Hd" = list(rd = 2e+05, vcmax = 2e+05, jmax = 2e+05, tpu = 2e+05, k_pepc = 2e+05),
           "Topt" = list(rd = 27.56, vcmax = 27.56, jmax = 19.89, tpu = 27.56, k_pepc = 27.56),
           "deltaS" = list(rd = 0, vcmax = 0, jmax = 0, tpu = 0, k_pepc = 0),
           "a_deltaS_t" = list(rd = 490, vcmax = 668, jmax = 660, tpu = 485),
           "b_deltaS_t" = list(rd = 0, vcmax = -1.07, jmax = -0.75, tpu = 0),
           "c_deltaS_t" = list(rd = 0, vcmax = 0, jmax = -0.52, tpu = 0),
           "q10" = list(rd = 2, vcmax = 2, jmax = 2, tpu = 2, k_pepc = 2, Kc = 2, Ko = 2, tau = 0.57),
           "tupp_cox" = list(vcmax = 36, rd = 45),
           "tlow_cox" = list(vcmax = 0, rd = 5),
           "exp_cox" = list(vcmax = 0.3, rd = 0.4),
           NULL
    )
  }
  
  # Update the input values when the user makes changes
  observeEvent(input$ui, {
    for (node_name in pars_leaf_node) {
      if (node_name %in% pars_grandchild_node_names) {
        updateTextInput(session, node_name, value = input[[node_name]])
      }
    }
  })
}



shinyApp(ui = ui, server = server)
