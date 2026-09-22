library(shiny)
library(xml2)
library(shinyWidgets)
library(magrittr)

# Read the XML file and extract the options
xml_file <- "/Users/fs8/Desktop/Project/MAAT/src/system_models/leaf/leaf_default.xml"
xml_data <- read_xml(xml_file)

# Find the "pars" parent node
pars_leaf_node <- xml_data %>% xml_find_first(".//pars/leaf")

# Find the grandchild node names under "pars"
pars_grandchild_node_names <- xml_children(pars_leaf_node) %>% xml_name()
pars_grandchild_node_values <- xml_children(pars_leaf_node) %>% xml_text()

pars_Topt <- xml_data %>% xml_find_first(".//pars/leaf/Topt")
pars_Topt_names <- xml_children(pars_Topt) %>% xml_name()
pars_Topt_values <- xml_children(pars_Topt) %>% xml_text()

# Define the UI
ui <- fluidPage(
  titlePanel(""),
  sidebarLayout(
    sidebarPanel(
      prettyRadioButtons(
        "parts",
        "Please Select",
        choices = "pars",
        inline = TRUE,
        status = "danger",
        fill = TRUE 
      ),
      width = '10%',
      hr(),
      uiOutput("ui")
    ),
    mainPanel(
      # Add the main panel content here
    )
  )
)

# Define the server
server <- function(input, output, session) {
  
  # Update the UI based on the selected parent node
  observeEvent(input$parts, {
    grandchild_nodes <- pars_grandchild_node_names
    output$ui <- renderUI({
      fluidRow(
        lapply(grandchild_nodes, function(node_name) {
          if (node_name %in% pars_grandchild_node_names) {
            # Replace the condition with the grandchild nodes you want as text inputs or dropdowns
            node_index <- match(node_name, pars_grandchild_node_names)
            default_value <- pars_grandchild_node_values[node_index]
            print(node_name)
            
            column(6, textInput(node_name, label = node_name, value = default_value)
                   )
          }  
        }) 
      )
    })
  })
  
  # Update the input values when the user makes changes
  observeEvent(input$ui, {
    for (node_name in grandchild_nodes) {
      if (node_name %in% pars_grandchild_node_names) {
        updateTextInput(session, node_name, value = input[[node_name]])
      }
    }
  })
}

shinyApp(ui = ui, server = server)