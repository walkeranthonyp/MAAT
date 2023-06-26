library(shiny)
library(xml2)
library(shinyWidgets)
library(magrittr)
library(shinyBS)

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
grandchild_nodes <- xml_children(pars_leaf_node)
pars_ggc_names <- xml_find_all(pars_leaf_node, ".//*[count(*) > 1]") %>% xml_name()

rt_leaf_node <- xml_data %>% xml_find_first(".//pars/leaf/reftemp")
print(rt_leaf_node)

env_leaf_node <- xml_data %>% xml_find_first(".//env/leaf")
env_grandchild_node_names <- xml_children(env_leaf_node) %>% xml_name()
env_grandchild_node_values <- xml_children(env_leaf_node) %>% xml_text()


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
           bsCollapse(id = "collapse_ex", open = "Panel 2", 
                           bsCollapsePanel("Child Parameters", uiOutput("ui"), style = "info"), 
                           bsCollapsePanel("Grand Child Parameters", 
                                      selectInput("grandchild_params", "Select Parameter", choices = pars_ggc_names),
                                      uiOutput("grandchild_textboxes"),
                                      style = "success")  
                      )
           )
    )
)

# Define the server
server <- function(input, output, session) {
  # Update the UI based on the selected parent node
  observeEvent(input$parts, {
    parent_node <- input$parts
    grandchild_nodes <- switch(parent_node,
                               "fnames" = fnames_grandchild_node_names,
                               "pars" = pars_grandchild_node_names,
                               "env" = env_grandchild_node_names,
                               NULL)
    # updateUI(session, grandchild_nodes)
  
  # Dynamically generate UI based on the selected parent node
  output$ui <- renderUI({
    parent_node <- input$parts
    all_grandchild_nodes <- switch(parent_node,
                                   "fnames" = fnames_grandchild_node_names,
                                   "pars" = pars_grandchild_node_names,
                                   "env" = env_grandchild_node_names,
                                   NULL)
    fnames_excluded_nodes <- c("tcor_asc", "tcor_des", "tcor_dep", "deltaS", "q10")
    pars_excluded_nodes <- c("reftemp", "atref", "Ha", "Hd", "Topt", "deltaS", "a_deltaS_t", "b_deltaS_t", "c_deltaS_t", "q10", "a_q10_t", "b_q10_t", "tupp_cox", "tlow_cox", "exp_cox")
    fluidRow(
      column(6,
             lapply(all_grandchild_nodes[1:ceiling(length(all_grandchild_nodes)/3)], function(node_name) {
               if (parent_node == "env") {
                 if (node_name %in% env_grandchild_node_names) {
                   # Replace the condition with the grandchild nodes you want as text inputs or dropdowns
                   node_index <- match(node_name, env_grandchild_node_names)
                   default_value <- env_grandchild_node_values[node_index]
                   if (node_name == "o2_conc") {
                     sliderInput(node_name, label = node_name, value = default_value, min = 0.0, max = 1.0, step = 0.1)
                   } else {
                     textInput(node_name, label = node_name, value = default_value)
                   }
                 } else {
                   # For other grandchild nodes, display as labels
                   tags$label(node_name)
                 }
               } else if (parent_node == "pars") {
                 if (node_name %in% pars_excluded_nodes) {
                   return(NULL)
                 }
                 if (node_name %in% pars_grandchild_node_names) {
                   # Replace the condition with the grandchild nodes you want as text inputs or dropdowns
                   node_index <- match(node_name, pars_grandchild_node_names)
                   default_value <- pars_grandchild_node_values[node_index]
                   
                   textInput(node_name, label = node_name, value = default_value)
                 } else {
                   # For other grandchild nodes, display as labels
                   tags$label(node_name)
                 }
               } else if (parent_node == "fnames") {
                 if (node_name %in% fnames_excluded_nodes) {
                   return(NULL)
                 }
                 if (node_name %in% fnames_grandchild_node_names) {
                   # Replace the condition with the grandchild nodes you want as text inputs or dropdowns
                   node_index <- match(node_name, fnames_grandchild_node_names)
                   default_value <- fnames_grandchild_node_values[node_index]
                   
                   if (node_name == "solver") {
                     choices <- leaf_xml_options %>% xml_find_first(".//fnames/leaf/solver")
                     # print(choices)
                     solver_options <- c('f_solver_analytical_leaf_c4_r0', 'f_solver_analytical_leaf_no_r', 
                                         'f_solver_analytical_leaf_quad', 'f_solver_analytical_leaf_quad_r0', 
                                         'f_solver_analytical_leaf_simple', 'f_solver_brent', 
                                         'f_solver_brent_diag', 'f_solver_semiana_leaf_Ar')
                     selectInput(node_name, label = node_name, choices = solver_options, selected = default_value)
                   } else {
                     textInput(node_name, label = node_name, value = default_value)
                   }
                 } else {
                   # For other grandchild nodes, display as labels
                   tags$label(node_name)
                 }
               }
             })
      ),
      column(6,
             lapply(all_grandchild_nodes[ceiling(length(all_grandchild_nodes)/2):length(all_grandchild_nodes)], function(node_name) {
               if (parent_node == "env") {
                 if (node_name %in% env_grandchild_node_names) {
                   # Replace the condition with the grandchild nodes you want as text inputs or dropdowns
                   node_index <- match(node_name, env_grandchild_node_names)
                   default_value <- env_grandchild_node_values[node_index]
                   if (node_name == "o2_conc") {
                     sliderInput(node_name, label = node_name, value = default_value, min = 0.0, max = 1.0, step = 0.1)
                   } else {
                     textInput(node_name, label = node_name, value = default_value)
                   }
                 } else {
                   # For other grandchild nodes, display as labels
                   tags$label(node_name)
                 }
               } else if (parent_node == "pars") {
                 if (node_name %in% pars_excluded_nodes) {
                   return(NULL)
                 }
                 if (node_name %in% pars_grandchild_node_names) {
                   # Replace the condition with the grandchild nodes you want as text inputs or dropdowns
                   node_index <- match(node_name, pars_grandchild_node_names)
                   default_value <- pars_grandchild_node_values[node_index]
                   
                   textInput(node_name, label = node_name, value = default_value)
                 } else {
                   # For other grandchild nodes, display as labels
                   tags$label(node_name)
                 }
               } else if (parent_node == "fnames") {
                 if (node_name %in% fnames_excluded_nodes) {
                   return(NULL)
                 }
                 if (node_name %in% fnames_grandchild_node_names) {
                   # Replace the condition with the grandchild nodes you want as text inputs or dropdowns
                   node_index <- match(node_name, fnames_grandchild_node_names)
                   default_value <- fnames_grandchild_node_values[node_index]
                   
                   if (node_name == "solver") {
                     choices <- leaf_xml_options %>% xml_find_first(".//fnames/leaf/solver")
                     # print(choices)
                     solver_options <- c('f_solver_analytical_leaf_c4_r0', 'f_solver_analytical_leaf_no_r', 
                                         'f_solver_analytical_leaf_quad', 'f_solver_analytical_leaf_quad_r0', 
                                         'f_solver_analytical_leaf_simple', 'f_solver_brent', 
                                         'f_solver_brent_diag', 'f_solver_semiana_leaf_Ar')
                     selectInput(node_name, label = node_name, choices = solver_options, selected = default_value)
                   } else {
                     textInput(node_name, label = node_name, value = default_value)
                   }
                 } else {
                   # For other grandchild nodes, display as labels
                   tags$label(node_name)
                 }
               }
             })
      )
    )
  })
  
  output$grandchild_textboxes <- renderUI({
    grandchild_param <- input$grandchild_params # grandchild_param gets the selected parameter from the dropdown menu
    
    grandchild_node_values <- pars_get_grandchild_node_values(grandchild_param)  # Function to get the values based on the selected parameter
    
    if (!is.null(grandchild_node_values)) {
      lapply(names(grandchild_node_values), function(node_name) {
        column(6, textInput(node_name, label = node_name, value = grandchild_node_values[[node_name]]))
      })
    } else {
      NULL
    }
  })
  })
  
  pars_get_grandchild_node_values <- function(param) {
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
  fnames_get_grandchild_node_values <- function(fnames) {
  switch(fnames,
         "tcor_asc" = list(vcmax = "f_tcor_asc_Arrhenius",
                           jmax = "f_tcor_asc_Arrhenius",
                           tpu = "f_tcor_asc_Arrhenius",
                           k_pepc = "f_tcor_asc_Arrhenius",
                           rd = "f_tcor_asc_Arrhenius",
                           gstar = "f_tcor_asc_quadratic_bf1985",
                           tau = "f_tcor_asc_Q10",
                           Kc = "f_tcor_asc_Arrhenius",
                           Ko = "f_tcor_asc_Arrhenius"),
         "tcor_des" = list(vcmax = "f_tcor_des_modArrhenius",
                           jmax = "f_tcor_des_modArrhenius",
                           tpu = "f_tcor_des_modArrhenius",
                           k_pepc = "f_tcor_des_modArrhenius",
                           rd = "f_scalar_none"),
         "tcor_dep" = list(tpu = "f_tcor_dep_independent",
                           rd = "f_tcor_dep_independent",
                           tau = "f_tcor_dep_independent"),
         "deltaS" = list(rd = "f_deltaS",
                         vcmax = "f_deltaS",
                         jmax = "f_deltaS",
                         tpu = "f_deltaS",
                         k_pepc = "f_deltaS"),
         "q10" = list(rd = "f_q10_constant",
                      vcmax = "f_q10_constant",
                      jmax = "f_q10_constant",
                      k_pepc = "f_q10_constant",
                      tau = "f_q10_constant",
                      Kc = "f_q10_constant",
                      Ko = "f_q10_constant"),
         NULL
  )
}

  
  # Update the input values when the user makes changes
  observeEvent(input$ui, {
    parent_node <- input$parts
    all_grandchild_nodes <- switch(parent_node,
                                   "fnames" = fnames_grandchild_node_names,
                                   "pars" = pars_grandchild_node_names,
                                   "env" = env_grandchild_node_names,
                                   NULL)
    for (node_name in all_grandchild_nodes) {
      if (node_name %in% env_grandchild_node_names) {
        updateTextInput(session, node_name, value = input[[node_name]])
      }
    }
  })
}


shinyApp(ui = ui, server = server)
