library(shiny)
library(xml2)
library(shinyWidgets)
library(magrittr)
library(shinyBS)
library(XML)
library(glue)
library(shinyjs)

# Read the XML file and extract the options
xml_data <- read_xml("/Users/fs8/Desktop/Project/MAAT/src/system_models/leaf/leaf_default.xml")
fnames_options <- read_xml("/Users/fs8/Desktop/Project/MAAT/src/system_models/leaf/leaf_options.xml")

# f_options <- xml_find_first(fnames_options, ".//fnames/leaf")
# f_o_names <- xml_children(f_options) %>% xml_name()
# f_o_values <- xml_children(f_options) %>% xml_text()
# f_o_values <- trimws(gsub("^\\s*'c\\(|\\)'\\s*$", "", f_o_values))
# 
# split_values <- strsplit(f_o_values, "', '")
# 
# nodes_with_multiple_values <- list()
# values_with_multiple_values <- list()
# 
# # Filter and store the nodes with more than one value
# for (i in seq_along(split_values)) {
#   if (length(split_values[[i]]) > 1) {
#     nodes_with_multiple_values[[length(nodes_with_multiple_values) + 1]] <- f_o_names[i]
#     formatted_values <- gsub("^'|'$", "", split_values[[i]])  # Remove single quotes from the beginning and end of each value
#     values_with_multiple_values[[length(values_with_multiple_values) + 1]] <- formatted_values
#   }
# }
# 
# f_o_final_names <- unlist(nodes_with_multiple_values)
# f_o_final_values <- unlist(values_with_multiple_values)
# 
# # The new code to read and populate f_o_final_values correctly
# # Create a list of choices for each node
# choices_list <- lapply(strsplit(f_o_values, "', '"), function(x) trimws(sub("^'|'$", "", x)))
# 
# # Filter and store the nodes with more than one value
# nodes_with_multiple_values <- sapply(choices_list, length) > 1
# f_o_final_names <- f_o_names[nodes_with_multiple_values]
# f_o_final_values <- choices_list[nodes_with_multiple_values]

# Find parent nodes
parent_nodes <- xml_data %>% xml_children() %>% xml_name() 

# Find grandchild node names for each parent node
fnames <- xml_data %>% xml_find_first(".//fnames/leaf")
fnames_names <- xml_children(fnames) %>% xml_name()
fnames_values <- xml_children(fnames)  %>% xml_text()
fnames_gc <- xml_find_all(fnames, ".//*[count(*) > 1]") %>% xml_name()

# print(fnames_gc)

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
  
  fnames_rv <- reactiveValues(names = NULL)
  pars_rv <- reactiveValues(names = NULL)
  env_rv <- reactiveValues(names = NULL)
  
  output$gc_params_ui <- renderUI({})
  output$gc_textboxes <- renderUI({})
  
  observeEvent(input$parts, {
    parent_node <- input$parts
    all_grandchild_nodes <- switch(parent_node, "fnames" = fnames_names, "pars" = pars_names, "env" = env_names, NULL)
    
    output$ui <- renderUI({
      fluidRow(
        lapply(all_grandchild_nodes, function(node_name) {
          if (parent_node == "fnames") {
            f_options <- xml_find_first(fnames_options, ".//fnames/leaf")
            f_o_names <- xml_children(f_options) %>% xml_name()
            f_o_values <- xml_children(f_options) %>% xml_text()
            f_o_values <- trimws(gsub("^\\s*'c\\(|\\)'\\s*$", "", f_o_values))
            
            split_values <- strsplit(f_o_values, "', '")

            nodes_with_multiple_values <- list()
            values_with_multiple_values <- list()

            # Filter and store the nodes with more than one value
            for (i in seq_along(split_values)) {
              if (length(split_values[[i]]) > 1) {
                nodes_with_multiple_values[[length(nodes_with_multiple_values) + 1]] <- f_o_names[i]
                formatted_values <- gsub("^'|'$", "", split_values[[i]])  # Remove single quotes from the beginning and end of each value
                values_with_multiple_values[[length(values_with_multiple_values) + 1]] <- formatted_values
              }
            }

            f_o_final_names <- unlist(nodes_with_multiple_values)
            f_o_final_values <- unlist(values_with_multiple_values)
           
            # The new code to read and populate f_o_final_values correctly
            # Create a list of choices for each node
            choices_list <- lapply(strsplit(f_o_values, "', '"), function(x) trimws(sub("^'|'$", "", x)))

            # Filter and store the nodes with more than one value
            nodes_with_multiple_values <- sapply(choices_list, length) > 1
            f_o_final_names <- f_o_names[nodes_with_multiple_values]
            f_o_final_values <- choices_list[nodes_with_multiple_values]

            if (node_name %in% fnames_gc) {
              return(NULL)
            }
            if (node_name %in% fnames_names) {
            node_index <- match(node_name, fnames_names)
            default_value <- fnames_values[node_index]

            if (node_name %in% f_o_final_names) {
              node_index <- match(node_name, f_o_final_names)
              dd_default_values <- f_o_final_values[[node_index]]

              default_value <- gsub("'", "", default_value)
              default_value_single <- toString(default_value)
              column(6, pickerInput(node_name, label = node_name, choices = dd_default_values, selected = default_value, options = list(`actions-box` = TRUE), multiple = TRUE))
            } else {
              disable("ui")
              column(6, textInput(node_name, label = node_name, value = if (is.null(fnames_rv[[node_name]]) || fnames_rv[[node_name]] == default_value) default_value else fnames_rv[[node_name]]))
              }
            }
          }
          else if (parent_node == "pars") {
            # print(input[[node_name]])
            if (node_name %in% pars_gc) {
              return(NULL)
            }
            # print(input[[node_name]])
            if (node_name %in% pars_names) {
              node_index <- match(node_name, pars_names)
              default_value <- pars_values[node_index]
              
              column(6, textInput(node_name, label = node_name, value = default_value))
            }
          } 
          else if (parent_node == "env") {
            if (node_name %in% env_names) {
              node_index <- match(node_name, env_names)
              default_value <- env_values[node_index]
              
              column(6, textInput(node_name, label = node_name, value = if (is.null(input[[node_name]]) || input[[node_name]] == default_value) default_value else input[[node_name]]))
            }
          }
        })
      )
    })
  })
  
  observeEvent(input$parts, {
    parent_node <- input$parts
    
    output$gc_params_ui <- renderUI({
      if (parent_node == "fnames") {
        selectInput("gc_params", "Select a Parameter", choices = fnames_gc)
      } else if (parent_node == "pars") {
        selectInput("gc_params", "Select a Parameter", choices = pars_gc)
      }
    })
    
    output$gc_textboxes <- renderUI({
      selected_option <- input$gc_params
      
      if (parent_node == "fnames" && !is.null(selected_option)) {
        fnames_gc_nodes <- xml_data %>% xml_find_first(paste0(".//fnames/leaf/", selected_option))
        fnames_gc_names$names <- xml_children(fnames_gc_nodes) %>% xml_name()
        fnames_gc_values <- xml_children(fnames_gc_nodes) %>% xml_text()
        
        f_options <- fnames_options %>% xml_find_first(paste0(".//fnames/leaf/", selected_option))
        f_o_names <- xml_text(f_options)
        f_o_values <- trimws(gsub("^\\s*'c\\(|\\)'\\s*$", "", f_o_names))
        
        choices_list <- strsplit(f_o_values, ", ")
        nodes_with_multiple_values <- sapply(choices_list, length) > 1
        f_o_final_names <- fnames_gc_names$names[nodes_with_multiple_values]
        f_o_final_values <- choices_list[nodes_with_multiple_values]
        
        # Ensure both lists have the same length
        max_len <- max(length(f_o_final_names), length(f_o_final_values))
        f_o_final_names <- rep(f_o_final_names, length.out = max_len)
        f_o_final_values <- rep(f_o_final_values, length.out = max_len)
        
        # Create a list to store the pickerInput elements
        input_list <- lapply(seq_along(f_o_final_names), function(i) {
          fnames_node_name <- f_o_final_names[i]
          xpath_expr <- paste0(".//*[local-name() = '", fnames_node_name, "']")
          default_value <- fnames_gc_values[fnames_gc_names$names == fnames_node_name]
          
          # Use fnames_node_name as the input ID for each pickerInput
          print(f_o_final_values[1])
          column(6, pickerInput(
            fnames_node_name, 
            label = fnames_node_name,
            choices = f_o_final_values[[i]], 
            selected = default_value, options = list(`actions-box` = TRUE),
            multiple = TRUE
          ))
        })
        
        # Combine all the pickerInput elements into a single output
        do.call(tagList, input_list)
      } else if (parent_node == "pars" && !is.null(selected_option)) {
        pars_gc_nodes <- xml_data %>% xml_find_first(paste0(".//pars/leaf/", selected_option))
        pars_gc_names$names <- xml_children(pars_gc_nodes) %>% xml_name()
        pars_gc_values <- xml_children(pars_gc_nodes) %>% xml_text()
        
        lapply(pars_gc_names$names, function(pars_nodes) {
          pars_node_index <- match(pars_nodes, pars_gc_names$names)
          pars_node_value <- pars_gc_values[pars_node_index]
          
          column(6, textInput(pars_nodes, label = pars_nodes, value = pars_node_value))
        })
      }
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
          # Check if the new_value is a string and add single quotes if necessary
          if (is.character(new_value) && !grepl("^'.*'$", new_value)) {
            new_value <- paste0("'", new_value, "'")
          }
          xml_text(node) <- new_value
        }
      }
    }
    
    if (parent_node == "pars") {
      for (node_name in pars_gc_names$names) {
        new_value <- input[[node_name]]
        
        node <- xml_find_first(updated_xml_data, paste0(".//pars/leaf/", selected_option, "/", node_name))
        
        if (!is.null(node)) {
          # Check if the new_value is a string and add single quotes if necessary
          if (is.character(new_value) && !grepl("^'.*'$", new_value)) {
            new_value <- paste0("'", new_value, "'")
          }
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
        default_value <- pars_values[match(node_name, pars_names)]
        
        if (!is.null(new_value) && new_value != default_value) { # new_value is check whether it's not NULL and it differes from default_value
          node <- xml_find_first(updated_xml_data, paste0(".//pars/leaf/", node_name))
          
          if (!is.null(node)) {
            xml_text(node) <- new_value
          }
        }
      }
    }
    
    if (parent_node == "fnames") {
      for (node_name in fnames_names) {
        if (node_name %in% fnames_gc) {next}
        new_value <- input[[node_name]]
        default_value <- fnames_values[match(node_name, fnames_names)]

        if (!is.null(new_value) && new_value != default_value) { # new_value is check whether it's not NULL and it differes from default_value
          node <- xml_find_first(updated_xml_data, paste0(".//fnames/leaf/", node_name))
          
          if (!is.null(node)) {
            xml_text(node) <- new_value
          }
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
    write_xml(updated_xml_data, new_xml_file, declaration = TRUE)
    
    showModal(modalDialog(
      title = "Update Successful",
      "The XML file has been updated.",
      easyClose = TRUE
    ))
  })
}



shinyApp(ui = ui, server = server)
