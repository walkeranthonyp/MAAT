library(shiny)
library(shinyFiles)
library(htmltools)
library(shinyhelper)
library(magrittr)
library(shinyjs)
library(xml2)
library(shinyWidgets)
library(shinyBS)
library(XML)
library(glue)

# # Get folder names within the system_models directory
# model_names <- read.csv("/Users/fs8/Desktop/Project/MAAT/src/system_models/model_names.csv")
# 
# # Read the XML file and extract the options
# xml_data <- read_xml("/Users/fs8/Desktop/Project/MAAT/src/system_models/leaf/leaf_default.xml")
# 
# # Find parent nodes
# # parent_nodes <- xml_data %>% xml_children() %>% xml_name() 
# 
# # Find grandchild node names for each parent node
# fnames <- xml_data %>% xml_find_first(".//fnames/leaf")
# fnames_names <- xml_children(fnames) %>% xml_name()
# fnames_values <- xml_children(fnames)  %>% xml_text()
# fnames_options <- read_xml("/Users/fs8/Desktop/Project/MAAT/src/system_models/leaf/leaf_options.xml")
# fnames_gc <- xml_find_all(fnames, ".//*[count(*) > 1]") %>% xml_name()
# 
# pars <- xml_data %>% xml_find_first(".//pars/leaf")
# pars_names <- xml_children(pars) %>% xml_name()
# pars_values <- xml_children(pars) %>% xml_text()
# pars_gc <- xml_find_all(pars, ".//*[count(*) > 1]") %>% xml_name()
# 
# env <- xml_data %>% xml_find_first(".//env/leaf")
# env_names <- xml_children(env) %>% xml_name()
# env_values <- xml_children(env) %>% xml_text()
# env_gc <- xml_find_all(env, ".//*[count(*) > 1]") %>% xml_name()

ui <- fluidPage(
  titlePanel('The Multi-Assumption Architecture and Testbed (MAAT) modelling system'),
  navbarPage(
    id = "maat_menu",
    tags$head(
      tags$style(HTML('.navbar-nav > li > a, .navbar-brand {
                            padding-top:10px !important; 
                            font-size: 23px;
                            padding-bottom:0 !important;
                            height: 45px;
                            }
                           .navbar {min-height:15px !important;}'))
    ),
    tabPanel("Welcome",
             mainPanel(
               style = "font-size: 20px",
               p("The MAAT framework is a", strong("tool"), "tool for comparing different", strong("models"),
                 "of processes in a simple and efficient way. It allows users to run multiple", strong("models"),
                 "with varying assumptions, parameters, and environmental conditions. With built-in", strong("sensitivity analysis"),
                 "and", strong("uncertainty quantification,"), "users can assess the impact of different process representations on 
                 model outputs. The framework also includes a", strong("parameter estimation method"), "called ", strong("Markov Chain 
                 Monte Carlo (MCMC)"), "for fine-tuning model configurations. Overall, MAAT simplifies", strong("model comparison,"),
                 " enables easy ", strong("simulation runs,"), " and provides valuable insights into model ", strong("variability"),
                 " and ", strong("parameter estimation.")), 
               hr(),
               tags$iframe(width = "560", height = "315",
                           src = "https://www.youtube.com/embed/MwDqLkZ90dk",
                           frameborder = "0", allowfullscreen = TRUE)
             )
    ),
    tabPanel("Model Object",
             fluidRow(column(
               5,
               wellPanel(
                 radioButtons(
                   'model_object_buttons',
                   h4('Select a Model Object'),
                   choices = c(""),
                   width = '37%'
                 ) %>% helper(
                   icon = "question",
                   colour = "green",
                   type = "markdown",
                   title = "Model Object Description",
                   content = "model_objs_desc"
                 ),
                 shinyDirButton('directory', 'Select Directory', 'Please Select a Directory')
                 %>% helper(
                   icon = "question",
                   colour = "green",
                   type = "markdown",
                   title = "Model Object Description",
                   content = "test_md"
                 ),
               )
             )),
             fluidRow(column(
               5,
               wellPanel(
                 radioButtons(
                   'change_pars_buttons',
                   h4('Would you like to run the model using the default parameter values?'),
                   choices = c("Yes", "No"),
                   width = '100%'
                 ) %>% helper(
                   icon = "question",
                   colour = "green",
                   type = "markdown",
                   title = "Leaf Parameters Description",
                   content = "run_model_quest_exp"
                 ),
                 conditionalPanel(
                   condition = "input.change_pars_buttons == 'No'",
                   actionButton('page_3', 'Next')
                 ), 
                 conditionalPanel(
                   condition = "input.change_pars_buttons == 'Yes'", 
                   actionButton('run_script', 'Run Script')
                 )
               ), 
             )),
             verbatimTextOutput("selected_directory"),
             verbatimTextOutput("model_object"),
             verbatimTextOutput("leaf_ex")
    ), 
    tabPanel("Parameters",
             titlePanel(""),
             fluidRow(
               column(2, 
                      wellPanel(
                        prettyRadioButtons(
                          "parts",
                          "Please Select",
                          choices = c("fnames", "pars", "env"),
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
  )
)

server <- function(input, output, session) {
  rv <- reactiveValues(directory_path = NULL, selected_model = NULL)
  
  # Get folder names within the system_models directory
  model_names <- read.csv("/Users/fs8/Desktop/Project/MAAT/src/system_models/model_names.csv")
  
  # selected_model <- reactive({input$model_object_buttons})
  # xml_file <- reactive({file.path("/Users/fs8/Desktop/Project/MAAT/src/system_models/", selected_model, "/", selected_model,"_default.xml")})
  # xml_data <- reactive({read_xml(xml_file())})
  
  # Read the XML file and extract the options
  xml_data <- read_xml("/Users/fs8/Desktop/Project/MAAT/src/system_models/leaf/leaf_default.xml")

  # Find parent nodes
  parent_nodes <- xml_data %>% xml_children() %>% xml_name()
  
  # Find grandchild node names for each parent node
  # fnames <- reactive({
  #   xml_data() %>% xml_find_first(".//fnames/leaf")
  # })
  # 
  # fnames_names <- reactive({
  #   xml_children(fnames()) %>% xml_name()
  # })
  # 
  # fnames_values <- reactive({
  #   xml_children(fnames()) %>% xml_text()
  # })
  # 
  # fnames_gc <- reactive({
  #   xml_find_all(fnames(), ".//*[count(*) > 1]") %>% xml_name()
  # })
  
  fnames <- xml_data %>% xml_find_first(".//fnames/leaf")
  fnames_names <- xml_children(fnames) %>% xml_name()
  fnames_values <- xml_children(fnames)  %>% xml_text()
  fnames_options <- read_xml("/Users/fs8/Desktop/Project/MAAT/src/system_models/leaf/leaf_options.xml")
  fnames_gc <- xml_find_all(fnames, ".//*[count(*) > 1]") %>% xml_name()
  
  pars <- xml_data %>% xml_find_first(".//pars/leaf")
  pars_names <- xml_children(pars) %>% xml_name()
  pars_values <- xml_children(pars) %>% xml_text()
  pars_gc <- xml_find_all(pars, ".//*[count(*) > 1]") %>% xml_name()

  env <- xml_data %>% xml_find_first(".//env/leaf")
  env_names <- xml_children(env) %>% xml_name()
  env_values <- xml_children(env) %>% xml_text()
  env_gc <- xml_find_all(env, ".//*[count(*) > 1]") %>% xml_name()

  observe_helpers(help_dir = "/Users/fs8/Desktop/Project/")
  
  pars_gc_names <- reactiveValues(names = NULL)
  fnames_gc_names <- reactiveValues(names = NULL)
  
  # roots <- c(home = '..')
  volumes = getVolumes()()
  shinyDirChoose(input, "directory", roots = volumes, filetypes = NULL)
  
  # Update choices for selectInput0
  updateRadioButtons(session, "model_object_buttons",
                     choices = setNames(model_names$original_name, model_names$display_name))
  
  observeEvent(input$directory, {
    rv$directory_path <- parseDirPath(volumes, input$directory)
    rv$directory_path <- gsub('^/Volumes/Macintosh HD', '', rv$directory_path)
    # print(rv$directory_path)
  })
  
  observeEvent(input$model_object_buttons, {
    rv$selected_model <- input$model_object_buttons
  })
  
  observeEvent(input$run_script, {
    if (!is.null(rv$directory_path) && !is.null(rv$selected_model)) {
      cmd <- sprintf("/Users/fs8/Desktop/Project/MAAT/run_scripts/setup_MAAT_project.bs '%s' '%s'", rv$selected_model, rv$directory_path)
      system(cmd)
    }
  })
  
  observeEvent(input$page_3, {
    updateTabsetPanel(session, "maat_menu", selected = "Parameters")
  })
  
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
              
              if (node_name == "solver") {
                solver <- fnames_options %>% xml_find_first(".//fnames/leaf/solver") %>% xml_text()
                solver <- gsub("[', )]|c\\(", "", unlist(strsplit(solver, ", ")))
                solver <- paste0("'", solver, "'")
                column(6, selectInput(node_name, label = node_name, choices = solver, selected = default_value))
              } else if (node_name == "residual_func") {
                residual_func <- fnames_options %>% xml_find_first(".//fnames/leaf/residual_func") %>% xml_text()
                residual_func <- gsub("[', )]|c\\(", "", unlist(strsplit(residual_func, ", ")))
                residual_func <- paste0("'", residual_func, "'")
                column(6, selectInput(node_name, label = node_name, choices = residual_func, selected = default_value))
              } else if (node_name == "assimilation") {
                assimilation <- fnames_options %>% xml_find_first(".//fnames/leaf/assimilation") %>% xml_text()
                assimilation <- gsub("[', )]|c\\(", "", unlist(strsplit(assimilation, ", ")))
                assimilation <- paste0("'", assimilation, "'")
                column(6, selectInput(node_name, label = node_name, choices = assimilation, selected = default_value))
              } else if (node_name == "semiana") {
                semiana <- fnames_options %>% xml_find_first(".//fnames/leaf/semiana") %>% xml_text()
                semiana <- gsub("[', )]|c\\(", "", unlist(strsplit(semiana, ", ")))
                semiana <- paste0("'", semiana, "'")
                column(6, selectInput(node_name, label = node_name, choices = semiana, selected = default_value))
              } else if (node_name == "Acg") {
                acg <- fnames_options %>% xml_find_first(".//fnames/leaf/Acg") %>% xml_text()
                acg <- gsub("[', )]|c\\(", "", unlist(strsplit(acg, ", ")))
                acg <- paste0("'", acg, "'")
                column(6, selectInput(node_name, label = node_name, choices = acg, selected = default_value))
              } else if (node_name == "Ajg") {
                ajg <- fnames_options %>% xml_find_first(".//fnames/leaf/Ajg") %>% xml_text()
                ajg <- gsub("[', )]|c\\(", "", unlist(strsplit(ajg, ", ")))
                ajg <- paste0("'", ajg, "'")
                column(6, selectInput(node_name, label = node_name, choices = ajg, selected = default_value))
              } else if (node_name == "Apg") {
                apg <- fnames_options %>% xml_find_first(".//fnames/leaf/Apg") %>% xml_text()
                apg <- gsub("[', )]|c\\(", "", unlist(strsplit(apg, ", ")))
                apg <- paste0("'", apg, "'")
                column(6, selectInput(node_name, label = node_name, choices = apg, selected = default_value))
              } else if (node_name == "etrans") {
                etrans <- fnames_options %>% xml_find_first(".//fnames/leaf/etrans") %>% xml_text()
                etrans <- gsub("[', )]|c\\(", "", unlist(strsplit(etrans, ", ")))
                etrans <- paste0("'", etrans, "'")
                column(6, selectInput(node_name, label = node_name, choices = etrans, selected = default_value))
              } else if (node_name == "gas_diff") {
                gas_diff <- fnames_options %>% xml_find_first(".//fnames/leaf/gas_diff") %>% xml_text()
                gas_diff <- gsub("[', )]|c\\(", "", unlist(strsplit(gas_diff, ", ")))
                gas_diff <- paste0("'", gas_diff, "'")
                column(6, selectInput(node_name, label = node_name, choices = gas_diff, selected = default_value))
              } else if (node_name == "Alim") {
                alim <- fnames_options %>% xml_find_first(".//fnames/leaf/Alim") %>% xml_text()
                alim <- gsub("[', )]|c\\(", "", unlist(strsplit(alim, ", ")))
                alim <- paste0("'", alim, "'")
                column(6, selectInput(node_name, label = node_name, choices = alim, selected = default_value))
              } else if (node_name == "vcmax") {
                vcmax <- fnames_options %>% xml_find_first(".//fnames/leaf/vcmax") %>% xml_text()
                vcmax <- gsub("[', )]|c\\(", "", unlist(strsplit(vcmax, ", ")))
                vcmax <- paste0("'", vcmax, "'")
                column(6, selectInput(node_name, label = node_name, choices = vcmax, selected = default_value))
              } else if (node_name == "jmax") {
                jmax <- fnames_options %>% xml_find_first(".//fnames/leaf/jmax") %>% xml_text()
                jmax <- gsub("[', )]|c\\(", "", unlist(strsplit(jmax, ", ")))
                jmax <- paste0("'", jmax, "'")
                column(6, selectInput(node_name, label = node_name, choices = jmax, selected = default_value))
              } else if (node_name == "tcor_jmax") {
                tcor_jmax <- fnames_options %>% xml_find_first(".//fnames/leaf/tcor_jmax") %>% xml_text()
                tcor_jmax <- gsub("[', )]|c\\(", "", unlist(strsplit(tcor_jmax, ", ")))
                tcor_jmax <- paste0("'", tcor_jmax, "'")
                column(6, selectInput(node_name, label = node_name, choices = tcor_jmax, selected = default_value))
              } else if (node_name == "tpu") {
                tpu <- fnames_options %>% xml_find_first(".//fnames/leaf/tpu") %>% xml_text()
                tpu <- gsub("[', )]|c\\(", "", unlist(strsplit(tpu, ", ")))
                tpu <- paste0("'", tpu, "'")
                column(6, selectInput(node_name, label = node_name, choices = tpu, selected = default_value))
              } else if (node_name == "k_pepc") {
                k_pepc <- fnames_options %>% xml_find_first(".//fnames/leaf/k_pepc") %>% xml_text()
                k_pepc <- gsub("[', )]|c\\(", "", unlist(strsplit(k_pepc, ", ")))
                k_pepc <- paste0("'", k_pepc, "'")
                column(6, selectInput(node_name, label = node_name, choices = k_pepc, selected = default_value))
              } else if (node_name == "rd") {
                rd <- fnames_options %>% xml_find_first(".//fnames/leaf/rd") %>% xml_text()
                rd <- gsub("[', )]|c\\(", "", unlist(strsplit(rd, ", ")))
                rd <- paste0("'", rd, "'")
                column(6, selectInput(node_name, label = node_name, choices = rd, selected = default_value))
              } else if (node_name == "rl_rd") {
                rl_rd <- fnames_options %>% xml_find_first(".//fnames/leaf/rl_rd") %>% xml_text()
                rl_rd <- gsub("[', )]|c\\(", "", unlist(strsplit(rl_rd, ", ")))
                rl_rd <- paste0("'", rl_rd, "'")
                column(6, selectInput(node_name, label = node_name, choices = rl_rd, selected = default_value))
              } else if (node_name == "gstar") {
                gstar <- fnames_options %>% xml_find_first(".//fnames/leaf/gstar") %>% xml_text()
                gstar <- gsub("[', )]|c\\(", "", unlist(strsplit(gstar, ", ")))
                gstar <- paste0("'", gstar, "'")
                column(6, selectInput(node_name, label = node_name, choices = gstar, selected = default_value))
              } else if (node_name == "ri") {
                ri <- fnames_options %>% xml_find_first(".//fnames/leaf/ri") %>% xml_text()
                ri <- gsub("[', )]|c\\(", "", unlist(strsplit(ri, ", ")))
                ri <- paste0("'", ri, "'")
                column(6, selectInput(node_name, label = node_name, choices = ri, selected = default_value))
              } else if (node_name == "rs") {
                rs <- fnames_options %>% xml_find_first(".//fnames/leaf/rs") %>% xml_text()
                rs <- gsub("[', )]|c\\(", "", unlist(strsplit(rs, ", ")))
                rs <- paste0("'", rs, "'")
                column(6, selectInput(node_name, label = node_name, choices = rs, selected = default_value))
              } else if (node_name == "rb") {
                rb <- fnames_options %>% xml_find_first(".//fnames/leaf/rb") %>% xml_text()
                rb <- gsub("[', )]|c\\(", "", unlist(strsplit(rb, ", ")))
                rb <- paste0("'", rb, "'")
                column(6, selectInput(node_name, label = node_name, choices = rb, selected = default_value))
              } else if (node_name == "d13c") {
                d13c <- fnames_options %>% xml_find_first(".//fnames/leaf/d13c") %>% xml_text()
                d13c <- gsub("[', )]|c\\(", "", unlist(strsplit(d13c, ", ")))
                d13c <- paste0("'", d13c, "'")
                column(6, selectInput(node_name, label = node_name, choices = d13c, selected = default_value))
              } else if (node_name == "tcor_asc") {
                tcor_asc <- fnames_options %>% xml_find_first(".//fnames/leaf/tcor_asc") %>% xml_text()
                tcor_asc <- gsub("[', )]|c\\(", "", unlist(strsplit(tcor_asc, ", ")))
                tcor_asc <- paste0("'", tcor_asc, "'")
                column(6, selectInput(node_name, label = node_name, choices = tcor_asc, selected = default_value))
              } else if (node_name == "tcor_des") {
                tcor_des <- fnames_options %>% xml_find_first(".//fnames/leaf/tcor_des") %>% xml_text()
                tcor_des <- gsub("[', )]|c\\(", "", unlist(strsplit(tcor_des, ", ")))
                tcor_des <- paste0("'", tcor_des, "'")
                column(6, selectInput(node_name, label = node_name, choices = tcor_des, selected = default_value))
              } else if (node_name == "tcor_dep") {
                tcor_dep <- fnames_options %>% xml_find_first(".//fnames/leaf/tcor_dep") %>% xml_text()
                tcor_dep <- gsub("[', )]|c\\(", "", unlist(strsplit(tcor_dep, ", ")))
                tcor_dep <- paste0("'", tcor_dep, "'")
                column(6, selectInput(node_name, label = node_name, choices = tcor_dep, selected = default_value))
              } else if (node_name == "deltaS") {
                deltaS <- fnames_options %>% xml_find_first(".//fnames/leaf/deltaS") %>% xml_text()
                deltaS <- gsub("[', )]|c\\(", "", unlist(strsplit(deltaS, ", ")))
                deltaS <- paste0("'", deltaS, "'")
                column(6, selectInput(node_name, label = node_name, choices = deltaS, selected = default_value))
              } else if (node_name == "q10") {
                q10 <- fnames_options %>% xml_find_first(".//fnames/leaf/q10") %>% xml_text()
                q10 <- gsub("[', )]|c\\(", "", unlist(strsplit(q10, ", ")))
                q10 <- paste0("'", q10, "'")
                column(6, selectInput(node_name, label = node_name, choices = q10, selected = default_value))
              } else {
                column(6, textInput(node_name, label = node_name, value = default_value))
              }
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

shinyApp(ui, server)

