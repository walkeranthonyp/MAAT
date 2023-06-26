library(shiny)
library(shinyFiles)
library(htmltools)
library(shinyhelper)
library(magrittr)

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
               p(
                 "The MAAT framework is a",
                 strong("tool"),
                 "tool for comparing different",
                 strong("models"),
                 "of processes in a simple and efficient way. It allows
          users to run multiple",
                 strong("models"),
                 "with varying assumptions,
          parameters, and environmental conditions. With built-in",
                 strong("sensitivity analysis"),
                 "and",
                 strong("uncertainty quantification,"),
                 "users can assess the impact of different process representations on model outputs.
          The framework also includes a",
                 strong("parameter estimation method"),
                 "called ",
                 strong("Markov Chain Monte Carlo (MCMC)"),
                 "for fine-tuning model configurations.
          Overall, MAAT simplifies",
                 strong("model comparison,"),
                 " enables easy ",
                 strong("simulation runs,"),
                 " and provides valuable insights into model ",
                 strong("variability"),
                 " and ",
                 strong("parameter estimation.")
               ), 
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
                   actionButton('next_button', 'Next')
                 ), 
                 conditionalPanel(
                   condition = "input.change_pars_buttons == 'Yes'", 
                   actionButton('run_script', 'Run Script')
                 )
               )
             )),
             verbatimTextOutput("selected_directory"),
             verbatimTextOutput("model_object"),
             verbatimTextOutput("leaf_ex")
    )
  )
)

server <- function(input, output, session) {
  observe_helpers(help_dir = "/Users/fs8/Desktop/Project/")
  
  rv <- reactiveValues(directory_path = NULL, selected_model = NULL)
  
  # roots <- c(home = '..')
  volumes = getVolumes()()
  shinyDirChoose(input, "directory", roots = volumes, filetypes = NULL)
  
  # Get folder names within the system_models directory
  model_names <- read.csv("/Users/fs8/Desktop/Project/MAAT/src/system_models/model_names.csv")
  
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
      cmd <- sprintf("/Users/fs8/Desktop/Project/MAAT/run_scripts/setup_MAAT_project.bs '%s' '%s'",
                     rv$selected_model, rv$directory_path)
      system(cmd)
    }
  })
  
  observeEvent(input$next_button, {
    if (!is.null(rv$directory_path) && !is.null(rv$selected_model)) {
      
    }
  })
}

shinyApp(ui, server)

