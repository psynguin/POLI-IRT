library(shiny)
library(shinyjs)
library(ggplot2)
library(dplyr)
library(tidyr)
library(bslib)

# Helper function for highly compressed synced inputs
param_input <- function(id, label, min, max, value) {
  tags$div(class = "param-container mb-2",
    tags$div(class = "d-flex justify-content-between align-items-end mb-3",
      tags$label(label, class = "form-label fw-bold mb-0", style = "font-size: 0.85rem;"),
      tags$div(style = "width: 70px;",
        numericInput(paste0(id, "_num"), NULL, value = value, min = min, max = max, step = 0.01)
      )
    ),
    sliderInput(id, NULL, min = min, max = max, value = value, step = 0.01)
  )
}

# UI Definition
ui <- page_navbar(
  title = "Polytomous IRT Item Category Response Curves",
  id = "model",
  theme = bs_theme(version = 5, preset = "shiny"),
  
  # Injecting CSS to remove scrollbar, compress margins, and shrink inputs
  tags$head(
    tags$style(HTML("
      .bslib-sidebar-layout > .sidebar { overflow-y: hidden !important; padding-top: 10px !important; padding-bottom: 10px !important; }
      .param-container .form-group { margin-bottom: 0 !important; }
      .param-container input[type='number'] { height: 26px; padding: 2px 6px; font-size: 0.85rem; }
      .param-container .irs { margin-top: -8px; margin-bottom: -10px; }
      .help-block { margin-bottom: 5px !important; font-size: 0.8rem; }
      hr { margin: 10px 0; }
    "))
  ),
  
  sidebar = sidebar(
    width = 350,
    useShinyjs(),
    
    helpText(textOutput("model_help"), class = "help-block"),
    hr(),
    
    param_input("a", "Discrimination (a)", 0.01, 4.0, 1.50),
    param_input("b1", "Threshold 1 (b₁)", -4.0, 4.0, -1.50),
    param_input("b2", "Threshold 2 (b₂)", -4.0, 4.0, 0.00),
    param_input("b3", "Threshold 3 (b₃)", -4.0, 4.0, 1.50)
  ),
  
  nav_panel("GRM", value = "GRM", card(plotOutput("irtPlot_grm", height = "550px"))),
  nav_panel("PCM", value = "PCM", card(plotOutput("irtPlot_pcm", height = "550px")))
)

# Server Logic
server <- function(input, output, session) {
  
  output$model_help <- renderText({
    if (input$model == "GRM") {
      "GRM enforces strictly ordered thresholds (b1 < b2 < b3)."
    } else {
      "PCM allows disordered step parameters (reversals) but assumes a fixed discrimination of 1."
    }
  })
  
  # Standardized bi-directional sync (avoids infinite loops and rounding overwrites)
  sync_param <- function(id) {
    num_id <- paste0(id, "_num")
    
    observeEvent(input[[id]], {
      req(is.numeric(input[[id]]), is.numeric(input[[num_id]]))
      if (input[[id]] != input[[num_id]]) {
        updateNumericInput(session, num_id, value = input[[id]])
      }
    }, ignoreInit = TRUE)
    
    observeEvent(input[[num_id]], {
      req(is.numeric(input[[num_id]]), is.numeric(input[[id]]))
      if (input[[num_id]] != input[[id]]) {
        updateSliderInput(session, id, value = input[[num_id]])
      }
    }, ignoreInit = TRUE)
  }
  
  sync_param("a")
  sync_param("b1")
  sync_param("b2")
  sync_param("b3")
  
  # Helper to update both inputs on constraint check
  update_val <- function(id, val) {
    updateSliderInput(session, id, value = val)
    updateNumericInput(session, paste0(id, "_num"), value = val)
  }
  
  # Handle model switching behavior
  observeEvent(input$model, {
    if (input$model == "PCM") {
      shinyjs::disable("a")
      shinyjs::disable("a_num")
      update_val("a", 1.0)
    } else {
      shinyjs::enable("a")
      shinyjs::enable("a_num")
      if (isTRUE(input$b1_num >= input$b2_num) || isTRUE(input$b2_num >= input$b3_num)) {
        update_val("b1", -1.5)
        update_val("b2", 0.0)
        update_val("b3", 1.5)
      }
    }
  })
  
  # Enforce GRM threshold order constraints using numeric inputs as the source of truth
  observeEvent(input$b1_num, {
    req(input$model == "GRM", is.numeric(input$b1_num), is.numeric(input$b2_num))
    if (input$b1_num >= input$b2_num) update_val("b1", input$b2_num - 0.01)
  })
  
  observeEvent(input$b2_num, {
    req(input$model == "GRM", is.numeric(input$b1_num), is.numeric(input$b2_num), is.numeric(input$b3_num))
    if (input$b2_num <= input$b1_num) update_val("b2", input$b1_num + 0.01)
    if (input$b2_num >= input$b3_num) update_val("b2", input$b3_num - 0.01)
  })
  
  observeEvent(input$b3_num, {
    req(input$model == "GRM", is.numeric(input$b2_num), is.numeric(input$b3_num))
    if (input$b3_num <= input$b2_num) update_val("b3", input$b2_num + 0.01)
  })
  
  # Central plot reactive expression
  plot_expr <- reactive({
    if (input$model == "GRM") {
      req(is.numeric(input$b1_num), is.numeric(input$b2_num), is.numeric(input$b3_num))
      req(input$b1_num < input$b2_num, input$b2_num < input$b3_num)
    }
    
    a <- if (input$model == "PCM") 1.0 else input$a_num
    b1 <- input$b1_num
    b2 <- input$b2_num
    b3 <- input$b3_num
    
    req(is.numeric(a), is.numeric(b1), is.numeric(b2), is.numeric(b3))
    
    df_plot <- tibble(theta = seq(-6, 6, by = 0.05)) %>%
      mutate(
        model_type = input$model,
        
        p_star2 = 1 / (1 + exp(-a * (theta - b1))),
        p_star3 = 1 / (1 + exp(-a * (theta - b2))),
        p_star4 = 1 / (1 + exp(-a * (theta - b3))),
        
        num1 = 1,
        num2 = exp(theta - b1),
        num3 = exp((theta - b1) + (theta - b2)),
        num4 = exp((theta - b1) + (theta - b2) + (theta - b3)),
        den  = num1 + num2 + num3 + num4,
        
        P1 = if_else(model_type == "GRM", pmax(0, 1 - p_star2), num1 / den),
        P2 = if_else(model_type == "GRM", pmax(0, p_star2 - p_star3), num2 / den),
        P3 = if_else(model_type == "GRM", pmax(0, p_star3 - p_star4), num3 / den),
        P4 = if_else(model_type == "GRM", pmax(0, p_star4), num4 / den)
      ) %>%
      select(theta, P1, P2, P3, P4) %>%
      pivot_longer(cols = P1:P4, names_to = "Category", values_to = "Probability")
    
    ggplot(df_plot, aes(x = theta, y = Probability, color = Category)) +
      geom_line(linewidth = 1) +
      geom_vline(xintercept = b1, color = "#475569", linetype = "dashed") +
      geom_vline(xintercept = b2, color = "#475569", linetype = "dashed") +
      geom_vline(xintercept = b3, color = "#475569", linetype = "dashed") +
      annotate("text", x = b1 - 0.15, y = 0.95, label = "b[1]", parse = TRUE, color = "#475569", angle = 90) +
      annotate("text", x = b2 - 0.15, y = 0.95, label = "b[2]", parse = TRUE, color = "#475569", angle = 90) +
      annotate("text", x = b3 - 0.15, y = 0.95, label = "b[3]", parse = TRUE, color = "#475569", angle = 90) +
      geom_hline(yintercept = 0.5, color = "#94a3b8", linetype = "dashed") +
      scale_color_manual(values = c(
        "P1" = "#2563EB", 
        "P2" = "#F59E0B", 
        "P3" = "#10B981", 
        "P4" = "#D97706"
      )) +
      scale_x_continuous(breaks = seq(-6, 6, 2)) +
      scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
      labs(
        x = expression("Latent Trait (" * theta * ")"),
        y = expression("Probability P(" * theta * ")")
      ) +
      theme_minimal(base_size = 14) +
      theme(
        legend.position = "top",
        legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(color = "#334155"),
        axis.title = element_text(color = "#0F172A", face = "bold")
      )
  })
  
  output$irtPlot_grm <- renderPlot({ plot_expr() })
  output$irtPlot_pcm <- renderPlot({ plot_expr() })
}

shinyApp(ui, server)