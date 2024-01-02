library(shiny)
library(plotly)
library(dplyr)
library(stringr)

data <- read.csv2("data.csv")
names(data) <- str_replace_all(names(data), '\\.', ' ')
Xs <- data %>% select(-Result)
Ys <- data$Result
rm(data)

ui <- fluidPage(theme = bslib::bs_theme(version = 4, bootswatch = "minty"), 
                tags$head(
                  tags$style(HTML("

     .multicol {

       -webkit-column-count: 3; /* Chrome, Safari, Opera */

       -moz-column-count: 3; /* Firefox */

       column-count: 3;

     }
     
   "))
                  
                ),
    title = "Phenotype-genotype association in patients with TMA and kidney involvement",
    
    fluidRow(align = "center",
      column(width = 8, 
             fluidRow( h3("Plot")),
             fluidRow(plotlyOutput(outputId = "pca", height = "600px"))
      ),
      column(width = 4, align = "left",
             fluidRow(h3("Parameters")),
             fluidRow(div("Use the following terms:")),
             fluidRow(uiOutput("checkbox")),
             fluidRow(
               selectInput(inputId = "lim", label = "Use terms occuring at least this time: ", choices = c(3:100), selected = 5),
               selectInput(inputId = "dim", label = "Number of dimensions: ", choices = c(2,3), selected = 3)
             )
      )
    ),
    
   HTML('<footer class="footer" align = "center"> 
          Patients with a variant explicative of the renal phenotype are depicted in red. Intended for exploratory purposes only. Source: <a href="">Lombardi et al., submitted.</a>
        </footer>')
    
    
)

#### SERVER SIDE ####
server <- function(input, output, session) 
{
  
  components <- reactive({
    req(input$checkbox)
    X <- Xs %>% select(all_of(input$checkbox))
    X_df <- X
    nr <- nrow(X_df)
    if(nr > 0)
    {
      X_df <- colSums(X_df)/nr
    }
    else
    {
      return(NULL)
    }
    X_idf <- log(1/X_df)
    
    X_tf <- X
    X_tf <- X_tf %>% mutate_all(as.integer)
    X_nterms <- rowSums(X_tf)
    for(i in 1:nrow(X_tf))
    {
      nt <- X_nterms[i]
      if(nt > 0)
      {
        X_tf[i,] <- X_tf[i,]/nt
      }
      else
      {
        X_tf[i,] <- 0
      }
    }
    
    for(i in 1:ncol(X_tf))
    {
      X_tf[,i] <- X_tf[,i]*X_idf[i]
    }
    string <- c()
    

    for(i in 1:nrow(X_tf))
    {
      curdf <- X_tf[i,]
      sums <- colSums(curdf)
      curnames <- paste(names(curdf[1,sums != 0]), collapse = "\n")
      string[i] <- paste(ifelse(!is.na(Ys[i]), paste('Result: variant in ', Ys[i], '\n', sep = ""), 'Result: no causal variant\n'), curnames, sep = "")
    }
    prin_comp <- prcomp(X_tf, rank. = as.numeric(input$dim[1]))
    
    components <- as.data.frame(prin_comp[["x"]])
    components$string <- string
    
    return(components)
  }) 
  
  output$checkbox <- renderUI({
    choice <- names(Xs[,colSums(Xs) >= as.numeric(input$lim[1])])
    wellPanel(
      
      tags$div(class = "multicol", checkboxGroupInput("checkbox", label = "", choices = choice, selected = choice))
      
    )
    
  }) 
  
  output$pca <- renderPlotly({
    components <- components() 
    
    if(ncol(components) == 4)
    {
      fig <- plot_ly(components, x = ~PC1, y = ~PC2, z = ~PC3, color = ~(!is.na(Ys)), colors = c('#636EFA','#EF553B') , type = "scatter3d", mode = "markers",
                     showlegend = F)
    }
    else if(ncol(components) == 3)
    {
      fig <- plot_ly(components, x = ~PC1, y = ~PC2, color = ~(!is.na(Ys)), colors = c('#636EFA','#EF553B'), type = "scatter", mode = "markers",
                     showlegend = F)
    }
    else
    {
      return(NULL)
    }
    
    fig <- fig %>%
      layout(
        scene = list(bgcolor = "#e5ecf6")
      ) %>%
      add_trace(
        text = components$string,
        hoverinfo = 'text'
      )
    
    fig
  })
}

shinyApp(ui = ui, server = server)