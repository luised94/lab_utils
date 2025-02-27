library(shiny)

# Create a simple plot
plot(1:10, 1:10)

# Create a shiny app using the plot
shinyApp(
  ui = basicPage(
    plotOutput("plot")
  ),
  server = function(input, output) {
    output$plot <- renderPlot({
      plot(1:10, 1:10)
    })
  }
)
