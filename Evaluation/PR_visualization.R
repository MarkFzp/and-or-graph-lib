args = commandArgs(trailingOnly=TRUE)
if(!require("plotly")){
  install.packages("plotly", repos="http://cran.rstudio.com/")
  library("plotly")
}

if(!require("gdata")){
  install.packages("gdata", repos="http://cran.rstudio.com/")
  library("gdata")
}

MCTS = read.xls (args[1], sheet = 1, header = TRUE)
Benchmark = read.xls(args[1], sheet = 2, header = TRUE)

trace1 <- list(
  x = MCTS$recall,
  y = MCTS$bleu,
  mode = "markers",
  showlegend = FALSE,
  text = MCTS$setting,
  type = "scatter",
  color = "black"
)

trace2 <- list(
  x = Benchmark$Recall,
  y = Benchmark$Bleu,
  mode = "markers",
  showlegend = FALSE,
  text = Benchmark$Setting,
  type = "scatter",
  color = "blue"
)

p <- plot_ly()
p <- add_trace(p, x=trace1$x, y=trace1$y, mode = trace1$mode, showlegend=trace1$showlegend, text=trace1$text, type=trace1$type) %>% 
add_lines(x = MCTS$recall, y = MCTS$bleu,
          line = list(color = '#07A4B5'), name = "MCTS", showlegend = TRUE)

p <- add_trace(p, x=trace2$x, y=trace2$y, mode = trace2$mode, showlegend=trace2$showlegend, text=trace2$text, type=trace2$type) %>% 
  add_lines(x = Benchmark$Recall, y = Benchmark$Bleu, 
            line = list(color = '#07A4B5'), name = "Benchmark", showlegend = TRUE)

htmlwidgets::saveWidget(as_widget(p), paste(args[1], ".html"))



