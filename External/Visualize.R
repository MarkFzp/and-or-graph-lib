args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
   stop("argument must be supplied (input file)", call.=FALSE)
} 

if(!require("visNetwork")){
  install.packages("visNetwork", repos="http://cran.rstudio.com/")
  library("visNetwork")
}

if(!require("shiny")){
  install.packages("shiny", repos = "http://cran.rstudio.com/")
  library("shinydashboard")
}

if(!require("shinydashboard")){
  install.packages("shinydashboard", repos = "http://cran.rstudio.com/")
}
#devtools::install_github("datastorm-open/visNetwork")

library(visNetwork)
example <- read.csv(args[1], header=FALSE, sep=",", stringsAsFactors = FALSE)
names(example) <- c("From", "To", "SourceContent", "ChildContent", "Weight", "Order", "Level")


Vertex <- unique(c((example[,1]), unique(example[,2])))

Contents <- rep(NA, length(Vertex))
for(i in 1:length(Vertex)){
  for(j in 1:length(example[,1])){
    if(Vertex[i] == example[,1][j]){
        Contents[i] <- example$SourceContent[j]
        break
    }
  }

  #if(!is.na(StateIDs[i]))
   # next

  for(k in 1:length(example[,2])){
    if(Vertex[i] == example[,2][k]){
        Contents[i] <- example$ChildContent[k]
    }
  }
}


Vertex_df <- data.frame(Vertex, Contents, stringsAsFactors = FALSE)
Link_df <- data.frame(example$From, example$To, example$Weight, example$Order, example$Level, stringsAsFactors = FALSE)
names(Link_df) <- c("From", "To", "Weight", "Order", "Level")

#Divide Nodes by Group
Vertex_df$Group <- NA
for(i in 1:length(Vertex_df$Vertex)){
  pos <- which(Link_df$From == Vertex_df$Vertex[i])[1]
  if(is.na(Link_df$Weight[pos])){
    Vertex_df$Group[i] <- "AndNode"
  }else{
    Vertex_df$Group[i] <- "OrNode"
  }
}

Vertex_df$Level <- NA
for(i in 1:length(Vertex_df$Vertex)){
  for(j in 1:length(Link_df$Level)){
    if(Vertex_df$Vertex[i] == Link_df$From[j])
      Vertex_df$Level[i] <- Link_df$Level[j]
    if(Vertex_df$Vertex[i] == Link_df$To[j])
      Vertex_df$Level[i] <- Link_df$Level[j] + 1
  }
}

Link_df$Order <- Link_df$Order
Vertex_df$Order <- NA
Vertex_df$Order[1] <- 1
for(i in 1:length(Vertex_df$Vertex)){
  for(j in 1:length(Link_df$Order)){
    if(Vertex_df$Vertex[i] == Link_df$To[j])
      Vertex_df$Order[i] <- Link_df$Order[j]
  }
}


# find leaf nodes
isLeafNode <- rep(NA, length(Vertex_df$Vertex))
for(i in 1:length(Vertex_df$Vertex)){
  if((!Vertex_df$Vertex[i] %in% Link_df$From) && (Vertex_df$Vertex[i] %in% Link_df$To)){
    isLeafNode[i] <- TRUE}else{
      isLeafNode[i] <- FALSE
    }
}
Vertex_df$Group[isLeafNode] <- "LeafNode"


#Divides edges by dash lines
Link_df$Dash <- NA
for(i in 1:length(Link_df$Weight)){
  if(!is.na(Link_df$Weight[i]))
    Link_df$Dash[i] <- TRUE
  else
    Link_df$Dash[i] <- FALSE
}

# Full Visualization
Vertex_df$label <- NA
for(i in 1:length(Vertex_df$Vertex)){
  if(Vertex_df$Group[i] != "LeafNode")
    Vertex_df$label[i] <- Vertex_df$Vertex[i]
  else{
    Vertex_df$label[i] <- Vertex_df$Contents[i]
  }
}

Vertex_df$hover_label <- NA
for(i in 1:length(Vertex_df$Vertex)){
  if(Vertex_df$Group[i] != "LeafNode")
    Vertex_df$hover_label[i] <- Vertex_df$Contents[i]
  else{
    Vertex_df$hover_label[i] <- Vertex_df$Vertex[i]
  }
}

nodes <- data.frame(id = Vertex_df[,1], group = Vertex_df$Group, label = Vertex_df$label, hover_label = Vertex_df$hover_label, level=Vertex_df$Level, x=Vertex_df$Order)
edges <- data.frame(from = Link_df$From, to = Link_df$To, label = Link_df$Weight,
    length=50,
    width=1,
    dashes = Link_df$Dash
    )


Visualize <- function(nodes, edges){
  visualization <- visNetwork(nodes, edges, width = "100%") %>%
    visGroups(groupname= "AndNode", color="darkblue", shape = "triangle") %>%
    visGroups(groupname= "OrNode", color="red", shape = "circle") %>%
    visGroups(groupname= "LeafNode", color="orange", shape = "square") %>%
    visInteraction(hover=TRUE, hoverConnectedEdge=TRUE) %>%
    visOptions(highlightNearest = list(enabled = T, degree = 1, hover = T)) %>%
    visEdges(arrows = "to") %>%
    visHierarchicalLayout(levelSeparation = 200) %>% 
    visInteraction(hover=T) %>% 
    visEvents(hoverNode  = "function(e){
            var label_info = this.body.data.nodes.get({
     fields: ['label', 'hover_label'],
     filter: function (item) {
     return item.id === e.node
     },
     returnType :'Array'
});
     this.body.data.nodes.update({id: e.node, label : label_info[0].hover_label, hover_label : label_info[0].label});
     }") %>% 
  visEvents(blurNode  = "function(e){
            var label_info = this.body.data.nodes.get({
            fields: ['label', 'hover_label'],
            filter: function (item) {
            return item.id === e.node
            },
            returnType :'Array'
            });
            this.body.data.nodes.update({id: e.node, label : label_info[0].hover_label, hover_label : label_info[0].label});
  }")

  visualization
  return (visualization)
}


visualization_full <- Visualize(nodes, edges)
visualization_full
length(Link_df$From)
sum(Vertex_df$Group == "LeafNode")


############################################################
### Construct visualization for each Dummy Nodes recursively:
############################################################
# rootVertexID <- Vertex_df$Vertex[1]
# RootToDummy <- Link_df[which(Link_df$From == rootVertexID),]
# DummyNodes <- data.frame(Vertex = integer(), Contents = character(), Group = character(), Level = integer(), Order=integer())
# for(i in 1:length(Vertex_df$Vertex)){
#   if(Vertex_df$Vertex[i] %in% RootToDummy$To){
#     DummyNodes <- rbind(DummyNodes, Vertex_df[i,])
#   }
# }

DummyNodesIDs <- Vertex_df[Vertex_df$Group == "AndNode" | Vertex_df$Group == "OrNode",]$Vertex


# FindSubVertex <- function(vertexID, SubGraph_Vertex_df){
#   SubGraph_Vertex_df <- rbind(SubGraph_Vertex_df, Vertex_df[Vertex_df$Vertex == vertexID,])
#   for(i in 1:length(Link_df$From)){
#     if(vertexID == Link_df$From[i]){
#       SubGraph_Vertex_df <- FindSubVertex(Link_df$To[i], SubGraph_Vertex_df)
#     }
#   }
#   return(SubGraph_Vertex_df)
# }
# 
# FindSubLinks <- function(subVertex, subGraph_Link_df){
#   for(vertex_from in subVertex$Vertex){
#     for(i in 1:length(Link_df$From)){
#       if(Link_df$From[i] == vertex_from){
#         for(vertex_to in subVertex$Vertex){
#           if(Link_df$To[i] == vertex_to){
#             subGraph_Link_df <- rbind(subGraph_Link_df, Link_df[i,])
#           }
#         }
#       }
#     }
#   }
# 
#   return (subGraph_Link_df)
# }


FindSubset <- function(vertexID, Subset){
  temp_Link_df <- data.frame(From = integer(), To = integer(), SourceContent = character(), ChildContent = character(), Weight = double(), Order= integer(), Level = integer())
  
  for(i in 1:length(example$From)){
    if(vertexID == example$From[i]){
      # check if we should add it to the subset:
      addOrNot <- TRUE
      if(length(Subset$From) != 0){
        for(j in 1:length(Subset$From)){
          if(Subset$From[j] == example$From[i] && Subset$To[j] == example$To[i]){
            addOrNot <- FALSE
          }
        }
      }
      
      if(addOrNot)
        temp_Link_df <- rbind(temp_Link_df, example[i,])
    }
  }
  Subset <- rbind(Subset, temp_Link_df)
  DeeperLinks <- temp_Link_df$To
  if(length(DeeperLinks) != 0){
    for(i in 1:length(DeeperLinks)){
      Subset <- FindSubset(DeeperLinks[i], Subset)
    }
  }
  return (Subset)
}


#################################################
#### Try Multiple Plots using Shiny App
#################################################
library(shiny)
library(visNetwork)
ui <- dashboardPage(
  dashboardHeader(
    title = "Visualizations",
    titleWidth = 350
  ),
  
  dashboardSidebar(
    width = 100,
    menuItem("Visualization", tabName= "Visualization", icon = icon("dashboard"))
  ),
  
  dashboardBody(
    tabItems(
      tabItem(tabName = "Visualization",
        fluidRow(
          box(width = 12, visNetworkOutput("full_model"))        
        ),

        
        fluidRow(
          box(width = 2,
            selectInput("var",
                        label = "Left",
                        choices = DummyNodesIDs),
            selectInput("var2",
                        label = "Right",
                        choices = DummyNodesIDs),
            align = "top"
          ),
          
          box(width=5,
            visNetworkOutput("part_vis"),
            align="top"
          ),
          
          box(width = 5,
            visNetworkOutput("part_vis2"),
            align="top"
          )
        )
      )
    )
)
)

#################################################################
### Data Preparation to be fed to server:
#################################################################
for(vertexID in DummyNodesIDs){
  Subset_ <- data.frame(From = integer(), To = integer(), SourceContent = character(), ChildContent = character(), Weight = double(), Order= integer(), Level = integer())
  Subset <- FindSubset(vertexID,Subset_)
  Subset <- unique(Subset)
  Subset_final <- data.frame(From = integer(), To = integer(), SourceContent = character(), ChildContent = character(), Weight = double(), Order= integer(), Level = integer())
  for(i in 1:length(Subset$From)){
    if(!(Subset[i,]$From %in% Subset_final$From) && !(Subset[i,]$To %in% Subset_final$To)){
      Subset_final <- rbind(Subset_final, Subset[i,])
    }
  }
  
  
  Vertex <- unique(c((Subset[,1]), unique(Subset[,2])))
  
  Contents <- rep(NA, length(Vertex))
  for(i in 1:length(Vertex)){
    for(j in 1:length(Subset[,1])){
      if(Vertex[i] == Subset[,1][j]){
        Contents[i] <- Subset$SourceContent[j]
        break
      }
    }
    
    #if(!is.na(StateIDs[i]))
    # next
    
    for(k in 1:length(Subset[,2])){
      if(Vertex[i] == Subset[,2][k]){
        Contents[i] <- Subset$ChildContent[k]
      }
    }
  }
  
  
  Sub_Vertex_df <- data.frame(Vertex, Contents,stringsAsFactors = FALSE)
  Sub_Link_df <- data.frame(Subset$From, Subset$To, Subset$Weight, Subset$Order, Subset$Level,stringsAsFactors = FALSE)
  names(Sub_Link_df) <- c("From", "To", "Weight", "Order", "Level")
  
  #Divide Nodes by Group
  Sub_Vertex_df$Group <- NA
  for(i in 1:length(Sub_Vertex_df$Vertex)){
    pos <- which(Sub_Link_df$From == Sub_Vertex_df$Vertex[i])[1]
    if(is.na(Sub_Link_df$Weight[pos])){
      Sub_Vertex_df$Group[i] <- "AndNode"
    }else{
      Sub_Vertex_df$Group[i] <- "OrNode"
    }
  }
  
  Sub_Vertex_df$Level <- NA
  for(i in 1:length(Sub_Vertex_df$Vertex)){
    if((Sub_Vertex_df$Vertex[i] %in% Sub_Link_df$From) && !(Sub_Vertex_df$Vertex[i] %in% Sub_Link_df$To)){
        Sub_Vertex_df$Level[i] <- Sub_Link_df[Sub_Link_df$From == Sub_Vertex_df$Vertex[i], ]$Level[1]
    }
  }
  
  for(i in 1:length(Sub_Vertex_df$Vertex)){
    if(is.na(Sub_Vertex_df[i,]$Level)){
      parentID <- Sub_Link_df[Sub_Link_df$To == Sub_Vertex_df$Vertex[i], ]$From[1]
      Sub_Vertex_df$Level[i] <- Sub_Vertex_df[Sub_Vertex_df$Vertex == parentID, ]$Level[1] + 1
    }
  }
    

  #Order
  Sub_Link_df$Order <- Sub_Link_df$Order
  Sub_Vertex_df$Order <- NA
  Sub_Vertex_df$Order[1] <- 1
  for(i in 1:length(Sub_Vertex_df$Vertex)){
    for(j in 1:length(Sub_Link_df$Order)){
      if(Sub_Vertex_df$Vertex[i] == Sub_Link_df$To[j])
        Sub_Vertex_df$Order[i] <- Sub_Link_df$Order[j]
    }
  }
  
  
  # find leaf nodes
  isLeafNode <- rep(NA, length(Sub_Vertex_df$Vertex))
  for(i in 1:length(Sub_Vertex_df$Vertex)){
    if((!Sub_Vertex_df$Vertex[i] %in% Sub_Link_df$From) && (Sub_Vertex_df$Vertex[i] %in% Sub_Link_df$To)){
      isLeafNode[i] <- TRUE}else{
        isLeafNode[i] <- FALSE
      }
  }
  Sub_Vertex_df$Group[isLeafNode] <- "LeafNode"
  
  
  #Divides edges by dash lines
  Sub_Link_df$Dash <- NA
  for(i in 1:length(Sub_Link_df$Weight)){
    if(!is.na(Sub_Link_df$Weight[i]))
      Sub_Link_df$Dash[i] <- TRUE
    else
      Sub_Link_df$Dash[i] <- FALSE
  }
  
  #Set labels and hover_labels
  Sub_Vertex_df$label <- NA
  for(i in 1:length(Sub_Vertex_df$Vertex)){
    if(Sub_Vertex_df$Group[i] != "LeafNode")
      Sub_Vertex_df$label[i] <- Sub_Vertex_df$Vertex[i]
    else{
      Sub_Vertex_df$label[i] <- Sub_Vertex_df$Contents[i]
    }
  }
  
  Sub_Vertex_df$hover_label <- NA
  for(i in 1:length(Sub_Vertex_df$Vertex)){
    if(Sub_Vertex_df$Group[i] != "LeafNode")
      Sub_Vertex_df$hover_label[i] <- Sub_Vertex_df$Contents[i]
    else{
      Sub_Vertex_df$hover_label[i] <- Sub_Vertex_df$Vertex[i]
    }
  }
  
  assign(paste("Sub_nodes_", vertexID, sep = ""), data.frame(id = Sub_Vertex_df[,1], group = Sub_Vertex_df$Group, label = Sub_Vertex_df$label, hover_label = Sub_Vertex_df$hover_label, level=Sub_Vertex_df$Level, x=Sub_Vertex_df$Order))
  assign(paste("Sub_edges_", vertexID, sep = ""), data.frame(from = Sub_Link_df$From, to = Sub_Link_df$To, label = Sub_Link_df$Weight, length=50, width=1, dashes = Sub_Link_df$Dash))
}

##################################################################
server <- function(input, output) {
   output$full_model <- renderVisNetwork({
     visualization_full
   })
   
   output$part_vis <- renderVisNetwork({
      node_name <- paste("Sub_nodes_", input$var,sep = "")
      Nodes <- get(node_name)
      
      edge_name <- paste("Sub_edges_", input$var, sep="")
      Edges <- get(edge_name)
      
      visNetwork(Nodes, Edges, width = "100%") %>%
        visGroups(groupname= "AndNode", color="darkblue", shape = "triangle") %>%
        visGroups(groupname= "OrNode", color="red", shape = "circle") %>%
        visGroups(groupname= "LeafNode", color="orange", shape = "square") %>%
        visHierarchicalLayout() %>%
        visInteraction(hover=TRUE, hoverConnectedEdge=TRUE) %>%
        visOptions(highlightNearest = list(enabled = T, degree = 1, hover = T)) %>%
        visEdges(arrows = "to") %>%
        visHierarchicalLayout(levelSeparation = 200) %>% 
      visInteraction(hover=T) %>% 
        visEvents(hoverNode  = "function(e){
         var label_info = this.body.data.nodes.get({
         fields: ['label', 'hover_label'],
         filter: function (item) {
         return item.id === e.node
         },
         returnType :'Array'
         });
         this.body.data.nodes.update({id: e.node, label : label_info[0].hover_label, hover_label : label_info[0].label});
   }") %>% 
  visEvents(blurNode  = "function(e){
            var label_info = this.body.data.nodes.get({
            fields: ['label', 'hover_label'],
            filter: function (item) {
            return item.id === e.node
            },
            returnType :'Array'
            });
            this.body.data.nodes.update({id: e.node, label : label_info[0].hover_label, hover_label : label_info[0].label});
  }")
   })
   
   output$part_vis2 <- renderVisNetwork({
     node_name <- paste("Sub_nodes_", input$var2,sep = "")
     Nodes <- get(node_name)
     
     edge_name <- paste("Sub_edges_", input$var2, sep="")
     Edges <- get(edge_name)
     
     visNetwork(Nodes, Edges, width = "100%") %>%
       visGroups(groupname= "AndNode", color="darkblue", shape = "triangle") %>%
       visGroups(groupname= "OrNode", color="red", shape = "circle") %>%
       visGroups(groupname= "LeafNode", color="orange", shape = "square") %>%
       visHierarchicalLayout() %>%
       visInteraction(hover=TRUE, hoverConnectedEdge=TRUE) %>%
       visOptions(highlightNearest = list(enabled = T, degree = 1, hover = T)) %>%
       visEdges(arrows = "to") %>%
       visHierarchicalLayout(levelSeparation = 200) %>% 
       visInteraction(hover=T) %>% 
       visEvents(hoverNode  = "function(e){
                 var label_info = this.body.data.nodes.get({
                 fields: ['label', 'hover_label'],
                 filter: function (item) {
                 return item.id === e.node
                 },
                 returnType :'Array'
                 });
                 this.body.data.nodes.update({id: e.node, label : label_info[0].hover_label, hover_label : label_info[0].label});
   }") %>% 
  visEvents(blurNode  = "function(e){
            var label_info = this.body.data.nodes.get({
            fields: ['label', 'hover_label'],
            filter: function (item) {
            return item.id === e.node
            },
            returnType :'Array'
            });
            this.body.data.nodes.update({id: e.node, label : label_info[0].hover_label, hover_label : label_info[0].label});
  }")
   })
}

shinyApp(ui = ui, server = server)
#############################################
##Further Improvements:
#############################################
#@@@ Change label of a node when hover: 
#https://github.com/datastorm-open/visNetwork/issues/146

#@@@ List all And Nodes and their children












