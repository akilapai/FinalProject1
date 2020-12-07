#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
load("FinalData.rda")

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  navbarPage("Final Project",
             tabPanel("Gabrielle Legaspi", 
                      fluidRow(selectInput("predictor", label = "Predictor", 
                                           choices = c("Time", "Treatment", "Gender")),
                               plotOutput("gabrielle"),
                               tableOutput("names"))),
             
             tabPanel("Girish Valluru", 
                      fluidRow(sliderInput("integer", "# of bins",
                                           min = 4, max = 500,
                                           value = 25),
                               plotOutput("girish"),
                               tableOutput("values"))),
             
             tabPanel("Zelong Liu",
                      fluidRow(sidebarPanel(radioButtons("Time","Please select the Time", 
                                                         choices = c("week0","week4"))),plotOutput("distPlot"))),
             
             tabPanel("Akila Pai",
                      fluidRow(sidebarLayout(position = "right",
                                             sidebarPanel(h3("Inputs for graph"), 
                                                          selectInput("Gender", "1. Select Gender", 
                                                                      choices = c("F" = "F","M" = "M"), selected = "F")),
                                             plotOutput("mygraph")))),
             
             tabPanel("Parker King", 
                      fluidRow(titlePanel("BMI & Treatment"),
                               sidebarLayout(sidebarPanel(sliderInput(inputId = "bins",
                                                                      label = "Number of bins:",
                                                                      min = 1, max = 50,value = 30)),
                                             splitLayout(cellWidths = c("50%", "50%"), 
                                                         plotOutput("plotgraph1"), plotOutput("plotgraph2"))))),
             tabPanel("Aislyn DiRisio", 
                      fluidRow(titlePanel("Comparing Expression of Proteins"),
                               sidebarLayout(position = "right", sidebarPanel(selectInput("ProteinName", "Protein", choices=unique(FinalData$ProteinName), selected = "TNFB")),
                                             mainPanel(plotOutput("aislyn_plot"), textOutput("aislyn_title"), textOutput("aislyn_drugA"), textOutput("aislyn_drugB"))))),
             tabPanel("Jane Brennan",
                      titlePanel(title="Protein vs. Expression"),
                      fluidRow(selectInput("Age", label = "Age",
                                           choices = c("Less than or Equal to 40", "Over 40")),
                               plotOutput("jane"))),
             tabPanel("Wendy Huang", 
                      fluidRow(selectInput("gt", label = "Gender and Treatment", 
                                           choices = c("Males on Drug A", "Males on Drug B","Females on Drug A","Females on Drug B")),
                               plotOutput("wendy"))),
             tabPanel("James Landefeld",
                      sidebarLayout(
                        sidebarPanel = sidebarPanel(
                          sliderInput(inputId = "SIZE.JL", label = "Number of models to include per protein",
                                      min = 1, max = 100, value = 1, step = 5),
                          
                          selectInput(inputId = 'protein_number.JL', label = 'Select a Protein',
                                      choices = unique(FinalData[['ProteinNumber']])),
                          
                          radioButtons(inputId = "MOD.JL", label = "Select Predicter",
                                       choices = c("Age" = "Age", "bmi" = "bmi", "Age*bmi" = "Age*bmi"),
                                       selected = "Age", inline = T)),
                        
                        mainPanel = mainPanel(tabsetPanel(tabPanel("All proteins", 
                                                                   plotOutput(outputId = "AllProt.JL"))),
                                              tabsetPanel(tabPanel(title='Selected Protein',
                                                                   plotOutput(outputId = "SingleProt.JL"))))))
             
             ))) 

