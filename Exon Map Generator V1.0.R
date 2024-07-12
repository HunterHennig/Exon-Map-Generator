# Install packages as needed
list_of_packages = c("shiny","xlsx","tidyverse","stringr","tidyr","dplyr","devEMF","plotrix","ggplot2","shinythemes","modelr", "colourpicker")

lapply(list_of_packages, 
       function(x) if(!require(x,character.only = TRUE)) install.packages(x))

# Load packages
library(shiny)
library(xlsx)
library(tidyverse)
library(stringr)
library(tidyr)
library(dplyr)
library(devEMF)
library(plotrix)
library(ggplot2)
library(colourpicker)
library(shinythemes)
library(modelr)

unlink("Files", recursive = TRUE)
dir.create("Files")
unlink("data.zip", recursive = TRUE)

# Define UI for dataset viewer app ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("Exon Map Generator"),
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Create user input options
      textInput(inputId = "app_gene_name",
                label = "Gene Name:",
                value = ""),
      
      textInput(inputId = "app_sequence",
                label = "Sequence:",
                value = ""),
      
      actionButton(inputId = "save",
                   label = "Save"),
      
      downloadButton(outputId = "downloadData",
                     label = "Download"),
      
      tags$h2("Customizations"),
      
      checkboxInput(inputId = "app_show_length",
                    label = "Show Exon Length",
                    value = TRUE),
      
      checkboxInput(inputId = "app_above_exon",
                    label = "Display Length Above The Exon",
                    value = FALSE),
      
      textInput(inputId = "app_intron_length",
                label = "Intron Length (default 0):",
                value = "0"),
      
      sliderInput(inputId = "app_exon_height",
                  label = "Exon Height:",
                  min = 0,
                  max = 100,
                  value = 50,
                  step = 0.1),
      
      sliderInput(inputId = "app_text_size",
                  label = "Text Size:",
                  min = 0.1,
                  max = 1.5,
                  value = 1,
                  step = 0.1),
      
      sliderInput(inputId = "app_phase_shift_size",
                  label = "Phase Shift Circle Size:",
                  min = 0,
                  max = 200,
                  value = 100,
                  step = 0.1),
      
      colourInput(inputId = "exon_color", 
                  label = "Select Exon Color:", 
                  value = "#47F6FF"),
      
      colourInput(inputId = "ps0_color", 
                  label = "Select Phase Shift 0 Color:", 
                  value = "#8EFF84"),
      
      colourInput(inputId = "ps1_color", 
                  label = "Select Phase Shift 1 Color:", 
                  value = "#FFEF35"),
      
      colourInput(inputId = "ps2_color", 
                  label = "Select Phase Shift 2 Color:", 
                  value = "#FF6055"),
      
      colourInput(inputId = "ps_bkgd_color",
                  label = "Select Phase Shift Background Color:",
                  value = "#353638"),
      
      actionButton(inputId = "resetcolor",
                   label = "Reset Color Options"),
      
      actionButton(inputId = "resetall",
                   label = "Reset All"),
      
    ),
    
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Plot preview + List of saved emf files ----
      plotOutput("plot"),
      textOutput("txt"),
      span(textOutput("txt"), style="font-size:15px; font-family:arial; font-style:bold"),
      
    )
  )
)

# Define server logic to summarize and view selected dataset ----
server <- function(input, output, session) {
  
  # Reset all selectable values
  observeEvent(input$resetall, 
               {
                 
                 updateTextInput(session, "app_gene_name", value = "")
                 updateTextInput(session, "app_sequence", value = "")
                 updateTextInput(session, "app_intron_length", value = "0")
                 updateSliderInput(session, "app_exon_height", value = 50)
                 updateSliderInput(session, "app_phase_shift_size", value = 100)
                 updateColourInput(session, "exon_color", value = "#47F6FF")
                 updateColourInput(session, "ps0_color", value = "#8EFF84")
                 updateColourInput(session, "ps1_color", value = "#FFEF35")
                 updateColourInput(session, "ps2_color", value = "#FF6055")
                 updateColourInput(session, "ps_bkgd_color", value = "#353638")
                 
               })
  
  # Reset color values
  observeEvent(input$resetcolor, 
               {
                 
                 updateColourInput(session, "exon_color", value = "#47F6FF")
                 updateColourInput(session, "ps0_color", value = "#8EFF84")
                 updateColourInput(session, "ps1_color", value = "#FFEF35")
                 updateColourInput(session, "ps2_color", value = "#FF6055")
                 updateColourInput(session, "ps_bkgd_color", value = "#353638")
                 
               })
  
  # Set up the text output
  txtout <- reactiveValues(outputText = "")
  txtout$txt <- "Maps Saved: "
  
  observeEvent(input$save,
               {
                 # Create an emf file in the "Files" folder for later download
                 plotSave(input$app_intron_length, input$app_sequence, input$ps_bkgd_color, input$exon_color, input$ps0_color, input$ps1_color, input$ps2_color)
                 
                 # Add name of saved plot to 'Maps Saved' list if it does not already exist (prevents a name being added multiple times to the list as there would only be one file output per name)
                 if(grepl(paste("^", as.character(input$app_gene_name)), txtout$txt) != TRUE){
                   
                   txtout$txt <- paste(txtout$txt, isolate(as.character(input$app_gene_name)), ", ", sep = "") 
                   
                 }
                 
               })
  
  # Create function to name the saved emf file
  filename <- function(){
    paste(isolate(as.character(input$app_gene_name)), " Exon Map ", ".emf", sep = "")
  }
  
  # Create function to retrieve and format the system date
  date <- function(){
    temp_date <- Sys.Date()
    return(format(temp_date, format = "%m-%d-%Y"))
  } 
  
  # Create function to save the drawn plot
  plotSave = function(x, y, b, e, in_ps0, in_ps1, in_ps2){
    
    intron_length <- as.numeric(x)
    
    upper_app <- gsub("\\s+", "", as.character(y))
    upper_app <- gsub("a|c|g|t", " ", upper_app)
    upper_app <- gsub("^\\s+|\\s+$", "", upper_app)
    upper_app <- gsub("\\s+", " ", upper_app)
    
    str_replace(upper_app, "\\s{2,}", "")
    
    df1_app <- data.frame(upper_app)
    temp_app <- separate_longer_delim(df1_app, cols = 1, delim = ' ')
    temp_app <- data.frame(NA, NA, str_count(temp_app$upper_app))
    colnames(temp_app) <- c("Start", "End", "Length")
    
    values <- y
    
    # Uses Length values of previous row to get Start values
    counter <- 2
    newdata <- c((intron_length))
    while(counter <= nrow(temp_app))
    {
      newdata <- append(newdata, (temp_app$Length[(counter - 1)] + newdata[(counter - 1)] + (intron_length)))
      counter <- counter + 1
    }
    
    temp_app$Start <- newdata
    
    # Use Start values to get End values of previous row
    counter2 <- 1
    newdata2 <- c()
    while(counter2 <= nrow(temp_app)){
      newdata2 <- append(newdata2, (temp_app$Start[(counter2)] + temp_app$Length[(counter2)]))
      counter2 <- counter2 + 1
    }
    
    temp_app$End <- newdata2
    
    # Add phase shift column
    ps <- c()
    counter4 <- 1
    
    # Calculate phase shifts and add to dataframe
    while(counter4 <= nrow(temp_app))
    {
      if (counter4 != 1 && ps[(counter4 - 1)] == 1)
      {
        b4 <- (temp_app$Length[counter4] - 2)%%3
        ps <- append(ps, (b4))
        counter4 <- counter4 + 1
      }
      else if (counter4 != 1 && ps[(counter4 - 1)] == 2)
      {
        b4 <- (temp_app$Length[counter4] - 1)%%3
        ps <- append(ps, (b4))
        counter4 <- counter4 + 1
      }
      else
      {
        ps <- append(ps, (temp_app$Length[counter4])%%3)
        counter4 <- counter4 + 1
      }
    }
    
    temp_app$PhaseShift <- ps
    
    species1_min = min(isolate(temp_app["Start"]))
    
    species1_max = max(isolate(temp_app["End"]))
    
    largestLength <- 3000
    
    emf(paste("Files/", paste(isolate(as.character(input$app_gene_name)), " Exon Map", ".emf", sep = "")), width=2560, height=180, family = "DejaVu Serif", emfPlus = FALSE, emfPlusFont=FALSE, emfPlusRaster =TRUE)
    par(mar=c(0, 0, 0, 0))
    
    plot.new()
    
    height1 = 0.5
    
    # Draw line under gene boxes
    end_of_intron <- (isolate(temp_app[nrow(temp_app),"End"])-species1_min)/largestLength
    polygon(x = c((end_of_intron), (end_of_intron), 0, 0) , y = c((height1 - 0.045), (height1 + 0.045), (height1 + 0.045), (height1 - 0.045)), col = "#333333")
    
    # Adding Gene Segments
    boxHeight = (as.numeric(input$app_exon_height))/700;
    
    for (i in 1:dim(temp_app)[1])
    {
      # Define the "left" and "right" x-coordinates of boxes
      left = (isolate(temp_app[i,"Start"])-species1_min)/largestLength
      right = (isolate(temp_app[i,"End"])-species1_min)/largestLength
      height = height1
      
      # Draw boxes
      boxStart = max (left, right)
      polygon (c(left, boxStart, right, boxStart, left), c(height-boxHeight/2, height-boxHeight/2, height, height+boxHeight/2, height+boxHeight/2), col = e, lwd = 0.75)
      
      # If the check box to show exon length is checked, display length above exon
      if(input$app_show_length && input$app_above_exon)
      {
        if(input$app_exon_height >= 40)
        {
          text((left + ((isolate(temp_app[i, "Length"]))/(2*largestLength))), (0.58 + boxHeight/2), font = 2, as.character(temp_app$Length[i]), col = "black", cex = 115*input$app_text_size)
        }
        else
        {
          text((left + ((isolate(temp_app[i, "Length"]))/(2*largestLength))), (0.58 + 20/700), font = 2, as.character(temp_app$Length[i]), col = "black", cex = 115*input$app_text_size)
        }
      }
      else if(input$app_show_length && input$app_above_exon == FALSE)
      {
        if(input$app_exon_height >= 40)
        {
          text((left + ((isolate(temp_app[i, "Length"]))/(2*largestLength))), (0.46 + boxHeight/2), font = 2, as.character(temp_app$Length[i]), col = "black", cex = 75*input$app_text_size)
        }
        else
        {
          text((left + ((isolate(temp_app[i, "Length"]))/(2*largestLength))), 0.4957, font = 2, as.character(temp_app$Length[i]), col = "black", cex = 75*input$app_text_size)
        }
      }
    }
    
    offset_value <- as.numeric(ifelse(as.numeric(input$app_intron_length) == 0, -44.75, -((44.75/as.numeric(input$app_intron_length)) + 44.75)))
    
    for (i in 1:dim(temp_app)[1])
    {
      # Define the "left" and "right" x-coordinates of boxes
      left = (isolate(temp_app[i,"Start"])-species1_min)/largestLength
      right = (isolate(temp_app[i,"End"])-species1_min)/largestLength
      height = height1
      height2 <- 0.49-(input$app_text_size/250)
      
      if(i < isolate(nrow(temp_app)))
      {
        if(isolate(temp_app$PhaseShift[i] == 0))
        {
          draw.circle((right + ((intron_length/2)/largestLength)), height1 , 0.005*((input$app_phase_shift_size)/100), col = b)
          text((right + ((intron_length/2)/largestLength)), height2, font = 2, as.character(temp_app$PhaseShift[i]), col = in_ps0, cex = 125*input$app_text_size, pos = 2, offset = offset_value, adj = c(0.5, 0.5))
        }
        else if(isolate(temp_app$PhaseShift[i] == 1))
        {
          draw.circle((right + ((intron_length/2)/largestLength)), height1 , 0.005*((input$app_phase_shift_size)/100), col = b)
          text((right + ((intron_length/2)/largestLength)), height2, font = 2, as.character(temp_app$PhaseShift[i]), col = in_ps1, cex = 125*input$app_text_size, pos = 2, offset = offset_value, adj = c(0.5, 0.5))
        }
        else
        {
          draw.circle((right + ((intron_length/2)/largestLength)), height1 , 0.005*((input$app_phase_shift_size)/100), col = b)
          text((right + ((intron_length/2)/largestLength)), height2, font = 2, as.character(temp_app$PhaseShift[i]), col = in_ps2, cex = 125*input$app_text_size, pos = 2, offset = offset_value, adj = c(0.5, 0.5))
        }
      }
    }
    dev.off() 
    
  }
  
  # Render a plot using the user inputs as a preview
  output$plot <- renderPlot({
    
    if(nchar(input$app_sequence) != 0){
    
      intron_length <- as.numeric(input$app_intron_length)
    
    upper_app <- gsub("\\s+", "", as.character(input$app_sequence))
    upper_app <- gsub("a|c|g|t", " ", upper_app)
    upper_app <- gsub("^\\s+|\\s+$", "", upper_app)
    upper_app <- gsub("\\s+", " ", upper_app)
    
    str_replace(upper_app, "\\s{2,}", "")
    
    df1_app <- data.frame(upper_app)
    temp_app <- separate_longer_delim(df1_app, cols = 1, delim = ' ')
    temp_app <- data.frame(NA, NA, str_count(temp_app$upper_app))
    colnames(temp_app) <- c("Start", "End", "Length")
    
    values <- input$app_sequence
    
    # Uses Length values of previous row to get Start values
    
    counter <- 2
    newdata <- c((intron_length))
    while(counter <= nrow(temp_app))
    {
      newdata <- append(newdata, (temp_app$Length[(counter - 1)] + newdata[(counter - 1)] + (intron_length)))
      counter <- counter + 1
    }
    
    temp_app$Start <- newdata
    
    # Use Start values to get End values of previous row
    counter2 <- 1
    newdata2 <- c()
    while(counter2 <= nrow(temp_app)){
      newdata2 <- append(newdata2, (temp_app$Start[(counter2)] + temp_app$Length[(counter2)]))
      counter2 <- counter2 + 1
    }
    
    temp_app$End <- newdata2
    
    # Add phase shift column
    ps <- c()
    counter4 <- 1
    
    while(counter4 <= nrow(temp_app))
    {
      if (counter4 != 1 && ps[(counter4 - 1)] == 1)
      {
        b4 <- (temp_app$Length[counter4] - 2)%%3
        ps <- append(ps, (b4))
        counter4 <- counter4 + 1
      }
      else if (counter4 != 1 && ps[(counter4 - 1)] == 2)
      {
        b4 <- (temp_app$Length[counter4] - 1)%%3
        ps <- append(ps, (b4))
        counter4 <- counter4 + 1
      }
      else
      {
        ps <- append(ps, (temp_app$Length[counter4])%%3)
        counter4 <- counter4 + 1
      }
    }
    
    temp_app$PhaseShift <- ps
    
    species1_min = min(isolate(temp_app["Start"]))
    
    species1_max = max(isolate(temp_app["End"]))
    
    largestLength <- 3000
    
    
    par(mar=c(0, 0, 0, 0), bg = "transparent")
    plot.new()
    
    height1 = 0.5
    
    # Draw line under gene boxes
    end_of_intron <- ((isolate(temp_app[nrow(temp_app),"End"])-species1_min)/largestLength)*1.25
    polygon(x = c((end_of_intron), (end_of_intron), 0, 0) , y = c((height1 - 0.005), (height1 + 0.005), (height1 + 0.005), (height1 - 0.005)), col = "#333333")
    
    # Adding Gene Segments
    boxHeight = (as.numeric(input$app_exon_height))/2200;
    
    for (i in 1:dim(temp_app)[1])
    {
      # Define the "left" and "right" x-coordinates of boWxes
      left = (isolate(temp_app[i,"Start"])-species1_min)/largestLength
      right = (isolate(temp_app[i,"End"])-species1_min)/largestLength
      height = height1
      
      # Draw boxes
      boxStart = max (left, right)
      polygon (c(left*1.25, boxStart*1.25, right*1.25, boxStart*1.25, left*1.25), c(height-boxHeight/2, height-boxHeight/2, height, height+boxHeight/2, height+boxHeight/2), col = input$exon_color, lwd = 0.75)
      
      if(input$app_show_length && input$app_above_exon)
      {
        if(input$app_exon_height >= 20)
        {
          text((left*1.25 + ((isolate(temp_app[i, "Length"]))/(2*largestLength))*1.25), (0.515 + boxHeight/2), font = 2, as.character(temp_app$Length[i]), col = "black", cex = 0.7*input$app_text_size)
        }
        else
        {
          text((left*1.25 + ((isolate(temp_app[i, "Length"]))/(2*largestLength))*1.25), (0.515 + 10/2200), font = 2, as.character(temp_app$Length[i]), col = "black", cex = 0.7*input$app_text_size)
        }
      }
      else if(input$app_show_length && input$app_above_exon == FALSE)
      {
        if(input$app_exon_height >= 20)
        {
          text((left*1.25 + ((isolate(temp_app[i, "Length"]))/(2*largestLength))*1.25), (0.4885 + boxHeight/2), font = 2, as.character(temp_app$Length[i]), col = "black", cex = 0.7*input$app_text_size)
        }
        else
        {
          text((left*1.25 + ((isolate(temp_app[i, "Length"]))/(2*largestLength))*1.25), 0.5, font = 2, as.character(temp_app$Length[i]), col = "black", cex = 0.7*input$app_text_size)
        }
      }
    }
    for (i in 1:dim(temp_app)[1])
    {
      # Define the "left" and "right" x-coordinates of boxes
      left = (isolate(temp_app[i,"Start"])-species1_min)/largestLength
      right = ((isolate(temp_app[i,"End"])-species1_min)/largestLength)
      height = height1
      height2 <- 0.504-(input$app_text_size/250)
      
      if(i < isolate(nrow(temp_app)))
      {
        if(isolate(temp_app$PhaseShift[i] == 0))
        {
          draw.circle((right*1.25 + ((intron_length/2)/largestLength)*1.25), height1 , 0.0125*((input$app_phase_shift_size)/100), col = input$ps_bkgd_color)
          text((right*1.25 + ((intron_length/2)/largestLength)*1.25), height2, font = 2, as.character(temp_app$PhaseShift[i]), col = input$ps0_color, cex = 0.75*input$app_text_size, pos = 3, offset = -0.22, adj = c(0.5, NA))
        }
        else if(isolate(temp_app$PhaseShift[i] == 1))
        {
          draw.circle((right*1.25 + ((intron_length/2)/largestLength)*1.25), height1 , 0.0125*((input$app_phase_shift_size)/100), col = input$ps_bkgd_color)
          text((right*1.25 + ((intron_length/2)/largestLength)*1.25), height2, font = 2, as.character(temp_app$PhaseShift[i]), col = input$ps1_color, cex = 0.75*input$app_text_size, pos = 3, offset = -0.22, adj = c(0.5, NA))
        }
        else
        {
          draw.circle((right*1.25 + ((intron_length/2)/largestLength)*1.25), height1 , 0.0125*((input$app_phase_shift_size)/100), col = input$ps_bkgd_color)
          text((right*1.25 + ((intron_length/2)/largestLength)*1.25), height2, font = 2, as.character(temp_app$PhaseShift[i]), col = input$ps2_color, cex = 0.75*input$app_text_size, pos = 3, offset = -0.22, adj = c(0.5, NA))
        }
      }
    }
  }
    
  })
  
  output$txt <- renderText({
    
   txtout$txt
    
  })
  
  # Download the currently displayed plot as an .emf file when the "Download" button is pressed  
  output$downloadData <- downloadHandler(
    filename = function(){paste("Exon Maps ", date(), ".zip", sep = "")},
    content = function(file){
      zip(zipfile = "data", files = "Files/")
      file.copy("data.zip", file)
    }
  )
}                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      