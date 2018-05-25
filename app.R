library(shiny)
library(shiny.semantic)
library(DT)
library(dplyr)
library(plotly)
library(webshot)
library(ggplot2)
library(reshape2)
library(htmlwidgets)

###############
#App Functions#
###############
dataOutput <- function(database, kinSelection){
  #Load correct Database and find indices of chosen kinases
  compoundTable = switch(database,
                         rbio = {read.csv('rbio_old_database.csv', header = TRUE, stringsAsFactors = FALSE)},
                         LINCS = {read.csv('LINCS_database.csv', header = TRUE, stringsAsFactors = FALSE)},
                         PKIS = {read.csv('PKIS_database.csv', header = TRUE, stringsAsFactors = FALSE)},
                         EMD = {read.csv('EMD_database.csv', header = TRUE, stringsAsFactors = FALSE)},
                         {read.csv('rbio_old_database.csv', header = TRUE, stringsAsFactors = FALSE)}
  )
  kinaseNames <- variable.names(compoundTable)[-c(1:3)]
  kinIndex <- which(kinaseNames %in% kinSelection)
  
  #Remove compounds not profiled for chosen set, begin making Table of Results with inhibition of chosen kinase(s)
  candidateTable <- compoundTable[!c(apply(as.matrix(compoundTable[,kinIndex+3]), 1, function(x) any(is.na(x)))),]
  resTable <- candidateTable[,c(1:3,kinIndex+3)]
  
  #Create scaling factors as maximum of each row (drug), table subsets with on/off/all target inhibitions only, and totals
  compoundScalings <- apply(candidateTable[,-(1:3)],1,max, na.rm = TRUE)
  totalTableScaled <- candidateTable[,-(1:3)] * 100 / compoundScalings
  onTargetTable <- as.matrix(resTable[,4:(length(kinSelection)+3)])
  onTargetTableScaled <- onTargetTable * 100 / compoundScalings
  offTargetTable <- as.matrix(candidateTable[,-c(1:3, kinIndex+3)])
  offTargetTableScaled <- offTargetTable * 100 / compoundScalings
  
  #Calculate On Target, Off Target, and Selectivity parameters and Scores
  resTable$OnGeomMean <- apply(onTargetTable, 1, function(x) exp(mean(log(x), na.rm = TRUE)))
  
  resTable$OnGeomMeanScaled <- apply(onTargetTableScaled, 1, function(x) exp(mean(log(x), na.rm = TRUE)))
  
  resTable$OffMean <- apply(offTargetTableScaled, 1, function(x) mean(x, na.rm = TRUE))
  resTable$OffRMS <- sqrt(apply(offTargetTableScaled, 1, function(x) mean(x^2, na.rm = TRUE)))
  
  resTable$OffVar <- apply(offTargetTableScaled, 1, function(x) var(x, na.rm = TRUE))
  resTable$AllVar <- apply(totalTableScaled, 1, function(x) var(x, na.rm = TRUE))
  #resTable$VarLoBound <- 5000 / rowSums(!is.na(totalTableScaled))
  resTable$VarLoBound <- apply(offTargetTableScaled, 1, max, na.rm = TRUE)^2 / (2*rowSums(!is.na(totalTableScaled)))
  
  resTable$VarRatioNorm <- (resTable$OffVar - resTable$VarLoBound) / (resTable$AllVar - resTable$VarLoBound)
  
  #resTable$Penalty <- 0.5*(resTable$OffMean + resTable$OffRMS)
  #resTable$Penalty <- 0.5*(resTable$OffMean + 100*resTable$VarRatioNorm)
  resTable$Penalty <- 0.5*(resTable$OffMean + 100*(resTable$VarRatioNorm^5))
  
  resTable$SS <- resTable$OnGeomMeanScaled - resTable$Penalty
  
  #Extract other helpful information about off target effects
  resTable$NumMoreInh <- rowSums(offTargetTable >= resTable$OnGeomMean, na.rm = TRUE)
  resTable$NumMoreHalf <- rowSums(offTargetTable >= 0.5*resTable$OnGeomMean, na.rm = TRUE)
  resTable$NumNAs <- rowSums(is.na(offTargetTable))
  
  row.names(resTable) <- NULL
  
  tf <- sweep(offTargetTable, 1, 0.5*resTable$OnGeomMean, ">=")
  res <- data.frame(
    row   = as.vector(row(offTargetTable)),
    Kinase  = colnames(offTargetTable)[col(offTargetTable)],
    Value = as.vector(offTargetTable),
    pass  = as.vector(tf)
  )
  res <- res[order(-res$Value), ]
  offTargets <- lapply(split(res, res$row), function(x) {
    x <- x[complete.cases(x), ]
    colnames(x)[3] <- '% Inhibition'
    x[x$pass, c("Kinase", "% Inhibition")]
  })
  
  #Format output table
  outputTable <- resTable
  outputTable <- outputTable[complete.cases(outputTable),]
  outputTable <- subset(outputTable, select = -c(OnGeomMean, OnGeomMeanScaled, OffMean, OffRMS, OffVar, AllVar, VarLoBound, VarRatioNorm, Penalty))
  
  colnames(outputTable) <- c(colnames(outputTable[1:3]),
                             paste0(colnames(outputTable)[4:(length(kinSelection)+3)], ' % inhibition'),
                             'KInhibition Selectivity Score', '# of Off Targets > Selection',
                             '# of Off Targets > Half Selection', '# of Uncharacterized Kinases')
  
  #Make full inhibition profiles for candidate drugs
  resultsProfile <- compoundTable[match(paste0(outputTable[,1],outputTable[,3]), paste0(compoundTable[,1],compoundTable[,3])),]
  
  #Output: table, database, kinase selection, compound profiles
  return(list(outputTable, compoundTable, kinSelection, resultsProfile, offTargets))
} 

hmOutput <- function(compoundTable, kinSelection, resultsProfile, inputRows){
  #Plotly Heatmap
  resultsProfile <- resultsProfile[inputRows,]
  hmData <- as.matrix(resultsProfile[,-c(1:3)])
  rownames(hmData) <- paste0(resultsProfile[,1], ' (', resultsProfile[,3], ')')
  hmMelt <- melt(hmData)
  colnames(hmMelt) <- c('Compound', 'Kinase', '% Inhibition')
  hmAnno <- subset(hmMelt, (Kinase %in% kinSelection & Compound == rownames(hmData)[1]), select = c(Compound, Kinase))
  valrange <- unique(sort(scales::rescale(c(as.matrix(compoundTable[,-c(1:3)])))))
  colorscale <- data.frame(seq(from=0, to = 1, by = 0.001), colorRampPalette(c('black', 'yellow', 'red'), bias = 1.00)(1001))
  p <- plot_ly(data = hmMelt, x = ~Kinase, y = ~Compound, z = ~`% Inhibition`, type = "heatmap",
               colorscale = colorscale, zauto = FALSE, zmin = 0, zmax = 100,
               text = ~paste0("Kinase: ", Kinase, '\nCompound: ', Compound, '\n% Inhibition: ', `% Inhibition`), 
               hoverinfo = "text", hoverlabel = list(bgcolor = toRGB('gray30'), bordercolor = toRGB('floralwhite'), font = list(color = toRGB('floralwhite'))),
               connectgaps = FALSE, width = 1300, height = 26*length(inputRows)+110)%>%
    config(displayModeBar = FALSE)%>% 
    layout(margin = list(b = 25, l = 8.3*max(nchar(rownames(hmData)))+10, t = 35),
           yaxis = list(autorange = "reversed", gridwidth = 0, gridcolor = "#FFFFFF", ticks = '', title = '', showline = FALSE,
                        tickfont = list(color = "#000000", size = 14)), 
           xaxis = list(gridwidth = 0, gridcolor = "#FFFFFF", ticks = '', title = 'All Kinases', titlefont = list(size = 22, color = 'black'),
                        showline = FALSE, showticklabels = FALSE))%>%
    add_annotations(x = hmAnno$Kinase, y = -0.5, text = hmAnno$Kinase, 
                    xref = "x", yref = "y",
                    showarrow = TRUE, arrowhead = 4, arrowsize = 0.5, arrowcolor = "#000000",
                    ax = 0, ay = -22)%>%
    colorbar(thickness = 30, thicknessmode = 'pixels', len = 250, lenmode = 'pixels',  ypad = 0, xpad = 0,
             y = 0.5,
             outlinecolor = "#000000", bordercolor = "#000000", outlinewidth = 2, borderwidth = 0,
             tickvals = c(0,10,20,30,40,50,60,70,80,90,100), ticklabels = c('0','10','20','30','40','50','60','70','80','90','100'),
             ticks = 'outside', tickcolor = '#000000', ticklen = 10, tickwidth = 2,
             title = '% Inhibition', titlefont = list(size = 16, color = '#000000'))
  #Output: heatmap
  return(p)
}
########
#App UI#
########
ui <- fluidPage(
  titlePanel(div(img(src='FredHutchLogo.png', align = "right", height = 75, width = 299),
             h1(strong('KInhibition: Kinase Inhibitor Selection Tool'), align = 'center')), windowTitle = 'KInhibition'),
  tags$head(tags$link(rel="icon", type="image/x-icon", href="data:image/x-icon;base64,AAABAAEAEBAAAAEAIABoBAAAFgAAACgAAAAQAAAAIAAAAAEAIAAAAAAAAAQAABILAAASCwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADwAAAH0AAABzAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAUQAAAOwAAAD/AAAAhAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAIQAAAD/AAAA/wAAAIQAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACEAAAA/wAAAP8AAACEAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAhAAAAP8AAAD/AAAAhAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAIQAAAD/AAAA/wAAAIQAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACEAAAA/wAAAP8AAACEAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAkwAAAP8AAAD/AAAAkwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAOgAAAPkAAAD/AAAA/wAAAPkAAAA6AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAFwAAAOMAAAD/AAAA/wAAAP8AAAD/AAAA4wAAABcAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAABAAAAL0AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAAC9AAAABAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAIgAAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAIgAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAFEAAAD9AAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD9AAAAUQAAAAAAAAAAAAAAAAAAACUAAADvAAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAO8AAAAlAAAAAAAAAAAAAADAAAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAAwAAAAAAAAAAAAAAA1wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAANcAAAAA/j8AAPw/AAD8PwAA/D8AAPw/AAD8PwAA/D8AAPw/AAD4HwAA8A8AAOAHAADAAwAAwAMAAIABAACAAQAAgAEAAA==")),
  sidebarLayout(
    sidebarPanel(
      selectizeInput("database", label = h3(strong("Database")),
                     choices = list("Reaction Biology" = 'rbio',
                                    "HMS LINCS" = 'LINCS',
                                    "GSK PKIS" = 'PKIS',
                                    "EMD Millipore" = "EMD")),
      uiOutput("kinaseSelection"),
      hr(), h4(strong('Instructions:')), p(strong('1. Select a Database'), 'to search (see \"Database Descriptions\" tab on the right).'), 
      p(strong('2. Select the kinase(s)'), 'to inhibit. Note that selecting more than one compound may result in poor selectivity for most compounds.'),
      p(strong('3. View the Table of Results'), 'tab to view the properties of the compounds. See below the table for descriptions of the column headers.'),
      p(strong('4. Customize the results'), 'shown in the table by using the column headers to sort the compounds by various properties, 
        the dropdown menu in the upper left to show more or fewer compounds, or the \"Column Visibility\" button to show or hide the columns.'),
      p(strong('5. Check off-target effects'), 'of the compounds in the table by clicking on the desired compound. Non-chosen kinases with inhibition higher 
        than half of the inhibition of the chosen kinase(s) will be shown below.'),
      p(strong('6. View the Heatmap'), 'tab to visualize the inhibition profiles of the compounds present in the \"Table of Results\".
        Mouse over cells in the heatmap for more information.'),
      p(strong('7. Download'), 'the tables or heatmaps using the buttons below them.'),
      hr(), p('Citation / Publication is'), p(a('coming soon...', href = 'http://www.pubmed.com')),
      hr(), p(em('Disclaimer: This application should not be considered, or used as a substitute for, medical advice, diagnosis or treatment. 
                 This site does not constitute the practice of any medical, nursing or other professional health care advice, diagnosis or treatment.'),
              style = 'color:grey')),
    
    mainPanel( tabsetPanel(
      tabPanel(strong("Table of Results"), uiOutput("tableUI"), uiOutput("offTargetUI"), uiOutput("descriptionUI"), hr()),
      tabPanel(strong("Heatmap"), uiOutput('heatmapUI'), hr(),
               span(em('Please note that heatmap may take a few seconds to fully load, update, or download.'), style = 'color:grey'), br(),
               span(em('To change or update the data being displayed, use the Table of Results tab.'), style = 'color:grey')),
      tabPanel(strong("Database Descriptions"), 
               h2('Reaction Biology'),
               p('This screen was performed through a collaboration between the Fox Chase Cancer Center and Reaction Biology Corporation. Details about the assay can be found in the original publication.'), 
               p('The data can also be explored', a("here.", href = 'http://reactionbiology.com/webapps/largedata/')),
               p(strong('Citation:'), 'Anastassiadis, T., Deacon, S. W., Devarajan, K., Ma, H., & Peterson, J. R. (2011). Comprehensive assay of kinase catalytic activity reveals features of kinase inhibitor selectivity. Nature Biotechnology, 29(11), 1039-1045. https://doi.org/10.1038/nbt.2017'),
               p(strong('Number of Compounds: '), '178'),
               p(strong('Number of compound/dose combinations: '), '178'),
               p(strong('Number of kinases: '), '300'),
               p(strong('Pairwise Coverage: '), '98.9%'),
               br(),
               h2('HMS LINCS'),
               p('The Library of Integrated Network-based Cellular Signatures (LINCS) is an NIH-funded program at Harvard Medical School. This dataset uses the', a('KINOMEscan Assay', href = 'https://www.discoverx.com/services/drug-discovery-development-services/kinase-profiling/kinomescan'), 'to profile kinases bound by inhibitors.'),
               p('The resulting datasets can be accessed ', a('here.', href = 'http://lincs.hms.harvard.edu/kinomescan/')),
               p(strong('Citation: '), 'Koleti, A., Terryn, R., Stathias, V., Chung, C., Cooper, D. J., Turner, J. P., . Schurer, S. C. (2017). Data Portal for the Library of Integrated Network-based Cellular Signatures (LINCS) program: integrated access to diverse large-scale cellular perturbation response data. Nucleic Acids Research. https://doi.org/10.1093/nar/gkx1063'),
               p(strong('Number of Compounds:'), '104'),
               p(strong('Number of compound/dose combinations:'), '112'),
               p(strong('Number of kinases: '), '494'),
               p(strong('Pairwise Coverage: '), '84.5%'),
               br(),
               h2('GSK PKIS'),
               p('A library of the GlaxoSmithKline (GSK) Published Kinase Inhibitor Set (PKIS) was assayed at the NIH using luciferase reporter enzymes. Details about the assay can be found in below publication.'),
               p('The published datasets can be accessed  by following the instructions found ', a('here.', href = 'https://www.ebi.ac.uk/chembldb/extra/PKIS/')),
               p(strong('Citation: '), 'Dranchak, P., MacArthur, R., Guha, R., Zuercher, W. J., Drewry, D. H., Auld, D. S., & Inglese, J. (2013). Profile of the GSK Published Protein Kinase Inhibitor Set Across ATP-Dependent and-Independent Luciferases: Implications for Reporter-Gene Assays. PLoS ONE, 8(3), e57888. https://doi.org/10.1371/journal.pone.0057888'),
               p(strong('Number of Compounds: '), '367'),
               p(strong('Number of compound/dose combinations: '), '734'),
               p(strong('Number of kinases: '), '224'),
               p(strong('Pairwise Coverage: '), '99.9%'),
               br(),
               h2('EMD Millipore'),
               p('Pairwise assays between purified recombinant human kinases and a collection of small molecule inhibitors were carried out by the EMD Millipore Corporation (now known as Millipore Sigma) using a filter binding radioactive ATP transferase assay.'),
               p('Raw data and details on the recombinant kinases and inhibitors used can be found in the supplemental files of the original publication ', a('here.', href = 'http://www.biochemj.org/content/451/2/313')),
               p(strong('Citation: '), 'Gao, Y., Davies, S. P., Augustin, M., Woodward, A., Patel, U. A., Kovelman, R., & Harvey, K. J. (2013). A broad activity screen in support of a chemogenomic map for kinase signalling research and drug discovery. Biochemical Journal, 451(2), 313-328. https://doi.org/10.1042/BJ20121418'),
               p(strong('Number of Compounds: '), '128'),
               p(strong('Number of compound/dose combinations: '), '255'),
               p(strong('Number of kinases: '), '234'),
               p(strong('Pairwise Coverage: '), '100%'),
               br(),
               hr(),
               p('If there is a database or compound screen that you would like to see included in this tool, please contact the authors of the publication.')      )
    ))
  )
)

############
#App Server#
############
server <- function(input, output, session) {
  
  output$kinaseSelection <- renderUI({
    validate(need(input$database, "Select a Database to search"))
    compoundTable = switch(input$database,
                           rbio = {read.csv('rbio_old_database.csv', header = TRUE, stringsAsFactors = FALSE)},
                           LINCS = {read.csv('LINCS_database.csv', header = TRUE, stringsAsFactors = FALSE)},
                           PKIS = {read.csv('PKIS_database.csv', header = TRUE, stringsAsFactors = FALSE)},
                           EMD = {read.csv('EMD_database.csv', header = TRUE, stringsAsFactors = FALSE)},
                           {read.csv('rbio_old_database.csv', header = TRUE, stringsAsFactors = FALSE)}
    )
    kinaseNames <- variable.names(compoundTable)[-c(1:3)]
    selectizeInput("kinSelection", label = h3(strong("Kinase(s) of Interest")),
                   choices = kinaseNames,
                   multiple = TRUE,
                   options = list(
                     placeholder = 'Select kinase(s) of interest...',
                     onInitialize = I('function() { this.setValue(""); }')))
  })
  
  output$nodrugs <- renderText("There are no compounds in the selected dataset characterized for this set of kinases.")
  output$nodrugs2 <- renderText("There are no compounds in the selected dataset characterized for this set of kinases.")
  
  observeEvent(input$database, updateTextInput(session, "kinSelection", value = ""), priority = 10)
  
  funcOutput <- eventReactive(c(input$kinSelection), {
    req(input$kinSelection)
    try(dataOutput(database = input$database, kinSelection = input$kinSelection), silent = FALSE)}, ignoreNULL = FALSE)
  
  heatmapPlot <- eventReactive(c(funcOutput(), input$table_rows_current), {
    req(input$kinSelection)
    hmOutput(isolate(funcOutput()[[2]]), isolate(funcOutput()[[3]]), isolate(funcOutput()[[4]]), isolate(input$table_rows_current))}, ignoreNULL = FALSE)
    
  output$tableUI <- renderUI({
    validate(need(input$kinSelection, message = "\nSelect kinase(s) of interest\n"))
    validate(need(funcOutput(), message = "\nThere are no compounds in the selected dataset characterized for this set of kinases.\n"))
    numKins <- length(isolate(input$kinSelection))
    output$table <- renderDT({
      req(funcOutput())
      DT::datatable(funcOutput()[[1]], filter = 'none', fillContainer = FALSE, selection = 'single', rownames = FALSE, 
                    options = list(dom = 'lfBrtip', autoWidth = TRUE, scrollX = TRUE,
                                   lengthMenu = list(c(10, 20, 50, 100, -1), c(c('10', '20', '50', '100', 'All'))),
                                   order = list(numKins+3, 'desc'),
                                   columnDefs = list(list(width = '200px', targets = 0),
                                                     list(width = '100px', targets = 1),
                                                     list(width = '10px', targets = 2),
                                                     list(width = '30px', targets = c(3:(numKins+2)))))) %>% 
        formatRound(columns = c((4:(numKins+3)), (4+numKins)), digits = 2) %>%
        formatStyle('KInhibition Selectivity Score', fontWeight = 'bold')
    })
    outputOptions(output, "table", suspendWhenHidden = FALSE)
    list(DTOutput("table"), downloadButton("downloadTable", 'Download Table of Results (.csv)'))
  })
  
  output$offTargetUI <- renderUI({
    req(funcOutput())
    req(input$table_rows_selected)
    output$offTarget <- renderTable({
      req(funcOutput())
      req(input$table_rows_selected)
      funcOutput()[[5]][[input$table_rows_selected]]
    })
    list(hr(),
         h3(strong('Off Target Effects of ', funcOutput()[[1]][input$table_rows_selected,1])),
         tableOutput("offTarget"),
         downloadButton("downloadOffTarget", 'Download Off Target Effects (.csv)'))
  })
  
  output$descriptionUI <- renderUI({
    validate(need(funcOutput(), ""))
    list(hr(), p(strong('Compound: '), 'The molecule or drug used in the study'), 
         p(strong('Identifier: '), 'The second column contains a unique identifier (CAS number, Library ID, etc.) for the given compound used in the study'), 
         p(strong('Dose: '), 'The concentration of the compound at which the given inhibition profile was observed'),
         p(strong('% Inhibition: '), 'The inhibition of the named kinase by the given compound, as percent compared to control.
                                       100 means total inhibition (no kinase activity detected), 0 means no inhibition (kinase activity is equal to or greater than control)'),
         p(strong('KInhibition Selectivity Score: '), 'a selectivity summary score, generally between -100 and 100, where 100 represents perfect on target inhibition, and -100 represents perfect off-target inhibition.
           For more details, refer to the manuscript (link in sidebar).'),
         p(strong('# of Off Targets > Selection: '), 'The number of non-chosen kinases with inhibition higher than the inhibition of the chosen kinase. For multiple chosen kinases, a geometric mean of all on-target inhibitions is used for comparisons.'),
         p(strong('# of Off Targets > Half Selection: '), 'The number of non-chosen kinases with inhibition higher than half of the inhibition of the chosen kinase. For multiple chosen kinases, a geometric mean of all on-target inhibitions is used for comparisons.'),
         p(strong('# of Uncharacterized Kinases: '), 'The number of non-chosen kinases for which inhibition data for this compound was not reported in the study'))
  })
  
  output$heatmapUI <- renderUI({
    validate(need(input$kinSelection, message = "\nSelect kinase(s) of interest\n"))
    validate(need(funcOutput(), message = "\nThere are no compounds in the selected dataset characterized for this set of kinases.\n"))
    req(input$table_rows_current)
    output$heatmap <- renderPlotly({heatmapPlot()})
    outputOptions(output, "heatmap", suspendWhenHidden = TRUE)
    list(plotlyOutput("heatmap", width = "1200", height = "100%"), 
         downloadButton("downloadHeatmap", "Download Static Heatmap (.png)"),
         downloadButton("downloadHeatmapInteractive", "Download Interactive Heatmap (.html)"))
  })
  
  output$downloadTable <- downloadHandler(
    filename = function() {paste0("KInhibition_", paste0(isolate(input$kinSelection), collapse = '+'), "_inhibitors.csv")},
    content = function(file) write.csv(funcOutput()[[1]][input$table_rows_all,], file = file, row.names = FALSE)
  )
  
  output$downloadOffTarget <- downloadHandler(
    filename = function() {paste("KInhibition",
                                paste0(isolate(input$kinSelection), collapse = '+'),
                                funcOutput()[[1]][input$table_rows_selected,1], 
                                "offTargets.csv", sep = '_')},
    content = function(file) {write.csv(funcOutput()[[5]][[input$table_rows_selected]], file = file, row.names = FALSE)}
  )
  
  output$downloadHeatmap <- downloadHandler(
    filename = function() {paste0("KInhibition_", paste0(input$kinSelection, collapse = '+'), "_heatmap.png")},
    content = function(file) {export(heatmapPlot(), file = file)}
  )
  
  output$downloadHeatmapInteractive <- downloadHandler(
    filename = function() {paste0("KInhibition_", paste0(isolate(input$kinSelection), collapse = '+'), "_heatmapInteractive.html")},
    content = function(file) {saveWidget(heatmapPlot(), file = file)}
  )
}

#####################
#Run The Application#
#####################
shinyApp(ui = ui, server = server)

