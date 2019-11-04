folderName = "/Users/Lab/FijiProgramming/output/p21 p6 1"
#Imports and Function Definitions#----
library(shiny)
library(shinyjs)
library(shinyBS)
knitr::opts_chunk$set(echo=FALSE,results='hide',fig.keep='all')
library(jsonlite)
savedParams=F
setwd(folderName)
if (savedParams){
  paramlist = as.list(unlist(read_json("params.json")))
  #This automatically assigns the read analysis parameters to R's global environment.
  list2env(paramlist, environment())
}
library(ggplot2)
th = theme_set(theme_bw() + theme(panel.grid.major = element_blank(), 
                                  panel.grid.minor = element_blank(),
                                  text = element_text(size=12),
                                  axis.title = element_text(face="bold")))
library(dplyr)
library(purrr)
library(autothresholdr)
library(rlang)
groupSet = character(0)
find_metadata = function(filename){
  strings = unlist(map(unlist(
    regmatches(filename, gregexpr("_[a-zA-Z1-9]+=[a-zA-Z0-9]+", filename))), 
    ~substring(., 2)))
  keys = unlist(map(unlist(
    regmatches(strings, gregexpr("^\\w+=", strings))), 
    ~substring(., 1, nchar(.)-1)))
  values = unlist(map(unlist(
    regmatches(strings, gregexpr("=\\w+$", strings))),
    ~substring(., 2)))
  setNames(values, keys)
}

bindOuterNames = function(folderName, pattern){
  setwd(folderName)
  retval = bind_rows(map(
    list.files(folderName, pattern), function(filename){
      df = read.csv(filename, stringsAsFactors=FALSE)
      pairs = find_metadata(filename)
      namelist = names(pairs)
      namebinder = function(y, z, x){
        i=1; while (i <= length(y)){
          if (!(y[i] %in% groupSet)){groupSet <<- append(groupSet, y[i])}
          x[y[i]] = as.character(z[i])
          i = i + 1}
        x}
      namebinder(namelist, pairs, df)
    }
  ))
}

findMaxDensity = function (intensity)
{  density(intensity)$x[which.max(density(intensity)$y)]}

#Data Imports#----
nuclei = bindOuterNames(folderName, "^nuclei_")
nucleilink = bindOuterNames(folderName, "^nucleilink_")
cm = bindOuterNames(folderName, "^cm_")
cyclenuclei = bindOuterNames(folderName, "^cyclenuclei_")
names(nuclei)[names(nuclei) == "X.1"] =  "NucleusID"
names(nucleilink)[names(nucleilink) == "X"] = "NucleusID"
names(cyclenuclei)[names(cyclenuclei) == "X"] = "NucleusID"
colnames(cyclenuclei)[colnames(cyclenuclei) == "Median"] = "cycleIntensity"
joins = character(0)
for (g in groupSet){joins[g] = g}
nuclei = inner_join(nuclei, cyclenuclei, by=c(c("NucleusID" = "NucleusID"), joins))
nuclei = inner_join(nuclei, nucleilink, by=c(c("NucleusID" = "NucleusID"), joins))

#Shiny Server#----
thisServer = shinyServer(function(input, output, session){
#Observed variables#----  
  observe({
    maxXVal = input$MaxX
    minXVal = input$MinX
    updateSliderInput(session, "MeanX", min = minXVal, max = maxXVal)
  })
  observe({
    maxXVal = input$AreaMaxX
    minXVal = input$AreaMinX
    updateSliderInput(session, "AreaMeanX", min = minXVal, max = maxXVal)
  })
  observe({
    maxXVal = input$CellCycleIntensityMaxX
    minXVal = input$CellCycleIntensityMinX
    updateSliderInput(session, "CellCycleIntensity", min = minXVal, max = maxXVal)
  })
  
#Plotting Code#----
  baseMean = ggplot(nuclei, aes(x=Mean)) + geom_density() + labs(x="Nuclear Mean Intensity", y="Density(AU)")
  baseArea = ggplot(nuclei, aes(x=Area)) + geom_density() + labs(x="Nuclear Area (pixels)", y="Density(AU)")
  baseFeret = ggplot(cm, aes(x=MinFeret)) + geom_density() + labs(x="CM Minimum Feret's Diameter  (pixels)", y="Density(AU)")
  output$selectMean = renderPlot(
    baseMean + geom_vline(xintercept=input$MeanX, color="blue")+ xlim(input$MinX, input$MaxX)
  )
  output$selectArea = renderPlot(
    baseArea + geom_vline(xintercept=input$AreaMeanX, color="red")+ xlim(input$AreaMinX, input$AreaMaxX)
  )
  output$selectFeret = renderPlot(
    baseFeret + geom_vline(xintercept=input$FeretMeanX, color="green")+ xlim(input$FeretMinX, input$FeretMaxX)
  )
  output$unnormalizedCellCycle = renderPlot(
    if(input$plotUnnormalized){
      
      ggplot(nuclei, aes(x=intensity, stat(count))) + geom_density() + labs(x="Integrated Intensity", y="Density (AU)")
    }
    else{
      ggplot()
    }
  )
  output$cycleIntensity = renderPlot(
    if(input$applyThresholds){
      
      ggplot(nuclei, aes(x=cycleIntensity, stat(count))) + geom_density() + labs(x="Cell Cycle Stain Intensity", y="Density (AU)")+ 
        geom_vline(xintercept=input$CellCycleIntensityMeanX, color="green") + 
        xlim(input$CellCycleIntensityMinX, input$CellCycleIntensityMaxX)
    }
    else{
      ggplot()
    }
  )
  output$normalizedCellCycle = renderPlot(
    if(input$plotNormalized){
      
      ggplot(nuclei, aes(x=intensity, stat(count))) + geom_density() + labs(x="Integrated Intensity", y="Density (AU)") +
        geom_vline(xintercept=input$Lower, color="red") +
        geom_vline(xintercept=input$Mid, color="green") +
        geom_vline(xintercept=input$Upper, color="blue")
    }
    else{
      ggplot()
    }
  )
  output$mononucleatedDistrubution = renderPlot(
    if (input$plotAndSave){
    ggplot(filter(cmnuclei, nucleiPerCm==1), aes(x=intensity)) + geom_histogram(bins=100, fill="black") + labs(y="Number of nuclei (mononucleated)", x="Estimated ploidy")} else {ggplot()}
  )
  output$binucleatedDistribution = renderPlot(
    if (input$plotAndSave){
    ggplot(filter(cmnuclei, nucleiPerCm==2), aes(x=intensity)) + geom_histogram(bins=100, fill="black") + labs(y="Number of nuclei (binucleated)", x="Estimated ploidy")} else {ggplot()}
  )
  output$ploidyDistributionByNucleation = renderPlot(
    if (input$plotAndSave){
    ggplot(filter(cmNucleusPortions, nucleiPerCm %in% c(1,2,3,4,5)), aes(x=as.factor(nucleiPerCm), y=portion)) + geom_col(fill="black") + labs(x="Nuclei per Cardiomyocyte", y="Cardiomyocyte population %") + ylim(0,100)} else {ggplot()}
  )
  # Insert the right number of plot output objects into the web page
  output$unnormPlots <- renderUI({
    plot_output_list <- lapply(groupSet, function(i) {
      plotname <- paste("plot", i, sep="")
      plotOutput(plotname)
    })
    
    # Convert the list to a tagList - this is necessary for the list of items
    # to display properly.
    do.call(tagList, plot_output_list)
  })
  
  output$normPlots <- renderUI({
    plot_output_list <- lapply(groupSet, function(i) {
      plotname <- paste("normPlot", i, sep="")
      plotOutput(plotname)
    })
    
    # Convert the list to a tagList - this is necessary for the list of items
    # to display properly.
    do.call(tagList, plot_output_list)
  })
#Event Observations#----
  observeEvent(input$applyCycleIntensityThresh,{
    nuclei <<- mutate(nuclei, cyclePositive = cycleIntensity > input$CellCycleIntensityMeanX)
  })
  
  observeEvent(input$applyThresholds, {
    cm <<- filter(cm, MinFeret < input$FeretMeanX)
    nuclei <<- filter(nuclei, Area > input$AreaMeanX, Mean > input$MeanX, Circ. > 0.4) %>% mutate(intensity = Mean * Area) %>%
    filter(intensity < 3 * quantile(intensity, 0.5))
  })
  observeEvent(input$normalize, {
    if (input$normalizeByWell) {
      nuclei <<- nuclei %>% group_by(!!!syms(groupSet)) %>% mutate(intensity = 2*intensity / findMaxDensity(intensity))
    } else {
      nuclei <<- nuclei %>% mutate(intensity = 2*intensity / findMaxDensity(intensity))    
    } 
  })
  observeEvent(input$ploidyNucleationCalculate, {
    nuclei <<- filter(nuclei, intensity > input$Lower & intensity < input$Upper) %>% mutate(FourN = intensity > input$Mid)
    cmnuclei <<- inner_join(filter(nuclei, Min>0), cm, by= c(c("Min" = "Mean"), joins))
    cmNucleusCounts <<- cmnuclei %>% group_by(!!!syms(c("Min", groupSet))) %>% summarize(nucleiPerCm = n())
    cmnuclei <<- inner_join(cmNucleusCounts, cmnuclei, by=c(c("Min"="Min"), joins))
    
    cmNucleusCyclePositivePortions <<- cmnuclei %>% group_by(nucleiPerCm) %>% summarize(portion = 100 * sum(cyclePositive == T) / n())
    ggplot(filter(cmNucleusCyclePositivePortions, nucleiPerCm %in% c(1,2,3,4,5)), aes(x=as.factor(nucleiPerCm), y=portion)) + geom_col(fill="black") + labs(x="Nuclei per Cardiomyocyte", y="Percentage Cycle Positive") + ylim(0,100)
    ggsave(paste(folderName,"cyclePositive.pdf"), width=3, height=3, device="pdf")
    cmnuclei = filter(cmnuclei, cyclePositive==T)
    cmNucleusPortions <<- cmNucleusCounts %>% group_by(nucleiPerCm) %>% summarize(portion = 100 * n()/nrow(cmNucleusCounts))
    cmNucleusFourPortions <<- cmnuclei %>% group_by(nucleiPerCm) %>% summarize(portion = 100 * sum(FourN == T) / n())
    write.csv(cmNucleusPortions, paste(folderName,"cmNucleusPortions.csv"))
    write.csv(cmNucleusFourPortions, paste(folderName,"cmNucleusFourPortions.csv"))
    write.csv(cmNucleusCyclePositivePortions, paste(folderName,"cmNucleusCyclePositivePortions"))
    
    th = theme_set(theme_bw() + theme(panel.grid.major = element_blank(), 
                                      panel.grid.minor = element_blank(),
                                      text = element_text(size=12),
                                      axis.title = element_text(face="bold")))
    
    ggplot(filter(cmnuclei, nucleiPerCm==1), aes(x=intensity)) + geom_histogram(bins=100, fill="black") + labs(y="Number of nuclei (mononucleated)", x="Estimated ploidy")
    ggsave(paste(folderName,"mononucPlot.pdf"), width=3, height=3, device="pdf")
    ggplot(filter(cmnuclei, nucleiPerCm==2), aes(x=intensity)) + geom_histogram(bins=100, fill="black") + labs(y="Number of nuclei (binucleated)", x="Estimated ploidy")
    ggsave(paste(folderName,"binucPlot.pdf"), width=3, height=3, device="pdf")
    ggplot(filter(cmnuclei, cyclePositive==T), aes(x=intensity)) + geom_histogram(bins=100, fill="black") + labs(y="Number of nuclei (cycle Positive)", x="Estimated ploidy")
    ggsave(paste(folderName,"cyclePositiveDNADistribution.pdf"), width=3, height=3, device="pdf")
    ggplot(filter(cmNucleusPortions, nucleiPerCm %in% c(1,2,3,4,5)), aes(x=as.factor(nucleiPerCm), y=portion)) + geom_col(fill="black") + labs(x="Nuclei per Cardiomyocyte", y="Cardiomyocyte population %") + ylim(0,100)
    ggsave(paste(folderName,"compnucPlot.pdf"), width=3, height=3, device="pdf")
  })
  
  # Call renderPlot for each one. Plots are only actually generated when they
  # are visible on the web page.
  for (g in groupSet) {
    # Need local so that each item gets its own number. Without it, the value
    # of i in the renderPlot() will be the same across all instances, because
    # of when the expression is evaluated.
    print(g)
    local({
      my_i <- g
      plotname <- paste("plot", my_i, sep="")
      plotname2 = paste("normPlot", my_i, sep="")
      
      output[[plotname]] <- renderPlot({
        if(input$plotUnnormalized){
        ggplot(nuclei, aes(x=intensity, color=as.factor(eval(parse(text=my_i))), stat(count))) + geom_density() + labs(x="Integrated Intensity",
                                                                                                                    y="Density (AU)",
                                                                                                                    color=my_i)} else {ggplot()}
      })
      output[[plotname2]] <- renderPlot({
        if(input$plotNormalized){
          ggplot(nuclei, aes(x=intensity, color=as.factor(eval(parse(text=my_i))), stat(count))) + geom_density() + labs(x="Ploidy",
                                                                                                                         y="Density (AU)",
                                                                                                                         color=my_i)} else {ggplot()}
      })
      
    })
  }
  
})



#UI#----
thisUI = fluidPage(
  titlePanel("Analyze Cardiomyocyte Nucleation and Ploidy"),
#Initial Thresholding Section#----
  fluidRow(
    column(4,
    sliderInput("MeanX", 
                "Minimum Nuclear Mean Intensity Threshold", 
                min = min(nuclei$Mean),
                max = max(nuclei$Mean), 
                value = 10000),
  sliderInput("MaxX", 
              "Upper X limit of Graph", 
              min = min(nuclei$Mean),
              max = max(nuclei$Mean), 
              value = 65535),
  sliderInput("MinX", 
            "Lower X limit of Graph", 
            min = min(nuclei$Mean),
            max = max(nuclei$Mean), 
            value = 1)),
  # Show a plot of the generated distribution
  column(8,
  mainPanel(
    plotOutput("selectMean")
  )
  )
  ),
  
  # Sidebar with a slider input for number of observations
  fluidRow(
    column(4,
           sliderInput("AreaMeanX", 
                       "Minimum Nuclear Area Threshold", 
                       min = 0,
                       max = 6*median(nuclei$Area), 
                       value = median(nuclei$Area)),
           sliderInput("AreaMaxX", 
                       "Upper X limit of Graph", 
                       min = 0,
                       max = 6*median(nuclei$Area), 
                       value = 6*median(nuclei$Area)),
           sliderInput("AreaMinX", 
                       "Lower X limit of Graph",
                       min = 0,
                       max = 6*median(nuclei$Area), 
                       value = 1)),
    # Show a plot of the generated distribution
    column(8,
           mainPanel(
             plotOutput("selectArea")
           )
           )
    ),
  fluidRow(
    column(4,
           sliderInput("FeretMeanX", 
                       "Maximum Cardiomyocyte Minimum Feret's Diameter Threshold", 
                       min = 0,
                       max = 6 * median(cm$MinFeret),
                       value = auto_thresh(as.integer(cm$MinFeret), "Triangle")[1]),
           sliderInput("FeretMaxX", 
                       "Upper X limit of Graph", 
                       min = 0,
                       max = 6 * median(cm$MinFeret),
                       value = 6 * median(cm$MinFeret)),
           sliderInput("FeretMinX", 
                       "Lower X limit of Graph",
                       min = 0,
                       max = 6 * median(cm$MinFeret),
                       value = 1)
           ),
    # Show a plot of the generated distribution
    column(8,
           mainPanel(
             plotOutput("selectFeret")
           )
    )
  ),
fluidRow(actionButton("applyThresholds", "Apply Selected Thresholds?")),
fluidRow(conditionalPanel(
  condition = "input.applyThresholds",
  bsCollapse(bsCollapsePanel("Cell Cycle Stain Intensity Distribution Plot",
                             column(4,            
                                    actionButton("applyCycleIntensityThresh", "Apply Cell Cycle Intensity Threshold"),
                                    sliderInput("CellCycleIntensityMeanX", 
                                                "Cell cycle stain positivity threshold", 
                                                min = min(nuclei$cycleIntensity),
                                                max = max(nuclei$cycleIntensity),
                                                value = max(nuclei$cycleIntensity) / 2)),
                             column(8,plotOutput("cycleIntensity"))
  )))),
fluidRow(conditionalPanel(
  condition = "input.applyCycleIntensityThresh",
  bsCollapse(bsCollapsePanel("DNA Integrated Intensity Distribution Plots",
                             column(4,            
                                    actionButton("plotUnnormalized", "Plot Intensity Distribution"),
                                    sliderInput("CellCycleIntensityMeanX", 
                                                "Cell cycle stain positivity threshold", 
                                                min = min(nuclei$cycleIntensity),
                                                max = max(nuclei$cycleIntensity),
                                                value = max(nuclei$cycleIntensity) / 2),
                                    sliderInput("CellCycleIntensityMaxX", 
                                                "Upper X limit of Graph", 
                                                min = min(nuclei$cycleIntensity),
                                                max = max(nuclei$cycleIntensity),
                                                value = max(nuclei$cycleIntensity)),
                                    sliderInput("CellCycleIntensityMinX", 
                                                "Lower X limit of Graph",
                                                min = min(nuclei$cycleIntensity),
                                                max = max(nuclei$cycleIntensity),
                                                value = min(nuclei$cycleIntensity))),
                             column(8,
                                    plotOutput("unnormalizedCellCycle"),
                                    uiOutput("unnormPlots"))
  )))),
fluidRow(conditionalPanel(
  condition = "input.plotUnnormalized",
  bsCollapse(bsCollapsePanel("Estimated Ploidy Distribution Plots",
                             column(4,            
                                    checkboxInput("normalizeByWell", "Normalize Separately by group?"),
                                    actionButton("normalize", "Calculate Ploidy"),
                                    actionButton("plotNormalized", "Plot Estimated Ploidy Distribution"),
                                    actionButton("ploidyNucleationCalculate", "Calculate Cardiomyocyte Ploidy and Nucleation Distributions"),
                                    sliderInput("Lower", 
                                                "Minimum diploid threshold", 
                                                min = 0,
                                                max = 7, 
                                                value = 1,
                                                step=0.1),
                                    sliderInput("Mid", 
                                                "Threshold between diploid and tetraploid", 
                                                min = 0,
                                                max = 7, 
                                                value = 3,
                                                step=0.1),
                                    sliderInput("Upper", 
                                                "Maximum tetraploid threshold",
                                                min = 0,
                                                max = 7, 
                                                value = 5.5,
                                                step=0.1)),
                                    
                             column(8,
                                    plotOutput("normalizedCellCycle"),
                                    uiOutput("normPlots"))
  )))),
fluidRow(conditionalPanel(
  condition = "input.ploidyNucleationCalculate",
  bsCollapse(bsCollapsePanel("Cell Cycle Stain Intensity Distribution Plot",
                             column(4,            
                                    actionButton("plotAndSave", "Plot and Save Into Output Folder")),
                             column(8,
                                    plotOutput("mononucleatedDistrubution"),
                                    plotOutput("binucleatedDistribution"),
                                    plotOutput("ploidyDistributionByNucleation"))
  ))))
)

shinyApp(ui = thisUI, server = thisServer)






