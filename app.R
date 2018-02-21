### DSF calculations
library(shiny)
library(plotrix)
library(pspline)
library(XLConnect)

##app
ui = fluidPage(
  includeCSS("styles.css"),
  headerPanel("ThermoFluor Protein Tm calculator"),
  sidebarLayout(
    sidebarPanel(
      wellPanel(
        h3("Upload data"),
        fileInput('dataset', label = 'Raw data'),
        fileInput('design', label = 'Experimental design'),
        hr(),
        downloadButton('design_template',"Experimental design template"),
        actionButton('sample_data', 'Use sample data')
      ),
      wellPanel(
        h3("Adjust parameters"),
        sliderInput('df', label = 'Df number for spline', step=1, min=2, max=200, value=50),
        sliderInput('treshold', label = 'Thresholds for Tm calculation', step=1, min=0, max=100, value=c(25,80))
      )
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Table", downloadButton("download_raw", 'Download CSV'),
                 tableOutput('table')),
        tabPanel("Evaluation Plot",
                 tableOutput('Sample_table'),
                 plotOutput('plot', height = '400px', hover = 'hover'),
                 textOutput('coords')),
        tabPanel("Means",
                 textInput('exclude', label = 'Exclude samples (sep by ",")', value = NULL),
                 downloadButton("download_means", 'Download CSV'), tableOutput('table_means'))
        )
      )
    )
)

server = function(input, output, session) {
  options(shiny.maxRequestSize=100^1024^2) 
  odd <- function(x) x[x%%2 != 0]
  
  wb = reactive({
    if(is.null(input$design) & input$sample_data==0) {
      return(NULL)
    } else if (!is.null(input$design) & input$sample_data==0){
      wb = loadWorkbook(input$design$datapath)
    } else {
      wb = loadWorkbook('exp_design_HopQ1.xlsx')
    }
    return(wb)
  })
  
  melts = reactive({
    if(is.null(input$dataset) & input$sample_data==0) {
      return(NULL)
    } else if (!is.null(input$dataset) & input$sample_data==0){
      melts = read.csv(input$dataset$datapath)
    } else {
      melts = read.csv('ca_ni_mg.csv')
    }
    return(melts)
  })
  
  structure = reactive({
      if(is.null(wb())) {
        return(NULL)
      } else {
        exp.design = readWorksheet(wb(), sheet = 1)
        structure = data.frame()
        r = 0
        for (i in 1:16) {
          for (j in 2:25) {
            r=r+1
            structure[r,"Sample"] = as.character(exp.design[i,j])
            structure[r, "col1"] = as.character(exp.design[i,26])
            structure[r, "col2"] = as.character(exp.design[i,27])
            structure[r, "col3"] = as.character(exp.design[i,28])
            structure[r, "col4"] = as.character(exp.design[i,29])
          }
        }
        colnames(structure) = c('Sample', colnames(exp.design)[26:29])
        structure$col1 = NULL
        structure$col2 = NULL
        structure$col3 = NULL
        structure$col4 = NULL
        return(structure) 
      }
  })

  melt_all = reactive({
    if(is.null(melts())) {
      return(NULL)
    } else {
      melts = melts()
      progress <- Progress$new(session, min=1, max=768)
      on.exit(progress$close())
      
      progress$set(message = 'Calculation in progress')
      melt_all = list()
      j=0
      for(i in odd(1:768)){
        j=j+1
        subset = cbind(melts[i+1], melts[i])
        colnames(subset) = c('fluo', 'temp')
        spline_model = sm.spline(subset$temp, subset$fluo, df=input$df)
        subset$fit = predict(spline_model, subset$temp)
        subset$deriv = predict(spline_model, subset$temp, 1)
        melt_all$temp[[j]] = subset$temp
        melt_all$fluo[[j]] = subset$fluo
        melt_all$deriv[[j]] = -subset$deriv
        melt_all$fit[[j]] = subset$fit
        progress$set(value = i)
      }
      return(melt_all)
    }
  })

  structure_Tm = reactive({
    if(is.null(melt_all()) | is.null(structure())) {
      return(NULL)
    } else {
      structure = structure()
      melt_all = melt_all()
      for (j in 1:length(melt_all$fluo)){
        subset = cbind(melt_all$temp[[j]], melt_all$deriv[[j]])
        subset = as.data.frame(subset); colnames(subset) = c('temp', 'deriv')
        Tm = mean(subset[subset$deriv == min(subset[subset$temp < input$treshold[2] &
                                                      subset$temp > input$treshold[1],'deriv']), 'temp'])
        structure[j, 'Tm'] = Tm
        structure[j, 'No'] = j
      }
      structure = structure[c(ncol(structure), 1:(ncol(structure)-1))]
      structure = structure[complete.cases(structure[,c('Sample', 'Tm')]),]
      rows = structure$No
      
      if(is.null(input$plot_Sample)){
        insertUI(
          selector='#Sample_table',
          where = 'beforeBegin',
          ui = selectInput('plot_Sample', label = 'Choose a sample for plotting',
                           choices = rows, selected=rows[1]))
      }
      return(structure)
    }
  })

  output$table = renderTable({
    if(is.null(structure_Tm())) {
      return(NULL)
    } else {
      return(structure_Tm())
    }
  })
    
  output$Sample_table = renderTable({
    if(is.null(input$plot_Sample)) {
      return(NULL)
    } else {
      structure_Tm = structure_Tm()
      structure_Tm[structure_Tm$No == input$plot_Sample,] 
    }
  })
  
  
  output$plot = renderPlot({
    if (is.null(structure_Tm())) {
      return(NULL)
    } else {
      tres_min = input$treshold[1]
      tres_max = input$treshold[2]
      
      j = as.numeric(input$plot_Sample)
      structure_Tm = structure_Tm()
      melt_all = melt_all()
      subset = cbind(melt_all$temp[[j]], melt_all$deriv[[j]], melt_all$fit[[j]], melt_all$fluo[[j]])
      subset = as.data.frame(subset)
      colnames(subset) = c('temp', 'deriv', 'fit', 'fluo')
      Tm = structure_Tm[structure_Tm$No==j, 'Tm']
      label = paste(as.character(structure_Tm[structure_Tm$No==j, 1:(ncol(structure_Tm)-1)]), collapse = ', ')
      
      twoord.stackplot(lx=subset$temp, rx=subset$temp, ldata=cbind(subset$fluo, subset$fit),
                       rdata=subset$deriv, rcol = 'blue', lcol=c('red','black'),
                       ltype=c('p', 'l'), rtype='l', xlab='Temperature [T]', rylab = '-dF/dT',
                       lylab = 'Fluorescence [F]', cex=0.6, pch=1,
                       # panel.first = rect(c(0,tres_max), -1e6, c(tres_min,1e6), 1e6, col='lightgrey', border=NA))                       panel.first = rect(c(0,tres_max), -1e6, c(tres_min,1e6), 1e6, col='lightgrey', border=NA))
                       panel.first = rect(tres_min, -1e6, tres_max, 1e6, col='#ffffb3', border=NA)) 

      abline(v=Tm, lty=4)
      abline(v=tres_min, lty=3)
      abline(v=tres_max, lty=3)

    }
  })
  
  output$coords = renderText({
    if (is.null(input$hover)) {
      return(" ")
    } else {
      temp = paste('Temperature: ', as.character(round(input$hover$x, digits=2)), ' deg. Celcius', sep='')
    }
    return(temp)
  })
  
  final_structure = reactive({
    if (is.null(input$exclude)) {
      return(structure_Tm())
    } else {
      string = input$exclude
      spl = as.numeric(strsplit(string, split=',')[[1]])
      structure_Tm = structure_Tm()
      structure_Tm = structure_Tm[!structure_Tm$No %in% spl,]
      return(structure_Tm)
    }
  })
  
  means = reactive({
    if(is.null(final_structure())){
      return(NULL)
    } else {
      structure = final_structure()
      var_list = structure[2:(ncol(structure)-1)]
      var_list  = as.list(var_list)
      means = aggregate(structure$Tm, by = var_list, FUN=mean)
      colnames(means)[ncol(means)]="mean"
      means$SD = aggregate(structure$Tm, by = var_list, FUN=sd)$x
      return(means)
    }
  })
    
  output$table_means = renderTable({
    if(is.null(means())) {
      return(NULL)
    } else {
      return(means())
    }
  })
  
  output$download_raw = downloadHandler(
    filename = function() {'Raw_Tm.csv'},
    content = function(file){
      write.csv(structure_Tm(), file)
    })
  
  output$download_means = downloadHandler(
    filename = function() {'Mean_Tm.csv'},
    content = function(file){
      write.csv(means(), file)
    })
  
  output$design_template = downloadHandler(
    filename = function() {'exp_design.xlsx'},
    content = function(file){
      file.copy('exp_design.xlsx', file)
    })
}
shinyApp(ui = ui, server = server)