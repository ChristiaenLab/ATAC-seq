# library(shiny)

# Define UI for application
shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel("Data"),
  
  do.call(sidebarPanel,append(
    list(
      textAreaInput('SQL',"SQL",value=dbGetQuery(con,"SELECT query FROM peaksets")[1,],height = '200px'),
      actionButton("update" ,"Show query", icon("refresh"), class = "btn btn-primary"),
      textInput('file.txt',"File name"),
      actionButton('write.tab','Write query to file'),
      actionButton('geneToPeak','Convert GeneID to PeakID'),
      actionButton('peakToGene','Convert PeakID to GeneID'),
      actionButton('writeGenePeak','Write gene-to-peak file'),
      actionButton('writePeakGene','Write peak-to-gene file')
      # actionButton("ngeneset" ,"Add gene set")
    ),
    # ignore this line
    # lapply(paste0('genesets',as.character(1:10)),uiOutput)
    list(
      h4("diamondplot"),
      textAreaInput('genes',"Gene sets",height='200px'),
      textAreaInput('peaks',"Peak sets",height='200px'),
      selectInput('rna','RNA-seq data',names(rna)),
      selectInput('atac','ATAC-seq data',names(atac)),
      selectInput("geneset","Gene set", c(names(scrna),names(bulkGS))),
      checkboxGroupInput('peakset',"Peak sets",names(peaksets)),
      checkboxInput(
        'intersect',
        'Plot only the intersect of the gene set and peak sets?'
      ),
      textInput('diamond.out',"File name",'diamond.eps'),
      actionButton('writePlot','Write diamondplot')
    )
  )),
  
  do.call(
    mainPanel,
    # append(
      list(
        h4("SQL"),
        tabsetPanel(
          tabPanel("query",dataTableOutput("SQL")),
          tabPanel("tables",dataTableOutput("tables")),
          tabPanel("Gene-to-peak",dataTableOutput("genes")),
          tabPanel("Peak-to-gene",dataTableOutput("peaks"))
        ),
    # show ATACseq and RNAseq as two tabsets with headers
        h4("RNA-seq result"),
        do.call(
          tabsetPanel,
          lapply(
            rna.atac.res$rnaseqRes,
            function(x) tabPanel(
              x,dataTableOutput(x)
            )
          )
        ),
        h4("ATAC-seq result"),
        do.call(
          tabsetPanel,
          lapply(
            rna.atac.res$atacRes,
            function(x) tabPanel(
              x,dataTableOutput(x)
            )
          )
        )
      )
  )
))
