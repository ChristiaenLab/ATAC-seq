# library(shiny)
# library(RSQLite)

# global variables here are (NOT!) visible in ui.R
# con <- dbConnect(SQLite(),'atacCiona.db')
# rna.atac.res <- dbReadTable(con,'rnaseq_atac_lib')
# 
# addInput <- function(i,...) renderUI({
#   
# })

shinyServer(function(input, output) {
  query <- reactive({
    input$SQL
  })
  sql <- reactive({
    dbGetQuery(con,input$SQL)
  })
  write.query <- reactive({
    write.table(sql(),input$file.txt,sep='\t',quote=F,row.names = F)
  })
  write.geneToPeak <- reactive({
    write.table(geneToPeak(con,isolate(sql())),input$file.txt,sep='\t',quote=F,row.names = F)
  })
  write.peakToGene <- reactive({
    write.table(peakToGene(con,isolate(sql())),input$file.txt,sep='\t',quote=F,row.names = F)
  })
  write.diamond <- reactive({
    dbDiamond(con,input$genes,input$peaks,input$diamond.out)
  })
  # tmp <- mapply(dbReadTable,rna.atac.res$rnaseqRes,MoreArgs = list(conn=con,row.names="GeneID"))
  # print(output)
  sapply(
    unlist(rna.atac.res),
    function(x) output[[x]] <<- renderDataTable({
      dbReadTable(con,x)
    })#,options = list(aLengthMenu = c(5, 30, 50), iDisplayLength = 5))
  )
  output$SQL <- renderDataTable({
    input$update
    isolate(sql())
  })
  output$genes <- renderDataTable({
    input$geneToPeak
    # tmp <- isolate(sql())
    geneToPeak(con,isolate(sql()))
  })
  output$peaks <- renderDataTable({
    input$peakToGene
    # tmp <- isolate(sql())
    peakToGene(con,isolate(sql()))
  })
  output$tables <- renderDataTable({
    data.frame(
      tables=dbListTables(con),
      fields=sapply(dbListTables(con),function(x) paste(dbListFields(con,x),collapse = ', '))
    )
  })
  # observe({
  #   if(input$ngeneset>0) {
  #     gsid <- paste0('geneset',as.character(ngeneset()))
  #     output[[gsid]] <<- renderUI({
  #       textAreaInput(gsid,gsid)
  #     })
  #   } else NULL
  # })
  observe({
    if(input$write.tab>0){
      write.query()
    } else NULL
  })
  observe({
    if(input$writeGenePeak>0){
      write.geneToPeak()
    } else NULL
  })
  observe({
    if(input$writePeakGene>0){
      write.peakToGene()
    } else NULL
  })
  observe({
    if(input$writePlot>0){
      write.diamond()
    } else NULL
  })
})
