setRowNames <- function(x,nm){
  if(is.data.frame(x)){ 
    row.names(x) <- nm
  }else rownames(x) <- nm
  return(x)
}

paste1 <- function(x,sep='',collapse=NULL){
  x <- lapply(x,as.character)
  x$sep <- sep
  x$collapse <- collapse
  return(do.call(paste,x))
} 
# simplifies Reduce syntax
sreduce <- function(fn,x,init=NULL, right = FALSE, accumulate = FALSE,...) Reduce(
  function(y,z) fn(y,z,...),x,init,right,accumulate)

lfold <- function(x,fn,...) {
  res <- x[[1]]
  lapply(2:length(x), function(y) res <<- fn(res,x[[y]],...))
  return(res)
}

# finds all peaks associated to all genes
pgrep <- function(genes,peaks) unique(unlist(sapply(genes,grep,peaks,fixed=T)))

# concatenates path and filename into output and writes a data.frame to csv.
# any folders in the path that do not exist are created.
dir.in <- function(filename, path = '.',...){
  filename <- paste(path, filename, sep = '/')
  read.delim(filename, ...)
}
mkdate <- function(filename,ext,path='.',append.date=T){
  if(append.date){
    if(grepl('^[~/\\.]',path)) path <- paste(path, Sys.Date(), sep = '/')
    else path <- paste(Sys.Date(), path, sep = '/')
  }
  filename <- paste0(path,'/' ,filename)
  path <- sub('(^.*\\/).*',"\\1",filename)
  if(!dir.exists(path)) dir.create(path,recursive = T)
  if(ext!='') filename <- paste0(filename,'.',ext)
  return(filename)
}
dir.out <- function(x,fn,filename,ext='txt',path='.',...,append.date=T){
  filename <- mkdate(filename,ext,path,append.date)
  fn(x,filename,...)
}
dir.img <- function(filename, fn,ext, path = '.', ...,append.date=T){
  filename <- mkdate(filename,ext,path,append.date)
  fn(filename, ...)
}
dir.tab <- function(x,filename, path = '.',ext='txt',quote=F,...,append.date=T){
  dir.out(x,write.table,filename,ext,path,sep='\t', quote=quote,...,append.date=append.date)
  # path <- paste(path, Sys.Date(), sep = '/')
  # # path <- get.filename(path)
  # if(!dir.exists(path)) dir.create(path,recursive = T)
  # filename <- paste(path, filename, sep = '/')
  # # filename <- get.filename(filename)
  # write.table(x,paste(filename,'txt',sep='.'),sep='\t', quote=F,...)
  # # if(is.data.frame(x)&summary) write.csv(summary(x),paste(filename,'summary.csv',sep = '.'))
}
dir.csv <- function(x,filename, path = '.', summary=F,quote=T,...,append.date=T){
  dir.out(x,write.csv,filename,'csv',path,quote=quote,...,append.date=append.date)
  # path <- paste(path, Sys.Date(), sep = '/')
  # if(!dir.exists(path)) dir.create(path,recursive = T)
  # filename <- paste(path, filename, sep = '/')
  # write.csv(x,paste(filename,'csv',sep='.'),quote=T,...)
}
dir.bed <- function(x,filename, path = '.',meta=NULL,...,append.date=T){
  sel <- !is.na(x)
  x <- x[sel]
  x <- cbind(do.call(rbind,strsplit(x,':')),x)
  x <- cbind(x,'.','*',meta[sel,],...)
  # path <- paste(path, Sys.Date(), sep = '/')
  dir.tab(x,filename,path,ext='bed',col.names = F,row.names = F,append.date=append.date)
  # if(!dir.exists(path)) dir.create(path,recursive = T)
  # filename <- paste(path, filename, sep = '/')
  # write.table(x,paste(filename,'bed',sep='.'),sep='\t', quote = F,col.names = F,row.names = F)
}

dir.export <- function(x,file,path='.',format='bed',...,append.date=T){ 
  require(rtracklayer)
  dir.out(x,export,file,ext=format,path,format,...,append.date=append.date)
}

# writes granges to bed file containing all metadata columns
dir.GRanges <- function(x,file,path='.',...,append.date=T){
  require(GenomicRanges)
  res <- data.frame(
    seqnames=seqnames(x),
    start=start(x),
    end=end(x),
    strand=strand(x),
    name=names(x)
  )
  res <- cbind(res,mcols(x))
  dir.tab(
    res,file,path=path,row.names=F,ext='bed',
    ...,append.date = append.date
  )
}

# export DNAString object
dir.DNAString <- function(x,file,path='.',...,append.date=T){
  mcols(x)$seq <- as.character(x)
  dir.GRanges(x,file,path = path,...,append.date = append.date)
}

dir.gsm <- function(x,filename,...,append.date=T) {
  dir.out(x,filename,...,append.date=append.date)
}
dir.png <- function(filename, path = '.', ...,append.date=T) dir.img(
  filename,png,'png',path,res=300,width=2000,height=2000,...,append.date=append.date
)
dir.pdf <- function(filename, path = '.', ...,append.date=T) dir.img(
  filename,pdf,'pdf',path,...,append.date=append.date
)
dir.svg <- function(filename, path = '.', ...,append.date=T) dir.img(
  filename,svg,'svg',path,...,append.date=append.date
)
dir.eps <- function(filename,path='.',...,append.date=T) {
  setEPS()
  dir.img(filename,postscript,'eps',path,...,append.date=append.date)
}

dir.gg <- function(x,filename,path='.',ext='pdf',...,append.date=T) {
  require(ggplot2)
  filename <- mkdate(filename,ext,path,append.date)
  ggsave(filename,x,ext,...)
}

resBed <- function(res,path)  lapply(names(res),function(x) {
  bed <- res[[x]][order(res[[x]]$log2FoldChange),]
  dir.bed(row.names(bed),x,path,bed$log2FoldChange)
})
read.bed <- function(x) {
  x <- read.table(x)
  row.names(x) <- do.call(paste,c(x[,1:3],list(sep=':')))#paste1(x[,1:3],':')
  # apply(read.table(x,colClasses = 'character')[,1:3],1,paste,collapse=':')
  if(ncol(x)>4) return(x[,-1:-3]) else return(row.names(x))
}

# dir.png <- function(filename, path = '~/ciona', ...) dir.img(
#   paste0(filename,'.png'),png,path,res=300,width=2000,height=2000
# )
# dir.pdf <- function(filename, path = '~/ciona', ...) dir.img(
#   paste0(filename,'.pdf'),pdf,path
# )
# dir.svg <- function(filename, path = '~/ciona', ...) dir.img(
#   paste0(filename,'.svg'),svg,path
# )
# dir.eps <- function(filename,path='~/ciona',...) {
#   setEPS()
#   dir.img(paste0(filename,'.eps'),postscript,path,...)
# }
dir.venn <- function(x,filename,path){
  require(VennDiagram)
  filename <- mkdate(filename,'tiff',path)
  venn.diagram(x,filename)
}
ann.venn <- function(x,file,path,feat='ann',rows='Row.names'){
  feats <- lapply(unique(x[,feat]),function(y) x[x[,feat]==y,rows])
  names(feats) <- unique(x[,feat])
  dir.venn(feats,file,path)
}
dir.apply <- function(x,path,fn=dir.tab,...) sapply(names(x),function(y) fn(x[[y]],y,path,...))

lrtab <- function(dir,fn=read.table,pattern=NULL,...) {
  res <- lapply(
    list.files(dir,pattern,full.names = T),
    function(x) if(file.size(x)>0) fn(x,...) else NULL)
  names(res) <- sub('\\.[A-Za-z0-9]*$','',list.files(dir,pattern))
  return(res)
}

# apply functions rewritten to use parallel processing
apply2 <- function(...,type = "FORK",threads=detectCores()-1) {
  require(BiocGenerics)
  cl <- makeCluster(threads,type = type,outfile='')
  result <- parApply(cl,...)
  stopCluster(cl)
  return(result)
}

lapply2 <- function(...,type = "FORK",threads=detectCores()-1) {
  require(BiocGenerics)
  cl <- makeCluster(threads,type = type,outfile='')
  result <- parLapply(cl,...)
  stopCluster(cl)
  return(result)
}

sapply2 <- function(...,type = "FORK",threads=detectCores()-1) {
  require(BiocGenerics)
  cl <- makeCluster(threads,type = type,outfile='')
  result <- parSapply(cl,...)
  stopCluster(cl)
  return(result)
}