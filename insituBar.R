library(reshape2)
library(lattice)
source('data/dirfns.R')

# define panel function
# This function will be used to plot each panel of the plot.
# lattice will automatically group x and y values by panel
# the panel function can also accept a subscripts argument, 
# which indicates which values from the original data.frame
# correspond to the panel values.
# It is not necessary to subset x and y in the panel function,
# but any additional column values passed to this function
# must be subset by subscripts.
# In this case we want the function to also accept a lower limit (l)
# and an upper limit (u).
panel.barchart.error <- function(
  x,y,subscripts,l,u,...
){
  # the barplot
  panel.barchart(
    x,y,
    ...,
    col=c('gray70','gray45','gray25')
  )
  # draw error bars on the barplot
  # we want to draw a vertical line from the lower limit to the upper limit
  panel.arrows(
  # The x-values will be the same as the bars.
    x0=x,
    x1=x,
  # The y-values should correspond to the confidence intervals of each variable.
    y0=l[subscripts],
    y1=u[subscripts],
    col='black',
  # To add a "T" to the end of the segments, we can plot an arrow with angle=90.
    angle = 90, 
    length = 1, unit = "mm", ends = "both"
  )
}

# define main function
# lattice functions accept a formula (x) corresponding to variables in a data.frame (dat)
# for our function,the formula should follow the syntax:
# percentage~response|experiment
# we define this function to also accept a vector of lower (l) and upper (u) values for error bars
# These values should correspond to the rows of dat
barchart.error <- function(
  x,dat,l,u,...
) barchart(
  x,
  data = dat,
  # the panel argument should be a function that accepts arguments 
  # x,y,subscripts,...
  # Here we define a function that calls panel.barchart.error. 
  # x,y,and subscripts can vary between panels
  # l and u are always the same values (which values are plotted 
  # depends on the subscripts argument).
  panel=function(
    x,y,
    subscripts,
    ...
  ) panel.barchart.error(
    x,y,
    subscripts,
    l=l,u=u,
    ...
  ),
  # Do not add a title for each panel
  strip=F,
  # draw bars vertically
  horizontal = F,
  # do not plot x-axis (but if it were drawn, plot it perpendicular to the axis)
  scales=list(
    x=list(rot=90,draw=F)
  ),
  # do not plot y-axis label
  ylab=NULL
)

# read data as a list of tables
insitu <- lapply(
  c('st24_ebf_enha.csv','st25_ebf_enha_2.csv',"st27_ebf_enha.csv"),
  read.csv,
  header=T,
  row.names=1
)

# proportion of row counts in each cell
insitu.p <- lapply(
  insitu,
  function(x) x/rowSums(x)
)
# standard deviation of proportion
insitu.sd <- mapply(
  function(x,y) sqrt(x*(1-x)/rowSums(y)),
  # x-values
  insitu.p,
  # y-values
  insitu,
  # return a list, not a data.frame
  SIMPLIFY = F
)
# lower limit of error bars
insitu.lower <- mapply(
  '-',
  insitu.p,
  insitu.sd,
  SIMPLIFY = F
)
# upper limit of error bars
insitu.upper <- mapply(
  '+',
  insitu.p,
  insitu.sd,
  SIMPLIFY = F
)

# convert data to "long" format and append error limits
insitu.dat <- mapply( 
  # define an anonymous function to merge the proportion, lower and upper limits
  function(p,l,u) Reduce(
    # define an anonymous function to merge on the combination of columns "Var1" and "Var2"
    # These names correspond to the default behavior of melt() for matrices
    function(x,y) merge(
      x,y,
      c("Var1","Var2")
    ),
    # melt tables before merging them
    lapply(
      list(p,l,u),
      function(x) melt(
        # converting a data.frame to a matrix before melting preserves the row.names
        as.matrix(x)
      )
    )
  ),
  insitu.p,
  insitu.lower,
  insitu.upper,
  SIMPLIFY = F
)

# apply barchart.error to each table
mapply(
  function(dat,file) {
    # open connection to image file
    dir.eps(
      file,
      width=3,
      height=2
    )
    # define plot
    res <- barchart.error(
      value.x*100 ~ Var2|Var1,
      dat,
      dat[,'value']*100,
      dat[,'value.y']*100,
      layout=c(6,1)
    )
    # write plot to connection
    plot(res)
    # close connection
    dev.off()
  },
  insitu.dat,
  # file names
  paste0(c('st24',"st25","st27"),'Bar')
)

