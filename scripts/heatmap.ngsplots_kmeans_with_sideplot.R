#clustering is kmeans then sorted within clusters



#source('jrb_R_scripts//heatmap.3-split.R')
#source('heatmap.2.2.R')
library('png')
library('reshape')
library('ggplot2')
library('gplots')




startRasterMode = function(width = 1, height = 1){
  png('tmp.png', units = 'in', width = width, height = height, res = 600)
  par(mai = rep(0,4))
}

stopRasterMode = function(mai = NA){
  dev.off()
  if(length(mai) > 1 && !is.na(mai)){
    par(mai = mai)
  }
  plot0()
  rasterImage(readPNG("tmp.png", native = FALSE), xleft = 0, xright = 1, ytop = 1, ybottom = 0,
              interpolate = FALSE)
}

heatmap.ngsplots = function(ngs_profiles, 
                            main_title = NULL, 
                            lmat_custom = NULL,
                            profiles_to_plot = NA, 
                            nclust = 6, 
                            labels_below = NA, 
                            labels_above = NA, 
                            fg_toPlot = character(), 
                            labels_right = NA, 
                            sortClustersByTotal = F, 
                            hmapColors = c('royalblue1', 'black', 'yellow'), 
                            labelWithCounts = F, 
                            fg_label = NA, 
                            label_clusters = T,
                            key.lab = 'log2 FE',
                            key.cex = 1.3,
                            cex.main = 2.5,
                            cex.row = 1.3,
                            cex.col = 2,
                            doSidePlot = T,
                            sidePlot_smoothing = 1,
                            extraData = NULL,
                            extraData_colors = NULL,
                            extraData_plotFunction = NA,
                            forPDF = T,
                            globalScale = 1){
  if(length(profiles_to_plot) == 1 && is.na(profiles_to_plot)){
    profiles_to_plot = names(ngs_profiles)
  }
  if(length(labels_below) < 2 && is.na(labels_below)){
    labels_below = profiles_to_plot
  }
  if(length(profiles_to_plot) %% length(labels_below) != 0){
    stop('length of labels_below must divide profiles_to_plot evenly!')
  }
  if(length(profiles_to_plot) %% length(labels_above) != 0){
    stop('length of labels_above must divide profiles_to_plot evenly!')
  }
  
  
  prof = matrix(0, ncol = 0, nrow = nrow(ngs_profiles[[1]]))#assemble single matrix by joining selected profiles (matrices)
  hidden = lapply(profiles_to_plot, function(x){
    prof <<- cbind(prof, ngs_profiles[[x]])
  })
  
  cr = colorRamp(hmapColors)
  tmp = 0:100 / 100
  colors = rgb(cr(tmp)/255)
  
  len = length(labels_below)
  len_major = length(labels_above)
  cseps = 1:(len-1) * ncol(prof) / len #column separators are between joined profiles
  if(len == 1) cseps = -1
  cseps_major = 1:(len_major-1) * ncol(prof) / len_major
  if(len_major == 1) cseps_major = -1
  
  labels.above = rep('', ncol(prof))#column labels are centered above each profile
  if(!is.na(labels_above[1])){
    labels_loc = round((1:length(labels_above) -.499) * ncol(prof) / length(labels_above))
    labels.above[labels_loc] = labels_above
  }
  
  
  labels.below = rep('', ncol(prof))
  labels_loc = round((1:length(labels_below)-.499) * ncol(prof) / length(labels_below))
  labels.below[labels_loc] = labels_below
  
  set.seed(1)
  kmclust = kmeans(prof, centers = nclust, iter.max = 10)
  
  
  
  if(length(fg_toPlot) > 0){#extract fg_plot as special cluster
    #rseps = c(0, rseps)#add an
    kmclust$cluster[fg_toPlot] = 0#change cluster of fg
    kmclust$cluster = kmclust$cluster + 1
    
    
    nclust = nclust + 1
    kmclust$size = sapply(1:nclust, function(x){
      return(sum(kmclust$cluster == x))
    })
    kmclust$centers = rbind(colMeans(prof[fg_toPlot,]), kmclust$centers)
    rownames(kmclust$centers) = 1:nclust
  }
  
  #   o = order(kmclust$cluster)#sort prof to be in cluster number order
  #   prof = prof[o,, drop = F]
  #   kmclust$cluster = kmclust$cluster[o]#keep cluster assignment sorted in same order as prof
  
  init_clust_data = function(clust, data){
    #initize cluster object
    #add data
    #sort by cluster id
    #add nclust
    o = order(clust$cluster)
    clust$data = data[o,, drop = F]
    names(clust$cluster) = NULL
    clust$cluster = clust$cluster[o]
    clust$nclust = length(clust$size)
    return(clust)
  }
  fetch_clust = function(clust, i){
    #returns the ith cluster
    if(length(clust$size) < i) stop('fetch_clust i greater than nclust')
    keep = clust$cluster == i
    return(clust$data[keep,, drop = F])
  }
  move_clust = function(clust, i_tomove, i_dest){
    #removes i_tomove and inserts it at i_dest
    #updates data, cluster, centers, and size
    if(clust$nclust < i_tomove) stop('move_clust i_tomove greater than nclust')
    if(clust$nclust < i_dest) stop('move_clust i_dest greater than nclust')
    curr_o = 1:clust$nclust
    new_o = curr_o[curr_o != i_tomove]
    new_o = c(new_o[1:length(new_o) < i_dest], i_tomove, new_o[1:length(new_o) > i_dest - 1]); new_o
    new_data = matrix(0, nrow = 0, ncol = ncol(clust$data))
    new_cluster = integer()
    hidden = sapply(new_o, function(x){
      keep = clust$cluster == x
      new_data <<- rbind(new_data, clust$data[keep,, drop = F])
      new_cluster <<- c(new_cluster, clust$cluster[keep])
    })
    clust$data = new_data
    clust$cluster = new_cluster
    clust$centers = clust$centers[new_o,]
    clust$size = clust$size[new_o]
    
    dict = 1:clust$nclust
    names(dict) = new_o
    
    clust$cluster = dict[as.character(clust$cluster)]
    rownames(clust$centers) = 1:clust$nclust
    return(clust)
  }
  replace_clust = function(clust, i, new_datai){
    #essentially removes cluster i and insert new_data in its place
    clust$size[i] = nrow(new_datai)
    clust$centers[i,] = colMeans(new_datai)
    before = clust$cluster < i
    after = clust$cluster > i
    data = rbind(clust$data[before,,drop = F], new_datai, clust$data[after,, drop = F])
    clust$data = data
    clusters = c(clust$cluster[before], rep(i, nrow(new_datai)) ,clust$cluster[after])
    clust$cluster = clusters
    return(clust)
  }
  
  kmclust = init_clust_data(kmclust, prof)
  
  for(i in 1:nclust){#sort within each cluster by rowSum
    dat = fetch_clust(kmclust, i)
    o = order(rowSums(dat), decreasing = T)
    dat = dat[o,,drop = F]
    kmclust = replace_clust(kmclust, i, dat)
  }
  
  prof_ordered = matrix(0, ncol = ncol(prof), nrow = 0)#prof_ordered contains profiles rearranged at cluster level
  get_kmclust = function(i){
    start = 1
    if(i > 1){
      start = sum(kmclust$size[1:(i-1)]) + 1
    }
    end = sum(kmclust$size[1:(i)])
    return(prof[start:end,,drop = F])
  }
  
  kmclust_order = 1:nclust
  if(length(fg_toPlot) > 0){
    sortByDist = F
    if(sortByDist){
      km_dist = as.matrix(dist(kmclust$centers))
      o = order(km_dist[,1])
      kmclust_order = o
    }else{
      kmclust_order = order(apply(kmclust$centers, 1, sum), decreasing = T)
      not_1 = kmclust_order != 1#move 1 to beginning
      kmclust_order = c(1, kmclust_order[not_1])
    }
    
    #    
  }else{
    if(sortClustersByTotal){#sort clusters by total of their centers
      kmclust_order = order(apply(kmclust$centers, 1, sum), decreasing = T)
    }else{
      hiclust = hclust(dist(kmclust$centers))#use hierarchical clustering on centers
      kmclust_order = hiclust$order
    }
  }
  
  rColorChoices = RColorBrewer::brewer.pal(nclust, 'Pastel2')#set cluster id colors
  if(length(fg_toPlot) > 0){#correct for special selection color
    rColorChoices = c('white', rColorChoices[2:length(rColorChoices)-1])
  }
  
  
  for(i in 1:length(kmclust_order)){#sort k means clusters 
    tomove = kmclust_order[i]
    dest = i
    kmclust = move_clust(kmclust, tomove, dest)
    kmclust_order = ifelse(kmclust_order <= i, kmclust_order + 1, kmclust_order)#position of clusters less than i must be increased by 1
    kmclust_order[1:i] = 1:i#up to i is sorted
    #     
    #     
    #     new_cluster = get_kmclust(kmclust_order[i])
    #     prof_ordered = rbind(prof_ordered, new_cluster)
    #     
    #     cluster_color = rColorChoices[i]
    #     rColors = c(rColors, rep(cluster_color, nrow(new_cluster)))
  }
  rColors = rColorChoices[kmclust$cluster]
  names(rColors) = rownames(kmclust$data)
  
  
  rseps = 1:(nclust-1)#cacluated row seperations
  #kmclust$size = kmclust$size[kmclust_order]
  
  for(i in 1:nclust){
    size = kmclust$size
    start = 1
    if(i > 1){
      start = sum(kmclust$size[1:(i-1)]) + 1
    }
    end = sum(kmclust$size[1:(i)])
    if(i < nclust){
      rseps[i] = end
    }
  }
  
  if(labelWithCounts){
    if(is.null(labels_right)){
      labels_right = 1:nclust
    } else {
      labels_right = paste(kmclust$size, labels_right)
    }
  }else{
    if(length(labels_right) == 1){
      labels_right = rep(labels_right, nclust)
    }else if(length(labels_right) != nclust){
      stop('length of labels_right must be 1 or = to # of clusters (+1 if fg_toPlot is specified)')
    }
  }
  if(length(fg_toPlot) == 0 && !is.na(fg_label)){
    stop('cannot add fg_label, fg_toPlot is empty')
  }
  RowSideLabels = character()
  if(label_clusters){
    RowSideLabels = 1:nclust#paste('cluster', 1:nclust)
    if(length(fg_toPlot > 0)){
      if(!is.na(fg_label)){
        RowSideLabels = c(fg_label, 1:(nclust-1))# paste("cluster", 1:(nclust-1)))#special label for fg_toPlot
      }
    }
  }else if(!is.na(fg_label)){
    RowSideLabels = c(fg_label, rep('', nclust-1))#special label for fg_toPlot
  }
  extra_spacer = 0
  if(length(fg_toPlot) > 0){
    extra_spacer = 2
  }
  args = list(
    x = kmclust$data,key.lab = key.lab, key.cex = key.cex, 
    RowSideLabels = RowSideLabels,
    RowSideColors = rColors,
    labels.below = labels.below, cexRow = cex.row, cexCol = cex.col, col = colors, 
    rowsep.major = c(rep(rseps[1], extra_spacer), rseps), 
    colsep.minor = cseps, 
    colsep.major = cseps_major,
    sepwidth.minor = .01, 
    sepwidth.major = .05, 
    labels.above = labels.above, 
    na.color = 'red', 
    labels_rowsep = c(labels_right[1], rep('', extra_spacer), labels_right[2:length(labels_right)]), 
    main = main_title, cex.main = cex.main
  )
  
  hidden = heatmap.2.2(args$x, key.lab = args$key.lab, args$key.cex, #key.par = list(mai = c(.5,0,0,0)),
                       clust = kmclust, globalScale = globalScale, forPDF = forPDF,
                       RowSideLabels = args$RowSideLabels,
                       RowSideColors = args$RowSideColors,
                       labels.below = args$labels.below, cexRow = args$cexRow, cexCol = args$cexCol, col = args$col, 
                       rowsep.major = args$rowsep.major, 
                       colsep.minor = args$colsep.minor, 
                       colsep.major = args$colsep.major,
                       sepwidth.minor = args$sepwidth.minor, 
                       sepwidth.major = args$sepwidth.major, 
                       labels.above = args$labels.above, 
                       na.color = args$na.colors, 
                       labels_rowsep = args$labels_rowsep,
                       main = args$main,
                       cex.main = args$cex.main,
                       doSidePlot = doSidePlot, 
                       sidePlot_smoothing = sidePlot_smoothing,
                       extraData = extraData,
                       extraData_plotFunction = extraData_plotFunction,
                       lmat_custom = lmat_custom, 
                       extraData_colors = extraData_colors)
  cluster_members = sapply(1:kmclust$nclust, function(x){
    return(rownames(kmclust$data)[kmclust$cluster == x])
  }) 
  ngs_res = list(class_sizes = sapply(cluster_members, length), cluster_members = cluster_members, colors = unique(rColors), as_plotted = kmclust$data, kmclust = kmclust, args = args)
  return(ngs_res)
}

heatmap.replot_ngsplots = function(ngs_profiles, hmap_res, labels_above, labels_below){
  
  args = hmap_res$args
  
  prof = matrix(0, ncol = 0, nrow = nrow(ngs_profiles[[1]]))#assemble single matrix by joining selected profiles (matrices)
  hidden = lapply(ngs_profiles, function(x){
    prof <<- cbind(prof, x)
  })
  
  len = length(ngs_profiles)
  
  cseps = 1:(len-1) * ncol(prof) / len #column separators are between joined profiles
  
  labels.above = rep('', ncol(prof))#column labels are centered above each profile
  if(!is.na(labels_above[1])){
    labels_loc = round((1:length(labels_above) -.499) * ncol(prof) / length(labels_above))
    labels.above[labels_loc] = labels_above
  }
  
  
  labels.below = rep('', ncol(prof))
  labels_loc = round((1:length(labels_below)-.499) * ncol(prof) / length(labels_below))
  labels.below[labels_loc] = labels_below
  
  prof = prof[rownames(args$x),]
  
  
  hidden = heatmap.2.2(x = prof, key.lab = args$key.lab, key.cex = args$key.cex, #key.par = list(mai = c(.5,0,0,0)),
                       clust = kmclust, globalScale = globalScale, forPDF = forPDF,
                       RowSideLabels = args$RowSideLabels,
                       RowSideColors = args$RowSideColors,
                       labels.below = labels.below,
                       labels.above = labels.above,
                       cexRow = args$cexRow,
                       cexCol = args$cexCol,
                       col = args$col, 
                       rowsep.major = args$rowsep.major, 
                       colsep.minor = cseps, 
                       sepwidth.minor = args$sepwidth.minor, 
                       sepwidth.major = args$sepwidth.major, 
                       
                       na.color = args$na.color, 
                       labels_rowsep = args$labels_rowsep, 
                       main = args$main,
                       cex.main = args$cex.main,
                       extraData_plotFunction = extraData_plotFunction,
                       doSidePlot = doSidePlot, sidePlot_smoothing = sidePlot_smoothing, extraData = extraData, lmat_custom = lmat_custom, extraData_colors = extraData_colors)
}

#dev.off()
#dummy out tracing, dummy out clustering
heatmap.2.2 = function (x,
                        clust,
                        col, 
                        lmat_custom = NULL,
                        globalScale = 1,
                        forPDF = T,
                        #dividing up the plot
                        colsep.minor = -5, 
                        colsep.major = -1, 
                        rowsep.minor = -1, 
                        rowsep.major = -1, 
                        sepwidth.minor = .01, 
                        sepwidth.major = .05, 
                        na.color = par("bg"), 
                        
                        #color and label for clusters
                        RowSideLabels = NA,
                        RowSideColors = NULL, 
                        
                        #title
                        main = NULL,
                        cex.main = 2,
                        
                        #labelling rows and columns
                        cexRow = 0.2 + 1/log10(nr), 
                        cexCol = 0.2 + 1/log10(nc),
                        labels.below = NULL, 
                        labels.above = NULL, 
                        labels_rowsep = NULL,
                        
                        #left margin size
                        left_mai = .2,
                        
                        #color key
                        key = T, 
                        key.height = 1.5,
                        key.lab = 'Color Key', 
                        key.cex = 1,
                        key.xtickfun = NULL, 
                        key.par = list(),
                        
                        #side plot of summary
                        doSidePlot = T,
                        sidePlot_smoothing = 1,
                        extraData = NULL,
                        extraData_colors = NULL,
                        extraData_plotFunction = NA) 
{
  x.original = x
  rowsep.minor = -1
  rowsep.minor.original = rowsep.minor
  rowsep.major.original = rowsep.major
  colsep.minor.original = colsep.minor
  colsep.major.original = colsep.major
  #sepwidth.minor = round(sepwidth.minor)
  #sepwidth.major = round(sepwidth.major)
  
  scale01 <- function(x, low = min(x), high = max(x)) {
    x <- (x - low)/(high - low)
    x
  }
  
  
  
  dev.width = par('din')[1]
  dev.height = par('din')[2]
  nclust = length(clust$size)
  
  body_height = dev.height
  body_width = dev.width
  body_iy = 1:nclust
  body_ix = 1:2
  lmat = matrix(rep(1, nclust), ncol = length(body_ix), nrow = nclust)
  key_width = min(4,dev.width/3)
  lwid = c(key_width, body_width - key_width)
  lhei = rep(body_height /nclust, nclust)
  
  label.height = .5
  main.height = .8
  extraDataSize = 1
  divLinesSize = .3
  sidePlotSize = 1
  labRowSize = 1
  
  add_lmat_left = function(added_width){
    lmat <<- cbind(rep(0, nrow(lmat)), lmat)
    added_width = added_width * globalScale
    lwid[max(body_ix)] <<- lwid[max(body_ix)] - added_width
    lwid <<- c(added_width, lwid)
    body_ix <<- body_ix + 1
  }
  add_lmat_right = function(added_width, solid = T){
    if(solid){
      lmat <<- cbind(lmat, c(rep(0, min(body_iy)-1), rep(max(lmat)+1, nclust), rep(0, nrow(lmat)-max(body_iy))))
    }else{
      lmat <<- lmat <- cbind(lmat, c(rep(0, min(body_iy)-1), (max(lmat)+1):(max(lmat) + nclust), rep(0, nrow(lmat)-max(body_iy))))
    }
    added_width = added_width * globalScale
    lwid[max(body_ix)] <<- lwid[max(body_ix)] - added_width
    lwid <<- c(lwid, added_width)
  }
  
  add_lmat_top = function(added_height){
    lmat <<- rbind(max(lmat) + 1, lmat)
    added_height = added_height * globalScale
    lhei[body_iy] <<- lhei[body_iy] - added_height / nclust
    lhei <<- c(added_height, lhei)
    body_iy <<- body_iy + 1
  }
  add_lmat_bottom = function(added_height, body_xpos = -1){
    if(body_xpos > 0){
      new_row = rep(0, ncol(lmat))
      new_row[body_ix[body_xpos]] =  max(lmat) + 1
      lmat <<- rbind(lmat, new_row)
    }else{
      lmat <<- rbind(lmat, max(lmat) + 1)
    }
    added_height = added_height * globalScale
    lhei[body_iy] <<- lhei[body_iy] - added_height / nclust
    lhei <<- c(lhei, added_height)
  }
  if(!is.null(labels.above)){
    add_lmat_top(label.height)
  }
  if(!is.null(main)){
    add_lmat_top(main.height)
  }
  if(!is.null(labels.below)){
    add_lmat_bottom(label.height)
  }
  if(key){
    add_lmat_bottom(key.height, body_xpos = 1)
  }
  RowSideColors_size = 1
  if (!is.null(RowSideColors)) {
    if (!is.character(RowSideColors) || length(RowSideColors) != nrow(x)) 
      stop("'RowSideColors' must be a character vector of length nrow(x)")
    add_lmat_right(RowSideColors_size, solid = T)
  }
  if(doSidePlot | !is.null(extraData)){
    add_lmat_right(divLinesSize, solid = T)
  }
  if(doSidePlot){
    add_lmat_right(sidePlotSize, solid = F)
  }
  if(!is.null(extraData)){
    add_lmat_right(extraDataSize, solid = F)
  }
  
  if (!is.null(labels_rowsep)) {
    add_lmat_right(labRowSize, solid = T)
  }
  if(left_mai > 0){
    add_lmat_left(left_mai)
    lmat <- cbind(rep(0, nrow(lmat)), lmat)
    lwid[body_ix] = lwid[body_ix] - left_mai
    lwid <- c(left_mai, lwid)
    body_ix = body_ix + 1
  }
  
  
  #   print(lmat)
  #   print(lmat_custom)
  if(!is.null(lmat_custom)){
    if(class(lmat_custom) != 'matrix'){
      stop('class of lmat_custom must be matrix')
    }
    if(any(dim(lmat_custom) != dim(lmat))){
      stop(paste("lmat_custom does not match lmat. dim was", 
                 paste(dim(lmat_custom), collapse = ', '), 
                 "dim should be", 
                 paste(dim(lmat), collapse = ', ')))
    }
    lmat_custom = ifelse(lmat_custom > 0, lmat_custom + max(lmat), 0)
    lmat = lmat + lmat_custom
  }
  #   print(lmat)
  #   print(lhei)
  #   print(lwid)
  nf = layout(lmat, heights = lhei, widths = lwid)
  
  #layout.show(nf)
  
  breaks <- length(col) + 1
  symbreaks = any(x < 0, na.rm = TRUE)
  symkey = F#any(x < 0, na.rm = TRUE) || symbreaks
  
  if (length(breaks) == 1) {
    if (!symbreaks) 
      breaks <- seq(min(x, na.rm = T), max(x, na.rm = T), 
                    length = breaks)
    else {
      extreme <- max(abs(x), na.rm = TRUE)
      breaks <- seq(-extreme, extreme, length = breaks)
    }
  }
  nbr <- length(breaks)
  ncol <- length(breaks) - 1
  if (class(col) == "function") 
    col <- col(ncol)
  min.breaks <- min(breaks)
  max.breaks <- max(breaks)
  x[x < min.breaks] <- min.breaks
  x[x > max.breaks] <- max.breaks
  
  x = x[nrow(x):1,]#reverse so heatmap top row is top row of input matrix
  colorClasses = unique(RowSideColors)
  RowSideColors = rev(RowSideColors)
  nr <- nrow(x)
  nc <- ncol(x)
  
  
  insert_col = function(n, i){#insert n column following index i
    x <<- cbind(x[,1:i], matrix(NA, nrow = nrow(x), ncol = n),  x[,(i+1):ncol(x)])
    affected = colsep.minor > i
    colsep.minor[affected] <<- colsep.minor[affected] + n
    affected = colsep.major > i
    colsep.major[affected] <<- colsep.major[affected] + n
    if(!is.null(labels.below)) labels.below <<- c(labels.below[1:i], rep('', n), labels.below[(i + 1):length(labels.below)])
    if(!is.null(labels.above)) labels.above <<- c(labels.above[1:i], rep('', n), labels.above[(i + 1):length(labels.above)])
  }
  
  insert_row = function(n, i){#insert n rows following index i
    x <<- rbind(x[1:i,], matrix(NA, ncol = ncol(x), nrow = n),  x[(i+1):nrow(x),])
    affected = rowsep.minor >= i
    rowsep.minor[affected] <<- rowsep.minor[affected] + n
    affected = rowsep.major >= i
    rowsep.major[affected] <<- rowsep.major[affected] + n
    RowSideColors <<- c(RowSideColors[1:i], rep(NA, n), RowSideColors[(i+1):length(RowSideColors)])
  }
  
  do_col_major = colsep.major.original[1] > -1#test before value are changed
  if (colsep.minor.original[1] > -1){ #insert  columns in x and color
    for (i in 1:length(colsep.minor)){ 
      csep = colsep.minor[i]
      insert_col(round(sepwidth.minor*nc), csep)
    }
  }
  if (do_col_major){ #insert  columns in x and color
    for (i in 1:length(colsep.major)){ 
      csep = colsep.major[i]
      insert_col(round(sepwidth.major*nc), csep)
    }
  }
  x = x[nrow(x):1,]#reverse rows temporarily
  RowSideColors = rev(RowSideColors)
  do_row_major = rowsep.major[1] > -1#test before values are changed
  if (rowsep.minor.original[1] > -1){
    for (i in 1:length(rowsep.minor)){ 
      rsep = rowsep.minor[i]
      insert_row(round(sepwidth.minor*nr), rsep)
    }
  }
  if (do_row_major){
    for (i in 1:length(rowsep.major)){ 
      rsep = rowsep.major[i]
      insert_row(round(sepwidth.major*nr), rsep)
    }
  }
  
  #unreverse rows
  x = x[nrow(x):1,]
  RowSideColors = rev(RowSideColors)
  nc = ncol(x)
  nr = nrow(x)
  
  #plots the body of the heatap
  if(forPDF) startRasterMode(width = lwid[body_ix], height = lhei[body_iy])
  if(!forPDF) par(mai = rep(0,4))
  image(0:nc, 0:nr, t(x), xlim = c(0, nc), ylim = c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, breaks = breaks)
  if(forPDF) stopRasterMode(mai = rep(0,4))
  
  #plots column labels above column
  srtCol = 0
  adjCol = c(.5,1)
  colCol = NULL
  
  apply_col_labels = function(col_labels, yshift = .5, yadj = .5){
    plot0(ncol(x))
    xpd.orig <- par("xpd")
    par(xpd = NA)
    xpos <- 1:nc-.5
    n = sum(col_labels != '')
    n_scale = 1 #/ n  #no longer necessary to resize text according to columns, plot gets wider
    
    text(x = xpos, y = yshift, labels = col_labels, adj = c(.5, yadj),
         cex = cexCol * globalScale * n_scale, srt = srtCol, col = colCol)
    par(xpd = xpd.orig)
  }
  if(!is.null(labels.above)){
    apply_col_labels(labels.above, .1, 0)
  }
  
  #draw title
  if(!is.null(main)){
    plot.new()
    par(xpd = NA)
    text(0,1, main, cex = cex.main * globalScale, adj = c(0,1))
  }
  
  #draws column labels below column
  if(!is.null(labels.below)){
    apply_col_labels(labels.below, .9, 1)
  }
  
  
  if (key) {
    mai <- c(.5, 0,.5,0) * globalScale
    par(mai = mai, cex = 0.75 * globalScale, mgp = c(2, 1, 0))
    if (length(key.par) > 0) 
      do.call(par, key.par)
    tmpbreaks <- breaks
    if (F) {
      max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
      min.raw <- -max.raw
      tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)-.001
      tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)+.001
    }
    else {
      min.raw <- min(x, na.rm = T)
      max.raw <- max(x, na.rm = T)
    }
    z <- seq(min.raw, max.raw, by = min(diff(breaks)/4))
    
    if(forPDF) startRasterMode(height = lhei[body_iy])
    #draw the key color gradient
    if(!forPDF) par(mai = mai)
    image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks, 
          xaxt = "n", yaxt = "n")
    if(forPDF) stopRasterMode(mai = mai)
    if (is.null(key.xtickfun)) {
      lv <- -100:100 * 2
      keep = lv > min.raw & lv < max.raw
      lv = lv[keep]
      xmin = par('usr')[1]
      xmax = par('usr')[2]
      lvmin = min(lv)
      lvmax = max(lv)
      xv <- (lv - min.raw) / (max.raw - min.raw) * (xmax - xmin) + xmin
      xargs <- list(at = xv, labels = lv)
    }
    else {
      xargs <- key.xtickfun()
    }
    xargs$side <- 1
    par(xpd = NA)
    do.call(axis, xargs)
    #     MAX = max(abs(c(min(x), max(x))))
    #     par(usr = c(-MAX, MAX, 0, 1))
    #     axis(side = 1, )
    if (!is.na(key.lab)) {
      mtext(side = 3, key.lab, line = 2, padj = 1, 
            cex = key.cex * globalScale)
    }
  }
  if(!is.null(RowSideColors)){#plot colored block for each cluster
    tmp = as.factor(RowSideColors)
    vals = matrix(as.numeric(tmp), ncol = 1)
    
    lev = levels(tmp)
    mai = c(0,.1,0,.1) * globalScale
    if(forPDF) startRasterMode()
    #par(mai = c(0,.1,0,.1))
    if(!forPDF) par(mai = mai)
    image(0:1, 0:nr, t(vals), xlim = c(0, 1), ylim = c(0, nr), axes = FALSE, xlab = "", ylab = "", col = lev)
    if(forPDF) stopRasterMode(mai = mai)
    cluster_levels = unique(vals)
    cluster_levels = cluster_levels[!is.na(cluster_levels)]
    cluster_ids = vals[,1]
    names(cluster_ids) = 1:length(cluster_ids)
    cluster_ids = cluster_ids[!is.na(cluster_ids)]
    par(xpd = NA)
    cluster_starts = 1:nclust
    cluster_ends = 1:nclust
    #print(par('usr'))
    for(i in 1:length(RowSideLabels)){
      keep = (cluster_ids == cluster_levels[i])
      cluster_starts[i] = min(as.numeric(names(keep[keep])))
      cluster_ends[i] = max(as.numeric(names(keep[keep])))
      center = mean(as.numeric(names(keep[keep])))
      if(forPDF) center = center / nr #stupid bandaid correction for start/stopRaster model
      rowLab = rev(RowSideLabels)[i]
      text(.5,center-.5, rowLab, adj = c(.5,.5), cex = cexRow * globalScale)
      
    }
    #par(xpd = NA)
    
  }
  if(doSidePlot | !is.null(extraData)){
    par(mai=rep(0,4))
    plot0()
    #     (x=c(0,1),frame.plot=FALSE, y=c(0,1), xaxt='n',yaxt='n', type="n", xlab="",
    #          ylab="",xlim=c(0,1),ylim=c(0,1))
    #no idea how these 'actual' margins are selected
    min=0#-.04
    max=1#1.04
    rng=max-min
    
    ### draw the dotted lines connecting line plots to heatmap
    #line at top and bottom
    lines(c(0,1),c(0,0),lty=2)
    lines(c(0,1),c(1,1),lty=2)
    for(i in 1:(nclust-1)){
      #calculate number of rows covered so far as fraction of total
      hFrac_a = (cluster_ends[i]) / (max(cluster_ends))
      hFrac_b = (cluster_starts[i+1]-1) / (max(cluster_ends))
      hFrac_mean = mean(c(hFrac_a, hFrac_b))
      lplotFrac = i/nclust
      #x is always from min to max
      #y goes from variable fraction on heatmap to constant fraction of plotting column
      meet_in_middle = 0
      if(meet_in_middle > 0){
        mim = 1 - meet_in_middle
        lines(c(min ,max - mim*rng),c(rng*hFrac_a,rng*hFrac_mean),lty=2)
        lines(c(min,max - mim*rng),c(rng*hFrac_b,rng*hFrac_mean),lty=2)
      }
      
      lines(c(min + meet_in_middle * rng, max), c(rng*hFrac_mean,rng*lplotFrac),lty=2)
    }
  }
  if(doSidePlot){
    avgA = matrix(0,nclust, ncol(clust$centers))
    for(i in 1:nrow(avgA)){
      avgA[i,] = clust$centers[i,]
    }
    
    par(mar = c(0,0, 0,0.5))
    
    #draw line chart represented each class from clustering
    nsplits = 1
    if(colsep.minor[1] > -1){
      nsplits = length(colsep.minor)+1
    }
    if(colsep.minor)
      win = ncol(avgA) / nsplits
    #     print(avgA)
    colorClasses = RColorBrewer::brewer.pal(max(nsplits,3), 'Dark2')
    colorClasses = colorClasses[1:nsplits]#stupid colorbrewer min classes warning workaround
    for (i in 1:nclust) { 
      xrange <- as.numeric(range(1:(ncol(avgA)/nsplits)))
      yrange <- range(avgA)
      xspace = (xrange[2] - xrange[1]) * .15
      yspace = (yrange[2] - yrange[1]) * .15
      # set up the plot 
      #op <- par(mar = rep(.01, 4))
      plot(xrange, yrange, xaxt='n',yaxt='n', type="n", xlab="",
           ylab="", ylim = c(min(yrange)-yspace, max(yrange)+yspace), xlim = c(min(xrange)-xspace, max(xrange)+xspace))  #c(minVal - .1*abs(minVal - maxVal),maxVal + .1*abs(minVal - maxVal))) 
      colors <- colorClasses[i] 
      linetype <- c(1:nclust) 
      plotchar <- seq(18,18+nclust,1)
      
      #   axis(side=1,tick=TRUE,at=days)
      vals <- avgA[i,]
      
      for(s in 1:nsplits){
        start = (s - 1) * win + 1
        end = s * win
        xs = 1:(ncol(avgA)/nsplits)#first half of profile
        
        line_dat = vals[start:end]
        if(sidePlot_smoothing > 1){
          line_dat = movingAverage(line_dat, n = sidePlot_smoothing)
        }
        lines(xs, line_dat, type="l", lwd=2.5,
              lty=1, col=colorClasses[s], pch=16) 
      }
      lines(rep(mean(xrange),2), yrange, lty = 3)
    }
  }
  if(!is.null(extraData)){
    par(mai = rep(0, 4), xpd = T)
    YMIN = 0
    YMAX = 20
    colorClasses = extraData_colors
    for(i in 1:nclust){
      #       plot0()
      #       box()
      subset = rownames(clust$data)[clust$cluster == i]
      dat = extraData[subset,, drop = F]# + 64
      colnames(dat) = NULL
      # extraData_asLines = T
      # extraData_asHmap = T
      if(is.function(extraData_plotFunction)){
        extraData_plotFunction(dat, c(YMIN,YMAX))
#         par(mai = rep(.1,4))
#         cr = colorRamp(c('green', 'red'))
#         colors = rgb(cr(0:100/100)/255)
#         sub_hclust = hclust(dist(dat))
#         dat = dat[sub_hclust$order,]
#         image(0:3, 0:nrow(dat), t(dat), xlim = c(0, 3), ylim = c(0, nrow(dat)), axes = FALSE, xlab = "", ylab = "", col = colors)
#       }else if(extraData_asLines){
#         xs = 1:ncol(dat)
#         plot(0, type = 'n', axes = F, xlim = c(.8,ncol(dat) + .2), ylim = c(0, max(dat)))
#         hidden = apply(dat, 1, function(x){
#           lines(xs, x, col = 'gray')
#         })
#         
#         lines(xs, apply(dat, 2, median), col = "#f62626", lwd = 3)
      }else{
        M = colMeans(dat)
        n = nrow(dat)
        s = apply(dat, 2, sd)
        SE = s / sqrt(n)
        E = 0
        if(n > 1) E = qt(.975, df=n - 1) * SE
        low = M - E
        mid = M
        high = M + E
        #       print(rbind(low, mid, high))
        if(any(is.na(low))) low = mid
        if(any(is.na(high))) high = mid
        
        low = ifelse(low <= 0, 1, low)
        barplot2(mid, ci.l = low, ci.u = high, col = colorClasses, plot.ci = T, ylim = c(YMIN,YMAX), axes = F, space = 0, bty = 'o')
        for(pos in c(4,8,12,16)){
          lines(c(0,nclust+1), c(pos,pos), lty = 2) 
        }
        
        #axis(side = 2, at = c(250,500,1000,2000))
        box()
      }
    }
  }
  par(xpd = NA)
  if(!is.null(labels_rowsep)){
    par(mai = rep(0,4))
    plot0(height = nr)
    for(i in 1:length(labels_rowsep)){
      if(!doSidePlot && is.null(extraData)){
        rsep_prev = 1
        if(i > 1){
          rsep_prev = rowsep.major[i - 1] + sepwidth.major
        }
        rsep_curr = rowsep.major.original[i]
        if(i > 1){
          rsep_curr = rsep_prev - 1 + rowsep.major.original[i] - rowsep.major.original[i - 1]
        }
        if(i > length(rowsep.major.original)){
          rsep_curr = nrow(x)
        }
        rsep_mid = (mean(c(rsep_prev, rsep_curr)))
        ypos = nr - rsep_mid + 1
      }else{
        ypos = nr - (i - .5) / nclust * nr
      }
      
      text(.5, ypos, labels_rowsep[i], adj = c(.5,.5), cex = cexRow * globalScale )
    }
  }
}

do_example = F
if(do_example){
  #   cl = 'MCF10A'
  #   hm = c('H3K4AC', 'H3K4ME3')#c('H3K27AC','H3K27ME3','H3K4AC', 'H3K4ME3', 'H4K20ME3')
  #   profiles_to_plot = paste(cl, hm, sep = '_')
  #   ngs_profiles = list()
  #   nr = 50
  #   nc = 49
  #   for(tp in profiles_to_plot){
  #     dat = matrix(runif(nr * nc, -3, 7), nrow = nr, ncol = nc)
  #     rownames(dat) = 1:nrow(dat)
  #     reduce = runif(nr/2, 1, nr)
  #     dat[reduce,] = dat[reduce,] - 4
  #     dat = ifelse(dat < -3, -3, dat)
  #     ngs_profiles[[tp]] = dat
  #   }
  
  
  
  load('norm_matched.save')
  load('promoter_wide_matched_prof.save')
  prof_normed = lapply(promoter_wide_matched_prof,  function(x){
    MIN = quantile(x, .1)
    MAX = quantile(x, .98)
    x = ifelse(x < MIN, MIN, x)
    x = ifelse(x > MAX, MAX, x)
    u = mean(x)
    sdev = sd(x)
    x = (x - u) / sdev
    return(x)
  })
  histone_mods = unique(sapply(strsplit(names(prof_normed), '_'), function(x)return(x[2])))
  cell_lines = c('MCF10A', 'MCF7', 'MDA231')
  pdf('compound_plots.pdf', width = 10, height = 8)
  nr = 11
  lmat_custom = matrix(0, ncol = 7, nrow = 11)
  lmat_custom[nr-1,5] = 1
  lmat_custom[nr-1,6] = 2
  lmat_custom[nr,5:6] = 3
  lmat_custom[1,7] = 4
  for(i in 0:5){
    main_title = paste(histone_mods[i+1], 'within 10kb of TSS')
    res = heatmap.ngsplots(ngs_profiles = prof_normed[c(1,7,13)+i], label_clusters = T, extraData = rna_norm_matched, labels_below = cell_lines, labels_above = main_title,labels_right = '', sortClustersByTotal = T, labelWithCounts = T, cex.col = 1.4, nclust = 8, lmat_custom = lmat_custom)
    plot0();text(.5,.5, 'average profile')
    plot0();text(.5,.5, 'log gene expression')
    plot0();legend('center', legend = c('MCF10A', 'MCF7', 'MDA231'), fill = RColorBrewer::brewer.pal(3, 'Set1'), horiz = T, bty = 'n')
    plot0();text(.5,.5, 'cluster size')
    
  }
  dev.off()
}