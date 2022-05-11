make_embedding_plot <- function (emb, cell.colors = NULL,
                                 do.par = F,  cell.border.alpha = 0.3) {
  
  if (do.par) {
    par(mfrow = c(1, 1), 
        mar = c(3.5, 3.5, 2.5, 1.5), 
        mgp = c(2,0.65, 0), 
        cex = 0.85)
  }
  
  celcol <- cell.colors[rownames(emb)]
  
  plot(emb, 
       bg = celcol, pch = 21, col = ac(1, alpha = cell.border.alpha))
  
  return(NULL)
}


my.add.arrows <- function (emb, vel, n = 100, cell.colors = NULL, 
                           arrow.color='black',
                           corr.sigma = 0.05, 
                           show.grid.flow = FALSE, grid.n = 20, grid.sd = NULL, min.grid.cell.mass = 1, 
                           min.arrow.size = NULL, arrow.scale = 1, max.grid.arrow.length = NULL, 
                           fixed.arrow.length = FALSE, plot.grid.points = FALSE, scale = "log", 
                           nPcs = NA, arrow.lwd = 1, xlab = "", ylab = "", n.cores = defaultNCores(), 
                           do.par = T, show.cell = NULL, cell.border.alpha = 0.3, cc = NULL, 
                           return.details = FALSE, expression.scaling = FALSE, ...) {
  
  if (do.par) 
    par(mfrow = c(1, 1), mar = c(3.5, 3.5, 2.5, 1.5), mgp = c(2, 
                                                              0.65, 0), cex = 0.85)
  
  # plot(emb, bg = celcol, pch = 21, col = ac(1, alpha = cell.border.alpha), 
  #      xlab = xlab, ylab = ylab, ...)
  em <- as.matrix(vel$current)
  ccells <- intersect(rownames(emb), colnames(em))
  em <- em[, ccells]
  emb <- emb[ccells, ]
  nd <- as.matrix(vel$deltaE[, ccells])
  cgenes <- intersect(rownames(em), rownames(nd))
  nd <- nd[cgenes, ]
  em <- em[cgenes, ]
  
  if (is.null(cc)) {
    cat("delta projections ... ")
    
    cat("sqrt ")
    cc <- colDeltaCorSqrt(em, (sqrt(abs(nd)) * sign(nd)), 
                          nthreads = n.cores)
    
    colnames(cc) <- rownames(cc) <- colnames(em)
    diag(cc) <- 0
  }
  cat("knn ... ")
  if (n > nrow(cc)) {
    n <- nrow(cc)
  }
  emb.knn <- balancedKNN(t(emb), k = n, maxl = nrow(emb), dist = "euclidean", 
                         n.threads = n.cores)
  diag(emb.knn) <- 1
  cat("transition probs ... ")
  tp <- exp(cc/corr.sigma) * emb.knn
  tp <- t(t(tp)/Matrix::colSums(tp))
  tp <- as(tp, "dgCMatrix")
  cat("done\n")
  
  cat("calculating arrows ... ")
  arsd <- data.frame(t(embArrows(emb, tp, arrow.scale, 
                                 n.cores)))
  rownames(arsd) <- rownames(emb)
  # if (expression.scaling) {
  #   tpb <- tp > 0
  #   tpb <- t(t(tpb)/colSums(tpb))
  #   es <- as.matrix(em %*% tp) - as.matrix(em %*% as.matrix(tpb))
  #   pl <- pmin(1, pmax(0, apply(as.matrix(vel$deltaE[, 
  #                                                    colnames(es)]) * es, 2, sum)/sqrt(colSums(es * 
  #                                                                                                es))))
  #   arsd <- arsd * pl
  # }
  ars <- data.frame(cbind(emb, emb + arsd))
  colnames(ars) <- c("x0", "y0", "x1", "y1")
  colnames(arsd) <- c("xd", "yd")
  rownames(ars) <- rownames(emb)
  cat("done\n")
  
  cat("grid estimates ... ")
  rx <- range(c(range(ars$x0), range(ars$x1)))
  ry <- range(c(range(ars$y0), range(ars$y1)))
  gx <- seq(rx[1], rx[2], length.out = grid.n)
  gy <- seq(ry[1], ry[2], length.out = grid.n)
  if (is.null(grid.sd)) {
    grid.sd <- sqrt((gx[2] - gx[1])^2 + (gy[2] - 
                                           gy[1])^2)/2
    cat("grid.sd=", grid.sd, " ")
  }
  if (is.null(min.arrow.size)) {
    min.arrow.size <- sqrt((gx[2] - gx[1])^2 + (gy[2] - 
                                                  gy[1])^2) * 0.01
    cat("min.arrow.size=", min.arrow.size, " ")
  }
  if (is.null(max.grid.arrow.length)) {
    max.grid.arrow.length <- sqrt(sum((par("pin")/c(length(gx), 
                                                    length(gy)))^2)) * 0.25
    cat("max.grid.arrow.length=", max.grid.arrow.length, 
        " ")
  }
  garrows <- do.call(rbind, lapply(gx, function(x) {
    cd <- sqrt(outer(emb[, 2], -gy, "+")^2 + (x - 
                                                emb[, 1])^2)
    cw <- dnorm(cd, sd = grid.sd)
    gw <- Matrix::colSums(cw)
    cws <- pmax(1, Matrix::colSums(cw))
    gxd <- Matrix::colSums(cw * arsd$xd)/cws
    gyd <- Matrix::colSums(cw * arsd$yd)/cws
    al <- sqrt(gxd^2 + gyd^2)
    vg <- gw >= min.grid.cell.mass & al >= min.arrow.size
    cbind(rep(x, sum(vg)), gy[vg], x + gxd[vg], gy[vg] + 
            gyd[vg])
  }))
  colnames(garrows) <- c("x0", "y0", "x1", "y1")
  
  
  alen <- pmin(max.grid.arrow.length, sqrt(((garrows[, 
                                                     3] - garrows[, 1]) * par("pin")[1]/diff(par("usr")[c(1, 
                                                                                                          2)]))^2 + ((garrows[, 4] - garrows[, 2]) * 
                                                                                                                       par("pin")[2]/diff(par("usr")[c(3, 4)]))^2))
  # if had made a plot, add the arrows to the plot
  suppressWarnings(lapply(1:nrow(garrows), 
                          function(i) arrows(garrows[i,1], garrows[i, 2], garrows[i, 3], garrows[i,4], 
                                             length = alen[i], 
                                             lwd = arrow.lwd,
                                             col=arrow.color)
  )
  )
  
  # return(list('garrows'=garrows,
  #             'alen'=alen,
  #             'arrow.lwd'=arrow.lwd))
  return()
  
}


colDeltaCorSqrt <- function(e, d, nthreads = 1L) {
  .Call('_velocyto_R_colDeltaCorSqrt', PACKAGE = 'velocyto.R', e, d, nthreads)
}


balancedKNN <- function(val,k,maxl=k,return.distance.values=FALSE,n.threads=1,dist='cor') {
  if(class(dist)=="dist") { # actual distance was passed
    if(!all(labels(dist)==colnames(val))) { stop("balancedKNN(): supplied distance doesn't match the columns of val") }
    cd <- as.matrix(dist);
  }  else {
    if(dist=='cor') {
      cd <- 1-cor(val);
    } else if(dist=='euclidean') {
      cd <- as.matrix(dist(t(val)))
    } else {
      stop(paste("unknown distance",dist,"specified"))
    }
  }
  z <-  balanced_knn(cd,k,maxl,return.distance.values,n.threads);
  rownames(z) <- colnames(z) <- colnames(val);
  z
}


balanced_knn <- function(d, k, maxl, returnDistanceValues = FALSE, nthreads = 1L) {
  .Call('_velocyto_R_balanced_knn', PACKAGE = 'velocyto.R', d, k, maxl, returnDistanceValues, nthreads)
}


embArrows <- function(emb, tp, arrowScale = 1.0, nthreads = 1L) {
  .Call('_velocyto_R_embArrows', PACKAGE = 'velocyto.R', emb, tp, arrowScale, nthreads)
}


color.vector <- function(seu, cl_annot){
  ## Obtain vector for cell coloring
  suppressPackageStartupMessages(require(scales, lib.loc = "/usr/local/lib/R/site-library/scales"))
  
  labels <- seu@meta.data[[cl_annot]]
  ncol <- length(unique(labels))
  x <- as.factor(labels)
  levels(x) <- 1:length(levels(x))
  x <- as.numeric(x)
  colvec <- hue_pal(l=80)(ncol)[x]
  
  names(colvec) <- colnames(seu)
  return(colvec)
}

