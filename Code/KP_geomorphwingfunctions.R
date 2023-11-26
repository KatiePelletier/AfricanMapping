
wing.spec <- gridPar(
  pt.bg = "black",
  pt.size = 1,
  link.col = "black",
  link.lwd = 2,
  link.lty = 1,
  out.col = "black",
  out.cex = 0.1,
  tar.pt.bg = "red",
  tar.pt.size = 1,
  tar.link.col = "red",
  tar.link.lwd = 2,
  tar.link.lty = 1,
  tar.out.col = "red",
  tar.out.cex = 0.1,
  n.col.cell = 20,
  grid.col = "black",
  grid.lwd = 1,
  grid.lty = 1,
  txt.adj = NULL,
  txt.pos = NULL,
  txt.cex = 0.8,
  txt.col = "black"
)


#Used for drawing wings in geomorph (when posterior x vein is missing)
cvl.links <- c(11, 12, 
               12, 13, 
               13, 14, 
               14, 15, 
               15, 16, 
               16, 1, 
               1, 17, 
               17, 18, 
               18, 19, 
               19, 2,
               2, 20, 
               20, 3, 
               3, 21, 
               21, 22, 
               22, 4, 
               4, 23, 
               23, 24, 
               24, 25, 
               25, 26, 
               26, 5, 
               5, 27, 
               27, 28, 
               5, 47, 
               47, 6, 
               4, 43, 
               43, 44, 
               44, 45, 
               45, 46, 
               46, 9, 
               9, 28, 
               8, 7, 
               1, 29, 
               29, 30, 
               30, 31, 
               31, 32, 
               32, 10,
               9, 42, 
               42, 8, 
               8, 41, 
               41, 40, 
               40, 39, 
               39, 38, 
               38, 3, 
               37, 7, 
               7, 36, 
               36, 35, 
               35, 34, 
               34, 33, 
               33, 2)
cvl.links <- matrix(cvl.links, ncol = 2, byrow = TRUE)


wing.links <- c(13, 14, 
                14, 15, 
                15, 16,
                16, 17,
                17, 18,
                18, 1, 
                1, 19,
                19, 20, 
                20, 21, 
                21, 2,
                2, 22, 
                22, 3, 
                3, 23, 
                23, 24, 
                24, 4, 
                4, 25, 
                25, 26, 
                26, 27, 
                27, 28, 
                28, 5, 
                5, 29, 
                29, 30, 
                5, 47, 
                47, 6, 
                4, 43, 
                43, 44, 
                44, 45, 
                45, 46, 
                46, 11, 
                11, 42,
                42, 10, 
                10, 9, 
                10, 41, 
                41, 40, 
                40, 39, 
                39, 38, 
                38, 3,
                37, 9, 
                9, 36, 
                36, 8, 
                8, 35, 
                35, 34, 
                34, 2, 
                8, 48, 
                48, 7, 
                12, 33, 
                33, 32, 
                32, 7, 
                7, 31, 
                31, 1)
wing.links <- matrix(wing.links, ncol = 2, byrow = TRUE)

####15 point shape links 

fifteen.links <- c(1, 7, 
                  7, 12, 
                  12, 13, 
                  13, 14,
                  14, 15,
                  15, 11, 
                  11, 10, 
                  10, 9,
                  9, 8, 
                  8, 6, 
                  6, 12, 
                  6, 2, 
                  8, 13, 
                  10, 14, 
                  9, 3, 
                  4, 5,
                  5, 3, 
                  5, 11)
fifteen.links <- matrix(fifteen.links, ncol = 2, byrow = TRUE)


# Maria's Altered Geomorph Functions

# plotAllSpecimens without the axes:
MP_altered_plotAllSpecimens <- function (A, mean = TRUE, links = NULL, label = FALSE, plot.param = list()) 
{
  if (length(dim(A)) != 3) {
    stop("Data matrix not a 3D array (see 'arrayspecs').")
  }
  if (any(is.na(A)) == T) {
    stop("Data matrix contains missing values. Estimate these first (see 'estimate.missing').")
  }
  k <- dim(A)[2]
  if (mean == TRUE) {
    mn <- mshape(A)
  }
  p.p <- plot.param
  if (is.null(p.p$pt.bg)) 
    p.p$pt.bg = "gray"
  if (is.null(p.p$pt.cex)) 
    p.p$pt.cex = 1
  if (is.null(p.p$mean.bg)) 
    p.p$mean.bg = "black"
  if (is.null(p.p$mean.cex)) 
    p.p$mean.cex = 2
  if (is.null(p.p$link.col)) 
    p.p$link.col = "black"
  if (is.null(p.p$link.lwd)) 
    p.p$link.lwd = 2
  if (is.null(p.p$link.lty)) 
    p.p$link.lty = 1
  if (is.null(p.p$txt.adj)) 
    p.p$txt.adj = c(-0.1, -0.1)
  if (is.null(p.p$txt.col)) 
    p.p$txt.col = "black"
  if (is.null(p.p$txt.cex)) 
    p.p$txt.cex = 0.8
  if (is.null(p.p$txt.pos)) 
    p.p$txt.pos = 1
  if (k == 2) {
    plot(A[, 1, ], A[, 2, ], asp = 1, pch = 21, bg = p.p$pt.bg, 
         cex = p.p$pt.cex * 1, axes = FALSE, xaxt='n', yaxt='n', ann=FALSE )
    if (mean == TRUE) {
      if (is.null(links) == FALSE) {
        linkcol <- rep(p.p$link.col, nrow(links))[1:nrow(links)]
        linklwd <- rep(p.p$link.lwd, nrow(links))[1:nrow(links)]
        linklty <- rep(p.p$link.lty, nrow(links))[1:nrow(links)]
        for (i in 1:nrow(links)) {
          segments(mn[links[i, 1], 1], mn[links[i, 1], 
                                          2], mn[links[i, 2], 1], mn[links[i, 2], 2], 
                   col = linkcol[i], lty = linklty[i], lwd = linklwd[i])
        }
      }
      points(mn, pch = 21, bg = p.p$mean.bg, cex = p.p$mean.cex)
      if (label == TRUE) {
        text(mn, label = paste(1:dim(mn)[1]), adj = (p.p$txt.adj + 
                                                       p.p$mean.cex), pos = p.p$txt.pos, cex = p.p$txt.cex, 
             col = p.p$txt.col)
      }
    }
  }
  if (k == 3) {
    A3d <- NULL
    for (i in 1:dim(A)[[3]]) {
      A3d <- rbind(A3d, A[, , i])
    }
    plot3d(A3d, type = "s", col = p.p$pt.bg, size = p.p$pt.cex * 1.5, 
           aspect = FALSE, axes = FALSE, xaxt='n', yaxt='n', ann=FALSE)
    if (mean == TRUE) {
      if (is.null(links) == FALSE) {
        linkcol <- rep(p.p$link.col, nrow(links))[1:nrow(links)]
        linklwd <- rep(p.p$link.lwd, nrow(links))[1:nrow(links)]
        linklty <- rep(p.p$link.lty, nrow(links))[1:nrow(links)]
        for (i in 1:nrow(links)) {
          segments3d(rbind(mn[links[i, 1], ], mn[links[i, 
                                                       2], ]), col = linkcol[i], lty = linklty[i], 
                     lwd = linklwd[i])
        }
      }
      points3d(mn, color = p.p$mean.bg, size = p.p$mean.cex * 
                 2)
      if (label == TRUE) {
        text3d(mn, texts = paste(1:dim(mn)[1]), adj = (p.p$txt.adj + 
                                                         p.p$mean.cex), pos = p.p$txt.pos, cex = p.p$txt.cex, 
               col = p.p$txt.col)
      }
    }
  }
}

# Changing plotRedToTarget to match my colour scheme

MP_altered_plotRefToTarget <- function (M1, M2, mesh = NULL, outline = NULL, method = c("TPS", 
                                                                                        "vector", "points", "surface"), mag = 1, links = NULL, label = FALSE, 
                                        axes = FALSE, gridPars = NULL, useRefPts = FALSE, ...) 
{
  method <- match.arg(method)
  if (any(is.na(M1))) {
    stop("Data contains missing values. Estimate these first (see 'estimate.missing').")
  }
  if (any(is.na(M2))) {
    stop("Data contains missing values. Estimate these first (see 'estimate.missing').")
  }
  if (is.null(gridPars)) 
    gP <- gridPar()
  else gP <- gridPars
  k <- dim(M1)[2]
  mag <- (mag - 1)
  M2 <- M2 + (M2 - M1) * mag
  limits <- function(x, s) {
    r <- range(x)
    rc <- scale(r, scale = FALSE)
    l <- mean(r) + s * rc
  }
  if (k == 2) {
    if (method == "TPS") {
      tps(M1, M2, gP$n.col.cell, sz = gP$tar.pt.size, pt.bg = gP$tar.pt.bg, 
          grid.col = gP$grid.col, grid.lwd = gP$grid.lwd, 
          grid.lty = gP$grid.lty, refpts = useRefPts)
      if (is.null(links) == FALSE) {
        linkcol <- rep(gP$tar.link.col, nrow(links))[1:nrow(links)]
        linklwd <- rep(gP$tar.link.lwd, nrow(links))[1:nrow(links)]
        linklty <- rep(gP$tar.link.lty, nrow(links))[1:nrow(links)]
        for (i in 1:nrow(links)) {
          segments(M2[links[i, 1], 1], M2[links[i, 1], 
                                          2], M2[links[i, 2], 1], M2[links[i, 2], 2], 
                   col = linkcol[i], lty = linklty[i], lwd = linklwd[i])
        }
      }
      if (label) {
        text(M2, label = paste(1:dim(M2)[1]), adj = gP$txt.adj, 
             pos = gP$txt.pos, cex = gP$txt.cex, col = gP$txt.col)
      }
      if (!is.null(outline)) {
        curve.warp <- xy.coords(tps2d(outline, M1, M2))
        plot.xy(curve.warp, type = "p", pch = 19, cex = gP$tar.out.cex, 
                col = gP$tar.out.col)
      }
      if (!useRefPts) {
        plot.xy(xy.coords(M2), type = "p", pch = 21, 
                cex = gP$tar.pt.size, bg = gP$tar.pt.bg)
      }
      else {
        plot.xy(xy.coords(M1), type = "p", pch = 21, 
                cex = gP$pt.size, bg = gP$pt.bg)
      }
    }
    if (method == "vector") {
      plot.new()
      if (axes) {
        plot.window(limits(M1[, 1], 1.25), limits(M1[, 
                                                     2], 1.25), xlab = "x", ylab = "y", asp = 1)
      }
      if (!axes) {
        plot.window(limits(M1[, 1], 1.25), limits(M1[, 
                                                     2], 1.25), xlab = "", ylab = "", asp = 1, xaxt = "n", 
                    yaxt = "n")
      }
      if (!is.null(links)) {
        linkcol <- rep(gP$link.col, nrow(links))[1:nrow(links)]
        linklwd <- rep(gP$link.lwd, nrow(links))[1:nrow(links)]
        linklty <- rep(gP$link.lty, nrow(links))[1:nrow(links)]
        for (i in 1:nrow(links)) {
          segments(M2[links[i, 1], 1], M2[links[i, 1], 
                                          2], M2[links[i, 2], 1], M2[links[i, 2], 2], 
                   col = linkcol[i], lty = linklty[i], lwd = linklwd[i])
        }
      }
      if (label) {
        text(M1, label = paste(1:dim(M1)[1]), adj = gP$txt.adj, 
             pos = gP$txt.pos, cex = gP$txt.cex, col = gP$txt.col)
      }
      arrows(M1[, 1], M1[, 2], M2[, 1], M2[, 2], length = 0.075, 
             lwd = 2)
      plot.xy(xy.coords(M1), type = "p", pch = 21, bg = gP$pt.bg, 
              cex = gP$pt.size)
    }
    if (method == "points") {
      plot.new()
      if (axes) {
        plot.window(limits(M1[, 1], 1.25), limits(M1[, 
                                                     2], 1.25), xlab = "x", ylab = "y", asp = 1)
      }
      if (!axes) {
        plot.window(limits(M1[, 1], 1.25), limits(M1[, 
                                                     2], 1.25), xlab = "", ylab = "", asp = 1, xaxt = "n", 
                    yaxt = "n")
      }
      if (label) {
        text(M1, label = paste(1:dim(M1)[1]), adj = gP$txt.adj, 
             pos = gP$txt.pos, cex = gP$txt.cex, col = gP$txt.col)
      }
      if (!is.null(outline)) {
        curve.warp <- tps2d(outline, M1, M2)
        plot.xy(xy.coords(outline), type = "p", pch = 19, 
                cex = gP$out.cex, col = gP$out.col)
        plot.xy(xy.coords(curve.warp), type = "p", pch = 19, 
                cex = gP$tar.out.cex, col = gP$tar.out.col)
      }
      if (!is.null(links)) {
        linkcol <- rep("black", nrow(links))[1:nrow(links)] # changed colour to my colour scheme
        linklwd <- rep(3, nrow(links))[1:nrow(links)]
        linklty <- rep(gP$link.lty, nrow(links))[1:nrow(links)]
        tarlinkcol <- rep("#E69F00", nrow(links))[1:nrow(links)] # changed colour to my colour scheme
        tarlinklwd <- rep(3, nrow(links))[1:nrow(links)]
        tarlinklty <- rep(gP$tar.link.lty, nrow(links))[1:nrow(links)]
        for (i in 1:nrow(links)) {
          segments(M1[links[i, 1], 1], M1[links[i, 1], 
                                          2], M1[links[i, 2], 1], M1[links[i, 2], 2], 
                   col = linkcol[i], lty = linklty[i], lwd = linklwd[i])
          segments(M2[links[i, 1], 1], M2[links[i, 1], 
                                          2], M2[links[i, 2], 1], M2[links[i, 2], 2], 
                   col = tarlinkcol[i], lty = tarlinklty[i], 
                   lwd = tarlinklwd[i])
        }
      }
      plot.xy(xy.coords(M2), type = "p", pch = 16, col = "#E69F00", bg = gP$tar.pt.bg, 
              cex = 2)
      plot.xy(xy.coords(M1), type = "p", pch = 16, col = "black", bg = gP$pt.bg, 
              cex = 2)
    }
    if (method == "surface") {
      stop("Surface plotting for 3D landmarks only.")
    }
  }
  if (k == 3) {
    if (method == "TPS" && class(M2) == "predshape.k2") {
      tps(M1[, 1:2], M2[, 1:2], gP$n.col.cell, sz = gP$tar.pt.size, 
          pt.bg = gP$tar.pt.bg, grid.col = gP$grid.col, 
          grid.lwd = gP$grid.lwd, grid.lty = gP$grid.lty, 
          refpts = useRefPts, k3 = TRUE)
      if (is.null(links) == FALSE) {
        linkcol <- rep(gP$tar.link.col, nrow(links))[1:nrow(links)]
        linklwd <- rep(gP$tar.link.lwd, nrow(links))[1:nrow(links)]
        linklty <- rep(gP$tar.link.lty, nrow(links))[1:nrow(links)]
        for (i in 1:nrow(links)) {
          segments3d(rbind(M2[links[i, 1], ], M2[links[i, 
                                                       2], ]), col = linkcol[i], lty = linklwd[i], 
                     lwd = linklty[i])
        }
      }
      if (label == TRUE) {
        text3d(M2, texts = paste(1:dim(M2)[1]), adj = (gP$txt.adj + 
                                                         gP$pt.size), cex = gP$txt.cex, col = gP$txt.col)
      }
      if (!is.null(outline)) {
        curve.warp <- tps2d(outline, M1[, 1:2], M2[, 
                                                   1:2])
        points3d(cbind(curve.warp, 0), size = gP$tar.out.cex, 
                 col = gP$tar.out.col)
      }
    }
    if (method == "TPS" && class(M2) == "predshape.k3") {
      layout3d(matrix(c(1, 2), 1, 2))
      tps(M1[, 1:2], M2[, 1:2], gP$n.col.cell, sz = gP$tar.pt.size, 
          pt.bg = gP$tar.pt.bg, grid.col = gP$grid.col, 
          grid.lwd = gP$grid.lwd, grid.lty = gP$grid.lty, 
          refpts = useRefPts, k3 = TRUE)
      title3d("X,Y tps grid")
      if (is.null(links) == FALSE) {
        linkcol <- rep(gP$tar.link.col, nrow(links))[1:nrow(links)]
        linklwd <- rep(gP$tar.link.lwd, nrow(links))[1:nrow(links)]
        linklty <- rep(gP$tar.link.lty, nrow(links))[1:nrow(links)]
        for (i in 1:nrow(links)) {
          segments3d(rbind(M2[links[i, 1], ], M2[links[i, 
                                                       2], ]), col = linkcol[i], lty = linklwd[i], 
                     lwd = linklty[i])
        }
      }
      if (label == TRUE) {
        text3d(M2, texts = paste(1:dim(M2)[1]), adj = (gP$txt.adj + 
                                                         gP$pt.size), cex = gP$txt.cex, col = gP$txt.col)
      }
      if (!is.null(outline)) {
        curve.warp <- tps2d(outline, M1[, 1:2], M2[, 
                                                   1:2])
        points3d(cbind(curve.warp, 0), size = gP$tar.out.cex, 
                 col = gP$tar.out.col)
      }
      b <- c(1, 3)
      tps(M1[, b], M2[, b], gP$n.col.cell, sz = gP$tar.pt.size, 
          pt.bg = gP$tar.pt.bg, grid.col = gP$grid.col, 
          grid.lwd = gP$grid.lwd, grid.lty = gP$grid.lty, 
          refpts = useRefPts, k3 = TRUE)
      title3d("X,Y tps grid")
      if (is.null(links) == FALSE) {
        linkcol <- rep(gP$tar.link.col, nrow(links))[1:nrow(links)]
        linklwd <- rep(gP$tar.link.lwd, nrow(links))[1:nrow(links)]
        linklty <- rep(gP$tar.link.lty, nrow(links))[1:nrow(links)]
        for (i in 1:nrow(links)) {
          segments3d(rbind(M2[links[i, 1], c(1, 3, 2)], 
                           M2[links[i, 2], c(1, 3, 2)]), col = linkcol[i], 
                     lty = linklwd[i], lwd = linklty[i])
        }
      }
      if (label == TRUE) {
        text3d(M2[, c(1, 3, 2)], texts = paste(1:dim(M2)[1]), 
               adj = (gP$txt.adj + gP$pt.size), cex = gP$txt.cex, 
               col = gP$txt.col)
      }
      if (!is.null(outline)) {
        curve.warp <- tps2d(outline, M1[, b], M2[, b])
        points3d(cbind(curve.warp, 0), size = gP$tar.out.cex, 
                 col = gP$tar.out.col)
      }
    }
    if (method == "TPS" && any(class(M2) == "matrix")) {
      old.par <- par(no.readonly = TRUE)
      layout(matrix(c(1, 2), 1, 2))
      par(mar = c(1, 1, 1, 1))
      tps(M1[, 1:2], M2[, 1:2], gP$n.col.cell, sz = gP$tar.pt.size, 
          pt.bg = gP$tar.pt.bg, grid.col = gP$grid.col, 
          grid.lwd = gP$grid.lwd, grid.lty = gP$grid.lty, 
          refpts = useRefPts)
      if (is.null(links) == FALSE) {
        for (i in 1:nrow(links)) {
          linkcol <- rep(gP$tar.link.col, nrow(links))[1:nrow(links)]
          linklwd <- rep(gP$tar.link.lwd, nrow(links))[1:nrow(links)]
          linklty <- rep(gP$tar.link.lty, nrow(links))[1:nrow(links)]
          segments(M2[links[i, 1], 1], M2[links[i, 1], 
                                          2], M2[links[i, 2], 1], M2[links[i, 2], 2], 
                   col = linkcol[i], lty = linklty[i], lwd = linklwd[i])
        }
      }
      title("X,Y tps grid")
      b <- c(1, 3)
      tps(M1[, b], M2[, b], gP$n.col.cell, sz = gP$tar.pt.size, 
          pt.bg = gP$tar.pt.bg, grid.col = gP$grid.col, 
          grid.lwd = gP$grid.lwd, grid.lty = gP$grid.lty, 
          refpts = useRefPts)
      if (is.null(links) == FALSE) {
        linkcol <- rep(gP$tar.link.col, nrow(links))[1:nrow(links)]
        linklwd <- rep(gP$tar.link.lwd, nrow(links))[1:nrow(links)]
        linklty <- rep(gP$tar.link.lty, nrow(links))[1:nrow(links)]
        for (i in 1:nrow(links)) {
          segments(M2[links[i, 1], 1], M2[links[i, 1], 
                                          3], M2[links[i, 2], 1], M2[links[i, 2], 3], 
                   col = linkcol[i], lty = linklty[i], lwd = linklwd[i])
        }
      }
      title("Y,Z tps grid")
      layout(1)
      on.exit(par(old.par))
    }
    if (method == "vector") {
      if (axes) {
        plot3d(M1, type = "s", col = gP$pt.bg, size = gP$pt.size, 
               aspect = FALSE, ...)
      }
      if (!axes) {
        plot3d(M1, type = "s", col = gP$pt.bg, size = gP$pt.size, 
               aspect = FALSE, xlab = "", ylab = "", zlab = "", 
               axes = F, ...)
      }
      if (label) {
        text3d(M1, texts = paste(1:dim(M1)[1]), adj = (gP$txt.adj + 
                                                         gP$pt.size), cex = gP$txt.cex, col = gP$txt.col)
      }
      for (i in 1:nrow(M1)) {
        segments3d(rbind(M1[i, ], M2[i, ]), lwd = 2)
      }
      if (!is.null(links)) {
        tarlinkcol <- rep(gP$tar.link.col, nrow(links))[1:nrow(links)]
        tarlinklwd <- rep(gP$tar.link.lwd, nrow(links))[1:nrow(links)]
        tarlinklty <- rep(gP$tar.link.lty, nrow(links))[1:nrow(links)]
        for (i in 1:nrow(links)) {
          segments3d(rbind(M2[links[i, 1], ], M2[links[i, 
                                                       2], ]), col = tarlinkcol[i], lty = tarlinklty[i], 
                     lwd = tarlinklwd[i])
        }
      }
    }
    if (method == "points") {
      if (axes) {
        plot3d(M1, type = "s", col = gP$pt.bg, size = gP$pt.size, 
               aspect = FALSE, ...)
        plot3d(M2, type = "s", col = gP$tar.pt.bg, size = gP$tar.pt.size, 
               add = TRUE)
      }
      if (!axes) {
        plot3d(M1, type = "s", col = gP$pt.bg, size = gP$pt.size, 
               aspect = FALSE, xlab = "", ylab = "", zlab = "", 
               axes = F, ...)
        plot3d(M2, type = "s", col = gP$tar.pt.bg, size = gP$tar.pt.size, 
               add = TRUE)
      }
      if (label) {
        text3d(M1, texts = paste(1:dim(M1)[1]), adj = (gP$txt.adj + 
                                                         gP$pt.size), cex = gP$txt.cex, col = gP$txt.col)
      }
      if (!is.null(links)) {
        linkcol <- rep(gP$link.col, nrow(links))[1:nrow(links)]
        linklwd <- rep(gP$link.lwd, nrow(links))[1:nrow(links)]
        linklty <- rep(gP$link.lty, nrow(links))[1:nrow(links)]
        tarlinkcol <- rep(gP$tar.link.col, nrow(links))[1:nrow(links)]
        tarlinklwd <- rep(gP$tar.link.lwd, nrow(links))[1:nrow(links)]
        tarlinklty <- rep(gP$tar.link.lty, nrow(links))[1:nrow(links)]
        for (i in 1:nrow(links)) {
          segments3d(rbind(M1[links[i, 1], ], M1[links[i, 
                                                       2], ]), col = linkcol[i], lty = linklty[i], 
                     lwd = linklwd[i])
          segments3d(rbind(M2[links[i, 1], ], M2[links[i, 
                                                       2], ]), col = tarlinkcol[i], lty = tarlinklty[i], 
                     lwd = tarlinklwd[i])
        }
      }
    }
    if (method == "surface") {
      if (is.null(mesh)) {
        stop("Surface plotting requires a template mesh3d object (see 'warpRefMesh').")
      }
      warp.PLY <- mesh
      vb <- as.matrix(t(mesh$vb)[, -4])
      cat("\nWarping mesh\n")
      warp <- tps2d3d(vb, M1, M2)
      warp.PLY$vb <- rbind(t(warp), 1)
      shade3d(warp.PLY, meshColor = "legacy", ...)
      return(warp.PLY)
    }
  }
}



wing15.links <- c(1, 7, 
                  2, 6,
                  6, 8,
                  6, 12,
                  12, 13, 
                  13, 14,
                  14, 15, 
                  15, 4, 
                  4, 3,
                  5, 11, 
                  11, 15, 
                  11, 10, 
                  7, 12, 
                  8, 13, 
                  8, 9, 
                  9, 10, 
                  10, 14)

wing15.links <- matrix(wing15.links, ncol = 2, byrow = TRUE)
