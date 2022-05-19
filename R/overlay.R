#' Build a grid
#'
#' Creates a regularly spaced grid covering a set of points
#'
#' @param data A matrix or data.frame containing 2 columns corresponding to x
#'     and y coordinates of points
#' @param stepsize numeric, the size of the step for the grid
#'
#' @return A data frame containing x and y coordinates for grid points

makeGrid = function(data, stepsize){
  maxstep = max(range(data)) + stepsize
  minstep = min(range(data)) - stepsize
  steps = seq(minstep, maxstep, by = stepsize)
  dots = as.data.frame(expand.grid(steps, steps))
  colnames(dots) = c("x", "y")
  return(dots)
}

#' Build a grid of squares
#'
#' Creates a grid of squares centered on an input grid
#'
#' @param grid output of \code{makeGrid}
#' @param stepsize numeric, the size of the step for the grid
#'
#' @return A data frame containing x and y coordinates for grid centers and for
#'     vertices of each square, named A, B, C, D and ordered clockwise from the
#'     top right one

makeSquares = function(grid, stepsize) {

  grid$Ax = grid$x + stepsize/2
  grid$Ay = grid$y + stepsize/2

  grid$Bx = grid$x + stepsize/2
  grid$By = grid$y - stepsize/2

  grid$Cx = grid$x - stepsize/2
  grid$Cy = grid$y - stepsize/2

  grid$Dx = grid$x - stepsize/2
  grid$Dy = grid$y + stepsize/2

  grid$square = seq_len(nrow(grid))
  return(grid)
}

#' Build a bounding box
#'
#' Creates a bounding box around a subset of square vertices on a grid
#'
#' @param squares output of \code{makeSquares}
#' @param stepsize numeric, the size of the step for the grid
#'
#' @return A data frame containing x and y coordinates for square vertices
#'    (with a value = 1) and the bounding box filled with value = 0. The box is
#'    larger than the set of vertices by \code{stepsize}.

makeBoundingBox <- function(squares,
                            stepsize = 1) {

  bounding_box <- c(range(squares[, 1]), range(squares[, 2])) + c(-stepsize, stepsize, -stepsize, stepsize)
  bounding_box_filled <- expand.grid(
    c(seq(bounding_box[1], bounding_box[2], by = stepsize)),
    c(seq(bounding_box[3], bounding_box[4], by = stepsize))
  )

  colnames(bounding_box_filled) <- c("x", "y")

  bounding_box_all <- rbind(bounding_box_filled, squares[, 1:2])

  bounding_box_all$value <- 1
  bounding_box_all$value[1:nrow(bounding_box_filled)] <- 0
  bounding_box_all[,1:2] = round(bounding_box_all[,1:2], digits = 4)
  bounding_box_all <- bounding_box_all[ !duplicated(bounding_box_all[, 1:2], fromLast = TRUE), ]

  return(bounding_box_all)
}

#' Make an isoband contour
#'
#' Creates an isoband around a cloud of points assigned to a square grid
#'
#' @param data a matrix or data frame containing a set of points
#' @param stepsize numeric, the spacing of the points
#' @param minsize numeric, the minimum size of the contour (in number of vertices)
#'    to be retained. Default is 8 (2 squares)
#' @param min_pts numeric, the minimum number of points to be included in a square
#'    for the contour to be drawn. Default is 1.
#'
#' @return A list of data frames containing the ordered set of coordinates to draw
#'     isoband polygons together with their clustering (disjointed polygon)
#'
#' @importFrom splancs inout
#' @importFrom reshape2 acast
#' @importFrom isoband isobands
#' @importFrom scales rescale

makeContour = function(data, stepsize, minsize = 8, min_pts = 1) {
  min_pts = as.integer(min_pts)
  if(min_pts < 1) stop("min_pts must be an integer strictly >= 1")

  gridpoints = makeGrid(data, stepsize)
  squares = makeSquares(gridpoints, stepsize)

  squares_list = split(squares[,3:10], squares$square)

  squares_list = lapply(squares_list, function(x) {
    square = as.data.frame(matrix(x, ncol = 2, byrow = TRUE))
    colnames(square) = c("x", "y")
    return(square)
  })

  squares_in = lapply(squares_list, function(x) which(inout(pts = data, poly = x, bound = TRUE) == TRUE))
  squares_in = squares_in[unlist(lapply(squares_in, function(x) length(x) >= min_pts))]

  squares_ok = as.data.frame(do.call(rbind, squares_list[names(squares_in)]))

  squares_empty = as.data.frame(do.call(rbind, squares_list[!names(squares_list) %in% names(squares_in)]))

  squares_ok$value = 1
  squares_empty$value = 0

  squares_all = rbind(squares_ok, squares_empty)
  colnames(squares_all) = c("x", "y", "z")
  squares_all$x = unlist(squares_all$x)
  squares_all$y = unlist(squares_all$y)
  squares_all$z = unlist(squares_all$z)
  squares_all = squares_all[!duplicated(squares_all[,1:2], fromLast = FALSE),]

  boxed = makeBoundingBox(squares_all[squares_all$z == 1,1:2], stepsize = stepsize)

  m = t(acast(boxed, formula = x ~ y, value.var = "value"))

  ib = isobands(x = seq_len(ncol(m)), y = seq_len(nrow(m)), m, levels_low = 0, levels_high = 1)

  ib_df = data.frame(
    "x" = ib[[1]]$x,
    "y" = ib[[1]]$y,
    "cluster" = ib[[1]]$id
  )

  ib_df <- ib_df[which(!ib_df$x %in% range(ib_df$x) & !ib_df$y %in% range(ib_df$y)), ]

  ib_df$x <- rescale(ib_df$x, to = range(squares_all[squares_all$z == 1, 1]))
  ib_df$y <- rescale(ib_df$y, to = range(squares_all[squares_all$z == 1, 2]))

  ib_list = split(ib_df, ib_df$cluster)
  if(!is.null(minsize)) {
    keep = unlist(lapply(ib_list, function(x) nrow(x) >= minsize))
    ib_list = ib_list[keep]
  }
  return(ib_list)
}


#' Find holes
#'
#' Finds which polygons represent holes within a set of polygons
#'
#' @param ib_df a data frame containing a series of isoband contours created by
#'
#' @return A data frame containing ordered coordinates for polygon vertices and
#'    three columns indicating whether they are a hole ("inner") or not ("outer"),
#'    to which cluster they belong, and a sub-clustering to allow \code{ggplot2}
#'    to draw them as holes.
#'
#' @importFrom splancs inout

find_holes <- function(ib_df){

  if (length(unique(ib_df$cluster)) > 1) {
    hole_list <- list()

    for (i in unique(ib_df$cluster)) {
      hole_list[[as.character(i)]] <- lapply(
        unique(ib_df$cluster)[unique(ib_df$cluster) != i],
        function(x) {
          any(inout(
            pts = ib_df[ib_df$cluster == x, 1:2],
            poly = ib_df[ib_df$cluster == i, 1:2]
          ))
        }
      )
    }

    for (i in seq_len(length(hole_list))) {
      names(hole_list[[i]]) <- unique(ib_df$cluster)[unique(ib_df$cluster) != names(hole_list)[i]]
    }

    hmat <- do.call(rbind, hole_list)
    hole_per_cluster <- apply(hmat, 1, function(x) colnames(hmat)[x == TRUE])
    hole_per_cluster <- hole_per_cluster[which(lengths(hole_per_cluster) > 0)]
    hole_per_cluster <- data.frame(
      "cluster_outer" = rep(names(hole_per_cluster), lengths(hole_per_cluster)),
      "cluster_inner" = unlist(hole_per_cluster)
    )

    hole_per_cluster$outer_size <- table(ib_df$cluster)[hole_per_cluster$cluster_outer]
    hole_per_cluster$inner_size <- table(ib_df$cluster)[hole_per_cluster$cluster_inner]
    hole_per_cluster <- hole_per_cluster[which(hole_per_cluster$inner_size < hole_per_cluster$outer_size), ]

    if(any(duplicated(hole_per_cluster$cluster_inner))) {
      dups <- unique(hole_per_cluster$cluster_inner[which(duplicated(hole_per_cluster$cluster_inner))])
      discard_list <- list()
      for(i in seq_len(length(dups))) {
        check_size <- hole_per_cluster[hole_per_cluster$cluster_inner == dups[[i]],]
        discard_list[[i]] <- rownames(check_size)[-which.max(check_size$outer_size)]
      }
      hole_per_cluster <- hole_per_cluster[setdiff(rownames(hole_per_cluster), unlist(discard_list)),]
    }
    ib_df$hole <- "outer"

    ib_df$cluster_hole <- ib_df$cluster

    for (i in unique(hole_per_cluster$cluster_inner)) {
      ib_df$hole[ib_df$cluster == i] <- "inner"
      ib_df$cluster_hole[ib_df$cluster == i] <- hole_per_cluster[hole_per_cluster$cluster_inner == i, 1]
    }

    ib_df$id_hole <- paste0(ib_df$structure, "_", ib_df$cluster_hole, "_", ib_df$cluster)
  } else {

    ib_df$hole <- "outer"
    ib_df$cluster_hole <- ib_df$cluster
    ib_df$id_hole <- paste0(ib_df$structure, "_", ib_df$cluster_hole, "_", ib_df$cluster)
  }
  return(ib_df)
}

#' Smooth polygons
#'
#' Uses kernel smoothing to smooth polygons
#'
#' @param poly a data frame containing ordered coordinates with polygon vertices
#' @param smoothness numeric, the extent of kernel smoothing. Higher means
#'     rounder shapes. Default is 3.
#' @param min_points numeric, the minimum number of vertices to smooth.
#'     Default is 8.
#'
#' @return A data frame containing ordered coordinates for polygon vertices and
#'    three columns indicating whether they are a hole ("inner") or not ("outer"),
#'    to which cluster they belong, and a sub-clustering to allow \code{ggplot2}
#'    to draw them as holes.
#'
#' @importFrom smoothr smooth_ksmooth
#' @importFrom utils installed.packages

poly_smooth <- function(poly, smoothness = 3, min_points = 8) {

  if (!"smoothr" %in% rownames(installed.packages())) stop("In order to use smoothing you must first install the package `smoothr`.")

  if (nrow(poly) < min_points) {
    poly_sm <- poly[, 1:2]
  } else {
    poly_sm <- as.data.frame(smooth_ksmooth(as.matrix(poly[, 1:2]),
                                            smoothness = smoothness,
                                            wrap = TRUE))
  }

  # Restore/add original columns
  colnames(poly_sm) <- c("x", "y")

  for (i in colnames(poly)[3:ncol(poly)]) {
    poly_sm[, i] <- unique(poly[, i])
  }

  return(poly_sm)
}

#' Make overlay contour
#'
#' Make an overlay contour by creating isobands of a square tessellation
#'
#' @param data a matrix ord data frame containing x and y coordinates of points
#' @param stepsize numeric, the relative step size used to build the square grid,
#'     expressed as proportion of the maximum span of the range of \code{data}
#' @param minsize numeric, the minimum size of an overlay cluster
#' @param min_pts numeric, the minimum number of points to be enclosed
#' @param offset_prop numeric, the proportion by which polygons are offset
#'     (inflated if outer, deflated if holes)
#' @param join_polys logical, should overlapping polygons be joined? Default is TRUE
#' @param smooth logical, perform kernel smoothing of polygons? Default is TRUE
#' @param smoothness numeric, smoothness parameter for kernel smoothing. Default is 3
#'
#' @return A data frame containing ordered coordinates for polygon vertices and
#'    three columns indicating whether they are a hole ("inner") or not ("outer"),
#'    to which cluster they belong, and a sub-clustering to allow \code{ggplot2}
#'    to draw them as holes.
#'
#' @importFrom polyclip polyoffset polyclip
#'
#' @examples
#'
#' # Normal usage
#' dat <- matrix(rnorm(1000), ncol = 2)
#' overlay <- makeOverlay(dat, min_pts = 1, stepsize = 0.06, minsize = 4)
#' plot(dat, pch = 16, cex = 0.5, xlim = range(overlay[,1:2]), ylim = range(overlay[,1:2]))
#' for(i in unique(overlay$cluster)) polygon(overlay[overlay$cluster == i, 1:2])
#'
#' # Decreasing step size (increasing granularity)
#' overlay <- makeOverlay(dat, min_pts = 1, stepsize = 0.02, minsize = 4)
#' plot(dat, pch = 16, cex = 0.5, xlim = range(overlay[,1:2]), ylim = range(overlay[,1:2]))
#' for(i in unique(overlay$cluster)) polygon(overlay[overlay$cluster == i, 1:2])
#'
#' # Increasing offset
#' overlay <- makeOverlay(dat,min_pts = 1, stepsize = 0.02, minsize = 4, offset_prop = 0.08)
#' plot(dat, pch = 16, cex = 0.5, xlim = range(overlay[,1:2]), ylim = range(overlay[,1:2]))
#' for(i in unique(overlay$cluster)) polygon(overlay[overlay$cluster == i, 1:2])
#'
#' # Increasing offset without joining polygons
#' overlay <- makeOverlay(dat,min_pts = 1, stepsize = 0.02, minsize = 4, offset_prop = 0.08, join_polys = FALSE)
#' plot(dat, pch = 16, cex = 0.5, xlim = range(overlay[,1:2]), ylim = range(overlay[,1:2]))
#' for(i in unique(overlay$cluster)) polygon(overlay[overlay$cluster == i, 1:2])
#'
#' @export

makeOverlay = function(data, stepsize, minsize, min_pts = 2, offset_prop = 0.01, join_polys = TRUE, smooth = TRUE, smoothness = 3) {

  stepsize = abs(diff(range(data))) * stepsize

  conts = makeContour(data, stepsize = stepsize, minsize = minsize, min_pts = min_pts)

  if(smooth) conts = lapply(conts, function(x) poly_smooth(x, smoothness = smoothness, min_points  = minsize))

  conts_holed = find_holes(do.call(rbind, conts))

  if(!is.null(offset_prop)) {
    conts_list = split(conts_holed, conts_holed$id_hole)

    outer_offset = unlist(lapply(conts_list, function(x) unique(x$hole) == "outer"))
    conts_outer_offset = conts_list[outer_offset]
    for(i in seq_len(length(conts_outer_offset))) {
      offset_points = polyoffset(conts_outer_offset[[i]][,1:2], delta = abs(diff(range(data[,1:2]))) * offset_prop)
      offset_coords = data.frame("x" = as.numeric(offset_points[[1]]$x), "y" = as.numeric(offset_points[[1]]$y))

      offset_coords$cluster = unique(conts_outer_offset[[i]]$cluster)
      offset_coords$hole = unique(conts_outer_offset[[i]]$hole)
      offset_coords$cluster_hole = unique(conts_outer_offset[[i]]$cluster_hole)
      offset_coords$id_hole = unique(conts_outer_offset[[i]]$id_hole)
      conts_outer_offset[[i]] = offset_coords
    }

    inner_offset = unlist(lapply(conts_list, function(x) unique(x$hole) == "inner"))

    if(sum(inner_offset) >= 1) {
      conts_inner_offset = conts_list[inner_offset]
      for(i in seq_len(length(conts_inner_offset))) {
        offset_points = polyoffset(conts_inner_offset[[i]][,1:2], delta = -abs(diff(range(data[,1:2]))) * offset_prop)
        if(length(offset_points) < 1) next
        offset_coords = data.frame("x" = as.numeric(offset_points[[1]]$x), "y" = as.numeric(offset_points[[1]]$y))

        offset_coords$cluster = unique(conts_inner_offset[[i]]$cluster)
        offset_coords$hole = unique(conts_inner_offset[[i]]$hole)
        offset_coords$cluster_hole = unique(conts_inner_offset[[i]]$cluster_hole)
        offset_coords$id_hole = unique(conts_inner_offset[[i]]$id_hole)
        conts_inner_offset[[i]] = offset_coords
      }
      conts_list[inner_offset] = conts_inner_offset
    }
    conts_list[outer_offset] = conts_outer_offset

    conts_holed = as.data.frame(do.call(rbind, conts_list))
  }

  if(join_polys & length(unique(conts_holed$cluster[conts_holed$hole == "outer"])) > 1) {
    conts_list_outer = split(conts_holed[conts_holed$hole == "outer", ], conts_holed[conts_holed$hole == "outer", "cluster_hole"])
    conts_list_inner = split(conts_holed[conts_holed$hole == "inner", ], conts_holed[conts_holed$hole == "inner", "cluster_hole"])
    conts_list_union = Reduce(x = conts_list_outer, f = function(x, y) polyclip(x, y, op = "union"))
    conts_list_union = lapply(seq_len(length(conts_list_union)), function(x) data.frame("x" = conts_list_union[[x]]$x, "y" = conts_list_union[[x]]$y, "cluster" = x))
    conts_holed = rbind(do.call(rbind, conts_list_union), do.call(rbind, conts_list_inner))
    conts_holed = find_holes(conts_holed)  
   }

  return(conts_holed)

}
