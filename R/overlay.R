#' Build a grid
#'
#' Creates a regularly spaced grid covering a set of points
#'
#' @param data A matrix or data.frame containing 2 columns corresponding to x
#'     and y coordinates of points
#' @param stepsize numeric, the size of the step for the grid
#'
#' @return A data frame containing x and y coordinates for grid points
#'
#' @author Giuseppe D'Agostino

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
#'
#' @return A data frame containing x and y coordinates for grid centers and for
#'     vertices of each square, named A, B, C, D and ordered clockwise from the
#'     top right one
#'
#' @author Giuseppe D'Agostino

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
#' @author Giuseppe D'Agostino
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
#' @author Giuseppe D'Agostino
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
#' @author Giuseppe D'Agostino
#'
#' @importFrom splancs inout

findHoles <- function(ib_df){

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
#' @param n_dense numeric, the number of points to add to the polygon for more
#'     smoothing. Default is 10.
#'
#' @return A data frame containing ordered coordinates for polygon vertices and
#'    three columns indicating whether they are a hole ("inner") or not ("outer"),
#'    to which cluster they belong, and a sub-clustering to allow \code{ggplot2}
#'    to draw them as holes.
#'
#' @details This is a refactoring of `smoothr::smooth_ksmooth()` to isolate the
#'    necessary code and avoid heavy GDAL-based dependencies. The code has been
#'    simplified assuming `wrap = TRUE` and adding some other bits to handle DF
#'    with other columns as input.
#'
#' @author Matthew Strimas-Mackey, modified by Giuseppe D'Agostino
#'
#' @importFrom stats ksmooth

smoothPolygon <- function(poly, smoothness = 3, min_points = 8, n_dense = 10) {

  if (nrow(poly) < min_points) {

    poly_sm <- poly

    } else {

    poly_coords <- as.matrix(poly[,1:2])

    d_poly = sqrt(rowSums(diff(as.matrix(poly_coords))^2))
    bandwidth = mean(d_poly) * smoothness
    dense_poly = addPoints(poly_coords, steps = n_dense)
    npt = nrow(dense_poly)
    wrapped <- rbind(dense_poly[-npt, ], dense_poly, dense_poly[-1, ])
    d_dense = sqrt(rowSums(diff(as.matrix(wrapped))^2))

    d_x = c(0, cumsum(d_dense))
    poly_sm <- NULL

    for (i in seq_len(ncol(wrapped))) {
      ks <- ksmooth(d_x, wrapped[, i], n.points = length(d_x),
                           kernel = "normal", bandwidth = bandwidth)
      poly_sm <- cbind(poly_sm, ks[["y"]])

      if (i == 1) {
       keep_rows <- (ks$x >= d_x[npt]) & (ks$x <= d_x[(2 * npt - 1)])
      }
    }

  poly_sm <- as.data.frame(poly_sm[keep_rows, ])
  poly_sm[nrow(poly_sm), ] <- poly_sm[1, ]

  #Restore/add original columns
  colnames(poly_sm) <- c("x", "y")
  poly_sm <- as.data.frame(poly_sm)

  for (i in colnames(poly)[3:ncol(poly)]) {
      poly_sm[, i] <- unique(poly[, i])
    }

  }
  return(poly_sm)
}


#' Make overlay contour
#'
#' Make an overlay contour by creating isobands of a square tessellation
#'
#' @param point_coordinates a matrix ord data frame containing x and y coordinates of points
#' @param stepsize numeric, the relative step size used to build the square grid,
#'     expressed as proportion of the maximum span of the range of \code{data}
#' @param minsize numeric, the minimum size of an overlay cluster.
#'     Default is 4 points (one square tile).
#' @param min_pts numeric, the minimum number of points to be enclosed.
#'     Default is 1 (all points)
#' @param offset_prop numeric, the proportion by which polygons are offset
#'     (inflated if outer, deflated if holes). Default is 0.05.
#' @param join_polys logical, should overlapping polygons be joined? Default is TRUE
#' @param smooth logical, perform kernel smoothing of polygons? Default is TRUE
#' @param smoothness numeric, smoothness parameter for kernel smoothing. Default is 3
#'
#' @return A data frame containing ordered coordinates for polygon vertices and
#'    three columns indicating whether they are a hole ("inner") or not ("outer"),
#'    to which cluster they belong, and a sub-clustering to allow \code{ggplot2}
#'    to draw them as holes.
#'
#' @details This function draws an overlay contour that nicely encases points from
#'    2D cloud at a user-selected level of granularity. The algorithm works through
#'    the following steps:
#' \enumerate{
#'    \item Tessellate the space occupied by the data by squares whose side is determined
#'      by the \code{stepsize} parameter. \code{stepsize} is a fraction of the maximum
#'      range in the data, i.e. it takes \code{abs(diff(range(point_coordinates)))} and multiplies it
#'      by the value of \code{stepsize} to obtain the length of each square.
#'    \item Check, for every square, how many points it contains. Discard empty squares
#'      and/or squares that do not contain a minimum of \code{min_pts}. If a point
#'      lies on the edge of a square, it is included.
#'    \item Join non-empty square vertices with isobands (through \code{isoband}) and
#'      obtain one or more encasing  polygons. These will be grouped by \code{cluster}.
#'    - Find which polygons constitute holes ("inner") and which are instead
#'      filled ("outer"), and label them separately.
#'    \item Optional: offset the encasing polygons by a proportional margin defined
#'      by \code{offset_prop}, using \code{polyclip}. Given the range of polygon
#'      coordinates, the total offset will be \code{abs(diff(range(coordinates))) * offset_prop}.
#'      "outer" polygons will be inflated (positive offset) and "inner" polygons
#'      will be deflated (negative offset).
#'    \item Optional: join offset polygons. If the offset increases the size of polygons
#'      so that they overlap, \code{polyclip} is used to calculate the union of
#'      polygons.
#'  }
#'    The resulting object is a \code{data.frame} containing the following columns:
#' \itemize{
#'    \item *x, y*: ordered coordinates for each polygon
#'    \item *cluster* the original cluster assigned by \code{isoband}
#'    \item *hole* either "outer" (not a hole) or "inner" (hole)
#'    \item *cluster_hole* group level for ggplot2 aesthetic
#'    \item *id_hole* subgroup level for ggplot2 aesthetic
#'  }
#'    Note that it is possible to use this output in base R graphics, but holes
#'    will be drawn as other polygons without removing space.
#'
#' @author Giuseppe D'Agostino
#'
#' @importFrom polyclip polyoffset polyclip
#' @importFrom scales rescale
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
#' overlay <- makeOverlay(dat,min_pts = 1, stepsize = 0.02,
#'    minsize = 4, offset_prop = 0.08, join_polys = FALSE)
#' plot(dat, pch = 16, cex = 0.5, xlim = range(overlay[,1:2]), ylim = range(overlay[,1:2]))
#' for(i in unique(overlay$cluster)) polygon(overlay[overlay$cluster == i, 1:2])
#'
#' \dontrun{
#' # With ggplot2, showcasing holes
#' library(ggplot2)
#' dat =  as.data.frame(matrix(rnorm(1000), ncol = 2))
#' colnames(dat) = c("x", "y")
#' dat = dat[abs(dat$x) > 0.5 | abs(dat$y) > 0.5,]

#' overlay <- makeOverlay(dat, min_pts = 1, stepsize = 0.02,
#'    minsize = 4, offset_prop = 0.01, join_polys = TRUE)

#' ggplot(data = dat, aes(x = x, y = y)) +
#'   geom_point() +
#'   geom_polygon(data = overlay,
#'                aes_string(x = "x", y = "y", group = "cluster_hole", subgroup = "id_hole"),
#'                color = "red", fill = "red", alpha = 0.3) +
#'   theme_bw()
#'}
#'
#' @export



makeOverlay = function(point_coordinates, stepsize = 0.1, minsize = 4,
                       min_pts = 2, offset_prop = 0.05,
                       join_polys = TRUE, smooth = TRUE,
                       smoothness = 3) {

  point_coordinates = as.matrix(point_coordinates)
  data_scaled = point_coordinates
  ranges = apply(point_coordinates, 2, range)
  range_diffs = apply(ranges, 2, diff)
  to_rescale = which.min(range_diffs)

  data_scaled[,to_rescale] = rescale(data_scaled[,to_rescale],
                                     from = ranges[,to_rescale],
                                     to = ranges[,-to_rescale])

  stepsize = abs(diff(range(data_scaled))) * stepsize

  conts = makeContour(data_scaled, stepsize = stepsize,
                      minsize = minsize, min_pts = min_pts)

  conts_holed = findHoles(do.call(rbind, conts))

  if(!is.null(offset_prop)) {
    conts_list = split(conts_holed, conts_holed$id_hole)

    outer_offset = unlist(lapply(conts_list, function(x) unique(x$hole) == "outer"))
    conts_outer_offset = conts_list[outer_offset]
    for(i in seq_len(length(conts_outer_offset))) {
      offset_points = polyoffset(conts_outer_offset[[i]][,1:2],
                                 delta = abs(diff(range(data_scaled[,1:2]))) * offset_prop)
      offset_coords = data.frame("x" = as.numeric(offset_points[[1]]$x),
                                 "y" = as.numeric(offset_points[[1]]$y))

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
        offset_points = polyoffset(conts_inner_offset[[i]][,1:2],
                                   delta = -abs(diff(range(data_scaled[,1:2]))) * offset_prop)
        if(length(offset_points) < 1) next
        offset_coords = data.frame("x" = as.numeric(offset_points[[1]]$x),
                                   "y" = as.numeric(offset_points[[1]]$y))

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

    conts_list_outer = split(conts_holed[conts_holed$hole == "outer", c("x", "y", "cluster")],
                             conts_holed[conts_holed$hole == "outer", "cluster"])
    conts_inner = conts_holed[conts_holed$hole == "inner", c("x", "y", "cluster")]
    conts_list_union = Reduce(x = conts_list_outer, f = function(x, y)
      polyclip(x, y, op = "union"))
    conts_list_union = lapply(seq_len(length(conts_list_union)), function(x)
      data.frame("x" = conts_list_union[[x]]$x, "y" = conts_list_union[[x]]$y, "cluster" = x))

    conts_joined = do.call(rbind, conts_list_union)

     if(any(unique(conts_joined$cluster) %in% conts_inner$cluster)) {
       conts_inner$cluster = conts_inner$cluster + max(conts_joined$cluster)
     }
     conts_all = rbind(conts_joined, conts_inner)
     conts_holed = findHoles(conts_all)
  }

  if(smooth)  {
    conts_holed = do.call(rbind,
                          lapply(split(conts_holed, conts_holed$id_hole),
                                 function(x) smoothPolygon(x,
                                                           smoothness = smoothness,
                                                           min_points  = minsize)))
  }


  conts_holed[,to_rescale] = rescale(conts_holed[,to_rescale],
                                     from = ranges[,-to_rescale],
                                     to = ranges[,to_rescale])

  return(conts_holed)
}


#' Add points
#'
#' Increases the number of points in a polygon while maintaining the shape
#'
#' @param poly a data frame containing ordered coordinates with polygon vertices
#' @param steps numeric, the number of points that should be added between
#'    each point. Default is 10.
#'
#' @return A data frame containing densified coordinates for polygon vertices and
#'    any other column (assuming original columns contained unique values).
#'
#' @details Internal use only.
#'
#' @author Giuseppe D'Agostino


addPoints <- function(poly, steps = 5) {

  colnames(poly) <- NULL
  polygon_coords = as.matrix(poly[,1:2])
  new_coords = poly[1,]

  for(i in 1:(nrow(poly)-1))
  {
    new_xy = cbind(seq(polygon_coords[i,1],
                       polygon_coords[i+1,1], length.out = steps+1),
                   seq(polygon_coords[i,2],
                       polygon_coords[i+1,2], length.out = steps+1))
    tmp_coords = polygon_coords[i,]
    new_coords = rbind(new_coords, tmp_coords, new_xy)
  }

  new_xy_last = cbind(seq(polygon_coords[nrow(polygon_coords),1],
                          polygon_coords[1,1], length.out = steps+1),
                      seq(polygon_coords[nrow(polygon_coords),2],
                          polygon_coords[1,2], length.out = steps+1))

  new_coords = new_coords[2:nrow(new_coords),]
  new_coords = new_coords[!duplicated(new_coords),]
  new_coords = as.data.frame(rbind(new_coords, new_xy_last))

  rownames(new_coords) = NULL
  colnames(new_coords) <- c("x", "y")

  for (i in colnames(poly)[3:ncol(poly)]) {
    new_coords[, i] <- unique(poly[, i])
  }

  return(new_coords)
}




