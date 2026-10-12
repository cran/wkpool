## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)
has_decido <- requireNamespace("decido", quietly = TRUE)
has_rtriangle <- requireNamespace("RTriangle", quietly = TRUE)
has_silicate <- requireNamespace("silicate", quietly = TRUE)

## ----setup--------------------------------------------------------------------
library(wkpool)

## ----earcut, eval = has_decido------------------------------------------------
x <- wk::as_wkb(c(
  "POLYGON ((0 0, 4 0, 4 4, 0 4, 0 0), (1 1, 1 3, 2 2, 3 3, 3 1, 1 1))",
  "POLYGON ((4 0, 6 0, 6 4, 4 4, 4 0))"
))
pool <- establish_topology(x)
tri <- pool_triangulate(pool)
tri
pool_triangles(tri)

## ----earcut-plot, eval = has_decido, fig.height = 4---------------------------
plot(tri, asp = 1)

## ----earcut-wk, eval = has_decido---------------------------------------------
wk::as_wkt(tri)[1:3]

## ----earcut-merge, eval = has_decido------------------------------------------
tri <- merge_coincident(tri)
edges <- tri_segments(tri)
find_shared_edges(edges)

## ----inlandwaters, eval = has_decido && has_silicate, message = FALSE---------
op <- par(mar = rep(0, 4L))
data("inlandwaters", package = "silicate")
# the geometry column as WKB
iw_geom <- wk::wk_handle(inlandwaters$geom, wk::wkb_writer())
iw <- merge_coincident(establish_topology(iw_geom))
iw_tri <- pool_triangulate(iw)
length(iw_tri)
nrow(pool_vertices(iw_tri)) == nrow(pool_vertices(iw))
plot(iw_tri, border = NA, asp = 1, axes = FALSE, xlab = "", ylab = "")
par(op)

## ----inlandwaters-trim, eval = has_decido && has_silicate---------------------
op <- par(mar = rep(0, 4L))
tt <- pool_triangles(iw_tri)
ev <- pool_vertices(iw_tri)
tx <- ev$x[match(tt$.vx0, ev$.vx)]
ty <- ev$y[match(tt$.vx0, ev$.vx)]
prov <- inlandwaters$Province
nsw <- tt$.feature == which(prov == "New South Wales")
islands <- c(
  "Macquarie Island" = tt$.part[which.min(ty)],
  "Lord Howe Island" = tt$.part[nsw][which.max(tx[nsw])],
  "New Caledonia" = unique(tt$.part[tt$.feature == which(prov == "New Caledonia")])
)
islands
plot(iw_tri[!tt$.part %in% islands], border = NA, asp = 1, axes = FALSE,
     xlab = "", ylab = "")
par(op)

## ----pslg---------------------------------------------------------------------
p <- as_pslg(establish_topology(wk::as_wkb(
  "POLYGON ((0 0, 4 0, 4 4, 0 4, 0 0), (1 1, 1 3, 2 2, 3 3, 3 1, 1 1))"
)))
str(p)

## ----rtriangle, eval = has_rtriangle, fig.height = 4--------------------------
h <- pool_holes(merge_coincident(establish_topology(wk::as_wkb(
  "POLYGON ((0 0, 4 0, 4 4, 0 4, 0 0), (1 1, 1 3, 2 2, 3 3, 3 1, 1 1))"
))))
h
rt <- RTriangle::triangulate(RTriangle::pslg(P = p$P, S = p$S, H = h),
                             q = 30, a = 0.1)
plot(rt$P, asp = 1, pch = ".", xlab = "", ylab = "")
polygon(rt$P[t(cbind(rt$T, NA)), ], col = "grey90")
points(h, col = "firebrick", pch = 19)

## ----back, eval = has_rtriangle-----------------------------------------------
v <- data.frame(.vx = seq_len(nrow(rt$P)), x = rt$P[, 1], y = rt$P[, 2])
rt_tri <- new_wkpool_tri(v, rt$T[, 1], rt$T[, 2], rt$T[, 3])
length(rt_tri)
rt_tri[1:3]

## ----act-rt, eval = has_rtriangle && has_decido && has_silicate---------------
act <- iw_tri[pool_feature(iw_tri) == 1]
act_p <- as_pslg(establish_topology(iw_geom[1]))
act_rt <- RTriangle::triangulate(RTriangle::pslg(P = act_p$P, S = act_p$S),
                                 q = 25)
act_v <- data.frame(.vx = seq_len(nrow(act_rt$P)),
                    x = act_rt$P[, 1], y = act_rt$P[, 2])
act_rt_tri <- new_wkpool_tri(act_v, act_rt$T[, 1], act_rt$T[, 2], act_rt$T[, 3])
c(earcut = length(act), RTriangle = length(act_rt_tri))
op <- par(mfrow = c(1, 2), mar = rep(0.5, 4))
plot(act, col = "grey90", border = "grey40", asp = 1, axes = FALSE,
     xlab = "", ylab = "")
plot(act_rt_tri, col = "grey90", border = "grey40", asp = 1, axes = FALSE,
     xlab = "", ylab = "")
par(op)

## ----compare, eval = has_rtriangle && has_decido && has_silicate--------------
tri_area <- function(P, T) {
  sum(abs((P[T[, 2], 1] - P[T[, 1], 1]) * (P[T[, 3], 2] - P[T[, 1], 2]) -
          (P[T[, 3], 1] - P[T[, 1], 1]) * (P[T[, 2], 2] - P[T[, 1], 2]))) / 2
}
iw_p <- as_pslg(iw)
iw_rt <- RTriangle::triangulate(
  RTriangle::pslg(P = iw_p$P, S = iw_p$S, H = pool_holes(iw))
)
ev <- pool_vertices(iw_tri)
c(
  earcut = tri_area(cbind(ev$x, ev$y), as.matrix(pool_triangles(iw_tri)[1:3])),
  RTriangle = tri_area(iw_rt$P, iw_rt$T)
) / 1e6

## ----islands, eval = has_decido && has_silicate, fig.height = 3.5-------------
op <- par(mfrow = c(1, 3), mar = c(0.5, 0.5, 2, 0.5))
for (nm in names(islands)) {
  plot(iw_tri[tt$.part == islands[nm]], border = "grey30", asp = 1,
       axes = FALSE, xlab = "", ylab = "", main = nm)
}
par(op)

## ----islands-pool, eval = has_decido && has_silicate--------------------------
mq <- iw_tri[tt$.part == islands["Macquarie Island"]]
c(triangles = length(mq), pool_vertices = nrow(pool_vertices(mq)))

## ----islands-count, eval = has_decido && has_silicate-------------------------
sapply(islands, function(k) {
  i <- tt$.part == k
  c(vertices = length(unique(unlist(tt[i, 1:3]))), triangles = sum(i))
})

