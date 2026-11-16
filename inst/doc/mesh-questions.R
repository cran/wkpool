## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 4.5
)
has_geos <- requireNamespace("geos", quietly = TRUE)
has_cdtr <- requireNamespace("cdtr", quietly = TRUE)
has_nanoarrow <- requireNamespace("nanoarrow", quietly = TRUE)

## ----setup--------------------------------------------------------------------
library(wkpool)

## ----grid---------------------------------------------------------------------
n <- 6L
g <- expand.grid(i = 0:n, j = 0:n)
v <- data.frame(.vx = seq_len(nrow(g)), x = as.double(g$i), y = as.double(g$j))
set.seed(1)
inner <- v$x > 0 & v$x < n & v$y > 0 & v$y < n
v$x[inner] <- v$x[inner] + runif(sum(inner), -0.3, 0.3)
v$y[inner] <- v$y[inner] + runif(sum(inner), -0.3, 0.3)
id <- function(i, j) j * (n + 1L) + i + 1L
cell <- expand.grid(i = seq_len(n) - 1L, j = seq_len(n) - 1L)
a <- id(cell$i, cell$j)
b <- id(cell$i + 1L, cell$j)
c <- id(cell$i + 1L, cell$j + 1L)
d <- id(cell$i, cell$j + 1L)
tri <- new_wkpool_tri(v, vx0 = c(a, a), vx1 = c(b, c), vx2 = c(c, d))
tri

## ----groups-------------------------------------------------------------------
region <- ifelse(cell$i < 3, ifelse(cell$j < 3, 1L, 2L), ifelse(cell$j < 3, 3L, 4L))
region[cell$i == 1 & cell$j == 1] <- NA   # a hole in region 1
region[cell$i == 5 & cell$j == 5] <- NA   # a notch out of region 4
group <- c(region, region)                # one value per triangle

## ----plot-groups--------------------------------------------------------------
plot(tri_dissolve(tri, group), col = hcl.colors(4, "Set 2"), asp = 1)

## ----measures-----------------------------------------------------------------
tri_measures(tri, group)

## ----relate-------------------------------------------------------------------
r <- tri_relate(tri, group)
r

## ----touches------------------------------------------------------------------
touches <- substr(r$relate, 1, 1) == "F"
r[touches & substr(r$relate, 5, 5) == "0", c("a", "b")]

## ----overlap------------------------------------------------------------------
centre <- which(cell$i %in% 2:3 & cell$j %in% 2:3)
gt <- data.frame(
  triangle = c(seq_along(tri), centre, centre + nrow(cell)),
  group = c(group, rep(5L, 2 * length(centre)))
)
gt <- gt[!is.na(gt$group), ]
r5 <- tri_relate(tri, gt)
r5[r5$b == 5, ]

## ----locate-------------------------------------------------------------------
set.seed(2)
pts <- wk::xy(runif(10, -0.5, 6.5), runif(10, -0.5, 6.5))
i <- tri_locate(tri, pts)
data.frame(x = round(wk::wk_coords(pts)$x, 2), y = round(wk::wk_coords(pts)$y, 2),
           triangle = i, region = group[i])

## ----structure----------------------------------------------------------------
tri_structure(tri, group)

## ----boundary-----------------------------------------------------------------
bnd <- tri_boundary(tri, group, once = TRUE)
plot(tri, col = NA, border = "grey85", asp = 1)
plot(bnd, col = "black", add = TRUE)

## ----clip---------------------------------------------------------------------
clipped <- tri_clip_rect(tri, c(1.5, 1.5, 4.5, 4.5))
clipped_group <- group[attr(clipped, "source")]
tri_measures(clipped, clipped_group)[, c("group", "area", "perimeter")]
wk::as_wkt(tri_convex_hull(tri, group))[1]

## ----prepare------------------------------------------------------------------
ptri <- tri_prepare(tri, group)
tri_is_prepared(ptri)
identical(tri_structure(ptri, group), tri_structure(tri, group))

## ----layer-b, eval = has_cdtr-------------------------------------------------
m <- 3L
gb <- expand.grid(i = 0:m, j = 0:m)
th <- 0.3
vb <- data.frame(.vx = seq_len(nrow(gb)),
                 x = 1.8 * (gb$i * cos(th) - gb$j * sin(th)) + 1.5,
                 y = 1.8 * (gb$i * sin(th) + gb$j * cos(th)) - 0.4)
idb <- function(i, j) j * (m + 1L) + i + 1L
cb <- expand.grid(i = seq_len(m) - 1L, j = seq_len(m) - 1L)
tri_b <- new_wkpool_tri(
  vb,
  vx0 = c(idb(cb$i, cb$j), idb(cb$i, cb$j)),
  vx1 = c(idb(cb$i + 1L, cb$j), idb(cb$i + 1L, cb$j + 1L)),
  vx2 = c(idb(cb$i + 1L, cb$j + 1L), idb(cb$i, cb$j + 1L))
)
group_b <- rep(c(1L, 1L, 2L, 3L, 3L, 2L, 3L, NA, 2L), 2)
ov <- tri_overlay(tri, tri_b, group, group_b)
ov

## ----plot-overlay, eval = has_cdtr--------------------------------------------
plot(ov$mesh, col = NA, border = "grey80", asp = 1)
plot(tri_dissolve(tri, group), col = NA, border = "black", lwd = 2, add = TRUE)
plot(tri_dissolve(tri_b, group_b), col = NA, border = "firebrick", lwd = 2, add = TRUE)

## ----overlay-relate, eval = has_cdtr------------------------------------------
tri_overlay_relate(ov)

## ----overlay-lines, eval = has_cdtr-------------------------------------------
lines <- wk::wkt(c("LINESTRING (-1 1.5, 7 1.5)", "LINESTRING (1 1, 2 2, 1 2, 1 1)",
                   "LINESTRING (3 0, 3 6)"))
ovl <- tri_overlay(tri, lines, group)
head(tri_overlay_relate(ovl))
ovp <- tri_overlay(lines, wk::xy(c(1, 3, 5), c(1.5, 3, 1.5)))
tri_overlay_relate(ovp)

## ----overlay-area, eval = has_cdtr--------------------------------------------
mem <- ov$members
g <- ov$groups
in_x <- mem[g$layer[mem$group] == "x", ]
in_y <- mem[g$layer[mem$group] == "y", ]
both <- merge(in_x, in_y, by = "item")
pieces <- data.frame(triangle = both$item,
                     group = as.integer(factor(paste(both$group.x, both$group.y))))
area <- tri_measures(ov$mesh, pieces)$area
key <- levels(factor(paste(both$group.x, both$group.y)))
data.frame(x = g$id[as.integer(sub(" .*", "", key))],
           y = g$id[as.integer(sub(".* ", "", key))], area = round(area, 3))

## ----arrow, eval = has_nanoarrow----------------------------------------------
arr <- tri_as_arrow(tri, group)
arr$length
back <- tri_from_arrow(arr)
all.equal(tri_measures(back$mesh, back$group), tri_measures(tri, group))

## ----geos, eval = has_geos----------------------------------------------------
polys <- geos::as_geos_geometry(tri_dissolve(tri, group))
all.equal(tri_measures(tri, group)$area, geos::geos_area(polys))
all(r$relate == geos::geos_relate(polys[r$a], polys[r$b]))

## ----geos-overlay, eval = has_geos && has_cdtr--------------------------------
polys_b <- geos::as_geos_geometry(tri_dissolve(tri_b, group_b))
ro <- tri_overlay_relate(ov, expand.grid(x = 1:4, y = 1:3))
all(ro$relate == geos::geos_relate(polys[ro$x], polys_b[ro$y]))
rl <- tri_overlay_relate(ovl, expand.grid(x = 1:4, y = 1:3))
all(rl$relate == geos::geos_relate(polys[rl$x], geos::as_geos_geometry(lines)[rl$y]))

