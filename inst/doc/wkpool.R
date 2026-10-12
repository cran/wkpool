## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(wkpool)

## ----basic-setup--------------------------------------------------------------
library(wk)
library(wkpool)

# Two adjacent unit squares
poly1 <- wkt("POLYGON ((0 0, 1 0, 1 1, 0 1, 0 0))")
poly2 <- wkt("POLYGON ((1 0, 2 0, 2 1, 1 1, 1 0))")
geoms <- c(poly1, poly2)

## ----establish----------------------------------------------------------------
pool <- establish_topology(geoms)
pool

## ----report-before------------------------------------------------------------
topology_report(pool)

## ----accessors----------------------------------------------------------------
pool_vertices(pool)
pool_segments(pool)
pool_feature(pool)

## ----merge--------------------------------------------------------------------
pool <- merge_coincident(pool)
topology_report(pool)

## ----shared-edges-------------------------------------------------------------
find_shared_edges(pool)

## ----internal-boundaries------------------------------------------------------
find_internal_boundaries(pool)

## ----neighbours---------------------------------------------------------------
find_neighbours(pool)

## ----cycles-------------------------------------------------------------------
cycles <- find_cycles(pool)
cycles

## ----winding------------------------------------------------------------------
# Area of the first cycle
cycle_signed_area(cycles[[1]], pool_vertices(pool))

# Classify all cycles
classify_cycles(pool)

## ----arc-node-----------------------------------------------------------------
vertex_degree(pool)
find_nodes(pool)
find_arcs(pool)
arc_node_summary(pool)

## ----arc-node-quotient--------------------------------------------------------
find_nodes(pool, quotient = TRUE)
find_arcs(pool, quotient = TRUE)

## ----round-trip---------------------------------------------------------------
# Arcs as linestrings
arcs_to_wkt(pool)

# Cycles as polygons
cycles_to_wkt(pool)

# Raw segments
segments_to_wkt(pool, type = "linestring")[1:3]

## ----period-x-----------------------------------------------------------------
hemi <- wk::as_wkb(c(
  "POLYGON ((-180 -60, 0 -60, 0 60, -180 60, -180 -60))",
  "POLYGON ((0 -60, 180 -60, 180 60, 0 60, 0 -60))"
))
hp <- establish_topology(hemi)
nrow(find_shared_edges(merge_coincident(hp)))
nrow(find_shared_edges(merge_coincident(hp, period_x = 360)))

## ----tri, eval = requireNamespace("decido", quietly = TRUE), fig.width = 6, fig.height = 4----
tri <- pool_triangulate(pool)
plot(tri, asp = 1)

## ----subset-------------------------------------------------------------------
pool[1:4]

## ----combine------------------------------------------------------------------
pool_a <- establish_topology(wk::wkt("POLYGON ((0 0, 1 0, 1 1, 0 1, 0 0))"))
pool_b <- establish_topology(wk::wkt("POLYGON ((5 5, 6 5, 6 6, 5 6, 5 5))"))
pool_combine(pool_a, pool_b)

## ----from-tables--------------------------------------------------------------
v <- data.frame(.vx = 1:4, x = c(0, 1, 1, 0), y = c(0, 0, 1, 1))
s <- data.frame(.vx0 = 1:4, .vx1 = c(2:4, 1L), .feature = 1L)
wkpool(v, s)

## ----plot, fig.width = 6, fig.height = 4--------------------------------------
plot(pool)

