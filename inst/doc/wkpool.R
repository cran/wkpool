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

## ----round-trip---------------------------------------------------------------
# Arcs as linestrings
arcs_to_wkt(pool)

# Cycles as polygons
cycles_to_wkt(pool)

# Raw segments
segments_to_wkt(pool, type = "linestring")[1:3]

## ----pslg---------------------------------------------------------------------
pslg <- as_pslg(pool)
str(pslg)

## ----triangulate, eval = requireNamespace("RTriangle", quietly = TRUE)--------
# A polygon with a hole
holed <- wk::as_wkb(
  "POLYGON ((0 0, 4 0, 4 4, 0 4, 0 0), (1 1, 3 1, 3 3, 1 3, 1 1))"
)
hp <- establish_topology(holed)
hp <- merge_coincident(hp)
pslg_h <- as_pslg(hp)
holes <- hole_points(hp)
tri <- RTriangle::triangulate(
  RTriangle::pslg(P = pslg_h$P, S = pslg_h$S, H = holes)
)

## ----decido, include = FALSE, eval = FALSE------------------------------------
# ### decido (earcut)
# 
# dec <- as_decido(pool)
# str(dec)

## ----subset-------------------------------------------------------------------
pool[1:4]

## ----combine------------------------------------------------------------------
pool_a <- establish_topology(wk::wkt("POLYGON ((0 0, 1 0, 1 1, 0 1, 0 0))"))
pool_b <- establish_topology(wk::wkt("POLYGON ((5 5, 6 5, 6 6, 5 6, 5 5))"))
pool_combine(pool_a, pool_b)

## ----plot, fig.width = 6, fig.height = 4--------------------------------------
plot(pool)

