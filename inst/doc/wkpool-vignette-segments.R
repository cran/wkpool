## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(wkpool)
library(wk)

two_squares <- c(
  wkt("POLYGON ((0 0, 1 0, 1 1, 0 1, 0 0))"),
  wkt("POLYGON ((1 0, 2 0, 2 1, 1 1, 1 0))")
)

## ----decompose----------------------------------------------------------------
pool <- establish_topology(two_squares)
pool

## ----report-raw---------------------------------------------------------------
topology_report(pool)

## ----merge--------------------------------------------------------------------
pool <- merge_coincident(pool)
topology_report(pool)

## ----shared-------------------------------------------------------------------
find_shared_edges(pool)

## ----neighbours---------------------------------------------------------------
find_neighbours(pool)

## ----internal-----------------------------------------------------------------
find_internal_boundaries(pool)

## ----degree-------------------------------------------------------------------
vertex_degree(pool)

## ----nodes--------------------------------------------------------------------
find_nodes(pool)

## ----arcs---------------------------------------------------------------------
find_arcs(pool)
arc_node_summary(pool)

## ----cycles-------------------------------------------------------------------
cycles <- find_cycles(pool)
cycles

## ----classify-----------------------------------------------------------------
classify_cycles(pool)

## ----roundtrip----------------------------------------------------------------
# Arcs as linestrings (the topological boundary representation)
arcs_to_wkt(pool)

# Cycles as polygons
cycles_to_wkt(pool)

## ----interch, include = FALSE-------------------------------------------------
pslg <- as_pslg(pool)
str(pslg)
#RTriangle::pslg(P = pslg$P, S = pslg$S)


## ----tracks, eval = requireNamespace("traipse", quietly = TRUE)---------------
library(traipse)

# 5 GPS fixes from a Southern Ocean track
x <- c(147.0, 147.5, 148.1, 148.3, 148.0)
y <- c(-42.0, -42.3, -42.5, -42.2, -41.9)

# Each function operates on the implicit segment between consecutive points
track_bearing(x, y)
track_distance(x, y)

