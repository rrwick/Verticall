#!/usr/bin/env Rscript

# This script take a PHYLIP distance matrix as input and outputs a number which represents how
# tree-like the distances are, with lower values indicating more tree-like distances. If the
# distances are perfectly compatible with a tree, then this script will return 0.0 (or at least a
# value close to zero due to floating point imprecision).

# Copyright 2022 Ryan Wick (rrwick@gmail.com)
# https://github.com/rrwick/Verticall

# This file is part of Verticall. Verticall is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by the Free Software Foundation,
# either version 3 of the License, or (at your option) any later version. Verticall is distributed
# in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details. You should have received a copy of the GNU General Public License along with Verticall.
# If not, see <https://www.gnu.org/licenses/>.

library(ape)
options(digits = 16) 


main <- function() {
    args = commandArgs(trailingOnly=TRUE)

    # Build a neighbour-joining tree from the distances in the PHYLIP matrix:
    matrix <- load_phylip_distances(args[1])
    tree <- bionj(matrix)
    tree$edge.length <- pmax(tree$edge.length, 0.0)  # set negative branch lengths to zero

    # Get all pairwise distances from the tree:
    matrix_from_tree <- cophenetic.phylo(tree)

    # Get non-redundant distances (lower triangle) from each matrix:
    distances_from_matrix <- matrix[lower.tri(matrix,diag = FALSE)]
    distances_from_tree <- matrix_from_tree[lower.tri(matrix_from_tree,diag = FALSE)]

    # Find the total relative distance between the two:
    relative_distances <- mapply(relative_difference, distances_from_matrix, distances_from_tree)
    cat(sum(relative_distances))
    cat("\n")
}


load_phylip_distances <- function(filename) {
    matrix <- read.table(filename, skip = 1, row.names = 1)
    colnames(matrix) <- rownames(matrix)
    return(as.matrix(matrix))
}


relative_difference <- function(a, b) {
    if (a == 0.0 && b == 0.0) {
        return(0.0)
    }
    if (a > b) {
        rel_diff <- (a / b) - 1
    } else {
        rel_diff <- (b / a) - 1
    }
    # To deal with infinities (when one value is 0), we put the result through a function with an
    # asymptote of one:
    return(1 - (1 / (rel_diff+1)))
}


if(!interactive()) {
    main()
}
