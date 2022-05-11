#!/usr/bin/env Rscript

# Copyright 2022 Ryan Wick (rrwick@gmail.com)
# https://github.com/rrwick/Verticall
#
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

    # Load in the PHYLIP distance matrix:
    matrix <- read.table(args[1], skip = 1, row.names = 1)
    samples <- rownames(matrix)
    colnames(matrix) <- samples
    matrix <- as.matrix(matrix)

    # Build a neighbour-joining tree:
    tree <- nj(matrix)
    tree$edge.length <- pmax(tree$edge.length, 0.0)  # set any negative branch lengths to zero

    # Get all pairwise distances from the tree:
    matrix_from_tree <- cophenetic.phylo(tree)

    # Get non-redundant distances (lower triangle) from each matrix:
    distances <- matrix[lower.tri(matrix,diag = FALSE)]
    distances_from_tree <- matrix_from_tree[lower.tri(matrix_from_tree,diag = FALSE)]

    # Find the total relative distance between the two:
    relative_distances <- mapply(relative_distance, distances, distances_from_tree)
    cat(sum(relative_distances))
    cat("\n")
}

relative_distance <- function(a, b) {
    if (a == 0.0 && b == 0.0) {
        return(0.0)
    } else if (a == 0.0) {
        return(1.0)
    } else if (b == 0.0) {
        return(1.0)
    } else if (a > b) {
        return(min(1, (a / b) - 1))
    } else {
        return(min(1, (b / a) - 1))
    }
}


if(!interactive()) {
    main()
}
