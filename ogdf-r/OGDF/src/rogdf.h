// -*- mode: C++ -*-
 
// OGDF R package
// Copyright (C) 2014  Gabor Csardi <csardi.gabor@gmail.com>
// 334 Harvard street, Cambridge, MA 02139 USA

// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
// 02110-1301 USA

#ifndef ROGDF_H
#define ROGDF_H

#include <Rcpp.h>
using namespace Rcpp;

#include <ogdf/misclayout/CircularLayout.h>
#include <ogdf/energybased/FMMMLayout.h>
using namespace ogdf;

#include "convert.h"

IntegerVector rogdf_version();

NumericMatrix rogdf_circular_layout(GraphAttributes graph, 
	 double minDistCircle, double minDistLevel,
	 double minDistSibling, double minDistCC, double pageRatio);

NumericMatrix rogdf_tree_layout(GraphAttributes graph,
	 double siblingDistance, double subtreeDistance,
	 double levelDistance, double treeDistance,
	 bool orthogonalLayout, Orientation orientation,
	 TreeLayout::RootSelectionType selectRoot);

#endif
