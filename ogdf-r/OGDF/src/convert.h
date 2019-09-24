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

#ifndef ROGDF_CONVERT_H
#define ROGDF_CONVERT_H

#include <RcppCommon.h>
#include <ogdf/basic/GraphAttributes.h>
#include <ogdf/tree/TreeLayout.h>

#include "rogdf.h"

struct ROGDF_FMMM_Options {

  // High-level options
  bool useHighLevelOptions;
  bool singleLevel;
  FMMMLayout::PageFormatType pageFormat;
  double unitEdgeLength;
  bool newInitialPlacement;
  FMMMLayout::QualityVsSpeed qualityVersusSpeed;

  // General low level options
  int randSeed;
  FMMMLayout::EdgeLengthMeasurement edgeLengthMeasurement;
  FMMMLayout::AllowedPositions allowedPositions;
  int maxIntPosExponent;

  // Options for the divide et impera step
  double pageRatio;
  int stepsForRotatingComponents;
  FMMMLayout::TipOver tipOverCCs;
  double minDistCC;
  FMMMLayout::PreSort presortCCs;

  // Options for the multilevel step
  int minGraphSize;
  FMMMLayout::GalaxyChoice galaxyChoice;
  int randomTries;
  FMMMLayout::MaxIterChange maxIterChange;
  int maxIterFactor;
  FMMMLayout::InitialPlacementMult initialPlacementMult;
  
  // Options for the force calculation step
  FMMMLayout::ForceModel forceModel;
  double springStrength;
  double repForcesStrength;
  FMMMLayout::RepulsiveForcesMethod repulsiveForcesCalculation;
  FMMMLayout::StopCriterion stopCriterion;
  double threshold;
  int fixedIterations;
  double forceScalingFactor;
  bool coolTemperature;
  double coolValue;
  FMMMLayout::InitialPlacementForces initialPlacementForces;

  // Options for the postprocessing step
  bool resizeDrawing;
  double resizingScalar;
  int fineTuningIterations;
  double fineTuneScalar;
  bool adjustPostRepStrengthDynamically;
  double postSpringStrength;
  double postStrengthOfRepForces;

  // Options for repulsice force approximation methods
  int frGridQuotient;
  FMMMLayout::ReducedTreeConstruction nmTreeConstruction;
  FMMMLayout::SmallestCellFinding nmSmallCell;
  int nmParticlesInLeaves;
  int nmPrecision;
};

namespace Rcpp {
  template <> ogdf::GraphAttributes as(SEXP igraph);
  template <> ogdf::Orientation as(SEXP orientation);
  template <> ogdf::TreeLayout::RootSelectionType as(SEXP selectRoot);
  template <> FMMMLayout::PageFormatType as(SEXP pageFormat);
  template <> FMMMLayout::EdgeLengthMeasurement 
              as(SEXP edgeLengthMeasurement);
  template <> FMMMLayout::AllowedPositions as(SEXP allowedPositions);
  template <> FMMMLayout::TipOver as(SEXP tipOverCCs);
  template <> FMMMLayout::PreSort as(SEXP presortCCs);
  template <> FMMMLayout::GalaxyChoice as(SEXP galaxyChoice);
  template <> FMMMLayout::MaxIterChange as(SEXP maxIterChange);
  template <> FMMMLayout::InitialPlacementMult as(SEXP initialPlacementMult);
  template <> FMMMLayout::ForceModel as(SEXP forceModel);
  template <> FMMMLayout::RepulsiveForcesMethod 
              as(SEXP repulsiveForcesCalculation);
  template <> FMMMLayout::StopCriterion as(SEXP stopCriterion);
  template <> FMMMLayout::InitialPlacementForces 
              as(SEXP initialPlacementForces);
  template <> FMMMLayout::ReducedTreeConstruction
              as(SEXP nmTreeConstruction);
  template <> FMMMLayout::SmallestCellFinding as(SEXP nmSmallCell);
  template <> ROGDF_FMMM_Options as(SEXP options);
}

#include <Rcpp.h>

#endif

