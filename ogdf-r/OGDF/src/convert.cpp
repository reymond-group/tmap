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

#include "convert.h"

#define USE_RINTERNALS
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

#undef AS_INTEGER
#define AS_INTEGER(x) Rf_coerceVector(x, INTSXP)
#undef PROTECT
#define PROTECT(x) Rf_protect(x)

namespace Rcpp {
  template <> ogdf::GraphAttributes as(SEXP igraph) {
    SEXP from, to;
    ogdf::Graph *G=new ogdf::Graph();
    ogdf::GraphAttributes GA(*G, ogdf::GraphAttributes::nodeGraphics | 
			     ogdf::GraphAttributes::edgeGraphics);
    int no_nodes=INTEGER(AS_INTEGER(VECTOR_ELT(igraph, 0)))[0];
    ogdf::Array<ogdf::node> nodemap(no_nodes);

    // TODO: check argument!

    PROTECT(from=AS_INTEGER(VECTOR_ELT(igraph, 2)));
    PROTECT(to=AS_INTEGER(VECTOR_ELT(igraph, 3)));

    int *ifrom=INTEGER(from);
    int *ito=INTEGER(to);
    int no_edges=Rf_length(from);
    
    for (int i=0; i<no_nodes; i++) { nodemap[i] = G->newNode(i); }
    for (int e=0; e<no_edges; e++) { 
      ogdf::node u=nodemap[ ifrom[e] ];
      ogdf::node v=nodemap[ ito[e] ];
      G->newEdge(u, v); 
    }
    
    UNPROTECT(2);
    return GA;
  }

  template <> ogdf::Orientation as(SEXP orientation) {
    CharacterVector ori=as<CharacterVector>(orientation);
    if (!strcmp(ori[0], "topToBottom")) {
      return ogdf::topToBottom;
    } else if (!strcmp(ori[0], "bottomToTop")) {
      return ogdf::bottomToTop;
    } else if (!strcmp(ori[0], "leftToRight")) {
      return ogdf::leftToRight;
    } else if (!strcmp(ori[0], "rightToLeft")) {
      return ogdf::rightToLeft;
    } else {
      // FIXME: error
      return ogdf::topToBottom;
    }
  }

  template <> ogdf::TreeLayout::RootSelectionType as(SEXP selectRoot) {
    ogdf::TreeLayout::RootSelectionType result;
    CharacterVector rs=as<CharacterVector>(selectRoot);
    if (!strcmp(rs[0], "rootIsSource")) {
      return ogdf::TreeLayout::rootIsSource;
    } else if (!strcmp(rs[0], "rootIsSink")) {
      return ogdf::TreeLayout::rootIsSink;
    } else if (!strcmp(rs[0], "rootByCoord")) {
      return ogdf::TreeLayout::rootByCoord;
    } else {
      // FIXME: error
      return ogdf::TreeLayout::rootIsSource;
    }
  }
  
  template <> FMMMLayout::PageFormatType as(SEXP pageFormat) {
    CharacterVector pf=as<CharacterVector>(pageFormat);
    if (!strcmp(pf[0], "Portrait")) {
      return ogdf::FMMMLayout::pfPortrait;      
    } else if (!strcmp(pf[0], "Landscape")) {
      return ogdf::FMMMLayout::pfLandscape;
    } else if (!strcmp(pf[0], "Square")) {
      return ogdf::FMMMLayout::pfSquare;
    } else {
      // FIXME: error
      return ogdf::FMMMLayout::pfSquare;
    }
  }

  template <> FMMMLayout::EdgeLengthMeasurement 
              as(SEXP edgeLengthMeasurement) {
    CharacterVector elm=as<CharacterVector>(edgeLengthMeasurement);
    if (!strcmp(elm[0], "midpoint")) {
      return ogdf::FMMMLayout::elmMidpoint;
    } else if (!strcmp(elm[0], "boundingcircle")) {
      return ogdf::FMMMLayout::elmBoundingCircle;
    } else {
      // FIXME: error
      return ogdf::FMMMLayout::elmMidpoint;
    }
  }

  template <> FMMMLayout::AllowedPositions as(SEXP allowedPositions) {
    CharacterVector ap=as<CharacterVector>(allowedPositions);
    if (!strcmp(ap[0], "all")) {
      return ogdf::FMMMLayout::apAll;
    } else if (!strcmp(ap[0], "integer")) {
      return ogdf::FMMMLayout::apInteger;      
    } else if (!strcmp(ap[0], "exponent")) {
      return ogdf::FMMMLayout::apExponent;      
    } else {
      // FIXME: error
      return ogdf::FMMMLayout::apAll;
    }
  }

  template <> FMMMLayout::TipOver as(SEXP tipOverCCs) {
    CharacterVector to=as<CharacterVector>(tipOverCCs);
    if (!strcmp(to[0], "none")) {
      return ogdf::FMMMLayout::toNone;
    } else if (!strcmp(to[0], "nogrowingrow")) {
      return ogdf::FMMMLayout::toNoGrowingRow;
    } else if (!strcmp(to[0], "always")) {
      return ogdf::FMMMLayout::toAlways;
    } else {
      // FIXME: error
      return ogdf::FMMMLayout::toNone;
    }
  }

  template <> FMMMLayout::PreSort as(SEXP presortCCs) {
    CharacterVector ps=as<CharacterVector>(presortCCs);
    if (!strcmp(ps[0], "none")) {
      return ogdf::FMMMLayout::psNone;
    } else if (!strcmp(ps[0], "decreasingheight")) {
      return ogdf::FMMMLayout::psDecreasingHeight;
    } else if (!strcmp(ps[0], "decreasingwidth")) {
      return ogdf::FMMMLayout::psDecreasingWidth;
    } else {
      // FIXME: error
      return ogdf::FMMMLayout::psNone;
    }
  }

  template <> FMMMLayout::GalaxyChoice as(SEXP galaxyChoice) {
    CharacterVector gc=as<CharacterVector>(galaxyChoice);
    if (!strcmp(gc[0], "uniformprob")) {
      return ogdf::FMMMLayout::gcUniformProb;
    } else if (!strcmp(gc[0], "nonuniformproblowermass")) {
      return ogdf::FMMMLayout::gcNonUniformProbLowerMass;
    } else if (!strcmp(gc[0], "nonuniformprobhighermass")) {
      return ogdf::FMMMLayout::gcNonUniformProbHigherMass;
    } else {
      // FIXME: error
      return ogdf::FMMMLayout::gcUniformProb;
    }
  }

  template <> FMMMLayout::MaxIterChange as(SEXP maxIterChange) {
    CharacterVector mic=as<CharacterVector>(maxIterChange);
    if (!strcmp(mic[0], "constant")) {
      return ogdf::FMMMLayout::micConstant;
    } else if (!strcmp(mic[0], "linearlydecreasing")) {
      return ogdf::FMMMLayout::micLinearlyDecreasing;
    } else if (!strcmp(mic[0], "rapidlydecreasing")) {
      return ogdf::FMMMLayout::micRapidlyDecreasing;
    } else {
      // FIXME: error
      return ogdf::FMMMLayout::micConstant;
    }    
  }

  template <> FMMMLayout::InitialPlacementMult 
              as(SEXP initialPlacementMult) {
    CharacterVector ipm=as<CharacterVector>(initialPlacementMult);
    if (!strcmp(ipm[0], "simple")) {
      return ogdf::FMMMLayout::ipmSimple;
    } else if (!strcmp(ipm[0], "advanced")) {
      return ogdf::FMMMLayout::ipmAdvanced;
    } else {
      // FIXME: error
      return ogdf::FMMMLayout::ipmSimple;
    }
  }

  template <> FMMMLayout::ForceModel as(SEXP forceModel) {
    CharacterVector fm=as<CharacterVector>(forceModel);
    if (!strcmp(fm[0], "fruchtermanreingold")) {
      return ogdf::FMMMLayout::fmFruchtermanReingold;
    } else if (!strcmp(fm[0], "eades")) {
      return ogdf::FMMMLayout::fmEades;
    } else if (!strcmp(fm[0], "new")) {
      return ogdf::FMMMLayout::fmNew;
    } else {
      // FIXME: error
      return ogdf::FMMMLayout::fmFruchtermanReingold;
    }    
  }

  template <> FMMMLayout::RepulsiveForcesMethod 
              as(SEXP repulsiveForcesCalculation) {
    CharacterVector rfc=as<CharacterVector>(repulsiveForcesCalculation);
    if (!strcmp(rfc[0], "exact")) {
      return ogdf::FMMMLayout::rfcExact;
    } else if (!strcmp(rfc[0], "gridapproximation")) {
      return ogdf::FMMMLayout::rfcGridApproximation;
    } else if (!strcmp(rfc[0], "nmm")) {
      return ogdf::FMMMLayout::rfcNMM;
    } else {
      // FIXME: error
      return ogdf::FMMMLayout::rfcExact;
    }    
  }

  template <> FMMMLayout::StopCriterion as(SEXP stopCriterion) {
    CharacterVector sc=as<CharacterVector>(stopCriterion);
    if (!strcmp(sc[0], "fixediterations")) {
      return ogdf::FMMMLayout::scFixedIterations;
    } else if (!strcmp(sc[0], "threshold")) {
      return ogdf::FMMMLayout::scThreshold;
    } else if (!strcmp(sc[0], "fixediterationsorthreshold")) {
      return ogdf::FMMMLayout::scFixedIterationsOrThreshold;
    } else {
      // FIXME: error
      return ogdf::FMMMLayout::scFixedIterations;
    }    
  }

  template <> FMMMLayout::InitialPlacementForces 
              as(SEXP initialPlacementForces) {
    CharacterVector ipf=as<CharacterVector>(initialPlacementForces);
    if (!strcmp(ipf[0], "uniformgrid")) {
      return ogdf::FMMMLayout::ipfUniformGrid;
    } else if (!strcmp(ipf[0], "randomtime")) {
      return ogdf::FMMMLayout::ipfRandomTime;
    } else if (!strcmp(ipf[0], "randomranditernr")) {
      return ogdf::FMMMLayout::ipfRandomRandIterNr;
    } else if (!strcmp(ipf[0], "keeppositions")) {
      return ogdf::FMMMLayout::ipfKeepPositions;
    } else {
      // FIXME: error
      return ogdf::FMMMLayout::ipfUniformGrid;
    }    
  }

  template <> FMMMLayout::ReducedTreeConstruction
              as(SEXP nmTreeConstruction) {
    CharacterVector rtc=as<CharacterVector>(nmTreeConstruction);
    if (!strcmp(rtc[0], "pathbypath")) {
      return ogdf::FMMMLayout::rtcPathByPath;
    } else if (!strcmp(rtc[0], "subtreebysubtree")) {
      return ogdf::FMMMLayout::rtcSubtreeBySubtree;
    } else {
      // FIXME: error
      return ogdf::FMMMLayout::rtcPathByPath;
    }    
  }

  template <> FMMMLayout::SmallestCellFinding as(SEXP nmSmallCell) {
    CharacterVector scf=as<CharacterVector>(nmSmallCell);
    if (!strcmp(scf[0], "iteratively")) {
      return ogdf::FMMMLayout::scfIteratively;
    } else if (!strcmp(scf[0], "aluru")) {
      return ogdf::FMMMLayout::scfAluru;
    } else {
      // FIXME: error
      return ogdf::FMMMLayout::scfIteratively;
    }    
  }
  
  template <> ROGDF_FMMM_Options as(SEXP options) {
    ROGDF_FMMM_Options result;
    List opt=as<List>(options);

    result.useHighLevelOptions = as<bool>(opt["useHighLevelOptions"]);
    result.singleLevel = as<bool>(opt["singleLevel"]);
    result.pageFormat = as<FMMMLayout::PageFormatType>(opt["pageFormat"]);
    result.unitEdgeLength = as<double>(opt["unitEdgeLength"]);
    result.newInitialPlacement = as<bool>(opt["newInitialPlacement"]);

    result.randSeed = as<int>(opt["randSeed"]);
    result.edgeLengthMeasurement = 
      as<FMMMLayout::EdgeLengthMeasurement>(opt["edgeLengthMeasurement"]);
    result.allowedPositions = 
      as<FMMMLayout::AllowedPositions>(opt["allowedPositions"]);
    result.maxIntPosExponent = as<int>(opt["maxIntPosExponent"]);

    result.pageRatio = as<double>(opt["pageRatio"]);
    result.stepsForRotatingComponents = 
      as<int>(opt["stepsForRotatingComponents"]);
    result.tipOverCCs = as<FMMMLayout::TipOver>(opt["tipOverCCs"]);
    result.minDistCC = as<double>(opt["minDistCC"]);
    result.presortCCs = as<FMMMLayout::PreSort>(opt["presortCCs"]);

    result.minGraphSize = as<int>(opt["minGraphSize"]);
    result.galaxyChoice = as<FMMMLayout::GalaxyChoice>(opt["galaxyChoice"]);
    result.randomTries = as<int>(opt["randomTries"]);
    result.maxIterChange = 
      as<FMMMLayout::MaxIterChange>(opt["maxIterChange"]);
    result.maxIterFactor = as<int>(opt["maxIterFactor"]);
    result.initialPlacementMult = 
      as<FMMMLayout::InitialPlacementMult>(opt["initialPlacementMult"]);

    result.forceModel = as<FMMMLayout::ForceModel>(opt["forceModel"]);
    result.springStrength = as<double>(opt["springStrength"]);
    result.repForcesStrength = as<double>(opt["repForcesStrength"]);
    result.repulsiveForcesCalculation = as<FMMMLayout::RepulsiveForcesMethod>
      (opt["repulsiveForcesCalculation"]);
    result.stopCriterion = 
      as<FMMMLayout::StopCriterion>(opt["stopCriterion"]);
    result.threshold = as<double>(opt["threshold"]);
    result.fixedIterations = as<int>(opt["fixedIterations"]);
    result.forceScalingFactor = as<double>(opt["forceScalingFactor"]);
    result.coolTemperature = as<bool>(opt["coolTemperature"]);
    result.coolValue = as<double>(opt["coolValue"]);
    result.initialPlacementForces = as<FMMMLayout::InitialPlacementForces>
      (opt["initialPlacementForces"]);
    
    result.resizeDrawing = as<bool>(opt["resizeDrawing"]);
    result.resizingScalar = as<double>(opt["resizingScalar"]);
    result.fineTuningIterations = as<int>(opt["fineTuningIterations"]);
    result.fineTuneScalar = as<double>(opt["fineTuneScalar"]);
    result.adjustPostRepStrengthDynamically = as<bool>
      (opt["adjustPostRepStrengthDynamically"]);
    result.postSpringStrength = as<double>(opt["postSpringStrength"]);
    result.postStrengthOfRepForces = as<double>
      (opt["postStrengthOfRepForces"]);
    
    result.frGridQuotient = as<int>(opt["frGridQuotient"]);
    result.nmTreeConstruction = as<FMMMLayout::ReducedTreeConstruction>
      (opt["nmTreeConstruction"]);
    result.nmSmallCell = as<FMMMLayout::SmallestCellFinding>
      (opt["nmSmallCell"]);
    result.nmParticlesInLeaves = as<int>(opt["nmParticlesInLeaves"]);
    result.nmPrecision = as<int>(opt["nmPrecision"]);

    return result;
  }

} // namespace Rcpp
