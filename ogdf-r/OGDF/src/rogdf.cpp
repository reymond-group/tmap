#include "rogdf.h"

// [[Rcpp::export]]
IntegerVector rogdf_version() {
  return IntegerVector::create(_["major"] = 2012,
			       _["minor"] = 7,
			       _["patch"] = 0);
}

NumericMatrix rogdf_get_layout(GraphAttributes &GA) {
  Graph G=GA.constGraph();
  NumericMatrix coords(G.numberOfNodes(), 2);
  int i; node v;
  for(i=0, v=G.firstNode(); v; i++, v=v->succ())  {
    coords(i, 0) = GA.x(v);
    coords(i, 1) = GA.y(v);
  }
  return coords;
}

// [[Rcpp::export("L_circular")]]  
NumericMatrix rogdf_circular_layout(GraphAttributes graph, 
	 double minDistCircle=20, double minDistLevel=20, 
	 double minDistSibling=10, double minDistCC=20, double pageRatio=1) {
  
  CircularLayout layout;
  layout.minDistCircle(minDistCircle);
  layout.minDistLevel(minDistLevel);
  layout.minDistSibling(minDistSibling);
  layout.minDistCC(minDistCC);
  layout.pageRatio(pageRatio);
  layout.call(graph);
  return rogdf_get_layout(graph);
}

// [[Rcpp::export("L_tree")]]
NumericMatrix rogdf_tree_layout(GraphAttributes graph,
	 double siblingDistance=20, double subtreeDistance=20,
	 double levelDistance=50, double treeDistance=50,
	 bool orthogonalLayout=false,
	 Orientation orientation=topToBottom,
	 TreeLayout::RootSelectionType selectRoot=TreeLayout::rootIsSource) {

  TreeLayout layout;
  layout.siblingDistance(siblingDistance);
  layout.subtreeDistance(subtreeDistance);
  layout.levelDistance(levelDistance);
  layout.treeDistance(treeDistance);
  layout.orthogonalLayout(orthogonalLayout);
  layout.orientation(orientation);
  layout.rootSelection(selectRoot);
  layout.call(graph);
  return rogdf_get_layout(graph);
}

// [[Rcpp::export("L_fmmm")]]
NumericMatrix rogdf_fmmm_layout(GraphAttributes graph,
	ROGDF_FMMM_Options options=ROGDF_FMMM_Options() 
				/* @R L_fmmm_options */) {
  FMMMLayout layout;

  // High-level options
  layout.useHighLevelOptions(options.useHighLevelOptions);
  layout.setSingleLevel(options.singleLevel);
  layout.pageFormat(options.pageFormat);
  layout.unitEdgeLength(options.unitEdgeLength);
  layout.newInitialPlacement(options.newInitialPlacement);
  layout.qualityVersusSpeed(options.qualityVersusSpeed);
  
  if (! options.useHighLevelOptions) {
    // General low level options
    layout.randSeed(options.randSeed);
    layout.edgeLengthMeasurement(options.edgeLengthMeasurement);
    layout.allowedPositions(options.allowedPositions);
    layout.maxIntPosExponent(options.maxIntPosExponent);
    
    // Options for the divide et impera step
    layout.pageRatio(options.pageRatio);
    layout.stepsForRotatingComponents(options.stepsForRotatingComponents);
    layout.tipOverCCs(options.tipOverCCs);
    layout.minDistCC(options.minDistCC);
    layout.presortCCs(options.presortCCs);
    
    // Options for the multilevel step
    layout.minGraphSize(options.minGraphSize);
    layout.galaxyChoice(options.galaxyChoice);
    layout.randomTries(options.randomTries);
    layout.maxIterChange(options.maxIterChange);
    layout.maxIterFactor(options.maxIterFactor);
    layout.initialPlacementMult(options.initialPlacementMult);
    
    // Options for the force calculation step
    layout.forceModel(options.forceModel);
    layout.springStrength(options.springStrength);
    layout.repForcesStrength(options.repForcesStrength);
    layout.repulsiveForcesCalculation(options.repulsiveForcesCalculation);
    layout.stopCriterion(options.stopCriterion);
    layout.threshold(options.threshold);
    layout.fixedIterations(options.fixedIterations);
    layout.forceScalingFactor(options.forceScalingFactor);
    layout.coolTemperature(options.coolTemperature);
    layout.coolValue(options.coolValue);
    layout.initialPlacementForces(options.initialPlacementForces);
    
    // Options for the postprocessing step
    layout.resizeDrawing(options.resizeDrawing);
    layout.resizingScalar(options.resizingScalar);
    layout.fineTuningIterations(options.fineTuningIterations);
    layout.fineTuneScalar(options.fineTuneScalar);
    layout.adjustPostRepStrengthDynamically
      (options.adjustPostRepStrengthDynamically);
    layout.postSpringStrength(options.postSpringStrength);
    layout.postStrengthOfRepForces(options.postStrengthOfRepForces);
    
    // Options for repulsice force approximation methods
    layout.frGridQuotient(options.frGridQuotient);
    layout.nmTreeConstruction(options.nmTreeConstruction);
    layout.nmSmallCell(options.nmSmallCell);
    layout.nmParticlesInLeaves(options.nmParticlesInLeaves);
    layout.nmPrecision(options.nmPrecision);        
  }

  layout.call(graph);
  return rogdf_get_layout(graph);
}
