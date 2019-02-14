#include "layout.hh"

using namespace ogdf;

template <class T>
static MultilevelBuilder *GetFactoredAdjustedMerger(double factor = 2.0, int adjustment = 0)
{
	T *merger = new T();
	merger->setFactor(factor);
	merger->setEdgeLengthAdjustment(adjustment);
	return merger;
}

template <class T>
static MultilevelBuilder *GetAdjustedMerger(int adjustment = 0)
{
	T *merger = new T();
	merger->setEdgeLengthAdjustment(adjustment);
	return merger;
}

static InitialPlacer *GetBarycenterPlacer()
{
	BarycenterPlacer *placer = new BarycenterPlacer();
	placer->weightedPositionPriority(true);
	return placer;
}

static InitialPlacer *GetSolarPlacer()
{
	SolarPlacer *placer = new SolarPlacer();
	return placer;
}

static InitialPlacer *GetCirclePlacer()
{
	CirclePlacer *placer = new CirclePlacer();
	return placer;
}

static InitialPlacer *GetMedianPlacer()
{
	MedianPlacer *placer = new MedianPlacer();
	return placer;
}

static InitialPlacer *GetRandomPlacer()
{
	RandomPlacer *placer = new RandomPlacer();
	return placer;
}

static InitialPlacer *GetZeroPlacer()
{
	ZeroPlacer *placer = new ZeroPlacer();
	return placer;
}

std::tuple<std::vector<double>, std::vector<double>>
Layout(uint32_t vertex_count, const std::vector<uint32_t> &from, const std::vector<uint32_t> &to,
	   LayoutConfiguration config, const std::vector<float> &weight)
{
	bool weighted = false;

	// from and to vectors have to have equal length
	if (from.size() != to.size())
	{
		throw std::length_error("From and to vectors are not of equal size");
	}

	// Also throw error if weight is set but not matching in length
	if (weight.size() > 0 && weight.size() != from.size())
	{
		throw std::length_error("Weight vector length does not match edge count");
	}
	else if (weight.size() == from.size())
	{
		weighted = true;
	}

	EdgeWeightedGraph<double> g;
	GraphAttributes ga(g);
	std::map<uint32_t, node> index_to_node;

	for (uint32_t i = 0; i < vertex_count; i++)
	{
		index_to_node[i] = g.newNode();
	}

	for (std::vector<uint32_t>::size_type i = 0; i != from.size(); i++)
	{
		if (weighted)
		{
			g.newEdge(index_to_node[from[i]], index_to_node[to[i]], weight[i]);
		}
		else 
		{
			g.newEdge(index_to_node[from[i]], index_to_node[to[i]], 1.0);
		}
	}

	for (node v : g.nodes)
	{
		ga.width(v) = ga.height(v) = 1.0;
	}

	// Check for isolated nodes. If there are isolated nodes,
	// call placement step later on.
	NodeArray<int> connected_components(g);
	List<node> isolated_nodes;
	int n_connected_components = connectedComponents(g, connected_components, &isolated_nodes);

	// Starting the layout
	MultilevelGraph mlg(ga);

	// The FastMultipoleEmbedder is used for the single level layout.
	FastMultipoleEmbedder *fme = new FastMultipoleEmbedder();
	fme->setNumIterations(config.fme_iterations);
	fme->setRandomize(config.fme_randomize);
	fme->setNumberOfThreads(config.fme_threads);

	// To minimize dispersion of the graph when more nodes are added, a
	// ScalingLayout can be used to scale up the graph on each level.
	ScalingLayout *sl = new ScalingLayout();
	sl->setLayoutRepeats(config.sl_repeats);
	sl->setSecondaryLayout(fme);

	// Used for the placement.
	InitialPlacer *placer = GetBarycenterPlacer();
	switch (config.placer)
	{
	case Placer::Barycenter:
		placer = GetBarycenterPlacer();
	case Placer::Circle:
		placer = GetCirclePlacer();
	case Placer::Median:
		placer = GetMedianPlacer();
	case Placer::Random:
		placer = GetRandomPlacer();
	case Placer::Solar:
		placer = GetSolarPlacer();
	case Placer::Zero:
		placer = GetZeroPlacer();
	}

	// Used for the coarsening phase.
	MultilevelBuilder *merger = GetFactoredAdjustedMerger<EdgeCoverMerger>();
	switch (config.merger)
	{
	case Merger::EdgeCover:
		merger = GetFactoredAdjustedMerger<EdgeCoverMerger>(config.merger_factor, config.merger_adjustment);
	case Merger::LocalBiconnected:
		merger = GetFactoredAdjustedMerger<LocalBiconnectedMerger>(config.merger_factor, config.merger_adjustment);
	case Merger::Solar:
		merger = GetAdjustedMerger<SolarMerger>(config.merger_adjustment);
	case Merger::IndependentSet:
		merger = GetAdjustedMerger<IndependentSetMerger>(config.merger_adjustment);
	}

	// Get the scaling type. As I do not want to expose any OGDF to Python,
	// there is this intermediate step.
	ScalingLayout::ScalingType scaling_type = ScalingLayout::ScalingType::RelativeToAvgLength;
	switch (config.sl_scaling_type)
	{
	case ScalingType::Absolute:
		scaling_type = ScalingLayout::ScalingType::Absolute;
	case ScalingType::RelativeToAvgLength:
		scaling_type = ScalingLayout::ScalingType::RelativeToAvgLength;
	case ScalingType::RelativeToDesiredLength:
		scaling_type = ScalingLayout::ScalingType::RelativeToDesiredLength;
	case ScalingType::RelativeToDrawing:
		scaling_type = ScalingLayout::ScalingType::RelativeToDrawing;
	}

	// Postprocessing is applied at each level after the single level layout.
	// In this example a FastMultipoleEmbedder with zero iterations is used for postprocessing.
	sl->setExtraScalingSteps(config.sl_extra_scaling_steps);
	sl->setScalingType(scaling_type);
	sl->setScaling(config.sl_scaling_x, config.sl_scaling_y);

	// Then the ModularMultilevelMixer is created.
	ModularMultilevelMixer *mmm = new ModularMultilevelMixer;
	mmm->setLayoutRepeats(1);
	// The single level layout, the placer and the merger are set.
	mmm->setLevelLayoutModule(sl);
	mmm->setInitialPlacer(placer);
	mmm->setMultilevelBuilder(merger);

	if (n_connected_components > 1)
	{
		// Since energybased algorithms are not doing well for disconnected
		// graphs, the ComponentSplitterLayout is used to split the graph and
		// computation is done separately for each connected component.
		ComponentSplitterLayout *csl = new ComponentSplitterLayout;
		// The TileToRowsPacker merges these connected components after computation.
		TileToRowsCCPacker *ttrccp = new TileToRowsCCPacker;
		csl->setPacker(ttrccp);
		csl->setLayoutModule(mmm);

		// At last the PreprocessorLayout removes double edges and loops.
		PreprocessorLayout ppl;
		ppl.setLayoutModule(csl);
		ppl.setRandomizePositions(true);

		ppl.call(mlg);
	}
	else
	{
		mmm->call(mlg);
	}

	mlg.exportAttributes(ga);

	std::vector<double> x(vertex_count);
	std::vector<double> y(vertex_count);

	int i = 0;
	for (node v : g.nodes)
	{
		x[i] = ga.x(v);
		y[i] = ga.y(v);
		i++;
	}

	return std::make_tuple(x, y);
}