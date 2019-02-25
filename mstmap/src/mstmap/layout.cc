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

std::vector<std::vector<uint32_t>> GetTreesFromForest(const Graph &g)
{
	NodeArray<int> connected_component_ids(g);
	List<node> isolated_nodes;
	int n_connected_components = connectedComponents(g, connected_component_ids, &isolated_nodes);

	std::vector<std::vector<uint32_t>> connected_components(n_connected_components);

	for (int i = 0; i < n_connected_components; i++)
		connected_components[i] = std::vector<uint32_t>();

	uint32_t i = 0;
	for (auto id : connected_component_ids)
		connected_components[id].emplace_back(i++);

	std::sort(connected_components.begin(), connected_components.end(), 
			  [](const std::vector<uint32_t> & a, const std::vector<uint32_t> & b){ return a.size() > b.size(); });

	return connected_components;
}

void ConnectGraph(Graph &g, std::map<uint32_t, node> &index_to_node, LSHForest &lsh_forest)
{
	auto trees = GetTreesFromForest(g);

	for (size_t i = 1; i < trees.size(); i++)
	{
		for (uint32_t v : trees[i])
		{
			auto nns = lsh_forest.QueryLinearScanExcludeById(v, 1, trees[i], 10);
			
			if (nns.size() == 0)
				continue;

			g.newEdge(index_to_node[v], index_to_node[std::get<1>(nns[0])], std::get<0>(nns[0]));
		}

		trees = GetTreesFromForest(g);
	}
}

std::tuple<std::vector<float>, std::vector<float>>
LayoutFromLSHForest(LSHForest &lsh_forest, LayoutConfiguration config, bool create_mst)
{
	uint32_t vertex_count = lsh_forest.size();
	auto edges = lsh_forest.GetKNNGraph(config.k);

	EdgeWeightedGraph<float> g;
	
	std::map<uint32_t, node> index_to_node;

	for (uint32_t i = 0; i < vertex_count; i++)
	{
		index_to_node[i] = g.newNode();
	}

	for (std::vector<uint32_t>::size_type i = 0; i != edges.size(); i++)
	{
		g.newEdge(index_to_node[std::get<0>(edges[i])], index_to_node[std::get<1>(edges[i])], std::get<2>(edges[i]));
	}

	if (create_mst)
	{
		ConnectGraph(g, index_to_node, lsh_forest);

		float weight = ogdf::makeMinimumSpanningTree(g, g.edgeWeights());
		std::cout << "Weight: " << weight << std::endl;

		auto connected_components = GetTreesFromForest(g);
	}

	return LayoutInternal(g, vertex_count, config);
}

std::tuple<std::vector<float>, std::vector<float>>
LayoutFromEdgeList(uint32_t vertex_count, const std::vector<std::tuple<uint32_t, uint32_t, float>> &edges,
       LayoutConfiguration config, bool create_mst)
{
	EdgeWeightedGraph<float> g;
	
	std::map<uint32_t, node> index_to_node;

	for (uint32_t i = 0; i < vertex_count; i++)
	{
		index_to_node[i] = g.newNode();
	}

	for (std::vector<uint32_t>::size_type i = 0; i != edges.size(); i++)
	{
		g.newEdge(index_to_node[std::get<0>(edges[i])], index_to_node[std::get<1>(edges[i])], std::get<2>(edges[i]));
	}

	if (create_mst)
	{
		float weight = ogdf::makeMinimumSpanningTree(g, g.edgeWeights());
		std::cout << "Weight: " << weight << std::endl;

		auto connected_components = GetTreesFromForest(g);
		
		for (auto cc : connected_components)
			std::cout << cc.size() << std::endl;
	}

	return LayoutInternal(g, vertex_count, config);
}

std::tuple<std::vector<float>, std::vector<float>>
LayoutInternal(Graph &g, uint32_t vertex_count, LayoutConfiguration config)
{
	GraphAttributes ga(g);

	for (node v : g.nodes)
	{
		ga.width(v) = ga.height(v) = 1.0;
	}

	// Check for isolated nodes. If there are isolated nodes,
	// call placement step later on.
	NodeArray<int> connected_components(g);
	List<node> isolated_nodes;
	int n_connected_components = connectedComponents(g, connected_components, &isolated_nodes);
	std::cout << "Connected components: " << n_connected_components << std::endl;
	std::cout << "Isolated nodes: " << isolated_nodes.size() << std::endl;

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
		break;
	case Placer::Circle:
		placer = GetCirclePlacer();
		break;
	case Placer::Median:
		placer = GetMedianPlacer();
		break;
	case Placer::Random:
		placer = GetRandomPlacer();
		break;
	case Placer::Solar:
		placer = GetSolarPlacer();
		break;
	case Placer::Zero:
		placer = GetZeroPlacer();
		break;
	}

	// Used for the coarsening phase.
	MultilevelBuilder *merger = GetFactoredAdjustedMerger<EdgeCoverMerger>();
	switch (config.merger)
	{
	case Merger::EdgeCover:
		merger = GetFactoredAdjustedMerger<EdgeCoverMerger>(config.merger_factor, config.merger_adjustment);
		break;
	case Merger::LocalBiconnected:
		merger = GetFactoredAdjustedMerger<LocalBiconnectedMerger>(config.merger_factor, config.merger_adjustment);
		break;
	case Merger::Solar:
		merger = GetAdjustedMerger<SolarMerger>(config.merger_adjustment);
		break;
	case Merger::IndependentSet:
		merger = GetAdjustedMerger<IndependentSetMerger>(config.merger_adjustment);
		break;
	}

	// Get the scaling type. As I do not want to expose any OGDF to Python,
	// there is this intermediate step.
	ScalingLayout::ScalingType scaling_type = ScalingLayout::ScalingType::RelativeToAvgLength;
	switch (config.sl_scaling_type)
	{
	case ScalingType::Absolute:
		scaling_type = ScalingLayout::ScalingType::Absolute;
		break;
	case ScalingType::RelativeToAvgLength:
		scaling_type = ScalingLayout::ScalingType::RelativeToAvgLength;
		break;
	case ScalingType::RelativeToDesiredLength:
		scaling_type = ScalingLayout::ScalingType::RelativeToDesiredLength;
		break;
	case ScalingType::RelativeToDrawing:
		scaling_type = ScalingLayout::ScalingType::RelativeToDrawing;
		break;
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

	std::vector<float> x(vertex_count);
	std::vector<float> y(vertex_count);

	int i = 0;
	for (node v : g.nodes)
	{
		x[i] = ga.x(v);
		y[i] = ga.y(v);
		i++;
	}

	return std::make_tuple(x, y);
}
