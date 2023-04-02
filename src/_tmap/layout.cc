/**
 * @file layout.cc
 * @author Daniel Probst (daenuprobst@gmail.com)
 * @brief Functions used for generating graph layouts from LSHForest instances
 * and edge lists.
 * @version 0.1
 * @date 2019-06-17
 *
 */

#include "layout.hh"
#include "lshforest.hh"

using namespace ogdf;

template <class T>
static MultilevelBuilder *
GetFactoredAdjustedMerger(double factor = 2.0, int adjustment = 0)
{
  T *merger = new T();
  merger->setFactor(factor);
  merger->setEdgeLengthAdjustment(adjustment);
  return merger;
}

template <class T>
static MultilevelBuilder *
GetAdjustedMerger(int adjustment = 0)
{
  T *merger = new T();
  merger->setEdgeLengthAdjustment(adjustment);
  return merger;
}

static InitialPlacer *
GetBarycenterPlacer()
{
  BarycenterPlacer *placer = new BarycenterPlacer();
  placer->weightedPositionPriority(true);
  return placer;
}

static InitialPlacer *
GetSolarPlacer()
{
  SolarPlacer *placer = new SolarPlacer();
  return placer;
}

static InitialPlacer *
GetCirclePlacer()
{
  CirclePlacer *placer = new CirclePlacer();
  return placer;
}

static InitialPlacer *
GetMedianPlacer()
{
  MedianPlacer *placer = new MedianPlacer();
  return placer;
}

static InitialPlacer *
GetRandomPlacer()
{
  RandomPlacer *placer = new RandomPlacer();
  return placer;
}

static InitialPlacer *
GetZeroPlacer()
{
  ZeroPlacer *placer = new ZeroPlacer();
  return placer;
}

std::vector<std::vector<uint32_t>>
GetTreesFromForest(const Graph &g)
{
  NodeArray<int> connected_component_ids(g);
  List<node> isolated_nodes;
  int n_connected_components =
      connectedComponents(g, connected_component_ids, &isolated_nodes);

  std::vector<std::vector<uint32_t>> connected_components(
      n_connected_components);

  for (int i = 0; i < n_connected_components; i++)
    connected_components[i] = std::vector<uint32_t>();

  uint32_t i = 0;
  for (auto id : connected_component_ids)
    connected_components[id].emplace_back(i++);

  std::sort(connected_components.begin(),
            connected_components.end(),
            [](const std::vector<uint32_t> &a, const std::vector<uint32_t> &b)
            {
              return a.size() > b.size();
            });

  return connected_components;
}

void RemoveDisconnectedComponents(Graph &g)
{
  GraphCopy t;
  Graph::CCsInfo info(g);

  int max = -1;
  int cc = 0;
  int num_ccs = info.numberOfCCs();

  for (int i = 0; i < num_ccs; i++)
  {
    int cc_size = info.numberOfNodes(i);
    if (cc_size > max)
    {
      max = cc_size;
      cc = i;
    }
  }

  for (int i = 0; i < num_ccs; i++)
  {
    if (i == cc)
      continue;

    for (int j = info.startNode(i); j < info.stopNode(i); ++j)
    {
      node v = info.v(j);
      g.delNode(v);
    }
  }
}

void ConnectGraph(Graph &g,
                  std::vector<node> &index_to_node,
                  tmap::LSHForest &lsh_forest)
{
  auto trees = GetTreesFromForest(g);

  for (size_t i = 1; i < trees.size(); i++)
  {
    for (uint32_t v : trees[i])
    {
      auto nns = lsh_forest.QueryLinearScanExcludeById(v, 1, trees[i], 10);

      if (nns.size() == 0)
        continue;

      g.newEdge(index_to_node[v],
                index_to_node[std::get<1>(nns[0])],
                std::get<0>(nns[0]));
    }

    trees = GetTreesFromForest(g);
  }
}

std::vector<std::pair<uint32_t, uint32_t>>
BFS(tmap::GraphProperties gp, uint32_t source,
    const std::vector<uint32_t> &targets)
{
  std::vector<std::pair<uint32_t, uint32_t>> result;
  std::vector<bool> visited(gp.adjacency_list.size(), false);
  std::queue<uint32_t> q;

  // Using max as the level delimiter
  uint32_t max = std::numeric_limits<uint32_t>::max();

  visited[source] = true;
  q.push(source);
  q.push(max);

  uint32_t depth = 0;
  while (q.size() > 1)
  {
    auto s = q.front();
    q.pop();

    if (s == max)
    {
      depth++;
      q.push(max);
      continue;
    }

    visited[s] = true;

    if (std::find(targets.begin(), targets.end(), s) != targets.end())
      result.push_back(std::make_pair(s, depth));

    if (result.size() == targets.size())
      break;

    for (size_t i = 0; i < gp.adjacency_list[s].size(); i++)
    {
      uint32_t t = gp.adjacency_list[s][i].first;

      if (visited[t])
        continue;

      q.push(t);
    }
  }

  return result;
}

std::vector<uint32_t>
BFSAll(tmap::GraphProperties &gp, uint32_t source)
{
  std::vector<uint32_t> result(gp.adjacency_list.size(), 0);
  std::vector<bool> visited(gp.adjacency_list.size(), false);
  std::queue<uint32_t> q;

  // Using max as the level delimiter
  uint32_t max = std::numeric_limits<uint32_t>::max();

  visited[source] = true;
  q.push(source);
  q.push(max);

  uint32_t depth = 0;
  while (q.size() > 1)
  {
    auto s = q.front();
    q.pop();
    if (s == max)
    {
      depth++;
      q.push(max);
      continue;
    }

    visited[s] = true;
    result[s] = depth;

    for (size_t i = 0; i < gp.adjacency_list[s].size(); i++)
    {
      uint32_t t = gp.adjacency_list[s][i].first;

      if (visited[t])
        continue;

      q.push(t);
    }
  }

  return result;
}

std::tuple<std::vector<uint32_t>, std::vector<uint32_t>, std::vector<float>>
tmap::MSTFromLSHForest(tmap::LSHForest &lsh_forest,
                       uint32_t k,
                       uint32_t kc)
{
  EdgeWeightedGraph<float> g;
  uint32_t vertex_count = lsh_forest.size();

  std::vector<uint32_t> from;
  std::vector<uint32_t> to;
  std::vector<float> weight;

  lsh_forest.GetKNNGraph(from, to, weight, k, kc);

  std::vector<node> index_to_node(lsh_forest.size());

  for (uint32_t i = 0; i < vertex_count; i++)
    index_to_node[i] = g.newNode();

  for (std::vector<uint32_t>::size_type i = 0; i != from.size(); i++)
    g.newEdge(index_to_node[from[i]], index_to_node[to[i]], weight[i]);

  ogdf::EdgeArray<float> edge_weights = g.edgeWeights();
  ogdf::makeMinimumSpanningTree(g, edge_weights);

  std::vector<uint32_t> x;
  std::vector<uint32_t> y;
  std::vector<float> w;

  for (edge e : g.edges)
  {
    x.emplace_back(e->source()->index());
    y.emplace_back(e->target()->index());
    w.emplace_back(edge_weights[e]);
  }

  return std::make_tuple(x, y, w);
}

std::tuple<std::vector<uint32_t>, std::vector<uint32_t>, std::vector<float>>
tmap::MSTFromEdgeList(
    uint32_t vertex_count,
    const std::vector<std::tuple<uint32_t, uint32_t, float>> &edges)
{
  tmap::GraphProperties gp;
  EdgeWeightedGraph<float> g;

  std::vector<std::vector<std::pair<uint32_t, float>>> adjacency_list(vertex_count);
  std::vector<node> index_to_node(vertex_count);

  for (uint32_t i = 0; i < vertex_count; i++)
    index_to_node[i] = g.newNode();

  // Normalize the edge weights, edge weights need to be positive
  float max_weight = 0.0f;
  for (size_t i = 0; i < edges.size(); i++)
    if (max_weight < std::get<2>(edges[i]))
      max_weight = std::get<2>(edges[i]);

  for (size_t i = 0; i < edges.size(); i++)
    g.newEdge(index_to_node[std::get<0>(edges[i])],
              index_to_node[std::get<1>(edges[i])],
              std::get<2>(edges[i]) / max_weight);

  ogdf::makeLoopFree(g);
  ogdf::makeParallelFreeUndirected(g);
  ogdf::EdgeArray<float> edge_weights = g.edgeWeights();
  ogdf::makeMinimumSpanningTree(g, edge_weights);

  std::vector<uint32_t> x;
  std::vector<uint32_t> y;
  std::vector<float> w;

  for (edge e : g.edges)
  {
    x.emplace_back(e->source()->index());
    y.emplace_back(e->target()->index());
    w.emplace_back(edge_weights[e]);
  }

  return std::make_tuple(x, y, w);
}

std::tuple<std::vector<float>,
           std::vector<float>,
           std::vector<uint32_t>,
           std::vector<uint32_t>,
           tmap::GraphProperties>
tmap::LayoutFromLSHForest(tmap::LSHForest &lsh_forest,
                          tmap::LayoutConfiguration config,
                          bool keep_knn,
                          bool create_mst,
                          bool clear_lsh_forest)
{
  tmap::GraphProperties gp;
  EdgeWeightedGraph<float> g;
  uint32_t vertex_count = lsh_forest.size();

  std::vector<uint32_t> from;
  std::vector<uint32_t> to;
  std::vector<float> weight;
  std::vector<uint32_t> degrees(vertex_count);
  std::vector<std::vector<std::pair<uint32_t, float>>> adjacency_list(vertex_count);

  lsh_forest.GetKNNGraph(from, to, weight, config.k, config.kc);

  if (clear_lsh_forest)
    lsh_forest.Clear();

  std::vector<node> index_to_node(lsh_forest.size());

  for (uint32_t i = 0; i < vertex_count; i++)
    index_to_node[i] = g.newNode();

  for (std::vector<uint32_t>::size_type i = 0; i != from.size(); i++)
    if (weight[i] >= 0.0f)
      g.newEdge(index_to_node[from[i]], index_to_node[to[i]], weight[i]);

  ogdf::makeLoopFree(g);
  ogdf::makeParallelFreeUndirected(g);
  ogdf::EdgeArray<float> edge_weights = g.edgeWeights();

  std::priority_queue<
      std::tuple<float, uint32_t, uint32_t>,
      std::vector<std::tuple<float, uint32_t, uint32_t>>,
      std::less<std::tuple<float, uint32_t, uint32_t>>>
      pq;

  // for (std::vector<uint32_t>::size_type i = 0; i != from.size(); i++)
  // {
  //   if (weight[i] < 0.0f)
  //     continue;

  //   pq.push(std::make_tuple(
  //       weight[i],
  //       from[i],
  //       to[i]));

  //   if (pq.size() > 999)
  //     pq.pop();
  // }

  uint32_t i = 0;

  // Get the adjancency list for the knn graph
  if (keep_knn)
  {
    std::vector<std::vector<std::pair<uint32_t, float>>> adjacency_list_knn(vertex_count);
    i = 0;
    for (node v : g.nodes)
    {
      adjacency_list_knn[i] = std::vector<std::pair<uint32_t, float>>(v->adjEntries.size());
      int j = 0;
      for (adjEntry adj : v->adjEntries)
      {
        uint32_t neighbor = adj->theEdge()->opposite(v)->index();
        float edge_weight = edge_weights[adj->theEdge()];
        adjacency_list_knn[i][j++] = std::make_pair(neighbor, edge_weight);
      }

      i++;
    }

    gp.adjacency_list_knn = adjacency_list_knn;
  }

  i = 0;
  for (node v : g.nodes)
    degrees[i++] = v->degree();

  gp.degrees = degrees;

  if (create_mst)
  {
    gp.mst_weight = ogdf::makeMinimumSpanningTree(g, edge_weights);
  }

  while (!pq.empty())
  {
    std::tuple<float, uint32_t, uint32_t> edge_tuple = pq.top();

    // std::cout << std::get<0>(edge_tuple) << ' ' << std::get<1>(edge_tuple) << ' ' << std::get<2>(edge_tuple) << std::endl;

    g.newEdge(index_to_node[std::get<0>(edge_tuple)],
              index_to_node[std::get<1>(edge_tuple)],
              std::get<2>(edge_tuple));

    pq.pop();
  }
  ogdf::makeLoopFree(g);
  ogdf::makeParallelFreeUndirected(g);

  i = 0;
  for (node v : g.nodes)
  {
    adjacency_list[i] = std::vector<std::pair<uint32_t, float>>(v->adjEntries.size());
    int j = 0;
    for (adjEntry adj : v->adjEntries)
    {
      uint32_t neighbor = adj->theEdge()->opposite(v)->index();
      float edge_weight = edge_weights[adj->theEdge()];
      adjacency_list[i][j++] = std::make_pair(neighbor, edge_weight);
    }

    i++;
  }

  gp.adjacency_list = adjacency_list;

  return LayoutInternal(g, vertex_count, config, gp);
}

std::tuple<std::vector<float>,
           std::vector<float>,
           std::vector<uint32_t>,
           std::vector<uint32_t>,
           tmap::GraphProperties>
tmap::LayoutFromEdgeList(
    uint32_t vertex_count,
    const std::vector<std::tuple<uint32_t, uint32_t, float>> &edges,
    tmap::LayoutConfiguration config,
    bool keep_knn,
    bool create_mst)
{
  tmap::GraphProperties gp;
  EdgeWeightedGraph<float> g;

  std::vector<std::vector<std::pair<uint32_t, float>>> adjacency_list(vertex_count);
  std::vector<uint32_t> degrees(vertex_count);
  std::vector<node> index_to_node(vertex_count);

  for (uint32_t i = 0; i < vertex_count; i++)
    index_to_node[i] = g.newNode();

  // Normalize the edge weights, edge weights need to be positive
  float max_weight = 0.0f;
  // std::priority_queue<
  //     std::tuple<float, uint32_t, uint32_t>,
  //     std::vector<std::tuple<float, uint32_t, uint32_t>>,
  //     std::less<std::tuple<float, uint32_t, uint32_t>>>
  //     pq;

  for (size_t i = 0; i < edges.size(); i++)
  {
    float edge_weight = std::get<2>(edges[i]);
    // if (edge_weight > 0.0f)
    // {
    //   pq.push(std::make_tuple(
    //       edge_weight,
    //       std::get<0>(edges[i]),
    //       std::get<1>(edges[i])));

    //   if (pq.size() > 9)
    //     pq.pop();
    // }

    if (max_weight < edge_weight)
      max_weight = edge_weight;
  }

  for (size_t i = 0; i < edges.size(); i++)
    g.newEdge(index_to_node[std::get<0>(edges[i])],
              index_to_node[std::get<1>(edges[i])],
              std::get<2>(edges[i]) / max_weight);

  ogdf::makeLoopFree(g);
  ogdf::makeParallelFreeUndirected(g);
  ogdf::EdgeArray<float> edge_weights = g.edgeWeights();

  // Get the adjancency list for the knn graph
  uint32_t i = 0;
  if (keep_knn)
  {
    std::vector<std::vector<std::pair<uint32_t, float>>> adjacency_list_knn(vertex_count);
    i = 0;
    for (node v : g.nodes)
    {
      adjacency_list_knn[i] = std::vector<std::pair<uint32_t, float>>(v->adjEntries.size());
      int j = 0;
      for (adjEntry adj : v->adjEntries)
      {
        uint32_t neighbor = adj->theEdge()->opposite(v)->index();
        float edge_weight = edge_weights[adj->theEdge()];
        adjacency_list_knn[i][j++] = std::make_pair(neighbor, edge_weight);
      }

      i++;
    }

    gp.adjacency_list_knn = adjacency_list_knn;
  }

  i = 0;
  for (node v : g.nodes)
    degrees[i++] = v->degree();

  gp.degrees = degrees;

  if (create_mst)
  {
    gp.mst_weight = ogdf::makeMinimumSpanningTree(g, edge_weights);
  }

  // while (!pq.empty())
  // {
  //   std::tuple<float, uint32_t, uint32_t> edge_tuple = pq.top();

  //   // std::cout << std::get<0>(edge_tuple) << ' ' << std::get<1>(edge_tuple) << ' ' << std::get<2>(edge_tuple) << std::endl;

  //   g.newEdge(index_to_node[std::get<0>(edge_tuple)],
  //             index_to_node[std::get<1>(edge_tuple)],
  //             std::get<2>(edge_tuple));

  //   pq.pop();
  // }
  // ogdf::makeLoopFree(g);
  // ogdf::makeParallelFreeUndirected(g);

  i = 0;
  for (node v : g.nodes)
  {
    adjacency_list[i] = std::vector<std::pair<uint32_t, float>>(v->adjEntries.size());
    int j = 0;
    for (adjEntry adj : v->adjEntries)
    {
      uint32_t neighbor = adj->theEdge()->opposite(v)->index();
      float edge_weight = edge_weights[adj->theEdge()];
      adjacency_list[i][j++] = std::make_pair(neighbor, edge_weight);
    }

    i++;
  }

  gp.adjacency_list = adjacency_list;

  GraphAttributes graph_attributes(g);
  MultilevelGraph mlg(g);

  return LayoutInternal(g, vertex_count, config, gp);
}

std::tuple<std::vector<float>,
           std::vector<float>,
           std::vector<uint32_t>,
           std::vector<uint32_t>,
           tmap::GraphProperties>
tmap::LayoutInternal(EdgeWeightedGraph<float> &g,
                     uint32_t vertex_count,
                     LayoutConfiguration config,
                     GraphProperties &gp)
{
  // Check for isolated nodes. If there are isolated nodes,
  // call placement step later on.
  NodeArray<int> connected_components(g);
  List<node> isolated_nodes;
  int n_connected_components =
      connectedComponents(g, connected_components, &isolated_nodes);
  gp.n_connected_components = n_connected_components;
  gp.n_isolated_vertices = isolated_nodes.size();

  GraphAttributes ga(g);

  ga.setAllHeight(config.node_size);
  ga.setAllWidth(config.node_size);

  // for (edge e : g.edges)
  //   g.setWeight(e, 1.0);

  // Starting the layout
  MultilevelGraph mlg(ga);

  // The FastMultipoleEmbedder is used for the single level layout.
  FastMultipoleEmbedder *fme = new FastMultipoleEmbedder();
  fme->setNumIterations(config.fme_iterations);
  fme->setRandomize(config.fme_randomize);
  // fme->setNumberOfThreads(config.fme_threads);
  fme->setMultipolePrec(config.fme_precision);
  fme->setDefaultEdgeLength(1);
  fme->setDefaultNodeSize(1);

  // FMMMLayout *fme = new FMMMLayout();

  // fme->useHighLevelOptions(true);
  // // fme->unitEdgeLength(15.0);
  // fme->newInitialPlacement(true);
  // fme->qualityVersusSpeed(FMMMOptions::QualityVsSpeed::GorgeousAndEfficient);

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
    merger = GetFactoredAdjustedMerger<EdgeCoverMerger>(
        config.merger_factor, config.merger_adjustment);
    break;
  case Merger::LocalBiconnected:
    merger = GetFactoredAdjustedMerger<LocalBiconnectedMerger>(
        config.merger_factor, config.merger_adjustment);
    break;
  case Merger::Solar:
    merger = GetAdjustedMerger<SolarMerger>(config.merger_adjustment);
    break;
  case Merger::IndependentSet:
    merger =
        GetAdjustedMerger<IndependentSetMerger>(config.merger_adjustment);
    break;
  }

  // Get the scaling type. As I do not want to expose any OGDF to Python,
  // there is this intermediate step.
  ScalingLayout::ScalingType scaling_type =
      ScalingLayout::ScalingType::RelativeToDrawing;
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
  // In this example a FastMultipoleEmbedder with zero iterations is used for
  // postprocessing.
  sl->setExtraScalingSteps(config.sl_extra_scaling_steps);
  sl->setScalingType(scaling_type);
  sl->setScaling(config.sl_scaling_min, config.sl_scaling_max);

  // Then the ModularMultilevelMixer is created.
  ModularMultilevelMixer *mmm = new ModularMultilevelMixer;
  mmm->setLayoutRepeats(config.mmm_repeats);
  // The single level layout, the placer and the merger are set.
  mmm->setLevelLayoutModule(sl);
  mmm->setInitialPlacer(placer);
  mmm->setMultilevelBuilder(merger);

  if (config.sl_scaling_type == ScalingType::Absolute)
    sl->setMMM(mmm);

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
    ppl.setRandomizePositions(false);

    ppl.call(mlg);
  }
  else
    mmm->call(mlg);

  mlg.exportAttributes(ga);

  // GraphAttributes ga(g);
  // for (node v : g.nodes)
  //   ga.width(v) = ga.height(v) = 0.1;

  // FMMMLayout fmmm;

  // fmmm.useHighLevelOptions(true);
  // fmmm.unitEdgeLength(15.0);
  // fmmm.newInitialPlacement(true);
  // fmmm.qualityVersusSpeed(FMMMOptions::QualityVsSpeed::GorgeousAndEfficient);

  // fmmm.call(ga);
  // SpringEmbedderFRExact layout;
  // layout.call(ga);

  std::vector<float>
      x(vertex_count);
  std::vector<float> y(vertex_count);

  std::vector<uint32_t> s(g.edges.size());
  std::vector<uint32_t> t(g.edges.size());

  int i = 0;
  for (node v : g.nodes)
  {
    x[i] = ga.x(v);
    y[i] = ga.y(v);
    i++;
  }

  // Also norm distances, as units are meaningless
  // Center on 0
  float max_x = *max_element(x.begin(), x.end());
  float max_y = *max_element(y.begin(), y.end());
  float min_x = *min_element(x.begin(), x.end());
  float min_y = *min_element(y.begin(), y.end());
  float diff_x = max_x - min_x;
  float diff_y = max_y - min_y;

  for (size_t i = 0; i < x.size(); i++)
  {
    x[i] = (x[i] - min_x) / diff_x - 0.5;
    y[i] = (y[i] - min_y) / diff_y - 0.5;
  }

  i = 0;
  for (edge e : g.edges)
  {
    s[i] = e->source()->index();
    t[i] = e->target()->index();
    i++;
  }

  return std::make_tuple(x, y, s, t, gp);
}

std::vector<std::tuple<uint32_t, float, uint32_t>>
tmap::VertexQuality(tmap::GraphProperties &gp, uint32_t v)
{
  if (gp.adjacency_list_knn.size() < 1)
    throw std::runtime_error("The GraphProperties object does not contain an adjancency list of the original knn graph. Run layout with 'keep_knn' set to true.");

  if (v >= gp.adjacency_list_knn.size())
    throw std::runtime_error("The argument v is out of range.");

  std::vector<std::tuple<uint32_t, float, uint32_t>> result(gp.adjacency_list_knn[v].size());
  std::vector<uint32_t> targets(gp.adjacency_list_knn[v].size());

  for (size_t i = 0; i < gp.adjacency_list_knn[v].size(); i++)
    targets[i] = gp.adjacency_list_knn[v][i].first;

  auto bfs = BFS(gp, v, targets);

  for (size_t i = 0; i < bfs.size(); i++)
  {
    result[i] = std::make_tuple(bfs[i].first, gp.adjacency_list_knn[v][i].second, bfs[i].second);
  }

  // Sort the targets by distance
  std::sort(result.begin(), result.end(), [](auto &left, auto &right)
            { return std::get<1>(left) < std::get<1>(right); });

  return result;
}

std::vector<float>
tmap::MeanQuality(GraphProperties &gp)
{
  // Get the max adjacency list knn size, as they can differ
  size_t max_size = 0;
  for (size_t i = 0; i < gp.adjacency_list_knn.size(); i++)
    if (gp.adjacency_list_knn[i].size() > max_size)
      max_size = gp.adjacency_list_knn[i].size();

  std::vector<float> result(max_size, 0.0);
  std::vector<uint32_t> counts(max_size, 0);

  for (size_t i = 0; i < gp.adjacency_list.size(); i++)
  {
    auto r = tmap::VertexQuality(gp, i);
    for (size_t j = 0; j < r.size(); j++)
    {
      result[j] += std::get<2>(r[j]);
      counts[j] += 1;
    }
  }

  for (size_t i = 0; i < max_size; i++)
    result[i] /= counts[i];

  return result;
}

std::vector<uint32_t>
tmap::GetTopologicalDistances(tmap::GraphProperties &gp, uint32_t v)
{
  return BFSAll(gp, v);
}

std::tuple<std::vector<float>,
           std::vector<float>,
           std::vector<float>,
           std::vector<float>>
tmap::MakeEdgeList(
    std::vector<float> x, std::vector<float> y,
    std::vector<uint32_t> s, std::vector<uint32_t> t)
{
  std::vector<float> x1(s.size());
  std::vector<float> y1(s.size());
  std::vector<float> x2(s.size());
  std::vector<float> y2(s.size());

  for (size_t i = 0; i < s.size(); i++)
  {
    x1[i] = x[s[i]];
    y1[i] = y[s[i]];
    x2[i] = x[t[i]];
    y2[i] = y[t[i]];
  }

  return std::make_tuple(x1, y1, x2, y2);
}