/**
 * @file analyse.cc
 * @author Daniel Probst (daenuprobst@gmail.com)
 * @brief Analysing TMAP graphs and trees.
 * @version 0.1
 * @date 2019-06-17
 *
 */

#include "analyse.hh"

typedef std::vector<std::vector<std::pair<uint32_t, float>>> AdjacencyList;
typedef std::tuple<uint32_t, uint32_t, float> Edge;
typedef std::vector<std::vector<uint32_t>> ConnectedComponents;

// Private helper functions
namespace {
  std::tuple<float, float> 
  mean_stdev(const AdjacencyList& adjacency_list) 
  {
    std::vector<float> weights;

    for (size_t i = 0; i < adjacency_list.size(); i++)
      for (size_t j = 0; j < adjacency_list[i].size(); j++)
        if (adjacency_list[i][j].first > i)
          weights.emplace_back(adjacency_list[i][j].second);

    float sum = std::accumulate(weights.begin(), weights.end(), 0.0);
    float mean = sum / weights.size();

    std::vector<float> diff(weights.size());
    std::transform(weights.begin(), weights.end(), diff.begin(),
                  std::bind(std::minus<float>(), std::placeholders::_1, mean));
    float sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
    float stdev = std::sqrt(sq_sum / weights.size());

    return std::make_tuple(mean, stdev);
  }

  float
  weighted_stdev(
    const AdjacencyList& adjacency_list,
    const ConnectedComponents& connected_components) 
  {
    std::vector<std::vector<float>> cc_weights(connected_components.size());
    std::vector<float> stdevs(connected_components.size());
    std::vector<uint32_t> cc_map(adjacency_list.size());
    
    for (size_t i = 0; i < connected_components.size(); i++)
      for (size_t j = 0; j < connected_components[i].size(); j++)
        cc_map[connected_components[i][j]] = i;

    for (size_t i = 0; i < adjacency_list.size(); i++)
      for (size_t j = 0; j < adjacency_list[i].size(); j++)
        if (adjacency_list[i][j].first > i)
          if (adjacency_list[i].size() > 1 || adjacency_list[adjacency_list[i][j].first].size() > 1)
            cc_weights[cc_map[i]].emplace_back(adjacency_list[i][j].second);

    for (size_t i = 0; i < cc_weights.size(); i++) {
      auto weights = cc_weights[i];

      float sum = std::accumulate(weights.begin(), weights.end(), 0.0);
      float mean = sum / weights.size();

      std::vector<float> diff(weights.size());
      std::transform(weights.begin(), weights.end(), diff.begin(),
                    std::bind(std::minus<float>(), std::placeholders::_1, mean));
      float sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
      float stdev = std::sqrt(sq_sum / weights.size());

      stdevs[i] = stdev;
    }

    
    float numerator = 0.0;
    float denominator = 0.0;

    for (size_t i = 0; i < connected_components.size(); i++) {
      if (cc_weights[i].size() > 0) {
        numerator += connected_components[i].size() * stdevs[i];
        denominator += connected_components[i].size();
      }
    }

    return numerator / denominator;
  }

  std::tuple<Edge, Edge>
  get_edge(const AdjacencyList& adjacency_list, 
    uint32_t u, 
    uint32_t v)
  {
    Edge a;
    Edge b;

    for (size_t i = 0; i < adjacency_list[u].size(); i++)
      if (adjacency_list[u][i].first == v)
        a = std::make_tuple(u, v, adjacency_list[u][i].second);

    for (size_t i = 0; i < adjacency_list[v].size(); i++)
      if (adjacency_list[v][i].first == u)
        b = std::make_tuple(v, u, adjacency_list[v][i].second);

    return std::make_tuple(a, b);
  }

  std::tuple<Edge, Edge>
  remove_edge(
    AdjacencyList& adjacency_list,
    uint32_t u,
    uint32_t v)
  {
    auto removed_edges = get_edge(adjacency_list, u, v);

    adjacency_list[u].erase(
      std::remove_if(
        adjacency_list[u].begin(), adjacency_list[u].end(),
        [v](std::pair<uint32_t, float> p) {
          return p.first == v;
        }
      ), adjacency_list[u].end()
    );

    adjacency_list[v].erase(
      std::remove_if(
        adjacency_list[v].begin(), adjacency_list[v].end(),
        [u](std::pair<uint32_t, float> p) {
          return p.first == u;
        }
      ), adjacency_list[v].end()
    );

    return removed_edges;
  }

  void
  add_edge(AdjacencyList& adjacency_list,
    Edge edge)
  {
    adjacency_list[std::get<0>(edge)].emplace_back(std::make_pair(std::get<1>(edge), std::get<2>(edge)));
    adjacency_list[std::get<1>(edge)].emplace_back(std::make_pair(std::get<0>(edge), std::get<2>(edge)));
  }

  ConnectedComponents
  get_connected_components(
    const AdjacencyList& adjacency_list)
  {
    ConnectedComponents connected_components;

    std::vector<uint8_t> visited(adjacency_list.size(), 0);

    for (size_t v = 0; v < adjacency_list.size(); v++) {
      if (visited[v] == 1)
        continue;
      
      std::vector<uint32_t> connected_component;

      // Do a BFS, so let's use a queue
      std::queue<uint32_t> q;
      q.push(v);
      visited[v] = 1;
      connected_component.emplace_back(v);

      while(!q.empty()) {
        uint32_t w = q.front();
        q.pop();

        for (size_t i = 0; i < adjacency_list[w].size(); i++) {
          auto k = adjacency_list[w][i].first;

          if (visited[k] == 0) {
            visited[k] = 1;
            q.push(k);
            connected_component.emplace_back(k);
          }
        }
      }

      connected_components.push_back(connected_component);
    }

    return connected_components;
  }
}

std::vector<std::tuple<uint32_t, std::vector<uint32_t>>>
tmap::GetClusters(
  const tmap::GraphProperties& gp,
  const std::vector<uint32_t>& classes) 
{
  std::vector<std::tuple<uint32_t, std::vector<uint32_t>>> result;
  std::vector<uint8_t> visited(gp.adjacency_list.size(), 0);

  for (size_t i = 0; i < gp.adjacency_list.size(); i++) {
    if (visited[i] == 1)
      continue;
    
    visited[i] = 1;
    // If this vertex hasn't been visited yet, create a new cluster with it's class
    auto cluster = std::make_tuple(i, std::vector<uint32_t>());
    std::get<1>(cluster).push_back(i);


    std::stack<uint32_t> stack;
    stack.push(i);

    // Kick off the discovery
    while (!stack.empty()) {
      uint32_t v = stack.top();
      stack.pop();

      for (size_t j = 0; j < gp.adjacency_list[v].size(); j++) {
        if (visited[gp.adjacency_list[v][j].first] == 1)
          continue;
        
        if (classes[i] == classes[gp.adjacency_list[v][j].first]) {
          std::get<1>(cluster).push_back(gp.adjacency_list[v][j].first);
          stack.push(gp.adjacency_list[v][j].first);
          visited[gp.adjacency_list[v][j].first] = 1;
        }
      }
    }
    
    result.push_back(cluster);
  }

  return result;
}

std::vector<std::vector<uint32_t>>
tmap::MSDR(tmap::GraphProperties gp)
{
  uint32_t u_remove, v_remove;
  float mean, stdev, temp;
  std::tie(mean, stdev) = mean_stdev(gp.adjacency_list);
  auto cc = get_connected_components(gp.adjacency_list);


  float epsilon = 0.0001;
  std::vector<float> d_stdev(1, 0.0);
  size_t i = 0;

  do {
    i += 1;
    temp = stdev;
    d_stdev.push_back(0.0);

    // Remove the edge which removal leads to the maximal stdev reduction
    for (size_t u = 0; u < gp.adjacency_list.size(); u++) {
      for (size_t j = 0; j < gp.adjacency_list[u].size(); j++) {
        uint32_t v = gp.adjacency_list[u][j].first;

        // Do not consider leafs
        if (gp.adjacency_list[u].size() == 1 || gp.adjacency_list[v].size() == 1)
          continue;

        auto tmp_edges = remove_edge(gp.adjacency_list, u, v);
        cc = get_connected_components(gp.adjacency_list);
        stdev = weighted_stdev(gp.adjacency_list, cc);

        add_edge(gp.adjacency_list, std::get<0>(tmp_edges));
        add_edge(gp.adjacency_list, std::get<1>(tmp_edges));

        if (d_stdev[i] < stdev - temp) {
          d_stdev[i] = stdev - temp;
          u_remove = u;
          v_remove = v;
        }
      }
    }

    remove_edge(gp.adjacency_list, u_remove, v_remove);
    stdev = temp - d_stdev[i];

    // std::cout << "stdev: " << stdev << std::endl;
    // std::cout << d_stdev[i] << ", " << d_stdev[i-1] << std::endl;
    std::cout << std::abs(d_stdev[i] - d_stdev[i-1]) << std::endl;
    std::cout << std::abs(epsilon * (d_stdev[i] + 1)) << std::endl;

  } while (i < 19);
  // } while(std::abs(d_stdev[i] - d_stdev[i-1]) < std::abs(epsilon * (d_stdev[i] + 1)));  

  return cc;
}
