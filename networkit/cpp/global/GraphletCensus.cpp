// networkit-format
#include <algorithm>
#include <iterator>
#include <stdexcept>
#include <vector>

#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/global/GraphletCensus.hpp>

namespace NetworKit {

GraphletCensus::GraphletCensus(const Graph &G) : G(&G) {
    if (G.isDirected()) {
        throw std::runtime_error("GraphletCensus: the graph must be undirected.");
    }
}

void GraphletCensus::run() {
    Aux::SignalHandler handler;
    const count z = G->upperNodeIdBound();

    std::vector<std::vector<node>> neighbors(z);
    G->forEdges([&](node u, node v) -> void {
        if (u != v) {
            neighbors[u].push_back(v);
            neighbors[v].push_back(u);
        }
    });

#pragma omp parallel for schedule(dynamic, 64)
    for (omp_index u = 0; u < static_cast<omp_index>(z); ++u) {
        auto &list = neighbors[u];
        std::sort(list.begin(), list.end());
        list.erase(std::unique(list.begin(), list.end()), list.end());
    }
    handler.assureRunning();

    // Occurrences as subgraphs, which do not have to be induced.
    count wedges = 0;
    count triangles = 0;
    count paths = 0;
    count stars = 0;
    count cycles = 0;
    count paws = 0;
    count diamonds = 0;
    count cliques = 0;

#pragma omp parallel reduction(+ : wedges, triangles, paths, stars, cycles, paws, diamonds, cliques)
    {
        std::vector<node> common;
        std::vector<node> commonAfterW;
        std::vector<count> wedgesTo(z, 0);
        std::vector<node> reached;

#pragma omp for schedule(dynamic, 16)
        for (omp_index i = 0; i < static_cast<omp_index>(z); ++i) {
            const node u = static_cast<node>(i);
            const count degreeU = neighbors[u].size();
            if (degreeU >= 2) {
                wedges += degreeU * (degreeU - 1) / 2;
            }
            if (degreeU >= 3) {
                // Dividing in two steps, as the product of all terms overflows more easily.
                const count pairs = degreeU * (degreeU - 1) / 2;
                stars += pairs % 3 == 0 ? pairs / 3 * (degreeU - 2) : pairs * ((degreeU - 2) / 3);
            }

            // Two nodes with k common neighbors are the opposite corners of (k choose 2) cycles.
            for (const node v : neighbors[u]) {
                for (const node w : neighbors[v]) {
                    if (w > u) {
                        if (wedgesTo[w]++ == 0) {
                            reached.push_back(w);
                        }
                    }
                }
            }
            for (const node w : reached) {
                cycles += wedgesTo[w] * (wedgesTo[w] - 1) / 2;
                wedgesTo[w] = 0;
            }
            reached.clear();

            for (const node v : neighbors[u]) {
                if (v < u) {
                    continue;
                }

                const count degreeV = neighbors[v].size();
                common.clear();
                std::set_intersection(neighbors[u].begin(), neighbors[u].end(),
                                      neighbors[v].begin(), neighbors[v].end(),
                                      std::back_inserter(common));

                paths += (degreeU - 1) * (degreeV - 1);
                if (!common.empty()) {
                    diamonds += common.size() * (common.size() - 1) / 2;
                }

                for (const node w : common) {
                    if (w < v) {
                        continue;
                    }
                    ++triangles;
                    paws += degreeU + degreeV + neighbors[w].size() - 6;

                    commonAfterW.clear();
                    std::set_intersection(common.begin(), common.end(), neighbors[w].begin(),
                                          neighbors[w].end(), std::back_inserter(commonAfterW));
                    cliques += static_cast<count>(
                        commonAfterW.end()
                        - std::upper_bound(commonAfterW.begin(), commonAfterW.end(), w));
                }
            }
        }
    }

    // Each cycle was found from both of its diagonals. Paths through a triangle are not paths.
    cycles /= 2;
    paths -= 3 * triangles;

    // Remove the occurrences within larger graphlets, by how many of them each of those contains.
    const count inducedDiamonds = diamonds - 6 * cliques;
    const count inducedCycles = cycles - inducedDiamonds - 3 * cliques;
    const count inducedPaws = paws - 4 * inducedDiamonds - 12 * cliques;

    counts[static_cast<index>(Graphlet::Wedge)] = wedges - 3 * triangles;
    counts[static_cast<index>(Graphlet::Triangle)] = triangles;
    counts[static_cast<index>(Graphlet::Path)] =
        paths - 2 * inducedPaws - 4 * inducedCycles - 6 * inducedDiamonds - 12 * cliques;
    counts[static_cast<index>(Graphlet::Star)] =
        stars - inducedPaws - 2 * inducedDiamonds - 4 * cliques;
    counts[static_cast<index>(Graphlet::Cycle)] = inducedCycles;
    counts[static_cast<index>(Graphlet::Paw)] = inducedPaws;
    counts[static_cast<index>(Graphlet::Diamond)] = inducedDiamonds;
    counts[static_cast<index>(Graphlet::Clique)] = cliques;

    hasRun = true;
}

} /* namespace NetworKit */
