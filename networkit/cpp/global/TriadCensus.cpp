// networkit-format
#include <algorithm>
#include <array>
#include <cstdint>
#include <iterator>
#include <stdexcept>
#include <vector>

#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/global/TriadCensus.hpp>

namespace NetworKit {

namespace {

constexpr uint8_t forward = 1;
constexpr uint8_t backward = 2;

// A node adjacent in either direction. Directions are seen from the node that stores the neighbor.
struct Neighbor {
    node id;
    uint8_t directions;
};

// Each pair of nodes contributes two bits to the code of a triple of nodes v, u and w: the lower
// ones for (v, u), then (v, w) and (u, w). Taken from the paper of Batagelj and Mrvar.
constexpr std::array<uint8_t, 64> triadOfCode{
    0, 1,  1,  2,  1, 3,  5, 7,  1,  5,  4,  6,  2, 7,  6,  10, 1,  5,  3,  7, 4, 8,
    8, 12, 5,  9,  8, 13, 6, 13, 11, 14, 1,  4,  5, 6,  5,  8,  9,  13, 3,  8, 8, 11,
    7, 12, 13, 14, 2, 6,  7, 10, 6,  11, 13, 14, 7, 13, 12, 14, 10, 14, 14, 15};

// Divides before multiplying, as the product overflows from about 2.6 million nodes.
count numberOfTriples(count n) {
    if (n < 3) {
        return 0;
    }
    const count pairs = n * (n - 1) / 2;
    return pairs % 3 == 0 ? pairs / 3 * (n - 2) : pairs * ((n - 2) / 3);
}

uint8_t directionsBetween(const std::vector<Neighbor> &neighbors, node v) {
    const auto it = std::lower_bound(
        neighbors.begin(), neighbors.end(), v,
        [](const Neighbor &neighbor, node id) -> bool { return neighbor.id < id; });
    return it != neighbors.end() && it->id == v ? it->directions : 0;
}

} // namespace

TriadCensus::TriadCensus(const Graph &G) : G(&G) {
    if (!G.isDirected()) {
        throw std::runtime_error("TriadCensus: the graph must be directed.");
    }
}

void TriadCensus::run() {
    Aux::SignalHandler handler;
    const count z = G->upperNodeIdBound();
    const count n = G->numberOfNodes();

    counts.fill(0);

    std::vector<std::vector<Neighbor>> neighbors(z);
    G->forEdges([&](node u, node v) -> void {
        if (u != v) {
            neighbors[u].push_back({v, forward});
            neighbors[v].push_back({u, backward});
        }
    });

#pragma omp parallel for schedule(dynamic, 64)
    for (omp_index i = 0; i < static_cast<omp_index>(z); ++i) {
        auto &list = neighbors[i];
        std::sort(list.begin(), list.end(),
                  [](const Neighbor &a, const Neighbor &b) -> bool { return a.id < b.id; });

        // Parallel edges and edges in both directions share an entry.
        auto last = list.begin();
        for (auto it = list.begin(); it != list.end(); ++it) {
            if (last != list.begin() && std::prev(last)->id == it->id) {
                std::prev(last)->directions |= it->directions;
            } else {
                *last++ = *it;
            }
        }
        list.erase(last, list.end());
    }
    handler.assureRunning();

    // Each connected triple is counted once, from its lowest pair of adjacent nodes.
#pragma omp parallel
    {
        std::array<count, numberOfTriads> localCounts{};
        std::vector<node> others;

#pragma omp for schedule(dynamic, 16)
        for (omp_index i = 0; i < static_cast<omp_index>(z); ++i) {
            const node v = static_cast<node>(i);
            for (const Neighbor &neighbor : neighbors[v]) {
                const node u = neighbor.id;
                if (u < v) {
                    continue;
                }

                others.clear();
                const auto addOther = [&](const Neighbor &other) -> void {
                    if (other.id != u && other.id != v) {
                        others.push_back(other.id);
                    }
                };
                const auto &ofV = neighbors[v];
                const auto &ofU = neighbors[u];
                auto itV = ofV.begin();
                auto itU = ofU.begin();
                while (itV != ofV.end() || itU != ofU.end()) {
                    if (itU == ofU.end() || (itV != ofV.end() && itV->id < itU->id)) {
                        addOther(*itV++);
                    } else if (itV == ofV.end() || itU->id < itV->id) {
                        addOther(*itU++);
                    } else {
                        addOther(*itV++);
                        ++itU;
                    }
                }

                // All remaining nodes are adjacent to neither of the two.
                const Triad dyad =
                    neighbor.directions == (forward | backward) ? Triad::T102 : Triad::T012;
                localCounts[static_cast<index>(dyad)] += n - others.size() - 2;

                for (const node w : others) {
                    const uint8_t betweenVW = directionsBetween(ofV, w);
                    if (u < w || (v < w && w < u && betweenVW == 0)) {
                        const index code =
                            neighbor.directions + 4 * betweenVW + 16 * directionsBetween(ofU, w);
                        ++localCounts[triadOfCode[code]];
                    }
                }
            }
        }

#pragma omp critical
        for (index triad = 0; triad < numberOfTriads; ++triad) {
            counts[triad] += localCounts[triad];
        }
    }

    count connectedTriples = 0;
    for (const count triples : counts) {
        connectedTriples += triples;
    }
    counts[static_cast<index>(Triad::T003)] = numberOfTriples(n) - connectedTriples;

    hasRun = true;
}

} /* namespace NetworKit */
