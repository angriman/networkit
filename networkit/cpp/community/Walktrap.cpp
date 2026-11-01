// networkit-format
#include <algorithm>
#include <functional>
#include <map>
#include <queue>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/community/Walktrap.hpp>
#include <networkit/structures/UnionFind.hpp>

namespace NetworKit {

namespace {

// Two adjacent communities: the variance increase of merging them, up to a constant factor, and the
// total weight of the edges between them.
struct Link {
    double delta;
    edgeweight weight;
};

} // namespace

Walktrap::Walktrap(const Graph &G, count walkLength)
    : CommunityDetectionAlgorithm(G), walkLength(walkLength) {
    if (walkLength == 0) {
        throw std::runtime_error("Walktrap: the walk length must be positive.");
    }
    G.forEdges([](node, node, edgeweight weight) -> void {
        if (weight < 0) {
            throw std::runtime_error("Walktrap: edge weights must not be negative.");
        }
    });
}

void Walktrap::run() {
    Aux::SignalHandler handler;
    const count z = G->upperNodeIdBound();

    // Node ids may have gaps, the walks are stored densely.
    std::vector<node> nodes;
    nodes.reserve(G->numberOfNodes());
    std::vector<index> positionOf(z, none);
    G->forNodes([&](node u) -> void {
        positionOf[u] = nodes.size();
        nodes.push_back(u);
    });
    const count n = nodes.size();

    // Each node gets a self-loop weighing as much as its average edge, which keeps the walks from
    // depending on the parity of the path lengths, as they would on bipartite graphs.
    std::vector<std::vector<std::pair<index, double>>> transitions(n);
    std::vector<double> degree(n, 0);
    for (index i = 0; i < n; ++i) {
        G->forEdgesOf(nodes[i], [&](node, node v, edgeweight weight) -> void {
            transitions[i].emplace_back(positionOf[v], weight);
            degree[i] += weight;
        });
        const double selfLoop = degree[i] > 0 ? degree[i] / transitions[i].size() : 1;
        transitions[i].emplace_back(i, selfLoop);
        degree[i] += selfLoop;
        for (auto &transition : transitions[i]) {
            transition.second /= degree[i];
        }
    }

    // Walk distribution of each community, indexed by community and then by destination.
    std::vector<std::vector<double>> walks(n);
#pragma omp parallel
    {
        std::vector<double> next(n);
#pragma omp for schedule(dynamic)
        for (omp_index start = 0; start < static_cast<omp_index>(n); ++start) {
            std::vector<double> current(n, 0);
            current[start] = 1;
            for (count step = 0; step < walkLength; ++step) {
                std::fill(next.begin(), next.end(), 0.);
                for (index i = 0; i < n; ++i) {
                    if (current[i] > 0) {
                        for (const auto &[to, probability] : transitions[i]) {
                            next[to] += current[i] * probability;
                        }
                    }
                }
                current.swap(next);
            }
            walks[start] = std::move(current);
        }
    }

    std::vector<count> size(n, 1);
    const auto mergeCost = [&](index a, index b) -> double {
        double distance = 0;
        for (index i = 0; i < n; ++i) {
            const double difference = walks[a][i] - walks[b][i];
            distance += difference * difference / degree[i];
        }
        return distance * static_cast<double>(size[a] * size[b])
               / static_cast<double>(size[a] + size[b]);
    };

    // Stale candidates are discarded when they are popped.
    using Candidate = std::tuple<double, index, index>;
    std::priority_queue<Candidate, std::vector<Candidate>, std::greater<>> candidates;
    std::vector<std::map<index, Link>> links(n);
    const auto updateCandidate = [&](index a, index b, const Link &link) -> void {
        links[a].insert_or_assign(b, link);
        links[b].insert_or_assign(a, link);
        candidates.emplace(link.delta, std::min(a, b), std::max(a, b));
    };

    // Zero-weight edges cannot be walked on and do not contribute to the modularity.
    G->forEdges([&](node u, node v, edgeweight weight) -> void {
        if (u != v && weight > 0) {
            links[positionOf[u]][positionOf[v]].weight += weight;
            links[positionOf[v]][positionOf[u]].weight += weight;
        }
    });
    for (index a = 0; a < n; ++a) {
        for (const auto &[b, link] : links[a]) {
            if (a < b) {
                updateCandidate(a, b, {mergeCost(a, b), link.weight});
            }
        }
    }

    std::vector<double> volume(n);
    for (index i = 0; i < n; ++i) {
        volume[i] = G->weightedDegree(nodes[i], true);
    }
    const double totalWeight = G->totalEdgeWeight();

    // Merges are replayed up to the best one, which is known only when all of them are done.
    std::vector<std::pair<index, index>> merges;
    double modularityGain = 0;
    double bestModularityGain = 0;
    count bestNumberOfMerges = 0;

    while (!candidates.empty()) {
        handler.assureRunning();
        const auto [delta, first, second] = candidates.top();
        candidates.pop();

        const auto candidate = links[first].find(second);
        if (candidate == links[first].end() || candidate->second.delta != delta) {
            continue;
        }

        // Keep the community with more links, so that fewer of them have to be moved.
        const index a = links[first].size() < links[second].size() ? second : first;
        const index b = a == first ? second : first;
        const count sizeA = size[a];
        const count sizeB = size[b];

        modularityGain += candidate->second.weight / totalWeight
                          - volume[a] * volume[b] / (2 * totalWeight * totalWeight);
        merges.emplace_back(a, b);
        if (modularityGain > bestModularityGain) {
            bestModularityGain = modularityGain;
            bestNumberOfMerges = merges.size();
        }

        links[a].erase(b);
        links[b].erase(a);
        for (index i = 0; i < n; ++i) {
            walks[a][i] = (static_cast<double>(sizeA) * walks[a][i]
                           + static_cast<double>(sizeB) * walks[b][i])
                          / static_cast<double>(sizeA + sizeB);
        }
        std::vector<double>().swap(walks[b]);
        size[a] += sizeB;
        volume[a] += volume[b];

        std::map<index, Link> linksOfB;
        linksOfB.swap(links[b]);
        std::vector<std::pair<index, Link>> updates;
        updates.reserve(links[a].size() + linksOfB.size());

        for (const auto &[c, link] : links[a]) {
            const auto shared = linksOfB.find(c);
            if (shared == linksOfB.end()) {
                updates.emplace_back(c, Link{mergeCost(a, c), link.weight});
            } else {
                // Lance-Williams update: the costs known for both communities are enough.
                const double sizeC = static_cast<double>(size[c]);
                const double cost =
                    ((static_cast<double>(sizeA) + sizeC) * link.delta
                     + (static_cast<double>(sizeB) + sizeC) * shared->second.delta - sizeC * delta)
                    / (static_cast<double>(sizeA + sizeB) + sizeC);
                updates.emplace_back(c, Link{cost, link.weight + shared->second.weight});
            }
        }
        for (const auto &[c, link] : linksOfB) {
            links[c].erase(b);
            if (links[a].find(c) == links[a].end()) {
                updates.emplace_back(c, Link{mergeCost(a, c), link.weight});
            }
        }

        for (const auto &[c, link] : updates) {
            updateCandidate(a, c, link);
        }
    }

    UnionFind communities(z);
    for (index i = 0; i < bestNumberOfMerges; ++i) {
        communities.merge(nodes[merges[i].first], nodes[merges[i].second]);
    }
    result = communities.toPartition();
    for (node u = 0; u < z; ++u) {
        if (!G->hasNode(u)) {
            result[u] = none;
        }
    }
    result.compact();

    hasRun = true;
}

} /* namespace NetworKit */
