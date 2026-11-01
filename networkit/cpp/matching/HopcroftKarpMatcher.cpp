// networkit-format
#include <numeric>
#include <queue>
#include <stdexcept>
#include <vector>

#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <networkit/matching/HopcroftKarpMatcher.hpp>

namespace NetworKit {

namespace {

// Returns the nodes of one side of the bipartite graph G.
std::vector<node> leftSide(const Graph &G) {
    enum class Side : unsigned char { UNKNOWN, LEFT, RIGHT };

    std::vector<Side> side(G.upperNodeIdBound(), Side::UNKNOWN);
    std::vector<node> left;
    std::queue<node> queue;

    G.forNodes([&](node root) -> void {
        if (side[root] != Side::UNKNOWN) {
            return;
        }
        side[root] = Side::LEFT;
        queue.push(root);
        while (!queue.empty()) {
            const node u = queue.front();
            queue.pop();
            if (side[u] == Side::LEFT) {
                left.push_back(u);
            }
            G.forNeighborsOf(u, [&](node v) -> void {
                if (side[v] == Side::UNKNOWN) {
                    side[v] = side[u] == Side::LEFT ? Side::RIGHT : Side::LEFT;
                    queue.push(v);
                }
            });
        }
    });

    return left;
}

} // namespace

HopcroftKarpMatcher::HopcroftKarpMatcher(const Graph &G) : Matcher(G) {
    if (G.isDirected()) {
        throw std::invalid_argument(
            "The input graph is directed, this algorithm only supports undirected graphs.");
    }
    if (!GraphTools::isBipartite(G)) {
        throw std::invalid_argument(
            "The input graph is not bipartite, this algorithm only supports bipartite graphs.");
    }
}

void HopcroftKarpMatcher::run() {
    Aux::SignalHandler handler;
    constexpr count unreachable = none;

    const count z = G->upperNodeIdBound();
    const std::vector<node> left = leftSide(*G);

    // Neighbors of the left nodes in compressed form, to resume the search of a node in O(1).
    std::vector<index> firstNeighbor(z + 1, 0);
    for (const node u : left) {
        firstNeighbor[u + 1] = G->degree(u);
    }
    std::partial_sum(firstNeighbor.begin(), firstNeighbor.end(), firstNeighbor.begin());
    std::vector<node> neighbors(firstNeighbor.back());
    for (const node u : left) {
        index position = firstNeighbor[u];
        G->forNeighborsOf(u, [&](node v) -> void { neighbors[position++] = v; });
    }

    std::vector<node> mate(z, none);
    std::vector<count> level(z);
    std::vector<index> nextNeighbor(z);
    std::vector<node> chosenRight(z);
    std::vector<node> queue, path;
    queue.reserve(left.size());

    while (true) {
        handler.assureRunning();

        // BFS from the free left nodes: compute the length of the shortest augmenting paths, and
        // the levels they go through.
        queue.clear();
        for (const node u : left) {
            if (mate[u] == none) {
                level[u] = 0;
                queue.push_back(u);
            } else {
                level[u] = unreachable;
            }
        }

        count lastLevel = unreachable;
        for (index i = 0; i < queue.size() && level[queue[i]] < lastLevel; ++i) {
            const node u = queue[i];
            for (index j = firstNeighbor[u]; j < firstNeighbor[u + 1]; ++j) {
                const node next = mate[neighbors[j]];
                if (next == none) {
                    lastLevel = level[u];
                } else if (level[next] == unreachable) {
                    level[next] = level[u] + 1;
                    queue.push_back(next);
                }
            }
        }

        if (lastLevel == unreachable) {
            break;
        }

        // Augment along a maximal set of disjoint shortest paths. Explicit stack, since paths can
        // be as long as the graph. A failed search makes its node unreachable, and all searches
        // resume from where the last one stopped: the total work of this phase is in O(m).
        for (const node u : left) {
            nextNeighbor[u] = firstNeighbor[u];
        }

        for (const node root : left) {
            if (mate[root] != none) {
                continue;
            }

            path.assign(1, root);
            while (!path.empty()) {
                const node u = path.back();
                if (nextNeighbor[u] == firstNeighbor[u + 1]) {
                    level[u] = unreachable;
                    path.pop_back();
                    continue;
                }

                const node right = neighbors[nextNeighbor[u]++];
                const node next = mate[right];
                chosenRight[u] = right;
                if (next == none) {
                    for (const node v : path) {
                        mate[v] = chosenRight[v];
                        mate[chosenRight[v]] = v;
                    }
                    path.clear();
                } else if (level[next] == level[u] + 1 && level[next] <= lastLevel) {
                    path.push_back(next);
                }
            }
        }
    }

    M = Matching(z);
    for (const node u : left) {
        if (mate[u] != none) {
            M.match(u, mate[u]);
        }
    }

    hasRun = true;
}

} /* namespace NetworKit */
