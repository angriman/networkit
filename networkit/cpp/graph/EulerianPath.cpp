// networkit-format
#include <algorithm>
#include <numeric>
#include <utility>
#include <vector>

#include <networkit/graph/EulerianPath.hpp>

namespace NetworKit {

namespace {

// Returns the node an Eulerian path has to start from, or none if the degrees rule out such a path.
// outDegree - inDegree must be 1 for the start, -1 for the end, and 0 for all the other nodes. In
// undirected graphs, there are at most two nodes with odd degree and any of them can be the start.
node findStart(const Graph &G, const std::vector<count> &outDegree,
               const std::vector<count> &inDegree) {
    node start = none, anyWithEdges = none;
    count numberOfStarts = 0, numberOfEnds = 0;
    bool degreesAreValid = true;

    G.forNodes([&](node u) -> void {
        if (anyWithEdges == none && outDegree[u] + inDegree[u] > 0) {
            anyWithEdges = u;
        }

        if (G.isDirected()) {
            if (outDegree[u] == inDegree[u] + 1) {
                ++numberOfStarts;
                start = u;
            } else if (inDegree[u] == outDegree[u] + 1) {
                ++numberOfEnds;
            } else if (outDegree[u] != inDegree[u]) {
                degreesAreValid = false;
            }
        } else if (outDegree[u] % 2 == 1) {
            ++numberOfStarts;
            if (start == none) {
                start = u;
            }
        }
    });

    if (!degreesAreValid) {
        return none;
    }

    if (G.isDirected() ? numberOfStarts != numberOfEnds || numberOfStarts > 1
                       : numberOfStarts != 0 && numberOfStarts != 2) {
        return none;
    }

    return start == none ? anyWithEdges : start;
}

} // namespace

EulerianPath::EulerianPath(const Graph &G) : G(&G) {}

void EulerianPath::run() {
    const count z = G->upperNodeIdBound();
    const count m = G->numberOfEdges();
    const bool directed = G->isDirected();

    path.clear();
    exists = true;

    if (m == 0) {
        hasRun = true;
        return;
    }

    // A self-loop adds two to the degree of an undirected node, but needs only one adjacency entry.
    std::vector<count> outDegree(z, 0), inDegree(z, 0);
    std::vector<index> firstEntry(z + 1, 0);
    G->forEdges([&](node u, node v) -> void {
        ++outDegree[u];
        ++firstEntry[u + 1];
        if (directed) {
            ++inDegree[v];
        } else {
            ++outDegree[v];
            if (u != v) {
                ++firstEntry[v + 1];
            }
        }
    });

    const node start = findStart(*G, outDegree, inDegree);
    if (start == none) {
        exists = false;
        hasRun = true;
        return;
    }

    std::partial_sum(firstEntry.begin(), firstEntry.end(), firstEntry.begin());

    // Each edge has an index, so that the two adjacency entries of an undirected edge share a flag.
    std::vector<std::pair<node, index>> entries(firstEntry.back());
    std::vector<index> nextEntry(firstEntry.begin(), firstEntry.end() - 1);
    index edgeIndex = 0;
    G->forEdges([&](node u, node v) -> void {
        entries[nextEntry[u]++] = {v, edgeIndex};
        if (!directed && u != v) {
            entries[nextEntry[v]++] = {u, edgeIndex};
        }
        ++edgeIndex;
    });
    std::copy(firstEntry.begin(), firstEntry.end() - 1, nextEntry.begin());

    std::vector<unsigned char> used(m, 0);
    std::vector<node> stack{start};
    path.reserve(m + 1);

    while (!stack.empty()) {
        const node u = stack.back();
        while (nextEntry[u] < firstEntry[u + 1] && used[entries[nextEntry[u]].second]) {
            ++nextEntry[u];
        }

        if (nextEntry[u] == firstEntry[u + 1]) {
            path.push_back(u);
            stack.pop_back();
        } else {
            const auto [v, edge] = entries[nextEntry[u]++];
            used[edge] = 1;
            stack.push_back(v);
        }
    }

    // The degrees allow a path, but if the edges are not connected some of them stay unreachable.
    if (path.size() != m + 1) {
        path.clear();
        exists = false;
    }

    std::reverse(path.begin(), path.end());
    hasRun = true;
}

} /* namespace NetworKit */
