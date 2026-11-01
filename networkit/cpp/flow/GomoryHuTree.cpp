// networkit-format
#include <algorithm>
#include <limits>
#include <map>
#include <queue>
#include <stdexcept>
#include <utility>
#include <vector>

#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/flow/EdmondsKarp.hpp>
#include <networkit/flow/GomoryHuTree.hpp>

namespace NetworKit {

namespace {

// Returns a copy of G with the same node ids, without self-loops and with parallel edges merged,
// which is what EdmondsKarp expects from undirected graphs.
Graph simpleIndexedCopy(const Graph &G) {
    std::map<std::pair<node, node>, edgeweight> weights;
    G.forEdges([&](node u, node v, edgeweight weight) -> void {
        if (u != v) {
            weights[std::minmax(u, v)] += weight;
        }
    });

    Graph simple(G.upperNodeIdBound(), true);
    for (node u = 0; u < G.upperNodeIdBound(); ++u) {
        if (!G.hasNode(u)) {
            simple.removeNode(u);
        }
    }
    for (const auto &[edge, weight] : weights) {
        simple.addEdge(edge.first, edge.second, weight);
    }
    simple.indexEdges();
    return simple;
}

} // namespace

GomoryHuTree::GomoryHuTree(const Graph &G) : G(&G) {
    if (G.isDirected()) {
        throw std::runtime_error("GomoryHuTree: the graph must be undirected.");
    }
    G.forEdges([&](node, node, edgeweight weight) -> void {
        if (weight < 0) {
            throw std::runtime_error("GomoryHuTree: edge weights must not be negative.");
        }
    });
}

void GomoryHuTree::run() {
    Aux::SignalHandler handler;
    const count z = G->upperNodeIdBound();

    const Graph simple = simpleIndexedCopy(*G);
    std::vector<node> nodes;
    nodes.reserve(simple.numberOfNodes());
    simple.forNodes([&](node u) -> void { nodes.push_back(u); });

    parent.assign(z, none);
    parentWeight.assign(z, 0.);
    depth.assign(z, 0);

    tree = Graph(z, true);
    for (node u = 0; u < z; ++u) {
        if (!G->hasNode(u)) {
            tree.removeNode(u);
        }
    }

    if (nodes.empty()) {
        hasRun = true;
        return;
    }

    const node root = nodes.front();
    for (const node u : nodes) {
        parent[u] = root;
    }

    // Gusfield's algorithm: cut each node off its current parent without contracting the graph.
    std::vector<unsigned char> onSourceSide(z);
    for (const node source : nodes) {
        if (source == root) {
            continue;
        }
        handler.assureRunning();

        const node target = parent[source];
        EdmondsKarp maxFlow(simple, source, target);
        maxFlow.run();

        const edgeweight cutValue = maxFlow.getMaxFlow();
        parentWeight[source] = cutValue;

        std::fill(onSourceSide.begin(), onSourceSide.end(), 0);
        for (const node u : maxFlow.getSourceSet()) {
            onSourceSide[u] = 1;
        }

        for (const node u : nodes) {
            if (u != source && onSourceSide[u] && parent[u] == target) {
                parent[u] = source;
            }
        }

        // The grandparent ended up on the side of source, which thus takes the place of target.
        if (onSourceSide[parent[target]]) {
            parent[source] = parent[target];
            parent[target] = source;
            parentWeight[source] = parentWeight[target];
            parentWeight[target] = cutValue;
        }
    }

    for (const node u : nodes) {
        if (u != root) {
            tree.addEdge(u, parent[u], parentWeight[u]);
        }
    }

    std::queue<node> queue;
    queue.push(root);
    while (!queue.empty()) {
        const node u = queue.front();
        queue.pop();
        tree.forNeighborsOf(u, [&](node v) -> void {
            if (v != parent[u]) {
                depth[v] = depth[u] + 1;
                queue.push(v);
            }
        });
    }

    hasRun = true;
}

edgeweight GomoryHuTree::getMinCutValue(node u, node v) const {
    assureFinished();
    if (!tree.hasNode(u) || !tree.hasNode(v)) {
        throw std::runtime_error("GomoryHuTree: both nodes must be in the graph.");
    }
    if (u == v) {
        throw std::runtime_error("GomoryHuTree: the two nodes must be different.");
    }

    edgeweight minCut = std::numeric_limits<edgeweight>::infinity();
    while (u != v) {
        if (depth[u] < depth[v]) {
            std::swap(u, v);
        }
        minCut = std::min(minCut, parentWeight[u]);
        u = parent[u];
    }
    return minCut;
}

} // namespace NetworKit
