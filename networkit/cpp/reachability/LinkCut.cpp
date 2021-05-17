#include <queue>
#include <set>
#include <stdexcept>
#include <unordered_set>

#include <networkit/auxiliary/Random.hpp>
#include <networkit/reachability/LinkCut.hpp>

namespace NetworKit {

LinkCut::LinkCut(const Graph &G) : G(&G) {
    if (!G.hasEdgeIds())
        throw std::runtime_error("Edges are not indexed.");

    pathLength.resize(G.upperNodeIdBound());
    edgeInST.resize(G.upperEdgeIdBound());

    parent.resize(G.upperNodeIdBound());
    nonSTEdges.resize(G.numberOfEdges() - G.numberOfNodes() + 1);
    seenNodes.reserve(G.numberOfNodes());

    edgeScores.resize(G.upperEdgeIdBound());
    root = 0;
    count maxDeg = 0;
    G.forNodes([&](node u) {
        if (G.degree(u) > maxDeg) {
            maxDeg = G.degree(u);
            root = u;
        }
    });

    buildNodeSequence();
    sampleUST();
}

void LinkCut::buildNodeSequence() {
    status.resize(G->upperNodeIdBound());
    sequence.clear();
    sequence.reserve(G->numberOfNodes());

    std::queue<node> q;
    q.push(root);
    status[root] = 1;
    sequence.push_back(root);

    do {
        const auto front = q.front();
        q.pop();
        G->forNeighborsOf(front, [&](node v) {
            if (!status[v]) {
                q.push(v);
                status[v] = 1;
                sequence.push_back(v);
            }
        });
    } while (!q.empty());
    assert(sequence.size() == G->numberOfNodes());
}

void LinkCut::sampleUST() {
    std::unordered_set<edgeid> stEdges;
    std::vector<std::pair<node, edgeid>> stack;
    G->parallelForNodes([&](const node u) { status[u] = 0; });
    status[root] = 1;
    parent[root] = none;

    for (count i = 1; i < sequence.size(); ++i) {
        const auto u = sequence[i];
        if (status[u] == 1) {
            continue;
        }
        status[u] = 2;
        stack.push_back({u, none});

        do {
            // Get a random neighbor
            std::pair<node, edgeid> randomNeighbor;
            randomNeighbor = G->randomNeighborWithId(stack.back().first);

            const auto v = randomNeighbor.first;

            // Spanning tree reached, adding the current random walk to the spanning tree
            if (status[v] == 1) {
                stack.push_back(randomNeighbor);
                break;
            }

            // Loop found, discarding it
            if (status[v] == 2) {
                auto rit = stack.rbegin();
                for (; rit->first != v; ++rit) {
                    status[rit->first] = 0;
                }
                stack.resize(stack.size() - (rit - stack.rbegin()));
            } else { // Non-visited node found, add it to the random walk
                status[v] = 2;
                stack.push_back(randomNeighbor);
            }
        } while (true);

        // Add the random walk to the spanning tree
        for (auto it = stack.begin() + 1; it != stack.end(); ++it) {
            status[it->first] = 1;
            parent[(it - 1)->first] = it->first;
            stEdges.insert(it->second);
        }
        status[stack[0].first] = 1;
        stack.clear();
    }

    count nonSTEdgeCount = 0;
    count stEdgeCount = 0;
    G->forEdges([&](node u, node v, edgeweight, edgeid eid) {
        if (stEdges.find(eid) == stEdges.end()) {
            nonSTEdges[nonSTEdgeCount++] = {u, v, eid};
            edgeInST[eid] = 0;
        } else {
            edgeInST[eid] = 1;
            ++stEdgeCount;
        }
    });
    assert(nonSTEdgeCount == nonSTEdges.size());
    assert(stEdgeCount == G->numberOfNodes() - 1);
}

void LinkCut::doLinkCut() {
    // Index of random edge outside spanning tree to be added
    const auto randInt = Aux::Random::integer(nonSTEdges.size() - 1);
    // Actual edge to be added
    const auto &randNonSTEdge = nonSTEdges[randInt];
    assert(!edgeInST[randNonSTEdge.eid]);
    node u = randNonSTEdge.u, u_ = randNonSTEdge.u;
    node v = randNonSTEdge.v, v_ = randNonSTEdge.v;
    seenNodes.push_back(u);
    seenNodes.push_back(v);

    auto processNode = [&](node &x, const char *s) {
        const auto xParent = parent[x];
        if (x == root) {
            return 1;
        }
        if (pathLength[xParent]) {
            return 2;
        }
        pathLength[xParent] = pathLength[x] + 1;
        seenNodes.push_back(xParent);
        x = xParent;
        return 0;
    };

    pathLength[u] = pathLength[v] = 1;

    // Finding lca
    while (processNode(u, "u") + processNode(v, "v") < 2) {
    }

    if (u == root || !pathLength[parent[u]]) {
        std::swap(u, v);
        std::swap(u_, v_);
    }

    const auto lca = parent[u];
    assert(u != root);
    assert(pathLength[parent[u]]);

    node pred = v_;

    // Here we omit the additional edge, the actual loop is 1 hop longer
    const auto loopLength = pathLength[u] + pathLength[lca] - 1;
    auto randomEdgeInLoop = Aux::Random::integer(loopLength - 1);
    if (randomEdgeInLoop >= pathLength[lca] - 1) {
        randomEdgeInLoop -= pathLength[lca] - 1;
        u = u_;
    } else {
        u = v_;
        pred = u_;
    }

    node tmp = u;
    do {
        assert(u != root);
        u = parent[u];
        parent[tmp] = pred;
        pred = tmp;
        tmp = u;
    } while (randomEdgeInLoop--);

    // Add edge to tree
    assert(!edgeInST[randNonSTEdge.eid]);
    edgeInST[randNonSTEdge.eid] = 1;
    // Update set of non-tree edges
    nonSTEdges[randInt] = {pred, u, G->edgeId(pred, u)}; // TODO get rid of edgeId
    // Remove edge from the tree
    assert(edgeInST[randNonSTEdge.eid]);
    edgeInST[randNonSTEdge.eid] = 0;

    for (auto seenNode : seenNodes) {
        pathLength[seenNode] = 0;
    }
    seenNodes.clear();
}

void LinkCut::doCuts(count n) {
    for (count i = 0; i < n; ++i)
        doLinkCut();
}

} // namespace NetworKit
