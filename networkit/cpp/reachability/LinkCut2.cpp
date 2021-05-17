#include <algorithm>
#include <iostream>
#include <queue>
#include <set>
#include <stdexcept>
#include <unordered_set>

#include <networkit/auxiliary/Random.hpp>
#include <networkit/reachability/LinkCut2.hpp>

namespace NetworKit {

LinkCut2::LinkCut2(const Graph &G) : G(&G) {
    if (!G.hasEdgeIds())
        throw std::runtime_error("Edges are not indexed.");

    pathLength.resize(G.upperNodeIdBound());
    edgeInST.resize(G.upperEdgeIdBound());

    parent.resize(G.upperNodeIdBound());
    nonSTEdges.resize(G.numberOfEdges() - G.numberOfNodes() + 1);
    seenNodes.reserve(G.numberOfNodes());

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

void LinkCut2::buildNodeSequence() {
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

void LinkCut2::sampleUST() {
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
            edgeInST[eid] = false;
        } else {
            edgeInST[eid] = true;
            ++stEdgeCount;
        }
    });
    assert(nonSTEdgeCount == nonSTEdges.size());
    assert(stEdgeCount == G->numberOfNodes() - 1);
}

void LinkCut2::doLinkCut() {
    // Remove a random edge
    auto edgeToRemove = Aux::Random::integer(G->numberOfNodes() - 2);
    index removedEdge = 0;
    for (index i = 0; i < G->numberOfEdges(); ++i) {
        if (edgeInST[i]) {
            if (edgeToRemove == 0) {
                removedEdge = i;
                break;
            } else
                --edgeToRemove;
        }
    }

    assert(edgeInST[removedEdge]);
    edgeInST[removedEdge] = false;

    // BFS from root
    std::vector<bool> reachableFromRoot(G->numberOfNodes());
    reachableFromRoot[root] = true;
    std::queue<node> queue;
    queue.push(root);

    do {
        const auto u = queue.front();
        queue.pop();

        G->forEdgesOf(u, [&](node, node v, edgeweight, edgeid eid) {
            if (!reachableFromRoot[v] && edgeInST[eid]) {
                reachableFromRoot[v] = true;
                queue.push(v);
            }
        });
    } while (!queue.empty());

    std::vector<edgeid> cutEdges;
    G->forEdges([&](node u, node v, edgeweight, edgeid eid) {
        if (reachableFromRoot[u] != reachableFromRoot[v])
            cutEdges.push_back(eid);
    });

    assert(!cutEdges.empty());
#ifndef NDEBUG
    for (edgeid eid : cutEdges)
        assert(!edgeInST[eid]);
#endif
    edgeInST[cutEdges[Aux::Random::integer(cutEdges.size() - 1)]] = true;
}

std::vector<double> LinkCut2::simulation(count reps, count cutsPerRep) {
    std::vector<double> result(G->numberOfEdges());
    for (count i = 0; i < reps; ++i) {
        for (count j = 0; j < cutsPerRep; ++j)
            doLinkCut();
        assert(std::count_if(edgeInST.begin(), edgeInST.end(), [](bool elem) { return elem; })
               == G->numberOfNodes() - 1);
        G->parallelForEdges([&](node, node, edgeweight, edgeid eid) {
            if (edgeInST[eid])
                ++result[eid];
        });
    }

    G->parallelForEdges(
        [&](node, node, edgeweight, edgeid eid) { result[eid] /= static_cast<double>(reps); });

    return result;
}

} // namespace NetworKit
