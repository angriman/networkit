/*
 * ApproxEffectiveResistance.cpp
 *
 *  Created on: 17.10.2019
 *      Author: Eugenio Angriman <angrimae@hu-berlin.de>
 */

// networkit-format

#include <omp.h>
#include <queue>

#include <networkit/auxiliary/Timer.hpp>
#include <networkit/centrality/ApproxEffectiveResistance.hpp>

namespace NetworKit {

ApproxEffectiveResistance::ApproxEffectiveResistance(const Graph &G, double inputEpsilon)
    : G(G), epsilon(0.7 * inputEpsilon), delta(1.0 / static_cast<double>(G.numberOfNodes())),
      bcc(new BiconnectedComponents(G)) {
    if (G.isDirected()) {
        throw std::runtime_error("Error: the input graph must be undirected.");
    }
    if (G.isWeighted()) {
        throw std::runtime_error("Error: the input graph must be unweighted.");
    }
    if (G.numberOfNodes() < 2) {
        throw std::runtime_error("Error: the graph should have at leasts two vertices");
    }

    const auto n = G.upperNodeIdBound();
    statusGlobal.resize(omp_get_max_threads(), std::vector<NodeStatus>(n, NodeStatus::NOT_VISITED));
    parentGlobal.resize(omp_get_max_threads(), std::vector<small_node>(n, inf));
    approxEffResistanceGlobal.resize(omp_get_max_threads(), std::vector<int>(n));
    tVisitGlobal.resize(omp_get_max_threads(), std::vector<uint32_t>(n));
    tFinishGlobal.resize(omp_get_max_threads(), std::vector<uint32_t>(n));
    generators.reserve(omp_get_max_threads());
    for (omp_index i = 0; i < omp_get_max_threads(); ++i) {
        generators.emplace_back((Aux::Random::integer(), Aux::Random::integer()));
        generators.back().seed(Aux::Random::integer());
    }

    ustChildPtrGlobal.resize(omp_get_max_threads(), std::vector<small_node>(n));
    ustSiblingPtrGlobal.resize(omp_get_max_threads(), std::vector<small_node>(n));
    bfsParent.resize(n, inf);
    samplingTime.resize(omp_get_max_threads(), 0);
    dfsTime.resize(omp_get_max_threads(), 0);
    aggregationTime.resize(omp_get_max_threads(), 0);
}

void ApproxEffectiveResistance::init() {
    computeNodeSequence();
    computeBFSTree();
    didInit = true;
}

small_node ApproxEffectiveResistance::approxMinEccNode() {
    auto &status = statusGlobal[0];
    std::vector<uint32_t> distance(G.upperNodeIdBound());
    std::vector<uint32_t> eccLowerBound(G.upperNodeIdBound());

    auto maxDegreeNode = [&]() {
        small_node maxDegNode = 0;
        count maxDeg = 0;
        G.forNodes([&](const small_node u) {
            const auto degU = G.degree(u);
            if (degU > maxDeg) {
                maxDeg = degU;
                maxDegNode = u;
            }
        });
        return maxDegNode;
    };

    auto doBFS = [&](const small_node source) {
        std::queue<small_node> q;
        q.push(source);
        status[source] = NodeStatus::VISITED;
        distance[source] = 0;
        small_node farthest = 0;

        do {
            const auto u = q.front();
            q.pop();
            eccLowerBound[u] = std::max(eccLowerBound[u], distance[u]);
            farthest = u;

            G.forNeighborsOf(u, [&](const small_node v) {
                if (status[v] == NodeStatus::NOT_VISITED) {
                    q.push(v);
                    status[v] = NodeStatus::VISITED;
                    distance[v] = distance[u] + 1;
                }
            });

        } while (!q.empty());

        std::fill(status.begin(), status.end(), NodeStatus::NOT_VISITED);
        return farthest;
    };

    small_node source = maxDegreeNode();

    for (int i = 0; i < sweeps; ++i) {
        source = doBFS(source);
    }

    // Return node with minimum ecc lower bound
    return static_cast<small_node>(std::min_element(eccLowerBound.begin(), eccLowerBound.end())
                                   - eccLowerBound.begin());
}

void ApproxEffectiveResistance::computeNodeSequence() {
    // We use thread 0's vector
    auto &status = statusGlobal[0];

    // Compute the biconnected components
    auto &bcc_ = *(bcc.get());
    bcc_.run();
    auto components = bcc_.getComponents();

    std::queue<small_node> queue;
    std::vector<small_node> curSequence;

    for (const auto &curComponent : components) {
        // Biconnected components with 3 nodes or less can be handled trivially.
        if (curComponent.size() == 2) {
            sequences.push_back({(small_node)curComponent[0], (small_node)curComponent[1]});
            continue;
        }

        // We take the node with highest degree in the component as source
        auto source = curComponent[0];
        for (const auto u : curComponent) {
            source = G.degree(u) > G.degree(source) ? u : source;
            status[u] = NodeStatus::VISITED;
        }

        // Start a BFS from source to determine the order of the nodes.
        // Later, we need to re-explore the graph again, so in the beginning we mark all nodes to be
        // visited as VISITED, and we set them as NOT_VISITED during the BFS. This avoids a call to
        // std::fill.
        queue.push(source);
        status[source] = NodeStatus::NOT_VISITED;

        do {
            const auto u = queue.front();
            queue.pop();
            curSequence.push_back(u);
            G.forNeighborsOf(u, [&](const small_node v) {
                if (status[v] == NodeStatus::VISITED) {
                    status[v] = NodeStatus::NOT_VISITED;
                    queue.push(v);
                }
            });
        } while (!queue.empty());

        sequences.push_back(std::move(curSequence));
    }

    // Set the root to the highest degree node within the largest biconnected component
    root = approxMinEccNode();

    biAnchor.resize(bcc_.numberOfComponents(), inf);
    biParent.resize(bcc_.numberOfComponents(), inf);

#ifndef NDEBUG
    G.forNodes([&](const small_node u) { assert(status[u] == NodeStatus::NOT_VISITED); });
#endif

    // Topological order of biconnected components: tree of biconnected components starting from the
    // root's biconencted component. If the root is in multiple biconencted components, we take one
    // of them arbitrarily select one of them.
    std::queue<std::pair<small_node, index>> q;
    const auto &rootComps = bcc_.getComponentsOfNode(root);
    q.push({root, *(rootComps.begin())});

    topOrder.reserve(bcc_.numberOfComponents());
    topOrder.insert(topOrder.begin(), rootComps.begin(), rootComps.end());

    std::vector<uint32_t> distance(G.upperNodeIdBound());
    status[root] = NodeStatus::VISITED;
    rootEcc = 0;
    do {
        const auto front = q.front();
        q.pop();
        G.forNeighborsOf(front.first, [&](const small_node v) {
            if (status[v] == NodeStatus::NOT_VISITED) {
                distance[v] = distance[front.first] + 1;
                rootEcc = std::max(rootEcc, distance[v]);
                const auto &vComps = bcc_.getComponentsOfNode(v);
                for (const auto vComponentIndex : vComps) {
                    // Check if a new biconnected components has been found.
                    if (vComponentIndex != front.second && biAnchor[vComponentIndex] == inf
                        && rootComps.find(vComponentIndex) == rootComps.end()) {
                        // The anchor cannot be the root, because the anchor does not have a parent.
                        // We handle biAnchor = inf cases later.
                        biAnchor[vComponentIndex] = (v == root) ? inf : v;
                        biParent[vComponentIndex] = front.second;
                        topOrder.push_back(vComponentIndex);
                    }
                }
                q.push({v, *(vComps.begin())});
                status[v] = NodeStatus::VISITED;
            }
        });
    } while (!q.empty());

    INFO("Root eccentricity = ", rootEcc);
    INFO("Root = ", root);

#ifndef NDEBUG
    G.forNodes([&](const small_node u) { assert(status[u] == NodeStatus::VISITED); });
    assert(topOrder.size() == bcc_.numberOfComponents());
    assert(std::unordered_set<index>(topOrder.begin(), topOrder.end()).size()
           == bcc_.numberOfComponents());
#endif
}

void ApproxEffectiveResistance::computeBFSTree() {
    // Using thread 0's vector
    auto &status = statusGlobal[0];
    std::fill(status.begin(), status.end(), NodeStatus::NOT_VISITED);

    std::queue<small_node> queue;
    queue.push(root);
    status[root] = NodeStatus::VISITED;

    do {
        const auto currentNode = queue.front();
        queue.pop();
        G.forNeighborsOf(currentNode, [&](const small_node v) {
            if (status[v] == NodeStatus::NOT_VISITED) {
                status[v] = NodeStatus::VISITED;
                queue.push(v);
                bfsParent[v] = currentNode;
            }
        });
    } while (!queue.empty());

#ifndef NDEBUG
    checkBFSTree();
#endif
}

void ApproxEffectiveResistance::sampleUST() {
    // Getting thread-local vectors
    auto &status = statusGlobal[omp_get_thread_num()];
    auto &parent = parentGlobal[omp_get_thread_num()];
    auto &childPtr = ustChildPtrGlobal[omp_get_thread_num()];
    auto &siblingPtr = ustSiblingPtrGlobal[omp_get_thread_num()];
    std::fill(status.begin(), status.end(), NodeStatus::NOT_IN_COMPONENT);
    std::fill(parent.begin(), parent.end(), inf);
    std::fill(childPtr.begin(), childPtr.end(), inf);
    std::fill(siblingPtr.begin(), siblingPtr.end(), inf);

    auto &generator = generators[omp_get_thread_num()];

    // Iterate over the biconnected components in their topological order.
    for (const auto componentIndex : topOrder) {
        // Current component, sorted by vertex degree.
        const auto &sequence = sequences[componentIndex];
        auto curAnchor = biAnchor[componentIndex];

        // Finds the parent of the current anchor node i.e., the anchor's neighbor that is in the
        // parent component.
        auto updateParentOfAnchor = [&]() {
            for (const auto v : G.neighborRange(curAnchor)) {
                const auto &vComps = bcc->getComponentsOfNode(v);
                if (vComps.find(biParent[componentIndex]) != vComps.end()) {
                    parent[curAnchor] = v;
                    break;
                }
            }
            assert(parent[curAnchor] != inf);
        };

        if (sequence.size() == 2) {
            // Happens when the current component is the root component in the topological
            // order. In this case, the root plays the anchor's role.
            if (curAnchor == inf) {
                const auto v = (sequence.front() == root) ? sequence[1] : sequence[0];
                assert(sequence.front() == root || sequence.back() == root);
                assert(v != root);
                parent[v] = root;
            } else {
                const auto v = (sequence[0] == curAnchor) ? sequence[1] : sequence[0];
                assert(v != curAnchor);
                assert(v != root);
                parent[v] = curAnchor;
                if (parent[curAnchor] == inf) {
                    updateParentOfAnchor();
                }
            }
#ifndef NDEBUG
            checkTwoNodesSequence(sequence);
#endif
            continue;
        }

        // We start building the spanning tree from the first node of the
        // sequence.
        // Root of the spanning tree
        status[sequence[0]] = NodeStatus::IN_SPANNING_TREE;
        const auto curAnchorParent = (curAnchor != inf) ? parent[curAnchor] : inf;

        // All the remaining nodes in the components need to be visited.
        std::for_each(sequence.begin() + 1, sequence.end(),
                      [&status](const small_node u) { status[u] = NodeStatus::NOT_VISITED; });

        uint32_t nodesInSpanningTree = 1;
        // Iterate over the remaining nodes to create the spanning tree.
        for (auto it = sequence.begin() + 1; it != sequence.end(); ++it) {
            const small_node walkStart = *it;
            if (status[walkStart] == NodeStatus::IN_SPANNING_TREE) {
                // Node already added to the spanning tree
                continue;
            }
            small_node currentNode = walkStart;

            // Start a new random walk from the current node.
            do {
                // Get a random neighbor within the component
                small_node randomNeighbor;
                do {
                    randomNeighbor =
                        G.getIthNeighbor(currentNode, generator.nextUInt(G.degree(currentNode)));
                } while (status[randomNeighbor] == NodeStatus::NOT_IN_COMPONENT);

                assert(randomNeighbor != inf);
                parent[currentNode] = randomNeighbor;
                currentNode = randomNeighbor;

            } while (status[currentNode] != NodeStatus::IN_SPANNING_TREE);

            // Last node encountered in the random walk (it is in the spanning tree);
            const auto walkEnd = currentNode;
            assert(status[walkEnd] == NodeStatus::IN_SPANNING_TREE);
            // Add the random walk to the spanning tree; eventually, reverse the path if the
            // anchor/root is encountered
            for (currentNode = walkStart; currentNode != walkEnd;
                 currentNode = parent[currentNode]) {

                status[currentNode] = NodeStatus::IN_SPANNING_TREE;
                ++nodesInSpanningTree;
                if (currentNode == curAnchor || currentNode == root) {

                    // Anchor of current component in the walk, we have to reverse the
                    // parent pointers
                    auto next = parent[currentNode];
                    auto nextParent = inf;
                    auto tmp = currentNode;
                    do {
                        status[next] = NodeStatus::IN_SPANNING_TREE;
                        nextParent = parent[next];
                        parent[next] = currentNode;
                        currentNode = next;
                        next = nextParent;
                    } while (next != inf);

                    if (tmp == root)
                        parent[root] = inf;
                    else if (parent[curAnchor] == inf) {
                        // Not inf if articulation node visited by the parent component before
                        updateParentOfAnchor();
                    } else
                        parent[curAnchor] = curAnchorParent;

                    break;
                }
            }

            if (nodesInSpanningTree == sequence.size())
                break;
        }

        for (const auto u : sequence)
            status[u] = NodeStatus::NOT_IN_COMPONENT;
    }

    uint32_t visitedNodes = 0;
    for (small_node u : G.nodeRange()) {
        while (status[u] == NodeStatus::NOT_IN_COMPONENT) {
            status[u] = NodeStatus::NOT_VISITED;
            ++visitedNodes;
            small_node parentU = parent[u];
            if (parent[u] != inf) {
                assert(siblingPtr[u] == inf);
                if(childPtr[parentU] != inf)
                    siblingPtr[u] = childPtr[parentU];
                childPtr[parentU] = u;
                u = parentU;
            } else
                break;
        }
        if (visitedNodes == G.numberOfNodes())
            break;
    }

#ifndef NDEBUG
    checkUST();
#endif
}

void ApproxEffectiveResistance::dfsUST() {
    auto &tVisit = tVisitGlobal[omp_get_thread_num()];
    auto &tFinish = tFinishGlobal[omp_get_thread_num()];
    const auto &childPtr = ustChildPtrGlobal[omp_get_thread_num()];
    const auto &siblingPtr = ustSiblingPtrGlobal[omp_get_thread_num()];

    auto &status = statusGlobal[omp_get_thread_num()];

    std::stack<std::pair<small_node, small_node>> stack;
    stack.push({root, childPtr[root]});

    uint32_t timestamp = 0;
    do {
        // v is a child of u that has not been visited yet.
        const small_node u = stack.top().first;
        const small_node v = stack.top().second;

        if (v == inf) {
            stack.pop();
            tFinish[u] = ++timestamp;
        } else {
            stack.top().second = siblingPtr[v];
            tVisit[v] = ++timestamp;
            stack.push({v, childPtr[v]});
            assert(parentGlobal[omp_get_thread_num()][v] == u);
        }
    } while (!stack.empty());
}

void ApproxEffectiveResistance::aggregateUST() {
#ifndef NDEBUG
    checkTimeStamps();
#endif

    auto &approxEffResistance = approxEffResistanceGlobal[omp_get_thread_num()];
    const auto &tVisit = tVisitGlobal[omp_get_thread_num()];
    const auto &tFinish = tFinishGlobal[omp_get_thread_num()];
    const auto &parent = parentGlobal[omp_get_thread_num()];

    // Doing aggregation
    G.forNodes([&](const small_node u) {
        small_node p = bfsParent[u];
        small_node c = u;
        auto goUp = [&]() {
            c = p;
            p = bfsParent[p];
        };
        while (p != inf) {
            // Edge in BSTree: e1 -> e2
            small_node e1 = p, e2 = c;
            bool reverse = false;
            if (e1 != parent[e2]) {
                if (e2 != parent[e1]) {
                    goUp();
                    continue;
                }
                std::swap(e1, e2);
                reverse = true;
            }

            if (tVisit[u] >= tVisit[e2] && tFinish[u] <= tFinish[e2])
                approxEffResistance[u] += reverse ? -1 : 1;

            goUp();
        }
    });
}

void ApproxEffectiveResistance::run() {
    numberOfUSTs = computeNumberOfUSTs();
#pragma omp parallel for
    for (omp_index i = 0; i < static_cast<omp_index>(numberOfUSTs); ++i) {
        Aux::Timer timer;
        timer.start();
        sampleUST();
        timer.stop();
        samplingTime[omp_get_thread_num()] += timer.elapsedNanoseconds();
        timer.start();
        dfsUST();
        timer.stop();
        dfsTime[omp_get_thread_num()] += timer.elapsedNanoseconds();
        // Update effective resistance values using sampled UST
        timer.start();
        aggregateUST();
        timer.stop();
        aggregationTime[omp_get_thread_num()] += timer.elapsedNanoseconds();
    }

    timeToSample = 0, timeDFS = 0, timeToAggregate = 0;
    for (index i = 0; i < omp_get_max_threads(); ++i) {
        timeToSample += ((double)samplingTime[i]) / 1e9;
        timeDFS += ((double)dfsTime[i]) / 1e9;
        timeToAggregate += ((double)aggregationTime[i]) / 1e9;
    }
    hasRun = true;
}

#ifndef NDEBUG
/*
 * Methods for sanity check.
 */
void ApproxEffectiveResistance::checkUST() const {
    std::vector<bool> visitedNodes(G.upperNodeIdBound());
    const auto &parent = parentGlobal[omp_get_thread_num()];
    // To debug
    G.forNodes([&](small_node u) {
        if (u == root) {
            assert(parent[u] == inf);
        } else {
            assert(parent[u] != inf);
            std::fill(visitedNodes.begin(), visitedNodes.end(), 0);
            visitedNodes[u] = 1;
            auto start = u;
            do {
                u = parent[u];
                assert(!visitedNodes[u]);
                visitedNodes[u] = 1;
            } while (u != root);
        }
    });
}

void ApproxEffectiveResistance::checkBFSTree() const {
    G.forNodes([&](small_node u) {
        if (u == root) {
            assert(bfsParent[u] == inf);
        } else {
            std::vector<bool> visited(G.upperNodeIdBound());
            visited[u] = true;
            do {
                u = bfsParent[u];
                if (visited[u]) {
                    WARN("Loop error");
                }
                assert(!visited[u]);
            } while (u != root);
        }
    });
}

void ApproxEffectiveResistance::checkTwoNodesSequence(
    const std::vector<small_node> &sequence) const {
    const std::vector<small_node> &parent = parentGlobal[omp_get_thread_num()];
    for (small_node u : sequence) {
        if (u == root) {
            assert(parent[u] == inf);
        } else {
            std::vector<bool> visited(G.upperNodeIdBound());
            visited[u] = true;
            do {
                u = parent[u];
                assert(!visited[u]);
                visited[u] = true;
            } while (u != root);
        }
    }
}

void ApproxEffectiveResistance::checkTimeStamps() const {
    const auto &tVisit = tVisitGlobal[omp_get_thread_num()];
    const auto &tFinish = tFinishGlobal[omp_get_thread_num()];
    G.forNodes([&](const small_node u) {
        assert(tVisit[u] < tFinish[u]);
        if (u == root)
            assert(tVisit[u] == 0);
        else
            assert(tVisit[u] > 0);
        assert(tVisit[u] < 2 * G.numberOfNodes());
    });
}

#endif // NDEBUG

} // namespace NetworKit
