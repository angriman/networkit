/*
 * ApproxEffectiveResistance.cpp
 *
 *  Created on: 17.10.2019
 *      Author: Eugenio Angriman <angrimae@hu-berlin.de>
 */

#include <omp.h>
#include <queue>

#include <networkit/auxiliary/Timer.hpp>
#include <networkit/centrality/ApproxEffectiveResistance.hpp>

static constexpr bool isDebug = false;

namespace NetworKit {

ApproxEffectiveResistance::ApproxEffectiveResistance(const Graph &G, double epsilon,
                                                     double tolerance)
    : G(G), epsilon(epsilon), delta(1.0 / static_cast<double>(G.numberOfNodes())),
      tolerance(tolerance), bcc(new BiconnectedComponents(G)), result(G.upperNodeIdBound()) {
    if (G.isDirected()) {
        throw std::runtime_error("Error: the input graph must be undirected.");
    }
    if (G.isWeighted()) {
        throw std::runtime_error("Error: the input graph must be unweighted.");
    }
    if (G.numberOfNodes() < 2) {
        throw std::runtime_error("Error: the graph should have at leasts two vertices");
    }

    Aux::Timer timer;
    timer.start();
    const auto n = G.upperNodeIdBound();
    statusGlobal.resize(omp_get_max_threads(), std::vector<NodeStatus>(n, NodeStatus::NOT_VISITED));
    parentGlobal.resize(omp_get_max_threads(), std::vector<node>(n, none));
    approxEffResistanceGlobal.resize(omp_get_max_threads(), std::vector<int>(n));
    tVisitGlobal.resize(omp_get_max_threads(), std::vector<uint32_t>(n));
    tFinishGlobal.resize(omp_get_max_threads(), std::vector<uint32_t>(n));
    generators.reserve(omp_get_max_threads());
    for (omp_index i = 0; i < omp_get_max_threads(); ++i) {
        generators.emplace_back((Aux::Random::integer(), Aux::Random::integer()));
        generators.back().seed(Aux::Random::integer());
    }

    bfsParent.resize(n, none);
    diagonal.resize(n);

    timer.stop();
    time["Initialization"] = static_cast<double>(timer.elapsedMicroseconds()) / 1e6;
}

void ApproxEffectiveResistance::init() {
    computeNodeSequence();
    computeBFSTree();
    didInit = true;
}

node ApproxEffectiveResistance::approxMinEccNode() {
    auto &status = statusGlobal[0];
    std::vector<uint32_t> distance(G.upperNodeIdBound()), eccLowerBound(G.upperNodeIdBound());

    auto maxDegreeNode = [&]() {
        node maxDegNode = 0;
        count maxDeg = 0;
        G.forNodes([&](const node u) {
            const auto degU = G.degree(u);
            if (degU > maxDeg) {
                maxDeg = degU;
                maxDegNode = u;
            }
        });
        return maxDegNode;
    };

    auto doBFS = [&](const node source) {
        std::queue<node> q;
        q.push(source);
        status[source] = NodeStatus::VISITED;
        distance[source] = 0;
        node farthest = 0;

        do {
            const auto u = q.front();
            q.pop();
            eccLowerBound[u] = std::max(eccLowerBound[u], distance[u]);
            farthest = u;

            G.forNeighborsOf(u, [&](const node v) {
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

    node source = maxDegreeNode();

    for (int i = 0; i < sweeps; ++i) {
        source = doBFS(source);
    }

    // Return node with minimum ecc lower bound
    return static_cast<node>(std::min_element(eccLowerBound.begin(), eccLowerBound.end())
                             - eccLowerBound.begin());
}

void ApproxEffectiveResistance::computeNodeSequence() {
    // We use thread 0's vector
    auto &status = statusGlobal[0];

    // Compute the biconnected components
    auto &bcc_ = *(bcc.get());
    bcc_.run();
    auto components = bcc_.getComponents();

    std::queue<node> queue;
    std::vector<node> curSequence;

    for (const auto &curComponent : components) {
        // Biconnected components with 3 nodes or less can be handled trivially.
        if (curComponent.size() <= 3) {
            sequences.push_back(std::move(curComponent));
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
            G.forNeighborsOf(u, [&](const node v) {
                if (status[v] == NodeStatus::VISITED) {
                    status[v] = NodeStatus::NOT_VISITED;
                    queue.push(v);
                }
            });
        } while (!queue.empty());

        sequences.push_back(std::move(curSequence));
    }

    // Set the root to the highest degree node within the largest biconnected component
    if (rootStrategy == RootStrategy::MaxDegree) {
        root = std::max_element(sequences.begin(), sequences.end(),
                                [](const std::vector<node> &c1, const std::vector<node> &c2) {
                                    return c1.size() < c2.size();
                                })
                   ->front();
       INFO("Using root strategy MaxDegree");
    } else if (rootStrategy == RootStrategy::Random) {
        auto &generator = generators.front();
        do {
            root = generator.nextUInt(G.upperNodeIdBound());
        } while (!G.hasNode(root));
       INFO("Using root strategy random");
    } else {
       INFO("Using root strategy MinApxEcc");
        root = approxMinEccNode();
    }

    biAnchor.resize(bcc_.numberOfComponents(), none);
    biParent.resize(bcc_.numberOfComponents(), none);

    if (isDebug) {
        G.forNodes([&](const node u) { assert(status[u] == NodeStatus::NOT_VISITED); });
    }

    // Topological order of biconnected components: tree of biconnected components starting from the
    // root's biconencted component. If the root is in multiple biconencted components, we take one
    // of them arbitrarily select one of them.
    std::queue<std::pair<node, index>> q;
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
        G.forNeighborsOf(front.first, [&](const node v) {
            if (status[v] == NodeStatus::NOT_VISITED) {
                distance[v] = distance[front.first] + 1;
                rootEcc = std::max(rootEcc, distance[v]);
                const auto &vComps = bcc_.getComponentsOfNode(v);
                for (const auto vComponentIndex : vComps) {
                    // Check if a new biconnected components has been found.
                    if (vComponentIndex != front.second && biAnchor[vComponentIndex] == none
                        && rootComps.find(vComponentIndex) == rootComps.end()) {
                        // The anchor cannot be the root, because the anchor does not have a parent.
                        // We handle biAnchor = none cases later.
                        biAnchor[vComponentIndex] = (v == root) ? none : v;
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

    if (isDebug) {
        G.forNodes([&](const node u) { assert(status[u] == NodeStatus::VISITED); });
        assert(topOrder.size() == bcc_.numberOfComponents());
        assert(std::unordered_set<index>(topOrder.begin(), topOrder.end()).size()
               == bcc_.numberOfComponents());
    }
}

void ApproxEffectiveResistance::computeBFSTree() {
    // Using thread 0's vector
    auto &status = statusGlobal[0];
    std::fill(status.begin(), status.end(), NodeStatus::NOT_VISITED);

    std::queue<node> queue;
    queue.push(root);
    status[root] = NodeStatus::VISITED;

    do {
        const auto currentNode = queue.front();
        queue.pop();
        G.forNeighborsOf(currentNode, [&](const node v) {
            if (status[v] == NodeStatus::NOT_VISITED) {
                status[v] = NodeStatus::VISITED;
                queue.push(v);
                bfsParent[v] = currentNode;
            }
        });
    } while (!queue.empty());

    if (isDebug) {
        checkBFSTree();
    }
}

void ApproxEffectiveResistance::sampleUST() {
    // Getting thread-local vectors
    auto &status = statusGlobal[omp_get_thread_num()];
    auto &parent = parentGlobal[omp_get_thread_num()];
    std::fill(status.begin(), status.end(), NodeStatus::NOT_IN_COMPONENT);
    std::fill(parent.begin(), parent.end(), none);

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
            assert(parent[curAnchor] != none);
        };

        if (sequence.size() == 2) {
            // Happens when the current component is the root component in the topological
            // order. In this case, the root plays the anchor's role.
            if (curAnchor == none) {
                const auto v = (sequence.front() == root) ? sequence.back() : sequence.front();
                assert(sequence.front() == root || sequence.back() == root);
                assert(v != root);
                parent[v] = root;
            } else {
                const auto v = (sequence[0] == curAnchor) ? sequence[1] : sequence[0];
                assert(v != curAnchor);
                assert(v != root);
                parent[v] = curAnchor;
                if (parent[curAnchor] == none) {
                    updateParentOfAnchor();
                }
            }
            if (isDebug) {
                checkTwoNodesSequence(sequence);
            }
            continue;
        }

        if (sequence.size() == 3) {
            // 3-vertex clique: we pick two edges at random
            std::pair<node, node> nonAnchors({none, none});
            if (curAnchor == none) {
                assert(!componentIndex);
                curAnchor = root;
            }

            for (const auto v : sequence) {
                if (v != curAnchor) {
                    if (nonAnchors.first == none) {
                        nonAnchors.first = v;
                    } else {
                        nonAnchors.second = v;
                        break;
                    }
                }
            }

            switch (generator.nextUInt(3)) {
            case 0:
                parent[nonAnchors.first] = curAnchor;
                parent[nonAnchors.second] = nonAnchors.first;
                break;
            case 1:
                parent[nonAnchors.first] = curAnchor;
                parent[nonAnchors.second] = curAnchor;
                break;
            case 2:
                parent[nonAnchors.first] = nonAnchors.second;
                parent[nonAnchors.second] = curAnchor;
                break;
            default:
                throw std::runtime_error("Unexpected random integer");
                break;
            }

            if (parent[curAnchor] == none) {
                updateParentOfAnchor();
            }
            continue;
        }

        // Components with > 3 vertices.

        // We start building the spanning tree from the first node of the
        // sequence.
        // Root of the spanning tree
        status[sequence[0]] = NodeStatus::IN_SPANNING_TREE;
        const auto curAnchorParent = (curAnchor != none) ? parent[curAnchor] : none;

        // All the remaining nodes in the components need to be visited.
        for (auto it = sequence.begin() + 1; it < sequence.end(); ++it) {
            status[*it] = NodeStatus::NOT_VISITED;
        }

        // Iterate over the remaining nodes to create the spanning tree.
        for (auto it = sequence.begin() + 1; it < sequence.end(); ++it) {
            auto currentNode = *it;

            if (status[currentNode] == NodeStatus::IN_SPANNING_TREE) {
                // Node already added to the spanning tree
                continue;
            }

            // Start a new random walk from the current node.
            do {
                // Get a random neighbor within the component
                node randomNeighbor;
                do {
                    randomNeighbor =
                        G.getIthNeighbor(currentNode, generator.nextUInt(G.degree(currentNode)));
                } while (status[randomNeighbor] == NodeStatus::NOT_IN_COMPONENT);

                assert(randomNeighbor != none);
                parent[currentNode] = randomNeighbor;
                currentNode = randomNeighbor;

            } while (status[currentNode] != NodeStatus::IN_SPANNING_TREE);

            // Last node encountered in the random walk (it is in the spanning tree);
            const auto walkEnd = currentNode;
            assert(status[walkEnd] == NodeStatus::IN_SPANNING_TREE);
            // Add the random walk to the spanning tree; eventually, reverse the path if the
            // anchor/root is encountered
            for (currentNode = *it; currentNode != walkEnd; currentNode = parent[currentNode]) {

                status[currentNode] = NodeStatus::IN_SPANNING_TREE;
                if (currentNode == curAnchor || currentNode == root) {

                    // Anchor of current component in the walk, we have to reverse the
                    // parent pointers
                    auto next = parent[currentNode];
                    auto nextParent = none;
                    auto tmp = currentNode;
                    do {
                        status[next] = NodeStatus::IN_SPANNING_TREE;
                        nextParent = parent[next];
                        parent[next] = currentNode;
                        currentNode = next;
                        next = nextParent;
                    } while (next != none);

                    if (tmp == root)
                        parent[root] = none;
                    else if (parent[curAnchor] == none) {
                        // Not none if articulation node visited by the parent component before
                        updateParentOfAnchor();
                    } else
                        parent[curAnchor] = curAnchorParent;

                    break;
                }
            }
        }

        for (const auto u : sequence) {
            assert(status[u] == NodeStatus::IN_SPANNING_TREE);
            status[u] = NodeStatus::NOT_IN_COMPONENT;
        }
    }

    // For debugging
    if (isDebug) {
        checkUST();
    }
}

void ApproxEffectiveResistance::dfsUST() {
    // Parent pointers in the UST
    const auto &parent = parentGlobal[omp_get_thread_num()];

    auto &tVisit = tVisitGlobal[omp_get_thread_num()];
    auto &tFinish = tFinishGlobal[omp_get_thread_num()];

    auto &status = statusGlobal[omp_get_thread_num()];
    std::fill(status.begin(), status.end(), NodeStatus::NOT_VISITED);

    std::stack<std::pair<node, Graph::NeighborIterator>> stack;
    stack.push({root, G.neighborRange(root).begin()});
    status[root] = NodeStatus::VISITED;

    uint32_t timestamp = 0;
    do {
        const auto u = stack.top().first;
        auto &iter = stack.top().second;
        assert(status[u] == NodeStatus::VISITED);
        const auto end = G.neighborRange(u).end();

        for (; iter != end; ++iter) {
            if (parent[*iter] == u && status[*iter] != NodeStatus::VISITED) {
                status[*iter] = NodeStatus::VISITED;
                tVisit[*iter] = ++timestamp;
                stack.push({*iter, G.neighborRange(*iter).begin()});
                break;
            }
        }

        if (iter == end) {
            stack.pop();
            tFinish[u] = ++timestamp;
        }
    } while (!stack.empty());

    if (isDebug) {
        G.forNodes([&](const node u) { assert(status[u] == NodeStatus::VISITED); });
        assert(!tVisit[root]);
    }
}

void ApproxEffectiveResistance::aggregateUST() {
    dfsUST();
    if (isDebug)
        checkTimeStamps();

    auto &approxEffResistance = approxEffResistanceGlobal[omp_get_thread_num()];
    const auto &tVisit = tVisitGlobal[omp_get_thread_num()];
    const auto &tFinish = tFinishGlobal[omp_get_thread_num()];
    const auto &parent = parentGlobal[omp_get_thread_num()];

    // Doing aggregation
    G.forNodes([&](const node u) {
        node p = bfsParent[u];
        node c = u;
        auto goUp = [&]() {
            c = p;
            p = bfsParent[p];
        };
        while (p != none) {
            // Edge in BSTree: e1 -> e2
            node e1 = p, e2 = c;
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
    if (!didInit)
        init();
    std::vector<count> ustSamplingTime(omp_get_max_threads());
    std::vector<count> ustAggregationTime(omp_get_max_threads());
    numberOfUSTs = static_cast<count>(std::ceil(rootEcc * computeNumberOfUSTs() / nProcessors));
#pragma omp parallel for schedule(dynamic)
    for (omp_index i = 0; i < static_cast<omp_index>(numberOfUSTs); ++i) {
        Aux::Timer timer;
        timer.start();
        sampleUST();
        timer.stop();
        ustSamplingTime[omp_get_thread_num()] += timer.elapsedNanoseconds();

        timer.start();
        // Update effective resistance values using sampled UST
        aggregateUST();
        timer.stop();
        ustAggregationTime[omp_get_thread_num()] += timer.elapsedNanoseconds();
    }

    for (omp_index i = 1; i < omp_get_max_threads(); ++i) {
        ustSamplingTime[0] += ustSamplingTime[i];
        ustAggregationTime[0] += ustAggregationTime[i];
    }

    time["UST sampling"] = static_cast<double>(ustSamplingTime[0]) / 1e9;
    time["UST aggregation"] = static_cast<double>(ustAggregationTime[0]) / 1e9;

    hasRun = true;

    computeDiagonal();
}

/*
 * Methods for sanity check.
 */
void ApproxEffectiveResistance::checkUST() const {
    std::vector<bool> visitedNodes(G.upperNodeIdBound());
    const auto &parent = parentGlobal[omp_get_thread_num()];
    // To debug
    G.forNodes([&](node u) {
        if (u == root) {
            assert(parent[u] == none);
        } else {
            assert(parent[u] != none);
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
    G.forNodes([&](node u) {
        if (u == root) {
            assert(bfsParent[u] == none);
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

void ApproxEffectiveResistance::checkTwoNodesSequence(const std::vector<node> &sequence) const {
    const std::vector<node> &parent = parentGlobal[omp_get_thread_num()];
    for (node u : sequence) {
        if (u == root) {
            assert(parent[u] == none);
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
    G.forNodes([&](const node u) {
        assert(tVisit[u] < tFinish[u]);
        if (u == root)
            assert(tVisit[u] == 0);
        else
            assert(tVisit[u] > 0);
        assert(tVisit[u] < 2 * G.numberOfNodes());
    });
}

} // namespace NetworKit
