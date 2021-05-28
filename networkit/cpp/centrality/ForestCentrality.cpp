/*
 * ForestCentrality.cpp
 *
 *  Created on: 17.10.2019
 *      Author: Eugenio Angriman <angrimae@hu-berlin.de>
 */

// networkit-format

#include <networkit/centrality/ForestCentrality.hpp>

namespace NetworKit {

ForestCentrality::ForestCentrality(const Graph &G, node root, double inputEpsilon)
    : G(G), root(root), epsilon(inputEpsilon),
      decDegree(G.nodeRange().begin(), G.nodeRange().end()) {

    if (G.isDirected())
        throw std::runtime_error("Error: the input graph must be undirected.");

    if (G.numberOfNodes() < 2)
        throw std::runtime_error("Error: the graph should have at leasts two vertices");

    if (G.degree(root) != G.numberOfNodes() - 1)
        throw std::runtime_error("Root is not connected to all the other nodes.");

    if (G.numberOfNodes() != G.upperNodeIdBound())
        throw std::runtime_error(
            "Error: while we develop this algorithm, the graph must be compact");

    const count n = G.upperNodeIdBound();
    parentsGlobal.resize(omp_get_max_threads(), std::vector<node>(n));
    statusGlobal.resize(omp_get_max_threads(), std::vector<uint8_t>(n));
    std::sort(decDegree.begin(), decDegree.end(), [&G](const node n1, const node n2) {
        return G.weightedDegree(n1) < G.weightedDegree(n2);
    });

    generators.reserve(omp_get_max_threads());
    for (omp_index i = 0; i < omp_get_max_threads(); ++i)
        generators.emplace_back(Aux::Random::getURNG());

    if (G.isWeighted()) {
        ustsTotalWeight = 0;
        approxEffResistanceWeightedGlobal.resize(omp_get_max_threads(),
                                                 std::vector<high_prec_float>(n));
        weightsToParentGlobal.resize(omp_get_max_threads(), std::vector<edgeweight>(n));

        std::vector<std::vector<edgeweight>> inverseEdgeWeights(n);
        G.parallelForNodes([&](node u) {
            inverseEdgeWeights[u].reserve(G.degree(u));
            G.forNeighborsOf(u, [&inverseEdgeWeights, u](node v, edgeweight ew) {
                assert(ew > 0);
                inverseEdgeWeights[u].push_back(1. / ew);
            });
        });

        G.forNodes([&weightedDistr = weightedDistr, &inverseEdgeWeights](node u) {
            weightedDistr.emplace_back(inverseEdgeWeights[u].begin(), inverseEdgeWeights[u].end());
        });
    } else {
        approxEffResistanceGlobal.resize(omp_get_max_threads(), std::vector<count>(n));
        G.forNodes([&G = G, &uniformDistr = uniformDistr](node u) {
            const auto degU = G.degree(u);
            if (degU == 0)
                uniformDistr.emplace_back();
            else
                uniformDistr.emplace_back(0, degU - 1);
        });
    }
}

void ForestCentrality::run() {
    if (G.isWeighted())
        sampleUSTsWeighted();
    else
        sampleUSTsUnweighted();
    hasRun = true;
}

void ForestCentrality::sampleUSTsUnweighted() {
#pragma omp parallel for schedule(dynamic)
    for (omp_index i = 0; i < static_cast<omp_index>(numberOfUSTs); ++i) {
        // Getting thread-local vectors
        const omp_index thread = omp_get_thread_num();
        auto &parent = parentsGlobal[thread];
        auto &r = approxEffResistanceGlobal[thread];
        auto &status = statusGlobal[thread];
        auto &generator = generators[thread];
        uint8_t rootStatus = ++status[root];

        count nodesInUST = 1;

        for (const node walkStart : decDegree) {
            if (status[walkStart] == rootStatus)
                continue;

            node currentNode = walkStart;
            do {
                auto idx = uniformDistr[currentNode](generator.get());
                assert(idx < G.degree(currentNode));
                const node randomNeighbor = G.getIthNeighbor(currentNode, idx);
                parent[currentNode] = randomNeighbor;
                currentNode = randomNeighbor;
            } while (status[currentNode] != rootStatus);

            const node walkEnd = currentNode;
            node next = parent[walkStart];
            currentNode = walkStart;
            do {
                status[currentNode] = rootStatus;
                r[currentNode] += static_cast<count>(next == root);
                currentNode = next;
                next = parent[next];
                ++nodesInUST;
            } while (currentNode != walkEnd);

            if (nodesInUST == G.numberOfNodes())
                break;
        }
    }
}

void ForestCentrality::sampleUSTsWeighted() {
    std::cout << "running weighted algo with " << numberOfUSTs << " usts\n";
#pragma omp parallel for
    for (omp_index i = 0; i < static_cast<omp_index>(numberOfUSTs); ++i) {

        // Getting thread-local vectors
        const omp_index thread = omp_get_thread_num();
        auto &parent = parentsGlobal[thread];
        auto &weightToParent = weightsToParentGlobal[thread];
        auto &r = approxEffResistanceWeightedGlobal[thread];
        auto &status = statusGlobal[thread];
        auto &generator = generators[thread];
        uint8_t rootStatus = ++status[root];
        high_prec_float ustWeight = 1;

        count nodesInUST = 1;

        for (const node walkStart : decDegree) {
            if (status[walkStart] == rootStatus)
                continue;

            node currentNode = walkStart;
            do {
                const index neighborIndex = weightedDistr[currentNode](generator.get());
                const auto randomNeighbor = G.getIthNeighborWithWeight(currentNode, neighborIndex);
                parent[currentNode] = randomNeighbor.first;
                weightToParent[currentNode] = randomNeighbor.second;
                currentNode = randomNeighbor.first;
            } while (status[currentNode] != rootStatus);

            const node walkEnd = currentNode;
            node next = parent[walkStart];
            currentNode = walkStart;
            do {
                status[currentNode] = rootStatus;
                if (next != root) // Edges adjacent to the root have weight 1
                    ustWeight *= static_cast<high_prec_float>(weightToParent[currentNode]);
                currentNode = next;
                next = parent[next];
                ++nodesInUST;
            } while (currentNode != walkEnd);

            if (nodesInUST == G.numberOfNodes())
                break;
        }

        for (node u = 0; u < G.upperNodeIdBound() - 1; ++u) {
            assert(u != root);
            if (parent[u] == root)
                r[u] += ustWeight;
        }

#pragma omp critical
        ustsTotalWeight += ustWeight;
    }
}

} // namespace NetworKit
