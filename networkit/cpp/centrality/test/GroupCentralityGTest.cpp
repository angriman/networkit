#include <gtest/gtest.h>

#include <networkit/auxiliary/Random.hpp>
#include <networkit/centrality/ApproxGroupBetweenness.hpp>
#include <networkit/centrality/GedWalk.hpp>
#include <networkit/centrality/GroupCloseness.hpp>
#include <networkit/centrality/GroupClosenessGrowShrink.hpp>
#include <networkit/centrality/GroupClosenessLocalSearch.hpp>
#include <networkit/centrality/GroupClosenessLocalSwaps.hpp>
#include <networkit/centrality/GroupDegree.hpp>
#include <networkit/centrality/GroupForestCloseness.hpp>
#include <networkit/centrality/GroupHarmonicCloseness.hpp>
#include <networkit/components/ConnectedComponents.hpp>
#include <networkit/distance/BFS.hpp>
#include <networkit/generators/ErdosRenyiGenerator.hpp>
#include <networkit/generators/HyperbolicGenerator.hpp>
#include <networkit/graph/BFS.hpp>
#include <networkit/graph/Dijkstra.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <networkit/io/EdgeListReader.hpp>

namespace NetworKit {

class GroupCentralityGTest : public testing::TestWithParam<std::pair<bool, bool>> {
protected:
    bool isDirected() const noexcept;
    bool isWeighted() const noexcept;
};

INSTANTIATE_TEST_SUITE_P(InstantiationName, GroupCentralityGTest,
                         testing::Values(std::make_pair(false, false), std::make_pair(true, false),
                                         std::make_pair(false, true), std::make_pair(true, true)));

bool GroupCentralityGTest::isWeighted() const noexcept {
    return GetParam().first;
}

bool GroupCentralityGTest::isDirected() const noexcept {
    return GetParam().second;
}

TEST_P(GroupCentralityGTest, testGroupDegree) {
    Aux::Random::setSeed(42, false);
    constexpr count nodes = 12;
    constexpr count k = 5;
    auto g = ErdosRenyiGenerator(nodes, 0.3, isDirected()).generate();

    auto computeGroupDegree = [&](const std::vector<bool> &curGroup, const Graph &g) {
        count result = 0;
        g.forNodes([&](node u) {
            if (!curGroup[u]) {
                bool neighborInGroup = false;
                g.forInNeighborsOf(u, [&](node v) {
                    if (!neighborInGroup && curGroup[v]) {
                        neighborInGroup = true;
                        ++result;
                    }
                });
            }
        });

        return result;
    };

    GroupDegree gd(g, k, false);
    gd.run();
    auto scoreNoGroup = gd.getScore();

    GroupDegree gdIncludeGroup(g, k, true);
    gdIncludeGroup.run();
    auto scorePlusGroup = gdIncludeGroup.getScore();

    std::vector<bool> reference(nodes, false);
    for (count i = nodes - k; i < nodes; ++i) {
        reference[i] = true;
    }

    count maxScore = 0;

    do {
        count curScore = computeGroupDegree(reference, g);
        if (curScore > maxScore) {
            maxScore = curScore;
        }
    } while (std::next_permutation(reference.begin(), reference.end()));

    EXPECT_TRUE(scoreNoGroup > 0.5 * maxScore);
    EXPECT_TRUE(scorePlusGroup > (1.0 - 1.0 / std::exp(1.0)) * static_cast<double>(maxScore + k));
    EXPECT_EQ(scoreNoGroup, gd.scoreOfGroup(gd.groupMaxDegree()));
    EXPECT_EQ(scorePlusGroup, gdIncludeGroup.scoreOfGroup(gdIncludeGroup.groupMaxDegree()));
}

TEST_F(GroupCentralityGTest, testGroupBetweennessScore) {
    /** 0           5
     *  | \       / |
     *  1 - 3 - 4 - 6
     *  | /       \ |
     *  2           7
     */
    Graph G(8, false, false);
    G.addEdge(0, 1);
    G.addEdge(1, 2);
    G.addEdge(0, 3);
    G.addEdge(1, 3);
    G.addEdge(2, 3);

    G.addEdge(3, 4);

    G.addEdge(4, 5);
    G.addEdge(4, 6);
    G.addEdge(4, 7);
    G.addEdge(5, 6);
    G.addEdge(6, 7);

    BFS bfs(G, 0, true, false);
    // Naively computes the group betweenness of a group of nodes S
    auto naiveGB = [&](const std::vector<node> &S) -> double {
        double score = 0;
        count n = G.upperNodeIdBound();
        std::vector<bool> inGroup(n);
        for (node u : S)
            inGroup[u] = true;
        for (node source = 0; source < n; ++source) {
            bfs.setSource(source);
            bfs.run();
            for (node target = 0; target < n; ++target) {
                if (target == source)
                    continue;

                auto paths = bfs.getPaths(target);
                if (paths.empty())
                    continue;
                double curScore = 0;
                for (auto &path : paths) {
                    for (node u : path) {
                        if (u != source && u != target && inGroup[u]) {
                            curScore += 1;
                            break;
                        }
                    }
                }

                score += curScore / static_cast<double>(paths.size());
            }
        }

        return score;
    };

    ApproxGroupBetweenness gb(G, 2, 0.1);
    count z = G.numberOfNodes();
    for (count k = 1; k <= z; ++k) {
        std::vector<bool> inGroup(z);
        for (index i = 0; i < k; ++i)
            inGroup[z - i - 1] = true;
        do {
            std::vector<node> group;
            group.reserve(k);
            for (node u = 0; u < z; ++u) {
                if (inGroup[u]) {
                    group.push_back(u);
                    if (group.size() == k)
                        break;
                }
            }
            double computedScore = gb.scoreOfGroup(group);
            double naiveScore = naiveGB(group);
            EXPECT_NEAR(computedScore, naiveScore, 1e-6);
        } while (std::next_permutation(inGroup.begin(), inGroup.end()));
    }
}

TEST_F(GroupCentralityGTest, runTestApproxGroupBetweennessSmallGraph) {

    Aux::Random::setSeed(42, false);

    count n = 8;
    Graph g(n, false, false);

    g.addEdge(0, 2);
    g.addEdge(1, 2);
    g.addEdge(2, 3);
    g.addEdge(2, 4);
    g.addEdge(3, 5);
    g.addEdge(4, 5);
    g.addEdge(5, 6);
    g.addEdge(5, 7);
    g.addEdge(0, 5);

    count k = 2;
    double eps = 0.1;
    ApproxGroupBetweenness gb(g, k, eps);
    gb.run();

    double maxScore = 0;
    std::vector<bool> inGroup(n);
    for (index i = 0; i < k; ++i)
        inGroup[n - i - 1] = true;
    do {
        std::vector<node> group;
        group.reserve(k);
        for (node u = 0; u < n; ++u) {
            if (inGroup[u]) {
                group.push_back(u);
                if (group.size() == k)
                    break;
            }
        }

        maxScore = std::max(maxScore, gb.scoreOfGroup(group));
    } while (std::next_permutation(inGroup.begin(), inGroup.end()));

    EXPECT_TRUE(gb.scoreOfGroup(gb.groupMaxBetweenness()) >= maxScore * eps);
}

TEST_F(GroupCentralityGTest, testGroupCloseness) {
    Aux::Random::setSeed(42, false);

    Graph g(8, false, false);

    g.addEdge(0, 2);
    g.addEdge(1, 2);
    g.addEdge(2, 3);
    g.addEdge(2, 4);
    g.addEdge(3, 5);
    g.addEdge(4, 5);
    g.addEdge(5, 6);
    g.addEdge(5, 7);
    g.addEdge(0, 5);

    count k = 3;

    GroupCloseness gc(g, k);
    gc.run();
    auto apx = gc.groupMaxCloseness();
    EXPECT_NEAR(gc.scoreOfGroup(apx), 1.0, 1e-5);
}

TEST_P(GroupCentralityGTest, testGroupClosenessGrowShrink) {
    if (isDirected()) { // directed graphs are not supported
        Graph G(10, isWeighted(), true);
        std::array<node, 1> group;
        EXPECT_THROW(GroupClosenessGrowShrink(G, group.begin(), group.end()), std::runtime_error);
        return;
    }

    { // Empty input groups are not supported
        Graph G(10, isWeighted(), false);
        std::vector<node> emptyGroup;
        EXPECT_THROW(GroupClosenessGrowShrink(G, emptyGroup.begin(), emptyGroup.end()),
                     std::runtime_error);
    }

    const count k = 5;
    auto G = EdgeListReader{'\t', 0, "#", false, false}.read("input/MIT8.edgelist");
    G = ConnectedComponents::extractLargestConnectedComponent(G);

    if (isWeighted()) {
        G = GraphTools::toWeighted(G);
        G.forEdges(
            [&G](const node u, const node v) { G.setWeight(u, v, Aux::Random::probability()); });
    }

    auto farnessOfGroup = [&](const Graph &G, const std::unordered_set<node> &group) -> edgeweight {
        edgeweight farness = 0;
        if (G.isWeighted()) {
            Traversal::DijkstraFrom(
                G, group.begin(), group.end(),
                [&farness](node, const edgeweight distance) { farness += distance; });
        } else {
            Traversal::BFSfrom(
                G, group.begin(), group.end(),
                [&farness](node, const edgeweight distance) { farness += distance; });
        }

        return farness;
    };

    for (int seed : {1, 2, 3}) {
        Aux::Random::setSeed(seed, true);

        std::unordered_set<node> group;

        do {
            group.insert(GraphTools::randomNode(G));
        } while (group.size() < k);

        edgeweight sumDist = farnessOfGroup(G, group);
        GroupClosenessGrowShrink gc(G, group.begin(), group.end());
        gc.run();
        auto groupMaxCC = gc.groupMaxCloseness();
        count nSwaps = gc.numberOfIterations();

        EXPECT_EQ(groupMaxCC.size(), k);

        edgeweight sumDistGroupMaxCC =
            farnessOfGroup(G, std::unordered_set<node>(groupMaxCC.begin(), groupMaxCC.end()));
        if (nSwaps > 0) {
            EXPECT_GE(sumDist, sumDistGroupMaxCC);
        } else {
            EXPECT_EQ(sumDist, sumDistGroupMaxCC);
            std::for_each(groupMaxCC.begin(), groupMaxCC.end(),
                          [&group](const node u) { EXPECT_NE(group.find(u), group.end()); });
        }
    }
}

TEST_P(GroupCentralityGTest, testGroupClosenessLocalSwaps) {
    if (isDirected()) { // directed graphs are not supported
        Graph G(10, isWeighted(), true);
        std::array<node, 1> group;
        EXPECT_THROW(GroupClosenessLocalSwaps(G, group.begin(), group.end()), std::runtime_error);
        return;
    }

    { // Empty input groups are not supported
        Graph G(10, isWeighted(), false);
        std::vector<node> emptyGroup;
        EXPECT_THROW(GroupClosenessLocalSwaps(G, emptyGroup.begin(), emptyGroup.end()),
                     std::runtime_error);
    }

    const count k = 5;
    auto G = EdgeListReader{'\t', 0, "#", false, false}.read("input/MIT8.edgelist");
    G = ConnectedComponents::extractLargestConnectedComponent(G);

    if (isWeighted()) {
        G = GraphTools::toWeighted(G);
        G.forEdges(
            [&G](const node u, const node v) { G.setWeight(u, v, Aux::Random::probability()); });
    }

    auto farnessOfGroup = [&](const Graph &G, const std::unordered_set<node> &group) -> count {
        count farness = 0;
        Traversal::BFSfrom(G, group.begin(), group.end(),
                           [&farness](node, count distance) { farness += distance; });

        return farness;
    };

    for (int seed : {1, 2, 3}) {
        Aux::Random::setSeed(seed, true);

        std::unordered_set<node> group;

        do {
            group.insert(GraphTools::randomNode(G));
        } while (group.size() < k);

        const count maxSwaps = 100;
        const count sumDist = farnessOfGroup(G, group);
        GroupClosenessLocalSwaps gc(G, group.begin(), group.end(), maxSwaps);
        gc.run();
        auto groupMaxCC = gc.groupMaxCloseness();
        const count nSwaps = gc.numberOfSwaps();

        EXPECT_LE(nSwaps, maxSwaps);
        EXPECT_EQ(groupMaxCC.size(), k);

        count sumDistGroupMaxCC =
            farnessOfGroup(G, std::unordered_set<node>(groupMaxCC.begin(), groupMaxCC.end()));

        if (nSwaps > 0) {
            EXPECT_GE(sumDist, sumDistGroupMaxCC);
        } else {
            EXPECT_EQ(sumDist, sumDistGroupMaxCC);
            std::for_each(groupMaxCC.begin(), groupMaxCC.end(),
                          [&group](const node u) { EXPECT_NE(group.find(u), group.end()); });
        }
    }
}

TEST_P(GroupCentralityGTest, testGroupHarmonicCloseness) {

    const auto computeOpt = [&](const Graph &G, count k) -> double {
        std::vector<bool> inGroup(G.upperNodeIdBound());
        std::fill(inGroup.begin(), inGroup.begin() + k, true);
        double opt = -std::numeric_limits<double>::max();
        std::vector<node> group;
        group.reserve(k);
        do {
            group.clear();
            G.forNodes([&](node u) {
                if (inGroup[u])
                    group.push_back(u);
            });
            opt =
                std::max(opt, GroupHarmonicCloseness::scoreOfGroup(G, group.begin(), group.end()));
        } while (std::prev_permutation(inGroup.begin(), inGroup.end()));

        return opt;
    };

    const double weightUB = 10;
    for (const int seed : {1, 2, 3}) {
        Aux::Random::setSeed(seed, true);
        auto G = ErdosRenyiGenerator(20, 0.1, isDirected()).generate();
        double lambda = 1;

        if (isWeighted()) {
            double maxWeight = 0, minWeight = weightUB;
            G = GraphTools::toWeighted(G);
            G.forEdges([&](node u, node v) {
                const double curWeight = Aux::Random::real(0.001, weightUB);
                maxWeight = std::max(maxWeight, curWeight);
                minWeight = std::min(minWeight, curWeight);
                G.setWeight(u, v, curWeight);
            });
            lambda = minWeight / maxWeight;
        }

        // Guaranted approximation ratio
        double approxRatio;
        if (isDirected())
            approxRatio = lambda * (1. - 1. / (2. * std::exp(1.)));
        else
            approxRatio = lambda * (1. - 1. / std::exp(1.)) / 2.;

        for (const count k : {3, 4, 5}) {
            GroupHarmonicCloseness ghc(G, k);
            ghc.run();
            const auto group = ghc.groupMaxHarmonicCloseness();

            // Test group
            EXPECT_EQ(group.size(), k);
            EXPECT_EQ(std::unordered_set<node>(group.begin(), group.end()).size(), k);

            // Test quality
            const double score =
                GroupHarmonicCloseness::scoreOfGroup(G, group.begin(), group.end());
            const double opt = computeOpt(G, k);
            EXPECT_GE(opt, score);
            EXPECT_GE(score / opt, approxRatio);
        }
    }
}

TEST_P(GroupCentralityGTest, testGroupClosenessLocalSearch) {
    { // Empty groups are not allowed
        std::vector<node> emptyVector;
        Graph G(10, isWeighted(), isDirected());
        EXPECT_THROW(GroupClosenessLocalSearch(G, emptyVector.begin(), emptyVector.end()),
                     std::runtime_error);
    }

    const auto groupCloseness = [&](const Graph &G, const std::vector<node> &group) -> edgeweight {
        edgeweight groupFarness = 0;
        Traversal::DijkstraFrom(G, group.begin(), group.end(),
                                [&groupFarness](node, edgeweight dist) { groupFarness += dist; });
        return groupFarness > 0 ? 1. / groupFarness : 0;
    };

    Aux::Random::setSeed(1, true);
    auto G = ErdosRenyiGenerator(100, 0.1, isDirected()).generate();
    if (isWeighted()) {
        G = GraphTools::toWeighted(G);
        G.forEdges([&](node u, node v) { G.setWeight(u, v, Aux::Random::real(10)); });
    }

    std::unordered_set<node> initGroup;
    const count k = 5;
    do {
        initGroup.insert(GraphTools::randomNode(G));
    } while (initGroup.size() < k);

    const auto initGC = groupCloseness(G, std::vector<node>{initGroup.begin(), initGroup.end()});

    GroupClosenessLocalSearch gcls(G, initGroup.begin(), initGroup.end(), !G.isDirected());
    gcls.run();

    const auto group = gcls.groupMaxCloseness();
    EXPECT_EQ(std::unordered_set<node>(group.begin(), group.end()).size(), k);
    EXPECT_GE(groupCloseness(G, group), initGC);

    GroupClosenessLocalSearch gcls2(G, group, false);
    gcls2.run();

    EXPECT_EQ(gcls2.numberOfIterations(), 0);
}

TEST_P(GroupCentralityGTest, testGedWalk) {
    Aux::Random::setSeed(42, true);
    constexpr count k = 3;
    constexpr double epsilon = 0.01;
    auto g = ErdosRenyiGenerator(20, 0.1, isDirected()).generate();

    for (const auto bs : {GedWalk::BoundStrategy::geometric, GedWalk::BoundStrategy::spectral}) {
        for (const auto gs : {GedWalk::GreedyStrategy::lazy, GedWalk::GreedyStrategy::stochastic}) {
            GedWalk gedWalk(g, k, epsilon, -1.0, bs, gs);
            gedWalk.run();

            const auto apxScore = gedWalk.getApproximateScore();
            EXPECT_GE(apxScore, 0);
            const auto apxGroup = gedWalk.groupMaxGedWalk();
            EXPECT_EQ(std::unordered_set<node>(apxGroup.begin(), apxGroup.end()).size(), k);

            double maxScore = 0.;
            std::vector<node> group;
            std::vector<node> max_group;
            std::vector<bool> reference(g.numberOfNodes(), false);
            std::fill(reference.end() - k, reference.end(), true);

            do {
                group.clear();
                for (count i = 0; i < reference.size(); ++i) {
                    if (reference[i]) {
                        group.push_back(i);
                    }
                }

                const auto curScore = gedWalk.scoreOfGroup(group.begin(), group.end(), epsilon);
                if (curScore > maxScore) {
                    maxScore = curScore;
                    max_group = group;
                }
            } while (std::next_permutation(reference.begin(), reference.end()));

            EXPECT_GE(apxScore, (1. - 1. / std::exp(1)) * maxScore - epsilon);
        }
    }
}

TEST_F(GroupCentralityGTest, testGroupForestCloseness) {
    Aux::Random::setSeed(42, true);
    auto G = HyperbolicGenerator(200, 5).generate();

    const auto root = GraphTools::augmentGraph(G);
    const count k = 5;
    GroupForestCloseness gfc(G, root, k);
    gfc.run();
}

} // namespace NetworKit
