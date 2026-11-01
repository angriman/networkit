// networkit-format
/*
 * ParallelScheduleBenchmark.cpp
 *
 * Compares OpenMP schedules for iterating over the nodes of a graph, on workloads whose cost per
 * node ranges from constant to strongly dependent on the degree. Run it in a Release build:
 *
 *   networkit_tests --gtest_filter='ParallelScheduleBenchmark.*'
 */

#include <algorithm>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <omp.h>
#include <string>
#include <utility>
#include <vector>
#include <gtest/gtest.h>

#include <networkit/Globals.hpp>
#include <networkit/auxiliary/Random.hpp>
#include <networkit/auxiliary/Timer.hpp>
#include <networkit/generators/BarabasiAlbertGenerator.hpp>
#include <networkit/generators/ErdosRenyiGenerator.hpp>
#include <networkit/generators/RmatGenerator.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/graph/GraphTools.hpp>

namespace NetworKit {

namespace {

constexpr count repetitions = 9;
constexpr count nodesLog2 = 18;
constexpr count numNodes = count{1} << nodesLog2;

struct Schedule {
    std::string name;
    omp_sched_t kind;
    int chunk;           // 0 selects the default of the schedule.
    int chunksPerThread; // If positive, overrides chunk with a size that scales with the graph.
};

struct NamedGraph {
    std::string name;
    Graph graph;
};

struct Measurement {
    double milliseconds;
    double imbalance; // Largest amount of work done by a thread, divided by the mean.
};

const std::vector<Schedule> &schedules() {
    static const std::vector<Schedule> all{
        {"static", omp_sched_static, 0, 0},
        {"static:64", omp_sched_static, 64, 0},
        {"guided", omp_sched_guided, 0, 0},
        {"guided:64", omp_sched_guided, 64, 0},
        {"dynamic:1", omp_sched_dynamic, 1, 0},
        {"dynamic:16", omp_sched_dynamic, 16, 0},
        {"dynamic:64", omp_sched_dynamic, 64, 0},
        {"dynamic:256", omp_sched_dynamic, 256, 0},
        {"dynamic:1024", omp_sched_dynamic, 1024, 0},
        {"dynamic:n/(32p)", omp_sched_dynamic, 0, 32},
        {"dynamic:n/(128p)", omp_sched_dynamic, 0, 128},
    };
    return all;
}

// The node ids of the graphs below are correlated with the degrees in different ways, which is
// what schedules that hand out large chunks first are sensitive to.
std::vector<NamedGraph> benchmarkGraphs() {
    std::vector<NamedGraph> graphs;

    // The hubs of a Barabasi-Albert graph are the nodes with the smallest ids.
    auto hubsFirst = BarabasiAlbertGenerator(8, numNodes).generate();

    std::vector<node> permutation(numNodes);
    std::iota(permutation.begin(), permutation.end(), node{0});
    std::shuffle(permutation.begin(), permutation.end(), Aux::Random::getURNG());
    graphs.push_back(
        {"ba-shuffled", GraphTools::getRemappedGraph(
                            hubsFirst, numNodes, [&](node u) -> node { return permutation[u]; })});
    graphs.push_back(
        {"ba-hubs-last", GraphTools::getRemappedGraph(hubsFirst, numNodes, [](node u) -> node {
             return numNodes - 1 - u;
         })});
    graphs.push_back({"ba-hubs-first", std::move(hubsFirst)});

    graphs.push_back({"rmat", RmatGenerator(nodesLog2, 8, 0.57, 0.19, 0.19, 0.05).generate()});

    // With few nodes, a chunk of fixed size can cover a large fraction of the graph.
    graphs.push_back({"rmat-small", RmatGenerator(12, 8, 0.57, 0.19, 0.19, 0.05).generate()});

    graphs.push_back({"erdos-renyi", ErdosRenyiGenerator(numNodes, 16. / numNodes).generate()});

    return graphs;
}

// Runs workload(u, scratch) on all nodes, which returns the amount of work it did on u.
template <class Workload>
Measurement measure(const Graph &G, int threads, const Schedule &schedule, Workload &&workload) {
    const auto z = static_cast<omp_index>(G.upperNodeIdBound());

    omp_set_num_threads(threads);
    omp_set_schedule(schedule.kind,
                     schedule.chunksPerThread > 0
                         ? std::max(1, static_cast<int>(z / (schedule.chunksPerThread * threads)))
                         : schedule.chunk);

    std::vector<std::vector<uint8_t>> scratch(threads, std::vector<uint8_t>(z));
    std::vector<uint64_t> work(threads);
    std::vector<double> times;
    double imbalance = 0;

    for (count repetition = 0; repetition < repetitions; ++repetition) {
        Aux::Timer timer;
        timer.start();
#pragma omp parallel
        {
            const int thread = omp_get_thread_num();
            uint64_t threadWork = 0;
#pragma omp for schedule(runtime) nowait
            for (omp_index u = 0; u < z; ++u) {
                if (G.hasNode(u)) {
                    threadWork += workload(static_cast<node>(u), scratch[thread]);
                }
            }
            work[thread] = threadWork;
        }
        timer.stop();
        times.push_back(static_cast<double>(timer.elapsedMicroseconds()) / 1000.);

        const auto total = std::accumulate(work.begin(), work.end(), uint64_t{0});
        imbalance = static_cast<double>(*std::max_element(work.begin(), work.end())) * threads
                    / static_cast<double>(total);
    }

    std::nth_element(times.begin(), times.begin() + times.size() / 2, times.end());
    return {times[times.size() / 2], imbalance};
}

template <class Workload>
void report(const NamedGraph &graph, const std::string &workloadName, Workload &&workload) {
    const auto sequential = measure(graph.graph, 1, schedules().front(), workload).milliseconds;
    for (const int threads : {2, 4}) {
        for (const auto &schedule : schedules()) {
            const auto [milliseconds, imbalance] =
                measure(graph.graph, threads, schedule, workload);
            std::cout << graph.name << ',' << workloadName << ',' << threads << ',' << schedule.name
                      << ',' << std::fixed << std::setprecision(3) << milliseconds << ','
                      << sequential / milliseconds << ',' << imbalance << '\n';
        }
    }
}

} // namespace

class ParallelScheduleBenchmark : public testing::Test {
protected:
    void SetUp() override {
        omp_get_schedule(&previousKind, &previousChunk);
        previousThreads = omp_get_max_threads();
    }

    void TearDown() override {
        omp_set_schedule(previousKind, previousChunk);
        omp_set_num_threads(previousThreads);
    }

private:
    omp_sched_t previousKind{};
    int previousChunk{};
    int previousThreads{};
};

TEST_F(ParallelScheduleBenchmark, nodeIteration) {
    Aux::Random::setSeed(42, true);
    const auto graphs = benchmarkGraphs();

    std::cout << "graph,workload,threads,schedule,milliseconds,speedup,imbalance\n";
    for (const auto &graph : graphs) {
        const auto &G = graph.graph;
        std::vector<double> input(G.upperNodeIdBound(), 1.);
        std::vector<double> output(G.upperNodeIdBound());

        // Constant work per node, where only the overhead of the schedule is visible.
        report(graph, "constant", [&](node u, std::vector<uint8_t> &) -> uint64_t {
            output[u] = static_cast<double>(G.degree(u));
            return 1;
        });

        // Work proportional to the degree, as in PageRank and in parallelForEdges.
        report(graph, "neighbors", [&](node u, std::vector<uint8_t> &) -> uint64_t {
            double sum = 0;
            G.forNeighborsOf(u, [&](node v) -> void { sum += input[v]; });
            output[u] = sum;
            return G.degree(u) + 1;
        });

        // Work proportional to the volume of the neighborhood, as in triangle counting.
        report(graph, "triangles", [&](node u, std::vector<uint8_t> &marked) -> uint64_t {
            uint64_t triangles = 0;
            uint64_t steps = G.degree(u) + 1;
            G.forNeighborsOf(u, [&](node v) -> void { marked[v] = 1; });
            G.forNeighborsOf(u, [&](node v) -> void {
                steps += G.degree(v);
                G.forNeighborsOf(v, [&](node w) -> void { triangles += marked[w]; });
            });
            G.forNeighborsOf(u, [&](node v) -> void { marked[v] = 0; });
            output[u] = static_cast<double>(triangles);
            return steps;
        });
    }
}

} /* namespace NetworKit */
