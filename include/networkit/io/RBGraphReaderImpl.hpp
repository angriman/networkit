/*
 *  RBGraphReaderImpl.hpp
 *
 *  Created on: 16.10.2024
 *      Author: bernlu
 */

#ifndef NETWORKIT_IO_RB_GRAPH_READER_IMPL_HPP_
#define NETWORKIT_IO_RB_GRAPH_READER_IMPL_HPP_

#include <cmath>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>

#include <networkit/io/RBGraphReader.hpp>
#include <networkit/io/RBMatrixReader.hpp>

namespace NetworKit {

template <class GraphT>
GraphT GenericRBGraphReader<GraphT>::read(std::string_view path) {
    using NodeT = typename GraphT::NodeT;
    using EdgeWeightT = typename GraphT::EdgeWeightT;

    std::ifstream in{std::string(path)};
    if (!in.is_open()) {
        throw std::runtime_error("could not open: " + std::string(path));
    }

    RBMatrixReader reader;

    reader.readToVectors(in);

    if (reader.nMatrixCols != reader.nMatrixRows)
        throw std::runtime_error(
            "File does not contain a square matrix - cannot parse this file into a graph!");

    const count n = reader.nMatrixCols;
    // The node ids 0, ..., n - 1 must be representable by NodeT and differ from NullNodeId<NodeT>,
    // which marks "no node". Since NullNodeId<NodeT> is the largest NodeT, this allows at most
    // max() nodes.
    static_assert(NullNodeId<NodeT> == std::numeric_limits<NodeT>::max());
    constexpr auto maxNumberOfNodes = static_cast<count>(std::numeric_limits<NodeT>::max());
    if (n > maxNumberOfNodes) {
        throw std::runtime_error("RB file " + std::string(path) + " has " + std::to_string(n)
                                 + " nodes, but the node type of the graph can hold at most "
                                 + std::to_string(maxNumberOfNodes) + " nodes");
    }

    // The values of the file are read as double.
    const auto toWeight = [&](double value) -> EdgeWeightT {
        if constexpr (std::is_integral_v<EdgeWeightT>) {
            // max() + 1 is a power of two and thus exact, unlike max() for 64-bit types.
            constexpr double lowest =
                static_cast<double>(std::numeric_limits<EdgeWeightT>::lowest());
            constexpr double maxPlusOne =
                2.0 * static_cast<double>(std::numeric_limits<EdgeWeightT>::max() / 2 + 1);
            if (std::trunc(value) != value || value < lowest || value >= maxPlusOne) {
                throw std::runtime_error("RB file " + std::string(path) + " has the value "
                                         + std::to_string(value)
                                         + ", which is not valid for the edge weight type");
            }
        }
        return static_cast<EdgeWeightT>(value);
    };

    GraphT graph(n, !reader.patternOnly, !reader.symmetric);

    if (reader.pointers.size() <= n) {
        throw std::runtime_error("RB file " + std::string(path) + " has "
                                 + std::to_string(reader.pointers.size())
                                 + " column pointers instead of " + std::to_string(n + 1));
    }

    // iterate values in csc matrix and add edges to graph
    for (index col = 0; col < n; ++col) {
        for (index idx = reader.pointers[col]; idx < reader.pointers[col + 1]; ++idx) {
            const index row = reader.rowindex[idx];
            if (row >= n) {
                throw std::runtime_error("RB file " + std::string(path) + " has row index "
                                         + std::to_string(row + 1)
                                         + ", but the valid ones are 1 to " + std::to_string(n));
            }
            // row and col are smaller than n, which fits into NodeT (checked above).
            if (reader.patternOnly)
                graph.addEdge(static_cast<NodeT>(row), static_cast<NodeT>(col));
            else
                graph.addEdge(static_cast<NodeT>(row), static_cast<NodeT>(col),
                              toWeight(reader.values[idx]));
        }
    }

    return graph;
}

} /* namespace NetworKit */

#endif // NETWORKIT_IO_RB_GRAPH_READER_IMPL_HPP_
