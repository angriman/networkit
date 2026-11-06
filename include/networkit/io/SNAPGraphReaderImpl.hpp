/*
 * SNAPGraphReaderImpl.hpp
 *
 *  Created on: 04.05.2018
 *      Author: Alexander van der Grinten
 */

#ifndef NETWORKIT_IO_SNAP_GRAPH_READER_IMPL_HPP_
#define NETWORKIT_IO_SNAP_GRAPH_READER_IMPL_HPP_

#include <charconv>
#include <stdexcept>
#include <string>
#include <system_error>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/io/MemoryMappedFile.hpp>
#include <networkit/io/SNAPGraphReader.hpp>

namespace NetworKit {

template <class GraphT>
GenericSNAPGraphReader<GraphT>::GenericSNAPGraphReader(bool directed, const bool &remapNodes,
                                                       const count &nodeCount)
    : directed(directed), nodeCount(nodeCount), remapNodes(remapNodes) {}

template <class GraphT>
GraphT GenericSNAPGraphReader<GraphT>::read(std::string_view path) {
    GraphT graph(0, false, directed);

    // The ids of a previously read file must not be reused.
    nodeIdMap.clear();

    // In the actual state this parameter has very little influence on the reader performance.
    // There can be a significant boost if it is possible to reserve space in the graph
    // initialization
    if (nodeCount != 0 && remapNodes)
        nodeIdMap.reserve(nodeCount);

    // Maps SNAP node IDs to consecutive NetworKit node IDs.
    auto mapNode = [&](count in) -> NodeT {
        if (remapNodes) {
            auto it = nodeIdMap.find(in);
            if (it != nodeIdMap.end())
                return it->second;
            if (graph.upperNodeIdBound() == NullNodeId<NodeT>) {
                throw std::runtime_error("The file has more nodes than the node type of the graph "
                                         "can hold");
            }
            auto result = nodeIdMap.insert({in, graph.addNode()});
            if (!result.second)
                throw std::runtime_error("Error in mapping nodes");
            return result.first->second;
        }
        // NullNodeId<NodeT> is the largest NodeT, so that the largest valid node id is smaller.
        if (in >= static_cast<count>(NullNodeId<NodeT>)) {
            throw std::runtime_error("Node id " + std::to_string(in)
                                     + " does not fit into the node type of the graph");
        }
        const auto upperNodeIdBound = static_cast<count>(graph.upperNodeIdBound());
        if (upperNodeIdBound <= in)
            graph.addNodes(in + 1 - upperNodeIdBound);
        return static_cast<NodeT>(in);
    };

    // This function modifies the graph on input.
    auto handleEdge = [&](NodeT source, NodeT target) -> void {
        if (!graph.hasEdge(source, target)) {
            graph.addEdge(source, target);
        } else {
            DEBUG("[" + std::to_string(source) + "->" + std::to_string(target)
                  + "] Multiple edges detected");
        }
    };

    MemoryMappedFile mmfile(path);
    auto it = mmfile.cbegin();
    auto end = mmfile.cend();

    // The following functions are helpers for parsing.
    auto skipWhitespace = [&]() -> void {
        while (it != end && (*it == ' ' || *it == '\t'))
            ++it;
    };

    auto scanId = [&]() -> count {
        count value{};
        const auto [past, ec] = std::from_chars(it, end, value);
        if (ec != std::errc{})
            throw std::runtime_error("Error in parsing file - looking for nodeId failed");
        it = past;
        return value;
    };

    // This loop does the actual parsing.
    while (it != end) {
        if (it >= end)
            throw std::runtime_error("Unexpected end of file");

        skipWhitespace();

        if (it == end)
            throw std::runtime_error("Unexpected end of file");

        if (*it == '\n') {
            // We ignore empty lines.
        } else if (*it == '#') {
            // Skip non-linebreak characters.
            while (it != end && *it != '\n')
                ++it;
        } else {
            auto sourceId = scanId();
            if (it == end)
                throw std::runtime_error("Unexpected end of file");
            if (!(*it == ' ' || *it == '\t'))
                throw std::runtime_error("Error in parsing file - pointer is whitespace");

            skipWhitespace();

            auto targetId = scanId();
            skipWhitespace();

            // Map the source first, so that remapped ids do not depend on the evaluation order of
            // function arguments.
            const NodeT source = mapNode(sourceId);
            handleEdge(source, mapNode(targetId));
        }

        if (it == end)
            throw std::runtime_error("Unexpected end of file");
        if (!(*it == '\n' || *it == '\r'))
            throw std::runtime_error("Line does not end with line break");
        ++it;
    }

    graph.shrinkToFit();
    return graph;
}

} // namespace NetworKit

#endif // NETWORKIT_IO_SNAP_GRAPH_READER_IMPL_HPP_
