/*
 * KONECTGraphReaderImpl.hpp
 *
 *  Created on: 11.05.2018
 *      Author: Roman Bange
 *
 * The KONECT format is described in detail in
 * http://konect.uni-koblenz.de/downloads/konect-handbook.pdf
 *
 */

#ifndef NETWORKIT_IO_KONECT_GRAPH_READER_IMPL_HPP_
#define NETWORKIT_IO_KONECT_GRAPH_READER_IMPL_HPP_

#include <charconv>
#include <cstdlib>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <system_error>
#include <type_traits>
#include <unordered_map>
#include <tlx/unused.hpp>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/io/KONECTGraphReader.hpp>
#include <networkit/io/MemoryMappedFile.hpp>

namespace NetworKit {

template <class GraphT>
GenericKONECTGraphReader<GraphT>::GenericKONECTGraphReader(
    bool remapNodes, GraphReaderOptions::MultipleEdgesHandling handlingmethod)
    : remapNodes(remapNodes), multipleEdgesHandlingMethod(handlingmethod) {}

template <class GraphT>
GraphT GenericKONECTGraphReader<GraphT>::read(std::string_view path) {
    using NodeT = typename GraphT::NodeT;
    using EdgeWeightT = typename GraphT::EdgeWeightT;
    using Options = GraphReaderOptions;

    std::string graphFormat = "";
    std::string graphType = "";
    count numberOfNodes = -1;
    count numberOfEdges = -1;

    bool directed = true;
    bool weighted = false;
    bool multiple = false;
    int valuesPerLine = 2;
    bool secondPropertyLine = false;
    std::unordered_map<count, NodeT> nodeIdMap;

    MemoryMappedFile mmfile(path);
    auto it = mmfile.cbegin();
    auto end = mmfile.cend();

    // Returns
    // 0 iff it does not point to a line ending
    // 1 iff it is a single-char ending ('\r' or '\n')
    // 2 iff it is two-char ending ("\r\n")
    auto detectLineEnding = [&]() -> size_t {
        if (*it == '\r') {
            if ((it + 1) != end && *(it + 1) == '\n')
                return 2;
            return 1;
        }
        if (*it == '\n')
            return 1;

        return 0;
    };

    // Returns true if a line ending was found and skipped
    auto skipLineEnding = [&](bool required) -> bool {
        auto length = detectLineEnding();
        if (length) {
            it += length;
            return true;
        }

        if (required)
            throw std::runtime_error("No break symbol after first property line");

        return false;
    };

    // The following functions are helpers for parsing.
    auto skipWhitespace = [&]() -> void {
        while (it != end && (*it == ' ' || *it == '\t'))
            ++it;
        if (it == end) { // proper error message if file ends unexpected
            throw std::runtime_error("Unexpected end of file. Whitespace at end of file found");
        }
    };

    // This function parses in whole words
    auto scanWord = [&]() -> std::string {
        std::string word = "";
        while (it != end && *it != ' ' && *it != '\t' && *it != '\n' && *it != '\r') {
            word += *it;
            ++it;
        }
        return word;
    };

    auto scanId = [&]() -> count {
        count value{};
        const auto [past, ec] = std::from_chars(it, end, value);
        if (ec != std::errc{}) {
            throw std::runtime_error("Scanning node failed. The file may be corrupt.");
        }
        it = past;
        return value;
    };

    // Parses the weight at it directly into EdgeWeightT, so that integral weights must be integers
    // representable by EdgeWeightT.
    auto scanWeight = [&]() -> EdgeWeightT {
        const char *tokenEnd = it;
        while (tokenEnd != end && *tokenEnd != ' ' && *tokenEnd != '\t' && *tokenEnd != '\n'
               && *tokenEnd != '\r') {
            ++tokenEnd;
        }
        EdgeWeightT value{};
        const char *past = it;
        if constexpr (std::is_integral_v<EdgeWeightT>) {
            const auto [ptr, ec] = std::from_chars(it, tokenEnd, value);
            if (ec == std::errc{}) {
                past = ptr;
            }
        } else {
            char *strtodPast;
            value = static_cast<EdgeWeightT>(std::strtod(it, &strtodPast));
            past = strtodPast;
        }
        if (past == it || past != tokenEnd) {
            throw std::runtime_error("Scanning weight " + std::string(it, tokenEnd)
                                     + " failed. It may not fit into the edge weight type.");
        }
        it = past;
        return value;
    };

    // parse graph properties
    if (*it == '%') {
        ++it;
        skipWhitespace();

        // Parse graph format directed / undirected
        graphFormat = scanWord();

        if (graphFormat == "sym" || graphFormat == "bip") {
            directed = false;
        } else if (graphFormat != "asym") {
            throw std::runtime_error("Graph format not tagged properly. Format is: " + graphFormat);
        }
        skipWhitespace();
        // Parse graph weighting
        graphType = scanWord();
        // NOTE: This parser disregards the edge weight classification of "interval scale" and
        // "ratio scale"
        if (graphType == "weighted" || graphType == "posweighted" || graphType == "signed") {
            weighted = true;
            valuesPerLine = 3;
        } else if (graphType == "positive") { // multiple edges
            if (multipleEdgesHandlingMethod == Options::SUM_WEIGHTS_UP) {
                weighted = true;
            }
            multiple = true; // weights will be added
        } else if (graphType == "multisigned" || graphType == "multiweighted"
                   || graphType == "multiposweighted") {
            weighted = true;
            multiple = true;
            valuesPerLine = 3;
        } else if (graphType == "dynamic") {
            throw std::runtime_error("Dynamic graphs are not supported yet");
        } else if (graphType != "unweighted") {
            throw std::runtime_error("Graph weight not tagged properly. Weight is: " + graphType);
        }
        DEBUG("First property line read in. Format: " + graphFormat + " / Type: " + graphType);
        if (graphFormat == "bip") {
            INFO("Your graph is flagged as a bipartite one. Keep in mind that"
                 "NetworKit cannot handle this kind of graph in its special cases."
                 "It is imported as a usual undirected graph and edges might be discarded.");
        }
        if (multiple) {
            DEBUG("Selected handling method for multiple edges is: "
                  + std::to_string(multipleEdgesHandlingMethod));
        }
    } else {
        throw std::runtime_error("No graph properties line found. This line is mandatory.");
    }

    skipWhitespace();
    skipLineEnding(true);
    skipWhitespace();
    // second optional property line
    if (*it == '%') {
        ++it;
        skipWhitespace();
        numberOfEdges = scanId();
        skipWhitespace();
        numberOfNodes = scanId();
        secondPropertyLine = true;
        while (it != end && !detectLineEnding())
            ++it;
        if (it >= end) { // proper error message if file ends unexpected
            throw std::runtime_error("Unexpected end of file");
        }

        skipLineEnding(true);
        DEBUG("Second property line read in. Edges: " + std::to_string(numberOfEdges)
              + " / Nodes: " + std::to_string(numberOfNodes));
        tlx::unused(numberOfEdges);

        // NullNodeId<NodeT> is the largest NodeT and marks "no node", so that at most max() nodes
        // fit.
        static_assert(NullNodeId<NodeT> == std::numeric_limits<NodeT>::max());
        constexpr auto maxNumberOfNodes = static_cast<count>(std::numeric_limits<NodeT>::max());
        if (numberOfNodes > maxNumberOfNodes) {
            throw std::runtime_error("KONECT file " + std::string(path) + " has "
                                     + std::to_string(numberOfNodes)
                                     + " nodes, but the node type of the graph can hold at most "
                                     + std::to_string(maxNumberOfNodes) + " nodes");
        }
    }

    GraphT graph((secondPropertyLine ? numberOfNodes : 0), weighted, directed);

    // Without remapping, the 1-based node id in becomes in - 1.
    const auto checkNodeId = [](count in) -> void {
        if (in == 0) {
            throw std::runtime_error("KONECT node ids start at 1, but the file has node id 0");
        }
        // NullNodeId<NodeT> is the largest NodeT, so that the largest valid node id is smaller.
        if (in - 1 >= static_cast<count>(NullNodeId<NodeT>)) {
            throw std::runtime_error("Node id " + std::to_string(in)
                                     + " does not fit into the node type of the graph");
        }
    };

    // Map nodes and increase graph size if no second property is defined
    std::function<NodeT(count)> mapNode = [&](count in) -> NodeT {
        if (secondPropertyLine) {
            if (in > numberOfNodes) { // if file is corrupted
                // secondPropertyLine is made useless
                ERROR("Given number of nodes by file does not match actual graph");
                secondPropertyLine = false;
                if (remapNodes) { // if remapNodes is selected true, the map has to be initalized
                                  // with the existing nodes
                    nodeIdMap.reserve(numberOfNodes);
                    graph.forNodes(
                        [&](NodeT u) -> void { nodeIdMap.insert({static_cast<count>(u) + 1, u}); });
                }
                return mapNode(in);
            } else {
                checkNodeId(in);
                return static_cast<NodeT>(in - 1); // minus firstNode
            }
        } else {
            if (remapNodes) {
                auto it = nodeIdMap.find(in);
                if (it != nodeIdMap.end())
                    return it->second;
                if (graph.upperNodeIdBound() == NullNodeId<NodeT>) {
                    throw std::runtime_error("The file has more nodes than the node type of the "
                                             "graph can hold");
                }
                auto result = nodeIdMap.insert({in, graph.addNode()});
                assert(result.second);
                return result.first->second;
            } else {
                checkNodeId(in);
                const auto upperNodeIdBound = static_cast<count>(graph.upperNodeIdBound());
                if (upperNodeIdBound < in)
                    graph.addNodes(in - upperNodeIdBound);
                return static_cast<NodeT>(in - 1);
            }
        }
    };

    // Helper function for handling edges
    auto handleEdge = [&](NodeT source, NodeT target, EdgeWeightT weight) -> void {
        if (!graph.hasEdge(source, target)) {
            if (!graph.addEdge(source, target, weight, true))
                WARN("Not adding edge ", source, "-", target, " since it is already present.");
        } else if (multiple) {
            switch (multipleEdgesHandlingMethod) {
            case Options::DISCARD_EDGES:
                break; // Do nothing
            case Options::SUM_WEIGHTS_UP:
                graph.increaseWeight(source, target, weight);
                break;
            case Options::KEEP_MINIMUM_WEIGHT:
                if (graph.weight(source, target) > weight) {
                    graph.setWeight(source, target, weight);
                }
                break;
            default:
                throw std::runtime_error("Invalid multipleEdgesHandlingMethod: "
                                         + std::to_string(multipleEdgesHandlingMethod));
            }
        } else {
            DEBUG("[" + std::to_string(source) + "->" + std::to_string(target)
                  + "] Multiple edges detected but declared as: " + graphFormat + " and "
                  + graphType);
        }
    };

    while (it != end) {
        skipWhitespace();
        if (skipLineEnding(false)) {
            // We ignore empty lines.
            continue;
        } else if (*it == '#' || *it == '%') {
            // Skip non-linebreak characters.
            while (it != end && !detectLineEnding())
                ++it;
        } else {
            // Normal case parsing
            auto sourceId = scanId();
            if (it >= end) { // proper error message if file ends unexpected
                throw std::runtime_error("Unexpected end of file");
            }
            if (!(*it == ' ' || *it == '\t')) {
                throw std::runtime_error("Target ID cannot be parsed, maybe the file is corrupt");
            }
            skipWhitespace();

            auto targetId = scanId();
            if (it >= end) { // proper error message if file ends unexpected
                throw std::runtime_error("Unexpected end of file");
            }
            skipWhitespace();
            // Map the source first, so that remapped ids do not depend on the evaluation order of
            // function arguments.
            const NodeT source = mapNode(sourceId);
            const NodeT target = mapNode(targetId);
            if (valuesPerLine > 2) {
                auto edgeWeight = scanWeight();
                handleEdge(source, target, edgeWeight);
            } else {
                handleEdge(source, target, EdgeWeightT{1});
            }
        }
        // break lines
        skipWhitespace();
        skipLineEnding(true);
    }

    graph.shrinkToFit();
    return graph;
}
} // namespace NetworKit

#endif // NETWORKIT_IO_KONECT_GRAPH_READER_IMPL_HPP_
