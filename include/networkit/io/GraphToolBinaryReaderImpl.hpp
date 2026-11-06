/*
 * GraphToolBinaryReaderImpl.hpp
 *
 *  Created on: 02.12.14
 *      Author: Maximilian Vogel
 */

#ifndef NETWORKIT_IO_GRAPH_TOOL_BINARY_READER_IMPL_HPP_
#define NETWORKIT_IO_GRAPH_TOOL_BINARY_READER_IMPL_HPP_

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <networkit/auxiliary/Enforce.hpp>
#include <networkit/auxiliary/Log.hpp>
#include <networkit/io/GraphToolBinaryReader.hpp>

namespace NetworKit {

template <class GraphT>
GraphT GenericGraphToolBinaryReader<GraphT>::read(std::string_view path) {
    using NodeT = typename GraphT::NodeT;

    std::ifstream file(path.data(), std::ios::binary | std::ios::in);
    Aux::enforceOpened(file);
    // if the header is ok, continue reading the file
    if (checkHeader(file)) {
        // TODO: write function to ignore the comment by skipping it.
        readComment(file);
        bool directed = getDirected(file);
        count n = getNumNodes(file);
        // NullNodeId<NodeT> is the largest NodeT and marks "no node", so that at most max() nodes
        // fit.
        static_assert(NullNodeId<NodeT> == std::numeric_limits<NodeT>::max());
        constexpr auto maxNumberOfNodes = static_cast<count>(std::numeric_limits<NodeT>::max());
        if (n > maxNumberOfNodes) {
            throw std::runtime_error("graph-tool file " + std::string(path) + " has "
                                     + std::to_string(n)
                                     + " nodes, but the node type of the graph can hold at most "
                                     + std::to_string(maxNumberOfNodes) + " nodes");
        }
        GraphT G(n, false, directed);
        addOutNeighbours(file, n, G);
        file.close();
        G.shrinkToFit();
        return G;
    } else {
        throw std::runtime_error("File header is broken");
    }
}

template <class GraphT>
bool GenericGraphToolBinaryReader<GraphT>::checkHeader(std::ifstream &file) {
    uint8_t header[8];
    uint8_t reference[6] = {0xe2, 0x9b, 0xbe, 0x20, 0x67, 0x74};
    file.read((char *)header, 8);
    // compare the first 6 bytes of the header
    if (!std::equal(std::begin(reference), std::end(reference), std::begin(header))) {
        return false;
    }

    if (header[6] != 0x01) {
        // version byte not valid
        return false;
    }

    // get information about endianness
    if (header[7] == 0x00) {
        this->littleEndianness = true;
    } else if (header[7] == 0x01) {
        this->littleEndianness = false;
    } else {
        // endianness byte is wrong
        return false;
    }
    return true;
}

template <class GraphT>
std::string GenericGraphToolBinaryReader<GraphT>::readComment(std::ifstream &file) {
    // get length of comment
    auto length = readType<uint64_t>(file, 8);
    char *comment = new char[length];
    // read comment to char array
    file.read((char *)comment, length);
    std::string s(comment, length);
    delete[] comment;
    return s;
}

template <class GraphT>
bool GenericGraphToolBinaryReader<GraphT>::getDirected(std::ifstream &file) {
    char c;
    bool directed = false;
    file.read(&c, sizeof(char));
    if (c == 0x00) {
        directed = false;
    } else if (c == 0x01) {
        directed = true;
    } // optional else case: error handling
    return directed;
}

template <class GraphT>
uint64_t GenericGraphToolBinaryReader<GraphT>::getNumNodes(std::ifstream &file) {
    return readType<uint64_t>(file, 8);
}

template <class GraphT>
uint8_t GenericGraphToolBinaryReader<GraphT>::getAdjacencyWidth(uint64_t n) {
    if (n < (uint64_t)1 << 8) {
        return 1;
    } else if (n < (uint64_t)1 << 16) {
        return 2;
    } else if (n < (uint64_t)1 << 32) {
        return 4;
    } else {
        return 8;
    } // error handling?
}

template <class GraphT>
std::vector<std::vector<uint64_t>>
GenericGraphToolBinaryReader<GraphT>::getOutNeighbours(std::ifstream &file, uint64_t numNodes) {
    std::vector<std::vector<uint64_t>> adjacencies(numNodes);
    // value of numNodes determines the size of the unsigned integer type storing the node ids
    int width = (int)getAdjacencyWidth(numNodes);
    // iterate over each node
    for (uint64_t u = 0; u < numNodes; ++u) {
        // get number of adjacencies for the current node
        uint64_t numOutNeighbours = readType<uint64_t>(file, 8);
        // iterate over adjacencies of a noede
        for (uint64_t v = 0; v < numOutNeighbours; ++v) {
            // read current adjacency and add it to the adjacency list
            uint64_t node = readType<uint64_t>(file, width);
            adjacencies[u].push_back(node);
        }
    }
    return adjacencies;
}

template <class GraphT>
void GenericGraphToolBinaryReader<GraphT>::addOutNeighbours(std::ifstream &file, uint64_t numNodes,
                                                            GraphT &G) {
    using NodeT = typename GraphT::NodeT;
    using EdgeWeightT = typename GraphT::EdgeWeightT;

    // value of numNodes determines the size of the unsigned integer type storing the node ids
    int width = (int)getAdjacencyWidth(numNodes);
    // iterate over each node
    for (uint64_t u = 0; u < numNodes; ++u) {
        // get number of adjacencies for the current node
        uint64_t numOutNeighbours = readType<uint64_t>(file, 8);
        for (uint64_t v = 0; v < numOutNeighbours; ++v) {
            // read current adjacency and add edge
            uint64_t node = readType<uint64_t>(file, width);
            if (node >= numNodes) {
                throw std::runtime_error("graph-tool file has node id " + std::to_string(node)
                                         + ", but the valid ids are 0 to "
                                         + std::to_string(numNodes - 1));
            }
            // u and node are smaller than numNodes, which fits into NodeT (checked in read).
            if (!G.addEdge(static_cast<NodeT>(u), static_cast<NodeT>(node), EdgeWeightT{1}, true))
                WARN("Not adding edge ", u, "-", node, " since it is already present.");
        }
    }
}

} /* namespace NetworKit */

#endif // NETWORKIT_IO_GRAPH_TOOL_BINARY_READER_IMPL_HPP_
