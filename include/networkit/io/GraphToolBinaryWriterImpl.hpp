/*
 * GraphToolBinaryWriterImpl.hpp
 *
 *  Created on: 02.12.14
 *      Author: Maximilian Vogel
 */

#ifndef NETWORKIT_IO_GRAPH_TOOL_BINARY_WRITER_IMPL_HPP_
#define NETWORKIT_IO_GRAPH_TOOL_BINARY_WRITER_IMPL_HPP_

#include <cstdint>
#include <fstream>
#include <string>
#include <vector>

#include <networkit/auxiliary/Enforce.hpp>
#include <networkit/io/GraphToolBinaryWriter.hpp>

namespace NetworKit {

template <class GraphT>
GenericGraphToolBinaryWriter<GraphT>::GenericGraphToolBinaryWriter(bool littleEndianness)
    : littleEndianness(littleEndianness) {}

template <class GraphT>
void GenericGraphToolBinaryWriter<GraphT>::write(const GraphT &G, std::string_view path) {
    std::ofstream file(path.data(), std::ios::binary | std::ios::out);
    Aux::enforceOpened(file);
    writeHeader(file);
    writeComment(file);
    // write directed byte
    char directed = (char)G.isDirected();
    file.write(&directed, 1);
    // write number of nodes to the file
    uint64_t n = G.numberOfNodes();
    writeType<uint64_t>(file, 8, n);
    // write the adjacencies to the file
    writeAdjacencies(file, G);
    // close the file
    file.close();
}

template <class GraphT>
uint8_t GenericGraphToolBinaryWriter<GraphT>::getAdjacencyWidth(uint64_t n) {
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
void GenericGraphToolBinaryWriter<GraphT>::writeHeader(std::ofstream &file) {
    uint8_t header[8] = {0xe2, 0x9b, 0xbe, 0x20, 0x67, 0x74, 0x01, 0x00};
    header[7] |= (uint8_t)!this->littleEndianness;
    file.write((char *)header, 8);
}

template <class GraphT>
void GenericGraphToolBinaryWriter<GraphT>::writeComment(std::ofstream &file) {
    std::string s = "";
    uint64_t size = (uint64_t)s.size();
    writeType<uint64_t>(file, 8, size);
    if (size > 0) {
        file.write(s.c_str(), size);
    }
}

template <class GraphT>
void GenericGraphToolBinaryWriter<GraphT>::writeAdjacencies(std::ofstream &file, const GraphT &G) {
    using NodeT = typename GraphT::NodeT;

    // value of numNodes determines the size of the unsigned integer type storing the node ids
    int width = (int)getAdjacencyWidth(G.numberOfNodes());
    // graph-tool node ids are continuous, so compact them (in node order) if the graph has deleted
    // nodes
    std::vector<uint64_t> continuousIds;
    const auto upperNodeIdBound = static_cast<count>(G.upperNodeIdBound());
    if (G.numberOfNodes() != upperNodeIdBound) {
        continuousIds.resize(upperNodeIdBound, none);
        uint64_t nextId = 0;
        G.forNodes([&](NodeT u) -> void { continuousIds[u] = nextId++; });
    }
    const auto fileId = [&](NodeT v) -> uint64_t {
        return continuousIds.empty() ? static_cast<uint64_t>(v) : continuousIds[v];
    };

    if (!G.isDirected()) {
        // iterate over each node
        G.forNodes([&](NodeT u) -> void {
            std::vector<NodeT> adjacencies;
            // collect number of adjacencies for undirected case
            G.forNeighborsOf(u, [&](NodeT v) -> void {
                if (v <= u)
                    adjacencies.push_back(v);
            });
            // write number of adjacencies to the file
            count deg = adjacencies.size();
            writeType<uint64_t>(file, 8, deg);
            // write each adjacent to the file
            for (auto v : adjacencies) {
                writeType<uint64_t>(file, width, fileId(v));
            }
        });
    } else {
        // iterate over each node
        G.forNodes([&](NodeT u) -> void {
            // write number of adjacencies to the file
            count deg = G.degree(u);
            writeType<uint64_t>(file, 8, deg);
            // write each adjacent to the file
            G.forNeighborsOf(u,
                             [&](NodeT v) -> void { writeType<uint64_t>(file, width, fileId(v)); });
        });
    }
}

} // namespace NetworKit

#endif // NETWORKIT_IO_GRAPH_TOOL_BINARY_WRITER_IMPL_HPP_
