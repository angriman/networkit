/*
 * GraphToolBinaryReader.hpp
 *
 *  Created on: 02.12.14
 *      Author: Maximilian Vogel
 */

#ifndef NETWORKIT_IO_GRAPH_TOOL_BINARY_READER_HPP_
#define NETWORKIT_IO_GRAPH_TOOL_BINARY_READER_HPP_

#include <fstream>
#include <string>

#include <networkit/io/GraphReader.hpp>

namespace NetworKit {

/**
 * Reads graphs from files in the binary format defined by graph-tool[1].
 * [1]: http://graph-tool.skewed.de/static/doc/gt_format.html
 *
 * Reading fails if the number of nodes does not fit into GraphT::NodeT.
 */
template <class GraphT>
class GenericGraphToolBinaryReader final : public GraphReaderBase<GraphT> {

public:
    GenericGraphToolBinaryReader() = default; // nullary constructor for Python shell

    /**
     * Given the path of an input file, read the graph contained.
     *
     * @param[in]  path  input file path
     */
    GraphT read(std::string_view path) override;

private:
    bool littleEndianness;

    void addOutNeighbours(std::ifstream &file, uint64_t numNodes, GraphT &G);

    std::vector<std::vector<uint64_t>> getOutNeighbours(std::ifstream &file, uint64_t numNodes);

    uint8_t getAdjacencyWidth(uint64_t n);

    uint64_t getNumNodes(std::ifstream &file);

    bool getDirected(std::ifstream &file);

    std::string readComment(std::ifstream &file);

    bool checkHeader(std::ifstream &file);

    template <typename Type>
    Type readType(std::ifstream &file, int width) {
        Type val = 0;
        uint8_t *bytes = new uint8_t[width];
        file.read((char *)bytes, width);
        if (this->littleEndianness) {
            for (int i = 0; i < width; ++i) {
                val |= ((Type)bytes[i] << (i * 8));
            }
        } else {
            for (int i = 0; i < width; ++i) {
                val |= ((Type)bytes[i] << (width - 1 - i) * 8);
            }
        }
        delete[] bytes;
        return val;
    }
};

using GraphToolBinaryReader = GenericGraphToolBinaryReader<Graph>;

} /* namespace NetworKit */

#include <networkit/io/GraphToolBinaryReaderImpl.hpp>

#endif // NETWORKIT_IO_GRAPH_TOOL_BINARY_READER_HPP_
