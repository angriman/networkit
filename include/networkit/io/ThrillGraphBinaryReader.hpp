/*
 * ThrillGraphBinaryReader.hpp
 *
 * @author Michael Hamann
 */

#ifndef NETWORKIT_IO_THRILL_GRAPH_BINARY_READER_HPP_
#define NETWORKIT_IO_THRILL_GRAPH_BINARY_READER_HPP_

#include <string>
#include <vector>

#include <networkit/io/GraphReader.hpp>

namespace NetworKit {

/**
 * Reads a graph format consisting of a serialized DIA of vector<uint32_t> from thrill.
 *
 * Reading fails if the number of nodes does not fit into GraphT::NodeT.
 */
template <class GraphT>
class GenericThrillGraphBinaryReader final : public GraphReaderBase<GraphT> {

public:
    /**
     * When the number of nodes is given, reading the graph is more efficient.
     *
     * @param[in]  n  The number of nodes
     */
    GenericThrillGraphBinaryReader(count n = 0);

    /**
     * Given the path of an input file, read the graph contained.
     *
     * @param[in]  path  input file path
     */
    GraphT read(std::string_view path) override;

    GraphT read(const std::vector<std::string> &paths);

private:
    const count n;
};

using ThrillGraphBinaryReader = GenericThrillGraphBinaryReader<Graph>;

} /* namespace NetworKit */

#include <networkit/io/ThrillGraphBinaryReaderImpl.hpp>

#endif // NETWORKIT_IO_THRILL_GRAPH_BINARY_READER_HPP_
