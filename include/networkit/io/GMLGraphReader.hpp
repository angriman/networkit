/*
 * GMLGraphReader.hpp
 *
 *  Created on: 18.09.2014
 *      Author: Maximilian Vogel (maximilian.vogel@student.kit.edu)
 */

#ifndef NETWORKIT_IO_GML_GRAPH_READER_HPP_
#define NETWORKIT_IO_GML_GRAPH_READER_HPP_

#include <networkit/io/GraphReader.hpp>

namespace NetworKit {

/**
 * @ingroup io
 *
 * Reading fails if the number of nodes does not fit into GraphT::NodeT.
 */
template <class GraphT>
class GenericGMLGraphReader final : public GraphReaderBase<GraphT> {
public:
    GenericGMLGraphReader() = default;

    /**
     * Takes a file path as parameter and returns a graph file.
     *
     * @param[in]  path  file path
     *
     * @param[out]  the graph read from file
     */
    GraphT read(std::string_view path) override;
};

using GMLGraphReader = GenericGMLGraphReader<Graph>;

} /* namespace NetworKit */

#include <networkit/io/GMLGraphReaderImpl.hpp>

#endif // NETWORKIT_IO_GML_GRAPH_READER_HPP_
