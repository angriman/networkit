/*
 * GMLGraphWriter.hpp
 *
 *  Created on: 30.01.2013
 *      Author: Christian Staudt
 */

#ifndef NETWORKIT_IO_GML_GRAPH_WRITER_HPP_
#define NETWORKIT_IO_GML_GRAPH_WRITER_HPP_

#include <networkit/io/GraphWriter.hpp>

namespace NetworKit {

/**
 * @ingroup io
 */
template <class GraphT>
class GenericGMLGraphWriter final : public GraphWriterBase<GraphT> {
public:
    /** Default constructor */
    GenericGMLGraphWriter() = default;

    /**
     * Write a graph @a G and its coordinates to a GML file.
     *
     * @param[in]  G     Graph of type NetworKit with 2D coordinates.
     * @param[in]  path  Path to file.
     */
    void write(const GraphT &G, std::string_view path) override;
};

using GMLGraphWriter = GenericGMLGraphWriter<Graph>;

} /* namespace NetworKit */

#include <networkit/io/GMLGraphWriterImpl.hpp>

#endif // NETWORKIT_IO_GML_GRAPH_WRITER_HPP_
