/*
 * DotGraphWriter.hpp
 *
 *  Created on: Jun 5, 2013
 *      Author: forigem
 */

#ifndef NETWORKIT_IO_DOT_GRAPH_WRITER_HPP_
#define NETWORKIT_IO_DOT_GRAPH_WRITER_HPP_

#include <networkit/graph/Graph.hpp>
#include <networkit/io/GraphWriter.hpp>

namespace NetworKit {

/**

 * @ingroup io
 *
 * This class turns a graph into a very basic GraphViz file as documented in the official manual
 [1].
 * If a more thorough support is desired, please contact the developers:
 https://github.com/networkit/networkit
 *
 * [1] https://graphviz.gitlab.io/_pages/pdf/dotguide.pdf
 */
template <class GraphT>
class GenericDotGraphWriter final : public GraphWriterBase<GraphT> {
public:
    /**
     * Write a graph as a GraphViz/file.
     *
     * @param[in]	graph	The graph object
     * @param[in]	path	The file path to be written to
     */
    void write(const GraphT &G, std::string_view path) override;
};

using DotGraphWriter = GenericDotGraphWriter<Graph>;

} /* namespace NetworKit */

#include <networkit/io/DotGraphWriterImpl.hpp>

#endif // NETWORKIT_IO_DOT_GRAPH_WRITER_HPP_
