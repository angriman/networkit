/*
 * DotWriter.cpp
 *
 *  Created on: Jun 5, 2013
 *      Author: forigem
 */

#include <fstream>
#include <string>
#include <networkit/Globals.hpp>
#include <networkit/graph/Graph.hpp>

#include <networkit/io/DotGraphWriter.hpp>

namespace NetworKit {

void DotGraphWriter::write(const Graph &G, const std::string &path) {
    std::ofstream file{path};

    file << "graph {\n";
    G.forEdges([&](node u, node v) { file << u << " -- " << v << ";\n"; });
    file << "}\n";
}

} /* namespace NetworKit */
