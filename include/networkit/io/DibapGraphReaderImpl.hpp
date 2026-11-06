/*
 * DibapGraphReaderImpl.hpp
 *
 *  Created on: Jun 12, 2013
 *      Author: Henning
 */

#ifndef NETWORKIT_IO_DIBAP_GRAPH_READER_IMPL_HPP_
#define NETWORKIT_IO_DIBAP_GRAPH_READER_IMPL_HPP_

#ifndef NETWORKIT_WINDOWS

#include <cstdio>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
#include <netinet/in.h>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/io/DibapGraphReader.hpp>

namespace NetworKit {

namespace DibapGraphReaderDetails {

// codes in file headers to distinguish type
enum class DibapIOType : short {
    XX = ('X' << 8) | 'X',
    GI = ('G' << 8) | 'I',
    GF = ('G' << 8) | 'F',
    HI = ('H' << 8) | 'I',
    HF = ('H' << 8) | 'F',
    P2 = ('P' << 8) | '2',
    P4 = ('P' << 8) | '4',
    AA = ('A' << 8) | 'A',
    T2 = ('T' << 8) | '2',
    TE = ('T' << 8) | 'E'
};

} // namespace DibapGraphReaderDetails

template <class GraphT>
GraphT GenericDibapGraphReader<GraphT>::read(std::string_view path) {
    using NodeT = typename GraphT::NodeT;
    using EdgeWeightT = typename GraphT::EdgeWeightT;

    int n, i;
    short type;
    FILE *file = NULL; // TODO: Use std::ifstream; this function won't reliably close the handle

    int V = 0;
    int dvw = 0;
    std::vector<int> vw;
    std::vector<int> of;
    std::vector<int> to;
    int dew = 0;
    std::vector<int> ew;
    int dxy = 0;
    int numE2 = 0;

    // init, try to open file
    DEBUG("reading graph in DibaP format... ");
    file = fopen(std::string(path).c_str(), "r");
    if (file == NULL) {
        throw std::runtime_error("cannot open file ");
    }

    // detect file type
    n = fread(&type, sizeof(short), 1, file);
    if (n != 1 || ntohs(type) != static_cast<short>(DibapGraphReaderDetails::DibapIOType::GI)) {
        throw std::runtime_error("bad file structure ");
    }

    // read number of vertices
    n = fread(&V, sizeof(int), 1, file);
    V = ntohl(V);
    TRACE("(V=%d %d) ", V, ", ", n);
    if (n != 1 || V < 0) {
        throw std::runtime_error("bad file structure ");
    }

    // read vertex weight dimension
    n = fread(&dvw, sizeof(int), 1, file);
    dvw = ntohl(dvw);
    TRACE("(dvw=%d %d) ", dvw, ", ", n);

    if (n != 1) {
        throw std::runtime_error("bad file structure ");
    }

    if (dvw > 0) {
        // read vertex weights
        vw.resize(V * dvw);
        n = fread(vw.data(), sizeof(int), V * dvw, file);
        TRACE("(vw %d) ", n);
        if (n != V * dvw) {
            throw std::runtime_error("bad file structure ");
        }
        for (i = 0; i < n; i++)
            vw[i] = ntohl(vw[i]);
    }

    of.resize(V + 1);
    n = fread(&of[0], sizeof(int), V + 1, file);
    TRACE("(of %d) ", n);
    if (n != V + 1) {
        throw std::runtime_error("bad file structure ");
    }
    for (i = 0; i < n; i++)
        of[i] = ntohl(of[i]);

    numE2 = of[V];
    to.resize(numE2);
    n = fread(&to[0], sizeof(int), numE2, file);
    TRACE("(to %d) ", n);
    if (n != numE2) {
        throw std::runtime_error("bad file structure ");
    }
    for (i = 0; i < n; i++)
        to[i] = ntohl(to[i]);

    n = fread(&dew, sizeof(int), 1, file);
    dew = ntohl(dew);
    TRACE("(dew=%d %d) ", dew, ", ", n);
    if (n != 1) {
        throw std::runtime_error("bad file structure ");
    }

    if (dew > 0) {
        int numWeights = numE2 * dew;
        ew.resize(numWeights);
        n = fread(&ew[0], sizeof(int), numWeights, file);
        TRACE("(ew %d) ", n);
        if (n != numWeights) {
            throw std::runtime_error("bad file structure ");
        }
        for (i = 0; i < n; i++)
            ew[i] = ntohl(ew[i]);
    }

    n = fread(&dxy, sizeof(int), 1, file);
    dxy = ntohl(dxy);
    TRACE("(dxy=%d %d) ", dxy, ", ", n);
    if (n != 1) {
        throw std::runtime_error("bad file structure ");
    }

    if (dxy > 0) {
        int numCoords = V * dxy;
        std::vector<float> xy;
        xy.resize(numCoords);
        n = fread(xy.data(), sizeof(float), numCoords, file);
        TRACE("(xy %d) ", n);
        if (n != numCoords) {
            throw std::runtime_error("bad file structure ");
        }

        coordinates.reserve(static_cast<size_t>(V));
        for (count u = 0; u < static_cast<count>(V); ++u) {
            coordinates.emplace_back(dxy);
            coordinates.back().copyFrom(xy.data() + (u * dxy));
        }
    }

    fclose(file);

    // The node ids 0, ..., V - 1 must be representable by NodeT and differ from NullNodeId<NodeT>,
    // which marks "no node". Since NullNodeId<NodeT> is the largest NodeT, this allows at most
    // max() nodes.
    static_assert(NullNodeId<NodeT> == std::numeric_limits<NodeT>::max());
    constexpr auto maxNumberOfNodes = static_cast<count>(std::numeric_limits<NodeT>::max());
    if (static_cast<count>(V) > maxNumberOfNodes) {
        throw std::runtime_error("DibaP file " + std::string(path) + " has " + std::to_string(V)
                                 + " nodes, but the node type of the graph can hold at most "
                                 + std::to_string(maxNumberOfNodes) + " nodes");
    }
    for (int e = 0; e < numE2; ++e) {
        if (to[e] < 0 || to[e] >= V) {
            throw std::runtime_error("DibaP file " + std::string(path) + " has node id "
                                     + std::to_string(to[e]) + ", but the valid ids are 0 to "
                                     + std::to_string(V - 1));
        }
    }

    // fill graph: FIXME: so far without node weights, extend to weights
    GraphT graph(static_cast<count>(V));

    if (dew > 0) {
        for (index v = 0; v < (count)V; ++v) {
            for (index e = of[v]; e < (index)of[v + 1]; ++e) {
                if (v <= (index)to[e]) {
                    // The graph is unweighted, so that the weight is ignored.
                    if (!graph.addEdge(static_cast<NodeT>(v), static_cast<NodeT>(to[e]),
                                       static_cast<EdgeWeightT>(ew[e]), true))
                        WARN("Not adding edge ", v, "-", to[e], " since it is already present.");
                }
            }
        }
    } else {
        for (index v = 0; v < (count)V; ++v) {
            for (index e = of[v]; e < (index)of[v + 1]; ++e) {
                if (v <= (index)to[e]) {
                    if (!graph.addEdge(static_cast<NodeT>(v), static_cast<NodeT>(to[e]),
                                       EdgeWeightT{1}, true))
                        WARN("Not adding edge ", v, "-", to[e], " since it is already present.");
                }
            }
        }
    }

    graph.shrinkToFit();
    return graph;
}

} /* namespace NetworKit */

#endif // NETWORKIT_WINDOWS

#endif // NETWORKIT_IO_DIBAP_GRAPH_READER_IMPL_HPP_
