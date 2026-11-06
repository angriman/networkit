/*
 * NetworkitBinaryWriterImpl.hpp
 *
 * @author Charmaine Ndolo <charmaine.ndolo@hu-berlin.de>
 */

#ifndef NETWORKIT_IO_NETWORKIT_BINARY_WRITER_IMPL_HPP_
#define NETWORKIT_IO_NETWORKIT_BINARY_WRITER_IMPL_HPP_

#include <cassert>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <vector>

#include <tlx/math/clz.hpp>

#include <networkit/auxiliary/Enforce.hpp>
#include <networkit/auxiliary/Log.hpp>
#include <networkit/io/NetworkitBinaryGraph.hpp>
#include <networkit/io/NetworkitBinaryWriter.hpp>

namespace NetworKit {

template <class GraphT>
GenericNetworkitBinaryWriter<GraphT>::GenericNetworkitBinaryWriter(
    uint64_t chunks, NetworkitBinaryWeights weightsType, NetworkitBinaryEdgeIDs edgeIndex)
    : chunks(chunks), weightsType(weightsType), edgeIndex(edgeIndex) {}

template <class GraphT>
void GenericNetworkitBinaryWriter<GraphT>::write(const GraphT &G, std::string_view path) {
    std::ofstream outfile(path.data(), std::ios::binary);
    Aux::enforceOpened(outfile);
    writeData(outfile, G);
    INFO("Written graph to ", path);
}

template <class GraphT>
std::vector<uint8_t> GenericNetworkitBinaryWriter<GraphT>::writeToBuffer(const GraphT &G) {
    std::stringstream outdata;
    writeData(outdata, G);
    std::string str = outdata.str();
    return std::vector<uint8_t>(str.begin(), str.end());
}

template <class GraphT>
template <class T>
void GenericNetworkitBinaryWriter<GraphT>::writeData(T &outStream, const GraphT &G) {
    using NodeT = typename GraphT::NodeT;
    using EdgeWeightT = typename GraphT::EdgeWeightT;

    nkbg::WEIGHT_FORMAT weightFormat;

    auto detectWeightsType = [&]() -> nkbg::WEIGHT_FORMAT {
        if (!G.isWeighted())
            return nkbg::WEIGHT_FORMAT::NONE;

        bool isUnsigned = true;
        bool fitsIntoInt64 = true;
        bool fitsIntoFloat = true;
        // The callback takes all four arguments, since shorter ones are ambiguous if EdgeWeightT
        // equals NodeT or edgeid.
        G.forEdges([&](NodeT, NodeT, EdgeWeightT w, edgeid) -> void {
            if constexpr (std::is_signed_v<EdgeWeightT>) {
                if (w < 0)
                    isUnsigned = false;
            }
            if constexpr (std::is_floating_point_v<EdgeWeightT>) {
                if (w != static_cast<int64_t>(w))
                    fitsIntoInt64 = false;
            }
            if (w != static_cast<float>(w))
                fitsIntoFloat = false;
        });
        if (fitsIntoInt64) {
            if (isUnsigned) {
                return nkbg::WEIGHT_FORMAT::VARINT;
            } else {
                return nkbg::WEIGHT_FORMAT::SIGNED_VARINT;
            }
        } else if (fitsIntoFloat) {
            return nkbg::WEIGHT_FORMAT::FLOAT;
        } else {
            return nkbg::WEIGHT_FORMAT::DOUBLE;
        }
    };

    switch (weightsType) {
    case NetworkitBinaryWeights::NONE:
        weightFormat = nkbg::WEIGHT_FORMAT::NONE;
        break;
    case NetworkitBinaryWeights::AUTO_DETECT:
        weightFormat = detectWeightsType();
        break;
    case NetworkitBinaryWeights::UNSIGNED_FORMAT:
        weightFormat = nkbg::WEIGHT_FORMAT::VARINT;
        break;
    case NetworkitBinaryWeights::SIGNED_FORMAT:
        weightFormat = nkbg::WEIGHT_FORMAT::SIGNED_VARINT;
        break;
    case NetworkitBinaryWeights::DOUBLE_FORMAT:
        weightFormat = nkbg::WEIGHT_FORMAT::DOUBLE;
        break;
    case NetworkitBinaryWeights::FLOAT_FORMAT:
        weightFormat = nkbg::WEIGHT_FORMAT::FLOAT;
        break;
    }

    switch (edgeIndex) {
    case NetworkitBinaryEdgeIDs::NO_EDGE_IDS:
        preserveEdgeIndex = false;
        break;
    case NetworkitBinaryEdgeIDs::AUTO_DETECT:
        preserveEdgeIndex = G.hasEdgeIds();
        break;
    case NetworkitBinaryEdgeIDs::WIRTE_EDGE_IDS:
        preserveEdgeIndex = G.hasEdgeIds();
        break;
    }

    nkbg::Header header;
    uint64_t adjListSize = 0;
    uint64_t adjTransposeSize = 0;

    auto setFeatures = [&]() -> void {
        header.features =
            (G.isDirected() & nkbg::DIR_MASK)
            | ((static_cast<uint64_t>(weightFormat) << nkbg::WGHT_SHIFT) & nkbg::WGHT_MASK)
            | ((preserveEdgeIndex << nkbg::INDEX_SHIFT) & nkbg::INDEX_MASK);
    };

    auto writeHeader = [&]() -> void {
        outStream.write(header.magic, 8);
        outStream.write(reinterpret_cast<char *>(&header.checksum), sizeof(uint64_t));
        outStream.write(reinterpret_cast<char *>(&header.features), sizeof(uint64_t));
        outStream.write(reinterpret_cast<char *>(&header.nodes), sizeof(uint64_t));
        outStream.write(reinterpret_cast<char *>(&header.chunks), sizeof(uint64_t));
        outStream.write(reinterpret_cast<char *>(&header.offsetBaseData), sizeof(uint64_t));
        outStream.write(reinterpret_cast<char *>(&header.offsetAdjLists), sizeof(uint64_t));
        outStream.write(reinterpret_cast<char *>(&header.offsetAdjTranspose), sizeof(uint64_t));
        outStream.write(reinterpret_cast<char *>(&header.offsetWeightLists), sizeof(uint64_t));
        outStream.write(reinterpret_cast<char *>(&header.offsetWeightTranspose), sizeof(uint64_t));
        outStream.write(reinterpret_cast<char *>(&header.offsetAdjIdLists), sizeof(uint64_t));
        outStream.write(reinterpret_cast<char *>(&header.offsetAdjIdTranspose), sizeof(uint64_t));
    };

    auto writeWeightsToFile = [&](EdgeWeightT w) -> void {
        uint8_t tmp[10];
        uint64_t weightSize;
        switch (weightFormat) {
        case nkbg::WEIGHT_FORMAT::VARINT:
            weightSize = nkbg::varIntEncode(static_cast<uint64_t>(w), tmp);
            outStream.write(reinterpret_cast<char *>(tmp), weightSize);
            break;
        case nkbg::WEIGHT_FORMAT::DOUBLE: {
            auto weight = static_cast<double>(w);
            outStream.write(reinterpret_cast<char *>(&weight), sizeof(double));
        } break;
        case nkbg::WEIGHT_FORMAT::SIGNED_VARINT: {
            uint64_t weight = nkbg::zigzagEncode(static_cast<int64_t>(w));
            weightSize = nkbg::varIntEncode(weight, tmp);
            outStream.write(reinterpret_cast<char *>(tmp), weightSize);
        } break;
        case nkbg::WEIGHT_FORMAT::FLOAT: {
            auto weight = static_cast<float>(w);
            outStream.write(reinterpret_cast<char *>(&weight), sizeof(float));
        } break;
        case nkbg::WEIGHT_FORMAT::NONE:
            break; // avoid compiler warning
        }
    };

    const auto nodes = static_cast<count>(G.upperNodeIdBound());
    if (nodes == 0) {
        strncpy(header.magic, FILE_FORMAT, 8);
        setFeatures();
        writeHeader();
        uint64_t zero = 0;
        outStream.write(reinterpret_cast<char *>(&zero), sizeof(zero));
        return;
    }

    if (nodes < chunks) {
        chunks = nodes;
        INFO("reducing chunks to ", chunks, " chunks");
    }

    // Compute first vertex of each chunk.
    std::vector<uint64_t> firstInChunk;
    firstInChunk.push_back(0);
    uint64_t firstNode = 0;
    for (uint64_t c = 1; c < chunks; c++) {
        firstNode += (nodes / chunks);
        firstInChunk.push_back(firstNode);
    }
    firstInChunk.push_back(nodes);

    // Compute encoded size of arrays and store in vector.
    uint64_t adjSize = 0;
    uint64_t transpSize = 0;
    uint64_t adjWeightSize = 0;
    uint64_t transpWeightSize = 0;
    uint64_t adjIndexSize = 0;
    uint64_t transpIndexSize = 0;
    std::unordered_map<count, uint64_t> nrOutNbrs;
    std::unordered_map<count, uint64_t> nrInNbrs;
    std::vector<size_t> adjOffsets;         // Prefix sum of size encoded adj arrays
    std::vector<size_t> transpOffsets;      // Prefix sum of encoded transposed adj arrays
    std::vector<size_t> adjWghtOffsets;     // Prefix sum of size encoded adj weights
    std::vector<size_t> transpWghtOffsets;  // Prefix sum of encoded transposed weights
    std::vector<size_t> adjIndexOffsets;    // Prefix sum of size encoded edge indexes
    std::vector<size_t> transpIndexOffsets; // Prefix sum of encoded transposed edge index

    auto computeWeightsOffsets = [&](EdgeWeightT w) -> uint64_t {
        uint64_t size = 0;
        uint8_t tmp[10];
        switch (weightFormat) {
        case nkbg::WEIGHT_FORMAT::VARINT:
            size = nkbg::varIntEncode(static_cast<uint64_t>(w), tmp);
            break;
        case nkbg::WEIGHT_FORMAT::DOUBLE:
            size = sizeof(double);
            break;
        case nkbg::WEIGHT_FORMAT::SIGNED_VARINT: {
            uint64_t weight = nkbg::zigzagEncode(static_cast<int64_t>(w));
            size += nkbg::varIntEncode(weight, tmp);
        } break;
        case nkbg::WEIGHT_FORMAT::FLOAT:
            size += sizeof(float);
            break;
        case nkbg::WEIGHT_FORMAT::NONE:
            break; // avoid compiler warning
        }

        return size;
    };

    auto computeIndexOffsets = [&](edgeid id) -> uint64_t {
        uint64_t size = 0;
        uint8_t tmp[10];
        size = nkbg::varIntEncode(id, tmp);
        return size;
    };

    for (uint64_t c = 0; c < chunks; c++) {
        for (uint64_t i = firstInChunk[c]; i < firstInChunk[c + 1]; i++) {
            uint64_t outNbrs = 0;
            uint64_t inNbrs = 0;
            uint8_t tmp[10];
            // i < nodes, so that it fits into NodeT.
            const auto n = static_cast<NodeT>(i);

            if (!G.hasNode(n)) {
                adjSize += nkbg::varIntEncode(0, tmp);
                transpSize += nkbg::varIntEncode(0, tmp);
                nrOutNbrs[i] = 0;
                nrInNbrs[i] = 0;
                continue;
            }
            if (!G.isDirected()) {
                G.forNeighborsOf(n, [&](NodeT, NodeT v, EdgeWeightT w, edgeid) -> void {
                    if (v <= n) {
                        outNbrs++;
                        adjSize += nkbg::varIntEncode(static_cast<uint64_t>(v), tmp);
                        adjWeightSize += computeWeightsOffsets(w);
                        if (preserveEdgeIndex)
                            adjIndexSize += computeIndexOffsets(G.edgeId(n, v));
                    }
                    if (v >= n) {
                        inNbrs++;
                        transpSize += nkbg::varIntEncode(static_cast<uint64_t>(v), tmp);
                        transpWeightSize += computeWeightsOffsets(w);
                        if (preserveEdgeIndex)
                            transpIndexSize += computeIndexOffsets(G.edgeId(n, v));
                    }
                });
            } else {
                G.forNeighborsOf(n, [&](NodeT, NodeT v, EdgeWeightT w, edgeid) -> void {
                    outNbrs++;
                    adjSize += nkbg::varIntEncode(static_cast<uint64_t>(v), tmp);
                    adjWeightSize += computeWeightsOffsets(w);
                    if (preserveEdgeIndex)
                        adjIndexSize += computeIndexOffsets(G.edgeId(n, v));
                });
                G.forInNeighborsOf(n, [&](NodeT, NodeT v, EdgeWeightT w, edgeid) -> void {
                    inNbrs++;
                    transpSize += nkbg::varIntEncode(static_cast<uint64_t>(v), tmp);
                    transpWeightSize += computeWeightsOffsets(w);
                    if (preserveEdgeIndex)
                        transpIndexSize += computeIndexOffsets(G.edgeId(v, n));
                });
            }
            adjListSize += outNbrs;
            adjSize += nkbg::varIntEncode(outNbrs, tmp);
            nrOutNbrs[i] = outNbrs;

            adjTransposeSize += inNbrs;
            transpSize += nkbg::varIntEncode(inNbrs, tmp);
            nrInNbrs[i] = inNbrs;
        }
        adjOffsets.push_back(adjSize);
        transpOffsets.push_back(transpSize);
        adjWghtOffsets.push_back(adjWeightSize);
        transpWghtOffsets.push_back(transpWeightSize);
        adjIndexOffsets.push_back(adjIndexSize);
        transpIndexOffsets.push_back(transpIndexSize);
    }
    // Write header.
    strncpy(header.magic, FILE_FORMAT, 8);
    header.checksum = 0;
    setFeatures();
    header.nodes = nodes;
    header.chunks = chunks;
    header.offsetBaseData = sizeof(nkbg::Header);
    auto nextOffset = header.offsetBaseData;
    nextOffset += nodes * sizeof(uint8_t)            // nodeFlags.
                  + (chunks - 1) * sizeof(uint64_t); // firstVertex.
    header.offsetAdjLists = nextOffset;
    nextOffset += (chunks - 1) * sizeof(uint64_t) // adjOffsets
                  + sizeof(uint64_t)              // adjListSize
                  + adjOffsets.back();            // Size of data
    header.offsetAdjTranspose = nextOffset;
    nextOffset += (chunks - 1) * sizeof(uint64_t) // transpOffsets
                  + sizeof(uint64_t)              // adjTranspSize
                  + transpOffsets.back();         // Size of data

    if (weightFormat != nkbg::WEIGHT_FORMAT::NONE) {
        header.offsetWeightLists = nextOffset;
        nextOffset += (chunks - 1) * sizeof(uint64_t) + adjWghtOffsets.back();
        header.offsetWeightTranspose = nextOffset;
        nextOffset += (chunks - 1) * sizeof(uint64_t) // transpWeightOffsets
                      + transpWghtOffsets.back();
    } else {
        header.offsetWeightLists = 0;
        header.offsetWeightTranspose = 0;
        nextOffset += (chunks - 1) * sizeof(uint64_t)    // weightOffsets
                      + (chunks - 1) * sizeof(uint64_t); // transpWeightOffsets
    }

    if (preserveEdgeIndex) {
        header.offsetAdjIdLists = nextOffset;
        nextOffset += (chunks - 1) * sizeof(uint64_t) + adjIndexOffsets.back();
        header.offsetAdjIdTranspose = nextOffset;
    } else {
        header.offsetAdjIdLists = 0;
        header.offsetAdjIdTranspose = 0;
    }
    writeHeader();
    // Write base data.
    for (count i = 0; i < nodes; ++i) {
        const auto u = static_cast<NodeT>(i);
        uint8_t nodeFlag = 0;
        if (!G.hasNode(u)) {
            nodeFlag |= nkbg::DELETED_BIT;
        }
        outStream.write(reinterpret_cast<char *>(&nodeFlag), sizeof(uint8_t));
    }

    assert(!firstInChunk[0]);
    for (uint64_t c = 1; c < chunks; c++) {
        outStream.write(reinterpret_cast<char *>(&firstInChunk[c]), sizeof(uint64_t));
    }

    // Write adjacency data.
    for (uint64_t c = 1; c < chunks; c++) {
        outStream.write(reinterpret_cast<char *>(&adjOffsets[c - 1]), sizeof(uint64_t));
    }
    // Write size of list
    outStream.write(reinterpret_cast<char *>(&adjListSize), sizeof(uint64_t));
    for (count i = 0; i < nodes; ++i) {
        const auto u = static_cast<NodeT>(i);
        uint8_t tmp[10];
        uint64_t nbrsSize = nkbg::varIntEncode(nrOutNbrs[i], tmp);
        outStream.write(reinterpret_cast<char *>(tmp), nbrsSize);
        if (!G.hasNode(u))
            continue; // deleted => no neighbors written
        G.forNeighborsOf(u, [&](NodeT v) -> void {
            uint64_t nodeSize;
            if (!G.isDirected()) {
                if (v <= u) {
                    nodeSize = nkbg::varIntEncode(static_cast<uint64_t>(v), tmp);
                    outStream.write(reinterpret_cast<char *>(tmp), nodeSize);
                }
            } else {
                nodeSize = nkbg::varIntEncode(static_cast<uint64_t>(v), tmp);
                outStream.write(reinterpret_cast<char *>(tmp), nodeSize);
            }
        });
    }

    // Write transpose data.
    for (uint64_t c = 1; c < chunks; c++) {
        outStream.write(reinterpret_cast<char *>(&transpOffsets[c - 1]), sizeof(uint64_t));
    }
    // Write size of transpose list.
    outStream.write(reinterpret_cast<char *>(&adjTransposeSize), sizeof(uint64_t));
    for (count i = 0; i < nodes; ++i) {
        const auto u = static_cast<NodeT>(i);
        uint8_t tmp[10];
        uint64_t nbrsSize = nkbg::varIntEncode(nrInNbrs[i], tmp);
        outStream.write(reinterpret_cast<char *>(tmp), nbrsSize);
        if (!G.hasNode(u))
            continue;
        uint64_t nodeSize;
        if (!G.isDirected()) {
            G.forNeighborsOf(u, [&](NodeT v) -> void {
                if (v >= u) {
                    nodeSize = nkbg::varIntEncode(static_cast<uint64_t>(v), tmp);
                    outStream.write(reinterpret_cast<char *>(tmp), nodeSize);
                }
            });
        } else {
            G.forInNeighborsOf(u, [&](NodeT v) -> void {
                nodeSize = nkbg::varIntEncode(static_cast<uint64_t>(v), tmp);
                outStream.write(reinterpret_cast<char *>(tmp), nodeSize);
            });
        }
    }

    // Write adj weights.
    for (uint64_t c = 1; c < chunks; c++) {
        outStream.write(reinterpret_cast<char *>(&adjWghtOffsets[c - 1]), sizeof(uint64_t));
    }
    for (count i = 0; i < nodes; ++i) {
        const auto u = static_cast<NodeT>(i);
        if (!G.hasNode(u))
            continue;
        G.forNeighborsOf(u, [&](NodeT, NodeT v, EdgeWeightT w, edgeid) -> void {
            if (!G.isDirected()) {
                if (v <= u) {
                    writeWeightsToFile(w);
                }
            } else {
                writeWeightsToFile(w);
            }
        });
    }

    // Write transpose weights.
    for (uint64_t c = 1; c < chunks; c++) {
        outStream.write(reinterpret_cast<char *>(&transpWghtOffsets[c - 1]), sizeof(uint64_t));
    }
    for (count i = 0; i < nodes; ++i) {
        const auto u = static_cast<NodeT>(i);
        if (!G.hasNode(u))
            continue;
        if (!G.isDirected()) {
            G.forNeighborsOf(u, [&](NodeT, NodeT v, EdgeWeightT w, edgeid) -> void {
                if (v >= u) {
                    writeWeightsToFile(w);
                }
            });
        } else {
            G.forInNeighborsOf(
                u, [&](NodeT, NodeT, EdgeWeightT w, edgeid) -> void { writeWeightsToFile(w); });
        }
    }

    // Write adj index.
    for (uint64_t c = 1; c < chunks; c++) {
        outStream.write(reinterpret_cast<char *>(&adjIndexOffsets[c - 1]), sizeof(uint64_t));
    }
    if (preserveEdgeIndex) {
        for (count i = 0; i < nodes; ++i) {
            const auto u = static_cast<NodeT>(i);
            if (!G.hasNode(u))
                continue;
            G.forNeighborsOf(u, [&](NodeT v) -> void {
                if (!G.isDirected()) {
                    if (v <= u) {
                        uint8_t tmp[10];
                        uint64_t edgeidSize = nkbg::varIntEncode(G.edgeId(u, v), tmp);
                        outStream.write(reinterpret_cast<char *>(tmp), edgeidSize);
                    }
                } else {
                    uint8_t tmp[10];
                    uint64_t edgeidSize = nkbg::varIntEncode(G.edgeId(u, v), tmp);
                    outStream.write(reinterpret_cast<char *>(tmp), edgeidSize);
                }
            });
        }
    }
    // Write transpose index.
    for (uint64_t c = 1; c < chunks; c++) {
        outStream.write(reinterpret_cast<char *>(&transpIndexOffsets[c - 1]), sizeof(uint64_t));
    }
    if (preserveEdgeIndex) {
        for (count i = 0; i < nodes; ++i) {
            const auto u = static_cast<NodeT>(i);
            if (!G.hasNode(u))
                continue;
            if (!G.isDirected()) {
                G.forNeighborsOf(u, [&](NodeT v) -> void {
                    if (v >= u) {
                        uint8_t tmp[10];
                        uint64_t edgeidSize = nkbg::varIntEncode(G.edgeId(u, v), tmp);
                        outStream.write(reinterpret_cast<char *>(tmp), edgeidSize);
                    }
                });
            } else {
                G.forInNeighborsOf(u, [&](NodeT v) -> void {
                    uint8_t tmp[10];
                    uint64_t edgeidSize = nkbg::varIntEncode(G.edgeId(v, u), tmp);
                    outStream.write(reinterpret_cast<char *>(tmp), edgeidSize);
                });
            }
        }
    }
}

} // namespace NetworKit

#endif // NETWORKIT_IO_NETWORKIT_BINARY_WRITER_IMPL_HPP_
