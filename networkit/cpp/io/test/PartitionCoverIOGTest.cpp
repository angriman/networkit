// networkit-format
/*
 * PartitionCoverIOGTest.cpp
 *
 *  Created on: 12.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include <fstream>
#include <string>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <networkit/community/ClusteringGenerator.hpp>
#include <networkit/community/GraphClusteringTools.hpp>
#include <networkit/generators/ErdosRenyiGenerator.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/io/BinaryEdgeListPartitionReader.hpp>
#include <networkit/io/BinaryEdgeListPartitionWriter.hpp>
#include <networkit/io/BinaryPartitionReader.hpp>
#include <networkit/io/BinaryPartitionWriter.hpp>
#include <networkit/io/CoverReader.hpp>
#include <networkit/io/CoverWriter.hpp>
#include <networkit/io/EdgeListCoverReader.hpp>
#include <networkit/io/EdgeListPartitionReader.hpp>
#include <networkit/io/EdgeListReader.hpp>
#include <networkit/io/PartitionReader.hpp>
#include <networkit/io/PartitionWriter.hpp>
#include <networkit/structures/Partition.hpp>

namespace NetworKit {

class PartitionCoverIOGTest : public testing::Test {};

TEST_F(PartitionCoverIOGTest, testPartitionWriterAndReader) {
    // write clustering first
    std::string path = "output/example.clust";

    count n = 100;
    count k = 3;
    ErdosRenyiGenerator graphGen(n, 0.1);
    Graph G = graphGen.generate();

    ClusteringGenerator clusteringGen;
    Partition zeta = clusteringGen.makeRandomClustering(G, k);

    PartitionWriter writer;
    writer.write(zeta, path);

    // check if file exists
    bool exists = false;
    std::ifstream file(path);
    if (file) {
        exists = true;
    }
    EXPECT_TRUE(exists) << "clustering file should have been written to: " << path;

    PartitionReader reader;
    Partition read = reader.read(path);

    EXPECT_EQ(n, read.numberOfElements()) << "read clustering should contain n nodes";
    EXPECT_TRUE(GraphClusteringTools::isProperClustering(G, read))
        << "read clustering should be proper clustering of G";
    EXPECT_TRUE(GraphClusteringTools::equalClusterings(read, zeta, G))
        << "read clustering should be identical to created clustering";
}

TEST_F(PartitionCoverIOGTest, testEdgeListPartitionReader) {
    EdgeListPartitionReader reader(1);

    Partition zeta = reader.read("input/community.dat");
    // EXPECT_EQ(10, zeta.size());
    EXPECT_EQ(1u, zeta[0]);
    EXPECT_EQ(3u, zeta[1]);
    EXPECT_EQ(2u, zeta[2]);
    EXPECT_EQ(10u, zeta.numberOfElements());
}

TEST_F(PartitionCoverIOGTest, testEdgeListCoverReader) {
    EdgeListCoverReader reader(1);
    EdgeListReader gReader('\t', 1);

    Graph G = gReader.read("input/network_overlapping.dat");
    Cover zeta = reader.read("input/community_overlapping.dat", G);
    EXPECT_EQ(9u, zeta.upperBound());
    EXPECT_EQ(10u, zeta.numberOfElements());
    EXPECT_EQ(1u, zeta[0].count(1));
    EXPECT_EQ(3u, zeta[0].size());
    EXPECT_EQ(1u, zeta[3].size());
}

TEST_F(PartitionCoverIOGTest, testCoverReader) {
    CoverReader reader;
    EdgeListReader gReader('\t', 1);

    Graph G = gReader.read("input/network_overlapping.dat");
    Cover zeta = reader.read("input/community_overlapping.cover", G);
    EXPECT_EQ(9u, zeta.upperBound());
    EXPECT_EQ(10u, zeta.numberOfElements());
    EXPECT_EQ(1u, zeta[0].count(1));
    EXPECT_EQ(3u, zeta[0].size());
    EXPECT_EQ(1u, zeta[3].size());
}

TEST_F(PartitionCoverIOGTest, testCoverWriter) {
    std::string outpath = "output/coverWriter_test.cover";
    CoverWriter writer;
    CoverReader reader;
    EdgeListReader gReader('\t', 1);

    Graph G = gReader.read("input/network_overlapping.dat");
    Cover zeta = reader.read("input/community_overlapping.cover", G);

    writer.write(zeta, outpath);

    Cover read = reader.read(outpath, G);
    EXPECT_EQ(9u, read.upperBound());
    EXPECT_EQ(10u, read.numberOfElements());
    EXPECT_EQ(1u, read[0].count(1));
    EXPECT_EQ(3u, read[0].size());
    EXPECT_EQ(1u, read[3].size());
}

TEST_F(PartitionCoverIOGTest, testBinaryPartitionWriterAndReader) {
    Partition P(5);
    P.setUpperBound((1ull << 32));
    P[0] = 0;
    P[1] = 2007;
    P[2] = none;
    P[3] = (1ull << 31 | 1ull << 28);
    P[4] = 3925932491;

    std::string path = "output/partition.bin";

    BinaryPartitionWriter writer(8);
    BinaryPartitionReader reader(8);

    writer.write(P, path);
    Partition Q(reader.read(path));

    EXPECT_EQ(P[0], Q[0]);
    EXPECT_EQ(P[1], Q[1]);
    EXPECT_EQ(P[2], Q[2]);
    EXPECT_EQ(P[3], Q[3]);
    EXPECT_EQ(P[4], Q[4]);
    EXPECT_EQ(Q.upperBound(), P[4] + 1);
}

TEST_F(PartitionCoverIOGTest, testBinaryEdgeListPartitionWriterAndReader) {
    Partition P(5);
    P.setUpperBound((1ull << 32));
    P[0] = 0;
    P[1] = 2007;
    P[2] = none;
    P[3] = (1ull << 31 | 1ull << 28);
    P[4] = 3925932491;

    std::string path = "output/partition.bin";

    BinaryEdgeListPartitionWriter writer(1, 8);
    BinaryEdgeListPartitionReader reader(1, 8);

    writer.write(P, path);
    Partition Q(reader.read(path));

    EXPECT_EQ(P[0], Q[0]);
    EXPECT_EQ(P[1], Q[1]);
    EXPECT_EQ(P[2], Q[2]);
    EXPECT_EQ(P[3], Q[3]);
    EXPECT_EQ(P[4], Q[4]);
    EXPECT_EQ(Q.upperBound(), P[4] + 1);
}

} /* namespace NetworKit */
