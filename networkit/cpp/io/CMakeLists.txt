networkit_add_module(io
    BinaryEdgeListPartitionReader.cpp
    BinaryEdgeListPartitionWriter.cpp
    BinaryPartitionReader.cpp
    BinaryPartitionWriter.cpp
    CoverReader.cpp
    CoverWriter.cpp
    DGSReader.cpp
    DibapGraphReader.cpp
    DotPartitionWriter.cpp
    EdgeListCoverReader.cpp
    EdgeListPartitionReader.cpp
    GraphIO.cpp
    MemoryMappedFile.cpp
    METISParser.cpp
    MTXParser.cpp
    MatrixMarketReader.cpp
    PartitionReader.cpp
    RasterReader.cpp
    RBMatrixReader.cpp
    SNAPEdgeListPartitionReader.cpp
)

networkit_module_link_modules(io
    algebraic auxiliary dynamics graph structures)

add_subdirectory(test)

