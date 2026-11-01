networkit_add_module(graph
    EulerianPath.cpp
    GraphBuilder.cpp
    GraphTools.cpp
    KruskalMSF.cpp
    PrimMSF.cpp
    RandomMaximumSpanningForest.cpp
    SpanningForest.cpp
    UnionMaximumSpanningForest.cpp
    DHBGraph.cpp
    )

networkit_module_link_modules(graph
    auxiliary structures)

add_subdirectory(test)

