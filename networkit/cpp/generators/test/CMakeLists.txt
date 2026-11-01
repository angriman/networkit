networkit_add_test(generators BarabasiAlbertGeneratorGTest
    auxiliary community dynamics global graph io viz)

networkit_add_test(generators ChungLuGeneratorGTest
    auxiliary community dynamics global graph io viz)

networkit_add_test(generators ClusteredRandomGraphGeneratorGTest
    auxiliary community dynamics global graph io viz)

networkit_add_test(generators DegreeSequenceGeneratorGTest
    auxiliary community dynamics global graph io viz)

networkit_add_test(generators GeneratorsGTest
    auxiliary community dynamics global graph io viz)

networkit_add_test(generators HyperbolicGeneratorGTest
    auxiliary community dynamics global graph io viz)

networkit_add_test(generators LFRGeneratorGTest
    auxiliary community dynamics global graph io viz)

networkit_add_test(generators RmatGeneratorGTest
    auxiliary community dynamics global graph io viz)

networkit_add_test(generators ErdosRenyiEnumeratorGTest)

networkit_add_benchmark(generators GeneratorsBenchmark
        auxiliary graph)
