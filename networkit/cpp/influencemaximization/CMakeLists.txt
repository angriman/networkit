networkit_add_module(influencemaximization
		IndependentCascade.cpp
		InfluenceMaximization.cpp
		LinearThreshold.cpp
	)

networkit_module_link_modules(influencemaximization
        auxiliary
        graph
    )

add_subdirectory(test)

