networkit_add_module(influencemaximization
		InfluenceMaximization.cpp
	)

networkit_module_link_modules(influencemaximization
        auxiliary
        graph
    )

add_subdirectory(test)

