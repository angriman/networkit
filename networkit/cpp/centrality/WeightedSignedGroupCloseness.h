#ifndef WEIGHTEDSIGNEDGROUPCLOSENESS_H_
#define WEIGHTEDSIGNEDGROUPCLOSENESS_H_

#include "WeightedGroupCloseness.h"

namespace NetworKit {

class WeightedSignedGroupCloseness : public WeightedGroupCloseness {

public:
	WeightedSignedGroupCloseness(const Graph &G, const count k);
	void run() override;

protected:
	void updateMarginalGain();
	std::vector<int> spSign;
	std::vector<std::vector<int>> contributionSign;
};

} // namespace NetworKit

#endif
