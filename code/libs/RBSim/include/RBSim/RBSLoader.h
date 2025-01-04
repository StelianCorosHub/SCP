#pragma once

#include <utils/logger.h>
#include <utils/utils.h>

#include <RBSim/RBUtils.h>
#include <RBSim/RBSRepo.h>

/**
 * Loader for rigid bodies and the joints that can be used to interconnect them.
 */
class RBSLoader : public BasicTextParser {
public:
    /**
     * Sometimes we only need physical entity of RBs (e.g. when we don't use
     * GUI). In that case, set loadVisuals=false.
     */
    bool loadGUIElements = true;

    int getRBSLineType();

    RB* createAndPopulateRB();
    RBJoint* createAndPopulateRBJoint(RBSRepo* rbsRepo);

public:
    RBSLoader(const char *fName, bool loadGUIElements = true);
    ~RBSLoader();

    //adds the rigid bodies / joints described in the input file
    // to the RB repo that is passed in as a parameter
    void populate(RBSRepo* rbsRepo);

private:

};
