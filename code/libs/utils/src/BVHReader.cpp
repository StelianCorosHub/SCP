
#include <utils/BVHReader.h>

bool BVHReader::load(BVHClip* bvhClip){
    bool parseSuccess = true;
    offset = 0;
    row = col = 1;

    bvhClip->skeleton = BVHSkeleton();
    bvhClip->frame_dt = -1;
    bvhClip->motionData.clear();

    if (buffer == nullptr) parseSuccess = false;

    // Hierarchy Data

    if (parseSuccess && !parseKeyword("HIERARCHY")) parseSuccess = false;
    if (parseSuccess && !parseNewline()) parseSuccess = false;

    if (parseSuccess && !parseKeyword("ROOT")) parseSuccess = false;
    BVHSkeletalFeature root;

    if (parseSuccess && !parseName(root.name)) parseSuccess = false;
    if (parseSuccess && !parseNewline()) parseSuccess = false;
    if (parseSuccess && !parseKeyword("{")) parseSuccess = false;
    if (parseSuccess && !parseNewline()) parseSuccess = false;
    if (parseSuccess && !parseOffset(root)) parseSuccess = false;
    if (parseSuccess && !parseChannels(root)) parseSuccess = false;

    bvhClip->skeleton.add(root);

    if (parseSuccess && !parseSkeletalFeatures(bvhClip->skeleton, 0)) parseSuccess = false;
    if (parseSuccess && !parseKeyword("}")) parseSuccess = false;
    if (parseSuccess && !parseNewline()) parseSuccess = false;

    bvhClip->skeleton.finalize();

    bvhClip->channels = bvhClip->skeleton.getChannels();

    // Motion Data

    if (parseSuccess && !parseKeyword("MOTION")) parseSuccess = false;
    if (parseSuccess && !parseNewline()) parseSuccess = false;
    int nFrames = 0;
    if (parseSuccess && !parseFrames(nFrames)) parseSuccess = false;

    if (parseSuccess && !parseFrameTime(bvhClip->frame_dt)) parseSuccess = false;

    if (parseSuccess && !parseMotionData(nFrames, bvhClip->skeleton.getChannelCount(), bvhClip->motionData)) parseSuccess = false;

    if (!parseSuccess && verbose == true) {
        Logger::print("Error: Could not parse BVH file %s", fName.c_str());
    } else {
        if (verbose) {
            Logger::print("Parsed '%s' successfully\n", fName.c_str());
            Logger::print("Loaded %d frames for a skeleton that has %d degrees of freedom. Motion data is available at %2.2lfHz.\n", bvhClip->motionData.size(), bvhClip->skeleton.getChannelCount(), 1.0 / bvhClip->frame_dt);
        }
    }

    return parseSuccess;
}
