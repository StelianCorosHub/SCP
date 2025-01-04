#pragma once

/*******************************************************************************************
*    BVHPlayer - A simple loader and container for .bvh motion capture files. Provides:
*
*     - A reader for the BVH file format
*     - A simple data structure to record joint / bone hierarchy
*     - A container for all the motion data stored in the BVH file
*     - Basic forward kinematics functions to access global positions for each body part
*
*******************************************************************************************/

#include <utils/mathUtils.h>
#include <utils/utils.h>

#define END_EFFECTOR_NAME "EndSite"

// Types of information "channels" that are possible in the BVH format
enum BVH_CHANNEL_TYPE {
    CHANNEL_X_POSITION = 0,
    CHANNEL_Y_POSITION = 1,
    CHANNEL_Z_POSITION = 2,
    CHANNEL_X_ROTATION = 3,
    CHANNEL_Y_ROTATION = 4,
    CHANNEL_Z_ROTATION = 5
};

// Types of information "channels" that are possible in the BVH format
enum SF_TYPE {
    SF_ROOT = 0,
    SF_JOINT = 1,
    SF_END_EFFECTOR = 2
};

// Data associated with skeletal features defined in the BVH format.
// These features can be a root, joints, or end effectors (called endSite).
// Bones are defined as the segments between a skeletal feature and its parent.
class BVHSkeletalFeature {
public:
    //store the index of the parent, if it's got any
    int parentIdx = -1;
    std::string name = "";
    // offset from parent to this skeletal feature, in local coordinates
    // (e.g. when all relative rotations are 0). It is used to define
    // bones of the skeletal structure
    V3D offset;
    // the configuration of each skeletal feature is defined through
    // channels of information. BVH motion data provides values for each
    // channel. These are the names
    Array<BVH_CHANNEL_TYPE> channels;

    bool isEndSite() const {
        return name == END_EFFECTOR_NAME;
    }

    bool isEndEffector() const {
        return isEndSite();
    }

    bool isRoot() const {
        return parentIdx == -1;
    }

    bool isJoint() const {
        return isEndSite() == false && isRoot() == false;
    }

    SF_TYPE getType() const {
        if (isRoot())
            return SF_ROOT;
        if (isEndSite())
            return SF_END_EFFECTOR;
        return SF_JOINT;
    }
};

typedef Array<double> BVHFrameData;
typedef Array<BVHFrameData> BVHMotionData;

//when t is 0, we'll return md1, when it's 1, we return md2.
//Between 0 and 1 we'll linearly interpolate
inline BVHFrameData interp(const Array<BVH_CHANNEL_TYPE>& channels, const BVHFrameData& md1, const BVHFrameData& md2, double t){
    clamp(t, 0, 1);
    BVHFrameData newFrame;

    for (uint i = 0; i < md1.size(); i++) {
        if (channels[i] == CHANNEL_X_ROTATION || channels[i] == CHANNEL_Y_ROTATION || channels[i] == CHANNEL_Z_ROTATION){
            //we have to worry about wrap-arounds and such...
            Quaternion q1(RAD(md1[i]), UP), q2(RAD(md2[i]), UP);
            Quaternion qInterp = q1.slerp(t, q2);
            newFrame.push_back(DEG(qInterp.getRotationAngle(UP)));
        } else
            newFrame.push_back((1-t) * md1[i] + t * md2[i]);
    }

    return newFrame;
}

// container that stores all positions of skeletal features, the relative orientations they
// have relative to their parent (i.e. joint angles), bone connectivity and the global
// orientation of every bone
class MocapSkeletonPose {
private:
    Array<P3D> sfPos;
    Array<SF_TYPE> sfType;
    Array<Quaternion> sfRelRot;         // store, in essence, joint angles / relative orientations
    Array<Quaternion> sfGlobalRot;      // store the global orientation of the bone who's got this SF as its parent

    //each bone is defined between a parent and child skeletal feature
    Array<int> bonePIdx;
    Array<int> boneCIdx;

public:
    MocapSkeletonPose(){

    }

    MocapSkeletonPose(int nSkeletalFeatures, int nBones) {
        sfPos.resize(nSkeletalFeatures, P3D());
        sfType.resize(nSkeletalFeatures, SF_JOINT);
        sfRelRot.resize(nSkeletalFeatures, Quaternion());
        sfGlobalRot.resize(nSkeletalFeatures, Quaternion());
        bonePIdx.resize(nBones, -1);
        boneCIdx.resize(nBones, -1);
    }

    int getSFCount() const {
        return sfPos.size();
    }

    int getBoneCount() const {
        return bonePIdx.size();
    }

    P3D getSFPos(int idx) const {
        return sfPos[idx];
    }

    void setSFPos(int idx, const P3D& p){
        sfPos[idx] = p;
    }

    void setSFRelRot(int idx, const Quaternion& qRel) {
        sfRelRot[idx] = qRel;
    }

    void setGlobalRotationForBoneWithSFIdxAsParent(int SFIdx, const Quaternion& q){
        sfGlobalRot[SFIdx] = q;
    }

    Quaternion getGlobalRotationForBoneWithSFIdxAsParent(int SFIdx) const {
        return sfGlobalRot[SFIdx];
    }

    void setGlobalBoneOrientation(int bIdx, const Quaternion& q) {
        setGlobalRotationForBoneWithSFIdxAsParent(bonePIdx[bIdx], q);
    }

    P3D getBoneStartPos(int idx) const {
        return sfPos[bonePIdx[idx]];
    }

    P3D getBoneEndPos(int idx) const {
        return sfPos[boneCIdx[idx]];
    }

    V3D getBoneVector(int idx) const {
        return V3D(getBoneStartPos(idx), getBoneEndPos(idx));
    }

    Quaternion getGlobalBoneRotation(int idx) const {
        return sfGlobalRot[bonePIdx[idx]];
    }

    Quaternion getBoneRotationRelativeToParent(int idx) const {
        return sfRelRot[bonePIdx[idx]];
    }

    SF_TYPE getSFType(int idx) const {
        return sfType[idx];
    }

    void setSFType(int idx, SF_TYPE t) {
        sfType[idx] = t;
    }

    void setBoneIndices(int bIdx, int pIdx, int cIdx){
        bonePIdx[bIdx] = pIdx;
        boneCIdx[bIdx] = cIdx;
    }

    double getBoneLength(int idx) {
        return V3D(getBoneStartPos(idx), getBoneEndPos(idx)).norm();
    }

    P3D getRootPosition() {
        //the first skeletal feature in a bvh skeleton needs to be the root...
        return sfPos[0];
    }

    Quaternion getRootOrientation() const {
        return sfRelRot[0];
    }

    //translates the entire skeleton pose by offset t
    void translateBy(const V3D& t){
        for (uint i = 0; i < sfPos.size(); i++){
            sfPos[i] += t;
        }
    }

    //translates the skeleton such that the root position is t
    void translateTo(const P3D& p){
        P3D rootP = getRootPosition();
        for (uint i = 0; i < sfPos.size(); i++){
            sfPos[i] -= rootP;
            sfPos[i] += p;
        }
    }

    void rotateInPlaceBy(const Quaternion& q){
        P3D rootP = getRootPosition();
        for (uint i = 0; i < sfPos.size(); i++)
            sfPos[i] = rootP + q * V3D(rootP, sfPos[i]);

        sfRelRot[0] *= q;
        for (int i = 0; i < getSFCount(); i++)
            sfGlobalRot[i] *= q;
    }

    void rotateBy(const Quaternion& q) {
        for (uint i = 0; i < sfPos.size(); i++)
            sfPos[i] = P3D() + q * V3D(P3D(0,0,0), sfPos[i]);
        sfRelRot[0] *= q;
        for (int i = 0; i < getSFCount(); i++)
            sfGlobalRot[i] *= q;
    }

    double getHeading(){
        return getRootOrientation().getHeadingAngle();
    }

    void setHeading(double h){
        rotateBy(Quaternion(h, UP) * Quaternion(-getHeading(), UP));
    }
};

class BVHSkeleton {
public:
    double scale = 0.1;

    // BVH skeletons are represented through location of joints. Note that end effectors are also stored as joints
    Array<BVHSkeletalFeature> skeletalFeatures;

    void add(const BVHSkeletalFeature& sf) {
        skeletalFeatures.push_back(sf);
    }

    //this method returns true if SF indexed by SFIdx has at least one child, otherwise it returns false
    bool hasChildren(int SFIdx) {
        //go through all the SFs, and see if any list this one as a parent
        for (auto sf : skeletalFeatures)
            if (sf.parentIdx == SFIdx)
                return true;
        return false;
    }

    void finalize(){
        //we want to make sure that the only skeletal features with no children are end effectors...
        //otherwise we might be missing some important bones / body part rotations
        for (uint i = 0; i < skeletalFeatures.size(); i++)
            if (hasChildren(i) == false && skeletalFeatures[i].isEndSite() == false){
                BVHSkeletalFeature newSF;
                newSF.parentIdx = i;
                newSF.name = END_EFFECTOR_NAME;
                add(newSF);
            }
    }

    //returns the total number of information channels that describe the configuration of the skeleton
    int getChannelCount() {
        int nCount = 0;
        for (auto sf : skeletalFeatures)
            nCount += sf.channels.size();
        return nCount;
    }

    Array<BVH_CHANNEL_TYPE> getChannels(){
        Array<BVH_CHANNEL_TYPE> res;
        for (auto sf : skeletalFeatures)
            for (uint i = 0; i < sf.channels.size(); i++)
                res.push_back(sf.channels[i]);
        return res;
    }

    //returns the number of bones in the skeleton. Bones are defined between every skeletal feature and its parent.
    int getNumberOfBones() {
        return skeletalFeatures.size() - 1;
    }

    //returns the number of joints in the skeleton.
    int getNumberOfJoints() {
        int n = 0;
        for (auto sf : skeletalFeatures)
            if (sf.isJoint())
                n++;
        return n;
    }

    //returns the number of end effectors
    int getNumberOfEndEffectors() {
        int n = 0;
        for (auto sf : skeletalFeatures)
            if (sf.isEndSite())
                n++;
        return n;
    }

    MocapSkeletonPose getDefaultPose(){
        MocapSkeletonPose pose(skeletalFeatures.size(), skeletalFeatures.size() - 1);

        int sfIdx = 0;
        int boneIdx = 0;
        for (auto sf : skeletalFeatures) {
            pose.setSFPos(sfIdx, P3D() + sf.offset * scale);
            pose.setSFType(sfIdx, sf.getType());
            if (sf.parentIdx != -1) {
                pose.setSFPos(sfIdx, pose.getSFPos(sfIdx) + pose.getSFPos(sf.parentIdx));
                pose.setBoneIndices(boneIdx, sf.parentIdx, sfIdx);
                boneIdx++;
            }
            sfIdx++;
        }

        return pose;
    }

    MocapSkeletonPose getPose(const BVHFrameData& frameData){
        MocapSkeletonPose pose = getDefaultPose();

        //first, make sure we've got enough channels of information here...
        if (frameData.size() != getChannelCount())
            return pose;

        //do all the forward kinematics work here...
        int channelIdx = 0;
        int sfIdx = 0;
        for (auto sf : skeletalFeatures) {
            V3D offset = sf.offset * scale;
            Quaternion q;
            for (uint j = 0; j < sf.channels.size(); j++) {
                double val = frameData[channelIdx + j];
                if (sf.channels[j] == CHANNEL_X_POSITION) offset += V3D(1,0,0) * val * scale;
                if (sf.channels[j] == CHANNEL_Y_POSITION) offset += V3D(0,1,0) * val * scale;
                if (sf.channels[j] == CHANNEL_Z_POSITION) offset += V3D(0,0,1) * val * scale;

                if (sf.channels[j] == CHANNEL_X_ROTATION) q = q * Quaternion(RAD(val), V3D(1,0,0));
                if (sf.channels[j] == CHANNEL_Y_ROTATION) q = q * Quaternion(RAD(val), V3D(0,1,0));
                if (sf.channels[j] == CHANNEL_Z_ROTATION) q = q * Quaternion(RAD(val), V3D(0,0,1));
            }

            P3D pos = P3D() + offset;

            pose.setSFRelRot(sfIdx, q);

            if (sf.parentIdx >= 0) {
                pos = pose.getSFPos(sf.parentIdx) + pose.getGlobalRotationForBoneWithSFIdxAsParent(sf.parentIdx) * (offset);
                q = pose.getGlobalRotationForBoneWithSFIdxAsParent(sf.parentIdx) * q;
            }

            pose.setGlobalRotationForBoneWithSFIdxAsParent(sfIdx, q);
            pose.setSFPos(sfIdx, pos);

            channelIdx += sf.channels.size();
            sfIdx++;
        }

        return pose;
    }

    double getUnscaledHeight() {
        MocapSkeletonPose pose = getDefaultPose();

        P3D pmin(INFINITY, INFINITY, INFINITY);
        P3D pmax(-INFINITY, -INFINITY, -INFINITY);

        for (int i = 0; i < pose.getSFCount(); i++){
            P3D sfPos = pose.getSFPos(i);
            pmin.x = MIN(pmin.x, sfPos.x);
            pmin.y = MIN(pmin.y, sfPos.y);
            pmin.z = MIN(pmin.z, sfPos.z);

            pmax.x = MAX(pmax.x, sfPos.x);
            pmax.y = MAX(pmax.y, sfPos.y);
            pmax.z = MAX(pmax.z, sfPos.z);
        }

        return V3D(pmin, pmax).dot(UP) / scale;
    }
};

class BVHClip;

//----------------------------------------------------------------------------------
// BVH Reader
//----------------------------------------------------------------------------------

class BVHReader : public BasicTextParser {
private:
    // Parse the "joint offset" part of the BVH File
    bool parseOffset(BVHSkeletalFeature &sf) {
        if (!parseKeyword("OFFSET")) { return false; }
        if (!parseDouble(sf.offset.x())) { return false; }
        if (!parseDouble(sf.offset.y())) { return false; }
        if (!parseDouble(sf.offset.z())) { return false; }
        if (!parseNewline()) { return false; }
        return true;
    }

    // Parse a channel type and return it in "channel"
    bool parseChannel(BVH_CHANNEL_TYPE &channel) {
        parseWhitespace();

        if (parseIfKeyword("Xposition")) {
            channel = CHANNEL_X_POSITION;
            return true;
        }

        if (parseIfKeyword("Yposition")) {
            channel = CHANNEL_Y_POSITION;
            return true;
        }

        if (parseIfKeyword("Zposition")) {
            channel = CHANNEL_Z_POSITION;
            return true;
        }

        if (parseIfKeyword("Xrotation")) {
            channel = CHANNEL_X_ROTATION;
            return true;
        }

        if (parseIfKeyword("Yrotation")) {
            channel = CHANNEL_Y_ROTATION;
            return true;
        }

        if (parseIfKeyword("Zrotation")) {
            channel = CHANNEL_Z_ROTATION;
            return true;
        }

        parserError("expected channel type");
        return false;
    }

// Parse the "channels" part of the BVH file format
    bool parseChannels(BVHSkeletalFeature &sf) {
        if (!parseKeyword("CHANNELS"))
            return false;
        int n = 0;
        if (!parseInt(n))
            return false;

        for (int i = 0; i < n; i++) {
            BVH_CHANNEL_TYPE ch;
            if (!parseChannel(ch))
                return false;
            sf.channels.push_back(ch);
        }

        return parseNewline();
    }

    // Parse skeletal features defined by the BVH file format
    bool parseSkeletalFeatures(BVHSkeleton &bvhSkeleton, int parentIndex) {
        while (isOneOf("JEje")) {
            BVHSkeletalFeature sf;
            sf.parentIdx = parentIndex;

            if (isOneOf("Jj")) {
                if (!parseKeyword("JOINT")) { return false; }
                if (!parseName(sf.name)) { return false; }
                if (!parseNewline()) { return false; }
                if (!parseKeyword("{")) { return false; }
                if (!parseNewline()) { return false; }
                if (!parseOffset(sf)) { return false; }
                if (!parseChannels(sf)) { return false; }

                bvhSkeleton.add(sf);

                if (!parseSkeletalFeatures(bvhSkeleton,
                                           (int) bvhSkeleton.skeletalFeatures.size() - 1)) { return false; }
                if (!parseKeyword("}")) { return false; }
                if (!parseNewline()) { return false; }
            } else if (isOneOf("Ee")) {
                if (!parseKeyword("End Site")) { return false; }
                sf.name = END_EFFECTOR_NAME;
                if (!parseNewline()) { return false; }
                if (!parseKeyword("{")) { return false; }
                if (!parseNewline()) { return false; }
                if (!parseOffset(sf)) { return false; }
                if (!parseKeyword("}")) { return false; }
                if (!parseNewline()) { return false; }
                bvhSkeleton.add(sf);
            } else
                return false;
        }
        return true;
    }

    // Parse the frame count
    bool parseFrames(int& frameCount) {
        if (!parseKeyword("Frames:")) { return false; }
        if (!parseInt(frameCount)) { return false; }
        if (!parseNewline()) { return false; }
        return true;
    }

    // Parse the frame time
    bool parseFrameTime(double &frameTime) {
        if (!parseKeyword("Frame Time:")) { return false; }
        if (!parseDouble(frameTime)) { return false; }
        if (!parseNewline()) { return false; }
        if (frameTime == 0.0f)
            frameTime = 1.0f / 60.0f;
        return true;
    }

    // Parse the motion data part of the BVH file format
    bool parseMotionData(int nFrames, int nChannels, Array<BVHFrameData>& motionData) {
        motionData.resize(nFrames);

        for (int i = 0; i < nFrames; i++) {
            motionData[i].resize(nChannels);

            for (int j = 0; j < nChannels; j++) {
                if (!parseDouble(motionData[i][j]))
                    return false;
            }

            if (!parseNewline())
                return false;
        }

        return true;
    }

public:

    BVHReader(const char *filename, bool verbose = true) {
        this->verbose = verbose;
        loadBufferFromFile(filename);
    }

    ~BVHReader(){
        freeBuffer();
    }

    bool load(BVHClip* bvhClip);
};

//----------------------------------------------------------------------------------
// BVH File Data
//----------------------------------------------------------------------------------

// Data structure matching what is present in the BVH file format
class BVHClip {
public:
    // Motion data applies to a BVH skeleton
    BVHSkeleton skeleton;

    //this is the time between frames
    double frame_dt = 0;

    //store the motion frames here
    Array<BVHFrameData> motionData;
    Array<BVH_CHANNEL_TYPE> channels;

    void loadFromFile(const char* fileName, bool verbose = true) {
        BVHReader bvhReader(fileName, verbose);
        bvhReader.load(this);
    }

    int getFrameCount() {
        return motionData.size();
    }

    int getChannelCount() {
        return channels.size();
    }

    //returns the length of the clip, measured in seconds
    double getLength(){
        return getFrameCount() * frame_dt;
    }

    BVHFrameData getZeroFrame(){
        BVHFrameData zFrame;
        zFrame.resize(skeleton.getChannelCount());
        for (uint i = 0; i < zFrame.size(); i++)
            zFrame[i] = 0;
        return zFrame;
    }

    MocapSkeletonPose getSkeletonPoseAtFrame(int idx) {
        clamp(idx, 0, getFrameCount()-1);

        return skeleton.getPose(motionData[idx]);
    }

    double getTimeForFrameIndex(int fIdx){
        return fIdx * frame_dt;
    }

    int getFrameIndexForTime(double t) {
        int idx = (int) (t / frame_dt);
        clamp(idx, 0, getFrameCount()-1);
        return idx;
    }

    MocapSkeletonPose getSkeletonPoseAtTime(double t) {
        int idx = getFrameIndexForTime(t);

        if (idx >= getFrameCount() - 1)
            return skeleton.getPose(motionData[getFrameCount() - 1]);

        double tLeft = t - getTimeForFrameIndex(idx);
        tLeft = mapTo01Range(tLeft, 0, frame_dt);

        return skeleton.getPose(interp(channels, motionData[idx], motionData[idx+1], tLeft));
    }

    MocapSkeletonPose getDefaultSkeletonPose(){
        return skeleton.getDefaultPose();
    }

};




