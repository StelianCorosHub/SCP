#pragma once

#include <imgui_widgets/implot.h>

#include <string>

struct PlotData {
    /* data */
    int maxSize = 1000;
    int offset = 0;
    ImVector<ImVec2> data;

    /* description and control */
    // bool enabled = true;
    std::string name = "";

    PlotData(int maxSize = 1000) : maxSize(maxSize) { data.reserve(maxSize); }

    void addData(float x, float y) {
        if (data.size() < maxSize)
            data.push_back(ImVec2(x, y));
        else {
            data[offset] = ImVec2(x, y);
            offset = (offset + 1) % maxSize;
        }
    }

    void erase() {
        if (data.size() > 0) {
            data.shrink(0);
            offset = 0;
        }
    }

    float xBegin() {
        if (data.size() < maxSize)
            return data[0].x;
        else {
            return data[offset].x;
        }
    }

    float xEnd() {
        if (data.size() < maxSize)
            return data[data.size() - 1].x;
        else {
            return data[offset - 1].x;
        }
    }
};
