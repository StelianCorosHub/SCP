#pragma once
#include "imgui.h"

namespace ImGui {

// create a double slider function
IMGUI_API bool SliderDouble(const char *label, double *v, double v_min,
                            double v_max, const char *format = NULL,
                            ImGuiSliderFlags flags = 0);
};  // namespace ImGui