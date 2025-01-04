#include <imgui_widgets/SliderDouble.h>

IMGUI_API bool ImGui::SliderDouble(const char *label, double *v, double v_min,
                                   double v_max, const char *format,
                                   ImGuiSliderFlags flags) {
    return ImGui::SliderScalar(label, ImGuiDataType_Double, v, &v_min, &v_max,
                               format, flags);
}