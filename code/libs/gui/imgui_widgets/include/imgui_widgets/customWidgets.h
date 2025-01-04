#pragma once

#include <imgui.h>

inline void drawCameraWidget(Camera& camera, double xpos, double ypos, double width, double height){
    glm::mat4 cameraview = camera.getViewMatrix();
    float* cameraviewptr = glm::value_ptr(cameraview);

    ImGuizmo::BeginFrame();
    ImGuizmo::ViewManipulate(cameraviewptr, 1.0,
                             ImVec2(xpos, ypos),
                             ImVec2(width, height), 0x10101010);

    camera.setRotation(glm::transpose(cameraview) );
}

inline bool ToggleButton(const char *str_id, bool *v) {
    bool clicked = false;
    ImVec2 p = ImGui::GetCursorScreenPos();
    ImDrawList *draw_list = ImGui::GetWindowDrawList();

    float height = ImGui::GetFrameHeight();
    float width = height * 1.55f;
    float radius = height * 0.50f;

    ImGui::InvisibleButton(str_id, ImVec2(width, height));
    if (ImGui::IsItemClicked()) {
        clicked = true;
        *v = !*v;
    }

    float t = *v ? 1.0f : 0.0f;

    ImGuiContext &g = *GImGui;
    float ANIM_SPEED = 0.08f;
    if (g.LastActiveId ==
        g.CurrentWindow->GetID(str_id))  // && g.LastActiveIdTimer < ANIM_SPEED)
    {
        float t_anim = ImSaturate(g.LastActiveIdTimer / ANIM_SPEED);
        t = *v ? (t_anim) : (1.0f - t_anim);
    }

    ImU32 col_bg;
    if (ImGui::IsItemHovered())
        col_bg =
                ImGui::GetColorU32(ImLerp(ImVec4(0.78f, 0.78f, 0.78f, 1.0f),
                                          ImVec4(0.64f, 0.83f, 0.34f, 1.0f), t));
    else
        col_bg =
                ImGui::GetColorU32(ImLerp(ImVec4(0.85f, 0.85f, 0.85f, 1.0f),
                                          ImVec4(0.56f, 0.83f, 0.26f, 1.0f), t));

    draw_list->AddRectFilled(p, ImVec2(p.x + width, p.y + height), col_bg,
                             height * 0.5f);
    draw_list->AddCircleFilled(
            ImVec2(p.x + radius + t * (width - radius * 2.0f), p.y + radius),
            radius - 1.5f, IM_COL32(255, 255, 255, 255));

    /*
            if (t < 0.1 || t > 0.9){
                    if (*v)
                            draw_list->AddText(ImVec2(p.x + radius / 2.0, p.y +
       radius / 3.0), ImGui::GetColorU32(ImVec4(0.0f, 0.0f, 0.0f, 1.0f)), "On");
                    else
                            draw_list->AddText(ImVec2(p.x + width - 2 * radius -
       radius / 2.0, p.y + radius / 3.0), ImGui::GetColorU32(ImVec4(0.0f, 0.0f,
       0.0f, 1.0f)), "Off");
            }
    */

//    ImGui::SameLine();
//    ImGui::Text("%s", str_id);
    return clicked;
}

inline bool DoubleBackButton(const char *str_id) {
    bool clicked = false;
    ImVec2 p = ImGui::GetCursorScreenPos();
    ImDrawList *draw_list = ImGui::GetWindowDrawList();

    float height = ImGui::GetFrameHeight();
    float width = height * 1.0f;

    ImGui::InvisibleButton(str_id, ImVec2(width, height));
    if (ImGui::IsItemClicked()) {
        clicked = true;
    }

    ImU32 col_bg;

    if (ImGui::IsItemHovered())
        col_bg = ImGui::GetColorU32(ImVec4(0.75f, 0.75f, 0.75f, 1.0f));
    else
        col_bg = ImGui::GetColorU32(ImVec4(0.85f, 0.85f, 0.85f, 1.0f));


    draw_list->AddRectFilled(p, ImVec2(p.x + width, p.y + height), col_bg,
                             height * 0.05f);

    draw_list->AddTriangleFilled(
            ImVec2(p.x + width * 0.76, p.y + height * 0.2),
            ImVec2(p.x + width * 0.44, p.y + height * 0.5),
            ImVec2(p.x + width * 0.76, p.y + height * 0.8),
            ImGui::GetColorU32(ImVec4(0.0f, 0.0f, 0.0f, 1.0f)));

    draw_list->AddTriangleFilled(
            ImVec2(p.x + width * 0.46, p.y + height * 0.2),
            ImVec2(p.x + width * 0.14, p.y + height * 0.5),
            ImVec2(p.x + width * 0.46, p.y + height * 0.8),
            ImGui::GetColorU32(ImVec4(0.0f, 0.0f, 0.0f, 1.0f)));

    return clicked;
}

inline bool PlayPauseButton(const char *str_id, bool *v) {
    bool clicked = false;
    ImVec2 p = ImGui::GetCursorScreenPos();
    ImDrawList *draw_list = ImGui::GetWindowDrawList();

    float height = ImGui::GetFrameHeight();
    float width = height * 1.0f;

    ImGui::InvisibleButton(str_id, ImVec2(width, height));
    if (ImGui::IsItemClicked()) {
        clicked = true;
        *v = !*v;
    }

    ImU32 col_bg;

    if (*v) {
        if (ImGui::IsItemHovered())
            col_bg = ImGui::GetColorU32(ImVec4(0.54f, 0.73f, 0.3f, 1.0f));
        else
            col_bg = ImGui::GetColorU32(ImVec4(0.64f, 0.83f, 0.34f, 1.0f));
    } else {
        if (ImGui::IsItemHovered())
            col_bg = ImGui::GetColorU32(ImVec4(0.75f, 0.75f, 0.75f, 1.0f));
        else
            col_bg = ImGui::GetColorU32(ImVec4(0.85f, 0.85f, 0.85f, 1.0f));
    }

    draw_list->AddRectFilled(p, ImVec2(p.x + width, p.y + height), col_bg,
                             height * 0.15f);

    if (!*v)
        draw_list->AddTriangleFilled(
                ImVec2(p.x + width * 0.25, p.y + height * 0.2),
                ImVec2(p.x + width * 0.75, p.y + height * 0.5),
                ImVec2(p.x + width * 0.25, p.y + height * 0.8),
                ImGui::GetColorU32(ImVec4(0.0f, 0.0f, 0.0f, 1.0f)));
    else {
        draw_list->AddRectFilled(
                ImVec2(p.x + width * 0.2, p.y + height * 0.2),
                ImVec2(p.x + width * 0.45, p.y + height * 0.8),
                ImGui::GetColorU32(ImVec4(0.0f, 0.0f, 0.0f, 1.0f)), 0);
        draw_list->AddRectFilled(
                ImVec2(p.x + width * 0.55, p.y + height * 0.2),
                ImVec2(p.x + width * 0.8, p.y + height * 0.8),
                ImGui::GetColorU32(ImVec4(0.0f, 0.0f, 0.0f, 1.0f)), 0);
    }

    return clicked;
}

