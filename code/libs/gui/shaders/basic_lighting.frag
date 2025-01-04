#version 330 core

out vec4 FragColor;
uniform vec3 objectColor;
uniform float alpha;
uniform sampler2D texture_diffuse1;
in vec2 TexCoords;
uniform bool use_textures;

#include "compute_shading.frag"

void main()
{
	if (use_textures == false){
		FragColor = vec4(computeBasicShading() * objectColor, alpha);
	}else{
		FragColor = vec4(computeBasicShading(), alpha) * texture(texture_diffuse1, TexCoords);
	}
} 
