#version 330 core
layout (location = 0) in vec3 aPos;
layout (location = 1) in vec3 aNormal;
layout (location = 2) in vec2 aTexCoords;

out vec3 FragPos;
out vec3 Normal;
out vec4 lightSpacePos;
out vec2 TexCoords;
out float zVal;
out float xVal;
out float yVal;

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;

uniform mat4 lightView;
uniform mat4 lightProjection;


void main()
{
	zVal = aPos.z;
	xVal = aPos.x;
	yVal = aPos.y;
	FragPos = vec3(model * vec4(aPos, 1.0));
	Normal = vec3(transpose(inverse(model)) * vec4(aNormal, 0));
    TexCoords = aTexCoords;
	
    gl_Position = projection * view * vec4(FragPos, 1.0);
	lightSpacePos = lightProjection * lightView * model * vec4(aPos, 1.0);
}

