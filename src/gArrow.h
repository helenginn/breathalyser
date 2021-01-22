#ifndef __breathalyser_gArrow__
#define __breathalyser_gArrow__

#include <string>

inline std::string Arrow_gsh()
{
	std::string str = 
	"#version 330 core\n"\
	"layout (lines) in;\n"\
	"layout (triangle_strip, max_vertices = 15) out;\n"\
	"\n"\
	"in vec4 vColor[];\n"\
	"in vec4 vPos[];\n"\
	"in vec2 vTex[];\n"\
	"\n"\
	"out vec4 fColor;\n"\
	"\n"\
	"void main()\n"\
	"{\n"\
	"    vec4 diff = gl_in[1].gl_Position - gl_in[0].gl_Position;\n"\
	"    vec2 axis = normalize(vec2(diff[0], diff[1]));\n"\
	"    mat2 rot = mat2(axis[0], axis[1],\n"\
	"                    axis[1], -axis[0]);\n"\
	"	 vec2 x = rot * vec2(1., 0.);\n"\
	"	 vec2 y = rot * vec2(0., 1.);\n"\
	"	 vec4 gl0 = gl_in[0].gl_Position;\n"\
	"	 vec4 gl1 = gl_in[1].gl_Position;\n"\
	"	 fColor = vec4(0., 0., 0., 1.);\n"\
	"\n"\
	"    gl_Position = gl0 + vec4(0.0, 0.0, 0.0, 0.0); \n"\
    "    EmitVertex();\n"\
	"    \n"\
    "    gl_Position = gl0 + vec4(5.*y[0], 5.*y[1], 0.0, 0.0);\n"\
    "    EmitVertex();\n"\
	"    gl_Position = gl1 + vec4(5.*y[0], 5.*y[1], 0.0, 0.0); \n"\
    "    EmitVertex();\n"\
	"    gl_Position = gl0 + vec4(0.0, 0.0, 0.0, 0.0); \n"\
    "    EmitVertex();\n"\
	"    gl_Position = gl1 + vec4(5.*y[0], 5.*y[1], 0.0, 0.0); \n"\
    "    EmitVertex();\n"\
	"    gl_Position = gl1 + vec4(-5.*y[0], -5.*y[1], 0., 0.); \n"\
    "    EmitVertex();\n"\
	"    gl_Position = gl0 + vec4(0.0, 0.0, 0.0, 0.0); \n"\
    "    EmitVertex();\n"\
	"    gl_Position = gl1 + vec4(-5.*y[0], -5.*y[1], 0., 0.); \n"\
    "    EmitVertex();\n"\
    "    gl_Position = gl0 + vec4(-5.*y[0], -5.*y[1], 0., 0.);\n"\
    "    EmitVertex();\n"\
	"    gl_Position = gl0 + vec4(0.0, 0.0, 0.0, 0.0); \n"\
    "    EmitVertex();\n"\
	"    gl_Position = gl0 + vec4(-10.*y[0], -10.*y[1], 0., 0.); \n"\
    "    EmitVertex();\n"\
    "    gl_Position = gl0 + vec4(-10.*x[0], -10.*x[1], 0., 0.);\n"\
    "    EmitVertex();\n"\
	"    gl_Position = gl0 + vec4(0.0, 0.0, 0.0, 0.0); \n"\
    "    EmitVertex();\n"\
	"    gl_Position = gl0 + vec4(10.*y[0], 10.*y[1], 0., 0.); \n"\
    "    EmitVertex();\n"\
    "    gl_Position = gl0 + vec4(-10.*x[0], -10.*x[1], 0., 0.);\n"\
    "    EmitVertex();\n"\
	"    \n"\
	"}";
	return str;
}

#endif
