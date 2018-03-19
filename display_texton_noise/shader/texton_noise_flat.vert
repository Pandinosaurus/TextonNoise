#version 150 compatibility
#extension GL_ARB_explicit_attrib_location : enable

layout (location = 0) in vec3 position;

out vec4 gl_Position;
out vec3 x_tex;

void main()
{
  gl_Position = vec4(position, 1.0f);
  x_tex = position;
}
