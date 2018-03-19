/*
 * Texton Noise 2D
 */

// #############################################################################

#include <cassert>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <fstream>
#include <iostream>
#include <GL/glew.h>
#include <GL/glut.h>
#include <AntTweakBar.h>
#include <getopt.h>

#include "io_png/io_png.h"

// #############################################################################

float camera_roi = 1.5f;
float camera_rotation_quat[4] = { 0.0f, 0.0f, 0.0f, 1.0f };

// Global variables for opengl program and shaders
GLuint program;
GLuint vertex_shader;
GLuint fragment_shader;
GLuint fbo;
GLuint fbo_color;
GLuint fbo_depth;
bool offscreen=false; // if offscreen=true, rendering is done in a Framebuffer object and output is saved as a .png image.



// Global variables for texton textures:
const int g_textonTextureUnit = 0; // Texture unit in Graphics card
GLuint g_textonTextureLocation;
GLuint g_textonTextureSampler;



// Procedural noise global parameters
float g_meannbofimpacts = 30.0f; //Can be changed within the interface
GLuint g_meannbofimpactsLocation;
float g_texture_mean[3];
GLuint g_noiseMeanLocation;
unsigned int g_offset=0;
GLuint g_offsetLocation;
float g_sumTexton[3]; // integral of kernel = sum texton for normalization
GLuint g_sumTextonLocation;
unsigned int g_noiseOff = 0u; // value to turn off noise
GLuint g_noiseOffLocation;
unsigned int g_filteringOn = 1u; // 1u: use filtering; 0u: do not use filtering, only local linear interpolation
unsigned int g_filteringOnOld = g_filteringOn;
float g_scale = 128.0f;
GLuint g_scaleLocation;

// Surface parameter for procedural computation of normal vectors
unsigned int g_surfaceType = 0u;
/* Surface Type
 * 0: square 
 * 1: cube
 * 2: sphere
 * 3: torus
 */
GLuint g_surfaceTypeLocation;
unsigned int Nvert[1]; /* number of vertcies */

GLuint vaoID[1];
GLuint vboID[1];


// int g_window_width = 128;
// int g_window_height = 128;
int g_window_width = 700;
int g_window_height = 700;
// int g_window_width = 1024;
// int g_window_height = 1024;
GLuint g_half_window_sizeLocation;

char texton_filename[200] = "../input/discharge_print_s128.texton";
char output_image_name[200] = "../output/out.png";
int waitingTimeBeforeOffscreenRendering = 2500;

// #############################################################################

void TwQuat4fToOpenGLMatrixf(const float* q, GLfloat* m)
{
  float yy2 = 2.0f * q[1] * q[1];
  float xy2 = 2.0f * q[0] * q[1];
  float xz2 = 2.0f * q[0] * q[2];
  float yz2 = 2.0f * q[1] * q[2];
  float zz2 = 2.0f * q[2] * q[2];
  float wz2 = 2.0f * q[3] * q[2];
  float wy2 = 2.0f * q[3] * q[1];
  float wx2 = 2.0f * q[3] * q[0];
  float xx2 = 2.0f * q[0] * q[0];
  m[ 0] = - yy2 - zz2 + 1.0f;
  m[ 1] = xy2 + wz2;
  m[ 2] = xz2 - wy2;
  m[ 3] = 0.0f;
  m[ 4] = xy2 - wz2;
  m[ 5] = - xx2 - zz2 + 1.0f;
  m[ 6] = yz2 + wx2;
  m[ 7] = 0.0f;
  m[ 8] = xz2 + wy2;
  m[ 9] = yz2 - wx2;
  m[10] = - xx2 - yy2 + 1.0f;
  m[11] = 0.0f;
  m[12] = m[13] = m[14] = 0.0f;
  m[15] = 1.0f;
}

GLchar* read_shader_source_from_file(const char* filename)
{
  FILE* fp = std::fopen(filename, "rb");
  if (fp == NULL) {
     std::cerr << " *** error: " << std::strerror(errno) << std::endl;
     std::exit(EXIT_FAILURE);
  }
  std::fseek(fp, 0, SEEK_END);
  std::size_t shader_source_size = std::ftell(fp);
  std::rewind(fp);
  GLchar* shader_source = new GLchar[shader_source_size + 1];
  shader_source[shader_source_size] = 0;
  if(std::fread(shader_source, 1, shader_source_size, fp));
  std::fclose(fp);
  return shader_source;
}

void shader_source_from_file(GLuint shader, const char* filename)
{
  const GLchar* shader_source = read_shader_source_from_file(filename);
  glShaderSource(shader, 1, &shader_source, NULL);
  delete [] shader_source;
}

void shader_info_log(GLuint shader)
{
  GLint info_log_length;
  glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &info_log_length);
  if (info_log_length > 0) {
    char* info_log = new char[info_log_length];
    glGetShaderInfoLog(shader, info_log_length, NULL, info_log);
    std::cout << info_log;
    delete [] info_log;
  }
}

void compile_shader(GLuint shader)
{
  glCompileShader(shader);
  GLint compile_status;
  glGetShaderiv(shader, GL_COMPILE_STATUS, &compile_status);
  if (compile_status == GL_FALSE) {
    std::cerr << " *** error:" << std::endl;
    shader_info_log(shader);
    std::exit(EXIT_FAILURE);
  }
}

void program_info_log(GLuint program)
{
  GLint info_log_length;
  glGetProgramiv(program, GL_INFO_LOG_LENGTH, &info_log_length);
  if (info_log_length > 0) {
    char* info_log = new char[info_log_length];
    glGetProgramInfoLog(program, info_log_length, NULL, info_log);
    std::cout << info_log;
    delete [] info_log;
  }
}

void link_program(GLuint program)
{
  glLinkProgram(program);
  GLint link_status;
  glGetProgramiv(program, GL_LINK_STATUS, &link_status);
  if (link_status == GL_FALSE) {
    std::cerr << " *** error:" << std::endl;
    program_info_log(program);
    std::exit(EXIT_FAILURE);
  }
}



//-------------------------------------------------------------------------
// Calculates the frames per second
// From: http://mycodelog.com/2010/04/16/fps/
//-------------------------------------------------------------------------
int frameCount = 0;
int currentTime = 0;
int previousTime = 0;
int totalTime = 0;
float fps = 0.0f;
float time_per_frame = 0.0f;
void calculateFPS()
{
    //  Increase frame count
    frameCount++;
 
    //  Get the number of milliseconds since glutInit called
    //  (or first call to glutGet(GLUT ELAPSED TIME)).
    currentTime = glutGet(GLUT_ELAPSED_TIME);
 
    //  Calculate time passed
    unsigned int timeInterval = currentTime - previousTime;
 
    if(timeInterval > 1000u)
    {
        //  calculate the number of frames per second
        fps = 1000.0f * frameCount / (float) timeInterval;
 
        time_per_frame = (float) timeInterval / frameCount;
        
        //  Set time
        previousTime = currentTime;
 
        //  Reset frame count
        frameCount = 0;
        printf("Framerate: %3.2f images/sec ; time per frame is: %3.2f ms\n", fps, time_per_frame);
    }
}


// #############################################################################
/**
 * This function check_black_frame check that the lines at the border of an interlaced RGB image is zero.
 * @return 1 if the border is (0,0,0), 0 otherwise.
 */
int check_black_frame(float *im, size_t w, size_t h)
{
  unsigned int i, j;
  
  for(i=0u;i<3u*w;i++)
  {
    if((im[i] != 0.0f)||(im[i+3u*w*(h-1u)] != 0.0f)) return(0);
  }
  for(j=0u;j<h;j++)
  {
    for(i=0u;i<3u;i++)
    {
      if((im[i+3*w*j] != 0.0f)||(im[i+3u*(w-1u)+3u*w*j] != 0.0f)) return(0);
    }
  }
  return(1);
}

// #############################################################################
/**
 *  This function load a .texton file
 */

float *load_texton_file(const char* filename, unsigned int* order, size_t* w, size_t* h, float texture_mean[3])
{
  
  float *texton=NULL;
  std::ifstream istream(filename);
  assert(istream.is_open());
  assert(istream);
  
  /* check filetype */
  std::string filetype;
  istream.width(6);
  istream >> filetype;
//   std::cout << filetype << std::endl;
  assert(filetype == "TEXTON");
  
  int version;
  istream >> version;
  if(version != 1) std::cout << "Error: version = " << version << " is not supported" << std::endl;
  
  
  istream >> *order;
//   std::cout << "Interpolation order is: " << *order << std::endl;
  
  istream >> texture_mean[0] >> texture_mean[1] >> texture_mean[2];
  std::cout << "TEXTON load info: Mean of texture is: " << texture_mean[0] << ", " << texture_mean[1] << ", " << texture_mean[2] << std::endl;
  istream >> *w >> *h;
  std::cout << "TEXTON load info: Texton size is: " << *w << " x " << *h << std::endl;
  
  /* Allocate and fill array of texton coefficients */
  texton = (float *) malloc(3*(*w)*(*h)*sizeof(float));
  for(unsigned int i=0; i < (*w)*(*h); i++)
    istream >> texton[3*i] >> texton[1+3*i] >> texton[2+3*i];
  return(texton);
}


// #############################################################################
/**
 * This functions transforms pixel data into an array of desinterlaced rgb image before writing on disc
 *  1) desinterlace RGB channels
 *  2) flip y-axis
 * 
 */
void from_gldata_to_desinterlace_rgb_image(unsigned char *im, size_t w, size_t h)
{
    unsigned char *temp;
    temp = (unsigned char *) malloc( 3*w*h*sizeof(unsigned char));
    size_t imsize = w*h;
    
    // 1) desinterlace RGB channels
    for(unsigned int i=0u; i<imsize; i++)
    {
        temp[i] = im[3*i];
        temp[i+imsize] = im[3*i+1];
        temp[i+2*imsize] = im[3*i+2];
    }
    
    // 2) flip y-axis
    for(unsigned int i=0u; i<w; i++)
    {
      for(unsigned int j=0u; j<h; j++)
      {  
        im[i+j*w] = temp[i+(h-1-j)*w];
        im[i+j*w+imsize] = temp[i+(h-1-j)*w+imsize];
        im[i+j*w+2*imsize] = temp[i+(h-1-j)*w+2*imsize];
      }
    }
    
    free(temp);
}



// #############################################################################


/**
 * This functions transforms pixel data into an array of desinterlaced rgb image before writing on disc
 *  1) desinterlace RGB channels
 *  2) flip y-axis
 * 
 */
void from_gldata_with_alpha_to_desinterlace_rgb_image(unsigned char *im, size_t w, size_t h)
{
    unsigned char *temp;
    temp = (unsigned char *) malloc( 3*w*h*sizeof(unsigned char));
    size_t imsize = w*h;
    
    // 1) desinterlace RGB channels
    for(unsigned int i=0u; i<imsize; i++)
    {
        temp[i] = im[4*i];
        temp[i+imsize] = im[4*i+1];
        temp[i+2*imsize] = im[4*i+2];
    }
    
    // 2) flip y-axis
    for(unsigned int i=0u; i<w; i++)
    {
      for(unsigned int j=0u; j<h; j++)
      {  
        im[i+j*w] = temp[i+(h-1-j)*w];
        im[i+j*w+imsize] = temp[i+(h-1-j)*w+imsize];
        im[i+j*w+2*imsize] = temp[i+(h-1-j)*w+2*imsize];
      }
    }
    
    free(temp);
}

// #############################################################################
void save_last_display()
{
 
  unsigned char *pixel_data; // pointer for pixel data to write on disc as png filename
  // Allocate memory for pixel_data
  pixel_data = (unsigned char *) malloc(3*g_window_width*g_window_height*sizeof(unsigned char));
  // Read pixel data
  glReadPixels(0, 0, g_window_width, g_window_height, GL_RGB, GL_UNSIGNED_BYTE, pixel_data);
  // Desinterlace image
  from_gldata_to_desinterlace_rgb_image(pixel_data, g_window_width, g_window_height);
  io_png_write_u8(output_image_name, pixel_data, g_window_width, g_window_height, 3);
}

// #############################################################################
/** 
 * 
 * This function save_fbo_image() save the content of the color texture of the framebuffer object fbo into a .png
 * image.
 * 
 */
int save_fbo_image()
{
  unsigned char *pixel_data; // pointer for pixel data to write on disc as png filename
  // Allocate memory for pixel_data
  pixel_data = (unsigned char *) malloc(4*g_window_width*g_window_height*sizeof(unsigned char));
  glBindTexture(GL_TEXTURE_2D, 0);
  glBindTexture(GL_TEXTURE_2D, fbo_color);
  glGetTexImage(GL_TEXTURE_2D, 0, GL_RGBA, GL_UNSIGNED_BYTE, pixel_data);
  // Desinterlace image
  from_gldata_with_alpha_to_desinterlace_rgb_image(pixel_data, g_window_width, g_window_height);
  io_png_write_u8(output_image_name, pixel_data, g_window_width, g_window_height, 3);

  return EXIT_SUCCESS;
}

// #############################################################################
void CHECK_FRAMEBUFFER_STATUS()
{                                                         
  GLenum status;
  status = glCheckFramebufferStatus(GL_DRAW_FRAMEBUFFER); 
  switch(status) {
  case GL_FRAMEBUFFER_COMPLETE:
    break;

  case GL_FRAMEBUFFER_UNSUPPORTED:
  /* choose different formats */
    break;

  default:
    /* programming error; will fail on all hardware */
    fputs("Framebuffer Error\n", stderr);
    exit(-1);
  }
}


// #############################################################################

void init()
{
  

  
  // ---------------------------------------------------------------------------
  /* loading shaders */
  
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_TEXTURE_2D); // Enable texturing so we can bind our frame buffer texture  
  program = glCreateProgram();
  assert(glIsProgram(program) == GL_TRUE);
  vertex_shader = glCreateShader(GL_VERTEX_SHADER);
  glAttachShader(program, vertex_shader);
  shader_source_from_file(vertex_shader, "texton_noise_2d.vert");
  compile_shader(vertex_shader);
  fragment_shader = glCreateShader(GL_FRAGMENT_SHADER);
  glAttachShader(program, fragment_shader);
  shader_source_from_file(fragment_shader, "texton_noise_2d.frag");
  compile_shader(fragment_shader);
  link_program(program);

  if(offscreen)
  {

    /* Generate Framebuffer Objects for offscreen rendering */
    glGenFramebuffers(1, &fbo);
    glGenTextures(1, &fbo_color);
    glGenRenderbuffers(1, &fbo_depth);
    
    glUseProgram(program);
    glBindFramebuffer(GL_FRAMEBUFFER, fbo);
    
    glBindTexture(GL_TEXTURE_2D, fbo_color);
    glTexImage2D(GL_TEXTURE_2D,0,GL_RGBA,g_window_width, g_window_height,0,GL_RGBA,GL_UNSIGNED_BYTE,NULL);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, fbo_color, 0);
    
    glBindRenderbuffer(GL_RENDERBUFFER, fbo_depth);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT24, g_window_width, g_window_height);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, fbo_depth);
    
    CHECK_FRAMEBUFFER_STATUS();

  }
  
  // ---------------------------------------------------------------------------
  /* create OpenGL texture and sampler */
  std::cout << "info: " << "TEST LOAD TEXTON FILE" << std::endl;
  float *texton;
  size_t wtexton, htexton;
  unsigned int interp_order;
  texton = load_texton_file(texton_filename, &interp_order, &wtexton, &htexton, g_texture_mean);
  std::cout << "info: " << "LOAD TEXTON FILE DONE" << std::endl << std::endl << std::endl;
  if(interp_order==1)
  {
    if(0==check_black_frame(texton, wtexton, htexton))
        std::cout << "WARNING: " << "FOR EXACT POISSON SAMPLING THE KERNEL MUST BE ZERO AT THE BORDER" << std::endl;
  }

  // Setting uniform variable for g_meannbofimpacts
  g_meannbofimpactsLocation = glGetUniformLocation(program, "meannbofimpacts");
  glUseProgram(program);
  glUniform1f(g_meannbofimpactsLocation, g_meannbofimpacts);
  
  // Setting uniform variable for mean value
  g_noiseMeanLocation = glGetUniformLocation(program, "noiseMean");
  glUseProgram(program);
  glUniform3f(g_noiseMeanLocation, g_texture_mean[0], g_texture_mean[1], g_texture_mean[2]);

  // Setting uniform variable for g_scale
  g_scaleLocation = glGetUniformLocation(program, "scale");
  glUseProgram(program);
  glUniform1f(g_scaleLocation, g_scale);
  
  // Setting offset
  g_offsetLocation = glGetUniformLocation(program, "offset");
  g_offset = (unsigned int) time((long int) 0);
  glUseProgram(program);
  glUniform1ui(g_offsetLocation, g_offset);
  
  // Setting uniform variable for sum of coefficients
  g_sumTexton[0] = 0.0f;
  g_sumTexton[1] = 0.0f;
  g_sumTexton[2] = 0.0f;
  for(size_t i=0; i<wtexton*htexton; ++i)
  {
    g_sumTexton[0] += texton[3*i];
    g_sumTexton[1] += texton[1+3*i];
    g_sumTexton[2] += texton[2+3*i];
  }
  g_sumTextonLocation = glGetUniformLocation(program, "sumTexton");
  glUseProgram(program);
  glUniform3f(g_sumTextonLocation, g_sumTexton[0], g_sumTexton[1], g_sumTexton[2]);
  
  // Setting uniform variable for g_noiseOff
  g_noiseOffLocation = glGetUniformLocation(program, "noiseOff");
  glUseProgram(program);
  glUniform1ui(g_noiseOffLocation, g_noiseOff);

  // Setting uniform texture variables
  GLuint textonTextureUniform = glGetUniformLocation(program, "textonTexture"); // textonTexture is the name of the sampler for the fragment shader
  glUseProgram(program);
  glUniform1i(textonTextureUniform, g_textonTextureUnit);
  glUseProgram(0);
  
  // Texture generation with filtering
  glGenTextures(1, &g_textonTextureLocation);
  glBindTexture(GL_TEXTURE_2D, g_textonTextureLocation);
  glTexImage2D(GL_TEXTURE_2D, 0,  GL_RGB32F, wtexton, htexton, 0, GL_RGB, GL_FLOAT, texton);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_BASE_LEVEL, 0u);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL, (unsigned int) floor(log2((double) wtexton)));
  glGenerateMipmap(GL_TEXTURE_2D);
  glBindTexture(GL_TEXTURE_2D, 0);
  
  
  // Sampler generation and option setting:
  glGenSamplers(1, &g_textonTextureSampler); // Generate sampler 
  if(interp_order==0) // Nearest neighbor interpolation
  {
    glSamplerParameteri(g_textonTextureSampler, GL_TEXTURE_MAG_FILTER, GL_NEAREST); // Setting evaluation for magnification: GL_NEAREST=Pointwise evaluation
    glSamplerParameteri(g_textonTextureSampler, GL_TEXTURE_MIN_FILTER, GL_NEAREST); // Setting evaluation for minification: GL_NEAREST=Pointwise evaluation
  }
  else // Bilinear interpolation
  {
    glSamplerParameteri(g_textonTextureSampler, GL_TEXTURE_MAG_FILTER, GL_LINEAR); // Setting evaluation for magnification: GL_NEAREST=Pointwise evaluation
    glSamplerParameteri(g_textonTextureSampler, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR); // Setting evaluation for minification: GL_NEAREST=Pointwise evaluation
    GLfloat maxAniso = 0.0f;
    glGetFloatv(GL_MAX_TEXTURE_MAX_ANISOTROPY_EXT, &maxAniso);
    std::cout << "info: Anisotropic filtering with maxAniso = " << maxAniso << std::endl;
    glSamplerParameterf(g_textonTextureSampler, GL_TEXTURE_MAX_ANISOTROPY_EXT, maxAniso);
  }
  // GL_CLAMP_TO_BORDER: set (0,0,0,0) outside of texture coordinates
  glSamplerParameteri(g_textonTextureSampler, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER); // Setting option for out of range evaluation for 1st variable S (2nd is T), GL_CLAMP_TO_EDGE, can also be GL_REPEAT for periodic repetition
  glSamplerParameteri(g_textonTextureSampler, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
  
  // ---------------------------------------------------------------------------
  // Load the plane mesh using VAO
  float geom_data[18] = {
    1.0f, 1.0f, 0.0f,
    -1.0f, 1.0f, 0.0f,
    1.0f, -1.0f, 0.0f,
    -1.0f, -1.0f, 0.0f,
    -1.0f, 1.0f, 0.0f,
    1.0f, -1.0f, 0.0f
  };
  Nvert[0] = 6;
  glGenVertexArrays(1, &vaoID[0]); // Create Vertex Array Object
  glBindVertexArray(vaoID[0]);     // Bind Vertex Array Object
  
  glGenBuffers(1, &vboID[0]); // Generate Vertex Buffer Object  
  glBindBuffer(GL_ARRAY_BUFFER, vboID[0]); // Bind Vertex Buffer Object  
  glBufferData(GL_ARRAY_BUFFER, 3 * sizeof(GLfloat) * Nvert[0], geom_data, GL_STATIC_DRAW); // load data
  
  // Enable the vertex array attributes.
  glEnableVertexAttribArray(0);  // Vertex position
  
  // Specify the location and format of the vertex position portion of the vertex buffer.
  glBindBuffer(GL_ARRAY_BUFFER, vboID[0]);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
  
  glBindVertexArray(0); // Disable Vertex Buffer Object
  
}

// -----------------------------------------------------------------------------

void reshape(int width, int height)
{
  TwWindowSize(width, height);
  glViewport(0, 0, width, height);
}

// -----------------------------------------------------------------------------
int frame=0;
int mytime=0;
int timebase=0;

void display()
{
  
  glUseProgram(program);

  //Before drawing: render into fbo instead of screen
  if(offscreen)
  {
    glBindTexture(GL_TEXTURE_2D, 0);
    glEnable(GL_TEXTURE_2D);
    glBindFramebuffer(GL_FRAMEBUFFER, fbo);
    glViewport(0,0, g_window_width, g_window_height);
  }
  else
  {
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
  }

  glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  float camera_fovy = 45.0f;
  float camera_aspect = (float) g_window_width / (float) g_window_height;
  float camera_eye_z = camera_roi * (1.0f / std::sin((camera_fovy * (M_PI / 180.0f)) / 2.0f));
  float camera_z_near = camera_eye_z - camera_roi;
  float camera_z_far = camera_eye_z + camera_roi;
  gluPerspective(camera_fovy, camera_aspect, camera_z_near, camera_z_far);  
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  gluLookAt(0.0f, 0.0f, camera_eye_z, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f);  

  GLfloat m[16];
  TwQuat4fToOpenGLMatrixf(camera_rotation_quat, m);
  glMultMatrixf(m);

  


  // Update meannbofimpacts
  glUniform1f(g_meannbofimpactsLocation, g_meannbofimpacts);
  
  // Update noiseOff
  glUniform1ui(g_noiseOffLocation, g_noiseOff);
  
  // Update scale
  glUseProgram(program);
  glUniform1f(g_scaleLocation, g_scale);
    
  // Update window size
  glUseProgram(program);
  glUniform2f(g_half_window_sizeLocation, (float) g_window_width / 2.0f, (float) g_window_height / 2.0f);
  
  // Activate texture depending on filtering: If the filtering switch is switched, the texture is reloaded.
  if(g_filteringOn!=g_filteringOnOld)
  {
    size_t wtexton, htexton;
    unsigned int interp_order;
    float *texton = load_texton_file(texton_filename, &interp_order, &wtexton, &htexton, g_texture_mean);
    if(texton==NULL) std::cout << "Probem while loading texton" << std::endl;
    if(g_filteringOn)
    { 
      std::cout << "Filtering ON" << std::endl;
      glBindTexture(GL_TEXTURE_2D, g_textonTextureLocation);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_BASE_LEVEL, 0u);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL, (unsigned int) floor(log2((double) wtexton)));
      glGenerateMipmap(GL_TEXTURE_2D);
      glBindTexture(GL_TEXTURE_2D, 0);
      /* sampler parameters */
      glSamplerParameteri(g_textonTextureSampler, GL_TEXTURE_MAG_FILTER, GL_LINEAR); // Setting evaluation for magnification: GL_NEAREST=Pointwise evaluation
      glSamplerParameteri(g_textonTextureSampler, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR); // Setting evaluation for minification: GL_NEAREST=Pointwise evaluation
      GLfloat maxAniso = 0.0f;
      glGetFloatv(GL_MAX_TEXTURE_MAX_ANISOTROPY_EXT, &maxAniso);
      glSamplerParameterf(g_textonTextureSampler, GL_TEXTURE_MAX_ANISOTROPY_EXT, maxAniso);
    }
    else
    {
      std::cout << "Filtering OFF" << std::endl;
      glBindTexture(GL_TEXTURE_2D, g_textonTextureLocation);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_BASE_LEVEL, 0u);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL, 0u);
      glBindTexture(GL_TEXTURE_2D, 0);
      glSamplerParameteri(g_textonTextureSampler, GL_TEXTURE_MAG_FILTER, GL_LINEAR); // Setting evaluation for magnification: GL_NEAREST=Pointwise evaluation
      glSamplerParameteri(g_textonTextureSampler, GL_TEXTURE_MIN_FILTER, GL_LINEAR); // Setting evaluation for minification: GL_NEAREST=Pointwise
      GLfloat maxAniso = 1.0f;
      glSamplerParameterf(g_textonTextureSampler, GL_TEXTURE_MAX_ANISOTROPY_EXT, maxAniso);     
    }
    g_filteringOnOld = g_filteringOn;
  }
  glActiveTexture(GL_TEXTURE0 + g_textonTextureUnit);  // Activate the texture
  glBindTexture(GL_TEXTURE_2D, g_textonTextureLocation);   // Bind the texture
  glBindSampler(g_textonTextureUnit, g_textonTextureSampler);  // Bind the sampler

  // rendering functions
  glBindVertexArray(vaoID[g_surfaceType]);
  glDrawArrays(GL_TRIANGLES, 0, Nvert[g_surfaceType]);
  
  glBindVertexArray(0); // Unbind our Vertex Array Object 

  // Texture cleanup
  glBindSampler(g_textonTextureUnit, 0); 
  glBindTexture(GL_TEXTURE_2D, 0);
  
  
  glUseProgram(0);
  
  if(offscreen)
  {
    glFinish();
    if(glutGet(GLUT_ELAPSED_TIME)>waitingTimeBeforeOffscreenRendering)
    {
      save_fbo_image();
      std::cout << "Image saved! Exit" << std::endl;
      std::exit(EXIT_SUCCESS);  
    }
    
  }
  else
  {
    TwDraw();
    glutSwapBuffers();
  }

}

// -----------------------------------------------------------------------------

void idle()
{ 
  calculateFPS();
  if(offscreen)
  {
      display();
  }
  else
  {
    g_window_width = glutGet(GLUT_WINDOW_WIDTH);
    g_window_height = glutGet(GLUT_WINDOW_HEIGHT);
    glutPostRedisplay();
  }

}

// -----------------------------------------------------------------------------

void mouse(int button, int state, int x, int y)
{
  TwEventMouseButtonGLUT(button, state, x, y);
}

// -----------------------------------------------------------------------------

void motion(int x, int y)
{
  TwEventMouseMotionGLUT(x, y);
}

// -----------------------------------------------------------------------------

void passive_motion(int x, int y)
{
  TwEventMouseMotionGLUT(x, y);
}

// -----------------------------------------------------------------------------

void keyboard(unsigned char key, int x, int y)
{
    if(!TwEventKeyboardGLUT(key, x, y)) {
      if (key == 27) // Escape = quit
      { 
        std::exit(EXIT_SUCCESS);
      }
      else if (key == 113) //Q = Write rendered image on disc and quit
      {
        save_last_display();
        std::exit(EXIT_SUCCESS);
      }
      else if (key == 115) //S = Write rendered image on disc and continu
      {
        save_last_display();
      }
    }
}

// -----------------------------------------------------------------------------

void special(int key, int x, int y)
{
  TwEventSpecialGLUT(key, x, y);
}

// -----------------------------------------------------------------------------

void exit()
{
  if(offscreen==false) TwTerminate();
}

// #############################################################################

void display_usage()
{
  std::cout << "Usage: ./texton_noise_2d" << std::endl;
  std::cout << "   or: ./texton_noise_2d [OPTIONS] -t inputtexton -o outputtexture\n" << std::endl;
  std::cout << "Optional arguments:" << std::endl;
  std::cout << "  -t texton   texton file" << std::endl;
  std::cout << "  -o output   output file" << std::endl;
  std::cout << "  -i mni      mean number of impacts" << std::endl;
  std::cout << "  -x sx       output width" << std::endl;
  std::cout << "  -y sy       output height" << std::endl;
  std::cout << "  -r          offscreen rendering" << std::endl;
  std::cout << "  -T time     waiting time before offscreen rendering " << std::endl;
  std::cout << "  -h          display usage\n" << std::endl;
}


int main(int argc, char* argv[])
{
  
  // ---------------------------------------------------------------------------
  std::cout << "\nTexton Noise 2D" << std::endl;
  glutInit(&argc, argv);  
  
  // ---------------------------------------------------------------------------
  /* Process input argument using getopt */
  int c; /* getopt flag */
  while ((c = getopt(argc, argv, "t:o:i:x:y:rT:h")) != -1) 
  {
    switch (c) 
    {
      case 't':
        /* texton file specified */
        if(strlen(optarg)>=200)
        {
          std::cout << "ERROR: file paths are limited to 200 characters" << std::endl;
          return(-1);
        }
        strcpy(texton_filename, optarg);
        break;
      case 'o':
        /* output filename for .png image specified */
        if(strlen(optarg)>=200)
        {
          std::cout << "ERROR: file paths are limited to 200 characters" << std::endl;
          return(-1);
        }
        strcpy(output_image_name, optarg);
        break;        
      case 'i':
        /* mean number of impact specified */
        g_meannbofimpacts = atof(optarg);
        break;
      case 'x':
        /* output width specified */
        g_window_width = atoi(optarg);
        break;
      case 'y':
        /* output height specified */
        g_window_height = atoi(optarg);
        break;
      case 'r':
        /* offscreen rendering */
        offscreen = true;
        break;
      case 'T':
        /* set waiting time before offscreen rendering (set 0 for immediate rendering) */
        waitingTimeBeforeOffscreenRendering = atoi(optarg);
        break;
      case 'h':
        /* display usage */
        display_usage();
        return (-1);
      default:
        abort();
    }
  }
  
  // ---------------------------------------------------------------------------
  /* Create window */
  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);  
  if(offscreen){
//     glutInitWindowSize(g_window_width, g_window_height);
    glutInitWindowSize(1, 1);
    glutCreateWindow("Hidden window");
  }
  else
  {
    /* use double buffer for display */
    glutInitWindowSize(g_window_width, g_window_height);
    glutCreateWindow("Texton Noise 2D");
  }
  // ---------------------------------------------------------------------------

  std::cout << "info: " << "GL_VENDOR: " << glGetString(GL_VENDOR) << std::endl;
  std::cout << "info: " << "GL_RENDERER: " << glGetString(GL_RENDERER) << std::endl;
  std::cout << "info: " << "GL_VERSION: " << glGetString(GL_VERSION) << std::endl;
  std::cout << "info: " << "GL_SHADING_LANGUAGE_VERSION: " << glGetString(GL_SHADING_LANGUAGE_VERSION) << std::endl;

  // ---------------------------------------------------------------------------

  GLenum return_value = glewInit();
  if (return_value != GLEW_OK) {
    std::cerr << " *** error: " << "glewInit() failed" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  std::cout << "info: " << "GLEW_VERSION: " << glewGetString(GLEW_VERSION) << std::endl;
  
  if(offscreen)
    glutHideWindow();
  
  // ---------------------------------------------------------------------------
  
  init();

  // ---------------------------------------------------------------------------
  /* Set functions for glutMainLoop(); */
  if(offscreen)
  {
    /* no interaction with offsreen rendering */
    glutDisplayFunc(display);
    glutIdleFunc(idle);
  }
  else
  {
    glutReshapeFunc(reshape);
    glutDisplayFunc(display);
    glutIdleFunc(idle);
    glutMouseFunc(mouse);  
    glutMotionFunc(motion);
    glutPassiveMotionFunc(passive_motion);
    glutKeyboardFunc(keyboard);
    glutSpecialFunc(special);
  }
  
  // ---------------------------------------------------------------------------

  std::atexit(exit);

  // ---------------------------------------------------------------------------
  // Set up AntTweakBar menu
  
  if(offscreen==false) /* no interaction with offsreen rendering */
  {
    TwInit(TW_OPENGL, NULL);
    TwWindowSize(glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT)); 
    TwGLUTModifiersFunc(glutGetModifiers);
    TwBar* bar = TwNewBar("Texton Noise 2D");

    TwAddVarRO(bar, "framerate", TW_TYPE_FLOAT, &fps, "group='Performance' label='framerate'");
    TwAddVarRO(bar, "time_per_frame", TW_TYPE_FLOAT, &time_per_frame, "group='Performance' label='ms per frame'");
    
    TwAddVarRW(bar, "camera_rotation", TW_TYPE_QUAT4F, &camera_rotation_quat, "group='Camera' label='tumble'");
    TwAddVarRW(bar, "camera_roi", TW_TYPE_FLOAT, &camera_roi, "group='Camera' label='dolly' min='1.0' max='4.00' step='0.01'");
    
    TwAddVarRW(bar, "meannbofimpacts", TW_TYPE_FLOAT, &g_meannbofimpacts, "group='Noise Parameters' label='mean nb of impacts' min='0.1' max='100' step='0.1'");
    TwAddVarRW(bar, "scale", TW_TYPE_FLOAT, &g_scale, "group='Noise Parameters' label='scale' min='10.0' max='1024.' step='0.5'");
    TwEnumVal noiseOff_enum_val[] = {
      { 0u, "on" },
      { 1u, "off" }
    };
    TwType noiseOff_type = TwDefineEnum("noiseOff", noiseOff_enum_val, 2);
    TwAddVarRW(bar, "noiseOff", noiseOff_type, &g_noiseOff, "group='Noise Parameters' label='noise switch'");
    TwEnumVal filteringOn_enum_val[] = {
      { 1u, "on" },
      { 0u, "off" }
    };
    TwType filteringOn_type = TwDefineEnum("filteringOn", filteringOn_enum_val, 2);
    TwAddVarRW(bar, "filteringOn", filteringOn_type, &g_filteringOn, "group='Noise Parameters' label='filtering switch'");
  }
  
  previousTime = glutGet(GLUT_ELAPSED_TIME);
  glutMainLoop();

  return EXIT_SUCCESS;
}

// #############################################################################
