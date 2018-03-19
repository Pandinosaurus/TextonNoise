
#version 150 compatibility
#extension GL_ARB_gpu_shader5 : enable
#extension GL_NV_shader_buffer_load : enable

const float pi = 3.14159265358979323846f;

/* 
 * From http://www.reedbeta.com/blog/2013/01/12/quick-and-easy-gpu-random-numbers-in-d3d11/
 * Same strategy as in Gabor noise by example
 * Apply hashtable to create cellseed
 * Use a linear congruential generator as fast PRNG
 */
uint wang_hash(uint seed)
{
  seed = (seed ^ 61u) ^ (seed >> 16u);
  seed *= 9u;
  seed = seed ^ (seed >> 4u);
  seed *= 668265261u;
  seed = seed ^ (seed >> 15u);
  return(seed);
}

uint cellseed(const in ivec2 c, const in uint offset)
{
  const uint period = 1024u;
  uint s = ((uint(c.y) % period) * period + (uint(c.x) % period))*period + offset;
  if (s == 0u) s = 1u;
  return(s);
}

struct noise_prng
{
  uint state;
};

void mysrand(inout noise_prng p, const in uint seed)
{
    uint s=seed;
    p.state = wang_hash(s);
}

uint myrand(inout noise_prng p)
{
// linear congruential generator: procudes correlated output. Similar patterns are visible
// p.state = 1664525u * p.state + 1013904223u;
// Xorshift algorithm from George Marsaglia's paper
  p.state ^= (p.state << 13u);
  p.state ^= (p.state >> 17u);
  p.state ^= (p.state << 5u);  
  return(p.state);
}


float myrand_uniform_0_1(inout noise_prng p)
{
    return(((float) myrand(p)) / ((float) 4294967295u));
}

float myrand_gaussian_0_1(inout noise_prng p)
{
    /* Box-Muller method for generating standard Gaussian variate */
    float u = myrand_uniform_0_1(p);
    float v = myrand_uniform_0_1(p);
    return( sqrt(-2.0 * log(u)) * cos(2.0 * pi * v) );
}

uint my_rand_poisson(inout noise_prng p, const in float lmbd)
{
      /* This is a crude approximation of Poisson distribution that should be used only for large lmbd */
      return (uint(floor(lmbd + 0.5 + (sqrt(lmbd) * myrand_gaussian_0_1(p)))));     
}

in vec3 x_tex;

out vec4 gl_FragColor;

uniform sampler2D textonTexture;
uniform float meannbofimpacts;
uniform vec3 noiseMean;
uniform uint offset;
uniform vec3 sumTexton;
uniform uint noiseOff;
uniform float scale;

void main()
{
  if(noiseOff==1u)
  {
    gl_FragColor = vec4(noiseMean, 1.0f);
    return;
  }

  vec2 textonTextureSize = textureSize(textonTexture,0);
  float r = textonTextureSize.x;
  
  float nsubdiv = 1.6f;

  vec2 xoverr =  scale/r*x_tex.xy;
  /* provide gradient for texture fetching */
  vec2 dx = scale/r*dFdx(x_tex.xy);
  vec2 dy = scale/r*dFdy(x_tex.xy);
  
  vec2 vx = vec2(0.5f,0.5f);
  ivec2 mink = ivec2(floor(nsubdiv*(xoverr - vx)));
  ivec2 maxk = ivec2(floor(nsubdiv*(xoverr + vx)));
  
  /* Initialize fragment color with (0,0,0,0) */
  gl_FragColor = vec4(0.0f,0.0f,0.0f,0.0f);
  
  /* Simulate Poisson process on the 4 neighborhood cells of (x,y) */
  for(int ncx = mink.x; ncx <= maxk.x; ncx++) /* x-cell number */
  {
    for(int ncy = mink.y; ncy <= maxk.y; ncy++) /* y-cell number */
    {
      /* seed cell = (x/w, y/w) */
      uint seed = cellseed(ivec2(ncx,ncy), offset);
      noise_prng p;
      mysrand(p, seed);
      /* Draw number of point in the cell */
      uint Ncell = my_rand_poisson(p, meannbofimpacts/(nsubdiv*nsubdiv));
      /* Draw the points of the cell */
      for(uint i=0u; i<Ncell; i++)
      {
        vec2 uv = xoverr + vec2(0.5f - (ncx+myrand_uniform_0_1(p))/nsubdiv, 0.5f - (ncy+myrand_uniform_0_1(p))/nsubdiv);
        if((uv.x>=0)&&(uv.x<=1)&&(uv.y>=0)&&(uv.y<=1))
        {
          gl_FragColor += texture2DGrad(textonTexture, uv, dx, dy);
        }
      }
    }
  }
  /* normalize and add mean */
  float lambda = meannbofimpacts/(r*r);
  gl_FragColor.xyz = 1.0f/sqrt(lambda) * (gl_FragColor.xyz - lambda * sumTexton) + noiseMean;
  gl_FragColor.w = 1.0f;
  
}



// -----------------------------------------------------------------------------    
