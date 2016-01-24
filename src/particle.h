#ifndef __PARTICLE_H__
#define __PARTICLE_H__

#include "vec4.h"

class particle
{
    public:
    vec4 position;
    vec4 velocity;
    double wx[8];
    double wy[8];
    double wz[8];
    particle(vec4 pos);
    ~particle(){};
};

#endif
