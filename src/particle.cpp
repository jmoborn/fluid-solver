#include "particle.h"

particle::particle(vec4 pos)
{
	this->position = pos;
	this->velocity = vec4(0,0,0);
}