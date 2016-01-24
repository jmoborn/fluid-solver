#ifndef __SIMULATION_H__
#define __SIMULATION_H__

#include "vec4.h"
#include "particle.h"

#include "pcg_solver.h"
#include "sparse_matrix.h"


#include <math.h>
#include <stdlib.h>
#include <algorithm>
#include <vector>
#include <string>
#include <fstream>
#include <omp.h>

// #include <openvdb/openvdb.h>

class cell
{
    public:
    cell();
    ~cell(){};
    vec4 velocity;
    vec4 old_velocity;
    vec4 velocity_diff;
    vec4 curl;
    vec4 weights;
    vec4 vector_fields[5];

    double mag_curl;
    double scalar_fields[4];

    int solid;
    int envelope;
    
    static const int DENSITY;
    static const int PRESSURE;
    static const int TEMPERATURE;
    static const int MAG_CURL;

    static const int VELOCITY;
    static const int OLD_VELOCITY;
    static const int DELTA_VELOCITY;
    static const int WEIGHTS;
    static const int CURL;

    double get_scalar_field(int field);
    void set_scalar_field(int field, double value);
    vec4 get_vector_field(int field);
    void set_vector_field(int field, vec4 value);
};

class source
{
    public:
    source(vec4 pos, double rad);
    ~source(){};
    vec4 position;
    double radius;
    double scale;
};

class grid
{
    public:
    grid(int x, int y, int z, vec4 size);
    ~grid(){};
    vec4 position_from_index(vec4 index);
    vec4 index_from_position(vec4 position);
    cell get_cell(int i, int j, int k, int &flag);
    cell* get_cell_ref(int i, int j, int k, int &flag);
    void get_neighbors(int, int, int, cell**);
    void transfer_to_grid();
    void transfer_to_particles();
    void calculate_timestep();
    void advance_timestep();
    void advect_fields();
    void advect_particles();
    void add_sources(bool particle_source=false);
    double blur_temperature(int i, int j, int k, int filter_width);
    void external_forces();
    void vorticity_confinement();
    void buoyancy_forces();
    void extrapolate_velocity(int);
    void update_grid();
    void solve_pressure();
    void enforce_boundaries();
    void enforce_boundaries_particles();
    vec4 multiply_weights_vector(vec4 index, double* weights, int field);
    double multiply_weights_scalar(vec4 cell, double* weights, int field);
    void get_interpolation_weights(double x, double y, double z, double* weights);
    vec4 get_interpolated_vector(vec4 position, int field);
    double get_interpolated_scalar(vec4 position, int field);
    void write_particles(std::string filename);
    void write_density(std::string filename);
    int dimx, dimy, dimz; // dimensions
    vec4 size;
    vec4 lower_bound;
    vec4 upper_bound;
    double cell_width;
    double half_cell_width;
    int resolution;
    double default_timestep;
    double timestep;
    std::vector<source> sources;
    std::vector<particle> particles;
    cell*** cells;
    cell*** next_cells;

    // simulation parameters
    double dissipation;
    double buoyancy;
    double ambient_temperature;
    double gravity;
    double epsilon;
    double flip;

    cell *empty;

};

#endif
