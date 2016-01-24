#include "simulation.h"
#include <sstream>
#include <sys/time.h>

double when()
{
    struct timeval tp;
    gettimeofday(&tp, NULL);
    return ((double) tp.tv_sec + (double) tp.tv_usec * 1e-6);
}

void checknans(grid sim)
{
    for(int i=0; i<sim.dimx; i++)
    {
        for(int j=0; j<sim.dimy; j++)
        {
            for(int k=0; k<sim.dimz; k++)
            {
                if(isnan(sim.cells[i][j][k].velocity.x) || isnan(sim.cells[i][j][k].velocity.y) || isnan(sim.cells[i][j][k].velocity.z))
                {
                    std::cout << "isnan" << std::endl;
                }
            }
        }
    }
}

bool cmp_float(double a, double b)
{
    double epsilon = 0.0000000000001;
    return (a<(b+epsilon) && a>(b-epsilon));
}

void solve_smoke(grid sim, int frames)
{
    sim.buoyancy = 3.0;
    sim.gravity = 0.0;
    double t = 0.0;
    double t_cur = 0.0; 
    sim.add_sources();
    for(int i=0; i<frames; i++)
    {
        double t = (i+1)*sim.default_timestep;
        std::ostringstream s;
        s << "out/volume_" << i+1 << ".vdb";
        std::string filename(s.str());
        sim.write_density(filename);

        while(t_cur < (t+sim.epsilon))
        {
            sim.calculate_timestep();
            sim.advect_fields();
            sim.enforce_boundaries();
            sim.external_forces();
            sim.buoyancy_forces();
            sim.vorticity_confinement();

            sim.enforce_boundaries();
            sim.solve_pressure();
            sim.enforce_boundaries();
            sim.add_sources();
            t_cur += sim.timestep;
            std::cout << "end substep" << std::endl;
        }
        std::cout << "============== frame " << i+1 << " ==============" << std::endl;
    }
}

void solve_flip(grid sim, int frames)
{
    sim.buoyancy = 0.0;
    sim.gravity = 3.0;
    double t = 0.0;
    double t_cur = 0.0;
    sim.add_sources(true);

    for(int i=0; i<frames; i++)
    {
        t = (i+1)*sim.default_timestep;
        std::ostringstream s;
        s << "out/particles_" << i+1 << ".txt";
        // s << "volume_" << i+1 << ".vdb";
        std::string filename(s.str());
        sim.write_particles(filename);
        // sim.write_density(filename);
        // if(i==9)
        // {
        //     std::cout << "TEST************************" << std::endl;
        //     std::cout << sim.particles[9].velocity.x << " " << sim.particles[9].velocity.y << " " << sim.particles[9].velocity.x << std::endl;
        //     sim.transfer_to_grid();
        //     sim.transfer_to_particles();
        //     std::cout << sim.particles[9].velocity.x << " " << sim.particles[9].velocity.y << " " << sim.particles[9].velocity.x << std::endl;
        //     break;
        // }
        while(t_cur < (t+sim.epsilon))
        {
            // checknans(sim);
            
            // checknans(sim);
            // std::cout << "advect_fields" << std::endl;
            // sim.advect_fields();
            // checknans(sim);
            // std::cout << "transfer_to_particles" << std::endl;
            // sim.transfer_to_particles();
            // checknans(sim);
            std::cout << "calculate_timestep" << std::endl;
            sim.calculate_timestep();
            std::cout << "advect_particles" << std::endl;
            sim.advect_particles();
            // checknans(sim);
            std::cout << "enforce_boundaries particles" << std::endl;
            sim.enforce_boundaries_particles();
            // sim.enforce_boundaries();
            std::cout << "transfer_to_grid" << std::endl;
            sim.transfer_to_grid();
            // checknans(sim);

            // std::cout << "update_grid" << std::endl;
            // sim.update_grid();

            std::cout << "external_forces" << std::endl;
            sim.external_forces();

            // checknans(sim);
            std::cout << "enforce_boundaries" << std::endl;
            sim.enforce_boundaries();
            std::cout << "extrapolate velocity" << std::endl;
            sim.extrapolate_velocity(2);
            std::cout << "solve_pressure" << std::endl;
            sim.solve_pressure();
            // checknans(sim);
            // std::cout << "enforce_boundaries" << std::endl;
            // sim.enforce_boundaries();
            std::cout << "extrapolate velocity" << std::endl;
            sim.extrapolate_velocity(2);
            std::cout << "transfer_to_particles" << std::endl;
            sim.transfer_to_particles();
            sim.add_sources(true);
            t_cur += sim.timestep;
            std::cout << "end substep" << std::endl;
        }
        std::cout << "============== frame " << i+1 << " ==============" << std::endl;
    }
}

/*
export OMP_NUM_THREADS=8
to compile: g++ *.cpp -O3 -o fluid_solver -L /usr/lib64/atlas -lcblas -fopenmp
(openvdb)
-I /opt/hfs14.0.291/toolkit/include/ -L/opt/hfs14.0.291/dsolib/ -lopenvdb_sesi -lHalf -ltbb
*/
int main(int argc, char* argv[])
{
    // grid sim(100, 100, 100, vec4(3.0,3.0,3.0));
    grid sim(50, 50, 50, vec4(3.0,3.0,3.0));
    std::cout << sim.resolution << " voxels" << std::endl;
    std::cout << "cell width: " << sim.cell_width << std::endl;
    std::cout << "half cell width: " << sim.half_cell_width << std::endl;
    // sim.sources.push_back(source(vec4(0.0, 0.0, 0.0), 0.3));
    sim.sources.push_back(source(vec4(-0.85, -0.2, 0.0), 0.28)); // best for particle fluid
    // sim.sources.push_back(source(vec4(-0.85, -0.75, 0.0), 0.3)); // best for smoke
    // sim.sources.push_back(source(vec4(-0.764, -0.73, 0.16), 0.15));
    // sim.sources.push_back(source(vec4(-0.85, -0.604, -0.118), 0.17));
    // sim.sources.push_back(source(vec4(-0.74, -0.783, -0.204), 0.18));
    // sim.sources.push_back(source(vec4(0.0, -0.75, 0.0), 0.5));
    // sim.sources.push_back(source(vec4(1.5,1.5,1.5), 0.5));

    int frames = 72;

    double start = when();
    // solve_smoke(sim, frames);
    solve_flip(sim, frames);
    double end = when();
    std::cout << "TIME: " << end - start << std::endl;
}
