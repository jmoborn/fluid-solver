#include "simulation.h"
// #include <GL/gl.h>
// #include <GL/glu.h>
// #include <GL/glx.h>
// #include <GL/glut.h>
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
        s << "volume_" << i+1 << ".vdb";
        std::string filename(s.str());
        sim.write_density(filename);

        while(t_cur < (t+sim.epsilon))
        {
            sim.calculate_timestep();
            sim.advect_fields();
            sim.enforce_boundaries();
            sim.external_forces();

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

// for(int i=0; i<sim.particles.size(); i++)
// {
//     vec4 myvel = sim.particles[i].position;
//     if(myvel.length()!=0)
//         myvel.normalize();
//     sim.particles[i].velocity = myvel;
// }
    for(int i=0; i<frames; i++)
    {
        t = (i+1)*sim.default_timestep;
        std::ostringstream s;
        s << "particles_" << i+1 << ".txt";
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
    // omp_set_num_threads(8);
    // grid sim(100, 100, 100, vec4(3.0,3.0,3.0));
    grid sim(50, 50, 50, vec4(3.0,3.0,3.0));
    std::cout << sim.resolution << " voxels" << std::endl;
    std::cout << "cell width: " << sim.cell_width << std::endl;
    std::cout << "half cell width: " << sim.half_cell_width << std::endl;
    // sim.sources.push_back(source(vec4(0.0, 0.0, 0.0), 0.3));
    sim.sources.push_back(source(vec4(-0.85, -0.4, 0.0), 0.28)); // best for particle fluid
    // sim.sources.push_back(source(vec4(-0.85, -0.75, 0.0), 0.3)); // best for smoke
    // sim.sources.push_back(source(vec4(-0.764, -0.73, 0.16), 0.15));
    // sim.sources.push_back(source(vec4(-0.85, -0.604, -0.118), 0.17));
    // sim.sources.push_back(source(vec4(-0.74, -0.783, -0.204), 0.18));
    // sim.sources.push_back(source(vec4(0.0, -0.75, 0.0), 0.5));
    // sim.sources.push_back(source(vec4(1.5,1.5,1.5), 0.5));

    int frames = 72;
    // sim.add_sources();
    // double t = 0.0;
    // double t_cur = 0.0;

    // vec4 epicenter(0.0,-0.75,0.0);
    // for(int i=0; i<sim.dimx; i++)
    // {
    //     for(int j=0; j<sim.dimy; j++)
    //     {
    //         for(int k=0; k<sim.dimz; k++)
    //         {
    //             vec4 pos = sim.position_from_index(vec4(i,j,k)) - epicenter;
    //             if(pos.length()!=0)
    //             {
    //                 pos.normalize();
    //                 pos *= 0.5;
    //                 sim.cells[i][j][k].velocity = pos;
    //                 // sim.cells[i][j][k].velocity = vec4(1,0,1);
    //             }       
    //         }
    //     }
    // }

    // for(int i=0; i<sim.dimx; i++)
    // {
    //     for(int j=0; j<sim.dimy; j++)
    //     {
    //         for(int k=0; k<sim.dimz; k++)
    //         {
    //             sim.cells[i][j][k].velocity = vec4(1,0,0);    
    //         }
    //     }
    // }

    // vec4 left_corner = sim.position_from_index(vec4(0,0,0));
    // vec4 right_corner = sim.position_from_index(vec4(sim.dimx-1, sim.dimy-1, sim.dimz-1));
    // std::cout << "lower corner: " << left_corner << std::endl;
    // std::cout << "upper corner: " << right_corner << std::endl;
    // std::cout << "lower_bound: " << sim.lower_bound << std::endl;
    // std::cout << "upper_bound: " << sim.upper_bound << std::endl;

    // TRANSFER TEST
    // sim.add_sources(true);
    // for(int i=0; i<sim.particles.size(); i++)
    // {
    //     sim.particles[i].velocity = vec4(1,1,1);
    // }
    // for(int t=0; t<5; t++)
    // {
    // sim.transfer_to_grid();
    // for(int i=0; i<sim.dimx; i++)
    // {
    //     for(int j=0; j<sim.dimy; j++)
    //     {
    //         for(int k=0; k<sim.dimz; k++)
    //         {
    //             sim.cells[i][j][k].velocity -= vec4(0.0,0.1,0.1);    
    //         }
    //     }
    // }
    // sim.transfer_to_particles();
    // }

    // vec4 xvel = vec4(0.9,0,0);
    // for(int i=0; i<sim.particles.size(); i++)
    // {
    //     if(!cmp_float(sim.particles[i].velocity.z,0.5) || !cmp_float(sim.particles[i].velocity.y,0.5))
    //     {
    //         std::cout << "incorrect: " << sim.particles[i].velocity.x << " " << sim.particles[i].velocity.y << " " << sim.particles[i].velocity.z << std::endl;
    //     }
    // }

    double start = when();
    // solve_smoke(sim, frames);
    solve_flip(sim, frames);
    //TODO: staggered MAC grid (fix interpolation functions)
    //TODO: add pressure envelope around fluid
    //TODO: fix vorticity confinement
    for(int i=0; i<frames; i++)
    {
        // t_cur = solve_smoke(sim, i, t_cur);
        // t = (i+1)*sim.default_timestep;
        // // t_cur = i*sim.default_timestep;
        // std::ostringstream s;
        // s << "volume_" << i+1 << ".vdb";
        // std::string filename(s.str());
        // sim.write_density(filename);
        // int frame = i%7==0 ? 1 : 0;
        
        // while(t_cur < (t+sim.epsilon))
        // {

        //     // std::cout << "============" << sim.cells[17][17][17].velocity.x << " " << sim.cells[17][17][17].velocity.y << " " << sim.cells[17][17][17].velocity.z << std::endl;
        //     // sim.calculate_timestep();
        //     // sim.advect_fields();
        //     // sim.enforce_boundaries();
        //     // sim.external_forces();
        //     // // if(frame)
        //     // // {
        //     // //     for(int i=0; i<sim.dimx; i++)
        //     // //     {
        //     // //         for(int j=0; j<sim.dimy; j++)
        //     // //         {
        //     // //             for(int k=0; k<sim.dimz; k++)
        //     // //             {
        //     // //                 for(int s=0; s<sim.sources.size(); s++)
        //     // //                 {
        //     // //                     vec4 pos = sim.position_from_index(vec4(i, j, k)) - sim.sources[s].position;
        //     // //                     if(pos.length() < sim.sources[s].radius) sim.cells[i][j][k].velocity += vec4(1.0, 0.0, 0.0);
        //     // //                 }
        //     // //             }
        //     // //         }
        //     // //     }
        //     // // }
        //     // frame = 0;
        //     // // sim.update_grid();
        //     // sim.enforce_boundaries();
        //     // sim.solve_pressure();
        //     // sim.enforce_boundaries();
        //     // sim.add_sources();

        //     solve_smoke(sim);
        //     t_cur += sim.timestep;
        //     std::cout << "end substep" << std::endl;
        // }
        
        
        // std::cout << "============== frame " << i+1 << " ==============" << std::endl;
    }
    double end = when();
    std::cout << "TIME: " << end - start << std::endl;
}
