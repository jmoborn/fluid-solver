#include "simulation.h"

cell::cell()
{
    velocity = vec4();
    curl = vec4();
    mag_curl = 0.0;
    temperature = 0.0;
    density = 0.0;
    pressure = 0.0;
    solid = 0;
    envelope = 0;
    weights = vec4(0,0,0);
}

source::source(vec4 pos, double rad)
{
    position = pos;
    radius = rad;
    scale = 1.6667;
}

grid::grid(int x, int y, int z, vec4 size)
{
    this->size = size;
    this->dimx = x;
    this->dimy = y;
    this->dimz = z;
    this->cell_width = size.x/x;
    this->half_cell_width = this->cell_width/2.0;
    this->upper_bound = size*0.5;
    this->lower_bound = this->upper_bound*-1;
    this->resolution = x*y*z;
    this->default_timestep = 1.0/24.0;
    this->timestep = this->default_timestep;
    if(cell_width != size.y/y || cell_width != size.z/z)
    {
        std::cout << "WARNING: cells are not square" << std::endl;
        size.y = cell_width*dimy;
        size.z = cell_width*dimz;
    }

    this->cells = new cell**[dimx+1];
    this->next_cells = new cell**[dimx+1];
    for(int i=0; i<dimx; i++)
    {
        cells[i] = new cell*[dimy+1];
        next_cells[i] = new cell*[dimy+1];
        for(int j=0; j<dimy; j++)
        {
            cells[i][j] = new cell[dimz+1];
            next_cells[i][j] = new cell[dimz+1];
            for(int k=0; k<dimz; k++)
            {
                cells[i][j][k] = cell();
                next_cells[i][j][k] = cell();
            }
        }
    }

    this->dissipation = 0.0;
    // this->buoyancy = 3.0;
    this->ambient_temperature = 0.1;
    this->epsilon = 0.0001;
    this->flip = 0.97;
    // this->flip = 0.0;

    empty = new cell();

    srand(resolution);
}

/*
returns the cell at the given index. If the index is outside of the
simulation grid, an empty cell is returned and the flag parm is set
to zero.
*/
cell grid::get_cell(int i, int j, int k, int &flag)
{
    if(i<0 || i>=dimx || j<0 || j>=dimy || k<0 || k>= dimz)
    {
        flag = 0;
        cell tmp;
        tmp.solid = 1;
        return tmp;
    }
    else
    {
        return cells[i][j][k];
    }
}

cell* grid::get_cell_ref(int i, int j, int k, int &flag)
{
    if(i<0 || i>=dimx || j<0 || j>=dimy || k<0 || k>= dimz)
    {
        flag = 0;
        return empty;
    }
    else
    {
        return &cells[i][j][k];
    }
}

void grid::get_neighbors(int x, int y, int z, cell** neighbor_list)
{
    int flag = 1;
    neighbor_list[0] = get_cell_ref(x-1,y,z,flag);
    neighbor_list[1] = get_cell_ref(x,y-1,z,flag);
    neighbor_list[2] = get_cell_ref(x,y,z-1,flag);
    neighbor_list[3] = get_cell_ref(x+1,y,z,flag);
    neighbor_list[4] = get_cell_ref(x,y+1,z,flag);
    neighbor_list[5] = get_cell_ref(x,y,z+1,flag);
}

void grid::transfer_to_grid()
{
    for(int i=0; i<dimx; i++)
    {
        for(int j=0; j<dimy; j++)
        {
            for(int k=0; k<dimz; k++)
            {
                cells[i][j][k].density = 0;
                cells[i][j][k].velocity = vec4(0,0,0);
                cells[i][j][k].old_velocity = vec4(0,0,0);
                cells[i][j][k].weights = vec4(0,0,0);

                cells[i][j][k].pressure = 0;
                cells[i][j][k].envelope = -1;
            }
        }
    }

    for(int i=0; i<particles.size(); i++)
    {
        // pos_norm.x = (position.x+size.x/2)/cell_width;
        // pos_norm.y = (position.y+size.y/2)/cell_width;
        // pos_norm.z = (position.z+size.z/2)/cell_width;

        // std::cout << "for each particle: " << i << std::endl;
        vec4 idx = index_from_position(particles[i].position);
        int x = idx.x; int y = idx.y; int z = idx.z;
        vec4 pos_norm = (particles[i].position + size*0.5)*(1/cell_width);
        double wx[8];
        vec4 pos_x(pos_norm); pos_x.y -= 0.5; pos_x.z -= 0.5; //pos_x.x -= 0.5;//
        get_interpolation_weights(pos_x.x, pos_x.y, pos_x.z, wx);
        double wy[8];
        vec4 pos_y(pos_norm); pos_y.x -= 0.5; pos_y.z -= 0.5; //pos_y.y -= 0.5;//
        get_interpolation_weights(pos_y.x, pos_y.y, pos_y.z, wy);
        double wz[8];
        vec4 pos_z(pos_norm); pos_z.x -= 0.5; pos_z.y -= 0.5; //pos_z.z -= 0.5;//
        get_interpolation_weights(pos_z.x, pos_z.y, pos_z.z, wz);
        // std::cout << "set density: " << x << ", " << y << ", " << z << std::endl;
        // std::cout << "position" << particles[i].position.x << ", " << particles[i].position.y << ", " << particles[i].position.z << std::endl;
        cells[x][y][z].density += 1;
// std::cout << "done setting density" << std::endl;
        int xx =floor(pos_x.x); int xy =floor(pos_x.y); int xz =floor(pos_x.z);
        int yx =floor(pos_y.x); int yy =floor(pos_y.y); int yz =floor(pos_y.z);
        int zx =floor(pos_z.x); int zy =floor(pos_z.y); int zz =floor(pos_z.z);

        // std::cout << "fuck this" << std::endl;
        // std::cout << "x: " << xx << " " << xy << " " << xz << std::endl;
        // std::cout << "y: " << yx << " " << yy << " " << yz << std::endl;
        // std::cout << "z: " << zx << " " << zy << " " << zz << std::endl;
        // cells[xx][xy][xz].velocity.x += particles[i].velocity.x*wx[0];
        // for(int q=0; q<8; q++)
        // {
        //     std::cout << wx[q] << " ";
        // }
        // std::cout << std::endl;

        // cells[x][y][z].velocity += particles[i].velocity;
        // cells[x][y][z].weights += vec4(1,1,1);

        int flag = 1;
        get_cell_ref(xx,xy,xz,flag)->velocity.x += particles[i].velocity.x*wx[0];
        // cells[xx][xy][xz].weights.x += wx[0];
        get_cell_ref(xx,xy,xz,flag)->weights.x += wx[0];
        // cells[yx][yy][yz].velocity.y += particles[i].velocity.y*wy[0];
        get_cell_ref(yx,yy,yz,flag)->velocity.y += particles[i].velocity.y*wy[0];
        // cells[yx][yy][yz].weights.y += wy[0];
        get_cell_ref(yx,yy,yz,flag)->weights.y += wy[0];
        // cells[zx][zy][zz].velocity.z += particles[i].velocity.z*wz[0];
        get_cell_ref(zx,zy,zz,flag)->velocity.z += particles[i].velocity.z*wz[0];
        // cells[zx][zy][zz].weights.z += wz[0];
        get_cell_ref(zx,zy,zz,flag)->weights.z += wz[0];

        // cells[xx+1][xy][xz].velocity.x += particles[i].velocity.x*wx[1];
        get_cell_ref(xx+1,xy,xz,flag)->velocity.x += particles[i].velocity.x*wx[1];
        // cells[xx+1][xy][xz].weights.x += wx[1];
        get_cell_ref(xx+1,xy,xz,flag)->weights.x += wx[1];
        // cells[yx+1][yy][yz].velocity.y += particles[i].velocity.y*wy[1];
        get_cell_ref(yx+1,yy,yz,flag)->velocity.y += particles[i].velocity.y*wy[1];
        // cells[yx+1][yy][yz].weights.y += wy[1];
        get_cell_ref(yx+1,yy,yz,flag)->weights.y += wy[1];
        // cells[zx+1][zy][zz].velocity.z += particles[i].velocity.z*wz[1];
        get_cell_ref(zx+1,zy,zz,flag)->velocity.z += particles[i].velocity.z*wz[1];
        // cells[zx+1][zy][zz].weights.z += wz[1];
        get_cell_ref(zx+1,zy,zz,flag)->weights.z += wz[1];

        // cells[xx][xy+1][xz].velocity.x += particles[i].velocity.x*wx[2];
        get_cell_ref(xx,xy+1,xz,flag)->velocity.x += particles[i].velocity.x*wx[2];
        // cells[xx][xy+1][xz].weights.x += wx[2];
        get_cell_ref(xx,xy+1,xz,flag)->weights.x += wx[2];
        // cells[yx][yy+1][yz].velocity.y += particles[i].velocity.y*wy[2];
        get_cell_ref(yx,yy+1,yz,flag)->velocity.y += particles[i].velocity.y*wy[2];
        // cells[yx][yy+1][yz].weights.y += wy[2];
        get_cell_ref(yx,yy+1,yz,flag)->weights.y += wy[2];
        // cells[zx][zy+1][zz].velocity.z += particles[i].velocity.z*wz[2];
        get_cell_ref(zx,zy+1,zz,flag)->velocity.z += particles[i].velocity.z*wz[2];
        // cells[zx][zy+1][zz].weights.z += wz[2];
        get_cell_ref(zx,zy+1,zz,flag)->weights.z += wz[2];

        // cells[xx+1][xy+1][xz].velocity.x += particles[i].velocity.x*wx[3];
        get_cell_ref(xx+1,xy+1,xz,flag)->velocity.x += particles[i].velocity.x*wx[3];
        // cells[xx+1][xy+1][xz].weights.x += wx[3];
        get_cell_ref(xx+1,xy+1,xz,flag)->weights.x += wx[3];
        // cells[yx+1][yy+1][yz].velocity.y += particles[i].velocity.y*wy[3];
        get_cell_ref(yx+1,yy+1,yz,flag)->velocity.y += particles[i].velocity.y*wy[3];
        // cells[yx+1][yy+1][yz].weights.y += wy[3];
        get_cell_ref(yx+1,yy+1,yz,flag)->weights.y += wy[3];
        // cells[zx+1][zy+1][zz].velocity.z += particles[i].velocity.z*wz[3];
        get_cell_ref(zx+1,zy+1,zz,flag)->velocity.z += particles[i].velocity.z*wz[3];
        // cells[zx+1][zy+1][zz].weights.z += wz[3];
        get_cell_ref(zx+1,zy+1,zz,flag)->weights.z += wz[3];

        // cells[xx][xy][xz+1].velocity.x += particles[i].velocity.x*wx[4];
        get_cell_ref(xx,xy,xz+1,flag)->velocity.x += particles[i].velocity.x*wx[4];
        // cells[xx][xy][xz+1].weights.x += wx[4];
        get_cell_ref(xx,xy,xz+1,flag)->weights.x += wx[4];
        // cells[yx][yy][yz+1].velocity.y += particles[i].velocity.y*wy[4];
        get_cell_ref(yx,yy,yz+1,flag)->velocity.y += particles[i].velocity.y*wy[4];
        // cells[yx][yy][yz+1].weights.y += wy[4];
        get_cell_ref(yx,yy,yz+1,flag)->weights.y += wy[4];
        // cells[zx][zy][zz+1].velocity.z += particles[i].velocity.z*wz[4];
        get_cell_ref(zx,zy,zz+1,flag)->velocity.z += particles[i].velocity.z*wz[4];
        // cells[zx][zy][zz+1].weights.z += wz[4];
        get_cell_ref(zx,zy,zz+1,flag)->weights.z += wz[4];

        // cells[xx+1][xy][xz+1].velocity.x += particles[i].velocity.x*wx[5];
        get_cell_ref(xx+1,xy,xz+1,flag)->velocity.x += particles[i].velocity.x*wx[5];
        // cells[xx+1][xy][xz+1].weights.x += wx[5];
        get_cell_ref(xx+1,xy,xz+1,flag)->weights.x += wx[5];
        // cells[yx+1][yy][yz+1].velocity.y += particles[i].velocity.y*wy[5];
        get_cell_ref(yx+1,yy,yz+1,flag)->velocity.y += particles[i].velocity.y*wy[5];
        // cells[yx+1][yy][yz+1].weights.y += wy[5];
        get_cell_ref(yx+1,yy,yz+1,flag)->weights.y += wy[5];
        // cells[zx+1][zy][zz+1].velocity.z += particles[i].velocity.z*wz[5];
        get_cell_ref(zx+1,zy,zz+1,flag)->velocity.z += particles[i].velocity.z*wz[5];
        // cells[zx+1][zy][zz+1].weights.z += wz[5];
        get_cell_ref(zx+1,zy,zz+1,flag)->weights.z += wz[5];

        // cells[xx][xy+1][xz+1].velocity.x += particles[i].velocity.x*wx[6];
        get_cell_ref(xx,xy+1,xz+1,flag)->velocity.x += particles[i].velocity.x*wx[6];
        // cells[xx][xy+1][xz+1].weights.x += wx[6];
        get_cell_ref(xx,xy+1,xz+1,flag)->weights.x += wx[6];
        // cells[yx][yy+1][yz+1].velocity.y += particles[i].velocity.y*wy[6];
        get_cell_ref(yx,yy+1,yz+1,flag)->velocity.y += particles[i].velocity.y*wy[6];
        // cells[yx][yy+1][yz+1].weights.y += wy[6];
        get_cell_ref(yx,yy+1,yz+1,flag)->weights.y += wy[6];
        // cells[zx][zy+1][zz+1].velocity.z += particles[i].velocity.z*wz[6];
        get_cell_ref(zx,zy+1,zz+1,flag)->velocity.z += particles[i].velocity.z*wz[6];
        // cells[zx][zy+1][zz+1].weights.z += wz[6];
        get_cell_ref(zx,zy+1,zz+1,flag)->weights.z += wz[6];

        // cells[xx+1][xy+1][xz+1].velocity.x += particles[i].velocity.x*wx[7];
        get_cell_ref(xx+1,xy+1,xz+1,flag)->velocity.x += particles[i].velocity.x*wx[7];
        // cells[xx+1][xy+1][xz+1].weights.x += wx[7];
        get_cell_ref(xx+1,xy+1,xz+1,flag)->weights.x += wx[7];
        // cells[yx+1][yy+1][yz+1].velocity.y += particles[i].velocity.y*wy[7];
        get_cell_ref(yx+1,yy+1,yz+1,flag)->velocity.y += particles[i].velocity.y*wy[7];
        // cells[yx+1][yy+1][yz+1].weights.y += wy[7];
        get_cell_ref(yx+1,yy+1,yz+1,flag)->weights.y += wy[7];
        // cells[zx+1][zy+1][zz+1].velocity.z += particles[i].velocity.z*wz[7];
        get_cell_ref(zx+1,zy+1,zz+1,flag)->velocity.z += particles[i].velocity.z*wz[7];
        // cells[zx+1][zy+1][zz+1].weights.z += wz[7];
        get_cell_ref(zx+1,zy+1,zz+1,flag)->weights.z += wz[7];

        // particles[i].wx = wx;
        // particles[i].wy = wy;
        // particles[i].wz = wz;

        // cells[(int)idx.x][(int)idx.y][(int)idx.z].geometry.push_back(&(particles[i]));
    }

    for(int i=0; i<dimx; i++)
    {
        for(int j=0; j<dimy; j++)
        {
            for(int k=0; k<dimz; k++)
            {
                // if(cells[i][j][k].density>0)
                if(cells[i][j][k].weights.length()>0)
                {
                    // std::cout << "particles contributing: " << cells[i][j][k].density << std::endl;
                    // std::cout << cells[i][j][k].velocity.x << " / " << cells[i][j][k].weights.x << std::endl;
                    // std::cout << cells[i][j][k].velocity.y << " / " << cells[i][j][k].weights.y << std::endl;
                    // std::cout << cells[i][j][k].velocity.z << " / " << cells[i][j][k].weights.z << std::endl;

                    cells[i][j][k].weights.x>0 ? cells[i][j][k].velocity.x /= cells[i][j][k].weights.x : cells[i][j][k].velocity.x = 0;
                    cells[i][j][k].weights.y>0 ? cells[i][j][k].velocity.y /= cells[i][j][k].weights.y : cells[i][j][k].velocity.y = 0;
                    cells[i][j][k].weights.z>0 ? cells[i][j][k].velocity.z /= cells[i][j][k].weights.z : cells[i][j][k].velocity.z = 0;
                    cells[i][j][k].old_velocity = cells[i][j][k].velocity;
                    // cells[i][j][k].density = 1;
                    // cells[i][j][k].updated = 1;
                    cells[i][j][k].envelope = 0;

                    // if(cells[i][j][k].old_velocity.x != cells[i][j][k].velocity.x) std::cout << "WTF???" << std::endl;
                    
                    // if(isnan(cells[i][j][k].velocity.x) || isnan(cells[i][j][k].velocity.y) || isnan(cells[i][j][k].velocity.z))
                    // {
                    //     std::cout << "grid transfer nan: " << cells[i][j][k].velocity.x << ", " << cells[i][j][k].velocity.y << ", " << cells[i][j][k].velocity.z << std::endl;
                    // }
                }
            }
        }
    }
}

void grid::transfer_to_particles()
{
    for(int i=0; i<dimx; i++)
    {
        for(int j=0; j<dimy; j++)
        {
            for(int k=0; k<dimz; k++)
            {
                cells[i][j][k].velocity_diff = cells[i][j][k].velocity - cells[i][j][k].old_velocity;
                // if(cells[i][j][k].velocity_diff.x!=0)
                // if(cells[i][j][k].velocity.x!=cells[i][j][k].old_velocity.x)
                // {
                //     // std::cout << cells[i][j][k].velocity_diff.x << " " << cells[i][j][k].velocity_diff.y << " " << cells[i][j][k].velocity_diff.z << std::endl;
                //     std::cout << "new " << cells[i][j][k].velocity.x << " " << cells[i][j][k].velocity.y << " " << cells[i][j][k].velocity.z << std::endl;
                //     std::cout << "old " << cells[i][j][k].old_velocity.x << " " << cells[i][j][k].old_velocity.y << " " << cells[i][j][k].old_velocity.z << std::endl;
                //     std::cout << "density: " << cells[i][j][k].density << std::endl;
                // }
                if(isnan(cells[i][j][k].velocity_diff.x) || isnan(cells[i][j][k].velocity_diff.y) || isnan(cells[i][j][k].velocity_diff.z))
                {
                    std::cout << "nan vel_diff due to: " << cells[i][j][k].old_velocity.x << ", " << cells[i][j][k].old_velocity.y << ", " << cells[i][j][k].old_velocity.z << std::endl;
                    std::cout << "    " << cells[i][j][k].velocity.x << ", " << cells[i][j][k].velocity.y << ", " << cells[i][j][k].velocity.z << std::endl;
                }
            }
        }
    }

    double pic = 1.0 - flip;
    for(int i=0; i<particles.size(); i++)
    {
        // vec4 p_x(particles[i].position); p_x.x -= half_cell_width;
        // particles[i].velocity.x += get_velocity_difference(p_x).x;
        // vec4 p_y(particles[i].position); p_y.y -= half_cell_width;
        // particles[i].velocity.y += get_velocity_difference(p_y).y;
        // vec4 p_z(particles[i].position); p_z.z -= half_cell_width;
        // particles[i].velocity.z += get_velocity_difference(p_z).z;
        // vec4 vel_diff = get_velocity_difference(particles[i].position);
        // std::cout << "vel_diff: " << vel_diff.x << " " << vel_diff.y << " " << vel_diff.z << std::endl;
        particles[i].velocity += get_velocity_difference(particles[i].position);
        particles[i].velocity *= flip;
        particles[i].velocity += get_velocity(particles[i].position)*pic;
    }
}

void grid::calculate_timestep()
{
    double vmax = 0.0;
    double vcur = 0.0;
    int ix, iy, iz;
    for(int i=0; i<dimx; i++)
    {
        for(int j=0; j<dimy; j++)
        {
            for(int k=0; k<dimz; k++)
            {
                vcur = cells[i][j][k].velocity.length();
                if(vcur > vmax && (cells[i][j][k].density > epsilon || cells[i][j][k].envelope==0)) vmax = vcur;
                // if(vcur > vmax)
                // {
                //     vmax = vcur;  
                //     ix = i; iy = j; iz = k;
                // }
            }
        }
    }
    double k_cfl = 2.0;
    if(vmax==0) timestep = default_timestep;
    timestep = cell_width*k_cfl/vmax;
    if(timestep>default_timestep) timestep = default_timestep;
    if(timestep<epsilon) timestep = epsilon;
    std::cout << "max_vel: " << vmax << std::endl;
    // std::cout << "at: " << ix << ", " << iy << ", " << iz << std::endl;
    std::cout << "timestep: " << timestep << std::endl;
}

void grid::advect_particles()
{
    for(int i=0; i<particles.size(); i++)
    {
        particles[i].position += particles[i].velocity*timestep; // TODO: RK3
        // if(isnan(particles[i].position.x) || isnan(particles[i].position.y) || isnan(particles[i].position.z))
        // {
        //     std::cout << "nan position due to vel: " << particles[i].velocity.x << ", " << particles[i].velocity.y << ", " << particles[i].velocity.z << std::endl;
        // }
    }
}

void grid::advect_fields()
{
    #pragma omp parallel for schedule(dynamic, 16)
    for(int i=0; i<dimx; i++)
    {
        for(int j=0; j<dimy; j++)
        {
            for(int k=0; k<dimz; k++)
            {
                // Backwards Particle Trace
                // vec4 v = this->cells[i][j][k].velocity;
// std::cout << v.x << ", " << v.y << ", " << v.z << std::endl;
                vec4 p = position_from_index(vec4(i,j,k));
                vec4 v = get_velocity(p);
                vec4 p_x(p); p_x.x -= half_cell_width;
                vec4 p_y(p); p_y.y -= half_cell_width;
                vec4 p_z(p); p_z.z -= half_cell_width;
                vec4 v_x = get_velocity(p_x);
                vec4 v_y = get_velocity(p_y);
                vec4 v_z = get_velocity(p_z);
                // if(isnan(p.x) || isnan(p.y) || isnan(p.z)) std::cout << i << " " << j << " " << k << std::endl;
// double w = 1;
// if((i==24 && j==15 && k==24) || (i==35 && j==15 && k==35)) w = 0.999;
// vec4 p_half(p.x-0.5*timestep*v.x, p.y-0.5*timestep*v.y, p.z-0.5*timestep*v.z);
// p_half.w = w;
// v = get_velocity(p_half);
                v = get_velocity(vec4(p.x-0.5*timestep*v.x, p.y-0.5*timestep*v.y, p.z-0.5*timestep*v.z));
                v_x = get_velocity(vec4(p_x.x-0.5*timestep*v_x.x, p_x.y-0.5*timestep*v_x.y, p_x.z-0.5*timestep*v_x.z));
                v_y = get_velocity(vec4(p_y.x-0.5*timestep*v_y.x, p_y.y-0.5*timestep*v_y.y, p_y.z-0.5*timestep*v_y.z));
                v_z = get_velocity(vec4(p_z.x-0.5*timestep*v_z.x, p_z.y-0.5*timestep*v_z.y, p_z.z-0.5*timestep*v_z.z));
                if(isnan(v.x) || isnan(v.y) || isnan(v.z)) std::cout << i << " " << j << " " << k << std::endl;
// std::cout << v.x << ", " << v.y << ", " << v.z << std::endl;
                p = p - v*timestep;
                p_x = p_x - v_x*timestep;
                p_y = p_y - v_y*timestep;
                p_z = p_z - v_z*timestep;
                // vec4 p_x(p); p_x.x -= half_cell_width;
                // vec4 p_y(p); p_y.y -= half_cell_width;
                // vec4 p_z(p); p_z.z -= half_cell_width;
                // this->next_cells[i][j][k].velocity = get_velocity(p);
                


                this->next_cells[i][j][k].velocity.x = get_velocity(p_x).x;
                this->next_cells[i][j][k].velocity.y = get_velocity(p_y).y;
                this->next_cells[i][j][k].velocity.z = get_velocity(p_z).z;
                // if(this->next_cells[i][j][k].velocity.length() > 1.0) std::cout << this->next_cells[i][j][k].velocity.x << " " << this->next_cells[i][j][k].velocity.y << " " << this->next_cells[i][j][k].velocity.z << std::endl;
// std::cout << i << " " << j << " " << k << std::endl;
// if(i==32 && j==32 && k==32)p.w = 0.999;
                this->next_cells[i][j][k].density = get_density(p);
                this->next_cells[i][j][k].temperature = get_temperature(p);
                this->next_cells[i][j][k].pressure = 0;
// if(this->next_cells[i][j][k].density < 0) std::cout << "ERROR: NEGATIVE DENSITY" << std::endl;
            }
        }
    }
    advance_timestep();
}

void grid::add_sources(bool particle_source)
{
    double fourth_cell_width = half_cell_width/2.0;
    
    for(int i=0; i<dimx; i++)
    {
        for(int j=0; j<dimy; j++)
        {
            for(int k=0; k<dimz; k++)
            {
                for(int s=0; s<sources.size(); s++)
                {
                    vec4 gpos = position_from_index(vec4(i, j, k));
                    vec4 pos = gpos - sources[s].position;
                    if (pos.length()<sources[s].radius)
                    {
                        if(particle_source)
                        {
                            vec4 fpos = gpos - half_cell_width;
                            vec4 cpos = gpos + half_cell_width;
                            double xp, yp, zp;
                            double rx, ry, rz;
                            for(int x=0; x<2; x++)
                            {
                                for(int y=0; y<2; y++)
                                {
                                    for(int z=0; z<2; z++)
                                    {
                                        // TODO: perfectly round source
                                        rx = rand() / (RAND_MAX + 1.)*half_cell_width - fourth_cell_width;
                                        ry = rand() / (RAND_MAX + 1.)*half_cell_width - fourth_cell_width;
                                        rz = rand() / (RAND_MAX + 1.)*half_cell_width - fourth_cell_width;
                                        x==0 ? xp = fpos.x + fourth_cell_width : xp = cpos.x - fourth_cell_width;
                                        y==0 ? yp = fpos.y + fourth_cell_width : yp = cpos.y - fourth_cell_width;
                                        z==0 ? zp = fpos.z + fourth_cell_width : zp = cpos.z - fourth_cell_width;
                                        particles.push_back(particle(vec4(xp,yp,zp)+vec4(rx,ry,rz)));
                                    }
                                }
                            }
                        }
                        else
                        {
                            cells[i][j][k].density += sources[s].scale*timestep;
                            cells[i][j][k].temperature += sources[s].scale*timestep;
                        }
                    }
                }
            }
        }
    }       
    
}

double grid::blur_temperature(int i, int j, int k, int filter_width)
{
    // temperature diffusion (the container absorbs heat . . . because i'm too lazy for boundary checking)
    // int filter_width = radius/cell_width;
    if(filter_width%2==0) filter_width--;
    if(filter_width<3) filter_width = 3;
    int offset = filter_width/2;
    double blur = 0;
    for(int x=-offset; x<=offset; x++)
    {
        for(int y=-offset; y<=offset; y++)
        {
            for(int z=-offset; z<=offset; z++)
            {
                int flag = 1; // in case we do boundary checking later . . .
                blur += get_cell(i+x, j+y, k+z, flag).temperature;
            }    
        }
    }
    blur /= 27.0;
    // if(blur <0.5 && blur > 0.0001) std::cout << blur << std::endl;
    return blur;
}

void grid::external_forces()
{
    // BE CAREFUL TO SET next_cells equal to
    for(int i=0; i<dimx; i++)
    {
        for(int j=0; j<dimy; j++)
        {
            for(int k=0; k<dimz; k++)
            {
                if( cells[i][j][k].density > epsilon ||  cells[i][j][k].temperature > epsilon)
                {

                    // dissipation
                    cells[i][j][k].density = std::max(0.0, cells[i][j][k].density - dissipation*timestep);
                    // next_cells[i][j][k].temperature = std::max(0.0, cells[i][j][k].temperature - dissipation*5.0*timestep);

                    // next_cells[i][j][k].temperature = blur_temperature(i, j, k, 3);
                    // // if(cells[i][j][k].temperature <0.5 && cells[i][j][k].temperature > 0.0001) std::cout << cells[i][j][k].temperature << std::endl;
                    // // if(isnan(cells[i][j][k].temperature)) std::cout << "nan temperature" << std::endl;

                    // vec4 corner = position_from_index(vec4(i,j,k)) - vec4(0,half_cell_width,0);
                    // double temp = get_temperature(corner);
                    // // double temp = cells[i][j][k].temperature;/
                    // double k_rise = temp*buoyancy;
                    // double k_fall = 0.0;//cells[i][j][k].density*-buoyancy/4;
                    // next_cells[i][j][k].velocity = cells[i][j][k].velocity + vec4(0.0, 1.0, 0.0)*((k_rise+k_fall)*timestep);


                    // pumps
                    for(int s=0; s<sources.size(); s++)
                    {
                        vec4 epicenter = sources[s].position - vec4(sources[s].radius, 0, 0);
                        vec4 pos = position_from_index(vec4(i, j, k));
                        vec4 dir = pos - epicenter;
                        if(dir.length()==0)dir = vec4(0,0,0);
                        else dir.normalize();
                        vec4 rad = pos - sources[s].position;
                        dir = vec4(1,0,0);
                        if(rad.length() < sources[s].radius) cells[i][j][k].velocity += dir*5.0*timestep;
                    }

                }
                // wind force
                // if( cells[i][j][k].density < epsilon)
                // {
                //     // TODO: make user option
                //     double k_wind = 0.5;
                //     // double k_wind = 0.0;
                //     vec4 wind_force(1,0,0);
                //     next_cells[i][j][k].velocity = cells[i][j][k].velocity*(1 - k_wind) + wind_force*(k_wind*timestep);
                // }
                //gravity
                if(cells[i][j][k].density > epsilon || cells[i][j][k].envelope==0)
                    cells[i][j][k].velocity +=(vec4(0,-1,0)*(gravity*timestep));
                    // cells[i][j][k].velocity +=(vec4(0,-.707,.707)*(gravity*timestep));

                if(isnan(cells[i][j][k].velocity.x) || isnan(cells[i][j][k].velocity.y) || isnan(cells[i][j][k].velocity.z))
                {
                    std::cout << "external_forces nan: " << cells[i][j][k].velocity.x << ", " << cells[i][j][k].velocity.y << ", " << cells[i][j][k].velocity.z << std::endl;
                }
                // if(cells[i][j][k].density < epsilon && !cells[i][j][k].envelope && cells[i][j][k].temperature < epsilon)
                //     next_cells[i][j][k].velocity = vec4(0,0,0);
            }
        }
    }

    // advance_timestep(); //TODO: figure out best way to to this without messing stuff up
    // vorticity confinement
    // TODO: vorticity confinement causes error when fluid reaches corners ???
    // for(int i=0; i<dimx; i++)
    // {
    //     for(int j=0; j<dimy; j++)
    //     {
    //         for(int k=0; k<dimz; k++)
    //         {
    //             // calculate curl
    //             int f = 0;
    //             vec4 pos = position_from_index(vec4(i,j,k));
    //             vec4 vx0 = get_velocity(vec4(pos.x-half_cell_width,pos.y,pos.z));
    //             vec4 vx1 = get_velocity(vec4(pos.x+half_cell_width,pos.y,pos.z));
    //             vec4 vy0 = get_velocity(vec4(pos.x,pos.y-half_cell_width,pos.z));
    //             vec4 vy1 = get_velocity(vec4(pos.x,pos.y+half_cell_width,pos.z));
    //             vec4 vz0 = get_velocity(vec4(pos.x,pos.y,pos.z-half_cell_width));
    //             vec4 vz1 = get_velocity(vec4(pos.x,pos.y,pos.z+half_cell_width));
                
    //             double dx = (vy1.z - vy0.z) - (vz1.y - vz0.y);
    //             double dy = (vz1.x - vz0.x) - (vx1.z - vx0.z);
    //             double dz = (vx1.y - vx0.y) - (vy1.x - vy0.x);

    //             // double dx = get_cell(i+1,j,k,f).velocity.x - get_cell(i,j,k,f).velocity.x;
    //             // double dy = get_cell(i,j+1,k,f).velocity.y - get_cell(i,j,k,f).velocity.y;
    //             // double dz = get_cell(i,j,k+1,f).velocity.z - get_cell(i,j,k,f).velocity.z;
    //             // cells[i][j][k].curl = vec4(dz-dy, dx-dz, dy-dx);
    //             cells[i][j][k].curl = vec4(dx,dy,dz);
    //             cells[i][j][k].mag_curl = cells[i][j][k].curl.length();
    //         }
    //     }
    // }

    // for(int i=0; i<dimx; i++)
    // {
    //     for(int j=0; j<dimy; j++)
    //     {
    //         for(int k=0; k<dimz; k++)
    //         {
    //             // vorticity confinement
    //             int f = 0;
    //             vec4 grad_curl = vec4(get_cell(i+1,j,k,f).mag_curl - get_cell(i-1,j,k,f).mag_curl,
    //                                   get_cell(i,j+1,k,f).mag_curl - get_cell(i,j-1,k,f).mag_curl,
    //                                   get_cell(i,j,k+1,f).mag_curl - get_cell(i,j,k-1,f).mag_curl);
    //             if(grad_curl.length()==0) grad_curl = vec4(0,0,0);
    //             else grad_curl.normalize();
    //             // grad_curl.normalize();
    //             vec4 vort = grad_curl.cross(cells[i][j][k].curl);
    //             vort *= 4.5;
    //             if(isnan(vort.x) || isnan(vort.y) || isnan(vort.z))
    //             {
    //                 std::cout << "grad_curl: " << grad_curl.x << " " << grad_curl.y << " " << grad_curl.z << std::endl;
    //             }
    //             // if(vort.x!=0) std::cout << "vort: " << vort.x << " " << vort.y << " " << vort.z << std::endl;
    //             if(cells[i][j][k].density > epsilon || cells[i][j][k].envelope)
    //                 cells[i][j][k].velocity += vort*timestep;
    //         }
    //     }
    // }
}

void grid::extrapolate_velocity(int num_cells)
{
    // initialize envelope field
    // for(int i=0; i<dimx; i++)
    // {
    //     for(int j=0; j<dimy; j++)
    //     {
    //         for(int k=0; k<dimz; k++)
    //         {
    //             if(cells[i][j][k].density>epsilon)
    //             {
    //                 cells[i][j][k].envelope = 0;
    //             }
    //             else
    //             {
    //                 cells[i][j][k].envelope = -1;
    //             }
    //         }
    //     }
    // }

    for(int layer=1; layer<=num_cells; layer++)
    {
        for(int i=0; i<dimx; i++)
        {
            for(int j=0; j<dimy; j++)
            {
                for(int k=0; k<dimz; k++)
                {
                    if(cells[i][j][k].envelope==-1)
                    {
                        cell* neighbors[6];
                        get_neighbors(i, j, k, neighbors);
                        // cells[i][j][k].velocity = vec4(0,0,0);
                        // vec4 count(0,0,0);
                        double count = 0.0;
                        vec4 old_vel = cells[i][j][k].velocity;
                        cells[i][j][k].velocity = vec4(0,0,0);
                        for(int c=0; c<6; c++)
                        {
                            // border[c] = neighbors[c]->envelope==layer-1;
                            if(neighbors[c]->envelope==layer-1)
                            {
                                cells[i][j][k].velocity += neighbors[c]->velocity;
                                count += 1;
                            }
                        }
                        if(count)
                        {
                            cells[i][j][k].velocity *= (1.0/count);
                            cells[i][j][k].envelope = layer;
                        }
                        // don't average components at fluid boundaries
                        if(neighbors[0]->envelope==0) cells[i][j][k].velocity.x = old_vel.x;
                        if(neighbors[1]->envelope==0) cells[i][j][k].velocity.y = old_vel.y;
                        if(neighbors[2]->envelope==0) cells[i][j][k].velocity.z = old_vel.z;
                    }
                }
            }
        }
    }
}

/*
identify cells that form a pressure envelope around the fluid
*/
void grid::update_grid()
{
    for(int i=0; i<dimx; i++)
    {
        for(int j=0; j<dimy; j++)
        {
            for(int k=0; k<dimz; k++)
            {
                cells[i][j][k].envelope = 0;
            }
        }
    }

    for(int i=0; i<dimx; i++)
    {
        for(int j=0; j<dimy; j++)
        {
            for(int k=0; k<dimz; k++)
            {
                if(cells[i][j][k].density > epsilon)
                {
                    int offset = 2;
                    for(int x=-offset; x<=offset; x++)
                    {
                        for(int y=-offset; y<=offset; y++)
                        {
                            for(int z=-offset; z<=offset; z++)
                            {
                                int flag = 1; // in case we do boundary checking later . . .
                                cell neighbor = get_cell(i+x, j+y, k+z, flag);
                                if(flag && (neighbor.density <= epsilon)) 
                                {
                                    cells[i+x][j+y][k+z].envelope = 1;
                                    // std::cout << "setting envelope" << std::endl;
                                }
                            }    
                        }
                    }
                }
            }
        }
    }
}

void grid::solve_pressure()
{
std::cout << "starting pressure solve" << std::endl;
    // allocate neighbor map
    int ***neighbor_map = new int**[dimx];
    for(int i=0; i<dimx; i++)
    {
        neighbor_map[i] = new int*[dimy];
        for(int j=0; j<dimy; j++)
        {
            neighbor_map[i][j] = new int[dimz];
            for(int k=0; k<dimz; k++)
            {
                neighbor_map[i][j][k] = -1;
            }
        }
    }

    // build neighbor map
    int A_idx = 0;
    for(int i=0; i<dimx; i++)
    {
        for(int j=0; j<dimy; j++)
        {
            for(int k=0; k<dimz; k++)
            {
                // if(cells[i][j][k].density < 3) cells[i][j][k].density == 0;
                if((cells[i][j][k].density > epsilon))// || (cells[i][j][k].envelope > 0))
                {
                    neighbor_map[i][j][k] = A_idx;
                    A_idx++;
                    // if(cells[i][j][k].density < epsilon) std::cout << "envelope added to neighbor map" << std::endl;
                }
            }
        }
    }
std::cout << "fluid cell count: " << A_idx << std::endl;
    // build matrix
    SparseMatrix<double> A(A_idx);
    std::vector<double> b;
    std::vector<double> p;
    int f = 0;
    for(int i=0; i<dimx; i++)
    {
        for(int j=0; j<dimy; j++)
        {
            for(int k=0; k<dimz; k++)
            {
                if(neighbor_map[i][j][k]>=0)
                {
                    int non_solids = 6;
                    cell c = get_cell(i-1,j,k,f);
                    if(c.density > epsilon)
                        A.set_element(neighbor_map[i][j][k], neighbor_map[i-1][j][k], -1);
                    else if(c.solid)
                        non_solids--;
                    c = get_cell(i+1,j,k,f);
                    if(c.density > epsilon)
                        A.set_element(neighbor_map[i][j][k], neighbor_map[i+1][j][k], -1);
                    else if(c.solid)
                        non_solids--;
                    c = get_cell(i,j-1,k,f);
                    if(c.density > epsilon)
                        A.set_element(neighbor_map[i][j][k], neighbor_map[i][j-1][k], -1);
                    else if(c.solid)
                        non_solids--;
                    c = get_cell(i,j+1,k,f);
                    if(c.density > epsilon)
                        A.set_element(neighbor_map[i][j][k], neighbor_map[i][j+1][k], -1);
                    else if(c.solid)
                        non_solids--;
                    c = get_cell(i,j,k-1,f);
                    if(c.density > epsilon)
                        A.set_element(neighbor_map[i][j][k], neighbor_map[i][j][k-1], -1);
                    else if(c.solid)
                        non_solids--;
                    c = get_cell(i,j,k+1,f);
                    if(c.density > epsilon)
                        A.set_element(neighbor_map[i][j][k], neighbor_map[i][j][k+1], -1);
                    else if(c.solid)
                        non_solids--;
                    A.set_element(neighbor_map[i][j][k], neighbor_map[i][j][k], non_solids);

                    // if(non_solids!=6) std::cout << "unique diagonal: " << i << ", " << j << ", " << k << ": " << non_solids << std::endl;
                    
                    c = get_cell(i,j,k,f);
                    cell x1 = get_cell(i+1,j,k,f);
                    cell y1 = get_cell(i,j+1,k,f);
                    cell z1 = get_cell(i,j,k+1,f);
                    double divx = x1.solid ? (0) : (x1.velocity.x - c.velocity.x);
                    double divy = y1.solid ? (0) : (y1.velocity.y - c.velocity.y);
                    double divz = z1.solid ? (0) : (z1.velocity.z - c.velocity.z);
                    // double divx = x1.solid ? (0) : (x1.velocity.x - get_cell(i-1,j,k,f).velocity.x);
                    // double divy = y1.solid ? (0) : (y1.velocity.y - get_cell(i,j-1,k,f).velocity.y);
                    // double divz = z1.solid ? (0) : (z1.velocity.z - get_cell(i,j,k-1,f).velocity.z);
                    double divergence = divx + divy + divz;
                    // if(isnan(divergence)) std::cout << "divergence is nan: " << x1.velocity.x << ", " << y1.velocity.y << ", " << z1.velocity.z << std::endl;
                    // double divergence = (get_cell(i+1,j,k,f).velocity.x - c.velocity.x) +
                    //                     (get_cell(i,j+1,k,f).velocity.y - c.velocity.y) +
                    //                     (get_cell(i,j,k+1,f).velocity.z - c.velocity.z);
                    // double b_i = (-(c.density+1)*cell_width*divergence)/timestep;
                    double b_i = -cell_width*divergence/timestep;
                    // std::cout << b_i << ",";
                    b.push_back(b_i);
                    // b.push_back(-0.5*cell_width*divergence);
                    p.push_back(0.0);
                    // if((i==24 && j==16 && k==24) || (i==35 && j==16 && k==35))
                    // {
                    //     std::cout << "  idx: " << i << std::endl;
                    //     std::cout << "  vel: " << cells[i][j][k].velocity.x << " " << cells[i][j][k].velocity.y << " " << cells[i][j][k].velocity.z << std::endl;
                    //     std::cout << "  vel+1: " << cells[i+1][j][k].velocity.x << " " << cells[i][j+1][k].velocity.y << " " << cells[i][j][k+1].velocity.z << std::endl;
                    //     std::cout << "  idx: " << A(neighbor_map[i][j][k],neighbor_map[i-1][j][k]) << " " << A(neighbor_map[i][j][k],neighbor_map[i][j-1][k]) << " " << A(neighbor_map[i][j][k],neighbor_map[i][j][k-1]) << " " << A(neighbor_map[i][j][k],neighbor_map[i+1][j][k]) << " " << A(neighbor_map[i][j][k],neighbor_map[i][j+1][k]) << " " << A(neighbor_map[i][j][k],neighbor_map[i][j][k+1]) << std::endl;
                    //     std::cout << "  diag: " << A(neighbor_map[i][j][k],neighbor_map[i][j][k]) << std::endl;
                    //     // std::cout << "  A(row)": << A(neighbor_map[i][j][k], neighbor_map[i-1][j][k] << " " << A(neighbor_map[i][j][k], neighbor_map[i][j][k])
                    //     std::cout << "  b_i: " << b_i << std::endl; 
                    //     std::cout << "  divergence: " << divergence << std::endl;

                    // }
                }
            }
        }
    }
    
    // std::ofstream ofs ("A.m", std::ofstream::out);
    // A.write_matlab(ofs, "A");
    // ofs << "lorem ipsum";

    // ofs.close();

    // matrix solve
    PCGSolver<double> solver;
    double residual;
    int iterations;
    std::cout << solver.solve(A, b, p, residual, iterations) << std::endl;

std::cout <<"pressure solved: " << p.size() << std::endl;
    // set pressure
    for(int i=0; i<dimx; i++)
    {
        for(int j=0; j<dimy; j++)
        {
            for(int k=0; k<dimz; k++)
            {
                if(neighbor_map[i][j][k]>=0)
                {
                    cells[i][j][k].pressure = p[neighbor_map[i][j][k]];
                }
                else
                {
                    cells[i][j][k].pressure = 0.0;
                }
            }
        }
    }
std::cout << "applying pressure" << std::endl;
    // apply pressure
    for(int i=0; i<dimx; i++)
    {
        for(int j=0; j<dimy; j++)
        {
            for(int k=0; k<dimz; k++)
            {
                // if(neighbor_map[i][j][k]>=0)
                // {
                double pressure = cells[i][j][k].pressure;
                vec4 gradient = vec4(pressure - get_cell(i-1,j,k,f).pressure,
                                     pressure - get_cell(i,j-1,k,f).pressure,
                                     pressure - get_cell(i,j,k-1,f).pressure);

                // vec4 gradient = vec4(get_cell(i+1,j,k,f).pressure - get_cell(i-1,j,k,f).pressure,
                //                      get_cell(i,j+1,k,f).pressure - get_cell(i,j-1,k,f).pressure,
                //                      get_cell(i,j,k+1,f).pressure - get_cell(i,j,k-1,f).pressure);
                // if(gradient.x>0) std::cout << gradient.x << " " << gradient.y << " " << gradient.z << std::endl;
                // cells[i][j][k].velocity = cells[i][j][k].velocity - ((gradient*timestep)*(1/((cells[i][j][k].density+1)*cell_width)));
                // if(cells[i][j][k].velocity.length() > 3) 
                // {
                //     std::cout << "vel: " << cells[i][j][k].velocity.x << " " << cells[i][j][k].velocity.y << " " << cells[i][j][k].velocity.z << std::endl;
                //     std::cout << "grad: " << gradient.x << " " << gradient.y << " " << gradient.z << std::endl;
                //     std::cout << "dens: " << cells[i][j][k].density+1 << std::endl;
                // }
                cells[i][j][k].velocity = cells[i][j][k].velocity - gradient*(timestep/(cell_width));
                // if((i==24 && j==15 && k==24) || (i==35 && j==15 && k==35))
                // {
                //     std::cout << "  idx: " << i << std::endl;
                //     std::cout << "  vel: " << cells[i][j][k].velocity.x << " " << cells[i][j][k].velocity.y << " " << cells[i][j][k].velocity.z << std::endl;
                //     std::cout << "  pressure: " << pressure << std::endl; 
                //     std::cout << "  gradient: " << gradient.x << " " << gradient.y << " " << gradient.z << std::endl;

                // }
                // if(isnan(gradient.x) || isnan(gradient.y) || isnan(gradient.z)) std::cout << "gradient is nan: " << i << ", " << j << ", " << k << std::endl;
                // }
            }
        }
    }

    // free neighbor map memory
    for(int i=0; i<dimx; i++)
    {
        for(int j=0; j<dimy; j++)
        {
            delete[] neighbor_map[i][j];
        }
        delete[] neighbor_map[i];
    }
    delete[] neighbor_map;
}

void grid::enforce_boundaries()
{
    for(int i=0; i<dimx; i++)
    {
        for(int j=0; j<dimy; j++)
        {
            for(int k=0; k<dimz; k++)
            {
                int f = 1;
                if((get_cell(i-1,j,k,f).solid == 1 && cells[i][j][k].velocity.x < 0.0) ||
                   (get_cell(i+1,j,k,f).solid == 1 && cells[i][j][k].velocity.x > 0.0))
                    cells[i][j][k].velocity.x = 0.0;
                if((get_cell(i,j-1,k,f).solid == 1 && cells[i][j][k].velocity.y < 0.0) ||
                   (get_cell(i,j+1,k,f).solid == 1 && cells[i][j][k].velocity.y > 0.0)) 
                    cells[i][j][k].velocity.y = 0.0;
                if((get_cell(i,j,k-1,f).solid == 1 && cells[i][j][k].velocity.z < 0.0) ||
                   (get_cell(i,j,k+1,f).solid == 1 && cells[i][j][k].velocity.z > 0.0)) 
                    cells[i][j][k].velocity.z = 0.0;
            }
        }
    }
}
void grid::enforce_boundaries_particles()
{
    // TODO: fix strange bouncy collisions
    for(int i=0; i<particles.size(); i++)
    {
        if(particles[i].position.x<lower_bound.x)
        { 
            particles[i].position.x = lower_bound.x+epsilon;
            particles[i].velocity.x = 0.0;
        }
        if(particles[i].position.y<lower_bound.y) 
        {
            particles[i].position.y = lower_bound.y+epsilon;
            particles[i].velocity.y = 0.0;
        }
        if(particles[i].position.z<lower_bound.z) 
        {
            particles[i].position.z = lower_bound.z+epsilon;
            particles[i].velocity.z = 0.0;
        }
        if(particles[i].position.x>upper_bound.x) 
        {
            particles[i].position.x = upper_bound.x-epsilon;
            particles[i].velocity.x = 0.0;
        }
        if(particles[i].position.y>upper_bound.y) 
        {
            particles[i].position.y = upper_bound.y-epsilon;
            particles[i].velocity.y = 0.0;
        }
        if(particles[i].position.z>upper_bound.z) 
        {
            particles[i].position.z = upper_bound.z-epsilon;
            particles[i].velocity.z = 0.0;
        }
    }
}

vec4 grid::multiply_weights_velocity(vec4 index, double* weights)
{
    int flags[8] = {1,1,1,1,1,1,1,1};
    vec4 val = get_cell(index.x  , index.y  , index.z  , flags[0]).velocity*weights[0] +
               get_cell(index.x+1, index.y  , index.z  , flags[1]).velocity*weights[1] +
               get_cell(index.x  , index.y+1, index.z  , flags[2]).velocity*weights[2] +
               get_cell(index.x+1, index.y+1, index.z  , flags[3]).velocity*weights[3] +
               get_cell(index.x  , index.y  , index.z+1, flags[4]).velocity*weights[4] +
               get_cell(index.x+1, index.y  , index.z+1, flags[5]).velocity*weights[5] +
               get_cell(index.x  , index.y+1, index.z+1, flags[6]).velocity*weights[6] +
               get_cell(index.x+1, index.y+1, index.z+1, flags[7]).velocity*weights[7];

    // vec4 val = get_cell(index.x  , index.y  , index.z  , flags[0]).velocity*weights[1] +
    //            get_cell(index.x+1, index.y  , index.z  , flags[1]).velocity*weights[0] +
    //            get_cell(index.x  , index.y+1, index.z  , flags[2]).velocity*weights[3] +
    //            get_cell(index.x+1, index.y+1, index.z  , flags[3]).velocity*weights[2] +
    //            get_cell(index.x  , index.y  , index.z+1, flags[4]).velocity*weights[5] +
    //            get_cell(index.x+1, index.y  , index.z+1, flags[5]).velocity*weights[4] +
    //            get_cell(index.x  , index.y+1, index.z+1, flags[6]).velocity*weights[7] +
    //            get_cell(index.x+1, index.y+1, index.z+1, flags[7]).velocity*weights[6];

    double weight_sum = 0;
    for(int i=0; i<8; i++)
    {
        weight_sum += weights[i]*flags[i];
    }
    val *= (1/weight_sum);
    if(weight_sum==0) val = vec4(0,0,0);
    // if(isnan(val.x) || isnan(val.y) || isnan(val.z))
    // {
    //     std::cout << "weight_sum: " << weight_sum << std::endl;
    //     std::cout << "index: " << index.x << " " << index.y << " " << index.z << std::endl;
    //     for(int i=0; i<8; i++)
    //     {
    //         std::cout << flags[i] << " ";
    //     }
    //     std::cout << std::endl;
    // }
    return val;
}

vec4 grid::multiply_weights_velocity_diff(vec4 index, double* weights)
{
    int flags[8] = {1,1,1,1,1,1,1,1};
    vec4 val = get_cell(index.x  , index.y  , index.z  , flags[0]).velocity_diff*weights[0] +
               get_cell(index.x+1, index.y  , index.z  , flags[1]).velocity_diff*weights[1] +
               get_cell(index.x  , index.y+1, index.z  , flags[2]).velocity_diff*weights[2] +
               get_cell(index.x+1, index.y+1, index.z  , flags[3]).velocity_diff*weights[3] +
               get_cell(index.x  , index.y  , index.z+1, flags[4]).velocity_diff*weights[4] +
               get_cell(index.x+1, index.y  , index.z+1, flags[5]).velocity_diff*weights[5] +
               get_cell(index.x  , index.y+1, index.z+1, flags[6]).velocity_diff*weights[6] +
               get_cell(index.x+1, index.y+1, index.z+1, flags[7]).velocity_diff*weights[7];

    double weight_sum = 0;
    for(int i=0; i<8; i++)
    {
        weight_sum += weights[i]*flags[i];
    }
    val *= (1/weight_sum);
    if(weight_sum==0) val = vec4(0,0,0);

    return val;
}

double grid::multiply_weights_density(vec4 index, double* weights)
{
    int flags[8] = {1,1,1,1,1,1,1,1};
    double val = get_cell(index.x  , index.y  , index.z  , flags[0]).density*weights[0] +
                 get_cell(index.x+1, index.y  , index.z  , flags[1]).density*weights[1] +
                 get_cell(index.x  , index.y+1, index.z  , flags[2]).density*weights[2] +
                 get_cell(index.x+1, index.y+1, index.z  , flags[3]).density*weights[3] +
                 get_cell(index.x  , index.y  , index.z+1, flags[4]).density*weights[4] +
                 get_cell(index.x+1, index.y  , index.z+1, flags[5]).density*weights[5] +
                 get_cell(index.x  , index.y+1, index.z+1, flags[6]).density*weights[6] +
                 get_cell(index.x+1, index.y+1, index.z+1, flags[7]).density*weights[7];

    // double val = get_cell(index.x  , index.y  , index.z  , flags[0]).density*weights[1] +
    //              get_cell(index.x+1, index.y  , index.z  , flags[1]).density*weights[0] +
    //              get_cell(index.x  , index.y+1, index.z  , flags[2]).density*weights[3] +
    //              get_cell(index.x+1, index.y+1, index.z  , flags[3]).density*weights[2] +
    //              get_cell(index.x  , index.y  , index.z+1, flags[4]).density*weights[5] +
    //              get_cell(index.x+1, index.y  , index.z+1, flags[5]).density*weights[4] +
    //              get_cell(index.x  , index.y+1, index.z+1, flags[6]).density*weights[7] +
    //              get_cell(index.x+1, index.y+1, index.z+1, flags[7]).density*weights[6];

    double weight_sum = 0;
    for(int i=0; i<8; i++)
    {
        weight_sum += weights[i]*flags[i];
    }
    val *= (1/weight_sum);
    if(weight_sum==0) val = 0.0;
    return val;
}

double grid::multiply_weights_temperature(vec4 index, double* weights)
{
    int flags[8] = {1,1,1,1,1,1,1,1};
    double val = get_cell(index.x  , index.y  , index.z  , flags[0]).temperature*weights[0] +
                 get_cell(index.x+1, index.y  , index.z  , flags[1]).temperature*weights[1] +
                 get_cell(index.x  , index.y+1, index.z  , flags[2]).temperature*weights[2] +
                 get_cell(index.x+1, index.y+1, index.z  , flags[3]).temperature*weights[3] +
                 get_cell(index.x  , index.y  , index.z+1, flags[4]).temperature*weights[4] +
                 get_cell(index.x+1, index.y  , index.z+1, flags[5]).temperature*weights[5] +
                 get_cell(index.x  , index.y+1, index.z+1, flags[6]).temperature*weights[6] +
                 get_cell(index.x+1, index.y+1, index.z+1, flags[7]).temperature*weights[7];

    // double val = get_cell(index.x  , index.y  , index.z  , flags[0]).temperature*weights[1] +
    //              get_cell(index.x+1, index.y  , index.z  , flags[1]).temperature*weights[0] +
    //              get_cell(index.x  , index.y+1, index.z  , flags[2]).temperature*weights[3] +
    //              get_cell(index.x+1, index.y+1, index.z  , flags[3]).temperature*weights[2] +
    //              get_cell(index.x  , index.y  , index.z+1, flags[4]).temperature*weights[5] +
    //              get_cell(index.x+1, index.y  , index.z+1, flags[5]).temperature*weights[4] +
    //              get_cell(index.x  , index.y+1, index.z+1, flags[6]).temperature*weights[7] +
    //              get_cell(index.x+1, index.y+1, index.z+1, flags[7]).temperature*weights[6];

    double weight_sum = 0;
    for(int i=0; i<8; i++)
    {
        weight_sum += weights[i]*flags[i];
    }
    val *= (1/weight_sum);
    if(weight_sum==0) val = 0.0;
    return val;
}

/*
Postion must be normalized by the cell_width
*/
void grid::get_interpolation_weights(double x, double y, double z, double* weights)
{
    int i = floor(x);
    int j = floor(y);
    int k = floor(z);
    weights[0] = (i+1-x)*(j+1-y)*(k+1-z);
    weights[1] = (x - i)*(j+1-y)*(k+1-z);
    weights[2] = (i+1-x)*(y - j)*(k+1-z);
    weights[3] = (x - i)*(y - j)*(k+1-z);
    weights[4] = (i+1-x)*(j+1-y)*(z - k);
    weights[5] = (x - i)*(j+1-y)*(z - k);
    weights[6] = (i+1-x)*(y - j)*(z - k);
    weights[7] = (x - i)*(y - j)*(z - k);
    // double sum = 0;
// for(int q=0; q<8; q++)
// {
//     sum += weights[q];
// }
// if(sum>1.01||sum<0.99) 
// {
// std::cout << "weights_sum: "<< sum << std::endl;
// }
}

/*
velocity is stored on the minimal faces of each cell
*/
vec4 grid::get_velocity(vec4 position)
{
    vec4 vel;
    vec4 pos_norm;// = position*(1.0/cell_width);
    pos_norm.x = (position.x+size.x/2)/cell_width;
    pos_norm.y = (position.y+size.y/2)/cell_width;
    pos_norm.z = (position.z+size.z/2)/cell_width;
    vec4 idx = index_from_position(position); 
    double weights[8];

    vec4 pos_x(pos_norm); pos_x.y -= 0.5; pos_x.z -= 0.5; //pos_x.x -= 0.5;
    // vec4 idx_x(idx); idx_x.y -= half_cell_width; idx_x.z -= half_cell_width; //idx_x.x -= cell_width/2;
    vec4 idx_x(floor(pos_x.x), floor(pos_x.y), floor(pos_x.z));
    // vec4 idx_x = index_from_position(vec4(position.x,position.y-half_cell_width,position.z-half_cell_width));
    // vec4 idx_x = index_from_position(vec4(position.x,position.y,position.z));
    // vec4 idx_x = index_from_position(vec4(position.x-half_cell_width,position.y-half_cell_width,position.z-half_cell_width));
    get_interpolation_weights(pos_x.x, pos_x.y, pos_x.z, weights);
    vel.x = multiply_weights_velocity(idx_x, weights).x;

    vec4 pos_y(pos_norm); pos_y.x -= 0.5; pos_y.z -= 0.5; //pos_y.y -= 0.5;
    // vec4 idx_y(idx); idx_y.x -= half_cell_width; idx_y.z -= half_cell_width; //idx_y.y -= cell_width/2;
    vec4 idx_y(floor(pos_y.x), floor(pos_y.y), floor(pos_y.z));
    // vec4 idx_y = index_from_position(vec4(position.x-half_cell_width,position.y,position.z-half_cell_width));
    // vec4 idx_y = index_from_position(vec4(position.x,position.y,position.z));
    // vec4 idx_y = index_from_position(vec4(position.x-half_cell_width,position.y-half_cell_width,position.z-half_cell_width));
    get_interpolation_weights(pos_y.x, pos_y.y, pos_y.z, weights);
    vel.y = multiply_weights_velocity(idx_y, weights).y;


    vec4 pos_z(pos_norm); pos_z.x -= 0.5; pos_z.y -= 0.5; //pos_z.z -= 0.5;
    // vec4 idx_z(idx); idx_z.x -= half_cell_width; idx_z.y -= neighbor_map[i][j][k]; //idx_z.z -= cell_width/2;
    vec4 idx_z(floor(pos_z.x), floor(pos_z.y), floor(pos_z.z));
    // vec4 idx_z = index_from_position(vec4(position.x-half_cell_width,position.y-half_cell_width,position.z));
    // vec4 idx_z = index_from_position(vec4(position.x,position.y,position.z));
    // vec4 idx_z = index_from_position(vec4(position.x-half_cell_width,position.y-half_cell_width,position.z-half_cell_width));
    get_interpolation_weights(pos_z.x, pos_z.y, pos_z.z, weights);
    vel.z = multiply_weights_velocity(idx_z, weights).z;

vec4 idx_real = index_from_position(position);
// if(idx_real.x != idx_x.x || idx_real.y != idx_y.y || idx_real.z != idx_z.z)
// {
//     std::cout << "WHOA!!!!!!!" << std::endl;
//     std::cout << "idx_real: " << idx_real.x << " " << idx_real.y << " " << idx_real.z << std::endl;
//     std::cout << "idx: " << idx_x.x << " " << idx_y.y << " " << idx_z.z << std::endl;
// }

if(position.w<1.0)
{
    std::cout << "position: " << position.x << " " << position.y << " " << position.z << std::endl;
    std::cout << "pos_norm: " << pos_norm.x << " " << pos_norm.y << " " << pos_norm.z << std::endl;
    
    std::cout << "idx_real: " << idx_real.x << " " << idx_real.y << " " << idx_real.z << std::endl;
    std::cout << "idx_x: " << idx_x.x << " " << idx_x.y << " " << idx_x.z << std::endl;
    std::cout << "idx_y: " << idx_y.x << " " << idx_y.y << " " << idx_y.z << std::endl;
    std::cout << "idx_z: " << idx_z.x << " " << idx_z.y << " " << idx_z.z << std::endl;
    // std::cout << "idx: " << idx_x.x << " " << idx_y.y << " " << idx_z.z << std::endl;
    // std::cout << "weights: ";
    // for(int i=0; i<8; i++)
    //     std::cout << weights[i] << " ";
    // std::cout << std::endl;
}

    return vel;
}

/*
get interpolated velocity increment for flip
*/
vec4 grid::get_velocity_difference(vec4 position)
{
    vec4 diff;
    vec4 pos_norm;// = position*(1.0/cell_width);
    pos_norm.x = (position.x+size.x/2)/cell_width;
    pos_norm.y = (position.y+size.y/2)/cell_width;
    pos_norm.z = (position.z+size.z/2)/cell_width;
    double weights[8];

    vec4 pos_x(pos_norm); pos_x.y -= 0.5; pos_x.z -= 0.5; //pos_x.x -= 0.5;
    vec4 idx_x(floor(pos_x.x), floor(pos_x.y), floor(pos_x.z));
    get_interpolation_weights(pos_x.x, pos_x.y, pos_x.z, weights);
    diff.x = multiply_weights_velocity_diff(idx_x, weights).x;

    vec4 pos_y(pos_norm); pos_y.x -= 0.5; pos_y.z -= 0.5; //pos_y.y -= 0.5;
    vec4 idx_y(floor(pos_y.x), floor(pos_y.y), floor(pos_y.z));
    get_interpolation_weights(pos_y.x, pos_y.y, pos_y.z, weights);
    diff.y = multiply_weights_velocity_diff(idx_y, weights).y;

    vec4 pos_z(pos_norm); pos_z.x -= 0.5; pos_z.y -= 0.5; //pos_z.z -= 0.5;
    vec4 idx_z(floor(pos_z.x), floor(pos_z.y), floor(pos_z.z));
    get_interpolation_weights(pos_z.x, pos_z.y, pos_z.z, weights);
    diff.z = multiply_weights_velocity_diff(idx_z, weights).z;

    return diff;
}

double grid::get_density(vec4 position)
{
    double density;
    vec4 pos_norm = position*(1.0/cell_width);
    // vec4 pos_norm_center(pos_norm);
    pos_norm -= 0.5;
// std::cout << pos_norm.x << ", " << pos_norm.y << ", " << pos_norm.z << std::endl;
    double weights[8];

    get_interpolation_weights(pos_norm.x, pos_norm.y, pos_norm.z, weights);
    // vec4 rl_idx = index_from_position(position);
    // vec4 position_center(position);
    position -= half_cell_width;
    vec4 idx = index_from_position(position);
// std::cout << idx.x << " " << idx.y << " " << idx.z << " . . ." << std::endl;



// if(position.w<1.0)
// {
//     std::cout << "position: " << position.x << " " << position.y << " " << position.z << std::endl;
//     std::cout << "position_center: " << position_center.x << " " << position_center.y << " " << position_center.z << std::endl;
//     std::cout << "pos_norm_center: " << pos_norm_center.x << " " << pos_norm_center.y << " " << pos_norm_center.z << std::endl;
//     std::cout << "pos_norm: " << pos_norm.x << " " << pos_norm.y << " " << pos_norm.z << std::endl;
//     std::cout << "idx: " << idx.x << " " << idx.y << " " << idx.z << std::endl;
//     std::cout << "rl_idx: " << rl_idx.x << " " << rl_idx.y << " " << rl_idx.z << std::endl;

//     std::cout << "weights: ";
//     for(int i=0; i<8; i++)
//         std::cout << weights[i] << " ";
//     std::cout << std::endl;
// }

    return multiply_weights_density(idx, weights);
}

double grid::get_temperature(vec4 position)
{
    vec4 pos_norm = position*(1.0/cell_width);
    pos_norm -= 0.5;
    double weights[8];

    get_interpolation_weights(pos_norm.x, pos_norm.y, pos_norm.z, weights);
    position -= half_cell_width;
    vec4 idx = index_from_position(position);

    return multiply_weights_temperature(idx, weights);
}

/*
Swap the previous and current grid pointers
*/
void grid::advance_timestep()
{
    cell*** tmp = next_cells;
    next_cells = cells;
    cells = tmp;
}

vec4 grid::position_from_index(vec4 index)
{
    double x = index.x*cell_width;
    double y = index.y*cell_width;
    double z = index.z*cell_width;
    return vec4(x-size.x/2+half_cell_width,y-size.y/2+half_cell_width,z-size.z/2+half_cell_width); 
    // return vec4(x+half_cell_width,y+half_cell_width,z+half_cell_width); 
}

vec4 grid::index_from_position(vec4 position)
{
    int i = floor((position.x+size.x/2)/cell_width);
    int j = floor((position.y+size.y/2)/cell_width);
    int k = floor((position.z+size.x/2)/cell_width);
    // int i = floor((position.x)/cell_width);
    // int j = floor((position.y)/cell_width);
    // int k = floor((position.z)/cell_width);
    return vec4(i,j,k);
}

void grid::write_particles(std::string filename)
{
    std::ofstream file;
    file.open(filename.c_str());

    for(int i=0; i<particles.size(); i++)
    {
        file << particles[i].position.x << " " << particles[i].position.y << " " << particles[i].position.z << " ";
    }

    file.close();
}

void grid::write_density(std::string filename)
{
    // openvdb::initialize();
    // openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create();
    // grid->setTransform(openvdb::math::Transform::createLinearTransform(this->cell_width));
    // grid->setGridClass(openvdb::GRID_FOG_VOLUME);
    // grid->setName("FluidDensity");
    // openvdb::FloatGrid::Accessor accessor = grid->getAccessor();
    // openvdb::Coord ijk;
    // int &i = ijk[0], &j = ijk[1], &k = ijk[2];
    // for(k=0; k<dimz; k++)
    // {
    //     for(j=0; j<dimy; j++)
    //     {
    //         for(i=0; i<dimx; i++)
    //         {
    //             accessor.setValue(ijk, cells[i][j][k].density);
    //         }
    //     }
    // }
    // openvdb::io::File file(filename);
    // openvdb::GridPtrVec container;
    // container.push_back(grid);
    // file.write(container);
    // file.close();


    std::ofstream file;
    file.open(filename.c_str());
    int count = 0;
    for(int k=0; k<dimz; k++)
    {
        for(int j=0; j<dimy; j++)
        {
            for(int i=0; i<dimx; i++)
            {
                std::string delim(",");
                // std::cout << "res: " << resolution << ", count: " << count << std::endl;
                if(count==resolution-1) delim = "";
                file << cells[i][j][k].density << delim;
                count++;
                // if(count==96){
                //     file << std::endl;
                //     count = 0;
                // }
            }
        }
    }
    file.close();
}
