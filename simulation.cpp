#include "simulation.h"

const int cell::DENSITY = 0;
const int cell::PRESSURE = 1;
const int cell::TEMPERATURE = 2;
const int cell::MAG_CURL = 3; // TODO

const int cell::VELOCITY = 0;
const int cell::OLD_VELOCITY = 1;
const int cell::DELTA_VELOCITY = 2;
const int cell::WEIGHTS = 3;
const int cell::CURL = 4; // TODO

cell::cell()
{
    solid = 0;
    envelope = 0;
    scalar_fields[DENSITY] = 0.0;
    scalar_fields[PRESSURE] = 0.0;
    scalar_fields[TEMPERATURE] = 0.0;
    scalar_fields[MAG_CURL] = 0.0;

    vector_fields[VELOCITY] = vec4();
    vector_fields[WEIGHTS] = vec4();
    vector_fields[CURL] = vec4();
}

double cell::get_scalar_field(int field)
{
    return this->scalar_fields[field];
}

void cell::set_scalar_field(int field, double value)
{
    this->scalar_fields[field] = value;
}

vec4 cell::get_vector_field(int field)
{
    return this->vector_fields[field];
}

void cell::set_vector_field(int field, vec4 value)
{
    this->vector_fields[field] = value;
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
                vcur = cells[i][j][k].get_vector_field(cell::VELOCITY).length();
                if(vcur > vmax && (cells[i][j][k].get_scalar_field(cell::DENSITY) > epsilon || cells[i][j][k].envelope==0)) vmax = vcur;
            }
        }
    }
    double k_cfl = 2.0;
    if(vmax==0) timestep = default_timestep;
    timestep = cell_width*k_cfl/vmax;
    if(timestep>default_timestep) timestep = default_timestep;
    if(timestep<epsilon) timestep = epsilon;
    std::cout << "max_vel: " << vmax << std::endl;
    std::cout << "timestep: " << timestep << std::endl;
}

void grid::advect_particles()
{
    for(int i=0; i<particles.size(); i++)
    {
        particles[i].position += particles[i].velocity*timestep; // TODO: RK3
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
                vec4 p = position_from_index(vec4(i,j,k));
                vec4 v = get_interpolated_vector(p, cell::VELOCITY);
                vec4 p_x(p); p_x.x -= half_cell_width;
                vec4 p_y(p); p_y.y -= half_cell_width;
                vec4 p_z(p); p_z.z -= half_cell_width;
                vec4 v_x = get_interpolated_vector(p_x, cell::VELOCITY);
                vec4 v_y = get_interpolated_vector(p_y, cell::VELOCITY);
                vec4 v_z = get_interpolated_vector(p_z, cell::VELOCITY);
                v = get_interpolated_vector(vec4(p.x-0.5*timestep*v.x, p.y-0.5*timestep*v.y, p.z-0.5*timestep*v.z), cell::VELOCITY);
                v_x = get_interpolated_vector(vec4(p_x.x-0.5*timestep*v_x.x, p_x.y-0.5*timestep*v_x.y, p_x.z-0.5*timestep*v_x.z), cell::VELOCITY);
                v_y = get_interpolated_vector(vec4(p_y.x-0.5*timestep*v_y.x, p_y.y-0.5*timestep*v_y.y, p_y.z-0.5*timestep*v_y.z), cell::VELOCITY);
                v_z = get_interpolated_vector(vec4(p_z.x-0.5*timestep*v_z.x, p_z.y-0.5*timestep*v_z.y, p_z.z-0.5*timestep*v_z.z), cell::VELOCITY);
                if(isnan(v.x) || isnan(v.y) || isnan(v.z)) std::cout << i << " " << j << " " << k << std::endl;

                p = p - v*timestep;
                p_x = p_x - v_x*timestep;
                p_y = p_y - v_y*timestep;
                p_z = p_z - v_z*timestep;

                this->next_cells[i][j][k].vector_fields[cell::VELOCITY].x = get_interpolated_vector(p_x, cell::VELOCITY).x;
                this->next_cells[i][j][k].vector_fields[cell::VELOCITY].y = get_interpolated_vector(p_y, cell::VELOCITY).y;
                this->next_cells[i][j][k].vector_fields[cell::VELOCITY].z = get_interpolated_vector(p_z, cell::VELOCITY).z;
                
                this->next_cells[i][j][k].set_scalar_field(cell::DENSITY, get_interpolated_scalar(p, cell::DENSITY));
                this->next_cells[i][j][k].set_scalar_field(cell::TEMPERATURE, get_interpolated_scalar(p, cell::TEMPERATURE));
                this->next_cells[i][j][k].set_scalar_field(cell::PRESSURE, 0.0);
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
                            double cur_density = cells[i][j][k].get_scalar_field(cell::DENSITY);
                            cells[i][j][k].set_scalar_field(cell::DENSITY, cur_density + sources[s].scale*timestep);
                            // cells[i][j][k].density += sources[s].scale*timestep;
                            double cur_temperature = cells[i][j][k].get_scalar_field(cell::TEMPERATURE);
                            cells[i][j][k].set_scalar_field(cell::TEMPERATURE, cur_temperature + sources[s].scale*timestep);
                            // cells[i][j][k].temperature += sources[s].scale*timestep;
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
                // blur += get_cell(i+x, j+y, k+z, flag).temperature;
                blur += get_cell(i+x, j+y, k+z, flag).get_scalar_field(cell::TEMPERATURE);
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
                double cur_density = cells[i][j][k].get_scalar_field(cell::DENSITY);
                double cur_temperature = cells[i][j][k].get_scalar_field(cell::TEMPERATURE);
                // if( cells[i][j][k].density > epsilon ||  cells[i][j][k].temperature > epsilon)
                if(cur_density > epsilon ||  cur_temperature > epsilon)
                {

                    // dissipation
                    // cells[i][j][k].density = std::max(0.0, cells[i][j][k].density - dissipation*timestep);
                    cells[i][j][k].set_scalar_field(cell::DENSITY, std::max(0.0, cur_density - dissipation*timestep));
                    // // next_cells[i][j][k].temperature = std::max(0.0, cells[i][j][k].temperature - dissipation*5.0*timestep);
                    // cells[i][j][k].set_scalar_field(cell::TEMPERATURE, std::max(0.0, cur_temperature - dissipation*5.0*timestep));

                    // // next_cells[i][j][k].temperature = blur_temperature(i, j, k, 3);
                    // next_cells[i][j][k].set_scalar_field(cell::TEMPERATURE, blur_temperature(i, j, k, 3));
                    // // if(cells[i][j][k].temperature <0.5 && cells[i][j][k].temperature > 0.0001) std::cout << cells[i][j][k].temperature << std::endl;
                    // // if(isnan(cells[i][j][k].temperature)) std::cout << "nan temperature" << std::endl;

                    // vec4 corner = position_from_index(vec4(i,j,k)) - vec4(0,half_cell_width,0);
                    // // double temp = get_temperature(corner);
                    // double temp = get_interpolated_scalar(corner, cell::TEMPERATURE);
                    // // double temp = cells[i][j][k].temperature;/
                    // double k_rise = temp*buoyancy;
                    // double k_fall = 0.0;//cells[i][j][k].density*-buoyancy/4;
                    // cells[i][j][k].vector_fields[cell::VELOCITY] = cells[i][j][k].vector_fields[cell::VELOCITY] + vec4(0.0, 1.0, 0.0)*((k_rise+k_fall)*timestep);


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
                        if(rad.length() < sources[s].radius) cells[i][j][k].vector_fields[cell::VELOCITY] += dir*5.0*timestep;
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
                if(cur_density > epsilon || cells[i][j][k].envelope==0)
                    cells[i][j][k].vector_fields[cell::VELOCITY] +=(vec4(0,-1,0)*(gravity*timestep));
                    // cells[i][j][k].velocity +=(vec4(0,-.707,.707)*(gravity*timestep));

                // if(isnan(cells[i][j][k].velocity.x) || isnan(cells[i][j][k].velocity.y) || isnan(cells[i][j][k].velocity.z))
                // {
                //     std::cout << "external_forces nan: " << cells[i][j][k].velocity.x << ", " << cells[i][j][k].velocity.y << ", " << cells[i][j][k].velocity.z << std::endl;
                // }
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
    //             // vec4 vx0 = get_velocity(vec4(pos.x-half_cell_width,pos.y,pos.z));
    //             vec4 vx0 = get_interpolated_vector(vec4(pos.x-half_cell_width,pos.y,pos.z), cell::VELOCITY);
    //             // vec4 vx1 = get_velocity(vec4(pos.x+half_cell_width,pos.y,pos.z));
    //             vec4 vx1 = get_interpolated_vector(vec4(pos.x+half_cell_width,pos.y,pos.z), cell::VELOCITY);
    //             // vec4 vy0 = get_velocity(vec4(pos.x,pos.y-half_cell_width,pos.z));
    //             vec4 vy0 = get_interpolated_vector(vec4(pos.x,pos.y-half_cell_width,pos.z), cell::VELOCITY);
    //             // vec4 vy1 = get_velocity(vec4(pos.x,pos.y+half_cell_width,pos.z));
    //             vec4 vy1 = get_interpolated_vector(vec4(pos.x,pos.y+half_cell_width,pos.z), cell::VELOCITY);
    //             // vec4 vz0 = get_velocity(vec4(pos.x,pos.y,pos.z-half_cell_width));
    //             vec4 vz0 = get_interpolated_vector(vec4(pos.x,pos.y,pos.z-half_cell_width), cell::VELOCITY);
    //             // vec4 vz1 = get_velocity(vec4(pos.x,pos.y,pos.z+half_cell_width));
    //             vec4 vz1 = get_interpolated_vector(vec4(pos.x,pos.y,pos.z+half_cell_width), cell::VELOCITY);
                
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
    //                 cells[i][j][k].vector_fields[cell::VELOCITY] += vort*timestep;
    //         }
    //     }
    // }
}

void grid::extrapolate_velocity(int num_cells)
{
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
                        double count = 0.0;
                        vec4 old_vel = cells[i][j][k].vector_fields[cell::VELOCITY];
                        cells[i][j][k].vector_fields[cell::VELOCITY] = vec4(0,0,0);
                        for(int c=0; c<6; c++)
                        {
                            if(neighbors[c]->envelope==layer-1)
                            {
                                cells[i][j][k].vector_fields[cell::VELOCITY] += neighbors[c]->vector_fields[cell::VELOCITY];
                                count += 1;
                            }
                        }
                        if(count)
                        {
                            cells[i][j][k].vector_fields[cell::VELOCITY] *= (1.0/count);
                            cells[i][j][k].envelope = layer;
                        }
                        // don't average components at fluid boundaries
                        if(neighbors[0]->envelope==0) cells[i][j][k].vector_fields[cell::VELOCITY].x = old_vel.x;
                        if(neighbors[1]->envelope==0) cells[i][j][k].vector_fields[cell::VELOCITY].y = old_vel.y;
                        if(neighbors[2]->envelope==0) cells[i][j][k].vector_fields[cell::VELOCITY].z = old_vel.z;
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
                if(cells[i][j][k].get_scalar_field(cell::DENSITY) > epsilon)
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
                                if(flag && (neighbor.get_scalar_field(cell::DENSITY) <= epsilon)) 
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
                if((cells[i][j][k].get_scalar_field(cell::DENSITY) > epsilon))// || (cells[i][j][k].envelope > 0))
                {
                    neighbor_map[i][j][k] = A_idx;
                    A_idx++;
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
                    if(c.get_scalar_field(cell::DENSITY) > epsilon)
                        A.set_element(neighbor_map[i][j][k], neighbor_map[i-1][j][k], -1);
                    else if(c.solid)
                        non_solids--;
                    c = get_cell(i+1,j,k,f);
                    if(c.get_scalar_field(cell::DENSITY) > epsilon)
                        A.set_element(neighbor_map[i][j][k], neighbor_map[i+1][j][k], -1);
                    else if(c.solid)
                        non_solids--;
                    c = get_cell(i,j-1,k,f);
                    if(c.get_scalar_field(cell::DENSITY) > epsilon)
                        A.set_element(neighbor_map[i][j][k], neighbor_map[i][j-1][k], -1);
                    else if(c.solid)
                        non_solids--;
                    c = get_cell(i,j+1,k,f);
                    if(c.get_scalar_field(cell::DENSITY) > epsilon)
                        A.set_element(neighbor_map[i][j][k], neighbor_map[i][j+1][k], -1);
                    else if(c.solid)
                        non_solids--;
                    c = get_cell(i,j,k-1,f);
                    if(c.get_scalar_field(cell::DENSITY) > epsilon)
                        A.set_element(neighbor_map[i][j][k], neighbor_map[i][j][k-1], -1);
                    else if(c.solid)
                        non_solids--;
                    c = get_cell(i,j,k+1,f);
                    if(c.get_scalar_field(cell::DENSITY) > epsilon)
                        A.set_element(neighbor_map[i][j][k], neighbor_map[i][j][k+1], -1);
                    else if(c.solid)
                        non_solids--;
                    A.set_element(neighbor_map[i][j][k], neighbor_map[i][j][k], non_solids);
                    
                    c = get_cell(i,j,k,f);
                    cell x1 = get_cell(i+1,j,k,f);
                    cell y1 = get_cell(i,j+1,k,f);
                    cell z1 = get_cell(i,j,k+1,f);
                    double divx = x1.solid ? (0) : (x1.vector_fields[cell::VELOCITY].x - c.vector_fields[cell::VELOCITY].x);
                    double divy = y1.solid ? (0) : (y1.vector_fields[cell::VELOCITY].y - c.vector_fields[cell::VELOCITY].y);
                    double divz = z1.solid ? (0) : (z1.vector_fields[cell::VELOCITY].z - c.vector_fields[cell::VELOCITY].z);

                    double divergence = divx + divy + divz;
                    double b_i = -cell_width*divergence/timestep;
                    b.push_back(b_i);
                    p.push_back(0.0);
                }
            }
        }
    }
    
    // std::ofstream ofs ("A.m", std::ofstream::out);
    // A.write_matlab(ofs, "A");

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
                    cells[i][j][k].set_scalar_field(cell::PRESSURE, p[neighbor_map[i][j][k]]);
                }
                else
                {
                    cells[i][j][k].set_scalar_field(cell::PRESSURE, 0.0);
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
                double pressure = cells[i][j][k].get_scalar_field(cell::PRESSURE);
                vec4 gradient = vec4(pressure - get_cell(i-1,j,k,f).get_scalar_field(cell::PRESSURE),
                                     pressure - get_cell(i,j-1,k,f).get_scalar_field(cell::PRESSURE),
                                     pressure - get_cell(i,j,k-1,f).get_scalar_field(cell::PRESSURE));
                // TODO: don't apply pressure on solid borders

                cells[i][j][k].vector_fields[cell::VELOCITY] -= gradient*(timestep/(cell_width));
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
                if((get_cell(i-1,j,k,f).solid == 1 && cells[i][j][k].vector_fields[cell::VELOCITY].x < 0.0) ||
                   (get_cell(i+1,j,k,f).solid == 1 && cells[i][j][k].vector_fields[cell::VELOCITY].x > 0.0))
                    cells[i][j][k].vector_fields[cell::VELOCITY].x = 0.0;
                if((get_cell(i,j-1,k,f).solid == 1 && cells[i][j][k].vector_fields[cell::VELOCITY].y < 0.0) ||
                   (get_cell(i,j+1,k,f).solid == 1 && cells[i][j][k].vector_fields[cell::VELOCITY].y > 0.0)) 
                    cells[i][j][k].vector_fields[cell::VELOCITY].y = 0.0;
                if((get_cell(i,j,k-1,f).solid == 1 && cells[i][j][k].vector_fields[cell::VELOCITY].z < 0.0) ||
                   (get_cell(i,j,k+1,f).solid == 1 && cells[i][j][k].vector_fields[cell::VELOCITY].z > 0.0)) 
                    cells[i][j][k].vector_fields[cell::VELOCITY].z = 0.0;
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

vec4 grid::multiply_weights_vector(vec4 index, double* weights, int field)
{
    int flags[8] = {1,1,1,1,1,1,1,1};
    vec4 val = get_cell(index.x  , index.y  , index.z  , flags[0]).get_vector_field(field)*weights[0] +
               get_cell(index.x+1, index.y  , index.z  , flags[1]).get_vector_field(field)*weights[1] +
               get_cell(index.x  , index.y+1, index.z  , flags[2]).get_vector_field(field)*weights[2] +
               get_cell(index.x+1, index.y+1, index.z  , flags[3]).get_vector_field(field)*weights[3] +
               get_cell(index.x  , index.y  , index.z+1, flags[4]).get_vector_field(field)*weights[4] +
               get_cell(index.x+1, index.y  , index.z+1, flags[5]).get_vector_field(field)*weights[5] +
               get_cell(index.x  , index.y+1, index.z+1, flags[6]).get_vector_field(field)*weights[6] +
               get_cell(index.x+1, index.y+1, index.z+1, flags[7]).get_vector_field(field)*weights[7];

    double weight_sum = 0;
    for(int i=0; i<8; i++)
    {
        weight_sum += weights[i]*flags[i];
    }
    val *= (1/weight_sum);
    if(weight_sum==0) val = vec4(0,0,0);

    return val;
}

double grid::multiply_weights_scalar(vec4 index, double* weights, int field)
{
    int flags[8] = {1,1,1,1,1,1,1,1};

    double val = get_cell(index.x  , index.y  , index.z  , flags[0]).get_scalar_field(field)*weights[0] +
                 get_cell(index.x+1, index.y  , index.z  , flags[1]).get_scalar_field(field)*weights[1] +
                 get_cell(index.x  , index.y+1, index.z  , flags[2]).get_scalar_field(field)*weights[2] +
                 get_cell(index.x+1, index.y+1, index.z  , flags[3]).get_scalar_field(field)*weights[3] +
                 get_cell(index.x  , index.y  , index.z+1, flags[4]).get_scalar_field(field)*weights[4] +
                 get_cell(index.x+1, index.y  , index.z+1, flags[5]).get_scalar_field(field)*weights[5] +
                 get_cell(index.x  , index.y+1, index.z+1, flags[6]).get_scalar_field(field)*weights[6] +
                 get_cell(index.x+1, index.y+1, index.z+1, flags[7]).get_scalar_field(field)*weights[7];

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
}

/*
get interpolated vector field
*/
vec4 grid::get_interpolated_vector(vec4 position, int field)
{
    vec4 val;
    vec4 pos_norm;
    pos_norm.x = (position.x+size.x/2)/cell_width;
    pos_norm.y = (position.y+size.y/2)/cell_width;
    pos_norm.z = (position.z+size.z/2)/cell_width;
    double weights[8];

    vec4 pos_x(pos_norm); pos_x.y -= 0.5; pos_x.z -= 0.5;
    vec4 idx_x(floor(pos_x.x), floor(pos_x.y), floor(pos_x.z));
    get_interpolation_weights(pos_x.x, pos_x.y, pos_x.z, weights);
    val.x = multiply_weights_vector(idx_x, weights, field).x;

    vec4 pos_y(pos_norm); pos_y.x -= 0.5; pos_y.z -= 0.5;
    vec4 idx_y(floor(pos_y.x), floor(pos_y.y), floor(pos_y.z));
    get_interpolation_weights(pos_y.x, pos_y.y, pos_y.z, weights);
    val.y = multiply_weights_vector(idx_y, weights, field).y;

    vec4 pos_z(pos_norm); pos_z.x -= 0.5; pos_z.y -= 0.5;
    vec4 idx_z(floor(pos_z.x), floor(pos_z.y), floor(pos_z.z));
    get_interpolation_weights(pos_z.x, pos_z.y, pos_z.z, weights);
    val.z = multiply_weights_vector(idx_z, weights, field).z;

    return val;
}

double grid::get_interpolated_scalar(vec4 position, int field)
{
    vec4 pos_norm = position*(1.0/cell_width);
    pos_norm -= 0.5;
    double weights[8];

    get_interpolation_weights(pos_norm.x, pos_norm.y, pos_norm.z, weights);
    position -= half_cell_width;
    vec4 idx = index_from_position(position);

    return multiply_weights_scalar(idx, weights, field);
}

void grid::transfer_to_grid()
{
    for(int i=0; i<dimx; i++)
    {
        for(int j=0; j<dimy; j++)
        {
            for(int k=0; k<dimz; k++)
            {
                cells[i][j][k].set_scalar_field(cell::DENSITY, 0.0);
                cells[i][j][k].set_vector_field(cell::VELOCITY, vec4());
                cells[i][j][k].set_vector_field(cell::OLD_VELOCITY, vec4());
                cells[i][j][k].set_vector_field(cell::WEIGHTS, vec4());
                cells[i][j][k].set_scalar_field(cell::PRESSURE, 0.0);
                cells[i][j][k].envelope = -1;
            }
        }
    }

    for(int i=0; i<particles.size(); i++)
    {
        vec4 idx = index_from_position(particles[i].position);
        int x = idx.x; int y = idx.y; int z = idx.z;
        vec4 pos_norm = (particles[i].position + size*0.5)*(1/cell_width);
        double wx[8];
        vec4 pos_x(pos_norm); pos_x.y -= 0.5; pos_x.z -= 0.5;
        get_interpolation_weights(pos_x.x, pos_x.y, pos_x.z, wx);
        double wy[8];
        vec4 pos_y(pos_norm); pos_y.x -= 0.5; pos_y.z -= 0.5;
        get_interpolation_weights(pos_y.x, pos_y.y, pos_y.z, wy);
        double wz[8];
        vec4 pos_z(pos_norm); pos_z.x -= 0.5; pos_z.y -= 0.5;
        get_interpolation_weights(pos_z.x, pos_z.y, pos_z.z, wz);

        double cur_density = cells[x][y][z].get_scalar_field(cell::DENSITY);
        cells[x][y][z].set_scalar_field(cell::DENSITY, cur_density+1);
        int xx =floor(pos_x.x); int xy =floor(pos_x.y); int xz =floor(pos_x.z);
        int yx =floor(pos_y.x); int yy =floor(pos_y.y); int yz =floor(pos_y.z);
        int zx =floor(pos_z.x); int zy =floor(pos_z.y); int zz =floor(pos_z.z);

        int flag = 1;
        get_cell_ref(xx,xy,xz,flag)->vector_fields[cell::VELOCITY].x += particles[i].velocity.x*wx[0];
        get_cell_ref(xx,xy,xz,flag)->vector_fields[cell::WEIGHTS].x += wx[0];
        get_cell_ref(yx,yy,yz,flag)->vector_fields[cell::VELOCITY].y += particles[i].velocity.y*wy[0];
        get_cell_ref(yx,yy,yz,flag)->vector_fields[cell::WEIGHTS].y += wy[0];
        get_cell_ref(zx,zy,zz,flag)->vector_fields[cell::VELOCITY].z += particles[i].velocity.z*wz[0];
        get_cell_ref(zx,zy,zz,flag)->vector_fields[cell::WEIGHTS].z += wz[0];

        get_cell_ref(xx+1,xy,xz,flag)->vector_fields[cell::VELOCITY].x += particles[i].velocity.x*wx[1];
        get_cell_ref(xx+1,xy,xz,flag)->vector_fields[cell::WEIGHTS].x += wx[1];
        get_cell_ref(yx+1,yy,yz,flag)->vector_fields[cell::VELOCITY].y += particles[i].velocity.y*wy[1];
        get_cell_ref(yx+1,yy,yz,flag)->vector_fields[cell::WEIGHTS].y += wy[1];
        get_cell_ref(zx+1,zy,zz,flag)->vector_fields[cell::VELOCITY].z += particles[i].velocity.z*wz[1];
        get_cell_ref(zx+1,zy,zz,flag)->vector_fields[cell::WEIGHTS].z += wz[1];

        get_cell_ref(xx,xy+1,xz,flag)->vector_fields[cell::VELOCITY].x += particles[i].velocity.x*wx[2];
        get_cell_ref(xx,xy+1,xz,flag)->vector_fields[cell::WEIGHTS].x += wx[2];
        get_cell_ref(yx,yy+1,yz,flag)->vector_fields[cell::VELOCITY].y += particles[i].velocity.y*wy[2];
        get_cell_ref(yx,yy+1,yz,flag)->vector_fields[cell::WEIGHTS].y += wy[2];
        get_cell_ref(zx,zy+1,zz,flag)->vector_fields[cell::VELOCITY].z += particles[i].velocity.z*wz[2];
        get_cell_ref(zx,zy+1,zz,flag)->vector_fields[cell::WEIGHTS].z += wz[2];

        get_cell_ref(xx+1,xy+1,xz,flag)->vector_fields[cell::VELOCITY].x += particles[i].velocity.x*wx[3];
        get_cell_ref(xx+1,xy+1,xz,flag)->vector_fields[cell::WEIGHTS].x += wx[3];
        get_cell_ref(yx+1,yy+1,yz,flag)->vector_fields[cell::VELOCITY].y += particles[i].velocity.y*wy[3];
        get_cell_ref(yx+1,yy+1,yz,flag)->vector_fields[cell::WEIGHTS].y += wy[3];
        get_cell_ref(zx+1,zy+1,zz,flag)->vector_fields[cell::VELOCITY].z += particles[i].velocity.z*wz[3];
        get_cell_ref(zx+1,zy+1,zz,flag)->vector_fields[cell::WEIGHTS].z += wz[3];

        get_cell_ref(xx,xy,xz+1,flag)->vector_fields[cell::VELOCITY].x += particles[i].velocity.x*wx[4];
        get_cell_ref(xx,xy,xz+1,flag)->vector_fields[cell::WEIGHTS].x += wx[4];
        get_cell_ref(yx,yy,yz+1,flag)->vector_fields[cell::VELOCITY].y += particles[i].velocity.y*wy[4];
        get_cell_ref(yx,yy,yz+1,flag)->vector_fields[cell::WEIGHTS].y += wy[4];
        get_cell_ref(zx,zy,zz+1,flag)->vector_fields[cell::VELOCITY].z += particles[i].velocity.z*wz[4];
        get_cell_ref(zx,zy,zz+1,flag)->vector_fields[cell::WEIGHTS].z += wz[4];

        get_cell_ref(xx+1,xy,xz+1,flag)->vector_fields[cell::VELOCITY].x += particles[i].velocity.x*wx[5];
        get_cell_ref(xx+1,xy,xz+1,flag)->vector_fields[cell::WEIGHTS].x += wx[5];
        get_cell_ref(yx+1,yy,yz+1,flag)->vector_fields[cell::VELOCITY].y += particles[i].velocity.y*wy[5];
        get_cell_ref(yx+1,yy,yz+1,flag)->vector_fields[cell::WEIGHTS].y += wy[5];
        get_cell_ref(zx+1,zy,zz+1,flag)->vector_fields[cell::VELOCITY].z += particles[i].velocity.z*wz[5];
        get_cell_ref(zx+1,zy,zz+1,flag)->vector_fields[cell::WEIGHTS].z += wz[5];

        get_cell_ref(xx,xy+1,xz+1,flag)->vector_fields[cell::VELOCITY].x += particles[i].velocity.x*wx[6];
        get_cell_ref(xx,xy+1,xz+1,flag)->vector_fields[cell::WEIGHTS].x += wx[6];
        get_cell_ref(yx,yy+1,yz+1,flag)->vector_fields[cell::VELOCITY].y += particles[i].velocity.y*wy[6];
        get_cell_ref(yx,yy+1,yz+1,flag)->vector_fields[cell::WEIGHTS].y += wy[6];
        get_cell_ref(zx,zy+1,zz+1,flag)->vector_fields[cell::VELOCITY].z += particles[i].velocity.z*wz[6];
        get_cell_ref(zx,zy+1,zz+1,flag)->vector_fields[cell::WEIGHTS].z += wz[6];

        get_cell_ref(xx+1,xy+1,xz+1,flag)->vector_fields[cell::VELOCITY].x += particles[i].velocity.x*wx[7];
        get_cell_ref(xx+1,xy+1,xz+1,flag)->vector_fields[cell::WEIGHTS].x += wx[7];
        get_cell_ref(yx+1,yy+1,yz+1,flag)->vector_fields[cell::VELOCITY].y += particles[i].velocity.y*wy[7];
        get_cell_ref(yx+1,yy+1,yz+1,flag)->vector_fields[cell::WEIGHTS].y += wy[7];
        get_cell_ref(zx+1,zy+1,zz+1,flag)->vector_fields[cell::VELOCITY].z += particles[i].velocity.z*wz[7];
        get_cell_ref(zx+1,zy+1,zz+1,flag)->vector_fields[cell::WEIGHTS].z += wz[7];
    }

    for(int i=0; i<dimx; i++)
    {
        for(int j=0; j<dimy; j++)
        {
            for(int k=0; k<dimz; k++)
            {
                if(cells[i][j][k].vector_fields[cell::WEIGHTS].length()>0)
                {
                    vec4 cell_weight = cells[i][j][k].get_vector_field(cell::WEIGHTS);
                    cell_weight.x>0 ? cells[i][j][k].vector_fields[cell::VELOCITY].x /= cell_weight.x : cells[i][j][k].vector_fields[cell::VELOCITY].x = 0;
                    cell_weight.y>0 ? cells[i][j][k].vector_fields[cell::VELOCITY].y /= cell_weight.y : cells[i][j][k].vector_fields[cell::VELOCITY].y = 0;
                    cell_weight.z>0 ? cells[i][j][k].vector_fields[cell::VELOCITY].z /= cell_weight.z : cells[i][j][k].vector_fields[cell::VELOCITY].z = 0;
                    cells[i][j][k].set_vector_field(cell::OLD_VELOCITY, cells[i][j][k].get_vector_field(cell::VELOCITY));

                    cells[i][j][k].envelope = 0;
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
                vec4 delta_vel = cells[i][j][k].get_vector_field(cell::VELOCITY) - cells[i][j][k].get_vector_field(cell::OLD_VELOCITY);
                cells[i][j][k].set_vector_field(cell::DELTA_VELOCITY, delta_vel);
            }
        }
    }

    double pic = 1.0 - flip;
    for(int i=0; i<particles.size(); i++)
    {
        particles[i].velocity += get_interpolated_vector(particles[i].position, cell::DELTA_VELOCITY);
        particles[i].velocity *= flip;
        particles[i].velocity += get_interpolated_vector(particles[i].position, cell::VELOCITY)*pic;
    }
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
}

vec4 grid::index_from_position(vec4 position)
{
    int i = floor((position.x+size.x/2)/cell_width);
    int j = floor((position.y+size.y/2)/cell_width);
    int k = floor((position.z+size.x/2)/cell_width);
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
                if(count==resolution-1) delim = "";
                file << cells[i][j][k].get_scalar_field(cell::DENSITY) << delim;
                count++;
            }
        }
    }
    file.close();
}
