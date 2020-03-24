void add_particle_to_density_field(double rho_field[N],/* the density field
				discretised on a grid of length N */
				int N,
				double cellsize, /* size of a single cell of the grid */
				double particle_pos, /* the particle position */
				double particle_mass) /* the particle mass */
{
/*change ........................................................................................................*/
	double xx = particle_pos / cellsize - 0.5;
/*...............................................................................................................*/
	int i = floor(xx); /* floor(x) truncates all decimal digits of the floating point
		number x, similar to ((int) x) in C (or: it does a strict rounding to the next
	lower integer number; floor(5.9) = 5) */
	double u = xx - i;
	int ii = i + 1;
	if(ii >= N){
	   ii = 0;
	}
	rho_field[i] += (1 - u) * particle_mass / cellsize;
	rho_field[ii] += u * particle_mass / cellsize;
/*change..........................................................................................................*/
	return rho_field[ii] - rho_field[i]
/*................................................................................................................*/
}
double interpolate_force_field_to_particle_position(double force_field[N],/* the force
							field discretised on a grid of length N */
							int N,
							double cellsize, /* size of a single cell of the grid */
							double particle_pos, /* the particle position */
							double particle_mass) /* the particle mass */
{
	double xx = particle_pos / cellsize;
	int i = floor(xx + 0.5);
	int ii = i + 1;
/*change.......see for subpart c...................................................................................*/
        int ii = i + 1;
	
	if(ii >= N) /* for boundary conditions */
		ii = 0;
/*................................................................................................................*/	
/*change.............................................................................................................*/

	double acceleration = u*force_field[i]+(1-u)*forcefield[ii]; /* In order to map the forces in a right way */
/*....................................................................................................................*/	
	return acceleration * particle_mass;
}
