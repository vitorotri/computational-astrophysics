/*
 * Copyright (c) 2025 Vito Romanelli Tricanico
 *
 * SPDX-License-Identifier: BSD-2-Clause (see LICENSE)
 */

#define N 1000
#define massa 1.0
#define G 1.0 // Gravitational constant, set to 1 for simplicity

// struct for single particle
typedef struct particle{
	double mass;
	double x, y, z;
	double vx, vy, vz;
	double soft;
	double pot;
	double r;
	double fx, fy, fz;
	int ID;
} particle;

typedef struct out{
	double *dens;
	double *radii;
} out;

double distance(double x1, double x2, double y1, double y2, double z1, double z2){
    double dx = x1 - x2;
    double dy = y1 - y2;
    double dz = z1 - z2;
    return sqrt(dx*dx + dy*dy + dz*dz);
}

// compute radii
// after array is sorted, no need to compute d_max or r_max, it is on the last position of the array
void compute_radii(particle ptcs[]){
	for (int i = 0; i < N; i++){
		ptcs[i].r = sqrt(ptcs[i].x*ptcs[i].x + ptcs[i].y*ptcs[i].y + ptcs[i].z*ptcs[i].z);
	}
}

// after sorting attribute IDs to particles
void attribute_IDs(particle ptcs[]){
	for (int i = 0; i < N; i++){
		ptcs[i].ID = i;
	}
}

// computes the absolute value of the force for given particle
double abs_force(particle p){
	return sqrt(p.fx*p.fx + p.fy*p.fy + p.fz*p.fz);
}

//heapsort functions for array (from Rosetta Code)
int maxheap(particle *a, int n, int i, int j, int k) {
    int m = i;
    if (j < n && a[j].r > a[m].r) {
        m = j;
    }
    if (k < n && a[k].r > a[m].r) {
        m = k;
    }
    return m;
}

void downheap(particle *a, int n, int i) {
    while (1) {
        int j = maxheap(a, n, i, 2 * i + 1, 2 * i + 2);
        if (j == i) {
            break;
        }
        particle t = a[i];
        a[i] = a[j];
        a[j] = t;
        i = j;
    }
}

void heapsort(particle *a, int n) {
    int i;
    for (i = (n - 2) / 2; i >= 0; i--) {
        downheap(a, n, i);
    }
    for (i = 0; i < n; i++) {
        particle t = a[n - i - 1];
        a[n - i - 1] = a[0];
        a[0] = t;
        downheap(a, n - i - 1, 0);
    }
}

// Function to generate a star cluster with N particles
void generate_cluster(particle *ptcs, double cluster_radius) {
    for (int i = 0; i < N; i++) {
        // Assign unique ID
        ptcs[i].ID = i;

        // Equal mass for all particles
        ptcs[i].mass = 1.0;

        // Assign a small softening length
        ptcs[i].soft = 0.01;

        // Generate positions within a sphere (Plummer-like distribution)
        double r = cluster_radius * pow(((double)rand() / RAND_MAX), 1.0 / 3.0); // Radial distance
        double theta = ((double)rand() / RAND_MAX) * M_PI;                      // Polar angle
        double phi = ((double)rand() / RAND_MAX) * (2.0 * M_PI);                // Azimuthal angle

        ptcs[i].x = r * sin(theta) * cos(phi);
        ptcs[i].y = r * sin(theta) * sin(phi);
        ptcs[i].z = r * cos(theta);

        // Store radial distance for convenience
        ptcs[i].r = sqrt(ptcs[i].x * ptcs[i].x + ptcs[i].y * ptcs[i].y + ptcs[i].z * ptcs[i].z);

        // Generate velocities assuming virial equilibrium (v ~ sqrt(GM/r))
        double speed = sqrt(G * N / r); // Approximate velocity magnitude
        double vtheta = ((double)rand() / RAND_MAX) * M_PI;
        double vphi = ((double)rand() / RAND_MAX) * (2.0 * M_PI);

        ptcs[i].vx = speed * sin(vtheta) * cos(vphi);
        ptcs[i].vy = speed * sin(vtheta) * sin(vphi);
        ptcs[i].vz = speed * cos(vtheta);

        // Initialize potential and forces to zero
        ptcs[i].pot = 0.0;
        ptcs[i].fx = 0.0;
        ptcs[i].fy = 0.0;
        ptcs[i].fz = 0.0;
    }
}

// returns array of densities for each shell (spherical bins), central bin radiis and cumulative mass
// ptcs: array of particle struct
// N_bins: number of shells (spherical bins)
// r_max: the radial position of the most distant particle, predetermined in the main program
out density(particle ptcs[],int N_bins,double r_max){
	// density = sum_mass_bin/volume_bin
	double ratio = r_max/(exp((double)N_bins) - 1);
	double sum, lower, upper;
	out outp;
	outp.dens = (double *) malloc(N_bins*sizeof(double));
	outp.radii = (double *) malloc(N_bins*sizeof(double));
	// innermost loops run first
	for (int j = 0; j < N_bins; j++){
		sum = 0.0;
		if (j == 0){
			lower = 0.0;
		}
		else {
			lower = exp(j)*ratio;
		}
		upper = exp(j+1)*ratio;
		for (int i = 0; i < N; i++){
			if (ptcs[i].r >= lower && ptcs[i].r < upper) {
				//sum += ptcs[i].mass*M_0/(4*M_PI*(pow(upper,3) - pow(lower,3))/3);
				sum++;
			}
		}
		outp.dens[j] = sum;
		outp.radii[j] = 0.5*(lower + upper); // central bin radii
	}
	return outp;
}

// initialize forces of particles to 0
void initialize_forces(particle ptcs[], int N_max){
	for (int i = 0; i < N_max; i++){
		ptcs[i].fx = 0.0;
		ptcs[i].fy = 0.0;
		ptcs[i].fz = 0.0;
		ptcs[i].pot = 0.0;
	}
}

// naive approach for compute forces by direct summation
void brute_forces_naive(particle ptcs[],double soft, int N_max){
	double Fx, Fy, Fz, dr;
	// compute forces for each particle. Remind that G = 1, and Newton's 3rd Law
	initialize_forces(ptcs, N_max);
	for (int i = 0; i < N; i++){
		for (int j = 0; j < N; j++){
			if (j != i){
				dr = sqrt(pow(ptcs[i].x - ptcs[j].x,2)
					+ pow(ptcs[i].y - ptcs[j].y,2)
					+ pow(ptcs[i].z - ptcs[j].z,2)
					+ soft*soft);
				Fx = ptcs[i].mass*ptcs[j].mass*(ptcs[j].x - ptcs[i].x)/pow(dr,3); // single force
				Fy = ptcs[i].mass*ptcs[j].mass*(ptcs[j].y - ptcs[i].y)/pow(dr,3); // single force
				Fz = ptcs[i].mass*ptcs[j].mass*(ptcs[j].z - ptcs[i].z)/pow(dr,3); // single force
				ptcs[i].fx += Fx;
				ptcs[i].fy += Fy;
				ptcs[i].fz += Fz;
			}
		}
	}
}

// routine to compute forces by direct summation (non-naive approach)
void brute_forces(particle ptcs[],double soft, int N_max){
	double Fx, Fy, Fz, dr;
	// compute forces for each particle. Remind that G = 1, and Newton's 3rd Law
	initialize_forces(ptcs, N_max);
	for (int i = 0; i < N; i++){
		for (int j = i+1; j < N; j++){
			dr = sqrt(pow(ptcs[i].x - ptcs[j].x,2)  + pow(ptcs[i].y - ptcs[j].y,2) + pow(ptcs[i].z - ptcs[j].z,2)
				+ soft*soft);
			Fx = ptcs[i].mass*ptcs[j].mass*(ptcs[j].x - ptcs[i].x)/pow(dr,3); // single force
			Fy = ptcs[i].mass*ptcs[j].mass*(ptcs[j].y - ptcs[i].y)/pow(dr,3); // single force
			Fz = ptcs[i].mass*ptcs[j].mass*(ptcs[j].z - ptcs[i].z)/pow(dr,3); // single force
			ptcs[i].fx += Fx;
			ptcs[j].fx -= Fx;
			ptcs[i].fy += Fy;
			ptcs[j].fy -= Fy;
			ptcs[i].fz += Fz;
			ptcs[j].fz -= Fz;
		}
	}
}

// routine to compute forces by direct summation (parallel) and potentials
void direct_forces_parallel(particle ptcs[], double soft, int N_max){
	double Fx, Fy, Fz, dr;
	// compute forces for each particle. Remind that G = 1, and Newton's 3rd Law
	initialize_forces(ptcs, N_max);
	#pragma acc parallel loop copy(ptcs[:N_max]) private(Fx, Fy, Fz, dr)
	for (int i = 0; i < N_max; i++){
		for (int j = i+1; j < N_max; j++){
			dr = sqrt(pow(ptcs[i].x - ptcs[j].x,2)  + pow(ptcs[i].y - ptcs[j].y,2) + pow(ptcs[i].z - ptcs[j].z,2)
				+ soft*soft);
			Fx = ptcs[i].mass*ptcs[j].mass*(ptcs[j].x - ptcs[i].x)/pow(dr,3); // single force
			Fy = ptcs[i].mass*ptcs[j].mass*(ptcs[j].y - ptcs[i].y)/pow(dr,3); // single force
			Fz = ptcs[i].mass*ptcs[j].mass*(ptcs[j].z - ptcs[i].z)/pow(dr,3); // single force
			#pragma acc atomic
			ptcs[i].fx += Fx;
			#pragma acc atomic
			ptcs[j].fx -= Fx;
			#pragma acc atomic
			ptcs[i].fy += Fy;
			#pragma acc atomic
			ptcs[j].fy -= Fy;
			#pragma acc atomic
			ptcs[i].fz += Fz;
			#pragma acc atomic
			ptcs[j].fz -= Fz;
			
			ptcs[i].pot += -massa/dr;
			ptcs[j].pot += -massa/dr;
		}
	}
}
/*
// estimates mean interparticle separation
double mean_is(particle ptcs[], int N_max){
	double mis = 0.0;
	int paircount = 0;
	// loop similar to the force for mis
	for (int i = 0; i < N_max; i++){
		for (int j = i+1; j < N_max; j++){
			mis += distance(ptcs[i].x, ptcs[j].x, ptcs[i].y, ptcs[j].y, ptcs[i].z, ptcs[j].z);
			paircount++;
		}
	}
	return mis/paircount;
}
*/
// computes mean interparticle separation
double mean_is(particle ptcs[], int N_p){
	return pow((4.0/3.0)*3.14159*pow(ptcs[N_p/2].r,3)/(0.5*N_p), 1.0/3.0);
}

// computes positions and velocities for a particle in one time-step
void leapfrog(particle *p, double dt){
	double v12_x, v12_y, v12_z;

	// Calculate half-step velocities
    v12_x = p->vx + 0.5 * (p->fx / p->mass) * dt;
    v12_y = p->vy + 0.5 * (p->fy / p->mass) * dt;
    v12_z = p->vz + 0.5 * (p->fz / p->mass) * dt;

    // Update positions
    p->x += v12_x * dt;
    p->y += v12_y * dt;
    p->z += v12_z * dt;

    // Update velocities with the new forces
    p->vx = v12_x + 0.5 * (p->fx / p->mass) * dt;
    p->vy = v12_y + 0.5 * (p->fy / p->mass) * dt;
    p->vz = v12_z + 0.5 * (p->fz / p->mass) * dt;
}

// for testing N-body model
void initialize_particles(particle *particles, int n_particles) {
    // Assuming you have an array of particles
    for (int i = 0; i < n_particles; ++i) {
        particles[i].mass = 1.0;
        particles[i].soft = 0.01;
        particles[i].pot = 0.0;
        particles[i].r = 0.0;
        particles[i].fx = particles[i].fy = particles[i].fz = 0.0;
    }

    particles[0].mass = 1.0;
    particles[0].x = 3.0;
    particles[0].y = 3.0;
    particles[0].z = 0.0;
    particles[0].vx = 0.1;
    particles[0].vy = 0.5;
    particles[0].vz = 0.0;
    particles[0].ID = 1;

    particles[1].mass = 1.0;
    particles[1].x = 2.0;
    particles[1].y = 5.0;
    particles[1].z = 0.0;
    particles[1].vx = 0.5;
    particles[1].vy = 0.5;
    particles[1].vz = 0.0;
    particles[1].ID = 2;
    
    particles[2].mass = 1.0;
    particles[2].x = 5.0;
    particles[2].y = 2.0;
    particles[2].z = 0.0;
    particles[2].vx = 0.5;
    particles[2].vy = 0.2;
    particles[2].vz = 0.0;
    particles[2].ID = 3;
    
    particles[3].mass = 1.0;
    particles[3].x = 0.6;
    particles[3].y = 2.5;
    particles[3].z = 0.0;
    particles[3].vx = 0.5;
    particles[3].vy = 0.5;
    particles[3].vz = 0.0;
    particles[3].ID = 4;
    
}

// computes mean squared velocity
double v_ms(particle ptcs[]){
	double v = 0.0;
	for (int i = 0; i < N; i++){
		v += ptcs[i].vx*ptcs[i].vx + ptcs[i].vy*ptcs[i].vy + ptcs[i].vz*ptcs[i].vz; // absolute value of velocity, squared
	}
	return v/N;
}

// computes the total energy of the system
double E_total(particle ptcs[]){
	double T, W, E;
	T = 0.0;
	W = 0.0;
	for (int i = 0; i < N; i++){
		T += massa*sqrt(ptcs[i].vx*ptcs[i].vx + ptcs[i].vy*ptcs[i].vy + ptcs[i].vz*ptcs[i].vz);
		W += massa*ptcs[i].pot;// eq 27 from Springel
	}
	E = 0.5*(T + W);
	return E;
}
