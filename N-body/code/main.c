/*
 * Copyright (c) 2025 Vito Romanelli Tricanico
 *
 * SPDX-License-Identifier: BSD-2-Clause (see LICENSE)
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include <string.h>
#include "treecode.h"
//#include <omp.h>

// usage: ./main [char] [N_bins] [softening_val]
//		  ./main [char] [N_bins] [softening_val] [opening_angle]
// Example:
// 
// 		  ./main T3 1 0.3 0.6

// For N = 1000, use T3 with softening = 0.3 and opening = 0.6
// Parameters such as dt could be adjusted according to the number of (super)particles.

int main(int argc,char *argv[]){
	
	clock_t tic, tac;
	double t_elapsed;
	double r_max, mis;
	char *eptr;
	int N_bins = atoi(argv[2]);
	int d = 0; // auxiliar index for while loop file scan
	double r_hm;
	double cluster_radius = 10.0; // Radius of the star cluster
	// declare output struct
	out dens_bins;
	dens_bins.dens = (double *) malloc(N_bins*sizeof(double));
	dens_bins.radii = (double *) malloc(N_bins*sizeof(double));
	
	// declare array of structs
	particle *ptcs;
    ptcs = (particle *)malloc(N * sizeof(particle));
    if (ptcs == NULL) {
        fprintf(stderr, "Memory allocation failed.\n");
        return EXIT_FAILURE;
    }    srand(time(NULL)); // Seed for random number generation
    generate_cluster(ptcs, cluster_radius);
	
	compute_radii(ptcs);
	heapsort(ptcs,N); // sort particles array by radial distance, once
	attribute_IDs(ptcs);
	r_max = ptcs[N-1].r;
	r_hm = ptcs[N/2].r;
	mis = mean_is(ptcs, N);
	printf("r_max = %lf * 10^-2 kpc   r_hm = %lf * 10^-2 kpc   mis = %lf * 10^-2 kpc\n", r_max, r_hm, mis);
	
	char start[100];
	printf("\n\nProgram start?\n\nThen type Yes\n\n");
	scanf("%s", &start);
	if (strcmp(start,"Yes") != 0){
		return 1;
	}
	printf("\n[OK] Program starting...\n\n");
	
	// density
	if (strcmp(argv[1],"D") == 0){
		dens_bins = density(ptcs,N_bins,r_max);
		// results for plotting densities
		for (int i = 0; i < N_bins; i++){
			printf("%lf %lf\n",dens_bins.radii[i],dens_bins.dens[i]);
		}
	}
	// compare elapsed times for computing forces
	else if (strcmp(argv[1],"F") == 0){
		// loops over different softening values
		for (int i = 3; i < argc; i++){
			/*
			tic = clock();
			brute_forces_naive(ptcs,strtod(argv[i],&eptr));
			tac = clock();
			t_elapsed = ((double) tac - tic) / CLOCKS_PER_SEC;
			printf("%lf %lf\n",sqrt(pow(ptcs[N-1].fx,2) + pow(ptcs[N-1].fy,2) + pow(ptcs[N-1].fz,2)), t_elapsed);
			*/
			/*
			tic = clock();
			brute_forces(ptcs,strtod(argv[i],&eptr));
			tac = clock();
			t_elapsed = ((double) tac - tic) / CLOCKS_PER_SEC;
			printf("%lf %lf\n",sqrt(pow(ptcs[N-1].fx,2) + pow(ptcs[N-1].fy,2) + pow(ptcs[N-1].fz,2)), t_elapsed);
			*/
			///*
			tic = clock();
			direct_forces_parallel(ptcs,strtod(argv[i],&eptr),N);
			tac = clock();
			t_elapsed = ((double) tac - tic) / CLOCKS_PER_SEC;
			printf("%lf %lf\n",abs_force(ptcs[N-1]), t_elapsed);
			//printf("%lf %lf\n",abs_force(ptcs[N-3]), t_elapsed);
			//printf("%lf %lf\n",abs_force(ptcs[0]), t_elapsed);
			//*/
		}
	}
	// tests treecode forces for different opening angles
	else if (strcmp(argv[1],"T1") == 0){
		// different opening-angles
		for (int j = 4; j < argc; j++){
			// create root node
			node *root = create_node(0.0, 0.0, 0.0, r_max);
			// insert particles, creating the octree
			for (int i = 0; i < N; i++){
				insert(&ptcs[i], root);
			}
			initialize_forces(ptcs, N); // zero the forces
			compute_attribs(root); // compute monopoles and CMs
			tic = clock();
			for (int i = 0; i < N; i++){
				barnes_hut(root, &ptcs[i], strtod(argv[3],&eptr), strtod(argv[j],&eptr));
			}
			free_octree(root);
			tac = clock();
			t_elapsed = ((double) tac - tic) / CLOCKS_PER_SEC;
			//printf("%lf %lf\n",abs_force(ptcs[N-1]), t_elapsed);
			printf("%lf %lf\n",abs_force(ptcs[0]), t_elapsed);
		}
	}
	// checks for numerical relaxation
	else if (strcmp(argv[1],"R") == 0){
		double t_cross, v_c, dt, t_max, soft, vms;
		int n_max; // maximum number of iterations
		v_c = sqrt(0.5*N*massa/r_hm);
		t_cross = r_hm/v_c;
		t_max = 2*t_cross;
		dt = 0.1*t_cross;
		n_max = t_max/dt;
		soft = 0.12*mis;
		FILE *file = fopen("output_R.txt", "w");
		printf("\nWill now do %d iterations...\n\n", n_max);
		for (int n = 0; n < n_max; n++){
			direct_forces_parallel(ptcs, soft, N);
			for (int i = 0; i < N; i++){
				leapfrog(&ptcs[i],dt);
			}
			vms = v_ms(ptcs);
			fprintf(file, "%lf %lf\n", ((double)n * dt)/t_cross, vms);
			// if softening is too small -> Miller's instability
			if (!(vms > 100*massa/soft)){
				break;
			}
			printf("[OK] from iteration %d\n",n);
		}
		fclose(file);
	}
	// compute using direct forces, but in parallel
	else if (strcmp(argv[1],"T2") == 0){
		double dt, t_max, soft;
		int n_max; // maximum number of iterations
		t_max = 50;
		dt = 0.1;
		n_max = t_max/dt;
		printf("t_max = %lf\n",t_max);
		
		int N_PARTICLES = 4;
		particle particles[N_PARTICLES];
    	initialize_particles(particles, N_PARTICLES);
		
		// The mean interparticle separation changes during evolution. That does not necessarily mean we need to
		// change the softening for each iteration.
		
		//printf("%f\n", soft);
		FILE *file = fopen("output_T2.txt", "w");
		// leapfrog time integration
		soft = mis;
		for (int n = 0; n < n_max; n++){
			//mis = mean_is(particles, N_PARTICLES);
			//soft = 0.5*mis;
			direct_forces_parallel(particles, soft, N_PARTICLES);
			// for each particle
			for (int i = 0; i < N_PARTICLES; i++){
				fprintf(file, "%d %lf %lf %lf\n", n, particles[i].x, particles[i].y, particles[i].z);
				leapfrog(&particles[i],dt);
			}
			printf("[OK] from iteration %d\n",n);
		}
		fclose(file);
	}
	// compute using Barnes-Hut
	else if (strcmp(argv[1],"T3") == 0){
		// t_max = 300*tcross
		// dt = 0.3*tcross
		// soft = 0.05
		// opening = 0.27
	
		double t_cross, v_c, dt, t_max, soft;
		int n_max, n_shift; // maximum number of iterations
		v_c = sqrt(0.5*N*massa/r_hm);
		t_cross = r_hm/v_c;
		t_max = 100*t_cross; // end time (page 197 - Binney & Tremaine)
		dt = 0.01*t_cross;
		n_max = t_max/dt;
		n_shift = (int)(0.005*n_max);
		// sum 1 to n_shift to avoid division by zero on conditional print
		if (n_shift == 0) n_shift += 1;
		//printf("%lf\n",0.02*sqrt(pow(r_hm,3)/(0.5*N*massa)));
		printf("dt = %lf\n",dt);
		printf("n_max = %d\n",n_max);
		printf("n_shift = %d\n",n_shift);
		// leapfrog time integration
		printf("[OK] Starting %d iterations\n",n_max + n_shift);
		mis = mean_is(ptcs, N);
		printf("[OK] mis computation\n");
		soft = mis;
		FILE *file = fopen("output.txt", "w");
		for (int n = 0; n < n_max + n_shift; n++){
		//for (int n = 0; n < 1; n++){
		
			initialize_forces(ptcs,N); // zero the forces
			//printf("[OK] initialized forces\n");
			// create root node
			node *root = create_node(0.0, 0.0, 0.0, r_max);
			//printf("[OK] created octree's root node\n");
			// insert particles, creating the octree
			for (int i = 0; i < N; i++){
				insert(&ptcs[i], root);
			}
			//printf("[OK] particles insertions\n");
			compute_attribs(root); // compute monopoles and CMs
			//printf("[OK] computed nodes attributes\n");
			for (int i = 0; i < N; i++){
				barnes_hut(root, &ptcs[i], soft, strtod(argv[3],&eptr));
			}
			//printf("[OK] Barnes-Hut\n");
			free_octree(root);
			//printf("[OK] freed octree\n");
			// for each particle
			for (int i = 0; i < N; i++){
				if (n % n_shift == 0){
					//printf("%d %lf %lf %lf\n", n, ptcs[i].x, ptcs[i].y, ptcs[i].z);
					fprintf(file, "%d %lf %lf %lf\n", n, ptcs[i].x, ptcs[i].y, ptcs[i].z);
				}
				leapfrog(&ptcs[i],dt);
			}
			
			compute_radii(ptcs);
			heapsort(ptcs,N); // sort particles array by radial distance
			attribute_IDs(ptcs);
			r_max = ptcs[N-1].r;
			
			printf("[OK] from iteration %d\n",n);
		}
		fclose(file);
	}
	else if (strcmp(argv[1],"BKP") == 0){
		// softening 0.05
		// opening 0.27
		d = 0;
		FILE *fptr = fopen("input_BKP.txt","r"); // pointer to file
		while (fscanf(fptr,"%lf %lf %lf", &ptcs[d].x, &ptcs[d].y, &ptcs[d].z, &ptcs[d].vx, &ptcs[d].vy, &ptcs[d].vz) == 6){
			d++; // increment d
		}
		fclose(fptr);
		// to implement if needed
	}
	free(dens_bins.dens);
	free(dens_bins.radii);
	free(ptcs); // maybe free right after creating tree, otherwise doubled data
	return 0;
}
