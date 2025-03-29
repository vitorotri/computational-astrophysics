/*
 * Copyright (c) 2025 Vito Romanelli Tricanico
 *
 * SPDX-License-Identifier: BSD-2-Clause (see LICENSE)
 */

#include <stdbool.h>
#include "nbodyheader.h"

typedef struct node{
	double cx, cy, cz; // center coordinates for the node
	double halfl; // halflength
	double CMx, CMy, CMz; // coordinates of the center of mass, if leaf == false
	double M; // monopole, if leaf == false
	bool leaf; // leaf attribute (all children are NULL)
	particle *ptc;
	struct node *child[8];
} node;

node *create_node(double cx, double cy, double cz, double halfl){
	node* newNode = (node *)malloc(sizeof(node));
    newNode->cx = cx;
    newNode->cy = cy;
    newNode->cz = cz;
    newNode->M = 0.0;
    newNode->CMx = 0.0;
    newNode->CMy = 0.0;
    newNode->CMz = 0.0;
    newNode->halfl = halfl;
    newNode->ptc = NULL;
    newNode->leaf = true;
    for (int i = 0; i < 8; ++i) {
        newNode->child[i] = NULL;
    }
    return newNode;
}

int compute_oct(particle *new_ptc, node *parent){
	int oct = 0;
	if (new_ptc->x > parent->cx) oct += 1;
    if (new_ptc->y > parent->cy) oct += 2;
    if (new_ptc->z > parent->cz) oct += 4;
    return oct;
}

int *sign_oct(int oct){
	int *aux = (int *)malloc(3*sizeof(int)); // aux[0] is for x, aux[1] is for y, aux[2] is for z
	aux[0] = (oct == 1 || oct == 3 || oct == 5 || oct == 7) ? 1 : -1;
	aux[1] = (oct == 2 || oct == 3 || oct == 6 || oct == 7) ? 1 : -1;
	aux[2] = (oct == 4 || oct == 5 || oct == 6 || oct == 7) ? 1 : -1;
	return aux;
}

void insert(particle *new_ptc, node *parent){
    if ((parent->leaf == true) && (parent->ptc == NULL)){
        parent->ptc = new_ptc;
    } else {
        int new_oct;
        int *sign;
        double new_halfl;
        new_halfl = 0.5 * parent->halfl;
        new_oct = compute_oct(new_ptc, parent);
        sign = sign_oct(new_oct);
        if (parent->child[new_oct] == NULL) {
            parent->child[new_oct] = create_node(parent->cx + (double)sign[0]*new_halfl,
                                                 parent->cy + (double)sign[1]*new_halfl,
                                                 parent->cz + (double)sign[2]*new_halfl,
                                                 new_halfl);
        }
        insert(new_ptc, parent->child[new_oct]);
        if (parent->ptc != NULL){
            int old_oct;
            old_oct = compute_oct(parent->ptc, parent);
            sign = sign_oct(old_oct);
            if (parent->child[old_oct] == NULL) {
                parent->child[old_oct] = create_node(parent->cx + (double)sign[0]*new_halfl,
                                                     parent->cy + (double)sign[1]*new_halfl,
                                                     parent->cz + (double)sign[2]*new_halfl,
                                                     new_halfl);
            }
            insert(parent->ptc, parent->child[old_oct]);
            parent->ptc = NULL; // Avoid freeing the particle memory, as it's owned by the array
        }
        if (parent->leaf == true) parent->leaf = false;
        free(sign);
    }
}

// depth-first traverse
void print_octree(node* root, int depth, int childNumber) {
    if (root != NULL) {
        for (int i = 0; i < depth - 1; ++i) {
            printf("    ");
        }
        
        if (depth > 0) {
            printf("╰── ");
        }

        if (root->ptc != NULL) {
            printf("* %d %.3lf\n", childNumber,root->M);
        } else {
            printf(". %d %.3lf\n", childNumber,root->M);
        }

        for (int i = 0; i < 8; ++i) {
            print_octree(root->child[i], depth + 1, i);
        }
    }
}

// compute node attributes such as monopole and CMs
// only after tree is created

void compute_attribs(node *current) {
    if (current != NULL) {
    	// after tree is constructed, nodes that contain a particle are leaves, so only check node->ptc != NULL
        if (current->ptc != NULL) {
            current->M = current->ptc->mass;
            current->CMx = current->ptc->x;
            current->CMy = current->ptc->y;
            current->CMz = current->ptc->z;
        }
        else {
        	double totalMass = 0.0;
            double totalCMx = 0.0, totalCMy = 0.0, totalCMz = 0.0;
            // If it's not a leaf, recursively compute attributes for children
            for (int i = 0; i < 8; ++i) {
                if (current->child[i] != NULL) {
                    compute_attribs(current->child[i]);
                    totalMass += current->child[i]->M;
                    totalCMx += current->child[i]->M * current->child[i]->CMx;
                    totalCMy += current->child[i]->M * current->child[i]->CMy;
                    totalCMz += current->child[i]->M * current->child[i]->CMz;
                }
            }
            current->M = totalMass;
            current->CMx = totalCMx / totalMass;
            current->CMy = totalCMy / totalMass;
            current->CMz = totalCMz / totalMass;
        }
    }
}

// computes theta
double theta(double halfl, double y){
	return halfl/y;
}

void quadrupole(node *current, double Q[]){
	//printf("[OK] Hello Q assignments\n");
	for (int i = 0; i < 8; i++){
		// need to iterate over all child particles, not only on leaf nodes
		if (current->child[i] != NULL){
			//printf("[OK] Hello inside Q 1st if\n");
			Q[0] += 2*current->child[i]->M*(current->CMx - current->child[i]->CMx)*(current->CMx - current->child[i]->CMx);
			Q[1] += 3*current->child[i]->M*(current->CMx - current->child[i]->CMx)*(current->CMy - current->child[i]->CMy);
			Q[2] += 3*current->child[i]->M*(current->CMx - current->child[i]->CMx)*(current->CMz - current->child[i]->CMz);
			Q[3] += 2*current->child[i]->M*(current->CMy - current->child[i]->CMy)*(current->CMy - current->child[i]->CMy);
			Q[4] += 3*current->child[i]->M*(current->CMy - current->child[i]->CMy)*(current->CMz - current->child[i]->CMz);
			Q[5] += 2*current->child[i]->M*(current->CMz - current->child[i]->CMz)*(current->CMz - current->child[i]->CMz);
		}
	}
}

// compute force using barnes-hut algorithm for given particle at one iteration
void barnes_hut(node *current, particle *p, double soft, double lim){
    double y[3], abs_y;
    y[0] = p->x - current->CMx;
    y[1] = p->y - current->CMy;
    y[2] = p->z - current->CMz;
    abs_y = sqrt(y[0]*y[0] + y[1]*y[1] + y[2]*y[2]);	
    if (theta(current->halfl, abs_y) < lim){
    	//printf("[OK] Hello inside barnes_hut if\n");
    	double Q[6];
    	for (int i = 0; i < 6; i++){
			Q[i] = 0.0;
		}
    	quadrupole(current,Q);
    	double y_TQy = (Q[0]*y[0]*y[0] + Q[3]*y[1]*y[1] + Q[5]*y[2]*y[2]
    					+ 2*(Q[1]*y[0]*y[1] + Q[2]*y[0]*y[2] + Q[4]*y[1]*y[2]));
    			
    	p->fx += massa*(-current->M*y[0]/pow(abs_y,3) + (Q[0]*y[0] + Q[1]*y[1] + Q[2]*y[2])/pow(abs_y,5)
    		- 2.5*y_TQy*y[0]/pow(abs_y,7));
    		
    	p->fy += massa*(-current->M*y[1]/pow(abs_y,3) + (Q[1]*y[0] + Q[3]*y[1] + Q[4]*y[2])/pow(abs_y,5)
    		- 2.5*y_TQy*y[1]/pow(abs_y,7));
    			
    	p->fz += massa*(-current->M*y[2]/pow(abs_y,3) + (Q[2]*y[0] + Q[4]*y[1] + Q[5]*y[2])/pow(abs_y,5)
    		- 2.5*y_TQy*y[2]/pow(abs_y,7));
    }
    else {
       	for (int i = 0; i < 8; ++i) {
           	if (current->child[i] != NULL){
           		if ((current->child[i]->ptc != NULL) && (current->child[i]->ptc->ID != p->ID)){
    				double dr;
    				dr = sqrt(pow(p->x - current->child[i]->ptc->x,2) + pow(p->y - current->child[i]->ptc->y,2)
    					+ pow(p->z - current->child[i]->ptc->z,2) + soft*soft);
    				p->fx += -massa*massa*(p->x - current->child[i]->ptc->x)/pow(dr,3);
    				p->fy += -massa*massa*(p->y - current->child[i]->ptc->y)/pow(dr,3);
    				p->fz += -massa*massa*(p->z - current->child[i]->ptc->z)/pow(dr,3);
           		}
           		else {
          			barnes_hut(current->child[i], p, soft, lim);
           		}
           	}
       	}
	}
}

void free_octree(node *root){
    if (root != NULL) {
        for (int i = 0; i < 8; ++i) {
            free_octree(root->child[i]);
        }
        free(root);
    }
}
