/*
 * Copyright (c) 2025 Vito Romanelli Tricanico
 *
 * SPDX-License-Identifier: BSD-2-Clause (see LICENSE)
 */
 
// Code that creates a simple quadtree, the 2D equivalent of an octree (3D.

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#define L 10.0

typedef struct particle{
	double x, y;
	double mass;
} particle;

typedef struct node{
	double cx, cy; // center coordinates for the node
	double halfl; // halflength
	bool leaf; // leaf attribute
	double M;
	particle *ptc;
	struct node *child[4];
} node;

particle create_particle(double x, double y, double mass){
	particle ptc;
	ptc.x = x;
	ptc.y = y;
	ptc.mass = mass;
	return ptc;
}

node *create_node(double cx, double cy, double halfl){
	node *newNode = (node *)malloc(sizeof(node));
    newNode->cx = cx;
    newNode->cy = cy;
    newNode->ptc = NULL;
    newNode->halfl = halfl;
    newNode->leaf = true;
    newNode->M = 0.0;
    for (int i = 0; i < 4; ++i) {
        newNode->child[i] = NULL;
    }
    return newNode;
}

int compute_quad(particle *new_ptc, node *parent){
	int quad = 0;
	if (new_ptc->x > parent->cx) quad += 1;
    if (new_ptc->y > parent->cy) quad += 2;
    //if (new_ptc->z > parent->cz) quad += 4; // octree case
    return quad;
}

int *sign_quad(int quad){
	int *aux = (int *)malloc(2*sizeof(int)); // aux[0] is for x, aux[1] is for y
	aux[0] = (quad == 1 || quad == 3) ? 1 : -1;
	aux[1] = (quad == 2 || quad == 3) ? 1 : -1;
	return aux;
}


void insert(particle *new_ptc, node *parent) {
    if ((parent->leaf == true) && (parent->ptc == NULL)) {
        parent->ptc = new_ptc;
        parent->leaf = false;
    } else {
        int new_quad;
        int *sign;
        double new_halfl;
        new_halfl = 0.5 * parent->halfl;
        new_quad = compute_quad(new_ptc, parent);
        sign = sign_quad(new_quad);
        if (parent->child[new_quad] == NULL) {
            parent->child[new_quad] = create_node(parent->cx + (double)sign[0] * new_halfl,
                                                  parent->cy + (double)sign[1] * new_halfl, new_halfl);
        }
        insert(new_ptc, parent->child[new_quad]);
        if (parent->ptc != NULL) {
            int old_quad;
            old_quad = compute_quad(parent->ptc, parent);
            sign = sign_quad(old_quad);
            if (parent->child[old_quad] == NULL) {
                parent->child[old_quad] = create_node(parent->cx + (double)sign[0] * new_halfl,
                                                      parent->cy + (double)sign[1] * new_halfl, new_halfl);
            }
            insert(parent->ptc, parent->child[old_quad]);
            parent->ptc = NULL; // Avoid freeing the particle memory, as it's owned by the array
        }
        free(sign);
    }
}

// depth-first traversal
void print_quadtree(node* root, int depth, int childNumber) {
    if (root != NULL) {
        for (int i = 0; i < depth - 1; ++i) {
            printf("    ");
        }

        if (depth > 0) {
            printf("╰── ");
        }

		// depth-first traversal
        if (root->ptc != NULL) {
            printf("* %d %.3lf\n", childNumber,root->M);
        } else {
            printf(". %d %.3lf\n", childNumber,root->M);
        }

        for (int i = 0; i < 4; ++i) {
            print_quadtree(root->child[i], depth + 1, i);
        }
    }
}

void free_quadtree(node *root) {
    if (root != NULL) {
        for (int i = 0; i < 4; ++i) {
            free_quadtree(root->child[i]);
        }
        free(root);
    }
}

// after tree is created only
void compute_mass(node *current) {
    if (current != NULL) {
    	// after tree is constructed, nodes that contain a particle are leaves, so only check node->ptc != NULL
        if (current->ptc != NULL) {
            current->M = current->ptc->mass;
        }
        else {
            // If it's not a leaf, recursively compute masses for children
            //current->M = 0.0;
            for (int i = 0; i < 4; ++i) {
                if (current->child[i] != NULL) {
                    compute_mass(current->child[i]);
                    current->M += current->child[i]->M;
                }
            }
        }
    }
}

int main() {
    node *root = create_node(0.0, 0.0, L);
    // carefull to level of recursion if particle are way to close
    particle ptc0 = create_particle(-8.0, -8.0, 1.0);
    particle ptc1 = create_particle(4.9, 5.9, 1.0);
    particle ptc2 = create_particle(1.7, 9.7, 1.0);
    particle ptc3 = create_particle(-2.0, -2.0, 1.0);
    particle ptc4 = create_particle(-8.0, 3.0, 1.0);
    insert(&ptc0, root);
    insert(&ptc1, root);
    insert(&ptc2, root);
    insert(&ptc3, root);
    insert(&ptc4, root);
    compute_mass(root);
    print_quadtree(root,0,0);
    free_quadtree(root);
    return 0;
}
