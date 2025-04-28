import matplotlib.pyplot as plt
import math
import numpy as np
from heapq import heappush, heapreplace
import random

class Particle:
    def __init__(self, r, mass=0.5):
        self.r = np.array(r)
        self.mass = mass
        self.rho = 0.0

class Cell:
    def __init__(self, rLow, rHigh, lower, upper):
        self.rLow = np.array(rLow)
        self.rHigh = np.array(rHigh)
        self.center = (self.rLow + self.rHigh)/2
        self.iLower = lower
        self.iUpper = upper
        self.pLower = None
        self.pUpper = None
        
    def isleaf(self):
        return self.pLower is None and self.pUpper is None

    def celldist2(self, r):
        d1 = r - self.rHigh
        d2 = self.rLow - r
        d1 = np.maximum(d1, d2)
        d1 = np.maximum(d1, 0)
        return d1.dot(d1)

class PriorityQueue:
    def __init__(self, k):
        self.heap = []
        sentinel = (-np.inf, None)
        for _ in range(k):
            heappush(self.heap, sentinel)
    
    @property
    def key(self):
        return -self.heap[0][0]
        
    def replace(self, d2, p):
        heapreplace(self.heap, (-d2, p))

def partition(A, i, j, v, d):
    while i < j:
        if A[i].r[d] < v:
            i += 1
        elif A[j].r[d] >= v:
            j -= 1
        else:
            A[i], A[j] = A[j], A[i]
    return i+1 if A[j].r[d] < v else i

def build_tree(A, root, dim):
    v = (root.rLow[dim] + root.rHigh[dim])/2
    s = partition(A, root.iLower, root.iUpper, v, dim)
    
    if dim == 0:
        rLow_upper = [v, root.rLow[1]]
        rHigh_lower = [v, root.rHigh[1]]
    else:
        rLow_upper = [root.rLow[0], v]
        rHigh_lower = [root.rHigh[0], v]

    if s > root.iLower:
        cLeft = Cell(root.rLow, rHigh_lower, root.iLower, s-1)
        root.pLower = cLeft
        if (s - root.iLower) > 8:
            build_tree(A, cLeft, 1-dim)
            
    if s <= root.iUpper:
        cRight = Cell(rLow_upper, root.rHigh, s, root.iUpper)
        root.pUpper = cRight
        if (root.iUpper - s + 1) > 8:
            build_tree(A, cRight, 1-dim)

def neighbor_search(pq, root, particles, r, rOffset):
    ri = r + rOffset
    if root.pLower and root.pUpper:
        if root.pLower.celldist2(ri) < root.pUpper.celldist2(ri):
            neighbor_search(pq, root.pLower, particles, r, rOffset)
            neighbor_search(pq, root.pUpper, particles, r, rOffset)
        else:
            neighbor_search(pq, root.pUpper, particles, r, rOffset)
            neighbor_search(pq, root.pLower, particles, r, rOffset)
    elif root.pLower:
        neighbor_search(pq, root.pLower, particles, r, rOffset)
    elif root.pUpper:
        neighbor_search(pq, root.pUpper, particles, r, rOffset)
    else:
        for j in range(root.iLower, root.iUpper+1):
            d2 = np.sum((particles[j].r - ri)**2)
            if d2 < pq.key:
                pq.replace(d2, particles[j])

def periodic_search(pq, root, particles, r, period):
    for y in [0, -period[1], period[1]]:
        for x in [0, -period[0], period[0]]:
            neighbor_search(pq, root, particles, r, np.array([x, y]))

def tophat_density(heap, R):
    total_mass = sum(p.mass for _, p in heap)
    return total_mass / ((4/3)*np.pi*R**3)

def monaghan_density(heap, h):
    total = 0.0
    coeff = 40 / (7*np.pi*h**2)
    for d2, p in heap:
        r = math.sqrt(abs(d2))
        q = r/h
        if q <= 0.5:
            W = 6*q**3 - 6*q**2 + 1
        elif q <= 1:
            W = 2*(1 - q)**3
        else:
            continue
        total += coeff * p.mass * W
    return total

# Create particles
Np = 1000
particles = []
for _ in range(Np//2):
    particles.append(Particle([random.triangular(0,1,0.2), 
                              random.triangular(0,1,0.2)]))
for _ in range(Np//2, Np):
    particles.append(Particle([random.triangular(0,1,0.8), 
                              random.triangular(0,1,0.8)]))

# Build tree
root = Cell([0,0], [1,1], 0, len(particles)-1)
build_tree(particles, root, 0)

# Calculate densities
k = 32
period = root.rHigh - root.rLow
c1, c2 = [], []
positions = np.array([p.r for p in particles])

for idx, p in enumerate(particles):
    pq = PriorityQueue(k)
    periodic_search(pq, root, particles, p.r, period)
    R = math.sqrt(-pq.heap[0][0])
    
    c1.append(tophat_density(pq.heap, R))
    c2.append(monaghan_density(pq.heap, R))

# Plotting
plt.style.use('dark_background')
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

def plot_cells(cell, ax):
    ax.add_patch(plt.Rectangle(cell.rLow, *(cell.rHigh - cell.rLow), 
                              fill=False, ec='white', lw=0.5))
    if cell.pLower:
        plot_cells(cell.pLower, ax)
    if cell.pUpper:
        plot_cells(cell.pUpper, ax)

for ax in (ax1, ax2):
    plot_cells(root, ax)
    ax.set_aspect('equal')

# Plot density fields
sc1 = ax1.scatter(positions[:,0], positions[:,1], c=c1, 
                cmap='magma', s=15, edgecolor='none')
sc2 = ax2.scatter(positions[:,0], positions[:,1], c=c2, 
                cmap='magma', s=15, edgecolor='none')

# Plot reference point and neighbors
ref_point = np.array([0.5, 0.5])
pq = PriorityQueue(k)
periodic_search(pq, root, particles, ref_point, period)
R = math.sqrt(-pq.heap[0][0])

for ax in (ax1, ax2):
    ax.scatter(*ref_point, c='cyan', s=50, edgecolor='white')
    ax.scatter([p.r[0] for _, p in pq.heap], [p.r[1] for _, p in pq.heap], 
              c='red', s=30, edgecolor='white')
    ax.add_patch(plt.Circle(ref_point, R, fill=False, ec='red', lw=1))

ax1.set_title('Tophat Kernel Density', pad=20)
ax2.set_title('Monaghan Kernel Density', pad=20)

plt.colorbar(sc1, ax=ax1, label='Density')
plt.colorbar(sc2, ax=ax2, label='Density')
plt.tight_layout()
plt.show()

