/* Umbrella Sampling (US) implementation for determining the free energy of a
 * nucleus of size n in a nucleating HS system
 * This script is based on my previous NPT MC scheme
 */






//    LIBRARIES
// ===============
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <time.h>
#include <math.h>       // need to compile with `-lm` flag
#include <string.h>
#include "mt19937ar.c"





//    PREPROCESSOR CONSTANTS
// ============================
#define INVFPI 0.07957747154594766788444188168625718101 // = 1 / (4 * pi)
#define MAX_PART_CELL 40
#define MAX_NEIGHBORS 50





//    PROTOTYPES
// ================
// structures
// ----------------
typedef struct Particle Particle;
typedef struct compl compl;
typedef struct bndT bndT;
typedef struct blistT blistT;
typedef struct nlistT nlistT;
typedef struct posnb posnb;


// functions
// ----------------
int overlapCheck(Particle* one, Particle* two);
int overlapCheckGlobal(void);
int overlapCheckGlobalCL(void);
int overlapCheckCL(int n);
void sanityCheck();
void readInit(char *filename);
void writeCoords(char *filename);
int particleMove(double particleStepTune);
int volumeMove(double volumeStepTune);
double tuneStepSize(int nSuccess, int nCycles, double acceptanceRate);
double measurePF(void);
int buildCL(int usecase);
int updateSingleCL(Particle* one);
int retrieveIndex(int u, int v, int w);
int snapshot(void);
int takeSnapshot(int n);
int findClusters(void);
int findClusters_init(void);
void init_nblist(void);
void update_nblistp_sann(int p);
int compare(const void * a, const void * b);
compl* calc_order(void);
double dotprod(compl *vec1, compl *vec2, int l);
double sqr(double x);
void order(int l, bndT *bnd, compl *res1, compl *res2);
float plgndr(int l, int m, double x);
double facs(int l, int m);
double gammln(double xx);
double minpow(int m);
int* calc_conn(compl* orderp);
int calc_clusters(int* conn, compl* orderp);
int calc_clusters_init(int* conn, compl* orderp);
void cycle(void);
double Ubias(int sizeClus);
void trajectory(void);






//    GLOBAL VARIABLES
// ======================
// hard spheres system
// -----------------------
int N = 0;                              // [READ FROM IMPUT] number of particles in the system
double L = 0.0f;                        // [READ FROM INPUT] reduced box size
double LOld = 0.0f;                     // backup copy of the above
double P = 15.0f;                       // [TUNABLE] reduced pressure --- only relevant for NPT
double sigma = 1.0f;                    // unit of length = 1 = max particle size
double PF = 0.0f;                       // [CALCULATED] system packing fraction
double rho = 0.0f;                      // [CALCULATED] reduced number density
Particle *particles = NULL;             // [READ FROM INPUT] pointer to the table of particles data
Particle *particlesCopy = NULL;         // backup copy of the above
// mc
// -----------------------
const int MCCycle = 400;            // [TUNABLE] number of cycles
int randseed = 848321;                  // [TUNABLE] seed for random number generator
int trajectoryLength = 20;              // [TUNABLE] length of a trajectory, in number of cycles, must be a divider of MCCycle
double targetedAcceptanceRate = 0.30f;  // [TUNABLE] targeted acceptance rate for particle moves and volume moves
// cell lists
// -----------------------
Particle **CLTable = NULL;
double sCell1D = 0.0f;
int nCell1D = 0;
int nCell = 0;
// snapshots
// -----------------------
int nSnapshot = MCCycle / 2.0f / 20.0f; // [TUNABLE] targeted number of snapshots
int cSnapshot = 0;                      // snapshot count
char snap_filename[100];                // name of the output file for snapshots
// bop
// -----------------------
int *clustersLog = NULL;                // log of the number of clusters of a given size
int *clustersLogCopy = NULL;            // backup copy of the above 
int countLog;                           // number of independent cluster sizes logged
const double bndLength = 1.4f;          // [TUNABLE] distance cutoff for bonds, if used
const double bndLengthSq = bndLength    // square of the distance cutoff for bonds
                     * bndLength;
double bnd_cutoff = 0.7f;               // [TUNABLE] order to be called a correlated bond
int nbnd_cutoff = 4;                    // [TUNABLE] number of correlated bonds for a crystalline particle
double obnd_cutoff = 0.0f;              // [TUNABLE] order to be in the same cluster (0.0 for all touching clusters 0.9 to see defects)
double maxr2 = 0.0f;
blistT *blist;
int maxClus = 0;                        // size of the largest cluster in the system at its current state
int maxClusOld = 0;                     // backup of the size of the largest cluster in the system at its former state
// umbrella sampling
// -----------------------
double k = 0.15;                        // [TUNABLE] coupling parameter for the bias potential
int n0 = 10;                            // [TUNABLE] targeted cluster size
double potBias = 0.0f;                  // value of the bias potential associated with the size of the largest cluster in the current state of the system
double potBiasOld = 0.0f;               // value of the bias potential associated with the size of the largest cluster in the previous state of the system
// input/output
// -----------------------
char* init_filename;                    // name of the initial configuration input file
char lastsnap_filename[200];            // name of the last snapshot output file
char movie_filename[200];               // name of the movie output file
char clusterlog_filename[200];          // name of the clusters log output file
int makeamovie = 1;                     // choice to make a movie of snapshots taken after each trajectory



//    STRUCTURES
// ================
struct Particle
{
        /*\
         *  Structure:  Particle
         *  --------------------
         *  Particle in the system, implemented for the use of cell lists (CL)
        \*/

	double x;	        // reduced x coordinate 
	double y;	        // reduced y coordinate
        double z;               // reduced z coordinate
	double r;	        // reduced radius --- redundant with sigma ... might remove later on
        char type;              // particle type --- for visualization purposes only
        Particle *next;         // pointer to the next particle in CL
        Particle *prev;         // pointer to the previous particle in CL
        int CLindex;            // index of CL in which particle is
        int index;              // index of the particle, in [0;N-1]
};



struct compl
{
        /*\
         *  Structure:  compl
         *  -----------------
         *  Complex number
        \*/

        double re;              // real part
        double im;              // imaginary part
};



struct bndT
{
        double nz;
        double si;
        double co;
        int n;
};



struct blistT
{
        int n;
        bndT * bnd;
};



struct nlistT
{
        int n;
        int * nb;
};



struct posnb
{
        Particle* part;
        double dist;
        double dx;
        double dy;
        double dz;
};






//    FUNCTIONS
// ===============
int overlapCheck(Particle* one, Particle* two)
{
        /*\
         *  Function:    overlapCheck
         *  -------------------------
         *  Checks for overlap between two particles
         *
         *  return:      1 if there is overlap, 0 if not
        \*/

        double dx = one->x - two->x,
               dy = one->y - two->y,
               dz = one->z - two->z;
        dx = dx - L * rint(dx / L);
        dy = dy - L * rint(dy / L);
        dz = dz - L * rint(dz / L);
        if ((dx * dx + dy * dy + dz * dz) < ((one->r + two->r) * (one->r + two->r)))
                return 1;       
        return 0;
}



int overlapCheckGlobal(void)
{
        /*\
         *  Function:     overlapCheckGlobal
         *  -----------------------------
         *  Checks for overlap in the whole system, i.e. for N(N-1)/2 pairs
         *
         *  return:     1 if there is overlap, 0 if not
        \*/

        int overlap = 0, n = 1, m = 0;
        while ((overlap == 0) && (n < N))
        {
                while ((overlap == 0) && (m < N))
                {
                        overlap = (n != m) * overlapCheck(particles+n, particles+m);
                        if (overlap)
                                printf("[OVERLAP] #%d & #%d\n", (particles+n)->index, (particles+m)->index);
                        m++;
                }
                m = n;
                n++;
        }
        return overlap;
}



int overlapCheckGlobalCL(void)
{
        /*\
         *  Function:     overlapCheckGlobalCL
         *  ----------------------------------
         *  Checks for overlap in the whole system using CLs; also checks pairs
         *  only once
         *
         *  return:     1 if there is overlap, 0 if not
        \*/

        for (int i = 0; i < N; i++)
        {
                if (overlapCheckCL(i))
                        return 1;
        }
        return 0;
}



int overlapCheckCL(int n)
{
        /*\
         *  Function:    overlapCheckCL
         *  ---------------------------
         *  Checks for overlap between a given particle and all other particles
         *  that can be found in neighbouring cells. Relies on the Cell Lists
         *  implementation.
         *
         *  n:           index of the particle to check overlap for
         *
         *  return:      1 if there's overlap, 0 if not
        \*/
        
        int a = particles[n].x / sCell1D;
        int b = particles[n].y / sCell1D;
        int c = particles[n].z / sCell1D;
        int index;
        Particle *current = NULL;
        for (int i = a-2; i < a+3; i++)
        {        
                for (int j = b-2; j < b+3; j++)
                {
                        for (int k = c-2; k < c+3; k++)
                        {
                                index = retrieveIndex(i,j,k);
                                current = CLTable[index];
                                while (current != NULL)
                                {
                                        if ((particles+n != current) && overlapCheck(particles+n, current))
                                                return 1;
                                        current = current->next;
                                }
                        }
                }
        }
        return 0;
}



void sanityCheck(void)
{
        /*\
         *  Function:    sanityCheck
         *  ------------------------
         *  Ends program execution if there is overlap in the system
        \*/
        
        if (overlapCheckGlobal())
        {
                writeCoords(movie_filename);
                printf("Something is VERY wrong... There's overlap...\n");
                exit(0);
        }
}



void readInit(char *filename)
{
        /*\
         *  Function:    readInit
         *  ---------------------
         *  Initializes the table of data of all N particles by reading
         *  it from a user supplied .sph file; particles coordinates
         *  are written in reduced units; also retrieves **cubic** box size
         *  NB: Reading is done only for files with typesetting indicated above
         *  NB: It is assumed that the supplied init file is written in standard units so
         *  conversion to reduced units is applied
         *
         *  *filename:   pointer to the name of the initialization .txt file
        \*/

        FILE *initfile = NULL;
        int i = 0;
        double rTemp = 0.0f;
        initfile = fopen(filename, "r");
        if (initfile != NULL)
        {
                // Read value of N and allocate memory accordingly
                if (fscanf(initfile, "%*c%d%*c", &N) != 1)
                {
                        printf("\nERROR --- Wrong input file\n\n");
                        exit(0);
                }
                particles = malloc(N * sizeof(*particles));
                particlesCopy = malloc(N * sizeof(*particlesCopy));

                // Read value of box size
                if (fscanf(initfile, "%lf %*f %*f%*c", &L) != 1)
                {
                        printf("\nERROR --- Wrong input file\n\n");
                        exit(0);
                }
                // Populate table of particles data in standard units
                for (i = 0; i < N; i++)
                {
                        if (fscanf  ( initfile,
                                  "%c %lf %lf %lf %lf%*c",
                                  &particles[i].type,
                                  &particles[i].x,
                                  &particles[i].y,
                                  &particles[i].z,
                                  &particles[i].r
                                )
                            != 5
                           )
                        {
                                printf("\nERROR --- Wrong input file\n\n");
                                exit(0);
                        }
                        // Find larger radius to rescale everything to obtain
                        // the expected sigma = 1 for the (larger) particle size
                        if (particles[i].r > rTemp)
                                rTemp = particles[i].r;
                }
                fclose(initfile);
                // Rescale everything w.r.t. larger particle size to obtain the
                // expected max sigma = 1, also converts to reduced units
                L /= (2.0f * rTemp);
                LOld = L;
                for (i = 0; i < N; i++)
                {
                        particles[i].x /= (2.0f * rTemp);
                        particles[i].y /= (2.0f * rTemp);
                        particles[i].z /= (2.0f * rTemp);
                        particles[i].r /= (2.0f * rTemp);
                        // Also attributes each particle a unique index
                        particles[i].index = i;
                        // Apply PBC to the read configuration --- necessary
                        // for inputs generated by Frank's EDMD where particles
                        // can lie outside the box
                        particles[i].x = fmod(particles[i].x + 2 * L, L);
                        particles[i].y = fmod(particles[i].y + 2 * L, L);
                        particles[i].z = fmod(particles[i].z + 2 * L, L);
                }
        }
        // Make a buffer copy of the particles table
        for (int j = 0; j < N; j++)
        {
                particlesCopy[j].x = particles[j].x;
                particlesCopy[j].y = particles[j].y;
                particlesCopy[j].z = particles[j].z;
                particlesCopy[j].r = particles[j].r;
                particlesCopy[j].type = particles[j].type;
                particlesCopy[j].next = particles[j].next;              // can cause problems
                particlesCopy[j].prev = particles[j].prev;              // can cause problems
                particlesCopy[j].CLindex = particles[j].CLindex;
                particlesCopy[j].index = particles[j].index;
        }

}



void writeCoords(char *filename)
{
        /*\
         *  Function:    writeCoords
         *  -------------------------
         *  Writes the position, radius, and type data for all N
         *  particles in a .sph file; standard units are used
         *  NB: Writing is done using the typesetting indicated above
         *
         *  *filename:   pointer to the name of the output .sph file
        \*/

        FILE *outfile = NULL;
        int i = 0;
        double boxSize = L * sigma;
        outfile = fopen(filename, "a");
        if (outfile != NULL)
        {
                fprintf(outfile, "&%d\n", N);
                fprintf(outfile, "%.12lf %.12lf %.12lf\n", boxSize, boxSize, boxSize);
                for (i = 0; i < N; i++)
                        fprintf(outfile,
                                "%c %.12lf %.12lf %.12lf %.12lf\n",
                                particles[i].type,
                                particles[i].x * sigma,
                                particles[i].y * sigma,
                                particles[i].z * sigma,
                                particles[i].r * sigma
                               );
                fclose(outfile);
        }
}



int particleMove(double particleStepTune)
{
        /*\
         *  Function: particleMove                
         *  ---------------------
         *  Tries to move a randomly selected particle in space by a small amount
         *  and checks for overlap
         *  If there is no overlap, overwrites the previous particle position
         *  with the new one
         *
         *  particleStepTune:    tuning factor for particle step size
         *
         *  return:              0 if attempt failed, 1 if attempt succeeded
        \*/ 
        
        // Generate random displacements in x, y, z
        double delta[3] = {(genrand() - 0.5) / sigma * (particleStepTune),
                           (genrand() - 0.5) / sigma * (particleStepTune),
                           (genrand() - 0.5) / sigma * (particleStepTune)
                          };
        // Randomly select a particle in the system
        int n = (int) (genrand() * N);
        Particle *selected = particles+n;
        // Set up a backup of the moved particle in case of overlap
        Particle bufferParticle = {     selected->x,
                                        selected->y,
                                        selected->z
                                  };
        // Move randomly selected particle by abiding to PBC and nearest image
        // convention
        selected->x = fmod((selected->x + delta[0]) + 2 * L, L);
        selected->y = fmod((selected->y + delta[1]) + 2 * L, L);
        selected->z = fmod((selected->z + delta[2]) + 2 * L, L);
        // Update Cell List
        updateSingleCL(particles+n);
        // Check for overlap and move back particle if need be
        if (overlapCheckCL(n))
        {
                selected->x = bufferParticle.x;
                selected->y = bufferParticle.y;
                selected->z = bufferParticle.z;
                updateSingleCL(particles+n);
                return 0;
        }
        else
                return 1;                
}



int volumeMove(double volumeStepTune)
{
        /*\
         *  Function:    volumeMove
         *  -----------------------
         *  Attempts to change the volume of the simulation box by a small
         *  amount according to the known acceptance rule
         * 
         *  volumeStepTune:      tuning factor for volume step size
         *
         *  return:              0 if attempt failed, 1 if attempt succeeded
        \*/
        
        // Generate random volume change, define new adequate quantities
        double  vol = L * L * L,
                delta[2] = {(genrand() - 0.5) / (sigma) * (volumeStepTune),
                            genrand()
                           },
                newVol = vol + delta[0],                  
                newL = cbrt(newVol),                      
                ratio = newL / L,
                rule = exp(- P * (newVol - vol) + N * log(newVol / vol));
        // Rejects move if rule is not respected
        if (delta[1] > rule)
                return 0;
        // Scale the system according to the volume change
        for (int i = 0; i < N; i++)
        {
                particles[i].x *= ratio;
                particles[i].y *= ratio;
                particles[i].z *= ratio;
        }
        L *= ratio;
        // Checks for overlap in case of box shrinking and rescale cells
        // accordingly in each case
        if (ratio < 1)
        { // Box shrinks
                //int overlap = 0;
                sCell1D *= ratio;
                if (sCell1D < sigma)
                // Cell size does not guarantee check for all possible
                // overlaps
                // Need to reduce cell size, thus number of cells
                        buildCL(3);
                else
                // Cell size is sufficient to check for all possible
                // overlaps
                // Only need to reduce cell size accordingly
                        buildCL(2);

                if (overlapCheckGlobalCL())
                { // Overlap
                        for (int i = 0; i < N; i++)
                        {
                                particles[i].x /= ratio;
                                particles[i].y /= ratio;
                                particles[i].z /= ratio;
                        }
                        L /= ratio;
                        sCell1D /= ratio;
                        buildCL(3);
                        return 0;
                }
                else
                { // No overlap
                        return 1;
                }
        }
        else
        { // Box expands = no overlap
                //sCell1D *= ratio;
                //if (sCell1D > (2.0f * sigma))
                // Cell size is such that we check for overlap where there
                // cannot be any
                // Need to reduce cell size, thus increase number of cells,
                // thus free existing CL table and re-allocate enough memory
                //        buildCL(1);
                //else
                // Cell size is necessary (?) to check for all possible overlaps
                // Only need to increase cell size accordingly
                //        buildCL(2);
                buildCL(3);
                return 1;
        }
}



double tuneStepSize(int nSuccess, int nCycles, double acceptanceRate)
{
        /*\
         *  Function:   tuneStepSize
         *  -------------------------
         *  Tunes the step size for volume or particle moves by +-2% depending
         *  on whether the success rate is above or below the specified
         *  acceptance rate
         *
         *  NSuccess:           number of successes over the last nCycles cyles
         *  acceptanceRate:     targeted acceptance rate
         *  nCycles:            number of cycles to average the number of successes over
         *
        \*/

        if (nSuccess - (int) (acceptanceRate * nCycles) > 0)
                return 1.02;
        else
                return 0.98; 
}



double measurePF(void)
{
        /*\
         *  Function:   measurePF
         *  -------------------
         *  Computes the reduced packing fraction for the given parameters
         *
         *  return:     reduced packing fraction
        \*/
        
        double PF = (double) N * sigma * sigma * sigma * M_PI / (6.0f * L * L * L);
        return PF;
}



int updateSingleCL(Particle *one)
{
        int index = one->CLindex;
        int u = one->x / sCell1D,
            v = one->y / sCell1D,
            w = one->z / sCell1D;
        
        // Remove particle from old CL
        if (CLTable[index] == one)
                CLTable[index] = one->next;
        if (one->next != NULL)
                one->next->prev = one->prev;
        if (one->prev != NULL)
                one->prev->next = one->next;
        
        // Fetches the index of the new CL
        index = retrieveIndex(u, v, w);
        one->CLindex = index;
        
        // Place particle on top of new CL
        one->prev = NULL;
        one->next = CLTable[index];
        if (CLTable[index] != NULL)
                CLTable[index]->prev = one;
        CLTable[index] = one;
        
        return 0;
}



int buildCL(int usecase)
{
        /*\
         *  Function:   buildCL
         *  --------------------
         *  Builds, updates, and expands cell lists and their table according to
         *  the new state of the system
         *
         *  usecase:    choice of what operation to perform on CL among:
         *              + 1:    initialization of the CL table and CLs
         *                      contents
         *              + 2:    update of the CLs contents only, no
         *                      redefinition of the table, its size, or the
         *                      size of the cells
         *              + 3:    redefinition of the size of the cell, which
         *                      decreases, and update of the CLs contents; less
         *                      of the previously allocated memory is used in
         *                      this particular case
         *
         * r eturn:     0 for normal termination
        \*/

        switch (usecase)
        {
                case 1:         // former initCL()
                // sCell goes back to ~sigma
                // nCell goes up, i.e. more cells
                // more memory is allocated to accomodate for increase in number of cells
                {
                        sCell1D = L / ((int) (L / sigma));
                        nCell1D = (int) (L / sCell1D);
                        nCell = nCell1D * nCell1D * nCell1D;
                        free(CLTable);
                        CLTable = malloc(nCell * sizeof(**CLTable));
                        if (CLTable == NULL)
                                exit(EXIT_FAILURE);
                        for (int i = 0; i < nCell; i++)
                                CLTable[i] = NULL;
                        break;
                }
                case 2:         // former resizeCL() 
                // sigma < sCell < 2*sigma
                // same number of cells
                // no more memory allocated
                {
                        for (int i = 0; i < nCell; i++)
                                CLTable[i] = NULL;
                        break;
                }
                case 3:         // former updateCL()
                // sCell goes back to ~sigma
                // nCell goes down, i.e. less cells
                // no more memory allocated, we juste use less of it
                {
                        for (int i = 0; i < nCell; i++)
                                CLTable[i] = NULL;
                        sCell1D = L / ((int) (L / sigma));
                        nCell1D = (int) (L / sCell1D);
                        nCell = nCell1D * nCell1D * nCell1D;
                        break;
                }
                default:
                        return 1;
        }
        for (int i = 0; i < N; i++)
        {
                Particle *one = particles+i;
                int u = one->x / sCell1D;
                int v = one->y / sCell1D;
                int w = one->z / sCell1D;
                int index = retrieveIndex(u, v, w);
                one->CLindex = index;
                one->prev = NULL;
                one->next = CLTable[index];
                if (CLTable[index] != NULL)
                        CLTable[index]->prev = one;
                CLTable[index] = one;
        }
        return 0;
}



int retrieveIndex(int u, int v, int w)
{
        /*\
         *  Function:   retrieveIndex
         *  --------------------------
         *  Retrieves index of the cell (from the table of CLs) to which a
         *  particle belongs based on the coordinates of the CL in the 3D
         *  coordinate system of the CL table; corresponds of a flattening
         *  operation of the CL index from a 3D view to standard 1D view
         *
         *  u:          x-coordinate of the cell
         *  v:          y-coordinate of the cell
         *  w:          z-coordinate of the cell
         *
         *  return:     index of the CL in the table of CLs
        \*/

        int index = 0;
        index = (w + 2 * nCell1D) % nCell1D * nCell1D * nCell1D
                + (v + 2 * nCell1D) % nCell1D * nCell1D
                + (u + 2 * nCell1D) % nCell1D;
        return index;

}



int snapshot(void)
{
        FILE *writeto = NULL;
         // coordinates
        //writeCoords("snap_coords.sph");
        writeCoords(snap_filename);
        // density, if needed
        writeto = fopen("snap_density.txt", "a");
        fprintf(writeto, "%.12lf\n", 6.0f * measurePF() / (M_PI * sigma * sigma * sigma));
        fclose(writeto);
        return 0;
}



int takeSnapshot(int n)
{
        if ((cSnapshot < n) && (genrand() < (2 * (double) (n) / (double) (MCCycle))))
        {
                snapshot(); // be careful, need to change writing method so it concatenates properly
                cSnapshot += 1;
        }
        return 0;
}



int findClusters(void)
{
        /*\
         *  Function:   findClusters
         *  ------------------------
         *  Wrapper for the identification of clusters from the local BOO
         *  analysis
         *  Updates the table accounting for the number of clusters of a given
         *  size
         *
         *  return:     size of the largest cluster in the current state of the
         *              system
        \*/

        compl* order = NULL;
        int* connections = NULL;
        int maxclus;

        init_nblist();
        order = calc_order();
        connections = calc_conn(order);
        maxclus = calc_clusters(connections, order);
        
        free(order);
        free(connections);
        for (int i = 0; i < N; i++)
                free(blist[i].bnd);
        free(blist);

        return maxclus;
}



int findClusters_init(void)
{
        /*\
         *  Function:   findClusters
         *  ------------------------
         *  Wrapper for the identification of clusters from the local BOO
         *  analysis
         *  Updates the table accounting for the number of clusters of a given
         *  size
         *
         *  return:     size of the largest cluster in the current state of the
         *              system
        \*/

        compl* order = NULL;
        int* connections = NULL;
        int maxclus;

        init_nblist();
        order = calc_order();
        connections = calc_conn(order);
        maxclus = calc_clusters_init(connections, order);
        
        free(order);
        free(connections);
        for (int i = 0; i < N; i++)
                free(blist[i].bnd);
        free(blist);

        return maxclus;
}



void init_nblist(void)
{
        /*\
         *  Function:   init_nblist
         *  -----------------------
         *  initializes the list of neighbors for each particle in the current
         *  state of the system
        \*/

        blist = malloc(N * sizeof(*blist));

        for (int p = N-1; p > -1; p--)
        {
                blist[p].bnd = malloc(MAX_NEIGHBORS * sizeof(bndT));
                update_nblistp_sann(p);
                //printf("Updated neighbors list of particle #%d: has %d nearest neighbours\n", p, blist[p].n);
        }
}



void update_nblistp_sann(int p)
{
        // local variables:
        int id;
        bndT* bnd;
        double d,
               dxy,
               dx,
               dy,
               dz,
               di;

        // finding to which cell belongs the current particle
        int a = particles[p].x / sCell1D;
        int b = particles[p].y / sCell1D;
        int c = particles[p].z / sCell1D;
        
        // initializes number of nearest neighbours for particle `p` to 0
        blist[p].n = 0;

        // defining another buffer particle
        int cellIndex;
        Particle *current = NULL;

        // saving info about neigbouring particles
        posnb neighbours[N];
        int numposnb = 0;

        // cycle over all particles in all neighboring CLs
        for (int i = a-2; i < a+3; i++)
        {
                for (int j = b-2; j < b+3; j++)
                {
                        for (int k = c-2; k < c+3; k++)
                        {
                                cellIndex = retrieveIndex(i,j,k);
                                current = CLTable[cellIndex];
                                while (current != NULL && numposnb <= N)
                                {
                                        int currentID = current->index;
                                        if (p == currentID)
                                        { // escapes if both particles are the same
                                                current = current->next;
                                                continue;
                                        }

                                        // computing distances
                                        dx = particles[p].x - current->x;
                                        dy = particles[p].y - current->y;
                                        dz = particles[p].z - current->z;

                                        // applying NIC
                                        dx = dx - L * rint(dx / L);
                                        dy = dy - L * rint(dy / L);
                                        dz = dz - L * rint(dz / L);

                                        // computes squared distance between two particles
                                        d = dx * dx + dy * dy + dz * dz;

                                        // saves neighbour data
	                                neighbours[numposnb].part = current;
                                        neighbours[numposnb].dist = sqrt(d);
	                                
                                        neighbours[numposnb].dx = dx;                   
	                                neighbours[numposnb].dy = dy;
	                                neighbours[numposnb].dz = dz;
	                                
                                        // continues with the next neighbouring particle
                                        numposnb ++;
	                                current = current->next;
                                }
                        }
                }
        }
        //printf("Particle #%d has %d neighbours\n", p, numposnb);

        
        qsort(neighbours, numposnb, sizeof(posnb), compare); // sorts the array 'neighbours' of neighboring particles
        
        // FINDING NEAREST NEIGHBORS USING CUTOFF RADIUS
        // ---------------------------------------------
        /*
        for (int i = 0; i < MAX_NEIGHBORS; i++) 
        {
                posnb* nb = neighbours+i;
                if (nb->dist < bndLength)
                {
                        bnd = &(blist[p].bnd[blist[p].n]);
                        blist[p].n++;
                        id = nb->part->index;
                        bnd->n = id;
                        if (id > p)
                        {
                                dx = nb->dx;
                                dy = nb->dy;
                                dz = nb->dz;
                                di = 1.0 / (nb->dist);
                                bnd->nz = dz * di;
                                dxy = 1.0 / sqrt(dx * dx + dy * dy);
                                bnd->si = dy * dxy;
                                bnd->co = dx * dxy;
                        }
                }
        }
        */

        // FINDING NEAREST NEIGHBORS USING SANN
        // ------------------------------------
        int m = 3; 
        int done = 0;
        while (!done)
        {
                double rim = 0;
                for (int i = 0; i < m; i++)
                        rim += neighbours[i].dist / (m - 2);
                if (rim > neighbours[m].dist)
                        m++;
                else 
                        done = 1;
                if (m > numposnb)
                { 
                        printf("NN algorithm did not converge!\n");
                        m--;
                        done = 1;
                }
        }

        // can redefine m at this point if we want to use only, e.g., the first 12 NN
        // m = 12;
        for (int i = 0; i < m; i++)
        {
                posnb* nb = neighbours+i;
                id = nb->part->index;
                dx = nb->dx;
                dy = nb->dy;
                dz = nb->dz;
                di = 1.0 / (nb->dist);
                bnd = &(blist[p].bnd[blist[p].n]);
                blist[p].n++;
                bnd->n = id;
                bnd->nz = dz * di;
                dxy = 1.0 / sqrt(dx * dx + dy * dy);
                bnd->si = dy * dxy;
                bnd->co = dx * dxy;
        }


}



int compare(const void * a, const void * b)
{
        /*\
         *  Function:   compare
         *  -------------------
         *  For use in qsort() function from stdlib.h
         *
         *  return:      1      if d > 0
         *               0      if d = 0
         *              -1      if d < 0
        \*/

        double d = ((posnb*) a)->dist - ((posnb*) b)->dist;

        return ((0 < d) - (d < 0));
}



compl* calc_order(void)
{
        /*\
         *  Function:   calc_order
         *  ----------------------
         *  Computes the normalized inner product of the local bond-order
         *  parameter (BOP)
         *
         *  return:     table of size 13*N of complex values of the inner
         *              product of the local BOP
        \*/

        compl *q1,
              *q2;
        const int l = 6;
        double temp;
        compl *orderp= (compl*) malloc(sizeof(compl) * N * (l * 2 + 1));
        memset(orderp, (int) 0.0, sizeof(compl) * N * (l * 2 + 1));
        
        for(int i = 0; i < N; i++)
        {
                q1 = (orderp+i * (2 * l + 1) + l);
                for(int j = 0; j < blist[i].n; j++)
                {
                        if(blist[i].bnd[j].n > i) 
                        {
                                q2 = (orderp + blist[i].bnd[j].n * (2 * l + 1) + l);
                                order(l, &(blist[i].bnd[j]), q1, q2);
                        }  
                }
        }  
    
        // normalize qs
        for(int i = 0; i < N; i++)
        {
                temp = sqrt(dotprod(orderp + i * (2 * l + 1), orderp + i * (2 * l + 1), l));
                temp = 1.0 / temp;
                for(int m = -l ; m <= l; m++) 
                {
                        (*(orderp+i * (2 * l + 1) + m + l)).re *= temp;
                        (*(orderp+i * (2 * l + 1) + m + l)).im *= temp;
                }
        }

        return orderp;
}



double dotprod(compl *vec1, compl *vec2, int l)
{
        /*\
         *  Function:   dotprod
         *  -------------------
         *  Computes the dot product of two complex vectors
         *
         *  *vec1:      pointer to a first complex vector
         *  *vec2:      pointer to a second complex vector
         *
         *  return:     dot product of two complex vectors
        \*/

        double res = 0.0f;
        for(int m = -l; m <= l; m++)
                res += (*(vec1 + m + l)).re * (*(vec2 + m + l)).re
                       + (*(vec1 + m + l)).im * (*(vec2 + m + l)).im;

        return res;
}



double sqr(double x)
{
        /*\
         *  Function:   sqr
         *  ---------------
         *  Computes the square of a double
         *
         *  x:          double to compute the square of
         *
         *  return:     x squared
        \*/

        return x*x;
}



void order(int l, bndT *bnd, compl *res1, compl *res2)
{
        /*\
         *  Function:   order
         *  -----------------
         *  
        \*/
        double fc,
               p,
               f,
               s,
               r,
               sp,
               spp,
               c,
               cp,
               cpp;
        double z;
        int m = 0;
        z = bnd->nz;

        // Not tested for l! = 6 

        p = plgndr(l,0,z);
        fc = facs(l,0);
        f = sqrt( (2*l+1) * INVFPI * fc );
        r = p*f;
        (res1+0)->re += r;
        (res1+0)->im += 0;
        (res2+0)->re += r * minpow(l);
        (res2+0)->im += 0;
        s=0;
        sp=0;
        c=1;
        cp=0;

        for(m = 1; m <= l; m++)
        {
                // positive m
                p = plgndr(l,m,z);
                fc = facs(l,m);
                f = sqrt((2 * l + 1) * INVFPI * fc);
                r = p * f;
                cpp = cp;
                cp = c;
                if(m == 1)
                        c = bnd->co;
                else
                        c = 2.0 * bnd->co * cp - cpp; //some cosine tricks

                spp = sp;
                sp = s;
                if(m == 1)
                        s = bnd->si;
                else
                        s = 2.0 * bnd->co * sp - spp; //some sine tricks

                (res1+m)->re += r*c;
                (res1+m)->im += r*s;
                (res2+m)->re += r*c;
                (res2+m)->im += r*s;

                //negative m
                r *= minpow(m);
                (res1-m)->re += r * c;
                (res1-m)->im += - r * s;
                (res2-m)->re += r * c;
                (res2-m)->im += - r * s;
        }
}



float plgndr(int l, int m, double x)
{
        /*\
         *  Function:   plgndr
         *  ------------------
        \*/

        // variables
        double fact,
               pll = 0.0,
               pmm,
               pmmp1,
               somx2;
        int i,
            ll;

        // checks for normal computation of Legendre polynoms
        if (m < 0 || m > l || fabs(x) > 1.0)
                printf("Bad arguments in routine plgndr %i %i %f\n", l, m, fabs(x));
        
        pmm = 1.0;
        
        if (m > 0)
        {
                somx2 = sqrt((1.0 - x) * (1.0 + x));
                fact = 1.0;
                for (i = 1; i <= m; i++)
                {
                        pmm *= -fact * somx2;
                        fact += 2.0;
                }
        }

        if (l == m)
                return pmm;
        else
        { 
                pmmp1 = x * (2 * m + 1) * pmm;
                if (l == (m + 1))
                        return pmmp1;
                else
                {
                        for (ll = m + 2; ll <= l; ll++)
                        {
                                pll = (x * (2 * ll - 1) * pmmp1 - (ll + m - 1) * pmm) / (ll - m);
                                pmm = pmmp1;
                                pmmp1 = pll;
                        }
                        return pll;
                }
        }
}



double facs(int l, int m)
{
        /*\
         *  Function:   facs
         *  ----------------
        \*/

        static double *fac_table[14] = {NULL,
                                        NULL,
                                        NULL,
                                        NULL,
                                        NULL,
                                        NULL,
                                        NULL,
                                        NULL,
                                        NULL,
                                        NULL,
                                        NULL,
                                        NULL,
                                        NULL,
                                        NULL
                                        }; //max l=10
        int a, b;
        if(fac_table[l] == NULL)
        {
                fac_table[l] = malloc(sizeof(double) * (2 * l + 1));
                for(a = 0; a < 2 * l + 1; a++)
                {
                        b = a - l;
                        fac_table[l][a]= exp(gammln(l - b + 1) - gammln(l + b + 1));
                }
        }
        return fac_table[l][m+l];
}



double gammln(double xx)
{
        /*\
         *  Function:   gammln
         *  ------------------
        \*/

        double x,
               y,
               tmp,
               ser;
        static double cof[6] = {76.18009172947146,
                                -86.50532032941677,
                                24.01409824083091,
                                -1.231739572450155,
                                0.1208650973866179e-2,
                                -0.5395239384953e-5
                                };
        int j;
        y = x = xx;
        tmp = x + 5.5;
        tmp -= (x + 0.5) * log(tmp);
        ser = 1.000000000190015;
        for (j = 0; j <= 5; j++) 
                ser += cof[j] / ++y;
        return -tmp + log(2.5066282746310005 * ser / x);
}



double minpow(int m)
{
        /*\
         *  Function:   minpow
         *  ------------------
        \*/

        if((m & 1) == 1)
                return -1.0;
        else
                return 1.0;
}



int* calc_conn(compl* orderp) 
{
        /*\
         *  Function:   calc_conn
         *  ---------------------
        \*/

        int z;
        const int l = 6;
        int *conn = malloc(sizeof(int) * N);
        for(int i = 0; i < N; i++)
        {
                z = 0;
                for(int j = 0; j < blist[i].n; j++)
                {
                        if(dotprod(orderp + i * (2 * l + 1), orderp + blist[i].bnd[j].n * (2 * l + 1), l) > bnd_cutoff)
                        { 
                                z++;
                        }
                }
                if(z >= nbnd_cutoff)
                {
                        conn[i]=1;
                } 
                else 
                {
                        conn[i]=0;
                }
        }

        return conn;
}



int calc_clusters(int* conn, compl* orderp)
{
        /*\
         *  Function:   calc_clusters
         *  -------------------------
         *
         *  return:     size of the largest cluster in the current state of the
         *              system
        \*/

        // variables
        // ---------
        int* cluss = malloc(sizeof(int) * N);
        int* size= malloc(sizeof(int) * N);
        int cs,
            cn = 1,
            big = 0;
        //int bc = -1;
        const int l = 6;
        
        // functions
        // ---------
        void setcluss(int pn)
        {
                cluss[pn] = cn;
                for(int jj = 0; jj < blist[pn].n; jj++)
                {
                        int tmp = blist[pn].bnd[jj].n;
                        if(conn[tmp] != 0
                           && cluss[tmp] == 0
                           && dotprod(orderp + pn * (2 * l + 1), orderp + tmp * (2 * l + 1), l) > obnd_cutoff
                           )
                        {
                        //obnd_cutoff = 0.9 gives nice results,
                        //obnd_cutoff =  0.6 gives all touching nuclei as one big nuclei
                                cs++;
                                setcluss(tmp);
                        }  
                }
        }

        // body
        // ----
        for(int i = 0; i < N; i++)
        {
                cluss[i] = 0;           
                size[i] = 0;
        }  

        for(int i = 0; i < N; i++)
        {
                cs = 0;
                if(conn[i] == 1 && cluss[i] == 0) 
                {
                        cs++;
                        setcluss(i); 
                        size[cn] = cs;
                        if(cs > big)
                        { 
                                big = cs;
                                //bc = cn;
                        }
                        cn++;
                }
        }
  
        int tcs = 0;                            // total number of particles in a cluster
        countLog = 0;
        for(int i = 0 ; i < cn ; i++)           // loop over all clusters
        {
                tcs += size[i];
                if (size[i] != 0)
                {
                        //printf("N = %d\n", N);
                        //printf("size[%d] = %d\n", i, size[i]);
                        clustersLog[size[i]]++;
                        if (clustersLog[size[i]] == 1)
                                countLog++;
                }
        }

        //printf("%i clusters, Max clustersize %i, percentage of crystalline particles %f\n", cn-1, big, tcs / (double) N);

        free(cluss);
        free(size);
        return big;
}



int calc_clusters_init(int* conn, compl* orderp)
{
        /*\
         *  Function:   calc_clusters
         *  -------------------------
         *
         *  return:     size of the largest cluster in the current state of the
         *              system
        \*/

        // variables
        // ---------
        int* cluss = malloc(sizeof(int) * N);
        int* size= malloc(sizeof(int) * N);
        int cs,
            cn = 1,
            big = 0;
        //int bc = -1;
        const int l = 6;
        
        // functions
        // ---------
        void setcluss(int pn)
        {
                cluss[pn] = cn;
                for(int jj = 0; jj < blist[pn].n; jj++)
                {
                        int tmp = blist[pn].bnd[jj].n;
                        if(conn[tmp] != 0
                           && cluss[tmp] == 0
                           && dotprod(orderp + pn * (2 * l + 1), orderp + tmp * (2 * l + 1), l) > obnd_cutoff
                           )
                        {
                        //obnd_cutoff = 0.9 gives nice results,
                        //obnd_cutoff =  0.6 gives all touching nuclei as one big nuclei
                                cs++;
                                setcluss(tmp);
                        }  
                }
        }

        // body
        // ----
        for(int i = 0; i < N; i++)
        {
                cluss[i] = 0;           
                size[i] = 0;
        }  

        for(int i = 0; i < N; i++)
        {
                cs = 0;
                if(conn[i] == 1 && cluss[i] == 0) 
                {
                        cs++;
                        setcluss(i); 
                        size[cn] = cs;
                        if(cs > big)
                        { 
                                big = cs;
                                //bc = cn;
                        }
                        cn++;
                }
        }
  
        int tcs = 0;                            // total number of particles in a cluster
        countLog = 0;
        for(int i = 0 ; i < cn ; i++)           // loop over all clusters
        {
                tcs += size[i];
                if (size[i] != 0)
                {
                        //printf("N = %d\n", N);
                        //printf("size[%d] = %d\n", i, size[i]);
                        clustersLog[size[i]]++;
                        if (clustersLog[size[i]] == 1)
                                countLog++;
                }
        }
        for (int i = 0; i < N; i++)
        {
                clustersLogCopy[i] = clustersLog[i];
        }

        //printf("%i clusters, Max clustersize %i, percentage of crystalline particles %f\n", cn-1, big, tcs / (double) N);

        free(cluss);
        free(size);
        return big;
}



void cycle(void)
{
        /*\
         *  Function:   cycle
         *  -----------------
         *  Performs a Monte Carlo cycle in {N,P,T}
        \*/

        double fixedStepSize = 1.0f;
        double proba = 0.0f;

        for (int j = 0; j < N; j++)
        {
                proba = genrand();
                if (proba < (1.0f / (float) (N + 1)))
                        volumeMove(fixedStepSize);
                else
                        particleMove(fixedStepSize);
        }

}



double Ubias(int sizeClus)
{
        /*\
         *  Function:   Ubias
         *  -----------------
         *  Computes the bias potential for a given cluster size
         *
         *  sizeClus:   cluster size to compute the bias potential for
         *
         *  return:     bias potential for the aforementioned cluster size
        \*/

        double w = 0.0f;
        w = 0.5f * k * ((double) sizeClus - (double) n0) * ((double) sizeClus - (double) n0);
        return w;

}


void trajectory(void)
{
        /*\
         *  Function: trajectory
         *  --------------------
         *  Performs a trajectory consisting of trajectoryLength MC cycles
        \*/
       
        double proba = genrand();
        double rule = 0.0f;
        clustersLog = malloc(N * sizeof(*clustersLog));
        memset(clustersLog, (int) 0.0, sizeof(*clustersLog) * N);
        //clustersLogCopy = malloc(N * sizeof(*clustersLogCopy));

        // computes new trajectory
        for (int i = 0; i < trajectoryLength; i++)
        {
                cycle();
                sanityCheck();          // SOMETHING IS HAPPENING HERE (SAME AS BELOW)
        }
        // compute table of clusters, find largest cluster & need to keep track of it
        maxClus = findClusters();
        potBias = Ubias(maxClus);
        // accept or reject trajectory based on biased potential
        rule = exp(potBiasOld - potBias);
        if (proba > rule)
        {// reject trajectory
                //printf("trajectory rejected!\n");
                // copies configuration backup
                for (int j = 0; j < N; j++)
                {
                        particles[j].x = particlesCopy[j].x;
                        particles[j].y = particlesCopy[j].y;
                        particles[j].z = particlesCopy[j].z;
                        particles[j].r = particlesCopy[j].r;
                        particles[j].type = particlesCopy[j].type;
                        particles[j].next = particlesCopy[j].next;              
                        particles[j].prev = particlesCopy[j].prev;              
                        particles[j].CLindex = particlesCopy[j].CLindex;
                        particles[j].index = particlesCopy[j].index;
                }
                // copies backup of clusters log
                for (int j = 0; j < N; j++)
                {
                        clustersLog[j] = clustersLogCopy[j];
                }
                // copies backup of box size
                L = LOld;
                buildCL(1);
        }
        else
        {// accept trajectory
                //printf("trajectory accepted!\n");
                maxClusOld = maxClus;
                potBiasOld = potBias;
                // backup of new configuration
                for (int j = 0; j < N; j++)
                {
                        particlesCopy[j].x = particles[j].x;
                        particlesCopy[j].y = particles[j].y;
                        particlesCopy[j].z = particles[j].z;
                        particlesCopy[j].r = particles[j].r;
                        particlesCopy[j].type = particles[j].type;
                        particlesCopy[j].next = particles[j].next;              
                        particlesCopy[j].prev = particles[j].prev;              
                        particlesCopy[j].CLindex = particles[j].CLindex;
                        particlesCopy[j].index = particles[j].index;
                }
                // backup of new clusters log
                for (int j = 0; j < N; j++)
                {
                        clustersLogCopy[j] = clustersLog[j];
                }
                // backup of new box size
                LOld = L;

        }

        // clustersLog can be saved at this point
        FILE *saveclusters = fopen(clusterlog_filename, "a");
        if (saveclusters != NULL)
        {
                fprintf(saveclusters, "&%d\n", countLog);
                for (int i = 0; i < N; i++)
                {
                        if (clustersLog[i] != 0)
                                fprintf(saveclusters, "%d\t%d\n", i, clustersLog[i]);
                }
        }

        free(clustersLog);
}



void relaxation(int nTrajectory)
{
        /*\
         *  Function: relaxation
         *  --------------------
         *  Relaxes the system without taking measurements for nTrajectory
         *  trajectories
         *
         *  nTrajectory:        number of trajectories to relax the system for
        \*/

        double proba[2] = {0.0f, genrand()};
        double rule = 0.0f;
        double volumeStepTune = 1.0f;                   // tuning factor for volume step size
        double particleStepTune = 1.0f;                 // tuning factor for particle step size
        int sV = 0;                                     // # of volumeMove() successes                       
        int sP = 0;                                     // # of particleMove() successes
        int rPstep = 500;                               // # of trajectories before tuning particle step size
        int rVstep = (N + 1) * rPstep;                  // # of trajectories before tuning volume step size

        // performs nTrajectory trajectories
        for (int i = 0; i < nTrajectory; i++)
        {
                //printf("traj. #%d\n", i);
                clustersLog = malloc(N * sizeof(*clustersLog));
                memset(clustersLog, (int) 0.0, sizeof(*clustersLog) * N);
                //clustersLogCopy = malloc(N * sizeof(*clustersLogCopy));

                // performs 1 trajectory
                for (int j = 0; j < trajectoryLength; j++)
                {
                        for (int k = 0; k < N; k++)
                        {
                                proba[0] = genrand();
                                if (proba[0] < (1.0f / (float) (N + 1)))
                                        sV += volumeMove(volumeStepTune);
                                else
                                        sP += particleMove(particleStepTune);
                        }
                }
                //sanityCheck();          // SOMETHING IS HAPPENING HERE

                // tunes acceptance rate when necessary
                if (i % rPstep == rPstep-1)
                {
                        sanityCheck();
                        particleStepTune *= tuneStepSize(sP, rPstep, targetedAcceptanceRate);
                        sP = 0;
                }
                if (i % rVstep == rVstep-1)
                {
                        volumeStepTune *= tuneStepSize(sV, rVstep, targetedAcceptanceRate);
                        sV = 0;
                }

                // compute table of clusters, find largest cluster & need to keep track of it
                maxClus = findClusters();
                potBias = Ubias(maxClus);

                // accept or reject trajectory based on biased potential
                rule = exp(potBiasOld - potBias);
                if (proba[1] > rule)
                {// reject trajectory
                        //printf("trajectory rejected!\n");
                        // copies configuration backup
                        for (int j = 0; j < N; j++)
                        {
                                particles[j].x = particlesCopy[j].x;
                                particles[j].y = particlesCopy[j].y;
                                particles[j].z = particlesCopy[j].z;
                                particles[j].r = particlesCopy[j].r;
                                particles[j].type = particlesCopy[j].type;
                                particles[j].next = particlesCopy[j].next;              // can cause problems
                                particles[j].prev = particlesCopy[j].prev;              // can cause problems
                                particles[j].CLindex = particlesCopy[j].CLindex;
                                particles[j].index = particlesCopy[j].index;
                        }
                        // copies backup of clusters log -- no need during relaxation: no outputting of clustersLog 
                        for (int j = 0; j < N; j++)
                        {
                                clustersLog[j] = clustersLogCopy[j];
                        }
                        // copies backup of box size
                        L = LOld;
                        buildCL(1);
                }
                else
                {// accept trajectory
                        //printf("trajectory accepted!\n");
                        maxClusOld = maxClus;
                        potBiasOld = potBias;
                        // backup of new configuration
                        for (int j = 0; j < N; j++)
                        {
                                particlesCopy[j].x = particles[j].x;
                                particlesCopy[j].y = particles[j].y;
                                particlesCopy[j].z = particles[j].z;
                                particlesCopy[j].r = particles[j].r;
                                particlesCopy[j].type = particles[j].type;
                                particlesCopy[j].next = particles[j].next;              // can cause problems
                                particlesCopy[j].prev = particles[j].prev;              // can cause problems
                                particlesCopy[j].CLindex = particles[j].CLindex;
                                particlesCopy[j].index = particles[j].index;
                        }
                        // backup of new clusters log
                        for (int j = 0; j < N; j++)
                        {
                                clustersLogCopy[j] = clustersLog[j];
                        }
                        // backup of new box size
                        LOld = L;
                }
                if (makeamovie)
                        writeCoords(movie_filename);
                sanityCheck();
                free(clustersLog);
        }
}



//    MAIN
// ==========
int main(int argc, char *argv[])
{
        // INPUT
        if (argc == 2)
                init_filename = argv[1];
        else
        {
                printf("ERROR: no input filename specified\n");
                exit(0);
        }

        // PARSING INPUT
        if (sscanf(init_filename, "init_configs/n%d_PF%lf.sph", &N, &PF) != 2)
        {
                printf("Error: could not parse string\n");
                exit(0);
        }

        // OUTPUT
        if (sprintf(clusterlog_filename, "./clusters_data/clusters_n%d_P%.2lf_k%.2lf.txt", N, P, k) < 0)
        {
                printf("[US] Could not write string. Exit.\n");
                exit(0);
        }
        if (sprintf(lastsnap_filename, "./snapshots/last_n%d_P%.2lf_k%.2lf.sph", N, P, k) < 0)
        {
                printf("[US] Could not write string. Exit.\n");
                exit(0);
        }
        if (sprintf(movie_filename, "./snapshots/movie_n%d_P%.2lf_k%.2lf.sph", N, P, k) < 0)
        {
                printf("[US] Could not write string. Exit.\n");
                exit(0);
        }

        // VARIABLES
        int halfSimtime = MCCycle / trajectoryLength / 2;
        //printf("halfSimtime = %d\n", halfSimtime);


        // OUTPUT FILES CLEANING
        FILE *writefile = fopen(clusterlog_filename, "w+");
        if (writefile != NULL)
                fclose(writefile);
        
        writefile = fopen(lastsnap_filename, "w+");
        if (writefile != NULL)
                fclose(writefile);

        if (makeamovie)
        {
                writefile = fopen(movie_filename, "w+");
                if (writefile != NULL)
                        fclose(writefile);
        }

        // SYSTEM INITIALIZATION        
        init_genrand(randseed);
        readInit(init_filename);
        if (makeamovie)
                writeCoords(movie_filename);
        buildCL(1);


        clustersLog = malloc(N * sizeof(*clustersLog));
        memset(clustersLog, (int) 0.0f, sizeof(*clustersLog) * N);
        clustersLogCopy = malloc(N * sizeof(*clustersLogCopy));

        maxClusOld = findClusters_init();
        //maxClusOld = findClusters();
        potBiasOld = Ubias(maxClusOld);

        // BODY
        //printf("Relaxation\n");
        relaxation(halfSimtime);
        //sanityCheck();
        //printf("passed sanity check\n\n Main loop\n");
        
        for (int i = 0; i < halfSimtime; i++)
        {
                trajectory();
                //sanityCheck();
                if (makeamovie)
                        writeCoords(movie_filename);
                //printf("Completed traj. #%d\n", i);
        }


        // END
        writeCoords(lastsnap_filename);

        free(CLTable);
        free(particles);
        free(particlesCopy);
        free(clustersLogCopy);

        printf("MC simulation terminated normally\n"); 
        return 0;
}
