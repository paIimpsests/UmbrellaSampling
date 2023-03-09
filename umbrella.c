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
#include <math.h>               // need to compile with `-lm` flag
#include <string.h>
#include <getopt.h>
#include "mt19937ar.c"





//    PREPROCESSOR CONSTANTS
// ============================
#define INVFPI 0.07957747154594766788444188168625718101 // = 1 / (4 * pi)
#define MAX_PART_CELL 40
#define MAX_NEIGHBORS 50





//    MACROS
// ============
#define OPTIONAL_ARGUMENT_IS_PRESENT \
    ((optind < argc-1 && optarg == NULL && optind < argc && argv[optind][0] != '-') \
     ? (uintptr_t) (optarg = argv[optind++]) \
     : (optarg != NULL))





//    PROTOTYPES
// ================
// structures
// ----------------
typedef struct Particle Particle;
typedef struct compl compl;
typedef struct bndT bndT;
typedef struct blistT blistT;
typedef struct posnb posnb;
typedef struct clusterBook clusterBook;
// functions
// ----------------
int parse_input(int argc, char* argv[]);
int overlapCheck(Particle* one, Particle* two);
int overlapCheckCL(int n);
int overlapCheckGlobalCL(void);
void sanityCheck(void);
void readInit(char *filename);
void writeCoords(char *filename, int fluidlike);
int particleMove(void);
int volumeMove(void);
void tuneStepSize(void);
int buildCL(int usecase);
int updateSingleCL(Particle* one);
int retrieveIndex(int u, int v, int w);
int findClusters(int init);
void build_nblist(int method);
int compare(const void * a, const void * b);
compl* calc_order(void);
double dotprod(compl *vec1, compl *vec2, int l);
double sqr(double x);
void compute_order(int l, bndT *bnd, compl *res1, compl *res2);
float plgndr(int l, int m, double x);
double facs(int l, int m);
double gammln(double xx);
double minpow(int m);
int* calc_conn(compl* orderp);
int calc_clusters(int* conn, compl* orderp, int init);
void cycle(void);
double Ubias(int sizeClus);
int concludeTrajectory(void);
int saveLogs(void);
int calcDGn(int c);






//    GLOBAL VARIABLES
// ======================
// hard spheres system
// ----------------------
int N = 4000;                           // [READ FROM IMPUT] [-N] number of particles in the system
double L = 0.0f;                        // [READ FROM INPUT] reduced box size
double LOld = 0.0f;                     // backup copy of the above
double P = 17.0f;                       // [TUNABLE] [-P] reduced pressure --- only relevant for NPT
double PD = 0.0f;                       // [TUNABLE] [-d] degree of polydispersity of the system
double sigma = 1.0f;                    // unit of length = 1 = max particle size
double PF = 0.0f;                       // [TUNABLE] [-p] system packing fraction
Particle *particles = NULL;             // [READ FROM INPUT] pointer to the table of particles data
Particle *particlesCopy = NULL;         // backup copy of the above
// mc
// ----------------------
int MCCycle = 12000;                    // [TUNABLE] [-c] number of cycles
int randseed = 1010;                    // [TUNABLE] seed for random number generator
int trajectoryLength = 20;              // [TUNABLE] length of a trajectory, in number of cycles, must be a divider of MCCycle
double relax = 0.5f;                    // [TUNABLE] [-r] proportion of simulation time the system spends relaxing
int relaxation;                         // number of cycles the system spends relaxing
int aT = 0, rT = 0;                     // total number of accepted, rejected trajectories
int aPm = 0, rPm = 0;                   // number of accepted, rejected particle moves since last step size tuning attempt
int aVm = 0, rVm = 0;                   // number of accepted, rejected volume moves since last step size tuning attempt
double pStepSize = 0.1f;                // particle move step size
double vStepSize = 0.4f;                // volume move step size
double arPm, arVm;                      // particles, volume moves acceptance rate
double tarP = 0.5;                      // [TUNABLE] particles moves targeted acceptance rate
double tarV = 0.3;                      // [TUNABLE] volume moves targeted acceptance rate
double tarTreshold = 0.02;              // [TUNABLE] targeted acceptance rate treshold
int tune = 0;                           // [TUNABLE] [-t] choice to tune the step size even after relaxation is finished
// cell lists
// ----------------------
Particle **CLTable = NULL;
double sCell1D = 0.0f;
int nCell1D = 0;
int nCell = 0;
// bop
// ----------------------
int* conn = NULL;
int* connections = NULL;
compl* order = NULL;
compl* orderp = NULL;
int* cluss = NULL;
int* size = NULL;
clusterBook* logBook = NULL;
clusterBook* logBookCopy = NULL;
const double bndLength = 1.4f;          // [TUNABLE] distance cutoff for bonds, if used
double bnd_cutoff = 0.7f;               // [TUNABLE] [-b] dot product cutoff to be call a connected particle
                                        // `c` in ten Wolde's work (`c=0.5`), `d_c` in Filion's work (`d_c=0.7`)
int nbnd_cutoff = 7;                    // [TUNABLE] [-c] minimum number of connected particles to be called a crystaline particle
                                        // `nc` in ten Wolde's work (`nc=?`), `xi_c` in Filion's work (`4<xi_c<10`)
double obnd_cutoff = 0.0f;              // [TUNABLE] [-o] dot product cutoff to be in the same cluster
                                        // additional criterion for cluster identification: <=bnd_cutoff for all touching clusters, =0.9 to see defects
double maxr2 = 0.0f;
blistT *blist;
int maxClus = 0;                        // size of the largest cluster in the system at its current state
int maxClusOld = 0;                     // backup of the size of the largest cluster in the system at its former state
int nbig = 0;                           // number of clusters of the largest recorded size
int NNidmethod = 0;                     // [TUNABLE] choice of method for NN identification, see build_nblist() doc
double degX;                            // crystallinity of the system
double degX_maxC;                       // crystallinity of the system the largest cluster contributes to
// umbrella sampling
// ----------------------
double k = 0.15;                        // [TUNABLE] [-k] coupling parameter for the bias potential
int n0 = 10;                            // [TUNABLE] [-n] targeted cluster size
double potBias = 0.0f;                  // value of the bias potential associated with the size of the largest cluster in the current state of the system
double potBiasOld = 0.0f;               // value of the bias potential associated with the size of the largest cluster in the previous state of the system
// ensemble average
// ----------------------
double* Nn_num = NULL;                  // numerator for DGn calculations running in background
double Nn_den = 0.0f;                   // denominator for DGn calculations running in background
// input/output
// ----------------------
char* init_filename;                    // name of the initial configuration input file
char lastsnap_filename[200];            // name of the last snapshot output file
char movie_filename[200];               // name of the movie output file
char clusterlog_filename[200];          // name of the clusters log output file
FILE* saveclusters;
char DGndata_filename[200];             // name of the DGn data output file
char tracking_filename[200];            // name of the system information tracking file
FILE* trackfile;
char snap_filename[200];                // name of the snapshot file used to launch the next window
int makeamovie = 0;                     // [TUNABLE] [-m] choice to make a movie of snapshots
int snapshots = 100;                    // [TUNABLE] [-m snapshots] number of snapshots for the movie
int hide_fluidlike = 0;                 // [TUNABLE] [-f] choice to reduce the size of fluidlike ('a'-type) particles in the movie snapshots
int save = 0;                           // [TUNABLE] [-l] choice to save cluster logs
int verbose = 0;                        // [TUNABLE] [-v] choiche to print more system info during simulation
int track = 0;                          // [TUNABLE] [-V] choice to save tracking of system info
int snap = 1;                           // tracker for if snapshot to launch next window was produced



//    STRUCTURES
// ================
struct Particle{
/* Structure:  Particle
 * --------------------
 * Particle in the system, implemented for the use of cell lists (CL)
 */
	double x;               // reduced x coordinate 
	double y;	        // reduced y coordinate
        double z;               // reduced z coordinate
	double r;	        // reduced radius --- redundant with sigma ... might remove later on
        char type;              // particle type --- for visualization purposes only
        Particle *next;         // pointer to the next particle in CL
        Particle *prev;         // pointer to the previous particle in CL
        int CLindex;            // index of CL in which particle is
        int index;              // index of the particle, in [0;N-1]
        int clusterindex;       // index of the cluster the particle belongs to; 0 if fluidlike
};



struct compl{
/* Structure:  compl
 * -----------------
 * Complex number
 */
        double re;              // real part
        double im;              // imaginary part
};



struct bndT{
        double nz;
        double si;
        double co;
        int n;
};



struct blistT{
/* Structure:  blistT
 * ------------------
 * List of particles connected to particles i
 */
        int n;                  // number of particles connected to particle i
        bndT * bnd;             // description of the bond connecting particle j to particle i
};



struct posnb{
        Particle* part;
        double dist;
        double dx;
        double dy;
        double dz;
};



struct clusterBook{
        int bookmark;
        int* size;
        int* qt;
};





//    FUNCTIONS
// ===============
int parse_input(int argc, char* argv[]){
/* Function:    parse_input
 * ------------------------
 * Parser for the programme execution
 *
 * argc:                ?
 * argv:                ?
 *
 * return:      0 for normal termination
 */
        // Parsing of command line arguments
	struct option longopt[] = {
		{"pressure",required_argument,NULL,'P'},
		{"packing-fraction",required_argument,NULL,'p'},
                {"cycles",required_argument,NULL,'c'},
                {"spring-constant",required_argument,NULL,'k'},
                {"cluster-size",required_argument, NULL,'n'},
                {"movie",optional_argument,NULL,'m'},
                {"particles",required_argument,NULL,'N'},
                {"polydispersity",required_argument,NULL,'d'},
                {"bnd-cutoff",required_argument,NULL,'b'},
                {"nbnd-cutoff",required_argument,NULL,'x'},
                {"obnd-cutoff",required_argument,NULL,'o'},
                {"hide-fluidlike",no_argument,NULL,'f'},
                {"randseed",required_argument,NULL,'g'},
                {"relaxation",required_argument,NULL,'r'},
                {"print-cluster-logs",no_argument,NULL,'l'},
                {"tune-step-size",no_argument,NULL,'t'},
                {"verbose",no_argument,NULL,'v'},
                {"save-verbose",no_argument,NULL,'V'},
		{"help",no_argument,NULL,'h'},
		{0,0,0,0}
        };

        int c;
        while ((c = getopt_long(argc, argv, "P:p:c:k:n:m::N:d:b:x:o:fg:r:ltvVh", longopt, NULL)) != -1){
                switch (c){
                        case 'P':
                                if (sscanf(optarg, "%lf", &P) != 1){
                                        printf("[umbrella] Could not parse pressure.\n");
                                        exit(0);
                                }
                                break;
                        case 'p':
                                if (sscanf(optarg, "%lf", &PF) != 1){
                                                printf("[umbrella] Could not parse packing fraction.\n");
                                        exit(0);
                                }
                                break;
                        case 'c':
                                if (sscanf(optarg, "%d", &MCCycle) != 1){
                                        printf("[umbrella] Could not parse number of cycles. Reverting to default: %d.\n", MCCycle);
                                }
                                break;
                        case 'k':
                                if (sscanf(optarg, "%lf", &k) != 1){
                                        printf("[umbrella] Could not parse spring constant.\n");
                                        exit(0);
                                }
                                break;
                        case 'n':
                                if (sscanf(optarg, "%d", &n0) != 1){
                                        printf("[umbrella] Could not parse cluster size.\n");
                                        exit(0);
                                }
                                break;
                        case 'm': // option with optional argument
                                makeamovie = 1;
                                if (OPTIONAL_ARGUMENT_IS_PRESENT){
                                // Handle is present
                                        if (sscanf(optarg, "%d", &snapshots) != 1){
                                        // Could also check for snapshots==0; in which case we can either revert back to makeamovie=1 or continue with default snapshots value
                                                printf("[umbrella] Could not parse number of snapshots. Reverting to default: %d.\n", snapshots);
                                        }
                                }
                                break;
                        case 'N':
                                if (sscanf(optarg, "%d", &N) != 1){
                                        printf("[umbrella] Could not parse number of particles. Reverting to default: %d. This is for naming convention only.\n", N);
                                }
                                break;
                        case 'd':
                                if (sscanf(optarg, "%lf", &PD) != 1){
                                        printf("[umbrella] Could not parse degree of polydispersity. Reverting to default: %.2lf. This is for naming convention only.\n", PD);
                                }
                                break;
                        case 'b':
                                if (sscanf(optarg, "%lf", &bnd_cutoff) != 1){
                                        printf("[umbrella] Could not parse bnd_cutoff value.\n");
                                        exit(0);
                                }
                                break;
                        case 'x':
                                if (sscanf(optarg, "%d", &nbnd_cutoff) != 1){
                                        printf("[umbrella] Could not parse nbnd_cutoff value.\n");
                                        exit(0);
                                }
                                break;
                        case 'o':
                                if (sscanf(optarg, "%lf", &obnd_cutoff) != 1){
                                        printf("[umbrella] Could not parse obnd_cutoff value.\n");
                                        exit(0);
                                }
                                break;
                        case 'f':
                                hide_fluidlike = 1;
                                break;
                        case 'g':
                                if (sscanf(optarg, "%d", &randseed) != 1){
                                        printf("[umbrella] Could not parse randseed. Falling back to default value: %d.\n", randseed);
                                }
                                break;
                        case 'r':
                                if (sscanf(optarg, "%lf", &relax) != 1)
                                        printf("[umbrella] Could not parse value for relaxation time proportion. Falling back to default value: %.2lf.\n", relax);
                                break;
                        case 'l':
                                save = 1;
                                break;
                        case 't':
                                tune = 1;
                                break;
                        case 'v':
                                verbose = 1;
                                break;
                        case 'V':
                                track = 1;
                                break;
                        case 'h':
                                printf("[umbrella]\n * Usage: ./umbrella [OPTIONS] SOURCE\n * Description: performs the Umbrella Sampling scheme based on the Monte-Carlo Metropolis method for a hard spheres system described in SOURCE and for a set of user specified simulation options.\n * Options:\n * -P [double]: fixed pressure of the system. In reduced units.\n * -p [double]: packing fraction of the system initial configuration.\n * -c [int]: number of cycles the simulation runs for.\n * -k [int]: spring constant for the bias potential.\n * -n [int]: targeted cluster size for the bias potential.\n * -m [int]: number of snapshots to take during the simulation.\n * -N [int]: number of particles in the system.\n * -d [double]: degree of polydispersity of the system.\n * -b [double]: value for bnd_cutoff in cluster identification method.\n * -x [int]: value for nbnd_cutoff in cluster identification method.\n * -o [double]: value for obnd_cutoff in cluster identification method.\n * -f: hides fluid-like particles when snapshotting.\n * -g [int]: seed for the pseudo-random number generation.\n * -r [double]: proportion of the total simulation time the system should relax for before snapshotting/measuring observables.\n * -l: save cluster logs to a file.\n * -t: continue tuning the acceptance rate even after relaxation.\n * -v: verbose.\n * -V: save verbose to a file.\n * -h: display this message and exit.\n");
                                exit(0);
                }
        }

        if(optind>argc-1){
		printf("[umbrella] Usage: ./umbrella [OPTIONS] SOURCE\n");
		exit(0);
	}
	
        init_filename = argv[optind];
        
        return 0;
}




int overlapCheck(Particle* one, Particle* two){
/* Function:    overlapCheck
 * -------------------------
 * Checks for overlap between two particles
 *
 * one:         pointer to the first particle
 * two:         pointer to the second particle
 *
 * return:      1 if there is overlap, 0 if not
 */
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



int overlapCheckCL(int n){
/* Function:    overlapCheckCL
 * ---------------------------
 * Checks for overlap between a given particle and all other particles
 * that can be found in neighbouring cells. Relies on the Cell Lists
 * implementation.
 *
 * n:           index of the particle to check overlap for
 *
 * return:      1 if there's overlap, 0 if not
 */
        int a = particles[n].x / sCell1D;
        int b = particles[n].y / sCell1D;
        int c = particles[n].z / sCell1D;
        int index;
        Particle *current = NULL;
        for (int i = a-1; i < a+2; i++){        
                for (int j = b-1; j < b+2; j++){
                        for (int k = c-1; k < c+2; k++){
                                index = retrieveIndex(i,j,k);
                                current = CLTable[index];
                                while (current != NULL){
                                        if ((particles+n != current) && overlapCheck(particles+n, current))
                                                return 1;
                                        current = current->next;
                                }
                        }
                }
        }
        return 0;
}



int overlapCheckGlobalCL(void){
/* Function:    overlapCheckGlobalCL
 * ---------------------------------
 * Checks for overlap in the whole system using CLs; also checks pairs
 * only once
 *
 * return:     1 if there is overlap, 0 if not
 */
        for (int i = 0; i < N; i++){
                if (overlapCheckCL(i))
                        return 1;
        }
        return 0;
}



void sanityCheck(void){
/* Function:    sanityCheck
 * ------------------------
 * Ends program execution if there is overlap in the system
 */ 
        if (overlapCheckGlobalCL()){
                printf("Something is VERY wrong... There's overlap...\n");
                exit(0);
        }
}



void readInit(char *filename){
/* Function:    readInit
 * ---------------------
 * Initializes the table of data of all N particles by reading
 * it from a user supplied .sph file; particles coordinates
 * are written in reduced units; also retrieves **cubic** box size
 * NB: Reading is done only for files with typesetting indicated above
 * NB: It is assumed that the supplied init file is written in standard units so
 * conversion to reduced units is applied
 *
 * *filename:   pointer to the name of the initialization .txt file
 */
        FILE *initfile = NULL;
        int i = 0;
        double rTemp = 0.0f;
        initfile = fopen(filename, "r");
        if (initfile != NULL){
                // Read value of N and allocate memory accordingly
                if (fscanf(initfile, "%*c%d%*c", &N) != 1){
                        printf("\nERROR --- Wrong input file\n\n");
                        exit(0);
                }
                particles = malloc(N * sizeof(*particles));
                particlesCopy = malloc(N * sizeof(*particlesCopy));
                Nn_num = malloc(N * sizeof(*Nn_num));
                memset(Nn_num, (double) 0.0f, N * sizeof(*Nn_num));
                conn = malloc(N * sizeof(*conn));
                connections = malloc(N * sizeof(*connections));
                order = malloc(N * (2 * 6 + 1) * sizeof(*order));
                orderp = malloc(N * (2 * 6 + 1) * sizeof(*orderp));
                cluss = malloc(N * sizeof(*cluss));
                size = malloc(N * sizeof(*size));
                logBook = malloc(sizeof(*logBook));
                logBook->size = malloc(N * sizeof(*logBook->size));
                logBook->qt = malloc(N * sizeof(*logBook->qt));
                logBookCopy = malloc(sizeof(*logBookCopy));
                logBookCopy->size = malloc(N * sizeof(*logBookCopy->size));
                logBookCopy->qt = malloc(N * sizeof(*logBookCopy->qt));
                blist = malloc(N * sizeof(*blist));
                for (i = 0; i < N; i++)
                        blist[i].bnd = malloc(MAX_NEIGHBORS * sizeof(bndT));
                // Read value of box size
                if (fscanf(initfile, "%lf %*f %*f%*c", &L) != 1){
                        printf("\nERROR --- Wrong input file\n\n");
                        exit(0);
                }
                // Populate table of particles data in standard units
                for (i = 0; i < N; i++){
                        if (fscanf  ( initfile,
                                  "%c %lf %lf %lf %lf%*c",
                                  &particles[i].type,
                                  &particles[i].x,
                                  &particles[i].y,
                                  &particles[i].z,
                                  &particles[i].r
                                )
                            != 5
                           ){
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
                for (i = 0; i < N; i++){
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
        for (int j = 0; j < N; j++){
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

}



void writeCoords(char *filename, int fluidlike){
/* Function:    writeCoords
 * ------------------------
 * Writes the position, radius, and type data for all N
 * particles in a .sph file; standard units are used
 * NB: Writing is done using the typesetting indicated above
 *
 * *filename:  pointer to the name of the output .sph file
 * fluidlike:  choice to print out fluid-like particles as small dots
 *             for visualization purposes
 */

        FILE *outfile = NULL;
        double boxSize = L * sigma;
        outfile = fopen(filename, "a");
        if (outfile != NULL){
                fprintf(outfile, "&%d\n", N);
                fprintf(outfile, "%.12lf %.12lf %.12lf\n", boxSize, boxSize, boxSize);
                switch (fluidlike){
                        case 0:
                                for (int i = 0; i < N; i++){
                                        fprintf(outfile,
                                        "%c %.12lf %.12lf %.12lf %.12lf\n",
                                        particles[i].type,
                                        particles[i].x * sigma,
                                        particles[i].y * sigma,
                                        particles[i].z * sigma,
                                        particles[i].r * sigma
                                        );
                                }
                                break;
                        case 1:
                                for (int i = 0; i < N; i++){
                                        if (particles[i].type == 'a'){
                                                fprintf(outfile,
                                                        "%c %.12lf %.12lf %.12lf %.12lf\n",
                                                        particles[i].type,
                                                        particles[i].x * sigma,
                                                        particles[i].y * sigma,
                                                        particles[i].z * sigma,
                                                        particles[i].r * sigma / 5.0f
                                                        );
                                        }
                                        else{
                                                fprintf(outfile,
                                                        "%c %.12lf %.12lf %.12lf %.12lf\n",
                                                        particles[i].type,
                                                        particles[i].x * sigma,
                                                        particles[i].y * sigma,
                                                        particles[i].z * sigma,
                                                        particles[i].r * sigma
                                                        );
                                        }
                                }
                                break;
                }
                fclose(outfile);
        }
}



int particleMove(void){
/* Function:    particleMove                
 * -------------------------
 * Tries to move a randomly selected particle in space by a small amount
 * and checks for overlap
 * If there is no overlap, overwrites the previous particle position
 * with the new one
 *
 * return:              0 if attempt failed, 1 if attempt succeeded
 */ 
        // Generate random displacements in x, y, z
        double delta[3] = {(genrand() - 0.5) / sigma * (pStepSize),
                           (genrand() - 0.5) / sigma * (pStepSize),
                           (genrand() - 0.5) / sigma * (pStepSize)
                          };
        // Randomly select a particle in the system
        int n = (int) (genrand() * N);
        Particle *selected = particles+n;
        // Set up a backup of the moved particle in case of overlap
        Particle bufferParticle = {selected->x, selected->y, selected->z};
        // Randomly move selected particle (PCB/NIC)
        selected->x = fmod((selected->x + delta[0]) + 2 * L, L);
        selected->y = fmod((selected->y + delta[1]) + 2 * L, L);
        selected->z = fmod((selected->z + delta[2]) + 2 * L, L);
        // Update Cell List
        updateSingleCL(particles+n);
        // Check for overlap and move back particle if need be
        if (overlapCheckCL(n)){
                selected->x = bufferParticle.x;
                selected->y = bufferParticle.y;
                selected->z = bufferParticle.z;
                updateSingleCL(particles+n);
                rPm++;
                return 0;
        }
        else {aPm++; return 1;}                
}



int volumeMove(void){
/* Function:    volumeMove
 * -----------------------
 * Attempts to change the volume of the simulation box by a small
 * amount according to the known acceptance rule
 * 
 * return:              0 if attempt failed, 1 if attempt succeeded
 */ 
        // Generate random volume change, define new adequate quantities
        double  vol = L * L * L,
                delta[2] = {(genrand() - 0.5) / (sigma) * (vStepSize),
                            genrand()
                           },
                newVol = vol + delta[0],                  
                newL = cbrt(newVol),                      
                ratio = newL / L,
                rule = exp(- P * (newVol - vol) + N * log(newVol / vol));
        // Reject move if rule is not respected
        if (delta[1] > rule) {rVm++; return 0;}
        // Scale the system according to the volume change
        for (int i = 0; i < N; i++){
                particles[i].x *= ratio;
                particles[i].y *= ratio;
                particles[i].z *= ratio;
        }
        L *= ratio;
        // Check for overlap in case of box shrinking and rescale cells
        // accordingly in each case
        if (ratio < 1){
                sCell1D *= ratio;
                if (sCell1D < sigma) {buildCL(3);}
                else {buildCL(2);}
                if (overlapCheckGlobalCL()){
                        for (int i = 0; i < N; i++){
                                particles[i].x /= ratio;
                                particles[i].y /= ratio;
                                particles[i].z /= ratio;
                        }
                        L /= ratio;
                        sCell1D /= ratio;
                        buildCL(3);
                        rVm++;
                        return 0;
                }
                else {aVm++; return 1;}
        }
        else {buildCL(3); aVm++; return 1;}
}



void tuneStepSize(void){
/*  Function:   tuneStepSize
 *  ------------------------
 *  Tunes the step size for volume and particle moves by +-5% depending
 *  on whether the success rate is above or below the acceptance rate
 *  (within a treshold)
 */
        // Volume moves
        arVm = (double) aVm / (double) (aVm + rVm);
        aVm = 0;
        rVm = 0;
        if (arVm > tarV + tarTreshold)
                vStepSize *= 1.05;
        else if (arVm < tarV - tarTreshold)
                vStepSize *= 0.95;

        // Particle moves
        arPm = (double) aPm / (double) (aPm + rPm);
        aPm = 0;
        rPm = 0;
        if (arPm > tarP + tarTreshold)
                pStepSize *= 1.05;
        else if (arPm < tarP - tarTreshold)
                pStepSize *= 0.95;
}



int updateSingleCL(Particle *one){
/* Function:    updatesingleCL
 * ---------------------------
 * Removes a particle from its CL, finds out to which CL it belongs based on
 * its coordinates and adds it to the top of the righ CL
 *
 * one:         pointer to the particle for which to update the CL
 *
 * return:      0 if terminated normally
 */
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



int buildCL(int usecase){
/* Function:   buildCL
 * --------------------
 * Builds, updates, and expands cell lists and their table according to
 * the new state of the system
 *
 * usecase:    choice of what operation to perform on CL among:
 *             + 1:    initialization of the CL table and CLs
 *                     contents
 *             + 2:    update of the CLs contents only, no
 *                     redefinition of the table, its size, or the
 *                     size of the cells
 *             + 3:    redefinition of the size of the cell, which
 *                     decreases, and update of the CLs contents; less
 *                     of the previously allocated memory is used in
 *                     this particular case
 *
 * return:      0 for normal termination
 */
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



int retrieveIndex(int u, int v, int w){
/* Function:   retrieveIndex
 * --------------------------
 * Retrieves index of the cell (from the table of CLs) to which a
 * particle belongs based on the coordinates of the CL in the 3D
 * coordinate system of the CL table; corresponds of a flattening
 * operation of the CL index from a 3D view to standard 1D view
 *
 * u:          x-coordinate of the cell
 * v:          y-coordinate of the cell
 * w:          z-coordinate of the cell
 *
 * return:     index of the CL in the table of CLs
 */
        int index = 0;
        index = (w + 2 * nCell1D) % nCell1D * nCell1D * nCell1D
                + (v + 2 * nCell1D) % nCell1D * nCell1D
                + (u + 2 * nCell1D) % nCell1D;
        return index;
}



int findClusters(int init){
/* Function:   findClusters
 * ------------------------
 * Wrapper for the identification of clusters from the local BOO
 * analysis
 * Updates the table accounting for the number of clusters of a given
 * size
 *
 * return:     size of the largest cluster in the current state of the
 *             system
 */
        int maxclus;     
        build_nblist(NNidmethod);
        order = calc_order();
        connections = calc_conn(order);
        maxclus = calc_clusters(connections, order, init);
        return maxclus;
}



void build_nblist(int method){
/* Function:    build_nblist
 * -------------------------
 * build the list of nearest neighbours (NN) for all particles in the system
 * based on the choice of identification method
 *
 * method:      choice of NN identification method among:
 *              + 0:    SANN
 *              + 1:    12 NN
 *              + 2:    cutoff radius rc
 */
        for (int p = N-1; p > -1; p--){
                // local variables:
                int id;
                bndT* bnd;
                double d,
                       dxy,
                       dx,
                       dy,
                       dz,
                       di;

                // find to which cell belongs the current particle
                int a = particles[p].x / sCell1D;
                int b = particles[p].y / sCell1D;
                int c = particles[p].z / sCell1D;
        
                // initialize number of nearest neighbours for particle `p` to 0
                blist[p].n = 0;

                // defining another buffer particle
                int cellIndex;
                Particle *current = NULL;

                // saving info about neigbouring particles
                posnb neighbours[N];
                int numposnb = 0;

                // cycle over all particles in all neighboring CLs
                for (int i = a-2; i < a+3; i++){
                        for (int j = b-2; j < b+3; j++){
                                for (int k = c-2; k < c+3; k++){
                                        cellIndex = retrieveIndex(i,j,k);
                                        current = CLTable[cellIndex];
                                        while (current != NULL && numposnb <= N){
                                                int currentID = current->index;
                                                if (p == currentID){
                                                // escapes if both particles are the same
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
                                                numposnb++;
                                                current = current->next;
                                        }
                                }
                        }
                } 
                // sorts neighbours by increasing distanc:
                qsort(neighbours, numposnb, sizeof(posnb), compare);
                // identify NN based on method of choice 
                int m = 0;
                switch (method){
                        case 0:
                                m = 3; 
                                int done = 0;
                                while (!done){
                                        double rim = 0;
                                        for (int i = 0; i < m; i++)
                                                rim += neighbours[i].dist / (m - 2);
                                        if (rim > neighbours[m].dist) {m++;}
                                        else {done = 1;}
                                        if (m > numposnb){ 
                                                printf("[umbrella] NN algorithm did not converge.\n");
                                                m--;
                                                done = 1;
                                        }
                                }
                                for (int i = 0; i < m; i++){
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
                                break;
                        case 1:
                                m = 12;
                                for (int i = 0; i < m; i++){
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
                                break;
                        case 2:
                                for (int i = 0; i < MAX_NEIGHBORS; i++){
                                        posnb* nb = neighbours+i;
                                        if (nb->dist < bndLength){
                                                bnd = &(blist[p].bnd[blist[p].n]);
                                                blist[p].n++;
                                                id = nb->part->index;
                                                bnd->n = id;
                                                if (id > p){
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
                                break;
                }

        }
}



int compare(const void * a, const void * b){
/* Function:    compare
 * --------------------
 * For use with qsort() function
 *
 * return:      1      if a > b
 *              0      if a = b
 *             -1      if a < b
 */
        double d = ((posnb*) a)->dist - ((posnb*) b)->dist;
        return ((0 < d) - (d < 0));
}



compl* calc_order(void){
/* Function:    calc_order
 * -----------------------
 * Computes the normalized qs (ten Wolde's local orientational order parameter) 
 *
 * return:     table of size 13*N of complex type of the normalized qs
 */
        compl* q1;
        compl* q2;
        const int l = 6;
        double temp;
        memset(orderp, (int) 0.0, sizeof(compl) * N * (l * 2 + 1));
        // NB: this memset is still mandatory to make sure orderp/order are
        // entirely overwritten

        // Compute qs
        for(int i = 0; i < N; i++){
                q1 = (orderp + i * (2 * l + 1) + l);
                for(int j = 0; j < blist[i].n; j++){
                        if(blist[i].bnd[j].n > i){
                                q2 = (orderp + blist[i].bnd[j].n * (2 * l + 1) + l);
                                compute_order(l, &(blist[i].bnd[j]), q1, q2);
                        }  
                }
        }  
        // Normalize qs
        for(int i = 0; i < N; i++){
                temp = sqrt(dotprod(orderp + i * (2 * l + 1), orderp + i * (2 * l + 1), l));
                temp = 1.0 / temp;
                for(int m = -l ; m <= l; m++){
                        (*(orderp+i * (2 * l + 1) + m + l)).re *= temp;
                        (*(orderp+i * (2 * l + 1) + m + l)).im *= temp;
                }
        }

        return orderp;
}



double dotprod(compl *vec1, compl *vec2, int l){
/* Function:    dotprod
 * --------------------
 * Computes the dot product of two complex vectors
 *
 * *vec1:      pointer to a first complex vector
 * *vec2:      pointer to a second complex vector
 *
 * return:     dot product of two complex vectors
 */
        double res = 0.0f;
        for(int m = -l; m <= l; m++)
                res += (*(vec1 + m + l)).re * (*(vec2 + m + l)).re
                       + (*(vec1 + m + l)).im * (*(vec2 + m + l)).im;

        return res;
}



double sqr(double x){
/* Function:    sqr
 * ----------------
 * Computes the square of a double
 *
 * x:          double to compute the square of
 *
 * return:     x squared
 */
        return x*x;
}



void compute_order(int l, bndT *bnd, compl *res1, compl *res2){
/* Function:    order
 * ------------------
 * Computes the spherical harmonics for two neighbouring particles and feeds
 * the sum of said spherical harmonics of all neighbouring particles for a
 * given particle in calculating the local orientational order parameter
 *
 * l:           = 6
 * bnd:         table of information about the bond
 * res1:        sum of the spherical harmonics for particle 1
 * res2:        sum of the spherical harmonics for particle 2
 */
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

        // Computes the spherical harmonics for m = 0
        p = plgndr(l,0,z);
        fc = facs(l,0);
        f = sqrt((2*l+1) * INVFPI * fc);
        r = p*f;
        (res1+0)->re += r;
        (res1+0)->im += 0;
        (res2+0)->re += r * minpow(l); // minpow(6)=1
        (res2+0)->im += 0;
        
        s=0;
        sp=0;
        c=1;
        cp=0;

        for(m = 1; m <= l; m++){
                // For m > 0
                p = plgndr(l,m,z);
                fc = facs(l,m);
                f = sqrt((2 * l + 1) * INVFPI * fc);
                r = p * f;
                // Chebyshev recursive method for computing cosine of multiple angles 
                cpp = cp;
                cp = c;
                if(m == 1){c = bnd->co;}
                else{c = 2.0 * bnd->co * cp - cpp;}
                // Chebyshev recursive method for computing sine of multiple angles 
                spp = sp;
                sp = s;
                if(m == 1){s = bnd->si;}
                else{s = 2.0 * bnd->co * sp - spp;}
                
                (res1+m)->re += r*c;
                (res1+m)->im += r*s;
                (res2+m)->re += r*c;
                (res2+m)->im += r*s;
                
                // For m < 0
                r *= minpow(m);
                (res1-m)->re += r*c;
                (res1-m)->im += -r*s;
                (res2-m)->re += r*c;
                (res2-m)->im += -r*s;
        }
}



float plgndr(int l, int m, double x){
/* Function:    plgndr
 * -------------------
 * Calculates the Legendre function P_{l,m}(x) = (1-x**2)**{m/2} (\frac{d}{dx})**m P_l(x)
 * where P_l(x) is the Legendre polynomial defined for x on [-1;1]
 *
 * l:           parameter
 * m:           parameter
 * x:           cos(theta), must be in [-1;1]
 *
 * return:      P_{l,m}(x)
 */
        // variables
        double fact,
               pll = 0.0f,
               pmm,
               pmmp1,
               somx2;
        int i,
            ll;

        // Check for normal computation of Legendre polynoms
        if (m < 0 || m > l || fabs(x) > 1.0)
                printf("Bad arguments in routine plgndr %i %i %f\n", l, m, fabs(x));
        
        pmm = 1.0;
        
        if (m > 0){
                somx2 = sqrt((1.0 - x) * (1.0 + x));
                fact = 1.0;
                for (i = 1; i <= m; i++){
                        pmm *= -fact * somx2;
                        fact += 2.0;
                }
        }

        if (l == m){return pmm;}
        else{ 
                pmmp1 = x * (2 * m + 1) * pmm;
                if (l == (m + 1)){return pmmp1;}
                else{
                        for (ll = m + 2; ll <= l; ll++){
                                pll = (x * (2 * ll - 1) * pmmp1 - (ll + m - 1) * pmm) / (ll - m);
                                pmm = pmmp1;
                                pmmp1 = pll;
                        }
                        return pll;
                }
        }
}



double facs(int l, int m){
/* Function:    facs
 * -----------------
 * Computes (l-m)!/(l+m)!
 *
 * l:           parameter
 * m:           parameter
 *
 * return:      (l-m)!/(l+m)!
 */
        static double* fac_table = NULL;
        int a, b;
        if(fac_table == NULL){
                fac_table = malloc((2*l+1) * sizeof(*fac_table));
                for(a = 0; a < 2*l+1; a++){
                        b = a - l;
                        fac_table[a]= exp(gammln(l - b + 1) - gammln(l + b + 1));
                }
        }
        return fac_table[m+l];
}



double gammln(double xx){
/* Function:    gammln
 * -------------------
 * Uses the Gamma function to compute factorials
 *
 * xx:          input number
 *
 * return:      factorial of xx-1
 */
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



double minpow(int m){
/* Function:    minpow
 * -------------------
 * Returns 1.0f if m is even, -1.0f if m is odd
 */
        if((m & 1) == 1){return -1.0;}
        else{return 1.0;}
}



int* calc_conn(compl* orderp){
/* Function:    calc_conn
 * ----------------------
 * Determines whether a particle is to be considered solid-like or not
 *
 * orderp:      table of normalized local oriental order parameters
 *
 * return:      table of int indicating whether a particle i is solid-like
 *              (conn[i]=1) or not (conn[i]=0)
 */
        int z;
        const int l = 6;
        for(int i = 0; i < N; i++){
                z = 0;
                for(int j = 0; j < blist[i].n; j++){
                        if(dotprod(orderp + i * (2 * l + 1), orderp + blist[i].bnd[j].n * (2 * l + 1), l) > bnd_cutoff)
                                z++;
                }
                if(z >= nbnd_cutoff){
                // particle is solid-like
                        conn[i]=1;
                } 
                else {
                // particle if fluid-like
                        conn[i]=0;
                }
                particles[i].type = 'a';
                particles[i].clusterindex = 0;
        }

        return conn;
}



int calc_clusters(int* conn, compl* orderp, int init){
/* Function:   calc_clusters
 * -------------------------
 *
 * return:     size of the largest cluster in the current state of the
 *             system
 */
        // VARIABLES
        int cs;         //cluster size
        int cn = 1;     //cluster id number
        int big = 0;    //largest cluster size
        int bc = -1;    //largest cluster id number
        int unread = 1;
        const int l = 6;
        
        // FUNCTIONS
        void setcluss(int pn){
                cluss[pn] = cn;
                for(int jj = 0; jj < blist[pn].n; jj++){
                        int tmp = blist[pn].bnd[jj].n;
                        if(conn[tmp] != 0
                           && cluss[tmp] == 0
                           && dotprod(orderp + pn * (2 * l + 1), orderp + tmp * (2 * l + 1), l) > obnd_cutoff
                           ){
                        //if(particle tmp is solid-like
                        //   && particle tmp does not yet belong to a cluster
                        //   && the dot product between q6(tmp) and q6(pn) > obnd_cutoff
                        //   )
                        //obnd_cutoff = 0.9 gives nice results
                        //obnd_cutoff <= bnd_cutoff gives all touching nuclei as one big nuclei
                                cs++; //cluster size goes up
                                particles[tmp].clusterindex = cn;
                                setcluss(tmp);
                        }  
                }
        }

        // BODY
        for(int i = 0; i < N; i++){
                cluss[i] = 0;           
                size[i] = 0;
        }  

        for(int i = 0; i < N; i++){
                cs = 0;
                if(conn[i] == 1 && cluss[i] == 0){
                //if(particle i is solid-like && particle i does not yet belong to a cluster)
                        cs++;
                        setcluss(i); 
                        size[cn] = cs;
                        if(cs > big){
                        //identifies largest cluster
                                big = cs;
                                bc = cn;
                        }
                        particles[i].clusterindex = cn;
                        cn++;
                }
        }
       
        // Only colour particles that belong to the largest cluster
        for (int i = 0; i < N; i++){
                if (particles[i].clusterindex == bc)
                        particles[i].type = 'b';
        }
  
        int tcs = 0;                            // total number of particles in a cluster
        nbig = 0;
        logBook->bookmark = 0;
        for(int i = 0 ; i < cn ; i++){
        // loop over all clusters
                tcs += size[i];
                if (size[i] != 0 && logBook->bookmark == 0){
                // first passage
                        logBook->size[logBook->bookmark] = size[i];
                        logBook->qt[logBook->bookmark] = 1;
                        logBook->bookmark++;
                }
                else if (size[i] != 0){
                        unread = 1;
                        for (int j = 0; j < logBook->bookmark; j++){
                        // sweep through book until bookmark-1
                                if (size[i] == logBook->size[j]){
                                // checks if page already exists
                                        logBook->qt[j]++;
                                        unread = 0;
                                }
                        }
                        if (unread){
                        // create new page if doesnt exist yet and moves bookmark
                                logBook->size[logBook->bookmark] = size[i];
                                logBook->qt[logBook->bookmark] = 1;
                                logBook->bookmark++;
                        }
                }
                if (size[i] == big){nbig++;}
        }

        if (init){
                logBookCopy->bookmark = logBook->bookmark;
                for (int i = 0; i < logBookCopy->bookmark; i++){
                        logBookCopy->size[i] = logBook->size[i];
                        logBookCopy->qt[i] = logBook->qt[i];
                }
        }
        degX = (double) tcs / (double) N;
        degX_maxC = (double) nbig * (double) big / (double) N;
        return big;
}



void cycle(void){
/*  Function:   cycle
 *  -----------------
 *  Performs a Monte Carlo cycle in {N,P,T}
 */
        for (int j = 0; j < N; j++){
                double proba = genrand();
                if (proba < (1.0f / (float) (N + 1)))
                        volumeMove();
                else
                        particleMove();
        }

}



double Ubias(int sizeClus){
/* Function:   Ubias
 * -----------------
 * Computes the bias potential for a given cluster size
 *
 * sizeClus:   cluster size to compute the bias potential for
 *
 * return:     bias potential for the aforementioned cluster size
 */
        double U = 0.0f;
        U = 0.5f * k * ((double) sizeClus - (double) n0) * ((double) sizeClus - (double) n0);
        return U;

}



int concludeTrajectory(void){
/* Function:    concludeTrajectory
 * -------------------------------
 * Core of the US scheme: accepts or rejects a trajectory based on the bias
 * potential
 *
 * return:      0 for normal termination
 */
        maxClus = findClusters(0);
        potBias = Ubias(maxClus);
        double rule = exp(potBiasOld - potBias);
        double proba = genrand();
        if (proba > rule){
        // reject trajectory
                rT++;
                maxClus = maxClusOld;
                potBias = potBiasOld;
                // copies configuration backup
                for (int j = 0; j < N; j++){
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
                logBook->bookmark = logBookCopy->bookmark;
                for (int i = 0; i < logBook->bookmark; i++){
                        logBook->size[i] = logBookCopy->size[i];
                        logBook->qt[i] = logBookCopy->qt[i];
                }
                // copies backup of box size
                L = LOld;
                buildCL(1); // really necessary?
        }
        else{
        // accept trajectory
                aT++;
                maxClusOld = maxClus;
                potBiasOld = potBias;
                // backup of new configuration
                for (int j = 0; j < N; j++){
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
                logBookCopy->bookmark = logBook->bookmark;
                for (int i = 0; i < logBookCopy->bookmark; i++){
                        logBookCopy->size[i] = logBook->size[i];
                        logBookCopy->qt[i] = logBook->qt[i];
                }
                // backup of new box size
                LOld = L;
        }
        return 0;
}



int saveLogs(void){
/* Function:    saveLogs
 * ---------------------
 * Writes cluster logs to a file
 *
 * return:      0 for normal termination
 */
        fprintf(saveclusters, "&%d\n", logBook->bookmark);
        for (int i = 0; i < logBook->bookmark; i++){
                // cluster sizes are not saved in increasing order
                fprintf(saveclusters, "%d\t%d\n", logBook->size[i], logBook->qt[i]);
        }
        return 0;
}



int calcDGn(int c){
/* Function:    calcDGn
 * --------------------
 * Runs the calculation for the free energy barrier DGn in the back
 *
 * c:           choice among:
 *              + 0:    feeds the calculation for DGn
 *              + 1:    writes calculated DGn values to a file
 * 
 * return:      0 for normal termination
 */ 
        switch (c){
                case 0:
                // Feeds DGn calculation
                        Nn_den += exp(potBias);
                        Nn_num[maxClus] += ((double) nbig / exp(- 1.0f * potBias));
                        break;
                case 1: ;
                // Writes DGn to a file
                        FILE* outfile = fopen(DGndata_filename, "a");
                        if (outfile != NULL){
                                for (int i = 0; i < N; i++){
                                        if (Nn_num[i] > 0.0f){
                                                fprintf(outfile, "%d\t%.12lf\n", i, -log(Nn_num[i] / Nn_den / (double) N));
                                        }
                                }
                                fclose(outfile);
                        }
                        break;
        }
        return 0;
}



void measurePF(void){
/* Function:    measurePF
 * ----------------------
 * Computes the reduced packing fraction for the given parameters
 *
 * return:     reduced packing fraction
 */
        PF = (double) N * sigma * sigma * sigma * M_PI / (6.0f * L * L * L);
}



void printVerbose(int i){
/* Function:    printVerbose
 * -------------------------
 * Prints out some information on the system state and simulation
 */
        printf("[umbrella] Cycle %d\n * acceptance rate\tlast 100 (1000 if relaxation is finished) cycles\n * particle moves\t%.3lf\n * volume moves\t\t%.3lf\n", i, arPm, arVm);
}



void saveVerbose(int i){
/* Function:    saveVerbose
 * ------------------------
 * Save some information on the system state and simulation
 */
        fprintf(trackfile, "%d %.3lf %.4lf %.3lf %.4lf %d %d %.4lf %.4lf %.3lf %.3lf\n", i, arPm, pStepSize, arVm, vStepSize, aT, rT, L, PF, degX, degX_maxC);
}





//    MAIN
// ==========
int main(int argc, char *argv[])
{
        // CPU TIME MEASUREMENT | START
        struct timespec beginc, endc;
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &beginc);
        // WALL TIME MEASUREMENT | START
        struct timespec beginw, endw;
        clock_gettime(CLOCK_REALTIME, &beginw);



        // PARSING INPUT
        parse_input(argc, argv);



        // OUTPUT
        if (sprintf(clusterlog_filename, "./data/clusters_n%d_P%.2lf_g%d_n0%d_k%.2lf.txt", N, P, randseed, n0, k) < 0){
                printf("[US] Could not write string. Exit.\n");
                exit(0);
        }
        if (sprintf(lastsnap_filename, "./snapshots/last_n%d_P%.2lf_g%d_n0%d_k%.2lf.sph", N, P, randseed, n0, k) < 0){
                printf("[US] Could not write string. Exit.\n");
                exit(0);
        }
        if (sprintf(movie_filename, "./snapshots/movie_n%d_P%.2lf_g%d_n0%d_k%.2lf.sph", N, P, randseed, n0, k) < 0){
                printf("[US] Could not write string. Exit.\n");
                exit(0);
        }
        if (sprintf(DGndata_filename, "./data/DGn_n%d_P%.2lf_g%d_n0%d_k%.2lf.txt", N, P, randseed, n0, k) < 0){
                printf("[US] Could not write string. Exit.\n");
                exit(0);
        }
        if (sprintf(tracking_filename, "./data/sysinfo_n%d_P%.2lf_g%d_n0%d_k%.2lf.txt", N, P, randseed, n0, k) < 0){
                printf("[US] Could not write string. Exit.\n");
                exit(0);
        }
        if (sprintf(snap_filename, "./snapshots/snap_n%d_P%.2lf_g%d_n0%d_k%.2lf.sph", N, P, randseed, n0, k) < 0){
                printf("[US] Could not write string. Exit.\n");
                exit(0);
        }
        FILE* writefile = fopen(lastsnap_filename, "w+");
        if (writefile != NULL){fclose(writefile);}
        writefile = fopen(DGndata_filename, "w+");
        if (writefile != NULL){fclose(writefile);}
        if (makeamovie){
                writefile = fopen(movie_filename, "w+");
                if (writefile != NULL){fclose(writefile);}
        }
        if (save){
                writefile = fopen(clusterlog_filename, "w+");
                if (writefile != NULL){fclose(writefile);}
                saveclusters = fopen(clusterlog_filename, "a");
                if (saveclusters == NULL) exit(0);
        }
        if (track){
                writefile = fopen(tracking_filename, "w+");
                if (writefile != NULL){
                        fprintf(writefile, "cycle arPm pStepSize arVm vStepSize aT rT L PF degC degC_maxC\n");
                        fclose(writefile);
                }
                trackfile = fopen(tracking_filename, "a");
                if (trackfile == NULL) exit(0);
        }



        // SYSTEM INITIALIZATION                i
        printf("[umbrella]\n * Starting US scheme with parameters n0 = %d, k = %.2lf, P = %.2lf\n * Supplied starting configuration: %s\n * Supplied starting configuration has N = %d\n * Simulation is to run for %d cycles\n * BOOP parameters are: bnd_cutoff = %.2lf, nbnd_cutoff = %d, obnd_cutoff = %.2lf\n * randseed = %d\n", n0, k, P, init_filename, N, MCCycle, bnd_cutoff, nbnd_cutoff, obnd_cutoff, randseed);
        if (makeamovie)
                printf(" * A movie of %d snapshots is to be taken\n", snapshots);
        if (hide_fluidlike)
                printf(" * Fluid-like particles will be reduced in size for easier visualization of clusters\n");
       
        relaxation = relax * MCCycle;
        if (makeamovie && snapshots > (MCCycle - relaxation) / trajectoryLength){
                printf("[umbrella] Error: targeted amount of snapshots suprpasses the duration of the simulation.\n[US] Exit.\n");
                exit(0);
        }
        init_genrand(randseed);
        readInit(init_filename);
        if (makeamovie)
                writeCoords(movie_filename, hide_fluidlike);
        buildCL(1);
        maxClusOld = findClusters(1);
        potBiasOld = Ubias(maxClusOld);



        // BODY
        for (int i = 0; i < MCCycle; i++){
                cycle();
                if (i%50 == 49){sanityCheck();}
                if (i%trajectoryLength == trajectoryLength-1){
                        concludeTrajectory();
                        if (snap && i > 0.1 * MCCycle && (maxClus > n0-2) && (maxClus < n0+2)) {
                        // Wait until 10% of simtime to take a snapshot to start next window
                                writefile = fopen(snap_filename, "w+");
                                if (writefile != NULL) {fclose(writefile);}
                                writeCoords(snap_filename, 0);
                                snap = 0;
                        }
                        if (i > relaxation){
                                if (save){saveLogs();}
                                calcDGn(0);
                                if (makeamovie && i%((MCCycle-relaxation)/snapshots) == ((MCCycle-relaxation)/snapshots)-1){
                                        writeCoords(movie_filename, hide_fluidlike);
                                }
                                if (i%1000 == 999){
                                        measurePF();
                                        arPm = (double) aPm / (double) (aPm + rPm);
                                        aPm = 0;
                                        rPm = 0;
                                        arVm = (double) aVm / (double) (aVm + rVm);
                                        aVm = 0;
                                        rVm = 0;
                                        if (tune)tuneStepSize();
                                        if (verbose) printVerbose(i);
                                        if (track) saveVerbose(i);
                                }
                        }
                }
                if (i < relaxation){
                        if (i%100 == 99){
                                tuneStepSize();
                                if (verbose){printVerbose(i);}
                                if (track){measurePF(); saveVerbose(i);}
                        }
                }
        }
        calcDGn(1);


        
        // END
        writeCoords(lastsnap_filename, 0);
        if (save) fclose(saveclusters);
        if (track) fclose(trackfile);
        printf("[umbrella] Program terminated normally\n[umbrella] Produced %s\n[umbrella] Produced %s\n[umbrella] Produced %s\n", lastsnap_filename, DGndata_filename, snap_filename);
        if (save){printf("[umbrella] Produced %s\n", clusterlog_filename);}
        if (track){printf("[umbrella] Produced %s\n", tracking_filename);}
        if (makeamovie){printf("[umbrella] Produced %s\n", movie_filename);}
      


        free(cluss);
        free(size);
        free(connections);
        free(order);
        free(Nn_num);
        free(CLTable);
        free(particles);
        free(particlesCopy);
        free(logBook->size);
        free(logBook->qt);
        free(logBook);
        free(logBookCopy->size);
        free(logBookCopy->qt);
        free(logBookCopy);
        for (int i = 0; i < N; i++)
                free(blist[i].bnd);
        free(blist);
        


        //  CPU TIME MEASUREMENT | END
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &endc);
        long secondsc = endc.tv_sec - beginc.tv_sec;       
        long nanosecondsc = endc.tv_nsec - beginc.tv_nsec;
        double elapsedc = secondsc + nanosecondsc * 1e-9;
        printf("[umbrella] Elapsed CPU time:\t%.3fs\n", elapsedc);
        // WALL TIME MEASUREMENT | END
        clock_gettime(CLOCK_REALTIME, &endw);
        long secondsw = endw.tv_sec - beginw.tv_sec;       
        long nanosecondsw = endw.tv_nsec - beginw.tv_nsec;
        double elapsedw = secondsw + nanosecondsw * 1e-9;
        printf("[umbrella] Elapsed wall time:\t%.3fs\n", elapsedw);

        return 0;
}
