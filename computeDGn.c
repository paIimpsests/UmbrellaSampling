/*
 * Subroutine for calculation of \Delta G(n) from the Umbrella Sampling scheme
 */





//    LIBRARIES
// ===============
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <time.h>
#include <math.h>       // need to compile with `-lm` flag
#include <string.h>





//    PREPROCESSOR CONSTANTS
// ============================
#define MAX_PART_CELL 40
#define MAX_NEIGHBORS 50
#define INVFPI 0.07957747154594766788444188168625718101 // = 1 / (4 * pi)





//    PROTOTYPES
// ================
// structutres
// ----------------
typedef struct Logs Logs;
// functions
// ----------------
void initReading(FILE* inputfile);
void readLogs(FILE* inputfile);
double Ubias(int sizeClus);
void calcNn(void);





//    GLOBAL VARIABLES
// ======================
// hard spheres system
// ----------------------
int N = 0;
double P = 0.0f;
// umbrella sampling
// ----------------------
int n0 = 10;
double k = 0.15;
// ensemble average
// ----------------------
Logs* Nn_ri = NULL;
double* Nn_num = NULL;
double Nn_den = 0.0f;
// input / output
// ----------------------
int count = 0;
int nLines;
char* input_filename;
long cursor_end;
long cursor_current;
char DGndata_filename[200];
char ndata_filename[200];




//    STRUCTURES
// ================
struct Logs
{
        int n;                  // cluster size
        int Nn;                 // number of clusters of cluster size n
};




//    FUNCTIONS
// ===============
void initReading(FILE* inputfile)
{
        fseek(inputfile, 0, SEEK_END);
        cursor_end = ftell(inputfile);
        rewind(inputfile);
        cursor_current = ftell(inputfile);
}



void readLogs(FILE* inputfile)
{
        if (fscanf(inputfile, "%*c%d%*c", &nLines) != 1)
        {
                printf("Error #1\n");
                exit(0);
        }
        for (int i = 0; i < nLines; i++)
        {
                if (fscanf(inputfile, "%d\t%d%*c", &Nn_ri[i].n, &Nn_ri[i].Nn) != 2)
                {
                        printf("Error #2\n");
                        exit(0);
                }
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

        double w =0.0f;
        w = 0.5f * k * ((double) sizeClus - (double) n0) * ((double) sizeClus - (double) n0);
        return w;

}



void calcNn(void)
{
        // Contribution to the denominator for one snapshot
        Nn_den += exp(Ubias(Nn_ri[nLines].Nn));

        // Contribution to the denomitnator of Nn for one snapshot
        for (int i = 0; i < nLines; i++)
                Nn_num[Nn_ri[i].n] += Nn_ri[i].Nn / exp(- Ubias(Nn_ri[nLines].Nn));
}



void saveLogs(void)
{
        FILE* writefile = fopen(ndata_filename, "a");
        if (writefile != NULL)
        {
                // Save all cluster sizes for a given trajectory iteration
                /*
                for (int i = 0; i < nLines; i++)
                {
                        fprintf(writefile, "%d\t%d\n", count, Nn_ri[i].n);
                }
                */
                
                // Save the largest cluster size only
                if (nLines == 0)
                        fprintf(writefile, "%d\t%d\n", count, 0);
                else
                        fprintf(writefile, "%d\t%d\n", count, Nn_ri[nLines-1].n);
                        
        }
        fclose(writefile);
}



//    MAIN
// ==========
int main(int argc, char* argv[])
{
        // INPUT
        if (argc == 2)
                input_filename = argv[1];
        else
        {
                printf("ERROR - no input file supplied\n");
                exit(0);
        }


        // PARSING INPUT
        if (sscanf(input_filename, "clusters_data/clusters_n%d_P%lf.txt", &N, &P) != 2)
        {
                printf("Error: Could not parse string\n");
                exit(0);
        }


        // OUTPUT n
        if (sprintf(ndata_filename, "./n_data/n_n%d_P%.2lf.txt", N, P) < 0)
        {
                printf("[US] Could not write string. Exit.\n");
                exit(0);
        }
        FILE *writefile2 = fopen(ndata_filename, "w+");
        if (writefile2 != NULL)
                fclose(writefile2);


        // OUTPUT DGn
        if (sprintf(DGndata_filename, "./DGn_data/DGn_n%d_P%.2lf.txt", N, P) < 0)
        {
                printf("[US] Could not write string. Exit.\n");
                exit(0);
        }
        FILE *writefile = fopen(DGndata_filename, "w+");
        if (writefile != NULL)
                fclose(writefile);


        // INITIALIZATION
        Nn_ri = malloc(N * sizeof(*Nn_ri));
        Nn_num = malloc(N * sizeof(*Nn_num));
        memset(Nn_num, (double) 0.0f, N * sizeof(*Nn_num));


        // BODY
        FILE *initfile = NULL;
        initfile = fopen(input_filename, "r");
        if (initfile != NULL)
        {
                initReading(initfile);
                while (cursor_current < cursor_end)
                {
                        count++;
                        readLogs(initfile);
                        calcNn();
                        saveLogs();
                        cursor_current = ftell(initfile);
                }
        }


        // SAVING DATA
        writefile = fopen(DGndata_filename, "a");
        if (writefile != NULL)
        {
                for (int i = 0; i < N; i++)
                {
                        if (Nn_num[i] > 0.0f)
                        {
                                fprintf(writefile, "%d\t%.12lf\n", i, -log(Nn_num[i] / Nn_den / N));
                        }
                }
        }
        fclose(writefile);


        // END
        free(Nn_ri);
        free(Nn_num);

        return 0;
}
