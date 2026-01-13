/*
 * read_data_ym1.c
 * This program reads data from a specific file. 
 * Data in define as to be specify fram log file.
 * 
 * Disclaimer: We specilize for .ms.dat files in order to 
 * get the decorrelated mesuere already. 
 *  
 */


#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#define t000 1
#define t010 6
#define t025 16
#define t100 62
#define t_0 6.2

#define FILE_PATH "/Users/ffara/Desktop/LFC/QCDanomalies/test_run/dat/run1_20.ms.dat"

static struct{
    int dn, nn, tmax;
    double eps;
} file_head;

static struct{
    int nt;
    double **Wsl, **Ysl, **Qsl;
} data;

static double *Wact, *Yact, *Qtop;

static FILE *fdat = NULL;

static void alloc_data(void)
{
    int nn, tmax;
    int in;
    double **pp, *p;

    nn= file_head.nn;
    tmax= file_head.tmax;

    pp = malloc(3*(nn + 1)*sizeof(*pp));
    p =  malloc(3*(nn + 1)*(tmax + 1)*sizeof(*p));

    if ((pp == NULL) || (p == NULL)) {
        fprintf(stderr, "alloc_data: memory allocation error\n");
    }

    data.Wsl = pp;
    data.Ysl = pp + nn + 1;
    data.Qsl = pp + 2*(nn + 1);

    for(in = 0; in < (3*(nn + 1)); in++) {
        *pp = p;
        pp++;
        p += tmax;
    }

    Wact = p;
    p += nn + 1;
    Yact = p;
    p += nn + 1;
    Qtop = p;
}   

static void  read_file_head(void)
{
    int ir;
    int istd[3];
    double dstd[1];

    ir = fread(istd, sizeof(int), 3, fdat);
    ir += fread(dstd, sizeof(double), 1, fdat);

    file_head.dn = istd[0];
    file_head.nn = istd[1];
    file_head.tmax = istd[2];
    file_head.eps = dstd[0];

    if (ir != 4) {
        fprintf(stderr, "read_file_head: incorrect read count\n");
    }
}

static int read_data(void)
{
    int ir, nn, tmax;
    int in, t;
    int istd[1];
    double dstd[1];

    ir = fread(istd, sizeof(int), 1, fdat);

    if (ir != 1) {
        return 0;
    }

    data.nt = (int)(istd[0]);

    nn = file_head.nn;
    tmax = file_head.tmax;

    for (in = 0; in <= nn; in++) {
        for(t = 0; t < tmax; t++) {
            ir += fread(dstd, sizeof(double), 1, fdat);
            data.Wsl[in][t] = dstd[0];
        }
    }

    for (in = 0; in <= nn; in++) {
        for(t = 0; t < tmax; t++) {
            ir += fread(dstd, sizeof(double), 1, fdat);
            data.Ysl[in][t] = dstd[0];
        }
    }

    for(in = 0; in <= nn; in++) {
        for(t = 0; t < tmax; t++) {
            ir += fread(dstd, sizeof(double), 1, fdat);
            data.Qsl[in][t] = dstd[0];
        }
    }

    if (ir != (1 + 3*(nn + 1)*tmax)) {
        fprintf(stderr, "read_data: read error\n");
    }
    
    return 1;
}

int main(int argc, char *argv[])
{
    FILE *fq, *fq_t000, *fq_t010, *fq_t025, *fq_t100;

    fdat = fopen(FILE_PATH, "r");
    if (!fdat) {
        perror("fopen of fdat impossible.");
        return 1;
    }

    fq = fopen("Q.txt", "w");
    if (!fq) {
        perror("fopen Q.txt");
        fclose(fdat);
        return 1;
    }
    fq_t000 = fopen("Qt000.txt", "w");
    fq_t010 = fopen("Qt010.txt", "w");
    fq_t025 = fopen("Qt025.txt", "w");
    fq_t100 = fopen("Qt100.txt", "w");
    if (!fq_t000 || !fq_t010 || !fq_t025 || !fq_t100) {
        perror("fopen Qt+*.txt impossible.");
        if (fq) fclose(fq);
        if (fq_t000) fclose(fq_t000);
        if (fq_t010) fclose(fq_t010);
        if (fq_t025) fclose(fq_t025);
        if (fq_t100) fclose(fq_t100);
        fclose(fdat);
        return 1;
    }

    read_file_head();

    alloc_data();

    double Q[file_head.nn+1];
    for (int i = 0; i <= file_head.nn; i++) {
        Q[i] = 0.0;
    }

    int i, t, n = 0;
    
    
    while (read_data() == 1) {

        for (i = 0; i <= file_head.nn; i++) {
            Q[i] = 0.0;
            for (t = 0; t < file_head.tmax; t++) {
                Q[i] += data.Qsl[i][t];
            }
        }

        fprintf(fq, "%d", n);
        for (i = 0; i <= file_head.nn; i++) {
            
            fprintf(fq, "\t%e", Q[i]);
        }
        fprintf(fq, "\n");

        fprintf(fq_t000, "%e\n", Q[t000]);
        fprintf(fq_t010, "%e\n", Q[t010]);
        fprintf(fq_t025, "%e\n", Q[t025]);
        fprintf(fq_t100, "%e\n", Q[t100]);

        n++;
    }

    fclose(fq);
    fclose(fq_t000);
    fclose(fq_t010);
    fclose(fq_t025);
    fclose(fq_t100);
    fclose(fdat);

    return 0;
}
