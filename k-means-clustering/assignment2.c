/*
FILENAME : assignment2.c
DESCRIPTION :
  Simple K-mean clustering with random picked centroid
NOTES :
  This C file is designed for assignment2 given by Dr. WONIL CHOI
AUTHOR :  Hong Geun Ji      START DATE : 11 Dec. 2021
CHANGES :
  NO   VERSION   DATE         WHO           DETAIL
  1    0.1      11 Dec. 2021  Hong Geun Ji  Initial version with datum setting and first clustering
  2    0.2      14 Dec. 2021  Hong Geun Ji  First datum clustering test
  3    0.3      15 Dec. 2021  Hong Geun Ji  Make the clustering function general (not only for first clustering)
  4    0.4      16 Dec. 2021  Hong Geun Ji  Input and Output file
  5    0.5      17 Dec. 2021  Hong Geun Ji  Code cleanup and add more comments
*/

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#define FILE_OPEN "test2"  // only the file name without .txt
                           // i.e. if the file full name is test.txt 
                           //      FILE_OPEN should be "test"
#define TRUE 1
#define FALSE 0
#define MAX_CHARS_LINE 200
#define MALLOC(p, t, s)  if( !((p) = (t*)malloc(s))  ){\
                            fprintf(stderr, "Insufficient memory!");\
                            exit(EXIT_FAILURE);\
                            }

typedef struct Data {
    int* elements;  // n-dimension info.
} Data;

typedef struct Cluster {
    Data** datum;
    Data* centroid;
} Cluster;

Data* get_datum(char* to_open, int* K, int* d, int* n)
{
    // to_open := txt file name to be opened
    // K := cluster number
    // d := data dimension
    // n := the number of data
    // returns the bunch of data with d dimension elements

    FILE* txt_file;  // read txt file
    if (!(txt_file = fopen(to_open, "r")))  // validation check
    {
        fprintf(stderr, "File cannot be read!\n");
        exit(EXIT_FAILURE);
    }

    char line[MAX_CHARS_LINE];
    if(fscanf(txt_file, "%d %d %d", K, d, n) == -1)  // get the first line
    {
        fclose(txt_file);
        return FALSE;
    }

    int i, j;
    Data* datum;
    MALLOC(datum, Data, sizeof(Data)*(*n)); // data bundle

    // n-dimension elements of each data
    for(i=0; i<*n; i++)
        MALLOC(datum[i].elements, int, sizeof(int)*(*d));

    // copy the values from the txt file
    for(i=0; i<*n; i++)
    {
        for(j=0; j<*d; j++) // running for d elements
        {
            // get each data element
            if(fscanf(txt_file, "%d ", &(datum[i].elements[j])) == -1)
                return NULL;
        }
    }

    // test code with rand()
    // for(i=0; i<*n; i++)
    // {
    //     for(j=0; j<*d; j++)
    //     {
    //         datum[i].elements[j] = rand() % 50 + 1;
    //     }
    // }

    fclose(txt_file);
    return datum;
}

double euclidean_distance(Data x, Data y, int d)
{
    // returns the distance between x and y

    int i;
    double distance = 0;
    for (i=0; i<d; i++)
        distance += pow(((x.elements)[i] - (y.elements)[i]), 2);

    return sqrt(distance);
}

void print_clusters(Cluster* clusters, int K, int d, int n, int iteration, FILE* out_file)
{
    // print out the clustered data to given output file

    int i, j, k;
    fprintf(out_file, "Iteration %d\n", iteration+1);
    for(i=0; i<K; i++)
    {
        fprintf(out_file, "Cluster %d:\n", i+1);
        for(j=0; j<n && clusters[i].datum[j]; j++)
        {
            for(k=0; k<d; k++)
                fprintf(out_file, "%d ", clusters[i].datum[j]->elements[k]);
            fprintf(out_file, "\n");
        }
        fprintf(out_file, "[Cluster %d has %d elements]\n\n", i+1, j);
    }
}

Data* new_centroid(Data** datum, int d, int n)
{
    // datum := array of data* (each row stands for the clustered data)
    // returns the new centroid with given datum

    Data* centroid;
    int i, j;
    int tmp;
    MALLOC(centroid, Data, sizeof(Data));
    MALLOC(centroid->elements, int, sizeof(int)*d);

    for (i=0; i<d; i++)
    {
        tmp = 0;
        for (j=0; j<n && datum[j]; j++)
        {
            tmp += datum[j]->elements[i];   // add up the dimensions element value
        }
        if (j > 0)
            centroid->elements[i] = tmp/j;  // divide into the number of each clustered data
    }
    
    return centroid;
}

int centroid_changed(Data** centroids1, Data** centroids2, int K, int d, int n)
{
    // centroids1 := new updated centroid
    // centroids2 := existing centroid
    // returns whether the centroid has changed or not

    int i, j, k;
    for(i=0; i<d; i++)
    {
        for(j=0; j<K && centroids1[j] && centroids2[j]; j++)
        {
            for(k=0; k<n && centroids1[j]->elements[k] && centroids2[j]->elements[k] > 0; k++)
            {
                // if something is not same, return TRUE right away
                if (centroids1[j]->elements[k] != centroids2[j]->elements[k])
                {
                    return TRUE;
                }
            }
        }
    }

    // if there was nothing different values, return FALSE
    return FALSE;
}

Cluster* clustering(Data** centroids, Data* datum, int K, int d, int n, char* to_out)
{
    // centroids := initial centroids (randomly picked)
    // datum := bunch of data to be clustered
    // K := the number of cluster
    // d := dimension
    // n := the number of data elements
    // to_out := output file name
    // returns the final clustered data

    Cluster* clusters;
    int i, j;
    double min_distance;    // with the data that has the minimum distance to centroid
    double distance;        // tmp calculated distance between the data and centroid
    int clusters_datum_runners[K];  // runner of each cluster
    int where_to_go;    // cluster number to be added
    int epoch = 0;  // how many iterations?
    int tmp;
    Data** current_centroids = centroids;   // existing centroids
    Data** updated_centroids;   // updated centroids
    FILE *out_file = fopen(to_out, "w"); // output file

    while(TRUE)
    {
        MALLOC(updated_centroids, Data*, sizeof(Data*)*K);
        MALLOC(clusters, Cluster, sizeof(Cluster)*K);

        // copy the each centroid to the cluster
        for (i=0; i<K; i++) {
            clusters[i].centroid = current_centroids[i];
            clusters_datum_runners[i] = 0;
            MALLOC(clusters[i].datum, Data*, sizeof(Data*)*n);
            for (j=0; j<n; j++) clusters[i].datum[j] = NULL;
        }

        // find the cluster number which makes the minimum distance
        for (i=0; i<n; i++)
        {
            min_distance = DBL_MAX;
            for (j=0; j<K; j++)
            {
                distance = euclidean_distance(datum[i], *(current_centroids[j]), d);
                if (min_distance > distance)
                {
                    min_distance = distance;
                    where_to_go = j;
                }
            }
            tmp = clusters_datum_runners[where_to_go]++;
            clusters[where_to_go].datum[tmp] = &datum[i];
        }
        
        print_clusters(clusters, K, d, n, epoch++, out_file);

        // update the centroid from the created clusters
        for (i=0; i<K; i++)
            updated_centroids[i] = new_centroid(clusters[i].datum, d, n);

        // if the centroid has not changed, there will be no more data shift so we stop here
        if(!centroid_changed(updated_centroids, current_centroids, K, d, n))
            break;
        
        // if the centroid has changed, use the updated centroids
        current_centroids = updated_centroids;
        
    }
    return clusters;
}

int main()
{
    int i, j, K, d, n;
    int centroid_idx;
    int cnt = 0;
    Data* datum;
    Data** centroids;
    char to_open[100];  // file name without extension
    char to_out[100];

     // file name set
    strcpy(to_open, FILE_OPEN);
    strcpy(to_out, to_open);
    strcat(to_out, "_out.txt");
    strcat(to_open, ".txt");

    if(!(datum = get_datum(to_open, &K, &d, &n)))
        exit(EXIT_FAILURE);
    int choosed[n];
    
    MALLOC(centroids, Data*, sizeof(Data*)*K);

    srand(time(NULL));

    // choose n random centroids
    for(i=0; i<n; i++)
        choosed[i] = FALSE;

    // randomly choose the initial centroids
    while(cnt < K)
    {
        centroid_idx = rand() % n;      // range would be [0, n)
        while (choosed[centroid_idx])   // find the centroid that has not chosen before
            centroid_idx = rand() % n;
        choosed[centroid_idx] = TRUE;
        centroids[cnt++] = &datum[centroid_idx];
    }
    
    // clustering start
    Cluster* clusters = clustering(centroids, datum, K, d, n, to_out);

    return 0;
}