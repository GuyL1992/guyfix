#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include "spkmeans.h"

enum Goal {wam, ddg, lnorm, jacobi, invalid};

#define EPSILON  pow(10,-5)
#define MAXITER 300

typedef struct EignValueObj{
    double value;
    int sourceIndex;
}EignValueObj;

/* implementation for K-means Calculation process and it's helper functions */

int eps_distance(double* vec1, double* vec2, int d){
    /**
     * @brief for two input vectors return 1 if are different more than EPSILON 
     *(according to the assigment defenition), else 0
     * 
     */
    double norm1;
    double norm2;
    int i;
    int result;
    double delta;

    norm1 = 0;
    norm2 = 0;

    for(i = 0; i < d; i ++){
        norm1 += pow(vec1[i],2);
        norm2 += pow(vec2[i], 2);
    }

    norm1 = pow(norm1,0.5);
    norm2 = pow(norm2,0.5);

    delta = fabs(norm1 - norm2);

    result = delta < EPSILON ? 1 : 0;
    return result;

}

double distance(double* a, double* b, int d){
    int i;
    double distance = 0;
    for (i = 0; i < d; i++)
    {
        distance += (a[i] - b[i]) * (a[i] - b[i]);
    }
    return distance;
}
int findMatch(double** centroids, double* vector,int k, int d){
    /**
     * @brief return the index of the closest centroid for a given vector
     * 
     */
    int i;
    double curr_distance;
    int closest_index = -1;
    double min_distance = 1000000;

    for (i = 0; i < k; i++){
        curr_distance = distance(vector, centroids[i],d);
        if (curr_distance < min_distance){
            min_distance = curr_distance;
            closest_index = i;
        }
    }
    return closest_index;
}

int checkNormCondition(double** centroids, double** centroids_in_progress, int k, int d){
    int i;
    int res;
    res = 1;

    for(i = 0; i < k; i++){
        res = res * eps_distance(centroids[i],centroids_in_progress[i],d);
    }

    return res;
}

void acumulate_sum(double* a, double* b, int d){
    /**
     * @brief add vector b to vector a
     * 
     */
    int i;
    for (i = 0; i < d; i ++){
        a[i] += b[i];
    }
}

void kmeansAlgorithm(double** observations, double** centroids, int n, int d,int k){
    int i;
    int j;
    int norm_condition = 0;
    int iterations = 0;
    int closest_index = 0;

    double* sizes = (double*) calloc(k, sizeof(double));
    double** centroids_in_progress = allocationMatrix(k,d);

    while(iterations < MAXITER && norm_condition == 0){
        for (i = 0; i < n; i ++){
            closest_index = findMatch(centroids, observations[i],k,d);
            sizes[closest_index] += 1;
            acumulate_sum(centroids_in_progress[closest_index], observations[i], d);
        }

        for(i = 0; i < k; i++){
            for (j = 0 ; j < d; j++){
                if (sizes[i] != 0){
                    centroids_in_progress[i][j] = centroids_in_progress[i][j] / sizes[i];
                }
            }
        }

        norm_condition = checkNormCondition(centroids, centroids_in_progress,k,d);

        if (norm_condition == 1){break;} /* After the centroids was updated check the E.Norm condition */
            
        for( i = 0; i < k; i++){
            for(j = 0; j < d; j ++){
                centroids[i][j] = centroids_in_progress[i][j];
                centroids_in_progress[i][j] = 0;
            }
            sizes[i] = 0;
        }

        iterations++;
    }

    freeMatrix(centroids_in_progress,k);
    free(sizes);
    printMatrix(centroids,k,k);    
}


/*General helper functions*/

/**
 * @brief method to convert input as a string into familiar enum
 * 
 * @param UserGoal user goal input
 * @return enum Goal 
 */

enum Goal converttoGoalEnum(char* UserGoal)
{   
        if (!strcmp(UserGoal,"wam")){
            return wam;
        }
        else if (!strcmp(UserGoal,"ddg"))
            return ddg;
        else if (!strcmp(UserGoal,"lnorm"))
            return lnorm;
        else if (!strcmp(UserGoal,"jacobi"))
            return jacobi;
        else
            return invalid;       
}
        
 
double** allocationMatrix(int n, int d){ /*allocate memory for new matrix - n X d*/ 

    int i = 0;
    double** allocatedMatrix =(double **) calloc (n,sizeof(double*));
    assert(allocatedMatrix!=NULL && "An Error Has Occured");

    for (i=0; i<n; i++){
        allocatedMatrix[i] = calloc (d,sizeof(double));
        assert(allocatedMatrix[i]!=NULL && "An Error Has Occured");
    }
    return allocatedMatrix;
}

void freeMatrix(double** matrix,int n){ /* for given 2D matrix - free all the allocated memory*/
    int i;
    for (i=0; i<n; i++){
        free(matrix[i]); 
    }
    free(matrix);
    return;
}

double* allocationVector(int n){
    double* vector;
    vector = calloc(n,sizeof(double));
    assert(vector!=NULL && "An Error Has Occured");
    return vector;
}

void copyMatrix(double** target, double** source,int n){
    /**
     * @brief copy the source matrix content into target matrix
     * 
     */
     int i;
     int j;
     for(i = 0; i < n; i++){
         for(j = 0 ; j < n; j ++){
             target[i][j] = source[i][j];
         }
     }
 }

void getDimentions(char arr[], int* n, int* d){ 
    /**
     * @brief for a given file name containing vectors data, return the vectors dimention
     * 
     */
    FILE *txt = NULL;
    char ch;
    int lines =0;
    int commas_num=0;
    txt = fopen(arr,"r");

    if(txt == NULL){
        printf("Invalid Input! size");
        exit(1);
    } 

    while (!feof(txt)){
        ch = fgetc(txt);
        if(lines ==0 && ch ==','){
            commas_num++;
        }
        if (ch == '\n'){
            lines++;
        }
    }
    fclose(txt);
    *(d) =  commas_num + 1;
    *(n) = lines;
}

double** getObservationsFromInput(char* input, int n, int d){
    /**
    * @brief read the input file and create the matrix represents the n vectors 
    * 
    * @param input file directory 
    * @param n number of vectors 
    * @param d the vector's dimention
    * @return double** - the matrix represents the n vectors 
    */

    int i = 0;
    int j = 0;
    int cnt = 0;
    FILE *txt = NULL;
    double **vectors = (double**) calloc (n,sizeof(double*));
    double *array =  calloc (n*d,sizeof(double));
    txt = fopen(input,"r");


    if (vectors == NULL || txt == NULL || array == NULL){
         printf("An Error Has Occurred");
         exit(1);
    }
    
    while (fscanf(txt, "%lf,", &array[i++])!=EOF){
        
    }
    fclose(txt);
    for(i = 0; i < n; i++){
        vectors[i]=(double*) calloc(d,sizeof(double));
        if (vectors[i] == NULL){
         printf("An Error Has Occurred");
         exit(1);
        }
    }
    for (i = 0; i < n; i++){
        for(j = 0; j < d; j++){
            vectors[i][j] = array[cnt];
            cnt ++;
        }
    }
    free(array); 
    return vectors;
}

void multiplyMatrix(double** result, double** A, double** B, int n){
    /**
     * @brief for given input matrix A and B, multiply A and B and save the result in result Matrix.
     * 
     */
    int k,i,j;
    for( i = 0; i < n; i++){
        for(j = 0; j < n; j ++){
            result[i][j] = 0;
        }
    }
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            for(k=0;k<n;k++)
            {
                result[i][j]+=A[i][k] * B[k][j];
            }
        }
    }
    return;
}

void printMatrix(double** matrix, int n, int d) { /* print a matrix of dimentions nXd as required format */
    int i,j;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j <  d - 1; j++)
        {
            printf("%.4f,", matrix[i][j]);
        }
        printf("%.4f", matrix[i][d-1]);
        putchar('\n');
    }    
}

void print_row_vector(double* vector, int n){
    int i;
    for(i=0;i<n;i++)
    {
        if (i<(n-1)){
            printf("%0.4f,",vector[i]);
        } else{
            printf("%0.4f\n",vector[i]);
        }
    }
}

/* helper functions for wam process */
double calcWeight(double* a, double* b, int d)/*compute the weight of two observations*/
{
    double sum = 0, norm=0, distance;
    int i=0;
    for (i=0;i<d;i++)
    {
        distance = (a[i]-b[i]) * (a[i]-b[i]);
        sum = sum + distance;
    }
    norm = pow(sum,0.5);
    norm = -norm/2;
    return exp(norm);
}

void formWeightedMatrix(double** weightedMatrix,double** vectorsMatrix, int n, int d){
    double weight;
    int row;
    int col;
    
    for(row=0;row<n;row++)
    {
        for(col=row+1;col<n;col++)
        {
            weight = calcWeight(vectorsMatrix[row],vectorsMatrix[col],d);
            weightedMatrix[row][col] = weight;
            weightedMatrix[col][row] = weight;
        }
    }
    return;
}
/**
 * @brief This function generates a degree matrix by all the '1' numbers of each row.
 * There are two cases: Regularmode - the values of the diagonal, notRegilarMode - 1/sqr(value on the diagonal)
 * 
 * 
 * @param weightedMatrix 
 * @param n - number of vertices
 * @param regularMode - If regularMode ==1 , it means that we want to compute the sum of each row, and if regularMode ==0 , it means 
 * that we in sqrt mode.
 * @return double** 
 */

/* helper functions for ddg process */

void formDegreeMatrix (double** degreeMatrix, double** weightedMatrix, int n, int regularMode){ 
    int row,col;
    double sum;

    for (row=0;row<n;row++){
        sum = 0;
        for (col=0;col<n;col++){
            sum += weightedMatrix[row][col];
        }
        if (regularMode == 1){
             degreeMatrix[row][row] =(sum);
        }
        else {
            degreeMatrix[row][row] = 1/sqrt(sum);
        }  
    }
    return;
}


/* helper functions for lnorm process */

/**
 * @brief form lNorm matrix out of the weighted adj matrix and the Degree Matrix
 * 
 * @param lNormMatrix - empty matrix, which the result is gonna save in
 * @param weightedMatrix 
 * @param degreeSqrtMatrix 
 * @param n - first dimention
 */
void formLnormMatrix (double** lNormMatrix, double** weightedMatrix,double** degreeSqrtMatrix, int n){ 
    int i=0,j=0;
    double** temp = allocationMatrix(n,n);
    multiplyMatrix(temp, degreeSqrtMatrix, weightedMatrix,n);
    multiplyMatrix(lNormMatrix, temp, degreeSqrtMatrix,n);
  
        for(i=0;i<n;i++)
        {
            for(j=0;j<n;j++)
            {
                if (i != j)
                {
                    lNormMatrix[i][j] = -lNormMatrix[i][j];
                    
                }
                else
                {
                    lNormMatrix[i][i] = 1 - lNormMatrix[i][i];
                }
            }
        }
    freeMatrix(temp,n);
    return;
}

/* helper functions for jacobi process */

void formRotaionAndAtagMatrix(double** P ,double** A, int n){
    /* find the largest element off the diaginal*/
    int i;
    int j;
    int r;
    int maxValueRow;
    int maxValueCol;
    double maxValue = 0;
    double s;
    double c; 
    double t;
    double theta, sign;
    double *maxRowcolumn, *maxColcolumn; 

    for (i = 0; i < n; i++){ /* the matrix is simetric so we can scan just the elemets above the diagonal*/
        for(j = i + 1; j < n; j++){
            if(fabs(A[i][j]) >= fabs(maxValue) && i != j ){
                maxValueRow = i;
                maxValueCol = j;
                maxValue = A[maxValueRow][maxValueCol];
            }
        }
    } 
    
    /*compute the rotation values - theta, s,t,c*/

    theta = (A[maxValueCol][maxValueCol]-A[maxValueRow][maxValueRow])/(2*A[maxValueRow][maxValueCol]);
    theta < 0 ? (sign = -1) : (sign = 1);
    /*tetha < 0 ? (absTetha = - tetha) : (absTetha = tetha);*/
    t = sign / (fabs(theta) + sqrt(pow(theta,2)+1));
    c = 1/(sqrt(pow(t,2)+1));
    s = t * c;

    for (i = 0; i < n; i++){ /* fill in the new  rotation matrix */
        for(j = 0; j < n; j++){
            if ((i == maxValueRow && j == maxValueRow) || (i == maxValueCol && j == maxValueCol)) /* should be complete*/
                P[i][j] = c;
            else if (i == j)
                P[i][j] = 1;
            else if (i == maxValueRow && j == maxValueCol) 
                P[i][j] = s;
            else if (i == maxValueCol && j == maxValueRow)
                P[i][j] = (-1) * s;
            else
                P[i][j] = 0;
        }
    }

    i = maxValueRow;
    j = maxValueCol;

    maxRowcolumn = allocationVector(n);
    maxColcolumn = allocationVector(n);

    for (r=0;r<n;r++)
    {
        maxRowcolumn[r] = A[r][i]; /* Icol = col maxValueRow */
        maxColcolumn[r] = A[r][j]; /* Jcol = col maxValueCol */
    }

    A[i][i] = pow(c,2)*maxRowcolumn[i] + pow(s,2)*maxColcolumn[j] - 2 * s * c * maxColcolumn[i];
    A[j][j] = pow(s,2) * maxRowcolumn[i] + pow(c,2) * maxColcolumn[j] + 2 * s * c * maxColcolumn[i];
    A[i][j] = 0;
    A[j][i] = 0;

    for(r=0; r<n;r++)
    {
        if(r!=i && r!=j)
        {
            A[r][i] = c * maxRowcolumn[r] - s * maxColcolumn[r];
            A[i][r] = A[r][i];
            A[r][j] = c * maxColcolumn[r] + s * maxRowcolumn[r];
            A[j][r] = A[r][j];
        }
    }

    free(maxRowcolumn);
    free(maxColcolumn);

    return;
}

void formIdentityMatrix(double** V, int n){ 
    int i;
    for( i = 0; i < n; i++){
        V[i][i] = 1;
    }
}

void fill_with_diagonal(double* eigenValues, double** A, int n){
    /**
     * @brief insert into eigenValues array the diagonal values of A 
     * 
     */
    int i;
    for(i=0; i < n; i++){
        eigenValues[i] = A[i][i];
    }
}

double getOff(double ** A, int n){
    /**
     * @brief return the off(A)^2 value
     * 
     */
    double offA = 0;
    int i = 0;
    int j = 0;

    for(i=0;i<n;i++)
    {
        for(j=i+1;j<n;j++)
        {
            offA += 2 * pow(A[i][j],2);
        }
    }

    return offA;
}

/**
 * @brief for an input matrix A (lNorm matrix in the spk process), perform jacobi algorithm and form it's eignVectors matrix, also find the eignValues
 * 
 * @param V - empty matrix which the eign vectors will be saved in the end of the process
 * @param A - the input matrix for jacobi algorithm
 * @param n - first dimentios
 */

 
void jaccobiAlgorithm (double** V, double** A, int n){ 

    int iteration = 0;
    double offA;
    double offAtag = getOff(A,n);
    int stopCondition = 1;
    double** P = allocationMatrix(n,n);
    double** temp = allocationMatrix(n,n);
    formIdentityMatrix(V,n); /* first: V equal to the Identity matrix */

    do{
        offA = offAtag;
        formRotaionAndAtagMatrix(P,A,n);
        multiplyMatrix(temp,V,P,n);
        copyMatrix(V,temp,n);

        offAtag = getOff(A,n);
        iteration ++;

       if ((offA-offAtag) <= EPSILON)
        {
            stopCondition = 0;
        }
        if (offAtag==0) {break;}
       
    } while(stopCondition && iteration < 100);
    
    freeMatrix(temp,n);
    freeMatrix(P,n);

    return;
 
}

/* helper functions for spk process */

int cmpfuncEignvalues (const void * a, const void * b) { /*Comperator for qsort of EignValues struct method*/
   EignValueObj first = *(EignValueObj*)a;
   EignValueObj second = *(EignValueObj*)b;

   if (first.value > second.value){return -1;}

   else if (first.value < second.value){return 1;}

   else {
       if (first.sourceIndex <= second.sourceIndex){return -1;}
       else {return 1;}
       }
}

/**
 * @brief sort the eigenVector deacrisingly and peek the greates K vectors and save them in T (kXk) Matrix  
 * @param T - empty matrix to save the result in
 * @param eigenVectors - source matrix 
 * @param eignValues - source eignValus array
 * @param n 
 * @param k 
 */
void sortEigenVectors(double**T, double** eignVectors, double* eignValues, int n,int k){

    int i;
    int j;
    int currIndex;
    EignValueObj* eignValuesObj = calloc(n,sizeof(EignValueObj));
    assert(eignValuesObj!=NULL && "An Error Has Occured");

    j = 0;
    for(i=0; i < n ; i++){ /*assign each eignValue to a new object with it's source index in the input matrix*/
        eignValuesObj[i].value = eignValues[i];
        eignValuesObj[i].sourceIndex = j;
        j ++;
        
    }

    qsort(eignValuesObj,n, sizeof(EignValueObj),cmpfuncEignvalues);
    
    for(i = 0; i < k; i++){ /* peek thw grater K vectors and form the T matrix */
        currIndex = eignValuesObj[i].sourceIndex;
        for (j = 0; j < n; j++){
            T[j][i] = eignVectors[j][currIndex];
        }  
    }

}

void normelize_matrix(double** matrix, int n, int d){ /*get "normal" valus of T matrix according to the defenition*/
    int i;
    int j;
    double s;

    double* rows_sums = allocationVector(n);

    for(i=0; i < n; i++){
        s = 0;
        for(j=0; j < d; j++){
            s+= matrix[i][j] * matrix[i][j];
        }
        rows_sums[i] = s;    
    }

    for(i=0; i < n; i++){
        for(j=0; j <d; j++){
            if (rows_sums[i] != 0){
                matrix[i][j] = matrix[i][j] / sqrt(rows_sums[i]);
            }
        }
    }

    free(rows_sums);
    return;

}

/**
 * @brief comparator using to  determine the number of cluskers k. k is the max gap between two 
 * following eingevalues, until half of the values.
 * 
 * @param eigenValues 
 * @param len len of eingevalues.
 * @return int - k 
 */

 int cmpfuncDouble (const void * a, const void * b) { /* required comparator for The Eigengap Heuristic */
     int result = 0;
     if (*(double*)a == *(double*)b)
        return result;

    result = *(double*)a < *(double*)b ? 1 : -1;
    return result;
 }

int TheEigengapHeuristic(double* eigenValues, int lenOfArr) {
    int index=0;
    double maxDelta = 0;
    double delta = 0;
    int i;

    qsort(eigenValues,lenOfArr, sizeof(double), cmpfuncDouble);
    for(i=1; i<=(lenOfArr/2);i++)
    {
        delta = fabs(eigenValues[i]-eigenValues[i-1]);
        if (delta > maxDelta){
            maxDelta = delta;
            index = i;
        }
    }
    
    return index;
}
/*Processes Functions - used for complete Process according to the User Input */

/* main process functions: 
    1. wam process
    2. ddg process 
    3.lnorm process
    4. jacobi process 
    5. get new data points for spk process (triggered by python)
*/

void wamProcess(double** observations, int n, int d){
    /**
     * @brief form and print the weighted adjacency matrix
     * 
     */
    
    double** weightedMatrix = allocationMatrix(n,n);
    formWeightedMatrix(weightedMatrix,observations,n,d);
    printMatrix(weightedMatrix,n,n);
    freeMatrix(weightedMatrix,n);
}

void ddgProcess (double** observations, int n, int d){
    /**
     * @brief form and print the  Diagonal Degree matrix
     * 
     */
    double** weightedMatrix = allocationMatrix(n,n);
    double** degreeMatrix = allocationMatrix(n,n);
    int regularMode = 1;
    formWeightedMatrix(weightedMatrix,observations,n,d);
    formDegreeMatrix(degreeMatrix,weightedMatrix,n,regularMode);
    freeMatrix(weightedMatrix,n);
    printMatrix(degreeMatrix,n,n);
    freeMatrix(degreeMatrix,n);

}

void lnormProcess (double** observations, int n, int d){
    /**
     * @brief form and print the Normalized Graph Laplacian
     * 
     */
    int regularMode = 0;
    double** weightedMatrix = allocationMatrix(n,n);
    double** degreeMatrix = allocationMatrix(n,n);
    double** lNormMatrix = allocationMatrix(n,n);
    formWeightedMatrix(weightedMatrix,observations,n,d);
    formDegreeMatrix(degreeMatrix,weightedMatrix,n,regularMode);
    formLnormMatrix(lNormMatrix, weightedMatrix,degreeMatrix,n);

    freeMatrix(weightedMatrix,n);
    freeMatrix(degreeMatrix,n);

    printMatrix(lNormMatrix,n,n);
    freeMatrix(lNormMatrix,n);

}

void jacobiProcess(double** observations, int n){
    /**
     * @brief form and print the eigenValues and Vectors produced by the jacobi algorithm 
     * 
     */

    double** eigenVectors = allocationMatrix(n,n);
    double* eigenValues = (double*) calloc (n,sizeof(double));


    jaccobiAlgorithm(eigenVectors,observations,n);
    fill_with_diagonal(eigenValues,observations,n);
    
    print_row_vector(eigenValues,n);
    printMatrix(eigenVectors,n,n);

    freeMatrix(eigenVectors,n);
    free(eigenValues);

    return;

}

double** getDataPoints(double** observations, int n, int d, int* k){
    /**
     * @brief return the T matrix ( n X k ) represent the new data point as the part of SPK proccess 
     * 
     */

    int regularMode = 0;
    double** weightedMatrix = allocationMatrix(n,n);
    double** degreeMatrix = allocationMatrix(n,n);
    double** lNormMatrix = allocationMatrix(n,n);
    double** T;
    double** eigenVectors = allocationMatrix(n,n);
    double* eignValues = (double*) calloc (n,sizeof(double));
    double* eignValuesForHuristic = (double*) calloc (n,sizeof(double));
    
    formWeightedMatrix(weightedMatrix,observations,n,d);
    formDegreeMatrix(degreeMatrix,weightedMatrix,n,regularMode);
    formLnormMatrix(lNormMatrix, weightedMatrix,degreeMatrix,n);
    freeMatrix(weightedMatrix,n);
    freeMatrix(degreeMatrix,n);

    

    jaccobiAlgorithm(eigenVectors,lNormMatrix,n);    
    fill_with_diagonal(eignValues,lNormMatrix,n);
    fill_with_diagonal(eignValuesForHuristic,lNormMatrix,n);
    freeMatrix(lNormMatrix,n);



    *(k) = (*(k) == 0) ? TheEigengapHeuristic(eignValuesForHuristic,n) : *(k);
    
    T = allocationMatrix(n,*k);
    if (n > *k)
        sortEigenVectors(T,eigenVectors,eignValues,n,*(k));

    freeMatrix(eigenVectors,n);
    normelize_matrix(T,n,*(k));
    return T;
}



int main(int argc, char *argv[]){ 

    char* input;
    char* flow;
    int d;
    int n;
    double ** observations;

    if (argc < 3){
        printf("Invalid Input!");
        exit(0);
    }


    flow = argv[1];
    input  = argv[2];
    
    getDimentions(input,&n,&d);
    observations= getObservationsFromInput(input,n,d);
    
    switch(converttoGoalEnum(flow)){ 

        case wam:
            wamProcess(observations,n,d);
            break;

        case ddg:
            ddgProcess(observations,n,d);
            break;

        case lnorm:
            lnormProcess(observations,n,d);
            break;

        case jacobi:
            jacobiProcess(observations,n);
            break;
        
        default:
            printf("Invalid Input!");
            freeMatrix(observations,n);
            exit(0);
    }
    freeMatrix(observations,n);
    return 0;
}




