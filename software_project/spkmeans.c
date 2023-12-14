#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>

#include "spkmeans.h"

#ifndef point_h
#define point_h

struct point{
    double (*cordinates);
    int length;
};

typedef struct point Point;

Point* newPoint(int i){
    int j;
    
    Point* p = (Point*) malloc(sizeof(Point));
    
    if(p==NULL){
        return NULL;
    }
    
    p->cordinates = (double*) malloc(i*(sizeof(double)));
    
    if(p->cordinates==NULL){
        free(p);
        return NULL;
    }
    
    
    p->length = i;
    for(j=0; j < (p->length); j++){
        p->cordinates[j] = 0;
    }
    return p;
}



void resetPoint(Point *p){
    int i;
    
    for(i=0; i < (p->length); i++){
        p->cordinates[i] = 0;
    }
}

void sumPoints(Point *p1, Point *p2){
    int i;
    for(i=0; i < (p1->length); i++){
        p1->cordinates[i] = p1->cordinates[i] + p2->cordinates[i];
    }
}



void copyCordinates(Point *p1, Point *p2){
    int i;
    for(i=0; i < (p1->length); i++){
        p1->cordinates[i] = p2->cordinates[i];
    }
}



void deletePoint(Point *p ){
    double *d;
    if(p==NULL){
        return;
    }
    d = (p)->cordinates;
    if(d!=NULL){
        free(d);
    }
    
    free(p);
    
    d=NULL;
    p=NULL;
}

Point** new_array_Points(int len, int size_Point){
    int i;
    int j;
    Point* p;
    Point** array_Points = (Point**) calloc(len, sizeof(Point*));
    if(array_Points==NULL){
        return NULL;
    }
    for(i=0; i<len; i++){
        p = newPoint(size_Point);
        if(p==NULL){
            for(j=0; j<i; j++){
                deletePoint(array_Points[j]);
            }
            free(array_Points);
            return NULL;
        }
        array_Points[i] = p;
    }
    return array_Points;
}

void delete_array_Points(Point** array, int len){
    int i;
    for(i=0; i<len; i++){
        deletePoint(array[i]);
    }
    free(array);
}

#endif /* point_h */


#ifndef list_h
#define list_h


typedef double Data;

struct node
{
    Data data;
    
    struct node *next;
    
    Point *point;
    
    int index;
    
};

typedef  struct node Node;

int push(Node** head, Data new_data){
    /* 1. allocate node */
    Node* new_node = (Node*) malloc(sizeof(Node));
    
    if(new_node==NULL){
        return 1;
    }
 
    /* 2. put in the data  */
    new_node->data  = new_data;
 
    /* 3. Make next of new node as head */
    new_node->next = (*head);
 
    /* 4. move the head to point to the new node */
    (*head)    = new_node;
    
    return 0;
    
}

int pushP(Node** head){
    /* 1. allocate node */
    Node* new_node = (Node*) malloc(sizeof(Node));
 
    if(new_node==NULL){
        return 1;
    }
    
    /* 3. Make next of new node as head */
    new_node->next = (*head);
 
    /* 4. move the head to point to the new node */
    (*head)    = new_node;
    
    return 0;
    
}





void deleteL( Node* head ) {
    if ( head != NULL ) {
        deleteL( head->next );
        free(head);
        head = NULL;
    }
}

void deletePlist( Node* head ) {
    if ( head != NULL ) {
        deletePlist(head->next);
        deletePoint(head->point);
        free(head);
        head = NULL;
    }
}




#endif /* list_h */

#ifndef Cluster_h
#define Cluster_h

struct Cluster{
    Point* centroid;
    
    Point* sumPoints;
    
    int numberOfPoints;
    
    int converge;

}
;

typedef struct Cluster cluster;

int newCluster(cluster* c, int i){

    c->centroid = newPoint(i);
    
    if(c->centroid==NULL){
        return 1;
    }
    
    c->sumPoints = newPoint(i);
    
    if(c->sumPoints==NULL){
        return 1;
    }
    
    
    c->numberOfPoints = 0;
    
    c->converge = 0;
    
    return 0;
}

void createArr(cluster** cl, int k){

    *cl = (cluster*) malloc(k*sizeof(cluster));
    
}


void deleteArrOfCluster(cluster* a, int k){
    int i;
    for(i=0; i<k; i++){
        if(a[i].centroid!=NULL){
            deletePoint(a[i].centroid);
        }
        if(a[i].sumPoints!=NULL){
            deletePoint(a[i].sumPoints);
        }
    }
}

#endif /* Cluster_h */


#ifndef matrix_h
#define matrix_h


typedef struct{
    int rows;
    int columns;
    double** vals;
}matrix;


matrix* initMatrix(int r, int c){
    int i;
    double** m;
    double* arr;
    matrix* mat;
    if(r<1||c==1){
        return NULL;
    }
    mat  = (matrix*) malloc(sizeof(matrix));
    if(mat==NULL) {return NULL;}
    
    mat->rows = r;
    mat->columns = c;
    
    arr = (double*) calloc(r*c, sizeof(double));
    if(arr==NULL) {
        free(mat);
        return NULL;
    }
    m = (double**) calloc(r, sizeof(double*));
    if(m==NULL) {
        free(mat);
        free(arr);
        return NULL;
    }
    mat->vals = m;
    for(i=0; i<r; i++){
        (mat->vals)[i] = &(arr[i*c]);
    }
    return mat;
}

matrix* multMatrices(matrix* m1, matrix* m2){
    int i, j, n;
    double sum;
    double** d1 = m1->vals;
    double** d2 = m2->vals;
    matrix* m3;
    m3 = initMatrix(m1->rows, m2->columns);
    if(m3==NULL) {return NULL;}
        
        for(i=0; i<m1->rows; i++){
            for(j=0; j<m2->columns; j++){
                sum=0;
                for(n=0; n < (m1->columns); n++){
                    sum+=(d1[i][n])*(d2[n][j]);
                }
                m3->vals[i][j] = sum;
            }
        }
    return m3;
}

void newmultMatrices(matrix* result, matrix* m1, matrix* m2){
    int i, j, n;
    double sum;
    double** d1 = m1->vals;
    double** d2 = m2->vals;
    
    for(i=0; i<m1->rows; i++){
        for(j=0; j<m2->columns; j++){
            sum=0;
            for(n=0; n < (m1->columns); n++){
                sum+=(d1[i][n])*(d2[n][j]);
            }
            result->vals[i][j] = sum;
        }
    }
    return;
}

void scalarMultMatrix(matrix* m, int scalar){
    int i;
    int j;
    for(i=0; i<m->rows; i++){
        for(j=0; j<m->columns; j++){
            m->vals[i][j] =  (m->vals[i][j])*scalar;
        }
    }
}

void subMatrices(matrix* m1, matrix* m2){
    int i;
    int j;
    for(i=0; i<m1->rows; i++){
        for(j=0; j<m1->columns; j++){
            m1->vals[i][j] =  (m1->vals[i][j]) - (m2->vals[i][j]);
        }
    }
}


void deleteMatrix(matrix* m){
    free(m->vals[0]);
    free(m->vals);
    free(m);
}


/* isDiagonalMatrix  return 0 if not and 1 if yes   */

int isDiagonalMatrix(matrix* m){
    int i;
    int j;
    for(i=0; i<m->rows; i++){
        for(j=0; j<m->columns; j++){
            if(i==j){
                continue;
            }
            else if( (m->vals[i][j]) != 0){
                return 0;
            }
        }
    }
    return 1;
}


void findPivot(matrix* m, int* pivot){
    int i;
    int j;
    double max = fabs(m->vals[0][1]);
    *pivot = 0;
    *(pivot+1) = 1;
    for(i=0; i<m->rows; i++){
        for(j=i+1; j < m->rows; j++){
            if(fabs(m->vals[i][j]) > max){
                max = fabs(m->vals[i][j]);
                *pivot = i;
                *(pivot+1) = j;
            }
        }
    }
}

matrix* ID_NewMatrix(int r, int c){
    int i;
    matrix* m;
    
    m = initMatrix(r, c);
    if(m == NULL) return NULL;
    
    for(i=0; i<m->rows; i++){
        m->vals[i][i] = 1;
    }
    return m;
}

void resetToIDMatrix(matrix* m){
    int i;
    int j;
    for(i=0; i<m->rows; i++){
        for(j=0; j<m->columns; j++){
            if(i==j){m->vals[i][j]=1; }
            else{ m->vals[i][j]=0; }
        }
    }
}

matrix* generate_Lnorm(matrix* d12, matrix* w){
    matrix* tmp;
    
    /* tmp = d12 * w but we create new matrix before */
    tmp = multMatrices(d12, w);
    if(tmp==NULL){
        return NULL;
    }
    
    /* w = tmp * d12 */
    newmultMatrices(w, tmp, d12);
    
    resetToIDMatrix(tmp);
    subMatrices(tmp, w);
    
    /* tmp = Lnorm now*/
    return tmp;
}

double off(matrix* m){
    int i;
    int j;
    double sum = 0;
    for(i=0; i<m->rows; i++){
        for(j=0; j<m->rows; j++){
            if(i!=j){ sum+= pow(m->vals[i][j], 2);}
            }
        }
    return sum;
}

int sign(double d){
    if(d<0) return -1;
    else return 1;
}

double calc_t(double d){
  return  sign(d) / (sqrt(pow(d, 2)+1)+fabs(d));
    
}

double calc_c(double t){
    return 1 / (sqrt(pow(t, 2)+1));
}

void matCopy_A_to_B(matrix* A, matrix* B){
    int i;
    int j;
    for(i=0; i<A->rows; i++){
        for(j=0; j<A->columns; j++){
            B->vals[i][j] = A->vals[i][j];
        }
    }
}


matrix** jacobi(matrix* A){
    matrix* newA;
    double isconverge;
    double epsilon = 1.0e-5;
    int iterNumb;
    int i,j;
    double s,c,t;
    double inter_val;
    matrix* tmpA;
    matrix** eigenvectors;
    matrix** doubleMatrices;
    int* pivot;
    matrix* p;
    matrix* pt;
    
    doubleMatrices  = (matrix**) calloc(2,sizeof(matrix*));
    if(doubleMatrices==NULL){ return NULL;}
    
    eigenvectors = doubleMatrices+1;
    
    *eigenvectors = ID_NewMatrix(A->rows, A->columns);
    if(*eigenvectors==NULL){
        deleteMatrix(A);
        free(doubleMatrices);
        return NULL;
    }
    if(isDiagonalMatrix(A)){
        *doubleMatrices = A;
        return doubleMatrices;
    }
    
    pivot = calloc(2, sizeof(int));
    if(pivot==NULL){
        deleteMatrix(A);
        deleteMatrix(*eigenvectors);
        free(doubleMatrices);
        return NULL;
    }
    
    p = ID_NewMatrix(A->rows, A->columns);
    if(p==NULL){
        deleteMatrix(A);
        free(pivot);
        deleteMatrix(*eigenvectors);
        free(doubleMatrices);
        return NULL;
    }
    
    pt = ID_NewMatrix(A->rows, A->columns);
    if(pt==NULL){
        deleteMatrix(A);
        deleteMatrix(p);
        free(pivot);
        deleteMatrix(*eigenvectors);
        free(doubleMatrices);
        return NULL;
    }
    
    tmpA = initMatrix(A->rows, A->columns);
    newA = initMatrix(A->rows, A->columns);
    
    for(iterNumb=0; iterNumb<100; iterNumb++){
        findPivot(A, pivot);
        i = *pivot;
        j = *(pivot+1);
        inter_val = ( A->vals[j][j] - A->vals[i][i] ) / ( 2*(A->vals[i][j]) );
        t = calc_t(inter_val);
        c = calc_c(t);
        s = t*c;
        p->vals[i][i] = (pt->vals[i][i] = c);
      
        p->vals[j][j] =  (pt->vals[j][j] = c);
        
        pt->vals[j][i] = (p->vals[i][j] = s);
        
        pt->vals[i][j] = (p->vals[j][i] = ((-1)*s));
        
        newmultMatrices(tmpA, pt, A);
        newmultMatrices(newA, tmpA, p);

        isconverge = off(A) - off(newA) - epsilon;
        matCopy_A_to_B(newA, A);
        matCopy_A_to_B(*eigenvectors, newA);
        
        newmultMatrices(*eigenvectors, newA, p);

        resetToIDMatrix(p);
        resetToIDMatrix(pt);
        if(isconverge<=0) break;
    }
    *doubleMatrices = A;
    deleteMatrix(p);
    deleteMatrix(pt);
    deleteMatrix(newA);
    deleteMatrix(tmpA);
    free(pivot);
    return doubleMatrices;
}

void printmat(matrix* mat){
    int i, j;
    for(i=0; i<mat->rows; i++){
        for(j=0; j<mat->columns; j++){
            printf("%0.4f", mat->vals[i][j]);
            if(j!=mat->columns-1){
                printf("%c", ',');
            }
        }
        printf("%c", '\n');
    }
}



#endif /* matrix_h */






        /* ################ start    declaration of functions    ################    */




/* ################ end   declaration of functions    ################    */





/* ################ start    declaration of functions    ################

 **************************************************************************
 **************************************************************************
 
################ new_kmean_fit FUNCTIONS DECLARATION  ################    */
int isNumber( const char number[]);

double Norm( Point* p1, Point* p2);

void matrix_array_Points(matrix* m, Point** arr);

void sort(int* array_indices, double* array_value, int len);

double* matrix_array_eigenvalues(matrix* m);

int find_k_from_arr(double* arr, int len);

void matrix_array_doubles(matrix* m, double* arr);

int* new_arr_of_indices(int n);

void generate_U_matrix(matrix* eighnvetors,  matrix* U,  int* array_indices);

void normalize_to_U(matrix* U);

void free_double_mat(matrix** two_mat);

void array_doubles_to_matrix(matrix* m, double* arr);

/* ############## END declaration of functions    ################  */


/* ################ START  kmean_new_fit   ################    */



int kmean_new_fit(double** plist_to_return ,double* megaListpoints, long* numberofpoints, long* pointSize, long* k, int goal){
    
    long i;
    
    long j;
    
    Point* p = NULL;
    
    Point** X_mat;
    
    matrix* W_mat, *D_mat, *Lnorm, **two_mat, *U_mat;
    
    double* arr_of_doubles;
    
    int* arr_of_indices;
    
    long numOfPoints = *numberofpoints;
    
    if(goal==4){
        Lnorm = initMatrix((int)*numberofpoints, (int)*numberofpoints);
        
        array_doubles_to_matrix(Lnorm, megaListpoints);
        
        two_mat = jacobi(Lnorm);
        
        *plist_to_return = calloc((two_mat[1]->columns + two_mat[1]->columns*two_mat[1]->rows), sizeof(double));
        
        if(*plist_to_return==NULL){
            free_double_mat(two_mat);
            return 1;
        }
        
        (*pointSize) = two_mat[1]->columns;
        *numberofpoints = numOfPoints + 1;
        matrix_array_doubles(two_mat[1], *plist_to_return + two_mat[1]->columns);
        
        for(i=0; i<two_mat[0]->rows; i++){
            /* make -0.0000 ----> 0 */
            if((two_mat[0]->vals[i][i] < 0) && ((two_mat[0]->vals[i][i] * 10) >= -0.0001)){
                (*plist_to_return)[i] = 0;
            }
            else (*plist_to_return)[i] = two_mat[0]->vals[i][i];
            
        }
        
        free_double_mat(two_mat);
        return 0;
    }
    
    else{
        X_mat  = new_array_Points((int)numOfPoints, (int)(*pointSize));
        if(X_mat==NULL){
            return 1;
        }
        
        for(i=0; i<numOfPoints; i++){
            p = X_mat[i];
            for(j=0; j<(*pointSize); j++){
                p->cordinates[j] = megaListpoints[(*pointSize)*i+j];
            }
        }
         /*make the W_matrix from the X_matrix*/
        W_mat = initMatrix((int)numOfPoints, (int)numOfPoints);

        if(W_mat==NULL){
            delete_array_Points(X_mat, (int)numOfPoints);
            return 1;
        }

        for(i=0; i<numOfPoints; i++){
            for(j=0; j<numOfPoints; j++){
                if(i==j) continue;
                (W_mat->vals)[i][j] = exp(-(sqrt(Norm(X_mat[i] , X_mat[j]))/2));
            }
        }
        /*we are done with X, free it*/
        delete_array_Points(X_mat, (int)numOfPoints);
        
        if(goal==1){
            *plist_to_return = calloc((W_mat->columns*W_mat->rows), sizeof(double));
            if(*plist_to_return==NULL){
                deleteMatrix(W_mat);
                return 1;
            }
            *pointSize = W_mat->columns;
            matrix_array_doubles(W_mat, *plist_to_return);
            deleteMatrix(W_mat);
            return 0;
        }

        /*make the D_matrix from the W_matrix*/
        D_mat = initMatrix((int)numOfPoints, (int)numOfPoints);
        
        if(D_mat==NULL){
            deleteMatrix(W_mat);
            return 1;
        }
        
        for(i=0; i<D_mat->columns; i++){
            double sum=0;
            for(j=0; j<D_mat->columns; j++){
                sum+=(W_mat->vals)[i][j];
            }
            D_mat->vals[i][i] = sum;
        }
        
        if(goal==2){
            
            *plist_to_return = calloc((D_mat->columns*D_mat->rows), sizeof(double));
            if(*plist_to_return==NULL){
                deleteMatrix(W_mat);
                deleteMatrix(D_mat);
                return 1;
            }
            *pointSize = D_mat->columns;
            matrix_array_doubles(D_mat, *plist_to_return);
            deleteMatrix(W_mat);
            deleteMatrix(D_mat);
            return 0;
        }
        
        /*   make the D^0.5_matrix from the D_matrix but we still use the D_mat pointer */
        for(i=0; i<D_mat->columns; i++){
            D_mat->vals[i][i] = pow((D_mat->vals[i][i]), -0.5);
        }
        
        /* generate new Lnorm matrix*/
        Lnorm =  generate_Lnorm(D_mat, W_mat);
        
        deleteMatrix(W_mat);
        deleteMatrix(D_mat);
        
        if(Lnorm==NULL){
            return 1;
        }
        
        
        if(goal==3){
            *plist_to_return = calloc((Lnorm->columns*Lnorm->rows), sizeof(double));
            
            if(*plist_to_return==NULL){
                deleteMatrix(Lnorm);
                return 1;
            }
            
            *pointSize = Lnorm->columns;
            
            matrix_array_doubles(Lnorm, *plist_to_return);
            deleteMatrix(Lnorm);
            return 0;
        }
    }
    
    two_mat = jacobi(Lnorm);
    
    if(two_mat==NULL){
       return 1;
    }
    
/* two_mat return 2 pointers to matrices. first is the eigenvalues of Lnorm and the 2nd is the  eigenvectors. */
    
    
    arr_of_doubles = matrix_array_eigenvalues(two_mat[0]);
    if(arr_of_doubles==NULL){
        free_double_mat(two_mat);
        return 1;
    }
    
    arr_of_indices = new_arr_of_indices(two_mat[0]->columns);
    
    if(arr_of_indices==NULL){
        free_double_mat(two_mat);
        free(arr_of_doubles);
        return 1;
    }
    
    /*each index represent a column in the eigenvectors matrix, we will sort it.*/
    
    sort(arr_of_indices, arr_of_doubles, two_mat[0]->columns);
    
    if((*k)==0){
        (*k) = find_k_from_arr(arr_of_doubles, (two_mat[0]->columns)/2);
        if(*k==1){
            free_double_mat(two_mat);
            free(arr_of_doubles);
            free(arr_of_indices);
            return 1;
        }
    }
    (*pointSize) = (*k);
    
    free(arr_of_doubles);
    
    U_mat = initMatrix(two_mat[1]->rows, (int)(*k));
    
    if(U_mat==NULL){
        free_double_mat(two_mat);
        free(arr_of_indices);
        return 1;
    }
    
    /* take the first k vectors and create the U matrix */
    generate_U_matrix(two_mat[1],  U_mat,  arr_of_indices);
    
    /*normalize each row*/
    normalize_to_U(U_mat);
    
    /*each row now represent a point, so we have n points at length k*/
    
    *plist_to_return = U_mat->vals[0];
    
    free(U_mat->vals);
    free(U_mat);
    /*
     deleteMatrix(U_mat);
    */
    free_double_mat(two_mat);
    free(arr_of_indices);
    
    return 0;

}




                            /*functions from here*/


double Norm( Point* p1, Point* p2){
    int i;
    double sum = 0;
    double* cord1 = p1->cordinates;
    double* cord2 = p2->cordinates;
    for(i=0; i < p1->length; i++){
        sum = sum + ((cord1[i] - cord2[i])*(cord1[i] - cord2[i]));
    }
    return sum;
}




double NormForD_mat( double* cord1, double* cord2, double len){
    unsigned int i;
    double sum=0;
    for(i=0; i < len; i++){
        sum = sum + ((cord1[i] - cord2[i])*(cord1[i] - cord2[i]));
    }
    return sum;
}

void matrix_array_Points(matrix* m, Point** arr){
    int i, j;
    Point* p;
    for(j=0; j<m->columns; j++){
        p = *(arr+j);
        for(i=0; i<m->rows; i++){
            p->cordinates[i] = m->vals[i][j];
            }
    }
}


double* matrix_array_eigenvalues(matrix* m){
    int i;
    double* arr = calloc(m->rows ,sizeof(double));
    if(arr==NULL){
        return NULL;
    }
    for(i=0; i<m->rows; i++){
        arr[i] = m->vals[i][i];
    }
    return arr;
}


void sort(int* array_indices, double* array_value, int len){
    int tmp_int;
    double tmp_double;
    int i, j;
    for(i=0; i<len; i++){
        for(j=1; j<len-i; j++){
            if(array_value[j]<array_value[j-1]){
                tmp_double = array_value[j];
                array_value[j] = array_value[j-1];
                array_value[j-1] = tmp_double;
                tmp_int = array_indices[j];
                array_indices[j] = array_indices[j-1];
                array_indices[j-1] = tmp_int;
            }
        }
    }
}



int find_k_from_arr(double* arr, int len){
    int i;
    int index=0;
    double sum = 0;
    double tmp;
    for(i=0; i<len; i++){
        tmp = fabs(arr[i+1]-arr[i]);
        if(tmp>sum){
            index = i+1;
            sum = tmp;
        }
    }
    return index;
}


int* new_arr_of_indices(int n){
    int i;
    int* arr_of_indices = (int*) calloc(n, sizeof(int));
    if(arr_of_indices==NULL) return NULL;
    for(i=0; i < n; i++){
        arr_of_indices[i] = i;
    }
    return arr_of_indices;
}

void generate_U_matrix(matrix* eighnvetors,  matrix* U,  int* array_indices){
    int j;
    int i;
    for(j=0; j<U->columns; j++){
        for(i=0; i<U->rows; i++){
            U->vals[i][j] = eighnvetors->vals[i][array_indices[j]];
        }
    }
}



void normalize_to_U(matrix* U){
    int i;
    int j;
    int m;
    double sum;
    double d;
    for(i=0; i< U->rows; i++){
        sum=0;
        
        for(j=0; j< U->columns; j++){
            sum+= pow(U->vals[i][j], 2);
        }
        
        if(sum==0) continue;
        
        else{
            d = sqrt(sum);
            for(m=0; m < U->columns; m++){
                U->vals[i][m] = U->vals[i][m] / d;
            }
        }
    }
}

void matrix_array_doubles(matrix* m, double* arr){
    int i;
    int j;
    for(i=0; i<m->rows; i++){
        for(j=0; j<m->columns; j++){
        arr[i*m->rows+j] = m->vals[i][j];
        }
    }
}
void free_double_mat(matrix** two_mat){
    deleteMatrix(two_mat[0]);
    deleteMatrix(two_mat[1]);
    free(two_mat);
}

void array_doubles_to_matrix(matrix* m, double* arr){
    int i;
    int j;
    for(i=0; i<m->rows; i++){
        for(j=0; j<m->columns; j++){
            m->vals[i][j]= arr[i*m->rows+j];
        }
    }
}

/*-----------kmean_fit from kmean++ project ----------------- */


        /*-------functions of kmeans_fit -------------*/


void closestCluster(cluster* arr, Point* p, long k){
    double min = Norm(p, arr[0].centroid);
    double tmp;
    long indexOfCluster=0;
    long i;
    for(i=0; i<k; i++){
        tmp = Norm((arr[i].centroid), p);
        if(tmp<min){
            indexOfCluster=i;
            min=tmp;
        }
    }
    arr[indexOfCluster].numberOfPoints++;
    sumPoints(arr[indexOfCluster].sumPoints, p);
}


void deleteAll(Node* head, cluster* arr, long k){
    
    deletePlist(head);
    
    deleteArrOfCluster(arr, k);
}

/*  -----------   end of functions -------------- */

int kmean_fit(double* megaListpoints, long numOfPoints, long pointSize, long k, int maxIter, double epsilon, double** plist) {

    double d;
    
    long countIter = 0;
    
    int check;
    
    long g;
    
    long i;
    
    long j;
    
    unsigned int psize;
    
    Node* head = NULL;
    
    Node* tmp;
    
    Point* pnt;
    
    Point* p = NULL;
    
    cluster* arrOfclusters = NULL;
    
    cluster* clstr;
    
    double* pplist;

    long countConverge=0;
    
    psize = (unsigned int)pointSize;
    
    for(i=0; i<numOfPoints; i++){
        check = pushP(&head);
        if(check==1){
            deletePlist(head);
            printf("An Error Has Occurred");
            return 1;
        }
        head->point = newPoint(psize);
        if(head->point==NULL){
            deletePlist(head);
            printf("An Error Has Occurred");
            return 1;
        }
        p = head->point;
        for(j=0; j<psize; j++){
            p->cordinates[j]=megaListpoints[i*psize+j];
        
    }
    
    }
    
    
    createArr(&(arrOfclusters), k);
    
    if(arrOfclusters==NULL){
        deletePlist(head);
        printf("An Error Has Occurred");
        return 1;
        
    }
    
    tmp= head;

    for(i=0; i<numOfPoints; i++){
        if( i >= numOfPoints-k){
            g = numOfPoints-i-1;
            clstr = &(arrOfclusters[g]);
            check = newCluster(clstr, psize);
            
            if(check==1){
            
                deleteAll(head, arrOfclusters, k);
                free(arrOfclusters);
                printf("An Error Has Occurred");
                return 1;
            }
          
            copyCordinates(clstr->centroid, tmp->point);
            clstr->numberOfPoints = 0;
        }
        tmp = tmp->next;
    }

    while (countConverge!=k&&(countIter<maxIter)) {
        
        countConverge=0;

        tmp = head;
        
        while(tmp!=NULL){
        
            closestCluster(arrOfclusters, tmp->point, k);
            tmp = tmp->next;

        }
    
        for(i=0; i<k; i++){
            for(j=0; j<p->length; j++){
                pnt = arrOfclusters[i].sumPoints;
                pnt->cordinates[j] =  (pnt->cordinates[j])/ (arrOfclusters[i].numberOfPoints);
            }
            
            
            arrOfclusters[i].numberOfPoints =0;
            d= Norm(arrOfclusters[i].centroid, arrOfclusters[i].sumPoints);
            if( d < epsilon ){
                
                countConverge++;
            }
            copyCordinates(arrOfclusters[i].centroid, arrOfclusters[i].sumPoints);
            resetPoint(arrOfclusters[i].sumPoints);
        }
        countIter++;
    }
    
    pplist = (double*)malloc((k*psize)*sizeof(double));
    if(pplist==NULL){
        deleteAll(head, arrOfclusters, k);
        free(arrOfclusters);
        printf("An Error Has Occurred");
        return 1;
    }
    *plist = pplist;
    for(i=0; i<k; i++){
        for(j=0; j<psize; j++){
            pnt = arrOfclusters[i].centroid;
            pplist[i*psize+j] = pnt->cordinates[j];
        }
    }
    
    deleteAll(head, arrOfclusters, k);
    free(arrOfclusters);

    return 0;
    
}

            /* ------ end of kmeans ++ --------*/



/* --------------------main-------------------------------------*/

int main(int argc, const char * argv[]) {
        
    int goal = 0;
    
    char firstLetter;
    
    char c;
    
    long k;
    
    FILE *ifp;
    
    double d;
    
    long countPoints = 0;
    
    unsigned int trackCordinates=0;
    
    long pointSize = 0;
    
    int check;
    
    long i;
    
    long j;
    
    Node* head = NULL;
    
    Node* tmp;
    
    Point* p;

    double* megaListpoints;
    
    double* arr_of_doubles;
    
    double* plist = NULL;
    
    ifp = fopen(argv[2] , "r");
   
    k = 1;
    sscanf(argv[1], "%s", &firstLetter);
    
    if(firstLetter == 'w') goal = 1;
    
    else if(firstLetter == 'd') goal = 2;
    
    else if(firstLetter == 'l') goal = 3;
    
    else if(firstLetter == 'j') goal = 4;
    
    else{
        printf("Invalid Input!\n");
        
        return 1;
    }
        
    if (argc!= 3) {
        printf("Invalid Input!\n");
        
        return 1;
    }
   
    
    while(fscanf(ifp, "%lf", &d)==1){
        check = push(&head, d);
        if(check==1){
            deleteL(head);
            fclose(ifp);
            printf("An Error Has Occurred");
            return 1;
        }
        ++pointSize;
        ++trackCordinates;
        c= getc(ifp);
        if(c=='\n'||c==EOF){
            break;
        }
    }
    
    p = newPoint((unsigned int)pointSize);
    if(p==NULL){
        deleteL(head);
        fclose(ifp);
        printf("An Error Has Occurred");
        return 1;
    }
    
    tmp = head;
    for(i=0; i<pointSize; i++){
        (p)->cordinates[pointSize-1-i] = tmp->data;
        tmp= tmp->next;
    }
    
    deleteL(head);
    
    head = NULL;
    
    check = pushP(&head);
    if(check==1){
        fclose(ifp);
        deletePoint(p);
        printf("An Error Has Occurred");
        return 1;
    }
    
    head->point = p;
    
    countPoints++;
    

    while(fscanf(ifp, "%lf", &d)==1){
        
        if(trackCordinates == pointSize){
            trackCordinates = 0;
            check = pushP(&head);
            if(check==1){
                deletePlist(head);
                fclose(ifp);
                printf("An Error Has Occurred");
                return 1;
            }
            head->point = newPoint((unsigned int)pointSize);
            if(head->point==NULL){
                deletePlist(head);
                fclose(ifp);
                printf("An Error Has Occurred");
                return 1;
            }
            p = head->point;
            countPoints++;
        }
        c = getc(ifp);
        p->cordinates[trackCordinates]= d;
        trackCordinates++;
    }
    fclose(ifp);
    
    megaListpoints = calloc((countPoints*pointSize), sizeof(double));

    if(megaListpoints==NULL){
        deletePlist(head);
        printf("An Error Has Occurred");
        return 1;
    }
    
    tmp=head;
    for(i=0; i<countPoints; i++){
        long tmpl;
        double tmpd;
        arr_of_doubles = (tmp->point)->cordinates;
        for(j=0; j<pointSize; j++){
            tmpl = countPoints*pointSize-pointSize*(i+1)+j;
            tmpd = arr_of_doubles[j];
            megaListpoints[tmpl] = tmpd;
        }
        tmp = tmp->next;
    }
    
    deletePlist(head);
    
    
    check = kmean_new_fit(&plist ,megaListpoints, &countPoints, &pointSize, &k, goal);
    
    free(megaListpoints);
    
    if(plist==NULL){
        printf("An Error Has Occurred");
        return 1;
    }
    if(check==1){
        free(plist);
        printf("An Error Has Occurred");
        return 1;
    }
    
    for(i=0; i<countPoints; i++){
        for(j=0; j<pointSize; j++){
            printf("%0.4f", plist[i*pointSize+j]);
            if(j != pointSize-1){
                printf("%c", ',');
            }
        }
            printf("%c", '\n');
    }
    free(plist);
    return 1;
}
    


/* functions of main */


int isNumber( const char number[]) {
    
    int i = 0;

    if (number[0] == '-')
        i = 1;
    for (; number[i] != 0; i++)
    {
        if (!(isdigit(number[i])))
 
            return 0;
    }
    return 1;
}




