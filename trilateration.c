#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>
#include<gsl/gsl_blas.h>
#include<lapacke.h>

// To compile: gcc trilateration.c -o trilateration -llapacke -llapack -lgsl -lgslcblas -lm

double epsl = 0.01;

typedef struct _vertex {
    int id;
    long double x, y, z;
    int *adjacents, adCount;
} Vertex;

typedef struct _edge{
    Vertex *v1, *v2;
    long double d;
} Edge;


typedef struct _graph{
    Vertex *V;
    Edge *E;
} Graph;

void usage(char *exec) {
    printf("%s -s <count_vertex>\n", exec);
    printf("%s -s <count_vertex> <tax_edges>\n", exec);
}

void updateFile(Graph G){
    FILE *file;
    file = fopen( "data.db" , "wb" );
    if (file == NULL) {
        printf("\nerror: Não foi possível abrir o arquivo");
    }else
    {
        fwrite(&G.E, sizeof(Edge), 5, file);
        fclose(file);    
    }
}

Graph createInstance(int n, int edges){
    Graph G;
    G.V = (Vertex*) calloc(n,sizeof(Vertex));
    G.E = (Edge*) calloc((n*n-n)/2,sizeof(Edge));
    int i, j;
    for (i = 0; i < n; i++)
    {
        (G.V+i)->id = i+1;
        (G.V+i)->x = rand() % 50;
        (G.V+i)->y = rand() % 50;
        (G.V+i)->z = rand() % 50;
    }
    int a = 0;
    for (i = n-1; i >= 0; i--)
    {   
        //printf("\neee %d %d\n", n,i);
        int ad[i], count =0;
        for (j = 0; j < i; j++)
        {
            //if (i == j) continue;
            (G.E+a)->v1 = (G.V+j);
            (G.E+a)->v2 = (G.V+i);
            if (j>4 && i>4 && rand() % 100>edges){
                (G.E+a++)->d = -1;    
                continue;
            } 
            (G.E+a++)->d = sqrt( pow(((G.V+i)->x - (G.V+j)->x),2) + pow(((G.V+i)->y - (G.V+j)->y),2) + pow(((G.V+i)->z - (G.V+j)->z),2));
            //if(count >= i) printf("\nooooo%d %d", count, i);else
            //printf("\n %d %d", j, count);
            ad[count++] = j;
        }
        //system("foi");
        (G.V+i)->adjacents = (int*)calloc(count, sizeof(int));
        (G.V+i)->adCount = count;
        for (int p = 0; p < count; p++)
        {
            //printf("\n %d %d", p, count);
            (G.V+i)->adjacents[p] = ad[p];
        }

    }
    return G;
}

long double buscaD(Edge *E, int v1, int v2, int n){
    for (int i = 0; i < (n*n-n)/2; i++)
    {   
        if ((((E+i)->v1 && (E+i)->v1->id == v1+1) && ((E+i)->v2 && (E+i)->v2->id == v2+1)) || (((E+i)->v1 && (E+i)->v1->id == v2+1) && ((E+i)->v2 && (E+i)->v2->id == v1+1)))
        {   
            return (E+i)->d;
        }
    }
    printf("\nERRORRR!!!!%d %d\n", v1+1, v2+1);
    return -1;
}

 /* Auxiliary routine: printing a matrix */
 void print_matrix_rowmajor( char* desc, lapack_int m, lapack_int n, double* mat, lapack_int ldm ) {
         lapack_int i, j;
         printf( "\n %s\n", desc );
  
         for( i = 0; i < m; i++ ) {
                 for( j = 0; j < n; j++ ) printf( " %6.2f", mat[i*ldm+j] );
                 printf( "\n" );
         }
 }

int main(int argc, char **argv){
    int n = 100, edges = 100;

    if (argc < 3) {
        //usage(argv[0]);
    }
    else {
        if (!strcmp(argv[1], "-s")) {
            int a, b;
            if ((a = atoi(argv[2])) == 0){
                perror("count_vertex param error: ");
                return 0;
            }
            n = a;
            if (argc > 3) {
                if ((b = atoi(argv[3])) <0){
                    perror("count_vertex param error: ");
                    return 0;
                }
                edges = b;
            }    
        } else {
            usage(argv[0]);
        }
    }

    int a[] = {5, 6, 7,10,15,20,30,40,50,100,200,400,500,600,700,800};
    for (int asd = 0; asd < 13; asd++)
    {
    n = a[asd];        
    
    //clock_t t; 
    //t = clock(); 

    Graph G = createInstance(n, edges);

    //t = clock() - t; 
    //double time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds 
  
    //printf("\nTime to preprocesor: %f\n", time_taken); 


     //for (int i = 0; i < n; i++)
    //{
       //printf("opa %d => (%Lf,%Lf,%Lf)\n", (G.V+i)->id,(G.V+i)->x,(G.V+i)->y,(G.V+i)->z);   
       //printf("(%Lg, & %Lg, & %Lg)\\\\\n", (G.V+i)->x,(G.V+i)->y,(G.V+i)->z);    
    //}
    //for (int i = 0; i < (n*n-n)/2; i++)
    //{  
       //if ((G.E+i)->v1) printf("opa (%d, %d) = %Lf\n", (G.E+i)->v1->id,(G.E+i)->v2->id, (G.E+i)->d );    
    //}
    /*
    
    long double *a = (long double *)calloc(n*n, sizeof(long double));

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            a[i*n +j] = (i == j) ? 0 : buscaD(G.E,i,j,n);   
        }
    }

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            printf(" %Lf &",a[i*n +j]);
            //*a[i,j] = sqrt( (long double)pow(((G.V+i)->x - (G.V+j)->x),2) + pow(((G.V+i)->y - (G.V+j)->y),2) + pow(((G.V+i)->z - (G.V+j)->z),2));   
        }
        printf("\\\\\n");
    } */

    Vertex vertex[n];
    
    for (int i = 0; i < 4; i++)
    {
        vertex[i].id = (G.V+i)->id;
        vertex[i].x = (G.V+i)->x;
        vertex[i].y = (G.V+i)->y;
        vertex[i].z = (G.V+i)->z;
    }

    clock_t t = clock(); 
    int countSearch = 0, countTrilateration = 0;

    for (int i = 4; i < n; i++)
    {      
        int W[4] = {0,0,0,0};
        /*int neigh = 0, aux = 0;
        while (neigh < 4)
        {        
            if (aux >= i)
            {
                printf("Erro de K-lateracao v: %d, buscas: %d, encontrados: %d",i+1, aux, neigh );
                exit(1);
            }
            if (buscaD(G.E,i,i-(++aux),n) > 0){
                W[neigh++] = i-aux;
            }

        }
        if ((G.V+i)->adCount < 4)
        {
            printf("Erro de K-lateracao v: %d, buscas: %d, encontrados: %d",i+1, aux, neigh );
            exit(1);
        }
        countSearch += aux;
        
        */
        W[0] = (G.V+i)->adjacents[(G.V+i)->adCount-1];
        W[1] = (G.V+i)->adjacents[(G.V+i)->adCount-2];
        W[2] = (G.V+i)->adjacents[(G.V+i)->adCount-3];
        W[3] = (G.V+i)->adjacents[(G.V+i)->adCount-4];
        

        //printf("%d %d %d %d\n", W[0], W[1], W[2], W[3]);

        lapack_int d, nrhs, lda, ldb, info;
        lapack_int *ipiv;
        d = 3; nrhs = 1;
        lda=d, ldb=nrhs;
        ipiv = (lapack_int *)calloc(d,sizeof(lapack_int)) ;
        if (ipiv==NULL){ printf("error of memory allocation\n"); exit(0); }
        double A[] = {2*(vertex[W[3]].x - vertex[W[0]].x),2*(vertex[W[3]].y - vertex[W[0]].y),2*(vertex[W[3]].z - vertex[W[0]].z),
                      2*(vertex[W[2]].x - vertex[W[0]].x),2*(vertex[W[2]].y - vertex[W[0]].y),2*(vertex[W[2]].z - vertex[W[0]].z),
                      2*(vertex[W[1]].x - vertex[W[0]].x),2*(vertex[W[1]].y - vertex[W[0]].y),2*(vertex[W[1]].z - vertex[W[0]].z)};

        double b[] = {  pow(vertex[W[3]].x,2)+pow(vertex[W[3]].y,2)+pow(vertex[W[3]].z,2) - (pow(vertex[W[0]].x,2)+pow(vertex[W[0]].y,2)+pow(vertex[W[0]].z,2)) -pow(buscaD(G.E,W[3],i,n),2) +pow(buscaD(G.E,W[0],i,n),2),
                        pow(vertex[W[2]].x,2)+pow(vertex[W[2]].y,2)+pow(vertex[W[2]].z,2) - (pow(vertex[W[0]].x,2)+pow(vertex[W[0]].y,2)+pow(vertex[W[0]].z,2)) -pow(buscaD(G.E,W[2],i,n),2) +pow(buscaD(G.E,W[0],i,n),2),
                        pow(vertex[W[1]].x,2)+pow(vertex[W[1]].y,2)+pow(vertex[W[1]].z,2) - (pow(vertex[W[0]].x,2)+pow(vertex[W[0]].y,2)+pow(vertex[W[0]].z,2)) -pow(buscaD(G.E,W[1],i,n),2) +pow(buscaD(G.E,W[0],i,n),2)};

        countTrilateration++;
        
        //print_matrix_rowmajor( "Details of LU factorization", d, d, A, lda );
         
        info = LAPACKE_dgesv( LAPACK_ROW_MAJOR, d, nrhs, A, lda, ipiv, b, ldb );

        if( info > 0 ) {
                printf( "The diagonal element of the triangular factor of A,\n" );
                printf( "U(%i,%i) is zero, so that A is singular;\n", info, info );
                printf( "the solution could not be computed.\n" );
                exit( 1 );
        }
        if (info <0) exit( 1 );

        vertex[i].id = i+1;
        vertex[i].x = b[0];
        vertex[i].y = b[1];
        vertex[i].z = b[2];

     //  printf("uhu %d => (%Lf,%Lf,%Lf)\n", vertex[i].id,vertex[i].x,vertex[i].y,vertex[i].z);   
     //  printf("uhu %d => (%Lf,%Lf,%Lf)\n", (G.V+i)->id,(G.V+i)->x,(G.V+i)->y,(G.V+i)->z);   
    }

    t = clock() - t; 
    double time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds 
  
    //printf("\nTime: %f\n", time_taken); 
    printf("%f, ", time_taken); 


    long double err = 0;
    for (int i = 0; i < n; i++)
    {
        err += sqrt(pow((G.V+i)->x - vertex[i].x,2)+pow((G.V+i)->y - vertex[i].y,2)+pow((G.V+i)->z - vertex[i].z,2)); 
    }    
    err = err/n;
    //printf("Erro médio: %Le\n", err);
    //printf("Quantidade de Trilaterações: %d\n", countTrilateration);
    //printf("Quantidade de iterações de busca em W: %d\n", countSearch);
    //printf("Média de iterações de busca: %g\n", (float)countSearch/countTrilateration);

/*
    long double *a = (long double *)calloc(n*n, sizeof(long double));
    long double *b = (long double *)calloc(n*n, sizeof(long double));

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            a[i*n +j] = sqrt(pow(((G.V+i)->x - (G.V+j)->x),2) + pow(((G.V+i)->y - (G.V+j)->y),2) + pow(((G.V+i)->z - (G.V+j)->z),2));   
            b[i*n +j] = sqrt(pow((vertex[i].x - vertex[j].x),2) + pow((vertex[i].y - vertex[j].y),2) + pow((vertex[i].z - vertex[j].z),2));   
        }
    }

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            printf(" %Lf ",a[i*n +j]);
            //*a[i,j] = sqrt( (long double)pow(((G.V+i)->x - (G.V+j)->x),2) + pow(((G.V+i)->y - (G.V+j)->y),2) + pow(((G.V+i)->z - (G.V+j)->z),2));   
        }
        printf("\n");
    }
        printf("\n\n");

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            printf(" %Lf ",b[i*n +j]);
            //*a[i,j] = sqrt( (long double)pow(((G.V+i)->x - (G.V+j)->x),2) + pow(((G.V+i)->y - (G.V+j)->y),2) + pow(((G.V+i)->z - (G.V+j)->z),2));   
        }
        printf("\n");
    }
    */

    int count = 0;
    long double dij = 0, mde;
    for (int i = 0; i < n; i++)
    {
        for (int j = i+1; j < n; j++)
        {
            long double dij2= (pow(((G.V+i)->x - (G.V+j)->x),2) + pow(((G.V+i)->y - (G.V+j)->y),2) + pow(((G.V+i)->z - (G.V+j)->z),2)); 
            if (dij2>0){
                dij = dij + pow((pow(vertex[i].x - vertex[j].x, 2) + pow(vertex[i].y - vertex[j].y, 2) + pow(vertex[i].z - vertex[j].z, 2) - dij2),2)/dij2;
                count++;
            }
        }
    }
    mde =  dij/count;

    //printf("%Le, ", mde); 
    //printf("Mde = %Le\n",mde);
    }
}