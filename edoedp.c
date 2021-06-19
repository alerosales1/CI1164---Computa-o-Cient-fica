
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>

#define MAXIT 50  // número máximo de iteração
#define ERRO  0.000003
#define L     6.0
#define W     8.0
#define PI    3.14159265359

//Matriz tridiagonal
typedef struct
 {
  double *x;          // pontos da malha
  double *D,*Di,*Ds;  // diagonais principal, inferior e superior
  double *B;          // vetor dos termos independentes
  double *Y;          // vetor solução
  int n;              // número de linhas e colunas da matriz
 } TRIDIAG;

// EDO:  y'' + p(x)y' + q(x)y = r(x)
typedef struct
 {
  int n;        // número de pontos internos da malha
  double a,b;   // intervalo [a,b]
  double ya,yb; // condições de contorno: y[0] = ya   y[n+1] = yb
  double (*p)(double x),(*q)(double x),(*r)(double x); // ponteiros sobre as funções p(x),q(x),r(x)
  double h;     // largura do passo da malha
  double erro;  // tolerância de erro
 } EDO;

//Matriz pentadiagonal
typedef struct
 {
  double *D,*Di1,*Di2,*Ds1,*Ds2;  // diagonais principal, inferiores e superiores
  double *B;          // vetor dos termos independentes
  double *Y;          // vetor solução
  int n,m;            // número de pontos internos da malha : n pontos em X, m pontos em Y 
  int size;           // tamanho da matriz size = n*m
 } PENTADIAG;

// EDP:  Txx + Tyy + kT = r(x,y)
typedef struct
 {
  int n,m;      // número de pontos internos da malha : n pontos em X e m pontos em Y
  double xa,xb; // intervalo [xa,xb]
  double ya,yb; // intervalo [ya,yb]
  double k;     // coeficiente constante de T(x,y) na EDP
  double (*r)(double x,double y); // ponteiro nas função r(x,y)
  double (*Txa)(double y); // condição de contorno T(xa,y)
  double (*Txb)(double y); // condição de contorno T(xb,y)
  double (*Tya)(double x); // condição de contorno T(x,ya)
  double (*Tyb)(double x); // condição de contorno T(x,yb)
  double hx,hy; // largura do passo da malha em X e em Y
  double erro;  // tolerância de erro
 } EDP; // condições de contorno:

double timestamp(void)
 {
  struct timeval tp;
  gettimeofday(&tp, NULL);
  return((double)(tp.tv_sec*1000.0 + tp.tv_usec/1000.0));
 }

static inline double P1(double x)
 {
  return 0;
 }

static inline double Q1(double x)
 {
  return 0;
 }

static inline double R1(double x)
 {
  return x*(6.0 - 0.5*x);
 }

static inline double P2(double x)
 {
  return 0;
 }

static inline double Q2(double x)
 {
  return 1;
 }

static inline double R2(double x)
 {
  return 0;
 }

static inline double P3(double x,double y)
 {
  return -1.0;
 }

static inline double R3(double x,double y)
 {
  return sin(x)*sin(x);
 }

static inline double R4(double x,double y)
 {
  return -2.0*cos(x)*cos(y);
 }

static inline double Txa3(double y)
 {
  return 20.0;
 }

static inline double Txb3(double y)
 {
  return 45.0;
 }

static inline double Tya3(double x)
 {
  return 0.0;
 }

static inline double Tyb3(double x)
 {
  return 100.0;
 }

static inline double Txa4(double y)
 {
  return cos(y);
 }

static inline double Txb4(double y)
 {
  return -cos(y);
 }

static inline double Tya4(double x)
 {
  return cos(x);
 }

static inline double Tyb4(double x)
 {
  return 0.0;
 }

void geraEDP(EDP *edp,int ind,int n,int m,int xa,int xb,int ya,int yb,double erro)
 {
  edp->n  = n;    // número de pontos dentro do intervalo em x
  edp->m  = m;    // número de pontos dentro do intervalo em y
  edp->xa = xa;   // limite inferior do intervalo em x: [xa,xb]
  edp->xb = xb;   // limite superior do intervalo em x: [xa,xb]
  edp->ya = ya;   // limite inferior do intervalo em y: [ya,yb]
  edp->yb = yb;   // limite superior do intervalo em y: [ya,yb]
  edp->hx  = (xb - xa) / (n + 1.0); // largura do passo da malha em x
  edp->hy  = (yb - ya) / (m + 1.0); // largura do passo da malha em y
  edp->erro = erro; // erro tolerado

  // Txx + Tyy + kT = r(x,y)

  if ( ind == 3 )
   {
    edp->k = -1.0;   // coeficiente constante de T(x,y) na EDP
    edp->r   = R3;   // função r(x)
    edp->Txa = Txa3; // condição de contorno T(xa,y)
    edp->Txb = Txb3; // condição de contorno T(xb,y)
    edp->Tya = Tya3; // condição de contorno T(x,ya)
    edp->Tyb = Tyb3; // condição de contorno T(x,yb)
   }
  else
   {
    edp->k   = 0.0;  // coeficiente constante de T(x,y) na EDP
    edp->r   = R4;   // função r(x)
    edp->Txa = Txa4; // condição de contorno T(xa,y)
    edp->Txb = Txb4; // condição de contorno T(xb,y)
    edp->Tya = Tya4; // condição de contorno T(x,ya)
    edp->Tyb = Tyb4; // condição de contorno T(x,yb)
   }
 }

void geraEDO(EDO *edo,int ind,int n,int a,int ya,int b,int yb,double erro)
 {
  edo->n  = n;    // número de pontos dentro do intervalo
  edo->a  = a;    // limite inferior do intervalo [a,b]
  edo->b  = b;    // limite superior do intervalo [a,b]
  edo->ya = ya;   // valor da função no limite inferior  y(x=a) = ya 
  edo->yb = yb;   // valor da função no limite superior  y(x=b) = yb
  edo->h  = (b - a) / (n + 1.0); // largura do passo da malha
  edo->erro = erro;

  // y'' + p(x)y' + q(x)y = r(x)

  if ( ind == 0 )
   {
    edo->p = P1; // função p(x)
    edo->q = Q1; // função q(x)
    edo->r = R1; // função r(x)
   }
  else
   {
    edo->p = P2; // função p(x)
    edo->q = Q2; // função q(x)
    edo->r = R2; // função r(x)
   }
 }
 
void geraPentadiag(EDP *edp,PENTADIAG *SL)
 {
  int i,n,m,size;
  double xa,xb,ya,yb,hx,hy;

  n = edp->n;
  m = edp->m;
  size = n*m;
  SL->n = n;
  SL->m = m;
  SL->size = m*n;

  SL->Di2 = (double*) calloc (size,sizeof(double)); // 2a diagonal inferior
  SL->Di1 = (double*) calloc (size,sizeof(double)); // 1a diagonal inferior
  SL->D   = (double*) calloc (size,sizeof(double)); // diagonal principal
  SL->Ds1 = (double*) calloc (size,sizeof(double)); // 1a diagonal superior
  SL->Ds2 = (double*) calloc (size,sizeof(double)); // 2a diagonal superior
  SL->B   = (double*) calloc (size,sizeof(double)); // vetor dos coeficientes independentes
  SL->Y   = (double*) calloc (size,sizeof(double)); // vetor solução

  xa = edp->xa;
  xb = edp->xb;
  ya = edp->xa;
  yb = edp->xb;
  hx = edp->hx;
  hy = edp->hy;

  for ( i = 0 ; i < size ; i++ )
   {
    SL->D[i]   = -2.0*(hx*hx+hy*hy) + hx*hx*hy*hy*edp->k; // diagonal principal
    SL->Ds1[i] = hy*hy; // 1a diagonal superior
    SL->Ds2[i] = hx*hx; // 2a diagonal superior
    SL->Di1[i] = hy*hy; // 1a diagonal inferior
    SL->Di2[i] = hx*hx; // 2a diagonal inferior
    SL->Y[i]   = 0.0;   // inicialização vetor solução
   }

  // esses termos são definidos pelas condições de contorno
  SL->Ds1[n-1]       = 0.0;
  SL->Ds1[(m-1)*n-1] = 0.0;
  SL->Di1[n]         = 0.0;
  SL->Di1[(m-1)*n]   = 0.0;
 
  // condições de contorno substraidas ao vetor dos coeficientes independentes
  SL->B[0] = hx*hx*hy*hy*edp->r(xa+hx,ya+hy) - hx*hx*edp->Tya(xa+hx) - hy*hy*edp->Txa(ya+hy);

  for ( i = 1 ; i < n-1 ; i++ ) SL->B[i] =  
    hx*hx*hy*hy*edp->r(xa+(i+1)*hx,ya+hy) - hx*hx*edp->Tya(xa+(i+1)*hx);

  SL->B[n-1] = hx*hx*hy*hy*edp->r(xa+n*hx,ya+hy) - hx*hx*edp->Tya(xa+n*hx) - hy*hy*edp->Txb(ya+hy);

  SL->B[n] = hx*hx*hy*hy*edp->r(xa+hx,ya+(m-1)*hy) - hy*hy*edp->Txa(ya+(m-1)*hy);

  for ( i = n+1 ; i < (m-1)*n-1 ; i++ ) SL->B[i] = hx*hx*hy*hy*edp->r(xa+(i-n+1)*hx,ya+(m-1)*hy);

  SL->B[(m-1)*n-1] = hx*hx*hy*hy*edp->r(xa+n*hx,ya+(m-1)*hy) - hy*hy*edp->Txb(ya+(m-1)*hy);  

  SL->B[(m-1)*n] = hx*hx*hy*hy*edp->r(xa+hx,ya+m*hy) - hx*hx*edp->Tyb(xa+hx)
                   - hy*hy*edp->Txa(ya+m*hy);

  for ( i = (m-1)*n+1 ; i < m*n-1 ; i++ ) 
    SL->B[i] = hx*hx*hy*hy*edp->r(xa+(i-(m-1)*n+1)*hx,ya+m*hy) - hx*hx*edp->Tyb(xa+(i-(m-1)*n+1)*hx);

  SL->B[m*n-1] = hx*hx*hy*hy*edp->r(xa+n*hx,ya+m*hy) - hx*hx*edp->Tyb(xa+n*hx) 
                 - hy*hy*edp->Txb(ya+m*hy);
 }

void geraTridiag(EDO *edo,TRIDIAG *SL)
 {
  double xi,h,a,b;
  int i,n;

  n = edo->n;
  a = edo->a;
  b = edo->b;
  h = edo->h;
  SL->n = edo->n;

  SL->x  = (double*) calloc (n,sizeof(double)); // valores de x nos pontos da malha
  SL->Di = (double*) calloc (n,sizeof(double)); // diagonal inferior
  SL->D  = (double*) calloc (n,sizeof(double)); // diagonal principal
  SL->Ds = (double*) calloc (n,sizeof(double)); // diagonal superior
  SL->B  = (double*) calloc (n,sizeof(double)); // vetor dos coeficientes independentes
  SL->Y  = (double*) calloc (n,sizeof(double)); // vetor solução
 
  // i = 0
  xi = a+h;
  SL->x[0]  =  xi;
  SL->D[0]  = -2.0 + h*h*edo->q(xi);   // diagonal principal
  SL->Ds[0] =  1.0 + 0.5*h*edo->p(xi); // diagonal superior
  SL->B[0]  =  h*h*edo->r(xi) - edo->ya*(1.0 - 0.5*h*edo->p(xi)); // termo independente  
  SL->Y[0]  =  0.0;  // inicialização do vetor solução

  for ( i = 1 ; i < n-1 ; i++ )
   {
    xi = a + (i+1)*h;  // ponto da malha
    SL->x[i]  =  xi;
    SL->Di[i] =  1.0 - 0.5*h*edo->p(xi); // diagonal inferior
    SL->D[i]  = -2.0 + h*h*edo->q(xi);   // diagonal principal
    SL->Ds[i] =  1.0 + 0.5*h*edo->p(xi); // diagonal superior
    SL->B[i]  =  h*h*edo->r(xi);         // termo independente  
    SL->Y[i]  =  0.0; // inicialização do vetor solução
   }

  i = n-1;
  xi = b-h;
  SL->x[i]  =  xi;
  SL->Di[i] =  1.0 - 0.5*h*edo->p(xi); // diagonal inferior
  SL->D[i]  = -2.0 + h*h*edo->q(xi);   // diagonal principal
  SL->B[i]  =  h*h*edo->r(xi) - edo->yb*(1.0 + 0.5*h*edo->p(xi)); // termo independente  
  SL->Y[i]  =  0.0; // inicialização do vetor solução
 }

void freeTridiagonal (TRIDIAG *SL)
 {
  free(SL->Di);
  free(SL->D);
  free(SL->Ds);
  free(SL->B);
  free(SL->Y);
 }

void freePentadiagonal (PENTADIAG *SL)
 {
  free(SL->Di1);
  free(SL->Di2);
  free(SL->D);
  free(SL->Ds1);
  free(SL->Ds2);
  free(SL->B);
  free(SL->Y);
 }

double norma (double *x0,double *x1, int n)
 {
  double max,diff;

  max = 0.0;

  for (int i = 0; i < n; i++)
   {
    diff = fabs(x1[i] - x0[i]);

    if (diff > max) max = diff; // máximum das diferenças
   }

  return max;
 }

// A*Y = B
// A = D + L + U
// D : matriz somente com a diagonal
// L : matriz triangular inferior
// U : matriz triangular superior
// (L + D)*YY = B - U*Y
// D*YY = B - L*YY - U*Y
// YY = (B - L*YY - U*Y)/D

int GaussSeidel_EDP (EDP *edp,PENTADIAG *SL,double *tTotal)
 {
  *tTotal = timestamp();
  int n,m,i,ni,continuar,size;
  double *Y,*YY,E;

  n = edp->n;
  m = edp->m;
  E = edp->erro;
  size = SL->size;
  Y = SL->Y;

  YY = (double *) malloc(sizeof(double)*size);

  for ( i = 0 ; i < size ; i++ )  // vetor solução inicial nulo
   {
    Y[i]  = 0.0; 
    YY[i] = 0.0;
   }

  continuar = 1;
  ni = 0;

  while (continuar && ni < MAXIT)
   {
    YY[0] = (SL->B[0] - SL->Ds1[0]*Y[1] - SL->Ds2[0]*Y[n])/SL->D[0];

    for ( i = 1; i < n; i++ )
      YY[i] = (SL->B[i] - SL->Di1[i]*YY[i-1] - SL->Ds1[i]*Y[i+1] - SL->Ds2[i]*Y[i+n])/SL->D[i];

    for ( i = n; i < 2*n; i++ )
      YY[i] = (SL->B[i] - SL->Di2[i]*YY[i-n] - SL->Di1[i]*YY[i-1]
                        - SL->Ds1[i]*Y[i+1]  - SL->Ds2[i]*Y[i+n])/SL->D[i];

    for ( i = 2*n; i < m*n-1; i++ )
      YY[i] = (SL->B[i] - SL->Di2[i]*YY[i-n] - SL->Di1[i]*YY[i-1] - SL->Ds1[i]*Y[i+1])/SL->D[i];

    YY[m*n-1] = (SL->B[m*n-1] - SL->Di2[m*n-1]*YY[m*n-1-n] - SL->Di1[m*n-1]*YY[m*n-2])/SL->D[m*n-1];

    // se | YY - Y | < Tolerância
    if ( norma(YY,Y,size) < E ) continuar = 0;  // temos uma solução
    else
     {
      for ( i = 0 ; i < size ; i++ ) Y[i] = YY[i];   
      ni++;
     }
   }
 
  free(YY);
  *tTotal = timestamp() - *tTotal;
  return ni;
 }

// A*Y = B
// A = D + L + U
// D : matriz somente com a diagonal
// L : matriz triangular inferior
// U : matriz triangular superior
// (L + D)*YY = B - U*Y
// D*YY = B - L*YY - U*Y
// YY = (B - L*YY - U*Y)/D
   
int GaussSeidel_EDO (EDO *edo,TRIDIAG *SL,double *tTotal)
 {
  *tTotal = timestamp();
  int n,i,ni,continuar;
  double *Y,*YY,E;

  n = edo->n;
  E = edo->erro;
  Y = SL->Y;

  YY = (double *) malloc(sizeof(double)*n);

  for ( i = 0 ; i < n ; i++ )  // vetor solução inicial nulo
   {
    Y[i]  = 0.0; 
    YY[i] = 0.0;
   }

  continuar = 1;
  ni = 0;

  while (continuar && ni < MAXIT)
   {
    YY[0] = (SL->B[0] - SL->Ds[0]*Y[1]) / SL->D[0];

    for ( i = 1; i < n-1; i++ )
      YY[i] = (SL->B[i] - SL->Di[i]*YY[i-1] - SL->Ds[i]*Y[i+1])/SL->D[i];

    YY[n-1] = (SL->B[n-1] - SL->Di[n-1]*YY[n-2])/SL->D[n-1];

    // se | YY - Y | < Tolerância
    if ( norma(YY,Y,n) < E ) continuar = 0;  // temos uma solução
    else
     {
      for ( i = 0 ; i < n ; i++ ) Y[i] = YY[i];   
      ni++;
     }
   }
 
  free(YY);
  *tTotal = timestamp() - *tTotal;
  return ni;
 }

double normaL2_EDO(TRIDIAG *SL)
 {
  int i,n;
  double norma;
  double *res;

  n = SL->n;
  res = (double *) malloc(sizeof(double)*n);

  res[0] = SL->B[0] - SL->D[0]*SL->Y[0] - SL->Ds[0]*SL->Y[1];

  for ( i = 1 ; i < n-1 ; i++ ) 
    res[i] = SL->B[i] - SL->Di[i]*SL->Y[i-1] - SL->D[i]*SL->Y[i] - SL->Ds[i]*SL->Y[i+1];

  res[n-1] = SL->B[n-1] - SL->Di[n-1]*SL->Y[n-2] - SL->D[n-1]*SL->Y[n-1];  

  norma = 0.0;
  for ( i = 0 ; i < n ; i++ ) norma += res[i]*res[i];

  free(res);
  return sqrt(norma); 
 }

double normaL2_EDP(PENTADIAG *SL)
 {
  int i,n,m,size;
  double norma;
  double *res;

  n = SL->n;
  m = SL->m;
  size = SL->size;

  res = (double *) malloc(sizeof(double)*size);

  res[0] = SL->B[0] - SL->D[0]*SL->Y[0] - SL->Ds1[0]*SL->Y[1] - SL->Ds2[0]*SL->Y[n];

  for ( i = 1 ; i < n ; i++ ) 
    res[i] = SL->B[i] - SL->Di1[i]*SL->Y[i-1] - SL->D[i]*SL->Y[i]
                      - SL->Ds1[i]*SL->Y[i+1] - SL->Ds2[i]*SL->Y[i+n];

  for ( i = n ; i < (m-1)*n ; i++ ) 
    res[i] = SL->B[i] - SL->Di2[i]*SL->Y[i-n] - SL->Di1[i]*SL->Y[i-1] - SL->D[i]*SL->Y[i]
                      - SL->Ds1[i]*SL->Y[i+1] - SL->Ds2[i]*SL->Y[i+n];

  for ( i = (m-1)*n ; i < m*n-1 ; i++ ) 
    res[i] = SL->B[i] - SL->Di2[i]*SL->Y[i-n] - SL->Di1[i]*SL->Y[i-1] - SL->D[i]*SL->Y[i]
                      - SL->Ds1[i]*SL->Y[i+1];

  res[m*n-1] = SL->B[m*n-1] - SL->Di2[m*n-1]*SL->Y[(m-1)*n-1] - SL->Di1[m*n-1]*SL->Y[m*n-2]
                            - SL->D[m*n-1]*SL->Y[m*n-1];  

  norma = 0.0;
  for ( i = 0 ; i < size ; i++ ) norma += res[i]*res[i];

  free(res);
  return sqrt(norma); 
 }

void resultadoEDO(EDO *edo,int ind,TRIDIAG *SL,int ni,double tempo)
 {
  int i,j,n;
  n = SL->n;

  if (ind == 0) printf("***** item (a) Y\" = 6x - 0.5 x²: ");
  else          printf("***** item (c) Y\" + Y = 0: ");
  printf("n = %d, H = %.6f\n",n,edo->h);

  printf("\nSistema Linear:\n\n");

  printf("%5.2f %5.2f ",SL->D[0],SL->Ds[0]);
  for ( j = 2 ; j < n ; j++ ) printf(" 0.00 ");
  printf("| %9.6f\n\n",SL->B[0]);
  
  for ( i = 1 ; i < n-1 ; i++ )
   {
    for ( j = 0 ; j < i-1 ; j++) printf(" 0.00 ");
    printf("%5.2f %5.2f %5.2f ",SL->Di[i],SL->D[i],SL->Ds[i]);
    for ( j = i+2 ; j < n ; j++ ) printf(" 0.00 ");
    printf("| %9.6f\n\n",SL->B[i]);
   }

  for ( j = 0 ; j < n-2 ; j++ ) printf(" 0.00 ");
  printf("%5.2f %5.2f ",SL->Di[n-1],SL->D[n-1]);
  printf("| %9.5f\n\n",SL->B[n-1]);
 
  printf("Y: ");
  for ( i = 0 ; i < n ; i++ ) printf("%9.5f ",SL->Y[i]);
  printf("\n\n");

  printf("Norma L2: %.7e, Tempo: %.8f ms\n\n",normaL2_EDO(SL),tempo);
 }

void resultadoEDP(EDP *edp,int ind,PENTADIAG *SL,int ni,double tempo)
 {
  int i,j,n,m,size;
  n = SL->n;
  m = SL->m;
  size = SL->size;

  if (ind == 3) printf("***** item (b) Txx + Tyy - T = sin²(x) : ");
  else          printf("***** item (d) Uxx + Uyy = -cos(x+y) - cos(x-y) : ");
  printf("L = %.2f, W = %.2f, n = %d, m = %d, Hx = %.2f , Hy = %.2f\n",
          edp->xb,edp->yb,n,m,edp->hx,edp->hy);

  printf("\nSistema Linear:\n\n");

  printf("2a diagonal superior\n");
  for ( i = 0 ; i < (m-1)*n ; i++ ) printf("%.2f ",SL->Ds2[i]); // 2a diagonal superior
  printf("\n\ndiagonal superior\n");
  for ( i = 0 ; i < size-1 ; i++ )  printf("%.2f ",SL->Ds1[i]); // 1a diagonal superior
  printf("\n\ndiagonal principal\n");
  for ( i = 0 ; i < size ; i++ )    printf("%.2f ",SL->D[i]);   // diagonal principal
  printf("\n\ndiagonal inferior\n");
  for ( i = 1 ; i < size ; i++ )    printf("%.2f ",SL->Di1[i]); // 1a diagonal inferior
  printf("\n\n2a diagonal inferior\n");
  for ( i = n ; i < size ; i++ )    printf("%.2f ",SL->Di2[i]); // 2a diagonal inferior

  printf("\n\nvetor dos coeficientes independentes\n");
  for ( i = 0 ; i < size ; i++ )    printf("%9.6f ",SL->B[i]); 

  printf("\n\nvetor solução\nT : ");
  for ( i = 0 ; i < size ; i++ )    printf("%9.6f ",SL->Y[i]); 
  printf("\n\n");

  printf("Norma L2: %.7e, Tempo: %.8f ms\n\n",normaL2_EDP(SL),tempo);
  printf("%d iterações\n\n",ni);

 }

int main()
 {
  EDO edo;
  EDP edp;
  TRIDIAG   SL;
  PENTADIAG SL2;
  double tempo;
  int result;
  
  // y'' = 6x -0.5*x²
  // equação 0 ; 5 pontos ; y(0) = ya = 0 ; y(12) = yb = 0
  geraEDO(&edo,0,5,0,0,12,0,ERRO);
  geraTridiag(&edo,&SL);
  result = GaussSeidel_EDO(&edo,&SL,&tempo);
  resultadoEDO(&edo,0,&SL,result,tempo);
  freeTridiagonal(&SL);
 
  // equação 0 ; 10 pontos ; y(0) = ya = 0 ; y(12) = yb = 0
  geraEDO(&edo,0,10,0,0,12,0,ERRO);
  geraTridiag(&edo,&SL);
  result = GaussSeidel_EDO(&edo,&SL,&tempo);
  resultadoEDO(&edo,0,&SL,result,tempo);
  freeTridiagonal(&SL);

  // y'' + y = 0
  // equação 1 ; 5 pontos ; y(0) = ya = 0 ; y(1) = yb = 1
  geraEDO(&edo,1,5,0,0,1,1,ERRO);
  geraTridiag(&edo,&SL);
  result = GaussSeidel_EDO(&edo,&SL,&tempo);
  resultadoEDO(&edo,1,&SL,result,tempo);
  freeTridiagonal(&SL);

  // equação 1 ; 10 pontos ; y(0) = ya = 0 ; y(1) = yb = 1
  geraEDO(&edo,1,10,0,0,1,1,ERRO);
  geraTridiag(&edo,&SL);
  result = GaussSeidel_EDO(&edo,&SL,&tempo);
  resultadoEDO(&edo,1,&SL,result,tempo);
  freeTridiagonal(&SL);

  // Uxx + Uyy - U = sin²(x)
  // equação 3 ; 5 pontos em x; 3 pontos em y ; x:[0,L] y:[0,W]
  geraEDP(&edp,3,5,3,0,L,0,W,ERRO);
  geraPentadiag(&edp,&SL2);
  result = GaussSeidel_EDP(&edp,&SL2,&tempo);
  resultadoEDP(&edp,3,&SL2,result,tempo);
  freePentadiagonal(&SL2);

  // Uxx + Uyy - U = sin²(x)
  // equação 3 ; 3 pontos em x; 3 pontos em y ; x:[0,L] y:[0,W]
  geraEDP(&edp,3,5,3,0,L,0,W,ERRO);
  geraPentadiag(&edp,&SL2);
  result = GaussSeidel_EDP(&edp,&SL2,&tempo);
  resultadoEDP(&edp,3,&SL2,result,tempo);
  freePentadiagonal(&SL2);

  // Uxx + Uyy - U = sin²(x)
  // equação 3 ; 10 pontos em x; 3 pontos em y ; x:[0,L] y:[0,W]
  geraEDP(&edp,3,10,3,0,L,0,W,ERRO);
  geraPentadiag(&edp,&SL2);
  result = GaussSeidel_EDP(&edp,&SL2,&tempo);
  resultadoEDP(&edp,3,&SL2,result,tempo);
  freePentadiagonal(&SL2);

  // Uxx + Uyy  = -cos(x+y) - cos(x-y)
  // equação 4 ; 5 pontos em x; 3 pontos em y ; x:[0,PI] y:[0,PI/2]
  geraEDP(&edp,4,5,3,0,PI,0,0.5*PI,ERRO);
  geraPentadiag(&edp,&SL2);
  result = GaussSeidel_EDP(&edp,&SL2,&tempo);
  resultadoEDP(&edp,4,&SL2,result,tempo);
  freePentadiagonal(&SL2);

  // Uxx + Uyy  = -cos(x+y) - cos(x-y)
  // equação 4 ; 10 pontos em x; 3 pontos em y ; x:[0,PI] y:[0,PI/2]
  geraEDP(&edp,4,10,3,0,PI,0,0.5*PI,ERRO);
  geraPentadiag(&edp,&SL2);
  result = GaussSeidel_EDP(&edp,&SL2,&tempo);
  resultadoEDP(&edp,4,&SL2,result,tempo);
  freePentadiagonal(&SL2);

  return 0;
 }
