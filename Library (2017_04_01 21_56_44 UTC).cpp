/*...........................................................................................................

// ������� ��� ���, ��� ������� ������ �������� ���������, ��������� ������.
FastFT,GetN,Log2,TriDec,TriSolv,InvMatr,SortEigVal,Conv,

// ������� ��� ���������� ����������� �������� � ��������.
FindEigenValVect,Hessen,Haus,TriMatrQRShift,CheckUnderDiagElem,FindMinUnderDiagElem,StepQRShift,
EigVectTriMatrRight,EigVectTriMatrLeft,SolvSystemTri,ReadString

// ������� ������ ����������.
GetFName,

...........................................................................................................*/

#include "stdafx.h"

#include <conio.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <omp.h>

#include "Headers.h"
#include "cmplx.h"
#include "Matrix.h"

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ���������.

#define flScanf_S 0 // ���� ������������� 'scanf_s'.

//-----------------------------------------------------------------------------------------------------------
// ������� �������������� �����.
// ������������� ��������� �� ����������� ������� 'w': 'i' �� 1 �� N/2-1,
// ������������� �� ����������� ������ ������� '-w': 'i' � �������� ������� �� N-1 �� N/2+1,
// w=0 ��� i=0, ������������ 'w' ��� i=N/2.

BYTE FastFT(complex *A,complex *B,BYTE M,SCHAR dir)
{
BYTE l; int i,j,k,n,le,leW,ip,N; double d,vd; complex U,W,T;

if(A==NULL) return 1; if(B==NULL) return 2; if(M==0) return 3; if(M>=31) return 4;
N=1<<M; n=N/2;
if(dir>=0) { d=1./(double)N; for(i=0;i<N;i++) B[i]=A[i]*complex(d,0.);} else for(i=0;i<N;i++) B[i]=A[i];
j=1; for(i=0;i<N-1;i++) { if(i<j-1) { T=B[j-1]; B[j-1]=B[i]; B[i]=T;}
k=n; b: if(k>=j) goto a; j-=k; k=k/2; goto b; a: j+=k;}
for(l=0;l<M;l++) { le=1<<(l+1); leW=le/2;
U=complex(1.,0.); vd=M_PI/(double)leW; d=sin(vd); if(dir<0) d=-d; W=complex(cos(vd),d);
for(j=0;j<leW;j++) { for(i=j;i<N;i+=le) { ip=i+leW; T=B[ip]*U; B[ip]=B[i]-T; B[i]+=T;} U*=W;}}
return 0;
}

//-----------------------------------------------------------------------------------------------------------
// ���������� 'N' - 2 � ������� 'M'.

BYTE GetN(int M,int *N)
{
if(N==NULL) return 1; if(M==0||M>=32) return 2; *N=1<<M; if(*N<2) return 3; return 0;
}

//-----------------------------------------------------------------------------------------------------------
// �������� �� ��������� 2.

double Log2(double n) { return log(n)/log(2.);}

//-----------------------------------------------------------------------------------------------------------
// ����������� ���������� ������� � ������ ���������� ������ � ������� �������� �������� �� ������� � ��������.

/*
//��������� ������.
BYTE TriDec(int N,const clMatrix &A,clMatrix &T,int *NumStr,int *NumRow)
{
BYTE err; int i,iEx,j,jSh,jShv,jShp,iv,jv,iM,jM,ip,jp,szCycle2; double v,vM; complex coe1,coe2,*pM;
double Time=0.,Time0=0.;

if(N<=0) return 1;
if(A.IsOK()!=0) return 2; if(A.Nx!=N) return 3; if(A.Ny!=N) return 4;
if(T.IsOK()!=0) return 5; if(T.Nx!=N) return 6; if(T.Ny!=N) return 7;
if(NumStr==NULL) return 8; if(NumRow==NULL) return 9;

T=A; pM=T.Matr; // ���������� �������������� ������� � ������� 'T', ������� ����� ����� ������� ��������� ������������ ���������� �� ������ ������.

// ����������� ���������� �������. ��������� ������������ � �� �� ������� 'T'.
for(j=0;j<N;j++) NumStr[j]=NumRow[j]=j; // ������������� �������� ������� ����������� ����� � ��������.

// ���� �� ������, ��� ���������� ���������� �����������. ����� ������� ������ - 'iEx'.
for(iEx=0;iEx<N-1;iEx++) { 
jSh=iEx*N;
// ������� ������������ ������� �� ������� � ��������.
vM=0.; ip=jp=-1;
for(i=iEx;i<N;i++) {
iv=NumStr[i]; if(iv<0||iv>=N) return 10; // ������ ������� 'i' ���������� ����������� 'NumStr' - 'iv'.
jShv=iv*N;
for(j=iEx;j<N;j++) {
jv=NumRow[j]; if(jv<0||jv>=N) return 11; // ������ ������� 'j' ���������� ����������� 'NumRow' - 'jv'.
v=abs(pM[jShv+jv]); if(v>vM) { vM=v; ip=i; jp=j;}}}
if(vM==0.||ip<0||jp<0) return messSingMatr;

iM=NumStr[ip]; if(iM<0||iM>=N) return 12; // ������ ������� 'ip' ���������� ����������� 'NumStr' - 'iM'.
jM=NumRow[jp]; if(jM<0||jM>=N) return 13; // ������ ������� 'jp' ���������� ����������� 'NumRow' - 'jM'.
iv=NumStr[iEx]; if(iM!=iv) { NumStr[iEx]=iM; NumStr[ip]=iv;} // ������ ������������ ������� ����� � ������� 'NumStr'.
jv=NumRow[iEx]; if(jM!=jv) { NumRow[iEx]=jM; NumRow[jp]=jv;} // ������ ������������ ������� �������� � ������� 'NumRow'.

jShv=iM*N; coe1=Cmplx_1/pM[jShv+jM];

Time=omp_get_wtime();
// �������� �� ����������� ����� ������ � ������� ��������� �� ���������, ���������� �� ����������� - ������� ������ � ������� ������� �������� ��������.
err=0; szCycle2=(N-iEx-1)*(N-iEx-1);
#pragma omp parallel for if(szCycle2>(int)MinNumMultParal/10) schedule(static) private(j,i,jv,iv,jShp,coe2)
for(j=iEx+1;j<N;j++) {
jv=NumRow[j]; if(jv<0||jv>=N) err=15; // ������ ������� 'j' ���������� ����������� 'NumRow' - 'jv'.
pM[jv+jShv]*=coe1; // ����� �� ������������ ������� ���������� ����� ������ � ������� ��������� �� ���������.
for(i=iEx+1;i<N;i++) {
iv=NumStr[i]; if(iv<0||iv>=N) err=16; // ������ ������� 'i' ���������� ����������� 'NumStr' - 'iv'.
jShp=iv*N; coe2=pM[jM+jShp];
pM[jv+jShp]-=pM[jv+jShv]*coe2;}
} // ����� ����� �� 'j'.
Time=omp_get_wtime()-Time; Time0+=Time;
} // ����� ����� �� 'iEx'.
printf("In TriDec Time sum : %lf\n",Time0);
return err;
}
*/

BYTE TriDec(int N,const clMatrix &A,clMatrix &T,int *NumStr,int *NumRow)
{
BYTE err; int i,iEx,j,jSh,jShv,jShp,iv,jv,iM,jM,szCycle2; int *ip,ipM,*jp,jpM,thr_num,n_th;
double v,*vM,vM_Max; complex coe1,coe2,*pM;

if(N<=0) return 1;
if(A.IsOK()!=0) return 2; if(A.Nx!=N) return 3; if(A.Ny!=N) return 4;
if(T.IsOK()!=0) return 5; if(T.Nx!=N) return 6; if(T.Ny!=N) return 7;
if(NumStr==NULL) return 8; if(NumRow==NULL) return 9;

n_th=omp_get_max_threads();
ip=jp=NULL; vM=NULL; err=0;

ip=new int[n_th]; if(ip==NULL) { err=11; goto end;} jp=new int[n_th]; if(jp==NULL) { err=12; goto end;}
vM=new double[n_th]; if(vM==NULL) { err=13; goto end;} // ��������� ������ ��� ��������������� ��������.

T=A; pM=T.Matr; // ���������� �������������� ������� � ������� 'T', ������� ����� ����� ������� ��������� ������������ ���������� �� ������ ������.

// ����������� ���������� �������. ��������� ������������ � �� �� ������� 'T'.
for(j=0;j<N;j++) NumStr[j]=NumRow[j]=j; // ������������� �������� ������� ����������� ����� � ��������.

// ������ ������������ �������.
#pragma omp parallel default(shared) private(iEx,jSh,thr_num)
{
thr_num=omp_get_thread_num();
// ���� �� ������, ��� ���������� ���������� �����������. ����� ������� ������ - 'iEx'.
for(iEx=0;iEx<N-1;iEx++) {
jSh=iEx*N;
// ������� ������������ ������� �� ������� � ��������.
vM[thr_num]=0.; ip[thr_num]=jp[thr_num]=-1;

#pragma omp barrier
#pragma omp for private(i,iv,j,jShv,jv,v) reduction(+:err)
for(i=iEx;i<N;i++) { if(err!=0) continue;
iv=NumStr[i]; if(iv<0||iv>=N) { err=10; continue;} // ������ ������� 'i' ���������� ����������� 'NumStr' - 'iv'.
jShv=iv*N;
for(j=iEx;j<N;j++) {
jv=NumRow[j]; if(jv<0||jv>=N) { err=11; continue;} // ������ ������� 'j' ���������� ����������� 'NumRow' - 'jv'.
v=abs(pM[jShv+jv]); if(v>vM[thr_num]) { vM[thr_num]=v; ip[thr_num]=i; jp[thr_num]=j;}}}
#pragma omp barrier

if(err==0) {
#pragma omp master
{
ipM=ip[0]; jpM=jp[0]; vM_Max=vM[0];
for(i=0;i<n_th;i++) if(vM[i]>vM_Max) { ipM=ip[i]; jpM=jp[i];}
iM=NumStr[ipM]; if(iM<0||iM>=N) err=12; // ������ ������� 'ip' ���������� ����������� 'NumStr' - 'iM'.
jM=NumRow[jpM]; if(jM<0||jM>=N) err=13; // ������ ������� 'jp' ���������� ����������� 'NumRow' - 'jM'.
iv=NumStr[iEx]; if(iM!=iv) { NumStr[iEx]=iM; NumStr[ipM]=iv;} // ������ ������������ ������� ����� � ������� 'NumStr'.
jv=NumRow[iEx]; if(jM!=jv) { NumRow[iEx]=jM; NumRow[jpM]=jv;} // ������ ������������ ������� �������� � ������� 'NumRow'.

jShv=iM*N; coe1=Cmplx_1/pM[jShv+jM];
szCycle2=(N-iEx-1)*(N-iEx-1);
} // ����� 'omp master'.
} // ����� 'if'

#pragma omp barrier
if(err==0) {
// �������� �� ����������� ����� ������ � ������� ��������� �� ���������, ���������� �� ����������� - ������� ������ � ������� ������� �������� ��������.
#pragma omp for schedule(static) private(j,i,jv,iv,jShp,coe2) reduction(+:err)
for(j=iEx+1;j<N;j++) { if(err!=0) continue;
jv=NumRow[j]; if(jv<0||jv>=N) { err=15; continue;} // ������ ������� 'j' ���������� ����������� 'NumRow' - 'jv'.
pM[jv+jShv]*=coe1; // ����� �� ������������ ������� ���������� ����� ������ � ������� ��������� �� ���������.
for(i=iEx+1;i<N;i++) {
iv=NumStr[i]; if(iv<0||iv>=N) { err=16; continue;} // ������ ������� 'i' ���������� ����������� 'NumStr' - 'iv'.
jShp=iv*N; coe2=pM[jM+jShp];
pM[jv+jShp]-=pM[jv+jShv]*coe2;}
} // ����� ����� �� 'j'.
} // ����� ����� �� 'iEx'.
#pragma omp barrier
} // ����� 'if'.
} // ����� ������������ �������.
end: SAFE_DELETE_ARR(ip); SAFE_DELETE_ARR(jp); SAFE_DELETE_ARR(vM); return err;
}

//-----------------------------------------------------------------------------------------------------------
// ������� ������� �������� ��������� �� ������������ ���������� ������� � ������ ���������� ������ � ������� �������� �������� �� ������� � ��������.

BYTE TriSolv(int N,const clMatrix &T,int *NumStr,int *NumRow,complex *F,complex *X)
{
BYTE err; int i,iEx,j,iv,jv,jSh,jShv,iM,jM; complex s,*G,*pM;

if(N<=0) return 1; if(T.IsOK()!=0) return 2; if(T.Nx!=N) return 3; if(T.Ny!=N) return 4; pM=T.Matr;
if(NumStr==NULL) return 5; if(NumRow==NULL) return 6; if(F==NULL) return 7; if(X==NULL) return 8;

G=NULL; err=0;
G=new complex[N]; if(G==NULL) { err=9; goto end;} // ������ ��� �������������� ������� ��� ������� ������� ������� ������.

// ���� ���������� ������� 'G'.
for(i=0;i<N;i++) G[i]=F[i]; // �������� ������� 'F' � 'G'.
for(iEx=0;iEx<N;iEx++) {
iM=NumStr[iEx]; if(iM<0||iM>=N) { err=10; goto end;} // ������ ������� 'iEx' ���������� ����������� 'NumStr' - 'iM'.
jM=NumRow[iEx]; if(jM<0||jM>=N) { err=11; goto end;} // ������ ������� 'iEx' ���������� ����������� 'NumRow' - 'jM'.
jSh=iM*N; G[iM]/=pM[jSh+jM]; if(iEx==N-1) break;
for(i=iEx+1;i<N;i++) {
iv=NumStr[i]; if(iv<0||iv>=N) { err=12; goto end;} // ������ ������� 'iEx' ���������� ����������� 'NumStr' - 'iv'.
jShv=iv*N; G[iv]-=G[iM]*pM[jM+jShv];}
} // ����� ����� �� 'iEx'.

// ���� ���������� ������� 'X'.
for(i=N-1;i>=0;i--) {
iv=NumStr[i]; if(iv<0||iv>=N) { err=13; goto end;} // ������ ������� 'i' ���������� ����������� 'NumStr' - 'iv'.
jM=NumRow[i]; if(jM<0||jM>=N) { err=14; goto end;} // ������ ������� 'i' ���������� ����������� 'NumRow' - 'jM'.
jSh=iv*N; s=0.;
if(i<N-1) for(j=i+1;j<N;j++) {
jv=NumRow[j]; if(jv<0||jv>=N) { err=15; goto end;} // ������ ������� 'j' ���������� ����������� 'NumRow' - 'jv'.
s+=pM[jv+jSh]*X[jv];}
X[jM]=G[iv]-s;}

end: SAFE_DELETE_ARR(G); return err;
}


//-----------------------------------------------------------------------------------------------------------
// ���������� �������� ������� ������� ���������� ������ � ������� �������� ��������.

BYTE InvMatr(const clMatrix &M,clMatrix &Inv)
{
BYTE TriDec(int N,const clMatrix &A,clMatrix &T,int *NumStr,int *NumRow); // ����������� ���������� ������� � ������ ���������� ������ � ������� �������� �������� �� ������� � ��������.
BYTE TriSolv(int N,const clMatrix &T,int *NumStr,int *NumRow,complex *F,complex *X); // ������� ������� �������� ��������� �� ������������ ���������� ������� � ������ ���������� ������ � ������� �������� �������� �� ������� � ��������.

BYTE be,err; int N,j,*NumStr,*NumRow; clMatrix T; double Time1,Time2;

if(M.IsOK()!=0) return 1; if(M.Nx!=M.Ny) return 2; N=M.Nx;
if(Inv.IsOK()!=0) return 3; if(Inv.Nx!=N) return 4; if(Inv.Ny!=N) return 5;
NumStr=NumRow=NULL; err=0;

if(T.Alloc(N,N)!=0) { err=6; goto end;}
NumStr=new int[N]; if(NumStr==NULL) { err=9; goto end;} // ������ ������� ����������� ����� ��� ������ ����������� ��������.
NumRow=new int[N]; if(NumRow==NULL) { err=10; goto end;} // ������ ������� ����������� �������� ��� ������ ����������� ��������.

Time1=omp_get_wtime();
be=TriDec(N,M,T,NumStr,NumRow); // ����������� ����������.
if(be!=0) { if(be==messSingMatr) err=messSingMatr; else err=11; goto end;}
Time1=omp_get_wtime()-Time1;

//printf("Time 1 = %lf\n",Time1);

Time2=omp_get_wtime();
// ���������� � ������� ������ ����� �� ������� ������� ��������� ������� �������� ������� �������� �������.
#pragma omp parallel default(shared)
{
complex *X; complex *F; F=X=NULL; int i;
X=new complex[N]; if(X==NULL) err=8; // ������ ����������� ��� ������� ������� ������� ������.
F=new complex[N]; if(F==NULL) { err=7;} // ������ ������ ����� ��� ������� ������� ������� ������.
for(i=0;i<N;i++) F[i]=Cmplx_0;
#pragma omp for schedule(static) private(j)
for(j=0;j<N;j++) { F[j]=Cmplx_1;
if(TriSolv(N,T,NumStr,NumRow,F,X)!=0) { err=12; continue;} // ������� �������.
for(i=0;i<N;i++) Inv.Matr[i*N+j]=X[i];
F[j]=Cmplx_0;
}
SAFE_DELETE_ARR(X); SAFE_DELETE_ARR(F);
}
Time2=omp_get_wtime()-Time2;
//printf("Time2 = %lf\n",Time2); getch();
end: SAFE_DELETE_ARR(NumStr); SAFE_DELETE_ARR(NumRow); return err;
}

//-----------------------------------------------------------------------------------------------------------
// ���������� ������������ ������� ����������� �������� ������� �������� � �������������� ������� �����������.

BYTE SortEigVal(int N,clRowMatr Arr,int *Num)
{
int i,j,num; complex ca,cb;

if(N<=0) return 1; if(N<2) return 0; if(Arr.IsOK()!=0) return 2; if(Arr.N!=N) return 3; if(Num==NULL) return 4;
for(j=0;j<N-1;j++) { 
for(i=0;i<N-j-1;i++) { ca=Arr.Vect[Num[i]]; cb=Arr.Vect[Num[i+1]];
if(real(ca)>real(cb)) { num=Num[i]; Num[i]=Num[i+1]; Num[i+1]=num;}}}
return 0;
}

//-----------------------------------------------------------------------------------------------------------
// ������� ������. � 'pTMatr' - ��������� �������.

BYTE Conv(int nTMatr,complex *pTMatr,int N,clRowMatr &Row,clRowMatr &Res)
{
int i,j,iv; complex s;

if(nTMatr<=0) return 1; if(pTMatr==NULL) return 2; if(N<=0) return 3; if(N*2-1!=nTMatr) return 4;
if(Row.IsOK()!=0) return 5; if(Row.N!=N) return 6;
if(Res.IsOK()!=0) return 7; if(Res.N!=N) return 8;
for(i=0;i<N;i++) {
s=Cmplx_0; for(j=0;j<N;j++) { if(i>=j) iv=i-j; else iv=j-i+N-1; if(iv<0||iv>=nTMatr) return 9;
s+=pTMatr[iv]*Row.Vect[j];}
Res.Vect[i]=s;}
return 0;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#define flDiag 255 // ��� ���������, � ���, ��� ���������� � ����� ����������� �� ��������� - ������� ��� ������������.

//-----------------------------------------------------------------------------------------------------------
// ���������� ����������� �������� � �������� ������� 'A'.
// ����� ������� QR �������� � ����������� ������� � ���� ����������� ����������� �����������.

// 'EigLam' - ������ ����������� �������� �������.
// 'P' - ������� ����������� �������� (�������������), ������� ������������� �� ��������.
// ������� 'A' ���������� � ������������� ���� A=P.D.P^-1.

BYTE FindEigenValVect(int N,const clMatrix &A,clRowMatr &EigLam,clMatrix &P,clMatrix &PInv,BYTE flInv)
{
BYTE Hessen(clMatrix &M,clMatrix &Q); // ���������� ������� � ���� ����������� ��������� ���������������, ��������� �� ��������� �����������.
BYTE TriMatrQRShift(clMatrix &M,clMatrix &Q); // ���������� ������� � ������������������ ���� ������� QR ���������� �� �������.
BYTE EigVectTriMatrRight(clMatrix &M,clMatrix &EV); // ���������� ����������� �������� ����������������� ������� ����� ����. ������� ����������� �� �������.
BYTE EigVectTriMatrLeft(clMatrix &T,clMatrix &EV); // ���������� ����� ����������� �������� (�����) ����������������� ������� ����� ����. ������� ����������� �� �������.
BYTE Transp(const clMatrix &M,clMatrix &MT); // ���������������� �������.

BYTE be; int i,j,jv,jSh; complex s; clMatrix M,Q,QH,PT,PR;

// ��������.
if(N<2) return 1; if(A.IsOK()!=0) return 2; if(A.Nx!=N) return 3; if(A.Ny!=N) return 4;
if(EigLam.IsOK()!=0) return 5; if(EigLam.N!=N) return 6;
if(P.IsOK()!=0) return 7; if(P.Nx!=N||P.Ny!=N) return 8;
if(flInv!=0) { if(PInv.IsOK()!=0) return 9; if(PInv.Nx!=N||PInv.Ny!=N) return 10;}

// ��������� ������ ��� ��������������� ������.
if(Q.Alloc(N,N)!=0) return 11; if(QH.Alloc(N,N)!=0) return 12;
if(PT.Alloc(N,N)!=0) return 13; if(PR.Alloc(N,N)!=0) return 14;

Q=Cmplx_1; // ������� � 'Q' ��������� �������, ����� �������� ��������� ������� � ����� ����.
M=A; // �������� ������� 'A' � ������� 'M' - ��� ����� ������� ������� ���� �����������, ����� ����������������� �������.
if(N>=2) {
be=Hessen(M,Q); // ���������� ������� � ���� ����������� ��������� ���������������, ��������� �� ��������� �����������.
if(be==flDiag) goto Tri; if(be!=0) { printf("Hessen err : %d\n",be); return 15;}}
be=TriMatrQRShift(M,Q); if(be!=0) { printf("TriMatrQRShift err : %d\n",be); return 16;} // ���������� ������� � ������������������ ���� ������� QR ���������� �� �������.
Tri: for(i=0;i<N;i++) EigLam.Vect[i]=M.Matr[i*(N+1)]; // ����������� ��������.
// ���������� ����������� �������� ����������������� ������� 'M'.
be=EigVectTriMatrRight(M,PR); if(be!=0) { printf("EigVectTriMatrRight err : %d\n",be); return 17;}
if(Transp(PR,PT)!=0) return 18; // ���������� ����������� ������� ������� 'M', ������������� �� �������, ������������ � �������, ��� ��� ����������� �� ��������.
if(Hermit(Q,QH)!=0) return 19; // ��� ��������� ����������� �������� ������� 'A' ����� �������� ���������� ������� ������� 'Q'.
P=QH*PT; // ����������� ������� ������� 'A' �������� ������� ������� 'QH' �� ������� ����������� �������� ����������������� ������� 'M'.

if(flInv==0) return 0;
// ���������� ����� ����������� �������� (�����) ����������������� 'M'.
be=EigVectTriMatrLeft(M,PT); if(be!=0) { printf("EigVectTriMatrLeft err : %d\n",be); return 20;}
// ��������� ����� ����������� �������, ����� ��� ������ 1 ��� ��������� �� ������ ����������� �������.
for(i=0;i<N;i++) { jSh=i*N;
s=Cmplx_0; for(j=0;j<N;j++) { jv=jSh+j; s+=PT.Matr[jv]*PR.Matr[jv];} if(abs(s)==0.) return 21;
s=Inv(s); for(j=0;j<N;j++) PT.Matr[jSh+j]*=s;}
PInv=PT*Q; // ����� ����������� ������� (������) ������� 'A' �������� ������� ������� ����� ����������� �������� (�����) ����������������� ������� 'M' �� ������� 'Q'.

return 0;
}

//-----------------------------------------------------------------------------------------------------------
// ���������� ������� � ���� ����������� ��������� ���������������, ��������� �� ��������� �����������.

BYTE Hessen(clMatrix &M,clMatrix &Q)
{
BYTE Haus(clMatrix &M,clMatrix &Q,int numJ); // ��������� ����������� ��� �������, ��������� ������ �� 'P'='I-uuH'.

BYTE fl,be; int j,N;

if(M.IsOK()!=0) return 1; if(M.Nx!=M.Ny) return 2; N=M.Nx;
if(Q.IsOK()!=0) return 3; if(Q.Nx!=N) return 4; if(Q.Ny!=N) return 5;
fl=flDiag; for(j=0;j<N-2;j++) { be=Haus(M,Q,j); if(be==flDiag) continue; fl=0;
if(be!=0) { printf("Haus err,j: %d %d\n",be,j); return 6;}} return fl;
}

//-----------------------------------------------------------------------------------------------------------
// ��������� ����������� ��� �������, ��������� ������ �� 'P'='I-uuH'.

BYTE Haus(clMatrix &M,clMatrix &Q,int numJ)
{
int i,j,jSh,jShv,numI,N; double s,coe,a; complex c,cs,alp,*pM,*pQ;

if(M.IsOK()!=0) return 1; if(M.Nx!=M.Ny) return 2; N=M.Nx; pM=M.Matr;
if(Q.IsOK()!=0) return 3; if(Q.Nx!=N) return 4; if(Q.Ny!=N) return 5; pQ=Q.Matr;
numI=numJ+1; if(numI>=N-1) return 6; if(numJ>=N-2) return 7;

// ���������� ������� 'u', �������� ��� � ������� ������� � 'numI' � ������� 'numJ'.
s=0.; for(i=numI;i<N;i++) { jSh=i*N+numJ; c=pM[jSh]; s+=real(c*conj(c));} if(s<=0.) return flDiag;
s=sqrt(s); // ����� ������� 'x'.
// ������� ����������� ��� ������� 'e1' ��� �������, ������������ � ������� ���������� �������.
jSh=numI*N+numJ; c=pM[jSh]; a=abs(c); if(a<SmCnst32_d) alp=complex(s,0.); else alp=c*complex(s/a,0.);
pM[jSh]-=alp; // �������� �� ������� �������� ������� ��������� ��������� ������� 'e1'.

// ���������� �������� ����� ������� 'u'.
s=0.; for(i=numI;i<N;i++) { jSh=i*N+numJ; c=pM[jSh]; s+=real(c*conj(c));} if(s<=0.) goto fin;
coe=2./s;

// ��������� �����������. ��������� ��������� ������� 'P' (P=PH) �� ��� ����������� ������� ������� � ����� �� 'j'.
for(j=numJ+1;j<N;j++) {
cs=Cmplx_0; for(i=numI;i<N;i++) { jSh=i*N; cs+=conj(pM[jSh+numJ])*pM[jSh+j];} cs*=coe;
for(i=numI;i<N;i++) { jSh=i*N; pM[jSh+j]-=cs*pM[jSh+numJ];}
} // ����� ����� �� j.

// ��������� ��������� ������� 'P' (P=PH) �� ������� 'Q' ��� ��������� �������� ������� �������������� 'Q'.
for(j=0;j<N;j++) {
cs=Cmplx_0; for(i=numI;i<N;i++) { jSh=i*N; cs+=conj(pM[jSh+numJ])*pQ[jSh+j];} cs*=coe;
for(i=numI;i<N;i++) { jSh=i*N; pQ[jSh+j]-=cs*pM[jSh+numJ];}
} // ����� ����� �� j.

// M1=P�.M.P; ��������� �������������� ������� M. �������� ������� ������ �� ������� 'P'='I-uuT'.
for(i=0;i<N;i++) { jSh=i*N;
cs=Cmplx_0; for(j=numJ+1;j<N;j++) { jShv=j*N; cs+=pM[jShv+numJ]*pM[jSh+j];} cs*=coe;
for(j=numJ+1;j<N;j++) { jShv=j*N; pM[jSh+j]-=cs*conj(pM[jShv+numJ]);}
} // ����� ����� �� j.

// ��������� ��������� �������: � ������ ������� ��� ��������� ��������, � ��������� �������� ����.
fin: pM[numI*N+numJ]=alp; for(i=numI+1;i<N;i++) pM[i*N+numJ]=Cmplx_0;

return 0;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ��������� ���������.

#define nCntMax 2000 // ����� ��������, ����� �������� �������� ��������� ������������.
#define parAccZEAbs 1.e-18 // �������� ������� ��������������� ��������� (����������).
#define parAccZEAbsMax 1.e-12 // ������������ �������� ������� (����������) ��������������� ���������, ��� ������� �������� ��������� ���������� �� �������� ������������ ����� �����.
#define parAccZER 1.e-12 // �������� ������� ��������������� ��������� (�������������).
#define parAccZERMax 1.e-8 // ������������ �������� ������� (�������������) ��������������� ���������, ��� ������� �������� ��������� ���������� �� �������� ������������ ����� �����.
#define parCoeDiagQR 0.5 // ������������� �������� ������� ������������ ��������� ��� ������� � �������������� ��������������� ���������� � 'StepQRShift'.

//-----------------------------------------------------------------------------------------------------------
// ���������� ������� � ������������������ ���� ������� QR ���������� �� �������.

BYTE TriMatrQRShift(clMatrix &M,clMatrix &Q)
{
BYTE StepQRShift(clMatrix &M,clMatrix &Q,int numI,int numF); // �������� ������� QR ���������� �� ������� ��� ���������� �� 'numI'-�� �� 'numF'-�� ������ ������������.
BYTE CheckUnderDiagElem(clMatrix &M,int N,int numI,int numF,int *Flags,BYTE &fl); // �������� ������� ��������������� ���������, ��������� ������.
BYTE FindMinUnderDiagElem(clMatrix &M,int N,int numI,int numF,BYTE &fl,int &iMin,double &Min,BYTE ind); // ������� ����������� ��������������� ������� �� �������� ������������� 'ind'=0 � ���������� 'ind'!=0 �������.

BYTE fl,ind,indA[2],be,ber; int i,N,jSh,cnt,numI,numF,nI,iMin,*Flags; double Min[2]; complex *pM;

if(M.IsOK()!=0) return 1; if(M.Nx!=M.Ny) return 2; N=M.Nx; pM=M.Matr;
if(Q.IsOK()!=0) return 3; if(Q.Nx!=N) return 4; if(Q.Ny!=N) return 5;

Flags=NULL; ber=0;
Flags=new int[N]; if(Flags==NULL) { ber=6; goto end;}
for(i=0;i<N;i++) Flags[i]=0;

numI=0; numF=N-1; // ��������� �������� �������� ����������.
be=CheckUnderDiagElem(M,N,numI,numF,Flags,fl); if(be!=0) { ber=7; goto end;} // ��������� �������� ������� ��������������� ���������, ��������� ������.
if(fl!=0) goto ChoiceLim;

// ���� ����� ������� QR �������� ��� ���������� �� ��� ���, ���� �����-���� ��������������� ������� �� ������ ����������� ������ ����.
repQR: cnt=0; while(1) {
be=StepQRShift(M,Q,numI,numF); 
if(be!=0) { printf("StepQRShift numI,numF,cnt: %d %d %d err: %d\n",numI,numF,cnt,be); ber=8; goto end;} // ��� ������� QR ��������.
be=CheckUnderDiagElem(M,N,numI,numF,Flags,fl); if(be!=0) { ber=9; goto end;} // �������� ������� ��������������� ���������, ��������� ������.
if(fl!=0) break; // ��������� ��������������� ������� ��� �������� ����� ������� ���� - ������� �� ����� QR �������� ��� ������ ����������.
cnt++;

// ����� �������� ��������� �������� ������. ��������� ���������� ��������� ���� ���������������� �������� � ����� ������ ��������� ��������.
if(cnt>nCntMax) { indA[0]=indA[1]=0;
// �������� ��������������� ������� �� �������� ������������� �������, ������� �������� ������ ����. ����� ���������� ��������.
ind=0; if(FindMinUnderDiagElem(M,N,numI,numF,fl,iMin,Min[ind],ind)!=0) { ber=10; goto end;}
if(fl!=0&&iMin>=numI+1&&iMin<=numF) {
if(Min[ind]<parAccZERMax) { i=iMin; jSh=i*N; pM[jSh+i-1]=Cmplx_0; Flags[i]=1; break;}
indA[ind]=1;}

// �������� ��������������� ������� �� �������� ���������� �������, ������� �������� ������ ����. ����� ���������� ��������.
ind=1; if(FindMinUnderDiagElem(M,N,numI,numF,fl,iMin,Min[ind],ind)!=0) { ber=11; goto end;}
if(fl!=0&&iMin>=numI+1&&iMin<=numF) {
if(Min[ind]<parAccZEAbsMax) { i=iMin; jSh=i*N; pM[jSh+i-1]=Cmplx_0; Flags[i]=1; break;}
indA[ind]=1;}

// QR �������� �� ������� �� �������� ����� �����, � ��� ���� ������������ ������������� ��� ���������� �������� ���������������� �������� ��������� �������� �������.
printf("QR iterations do not converged: cnt: %d numI,numF: %d %d\n",cnt,numI,numF);
for(ind=0;ind<2;ind++) { if(indA[ind]==0) continue; printf("ind: %d Min: %16.9le\n",ind,Min[ind]);}
ber=12; goto end;}
} // ����� ����� while.

// ����� ��������� ���������� ��� ���������� ������ ���������, ���� ��� ��������������� �������� ����� ������.
ChoiceLim: nI=1;
numI=INT_MAX; for(i=nI;i<N;i++) if(Flags[i]==0) { numI=i-1; break;} if(numI==INT_MAX) goto end; // ��� ��������������� �������� ���� - ��������� ��������� ���������� � ������������� ����.
numF=INT_MAX; for(i=numI+1;i<N;i++) if(Flags[i]!=0) { numF=i-1; break;} if(numF==INT_MAX) numF=N-1; // ����� �� ����� �������.
goto repQR;

end: SAFE_DELETE_ARR(Flags); return ber;
}

//-----------------------------------------------------------------------------------------------------------
// �������� ������� ��������������� ���������, ��������� ������.

BYTE CheckUnderDiagElem(clMatrix &M,int N,int numI,int numF,int *Flags,BYTE &fl)
{
int i,jSh; double vd,ve; complex c,*pM;

if(M.IsOK()!=0) return 1; if(M.Nx!=N||M.Ny!=N) return 2; pM=M.Matr;
fl=0; for(i=numI+1;i<=numF;i++) { jSh=i*N; c=pM[jSh+i]; vd=abs(c); c=pM[jSh+i-1]; ve=abs(c);
if(vd<=0.) { if(ve>parAccZEAbs) continue;} else { if(ve>vd*parAccZER) continue;}
pM[jSh+i-1]=Cmplx_0; Flags[i]=1; fl=1;}
return 0;
}

//-----------------------------------------------------------------------------------------------------------
// ������� ����������� ��������������� ������� �� �������� ������������� 'ind'=0 � ���������� 'ind'!=0 �������.

BYTE FindMinUnderDiagElem(clMatrix &M,int N,int numI,int numF,BYTE &fl,int &iMin,double &Min,BYTE ind)
{
int i,jSh; double vd,ve,r; complex c,*pM;

if(M.IsOK()!=0) return 1; if(M.Nx!=N||M.Ny!=N) return 2; pM=M.Matr; Min=LrgCnst64_d; iMin=INT_MAX;
fl=0; for(i=numI+1;i<=numF;i++) { jSh=i*N; c=pM[jSh+i]; vd=abs(c);
if(ind==0) { if(vd<=0.) continue;} c=pM[jSh+i-1]; ve=abs(c);
if(ind==0) { r=ve/vd; if(r>Min) continue; Min=r;} else { if(ve>Min) continue; Min=ve;} iMin=i; fl=1;}
return 0;
}

//-----------------------------------------------------------------------------------------------------------
// �������� ������� QR ���������� �� ������� ��� ���������� �� 'numI'-�� �� 'numF'-�� ������ ������������.

BYTE StepQRShift(clMatrix &M,clMatrix &QS,int numI,int numF)
{
BYTE k; int nt,i,j,iF,jI,N,jSh,jShv,jShA,jShB; double v; complex *pM,*pQ,c,s,sig,Q[2][2],QH[2][2],Arr[2],a01,a10;

if(M.IsOK()!=0) return 1; if(M.Nx!=M.Ny) return 2; N=M.Nx; pM=M.Matr;
if(QS.IsOK()!=0) return 3; if(QS.Nx!=N) return 4; if(QS.Ny!=N) return 5; pQ=QS.Matr;
if(numI<0) return 6; if(numF>=N) return 7; if(numI>=numF) return 8;

// ������ �������������� �������������� � ������ ���� ������� � �������� ����������
// ��� ��������� ���������� � ������ ����� ����� ������������ ���������.
v=parCoeDiagQR; jSh=numI*N; jShv=jSh+N; jShA=jSh+numI; jShB=jShv+numI; a01=pM[jShA+1]; a10=pM[jShB];
if(abs(pM[jShA])<abs(a10)*v&&abs(pM[jShB+1])<abs(a01)*v) {
c=sqrt(a01); s=conj(sqrt(a10));
v=real(c*conj(c)+s*conj(s)); if(v<=0.) return 0; v=1./sqrt(v); c*=v; s*=v; // ����������.
Q[0][0]=c; Q[0][1]=-s; Q[1][0]=conj(s); Q[1][1]=conj(c); // ������� 'Q'.
QH[0][0]=conj(c); QH[0][1]=s; QH[1][0]=-conj(s); QH[1][1]=c; // �������� ���������� ������� 'QH'.

// �������� �������� ���������� ������� 'QH' �� 'M'. ���� �� k - �� ���� ������� ������� 'QH'. ������������ �� �������� ������� 'Q' � ������� ������� 'M'.
jI=numI;
for(j=jI;j<N;j++) { jShA=jSh+j; jShB=jShv+j;
for(k=0;k<2;k++) Arr[k]=pM[jShA]*QH[k][0]+pM[jShB]*QH[k][1]; pM[jShA]=Arr[0]; pM[jShB]=Arr[1];}

// �������� ������� 'M' �� 'Q'. ���� �� k - �� ���� ������� ������� 'M'. ������������ �� �������� ������� 'M' � ������� ������� 'Q'.
iF=numI+3; if(iF>numF) iF=numF;
for(i=0;i<=iF;i++) { jShA=i*N+numI; for(k=0;k<2;k++) Arr[k]=pM[jShA]*Q[0][k]+pM[jShA+1]*Q[1][k];
for(k=0;k<2;k++) pM[jShA+(int)k]=Arr[k];}

// �������� �������� ���������� ������� 'QH' �� 'QS'. ���� �� k - �� ���� ������� ������� 'QH'. ������������ �� �������� ������� 'Q' � ������� ������� 'QS'.
for(j=0;j<N;j++) { jShA=jSh+j; jShB=jShv+j;
for(k=0;k<2;k++) Arr[k]=pQ[jShA]*QH[k][0]+pQ[jShB]*QH[k][1]; pQ[jShA]=Arr[0]; pQ[jShB]=Arr[1];}
} // ����� �������������� ��������������.

sig=pM[numF*N+numF]; // ����� ������ - ��������� ������������ ������� ����������.

// ���� ����������� �������� �������.
for(nt=numI;nt<numF;nt++) {

// ����� ��������� ������� �������� Q.
// ��� ������ ������� Q1 �������� �������� ���, ����� ������ ������� Q (c,s^) ��� �������������� ������� ������� 'A-sig.I'.
if(nt==numI) { jSh=numI*N+numI; c=pM[jSh]-sig; s=conj(pM[jSh+N]);} // ������� ������� ������� 'A-sig.I'.
// ��� ����������� ������ Q �������� � �������� ���, ����� ������� ���� '�����������' �������
// ���� ������������ (��� ���������� '����'). ������ ����: -conj(s)*a10+c*a20=0.
else { jSh=nt*(N+1)-1; c=pM[jSh]; s=conj(pM[jSh+N]);} // c=a10 - M[nt;nt-1], s=a20^ - M[nt+1;nt-1]^.
v=real(c*conj(c)+s*conj(s)); if(v<=0.) return 0; v=1./sqrt(v); c*=v; s*=v; // ����������.
Q[0][0]=c; Q[0][1]=-s; Q[1][0]=conj(s); Q[1][1]=conj(c); // ������� 'Q'.
QH[0][0]=conj(c); QH[0][1]=s; QH[1][0]=-conj(s); QH[1][1]=c; // �������� ���������� ������� 'QH'.

// �������� �������� ���������� ������� 'QH' �� 'M'. ���� �� k - �� ���� ������� ������� 'QH'. ������������ �� �������� ������� 'Q' � ������� ������� 'M'.
jI=nt-1; if(jI<numI) jI=numI; jSh=nt*N; jShv=(nt+1)*N;
for(j=jI;j<N;j++) { jShA=jSh+j; jShB=jShv+j;
for(k=0;k<2;k++) Arr[k]=pM[jShA]*QH[k][0]+pM[jShB]*QH[k][1]; pM[jShA]=Arr[0]; pM[jShB]=Arr[1];}

// �������� ������� 'M' �� 'Q'. ���� �� k - �� ���� ������� ������� 'M'. ������������ �� �������� ������� 'M' � ������� ������� 'Q'.
iF=nt+3; if(iF>numF) iF=numF;
for(i=0;i<=iF;i++) { jSh=i*N+nt; for(k=0;k<2;k++) Arr[k]=pM[jSh]*Q[0][k]+pM[jSh+1]*Q[1][k];
for(k=0;k<2;k++) pM[jSh+(int)k]=Arr[k];}

// �������� �������� ���������� ������� 'QH' �� 'QS'. ���� �� k - �� ���� ������� ������� 'QH'. ������������ �� �������� ������� 'Q' � ������� ������� 'QS'.
jSh=nt*N; jShv=(nt+1)*N;
for(j=0;j<N;j++) { jShA=jSh+j; jShB=jShv+j;
for(k=0;k<2;k++) Arr[k]=pQ[jShA]*QH[k][0]+pQ[jShB]*QH[k][1]; pQ[jShA]=Arr[0]; pQ[jShB]=Arr[1];}

} // ����� ����� �� ��������� �������.

return 0;
}

/*
//-----------------------------------------------------------------------------------------------------------
// �������� ������� QR ���������� �� ������� ��� ���������� �� 'numI'-�� �� 'numF'-�� ������ ������������.

BYTE StepQRShift(clMatrix &M,clMatrix &QS,int numI,int numF)
{
BYTE k,flOut; int nt,i,j,iF,jI,N,jSh,jShv,jShA,jShB; double v;
complex *pM,*pQ,c,s,sig,Q[2][2],QH[2][2],Arr[2],a01,a10;

if(M.IsOK()!=0) return 1; if(M.Nx!=M.Ny) return 2; N=M.Nx; pM=M.Matr;
if(QS.IsOK()!=0) return 3; if(QS.Nx!=N) return 4; if(QS.Ny!=N) return 5; pQ=QS.Matr;
if(numI<0) return 6; if(numF>=N) return 7; if(numI>=numF) return 8;

// ������ �������������� �������������� � ������ ���� ������� � �������� ����������
// ��� ��������� ���������� � ������ ����� ����� ������������ ���������.
v=parCoeDiagQR; jSh=numI*N; jShv=jSh+N; jShA=jSh+numI; jShB=jShv+numI; a01=pM[jShA+1]; a10=pM[jShB];
if(abs(pM[jShA])<abs(a10)*v&&abs(pM[jShB+1])<abs(a01)*v) {
c=sqrt(a01); s=conj(sqrt(a10));
v=real(c*conj(c)+s*conj(s)); if(v<=0.) return 0; v=1./sqrt(v); c*=v; s*=v; // ����������.
Q[0][0]=c; Q[0][1]=-s; Q[1][0]=conj(s); Q[1][1]=conj(c); // ������� 'Q'.
QH[0][0]=conj(c); QH[0][1]=s; QH[1][0]=-conj(s); QH[1][1]=c; // �������� ���������� ������� 'QH'.

// �������� �������� ���������� ������� 'QH' �� 'M'. ���� �� k - �� ���� ������� ������� 'QH'. ������������ �� �������� ������� 'Q' � ������� ������� 'M'.
jI=numI;
for(j=jI;j<N;j++) { jShA=jSh+j; jShB=jShv+j;
for(k=0;k<2;k++) Arr[k]=pM[jShA]*QH[k][0]+pM[jShB]*QH[k][1]; pM[jShA]=Arr[0]; pM[jShB]=Arr[1];}

// �������� ������� 'M' �� 'Q'. ���� �� k - �� ���� ������� ������� 'M'. ������������ �� �������� ������� 'M' � ������� ������� 'Q'.
iF=numI+3; if(iF>numF) iF=numF;
for(i=0;i<=iF;i++) { jShA=i*N+numI; for(k=0;k<2;k++) Arr[k]=pM[jShA]*Q[0][k]+pM[jShA+1]*Q[1][k];
for(k=0;k<2;k++) pM[jShA+(int)k]=Arr[k];}

// �������� �������� ���������� ������� 'QH' �� 'QS'. ���� �� k - �� ���� ������� ������� 'QH'. ������������ �� �������� ������� 'Q' � ������� ������� 'QS'.
for(j=0;j<N;j++) { jShA=jSh+j; jShB=jShv+j;
for(k=0;k<2;k++) Arr[k]=pQ[jShA]*QH[k][0]+pQ[jShB]*QH[k][1]; pQ[jShA]=Arr[0]; pQ[jShB]=Arr[1];}
} // ����� �������������� ��������������.

sig=pM[numF*N+numF]; // ����� ������ - ��������� ������������ ������� ����������.

flOut=0;
// ������ ������������ �������.
#pragma omp parallel default(shared) private(nt,c,s,v,jI,jSh,jShv,iF) firstprivate(flOut)
{
complex Q_[2][2],QH_[2][2],Arr_[2];

// ���� ����������� �������� �������.
for(nt=numI;nt<numF;nt++) {
if(flOut!=0) continue;

// ����� ��������� ������� �������� Q.
// ��� ������ ������� Q1 �������� �������� ���, ����� ������ ������� Q (c,s^) ��� �������������� ������� ������� 'A-sig.I'.
if(nt==numI) { jSh=numI*N+numI; c=pM[jSh]-sig; s=conj(pM[jSh+N]);} // ������� ������� ������� 'A-sig.I'.
// ��� ����������� ������ Q �������� � �������� ���, ����� ������� ���� '�����������' �������
// ���� ������������ (��� ���������� '����'). ������ ����: -conj(s)*a10+c*a20=0.
else { jSh=nt*(N+1)-1; c=pM[jSh]; s=conj(pM[jSh+N]);} // c=a10 - M[nt;nt-1], s=a20^ - M[nt+1;nt-1]^.
v=real(c*conj(c)+s*conj(s)); if(v<=0.) { flOut=1; continue;} v=1./sqrt(v); c*=v; s*=v; // ����������.
Q_[0][0]=c; Q_[0][1]=-s; Q_[1][0]=conj(s); Q_[1][1]=conj(c); // ������� 'Q'.
QH_[0][0]=conj(c); QH_[0][1]=s; QH_[1][0]=-conj(s); QH_[1][1]=c; // �������� ���������� ������� 'QH'.

// �������� �������� ���������� ������� 'QH' �� 'M'. ���� �� k - �� ���� ������� ������� 'QH'. ������������ �� �������� ������� 'Q' � ������� ������� 'M'.
jI=nt-1; if(jI<numI) jI=numI; jSh=nt*N; jShv=(nt+1)*N;
#pragma omp barrier
#pragma omp for schedule(static) private(j,jShA,jShB,k)
for(j=jI;j<N;j++) { jShA=jSh+j; jShB=jShv+j;
for(k=0;k<2;k++) Arr_[k]=pM[jShA]*QH_[k][0]+pM[jShB]*QH_[k][1]; pM[jShA]=Arr_[0]; pM[jShB]=Arr_[1];}
#pragma omp barrier

// �������� ������� 'M' �� 'Q'. ���� �� k - �� ���� ������� ������� 'M'. ������������ �� �������� ������� 'M' � ������� ������� 'Q'.
iF=nt+3; if(iF>numF) iF=numF;
#pragma omp for schedule(static) private(i,jSh,k)
for(i=0;i<=iF;i++) { jSh=i*N+nt; for(k=0;k<2;k++) Arr_[k]=pM[jSh]*Q_[0][k]+pM[jSh+1]*Q_[1][k];
for(k=0;k<2;k++) pM[jSh+(int)k]=Arr_[k];}
#pragma omp barrier

// �������� �������� ���������� ������� 'QH' �� 'QS'. ���� �� k - �� ���� ������� ������� 'QH'. ������������ �� �������� ������� 'Q' � ������� ������� 'QS'.
jSh=nt*N; jShv=(nt+1)*N;
#pragma omp for schedule(static) private(j,jShA,jShB,k)
for(j=0;j<N;j++) { jShA=jSh+j; jShB=jShv+j;
for(k=0;k<2;k++) Arr_[k]=pQ[jShA]*QH_[k][0]+pQ[jShB]*QH_[k][1]; pQ[jShA]=Arr_[0]; pQ[jShB]=Arr_[1];}
#pragma omp barrier

} // ����� ����� �� ��������� �������.
} // ����� ������������ �������.

return 0;
}
*/

//-----------------------------------------------------------------------------------------------------------
// ���������� ����������� �������� ����������������� ������� ����� ����. ������� ����������� �� �������.

BYTE EigVectTriMatrRight(clMatrix &T,clMatrix &EV)
{
BYTE SolvSystemTri(clMatrix &T,complex Lam,clRowMatr &F,clRowMatr &X,int numI,int numF); // ���������� ������� ������� � ����������������� �������� (T-Lam.I).X=F. ��������������� ���������� �� 'numI' �� 'NumF'.

int i,j,jSh,k,N,numF; double s,v; complex *pT,*pEV,Lam,c; clRowMatr F,X;

if(T.IsOK()!=0) return 1; if(T.Nx!=T.Ny) return 2; N=T.Nx; pT=T.Matr;
if(EV.IsOK()!=0) return 3; if(EV.Nx!=N) return 4; if(EV.Ny!=N) return 5; pEV=EV.Matr;
if(F.Alloc(N)!=0) return 6; if(X.Alloc(N)!=0) return 7;

pEV[0]=Cmplx_1; for(j=1;j<N;j++) pEV[j]=Cmplx_0; // ������ ����������� ������.

// ���� �� ����������� ������ � ����������� ��������.
for(k=1;k<N;k++) { jSh=k*N; Lam=pT[jSh+k]; numF=k-1;
for(i=0;i<=numF;i++) F.Vect[i]=pT[i*N+k]; // ���������� ������� 'F' ��� ������� ������� � ����������������� ��������.
if(SolvSystemTri(T,Lam,F,X,0,numF)!=0) return 8; // ���������� ��������� (�������) ����� ������������ �������.
for(j=0;j<=numF;j++) pEV[j+jSh]=-X.Vect[j]; pEV[k+jSh]=Cmplx_1; // ���������� ��������� (�������) ����� ������������ �������.
if(k+1<N) for(j=k+1;j<N;j++) pEV[j+jSh]=Cmplx_0; // ������ ����� ������������ ������� ��������� ������.
s=0.; for(j=0;j<=k;j++) { c=pEV[j+jSh]; s+=real(c*conj(c));} if(s<=0.) return 9; s=sqrt(s); // ����� ������������ �������.
v=1./s; for(j=0;j<=k;j++) pEV[j+jSh]*=v; // ���������� ������������ �������.
}
return 0;
}

//-----------------------------------------------------------------------------------------------------------
// ���������� ����� ����������� �������� (�����) ����������������� ������� ����� ����. ������� ����������� �� �������.

BYTE EigVectTriMatrLeft(clMatrix &T,clMatrix &EV)
{
BYTE SolvSystemTri(clMatrix &T,complex Lam,clRowMatr &F,clRowMatr &X,int numI,int numF); // ���������� ������� ������� � ����������������� �������� (T-Lam.I).X=F. ��������������� ���������� �� 'numI' �� 'NumF'.

int i,j,jSh,jShv,k,kv,N,numF; double s,v; complex *pT,*pEV,Lam,c; clMatrix TS; clRowMatr F,X;

if(T.IsOK()!=0) return 1; if(T.Nx!=T.Ny) return 2; N=T.Nx;
if(EV.IsOK()!=0) return 3; if(EV.Nx!=N) return 4; if(EV.Ny!=N) return 5; pEV=EV.Matr;
if(F.Alloc(N)!=0) return 6; if(X.Alloc(N)!=0) return 7;
if(TS.Alloc(N,N)!=0) return 8; pT=TS.Matr;

// ������� ����������� ������� 'TS', ������������ ������� 'T' ������������ ��������� 'i=N-1-j'.
for(i=0;i<N;i++) { jSh=i*N; for(j=i;j<N;j++) { jShv=(N-1-j)*N; TS.Matr[jShv+N-1-i]=T.Matr[jSh+j];}}
for(i=1;i<N;i++) { jSh=i*N; for(j=0;j<i;j++) TS.Matr[jSh+j]=Cmplx_0;}

jSh=(N-1)*N; pEV[jSh+N-1]=Cmplx_1; for(j=0;j<N-1;j++) pEV[jSh+j]=Cmplx_0; // ��������� ����������� ������.

// ���� �� ����������� ������ � ����������� ��������.
for(k=1;k<N;k++) { jSh=k*N; kv=N-1-k; jShv=kv*N; Lam=pT[jSh+k]; numF=k-1;
for(i=0;i<=numF;i++) F.Vect[i]=pT[i*N+k]; // ���������� ������� 'F' ��� ������� ������� � ����������������� ��������.
if(SolvSystemTri(TS,Lam,F,X,0,numF)!=0) return 8; // ���������� �������� ����� ������������ �������.
for(j=0;j<=numF;j++) pEV[N-1-j+jShv]=-X.Vect[j]; pEV[kv+jShv]=Cmplx_1; // ���������� �������� ����� ������������ �������.
if(k+1<N) for(j=k+1;j<N;j++) pEV[N-1-j+jShv]=Cmplx_0; // ��������� ����� ������������ ������� ��������� ������.
s=0.; for(j=0;j<=k;j++) { c=pEV[N-1-j+jShv]; s+=real(c*conj(c));} if(s<=0.) return 9; s=sqrt(s); // ����� ������������ �������.
v=1./s; for(j=0;j<=k;j++) pEV[N-1-j+jShv]*=v; // ���������� ������������ �������.
}
return 0;
}

//-----------------------------------------------------------------------------------------------------------
// ���������� ������� ������� � ����������������� �������� (T-Lam.I).X=F. ��������������� ���������� �� 'numI' �� 'NumF'.

BYTE SolvSystemTri(clMatrix &T,complex Lam,clRowMatr &F,clRowMatr &X,int numI,int numF)
{
int i,j,jSh,N; complex s,*pT,*pF,*pX;

if(T.IsOK()!=0) return 1; if(T.Nx!=T.Ny) return 2; N=T.Nx; pT=T.Matr;
if(F.IsOK()!=0) return 3; if(F.N!=N) return 4; pF=F.Vect;
if(X.IsOK()!=0) return 5; if(X.N!=N) return 6; pX=X.Vect;
if(numI<0) return 7; if(numF>=N) return 8; if(numI>numF) return 9;
for(i=numF;i>=numI;i--) { jSh=i*N;
s=Cmplx_0; if(i<numF) for(j=i+1;j<=numF;j++) s+=pT[j+jSh]*pX[j]; pX[i]=(pF[i]-s)/(pT[i+jSh]-Lam);}
return 0;
}

//-----------------------------------------------------------------------------------------------------------
// ������ ������ �� ���������� �����.

BYTE ReadString(FILE *file,char *pStr,int lnStr)
{
if(file==NULL) return 1; if(pStr==NULL) return 3; if(lnStr<=0) return 4;
#if flScanf_S!=0
if(fscanf_s(file,"%s\n",pStr,lnStr)==EOF) return 4;
#else
if(fscanf(file,"%s\n",pStr)==EOF) return 5;
#endif
return 0;
}

//-----------------------------------------------------------------------------------------------------------
// ��������� ����� ����� � �����������.

char *GetFName(char *FName,char *ext)
{
size_t ln; char *cp,*fn;

if(FName==NULL) return NULL; ln=strlen(FName); if(ln==0) return NULL;
cp=strchr(FName,'.'); if(cp!=NULL) return NULL;
if(ext==NULL) return NULL; cp=strchr(ext,'.'); if(cp==NULL) return NULL; if(cp!=ext) return NULL;
ln+=strlen(ext); fn=new char[ln+1]; if(fn==NULL) return NULL; fn[ln]='\0';
strcpy(fn,FName); strcat(fn,ext); return fn;
}