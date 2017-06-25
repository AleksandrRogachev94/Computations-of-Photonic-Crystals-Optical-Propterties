/*...........................................................................................................

// ����� ������� �������.
class clRowMatr { clRowMatr,clRowMatr,clRowMatr,=,Zero,Free,Alloc,IsOK,-,+=R,-=R,+=c,-=c,*=c,/=c,*=d,/=d,
SetZero,SetRand,GetAbsMax,GetNorm,Norm},

// ��������� � ������� ��� ������ �������.
+,-,*c,/c,c*,R+c,c+R,R-c,R*R,!,Conj,

// ����� ������� �������� 'clRowMatr'.
class clRows { clRows,clRows,~clRows,=,Zero,Free,Alloc_,Alloc,IsOK,IsOK_All,Get,Write},

// ����� ��������� �������.
class clMatrix { clMatrix,clMatrix,~clMatrix,=,Zero,Free,Alloc,IsOK,=c,<<c,=R,-,+=c,-=c,*=c,/=c,*=d,/=d,
+=R,-=R},

// ��������� � ������� ��� ������ ��������� �������.
+,-,*,*c,~,+c,c+,-c,c-,*c,c*,Transp,Hermit,

// ��������� ��� ������� ��������� ������� � �������.
M*R,M+R,R+M,M-R,R-M,M%R,R%M

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

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// ����� ������� �������.

//-----------------------------------------------------------------------------------------------------------
// �����������.

clRowMatr::clRowMatr(void)
{
Zero();
}

//-----------------------------------------------------------------------------------------------------------
// ���������� �����������.

clRowMatr::clRowMatr(const clRowMatr &RM_)
{
int i;

Zero(); if(RM_.IsOK()!=0) return; if(Alloc(RM_.N)!=0) return; for(i=0;i<N;i++) Vect[i]=RM_.Vect[i];
}

//-----------------------------------------------------------------------------------------------------------
// ����������.

clRowMatr::~clRowMatr(void)
{
Free();
}

//-----------------------------------------------------------------------------------------------------------
// �������� �����������. 

BYTE clRowMatr::operator =(const clRowMatr &RM_)
{
BYTE fl; int i;

if(RM_.IsOK()!=0) return 1; fl=0; if(IsOK()!=0) fl=1; else if(N!=RM_.N) fl=1;
if(fl!=0) { if(Alloc(RM_.N)!=0) return 2;} for(i=0;i<N;i++) Vect[i]=RM_.Vect[i]; return 0;
}

//-----------------------------------------------------------------------------------------------------------
// ��������� ����������.

void clRowMatr::Zero(void)
{
Vect=NULL; N=0;
}

//-----------------------------------------------------------------------------------------------------------
// ������������ ������.

void clRowMatr::Free(void)
{
SAFE_DELETE_ARR(Vect); N=0;
}

//-----------------------------------------------------------------------------------------------------------
// ��������� ������ ��� �������.

BYTE clRowMatr::Alloc(int N_)
{
Free(); if(N_<=0) return 1; N=N_; Vect=new complex[N]; if(Vect==NULL) return 2; return 0;
} 

//-----------------------------------------------------------------------------------------------------------
// ��������.

BYTE clRowMatr::IsOK(void) const
{
if(N<=0) return 1; if(Vect==NULL) return 2; return 0;
}

//-----------------------------------------------------------------------------------------------------------
// �������� ���������� ���������� ������.

BYTE clRowMatr::operator =(complex c)
{
int i;

if(IsOK()!=0) return 1; for(i=0;i<N;i++) Vect[i]=c; return 0;
}

//-----------------------------------------------------------------------------------------------------------
// ���� �����.

clRowMatr clRowMatr::operator -()
{
int i; clRowMatr RM;

if(IsOK()!=0) goto end; if(RM.Alloc(N)!=0) goto end; for(i=0;i<N;i++) RM.Vect[i]=-Vect[i];
end: return RM;
}

//-----------------------------------------------------------------------------------------------------------
// ����������� �������.

clRowMatr &clRowMatr::operator +=(clRowMatr const &RM)
{
int i;

if(IsOK()!=0) goto end; if(RM.IsOK()!=0) goto end; if(N!=RM.N) goto end;
for(i=0;i<N;i++) Vect[i]+=RM.Vect[i];
end: return *this;
}

//-----------------------------------------------------------------------------------------------------------
// ��������� �������.

clRowMatr &clRowMatr::operator -=(clRowMatr const &RM)
{
int i;

if(IsOK()!=0) goto end; if(RM.IsOK()!=0) goto end; if(N!=RM.N) goto end;
for(i=0;i<N;i++) Vect[i]-=RM.Vect[i];
end: return *this;
}

//-----------------------------------------------------------------------------------------------------------
// ����������� ��������� ������� � ���������� 'c'.

clRowMatr &clRowMatr::operator +=(complex c)
{
int i;

if(IsOK()!=0) return *this; for(i=0;i<N;i++) Vect[i]+=c; return *this;
}

//-----------------------------------------------------------------------------------------------------------
// ��������� ��������� ������� � ���������� 'c'.

clRowMatr &clRowMatr::operator -=(complex c)
{
int i;

if(IsOK()!=0) return *this; for(i=0;i<N;i++) Vect[i]-=c; return *this;
}

//-----------------------------------------------------------------------------------------------------------
// ��������� �� ����������� �����.

clRowMatr &clRowMatr::operator *=(complex c)
{
int i;

if(IsOK()!=0) return *this; for(i=0;i<N;i++) Vect[i]*=c; return *this;
}

//-----------------------------------------------------------------------------------------------------------
// ������� �� ����������� �����.

clRowMatr &clRowMatr::operator /=(complex c)
{
int i; complex cInv;

if(IsOK()!=0) return *this; cInv=Inv(c); for(i=0;i<N;i++) Vect[i]*=cInv; return *this;
}

//-----------------------------------------------------------------------------------------------------------
// ��������� �� �������������� �����.

clRowMatr &clRowMatr::operator *=(double d)
{
int i;

if(IsOK()!=0) return *this; for(i=0;i<N;i++) Vect[i]*=d; return *this;
}

//-----------------------------------------------------------------------------------------------------------
// ������� �� �������������� �����.

clRowMatr &clRowMatr::operator /=(double d)
{
int i; double v;

if(IsOK()!=0) return *this; if(d==0.) return *this; v=1./d; for(i=0;i<N;i++) Vect[i]*=v; return *this;
}

//-----------------------------------------------------------------------------------------------------------
// ��������� �������.

BYTE clRowMatr::SetZero(void)
{
int i;

if(IsOK()!=0) return 1; for(i=0;i<N;i++) Vect[i]=Cmplx_0; return 0;
}

//-----------------------------------------------------------------------------------------------------------
// ���������� ���������� ������� � �������� Re,Im �� -Lim �� +Lim.

BYTE clRowMatr::SetRand(double Lim)
{
int i; double coe,Re,Im;

if(IsOK()!=0) return 1; coe=2./(double)RAND_MAX;
for(i=0;i<N;i++) { Re=((double)rand()*coe-1.)*Lim; Im=((double)rand()*coe-1.)*Lim; Vect[i]=complex(Re,Im);}
return 0;
}

//-----------------------------------------------------------------------------------------------------------
// ��������� ������������� ��������.

double clRowMatr::GetAbsMax(void)
{
int i; double v,vM;

if(IsOK()!=0) return 0.; vM=0.; for(i=0;i<N;i++) { v=abs(Vect[i]); if(v>vM) vM=v;} return vM;
}

//-----------------------------------------------------------------------------------------------------------
// ��������� ����� �������.

double clRowMatr::GetNorm(void)
{
int i; double s; complex c;

if(IsOK()!=0) return 0.; s=0.; for(i=0;i<N;i++) { c=Vect[i]; s+=real(c*conj(c));} if(s<=0.) return 0.;
return sqrt(s);
}

//-----------------------------------------------------------------------------------------------------------
// ���������� ������� �������.

BYTE clRowMatr::Norm(void)
{
int i; double s; complex c;

if(IsOK()!=0) return 1; s=0.; for(i=0;i<N;i++) { c=Vect[i]; s+=real(c*conj(c));} if(s<=0.) return 2;
s=1./sqrt(s); for(i=0;i<N;i++) Vect[i]*=s; return 0;
}

//-----------------------------------------------------------------------------------------------------------
// ������ ������� � ����.

BYTE clRowMatr::Write(FILE *fp)
{
int i;

if(fp==NULL) return 1; for(i=0;i<N;i++) fprintf(fp,"%lf %lf   ",Vect[i].re,Vect[i].im); return 0;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// ��������� � ������� ��� ������ �������.

//-----------------------------------------------------------------------------------------------------------
// ��������.

clRowMatr operator +(const clRowMatr &RM1,const clRowMatr &RM2)
{
int i; clRowMatr RM;

if(RM1.IsOK()!=0||RM2.IsOK()!=0) goto end; if(RM1.N!=RM2.N) goto end; if(RM.Alloc(RM1.N)!=0) goto end;
for(i=0;i<RM.N;i++) RM.Vect[i]=RM1.Vect[i]+RM2.Vect[i]; end: return RM;
}

//-----------------------------------------------------------------------------------------------------------
// ���������.

clRowMatr operator -(const clRowMatr &RM1,const clRowMatr &RM2)
{
int i; clRowMatr RM;

if(RM1.IsOK()!=0||RM2.IsOK()!=0) goto end; if(RM1.N!=RM2.N) goto end; if(RM.Alloc(RM1.N)!=0) goto end;
for(i=0;i<RM.N;i++) RM.Vect[i]=RM1.Vect[i]-RM2.Vect[i]; end: return RM;
}

//-----------------------------------------------------------------------------------------------------------
// ��������� �� �����.

clRowMatr operator *(const clRowMatr &RM_,complex c)
{
int i; clRowMatr RM;

if(RM_.IsOK()!=0) goto end; if(RM.Alloc(RM_.N)!=0) goto end;
for(i=0;i<RM.N;i++) RM.Vect[i]=RM_.Vect[i]*c; end: return RM;
}

//-----------------------------------------------------------------------------------------------------------
// ������� �� �����.

clRowMatr operator /(const clRowMatr &RM_,complex c)
{
int i; complex cInv; clRowMatr RM;

if(RM_.IsOK()!=0) goto end; if(RM.Alloc(RM_.N)!=0) goto end; cInv=Inv(c);
for(i=0;i<RM.N;i++) RM.Vect[i]=RM_.Vect[i]*cInv; end: return RM;
}

//-----------------------------------------------------------------------------------------------------------
// ��������� �� �����.

clRowMatr operator *(complex c,const clRowMatr &RM_)
{
int i; clRowMatr RM;

if(RM_.IsOK()!=0) goto end; if(RM.Alloc(RM_.N)!=0) goto end;
for(i=0;i<RM.N;i++) RM.Vect[i]=RM_.Vect[i]*c; end: return RM;
}

//-----------------------------------------------------------------------------------------------------------
// ����������� ��������� ������� � ���������� 'c'.

clRowMatr operator +(const clRowMatr &RM_,complex c)
{
int i; clRowMatr RM;

if(RM_.IsOK()!=0) goto end; if(RM.Alloc(RM_.N)!=0) goto end;
for(i=0;i<RM.N;i++) RM.Vect[i]=RM_.Vect[i]+c; end: return RM;
}

//-----------------------------------------------------------------------------------------------------------
// ����������� ��������� ������� � ���������� 'c'.

clRowMatr operator +(complex c,const clRowMatr &RM_)
{
int i; clRowMatr RM;

if(RM_.IsOK()!=0) goto end; if(RM.Alloc(RM_.N)!=0) goto end;
for(i=0;i<RM.N;i++) RM.Vect[i]=RM_.Vect[i]+c; end: return RM;
}

//-----------------------------------------------------------------------------------------------------------
// ��������� ��������� ������� � ���������� 'c'.

clRowMatr operator -(const clRowMatr &RM_,complex c)
{
int i; clRowMatr RM;

if(RM_.IsOK()!=0) goto end; if(RM.Alloc(RM_.N)!=0) goto end;
for(i=0;i<RM.N;i++) RM.Vect[i]=RM_.Vect[i]-c; end: return RM;
}

//-----------------------------------------------------------------------------------------------------------
// ��������� ������������ ������.

clRowMatr operator *(const clRowMatr &RM1,const clRowMatr &RM2)
{
int i,N; clRowMatr RM;

if(RM1.IsOK()!=0) goto end; if(RM2.IsOK()!=0) goto end; N=RM1.N; if(RM2.N!=N) goto end;
if(RM.Alloc(N)!=0) goto end; for(i=0;i<N;i++) RM.Vect[i]=RM1.Vect[i]*RM2.Vect[i];
end: return RM;
}

//-----------------------------------------------------------------------------------------------------------
// �������� ������������ �������.

clRowMatr operator !(const clRowMatr &RM_)
{
int i; clRowMatr RM;

if(RM_.IsOK()!=0) goto end; if(RM.Alloc(RM_.N)!=0) goto end;
for(i=0;i<RM.N;i++) RM.Vect[i]=Inv(RM_.Vect[i]); end: return RM;
}

//-----------------------------------------------------------------------------------------------------------
// ����������� ����������.

clRowMatr Conj(clRowMatr const &RM_)
{
int i; clRowMatr RM;

if(RM_.IsOK()!=0) goto end; if(RM.Alloc(RM_.N)!=0) goto end;
for(i=0;i<RM.N;i++) RM.Vect[i]=conj(RM_.Vect[i]); end: return RM;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// ����� ������� �������� 'clRowMatr'.

//-----------------------------------------------------------------------------------------------------------
// �����������.

clRows::clRows(void)
{
Zero();
}

//-----------------------------------------------------------------------------------------------------------
// ���������� �����������.

clRows::clRows(const clRows &R_)
{
int i;

Zero(); if(R_.IsOK()!=0) return; if(Alloc_(R_.nRows)!=0) return;
for(i=0;i<nRows;i++) Rows[i]=R_.Rows[i];
}

//-----------------------------------------------------------------------------------------------------------
// ����������.

clRows::~clRows(void)
{
Free();
}

//-----------------------------------------------------------------------------------------------------------
// �������� �����������. 

BYTE clRows::operator =(const clRows &R_)
{
BYTE fl; int i;

if(R_.IsOK()!=0) return 1; fl=0; if(IsOK()!=0) fl=1; else if(nRows!=R_.nRows) fl=1;
if(fl!=0) { if(Alloc_(R_.nRows)!=0) return 2;} N=R_.N;
for(i=0;i<nRows;i++) Rows[i]=R_.Rows[i]; return 0;
}

//-----------------------------------------------------------------------------------------------------------
// ��������� ����������.

void clRows::Zero(void)
{
Rows=NULL; nRows=N=0;
}

//-----------------------------------------------------------------------------------------------------------
// ������������ ������.

void clRows::Free(void)
{
int i;

if(IsOK()!=0) for(i=0;i<nRows;i++) Rows[i].Free(); SAFE_DELETE_ARR(Rows); nRows=N=0;
}

//-----------------------------------------------------------------------------------------------------------
// ��������� ������ ��� �������� 'clRowMatr'.

BYTE clRows::Alloc_(int nRows_)
{
Free(); if(nRows_<=0) return 1; Rows=new clRowMatr[nRows_]; if(Rows==NULL) return 2; nRows=nRows_; return 0;
}

//-----------------------------------------------------------------------------------------------------------
// ��������� ������ ��� ���� ��������.

BYTE clRows::Alloc(int nRows_,int N_)
{
BYTE ber; int i;

Free(); if(nRows_<=0) return 1; if(N_<=0) return 2; if(Alloc_(nRows_)!=0) return 3;
N=N_; ber=0;
for(i=0;i<nRows;i++) if(Rows[i].Alloc(N)!=0) { ber=4; goto end;}
end: if(ber!=0) Free(); return ber;
} 

//-----------------------------------------------------------------------------------------------------------
// ��������.

BYTE clRows::IsOK(void) const
{
if(nRows<=0) return 1; if(Rows==NULL) return 2; return 0;
}

//-----------------------------------------------------------------------------------------------------------
// ��������, ������� �������.

BYTE clRows::IsOK_All(void) const
{
int i;

if(IsOK()!=0) return 1; for(i=0;i<nRows;i++) if(Rows[i].IsOK()!=0) return 2; return 0;
}

//-----------------------------------------------------------------------------------------------------------
// ��������� ��������� �� ������� � ������� 'num'.

clRowMatr *clRows::Get(int num)
{
if(IsOK()!=0) return NULL; if(num<0||num>=nRows) return NULL; return Rows+num;
}

//-----------------------------------------------------------------------------------------------------------
// ������ ������� � ����.

BYTE clRows::Write(FILE *fp)
{
int i;

if(fp==NULL) return 1; if(IsOK_All()!=0) return 2;
fprintf(fp,"%d %d\n",N,nRows);
for(i=0;i<nRows;i++) { if(Rows[i].Write(fp)!=0) return 3; fprintf(fp,"\n");}
return 0;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// ����� ��������� �������.

//-----------------------------------------------------------------------------------------------------------
// �����������.

clMatrix::clMatrix(void)
{
Zero();
}

//-----------------------------------------------------------------------------------------------------------
// ���������� �����������.

clMatrix::clMatrix(const clMatrix &M_)
{
int i,j,jSh;

Zero(); if(M_.IsOK()!=0) return; if(Alloc(M_.Nx,M_.Ny)!=0) return;
for(i=0;i<Nx;i++) { jSh=i*Ny; for(j=0;j<Ny;j++) Matr[j+jSh]=M_.Matr[j+jSh];}
}

//-----------------------------------------------------------------------------------------------------------
// ����������.

clMatrix::~clMatrix(void)
{
Free();
}

//-----------------------------------------------------------------------------------------------------------
// �������� �����������. 

BYTE clMatrix::operator =(const clMatrix &M_)
{
BYTE fl; int i,j,jSh;

if(M_.IsOK()!=0) return 1; fl=0; if(IsOK()!=0) fl=1; else if(Nx!=M_.Nx||Ny!=M_.Ny) fl=1;
if(fl!=0) { if(Alloc(M_.Nx,M_.Ny)!=0) return 2;}
for(i=0;i<Nx;i++) { jSh=i*Ny; for(j=0;j<Ny;j++) Matr[j+jSh]=M_.Matr[j+jSh];}
return 0;
}

//-----------------------------------------------------------------------------------------------------------
// ��������� ����������.

void clMatrix::Zero(void)
{
Matr=NULL; Nx=Ny=0;
}

//-----------------------------------------------------------------------------------------------------------
// ������������ ������.

void clMatrix::Free(void)
{
SAFE_DELETE_ARR(Matr); Nx=Ny=0;
}

//-----------------------------------------------------------------------------------------------------------
// ��������� ������ ��� �������.

BYTE clMatrix::Alloc(int Nx_,int Ny_)
{
Free(); if(Nx_<=0||Ny_<=0) return 1; if(INT_MAX/Ny_<Nx_) return 2; Nx=Nx_; Ny=Ny_;
Matr=new complex[Nx*Ny]; if(Matr==NULL) return 3; return 0;
} 

//-----------------------------------------------------------------------------------------------------------
// ��������.

BYTE clMatrix::IsOK(void) const
{
if(Nx<=0||Ny<=0) return 1; if(Matr==NULL) return 2; return 0;
}

//-----------------------------------------------------------------------------------------------------------
// ��������� ������� � ����������� ���������� �� ���������.

BYTE clMatrix::operator =(complex c)
{
int i,j,jSh;

if(IsOK()!=0) return 1; if(Nx!=Ny) return 2;
for(i=0;i<Nx;i++) { jSh=i*Ny; for(j=0;j<Ny;j++) if(i!=j) Matr[j+jSh]=Cmplx_0; else Matr[j+jSh]=c;}
return 0;
}

//-----------------------------------------------------------------------------------------------------------
// ��������� ������� �� ����� ����������� ����������.

BYTE clMatrix::operator <<(complex c)
{
int i,j,jSh;

if(IsOK()!=0) return 1; for(i=0;i<Nx;i++) { jSh=i*Ny; for(j=0;j<Ny;j++) Matr[j+jSh]=c;} return 0;
}

//-----------------------------------------------------------------------------------------------------------
// ���������� ������� � ������������� ����������.

BYTE clMatrix::operator =(clRowMatr const &RM)
{
int i,j,jSh;

if(IsOK()!=0) return 1; if(Nx!=Ny) return 2; if(RM.IsOK()!=0) return 3; if(RM.N!=Nx) return 4;
for(i=0;i<Nx;i++) { jSh=i*Ny; for(j=0;j<Ny;j++) if(i!=j) Matr[j+jSh]=Cmplx_0;}
for(i=0;i<Nx;i++) { jSh=i*(Ny+1); Matr[jSh]=RM.Vect[i];} return 0;
}

//-----------------------------------------------------------------------------------------------------------
// ���� �����.

clMatrix clMatrix::operator -()
{
int i,j,jSh; clMatrix M;

if(IsOK()!=0) return M; if(M.Alloc(Nx,Ny)!=0) return M;
for(i=0;i<Nx;i++) { jSh=i*Ny; for(j=0;j<Ny;j++) M.Matr[j+jSh]=-Matr[j+jSh];} return M;
}

//-----------------------------------------------------------------------------------------------------------
// ����������� ��������� ������� � ���������� 'c'.

clMatrix &clMatrix::operator +=(complex c)
{
int i,jSh,N;

N=MIN(Nx,Ny); if(N<=0) goto end; if(IsOK()!=0) goto end;
for(i=0;i<N;i++) { jSh=i*(Ny+1); Matr[jSh]+=c;}
end: return *this;
}

//-----------------------------------------------------------------------------------------------------------
// ��������� ��������� ������� � ���������� 'c'.

clMatrix &clMatrix::operator -=(complex c)
{
int i,jSh,N;

N=MIN(Nx,Ny); if(N<=0) goto end; if(IsOK()!=0) goto end;
for(i=0;i<N;i++) { jSh=i*(Ny+1); Matr[jSh]-=c;}
end: return *this;
}

//-----------------------------------------------------------------------------------------------------------
// ��������� �� ����������� �����.

clMatrix &clMatrix::operator *=(complex c)
{
int i,j,jSh;

if(IsOK()!=0) return *this;
#pragma omp parallel for if(Nx*Ny>MinNumMultParal/6) schedule(static) default(shared) private(i,j,jSh)
for(i=0;i<Nx;i++) { jSh=i*Ny; for(j=0;j<Ny;j++) Matr[j+jSh]*=c;} return *this;
}

//-----------------------------------------------------------------------------------------------------------
// ������� �� ����������� �����.

clMatrix &clMatrix::operator /=(complex c)
{
int i,j,jSh; complex cInv;

if(IsOK()!=0) return *this; cInv=Inv(c);
#pragma omp parallel for if(Nx*Ny>MinNumMultParal/6) schedule(static) default(shared) private(i,j,jSh)
for(i=0;i<Nx;i++) { jSh=i*Ny; for(j=0;j<Ny;j++) Matr[j+jSh]*=cInv;} return *this;
}

//-----------------------------------------------------------------------------------------------------------
// ��������� �� �������������� �����.

clMatrix &clMatrix::operator *=(double d)
{
int i,j,jSh;

if(IsOK()!=0) return *this;
#pragma omp parallel for if(Nx*Ny>MinNumMultParal/2) schedule(static) default(shared) private(i,j,jSh)
for(i=0;i<Nx;i++) { jSh=i*Ny; for(j=0;j<Ny;j++) Matr[j+jSh]*=d;} return *this;
}

//-----------------------------------------------------------------------------------------------------------
// ������� �� �������������� �����.

clMatrix &clMatrix::operator /=(double d)
{
int i,j,jSh; double v;

if(IsOK()!=0) return *this; if(d==0.) return *this; v=1./d;
#pragma omp parallel for if(Nx*Ny>MinNumMultParal/2) schedule(static) default(shared) private(i,j,jSh)
for(i=0;i<Nx;i++) { jSh=i*Ny; for(j=0;j<Ny;j++) Matr[j+jSh]*=v;} return *this;
}

//-----------------------------------------------------------------------------------------------------------
// ����������� ������������ �������.

clMatrix &clMatrix::operator +=(const clRowMatr &RM)
{
int i,jSh,N;

N=MIN(Nx,Ny); N=MIN(N,RM.N); if(N<=0) goto end; if(IsOK()!=0) goto end; if(RM.IsOK()!=0) goto end;
for(i=0;i<N;i++) { jSh=i*(Ny+1); Matr[jSh]+=RM.Vect[i];}
end: return *this;
}

//-----------------------------------------------------------------------------------------------------------
// ��������� ������������ �������.

clMatrix &clMatrix::operator -=(const clRowMatr &RM)
{
int i,jSh,N;

N=MIN(Nx,Ny); N=MIN(N,RM.N); if(N<=0) goto end; if(IsOK()!=0) goto end; if(RM.IsOK()!=0) goto end;
for(i=0;i<N;i++) { jSh=i*(Ny+1); Matr[jSh]-=RM.Vect[i];}
end: return *this;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// ��������� � ������� ��� ������ ��������� �������.

//-----------------------------------------------------------------------------------------------------------
// ��������.

clMatrix operator +(const clMatrix &M1,const clMatrix &M2)
{
int i,j,jSh,jv; clMatrix M3;

if(M1.IsOK()!=0||M2.IsOK()!=0||M1.Nx!=M2.Nx||M1.Ny!=M2.Ny) goto end;
if(M3.Alloc(M1.Nx,M1.Ny)!=0) goto end;
#pragma omp parallel for if(M1.Nx*M1.Ny>MinNumMultParal/3) schedule(static) default(shared) private(i,j,jv,jSh)
for(i=0;i<M3.Nx;i++) { jSh=i*M3.Ny; for(j=0;j<M3.Ny;j++) { jv=j+jSh; M3.Matr[jv]=M1.Matr[jv]+M2.Matr[jv];}} 
end: return M3;
}

//-----------------------------------------------------------------------------------------------------------
// ���������.

clMatrix operator -(const clMatrix &M1,const clMatrix &M2)
{
int i,j,jSh,jv; clMatrix M3;

if(M1.IsOK()!=0||M2.IsOK()!=0||M1.Nx!=M2.Nx||M1.Ny!=M2.Ny) goto end;
if(M3.Alloc(M1.Nx,M1.Ny)!=0) goto end;
#pragma omp parallel for if(M1.Nx*M1.Ny>MinNumMultParal/3) schedule(static) default(shared) private(i,j,jv,jSh)
for(i=0;i<M3.Nx;i++) { jSh=i*M3.Ny; for(j=0;j<M3.Ny;j++) { jv=j+jSh; M3.Matr[jv]=M1.Matr[jv]-M2.Matr[jv];}}
end: return M3;
}

//-----------------------------------------------------------------------------------------------------------
// �������� ��������� ������.

clMatrix operator *(const clMatrix &M1,const clMatrix &M2)
{
int i,j,k,jSh1,jSh2,jSh3; complex s; clMatrix M3;

if(M1.IsOK()!=0) goto end; if(M2.IsOK()!=0) goto end; if(M1.Ny!=M2.Nx) goto end;
if(M3.Alloc(M1.Nx,M2.Ny)!=0) goto end;

#pragma omp parallel for schedule(static) default(shared) private(i,j,k,jSh1,jSh2,jSh3,s)
for(i=0;i<M1.Nx;i++) { jSh1=i*M1.Ny; jSh3=i*M3.Ny;
for(j=0;j<M2.Ny;j++) {
s=Cmplx_0; for(k=0;k<M1.Ny;k++) { jSh2=k*M2.Ny; s+=M1.Matr[k+jSh1]*M2.Matr[j+jSh2];}
M3.Matr[j+jSh3]=s;}}
end: return M3;
}

//-----------------------------------------------------------------------------------------------------------
// ����������������� �������.

clMatrix operator ~(const clMatrix &M_)
{
int i,j,iSh,jSh,N; clMatrix M;

if(M_.IsOK()!=0) goto end; if(M_.Nx!=M_.Ny) goto end; N=M_.Nx; if(M.Alloc(N,N)!=0) goto end;
for(i=0;i<N;i++) { jSh=i*N; for(j=0;j<N;j++) { iSh=j*N; M.Matr[i+iSh]=M_.Matr[j+jSh];}}
end: return M;
}

//-----------------------------------------------------------------------------------------------------------
// ����������� ��������� ������� � ���������� 'c'.

clMatrix operator +(const clMatrix &M_,complex c)
{
int i,j,jv,jSh,N; clMatrix M;

if(M_.IsOK()!=0) goto end; if(M_.Nx!=M_.Ny) goto end; N=M_.Nx; if(M.Alloc(N,N)!=0) goto end;
for(i=0;i<N;i++) { jSh=i*N; for(j=0;j<N;j++) { jv=j+jSh; M.Matr[jv]=M_.Matr[jv]; if(j==i) M.Matr[jv]+=c;}}
end: return M;
}

//-----------------------------------------------------------------------------------------------------------
// ����������� ��������� ������� � ���������� 'c'.

clMatrix operator +(complex c,const clMatrix &M_)
{
int i,j,jv,jSh,N; clMatrix M;

if(M_.IsOK()!=0) goto end; if(M_.Nx!=M_.Ny) goto end; N=M_.Nx; if(M.Alloc(N,N)!=0) goto end;
for(i=0;i<N;i++) { jSh=i*N; for(j=0;j<N;j++) { jv=j+jSh; M.Matr[jv]=M_.Matr[jv]; if(j==i) M.Matr[jv]+=c;}}
end: return M;
}

//-----------------------------------------------------------------------------------------------------------
// ��������� ��������� ������� � ���������� 'c'.

clMatrix operator -(const clMatrix &M_,complex c)
{
int i,j,jv,jSh,N; clMatrix M;

if(M_.IsOK()!=0) goto end; if(M_.Nx!=M_.Ny) goto end; N=M_.Nx; if(M.Alloc(N,N)!=0) goto end;
for(i=0;i<N;i++) { jSh=i*N; for(j=0;j<N;j++) { jv=j+jSh; M.Matr[jv]=M_.Matr[jv]; if(j==i) M.Matr[jv]-=c;}}
end: return M;
}

//-----------------------------------------------------------------------------------------------------------
// ��������� �� ��������� ������� � ���������� 'c'.

clMatrix operator -(complex c,const clMatrix &M_)
{
int i,j,jv,jSh,N; clMatrix M;

if(M_.IsOK()!=0) goto end; if(M_.Nx!=M_.Ny) goto end; N=M_.Nx; if(M.Alloc(N,N)!=0) goto end;
for(i=0;i<N;i++) { jSh=i*N; for(j=0;j<N;j++) { jv=j+jSh; M.Matr[jv]=-M_.Matr[jv]; if(j==i) M.Matr[jv]+=c;}}
end: return M;
}

//-----------------------------------------------------------------------------------------------------------
// ��������� ������� �� �����.

clMatrix operator *(const clMatrix &M_,complex c)
{
int i,j,jSh; clMatrix M;

if(M_.IsOK()!=0) goto end; if(M.Alloc(M_.Nx,M_.Ny)!=0) goto end;
#pragma omp parallel for if(M.Nx*M.Ny>MinNumMultParal/6) schedule(static) default(shared) private(i,j,jSh)
for(i=0;i<M.Nx;i++) { jSh=i*M.Ny; for(j=0;j<M.Ny;j++) M.Matr[j+jSh]=M_.Matr[j+jSh]*c;}
end: return M;
}

//-----------------------------------------------------------------------------------------------------------
// ��������� ������� �� �����.

clMatrix operator *(complex c,const clMatrix &M_)
{
int i,j,jSh; clMatrix M;

if(M_.IsOK()!=0) goto end; if(M.Alloc(M_.Nx,M_.Ny)!=0) goto end;
#pragma omp parallel for if(M.Nx*M.Ny>MinNumMultParal/6) schedule(static) default(shared) private(i,j,jSh)
for(i=0;i<M.Nx;i++) { jSh=i*M.Ny; for(j=0;j<M.Ny;j++) M.Matr[j+jSh]=M_.Matr[j+jSh]*c;}
end: return M;
}

//-----------------------------------------------------------------------------------------------------------
// ���������������� �������.

BYTE Transp(const clMatrix &M,clMatrix &MT)
{
int i,j,iSh,jSh,N;

if(M.IsOK()!=0) return 1; if(M.Nx!=M.Ny) return 2; N=M.Nx;
if(MT.IsOK()!=0) return 3; if(MT.Nx!=N) return 4; if(MT.Ny!=N) return 5;
for(i=0;i<N;i++) { jSh=i*N; for(j=0;j<N;j++) { iSh=j*N; MT.Matr[i+iSh]=M.Matr[j+jSh];}} return 0;
}

//-----------------------------------------------------------------------------------------------------------
// ��������� ���������� �������.

BYTE Hermit(const clMatrix &M,clMatrix &MH)
{
int i,j,iSh,jSh,N;

if(M.IsOK()!=0) return 1; if(M.Nx!=M.Ny) return 2; N=M.Nx;
if(MH.IsOK()!=0) return 3; if(MH.Nx!=N) return 4; if(MH.Ny!=N) return 5;
for(i=0;i<N;i++) { jSh=i*N; for(j=0;j<N;j++) { iSh=j*N; MH.Matr[i+iSh]=conj(M.Matr[j+jSh]);}} return 0;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// ��������� ��� ������� ��������� ������� � �������.

//-----------------------------------------------------------------------------------------------------------
// ��������� ������� �� �������.

clRowMatr operator *(const clMatrix &M,const clRowMatr &R)
{
int i,j,jSh; complex s; clRowMatr RM;

if(M.IsOK()!=0) goto end; if(R.IsOK()!=0) goto end; if(M.Ny!=R.N) goto end; if(RM.Alloc(M.Nx)!=0) goto end;

#pragma omp parallel for if(M.Nx*M.Ny>MinNumMultParal/8) schedule(static) default(shared) private(i,j,jSh,s)
for(i=0;i<M.Nx;i++) { jSh=i*M.Ny;
s=Cmplx_0; for(j=0;j<M.Ny;j++) s+=M.Matr[j+jSh]*R.Vect[j]; RM.Vect[i]=s;}
end: return RM;
}

//-----------------------------------------------------------------------------------------------------------
// �������� ������� � ������������ ��������.

clMatrix operator +(const clMatrix &M,const clRowMatr &R)
{
int i,j,jSh,N; clMatrix MF;

if(M.IsOK()!=0) goto end; if(M.Nx!=M.Ny) goto end; N=M.Nx; if(R.IsOK()!=0) goto end; if(R.N!=N) goto end;
if(MF.Alloc(N,N)!=0) goto end;
for(i=0;i<N;i++) { jSh=i*N;
for(j=0;j<N;j++) { MF.Matr[j+jSh]=M.Matr[j+jSh]; if(i==j) MF.Matr[j+jSh]+=R.Vect[i];}}
end: return MF;
}

//-----------------------------------------------------------------------------------------------------------
// �������� ������������ ������� � ��������.

clMatrix operator +(const clRowMatr &R,const clMatrix &M)
{
int i,j,jSh,N; clMatrix MF;

if(M.IsOK()!=0) goto end; if(M.Nx!=M.Ny) goto end; N=M.Nx; if(R.IsOK()!=0) goto end; if(R.N!=N) goto end;
if(MF.Alloc(N,N)!=0) goto end;
for(i=0;i<N;i++) { jSh=i*N;
for(j=0;j<N;j++) { MF.Matr[j+jSh]=M.Matr[j+jSh]; if(i==j) MF.Matr[j+jSh]+=R.Vect[i];}}
end: return MF;
}

//-----------------------------------------------------------------------------------------------------------
// ��������� ������������ ������� �� �������.

clMatrix operator -(const clMatrix &M,const clRowMatr &R)
{
int i,j,jSh,N; clMatrix MF;

if(M.IsOK()!=0) goto end; if(M.Nx!=M.Ny) goto end; N=M.Nx; if(R.IsOK()!=0) goto end; if(R.N!=N) goto end;
if(MF.Alloc(N,N)!=0) goto end;
for(i=0;i<N;i++) { jSh=i*N;
for(j=0;j<N;j++) { MF.Matr[j+jSh]=M.Matr[j+jSh]; if(i==j) MF.Matr[j+jSh]-=R.Vect[i];}}
end: return MF;
}

//-----------------------------------------------------------------------------------------------------------
// ��������� ������� �� ������������ �������.

clMatrix operator -(const clRowMatr &R,const clMatrix &M)
{
int i,j,jSh,N; clMatrix MF;

if(M.IsOK()!=0) goto end; if(M.Nx!=M.Ny) goto end; N=M.Nx; if(R.IsOK()!=0) goto end; if(R.N!=N) goto end;
if(MF.Alloc(N,N)!=0) goto end;
for(i=0;i<N;i++) { jSh=i*N;
for(j=0;j<N;j++) { MF.Matr[j+jSh]=-M.Matr[j+jSh]; if(i==j) MF.Matr[j+jSh]+=R.Vect[i];}}
end: return MF;
}

//-----------------------------------------------------------------------------------------------------------
// ��������� ���������� ������� �� ������������ �������.

clMatrix operator %(const clMatrix &M,const clRowMatr &R)
{
int i,j,jSh,N; clMatrix MF;

if(M.IsOK()!=0) goto end; if(M.Nx!=M.Ny) goto end; N=M.Nx; if(R.IsOK()!=0) goto end; if(R.N!=N) goto end;
if(MF.Alloc(N,N)!=0) goto end;
#pragma omp parallel for if(M.Nx*M.Ny>MinNumMultParal/7) schedule(static) default(shared) private(i,j,jSh)
for(j=0;j<N;j++) for(i=0;i<N;i++) { jSh=i*N; MF.Matr[j+jSh]=M.Matr[j+jSh]*R.Vect[j];}
end: return MF;
}

//-----------------------------------------------------------------------------------------------------------
// ��������� ������������ ������� �� ���������� �������.

clMatrix operator %(const clRowMatr &R,const clMatrix &M)
{
int i,j,jSh,N; clMatrix MF;

if(M.IsOK()!=0) goto end; if(M.Nx!=M.Ny) goto end; N=M.Nx; if(R.IsOK()!=0) goto end; if(R.N!=N) goto end;
if(MF.Alloc(N,N)!=0) goto end;
#pragma omp parallel for if(M.Nx*M.Ny>MinNumMultParal/7) schedule(static) default(shared) private(i,j,jSh)
for(i=0;i<N;i++) { jSh=i*N; for(j=0;j<N;j++) MF.Matr[j+jSh]=M.Matr[j+jSh]*R.Vect[i];}
end: return MF;
}
