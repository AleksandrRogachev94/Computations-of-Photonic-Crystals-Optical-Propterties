/*............................................................................................................

clTask,strMaterial,clMaterialAll,clLayGeom,strLay,clCrystal

............................................................................................................*/

#include <conio.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "Headers.h"
#include "Matrix.h"

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ����������� �����.

typedef BYTE (*funWinSmooth)(float *Weight,int nPoints); // ������� ���� ����������� �������.


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ��������� ��������������� ������������� ��� ������ ����� �����.

struct strEps
{
double eps; double wlen; 
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ��������� ���������.

struct strMaterial
{
char *name; // �������� ���������.
BYTE flMeth; // ���� ������� ������� eps. flMeth=0 - ������� ��� ���������� ����� �����. flMeth=1 - ������� � ������� ������ �����. flMeth=2,3 - ������� � ������� ������������. 2-n, 3-eps.
complex eps; // ��������������� ������������� ������� ���������.
double WaveLength; // ����� �����, ��� ������� ������ ��������������� �������������.
char chr; // ������, ������������ ��� ����������� ������� ���������.
struct strEps *nReal; // ���������� ����������� ���������. �������������� �����. (������������� ��� flMeth=2;3)
struct strEps *nImag; // ���������� ����������� ���������. ������ �����. (������������� ��� flMeth=2;3)
int Size_nReal,Size_nImag; // ������� ��������.

strMaterial(void); // �����������.
~strMaterial(void); // ����������.
void ZeroDisp(void); // ��������� �������� nImag,nReal.
void FreeDisp(void); // ������������ ������.
BYTE Read(FILE *fp); // ������ ������ � ��������� �� �����.
BYTE WriteInfo(FILE *fp); // ������ � ���� ���������� � ���������.
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ����� ������� ���� ����������, ������������ � ��������� (���� ������ ����������).

class clMaterialAll
{
void Zero(void); // ��������� ����������.
void Free(void); // ������������ ������.

public:

int N; // ������ �������.
struct strMaterial *Mat; // ������ ���� ����������.

clMaterialAll(void); // �����������.
~clMaterialAll(void); // ����������.
BYTE Alloc(int N_); // ��������� ������ ��� �������.
BYTE IsOK(void) const; // ��������.
BYTE Read(char *File_Name); // ������ ������ � ���������� �� �����.
int GetNum(char chr,double WaveLength) const; // ����� ������� ��������� � ����������� ��� ������ � ������ ����.
BYTE SetEps(double wLength); // ��������� eps ���������� �� ������ ����� ����� ��������� �����.
BYTE WriteInfo(FILE *fp); // ������ � ���� ���������� � ����������.
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ����� ���������� ������.

class clTask
{
void Zero(void); // ��������� ���������� � ������������� ����������.
void Free(void); // ������������ ������.
void ZeroFNames(void); // ��������� ���������� �� ������ - �������� ������.
void FreeFNames(void); // �������� ����� - �������� ������ �� ������.
BYTE ReadFNames(FILE *file); // ������ ��� ������.

public:

int L1,M1; // ������� ��� �������������� �����.
BYTE flTypComp; // ���� ���� �������. 0 - ������ ������������� ���� ��� �������� ����� �����. 1 - ������ ��������. 2 - ������ ������� ���������.
double wLength,wLength_St,wLength_Fin,dwLength; // ����� ����� ��������� �����.
double Theta,Theta_St,Theta_Fin,dTheta; // ���� ������� ����� (� �������).
BYTE flPol; // ���� ���� ����������� �������� �����. 0 - p, 1 - s.
complex EpsInc,EpsOut; // ��������������� ������������� ����, ������ ������, � ���� ������ ����.
char chEpsInc,chEpsOut; // ��������������� ������������� ����, ������ ������, � ���� ������ ����.
double dz; // ������� ��� ���������� ������������� ����� �� ��� z (���������).
char *FName_Cryst; // ����� ��� ��� ���� ������ ������������.
char *FName_Mat; // ��� ����� ����������.

clTask(void); // �����������.
~clTask(void); // ����������.
BYTE Read(char *File_Name); // ������ ������ �� �����.
BYTE SetEpsIncOut(clMaterialAll &Matter); // ��������� �������� ��������������� ������������� ����, ������ ������, � ���� ������ ����. 
BYTE WriteInfo(FILE *fp); // ������ ���������� ������ � ����.
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ����� ��������� ������� ���� ��������� ���������.

class clLaySmplGeom
{
void Zero(void); // ��������� ����������.
void Free(void); // ������������ ������.

public:

int laySize; // ������ ������� ������� �� ��� x.
double Period; // ������ ��������� �� ��� x � ����������.
char *GeomSymbs; // ������ ��������� ��������� � �������� ����������.
int *Geom; // ������ ��������� �������.

clLaySmplGeom(void); // �����������.
~clLaySmplGeom(void); // ����������.
BYTE Alloc(int laySize_); // ��������� ������ ��� �������.
BYTE IsOK(void) const; // ��������.
BYTE Read(FILE *fp,const clMaterialAll &Matter,double WaveLength,int numLay); // ������ ��������� ������� �� �����.
BYTE IsHomogen(void); // �������� ������������ �������.
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ����������.

class clCrystalSmpl; // ����� �������� ����� ��������� ���������.

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ��������� ������ ������� ���� �� ��� z ��������� ���������.

struct strLaySmpl
{
int numLay; // ����� �������.
clLaySmplGeom LayGeom; // ��������� �������.
int L1,M1; // L1 - ����� ������������� ��������, M1 - �������������.
int M; // ������� ������. ��� ����� ��������������.
int N; // ������ �������� ��������� ����� (L1+M1+1).
int N1; // ������ �������� 'Eps_Four', 'EpsInv_Four' 2*(L1+M1)+1.
int NF; // ������ �������� ��� �������� �������������� �����.
complex *Eps_Four,*EpsInv_Four; // ��������� ����� ������������� ��� 'eps' � '1/eps'.
BYTE flHomogen; // ���� ������������ ����.
complex EpsHomogen; // ��������������� ������������� ��� ����������� ����.
clMatrix Ek,Ak,EkInv,AkInv;
clMatrix Mk; // ������� ���������� ��������� ��������� ����� �������������� �����.
clRowMatr EigLamMk; // ������ ��� ����������� �������� ������� 'Mk' � ���� �������.
clMatrix EigVectMk,EigVectMkInv; // ������� ��� �������� ����������� �������� ������� 'Mk'.
clMatrix GamP[2][2],GamPInv[2][2]; // ����� ������ 'EigVectMk','EigVectMkInv'.

strLaySmpl(void); // �����������.
~strLaySmpl(void); // ����������.
void Zero(void); // ��������� ����������.
void Free(void); // ������������ ������.
BYTE IsOK(void) const; // ��������.
BYTE IsOKFourier(void); // �������� ������� ����� ��������������.
BYTE IsOKGamma(void) const; // �������� ������ ������ 'Gamma'.
BYTE Read(FILE *fp,const clCrystalSmpl &Samples); // ������ ������ � ���� �� �����.
BYTE IsHomogen(void); // ��������� ����� ������������ ����.
BYTE FastFTLay(int L1_,int M1_,const clMaterialAll &Matter); // ������� �������������� ����� ��� ����.
BYTE SetParHomogen(int L1_,int M1_,const clMaterialAll &Matter); // ��������� ���������� ��� ����������� ����.
BYTE CompMk(const clTask &Task); // ���������� ������� Mk.
BYTE FindEigVectMk(void); // ���������� ����������� �������� ������� Mk.
BYTE FindEigVectMkHomogen(const clTask &Task); // ���������� ����������� �������� � �������� ������� Mk (������ ����������� ����).
BYTE FixFormEig(void); // ���������� ������� ����������� �������� � ������� � ������ ������ ���� (�� ����������� �������������� �����).
BYTE FindGamma(void); // ���������� ������ ������ 'Gamma'.
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ����� �������� ����� ��������� ���������.

class clCrystalSmpl
{
void Zero(void); // ��������� ����������.
void Free(void); // ������������ ������.

public:

int laysNum; // ����� ��������.
double Period; // ������ ��������� �� ��� x � ����������.
struct strLaySmpl *Lays; // ������ ���� ��������.
clMaterialAll Matter; // ���������, �� ������� ����������� �������.
clTask Task; // ��������� ������.

clCrystalSmpl(void); // �����������.
~clCrystalSmpl(void); // ����������.
BYTE Alloc(int laysNum_); // ��������� ������ ��� �������.
BYTE IsOK(void) const; // ��������.
BYTE IsOKGamma(void); // �������� ������ ������ 'Gamma' �� ���� �����.
BYTE ReadInpData(void); // ������ ������ �� �������� �� �����.
BYTE FastFTAll(void); // ������� �������������� ����� ��� ���� ��������.
BYTE CompMk(void); // ���������� ������ 'Mk' ��� ���� ��������.
BYTE FindEigVectMk(void); // ���������� ����������� �������� ������ 'Mk' ��� ���� ��������.
strLaySmpl *Get(int num); // ��������� ��������� �� ������� ���� � ������� 'num'.
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ��������� ������ ���� �� ��� z ��������� ���������.

struct strLay
{
int numLay; // ����� ���� � ���������.
int numSmpl; // ����� �������, �������� ������������� ������ ����.
double depth; // ������� ���� � ����������.
int N; // ������ �������� ��������� ����� (L1+M1+1). ������ �������� DExp00Inv,DExp11.
clRowMatr DExp00Inv,DExp11; // ������������ ������� � ������������: �������� � 'Gamma00' � 'Gamma11'.
clMatrix GK[4]; // ������� 'Gk' ��� ������� ����� ������ ��������� (�������������� ����� ���������� � 'P^-1','DExp','P' � ����� �������� � 'P').
clRowMatr EFour[4],HFour[4]; // ������� ������������� ����� ���������� ����� � � H �������������� �� �������� �������, ����� ���������� �� 'P', ����� ��������� �� 'P^-1' � �� ������ ������� ����.
double dz; // ������� ������� ���� ��� ���������� ���� ������ ����.
int nDiv; // ����� ��������� ���� ��� ���������� ������������� �����.
clRows EDistrFour,HDistrFour; // ����� ������������ ����� E � H. ��� 'p' ����� Ex,Hy, ��� 's' ����� Ey,Hx.
clRows EDistr,HDistr,ZDistr; // ������������� ��������� ����� E,H � Ez ��� Hz � ����������� �� �����������. ��� 'p' ����� Ex,Hy,Ez, ��� 's' ����� Ey,Hx,Hz.

strLay(void); // �����������.
void Zero(void); // ��������� ����������.
BYTE IsOK(void); // ��������.
BYTE FindGammaExp(const strLaySmpl &Smpl); // ���������� ������ DExp00Inv, DExp11.
BYTE IsOKGammaExp(void); // �������� ������ ������ 'Gamma'.
BYTE FindFullFieldDistr(clCrystalSmpl &Samples); // ���������� ������������� ���� ������ ���� ���������.
BYTE CompDz(double dz_); // ���������� 'dz' �� ����������, ������������ 'dz' � ��������� ������ ��� �������� �������� ������������� ����� ����� E � H.
BYTE FindDistrFourEH(const strLaySmpl &Smpl); // // ���������� ������������� ����� ���������� ����� � � H � ����.
BYTE FindDistrEH(const strLaySmpl &Smpl,clTask &Task); // ���������� ������������� ����� � � H � ���� (���������� ��������� �������������� �����).
BYTE WriteInfo(FILE *fp,const strLaySmpl &Smpl); // ������ � ���� ���������� � ����.
BYTE WriteFieldDistr(FILE *fp); // ������ � ���� ������������� �����.
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ����� ��������� ���������.

class clCrystal
{
void Zero(void); // ��������� ����������.
void Free(void); // ������������ ������.

public:

clCrystalSmpl Samples; // ���������� ��� ���� �������� ����� ���������.
int laysNum; // ����� �����.
struct strLay *Lays; // ������ ���� �����.
clRowMatr CInc,CIncInv,COut; // ����� ������ '�11','�11^-1' ��� �����, ������ ������ ����, � '�11' - ��� �����, ���� �������� ����.
clRowMatr J,T,R; // �������, ��������������� �������� ����� � �������� ��������� � ��������� ������.

clCrystal(void); // �����������.
~clCrystal(void); // ����������.
BYTE Alloc(int laysNum_); // ��������� ������ ��� �������.
BYTE IsOK(void) const; // ��������.
BYTE ReadInpData(char *FileName_Task); // ������ �������� ������.
BYTE ReadGeom(char *FileName_Cryst); // ������ ���������� � ��������� �� �����.
BYTE FindGammaExp(void); // ���������� ������ 'DExp00Inv', 'DExp11' ��� ���� ����.
BYTE IsOKGamma(void); // �������� ������ ������ 'Gamma' �� ���� �����.
BYTE CompC(void); // ���������� ������ 'C'.
BYTE IsOK_C(void); // �������� ������ ������ 'C'.
BYTE FindRTW(void); // ���������� �������� ��������� � ��������� ����.
BYTE FindBoundsFieldDistr(void); // ���������� �������� �������� ���������� ����� ������������� ���� �� �������� ����.
BYTE FindFullFieldDistr(void); // ���������� ������������� ���� ������ ���������.
BYTE WriteInfo(FILE *fp); // ������ � ���� ���������� � ��������� � ������.
BYTE AngDistrOutput(void); // ����� �������� ������������� �� �������� ���������.
BYTE CompPartRefTr(double *pRef,double *pTr); // ����� �������� ������������� �� �������� ���������.
BYTE SpectrOutput(int n,double *T_Sp,double *R_Sp,double *A_Sp); // ����� � ���� ���� ��������.
BYTE AnglDispOutput(int nRows,int n,clRows &T_Sp,clRows &R_Sp,clRows &A_Sp); // ����� � ���� ���� ��������.
BYTE FieldOutput(void); // ����� � ���� ������������� ����� � ���������.
BYTE Compute(double *Time); // ������ ��� ��������� ���������.
BYTE ComputeSpectr(double *Time); // ������ �������� �����������, ��������� � ���������� ��� ��������� ���������.
BYTE ComputeAnglDisp(double *Time); // ������ ������� ���������� ��� ��������� ���������.
};
