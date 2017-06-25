#include "stdafx.h"

#include <conio.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#include "Headers.h"
#include "cmplx.h"
#include "Properties.h"

//------------------------------------------------------------------------------------------------------------
// Вывод в файл данных формулы Друде для золота.

BYTE WriteDrudeData(double st,double fin,double step)
{
complex GetEps_Drude(double Lambda,double Eps0,double gamma,double wp);

BYTE err; double v,lambda; FILE *fp_r,*fp_i; complex eps;

if(st<=0.) return 1; if(fin<=0.) return 2; if(step<=0.) return 3;
if(st>fin) { v=st; st=fin; fin=v;} if(step>fin-st) return 4;

fp_r=fp_i=NULL; err=0; lambda=st;
fp_r=fopen("DrudeData_real.dr","w"); if(fp_r==NULL) { err=5; goto end;}
fp_i=fopen("DrudeData_imag.dr","w"); if(fp_i==NULL) { err=6; goto end;}

fprintf(fp_r,"permittivity(real) wavelength\n"); fprintf(fp_i,"permittivity(imag) wavelength\n");
while(lambda<=fin) {
eps=GetEps_Drude(lambda,Eps0_Gold,gamma_Gold,wp_Gold);
fprintf(fp_r,"%lf %lf\n",lambda,eps.re); fprintf(fp_i,"%lf %lf\n",lambda,eps.im);
lambda+=step;
}

end: SAFE_CLOSE(fp_r); SAFE_CLOSE(fp_i); return err;
}

//------------------------------------------------------------------------------------------------------------
// Вывод в файл экспериментальных данных диэлектрической проницаемости для золота.

BYTE WriteExpData(const clMaterialAll &Matter,double st,double fin,double step)
{
complex EpsInterpol_lin(int Size_nReal_,strEps *nReal_,int Size_nImag_,strEps *nImag_,double wLength,BYTE flMeth_); // Интерполяции экспериментальных данных коэффициента преломления.

BYTE err; int i; double v,lambda; FILE *fp_r,*fp_i; complex eps; strMaterial *pM;

if(st<=0.) return 1; if(fin<=0.) return 2; if(step<=0.) return 3;
if(st>fin) { v=st; st=fin; fin=v;} if(step>fin-st) return 4;
if(Matter.IsOK()!=0) return 5;

fp_r=fp_i=NULL; pM=NULL; err=0; lambda=st;
fp_r=fopen("ExpData_real.exper","w"); if(fp_r==NULL) { err=6; goto end;}
fp_i=fopen("ExpData_imag.exper","w"); if(fp_i==NULL) { err=7; goto end;}

for(i=0;i<Matter.N;i++) if(Matter.Mat[i].flMeth>=2) { pM=Matter.Mat+i; break;} // ~выбор золота.
if(pM==NULL) { err=8; goto end;}

fprintf(fp_r,"permittivity(real) wavelength\n"); fprintf(fp_i,"permittivity(imag) wavelength\n");
while(lambda<=fin) {
eps=EpsInterpol_lin(pM->Size_nReal,pM->nReal,pM->Size_nImag,pM->nImag,lambda,pM->flMeth);
fprintf(fp_r,"%lf %lf\n",lambda,eps.re); fprintf(fp_i,"%lf %lf\n",lambda,eps.im);
lambda+=step;
}

end: SAFE_CLOSE(fp_r); SAFE_CLOSE(fp_i); return err;
}

//------------------------------------------------------------------------------------------------------------
// Главная программа.

int main(void)
{
BYTE err=0; int n_Th; double E,Time; clCrystal Cryst;
complex eps;

// Чтение исходных данных.
if(Cryst.ReadInpData(File_Name_Task)!=0) { err=1; printf("Error in reading Input Data !\n"); goto end;}

//WriteDrudeData(190,2480,5); WriteExpData(Cryst.Samples.Matter,190,2480,5);

E=hPlank*cVel*2.*M_PI/(Cryst.Samples.Task.wLength*1.e-9*eChrg);
printf("Energy of photon (eV): %lf\n\n",E);

printf("This computer has %d cores. How many threads do you want to use?\n",omp_get_num_procs());
scanf("%d",&n_Th); if(n_Th<=0) { err=2; goto end;} omp_set_num_threads(n_Th);
omp_set_nested(1); // Включение вложенного параллелизма.

// Расчёт для фотонного кристалла.
switch(Cryst.Samples.Task.flTypComp) {
default: break;
case 0: err=Cryst.Compute(&Time); break;
case 1: err=Cryst.ComputeSpectr(&Time); break;
case 2: err=Cryst.ComputeAnglDisp(&Time); break;
}
end: if(err==0) printf("\nAll done ! Time of computations : %lf seconds\nResults are in output files",Time);
else printf("error=%d",err);
_getch(); return 0;
}