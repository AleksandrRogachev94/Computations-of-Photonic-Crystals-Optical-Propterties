/*...........................................................................................................
// Функции для вычисления физических параметров.
GetEps_Drude,EpsInterpol_lin,

// Структура материала.
struct strMaterial { strMaterial,~strMaterial,Read,WriteInfo},

// Класс массива всех материалов, используемых в кристалле (база данных материалов).
class clMaterialAll { clMaterialAll,~clMaterialAll,Zero,Free,Alloc,IsOK,Read,GetNum,WriteInfo},

// Класс параметров задачи.
class clTask { clTask,~clTask,Zero,ZeroFNames,FreeFNames,ReadFNames,Read,WriteInfo},

// Класс геометрии одного слоя фотонного кристалла.
class clLayGeom { clLayGeom,~clLayGeom,Zero,Free,Alloc,IsOK,Read,IsHomogen},

// Структура одного образца по оси z фотонного кристалла.
struct strLaySmpl { strLaySmpl,~strLaySmpl,Zero,Free,IsOK,IsOKFourier,IsOKGamma,Read,IsHomogen,FastFTLay,
SetParHomogen,CompMk,FindEigVectMk,FindEigVectMkHomogen,FixFormEig,FindGamma,WriteInfo},

// Класс фотонного кристалла.
class clCrystalSmpl { clCrystalSmpl,~clCrystalSmpl,Zero,Free,Alloc,IsOK,IsOKGamma,ReadInpData,
FastFTAll,CompMk,FindEigVectMk,Get},

// Структура одного слоя по оси z фотонного кристалла.
struct strLay { strLay,Zero,IsOK,FindGammaExp,IsOKGammaExp,FindFullFieldDistr,CompDz,FindDistrFourEH,
FindDistrEH,WriteInfo,WriteFieldDistr},
FillArrInvFour,SmoothSpectrum,GetWeightWinCos,GetWeightKaisBess,GetBesselI0,

// Класс фотонного кристалла.
class clCrystal { clCrystal,~clCrystal,Zero,Free,Alloc,IsOK,ReadInpData,ReadGeom,FindGammaExp,IsOKGamma,
CompC,IsOK_C,FindRTW,FindBoundsFieldDistr,FindFullFieldDistr,WriteInfo,AngDistrOutput,FieldOutput,Compute},
MessageSmallEps,MessageSingLam

...........................................................................................................*/

#include "stdafx.h"

#include <conio.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

#include "Headers.h"
#include "cmplx.h"
#include "Properties.h"

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Константы.

#define flFindPInv 1 // Флаг нахождения обратной матрицы 'PInv' матрицы собственных векторов как матрицы левых собственных векторов.

#define flSmooth 1 // Флаг сглаживания пространственного спектра при обратном преобразовании Фурье для нахождения полей.
#define CoeWinSmooth 0.7F // Относительная ширина интервала сглаживания на больших по модулю гармониках пространственного спектра.
#define coeWinKaisBess 4.54 // Коэффициент окна Кайзера-Бесселя.

// Типы окнных функций сглаживания спектра.
#define typFunSmoothCos 1 // Косинусоидальное сглаживание.
#define typFunKaisBess 2 // Окно Кайзера-Бесселя.

#define typFunSmooth typFunKaisBess // Задание типа окнной функции сглаживания спектра.

#define smValEps 1.e-6 // Порядок малости для диэлектрической проницаемости.
#define smValLam 1.e-8 // Порядок малости для собственных значений (компоненты 'Kz' волнового вектора).

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Функции.

char *GetFName(char *FName,char *ext); // Получение имени файла с расширением.

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Строки.

char *extCrystSamples=".smpl"; // Расширение файла для образцов слоёв кристалла.
char *extCrystal=".crys"; // Расширение файла для геометрии кристалла.
char *extOutAngle=".outa"; // Расширение файла с угловым распределением.
char *extOutFieldDistr=".outd"; // Расширение файла с распределением поля.
char *extOutR_Sp=".outrs";
char *extOutT_Sp=".outts";
char *extOutA_Sp=".outas"; // Расширения файлов со спектрами.
char *extOutR_AnglDisp=".outrad";
char *extOutT_AnglDisp=".outtad";
char *extOutA_AnglDisp=".outaad";// Расширения файлов с угловыми дисперсиями.

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Функции для вычисления физических параметров.

//-----------------------------------------------------------------------------------------------------------
// Вычисление диэлектрической проницаемости материала по формулам Друде.
// wp и gamma в терагерцах (10^12). Lambda в нанометрах.
// Формула верна для длин волн >650 нм. !

complex GetEps_Drude(double Lambda,double Eps0,double gamma,double wp)
{
double w;
if(Lambda<=0.) return Cmplx_0; if(wp<=0.) return Cmplx_0;
w=2*M_PI*cVel*(1.e-3)/Lambda; 

return complex(Eps0,0.)-wp*wp/complex(w*w,gamma*w);
}

//-----------------------------------------------------------------------------------------------------------
// Интерполяции экспериментальных данных коэффициента преломления.
// flMeth=2 --> в массивах коэффициент преломления. flMeth=3 --> в массивах диэлектрическая проницаемость.

complex EpsInterpol_lin(int Size_nReal,strEps *nReal,int Size_nImag,strEps *nImag,double wLength,BYTE flMeth)
{
int j,k,Size; double coe,b,Real,Imag; complex eps; strEps *n;

if(flMeth!=2) return Cmplx_0; if(nReal==NULL) return Cmplx_0; if(nImag==NULL) return Cmplx_0;
if(Size_nReal<=0) return Cmplx_0; if(Size_nImag<=0) return Cmplx_0; if(wLength<=0.) return Cmplx_0;

// Интерполяция действительной и мнимой частей отдельно.
for(k=0;k<2;k++) { if(k==0) { Size=Size_nReal; n=nReal;} else { Size=Size_nImag; n=nImag;}
if(wLength>n[Size-1].wlen||wLength<n[0].wlen) return Cmplx_0; // Проверка.
for(j=1;j<Size;j++) if(n[j].wlen>wLength) break;
coe=(n[j].eps-n[j-1].eps)/(n[j].wlen-n[j-1].wlen); b=n[j].eps-coe*n[j].wlen;
if(k==0) Real=coe*wLength+b; else Imag=coe*wLength+b;
}

switch(flMeth) { default: return Cmplx_0;
case 2: eps=complex(Real,Imag)*complex(Real,Imag); break;
case 3: eps=complex(Real,Imag); break;
}
return eps;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Структура материала.

//-----------------------------------------------------------------------------------------------------------
// Конструктор.

strMaterial::strMaterial(void)
{
flMeth=1; eps=Cmplx_1; WaveLength=0.; chr=' '; name=NULL; ZeroDisp();
name=new char[szNameSubst];
}

//-----------------------------------------------------------------------------------------------------------
// Деструктор.

strMaterial::~strMaterial(void)
{
SAFE_DELETE_ARR(name); FreeDisp();
}

//-----------------------------------------------------------------------------------------------------------
// Обнуление массивов nImag,nReal.

void strMaterial::ZeroDisp(void)
{
Size_nReal=Size_nReal=0; nReal=nImag=NULL;
}

//-----------------------------------------------------------------------------------------------------------
// Освобождение памяти.

void strMaterial::FreeDisp(void)
{
SAFE_DELETE_ARR(nReal); SAFE_DELETE_ARR(nImag);
}

//-----------------------------------------------------------------------------------------------------------
// Чтение данных о материале из файла.

BYTE strMaterial::Read(FILE *fp)
{
BYTE err; int i,k,*Size; strEps **n;

if(fp==NULL) return 1; if(name==NULL) return 2;
if(fscanf(fp,"%s",name)==EOF) return 3; if(fscanf(fp,"%d\n",&flMeth)==EOF) return 3; if(flMeth>2||flMeth<0) return 4;

err=0;
switch(flMeth) {
default: break;
case 0: if(fscanf(fp,"%lf\n%lf%lf",&WaveLength,&eps.re,&eps.im)==EOF) return 5; break;
case 1: break;
case 2: FreeDisp();
for(k=0;k<2;k++) { if(k==0) { Size=&Size_nReal; n=&nReal;} else { Size=&Size_nImag; n=&nImag;}
if(fscanf(fp,"%d\n",Size)==EOF) { err=6; goto end;}
(*n)=new strEps[*Size]; if((*n)==NULL) { err=7; goto end;}
for(i=0;i<(*Size);i++) if(fscanf(fp,"%lf %lf\n",&(*n)[i].wlen,&(*n)[i].eps)==EOF) { err=8; goto end;}
}
break;
}

if(fscanf(fp,"\n%c",&chr)==EOF) goto end;
end: return err;
}

//-----------------------------------------------------------------------------------------------------------
// Запись в файл информации о материале.

BYTE strMaterial::WriteInfo(FILE *fp)
{
if(fp==NULL) return 1;
if(fprintf(fp,"%s\n",name)<0) return 2;
//xxx mmm aaa if(fprintf(fp,"%lf %lf %d %d %d %c\n",eps.re,eps.im,R,G,B,chr)<0) return 3;
if(fprintf(fp,"%lf %lf %c\n",eps.re,eps.im,chr)<0) return 3;
return 0;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Класс массива всех материалов, используемых в кристалле.

//-----------------------------------------------------------------------------------------------------------
// Конструктор.

clMaterialAll::clMaterialAll(void)
{
Zero();
}

//-----------------------------------------------------------------------------------------------------------
// Деструктор.

clMaterialAll::~clMaterialAll(void)
{
Free();
}

//-----------------------------------------------------------------------------------------------------------
// Обнуление переменных.

void clMaterialAll::Zero(void)
{
Mat=NULL; N=0;
}

//-----------------------------------------------------------------------------------------------------------
// Освобождение памяти.

void clMaterialAll::Free(void)
{
SAFE_DELETE_ARR(Mat); N=0;
}

//-----------------------------------------------------------------------------------------------------------
// Выделение памяти для массива.

BYTE clMaterialAll::Alloc(int N_)
{
Free(); if(N_<=0) return 1; N=N_;
Mat=new strMaterial[N]; if(Mat==NULL) return 2; return 0;
} 

//-----------------------------------------------------------------------------------------------------------
// Проверки.

BYTE clMaterialAll::IsOK(void) const
{
if(N<=0) return 1; if(Mat==NULL) return 2; return 0;
}

//-----------------------------------------------------------------------------------------------------------
// Чтение данных о материалах из файла.

BYTE clMaterialAll::Read(char *File_Name)
{
BYTE err; int i,N_; FILE *fp;

fp=NULL; err=0;
fp=fopen(File_Name,"r"); if(fp==NULL) { printf("Can not open file: %s !\n",File_Name); err=1; goto end;}
if(fscanf(fp,"%d",&N_)==EOF) { err=2; goto end;} if(Alloc(N_)!=0) { err=3; goto end;}
for(i=0;i<N;i++) { if(fscanf(fp,"\n")==EOF) { err=4; goto end;} if(Mat[i].Read(fp)!=0) { err=5; goto end;}}
end: SAFE_CLOSE(fp); return err;
}

//-----------------------------------------------------------------------------------------------------------
// Запись в файл информации о материалах.

BYTE clMaterialAll::WriteInfo(FILE *fp)
{
int i;

if(fp==NULL) return 1; if(IsOK()!=0) return 2;
if(fprintf(fp,"%d\n",N)<0) return 3;
for(i=0;i<N;i++) if(Mat[i].WriteInfo(fp)!=0) return 4; fprintf(fp,"\n"); return 0;
}

//-----------------------------------------------------------------------------------------------------------
// Поиск символа материала и возвращение его номера в списке сред.

int clMaterialAll::GetNum(char chr,double WaveLength) const
{
int i,num; strMaterial *pM;

if(IsOK()!=0) return -1; num=-2;
for(i=0;i<N;i++) { pM=Mat+i; if(pM->chr!=chr) continue;
if(pM->flMeth==0) if(pM->WaveLength>0.) { if(fabs(pM->WaveLength-WaveLength)>1.e-3) continue;}
num=i; break;} return num;
}

//-----------------------------------------------------------------------------------------------------------
// Установка eps материалов на данной длине волны.

BYTE clMaterialAll::SetEps(double wLength)
{
complex GetEps_Drude(double Lambda,double Eps0,double gamma,double wp); // Вычисление диэлектрической проницаемости материала по формулам Друде.
complex EpsInterpol_lin(int Size_nReal_,strEps *nReal_,int Size_nImag_,strEps *nImag_,double wLength,BYTE flMeth_); // Интерполяции экспериментальных данных коэффициента преломления.

int i; strMaterial *pM;

if(wLength<=0.) return 1;

for(i=0;i<N;i++) { pM=Mat+i; 

switch(pM->flMeth) {
default: break;
case 1: pM->eps=GetEps_Drude(wLength,Eps0_Gold,gamma_Gold,wp_Gold);
if(pM->eps.re==0.&&pM->eps.im==0.) return 2; pM->WaveLength=wLength; break;
case 2: pM->eps=EpsInterpol_lin(pM->Size_nReal,pM->nReal,pM->Size_nImag,pM->nImag,wLength,pM->flMeth);
if(pM->eps.re==0.&&pM->eps.im==0.) { printf("Error in the interpolation of permittivities !\n"); getch(); return 3;}
pM->WaveLength=wLength; break;
}
}
return 0;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Класс параметров задачи.

//-----------------------------------------------------------------------------------------------------------
// Конструктор.

clTask::clTask(void)
{
Zero();
}

//-----------------------------------------------------------------------------------------------------------
// Деструктор.

clTask::~clTask(void)
{
FreeFNames();
}

//-----------------------------------------------------------------------------------------------------------
// Обнуление указателей и инициализация переменных.

void clTask::Zero(void)
{
ZeroFNames(); L1=50; M1=50; flPol=0; EpsInc=Cmplx_1; EpsOut=Cmplx_1; dz=0.;
wLength=wLength_St=wLength_Fin=dwLength=0.; flTypComp=0;
}

//-----------------------------------------------------------------------------------------------------------
// Обнуление указателей на строки - названия файлов.

void clTask::ZeroFNames(void)
{
FName_Cryst=FName_Mat=NULL;
}

//-----------------------------------------------------------------------------------------------------------
// Удаление строк - названий файлов из памяти.

void clTask::FreeFNames(void)
{
SAFE_DELETE_ARR(FName_Mat); SAFE_DELETE_ARR(FName_Cryst);
}

//-----------------------------------------------------------------------------------------------------------
// Чтение имён файлов.

BYTE clTask::ReadFNames(FILE *file)
{
BYTE ReadString(FILE *file,char *pStr,int lnStr); // Чтение строки из текстового файла.

BYTE i,err; size_t ln; char **pStr,*str;

FreeFNames(); if(file==NULL) return 1;
str=NULL; err=0;
str=new char[laySize_Max]; if(str==NULL) { err=2; goto end;} str[laySize_Max-1]='\0';

for(i=0;i<2;i++) {
switch(i) { default: continue;
case 0: pStr=&FName_Mat; break; // Название файла материалов.
case 1: pStr=&FName_Cryst; break; // Общее название файлов для кристалла.
}
if(ReadString(file,str,laySize_Max-1)!=0) { err=3; goto end;}
ln=strlen(str); if(ln==0) { err=4; goto end;}
*pStr=new char[ln+1]; (*pStr)[ln]='\0'; strcpy(*pStr,str);}

if(strchr(FName_Cryst,'.')!=0) { printf("Error: Common File Name should not have '.' !\n"); err=5; goto end;}

end: SAFE_DELETE_ARR(str); if(err!=0) FreeFNames(); return err;
}

//-----------------------------------------------------------------------------------------------------------
// Чтение данных из файла.

BYTE clTask::Read(char *File_Name)
{
BYTE err; FILE *fp; double v;

fp=NULL; err=0;
fp=fopen(File_Name,"r"); if(fp==NULL) { printf("Can not open file: %s !\n",File_Name); err=1; goto end;}
if(fscanf(fp,"%d\n",&flTypComp)==EOF) { err=2; goto end;} if(flTypComp<0||flTypComp>2) flTypComp=0;

switch(flTypComp) {
default: err=3; goto end;
case 0: if(fscanf(fp,"%lf\n",&wLength)==EOF) { err=4; goto end;} // Распределение.
if(wLength<=0.) { err=5; goto end;} break;
case 1: // Спектр.
if(fscanf(fp,"%lf %lf %lf\n",&wLength_St,&wLength_Fin,&dwLength)==EOF) { err=6; goto end;}
if(wLength_St<=0.||wLength_Fin<=0.||dwLength<=0.) { err=7; goto end;} // Проверки.
if(wLength_St>wLength_Fin) { v=wLength_St; wLength_St=wLength_Fin; wLength_Fin=v;} // Проверки.
if(wLength_St==wLength_Fin) { flTypComp=0; wLength=wLength_Fin;} // Проверки.
if(dwLength>wLength_Fin-wLength_St) dwLength=wLength_Fin-wLength_St; // Проверки.
wLength=wLength_St;
break;
case 2: // Угловое распределение.
if(fscanf(fp,"%lf %lf %lf\n",&wLength_St,&wLength_Fin,&dwLength)==EOF) { err=8; goto end;}
if(fscanf(fp,"%lf %lf %lf\n",&Theta_St,&Theta_Fin,&dTheta)==EOF) { err=9; goto end;}
if(wLength_St<=0.||wLength_Fin<=0.) { err=10; goto end;} if(wLength_Fin-wLength_St!=0&&dwLength==0) { err=11; goto end;} if(Theta_Fin-Theta_St!=0&&dTheta==0) { err=11; goto end;} // Проверки.
if(wLength_St>wLength_Fin) { v=wLength_St; wLength_St=wLength_Fin; wLength_Fin=v;}
if(wLength_St==wLength_Fin) wLength=wLength_Fin; if(Theta_St==Theta_Fin) Theta=Theta_Fin; // Проверки.
if(dwLength>wLength_Fin-wLength_St) dwLength=wLength_Fin-wLength_St; if(dTheta>Theta_Fin-Theta_St) dTheta=Theta_Fin-Theta_St; // Проверки.
wLength=wLength_St; Theta=Theta_St;
break;
}

if(flTypComp!=2) { if(fscanf(fp,"%lf\n",&Theta)==EOF) { err=12; goto end;} if(Theta<-90.||Theta>90.) { err=13; goto end;}}
if(fscanf(fp,"%d\n",&flPol)==EOF) { err=14; goto end;} if(flPol>=nPol) flPol=0;
if(fscanf(fp,"%d%d\n",&L1,&M1)==EOF) { err=15; goto end;} if(L1<=0||M1<=0) { err=16; goto end;}
if(fscanf(fp,"%c ",&chEpsInc)==EOF) { err=17; goto end;}
if(fscanf(fp,"%c\n",&chEpsOut)==EOF) { err=18; goto end;}
if(fscanf(fp,"%lf\n",&dz)==EOF) { err=19; goto end;} if(dz<=0.) { err=20; goto end;}

if(ReadFNames(fp)!=0) { err=21; goto end;}
end: SAFE_CLOSE(fp); return err;
}

//-----------------------------------------------------------------------------------------------------------
// Запись параметров задачи в файл.

BYTE clTask::WriteInfo(FILE *fp)
{
if(fp==NULL) return 1;
fprintf(fp,"%lf\n",wLength); fprintf(fp,"%lf\n",Theta);
fprintf(fp,"%d\n",flPol);
fprintf(fp,"%lf %lf\n",EpsInc); fprintf(fp,"%lf %lf\n",EpsOut);
fprintf(fp,"%d %d\n",L1,M1); fprintf(fp,"%lf\n",dz);
return 0;
}

//-----------------------------------------------------------------------------------------------------------
// Установка значений диэлектрической проницаемости сред, откуда падает, и куда уходит свет.

BYTE clTask::SetEpsIncOut(clMaterialAll &Matter)
{
int num,ind; char c,*cp; complex *pEps;

if(Matter.IsOK()!=0) return 1;

// Установка диэлектрической проницаемости сред откуда падает, и куда выходит свет.
for(ind=0;ind<2;ind++) {
if(ind==0) { c=chEpsInc; cp="Medium of Incidence"; pEps=&EpsInc;}
else { c=chEpsOut; cp="Medium of Output"; pEps=&EpsOut;}
num=Matter.GetNum(c,wLength);
if(num<0||num>=Matter.N) { printf("Error: Symbol: '%c' in the %s does not correspond to a Substance !\n",c,cp); return 2;}
*pEps=Matter.Mat[num].eps;}

return 0;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Класс геометрии одного слоя фотонного кристалла.

//-----------------------------------------------------------------------------------------------------------
// Конструктор.

clLaySmplGeom::clLaySmplGeom(void)
{
Zero();
}

//-----------------------------------------------------------------------------------------------------------
// Деструктор.

clLaySmplGeom::~clLaySmplGeom(void)
{
Free();
}

//-----------------------------------------------------------------------------------------------------------
// Обнуление переменных.

void clLaySmplGeom::Zero(void)
{
Geom=NULL; laySize=0; Period=0.; GeomSymbs=NULL;
}

//-----------------------------------------------------------------------------------------------------------
// Освобождение памяти.

void clLaySmplGeom::Free(void)
{
SAFE_DELETE_ARR(Geom); SAFE_DELETE_ARR(GeomSymbs); laySize=0;
}

//-----------------------------------------------------------------------------------------------------------
// Выделение памяти для массива.

BYTE clLaySmplGeom::Alloc(int laySize_)
{
Free(); if(laySize_<=0) return 1; laySize=laySize_;
Geom=new int[laySize]; if(Geom==NULL) return 2;
GeomSymbs=new char[laySize+1]; if(GeomSymbs==NULL) return 3; GeomSymbs[laySize]='\0';
return 0;
} 

//-----------------------------------------------------------------------------------------------------------
// Проверки.

BYTE clLaySmplGeom::IsOK(void) const
{
if(laySize<=0) return 1; if(Geom==NULL) return 2; if(GeomSymbs==NULL) return 3; return 0;
}

//-----------------------------------------------------------------------------------------------------------
// Чтение геометрии образца из файла.

BYTE clLaySmplGeom::Read(FILE *fp,const clMaterialAll &Matter,double WaveLength,int numLay)
{
BYTE ReadString(FILE *file,char *pStr,int lnStr); // Чтение строки из текстового файла.

BYTE err; int i,num; char *str; 

if(fp==NULL) return 1; if(Matter.IsOK()!=0) return 2;
Free(); str=NULL; err=0;
str=new char[laySize_Max]; if(str==NULL) { err=3; goto end;} str[laySize_Max-1]='\0';
if(ReadString(fp,str,laySize_Max-1)!=0) { err=4; goto end;} // Чтение строки из текстового файла.
laySize=(int)strlen(str); if(laySize<=0) { err=5; goto end;}
if(Alloc(laySize)!=0) { err=6; goto end;}
for(i=0;i<laySize;i++) GeomSymbs[i]=str[i]; // Копирование.
for(i=0;i<laySize;i++) { num=Matter.GetNum(GeomSymbs[i],WaveLength); // Поиск совпадения в классе материалов кристалла.
if(num<0) {
printf("Error: Symbol: '%c' in the sample %d does not correspond to a Substance !\n",GeomSymbs[i],numLay);
err=7; goto end;}
Geom[i]=num;}

end: SAFE_DELETE_ARR(str); return err;
}

//-----------------------------------------------------------------------------------------------------------
// Проверка однородности образца.

BYTE clLaySmplGeom::IsHomogen(void)
{
BYTE fl; int i,nm,ng;

if(IsOK()!=0) return 0;
fl=1; for(i=0;i<laySize;i++) { ng=Geom[i]; if(i==0) nm=ng; else if(nm!=ng) { fl=0; break;}} return fl;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Структура одного образца по оси z фотонного кристалла.

//-----------------------------------------------------------------------------------------------------------
// Конструктор.

strLaySmpl::strLaySmpl(void)
{
Zero();
}

//-----------------------------------------------------------------------------------------------------------
// Деструктор.

strLaySmpl::~strLaySmpl(void)
{
Free();
}

//-----------------------------------------------------------------------------------------------------------
// Обнуление переменных.

void strLaySmpl::Zero(void)
{
numLay=0; Eps_Four=EpsInv_Four=NULL; N=L1=M1=0; N1=NF=0; M=0; flHomogen=0; EpsHomogen=Cmplx_1;
}

//-----------------------------------------------------------------------------------------------------------
// Освобождение памяти.

void strLaySmpl::Free(void)
{
SAFE_DELETE_ARR(Eps_Four); SAFE_DELETE_ARR(EpsInv_Four); N=L1=M1=0; N1=NF=0; M=0;
}

//-----------------------------------------------------------------------------------------------------------
// Проверки.

BYTE strLaySmpl::IsOK(void) const
{
if(LayGeom.IsOK()!=0) return 1; return 0;
}

//-----------------------------------------------------------------------------------------------------------
// Проверки массива Фурье преобразования.

BYTE strLaySmpl::IsOKFourier(void)
{
if(L1<=0) return 1; if(M1<=0) return 2; if(N1<=0) return 3; if(N<=0) return 4; if(NF<=0) return 5;
if(Eps_Four==NULL) return 6; if(EpsInv_Four==NULL) return 7; return 0;
}

//-----------------------------------------------------------------------------------------------------------
// Проверка блоков матриц 'Gamma'.

BYTE strLaySmpl::IsOKGamma(void) const
{
BYTE ki,kj,ind; const clMatrix *pM;

if(N<=0) return 1;
// Матрицы - блоки матриц 'Gamma' для 'P','P^-1'.
for(ind=0;ind<2;ind++) for(ki=0;ki<2;ki++) for(kj=0;kj<2;kj++) {
if(ind==0) pM=&GamP[ki][kj]; else pM=&GamPInv[ki][kj];
if(pM->IsOK()!=0) return 2; if(pM->Nx!=N||pM->Ny!=N) return 3;}
return 0;
}

//-----------------------------------------------------------------------------------------------------------
// Чтение данных об образце из файла.

BYTE strLaySmpl::Read(FILE *fp,const clCrystalSmpl &Samples)
{
if(fp==NULL) return 1; if(Samples.Matter.IsOK()!=0) return 2;
if(LayGeom.Read(fp,Samples.Matter,Samples.Task.wLength,numLay)!=0) return 3;
if(IsHomogen()!=0) return 4; return 0;
}

//-----------------------------------------------------------------------------------------------------------
// Установка флага однородности образца.

BYTE strLaySmpl::IsHomogen(void)
{
flHomogen=0; if(LayGeom.IsOK()!=0) return 1; flHomogen=LayGeom.IsHomogen(); return 0;
}

//-----------------------------------------------------------------------------------------------------------
// Быстрое преобразование Фурье для образца.

BYTE strLaySmpl::FastFTLay(int L1_,int M1_,const clMaterialAll &Matter)
{
BYTE GetEpsDistr(clRowMatr &Distr,const clLaySmplGeom &LayGeom,const clMaterialAll &Matter,BYTE flInv); // Составление массива диэлектрической проницаемости или её обратной величины из массива геометрии слоя.
BYTE GetN(int M,int *N); // Нахождение 'N' - 2 в степени 'M'.
BYTE FastFT(complex *A,complex *B,BYTE M,SCHAR dir); // Быстрое преобразование Фурье.
double Log2(double n); // Логарифм по основанию 2.

BYTE err; int i,iv,j; complex *Eps_Four_,*EpsInv_Four_; clRowMatr Distr;

if(flHomogen!=0) return 1; Free();
if(L1_<=0) return 2; if(M1_<=0) return 3; if(IsOK()!=0) return 4; if(Matter.IsOK()!=0) return 5;
L1=L1_; M1=M1_; N=L1+M1+1; // Установка числа гармоник.
N1=(L1+M1)*2+1; M=(int)Log2(N1)+1; if(GetN(M,&NF)!=0) return 6; if(NF<N1) return 7;
if(NF<LayGeom.laySize) { printf("\nNumber of Fourier Harmonics is Too Small !\n"); return 8;}

Eps_Four_=EpsInv_Four_=NULL; err=0;
// Выделение памяти для вспомогательных массивов.
Eps_Four_=new complex[NF]; if(Eps_Four_==NULL) { err=9; goto end;}
EpsInv_Four_=new complex[NF]; if(EpsInv_Four_==NULL) { err=10; goto end;}
Eps_Four=new complex[N1]; if(Eps_Four==NULL) { err=11; goto end;}
EpsInv_Four=new complex[N1]; if(EpsInv_Four==NULL) { err=12; goto end;}

// Выделение памяти для массивов распределения диэлектрической проницаемости слоя и её обратной величины по оси 'X'.
if(Distr.Alloc(NF)!=0) { err=13; goto end;}

// Заполнение массива диэлектрической проницаемости из массива геометрии слоя.
if(GetEpsDistr(Distr,LayGeom,Matter,0)!=0) { err=14; goto end;}
if(FastFT(Distr.Vect,Eps_Four_,M,1)!=0) { err=15; goto end;} // FFT для распределения 'eps'.

// Заполнение массива обратных значений диэлектрической проницаемости из массива геометрии слоя.
if(GetEpsDistr(Distr,LayGeom,Matter,1)!=0) { err=16; goto end;}
if(FastFT(Distr.Vect,EpsInv_Four_,M,1)!=0) { err=17; goto end;} // FFT для распределения '1/eps'.

// Урезание Фурье ряда. В массивах разложения Фурье 'eps', '1/eps' сначала записаны положительные гармоники, затем начиная с 'L1+M1', отрицательные.
for(i=0;i<=L1+M1;i++) { Eps_Four[i]=Eps_Four_[i]; EpsInv_Four[i]=EpsInv_Four_[i];} // Заполнение положительных гармоник.
for(i=0;i<L1+M1;i++) { iv=L1+M1+1+i; j=NF-1-i; Eps_Four[iv]=Eps_Four_[j]; EpsInv_Four[iv]=EpsInv_Four_[j];} // Заполнение отрицательных гармоник.

end: SAFE_DELETE_ARR(Eps_Four_); SAFE_DELETE_ARR(EpsInv_Four_); return err;
}

//-----------------------------------------------------------------------------------------------------------
// Установка параметров для однородного образца.

BYTE strLaySmpl::SetParHomogen(int L1_,int M1_,const clMaterialAll &Matter)
{
BYTE GetN(int M,int *N); // Нахождение 'N' - 2 в степени 'M'.
double Log2(double n); // Логарифм по основанию 2.

int nMat;

if(flHomogen==0) return 1;
if(L1_<=0) return 2; if(M1_<=0) return 3; if(IsOK()!=0) return 4; if(Matter.IsOK()!=0) return 5;
L1=L1_; M1=M1_; N=L1+M1+1; // Установка числа гармоник.
N1=(L1+M1)*2+1; M=(int)Log2(N1)+1; if(GetN(M,&NF)!=0) return 6; if(NF<N1) return 7;
nMat=LayGeom.Geom[0]; EpsHomogen=Matter.Mat[nMat].eps; // Установка диэлектрической проницаемости однородного слоя.
return 0;
}


//-----------------------------------------------------------------------------------------------------------
// Вычисление матрицы Mk для расчёта Фурье коэффициентов компонент поля.

BYTE strLaySmpl::CompMk(const clTask &Task)
{
BYTE InvMatr(const clMatrix &M,clMatrix &Inv); // Нахождение обратной матрицы методом исключения Гаусса с выбором главного элемента.

int i,j,iv,jv,jSh,jShAk,jShV,ND; double dK,k0X,k0,v; complex kz,eps,c; clMatrix V;
clRowMatr Kx,Kx2;

if(flHomogen!=0) return 1; if(N<=0) return 2; if(IsOKFourier()!=0) return 3; ND=N*2;
if(Task.wLength<=0.) return 4; if(LayGeom.Period<=0.) return 5;

// Расчёт волнового вектора и его компоненты по оси X.
k0=2.*M_PI/Task.wLength; k0X=k0*sin(Task.Theta*M_PI/180.); dK=2.*M_PI/LayGeom.Period;

// Выделение памяти для матриц операторов уравнения Максвелла после преобразования Фурье.
if(Mk.Alloc(ND,ND)!=0) return 6;

// Выделение памяти для вспомогательных матриц.
if(Ek.Alloc(N,N)!=0) return 7; if(Ak.Alloc(N,N)!=0) return 8;
if(Kx.Alloc(N)!=0) return 9; if(Kx2.Alloc(N)!=0) return 10;
if(EkInv.Alloc(N,N)!=0) return 11;
if(AkInv.Alloc(N,N)!=0) return 12;

// Заполнение матриц Ek и Ak.
for(i=0;i<N;i++) { jSh=i*N;
for(j=0;j<N;j++) { jv=j+jSh; if(i>=j) iv=i-j; else iv=j-i+L1+M1; if(iv<0||iv>=N1) return 13;
Ek.Matr[jv]=Eps_Four[iv]; Ak.Matr[jv]=EpsInv_Four[iv];}}

// Вычисление матриц Kx,Kx2.
for(i=0;i<N;i++) { v=((double)(-L1+i)*dK+k0X)/k0; Kx.Vect[i]=complex(v,0.); Kx2.Vect[i]=complex(v*v,0.);}

// Заполнение матрицы Mk в зависимости от типа поляризации волны.
if((Mk<<Cmplx_0)!=0) return 14; // Начальное заполнение нулями.

switch(Task.flPol) { default: return 15;
// Вектор E в плоскости падения.
case p_wave_Pol:
if(InvMatr(Ek,EkInv)!=0) return 16;
V=Kx%EkInv%Kx; V-=Cmplx_1; if(V.IsOK()!=0) return 17; if(V.Nx!=N||V.Ny!=N) return 18;
if(InvMatr(Ak,AkInv)!=0) return 19;
for(i=0;i<N;i++) { jSh=i*ND; jShV=i*N; for(j=N;j<ND;j++) Mk.Matr[j+jSh]=V.Matr[j-N+jShV];}
for(i=N;i<ND;i++) { jSh=i*ND; jShAk=(i-N)*N; for(j=0;j<N;j++) Mk.Matr[j+jSh]=AkInv.Matr[j+jShAk];}
break;

// Вектор E перпендикулярен плоскости падения.
case s_wave_Pol:
V=Kx2-Ek; if(V.IsOK()!=0) return 20; if(V.Nx!=N||V.Ny!=N) return 21;
for(i=0;i<N;i++) { jSh=i*ND; for(j=0;j<N;j++) if(i==j) Mk.Matr[j+N+jSh]=Cmplx_1;}
for(i=N;i<ND;i++) { jSh=i*ND; jShAk=(i-N)*N; for(j=0;j<N;j++) Mk.Matr[j+jSh]=V.Matr[j+jShAk];}
break;
}
Mk*=k0; // Множитель для матрицы 'Mk'.

return 0;
}

//-----------------------------------------------------------------------------------------------------------
// Нахождение собственных значений и векторов матрицы Mk.

BYTE strLaySmpl::FindEigVectMk(void)
{
BYTE FindEigenValVect(int N,const clMatrix &A,clRowMatr &EigLam,clMatrix &P,clMatrix &PInv,BYTE flInv); // Нахождение собственных значений и векторов матрицы 'A'. Метод неявных QR итераций с приведением матрицы к виду Хессенберга отражениями Хаусхолдера.
BYTE InvMatr(const clMatrix &M,clMatrix &Inv); // Нахождение обратной матрицы методом исключения Гаусса с выбором главного элемента.

BYTE fl,err; int ND;

if(flHomogen!=0) return 1; if(N<=0) return 2; ND=N*2;
if(Mk.IsOK()!=0) return 3; if(Mk.Nx!=ND||Mk.Ny!=ND) return 4;
EigLamMk.Free(); EigVectMk.Free(); EigVectMkInv.Free();

// Выделение памяти для собственных значений и собственных векторов матрицы 'Mk'.
if(EigLamMk.Alloc(ND)!=0) return 5;
if(EigVectMk.Alloc(ND,ND)!=0) return 6;
if(EigVectMkInv.Alloc(ND,ND)!=0) return 7;

fl=flFindPInv;
err=FindEigenValVect(ND,Mk,EigLamMk,EigVectMk,EigVectMkInv,fl); // Нахождение собственных значений и векторов матрицы 'Mk'.
if(err!=0) { printf("FindEigenValVect err: %d\n",err); return 8;}
if(fl==0) { if(InvMatr(EigVectMk,EigVectMkInv)!=0) return 9;} // Нахождение обратной матрицы 'P^-1'.

if(FixFormEig()!=0) return 10; // Приведение столбца собственных значений к нужному в данной задаче виде (по возрастанию действительной части).
if(FindGamma()!=0) return 11; // Нахождение блоков матриц 'Gamma'.
return 0;
}

//-----------------------------------------------------------------------------------------------------------
// Нахождение собственных значений и векторов матрицы Mk (случай однородного слоя).

BYTE strLaySmpl::FindEigVectMkHomogen(const clTask &Task)
{
void MessageSmallEps_(int numLay); // Сообщение о близком к нулю значении диэлектрической проницаемости.
void MessageSingLam_(int numLay); // Сообщение о близком к нулю собственном значении.

int i,j,jv,jSh,jShv,ND; double k0,k0X,dK,kx,k2,s,v; complex lam,c; clRowMatr D,DInv; clMatrix *pM;

if(flHomogen==0) return 1; if(Task.wLength<=0.) return 2; if(LayGeom.Period<=0.) return 3;
if(N<=0) return 4; ND=N*2;
EigLamMk.Free(); EigVectMk.Free(); EigVectMkInv.Free();

dK=2.*M_PI/LayGeom.Period; k0=2.*M_PI/Task.wLength; k2=k0*k0; k0X=k0*sin(Task.Theta*M_PI/180.);

// Выделение памяти для собственных значений и собственных векторов матрицы 'Mk'.
if(EigLamMk.Alloc(ND)!=0) return 5;
if(EigVectMk.Alloc(ND,ND)!=0) return 6;
if(EigVectMkInv.Alloc(ND,ND)!=0) return 7;

// Нахождение собственных значений.
for(i=0;i<N;i++) { kx=k0X+(double)(i-L1)*dK; lam=Cmplx_I*sqrt(complex(k2,0.)*EpsHomogen-complex(kx*kx,0.));
EigLamMk.Vect[i]=lam; EigLamMk.Vect[i+N]=-lam;}

EigVectMk=Cmplx_0; EigVectMkInv=Cmplx_0; // Начальное заполнение нулями.

// Нахождение матриц собственных векторов.
switch(Task.flPol) { default: return 8;
// Вектор E в плоскости падения.
case p_wave_Pol:
if(abs(EpsHomogen)<smValEps) { MessageSmallEps_(numLay); return messSmallEps;}
for(i=0;i<N;i++) { lam=EigLamMk.Vect[i]; j=i; jv=i+N; jSh=i*ND; jShv=(i+N)*ND;
if(abs(lam)<smValLam) { MessageSingLam_(numLay); return messSingLam;}
c=lam/(k0*EpsHomogen); s=1.+real(c*conj(c)); if(s<=0.) return 9; s=1./sqrt(s);
pM=&EigVectMk; // Матрица 'P'.
pM->Matr[j+jSh]=c*s; pM->Matr[j+jShv]=complex(s,0.);
pM->Matr[jv+jSh]=-c*s; pM->Matr[jv+jShv]=complex(s,0.);
pM=&EigVectMkInv; // Обратная матрица 'P^-1'.
c=0.5/(c*s); pM->Matr[j+jSh]=c; pM->Matr[j+jShv]=-c;
v=0.5/s; pM->Matr[jv+jSh]=pM->Matr[jv+jShv]=complex(v,0.);}
break;

// Вектор E перпендикулярен плоскости падения.
case s_wave_Pol:
for(i=0;i<N;i++) { lam=EigLamMk.Vect[i]; j=i; jv=i+N; jSh=i*ND; jShv=(i+N)*ND;
if(abs(lam)<smValLam) { MessageSingLam_(numLay); return messSingLam;}
c=lam/k0; s=1.+real(c*conj(c)); if(s<=0.) return 10; s=1./sqrt(s);
pM=&EigVectMk; // Матрица 'P'.
pM->Matr[j+jSh]=complex(s,0.); pM->Matr[j+jShv]=c*s;
pM->Matr[jv+jSh]=complex(s,0.); pM->Matr[jv+jShv]=-c*s;
pM=&EigVectMkInv; // Обратная матрица 'P^-1'.
v=0.5/s; pM->Matr[j+jSh]=pM->Matr[j+jShv]=complex(v,0.);
c=0.5/(c*s); pM->Matr[jv+jSh]=c; pM->Matr[jv+jShv]=-c;}
break;
}

if(FixFormEig()!=0) return 11; // Приведение столбца собственных значений к нужному в данной задаче виде (по возрастанию действительной части).
if(FindGamma()!=0) return 12; // Нахождение блоков матриц 'Gamma'.

return 0;
}

//-----------------------------------------------------------------------------------------------------------
// Приведение столбца собственных значений к нужному в данной задаче виде (по возрастанию действительной части).

BYTE strLaySmpl::FixFormEig(void)
{
BYTE SortEigVal(int N,clRowMatr Arr,int *Num); // Сортировка комплексного массива собственных значений методом пузырька с использованием массива подстановок.

BYTE err; int i,j,jSh,jShv,ND,*Num,*NumV; clRowMatr EigLam; clMatrix EigVect,EigVectInv,*pM;

if(N<=0) return 1; ND=N*2;
if(EigLamMk.IsOK()!=0) return 2; if(EigLamMk.N!=ND) return 3;
pM=&EigVectMk; if(pM->IsOK()!=0) return 4; if(pM->Nx!=ND||pM->Ny!=ND) return 5;
pM=&EigVectMkInv; if(pM->IsOK()!=0) return 6; if(pM->Nx!=ND||pM->Ny!=ND) return 7;

Num=NumV=NULL; err=0;

// Выделение памяти.
Num=new int[ND]; if(Num==NULL) { err=8; goto end;}
NumV=new int[ND]; if(NumV==NULL) { err=9; goto end;}

EigLam=EigLamMk; EigVect=EigVectMk; EigVectInv=EigVectMkInv; // Инициализация столбца 'EigLam' и матриц 'EigVect','EigVectInv'.
for(i=0;i<ND;i++) Num[i]=ND-1-i; // Инициализация массива подстановок.
if(SortEigVal(ND,EigLam,Num)!=0) { err=10; goto end;} // Сортировка по возрастанию действительной части.

// Заполнение итоговых столбца собственных значений и матриц собственных векторов.
for(i=0;i<ND;i++) EigLamMk.Vect[i]=EigLam.Vect[Num[i]];
for(i=0;i<ND;i++) { jSh=i*ND; for(j=0;j<ND;j++) EigVectMk.Matr[j+jSh]=EigVect.Matr[Num[j]+jSh];} // Сортировка матрицы P.
for(i=0;i<ND;i++) { jSh=i*ND; jShv=Num[i]*ND; for(j=0;j<ND;j++) EigVectMkInv.Matr[j+jSh]=EigVectInv.Matr[j+jShv];} // Сортировка матрицы P^-1.

end: SAFE_DELETE_ARR(NumV); SAFE_DELETE_ARR(Num); return err;
}

//-----------------------------------------------------------------------------------------------------------
// Нахождение блоков матриц 'Gamma'.

BYTE strLaySmpl::FindGamma(void)
{
BYTE ki,kj,ind; int i,j,jB,jShB,jSh0,jSh1,ND; clRowMatr DiagExp; clMatrix *pM,*pMG,(*pGam)[2][2];

if(N<=0) return 1; ND=N*2;
if(EigLamMk.IsOK()!=0) return 2; if(EigLamMk.N!=ND) return 3;
pM=&EigVectMk; if(pM->IsOK()!=0) return 4; if(pM->Nx!=ND||pM->Ny!=ND) return 5;
pM=&EigVectMkInv; if(pM->IsOK()!=0) return 6; if(pM->Nx!=ND||pM->Ny!=ND) return 7;

// Матрицы - блоки матриц 'Gamma' для 'P','P^-1'.
for(ind=0;ind<2;ind++) for(ki=0;ki<2;ki++) for(kj=0;kj<2;kj++) {
if(ind==0) pM=&GamP[ki][kj]; else pM=&GamPInv[ki][kj]; if(pM->Alloc(N,N)!=0) return 10;}

// Нахождение четырех блоков матриц 'Gamma' для матриц 'P','P^-1'.
for(ind=0;ind<2;ind++) { if(ind==0) { pMG=&EigVectMk; pGam=&GamP;} else { pMG=&EigVectMkInv; pGam=&GamPInv;}
for(i=0;i<N;i++) { jShB=i*N; jSh0=i*ND; jSh1=(i+N)*ND; 
for(j=0;j<N;j++) { jB=j+jShB;
(*pGam)[0][0].Matr[jB]=pMG->Matr[j+jSh0]; (*pGam)[0][1].Matr[jB]=pMG->Matr[j+N+jSh0];
(*pGam)[1][0].Matr[jB]=pMG->Matr[j+jSh1]; (*pGam)[1][1].Matr[jB]=pMG->Matr[j+N+jSh1];
}} // Конец циклов по 'i','j' - заполнение матриц.
} // Конец цикла по 'ind' - тип матрицы.
return 0;
}

//-----------------------------------------------------------------------------------------------------------
// Составление массива диэлектрической проницаемости или её обратной величины из массива геометрии слоя.

BYTE GetEpsDistr(clRowMatr &Distr,const clLaySmplGeom &LayGeom,const clMaterialAll &Matter,BYTE flInv)
{
BYTE fl[2]; int i,k,nMat[3],nM,nMv,ls,n; double st,x,wei; complex v,Val;

if(Distr.IsOK()!=0) return 1; if(LayGeom.IsOK()!=0) return 2; if(Distr.N<LayGeom.laySize) return 3;
if(Matter.IsOK()!=0) return 4;
ls=LayGeom.laySize; st=(double)ls/(double)Distr.N;
// Цикл по отрезкам разбиения по числу N для последующего Фурье преобразования.
for(i=0;i<Distr.N;i++) { x=((double)i+0.5)*st;
for(k=0;k<3;k++) { n=(int)(x+st*0.5*(double)(k-1));
if(n<0) n+=ls; if(n>=ls) n-=ls; nMat[k]=LayGeom.Geom[n];}
for(k=0;k<2;k++) if(nMat[k]==nMat[k+1]) fl[k]=0; else fl[k]=1;
if(fl[0]!=0&&fl[1]!=0) return 5; // На начальной и конечной границах отрезка вещество отличается от вещества в центре.
// Вещество одинаково на всём отрезке.
nM=nMat[1]; // Номер вещества в центре отрезка.
if(fl[0]==0&&fl[1]==0) { v=Matter.Mat[nM].eps; if(flInv==0) Val=v; else Val=Inv(v);}
// Вещество неодинаково на протяжении отрезка.
else {
if(fl[0]!=0) { nMv=nMat[0]; wei=(st*0.5-(x-(double)(int)x))/st;}
else { nMv=nMat[2]; wei=(st*0.5-((double)((int)x+1)-x))/st;}
v=Matter.Mat[nM].eps; if(flInv!=0) v=Inv(v); Val=v*(1.-wei);
v=Matter.Mat[nMv].eps; if(flInv!=0) v=Inv(v); Val+=v*wei;}
Distr.Vect[i]=Val;
} // Конец цикла по отрезкам разбиения 'i'.
return 0;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Класс образцов слоев фотонного кристалла.

//-----------------------------------------------------------------------------------------------------------
// Конструктор.

clCrystalSmpl::clCrystalSmpl(void)
{
Zero();
}

//-----------------------------------------------------------------------------------------------------------
// Деструктор.

clCrystalSmpl::~clCrystalSmpl(void)
{
Free();
}

//-----------------------------------------------------------------------------------------------------------
// Обнуление переменных.

void clCrystalSmpl::Zero(void)
{
Lays=NULL; laysNum=0; Period=0.;
}

//-----------------------------------------------------------------------------------------------------------
// Освобождение памяти.

void clCrystalSmpl::Free(void)
{
SAFE_DELETE_ARR(Lays); laysNum=0;
}

//-----------------------------------------------------------------------------------------------------------
// Выделение памяти для массива.

BYTE clCrystalSmpl::Alloc(int laysNum_)
{
Free(); if(laysNum_<=0) return 1; laysNum=laysNum_;
Lays=new strLaySmpl[laysNum]; if(Lays==NULL) return 2; return 0;
} 

//-----------------------------------------------------------------------------------------------------------
// Проверки.

BYTE clCrystalSmpl::IsOK(void) const
{
int i;

if(laysNum<=0) return 1; if(Lays==NULL) return 2; 
for(i=0;i<laysNum;i++) if(Lays[i].IsOK()!=0) return 3; return 0;
}

//-----------------------------------------------------------------------------------------------------------
// Проверка блоков матриц 'Gamma' во всех образцах.

BYTE clCrystalSmpl::IsOKGamma(void)
{
int i,N;

if(IsOK()!=0) return 1; N=Lays[0].N; if(N<=0) return 2;
for(i=0;i<laysNum;i++) { if(Lays[i].IsOKGamma()!=0) return 3; if(Lays[i].N!=N) return 4;} return 0;
}

//-----------------------------------------------------------------------------------------------------------
// Чтение данных об образцах из файла.

BYTE clCrystalSmpl::ReadInpData(void)
{
BYTE err; int i,LaysNum_; char *fname; FILE *fp; struct strLaySmpl *pL;

// Чтение списка материалов.
if(Matter.Read(Task.FName_Mat)!=0) { printf("Error in reading of file of materials !\n"); _getch(); return 1;}
if(Matter.SetEps(Task.wLength)!=0) return 2;
if(Task.SetEpsIncOut(Matter)!=0) return 3;

fp=NULL; fname=NULL; err=0;
fname=GetFName(Task.FName_Cryst,extCrystSamples); if(fname==NULL) { err=4; goto end;}
//ccc
// Послойное чтение информации об образцах.
fp=fopen(fname,"r"); if(fp==NULL) { printf("Can not open file: %s !\n",fname); err=5; goto end;}
if(fscanf(fp,"%d",&LaysNum_)==EOF) { err=6; goto end;}
if(Alloc(LaysNum_)!=0) { err=7; goto end;}
if(fscanf(fp,"%lf",&Period)==EOF) { err=8; goto end;}
if(fscanf(fp,"\n")==EOF) { err=9; goto end;}
for(i=0;i<laysNum;i++) { pL=Lays+i; pL->LayGeom.Period=Period; pL->numLay=i;
if(fscanf(fp,"\n")==EOF) { err=10; goto end;} 
if(pL->Read(fp,*this)!=0) { err=11; goto end;}}

end: SAFE_CLOSE(fp); SAFE_DELETE_ARR(fname); return err;
}

//-----------------------------------------------------------------------------------------------------------
// Быстрое преобразование Фурье для всех образцов.

BYTE clCrystalSmpl::FastFTAll(void)
{
int i; struct strLaySmpl *pL;

if(Task.L1<=0) return 1; if(Task.M1<=0) return 2; if(IsOK()!=0) return 3;
for(i=0;i<laysNum;i++) { pL=Lays+i; if(Task.flTypComp==0) printf("Computing of the Sample %d",i);
if(pL->flHomogen==0) { if(pL->FastFTLay(Task.L1,Task.M1,Matter)!=0) return 4;}
else { if(pL->SetParHomogen(Task.L1,Task.M1,Matter)!=0) return 5;}
if(Task.flTypComp==0) printf(" - Completed\n");}
return 0;
} 

//-----------------------------------------------------------------------------------------------------------
// Вычисление матрицы 'Mk' для всех образцов слоёв.

BYTE clCrystalSmpl::CompMk(void)
{
BYTE flParallel,numHomogen,err,flBreak; int i;

if(IsOK()!=0) return 1;

// Проверка условия распараллеливания.
flParallel=1;
if(laysNum<=1) flParallel=0;
numHomogen=0; for(i=0;i<laysNum;i++) if(Lays[i].flHomogen==1) numHomogen++; // Подсчет числа однородных слоев.
if(numHomogen>=laysNum-1) flParallel=0;

err=0;
#pragma omp parallel if(flParallel==1) default(shared) private(flBreak) // Параллельная область.
{ 
struct strLaySmpl *pL; flBreak=0;

#pragma omp for schedule(static) private(i)
for(i=0;i<laysNum;i++) { if(flBreak==1) continue;
pL=Lays+i; if(pL->flHomogen!=0) continue; 
if(Task.flTypComp==0) if(flParallel==0) printf("Computing of the Sample %d",i);
if(pL->CompMk(Task)!=0) { flBreak=1; err=2; continue; } 
if(Task.flTypComp==0) if(flParallel==0) printf(" - Completed\n");}
} // Конец параллельной области.
return err;
}

//-----------------------------------------------------------------------------------------------------------
// Нахождение собственных значений матриц 'Mk' для всех слоёв.

BYTE clCrystalSmpl::FindEigVectMk(void)
{
BYTE err,flParallel,flBreak,numHomogen; int i;

if(IsOK()!=0) return 1;

// Проверка условия распараллеливания.
flParallel=1; if(laysNum<=1) flParallel=0;
numHomogen=0; for(i=0;i<laysNum;i++) if(Lays[i].flHomogen==1) numHomogen++; // Подсчет числа однородных слоев.
if(numHomogen>=laysNum-1) flParallel=0;

err=0;
#pragma omp parallel if(flParallel==1) default(shared) private(flBreak) // Параллельная область.
{ 
struct strLaySmpl *pL; flBreak=0;
#pragma omp for schedule(static) private(i)
for(i=0;i<laysNum;i++) { if(flBreak==1) continue;
pL=Lays+i; if(Task.flTypComp==0) if(flParallel==0) printf("Computing of the Sample %d",i);
if(pL->flHomogen==0) { if(pL->FindEigVectMk()!=0) { flBreak=1; err=2; continue;}}
else { if(pL->FindEigVectMkHomogen(Task)!=0) { flBreak=1; err=3; continue;}}
if(Task.flTypComp==0) if(flParallel==0) printf(" - Completed\n");}
}
return err;
}

//-----------------------------------------------------------------------------------------------------------
// Получение указателя на образец слоя с номером 'num'.

strLaySmpl *clCrystalSmpl::Get(int num)
{
if(IsOK()!=0) return NULL; if(num<0||num>=laysNum) return NULL; return Lays+num;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Структура одного слоя по оси z фотонного кристалла.

//-----------------------------------------------------------------------------------------------------------
// Конструктор.

strLay::strLay(void)
{
Zero();
}

//-----------------------------------------------------------------------------------------------------------
// Обнуление переменных.

void strLay::Zero(void)
{
numLay=numSmpl=-1; depth=dz=0.; N=0; nDiv=0;
}

//-----------------------------------------------------------------------------------------------------------
// Проверки.

BYTE strLay::IsOK(void)
{
if(numLay<0) return 1; if(numSmpl<0) return 2; if(depth<=0.) return 3; return 0;
}

//-----------------------------------------------------------------------------------------------------------
// Вычисление матриц 'DExp00Inv', 'DExp11'.

BYTE strLay::FindGammaExp(const strLaySmpl &Smpl)
{
int i;

if(IsOK()!=0) return 1; if(Smpl.IsOK()!=0) return 2; N=Smpl.N;
if(Smpl.EigLamMk.IsOK()!=0) return 3; if(Smpl.EigLamMk.N!=N*2) return 4;

// Выделение памяти.
if(DExp00Inv.Alloc(N)!=0) return 5; // Диагональная матрица, обратная к 'Gamma00'.
if(DExp11.Alloc(N)!=0) return 6; // Диагональная матрица 'Gamma11'.
// Нахождение диагональных матриц, составленных из экспонент.
for(i=0;i<N;i++) {
DExp00Inv.Vect[i]=exp(Smpl.EigLamMk.Vect[i]*depth);
DExp11.Vect[i]=exp(-Smpl.EigLamMk.Vect[i+N]*depth);}
return 0;
}

//-----------------------------------------------------------------------------------------------------------
// Проверка блоков матриц 'Gamma'.

BYTE strLay::IsOKGammaExp(void)
{
if(N<=0) return 1;
if(DExp00Inv.IsOK()!=0) return 2; if(DExp00Inv.N!=N) return 3; // Диагональная матрица, обратная к 'Gamma00'.
if(DExp11.IsOK()!=0) return 4; if(DExp11.N!=N) return 5; // Диагональная матрица 'Gamma11'.
return 0;
}

//-----------------------------------------------------------------------------------------------------------
// Нахождение распределения поля внутри слоя кристалла.

BYTE strLay::FindFullFieldDistr(clCrystalSmpl &Samples)
{
struct strLaySmpl *pLS;

if(CompDz(Samples.Task.dz)!=0) return 1; // Расчёт 'dz'.
pLS=Samples.Get(numSmpl); if(pLS==NULL) return 2;
if(FindDistrFourEH(*pLS)!=0) return 3; // Нахождение коэффициентов Фурье разложения полей Е и H.
if(FindDistrEH(*pLS,Samples.Task)!=0) return 4; // Нахождение всего распределения поля.
return 0;
}

//-----------------------------------------------------------------------------------------------------------
// Нахождение 'dz' из начального, приближённого 'dz' и выделение памяти для массивов столбцов коэффициентов Фурье полей E и H.

BYTE strLay::CompDz(double dz_)
{
if(dz_<=0.) return 1; nDiv=int(depth/dz_+0.5); if(nDiv<1) nDiv=1; dz=depth/(double)nDiv; // Определение 'dz'.
if(EDistrFour.Alloc(nDiv+1,N)!=0) return 2;
if(HDistrFour.Alloc(nDiv+1,N)!=0) return 3;
return 0;
}

//-----------------------------------------------------------------------------------------------------------
// Нахождение коэффициентов Фурье разложения полей Е и H в слое.

BYTE strLay::FindDistrFourEH(const strLaySmpl &Smpl)
{
BYTE ki,kj,ind; int i,j; double zX,zY; clRowMatr X,Y,*pRM; const clMatrix *pM;

if(nDiv<1) return 1; if(dz<=0.) return 2; if(IsOK()!=0) return 3; if(IsOKGammaExp()!=0) return 4;
if(Smpl.EigLamMk.IsOK()!=0) return 5; if(Smpl.EigLamMk.N!=N*2) return 6;

// Проверка блоков матрицы 'P' для образца слоя.
for(ki=0;ki<2;ki++) for(kj=0;kj<2;kj++) { pM=&Smpl.GamP[ki][kj];
if(pM->IsOK()!=0) return 7; if(pM->Nx!=N||pM->Ny!=N) return 8;}

// Проверка матриц 'EFour', 'HFour' для слоя.
for(ind=0;ind<6;ind++) {
switch(ind) { default: continue;
case 0: pRM=&EFour[0]; break;
case 1: pRM=&HFour[0]; break;
case 2: pRM=&EFour[3]; break;
case 3: pRM=&HFour[3]; break;
case 4: pRM=&EFour[1]; break;
case 5: pRM=&HFour[2]; break;}
if(pRM->IsOK()!=0) return 9; if(pRM->N!=N) return 10;}

// Копирование уже найденных столбцов Фурье разложения полей E,H на границах слоя.
pRM=EDistrFour.Get(0); if(pRM==NULL) return 11; *pRM=EFour[0];
pRM=HDistrFour.Get(0); if(pRM==NULL) return 12; *pRM=HFour[0];
pRM=EDistrFour.Get(nDiv); if(pRM==NULL) return 13; *pRM=EFour[3];
pRM=HDistrFour.Get(nDiv); if(pRM==NULL) return 14; *pRM=HFour[3];

// Вычисляем столбец 'X' от начала слоя к концу, а столбец 'Y' с конца слоя к началу, чтобы избежать роста амплитуд собственных волн.
for(j=1;j<nDiv;j++) { zX=dz*(double)j; zY=dz*(double)(nDiv-j);
// Нахождение диагональных матриц, составленных из экспонент.
for(i=0;i<N;i++) {
DExp00Inv.Vect[i]=exp(Smpl.EigLamMk.Vect[i]*zX);
DExp11.Vect[i]=exp(-Smpl.EigLamMk.Vect[i+N]*zY);}

// Находим коэффициенты Фурье разложения полей Е и H.
X=DExp00Inv*EFour[1]; Y=DExp11*HFour[2];
pRM=EDistrFour.Get(j); if(pRM==NULL) return 15; *pRM=Smpl.GamP[0][0]*X+Smpl.GamP[0][1]*Y;
pRM=HDistrFour.Get(j); if(pRM==NULL) return 16; *pRM=Smpl.GamP[1][0]*X+Smpl.GamP[1][1]*Y;
} // Конец цикла по 'j'.

return 0;
}

#define flLiFactEz 1
#define flCompEzConv 0
//-----------------------------------------------------------------------------------------------------------
// Нахождение распределения полей Е и H в слое (выполнение обратного преобразования Фурье).

BYTE strLay::FindDistrEH(const strLaySmpl &Smpl,clTask &Task)
{
BYTE GetN(int M,int *N); // Нахождение 'N' - 2 в степени 'M'.
BYTE GetWeightWinCos(float *Weight,int nPoints); // Нахождение функции для косинусоидального сглаживания.
BYTE GetWeightKaisBess(float *Weight,int nPoints); // Окно Кайзера-Бесселя.
BYTE FillArrInvFour(int NF,complex *A,const clRowMatr &Row,int L1,int M1); // Заполнение массива 'A' для обратного преобразования Фурье.
BYTE SmoothSpectrum(int NF,complex *A,int NP,float *WeiP,int nWeiP,int NN,float *WeiN,int nWeiN); // Сглаживание спектра.
BYTE FastFT(complex *A,complex *B,BYTE M,SCHAR dir); // Быстрое преобразование Фурье.
void MessageSmallEps_(int numLay); // Сообщение о близком к нулю значении диэлектрической проницаемости.

BYTE ind,flSm,err,flBreak; int NF,j,i,M,L1,M1,nWeiP,nWeiN; float *WeiP,*WeiN; double dK,k0,k0X,kx;
clRows *pRW; clRowMatr ZFourV,ZFour; funWinSmooth pFunSm;

if(N<=0) return 1; if(nDiv<=0) return 2; if(Task.wLength<=0.) return 3;
if(Smpl.LayGeom.Period<=0.) return 4; if(Smpl.M<=0||Smpl.NF<=0) return 5;
// Проверка массивов гармоник Фурье полей.
for(ind=0;ind<2;ind++) { if(ind==0) pRW=&EDistrFour; else pRW=&HDistrFour;
if(pRW->IsOK_All()!=0) return 6; if(pRW->nRows!=nDiv+1) return 7; if(pRW->N!=N) return 8;}

// Установка числа точек по 'X' для распределения полей (числа 'N'=2^M для быстрого преобразования Фурье).
M=Smpl.M; if(GetN(M,&NF)!=0) return 9; if(NF!=Smpl.NF) return 10;
L1=Smpl.L1; M1=Smpl.M1; if(N!=L1+M1+1) return 11;
if(Task.flPol==p_wave_Pol&&Smpl.flHomogen==0) { if(Smpl.EpsInv_Four==NULL) return 12;}

// Вспомогательные величины.
dK=2.*M_PI/Smpl.LayGeom.Period; k0=2.*M_PI/Task.wLength; k0X=k0*sin(Task.Theta*M_PI/180.);

WeiP=WeiN=NULL; err=0;

// Выделение памяти для функций сглаживания пространственного спектра.
flSm=flSmooth; if(flSm!=0) {
switch(typFunSmooth) { default: err=14; goto end;
// Косинусоидальное сглаживание.
case typFunSmoothCos:
nWeiP=(int)((float)M1*CoeWinSmooth); if(nWeiP>M1) nWeiP=M1;
nWeiN=(int)((float)L1*CoeWinSmooth); if(nWeiN>L1) nWeiN=L1;
pFunSm=GetWeightWinCos; break;
// Окно Кайзера-Бесселя.
case typFunKaisBess: nWeiP=M1; nWeiN=L1; pFunSm=GetWeightKaisBess; break;
}
if(nWeiP>0) { WeiP=new float[nWeiP]; if(WeiP==NULL) { err=15; goto end;}}
if(nWeiN>0) { WeiN=new float[nWeiN]; if(WeiN==NULL) { err=16; goto end;}}
}

// Выделение памяти для массивов распределений.
if(EDistr.Alloc(nDiv+1,NF)!=0) { err=17; goto end;}
if(HDistr.Alloc(nDiv+1,NF)!=0) { err=18; goto end;}
if(ZDistr.Alloc(nDiv+1,NF)!=0) { err=19; goto end;}
if(ZFour.Alloc(N)!=0) { err=20;} // Вспомогательный массив.

// Нахождение функций сглаживания пространственного спектра.
if(flSm!=0) { if(pFunSm==NULL) { err=21; goto end;}
if((*pFunSm)(WeiP,nWeiP)!=0) { err=22; goto end;}
if((*pFunSm)(WeiN,nWeiN)!=0) { err=23; goto end;}}

// Цикл по слоям для расчёта обратного преобразования Фурье полей.
// Для 'p' волны EDistr=Ex,HDistr=Hy,ZDistr=Ez; Для 's' волны EDistr=Ey,HDistr=Hx,ZDistr=Hz;
if(Task.flPol==p_wave_Pol) pRW=&HDistrFour; else pRW=&EDistrFour;

#pragma omp parallel default(shared) // Параллельная область.
{
// Выделение памяти.
complex *A; clRowMatr *pR,*pRM,*pZF;
A=NULL; A=new complex[NF]; if(A==NULL) { err=13;} // Массив, используемый для обратного преобразования Фурье.
flBreak=0;

#pragma omp for schedule(static) private(kx,j,i,ZFourV) firstprivate(ZFour)
for(j=0;j<=nDiv;j++) { if(flBreak==1) continue;

pR=EDistr.Get(j); if(pR==NULL) { flBreak=1; err=24; continue;} if(pR->N!=NF) { flBreak=1; err=25; continue;}
if(FillArrInvFour(NF,A,EDistrFour.Rows[j],L1,M1)!=0) { flBreak=1; err=26; continue;} // Заполнение массива 'A' для обратного преобразования Фурье для 'EDistr'.
if(flSm!=0) { if(SmoothSpectrum(NF,A,M1,WeiP,nWeiP,L1,WeiN,nWeiN)!=0) { flBreak=1; err=27; continue;}} // Сглаживание спектра.
if(FastFT(A,pR->Vect,M,-1)!=0) { flBreak=1; err=28; continue;} // Обратное преобразование Фурье для полей 'Е'.

pR=HDistr.Get(j); if(pR==NULL) { flBreak=1; err=29; continue;}
if(pR->N!=NF) { flBreak=1; err=30; continue;}
if(FillArrInvFour(NF,A,HDistrFour.Rows[j],L1,M1)!=0) { flBreak=1; err=31; continue;} // Заполнение массива 'A' для обратного преобразования Фурье для 'HDistr'.
if(flSm!=0) { if(SmoothSpectrum(NF,A,M1,WeiP,nWeiP,L1,WeiN,nWeiN)!=0) { flBreak=1; err=32; continue;}} // Сглаживание спектра.
if(FastFT(A,pR->Vect,M,-1)!=0) { flBreak=1; err=33; continue;} // Обратное преобразование Фурье для полей 'H'.

// Вычисление 'Ez' или 'Hz', в зависимости от поляризации.
pR=ZDistr.Get(j); if(pR==NULL) { flBreak=1; err=34; continue;} if(pR->N!=NF) { flBreak=1; err=35; continue;}
pRM=pRW->Get(j); if(pRM==NULL) { flBreak=1; err=36; continue;} if(pRM->N!=N) { flBreak=1; err=37; continue;}
for(i=0;i<N;i++) { kx=(k0X+(double)(i-L1)*dK)/k0; ZFour.Vect[i]=-pRM->Vect[i]*complex(0.,kx);}
// Для p-волны при расчёте 'Ez' неоднородного слоя ищем свёртку с результатом Фурье преобразования от '1/eps'.
if(Task.flPol==p_wave_Pol&&Smpl.flHomogen==0) {
if(ZFourV.Alloc(N)!=0) { flBreak=1; err=38; continue;} // Вспомогательный массив.
ZFourV=Smpl.EkInv*ZFour; pZF=&ZFourV;}
else pZF=&ZFour;

if(FillArrInvFour(NF,A,*pZF,L1,M1)!=0) { flBreak=1; err=40; continue;} // Заполнение матрицы 'A' для обратного преобразования Фурье..
if(flSm!=0) { if(SmoothSpectrum(NF,A,M1,WeiP,nWeiP,L1,WeiN,nWeiN)!=0) { flBreak=1; err=41; continue;}} // Сглаживание спектра.
if(FastFT(A,pR->Vect,M,-1)!=0) { flBreak=1; err=42; continue;} // Обратное преобразование Фурье для полей 'Ez' или 'Hz'.
// Для p-волны при расчёте 'Ez' однородного слоя результат обратного преобразования Фурье ('Ez*eps') делим на диэлектрическую проницаемость.
if(Task.flPol==p_wave_Pol&&Smpl.flHomogen!=0) {
if(abs(Smpl.EpsHomogen)<smValEps) { MessageSmallEps_(numLay); err=messSmallEps;}
*pR/=Smpl.EpsHomogen;}
} // Конец цикла по 'j' (по уровням 'z' внутри слоя).
SAFE_DELETE_ARR(A);
} 
end: SAFE_DELETE_ARR(WeiP); SAFE_DELETE_ARR(WeiN); return err;
}

//-----------------------------------------------------------------------------------------------------------
// Запись в файл информацию о слое.

BYTE strLay::WriteInfo(FILE *fp,const strLaySmpl &Smpl)
{
if(fp==NULL) return 1; if(IsOK()!=0) return 2; if(Smpl.IsOK()!=0) return 3;
if(fprintf(fp,"%s %lf %lf %d\n",Smpl.LayGeom.GeomSymbs,depth,dz,nDiv)<0) return 4; return 0;
}

//-----------------------------------------------------------------------------------------------------------
// Запись в файл распределения полей.

BYTE strLay::WriteFieldDistr(FILE *fp)
{
if(fp==NULL) return 1; if(IsOK()!=0) return 2;
if(EDistr.Write(fp)!=0) return 3;
if(HDistr.Write(fp)!=0) return 4;
if(ZDistr.Write(fp)!=0) return 5;
return 0;
}

//-----------------------------------------------------------------------------------------------------------
// Заполнение массива 'A' для обратного преобразования Фурье.

BYTE FillArrInvFour(int NF,complex *A,const clRowMatr &Row,int L1,int M1)
{
int i,NFh;

if(NF<=0) return 1; if(A==NULL) return 2; if(Row.IsOK()!=0) return 3; if(L1<=0||M1<=0) return 4;
NFh=NF/2; if(NFh<L1+M1+1) return 5;
for(i=0;i<=M1;i++) A[i]=Row.Vect[i+L1]; // Заполнение положительных гармоник (и нулевой).
for(i=0;i<L1;i++) A[NF-1-i]=Row.Vect[L1-1-i]; // Заполнение отрицательных гармоник.
for(i=M1+1;i<NF-L1;i++) A[i]=Cmplx_0; return 0;
}

//-----------------------------------------------------------------------------------------------------------
// Сглаживание спектра.

BYTE SmoothSpectrum(int NF,complex *A,int NP,float *WeiP,int nWeiP,int NN,float *WeiN,int nWeiN)
{
int i,j,NFh;

if(NF<=0) return 1; if(A==NULL) return 2; if(NP<=0) return 3; if(NN<=0) return 4;
NFh=NF/2; if(NFh<NN+NP+1) return 5;
// Сглаживание спектра на положительных гармониках.
if(WeiP!=NULL&&nWeiP>0) { if(NP<nWeiP) return 6;
for(i=0;i<nWeiP;i++) { j=NP-i; A[j]*=(double)WeiP[i];}}
// Сглаживание спектра на отрицательных гармониках.
if(WeiN!=NULL&&nWeiN>0) { if(NN<nWeiN) return 7;
for(i=0;i<nWeiN;i++) { j=NF-NN+i; A[j]*=(double)WeiN[i];}}
return 0;
}

//-----------------------------------------------------------------------------------------------------------
// Нахождение функции для косинусоидального сглаживания.

BYTE GetWeightWinCos(float *Weight,int nPoints)
{
int i; double x,coe;

if(Weight==NULL) return 1; if(nPoints<=0) return 2; coe=1./(double)nPoints;
for(i=0;i<nPoints;i++) { x=(double)(i+1)*coe; Weight[i]=(float)(0.5*(1.-cos(x*M_PI)));}
return 0;
}

//-----------------------------------------------------------------------------------------------------------
// Окно Кайзера-Бесселя.

BYTE GetWeightKaisBess(float *Weight,int nPoints)
{
double GetBesselI0(double x); // Расчёт модифицированной функции Бесселя 'I0'.

int i; double x,v,coe,coeF,alp;

if(Weight==NULL) return 1; if(nPoints<=0) return 2;
coe=1.0F/(double)nPoints; alp=M_PI*coeWinKaisBess; v=GetBesselI0(alp); if(v<1.e-12) return 3; coeF=1./v;
for(i=0;i<nPoints;i++) { x=(double)(nPoints-1-i)*coe; v=1.-x*x; if(v>0.) v=sqrt(v); else v=0.; v*=alp;
v=GetBesselI0(v)*coeF; Weight[i]=(float)v;}
return 0;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#define nSumBessI0 36 // Число членов в сумме при вычислении ряда для нахождении модифицированной функции Бесселя I0.

//-----------------------------------------------------------------------------------------------------------
// Расчёт модифицированной функции Бесселя 'I0'.

double GetBesselI0(double x)
{
BYTE i; double s,xh,p;

xh=x*0.5; s=1.; p=1.; for(i=1;i<nSumBessI0;i++) { p*=xh/(double)i; s+=p*p;} return s;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Класс фотонного кристалла.

//-----------------------------------------------------------------------------------------------------------
// Конструктор.

clCrystal::clCrystal(void)
{
Zero();
}

//-----------------------------------------------------------------------------------------------------------
// Деструктор.

clCrystal::~clCrystal(void)
{
Free();
}

//-----------------------------------------------------------------------------------------------------------
// Обнуление переменных.

void clCrystal::Zero(void)
{
Lays=NULL; laysNum=0;
}

//-----------------------------------------------------------------------------------------------------------
// Освобождение памяти.

void clCrystal::Free(void)
{
SAFE_DELETE_ARR(Lays); laysNum=0;
}

//-----------------------------------------------------------------------------------------------------------
// Выделение памяти для массива.

BYTE clCrystal::Alloc(int laysNum_)
{
Free(); if(laysNum_<=0) return 1; laysNum=laysNum_;
Lays=new strLay[laysNum]; if(Lays==NULL) return 2; return 0;
} 

//-----------------------------------------------------------------------------------------------------------
// Проверки.

BYTE clCrystal::IsOK(void) const
{
int i;

if(laysNum<=0) return 1; if(Lays==NULL) return 2; for(i=0;i<laysNum;i++) if(Lays[i].IsOK()!=0) return 3;
if(Samples.IsOK()!=0) return 4; return 0;
}

//-----------------------------------------------------------------------------------------------------------
// Чтение исходных данных.

BYTE clCrystal::ReadInpData(char *FileName_Task)
{
if(File_Name_Task==NULL) return 1;
if(Samples.Task.Read(FileName_Task)!=0) return 2; // Чтение параметров задачи из файла.
if(Samples.ReadInpData()!=0) return 3; // Чтение данных о кристалле из файлов.
if(ReadGeom(Samples.Task.FName_Cryst)!=0) return 5; // Чтение данных о кристалле из файлов.
return 0;
}

//-----------------------------------------------------------------------------------------------------------
// Чтение данных о кристалле из файла.

BYTE clCrystal::ReadGeom(char *FileName_Cryst)
{
BYTE ReadString(FILE *file,char *pStr,int lnStr); // Чтение строки из текстового файла.

BYTE err,fl; int i,j,laysNum_,num; char *str,*fname; FILE *fp;
str=NULL; fp=NULL; err=0;

// Выделение памяти для вспомогательной строки.
str=new char[laySize_Max]; if(str==NULL) { err=1; goto end;} str[laySize_Max-1]='\0';
fname=GetFName(Samples.Task.FName_Cryst,extCrystal); if(fname==NULL) { err=2; goto end;}
fp=fopen(fname,"r"); if(fp==NULL) { printf("Can not open file: %s !\n",fname); err=3; goto end;}
if(fscanf(fp,"%d",&laysNum_)==EOF) { err=4; goto end;} 
if(Alloc(laysNum_)!=0) { err=5; goto end;} // Выделение памяти для массива слоев. 
if(fscanf(fp,"\n")==EOF) { err=6; goto end;}

// Послойное чтение.
for(i=0;i<laysNum;i++) { Lays[i].numLay=i;
if(ReadString(fp,str,laySize_Max-1)!=0) { err=7; goto end;} // Чтение строки из текстового файла.
num=(int)strlen(str); if(num<=0) { err=8; goto end;} 

fl=1; for(j=0;j<Samples.laysNum;j++) // Поиск соответствия в классе образцов слоёв.
if(strcmp(Samples.Lays[j].LayGeom.GeomSymbs,str)==0) { Lays[i].numSmpl=j; fl=0; break;}
if(fl==1) { printf("Error: There is no match between samples and Lay %d from file\n",i); err=9; goto end;}
if(fscanf(fp,"%lf\n",&Lays[i].depth)==EOF) { err=10; goto end;} // Чтение толщины слоя.
} // Конец цикла по 'i'.

end: SAFE_DELETE_ARR(str); SAFE_DELETE_ARR(fname); SAFE_CLOSE(fp); return err;
}

//-----------------------------------------------------------------------------------------------------------
// Вычисление матриц DExp00Inv, DExp11 для всех слоев.

BYTE clCrystal::FindGammaExp(void)
{
int i;

if(IsOK()!=0) return 1; if(Samples.IsOK()!=0) return 2;
for(i=0;i<laysNum;i++) if(Lays[i].FindGammaExp(Samples.Lays[Lays[i].numSmpl])!=0) return 3; return 0;
}

//-----------------------------------------------------------------------------------------------------------
// Проверка блоков матриц 'Gamma' во всех слоях.

BYTE clCrystal::IsOKGamma(void)
{
int i;

if(IsOK()!=0) return 1; if(Samples.IsOKGamma()!=0) return 2;
for(i=0;i<laysNum;i++) if(Lays[i].IsOKGammaExp()!=0) return 3; return 0;
}

//-----------------------------------------------------------------------------------------------------------
// Вычисление матриц 'C'.

BYTE clCrystal::CompC(void)
{
void MessageSmallEps(char *str); // Сообщение о близком к нулю значении диэлектрической проницаемости.
void MessageSingLam(char *str); // Сообщение о близком к нулю собственном значении.

BYTE ind; int i,N,L1,M1; double dK,k0X,k0,k2,ki; complex kZ,kZI,eps,c,cI,epsInc,epsOut; char *cp;
clRowMatr *pRM,*pRMI;

L1=Samples.Task.L1; if(L1<=0) return 1;
M1=Samples.Task.M1; if(M1<=0) return 2;
N=L1+M1+1;
if(Samples.Task.wLength<=0.) return 3; if(Samples.Period<=0.) return 4;
epsInc=Samples.Task.EpsInc; epsOut=Samples.Task.EpsOut;

// Расчёт волнового вектора и его компоненты по оси X.
k0=2.*M_PI/Samples.Task.wLength; k2=k0*k0; k0X=k0*sin(Samples.Task.Theta*M_PI/180.);
dK=2.*M_PI/Samples.Period;

// Нахождение диагональных матриц 'С11Inc', 'С11Out'.
// Выделение памяти для матриц.
for(ind=0;ind<3;ind++) {
switch(ind) { default: continue;
case 0: pRM=&CInc; break;
case 1: pRM=&CIncInv; break;
case 2: pRM=&COut; break;}
if(pRM->Alloc(N)!=0) return 5;}

switch(Samples.Task.flPol) { default: return 6;
// Вектор E в плоскости падения.
case p_wave_Pol:
for(ind=0;ind<2;ind++) {
if(ind==0) { eps=epsInc; pRM=&CInc; pRMI=&CIncInv; cp="Medium of Incidence";}
else { eps=epsOut; pRM=&COut; pRMI=NULL; cp="Medium of Output";}
if(abs(eps)<smValEps) { MessageSmallEps(cp); return messSmallEps;}
c=complex(0.,k0)*eps; cI=Inv(c);
for(i=0;i<N;i++) { ki=k0X+(double)(i-L1)*dK; kZ=sqrt(complex(k2,0.)*eps-complex(ki*ki,0.));
if(abs(kZ)<smValLam) { MessageSingLam(cp); return messSingLam;}
kZI=Inv(kZ); pRM->Vect[i]=c*kZI; if(pRMI!=NULL) pRMI->Vect[i]=kZ*cI;}}
break;

// Вектор E перпендикулярен плоскости падения.
case s_wave_Pol:
for(ind=0;ind<2;ind++) {
if(ind==0) { eps=epsInc; pRM=&CInc; pRMI=&CIncInv; cp="Medium of Incidence";}
else { eps=epsOut; pRM=&COut; pRMI=NULL; cp="Medium of Output";}
if(abs(eps)<smValEps) { MessageSmallEps(cp); return messSmallEps;}
c=complex(0.,-1./k0); cI=Inv(c);
for(i=0;i<N;i++) { ki=k0X+(double)(i-L1)*dK; kZ=sqrt(complex(k2,0.)*eps-complex(ki*ki,0.));
if(abs(kZ)<smValLam) { MessageSingLam(cp); return messSingLam;}
kZI=Inv(kZ); pRM->Vect[i]=c*kZ; if(pRMI!=NULL) pRMI->Vect[i]=kZI*cI;}}
break;
}

return 0;
}

//-----------------------------------------------------------------------------------------------------------
// Проверка блоков матриц 'C'.

BYTE clCrystal::IsOK_C(void)
{
BYTE ind; int N; clRowMatr *pRM;

N=Samples.Task.L1+Samples.Task.M1+1; if(N<=0) return 1;
for(ind=0;ind<3;ind++) {
switch(ind) { default: continue;
case 0: pRM=&CInc; break;
case 1: pRM=&CIncInv; break;
case 2: pRM=&COut; break;}
if(pRM->IsOK()!=0) return 2; if(pRM->N!=N) return 3;}
return 0;
}

//-----------------------------------------------------------------------------------------------------------
// Нахождение амплитуд прошедших и отражённых волн.

BYTE clCrystal::FindRTW(void)
{
BYTE InvMatr(const clMatrix &M,clMatrix &Inv); // Нахождение обратной матрицы методом исключения Гаусса с выбором главного элемента.

BYTE ind,cnt,be; int i,N,k,L1,M1; double v; clMatrix *pGC,*pGP,*pZC,*pZP,*pM,(*pGam)[2][2],GM[2],ZM[2],V,Inv;
struct strLaySmpl *pLS; struct strLay *pL;

if(IsOKGamma()!=0) return 1; if(IsOK_C()!=0) return 2;
L1=Samples.Task.L1; if(L1<=0) return 3;
M1=Samples.Task.M1; if(M1<=0) return 4;
N=L1+M1+1; if(N<=0) return 5;

// Выделение памяти.
for(k=0;k<2;k++) for(ind=0;ind<2;ind++) { if(ind==0) pM=GM+k; else pM=ZM+k; if(pM->Alloc(N,N)!=0) return 6;}
if(R.Alloc(N)!=0) return 7;
if(T.Alloc(N)!=0) return 8;
if(J.Alloc(N)!=0) return 9;
if(Inv.Alloc(N,N)!=0) return 10;

cnt=0; pGP=GM+cnt; pZP=ZM+cnt; pGC=GM+1-cnt; pZC=ZM+1-cnt; // Установка указателей на матрицы G,Z.

*pGP=-COut; *pZP=Cmplx_1; // Инициализация начальных матриц G и Z после умножения на 'C'.
if(Samples.Task.flTypComp==0) printf("."); // Полоса прогресса.

// Расчёт G0,Z0 по реккурентным формулам.
// Идём по слоям от последнего слоя к первому, перебирая матрицы в порядке: 'PInv','DiagExp','P'.
for(i=laysNum-1;i>=0;i--) { pL=Lays+i; if(pL->N!=N) return 11;
if(pL->DExp00Inv.N!=N||pL->DExp11.N!=N) return 12;
pLS=Samples.Get(pL->numSmpl); if(pLS==NULL) return 13; if(pLS->N!=N) return 14;

// Итерации для матриц 'Gk', 'Zk' с матрицами 'Gamma' в порядке: 'PInv','DiagExp','P'.
for(k=0;k<3;k++) {
pL->GK[k]=*pGP; // Запоминаем матрицы 'Gk' в классе слоя.

// Матрица 'Gamma' - диагональная 'DiagExp' с экспонентами на диагонали.
if(k==1) { *pGC=pL->DExp11%(*pGP)%pL->DExp00Inv; *pZC=*pZP%pL->DExp00Inv;}
// Матрица 'Gamma' - 'PInv' или 'P'.
else { if(k==0) pGam=&pLS->GamPInv; else pGam=&pLS->GamP;
V=(*pGam)[0][0]+(*pGam)[0][1]*(*pGP);
be=InvMatr(V,Inv); if(be!=0) { printf("\nError in 'InvMatr' %d Layer: %d k: %d\n",be,i,k); return 15;}
*pGC=((*pGam)[1][0]+(*pGam)[1][1]*(*pGP))*Inv; *pZC=*pZP*Inv;}

cnt=1-cnt; pGP=GM+cnt; pZP=ZM+cnt; pGC=GM+1-cnt; pZC=ZM+1-cnt; // Установка указателей. Меняем местами 'GP','GC' и 'ZP','ZC'.
if(Samples.Task.flTypComp==0) printf("."); // Полоса прогресса.
} // Конец цикла по 'k' - по матрицам 'PInv','DiagExp','P' внутри каждого слоя кристалла.
pL->GK[3]=*pGP; // Запоминаем матрицу 'Gk' после последней итерации в классе слоя.
} // Конец цикла по слоям кристалла.

// Расчёт по реккурентным формулам с матрицей 'C^-1' среды, откуда падает свет.
V=Cmplx_1-CIncInv%(*pGP);
be=InvMatr(V,Inv); if(be!=0) { printf("\nError in Final 'InvMatr' %d\n",be); return 16;}
*pGC=(Cmplx_1+CIncInv%(*pGP))*Inv; *pZC=*pZP*Inv*complex(2.,0.);
cnt=1-cnt; pGP=GM+cnt; pZP=ZM+cnt; // Установка указателей для 'GP', 'ZP': G0=*pGP, Z0=*pZP.
if(Samples.Task.flTypComp==0) printf("."); // Полоса прогресса.

J=Cmplx_0; // Начальное обнуление столбца амплитуд падающих волн.
switch(Samples.Task.flPol) { default: return 17;
case p_wave_Pol: v=cos(Samples.Task.Theta/180.*M_PI); break;
case s_wave_Pol: v=1.; break;}
J.Vect[L1]=complex(v,0.);

R=*pGP*J; T=*pZP*J; // Нахождение амплитуд прошедшей и отражённой волн.
if(Samples.Task.flTypComp==0) printf("."); // Полоса прогресса.
return 0;
}

//-----------------------------------------------------------------------------------------------------------
// Нахождение столбцов амплитуд разложения Фурье распределения поля на границах слоёв.

BYTE clCrystal::FindBoundsFieldDistr(void)
{
BYTE InvMatr(const clMatrix &M,clMatrix &Inv); // Нахождение обратной матрицы методом исключения Гаусса с выбором главного элемента.

BYTE k,cnt,be; int i,N; struct strLay *pL; struct strLaySmpl *pLS; clRowMatr *pC,*pP,*pD,RV[2],RMV;
clMatrix *pGK,(*pGam)[2][2],MV,Inv;

// Проверки.
N=Samples.Task.L1+Samples.Task.M1+1; if(N<=0) return 1; if(IsOK()!=0) return 2;
if(CInc.IsOK()!=0) return 3; if(CInc.N!=N) return 4;

for(k=0;k<3;k++) {
switch(k) { default: continue;
case 0: pC=&J; break; // Падающая волна.
case 1: pC=&T; break; // Прошедшая волна.
case 2: pC=&R; break;} // Отражённая волна.
if(pC->IsOK()!=0) return 5; if(pC->N!=N) return 6;}

// Выделение памяти для вспомогательной матрицы.
if(Inv.Alloc(N,N)!=0) return 7;

//...........................................................................................................
// Считаем столбцы полей 'E' от начала кристалла в сторону конца. Столбцы полей 'H' получаем умножая матрицу 'Gk' на столбцы полей 'E'.
// Итерации по слоям фотонного кристалла и по матрицам 'P', 'DExp', 'PInv' внутри слоя кристалла.
cnt=0; pP=RV+cnt; pC=RV+1-cnt; // Установка указателей.
*pP=J+R; // Начальный столбец коэффициентов Фурье разложения полей E на передней границе первого слоя.
printf("."); // Полоса прогресса.

// Цикл итераций по слоям кристалла.
for(i=0;i<laysNum;i++) { pL=Lays+i; if(pL->N!=N) return 8;
pD=&pL->DExp00Inv; if(pD->IsOK()!=0) return 9; if(pD->N!=N) return 10; // Проверка диагональной матрицы 'DExp00Inv'.
pLS=Samples.Get(pL->numSmpl); if(pLS==NULL) return 11; if(pLS->N!=N) return 12;
if(pLS->IsOKGamma()!=0) return 13; // Проверка матриц 'Gamma'.

pGK=pL->GK+3; if(pGK->IsOK()!=0) return 14; if(pGK->Nx!=N||pGK->Ny!=N) return 15; // Выбираем соответствующую матрицу 'Gk'.
pL->EFour[0]=*pP; // Запоминаем столбец коэффициентов Фурье разложения полей E на передней границе i-го слоя.
pL->HFour[0]=*pGK*(*pP); // Столбец коэффициентов Фурье разложения полей H получаем умножая матрицу 'Gk' на столбец Фурье разложения полей E.

// Итерации по матрицам 'P', 'DExp', 'PInv' внутри слоя кристалла.
for(k=0;k<3;k++) {
pGK=pL->GK+2-k; if(pGK->IsOK()!=0) return 16; if(pGK->Nx!=N||pGK->Ny!=N) return 17; // Выбираем соответствующую матрицу 'Gk'.
if(k==1) *pC=*pD*(*pP); // Матрица 'Gamma' - диагональная 'DiagExp' с экспонентами на диагонали.
else { if(k==0) pGam=&pLS->GamP; else pGam=&pLS->GamPInv; // Матрица 'Gamma' - 'P' или 'PInv'.
be=InvMatr((*pGam)[0][0],Inv);
if(be!=0) { printf("\nError in 'InvMatr' for 'Gam00' %d Layer: %d k: %d\n",be,i,k); return 18;}
RMV=Inv*(*pP);
MV=(*pGam)[0][1]*(*pGK); MV=Inv*MV; MV+=Cmplx_1;
be=InvMatr(MV,Inv);
if(be!=0) { printf("\nError in 'InvMatr' for 'MV' %d Layer: %d k: %d\n",be,i,k); return 19;}
*pC=Inv*RMV;} // Находим следующий столбец 'X'.

cnt=1-cnt; pP=RV+cnt; pC=RV+1-cnt; // Установка указателей.
pL->EFour[k+1]=*pP; // Запоминаем столбец 'X' (или коэффициентов Фурье разложения полей E) после очередной итерации для i-го слоя.
pL->HFour[k+1]=*pGK*(*pP); // Столбец 'Y' (или коэффициентов Фурье разложения полей H) получаем умножая матрицу 'Gk' на столбец 'X'.
printf("."); // Полоса прогресса.
} // Конец цикла по 'k' - по матрицам 'P','DiagExp','PInv' внутри каждого слоя кристалла.
} // Конец цикла по слоям кристалла.

return 0;
}

//-----------------------------------------------------------------------------------------------------------
// Нахождение распределения поля внутри кристалла.

BYTE clCrystal::FindFullFieldDistr(void)
{
int i; struct strLay *pL;

if(IsOK()!=0) return 1;
// Нахождение столбцов амплитуд разложения Фурье распределения поля на границах слоёв.
printf("Computing of the Fourier Harmonics on Boundaries of Layers\n");
if(FindBoundsFieldDistr()!=0) return 2; printf(" - Completed\n");
// Нахождение коэффициентов Фурье разложения полей Е и H и звтем нахождение распределения поля в каждом слое.
for(i=0;i<laysNum;i++) { printf("Computing of the Layer %d",i);
pL=Lays+i; if(pL->FindFullFieldDistr(Samples)!=0) return 3;
printf(" - Completed\n");}
return 0;
}

//-----------------------------------------------------------------------------------------------------------
// Запись в файл информации о кристалле и задаче.

BYTE clCrystal::WriteInfo(FILE *fp)
{
int i;

if(fp==NULL) return 1; if(IsOK()!=0) return 2;
fprintf(fp,"%d\n",laysNum);
fprintf(fp,"%lf\n",Samples.Period);
for(i=0;i<laysNum;i++) if(Lays[i].WriteInfo(fp,Samples.Lays[Lays[i].numSmpl])!=0) return 3;
fprintf(fp,"\n"); return 0;
}

//-----------------------------------------------------------------------------------------------------------
// Вывод углового распределения по порядкам дифракции.

BYTE clCrystal::AngDistrOutput(void)
{
BYTE typPol,err; int i,n,L1;
double ThetaT,ThetaR,coeA,tr,ref,sumR,sumT,pInc,pRef,pTr,pLoss,k0,k2,k0X,dK,kx,kz,ct,v;
complex c,KzR,KzT,epsInc,epsOut; char *fname; FILE *fp;

if(Samples.Period<=0.) return 1; if(Samples.Task.wLength<=0.) return 2; if(T.IsOK()!=0) return 3;
fname=NULL; fp=NULL; err=0;
fname=GetFName(Samples.Task.FName_Cryst,extOutAngle); if(fname==NULL) { err=4; goto end;}
fp=fopen(fname,"w"); if(fp==NULL) { printf("Can not open file: %s !\n",fname); err=5; goto end;}
typPol=Samples.Task.flPol;

// Печать информации о задаче.
if(Samples.Matter.WriteInfo(fp)!=0) { err=6; goto end;}
if(WriteInfo(fp)!=0) { err=7; goto end;}
if(Samples.Task.WriteInfo(fp)!=0) { err=8; goto end;}
coeA=180./M_PI; L1=Samples.Task.L1;

fprintf(fp,"n   ThetaR R   ThetaT T   KzR KzT\n");
sumR=sumT=0.;
dK=2.*M_PI/Samples.Period; k0=2.*M_PI/Samples.Task.wLength; k2=k0*k0; k0X=k0*sin(Samples.Task.Theta*M_PI/180.);
epsInc=Samples.Task.EpsInc;
epsOut=Samples.Task.EpsOut;

// Печать амплитуд прошедшей и отражённой волн в различных порядках дифракционных спектров.
for(i=0;i<T.N;i++) { n=i-L1; tr=abs(T.Vect[i]); ref=abs(R.Vect[i]);
kx=k0X+(double)n*dK;
c=complex(k2,0.)*epsInc-complex(kx*kx,0.); KzR=sqrt(c); kz=real(KzR); ThetaR=atan2(kx,kz);
c=complex(k2,0.)*epsOut-complex(kx*kx,0.); KzT=sqrt(c); kz=real(KzT); ThetaT=atan2(kx,kz);
if(fprintf(fp,"%d   %lf %lf   %lf %lf   %lf %lf   %lf %lf\n",n,ThetaR*coeA,ref,ThetaT*coeA,tr,
real(KzR),imag(KzR),real(KzT),imag(KzT))<0) { err=9; goto end;}
// Суммируем мощности прошедших и отражённых волн.
if(real(KzR)>imag(KzR)) { ct=cos(ThetaR); v=ref; if(typPol==p_wave_Pol) v/=ct; sumR+=v*v*ct;}
if(real(KzT)>imag(KzT)) { ct=cos(ThetaT); v=tr; if(typPol==p_wave_Pol) v/=ct; sumT+=v*v*ct;}
}

// Умножаем на множитель, зависящий от диэлектрической проницаемости среды, нормируем на мощность падающей волны.
ct=cos(Samples.Task.Theta/coeA);
v=real(sqrt(epsInc)); pInc=v*ct; pRef=sumR*v; pTr=sumT*real(sqrt(epsOut));
pRef/=pInc; pTr/=pInc; pLoss=1.-pRef-pTr;

// Печать суммарных мощностей для прошедших и отражённых волн и мощности потерь.
if(fprintf(fp,"Relative Power - Transm.,Refl.: %lf %lf Losses: %lf\n",pTr,pRef,pLoss)<0) { err=11; goto end;}

end: SAFE_DELETE_ARR(fname); SAFE_CLOSE(fp); return err;
}

//-----------------------------------------------------------------------------------------------------------
// Расчет долей прошедшего и отраженного света.

BYTE clCrystal::CompPartRefTr(double *pRef,double *pTr)
{
BYTE typPol,err; int i,n,L1; double ThetaT,ThetaR,coeA,tr,ref,sumR,sumT,pInc,k0,k2,k0X,dK,kx,kz,ct,v;
complex c,KzR,KzT,epsInc,epsOut;

if(pRef==NULL||pTr==NULL) return 1;
if(Samples.Period<=0.) return 1; if(Samples.Task.wLength<=0.) return 2; if(T.IsOK()!=0) return 3;
err=0;
typPol=Samples.Task.flPol;

coeA=180./M_PI; L1=Samples.Task.L1;

sumR=sumT=0.;
dK=2.*M_PI/Samples.Period; k0=2.*M_PI/Samples.Task.wLength; k2=k0*k0; k0X=k0*sin(Samples.Task.Theta*M_PI/180.);
epsInc=Samples.Task.EpsInc;
epsOut=Samples.Task.EpsOut;

// Печать амплитуд прошедшей и отражённой волн в различных порядках дифракционных спектров.
for(i=0;i<T.N;i++) { n=i-L1; tr=abs(T.Vect[i]); ref=abs(R.Vect[i]);
kx=k0X+(double)n*dK;
c=complex(k2,0.)*epsInc-complex(kx*kx,0.); KzR=sqrt(c); kz=real(KzR); ThetaR=atan2(kx,kz);
c=complex(k2,0.)*epsOut-complex(kx*kx,0.); KzT=sqrt(c); kz=real(KzT); ThetaT=atan2(kx,kz);
// Суммируем мощности прошедших и отражённых волн.
if(real(KzR)>imag(KzR)) { ct=cos(ThetaR); v=ref; if(typPol==p_wave_Pol) v/=ct; sumR+=v*v*ct;}
if(real(KzT)>imag(KzT)) { ct=cos(ThetaT); v=tr; if(typPol==p_wave_Pol) v/=ct; sumT+=v*v*ct;}
}

// Умножаем на множитель, зависящий от диэлектрической проницаемости среды, нормируем на мощность падающей волны.
ct=cos(Samples.Task.Theta/coeA);
v=real(sqrt(epsInc)); pInc=v*ct; *pRef=sumR*v; *pTr=sumT*real(sqrt(epsOut));
*pRef/=pInc; *pTr/=pInc;
return 0;
}

//-----------------------------------------------------------------------------------------------------------
// Вывод в файл всех спектров.

BYTE clCrystal::SpectrOutput(int n,double *T_Sp,double *R_Sp,double *A_Sp)
{
BYTE err; int i,j; double wL_St,dwL,*Row; char *fname; FILE *fp;

if(T_Sp==NULL) return 1; if(R_Sp==NULL) return 2; if(A_Sp==NULL) return 3;
if(n<=0) return 4; if(Samples.Task.flTypComp!=1) return 5;

err=0; wL_St=Samples.Task.wLength_St; dwL=Samples.Task.dwLength;

// Цикл по файлам.
for(i=0;i<3;i++) {

fname=NULL;
switch(i) { default: break;
case 0: fname=GetFName(Samples.Task.FName_Cryst,extOutR_Sp); Row=R_Sp; break;
case 1: fname=GetFName(Samples.Task.FName_Cryst,extOutT_Sp); Row=T_Sp; break;
case 2: fname=GetFName(Samples.Task.FName_Cryst,extOutA_Sp); Row=A_Sp; break;
}
if(fname==NULL) { err=6; goto end;}
fp=NULL; fp=fopen(fname,"w"); if(fp==NULL) { printf("Can not open file: %s !\n",fname); err=7; goto end;}

switch(i) { default: break;
case 0: if(fprintf(fp,"wLength R\n")<0) { err=8; goto end;} break;
case 1: if(fprintf(fp,"wLength T\n")<0) { err=8; goto end;} break;
case 2: if(fprintf(fp,"wLength A\n")<0) { err=8; goto end;} break;
}
// Цикл вывода результатов в файл.
for(j=0;j<n;j++) if(fprintf(fp,"%lf %.8lf\n",wL_St+j*dwL,Row[j])<0) { err=9; goto end;}
SAFE_CLOSE(fp);
} // Конец цикла по 'i'.

end: SAFE_CLOSE(fp); SAFE_DELETE_ARR(fname); return err;
}

//-----------------------------------------------------------------------------------------------------------
// Вывод в файл всех спектров.

BYTE clCrystal::AnglDispOutput(int nRows,int n,clRows &T_Sp,clRows &R_Sp,clRows &A_Sp)
{
int i,j,k; double wL_St,dwL; char *fname; FILE *fp; clRows Row; clTask *Tsk;

if(T_Sp.IsOK()!=0) return 1; if(R_Sp.IsOK()!=0) return 2; if(A_Sp.IsOK()!=0) return 3;
if(n<=0) return 4; if(nRows<=0) return 5; if(Samples.Task.flTypComp!=2) return 6;

wL_St=Samples.Task.wLength_St; dwL=Samples.Task.dwLength;

// Цикл по файлам.
for(i=0;i<3;i++) {

fname=NULL;
switch(i) { default: break;
case 0: fname=GetFName(Samples.Task.FName_Cryst,extOutR_AnglDisp); Row=R_Sp; break;
case 1: fname=GetFName(Samples.Task.FName_Cryst,extOutT_AnglDisp); Row=T_Sp; break;
case 2: fname=GetFName(Samples.Task.FName_Cryst,extOutA_AnglDisp); Row=A_Sp; break;
}
if(fname==NULL) return 7;
fp=NULL; fp=fopen(fname,"w"); if(fp==NULL) { printf("Can not open file: %s !\n",fname); return 8;}
Tsk=&Samples.Task;
if(fprintf(fp,"%lf %lf %lf %lf %lf %lf %d %d\n",Tsk->wLength_St,Tsk->wLength_Fin,Tsk->dwLength,Tsk->Theta_St,Tsk->Theta_Fin,Tsk->dTheta,nRows,n)<0) return 9;

// Цикл вывода результатов в файл.
for(j=0;j<nRows;j++) { for(k=0;k<n;k++) if(fprintf(fp,"%lf ",Row.Rows[j].Vect[k].re)<0) return 10; if(fprintf(fp,"\n")<0) return 11;}
SAFE_CLOSE(fp);
} // Конец цикла по 'i'.

SAFE_CLOSE(fp); SAFE_DELETE_ARR(fname); return 0;
}

//-----------------------------------------------------------------------------------------------------------
// Вывод в файл распределения полей в кристалле.

BYTE clCrystal::FieldOutput(void)
{
BYTE err; int i; char *fname; FILE *fp;

if(IsOK()!=0) return 1;
fname=NULL; fp=NULL; err=0;
fname=GetFName(Samples.Task.FName_Cryst,extOutFieldDistr); if(fname==NULL) { err=2; goto end;}
fp=fopen(fname,"w"); if(fp==NULL) { printf("Can not open file: %s !\n",fname); err=3; goto end;}
// Печать информации о задаче.
if(Samples.Task.WriteInfo(fp)!=0) { err=4; goto end;}
if(Samples.Matter.WriteInfo(fp)!=0) { err=5; goto end;}
if(WriteInfo(fp)!=0) { err=6; goto end;}
// Запись в файл распределения полей.
//xxx for(i=0;i<laysNum;i++) if(Lays[i].WriteFieldDistr(fp)!=0) { err=7; goto end;}

for(i=0;i<laysNum;i++) if(Lays[i].EDistr.Write(fp)!=0) { err=7; goto end;}
for(i=0;i<laysNum;i++) if(Lays[i].HDistr.Write(fp)!=0) { err=8; goto end;}
for(i=0;i<laysNum;i++) if(Lays[i].ZDistr.Write(fp)!=0) { err=9; goto end;}

end: SAFE_DELETE_ARR(fname); SAFE_CLOSE(fp); return err;
}

//-----------------------------------------------------------------------------------------------------------
// Расчёт для фотонного кристалла.

BYTE clCrystal::Compute(double *Time)
{
if(Time==NULL) return 1; 

double Time1,Time2,Time3;

*Time=Time1=omp_get_wtime(); 

printf("Computing of Fast Fourier Transform...\n");
if(Samples.FastFTAll()!=0) return 1; // Быстрое преобразование Фурье для всего кристалла.
printf("Completed\n");

printf("Creating of Matrix 'Mk'...\n");
if(Samples.CompMk()!=0) return 2; // Вычисление матрицы Mk.
printf("Completed\n");

printf("Computing of Eigen Vectors and Values of 'Mk'...\n");
if(Samples.FindEigVectMk()!=0) return 3; // Нахождение собственных значений матрицы Mk.
printf("Completed\n");

Time1=omp_get_wtime()-Time1; Time2=omp_get_wtime();

printf("Computing of Matrices C...");
if(CompC()!=0) return 4; // Вычисление матриц 'C'.
printf(" - Completed\n");

if(FindGammaExp()!=0) return 5; // Вычисление матриц DExp00Inv, DExp11.

printf("Computing of Transmitted and Reflected Wawes\n");
if(FindRTW()!=0) return 6; // Нахождение амплитуд прошедших и отражённых волн.
printf(" - Completed\n");

Time2=omp_get_wtime()-Time2; Time3=omp_get_wtime();

printf("Computing of the field distribution inside of the crystal...\n");
if(FindFullFieldDistr()!=0) return 7; // Нахождение распределения поля внутри кристалла.
printf("Completed\n");

Time3=omp_get_wtime()-Time3; *Time=omp_get_wtime()-*Time;

printf("Writing of the angle distribution...");
if(AngDistrOutput()!=0) return 8; // Вывод углового распределения по порядкам дифракции.
printf(" - Completed\n");

printf("Writing of the field distributions...");
if(FieldOutput()!=0) return 9; // Вывод в файл распределения полей в кристалле.
printf(" - Completed\n");

printf("\nTime of the first section %lf\n",Time1); printf("Time of the second section %lf\n",Time2);
printf("Time of the third section %lf\n",Time3);
return 0;
}

//-----------------------------------------------------------------------------------------------------------
// Расчёт спектров прохождения, отражения и поглощения для фотонного кристалла.

BYTE clCrystal::ComputeSpectr(double *Time)
{
BYTE err; int n,i; double dwL,*wL; double *T_Sp,*R_Sp,*A_Sp;

if(Time==NULL) return 1; *Time=omp_get_wtime();
if(Samples.Task.flTypComp!=1) return 2;

T_Sp=R_Sp=A_Sp=NULL;

dwL=Samples.Task.dwLength; n=int((Samples.Task.wLength_Fin-Samples.Task.wLength_St)/dwL)+1;

err=0;
T_Sp=new double[n]; if(T_Sp==NULL) { err=3; goto end;}
R_Sp=new double[n]; if(T_Sp==NULL) { err=4; goto end;}
A_Sp=new double[n]; if(T_Sp==NULL) { err=5; goto end;}

printf("\nComputing of spectrum...\n");
printf("Number of points %d\n",n);

// Расчет спектра.
wL=&Samples.Task.wLength;
*wL=Samples.Task.wLength_St;
for(i=0;i<n;i++) { 
if(Samples.FastFTAll()!=0) { err=6; goto end;} // Быстрое преобразование Фурье для всего кристалла.
if(Samples.CompMk()!=0) { err=7; goto end;} // Вычисление матрицы Mk.
if(Samples.FindEigVectMk()!=0) { err=8; goto end;} // Нахождение собственных значений матрицы Mk.
if(CompC()!=0) { err=9; goto end;} // Вычисление матриц 'C'.
if(FindGammaExp()!=0) { err=10; goto end;} // Вычисление матриц DExp00Inv, DExp11.
if(FindRTW()!=0) { err=11; goto end;} // Нахождение амплитуд прошедших и отражённых волн.

if(CompPartRefTr(R_Sp+i,T_Sp+i)!=0) { err=12; goto end;} A_Sp[i]=1.-R_Sp[i]-T_Sp[i];
*wL+=dwL; printf(".");
if(Samples.Matter.SetEps(*wL)!=0) { err=13; goto end;}
if(Samples.Task.SetEpsIncOut(Samples.Matter)!=0) { err=14; goto end;}
}

*Time=omp_get_wtime()-*Time;
printf("Completed\n");

printf("Writing spectrum in files...");
if(SpectrOutput(n,T_Sp,R_Sp,A_Sp)!=0) { err=15; goto end;} // Вывод в файл всех спектров.
printf(" - Completed\n");

end: SAFE_DELETE_ARR(T_Sp); SAFE_DELETE_ARR(R_Sp); SAFE_DELETE_ARR(A_Sp); return err;
}

//-----------------------------------------------------------------------------------------------------------
// Расчёт угловой диспперсии для фотонного кристалла.

BYTE clCrystal::ComputeAnglDisp(double *Time)
{
BYTE error; int n,nRows,i,j,cntErr; double dwL,*wL,dTh,*Th; clRows T_Sp,R_Sp,A_Sp;

if(Time==NULL) return 1; *Time=omp_get_wtime();
if(Samples.Task.flTypComp!=2) return 2;

dwL=Samples.Task.dwLength; dTh=Samples.Task.dTheta; 
if(dwL==0) nRows=1; else nRows=int((Samples.Task.wLength_Fin-Samples.Task.wLength_St)/dwL)+1; if(dTh==0) n=1; else n=int((Samples.Task.Theta_Fin-Samples.Task.Theta_St)/dTh)+1;
printf("\nn nrows %d %d\n",n,nRows);
if(n<=0||nRows<=0) return 3;

if(T_Sp.Alloc(nRows,n)!=0) return 4; if(R_Sp.Alloc(nRows,n)!=0) return 5; if(A_Sp.Alloc(nRows,n)!=0) return 6;

printf("\nComputing of angle dispersion...\n"); printf("Number of points %d * %d = %d\n",n,nRows,n*nRows);

// Расчет спектра.
wL=&Samples.Task.wLength; Th=&Samples.Task.Theta;
*wL=Samples.Task.wLength_St; cntErr=0;
for(j=0;j<nRows;j++) { *Th=Samples.Task.Theta_St;
for(i=0;i<n;i++) { error=0;
if(Samples.FastFTAll()!=0) { error=7; goto err;} // Быстрое преобразование Фурье для всего кристалла.
if(Samples.CompMk()!=0) { error=8; goto err;} // Вычисление матрицы Mk.
if(Samples.FindEigVectMk()!=0) { error=9; goto err;} // Нахождение собственных значений матрицы Mk.
if(CompC()!=0) { error=10; goto err;} // Вычисление матриц 'C'.
if(FindGammaExp()!=0) { error=11; goto err;} // Вычисление матриц DExp00Inv, DExp11.
if(FindRTW()!=0) { error=12; goto err;} // Нахождение амплитуд прошедших и отражённых волн.

// Прием: нужна действительная матрица. Используем только действительную часть комплексной матрицы.
if(CompPartRefTr(&R_Sp.Rows[j].Vect[i].re,&T_Sp.Rows[j].Vect[i].re)!=0) { error=13; goto err;}
err: if(error!=0) { cntErr++; R_Sp.Rows[j].Vect[i].re=T_Sp.Rows[j].Vect[i].re=A_Sp.Rows[j].Vect[i].re=0.; continue;}
A_Sp.Rows[j].Vect[i].re=1.-R_Sp.Rows[j].Vect[i].re-T_Sp.Rows[j].Vect[i].re;
*Th+=dTh; printf(".");
}
*wL+=dwL;
if(Samples.Matter.SetEps(*wL)!=0) return 14;
if(Samples.Task.SetEpsIncOut(Samples.Matter)!=0) return 15;
}

*Time=omp_get_wtime()-*Time;
printf("Completed\n"); printf("Number of errors in distribution: %d\n",cntErr);

printf("Writing spectrum in files...");
if(AnglDispOutput(nRows,n,T_Sp,R_Sp,A_Sp)!=0) return 16; // Вывод в файл всех спектров.
printf(" - Completed\n");

return 0;
}

//-----------------------------------------------------------------------------------------------------------
// Сообщение о близком к нулю значении диэлектрической проницаемости.

void MessageSmallEps(char *str)
{
printf("\nValue of Dielectric Permittivity 'eps' is too small");
if(str!=NULL) printf("\nIn the %s !\n",str); else printf(" !\n");
}

//-----------------------------------------------------------------------------------------------------------
// Сообщение о близком к нулю значении диэлектрической проницаемости.

void MessageSmallEps_(int numLay)
{
printf("\nValue of Dielectric Permittivity in the Layer %d is too small !\n",numLay);
}

//-----------------------------------------------------------------------------------------------------------
// Сообщение о близком к нулю собственном значении.

void MessageSingLam(char *str)
{
printf("\nValue of 'Kz' - z-component of Wave Vector is too small");
if(str!=NULL) printf("\nIn the %s !\n",str); else printf(" !\n");
printf("You should put Attenuation (small eps''>0) into the Substance\n");
printf("or choose another Incident Angle or Period of the Crystal.\n");
}

//-----------------------------------------------------------------------------------------------------------
// Сообщение о близком к нулю собственном значении.

void MessageSingLam_(int numLay)
{
printf("\nValue of Eigen Value in the Layer %d is too small !\n",numLay);
printf("You should put Attenuation (small eps''>0) into the Substance\n");
printf("or choose another Incident Angle or Period of the Crystal.\n");
}
