//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Программы.

BYTE InvMatr(const Matrix &A_Ini,Matrix &A_Inv); // Нахождение обратной матрицы методом исключения Гаусса с выбором столбца с главным элементом.

//-----------------------------------------------------------------------------------------------------------
// Число 1/c.

complex Inv(complex c)
{
double x,y; complex c2;

x=c.re; y=c.im; c2.re=x/(x*x+y*y); c2.im=-y/(x*x+y*y);
return c2;
}

//-----------------------------------------------------------------------------------------------------------
// Деление на комплексное число.

complex complex::operator /=(complex &a)
{
complex c;

c.re=re; c.im=im; c=c/a; re=c.re; im=c.im; return *this;
}

//-----------------------------------------------------------------------------------------------------------
// Нахождение разности фаз на концах отрезка.
// fl=0 -> вдоль оси x, fl=1 -> вдоль оси y.

BYTE Find_dPhase(complex R0,double R,BYTE fl,int nDiv,complex *P,int N,double *dPhase)
{
complex Lambda_Fun(complex *P,int N,complex L); // Функция-полином от Лямбда.

double dr,FP,FC; complex r;

if(R<=0) return 1; if(nDiv<=0) return 2; if(N<=0) return 3; if(P==NULL) return 4; if(dPhase==NULL) return 5;

r=R0; dr=R/nDiv;

FP=FC=arg(Lambda_Fun(P,N,r));
num=0; for(i=0;i<nDiv;i++) { FP=FC; if(fl==0) r.re+=dr; else r.im+=dr; FC=arg(Lambda_Fun(P,N,r)); 
if(FC-FP>M_PI) num--; if(FC-FP<-M_PI) num++;
}
FC+=num*2.*M_PI; *dPhase=FC-arg(Lambda_Fun(P,N,R0));
return 0;
}

//-----------------------------------------------------------------------------------------------------------
//
// R0 - левый верхний угол.

BYTE Div&FindNumSol(complex R0,double R,int nDiv,int *NumSol,complex *P,int N)
{
int i,j,jSh,jShN; double dr,*dPhaseX,*dPhaseY,F;
complex r;

if(R<=0.) return 1; if(nDiv<=0) return 2; if(N<=0) return 3; if(P==NULL) return 4;

err=0;
dPhaseX=NULL; dPhaseX=new double[nDiv*(nDiv+1)]; if(dPhaseX==NULL) { err=5; goto end;}
dPhaseY=NULL; dPhaseY=new double[nDiv*(nDiv+1)]; if(dPhaseY==NULL) { err=6; goto end;}
SAFE_DELETE_ARR(NumSol); NumSol=new int[nDiv*nDiv]; if(NumSol==NULL) { err=7; goto end;} 

dr=R/nDiv; // Деление квадрата.

// Заполняем массивы dPhase.
// Идем параллельно оси x.
r=R0;
for(i=0;i<nDiv+1;i++) { jSh=i*nDiv;
for(j=0;j<nDiv;j++) { if(Find_dPhase(r,dr,0,num_div_dPhase,P,N,&dPhaseX[j+jSh])!=0) { err=3; goto end;}
r.re+=dr;
}
r.re=R0.re; r.im+=dr;
}

// Идем параллельно оси y.
r=R0;
for(j=0;j<nDiv+1;j++) {
for(i=0;i<nDiv;i++) { jSh=i*nDiv;
if(Find_dPhase(r,dr,1,num_div_dPhase,&dPhaseН[j+jSh])!=0) { err=3; goto end;} r.im+=dr;
}
r.im=R0.im; r.re+=dr;
}

// Идем по всем маленьким квадратам и вычисляем число нулей внутри них.
for(i=0;i<nDiv;i++) { jSh=i*nDiv; jShN=(i+1)*nDiv; 
for(j=0;j<nDiv;j++) {
F=-dPhaseX[j+jSh]+dPhaseY[j+jSh]+dPhaseX[j+jShN]-dPhaseY[j+1+jSh]; if(F<0) { err=3; goto end;}
F/=2.*M_PI; NumSol[j+jSh]=(int)(F+0.5); if(abs((double)NumSol[j+jSh]-F)>eps_NumSol) { err=3; goto end;}
}}

end: SAFE_DELETE_ARR(dPhaseX); SAFE_DELETE_ARR(dPhaseY);
return err;
}

/*
//-----------------------------------------------------------------------------------------------------------
// Нахождение всех собственных значений с помощью принципа аргумента ТФКП.
// В L[0] уже хранится максимальный по модулю собственный вектор !.

BYTE FindLAll(complex *L,complex *P,int N)
{
BYTE err;
int ndiv,num,nSol,i,j;
double *dPhaseX,*dPhaseY;
double dx,dy,FC,FP;
complex r,r0;

if(L==NULL) return 1; if(P==NULL) return 2; if(N<=0) return 3;

ndiv=100*N; dx=2.*L[0].re/ndiv; dy=2.*L[0].im/ndiv; // Деление квадрата.

err=0;
dPhaseX=new double[ndiv]; if(dPhaseX==NULL) { err=4; goto end;}
dPhaseY=new double[ndiv]; if(dPhaseY==NULL) { err=5; goto end;}

r.re=-fabs(L[0].re); r.im=-fabs(L[0].im); r0=r;
for(i=0;i<ndiv;i++) {
FP=FC=arg(Lambda_Fun(P,N,r));
num=0; for(j=0;j<ndiv;j++) { FP=FC; r.im+=dy; FC=arg(Lambda_Fun(P,N,r)); 
if(FC-FP>M_PI) { FC-=2*M_PI; num--;} if(FC-FP<-M_PI) { FC+=2*M_PI; num++;}
}
FC+=num*2*M_PI; dPhaseX[i]=FC-arg(Lambda_Fun(P,N,r0));
r.re+=dx; r.im=-fabs(L[0].im); r0=r;
}

r.re=-fabs(L[0].re); r.im=-fabs(L[0].im); r0=r;
// Идем по отрезкам, параллельным оси x.
for(i=0;i<ndiv;i++) {
FP=FC=arg(Lambda_Fun(P,N,r));
num=0; for(j=0;j<ndiv;j++) { FP=FC; r.re+=dx; FC=arg(Lambda_Fun(P,N,r)); 
if(FC-FP>M_PI) num--; if(FC-FP<-M_PI) num++;
}
FC+=num*2*M_PI; dPhaseY[i]=FC-arg(Lambda_Fun(P,N,r0));
r.im+=dy; r.re=-fabs(L[0].re); r0=r;
}

end: SAFE_DELETE_ARR(dPhaseX); SAFE_DELETE_ARR(dPhaseY);
return err;
}
*/

//-----------------------------------------------------------------------------------------------------------
// Функция-полином от Лямбда.

complex Lambda_Fun(complex *P,int N,complex L)
{
int i; complex Res;

Res=Cmplx_0; for(i=N-1;i>=0;i--) Res+=P[i]*pow(L,N-i-1); return Res+pow(L,N);
}


clMatrix Mk; // Матрица для расчета Фурье коэффициентов компонент поля.

case p_wave_Pol:
V=((Kx*EkInv)*Kx)-E; 
for(i=0;i<N;i++) { jShMk=i*Mk.Ny; for(j=0;j<N;j++) Mk.Matr[j+jShMk]=Cmplx_0;}
for(i=0;i<V.Nx;i++) { jSh=i*V.Ny; jShMk=i*Mk.Ny; for(j=0;j<V.Ny;j++) Mk.Matr[j+jShMk+N]=V.Matr[j+jSh];}
for(i=0;i<AkInv.Nx;i++) { jSh=i*AkInv.Ny; jShMk=(i+N)*Mk.Ny;
for(j=0;j<AkInv.Ny;j++) Mk.Matr[j+jShMk]=AkInv.Matr[j+jSh];}
for(i=0;i<N;i++) { jShMk=(i+N)*Mk.Ny; for(j=0;j<N;j++) Mk.Matr[j+jShMk+N]=Cmplx_0;}
break;

// Единичная матрица.
//xxx aaa mmm
for(i=0;i<E.Nx;i++) { jSh=i*E.Ny; for(j=0;j<E.Ny;j++) if(i==j) E.Matr[j+jSh]=Cmplx_1; else E.Matr[j+jSh]=Cmplx_0;}

if(Mk.Alloc(2*N,2*N)!=0) return 2;
V=A*VI; cv=VI*VI; if(abs(cv)<1.e-24) return 13; c=(V*VI)/cv;

//d=abs(VI*VI); d/=(double)N; d=sqrt(d);
//printf("i: %d VI Norm: %lf\n",i,d); _getch();

//-----------------------------------------------------------------------------------------------------------
// Нахождение собственных значений матрицы методом Крылова.
// P - Массив коэффициентов полинома.
// L - массив собственных значений матрицы.
// V - Массив собственных векторов матриц.

BYTE Krylov(const clMatrix &A_Ini,complex *P,complex *L,clMatrix &V)
{
BYTE Gauss(const clMatrix &A_Ini,complex *F,complex *X); // Метод исключения Гаусса с выбором столбца с главным элементом.
BYTE FindMaxLambda(const clMatrix &A_Ini,complex *L,double eps); // Нахождение максимального собственного значения.

BYTE err; int i,j,k,jSh; double coe; complex *F,*Y,*YV,*q; clMatrix A;

err=0; F=NULL;
if(P==NULL||A_Ini.IsOK()!=0||L==NULL||V.IsOK()!=0) { err=1; goto end;}

// Выделение памяти.
if(A.Alloc(A_Ini.Nx,A_Ini.Nx)!=0) { err=2; goto end;} // Выделение памяти для матрицы A.
F=new complex[A_Ini.Nx]; if(F==NULL) { err=3; goto end;} // Массив правой части системы.
Y=new complex[A_Ini.Nx]; if(Y==NULL) { err=6; goto end;}
YV=new complex[A_Ini.Nx]; if(YV==NULL) { err=6; goto end;} // Вектор из алгоритма.
q=new complex[A_Ini.Nx]; if(q==NULL) { err=7; goto end;} // Массив коэффициентов для нахождения собственных векторов.

// Выбор начального вектора.
srand(1); coe=1./RAND_MAX;
for(i=0;i<A.Nx;i++) Y[i]=complex(rand()*coe,rand()*coe); 

// Заполнение матрицы A левой части системы.
for(j=A.Nx-1;j>=0;j--) { for(i=0;i<A.Nx;i++) { jSh=i*A.Ny;
if(j==A.Nx-1) A.Matr[j+jSh]=Y[i]; 
else { A.Matr[j+jSh]=0.;
for(k=0;k<A.Nx;k++) A.Matr[j+jSh]+=A_Ini.Matr[jSh+k]*Y[k]; YV[i]=A.Matr[j+jSh];
} // Конец else.
} // Конец цикла по i.
if(j!=A.Nx-1) for(k=0;k<A.Nx;k++) Y[k]=YV[k];
} // Конец цикла по j.

//printf("\n");
//printf(" Matrix A :\n");
//for(i=0;i<A.Nx;i++){ jSh=i*A.Ny; printf("\n"); for(j=0;j<A.Nx;j++) printf("%lf ",A.Matr[j+jSh]);}

// Заполнение правой части F системы.
for(i=0;i<A.Nx;i++) { jSh=i*A.Ny; F[i]=0.; for(j=0;j<A.Nx;j++) F[i]-=A_Ini.Matr[j+jSh]*Y[j];}
//printf("\n\n");
//printf(" Row F :\n\n");
//for(i=0;i<A.Nx;i++) printf("%lf ",F[i]);

// Решение системы линейных уравнений методом Гаусса для получения коэффициентов полинома.
if(Gauss(A,F,P)!=0) { err=5; goto end;}

// Вычисление собственных значений
if(FindMaxLambda(A_Ini,&L[0],0.0000001)!=0) { err=6; goto end;}
//if(FindLAll(L,P,A_Ini.Nx)!=0) { err=7; goto end;}
//printf("\nLambdas:\n");
//for(i=0;i<A_Ini.Nx;i++) printf("%lf ",L[i]);

// Получение собственных векторов из собственных значений методом Крылова.
// Тут собственные векторы хранятся в матрице V по столбцам.

q[0]=Cmplx_1;
for(j=0;j<V.Nx;j++) { // Цикл по столбцам.
for(k=1;k<V.Nx;k++) q[k]=L[j]*q[k-1]+P[k-1]; // Цикл заполнения массива коэффициентов q (из схемы Горнера).
for(i=0;i<V.Nx;i++) { jSh=i*A.Ny; V.Matr[j+jSh]=0.;  // Цикл по строкам.
for(k=0;k<V.Nx;k++) V.Matr[j+jSh]+=A.Matr[k+jSh]*q[k];
} // Конец цикла по j.
} // Конец цикла по i.

// Приведение матрицы V к естественному виду(векторы хранятся по строкам)(транспонируем).
// (Не сделано) 

// Освобождение памяти.
end: SAFE_DELETE_ARR(F); SAFE_DELETE_ARR(Y); SAFE_DELETE_ARR(YV); SAFE_DELETE_ARR(q);
return err;
}

//-----------------------------------------------------------------------------------------------------------
// Нахождение максимального собственного значения.

BYTE FindMaxLambda(const clMatrix &A_Ini,complex *L,double eps)
{
BYTE err; int i,j,jSh,k; double coe; complex *Y,*YP,L1,L1P;

if(A_Ini.IsOK()!=0) return 1; if(L==NULL) return 2; if(eps<=0.) return 3;

err=0; Y=new complex[A_Ini.Nx]; if(Y==NULL) { err=4; goto end;} 
YP=new complex[A_Ini.Nx]; if(YP==NULL) { err=5; goto end;} // Произвольный вектор.

// Выбор начального вектора.
srand(1); coe=1./RAND_MAX;
for(i=0;i<A_Ini.Nx;i++) YP[i]=complex(rand()*coe,rand()*coe);

// Алгоритм нахождения L.
k=0; while(1) { for(i=0;i<A_Ini.Nx;i++) YP[i]=Y[i];
for(i=0;i<A_Ini.Nx;i++) { jSh=i*A_Ini.Ny; Y[i]=Cmplx_0; for(j=0;j<A_Ini.Nx;j++) Y[i]+=A_Ini.Matr[j+jSh]*YP[j];}
if(k!=0) L1P=L1; L1=Y[0]/YP[0]; if(k!=0) if(fabs(abs(L1)-abs(L1P))<=eps) break; 
if(k>10000) for(i=0;i<A_Ini.Nx;i++) YP[i]=complex(rand()*coe,rand()*coe); // Выбор другого начального вектора.
if(k>30000) break;
k++;}

*L=L1;
end: SAFE_DELETE_ARR(Y); SAFE_DELETE_ARR(YP); if(err!=0) return err;
return 0;
}


//xxx
#define SigFindLam 1.e-6 // Относительная погрешность при нахождении собственного значения.
#define SigDiffVect 1.e-4 // Относительная погрешность при нахождении собственного вектора.

//-----------------------------------------------------------------------------------------------------------
// Нахождение собственных значений матрицы методом QR итераций с приведением матрицы к виду Хесенберга отражениями Хаусхолдера.

// 'EigLam' - массив собственных значений матрицы.
// 'EigVect' - Массив собственных векторов матриц в виде матрицы.

BYTE FindEigenVectors(int N,const clMatrix &A,clRowMatr &EigLam,clMatrix &EigVect)
{
BYTE GetCoeff(clRowMatr &VN,clRowMatr &VD,complex *pRes); // Нахождение коэффициента пропорциональности между столбцами.

int k,np; double d,v,vn; complex c,cp,cv,lam; clRowMatr VI,V,VV; clMatrix M;

if(N<=0) return 1; if(A.IsOK()!=0) return 2; if(A.Nx!=N) return 3; if(A.Ny!=N) return 4;
if(EigLam.IsOK()!=0) return 5; if(EigLam.N!=N) return 6;
if(EigVect.IsOK()!=0) return 7; if(EigVect.Nx!=N) return 8; if(EigVect.Ny!=N) return 9;

// Выделение памяти.
if(VI.Alloc(N)!=0) return 10;
if(V.Alloc(N)!=0) return 11;

srand(2014); // Инициализация последовательности случайных чисел.
if(VI.SetRand(1.)!=0) return 12; // Случайно инициализируем начальный вектор.

M=A;
M+=complex(5.,0.);

np=4;
M=M*M; M=M*M; // M=M*M;

for(k=0;k<300;k++) { //VI=M*VI;
//d=VI.GetAbsMax(); if(d==0.) return 13; VI/=d;
if(VI.Norm()!=0) return 13;
//xxx V=A*VI;
V=M*VI; vn=V.GetNorm(); if(vn==0.) return 14;
if(GetCoeff(V,VI,&c)!=0) return 15; if(k==0) goto e;
v=MAX(abs(c),abs(cp)); d=abs(c-cp)/v;
if(d<SigFindLam) { VI=V; if(VI.Norm()!=0) return 16; break;}
VV=V-VI*c; v=VV.GetNorm()/vn;
//xxx mmm aaa
printf("k: %d Diff: %lf Lambda (re,im): %lf %lf Diff Vect: %lf\n",k,d,real(c),imag(c),v); _getch();
e: cp=c; VI=V;}


lam=pow(c,1./(double)np);
//xxx mmm aaa
printf("lam: %lf %lf\n",real(lam),imag(lam)); _getch();

VV=A*VI;
for(k=0;k<np;k++) { v=(double)k*2.*M_PI/(double)np; c=lam*polar(1.,v)-complex(5.,0.);
V=(VI*c)-VV; d=V.GetAbsMax()/abs(c);
//xxx mmm aaa
printf("num: %d c: %lf %lf d: %lf\n",k,real(c),imag(c),d); _getch();

}


return 0;
}

//BYTE FindMaxEigValue(int N,const clMatrix &A); // Нахождение максимального собственного значения матрицы.

//if(FindMaxEigValue(N,A)!=0) return 10; // Нахождение максимального собственного значения матрицы.

//xxx
#define SigFindLam 1.e-6 // Относительная погрешность при нахождении собственного значения.
#define SigDiffVect 1.e-4 // Относительная погрешность при нахождении собственного вектора.

//-----------------------------------------------------------------------------------------------------------
// Нахождение максимального собственного значения матрицы.

BYTE FindMaxEigValue(int N,const clMatrix &A)
{
BYTE GetCoeff(clRowMatr &VN,clRowMatr &VD,complex *pRes); // Нахождение коэффициента пропорциональности между столбцами.

int k,np; double d,v,vn; complex c,cp,cv,lam,shift; clRowMatr VI,V,VV; clMatrix M;

if(N<=0) return 1; if(A.IsOK()!=0) return 2; if(A.Nx!=N) return 3; if(A.Ny!=N) return 4;

// Выделение памяти.
if(VI.Alloc(N)!=0) return 5;
if(V.Alloc(N)!=0) return 6;

srand(2014); // Инициализация последовательности случайных чисел.
if(VI.SetRand(1.)!=0) return 7; // Случайно инициализируем начальный вектор.

shift=complex(0.,0.); // 5.

M=A;
M+=shift;

np=4;
M=M*M; M=M*M; // M=M*M;

for(k=0;k<700;k++) { //VI=M*VI;
//d=VI.GetAbsMax(); if(d==0.) return 8; VI/=d;
if(VI.Norm()!=0) return 9;
//xxx V=A*VI;
V=M*VI; vn=V.GetNorm(); if(vn==0.) return 10;
if(GetCoeff(V,VI,&c)!=0) return 11; if(k==0) goto e;
v=MAX(abs(c),abs(cp)); d=abs(c-cp)/v;
VV=V-VI*c; v=VV.GetNorm()/vn;
if(d<SigFindLam||v<SigDiffVect) { VI=V; if(VI.Norm()!=0) return 12; break;}
//xxx mmm aaa
printf("k: %d Diff: %lf Lambda (re,im): %lf %lf Diff Vect: %lf\n",k,d,real(c),imag(c),v); _getch();
e: cp=c; VI=V;}


lam=pow(c,1./(double)np);
//xxx mmm aaa
printf("lam: %lf %lf\n",real(lam),imag(lam)); _getch();

VV=A*VI;
for(k=0;k<np;k++) { v=(double)k*2.*M_PI/(double)np; c=lam*polar(1.,v)-shift;
V=(VI*c)-VV; d=V.GetAbsMax()/abs(c);
//xxx mmm aaa
printf("num: %d c: %lf %lf d: %lf\n",k,real(c),imag(c),d); _getch();

}


return 0;
}

//-----------------------------------------------------------------------------------------------------------
// Нахождение коэффициента пропорциональности между столбцами.

BYTE GetCoeff(clRowMatr &VN,clRowMatr &VD,complex *pRes)
{
int i; double sw,v; complex s,rat,cN,cD;

if(VN.IsOK()!=0) return 1; if(VD.IsOK()!=0) return 2; if(VN.N!=VD.N) return 3; if(pRes==NULL) return 4;
s=Cmplx_0; sw=0.; for(i=0;i<VN.N;i++) { cN=VN.Vect[i]; cD=VD.Vect[i]; v=abs(cD); if(v<SmCnst12_d) continue;
sw+=v; rat=cN/cD; s+=rat*complex(v,0.);}
if(sw==0.) return 5; *pRes=s/sw; return 0;
}

//-----------------------------------------------------------------------------------------------------------
// Отражение Хаусхолдера для матрицы, умножение справа на 'I-uuT'.

BYTE clMatrix::Haus(int numJ)
{
int i,j,jSh,jShv,numI; double s; complex c,cs,coe;

if(IsOK()!=0) return 1; if(Nx!=Ny) return 2; numI=numJ+1; if(numI>=Nx-1) return 3; if(numJ>=Ny-1) return 4;
// Нахождение столбца 'u', хранение его в матрице начиная с 'numI', 'numJ'.
s=0.; for(i=numI;i<Nx;i++) { jSh=i*Ny+numJ; c=Matr[jSh]; s+=real(c*conj(c));}
if(s<0.) s=0.; else s=sqrt(s);
jSh=numI*Ny+numJ; if(real(Matr[jSh])<0.) s=-s; Matr[jSh]+=complex(s,0.);

// Нахождение нормы вектора 'u'.
s=0.; for(i=numI;i<Nx;i++) { jSh=i*Ny+numJ; c=Matr[jSh]; s+=real(c*conj(c));} if(s<=0.) return 5;
coe=complex(2.,0.)/s;

// Отражение Хаусхолдера. Умножение матрицы 'P' на все последующие столбцы матрицы в цикле по 'j'.
for(j=numJ+1;j<Ny;j++) {
cs=Cmplx_0; for(i=numI;i<Nx;i++) { jSh=i*Ny; cs+=Matr[jSh+numJ]*Matr[jSh+j];} cs*=coe;
for(i=numI;i<Nx;i++) { jSh=i*Ny; Matr[jSh+j]-=cs*Matr[jSh+numJ];}
} // Конец цикла по j.

// A1=PTAP; Унитарное преобразование матрицы A. Умножаем матрицу справа на матрицу 'P'='I-uuT'. 
for(i=numI;i<Nx;i++) { jSh=i*Ny;
cs=Cmplx_0; for(j=numJ+1;j<Ny;j++) { jShv=j*Ny; cs+=Matr[jShv+numJ]*Matr[jSh+j];} cs*=coe;
for(j=numJ+1;j<Ny;j++) { jShv=j*Ny; Matr[jSh+j]-=cs*Matr[jShv+numJ];}
} // Конец цикла по j.

return 0;
}
BYTE Haus(int numJ); // Отражение Хаусхолдера для матрицы, умножение справа на 'I-uuT'.

//xxx mmm aaa проверка выступающего элемента.
if(nt>numI) { jSh=nt*(N+1)-1; c=pM[jSh+N];
//xxx mmm aaa
printf("nt: %d numI,numF: %d %d Elem: %lf\n",nt,numI,numF,abs(c)); _getch();
}

//-----------------------------------------------------------------------------------------------------------
// Печать подматрицы.

BYTE PrintMatr(clMatrix &M,int numI,int numF)
{
int i,j,jSh,N; complex c;

if(M.IsOK()!=0) return 1; if(M.Nx!=M.Ny) return 2; N=M.Nx;
if(numI>=N||numF>=N) return 3; if(numI>=numF) return 4;
for(i=numI;i<=numF;i++) { jSh=i*N; printf("i: %d --> ",i);
for(j=numI;j<=numF;j++) { c=M.Matr[jSh+j]; printf("[%lf %lf] ",real(c),imag(c));}
printf("\n");}
return 0;
}

//xxx mmm aaa
//printf("--->Initial Matrix\n");
//PrintMatr(M,numI,numF); _getch();

//xxx mmm aaa
//printf("--->  nt-numI: %d numF-numI: %d\n",nt-numI,numF-numI);
//PrintMatr(M,numI,numF); _getch();


//xxx mmm aaa
BYTE PrintMatr(clMatrix &M,int numI,int numF); // Печать подматрицы.

//-----------------------------------------------------------------------------------------------------------
// Метод исключения Гаусса с выбором столбца с главным элементом.

BYTE Gauss(const clMatrix &A_Ini,complex *F,complex *X)
{
BYTE CompMaxDiff(int sz,complex *r1,complex *r2,double *MaxDif); // Расчет максимального отклонения двух столбцов.

BYTE err; int i,iEx,j,jSh,jShv,jv,jM,jp,*NumRow; double v,vM,MaxDif;
complex coe,s,*G,*F_Ch; clMatrix A;

err=0;
if(F==NULL||X==NULL||A_Ini.IsOK()!=0) { err=1; goto end;} 

G=F_Ch=NULL; NumRow=NULL; 

G=new complex[A_Ini.Nx]; if(G==NULL) { err=2; goto end;} // Массив для промежуточного столбца при решении системы методом Гаусса.
F_Ch=new complex[A_Ini.Nx]; if(F_Ch==NULL) { err=3; goto end;} // Массив, где будет результат умножения матрицы на решение X, для сравнения с исходным массивом F.
NumRow=new int[A_Ini.Nx]; if(NumRow==NULL) { err=4; goto end;} // Массив номеров подстановки столбцов при выборе наибольшего элемента.

// В матрицу A записываются результаты треугольного разложения.

A=A_Ini; // Записываем первоначальную матрицу в матрицу A, которая затем будет хранить результат треугольного разложения по методу Гаусса.

// Треугольное разложение матрицы. Результат записывается в ту же матрицу A.
for(j=0;j<A.Ny;j++) NumRow[j]=j; // Инициализация номеров подстановок столбцов.

// Цикл по строке, где проводится исключение неизвестной. Номер позиции строки - iEx.
for(iEx=0;iEx<A.Nx-1;iEx++) { jSh=iEx*A.Ny;

// Находим максимальный элемент по столбцам.
vM=0.; jp=-1; for(j=iEx;j<A.Ny;j++) {
jv=NumRow[j]; if(jv<0||jv>=A.Ny) { err=5; goto end;} // Вместо индекса j используем подстановку NumRow - jv.
v=abs(A.Matr[jSh+jv]); if(v>vM) { vM=v; jp=j;}}
if(vM==0.||jp<0) { err=6; goto end;}
jM=NumRow[jp]; if(jM<0||jM>=A.Ny) { err=7; goto end;} // Вместо индекса jp используем подстановку NumRow - jM.
jv=NumRow[iEx]; if(jM!=jv) { NumRow[iEx]=jM; NumRow[jp]=jv;} // Делаем перестановку номеров столбцов в массиве NuRow.

coe=1./A.Matr[jSh+jM];
for(j=iEx+1;j<A.Ny;j++) {
jv=NumRow[j]; if(jv<0||jv>=A.Ny) { err=8; goto end;} // Вместо индекса j используем подстановку NumRow - jv.
A.Matr[jv+jSh]*=coe;}

for(i=iEx+1;i<A.Nx;i++) { jShv=i*A.Ny; coe=A.Matr[jM+jShv];
for(j=iEx+1;j<A.Ny;j++) {
jv=NumRow[j]; if(jv<0||jv>=A.Ny) { err=9; goto end;} // Вместо индекса j используем подстановку NumRow - jv.
A.Matr[jv+jShv]-=A.Matr[jv+jSh]*coe;}
} // Конец цикла по i.
} // Конец цикла по iEx.

// Цикл вычисления столбца G.
for(i=0;i<A.Nx;i++) G[i]=F[i]; // Копируем столбец F в G.
for(iEx=0;iEx<A.Nx;iEx++) { jSh=iEx*A.Ny;
jM=NumRow[iEx]; if(jM<0||jM>=A.Ny) { err=10; goto end;} // Вместо индекса iEx используем подстановку NumRow - jM.
G[iEx]/=A.Matr[jSh+jM];
if(iEx==A.Nx-1) break;
for(i=iEx+1;i<A.Nx;i++) { jShv=i*A.Ny; G[i]-=G[iEx]*A.Matr[jM+jShv];}
} // Конец цикла по iEx.

// Цикл вычисление столбца X.
for(i=A.Nx-1;i>=0;i--) { jSh=i*A.Ny; s=0.;
jM=NumRow[i]; if(jM<0||jM>=A.Ny) { err=11; goto end;} // Вместо индекса i используем подстановку NumRow - jM.
if(i<A.Nx-1) for(j=i+1;j<A.Ny;j++) {
jv=NumRow[j]; if(jv<0||jv>=A.Ny) { err=12; goto end;} // Вместо индекса j используем подстановку NumRow - jv.
s+=A.Matr[jv+jSh]*X[jv];}
X[jM]=G[i]-s;}

// Проверяем решение - умножаем исходную матрицу на столбец решения.
for(i=0;i<A_Ini.Nx;i++) { jSh=i*A_Ini.Ny; s=0.;
for(j=0;j<A_Ini.Ny;j++) s+=A_Ini.Matr[j+jSh]*X[j]; F_Ch[i]=s;}

// Находим максимальное отклонение исходной правой части от столбца, полученного умножением матрицы системы на столбец неизвестных.
if(CompMaxDiff(A_Ini.Nx,F,F_Ch,&MaxDif)!=0) { err=13; goto end;}
printf("\n\nGauss_Meth : maximum deviation = %12.5le\n",MaxDif);

//for(i=0;i<A_Ini.Nx;i++) printf("%lf ",X[i]); printf("\n");
end: SAFE_DELETE_ARR(G); SAFE_DELETE_ARR(F_Ch); SAFE_DELETE_ARR(NumRow);
return err;
}

//-----------------------------------------------------------------------------------------------------------
// Печать диагональных элементов.

BYTE PrintDiag(clMatrix &M)
{
int i,N; complex c;

if(M.IsOK()!=0) return 1; if(M.Nx!=M.Ny) return 2; N=M.Nx;
for(i=0;i<N;i++) { c=M.Matr[i*(N+1)];
printf("n: %d   %lf %lf  %lf %lf\n",i,abs(c),arg(c)/M_PI*180.,real(c),imag(c));}
return 0;
}


CompMaxDiff,
//-----------------------------------------------------------------------------------------------------------
// Расчет максимального отклонения двух столбцов.

BYTE CompMaxDiff(int sz,complex *r1,complex *r2,double *MaxDif)
{
double v; double MD; int i;

if(sz<=0) return 1; if(r1==NULL) return 2; if(r2==NULL) return 3; if(MaxDif==NULL) return 4;
MD=fabs(abs(r1[0])-abs(r2[0])); 
for(i=0;i<sz;i++) { v=fabs(abs(r1[i])-abs(r2[i])); if(v>MD) MD=v;}
*MaxDif=MD; return 0;
}


//-----------------------------------------------------------------------------------------------------------
// Нахождение обратной матрицы методом исключения Гаусса с выбором столбца с главным элементом.

BYTE InvMatr(const clMatrix &A_Ini,clMatrix &A_Inv)
{
BYTE err; int i,iEx,j,jG,jSh,jShv,jv,jM,jp,*NumRow; double v,vM; complex coe,s;
clMatrix A,E,G;

err=0;
if(A_Ini.IsOK()!=0) { err=1; goto end;} 

NumRow=NULL; 
NumRow=new int[A_Ini.Nx]; if(NumRow==NULL) { err=4; goto end;} // Массив номеров подстановки столбцов при выборе наибольшего элемента.

// В матрицу A записываются результаты треугольного разложения.

A=A_Ini; // Записываем первоначальную матрицу в матрицу A, которая затем будет хранить результат треугольного разложения по методу Гаусса.
if(E.Alloc(A.Nx,A.Ny)!=0) { err=4; goto end;} // Единичная матрица.
if(G.Alloc(A.Nx,A.Ny)!=0) { err=5; goto end;} // результат преобразования единичной матрицы.
for(i=0;i<E.Nx;i++) { jSh=i*E.Ny; for(j=0;j<E.Ny;j++) if(i==j) E.Matr[j+jSh]=1.; else E.Matr[j+jSh]=0.;}

// Треугольное разложение матрицы. Результат записывается в ту же матрицу A.
for(j=0;j<A.Ny;j++) NumRow[j]=j; // Инициализация номеров подстановок столбцов.

// Цикл по строке, где проводится исключение неизвестной. Номер позиции строки - iEx.
for(iEx=0;iEx<A.Nx-1;iEx++) { jSh=iEx*A.Ny;

// Находим максимальный элемент по столбцам.
vM=0.; jp=-1; for(j=iEx;j<A.Ny;j++) {
jv=NumRow[j]; if(jv<0||jv>=A.Ny) { err=5; goto end;} // Вместо индекса j используем подстановку NumRow - jv.
v=abs(A.Matr[jSh+jv]); if(v>vM) { vM=v; jp=j;}}
if(vM==0.||jp<0) { err=6; goto end;}
jM=NumRow[jp]; if(jM<0||jM>=A.Ny) { err=7; goto end;} // Вместо индекса jp используем подстановку NumRow - jM.
jv=NumRow[iEx]; if(jM!=jv) { NumRow[iEx]=jM; NumRow[jp]=jv;} // Делаем перестановку номеров столбцов в массиве NuRow.

coe=1./A.Matr[jSh+jM];
for(j=iEx+1;j<A.Ny;j++) {
jv=NumRow[j]; if(jv<0||jv>=A.Ny) { err=8; goto end;} // Вместо индекса j используем подстановку NumRow - jv.
A.Matr[jv+jSh]*=coe;}

for(i=iEx+1;i<A.Nx;i++) { jShv=i*A.Ny; coe=A.Matr[jM+jShv];
for(j=iEx+1;j<A.Ny;j++) {
jv=NumRow[j]; if(jv<0||jv>=A.Ny) { err=9; goto end;} // Вместо индекса j используем подстановку NumRow - jv.
A.Matr[jv+jShv]-=A.Matr[jv+jSh]*coe;}
} // Конец цикла по i.
} // Конец цикла по iEx.

// Цикл вычисления столбцов G.
for(i=0;i<G.Nx;i++){ jSh=i*G.Ny; for(j=0;j<G.Ny;j++) G.Matr[j+jSh]=E.Matr[j+jSh];} // Копируем столбцы E в G.

for(jG=0;jG<G.Ny;jG++)
for(iEx=0;iEx<A.Nx;iEx++) { jSh=iEx*A.Ny;
jM=NumRow[iEx]; if(jM<0||jM>=A.Ny) { err=10; goto end;} // Вместо индекса iEx используем подстановку NumRow - jM.
G.Matr[iEx*G.Ny+jG]/=A.Matr[jSh+jM];
if(iEx==A.Nx-1) break;
for(i=iEx+1;i<A.Nx;i++) { jShv=i*A.Ny; G.Matr[i*G.Ny+jG]-=G.Matr[iEx*G.Ny+jG]*A.Matr[jM+jShv];}
} // Конец цикла по iEx.

// Цикл вычисление столбца самой обратной матрицы.
for(jG=0;jG<G.Ny;jG++) 
for(i=A.Nx-1;i>=0;i--) { jSh=i*A.Ny; s=0.;
jM=NumRow[i]; if(jM<0||jM>=A.Ny) { err=11; goto end;} // Вместо индекса i используем подстановку NumRow - jM.
if(i<A.Nx-1) for(j=i+1;j<A.Ny;j++) {
jv=NumRow[j]; if(jv<0||jv>=A.Ny) { err=12; goto end;} // Вместо индекса j используем подстановку NumRow - jv.
s+=A.Matr[jv+jSh]*A_Inv.Matr[jG+jv*G.Ny];}
A_Inv.Matr[jG+jM*G.Ny]=G.Matr[i*G.Ny+jG]-s;
}

end: SAFE_DELETE_ARR(NumRow);
return err;
}


//-----------------------------------------------------------------------------------------------------------
// Установка матрицы с одинаковыми элементами на диагонали.

BYTE clMatrix::SetDiag(complex Val)
{
int i,j,jSh;

if(IsOK()!=0) return 1; if(Nx!=Ny) return 2;
for(i=0;i<Nx;i++) { jSh=i*Ny; for(j=0;j<Ny;j++) if(i!=j) Matr[j+jSh]=Cmplx_0; else Matr[j+jSh]=Val;}
return 0;
}

//-----------------------------------------------------------------------------------------------------------
// Заполнение матрицы с диагональными элементами.

BYTE clMatrix::SetDiag(clRowMatr const &RM)
{
int i,j,jSh;

if(IsOK()!=0) return 1; if(Nx!=Ny) return 2; if(RM.IsOK()!=0) return 3; if(RM.N!=Nx) return 4;
for(i=0;i<Nx;i++) { jSh=i*Ny; for(j=0;j<Ny;j++) if(i!=j) Matr[j+jSh]=Cmplx_0;}
for(i=0;i<Nx;i++) { jSh=i*(Ny+1); Matr[jSh]=RM.Vect[i];} return 0;
}

//-----------------------------------------------------------------------------------------------------------
// Комплексное сопряжение матрицы.

BYTE Conj(clMatrix const &M)
{
int i,j,jSh; complex *pc;

if(M.IsOK()!=0) return 1;
for(i=0;i<M.Nx;i++) { jSh=i*M.Ny; for(j=0;j<M.Ny;j++) { pc=M.Matr+j+jSh; *pc=conj(*pc);}} return 0;
}

/*
DM=A*EV;
for(j=0;j<N;j++) { jShJ=j*N;
s=Cmplx_0; for(i=0;i<N;i++) { jShI=i*N; s+=EVH.Matr[jShJ+i]*DM.Matr[jShI+j];}
s/=M.Matr[jShJ+j];
printf("j: %d s: %lf %lf abs(s): %lf\n",j,real(s),imag(s),abs(s));
}
*/


/*
//printf("M x Eigen vectors:\n");
DM=EVH*EV;
PrintMatr(DM,0,5); // Печать подматрицы.
_getch();
*/

/*
for(ki=0;ki<10;ki++) { jShI=ki*N;
for(kj=0;kj<10;kj++) { jShJ=kj*N;
s=Cmplx_0; for(j=0;j<N;j++) { s+=conj(EVH.Matr[j+jShI])*EVH.Matr[j+jShJ];}
printf("i,j: %d %d s: %lf %lf abs(s): %lf\n",ki,kj,real(s),imag(s),abs(s));
}}
*/

//xxx mmm aaa
//if(Hermit(Q,QV)!=0) return 100;


DM=QV*Q;
printf("E matrix QV*Q ?:\n");
PrintMatr(DM,0,5); // Печать подматрицы.
_getch();


//printf("--->Matrix A\n");
//FindMaxEigValue(N,A); // Нахождение максимального собственного значения матрицы.

//printf("--->Matrix M\n");
//FindMaxEigValue(N,M); // Нахождение максимального собственного значения матрицы.
//_getch();

//xxx mmm aaa
Q=Cmplx_1; MV=M;

//xxx mmm aaa
// Test.
j=numJ; cs=Cmplx_0;
for(i=numI;i<N;i++) { jSh=i*N; c=pM[jSh+j]; if(i==numI) c+=alp; cs+=conj(pM[jSh+numJ])*c;} cs*=coe;
for(i=numI;i<N;i++) { jSh=i*N; c=pM[jSh+j]; if(i==numI) c+=alp; c-=cs*pM[jSh+numJ];
printf("Test numJ,i: %d %d s: %lf %lf abs(s): %lf\n",numJ,i,real(c),imag(c),abs(c));}
_getch();

//xxx mmm aaa
//printf("QR iterations - numI,numF: %d %d\n",numI,numF); // _getch();

//xxx mmm aaa
//printf("QR iterations succeed - numI,numF: %d %d cnt: %d\n",numI,numF,cnt); // _getch();

//xxx mmm aaa
//printf("All OK\n"); _getch();
//PrintDiag(M);

//xxx aaa mmm
if(Hermit(Q,QV)!=0) return 100;

//xxx mmm aaa
/*
DM=Q*QV;
printf("E matrix Q*QV after Hessen ?:\n");
PrintMatr(DM,0,5,0,5); // Печать подматрицы.
_getch();
*/

//xxx mmm aaa
printf("Hessen(A):\n");
PrintMatr(M,0,5,0,5); // Печать подматрицы.
_getch();


M=Q*A*QV; // Должна быть та же матрица Хессенберга.

//xxx mmm aaa
printf("Hessen Q*A*QV:\n");
PrintMatr(M,0,5,0,5); // Печать подматрицы.
_getch();

//xxx mmm aaa
//Q=Cmplx_1; MV=M;

//xxx mmm aaa
DM=QV*Q;
printf("E matrix QV*Q ?:\n");
PrintMatr(DM,0,5,0,5); // Печать подматрицы.
_getch();


for(j=0;j<GP.Ny;j++) { GP.Matr[j+jSh]=GC.Matr[j+jSh]=Cmplx_0; 
if(i==j) ZP.Matr[j+jSh]=ZC.Matr[j+jSh]=Cmplx_1; else ZP.Matr[j+jSh]=ZC.Matr[j+jSh]=Cmplx_0;}
}

 k[1]=k0*cos(Task.Theta*M_PI/180.);
int ND; // ND=N*2 - размер матриц 'Mk', 'EigVectMk' и число собственных значений в массиве 'EigLamMk'.
for(i=0;i<ND;i++) { jSh=i*ND; for(j=0;j<ND;j++) Mk.Matr[j+jSh]=Cmplx_0;}

// Инициализация матриц EigLam,EigVect.
for(i=0;i<ND;i++) EigLam.Vect[i]=EigLamMk.Vect[i];
for(i=0;i<EigVect.Nx;i++) { jSh=i*EigVect.Ny;
for(j=0;j<EigVect.Ny;j++) EigVect.Matr[j+jSh]=EigVectMk.Matr[j+jSh];

// Вычисление матрицы Kx.
//xxx Убрать матрицу Kx - сделать через столбец.

K=2.*M_PI/LayGeom.Period;
for(i=0;i<N;i++) { jSh=i*N;
for(j=0;j<N;j++) if(i==j) Kx.Matr[j+jSh]=complex(((double)(-L1+i)*K+k0X)/k0,0.);
else Kx.Matr[j+jSh]=Cmplx_0;}

// Инициализация матриц EigLam,EigVect.
EigLam=EigLamMk; EigVect=EigVectMk;

//xxx mmm aaa
for(i=0;i<ND;i++) EigLam.Vect[i]=EigLamMk.Vect[i];
for(i=0;i<EigVect.Nx;i++) { jSh=i*EigVect.Ny;
for(j=0;j<EigVect.Ny;j++) EigVect.Matr[j+jSh]=EigVectMk.Matr[j+jSh];
}


clMatrix Gamma00,Gamma01,Gamma10,Gamma11; // Блоки матрицы Gamma=EigVectMk*DiagExp*!EigVectMk.

//-----------------------------------------------------------------------------------------------------------
// Нахождение блоков матриц 'Gamma'.

BYTE strLay::FindGamma(void)
{
int i,j,jShB,jSh0,jSh1,ND; clRowMatr DiagExp; clMatrix Gamma; 

if(N<=0) return 1; ND=N*2;
if(EigLamMk.IsOK()!=0) return 2; if(EigLamMk.N!=ND) return 3;
if(EigVectMk.IsOK()!=0) return 4; if(EigVectMk.Nx!=ND||EigVectMk.Ny!=ND) return 5;

// Выделение памяти.
if(DiagExp.Alloc(ND)!=0) return 6;
if(Gamma.Alloc(ND,ND)!=0) return 7;
if(Gamma00.Alloc(N,N)!=0) return 8;
if(Gamma01.Alloc(N,N)!=0) return 9;
if(Gamma10.Alloc(N,N)!=0) return 10;
if(Gamma11.Alloc(N,N)!=0) return 11;

for(i=0;i<ND;i++) DiagExp.Vect[i]=exp(-EigLamMk.Vect[i]*depth); // Нахождение диагональной матрицы, составленной из экспонент.

Gamma=(EigVectMk%DiagExp)*!EigVectMk; // Нахождение матрицы Gamma.

//xxx mmm aaa

// Нахождение четырех блоков матрицы Gamma.
for(i=0;i<N;i++) { jShB=i*N; jSh0=i*ND; jSh1=(i+N)*ND; 
for(j=0;j<N;j++) {
Gamma00.Matr[j+jShB]=Gamma.Matr[j+jSh0]; Gamma01.Matr[j+jShB]=Gamma.Matr[j+N+jSh0];
Gamma10.Matr[j+jShB]=Gamma.Matr[j+jSh1]; Gamma11.Matr[j+jShB]=Gamma.Matr[j+N+jSh1];
}}
return 0;
}

// Расчёт G0 по реккурентной формуле. Идём по слоям.
for(i=laysNum-1;i>=0;i--) {
GP=GC; ZP=ZC; pL=&Lays[i];
GC=(pL->Gamma10+pL->Gamma11*GP)*!(pL->Gamma00+pL->Gamma01*GP);
ZC=ZP*!(pL->Gamma00+pL->Gamma01*GP);
}
// G0=GC.

//xxx mmm aaa
printf("StepQRShift err: nt,nimI,numF: %d %d %d c: %lf %lf s: %lf %lf\n",
nt,numI,numF,real(c),imag(c),real(s),imag(s));
if(nt==numI) {
//xxx mmm aaa
printf("pM[jSh]: %lf %lf sig: %lf %lf\n",nt,numI,numF,real(pM[jSh]),imag(pM[jSh]),real(sig),imag(sig));
}
getch();

// Разделение на волны, распространяющиеся в противоположные стороны.
s=0; for(i=0;i<EigLam.N;i++) { if(EigLam.Vect[Num[i]].im>0) {
if(s>=EigLam.N/2) { err=10; goto end;}
for(k=i;k>s;k--) { v=Num[k-1]; Num[k-1]=Num[k]; Num[k]=v;} // Пузырьковое смещение.
s++;}
}

//-----------------------------------------------------------------------------------------------------------
// Проверка точности алгоритма нахождения собственных векторов, печать диагонализованной матрицы.

BYTE CheckEigVect(FILE *file,clMatrix &A,clMatrix &P,BYTE numDiag,BYTE flAbs)
{
BYTE InvMatr(const clMatrix &M,clMatrix &Inv); // Нахождение обратной матрицы методом исключения Гаусса с выбором столбца с главным элементом.
BYTE PrintMatr(FILE *file,clMatrix &M,int iI,int iF,int jI,int jF); // Печать подматрицы.
BYTE PrintMatrDiag(FILE *file,clMatrix &M,BYTE numDiag,BYTE flAbs); // Печать диагоналей матрицы.

int N; complex s; clMatrix PInv,DM; char *Form;

if(A.IsOK()!=0) return 1; if(A.Nx!=A.Ny) return 2; N=A.Nx;
if(P.IsOK()!=0) return 3; if(P.Nx!=P.Ny) return 4; if(P.Nx!=N) return 5;
if(PInv.Alloc(N,N)!=0) return 6;

// Проверяем правильность алгоритма. Умножаем матрицу M слева на PInv, справа на P. Должны получить диагональную.
Form="Diagonal matrix PInv*A*P ?:\n"; if(file==NULL) printf(Form); else fprintf(file,Form);
if(InvMatr(P,PInv)!=0) return 7;
DM=PInv*A*P;
//xxx mmm aaa
PrintMatrDiag(file,DM,3,1); // Печать диагоналей матрицы.

//Form="numI,numF: %d %d\n"; if(file==NULL) printf(Form,numI,numF); else fprintf(file,Form,numI,numF);
//PrintMatr(file,DM,numI,numF,numI,numF); // Печать подматрицы.
if(file==NULL) _getch(); return 0;
}

//xxx mmm aaa
BYTE PrintMatrDiag(FILE *file,clMatrix &M,BYTE numDiag,BYTE flAbs); // Печать диагоналей матрицы.
//xxx mmm aaa
clMatrix DM;
FILE *fp;

fp=fopen("ChUnitaryEigVect.out","w"); if(fp!=NULL) {
DM=QH*Q;
fprintf(fp,"Unity matrix QH*Q ?:\n");
PrintMatrDiag(fp,DM,3,1); // Печать диагоналей матрицы.
SAFE_CLOSE(fp);}

//-----------------------------------------------------------------------------------------------------------
// Треугольное разложение матрицы в методе исключения Гаусса с выбором столбца с главным элементом.

BYTE TriDecRow(int N,const clMatrix &A,clMatrix &T,int *NumRow)
{
int i,iEx,j,jSh,jShv,jv,jM,jp; double v,vM; complex coe,*pM;

if(N<=0) return 1;
if(A.IsOK()!=0) return 2; if(A.Nx!=N) return 3; if(A.Ny!=N) return 4;
if(T.IsOK()!=0) return 5; if(T.Nx!=N) return 6; if(T.Ny!=N) return 7;
if(NumRow==NULL) return 8;

T=A; pM=T.Matr; // Записываем первоначальную матрицу в матрицу 'T', которая затем будет хранить результат треугольного разложения по методу Гаусса.

// Треугольное разложение матрицы. Результат записывается в ту же матрицу 'T'.
for(j=0;j<N;j++) NumRow[j]=j; // Инициализация номеров подстановок столбцов.

// Цикл по строке, где проводится исключение неизвестной. Номер позиции строки - 'iEx'.
for(iEx=0;iEx<N-1;iEx++) { jSh=iEx*N;

// Находим максимальный элемент по столбцам.
vM=0.; jp=-1; for(j=iEx;j<N;j++) {
jv=NumRow[j]; if(jv<0||jv>=N) return 9; // Вместо индекса 'j' используем подстановку 'NumRow' - 'jv'.
v=abs(pM[jSh+jv]); if(v>vM) { vM=v; jp=j;}}
if(vM==0.||jp<0) return 10;
jM=NumRow[jp]; if(jM<0||jM>=N) return 11; // Вместо индекса 'jp' используем подстановку 'NumRow' - 'jM'.
jv=NumRow[iEx]; if(jM!=jv) { NumRow[iEx]=jM; NumRow[jp]=jv;} // Делаем перестановку номеров столбцов в массиве 'NumRow'.

coe=Cmplx_1/pM[jSh+jM];
for(j=iEx+1;j<N;j++) {
jv=NumRow[j]; if(jv<0||jv>=N) return 12; // Вместо индекса 'j' используем подстановку 'NumRow' - 'jv'.
pM[jv+jSh]*=coe;}

for(i=iEx+1;i<N;i++) { jShv=i*N; coe=pM[jM+jShv];
for(j=iEx+1;j<N;j++) {
jv=NumRow[j]; if(jv<0||jv>=N) return 13; // Вместо индекса 'j' используем подстановку 'NumRow' - 'jv'.
pM[jv+jShv]-=pM[jv+jSh]*coe;}
} // Конец цикла по 'i'.
} // Конец цикла по 'iEx'.

return 0;
}

//-----------------------------------------------------------------------------------------------------------
// Решение системы линейных уравнений по треугольному разложению матрицы в методе исключения Гаусса с выбором столбца с главным элементом.

BYTE TriSolvRow(int N,const clMatrix &T,int *NumRow,complex *F,complex *X)
{
BYTE err; int i,iEx,j,jv,jSh,jShv,jM; complex s,*G,*pM;

if(N<=0) return 1; if(T.IsOK()!=0) return 2; if(T.Nx!=N) return 3; if(T.Ny!=N) return 4; pM=T.Matr;
if(NumRow==NULL) return 5; if(F==NULL) return 6; if(X==NULL) return 7;

G=NULL; err=0;
G=new complex[N]; if(G==NULL) { err=8; goto end;} // Массив для промежуточного столбца при решении системы методом Гаусса.

// Цикл вычисления столбца 'G'.
for(i=0;i<N;i++) G[i]=F[i]; // Копируем столбец 'F' в 'G'.
for(iEx=0;iEx<N;iEx++) { jSh=iEx*N;
jM=NumRow[iEx]; if(jM<0||jM>=N) { err=9; goto end;} // Вместо индекса 'iEx' используем подстановку 'NumRow' - 'jM'.
G[iEx]/=pM[jSh+jM]; if(iEx==N-1) break;
for(i=iEx+1;i<N;i++) { jShv=i*N; G[i]-=G[iEx]*pM[jM+jShv];}
} // Конец цикла по 'iEx'.

// Цикл вычисления столбца 'X'.
for(i=N-1;i>=0;i--) { jSh=i*N; s=0.;
jM=NumRow[i]; if(jM<0||jM>=N) { err=10; goto end;} // Вместо индекса 'i' используем подстановку 'NumRow' - 'jM'.
if(i<N-1) for(j=i+1;j<N;j++) {
jv=NumRow[j]; if(jv<0||jv>=N) { err=11; goto end;} // Вместо индекса 'j' используем подстановку 'NumRow' - 'jv'.
s+=pM[jv+jSh]*X[jv];}
X[jM]=G[i]-s;}

end: SAFE_DELETE_ARR(G); return err;
}

//-----------------------------------------------------------------------------------------------------------
// Треугольное разложение матрицы в методе исключения Гаусса с выбором строки с главным элементом.

BYTE TriDecStr(int N,const clMatrix &A,clMatrix &T,int *NumStr)
{
int i,iEx,j,jSh,jShv,iv,iM,ip; double v,vM; complex coe,*pM;

if(N<=0) return 1;
if(A.IsOK()!=0) return 2; if(A.Nx!=N) return 3; if(A.Ny!=N) return 4;
if(T.IsOK()!=0) return 5; if(T.Nx!=N) return 6; if(T.Ny!=N) return 7;
if(NumStr==NULL) return 8;

T=A; pM=T.Matr; // Записываем первоначальную матрицу в матрицу 'T', которая затем будет хранить результат треугольного разложения по методу Гаусса.

// Треугольное разложение матрицы. Результат записывается в ту же матрицу 'T'.
for(j=0;j<N;j++) NumStr[j]=j; // Инициализация номеров подстановок столбцов.

// Цикл по строке, где проводится исключение неизвестной. Номер позиции строки без учёта перестановок - iEx.
for(iEx=0;iEx<N-1;iEx++) {

// Находим максимальный элемент по строкам.
vM=0.; ip=-1; for(i=iEx;i<N;i++) {
iv=NumStr[i]; if(iv<0||iv>=N) return 9; // Вместо индекса i используем подстановку NumStr - iv.
jShv=iv*N; v=abs(pM[jShv+iEx]); if(v>vM) { vM=v; ip=i;}}
if(vM==0.||ip<0) return 10;
iM=NumStr[ip]; if(iM<0||iM>=N) return 11; // Вместо индекса ip используем подстановку NumStr - iM.
iv=NumStr[iEx]; if(iM!=iv) { NumStr[iEx]=iM; NumStr[ip]=iv;} // Делаем перестановку номеров строк в массиве NumStr.

jSh=iM*N; coe=1./pM[jSh+iEx];
for(j=iEx+1;j<N;j++) pM[j+jSh]*=coe;

for(i=iEx+1;i<N;i++) {
iv=NumStr[i]; if(iv<0||iv>=N) return 12; // Вместо индекса i используем подстановку NumStr - iv.
jShv=iv*N; coe=pM[iEx+jShv];
for(j=iEx+1;j<N;j++) pM[j+jShv]-=pM[j+jSh]*coe;
} // Конец цикла по i.
} // Конец цикла по iEx.

return 0;
}

//-----------------------------------------------------------------------------------------------------------
// Решение системы линейных уравнений по треугольному разложению матрицы в методе исключения Гаусса с выбором строки с главным элементом.

BYTE TriSolvStr(int N,const clMatrix &T,int *NumStr,complex *F,complex *X)
{
BYTE err; int i,iEx,iv,j,jSh,jShv,iM; complex s,*G,*pM;

if(N<=0) return 1; if(T.IsOK()!=0) return 2; if(T.Nx!=N) return 3; if(T.Ny!=N) return 4; pM=T.Matr;
if(NumStr==NULL) return 5; if(F==NULL) return 6; if(X==NULL) return 7;

G=NULL; err=0;
G=new complex[N]; if(G==NULL) { err=8; goto end;} // Массив для промежуточного столбца при решении системы методом Гаусса.

// Цикл вычисления столбца G.
for(i=0;i<N;i++) G[i]=F[i]; // Копируем столбец F в G.
for(iEx=0;iEx<N;iEx++) {
iM=NumStr[iEx]; if(iM<0||iM>=N) { err=9; goto end;} // Вместо индекса iEx используем подстановку NumStr - iM.
jSh=iM*N; G[iM]/=pM[jSh+iEx]; if(iEx==N-1) break;
for(i=iEx+1;i<N;i++) {
iv=NumStr[i]; if(iv<0||iv>=N) { err=10; goto end;} // Вместо индекса iEx используем подстановку NumStr - iv.
jShv=iv*N; G[iv]-=G[iM]*pM[iEx+jShv];}
} // Конец цикла по iEx.

// Цикл вычисление столбца X.
for(i=N-1;i>=0;i--) {
iv=NumStr[i]; if(iv<0||iv>=N) { err=11; goto end;} // Вместо индекса i используем подстановку NumStr - iv.
jSh=iv*N; s=0.;
if(i<N-1) for(j=i+1;j<N;j++) s+=pM[j+jSh]*X[j];
X[i]=G[iv]-s;}

end: SAFE_DELETE_ARR(G); return err;
}

//-----------------------------------------------------------------------------------------------------------
// Сортировка элементов комплексного массива методом пузырька с использованием массива подстановок.
// Сортировка комплексного массива собственных значений методом пузырька с использованием массива подстановок.

BYTE SortRowMatr(int N,clRowMatr Arr,int *Num)
{
int i,j,v;

if(N==0) return 0; if(N<0) return 1; if(Arr.IsOK()!=0) return 2; if(Num==NULL) return 3;
for(j=0;j<N-1;j++) { 
for(i=0;i<N-j-1;i++) if(abs(Arr.Vect[Num[i]])>=abs(Arr.Vect[Num[i+1]])) { v=Num[i]; Num[i]=Num[i+1]; Num[i+1]=v;}
}
return 0;
}

//xxx mmm aaa
// Разделение на волны, распространяющиеся в противоположные стороны.
s=0; for(i=0;i<EigLam.N;i++) if(EigLam.Vect[Num[i]].im>0) {
if(s>=EigLam.N/2) { err=10; goto end;}
for(k=i;k>s;k--) { v=Num[k-1]; Num[k-1]=Num[k]; Num[k]=v;} // Пузырьковое смещение.
s++;}


//xxx mmm aaa
//-----------------------------------------------------------------------------------------------------------
// Сортировка комплексного массива собственных значений методом пузырька с использованием массива подстановок.

BYTE SortEigVal_(int N,clRowMatr Arr,int *Num)
{
int i,j,v;

if(N==0) return 0; if(N<0) return 1; if(Arr.IsOK()!=0) return 2; if(Num==NULL) return 3;
for(j=0;j<N-1;j++) { 
for(i=0;i<N-j-1;i++) if(real(Arr.Vect[Num[i]])>=real(Arr.Vect[Num[i+1]])) { v=Num[i]; Num[i]=Num[i+1]; Num[i+1]=v;}}
return 0;
}

//xxx mmm aaa
//-----------------------------------------------------------------------------------------------------------
// Сортировка комплексного массива собственных значений методом пузырька с использованием массива подстановок.

BYTE SortEigVal__(int N,clRowMatr Arr,int *Num)
{
int i,j,v;
//xxx mmm aaa
int s,k;

if(N==0) return 0; if(N<0) return 1; if(Arr.IsOK()!=0) return 2; if(Num==NULL) return 3;
for(j=0;j<N-1;j++) { 
for(i=0;i<N-j-1;i++) if(abs(Arr.Vect[Num[i]])>=abs(Arr.Vect[Num[i+1]])) { v=Num[i]; Num[i]=Num[i+1]; Num[i+1]=v;}}

// Сортировка по действительной части.
s=0; for(i=0;i<Arr.N;i++) if(Arr.Vect[Num[i]].re<0) {
if(s>=Arr.N/2) return 4;
for(k=i;k>s;k--) { v=Num[k-1]; Num[k-1]=Num[k]; Num[k]=v;} // Пузырьковое смещение.
s++;}

return 0;
}


/*
complex operator *(const clRowMatr &,const clRowMatr &); // Умножение элементов столбцов.
complex operator &(const clRowMatr &,const clRowMatr &); // Скалярное произведение векторов столбцов.

//-----------------------------------------------------------------------------------------------------------
// Умножение элементов столбцов.

complex operator *(const clRowMatr &RM1,const clRowMatr &RM2)
{
int i,N; complex s;

s=Cmplx_0; if(RM1.IsOK()!=0) return s; if(RM2.IsOK()!=0) return s;
N=MIN(RM1.N,RM2.N); for(i=0;i<N;i++) s+=RM1.Vect[i]*RM2.Vect[i]; return s;
}

//-----------------------------------------------------------------------------------------------------------
// Скалярное произведение векторов столбцов.

complex operator &(const clRowMatr &RM1,const clRowMatr &RM2)
{
int i,N; complex s;

s=Cmplx_0; if(RM1.IsOK()!=0) return s; if(RM2.IsOK()!=0) return s;
N=MIN(RM1.N,RM2.N); for(i=0;i<N;i++) s+=RM1.Vect[i]*conj(RM2.Vect[i]); return s;
}
*/

n=Task.L1+Task.nInc; if(n<0||n>=J.N) return 10;

// Проверки и расчёт nInc.
if(Task.wLength<=0.) return 1;
Task.nInc=(int)(Period*2.*M_PI*sin(Task.Theta*M_PI/180.)/Task.wLength);
if(Task.nInc>=Task.M1||Task.nInc<-Task.L1) return 2;
v=Task.nInc*Task.wLength/(Period*2.*M_PI); Task.Theta=asin(v)*180./M_PI;
printf("Theta corrected: %lf nInc: %d\n",Task.Theta,Task.nInc);

int nInc; // Номер гармоники для падающей волны.

//xxx mmm aaa

fl=1; for(j=0;j<Matter.N;j++) if(c==Matter.Mat[j].chr) {  fl=0; break;} if(fl==1) return 5;
 v=(double)(i-Task.L1)*coe;

//xxx mmm aaa
BYTE PrintMatr(FILE *file,clMatrix &M,int iI,int iF,int jI,int jF); // Печать подматрицы.
//..........................................................
FILE *fp;
if(flHomogen!=0) {
fp=fopen("MatricesEps.out","w"); if(fp==NULL) return 100;
fprintf(fp,"Ek Matrix -------->");
PrintMatr(fp,Ek,0,6,0,6); // Печать подматрицы.
fprintf(fp,"Ak Matrix -------->");
PrintMatr(fp,Ak,0,6,0,6); // Печать подматрицы.
SAFE_CLOSE(fp);
}
//..........................................................

// Нахождение обратной матрицы 'P^-1'.
if(InvMatr(EigVectMk,EigInv)!=0) return 10;
//xxx mmm aaa
if(EigInv.Alloc(ND,ND)!=0) return 9; // Матрица для обратной матрицы 'EigVectMk'.

if(EigLam.Alloc(ND)!=0) { err=8; goto end;}
if(EigVect.Alloc(ND,ND)!=0) { err=9; goto end;}
if(EigVectInv.Alloc(ND,ND)!=0) { err=10; goto end;}

//xxx mmm aaa
/*
//..........................................................
//xxx mmm aaa
BYTE PrintEigVal(FILE *file,clRowMatr &EigLam); // Печать собственных значений.

FILE *fp;

fp=fopen("FindEigVectMkHomogen.out","w"); if(fp==NULL) return 100;
fprintf(fp,"Eigen Values -------->");
PrintEigVal(fp,EigLamMk); // Печать собственных значений.
SAFE_CLOSE(fp);
//..........................................................
*/

else { if(pL->CompMkHomogen(Task)!=0) return 3;}
//-----------------------------------------------------------------------------------------------------------
// Вычисление матрицы Mk для расчёта Фурье коэффициентов компонент поля (случай однородного слоя).

BYTE strLay::CompMkHomogen(const clTask &Task)
{
int i; double dK,k0X,k0,v; clRowMatr Kx2;

if(flHomogen==0) return 1; if(N<=0) return 2;
if(Task.wLength<=0.) return 3; if(LayGeom.Period<=0.) return 4;

// Расчёт волнового вектора и его компоненты по оси X.
k0=2.*M_PI/Task.wLength; k0X=k0*sin(Task.Theta*M_PI/180.); dK=2.*M_PI/LayGeom.Period;

// Выделение памяти для вспомогательных матриц.
if(Kx2.Alloc(N)!=0) return 5;

// Вычисление матрицы Kx2.
for(i=0;i<N;i++) { v=((double)(-L1+i)*dK+k0X)/k0; Kx2.Vect[i]=complex(v*v,0.);}

switch(Task.flPol) { default: return 6;
// Вектор E в плоскости падения.
case p_wave_Pol: MkDiag=Kx2/EpsHomogen-Cmplx_1; break;
// Вектор E перпендикулярен плоскости падения.
case s_wave_Pol: MkDiag=Kx2-EpsHomogen; break;
}

return 0;
}

BYTE CompMkHomogen(const clTask &Task); // Вычисление матрицы Mk для расчёта Фурье коэффициентов компонент поля (случай однородного слоя).
if(MkDiag.IsOK()!=0) return 5; if(MkDiag.N!=N) return 6;
clRowMatr MkDiag; // Диагональная матрица в матрице операторов уравнения Максвелла для однородного слоя.

//-----------------------------------------------------------------------------------------------------------
// Нахождение распределения поля на границах слоёв.

BYTE clCrystal::FindBoundsFieldDistr(void)
{
int i,N; struct strLay *pL; clMatrix G_T,Z_T;

if(IsOK()!=0) return 1; if(T.IsOK()!=0) return 2; // Проверки.
N=Lays[0].N;

// Нулевая итерация.

//xxx mmm aaa
// Цикл по слоям.
//for(i=laysNum-1;i>=0;i--) { pL=Lays+i;

return 0;
}

//-----------------------------------------------------------------------------------------------------------
// Нахождение распределения поля внутри кристалла.

BYTE clCrystal::FindFullFieldDistr(void)
{
BYTE be;

if(IsOK()!=0) return 1; if(T.IsOK()!=0) return 2; BYTE be;
be=FindBoundsFieldDistr(); if(be!=0) { printf("\nbe=%d\n",be); return 3;}
//xxx mmm aaa
return 0;
}

clRowMatr EFour,HFour; // Столбцы коэффициентов Фурье разложения полей Е и H на передней границе первого слоя, куда падает волна.
clMatrix GKInc; // Матрица 'Gk' на передней границе первого слоя, куда падает волна.
GKInc=*pGP; // Запоминаем матрицу 'Gk' на передней границе первого слоя, куда падает волна.


pL->GK[3]=*pGP; // Запоминаем последнюю матрицу 'Gk' на выходе слоя.
clMatrix GK[4]; // Матрицы 'Gk' для расчёта полей внутри кристалла (соответственно перед итерациями с 'P^-1','DExp','P' и после итераций).
if(T.Alloc(N)!=0) return 6;
if(J.Alloc(N)!=0) return 7;

//xxx mmm aaa
/*
BYTE PrintRowMatr(FILE *file,clRowMatr &RM,int iI,int iF); // Печать части столбца.
for(i=0;i<laysNum;i++) { pL=Lays+i;
for(k=0;k<4;k++) {
printf("i,k: %d %d EFields----->\n",i,k); pP=pL->EFour+k; PrintRowMatr(NULL,*pP,0,N-1); _getch();
printf("i,k: %d %d HFields----->\n",i,k); pP=pL->HFour+k; PrintRowMatr(NULL,*pP,0,N-1); _getch();
}}
*/
// Для 'p' волны EDistrFour=ExFour,HDistrFour=HyFour; Для 's' волны EDistrFour=EyFour,HDistrFour=HxFour;
clRows EDistr,HDistr,ZDistr; // Распределение компонент полей в кристалле (в зависимости от типа поляризации). 
// Для 'p' волны EDistr=Ex,HDistr=Hy,ZDistr=Ez; Для 's' волны EDistr=Ey,HDistr=Hx,ZDistr=Hz;

// Выделение памяти.
if(DExp00Inv.Alloc(N)!=0) return 3; // Диагональная матрица, обратная к 'Gamma00'.
if(DExp11.Alloc(N)!=0) return 4; // Диагональная матрица 'Gamma11'.



//-----------------------------------------------------------------------------------------------------------
// Нахождение коэффициентов Фурье разложения полей Е и H в слое.

BYTE strLay::FindDistrFourEH(const strLaySmpl &Smpl)
{
int i,j; double z; clRowMatr X,Y,*pRH,*pRE; 

if(nDiv<1) return 1; if(dz<=0.) return 2; if(IsOK()!=0) return 3; if(IsOKGammaExp()!=0) return 4;
if(Smpl.EigLamMk.IsOK()!=0) return 5; if(Smpl.EigLamMk.N!=N*2) return 6;
//xxx mmm aaa
//if(Smpl.EigVectMk.IsOK()!=0) return 5;

// Вычисляем столбец 'Y' с конца слоя к началу, а столбец 'X' с начала слоя к концу, для избежания роста амплитуд собственных волн.
for(j=1;j<nDiv;j++) { z=dz*(double)j;
// Нахождение диагональных матриц, составленных из экспонент.
for(i=0;i<N;i++) {
DExp00Inv.Vect[i]=exp(Smpl.EigLamMk.Vect[i]*z);
DExp11.Vect[i]=exp(-Smpl.EigLamMk.Vect[i+N]*z);}

// Находим коэффициенты Фурье разложения полей Е и H.
X=DExp00Inv*HFour[2]; Y=DExp11*HFour[1];
pRE=EDistrFour.Get(j-1); *pRE=(Smpl.GamP[0][0]*X)+(Smpl.GamP[0][1]*Y);
pRH=HDistrFour.Get(nDiv-j-1); *pRH=(Smpl.GamP[1][0]*X)+(Smpl.GamP[1][1]*Y);
} // Конец цикла по 'j'.

return 0;
}


//-----------------------------------------------------------------------------------------------------------
// Нахождение распределения полей Е и H в слое (выполнение обратного преобразования Фурье для нужных компонент полей).

BYTE strLay::FindDistrEH(const strLaySmpl &Smpl,clTask &Task)
{
BYTE GetN(BYTE M,int *N); // Нахождение 'N' - 2 в степени 'M'.
BYTE FastFT(complex *A,complex *B,BYTE M,SCHAR dir); // Быстрое преобразование Фурье.
double Log2(double n); // Логарифм по основанию 2.
BYTE FillA(int NF,complex *A,const clRowMatr &Row,int L1,int M1); // Заполнение матрицы A для обратного преобразования Фурье.

BYTE err; int NF,j,i,M; double dK,k0,k0X,kx; complex *A; clRowMatr *pR,ZFour;

if(N<=0) return 1; if(Smpl.LayGeom.Period<=0.) return 2; if(Task.wLength<=0.) return 3;
if(EDistrFour.IsOK_All()!=0) return 2; if(EDistrFour.N!=N) return 3;
if(HDistrFour.IsOK_All()!=0) return 4; if(HDistrFour.N!=N) return 5;

A=NULL; err=0;

// Вспомогательные величины.
dK=2.*M_PI/Smpl.LayGeom.Period; k0=2.*M_PI/Task.wLength; k0X=k0*sin(Task.Theta*M_PI/180.);

// Установка числа гармоник.
M=(int)Log2(N)+1; if(GetN(M,&NF)!=0) return 4; if(NF<N) return 5;

// Выделение памяти. 

A=new complex[NF]; if(A==NULL) { err=6; goto end;} // Матрица, используемая для обратного преобразования Фурье.
if(EDistr.Alloc(NF,nDiv+1)!=0) { err=7; goto end;}
if(HDistr.Alloc(NF,nDiv+1)!=0) { err=8; goto end;}
if(ZDistr.Alloc(NF,nDiv+1)!=0) { err=9; goto end;} // Массивы распределений.
if(ZFour.Alloc(N)!=0) { err=10; goto end;} // Вспомогательный массив.

//---------------------------------------------------------------------------------------
// Цикл по слоям для расчета обратного преобразования Фурье полей.
// Для 'p' волны EDistr=Ex,HDistr=Hy,ZDistr=Ez; Для 's' волны EDistr=Ey,HDistr=Hx,ZDistr=Hz;

// Отдельный просчёт на первой границе слоя.
pR=EDistr.Rows;
if(FillA(NF,A,EFour[3],Smpl.L1,Smpl.M1)!=0) { err=11; goto end;} // Заполнение матрицы 'A' для обратного преобразования Фурье.
if(FastFT(A,pR->Vect,M,-1)!=0) { err=12; goto end;} // Обратное преобразование Фурье для полей 'Е'.
pR=HDistr.Rows;
if(FillA(NF,A,HFour[3],Smpl.L1,Smpl.M1)!=0) { err=13; goto end;} // Заполнение матрицы 'A' для обратного преобразования Фурье..
if(FastFT(A,pR->Vect,M,-1)!=0) { err=14; goto end;} // Обратное преобразование Фурье для полей 'H'.
// Вычисление Ez или Hz, в зависимости от поляризации.
pR=ZDistr.Rows;
for(i=0;i<N;i++) { kx=k0X+(double)(i-Smpl.L1)*dK; 
if(Task.flPol==p_wave_Pol) ZFour.Vect[i]=-HFour[3].Vect[i]*kx*Cmplx_I/k0;
else ZFour.Vect[i]=-EFour[3].Vect[i]*kx*Cmplx_I/k0;
}
if(FillA(NF,A,ZFour,Smpl.L1,Smpl.M1)!=0) { err=15; goto end;} // Заполнение матрицы 'A' для обратного преобразования Фурье..
if(FastFT(A,pR->Vect,M,-1)!=0) { err=16; goto end;} // Обратное преобразование Фурье для полей 'H'.

// Сам цикл.
for(j=1;j<nDiv;j++) { 

pR=EDistr.Rows+j;
if(FillA(NF,A,EDistrFour.Rows[j-1],Smpl.L1,Smpl.M1)!=0) { err=17; goto end;} // Заполнение матрицы 'A' для обратного преобразования Фурье для 'EDistr'.
if(FastFT(A,pR->Vect,M,-1)!=0) { err=18; goto end;} // Обратное преобразование Фурье для полей 'Е'.

pR=HDistr.Rows+j;
if(FillA(NF,A,HDistrFour.Rows[j-1],Smpl.L1,Smpl.M1)!=0) { err=19; goto end;} // Заполнение матрицы 'A' для обратного преобразования Фурье для 'EDistr'.
if(FastFT(A,pR->Vect,M,-1)!=0) { err=20; goto end;} // Обратное преобразование Фурье для полей 'H'.

// Вычисление Ez или Hz, в зависимости от поляризации.
pR=ZDistr.Rows+j; 
for(i=0;i<N;i++) { kx=k0X+(double)(i-Smpl.L1)*dK; 
if(Task.flPol==p_wave_Pol) ZFour.Vect[i]=-HDistrFour.Rows[j-1].Vect[i]*kx*Cmplx_I/k0;
else ZFour.Vect[i]=-EDistrFour.Rows[j-1].Vect[i]*kx*Cmplx_I/k0;
}
if(FillA(NF,A,ZFour,Smpl.L1,Smpl.M1)!=0) { err=21; goto end;} // Заполнение матрицы 'A' для обратного преобразования Фурье..
if(FastFT(A,pR->Vect,M,-1)!=0) { err=22; goto end;} // Обратное преобразование Фурье для полей 'H'.

} // Конец цикла по 'j'.

// Отдельный просчёт на последней границе слоя.
pR=EDistr.Rows+nDiv;
if(FillA(NF,A,EFour[0],Smpl.L1,Smpl.M1)!=0) { err=23; goto end;} // Заполнение матрицы 'A' для обратного преобразования Фурье для 'EDistr'.
if(FastFT(A,pR->Vect,M,-1)!=0) { err=24; goto end;} // Обратное преобразование Фурье для полей 'Е'.
pR=HDistr.Rows+nDiv;
if(FillA(NF,A,HFour[0],Smpl.L1,Smpl.M1)!=0) { err=25; goto end;} // Заполнение матрицы 'A' для обратного преобразования Фурье для 'HDistr'.
if(FastFT(A,pR->Vect,M,-1)!=0) { err=26; goto end;} // Обратное преобразование Фурье для полей 'H'.
// Вычисление Ez или Hz, в зависимости от поляризации.
pR=ZDistr.Rows+nDiv;
for(i=0;i<N;i++) { kx=k0X+(double)(i-Smpl.L1)*dK; 
if(Task.flPol==p_wave_Pol) ZFour.Vect[i]=-HFour[0].Vect[i]*kx*Cmplx_I/k0;
else ZFour.Vect[i]=-EFour[0].Vect[i]*kx*Cmplx_I/k0;
}
if(FillA(NF,A,ZFour,Smpl.L1,Smpl.M1)!=0) { err=27; goto end;} // Заполнение матрицы 'A' для обратного преобразования Фурье..
if(FastFT(A,pR->Vect,M,-1)!=0) { err=28; goto end;} // Обратное преобразование Фурье для полей 'H'.

end: SAFE_DELETE_ARR(A); return err;
}


//-----------------------------------------------------------------------------------------------------------
// Заполнение массива 'A' для обратного преобразования Фурье.

BYTE FillArrInvFour(int NF,complex *A,const clRowMatr &Row,int L1,int M1)
{
int i;

if(NF<=0) return 1; if(A==NULL) return 2; if(Row.IsOK()!=0) return 3; if(L1<=0||M1<=0) return 4;
if(NF<L1+M1+1) return 5;

for(i=0;i<NF;i++) A[i]=Cmplx_0;
for(i=0;i<=M1;i++) A[i]=Row.Vect[i]; // Заполнение положительных гармоник.
for(i=0;i<L1;i++) A[NF-1-i]=Row.Vect[i+M1+1]; // Заполнение отрицательных гармоник.
return 0;
}
