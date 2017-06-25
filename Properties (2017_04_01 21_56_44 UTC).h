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
// Определения типов.

typedef BYTE (*funWinSmooth)(float *Weight,int nPoints); // Функция окна сглаживания спектра.


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Структура диэлектрической проницаемости при данной длине волны.

struct strEps
{
double eps; double wlen; 
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Структура материала.

struct strMaterial
{
char *name; // Название материала.
BYTE flMeth; // Флаг способо задания eps. flMeth=0 - задание при постоянной длине волны. flMeth=1 - задание с помощью формул Друде. flMeth=2,3 - задание с помощью интерполяции. 2-n, 3-eps.
complex eps; // Диэлектрическая проницаемость данного материала.
double WaveLength; // Длина волны, для которой задана диэлектрическая проницаемость.
char chr; // Символ, используемый при кодировании данного материала.
struct strEps *nReal; // Показатель преломления материала. Действительная часть. (использование при flMeth=2;3)
struct strEps *nImag; // Показатель преломления материала. Мнимая часть. (использование при flMeth=2;3)
int Size_nReal,Size_nImag; // Размеры массивов.

strMaterial(void); // Конструктор.
~strMaterial(void); // Деструктор.
void ZeroDisp(void); // Обнуление массивов nImag,nReal.
void FreeDisp(void); // Освобождение памяти.
BYTE Read(FILE *fp); // Чтение данных о материале из файла.
BYTE WriteInfo(FILE *fp); // Запись в файл информации о материале.
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Класс массива всех материалов, используемых в кристалле (база данных материалов).

class clMaterialAll
{
void Zero(void); // Обнуление переменных.
void Free(void); // Освобождение памяти.

public:

int N; // Размер массива.
struct strMaterial *Mat; // Массив всех материалов.

clMaterialAll(void); // Конструктор.
~clMaterialAll(void); // Деструктор.
BYTE Alloc(int N_); // Выделение памяти для массива.
BYTE IsOK(void) const; // Проверки.
BYTE Read(char *File_Name); // Чтение данных о материалах из файла.
int GetNum(char chr,double WaveLength) const; // Поиск символа материала и возвращение его номера в списке сред.
BYTE SetEps(double wLength); // Установка eps материалов на данной длине волны формулами Друде.
BYTE WriteInfo(FILE *fp); // Запись в файл информации о материалах.
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Класс параметров задачи.

class clTask
{
void Zero(void); // Обнуление указателей и инициализация переменных.
void Free(void); // Освобождение памяти.
void ZeroFNames(void); // Обнуление указателей на строки - названия файлов.
void FreeFNames(void); // Удаление строк - названий файлов из памяти.
BYTE ReadFNames(FILE *file); // Чтение имён файлов.

public:

int L1,M1; // Размеры для преобразования Фурье.
BYTE flTypComp; // Флаг типа расчета. 0 - расчет распределения поля при заданной длине волны. 1 - расчет спектров. 2 - расчет угловой дисперсии.
double wLength,wLength_St,wLength_Fin,dwLength; // Длина волны падающего света.
double Theta,Theta_St,Theta_Fin,dTheta; // Угол падения волны (к нормали).
BYTE flPol; // Флаг типа поляризации падающей волны. 0 - p, 1 - s.
complex EpsInc,EpsOut; // Диэлектрическая проницаемость сред, откуда падает, и куда уходит свет.
char chEpsInc,chEpsOut; // Диэлектрическая проницаемость сред, откуда падает, и куда уходит свет.
double dz; // Дискрет для нахождения распределения полей по оси z (нанометры).
char *FName_Cryst; // Общее имя для всех файлов конфигурации.
char *FName_Mat; // Имя файла материалов.

clTask(void); // Конструктор.
~clTask(void); // Деструктор.
BYTE Read(char *File_Name); // Чтение данных из файла.
BYTE SetEpsIncOut(clMaterialAll &Matter); // Установка значений диэлектрической проницаемости сред, откуда падает, и куда уходит свет. 
BYTE WriteInfo(FILE *fp); // Запись параметров задачи в файл.
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Класс геометрии образца слоя фотонного кристалла.

class clLaySmplGeom
{
void Zero(void); // Обнуление переменных.
void Free(void); // Освобождение памяти.

public:

int laySize; // Размер периода образца по оси x.
double Period; // Период кристалла по оси x в нанометрах.
char *GeomSymbs; // Массив геометрии кристалла в символах материалов.
int *Geom; // Массив геометрии образца.

clLaySmplGeom(void); // Конструктор.
~clLaySmplGeom(void); // Деструктор.
BYTE Alloc(int laySize_); // Выделение памяти для массива.
BYTE IsOK(void) const; // Проверки.
BYTE Read(FILE *fp,const clMaterialAll &Matter,double WaveLength,int numLay); // Чтение геометрии образца из файла.
BYTE IsHomogen(void); // Проверка однородности образца.
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Декларации.

class clCrystalSmpl; // Класс образцов слоев фотонного кристалла.

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Структура одного образца слоя по оси z фотонного кристалла.

struct strLaySmpl
{
int numLay; // Номер образца.
clLaySmplGeom LayGeom; // Геометрия образца.
int L1,M1; // L1 - число отрицательных гармоник, M1 - положительных.
int M; // Степень двойки. Для Фурье преобразования.
int N; // Размер массивов компонент полей (L1+M1+1).
int N1; // Размер массивов 'Eps_Four', 'EpsInv_Four' 2*(L1+M1)+1.
int NF; // Размер массивов для быстрого преобразования Фурье.
complex *Eps_Four,*EpsInv_Four; // Результат Фурье пробразования для 'eps' и '1/eps'.
BYTE flHomogen; // Флаг однородности слоя.
complex EpsHomogen; // Диэлектрическая проницаемость для однородного слоя.
clMatrix Ek,Ak,EkInv,AkInv;
clMatrix Mk; // Матрица операторов уравнения Максвелла после преобразования Фурье.
clRowMatr EigLamMk; // Массив для собственных значений матрицы 'Mk' в виде столбца.
clMatrix EigVectMk,EigVectMkInv; // Матрицы для хранения собственных векторов матрицы 'Mk'.
clMatrix GamP[2][2],GamPInv[2][2]; // Блоки матриц 'EigVectMk','EigVectMkInv'.

strLaySmpl(void); // Конструктор.
~strLaySmpl(void); // Деструктор.
void Zero(void); // Обнуление переменных.
void Free(void); // Освобождение памяти.
BYTE IsOK(void) const; // Проверки.
BYTE IsOKFourier(void); // Проверки массива Фурье преобразования.
BYTE IsOKGamma(void) const; // Проверка блоков матриц 'Gamma'.
BYTE Read(FILE *fp,const clCrystalSmpl &Samples); // Чтение данных о слое из файла.
BYTE IsHomogen(void); // Установка флага однородности слоя.
BYTE FastFTLay(int L1_,int M1_,const clMaterialAll &Matter); // Быстрое преобразование Фурье для слоя.
BYTE SetParHomogen(int L1_,int M1_,const clMaterialAll &Matter); // Установка параметров для однородного слоя.
BYTE CompMk(const clTask &Task); // Вычисление матрицы Mk.
BYTE FindEigVectMk(void); // Нахождение собственных значений матрицы Mk.
BYTE FindEigVectMkHomogen(const clTask &Task); // Нахождение собственных значений и векторов матрицы Mk (случай однородного слоя).
BYTE FixFormEig(void); // Приведение столбца собственных значений к нужному в данной задаче виде (по возрастанию действительной части).
BYTE FindGamma(void); // Нахождение блоков матриц 'Gamma'.
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Класс образцов слоев фотонного кристалла.

class clCrystalSmpl
{
void Zero(void); // Обнуление переменных.
void Free(void); // Освобождение памяти.

public:

int laysNum; // Число образцов.
double Period; // Период кристалла по оси x в нанометрах.
struct strLaySmpl *Lays; // Массив всех образцов.
clMaterialAll Matter; // Материалы, из которых изготовлены образцы.
clTask Task; // Параметры задачи.

clCrystalSmpl(void); // Конструктор.
~clCrystalSmpl(void); // Деструктор.
BYTE Alloc(int laysNum_); // Выделение памяти для массива.
BYTE IsOK(void) const; // Проверки.
BYTE IsOKGamma(void); // Проверка блоков матриц 'Gamma' во всех слоях.
BYTE ReadInpData(void); // Чтение данных об образцах из файла.
BYTE FastFTAll(void); // Быстрое преобразование Фурье для всех образцов.
BYTE CompMk(void); // Вычисление матриц 'Mk' для всех образцов.
BYTE FindEigVectMk(void); // Нахождение собственных значений матриц 'Mk' для всех образцов.
strLaySmpl *Get(int num); // Получение указателя на образец слоя с номером 'num'.
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Структура одного слоя по оси z фотонного кристалла.

struct strLay
{
int numLay; // Номер слоя в кристалле.
int numSmpl; // Номер образца, которому соответствует данный слой.
double depth; // Толщина слоя в нанометрах.
int N; // Размер массивов компонент полей (L1+M1+1). Размер массивов DExp00Inv,DExp11.
clRowMatr DExp00Inv,DExp11; // Диагональные матрицы с экспонентами: обратная к 'Gamma00' и 'Gamma11'.
clMatrix GK[4]; // Матрицы 'Gk' для расчёта полей внутри кристалла (соответственно перед итерациями с 'P^-1','DExp','P' и после итерации с 'P').
clRowMatr EFour[4],HFour[4]; // Столбцы коэффициентов Фурье разложения полей Е и H соответственно на передней границе, перед умножением на 'P', после умножения на 'P^-1' и на задней границе слоя.
double dz; // Дискрет деления слоя для нахождения поля внутри слоя.
int nDiv; // Число разбиений слоя для нахождения распределения полей.
clRows EDistrFour,HDistrFour; // Фурье коэффициенты полей E и H. Для 'p' волны Ex,Hy, для 's' волны Ey,Hx.
clRows EDistr,HDistr,ZDistr; // Распределение компонент полей E,H и Ez или Hz в зависимости от поляризации. Для 'p' волны Ex,Hy,Ez, для 's' волны Ey,Hx,Hz.

strLay(void); // Конструктор.
void Zero(void); // Обнуление переменных.
BYTE IsOK(void); // Проверки.
BYTE FindGammaExp(const strLaySmpl &Smpl); // Вычисление матриц DExp00Inv, DExp11.
BYTE IsOKGammaExp(void); // Проверка блоков матриц 'Gamma'.
BYTE FindFullFieldDistr(clCrystalSmpl &Samples); // Нахождение распределения поля внутри слоя кристалла.
BYTE CompDz(double dz_); // Нахождение 'dz' из начального, приближённого 'dz' и выделение памяти для массивов столбцов коэффициентов Фурье полей E и H.
BYTE FindDistrFourEH(const strLaySmpl &Smpl); // // Нахождение коэффициентов Фурье разложения полей Е и H в слое.
BYTE FindDistrEH(const strLaySmpl &Smpl,clTask &Task); // Нахождение распределения полей Е и H в слое (выполнение обратного преобразования Фурье).
BYTE WriteInfo(FILE *fp,const strLaySmpl &Smpl); // Запись в файл информации о слое.
BYTE WriteFieldDistr(FILE *fp); // Запись в файл распределения полей.
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Класс фотонного кристалла.

class clCrystal
{
void Zero(void); // Обнуление переменных.
void Free(void); // Освобождение памяти.

public:

clCrystalSmpl Samples; // Информация обо всех образцах слоев кристалла.
int laysNum; // Число слоев.
struct strLay *Lays; // Массив всех слоев.
clRowMatr CInc,CIncInv,COut; // Блоки матриц 'С11','С11^-1' для среды, откуда падает свет, и 'С11' - для среды, куда проходит свет.
clRowMatr J,T,R; // Столбцы, соответствующие падающей волне и итоговым прошедшим и отражённым волнам.

clCrystal(void); // Конструктор.
~clCrystal(void); // Деструктор.
BYTE Alloc(int laysNum_); // Выделение памяти для массива.
BYTE IsOK(void) const; // Проверки.
BYTE ReadInpData(char *FileName_Task); // Чтение исходных данных.
BYTE ReadGeom(char *FileName_Cryst); // Чтение информации о кристалле из файла.
BYTE FindGammaExp(void); // Вычисление матриц 'DExp00Inv', 'DExp11' для всех слоёв.
BYTE IsOKGamma(void); // Проверка блоков матриц 'Gamma' во всех слоях.
BYTE CompC(void); // Вычисление матриц 'C'.
BYTE IsOK_C(void); // Проверка блоков матриц 'C'.
BYTE FindRTW(void); // Нахождение амплитуд прошедших и отражённых волн.
BYTE FindBoundsFieldDistr(void); // Нахождение столбцов амплитуд разложения Фурье распределения поля на границах слоёв.
BYTE FindFullFieldDistr(void); // Нахождение распределения поля внутри кристалла.
BYTE WriteInfo(FILE *fp); // Запись в файл информации о кристалле и задаче.
BYTE AngDistrOutput(void); // Вывод углового распределения по порядкам дифракции.
BYTE CompPartRefTr(double *pRef,double *pTr); // Вывод углового распределения по порядкам дифракции.
BYTE SpectrOutput(int n,double *T_Sp,double *R_Sp,double *A_Sp); // Вывод в файл всех спектров.
BYTE AnglDispOutput(int nRows,int n,clRows &T_Sp,clRows &R_Sp,clRows &A_Sp); // Вывод в файл всех спектров.
BYTE FieldOutput(void); // Вывод в файл распределения полей в кристалле.
BYTE Compute(double *Time); // Расчёт для фотонного кристалла.
BYTE ComputeSpectr(double *Time); // Расчёт спектров прохождения, отражения и поглощения для фотонного кристалла.
BYTE ComputeAnglDisp(double *Time); // Расчёт угловой диспперсии для фотонного кристалла.
};
