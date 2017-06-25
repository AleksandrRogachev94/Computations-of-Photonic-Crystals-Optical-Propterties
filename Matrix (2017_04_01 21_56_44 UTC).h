/*...........................................................................................................

clRowMatr,clMatrix

...........................................................................................................*/

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Класс столбца матрицы.

class clRowMatr {

public:
int N; complex *Vect;

clRowMatr(void); // Конструктор.
clRowMatr(const clRowMatr &); // Копирующий конструктор.
~clRowMatr(void); // Деструктор.
BYTE operator =(const clRowMatr &); // Оператор копирования.
void Zero(void); // Обнуление переменных.
void Free(void); // Освобождение памяти.
BYTE Alloc(int N_); // Выделение памяти для массива.
BYTE IsOK(void) const; // Проверки.
BYTE operator =(complex); // Оператор заполнения одинаковым числом.
clRowMatr operator -(); // Знак минус.
clRowMatr &operator +=(clRowMatr const &); // Прибавление столбца.
clRowMatr &operator -=(clRowMatr const &); // Вычитание столбца.
clRowMatr &operator +=(complex); // Прибавление единичной матрицы с множителем 'c'.
clRowMatr &operator -=(complex); // Вычитание единичной матрицы с множителем 'c'.
clRowMatr &operator *=(complex); // Умножение на комплексное число.
clRowMatr &operator /=(complex); // Деление на комплексное число.
clRowMatr &operator *=(double); // Умножение на действительное число.
clRowMatr &operator /=(double); // Деление на действительное число.
BYTE SetZero(void); // Обнуление вектора.
BYTE SetRand(double Lim); // Заполнение случайными числами в пределах Re,Im от -Lim до +Lim.
double GetAbsMax(void); // Получение максимального значения.
double GetNorm(void); // Получение нормы вектора.
BYTE Norm(void); // Нормировка вектора столбца.
BYTE Write(FILE *fp); // Запись столбца в файл.
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Операторы и функции для класса столбца.

clRowMatr operator +(const clRowMatr &,const clRowMatr &); // Сложение.
clRowMatr operator -(const clRowMatr &,const clRowMatr &); // Вычитание.
clRowMatr operator *(const clRowMatr &,complex); // Умножение на число.
clRowMatr operator /(const clRowMatr &,complex); // Деление на число.
clRowMatr operator *(complex,const clRowMatr &); // Умножение на число.
clRowMatr operator +(const clRowMatr &RM_,complex c); // Прибавление единичной матрицы с множителем 'c'.
clRowMatr operator +(complex c,const clRowMatr &RM_); // Прибавление единичной матрицы с множителем 'c'.
clRowMatr operator -(const clRowMatr &RM_,complex c); // Вычитание единичной матрицы с множителем 'c'.
clRowMatr operator *(const clRowMatr &,const clRowMatr &); // Умножение диагональных матриц.
clRowMatr operator !(const clRowMatr &); // Обратная диагональная матрица.
clRowMatr conj(clRowMatr const &); // Комплексное сопряжение.

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Класс массива столбцов 'clRowMatr'.

class clRows
{
public:

int nRows,N; clRowMatr *Rows;

clRows(void); // Конструктор.
clRows(const clRows &); // Копирующий конструктор.
~clRows(void); // Деструктор.
BYTE operator =(const clRows &); // Оператор копирования. 
void Zero(void); // Обнуление переменных.
void Free(void); // Освобождение памяти.
BYTE Alloc_(int nRows_); // Выделение памяти для объектов 'clRowMatr'.
BYTE Alloc(int nRows_,int N_); // Выделение памяти для всех массивов.
BYTE IsOK(void) const; // Проверки.
BYTE IsOK_All(void) const; // Проверки, включая столбцы.
clRowMatr *Get(int num); // Получение указателя на столбец с номером 'num'.
BYTE Write(FILE *fp); // Запись матрицы в файл.
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Класс двумерной матрицы.

class clMatrix {

public:

int Nx,Ny; complex *Matr;

clMatrix(void); // Конструктор.
clMatrix(const clMatrix &); // Копирующий конструктор.
~clMatrix(void); // Деструктор.
BYTE operator =(const clMatrix &); // Оператор копирования. 
void Zero(void); // Обнуление переменных.
void Free(void); // Освобождение памяти.
BYTE Alloc(int Nx_,int Ny_); // Выделение памяти для массива.
BYTE IsOK(void) const; // Проверки.
BYTE operator =(complex); // Заполнение матрицы с одинаковыми элементами на диагонали.
BYTE operator <<(complex); // Установка матрицы со всеми одинаковыми элементами.
BYTE operator =(clRowMatr const &); // Заполнение матрицы с диагональными элементами.
clMatrix operator -(); // Знак минус.
clMatrix &operator +=(complex); // Прибавление единичной матрицы с множителем 'c'.
clMatrix &operator -=(complex); // Вычитание единичной матрицы с множителем 'c'.
clMatrix &operator *=(complex); // Умножение на комплексное число.
clMatrix &operator /=(complex); // Деление на комплексное число.
clMatrix &operator *=(double); // Умножение на действительное число.
clMatrix &operator /=(double); // Деление на действительное число.
clMatrix &operator +=(const clRowMatr &); // Прибавление диагональной матрицы.
clMatrix &operator -=(const clRowMatr &); // Вычитание диагональной матрицы.
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Операторы для класса двумерной матрицы.

clMatrix operator +(const clMatrix &,const clMatrix &); // Сложение.
clMatrix operator -(const clMatrix &,const clMatrix &); // Вычитание.
clMatrix operator *(const clMatrix &,const clMatrix &); // Умножение.
clMatrix operator ~(const clMatrix &); // Транспонированная матрица.
clMatrix operator +(const clMatrix &,complex); // Прибавление единичной матрицы с множителем 'c'.
clMatrix operator +(complex,const clMatrix &); // Прибавление единичной матрицы с множителем 'c'.
clMatrix operator -(const clMatrix &,complex); // Вычитание единичной матрицы с множителем 'c'.
clMatrix operator -(complex,const clMatrix &); // Вычитание из единичной матрицы с множителем 'c'.
clMatrix operator *(const clMatrix &,complex); // Умножение матрицы на число.
clMatrix operator *(complex,const clMatrix &); // Умножение матрицы на число.
BYTE Transp(const clMatrix &M,clMatrix &MT); // Транспонирование матрицы.
BYTE Hermit(const clMatrix &M,clMatrix &MH); // Эрмитовое сопряжение матрицы.

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Операторы для классов двумерной матрицы и столбца.

clRowMatr operator *(const clMatrix &,const clRowMatr &); // Умножение матрицы на столбец.
clMatrix operator +(const clMatrix &,const clRowMatr &); // Сложение матрицы с диагональной матрицей.
clMatrix operator +(const clRowMatr &,const clMatrix &); // Сложение диагональной матрицы с матрицей.
clMatrix operator -(const clMatrix &,const clRowMatr &); // Вычитание диагональной матрицы из матрицы.
clMatrix operator -(const clRowMatr &,const clMatrix &); // Вычитание матрицы из диагональной матрицы.
clMatrix operator %(const clMatrix &,const clRowMatr &); // Умножение квадратной матрицы на диагональную матрицу.
clMatrix operator %(const clRowMatr &,const clMatrix &); // Умножение диагональной матрицы на квадратную матрицу.
