// Определение типов.
typedef unsigned char BYTE;
typedef unsigned int UINT;
typedef signed char SCHAR;
typedef short int SHRT;
typedef unsigned short int USHRT;
typedef unsigned long U_LONG;

// Макросы.
#define SAFE_DELETE_ARR(Arr) { if(Arr!=NULL) { delete[] Arr; Arr=NULL;}} // Безопасное удаление массива.
#define SAFE_CLOSE(file) { if(file!=NULL) { fclose(file); file=NULL;}} // Безопасное закрытие файла.
#define MAX(x,y) (x>y?x:y) // Нахождение максимума.
#define MIN(x,y) (x<y?x:y) // Нахождение минимума.

// Константы.
#define szNameSubst 64 // Размер строки для названия вещества.
#define szFileName  32 // Размер строки для имени файла.
#define laySize_Max 1000 // Максимально допустимый размер слоя по оси x.
#define MinNumMultParal 8000 // Минимальное число операций, по времени эквивалентных умножению, для распараллеливания

// Сообщения.
#define messSingMatr 255 // Сингулярная матрица.
#define messSingLam  254 // Собственное значение близко к нулю.
#define messSmallEps 253 // Малое значение 'eps'.

// Названия файлов.
#define File_Name_Task "Task.txt" // Название файла параметров задачи.

// Типы поляризаций падающей волны.
#define nPol 2 // Число поляризаций.
#define p_wave_Pol 0 // Поле E лежит в плоскости падения.
#define s_wave_Pol 1 // Поле E перпендикулярно плоскости падения.

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Фундаментальные константы.

#define M_PI 3.1415926535897932384626433832795 // Число "Пи".
#define cVel 2.99792458e+8 // Скорость света в вакууме (м/с).
#define hPlank 1.054571726e-34 // Постоянная Планка (с планкой) (Дж.с). 
#define eChrg 1.60217657e-19 // Заряд электрона (Кл).

// Константы для формул Друде.
// Золото.
#define Eps0_Gold 7.9
#define gamma_Gold 113
#define wp_Gold 13400

// Малые константы.
#define SmCnst12_d 1.e-12
#define SmCnst32_d 1.e-32
#define LrgCnst64_d 1.e+64
