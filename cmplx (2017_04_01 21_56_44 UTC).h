//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Класс комплексных чисел.

class complex {

public:

double re,im;

complex(); // Конструктор.
complex(double re_,double im_); // Конструктор.
complex operator +(); // Знак плюс.
complex operator -(); // Знак минус.
complex operator +=(const complex &); // Прибавление комплексного числа.
complex operator -=(const complex &); // Вычитание комплексного числа.
complex operator *=(const complex &); // Умножение на комплексное число.
complex operator /=(const complex &); // Деление на комплексное число.
complex operator =(double); // Присваивание действительного числа.
complex operator *=(double); // Умножение на число.
complex operator /=(double); // Деление на число.
};

complex operator +(const complex &,const complex &); // Сложение.
complex operator -(const complex &,const complex &); // Вычитание.
complex operator *(const complex &,const complex &); // Умножение.
complex operator /(const complex &,const complex &); // Деление.
complex operator *(const complex &,double); // Умножение на число.
complex operator *(double,const complex &); // Умножение на число.
complex operator /(const complex &,double); // Деление на число.
complex operator /(double,const complex &); // Деление на число.

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Функции для комплексных чисел.

complex polar(double Abs,double Phase); // Создание комплексного числа по амплитуде и фазе.
complex conj(const complex &c); // Комплексное сопряжение.
double real(const complex &c); // Действительная часть.
double imag(const complex &c); // Мнимая часть.
double abs(const complex &c); // Модуль комплексного числа.
double arg(const complex &c); // Аргумент комплексного числа.
complex sqrt(const complex &c); // Квадратный корень.
complex pow(const complex &c,double p); // Возведение в степень.
complex exp(const complex &c); // Экспонента.
complex log(const complex &c); // Логарифм.
complex Inv(const complex &c); // Обратное число.

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Константы.

#define Cmplx_0 complex(0.,0.)
#define Cmplx_1 complex(1.,0.)
#define Cmplx_I complex(0.,1.)
