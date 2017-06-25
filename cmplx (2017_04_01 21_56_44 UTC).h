//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ����� ����������� �����.

class complex {

public:

double re,im;

complex(); // �����������.
complex(double re_,double im_); // �����������.
complex operator +(); // ���� ����.
complex operator -(); // ���� �����.
complex operator +=(const complex &); // ����������� ������������ �����.
complex operator -=(const complex &); // ��������� ������������ �����.
complex operator *=(const complex &); // ��������� �� ����������� �����.
complex operator /=(const complex &); // ������� �� ����������� �����.
complex operator =(double); // ������������ ��������������� �����.
complex operator *=(double); // ��������� �� �����.
complex operator /=(double); // ������� �� �����.
};

complex operator +(const complex &,const complex &); // ��������.
complex operator -(const complex &,const complex &); // ���������.
complex operator *(const complex &,const complex &); // ���������.
complex operator /(const complex &,const complex &); // �������.
complex operator *(const complex &,double); // ��������� �� �����.
complex operator *(double,const complex &); // ��������� �� �����.
complex operator /(const complex &,double); // ������� �� �����.
complex operator /(double,const complex &); // ������� �� �����.

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ������� ��� ����������� �����.

complex polar(double Abs,double Phase); // �������� ������������ ����� �� ��������� � ����.
complex conj(const complex &c); // ����������� ����������.
double real(const complex &c); // �������������� �����.
double imag(const complex &c); // ������ �����.
double abs(const complex &c); // ������ ������������ �����.
double arg(const complex &c); // �������� ������������ �����.
complex sqrt(const complex &c); // ���������� ������.
complex pow(const complex &c,double p); // ���������� � �������.
complex exp(const complex &c); // ����������.
complex log(const complex &c); // ��������.
complex Inv(const complex &c); // �������� �����.

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ���������.

#define Cmplx_0 complex(0.,0.)
#define Cmplx_1 complex(1.,0.)
#define Cmplx_I complex(0.,1.)
