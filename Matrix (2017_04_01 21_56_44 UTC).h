/*...........................................................................................................

clRowMatr,clMatrix

...........................................................................................................*/

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ����� ������� �������.

class clRowMatr {

public:
int N; complex *Vect;

clRowMatr(void); // �����������.
clRowMatr(const clRowMatr &); // ���������� �����������.
~clRowMatr(void); // ����������.
BYTE operator =(const clRowMatr &); // �������� �����������.
void Zero(void); // ��������� ����������.
void Free(void); // ������������ ������.
BYTE Alloc(int N_); // ��������� ������ ��� �������.
BYTE IsOK(void) const; // ��������.
BYTE operator =(complex); // �������� ���������� ���������� ������.
clRowMatr operator -(); // ���� �����.
clRowMatr &operator +=(clRowMatr const &); // ����������� �������.
clRowMatr &operator -=(clRowMatr const &); // ��������� �������.
clRowMatr &operator +=(complex); // ����������� ��������� ������� � ���������� 'c'.
clRowMatr &operator -=(complex); // ��������� ��������� ������� � ���������� 'c'.
clRowMatr &operator *=(complex); // ��������� �� ����������� �����.
clRowMatr &operator /=(complex); // ������� �� ����������� �����.
clRowMatr &operator *=(double); // ��������� �� �������������� �����.
clRowMatr &operator /=(double); // ������� �� �������������� �����.
BYTE SetZero(void); // ��������� �������.
BYTE SetRand(double Lim); // ���������� ���������� ������� � �������� Re,Im �� -Lim �� +Lim.
double GetAbsMax(void); // ��������� ������������� ��������.
double GetNorm(void); // ��������� ����� �������.
BYTE Norm(void); // ���������� ������� �������.
BYTE Write(FILE *fp); // ������ ������� � ����.
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ��������� � ������� ��� ������ �������.

clRowMatr operator +(const clRowMatr &,const clRowMatr &); // ��������.
clRowMatr operator -(const clRowMatr &,const clRowMatr &); // ���������.
clRowMatr operator *(const clRowMatr &,complex); // ��������� �� �����.
clRowMatr operator /(const clRowMatr &,complex); // ������� �� �����.
clRowMatr operator *(complex,const clRowMatr &); // ��������� �� �����.
clRowMatr operator +(const clRowMatr &RM_,complex c); // ����������� ��������� ������� � ���������� 'c'.
clRowMatr operator +(complex c,const clRowMatr &RM_); // ����������� ��������� ������� � ���������� 'c'.
clRowMatr operator -(const clRowMatr &RM_,complex c); // ��������� ��������� ������� � ���������� 'c'.
clRowMatr operator *(const clRowMatr &,const clRowMatr &); // ��������� ������������ ������.
clRowMatr operator !(const clRowMatr &); // �������� ������������ �������.
clRowMatr conj(clRowMatr const &); // ����������� ����������.

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ����� ������� �������� 'clRowMatr'.

class clRows
{
public:

int nRows,N; clRowMatr *Rows;

clRows(void); // �����������.
clRows(const clRows &); // ���������� �����������.
~clRows(void); // ����������.
BYTE operator =(const clRows &); // �������� �����������. 
void Zero(void); // ��������� ����������.
void Free(void); // ������������ ������.
BYTE Alloc_(int nRows_); // ��������� ������ ��� �������� 'clRowMatr'.
BYTE Alloc(int nRows_,int N_); // ��������� ������ ��� ���� ��������.
BYTE IsOK(void) const; // ��������.
BYTE IsOK_All(void) const; // ��������, ������� �������.
clRowMatr *Get(int num); // ��������� ��������� �� ������� � ������� 'num'.
BYTE Write(FILE *fp); // ������ ������� � ����.
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ����� ��������� �������.

class clMatrix {

public:

int Nx,Ny; complex *Matr;

clMatrix(void); // �����������.
clMatrix(const clMatrix &); // ���������� �����������.
~clMatrix(void); // ����������.
BYTE operator =(const clMatrix &); // �������� �����������. 
void Zero(void); // ��������� ����������.
void Free(void); // ������������ ������.
BYTE Alloc(int Nx_,int Ny_); // ��������� ������ ��� �������.
BYTE IsOK(void) const; // ��������.
BYTE operator =(complex); // ���������� ������� � ����������� ���������� �� ���������.
BYTE operator <<(complex); // ��������� ������� �� ����� ����������� ����������.
BYTE operator =(clRowMatr const &); // ���������� ������� � ������������� ����������.
clMatrix operator -(); // ���� �����.
clMatrix &operator +=(complex); // ����������� ��������� ������� � ���������� 'c'.
clMatrix &operator -=(complex); // ��������� ��������� ������� � ���������� 'c'.
clMatrix &operator *=(complex); // ��������� �� ����������� �����.
clMatrix &operator /=(complex); // ������� �� ����������� �����.
clMatrix &operator *=(double); // ��������� �� �������������� �����.
clMatrix &operator /=(double); // ������� �� �������������� �����.
clMatrix &operator +=(const clRowMatr &); // ����������� ������������ �������.
clMatrix &operator -=(const clRowMatr &); // ��������� ������������ �������.
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ��������� ��� ������ ��������� �������.

clMatrix operator +(const clMatrix &,const clMatrix &); // ��������.
clMatrix operator -(const clMatrix &,const clMatrix &); // ���������.
clMatrix operator *(const clMatrix &,const clMatrix &); // ���������.
clMatrix operator ~(const clMatrix &); // ����������������� �������.
clMatrix operator +(const clMatrix &,complex); // ����������� ��������� ������� � ���������� 'c'.
clMatrix operator +(complex,const clMatrix &); // ����������� ��������� ������� � ���������� 'c'.
clMatrix operator -(const clMatrix &,complex); // ��������� ��������� ������� � ���������� 'c'.
clMatrix operator -(complex,const clMatrix &); // ��������� �� ��������� ������� � ���������� 'c'.
clMatrix operator *(const clMatrix &,complex); // ��������� ������� �� �����.
clMatrix operator *(complex,const clMatrix &); // ��������� ������� �� �����.
BYTE Transp(const clMatrix &M,clMatrix &MT); // ���������������� �������.
BYTE Hermit(const clMatrix &M,clMatrix &MH); // ��������� ���������� �������.

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ��������� ��� ������� ��������� ������� � �������.

clRowMatr operator *(const clMatrix &,const clRowMatr &); // ��������� ������� �� �������.
clMatrix operator +(const clMatrix &,const clRowMatr &); // �������� ������� � ������������ ��������.
clMatrix operator +(const clRowMatr &,const clMatrix &); // �������� ������������ ������� � ��������.
clMatrix operator -(const clMatrix &,const clRowMatr &); // ��������� ������������ ������� �� �������.
clMatrix operator -(const clRowMatr &,const clMatrix &); // ��������� ������� �� ������������ �������.
clMatrix operator %(const clMatrix &,const clRowMatr &); // ��������� ���������� ������� �� ������������ �������.
clMatrix operator %(const clRowMatr &,const clMatrix &); // ��������� ������������ ������� �� ���������� �������.
