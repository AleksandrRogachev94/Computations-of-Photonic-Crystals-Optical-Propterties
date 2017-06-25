// ����������� �����.
typedef unsigned char BYTE;
typedef unsigned int UINT;
typedef signed char SCHAR;
typedef short int SHRT;
typedef unsigned short int USHRT;
typedef unsigned long U_LONG;

// �������.
#define SAFE_DELETE_ARR(Arr) { if(Arr!=NULL) { delete[] Arr; Arr=NULL;}} // ���������� �������� �������.
#define SAFE_CLOSE(file) { if(file!=NULL) { fclose(file); file=NULL;}} // ���������� �������� �����.
#define MAX(x,y) (x>y?x:y) // ���������� ���������.
#define MIN(x,y) (x<y?x:y) // ���������� ��������.

// ���������.
#define szNameSubst 64 // ������ ������ ��� �������� ��������.
#define szFileName  32 // ������ ������ ��� ����� �����.
#define laySize_Max 1000 // ����������� ���������� ������ ���� �� ��� x.
#define MinNumMultParal 8000 // ����������� ����� ��������, �� ������� ������������� ���������, ��� �����������������

// ���������.
#define messSingMatr 255 // ����������� �������.
#define messSingLam  254 // ����������� �������� ������ � ����.
#define messSmallEps 253 // ����� �������� 'eps'.

// �������� ������.
#define File_Name_Task "Task.txt" // �������� ����� ���������� ������.

// ���� ����������� �������� �����.
#define nPol 2 // ����� �����������.
#define p_wave_Pol 0 // ���� E ����� � ��������� �������.
#define s_wave_Pol 1 // ���� E ��������������� ��������� �������.

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ��������������� ���������.

#define M_PI 3.1415926535897932384626433832795 // ����� "��".
#define cVel 2.99792458e+8 // �������� ����� � ������� (�/�).
#define hPlank 1.054571726e-34 // ���������� ������ (� �������) (��.�). 
#define eChrg 1.60217657e-19 // ����� ��������� (��).

// ��������� ��� ������ �����.
// ������.
#define Eps0_Gold 7.9
#define gamma_Gold 113
#define wp_Gold 13400

// ����� ���������.
#define SmCnst12_d 1.e-12
#define SmCnst32_d 1.e-32
#define LrgCnst64_d 1.e+64
