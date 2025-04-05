/*
   Boneh & Franklin's Identity Based Encryption
   
   Set-up phase

   After this program has run the file common.ibe contains

   <Size of prime modulus in bits>
   <Prime p>
   <Prime q (divides p+1) >
   <Point P - x coordinate>
   <Point P - y coordinate>
   <Point Ppub - x coordinate>
   <Point Ppub - y coordinate>
   <Cube root of unity in Fp2 - x component >
   <Cube root of unity in Fp2 - y component >

   The file master.ibe contains

   <The master secret s>

   Requires: zzn2.cpp big.cpp zzn.cpp ecn.cpp

 */
#include <time.h>
#include <iostream>
#include <fstream>
#include <cstring>
using namespace std;

#include "ecn.h"
#include "zzn.h"
#include "zzn2.h"

//C++���������߱������ⲿ�ִ�����C���Ա�д�ģ����������ʱӦ��ʹ��C���Ե�����Լ��
extern "C"
{
   #include"miracl.h"
   #include"mirdef.h"
}

//2015�����
//�ṩһ�������Բ�,Ϊ���ڲ�ͬ�汾��MSVC����ʱ��֮���ṩһ�ּ����Ի��ƣ�ʹ�þɴ����ܹ����°汾�Ŀ����������У�������Ҫ�Դ�����д����޸ġ�
FILE* __cdecl __iob_func(unsigned i) {
    return __acrt_iob_func(i);
}

#ifdef __cplusplus
extern "C" 
#endif
//extern "C" { FILE __iob_func[3] = { *stdin,*stdout,*stderr }; }
FILE _iob[3] = {__iob_func(0)[0], __iob_func(1)[1], __iob_func(2)[2]}; 

//һ���ض���Microsoft Visual C++��������ָ����ڿ�������������Ϊ�����������Ǹ��������������ӹ����в�ҪĬ�����ӵ�libc.lib����⡣
#pragma comment(linker, "/NODEFAULTLIB:libc.lib")


#define renum 2500

#define HASH_LEN 32
#define HASH_LEN1 20   //������H2����Ϊ��������q��160λ�Ķ�����������160/8=20
                                        //H2�в���sha256Ҫ��HASH_LEN1 ������32�ı�������ˣ��Լ���H2�ڲ��������Ϊsha-1


#define PBITS 512
#define QBITS 160

// Using SHA-256 as basic hash algorithm

//
// Define one or the other of these
//
// Which is faster depends on the I/M ratio - See imratio.c
// Roughly if I/M ratio > 16 use PROJECTIVE, otherwise use AFFINE
//

// #define AFFINE
#define PROJECTIVE

// Define this to use this idea ftp://ftp.computing.dcu.ie/pub/resources/crypto/short.pdf
// which enables denominator elimination
#define SCOTT

Miracl precision(16,0);  // increase if PBITS increases. (32,0) for 1024 bit p

/*----------------------------------------------------------------------------Tate Paring ��������Ҫ�ĺ���-----------------------------------------------------*/
void extract(ECn& A,ZZn& x,ZZn& y)  //��������
{ 
    x=(A.get_point())->X;
    y=(A.get_point())->Y;
}

void extract(ECn& A,ZZn& x,ZZn& y,ZZn& z)  //��Ӱ����
{ 
    big t;
    x=(A.get_point())->X;
    y=(A.get_point())->Y;
    t=(A.get_point())->Z;
    if (A.get_status()!=MR_EPOINT_GENERAL) 
        z=1;
    else                                   
        z=t;
}

//
// Line from A to destination C. Let A=(x,y)
// Line Y-slope.X-c=0, through A, so intercept c=y-slope.x
// Line Y-slope.X-y+slope.x = (Y-y)-slope.(X-x) = 0
// Now evaluate at Q -> return (Qy-y)-slope.(Qx-x)
//

ZZn2 line(ECn& A, ECn& C, ZZn& slope, ZZn2& Qx, ZZn2& Qy)  //������Բ�����ϵĵ�Q��A��C��ֱ��
{ 
    ZZn2 n=Qx,w=Qy;
    ZZn x,y,z,t;
#ifdef AFFINE
    extract(A,x,y);
    n-=x; n*=slope;            // 2 ZZn muls
    w-=y; n-=w;
#endif
#ifdef PROJECTIVE
    extract(A,x,y,z);
    x*=z; t=z; z*=z; z*=t;          
    n*=z; n-=x;                // 9 ZZn muls
    w*=z; w-=y; 
    extract(C,x,y,z);
    w*=z; n*=slope; n-=w;                     
#endif
    return n;
}

#ifndef SCOTT

//
// Vertical line through point A
//

ZZn2 vertical(ECn& A,ZZn2& Qx)
{
    ZZn2 n=Qx;
    ZZn x,y,z;
#ifdef AFFINE
    extract(A,x,y);
    n-=x;
#endif
#ifdef PROJECTIVE
    extract(A,x,y,z);
    z*=z;                    
    n*=z; n-=x;                // 3 ZZn muls
#endif
    return n;
}

#endif

//
// Add A=A+B  (or A=A+A) 
// Bump up num and denom
//
// AFFINE doubling     - 12 ZZn muls, plus 1 inversion
// AFFINE adding       - 11 ZZn muls, plus 1 inversion
//
// PROJECTIVE doubling - 26 ZZn muls
// PROJECTIVE adding   - 34 ZZn muls
//


void g(ECn& A,ECn& B,ZZn2& Qx,ZZn2& Qy,ZZn2& num) 
{
    ZZn  lam,mQy;
    ZZn2 d,u;
    big ptr;
    ECn P=A;

// Evaluate line from A
    ptr=A.add(B);

#ifndef SCOTT
    if (A.iszero())   { u=vertical(P,Qx); d=1; }
    else
    {
#endif
        if (ptr==NULL)
            u=1;
        else 
        {
            lam=ptr;
            u=line(P,A,lam,Qx,Qy);
        }
#ifndef SCOTT
        d=vertical(A,Qx);
    }

    num*=(u*conj(d));    // 6 ZZn muls  
#else
// denominator elimination!
    num*=u;
#endif
}

//
// Tate Pairing 
//

BOOL fast_tate_pairing(ECn& P, ZZn2& Qx, ZZn2& Qy, Big& q, ZZn2& res) //P:����Ԫ��Qx,Qy:��Բ�����ϵĵ�Q�����꣬q:�����ף�res:˫���ԶԽ��
{ 
    int i,nb;
    Big n,p;
    ECn A;


// q.P = 2^17*(2^142.P +P) + P

    res=1;
    A=P;    // reset A

#ifdef SCOTT
// we can avoid last iteration..
    n=q-1;
#else
    n=q;
#endif
    nb=bits(n);

    for (i=nb-2;i>=0;i--)
    {
        res*=res;         
        g(A,A,Qx,Qy,res); 
        if (bit(n,i))
            g(A,P,Qx,Qy,res);       
    }

#ifdef SCOTT
    if (A!=-P || res.iszero()) return FALSE;
#else
    if (!A.iszero()) return FALSE;
#endif

    p=get_modulus();         // get p
    res= pow(res,(p+1)/q);   // raise to power of (p^2-1)/q
    res=conj(res)/res;
    if (res.isunity()) return FALSE;
    return TRUE;   
}
BOOL ecap(ECn& P,ECn& Q,Big& order,ZZn2& cube,ZZn2& res)  //P:����Ԫ��Q:����㣬order:�����ף�cube:��������res:˫���ԶԽ��
{
     ZZn2 Qx,Qy;
     Big xx,yy;
#ifdef SCOTT
     ZZn a,b,x,y,ib,w,t1,y2,ib2;
#else
     ZZn2 lambda,ox;
#endif
     Q.get(xx,yy);
     Qx=(ZZn)xx*cube;
     Qy=(ZZn)yy;

#ifndef SCOTT
// point doubling
     lambda=(3*Qx*Qx)/(Qy+Qy);
     ox=Qx;
     Qx=lambda*lambda-(Qx+Qx);
     Qy=lambda*(ox-Qx)-Qy;
#else
 //explicit point subtraction
     Qx.get(a,b);
     y=yy;
     ib=(ZZn)1/b;

     t1=a*b*b;
     y2=y*y;
     ib2=ib*ib;
     w=y2+2*t1;
     x=-w*ib2;
     y=-y*(w+t1)*(ib2*ib);
     Qx.set(x); 
     Qy.set((ZZn)0,y);

#endif

     if (fast_tate_pairing(P,Qx,Qy,order,res)) return TRUE;
     return FALSE;
}


//
// ecap(.) function - apply distortion map
//
// Qx is in ZZn if SCOTT is defined. Qy is in ZZn if SCOTT is not defined. 
// This can be exploited for some further optimisations. 
/*----------------------------------------------------------------------------Tate Paring ��������Ҫ�ĺ���-----------------------------------------------------*/


/*----------------------------------------------------------------------------���Hash��������ĺ���-----------------------------------------------------*/
// ʵ����һ�����ַ�����ϣ��С��ģ�� p �Ĵ����� Big �ĺ��� H1
Big H1(char *string)
{ // Hash a zero-terminated string to a number < modulus
    Big h,p;
    char s[HASH_LEN];
    int i,j; 
    sha256 sh;

    shs256_init(&sh);

    for (i=0;;i++)
    {
        if (string[i]==0) 
            break;
        shs256_process(&sh,string[i]);
    }
    shs256_hash(&sh,s);
    p=get_modulus();
	//cout<<"modulus"<<p<<endl;//�Լ��ӵĲ鿴pֵ����䣬ͨ��pֵ��֪get_modulus()������get_mip()������
	//��get_mip()�õ����ǵ�ǰ��������Ⱥ�Ľ�ֵq.
    h=1; j=0; i=1;
    forever
    {
        h*=256; 
        if (j==HASH_LEN)  
        {h+=i++; j=0;}
        else        
            h+=s[j++];
        if (h>=p)
            break;
    }
    h%=p;
    return h;
}

//��δ���ʵ����һ���� Fp2 Ԫ�ع�ϣ��һ�� n �ֽ��ַ����ĺ��� H2�����������Ŀ���ǽ� Fp2 ���е�Ԫ��ӳ�䵽һ���̶����ȵ��ַ���
//�����������Ҫ�����ǽ� Fp2 ���е�Ԫ��ͨ�� SHA-1 ��ϣ�㷨ת��Ϊһ���̶����ȵ��ַ���
int H2(ZZn2 x,char *s)
{ // Hash an Fp2 to an n-byte string s[.]. Return n
    sha sh;
    Big a,b;
    int m;  

    shs_init(&sh);
    x.get(a,b);

    while (a>0)
    {
        m=a%160;
        shs_process(&sh,m);
        a/=160;
    }
    while (b>0)
    {
        m=b%160;
        shs_process(&sh,m);
        b/=160;
    }
    shs_hash(&sh,s);

  return HASH_LEN1;
	/*sha256 sh;
    Big a,b;
    int m;

    shs256_init(&sh);
    x.get(a,b);

    while (a>0)
    {
        m=a%256;
        shs256_process(&sh,m);
        a/=256;
    }
    while (b>0)
    {
        m=b%256;
        shs256_process(&sh,m);
        b/=256;
    }
    shs256_hash(&sh,s);

  return HASH_LEN1;
	// return 20;*/
}

// ��δ���ʵ����һ��������ֹ�ַ�����ϣ��С�ڸ���ģ�� qm �Ĵ����� Big �ĺ��� H3��
// �����������Ҫ�����ǽ�������ַ���ͨ�� SHA-1 ��ϣ�㷨ת��Ϊһ���̶����ȵĹ�ϣֵ��
// Ȼ�������ϣֵӳ�䵽һ������������ȷ�����������С�ڸ�����ģ�� qm
Big H3(char *string,Big qm)
{ // Hash a zero-terminated string to a number < modulus q
    Big h;
    char s[HASH_LEN1];
    int i,j; 
    sha sh;

    shs_init(&sh);

    for (i=0;;i++)
    {
        if (string[i]==0) break;
        shs_process(&sh,string[i]);
    }
    shs_hash(&sh,s);
    //q=get_modulus();
	//cout<<"modulus"<<p<<endl;//�Լ��ӵĲ鿴pֵ����䣬ͨ��pֵ��֪get_modulus()�õ�����Բ�������������������P
    h=1; j=0; i=1;
    forever
    {
        h*=160; 
        if (j==HASH_LEN1)  
        {h+=i++; j=0;}
        else        
            h+=s[j++];
        if (h>=qm) break;
    }
    h%=qm;
    return h;
}
   
//
// Given y, get x=(y^2-1)^(1/3) mod p (from curve equation)
// �ڸ��� (y) ֵ������£�
// �ҵ�������Բ���߷��� (y^2 = x^3 + 1) �� (x) ֵ
//

Big getx(Big y)
{
    Big p=get_modulus();
    Big t=modmult(y+1,y-1,p);   // avoids overflow
    return pow(t,(2*p-1)/3,p);
}
 
// MapToPoint
//��һ���ַ�����ʶ����ID��ӳ�䵽��Բ�����ϵ�һ����Ĺ���
ECn map_to_point(char *ID)
{
    ECn Q;
    Big x0,y0=H1(ID);
 
    x0=getx(y0);

    Q.set(x0,y0);

    return Q;
}

/*-------------------------------------------------------------------�ṹ������-------------------------------------------------------------*/
// ��������
typedef struct params{
    ECn P, P1, P2, P3, P4;
    ECn u[257];
    ECn v[257];
    ZZn2 e_g1_g2;
    ZZn2 e_g3_g4;
}params;

// PKG˽Կ
typedef struct sk_pkg {
    Big alpha;
    ECn sk_pkg;
}sk_pkg;

// TimeServer˽Կ
typedef struct sk_ts {
    Big beta;
    ECn sk_ts;
}sk_ts;

// ˽Կ�ṹ��
typedef struct PrivateKey {
    ECn d1;
    ECn d2;
}PrivateKey;

// ʱ�����Žṹ��
typedef struct TimetrapDoor {
    ECn st1;
    ECn st2;
}Timetrapdoor;

typedef struct Ciphtertext {
    ZZn2 c1;
    ECn c2, c3, c4, c5;
}Ciphtertext;



/*------------------------------------�Զ��庯��--------------------------------------------*/
void SkGen(params&params, sk_pkg& sk_pkg, Big q, ECn& Qid, PrivateKey& userkey) {
    Big r1 = rand(q);
    userkey.d1 = sk_pkg.sk_pkg + ( r1 * Qid );
    userkey.d2 = r1 * params.P;

	//cout << "489-userkey.d1:" << userkey.d1 << endl;
	//cout << "490-userkey.d2:" << userkey.d2 << endl;
    
}

void StGen(params& params, sk_ts& sk_ts, Big q, ECn& T, Timetrapdoor& st) {
    Big r4 = rand(q);
    st.st1 = sk_ts.sk_ts + ( r4 * T) ;
    st.st2 = r4 * params.P;
}

void Enc(ZZn2& g0, Big& m, params& params, ECn& QID, ECn& QT, Big q, Ciphtertext& CT) {
    Big r2;
    Big r3;
    r2 = rand(q);
    r3 = rand(q);

    ZZn2 Q;
    Q = pow(g0, m);

    //cout << "m:" << m << endl;
    //cout << "Q:" << Q << endl;
    //cout << "510:QID" << QID << endl;

    
    //cout << "511:m:" << m << endl;
    CT.c1 = Q *  pow(params.e_g1_g2, r2) * pow(params.e_g3_g4, r3);
	//CT.c1 *= g0;
	//cout << "515-CT.c1:" << CT.c1 << endl;
    //cout << "514: m:" << m << endl;
    //cout << "m after dec:" << CT.c1 / pow(params.e_g1_g2, r2) << endl;
    //cout << "518-r2" << r2 << endl;
    CT.c2 = r2 * params.P;
    //cout << "params.P:" << params.P << endl;
    //cout << "521:" << CT.c2 << endl;
    //cout << "520-r2:" << CT.c2 - (r2 - 1) * params.P << endl;

    //cout << "524:QID" << QID << endl;
    CT.c3 = r2 * QID;
	//cout << "526:" << CT.c3 << endl;
	//cout << "527:" << CT.c3 - (r2 - 1) * QID << endl;

    CT.c4 = r3 * params.P;

    CT.c5 = r3 * QT;
}

void Dec(Ciphtertext& CT, PrivateKey& privatekey, TimetrapDoor& st,Big& m1dec, ZZn2& cube, ZZn2& g0, Big q) {
	

    ZZn2 e1, e2, e3, e4, temp, Q;
    
    ecap(privatekey.d2, CT.c3, q, cube, e1);
    
	if (ecap(privatekey.d2, CT.c3, q, cube, e1)) {
		cout << "e1 success" << endl;
	}
	else {
		cout << "e1 error" << endl;
	}

    ecap(st.st2, CT.c5, q, cube ,e2);
    if (ecap(st.st2, CT.c5, q, cube, e2)) {
        cout << "e2 success" << endl;
    }
    else {
        cout << "e2 error" << endl;
    }


    ecap(CT.c2, privatekey.d1, q, cube, e3);
    if (ecap(CT.c2, privatekey.d1, q, cube, e3)) {
        cout << "e3 success" << endl;
    }
    else {
        cout << "e3 error" << endl;
    }

    ecap(CT.c4, st.st1, q, cube , e4);
    if (ecap(CT.c4, st.st1, q, cube, e4)) {
        cout << "e4 success" << endl;
    }
    else {
        cout << "e4 error" << endl;
    }



	//mdec = CT.c1 * e1 / e3;
	//cout << "550-mdec" << mdec << endl;
    temp = CT.c1;
    //cout << "temp0:" << temp << endl;
	temp *= e1;
    //cout << "temp1:" << temp << endl;
	temp *= e2;
    //cout << "temp2:" << temp << endl;
	temp /= e3;
    //cout << "temp3:" << temp << endl;
	temp /= e4;
    //cout << "temp4:" << temp << endl;
    
    
    
	for (Big i = 1; i <= 99999; i += 1) {
		Q = pow(g0, i);
		if (Q == temp)  {
            m1dec = i;
            break;
		}
    }
	
}

void Eval(Ciphtertext& CT1, Ciphtertext& CT2, Ciphtertext& CT3) {
	CT3.c1 = CT1.c1 * CT2.c1;
	CT3.c2 = CT1.c2 + CT2.c2;
	CT3.c3 = CT1.c3 + CT2.c3;
	CT3.c4 = CT1.c4 + CT2.c4;
    CT3.c5 = CT1.c5 + CT2.c5;
}


void print_centered_divider(const std::string& content, int total_length = 80) {
    int content_length = content.length();

    if (content_length > total_length) {
        std::cout << content << std::endl;
        return;
    }

    int left_padding = (total_length - content_length) / 2;
    int right_padding = total_length - content_length - left_padding;

    std::string divider;
    divider.append(left_padding, '-');
    divider.append(content);
    divider.append(right_padding, '-');

    std::cout << divider << std::endl;
}





/*---------------------------------------------------------���Hash��������ĺ���-----------------------------------------------------*/
/*---------------------------------------------------------------------------------------------------------------------������---------------------------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------------------------------------------������---------------------------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------------------------------------------������---------------------------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------------------------------------------������---------------------------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------------------------------------------������---------------------------------------------------------------------------------------------------*/
int main()
{
    //ofstream common("common.ibe");
    //ofstream master("master.ibe");
    //ECn P,Ppub;
    ECn P, Q, R, Ppub, Qid;                    //Ppub:Master Public Key, Qid:Identity Public Key P:����Ԫ��Q:����㣬R:�м����
	ECn QidAlice, QidBob, QTimeToBeDec, QTimeNow;                       //QidAlice:Identity Public Key Alice, QidBob:Identity Public Key Bob
    //ZZn px,py;//�Լӣ�������ǵ��x,y����
    ZZn2 g0, Qx, Qy, gid, gid1, cube, w;//�Լ�     //gid:˫���ԶԽ����gid1:˫���ԶԽ����cube:��������w:�������� Qx,Qy:��Բ�����ϵĵ�Q������
   // Big xx,yy,ab,r1;//�Լ�
    Big px, py, qx, qy, ab, r1;//�Լ�         //px,py:��Բ�����ϵĵ�P�����꣬qx,qy:��Բ�����ϵĵ�Q�����꣬ab:��Բ�����ϵĵ�Qid�����꣬r1:��ϣֵ
    //int i, tag;//�Լ�                         //i:ѭ��������tag:��־λ
   // float t1,t_sub,t2,t_div,t3,t4,t5,t_ecsub,t6,t7,t8,t9,t10,t11,t_div_G2,t_comp;//�Լ�
    //float t1, t_sub, t2, t_div, t3, t4, t5, t_ecsub, t6, t7, t8, t88, t9, t10, t11, t_div_G2, t_comp, t_G1_xor;//�Լ�  //:�����ף�t1:��ϣֵ����ʱ�䣬t_sub:�������ʱ�䣬t2:�������ʱ�䣬t_div:ģ������ʱ�䣬t3:ģ������ʱ�䣬t4:�������ʱ�䣬t5:�������ʱ�䣬t_ecsub:�������ʱ�䣬t6:˫���Զ�����ʱ�䣬t7:������ʱ�䣬t8:map_to_point,��H2:G2--->G1����ʱ�䣬t9:H2:G2--->{0,1}^logp����ʱ�䣬t10:H:{0,1}^*--->Z^*_q����ʱ�䣬t11:ģ������in G2��gid*gid1 ����ʱ�䣬t_div_G2:ģ������in G2��gid1/gid ����ʱ�䣬t_comp:����ZZn2�ϵ������Ƿ���ģ����
   // ZZn2 cube;
   // Big s,p,q,t,n,cof,x,y;
    Big a, b, c, d1, d2, p, q, t, n, cof, x, y;  //a:˽Կ��b:˽Կ��c:˽Կ��d1:˽Կ��d2:˽Կ��p:����p��q:����q��t:�����ף�n:�������cof:ϵ����x:��Բ�����ϵĵ�����꣬y:��Բ�����ϵĵ������
    big x1, y1, x2, y2;                          //x1,y1:��Բ�����ϵĵ�����꣬x2,y2:��Բ�����ϵĵ������
    long seed;                                   //seed:���������
    char pad[HASH_LEN1];//�Լ�                   //pad:��ϣֵ
    //char pad[20]={0};

    Big m000; ZZn2 g111; ZZn2 g222;

    Big m1 = 14562;
    ZZn2 m1test, m1testdec;
    Big m2 = 11111;
	Big m1dec, m2dec, m3dec;
    params params;
	sk_pkg sk_pkg;
	sk_ts sk_ts;

	char Alice[9] = "Alice";
	char Bob[100] = "Bob@qq.com";
	char TimeToBeDec[9] = "22222222";
    char TimeNow[9] = "22222222";
    Ciphtertext CT1;
    Ciphtertext CT2;
    Ciphtertext CT3;

    PrivateKey sk_Alice;
    PrivateKey sk_Bob;
	TimetrapDoor st_TimeNow;

    miracl* mip = &precision;                    //miracl* mip:����

    cout << "������Щ����������ʱ����1���룬���л������������ظ�ִ��" << renum << "�� " << endl;
    cout << "Enter 9 digit random number seed  = ";
    cin >> seed;
    irand(seed);

    cout << "\n" << endl;
    print_centered_divider("System Initialization");

    // SET-UP
    /*-------------------------------------------------------------����������q-------------------------------------------------------------------------*/
    q = pow((Big)2, 159) + pow((Big)2, 17) + 1; 
    cout<<"����������q��"<<endl;
    cout << "q= " << q << endl;

/*--------------------------------------------------------------��������p-------------------------------------------------------------------------*/
    t = (pow((Big)2, PBITS) - 1) / (2 * q);
    a = (pow((Big)2, PBITS - 1) - 1) / (2 * q);
    forever
    {
        n = rand(t);
        if (n < a) continue;
        p = 2 * n * q - 1;
        if (p % 24 != 11) continue;  // must be 2 mod 3, also 3 mod 8
        if (prime(p)) break;
    }
	cout << "����ģ��p:" << endl;
    cout << "p= " << p << endl;


    cof = 2 * n;  //��Բ����������

    ecurve(0, 1, p, MR_PROJECTIVE);    // elliptic curve y^2=x^3+1 mod p����Ӱ����ϵͳ

   
    forever
    {
            cube = pow(randn2(),(p + 1) / 3); 
            cube = pow(cube,p - 1);
            if (!cube.isunity()) break;  
    }

    cout << "Cube root of unity= " << cube << endl;

    if (!(cube * cube * cube).isunity())
    {
        cout << "sanity check failed" << endl;
        exit(0);
    }

    sk_pkg.alpha = rand(q);
    sk_ts.beta = rand(q);
    cout<< "����PKG˽Կalpha:" << endl;
	cout << "alpha= " << sk_pkg.alpha << endl;
    cout << "����TS˽Կbeta:" << endl;
	cout << "beta= " << sk_ts.beta << endl;
    /*---------------------------------------------------------���ɹ�������params----------------------------------------------*/
    print_centered_divider("��������params����");
  
    //������Բ����������ԪP(g)
    forever
    {
        while (!params.P.set(randn()));
        params.P *= cof;
        if (!params.P.iszero()) break;
    }
        cout<< "������Բ����������Ԫg:" << endl;
        cout<< "g:" << params.P << endl;

    //������Բ����������ԪP1(g1)
	params.P1 = sk_pkg.alpha * params.P;
    cout << "������Բ����������Ԫg1:" << endl;
    cout<<"g1:" << params.P1 << endl;

    //������Բ����������ԪP2(g2)
    forever
    {
        while (!params.P2.set(randn()));
        params.P2 *= cof;
        if (!params.P2.iszero()) {
			cout << "685" << params.P2 << endl;
            break;
        }
    }
        cout << "������Բ����������Ԫg2:" << endl;
        cout << "g2:" << params.P2 << endl;

    //������Բ����������ԪP3(g3)
    params.P3 = sk_ts.beta * params.P;
    cout << "������Բ����������Ԫg3:" << endl;
    cout << "g3:" << params.P3 << endl;

    //������Բ����������ԪP4(g4)
    forever
    {
        while (!params.P4.set(randn()));
        params.P4 *= cof;
        if (!params.P4.iszero()) break;
    }
        cout << "������Բ����������Ԫg4:" << endl;
		cout << "g4:" << params.P4 << endl;

    //������Բ������u[257]
    for (int i = 0; i < 257; i++) {
        forever
        {
            while (!params.u[i].set(randn()));
            params.u[i] *= cof;
            if (!params.u[i].iszero()) break;
        }
    }
    cout<< "��Բ������u[257]�����ɹ�" << endl;

    //������Բ������v[257]
    for (int i = 0; i < 257; i++) {
        forever
        {
            while (!params.v[i].set(randn()));
            params.v[i] *= cof;
            if (!params.v[i].iszero()) break;
        }
    }
    cout<< "��Բ������v[257]�����ɹ�" << endl;

    //������Բ������˫����ӳ�� e(g1,g2) �� e(g3,g4)
    ecap(params.P1, params.P2, q, cube, params.e_g1_g2);
    if (ecap(params.P1, params.P2, q, cube, params.e_g1_g2)) {
        cout << "ecap:e_g1_g2 success" << endl;
    }
    else {
        cout << "ecap:e_g1_g2 error" << endl;
    }
    cout << "��Բ������˫����ӳ�� e(g1, g2)�����ɹ�" << endl;
    cout << "e(g1,g2):" << params.e_g1_g2 << endl;
    


    //ecap(params.P3, params.P4, q, cube, params.e_g3_g4);
    if (ecap(params.P3, params.P4, q, cube, params.e_g3_g4)) {
        cout << "ecap:e_g3_g4 success" << endl;
	}
	else {
		cout << "ecap:e_g3_g4 error" << endl;
	}
    cout<< "��Բ������˫����ӳ�� e(g3, g4)�����ɹ�" << endl;
    cout << "e(g3,g4):" << params.e_g3_g4 << endl;
    
	/*------------------------------------------------------��������params���� ���-------------------------------------------*/


    /*
    m000 = 123; 
    g111 = randn2(); 
    g222 = pow(g111, m000);
    for (int i = 1; i < 1000; i + 1) {
        if (g222 == pow(g111, m000)) {
			cout << "795" << i << endl;
            break;
        }
    }
    */




    //g0 = randn2();
	//ZZn2 g00 = randn2();
    //ZZn2 g1 = randn2();
    //ZZn2 g2 = g0 * g1 / g00;
	//cout << "g0:" << g0 << endl;
    //cout << "g00:" << g00 << endl;
    //cout << "g1:" << g1 << endl;
    //cout << "g2:" << g2 << endl;
    //cout << "g0:" <<  g2 / g1 * g00 <<endl;
    //cout << "g0:" <<  g2 * g00 / g1 << endl;


    /*-------------------------------------------���ɣ�PKG��˽Կ��TS��˽Կ��Alice ˽Կ��Bob˽Կ��ʱ������ST----------------------------------------------*/
    print_centered_divider("˽Կ����");
   
    sk_pkg.sk_pkg = sk_pkg.alpha * params.P2;
    cout<< "PKG��˽Կ���ɳɹ�" << endl;
    cout<< "PKG��˽Կsk_pkg:" << sk_pkg.sk_pkg << endl;
	sk_ts.sk_ts = sk_ts.beta * params.P4;
    cout<< "TS��˽Կ���ɳɹ�" << endl;
	cout<< "TS��˽Կsk_ts:" << sk_ts.sk_ts << endl;

	QidAlice = map_to_point(Alice);
	cout << "Alice�����˽Կ���ɳɹ�" << endl;
    cout << "Alice�����˽ԿQidAlice:" << QidAlice << endl;
	QidBob = map_to_point(Bob);
    cout << "Bob�����˽Կ���ɳɹ�" << endl;
    cout<<"Bob�����˽ԿQidBob:" << QidBob << endl;

	//QTimeNow = map_to_point(TimeNow);
	QTimeToBeDec = map_to_point(TimeToBeDec);
    cout << "ʱ������ST���˽Կ���ɳɹ�" << endl;
    cout << "ʱ������ST�����˽ԿQTimeToBeDec:" << QTimeToBeDec << endl;

    //SkGen(params, sk_pkg, p, QidAlice, sk_Alice);
    SkGen(params, sk_pkg, q, QidBob, sk_Bob);

	StGen(params, sk_ts, q, QTimeToBeDec, st_TimeNow);

    g0 = randn2();
    //m1test = randn2();
	//cout << "772-m1test:" << m1test << endl;


    //cout << "766-QID" << QidBob << endl;
    print_centered_divider("�ӽ�����֤");
    
    cout << "841:m1:::" << m1 << endl;
    cout << "842:m2:::" << m2 << endl;
    Enc(g0, m1, params, QidBob, QTimeToBeDec, q, CT1);
    Enc(g0, m2, params, QidBob, QTimeToBeDec, q, CT2);
    //Enc(m1test, params, QidBob, QTimeToBeDec, q, CT1
    //Enc(m1test, params, QidBob, q, CT1);
    cout << "���ܳɹ�!" << endl;

    Dec(CT1, sk_Bob, st_TimeNow, m1dec, cube, g0 ,q);
	cout << "848:m1dec" << m1dec << endl;
    Dec(CT2, sk_Bob, st_TimeNow, m2dec, cube, g0, q);
    cout << "850:m2dec" << m2dec << endl;
    cout<< "���ܳɹ�!" << endl;
  
    print_centered_divider("̬ͬ����֤");
    Eval(CT1, CT2, CT3);
    Dec(CT3, sk_Bob, st_TimeNow, m3dec, cube, g0, q);
    cout << "861:(m1+m2)dec" << m3dec << endl;
    cout<<"̬ͬ����֤�ɹ�"<<endl;
    //Dec(CT1, sk_Bob, m1testdec, cube, q);

	//cout << "m1testdec:" << m1testdec << endl;
    /*
    if (ZZn2_compare(m1testdec, m1test)) {
        cout << "sucess" << endl;
    }
    else {
        cout << "error" << endl;
    }
    */
    //cout << "m2����ǰ�� " << m1 << endl;
    //Dec(CT2, sk_Bob, st_TimeNow, m2dec, q, cube, g0);

    return 0;

    //int Readkey();
    //int aa;
    //cin >> aa;
}


	