/*
 *    MIRACL  C++  functions ecn.cpp
 *
 *    AUTHOR  :    M. Scott
 *             
 *    PURPOSE :    Implementation of class ECn functions using Montgomery
 *                 representation
 *    NOTE    :    Must be used in conjunction with big.h and big.cpp
 *
 *    Copyright (c) 1988-2004 Shamus Software Ltd.
 */

#include "ecn.h"

int ECn::get(Big& x,Big& y) const 
        {return epoint_get(p,x.getbig(),y.getbig());}
int ECn::get(Big& x) const   
        {return epoint_get(p,x.getbig(),x.getbig());}
void ECn::getx(Big &x) const
        {epoint_getxyz(p,x.getbig(),NULL,NULL);}
void ECn::getxy(Big &x,Big &y) const
        {epoint_getxyz(p,x.getbig(),y.getbig(),NULL);}
void ECn::getxyz(Big &x,Big &y, Big &z) const
        {epoint_getxyz(p,x.getbig(),y.getbig(),z.getbig());}


BOOL ECn::iszero() const
        {if (p->marker==MR_EPOINT_INFINITY) return TRUE; return FALSE;}

epoint * ECn::get_point() const
{ return p; }

ECn operator-(const ECn& e)
{ ECn t=e; epoint_negate(t.p); return t;}

ECn mul(const Big& e1,const ECn& p1,const Big& e2,const ECn& p2)
{ 
    ECn t; 
    ecurve_mult2(e1.getbig(),p1.get_point(),e2.getbig(),p2.get_point(),t.get_point());
    return t;
}

ECn operator*(const Big& e,const ECn& b)
{
    ECn t;
    ecurve_mult(e.getbig(),b.p,t.p);
    return t;
}
/*---------yk���ӵ�������������:���ڵ��������Ч�ʲ��ߣ��Ҵ���ʵ�ֱ�ʾ�Ĳ����㣨��������ĵ�ܿ��ܲ�����Բ�����ϣ�����ϵͳ����ֵ��Ϊ��0,0�������޸ķ����е����Ϊ�������,����ʹ�øò����������ǲ�ǿ���ã�0,0���Ĵ������⣬PBC���Ƿ����㣿-------------------------------------------------*/
//ECn operator^(const ECn& a,const ECn& b)
//void ecn_xor(const ECn& a,const ECn& b,Big& c, Big& d)
void ecn_xor(Big& px,Big& py,Big& qx,Big& qy)
{
	int i;
/*����������������������������1��ʹ�ý϶�ռ�,��ʱ����0.00179��������������������������������������*/
	char ptrpx[64]={0},ptrpy[64]={0},ptrqx[64]={0},ptrqy[64]={0};
	BOOL justify=1;
	big_to_bytes(64,px.getbig(),ptrpx,justify);//������ת��Ϊ������
	big_to_bytes(64,py.getbig(),ptrpy,justify);
	big_to_bytes(64,qx.getbig(),ptrqx,justify);
	big_to_bytes(64,qy.getbig(),ptrqy,justify);

	for(i=0;i<64;i++)//ʵ�����
		{
			ptrqx[i]=ptrpx[i]^ptrqx[i];//�����Ƿ�ʵ������򣬾����ԣ��ɹ�ʵ�������
			ptrqy[i]=ptrpy[i]^ptrqy[i];
	    }
	bytes_to_big(64,ptrqx,qx.getbig());//ʵ�ֶ�����ת��Ϊ������
	bytes_to_big(64,ptrqy,qy.getbig());//ʵ�ֶ�����ת��Ϊ������*/
/*����������������������������2��ʹ�ý��ٿռ�,��ʱ����0.00179��������������������������������������*/
	/*char ptr1[64]={0},ptr2[64]={0};
	BOOL justify=1;
	big_to_bytes(64,px.getbig(),ptr1,justify);//������ת��Ϊ������
	big_to_bytes(64,qx.getbig(),ptr2,justify);
	for(i=0;i<64;i++)//ʵ�����
		ptr2[i]=ptr1[i]^ptr2[i];//�����Ƿ�ʵ������򣬾����ԣ��ɹ�ʵ�������
	bytes_to_big(64,ptr2,qx.getbig());//ʵ�ֶ�����ת��Ϊ������

    big_to_bytes(64,py.getbig(),ptr1,justify);//������ת��Ϊ������
	big_to_bytes(64,qy.getbig(),ptr2,justify);
	for(i=0;i<64;i++)//ʵ�����
		ptr2[i]=ptr1[i]^ptr2[i];//�����Ƿ�ʵ������򣬾����ԣ��ɹ�ʵ�������
	bytes_to_big(64,ptr2,qy.getbig());//ʵ�ֶ�����ת��Ϊ������*/







	
	/*int i;
	big px,py,qx,qy,tx,ty;
	char ptrpx[64]={0},ptrpy[64]={0},ptrqx[64]={0},ptrqy[64]={0},ptrtx[64]={0},ptrty[64]={0};
	BOOL justify=1;
	px=mirvar(0);  // all big variables need to be "mirvar"ed 
    py=mirvar(0);
	qx=mirvar(0);  // all big variables need to be "mirvar"ed 
    qy=mirvar(0);
	tx=mirvar(0);  // all big variables need to be "mirvar"ed 
    ty=mirvar(0);

	epoint_get(a.get_point(),px,py);//ȡ����
	epoint_get(b.get_point(),qx,qy);

	//cotnum(px,stdout);//�������������
	//cotnum(py,stdout);//
	//cotnum(qx,stdout);//�������������
	//cotnum(qy,stdout);//

	big_to_bytes(64,px,ptrpx,justify);//������ת��Ϊ������
	big_to_bytes(64,py,ptrpy,justify);
	big_to_bytes(64,qx,ptrqx,justify);
	big_to_bytes(64,qy,ptrqy,justify);

	/*bytes_to_big(64,ptrpx,px);//ʵ�ֶ�����ת��Ϊ������
	bytes_to_big(64,ptrpy,py);
	bytes_to_big(64,ptrqx,qx);//ʵ�ֶ�����ת��Ϊ������
	bytes_to_big(64,ptrqy,qy);

	cotnum(px,stdout);//�������������
	cotnum(py,stdout);//
	cotnum(qx,stdout);//�������������
	cotnum(qy,stdout);//*/

	
	//bytes_to_big(64,ptrpx,px);//�����Ƿ�ʵ������򣬾����ԣ��ɹ�ʵ�������
	//bytes_to_big(64,ptrpy,py);//�����Ƿ�ʵ������򣬾����ԣ��ɹ�ʵ�������

	//cotnum(tx,stdout);//�������������
	//cotnum(ty,stdout);//
	//cotnum(px,stdout);//�����Ƿ�ʵ������򣬾����ԣ��ɹ�ʵ�������
	//cotnum(py,stdout);//�����Ƿ�ʵ������򣬾����ԣ��ɹ�ʵ�������

	//c=tx;
	//d=ty;
	///cotnum(c.getbig(),stdout);//�������������
	//cotnum(d.getbig(),stdout);//

	//epoint_set(tx,ty,1,t.get_point());//ʵ��Ϊ�����긳ֵ
	//epoint_get(t.get_point(),tx,ty);
	//big_to_bytes(64,tx,ptrtx,justify);
	//big_to_bytes(64,ty,ptrty,justify);
	//for (i=0;i<64;i++) printf("%02x",ptrtx[i]);
	//printf("\n");
	//for (i=0;i<64;i++) printf("%02x",ptrty[i]);
	//printf("\n");
	//return t;//����t�������������Ľ������ˣ�t���ܲ�����Բ�����ϣ��������Ϊ��Infinity����Ϊ��������Զ�㡣
 
}
/*---------------------------------------------------------------------�Լ����ӵ�������������-------------------------------------------------------------------------------------------------*/
/*ECn operator^(const ECn& a,const ECn& b)
//void ecn_xor(const ECn& a,const ECn& b,Big& c, Big& d)
{
	ECn t;
    int i;
	big px,py,qx,qy,tx,ty;
	char ptrpx[64]={0},ptrpy[64]={0},ptrqx[64]={0},ptrqy[64]={0},ptrtx[64]={0},ptrty[64]={0};
	BOOL justify=1;
	px=mirvar(0);  // all big variables need to be "mirvar"ed 
    py=mirvar(0);
	qx=mirvar(0);  // all big variables need to be "mirvar"ed 
    qy=mirvar(0);
	tx=mirvar(0);  // all big variables need to be "mirvar"ed 
    ty=mirvar(0);

	epoint_get(a.get_point(),px,py);//ȡ����
	epoint_get(b.get_point(),qx,qy);

	cotnum(px,stdout);//�������������
	cotnum(py,stdout);//
	cotnum(qx,stdout);//�������������
	cotnum(qy,stdout);//

	big_to_bytes(64,px,ptrpx,justify);//������ת��Ϊ������
	big_to_bytes(64,py,ptrpy,justify);
	big_to_bytes(64,qx,ptrqx,justify);
	big_to_bytes(64,qy,ptrqy,justify);

	/*bytes_to_big(64,ptrpx,px);//ʵ�ֶ�����ת��Ϊ������
	bytes_to_big(64,ptrpy,py);
	bytes_to_big(64,ptrqx,qx);//ʵ�ֶ�����ת��Ϊ������
	bytes_to_big(64,ptrqy,qy);

	cotnum(px,stdout);//�������������
	cotnum(py,stdout);//
	cotnum(qx,stdout);//�������������
	cotnum(qy,stdout);//*/

	/*for(i=0;i<64;i++)//ʵ�����
		{
			ptrtx[i]=ptrpx[i]^ptrqx[i];
			ptrty[i]=ptrpy[i]^ptrqy[i];
			//ptrtx[i]=ptrpx[i]^0;
			//ptrty[i]=ptrpy[i]^0;
	    }
	/*for (i=0;i<64;i++) printf("%02x",ptrtx[i]);
	printf("\n");
	for (i=0;i<64;i++) printf("%02x",ptrty[i]);
	printf("\n");
	bytes_to_big(64,ptrtx,tx);//ʵ�ֶ�����ת��Ϊ������
	bytes_to_big(64,ptrty,ty);
	big_to_bytes(64,tx,ptrtx,justify);
	big_to_bytes(64,ty,ptrty,justify);
	for (i=0;i<64;i++) printf("%02x",ptrtx[i]);
	printf("\n");
	for (i=0;i<64;i++) printf("%02x",ptrty[i]);
	printf("\n");*/

	/*bytes_to_big(64,ptrtx,tx);//ʵ�ֶ�����ת��Ϊ������
	bytes_to_big(64,ptrty,ty);
	

	cotnum(tx,stdout);//�������������
	cotnum(ty,stdout);//

	epoint_set(tx,ty,1,t.get_point());//ʵ��Ϊ�����긳ֵ
	//epoint_get(t.get_point(),tx,ty);
	//big_to_bytes(64,tx,ptrtx,justify);
	//big_to_bytes(64,ty,ptrty,justify);
	//for (i=0;i<64;i++) printf("%02x",ptrtx[i]);
	//printf("\n");
	//for (i=0;i<64;i++) printf("%02x",ptrty[i]);
	//printf("\n");
	return t;//����t�������������Ľ������ˣ�t���ܲ�����Բ�����ϣ��������Ϊ��Infinity����Ϊ��������Զ�㡣
	
	
	
	
	/* ECn t;
	int i;
	big px,py,qx,qy,tx,ty;
	char ptrpx[64],ptrpy[64],ptrqx[64],ptrqy[64],ptrtx[64],ptrty[64];
	BOOL justify=1;
	px=mirvar(0);  // all big variables need to be "mirvar"ed 
    py=mirvar(0);
	qx=mirvar(0);  // all big variables need to be "mirvar"ed 
    qy=mirvar(0);
	tx=mirvar(0);  // all big variables need to be "mirvar"ed 
    ty=mirvar(0);

	epoint_get(a.get_point(),px,py);//ȡ����
	epoint_get(b.get_point(),qx,qy);

	big_to_bytes(64,px,ptrpx,justify);//������ת��Ϊ������
	big_to_bytes(64,py,ptrpy,justify);
	big_to_bytes(64,qx,ptrqx,justify);
	big_to_bytes(64,qy,ptrqy,justify);

	for(i=0;i<64;i++)//ʵ�����
		{
			ptrtx[i]=ptrpx[i]^ptrqx[i];
			ptrty[i]=ptrpy[i]^ptrqy[i];
	    }
	for (i=0;i<64;i++) printf("%02x",ptrtx[i]);
	printf("\n");
	for (i=0;i<64;i++) printf("%02x",ptrty[i]);
	printf("\n");
	bytes_to_big(64,ptrtx,tx);//ʵ�ֶ�����ת��Ϊ������
	bytes_to_big(64,ptrty,ty);
	big_to_bytes(64,tx,ptrtx,justify);
	big_to_bytes(64,ty,ptrty,justify);
	for (i=0;i<64;i++) printf("%02x",ptrtx[i]);
	printf("\n");
	for (i=0;i<64;i++) printf("%02x",ptrty[i]);
	printf("\n");
	bytes_to_big(64,ptrtx,tx);//ʵ�ֶ�����ת��Ϊ������
	bytes_to_big(64,ptrty,ty);

	epoint_set(tx,ty,1,t.get_point());//ʵ��Ϊ�����긳ֵ
	epoint_get(t.get_point(),tx,ty);
	big_to_bytes(64,tx,ptrtx,justify);
	big_to_bytes(64,ty,ptrty,justify);
	for (i=0;i<64;i++) printf("%02x",ptrtx[i]);
	printf("\n");
	for (i=0;i<64;i++) printf("%02x",ptrty[i]);
	printf("\n");
	return t;*/
	/*epoint_get(a.get_point(),px,py);//������Ա�����epoint_get(a.get_point(),px,py)��ʵ����ȡ������ĺ������������ǣ�px=(a.get_point())->X;py=(a.get_point())->Y;
	//px=(a.get_point())->X;//����ʵ��ȡ��������
	//py=(a.get_point())->Y;//����ʵ��ȡ��������
	big_to_bytes(64,px,ptrpx,justify);
	big_to_bytes(64,py,ptrpy,justify);
	for (i=0;i<64;i++) printf("%02x",ptrpx[i]);
	printf("\n");
	for (i=0;i<64;i++) printf("%02x",ptrpy[i]);
	printf("\n");
	epoint_set(px,py,1,a.get_point());*/

    /*px=(a.get_point())->X;
	py=(a.get_point())->Y;
	qx=(b.get_point())->X;
	qy=(b.get_point())->Y;
	big_to_bytes(64,px,ptrpx,justify);
	//for (i=0;i<64;i++) printf("%02x",ptrpx[i]);
   // printf("\n");
	big_to_bytes(64,py,ptrpy,justify);
	big_to_bytes(64,qx,ptrqx,justify);
	big_to_bytes(64,qy,ptrqy,justify);

	for (i=0;i<64;i++) printf("%02x",ptrpx[i]);
	printf("\n");
	for (i=0;i<64;i++) printf("%02x",ptrpy[i]);
	printf("\n");
	for (i=0;i<64;i++) printf("%02x",ptrqx[i]);
	printf("\n");
	for (i=0;i<64;i++) printf("%02x",ptrqy[i]);
	printf("\n");


		tx=(t.get_point())->X;
		ty=(t.get_point())->Y;
	big_to_bytes(64,tx,ptrtx,justify);
	big_to_bytes(64,ty,ptrty,justify);
	for (i=0;i<64;i++) printf("%02x",ptrtx[i]);
	printf("\n");
	for (i=0;i<64;i++) printf("%02x",ptrty[i]);
	printf("\n");


	for(i=0;i<64;i++)
		{
			*(ptrtx+i)=*(ptrpx+i)^*(ptrqx+i);
			*(ptrty+i)=*(ptrpy+i)^*(ptrqy+i);
	    }
	for (i=0;i<64;i++) printf("%02x",ptrtx[i]);
    printf("\n");
	bytes_to_big(64,ptrtx,tx);
	bytes_to_big(64,ptrty,ty);
	//(t.get_point())->X=tx;
	//(t.get_point())->Y=ty;
	//t.set(tx,ty);
	tx=(t.get_point())->X;
	big_to_bytes(64,tx,ptrtx,justify);
	for (i=0;i<64;i++) printf("%02x",ptrtx[i]);//
	printf("\n");
	ty=(t.get_point())->Y;
	big_to_bytes(64,ty,ptrty,justify);
	for (i=0;i<64;i++) printf("%02x",ptrty[i]);//

    printf("\n");
   
}*/
/*---------------------------------------------------------------------yk���ӵ�������������-------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------------------  ---yk���ӵĵ�Ӻ���-------------------------------------------------------------------------------------------------*/
 ECn operator+(const ECn& b1, const ECn& b2)
{ ecurve_add(b1.p,b2.p); return  b2;}
/*------------------------------------------------------------------  ---yk���ӵĵ�Ӻ���-------------------------------------------------------------------------------------------------*/
 /*------------------------------------------------------------------  ---yk���ӵĵ������-------------------------------------------------------------------------------------------------*/
 ECn operator-(const ECn& b1, const ECn& b2)
{ ecurve_sub(b2.p,b1.p); return  b1;}
/*------------------------------------------------------------------  ---yk���ӵĵ������-------------------------------------------------------------------------------------------------*/

#ifndef MR_STATIC

ECn mul(int n,const Big *y,ECn *x)
{
    ECn w;
    int i;
    big *a=(big *)mr_alloc(n,sizeof(big));
    epoint **b=(epoint **)mr_alloc(n,sizeof(epoint *));
    for (i=0;i<n;i++)
    {
        a[i]=y[i].getbig();
        b[i]=x[i].p;
    }
    ecurve_multn(n,a,b,w.p);

    mr_free(b);
    mr_free(a);
    return w;  
}

void multi_add(int m,ECn *x, ECn *w)
{
    int i;
    epoint **xp=(epoint **)mr_alloc(m,sizeof(epoint *));
    epoint **wp=(epoint **)mr_alloc(m,sizeof(epoint *));
    for (i=0;i<m;i++)
    {
        xp[i]=x[i].p;
        wp[i]=w[i].p;
    }
    ecurve_multi_add(m,xp,wp);
    mr_free(wp);
    mr_free(xp);
}

#endif

void double_add(ECn& A,ECn& B,ECn& C,ECn& D,big& s1,big& s2)
{
    ecurve_double_add(A.p,B.p,C.p,D.p,&s1,&s2);
}

#ifndef MR_NO_STANDARD_IO

ostream& operator<<(ostream& s,const ECn& b)
{
    Big x,y;
    if (b.iszero())
        s << "(Infinity)";
    else
    {
        b.get(x,y);
        s << "(" << x << "," << y << ")";
    }
    return s;
}

#endif

void ecurve(const Big& a,const Big& b,const Big& p,int t)
{
    ecurve_init(a.fn,b.fn,p.fn,t);
}

#ifndef MR_NOSUPPORT_COMPRESSION
#ifndef MR_NOTESTXONCURVE


BOOL is_on_curve(const Big& a)
{ return epoint_x(a.fn);}

#endif
#endif
