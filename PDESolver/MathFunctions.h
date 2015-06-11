#ifndef MATH1
#define MATH1

#include "knotSet.h"
#include "matrix3D.h"
#include "fmatrix.h"


class Math {

public:

	// number of combinations of k from n

static double Combin(int n,int k)
	{
		int j;
 
		double comb;
		if(k>(n/2))
			j=n-k;
		else
			j=k;
		if(j<0)
			comb=0; // changed from 1
		else if(j==0)
			comb=1;
		else 
		{
			comb=1;
			for(int i=0;i<=j-1;i++)
				comb=(comb*(double)(n-i))/(double)(j-i);
		}
		return comb;
	}

	// cyles through all permutations of r things from n, returning permutation vector each time

	 static int Combin1(int n,int r,Vector<int>& iarray,int first)
	{
		if(first==1)
		{
			for(int j=0;j<r;j++)
				iarray[j] = j+1;
			first=0;
			return first;
		}
		if(iarray[r-1]<n)
		{
			iarray[r-1] = iarray[r-1]+1;
			return first;
		}
		for(int j=r-1;j>=1;j--)
		{
			if(iarray[j-1]<n-r+j)
			{
				iarray[j-1] = iarray[j-1]+1;
				for(int k=j;k<r;k++)
					iarray[k] = iarray[j-1]+k-j+1;
				return first;
			}
		}
		first=1;
		return first;
	}               

	// factorial n

	 static double Fact(int n)
	{
		double fact;
   
		if(n<0)
			fact=0; // changed from 1
		else if(n==0)
			fact=1;
		else 
		{
			fact=1;
			for(int i=0;i<=n-1;i++)
				fact=fact*(n-i);
		}
		return fact;
	}
 

/************************* <<  gasdev  >> ***********************************/
/* Returns a Gaussian distributed deviate with zero mean and unit variance, */
/* using RANG(IDUM) as the source of uniform deviates.                      */
/*                                                                          */
/* (From 'Numerical Recipes' by Press et al. 1987)                          */
/****************************************************************************/

	 static double gasdev(int *idum, int iset)
	{
		double gset;
		double v1, v2, r, fac;

		if (iset == 0)
		{
			do
			{
				v1 = 2.0 * rang(idum,iset) - 1.0;
				v2 = 2.0 * rang(idum,iset) - 1.0;

				r = pow( v1, 2 ) + pow( v2, 2 );
			}
			while (r >= 1.0);

			fac = sqrt( (double) -2.0 * log( r ) / r );
			gset = v1 * fac;

			iset = 1;
			return( v2 * fac );
		} else {
			iset = 0;
			return( gset );
		}

	}


/************************ <<  rang  >> **************************************/
/* - This function returns a uniform random number between 0. and 1.        */
/* - Set IDUM to any negative value to initiate or reinitiate the sequence. */
/*                                                                          */
/*   (From 'Numerical Recipes' by Press et al. 1987)                        */
/****************************************************************************/

	 static double rang(int *idum, int iff)
	{
		double r[97 + 1];
		long ix1, ix2, ix3;
		int i;
		double ri;

		if ((*idum < 0) || (iff == 0))
		{
			iff = 1;

			ix1 = (IC1 - *idum) % M1;
			ix1 = (IA1 * ix1 + IC1) % M1;
			ix2 = ix1 % M2;
			ix1 = (IA1 * ix1 + IC1) % M1;
			ix3 = ix1 % M3;

			for ( i=1; i<=97; i++ )
			{
				ix1 = (IA1 * ix1 + IC1) % M1;
				ix2 = (IA2 * ix2 + IC2) % M2;

				r[i] = ((double) ix1 + (double) ix2 * RM2) * RM1;
			}

		   *idum = 1;
		}

		ix1 = (IA1 * ix1 + IC1) % M1;
		ix2 = (IA2 * ix2 + IC2) % M2;
		ix3 = (IA3 * ix3 + IC3) % M3;

		i = 1 + (97 * ix3) / M3;

		if ((i > 97) || (i < 1))
		{
			printf( "GNOISE: Index outside range!\n" );
			exit( 1 );
		}

		ri = r[i];
		r[i] = ((double) ix1 + (double) ix2 * RM2) * RM1;

		return( ri );
	}


	template<class T>
	 static void swap(T& x, T& y)
	{
		T temp;

		temp = x;
		x = y;
		y = temp;
	}

	// max of two ints
	template<class T>
	static T max1(T a, T b)
	{ 
		if (a > b) return a;
		else return b;
	}

	// min of two ints
	template<class T>
	static T min1(T a, T b)
	{
		if (a > b) return b;
		else return a;
	}

	
	 static double Poisson3D(double x, double y, double z, int num)
	{
		double sum=0.0;
		for (int i=1; i<num; i+=2)
			for (int j=1; j<num; j+=2)
				for (int k=1; k<num; k+=2)
					sum = sum + sin(i*pi*x)*sin(j*pi*y)*sin(k*pi*z)/(double)(i*j*k*(i*i+j*j+k*k));

		return sum*64.0/pow(pi,5.0);
	}


	template<class T>
	 static double Integral(Vol<T>* f)
	{
		FnVolObject<T>* vol = new FnVolObject<T>(f);
		FVector<double> intgrlv(1,3,0.0);
		FVector<double> wgt(1,5);
		FVector<double> snode(1,5), tnode(1,5), unode(1,5);
		Orthonormal o;
		int rndiv = 5;
		int ndiv = 5;
		double ridiv, rjdiv, rkdiv;
		double sstart = f->GetLeftLimitU();
		double send = f->GetRightLimitU();
		double tstart = f->GetLeftLimitV();
		double tend = f->GetRightLimitV();
		double ustart = f->GetLeftLimitW();
		double uend = f->GetRightLimitW();
		double sglbal, tglbal, uglbal;
		o.InitialiseQuadrature();
	
		
		for (int idiv=1; idiv<=ndiv; idiv++) {
			ridiv = (double)(idiv-1);
			for (int inode=1; inode<=5; inode++) {
				tnode[inode] = (ridiv/rndiv)+(o.nnode[inode]/rndiv);
			}
		}
	
		for (int jdiv=1; jdiv<=ndiv; jdiv++) {
			rjdiv = (double)(jdiv-1);

			for (int jnode=1; jnode<=5; jnode++) {
				snode[jnode] = (rjdiv/rndiv)+(o.nnode[jnode]/rndiv);
			}
		}


		for (int kdiv=1; kdiv<=ndiv; kdiv++) {
			rkdiv = (double)(kdiv-1);

			for (int knode=1; knode<=5; knode++) {
				unode[knode] = (rjdiv/rndiv)+(o.nnode[knode]/rndiv);
			}
		}
		
		for (int i=1; i<=5; i++) wgt[i] = 10.0*o.wwgt[i]/rndiv;
	
		for (int knode=1; knode<=5; knode++) {
			uglbal = ustart + o.nnode[knode]*(uend-ustart);
			for (int inode=1; inode<=5; inode++) {
				tglbal = tstart + o.nnode[inode]*(tend-tstart);
				for (int jnode=1; jnode<=5; jnode++) {
					sglbal = sstart+o.nnode[jnode]*(send-sstart);
			
					T fmod = (*vol)(sglbal,tglbal,uglbal);
			
					FVector<double> fmodv = fmod;
				
					for (int i=1; i<=3; i++) { 
						intgrlv[i] = intgrlv[i]+fmodv[i]*wgt[inode]*wgt[jnode]*wgt[knode];
			
					}
				}
			}
		}
	
		double scale = (send-sstart)*(tend-tstart)*(uend-ustart)/8.0;

		return scale*sqrt(intgrlv[1]*intgrlv[1]+intgrlv[2]*intgrlv[2]+intgrlv[3]*intgrlv[3]);
	}



	template<class T>
	 static double Integral(Surf<T>* f)
	{
		FnSurfObject<T>* surf = new FnSurfObject<T>(f);
		FVector<double> intgrlv(1,3,0.0);
		FVector<double> wgt(1,5);
		FVector<double> snode(1,5), tnode(1,5);
		Orthonormal o;
		int rndiv = 5;
		int ndiv = 5;
		double ridiv, rjdiv;
		double sstart = f->GetLeftLimitU();
		double send = f->GetRightLimitU();
		double tstart = f->GetLeftLimitV();
		double tend = f->GetRightLimitV();
		double sglbal, tglbal;
		o.InitialiseQuadrature();


		for (int idiv=1; idiv<=ndiv; idiv++) {
			ridiv = (double)(idiv-1);
			for (int inode=1; inode<=5; inode++) {
				tnode[inode] = (ridiv/rndiv)+(o.nnode[inode]/rndiv);
			}
		}

		for (int jdiv=1; jdiv<=ndiv; jdiv++) {
			rjdiv = (double)(jdiv-1);

			for (int jnode=1; jnode<=5; jnode++) {
				snode[jnode] = (rjdiv/rndiv)+(o.nnode[jnode]/rndiv);
			}
		}
	

		for (int i=1; i<=5; i++) wgt[i] = 10.0*o.wwgt[i]/rndiv;
	

		for (int inode=1; inode<=5; inode++) {
			tglbal = tstart + o.nnode[inode]*(tend-tstart);
			for (int jnode=1; jnode<=5; jnode++) {
				sglbal = sstart+o.nnode[jnode]*(send-sstart);
		
				T fmod = (*surf)(sglbal,tglbal);
			
				FVector<double> fmodv = fmod;
			
				for (int i=1; i<=3; i++) { 
					intgrlv[i] = intgrlv[i]+fmodv[i]*wgt[inode]*wgt[jnode];
			
				}
			}
		}
	
	
		double scale = (send-sstart)*(tend-tstart)/4.0;

		return scale*sqrt(intgrlv[1]*intgrlv[1]+intgrlv[2]*intgrlv[2]+intgrlv[3]*intgrlv[3]);
	}



	template<class T>
	 static double Integral(Curve<T>* c)
	{
		FnCurvObject<T>* curv = new FnCurvObject<T>(c);
		FVector<double> intgrlv(1,3,0.0);
		FVector<double> wgt(1,5);
		FVector<double> snode(1,5);
		Orthonormal o;
		int rndiv = 5;
		int ndiv = 5;
		double ridiv, rjdiv;
		double sstart = c->GetLeftLimit();
		double send = c->GetRightLimit();
		double sglbal;
		o.InitialiseQuadrature();

		for (int jdiv=1; jdiv<=ndiv; jdiv++) {
			rjdiv = (double)(jdiv-1);

			for (int jnode=1; jnode<=5; jnode++) {
				snode[jnode] = (rjdiv/rndiv)+(o.nnode[jnode]/rndiv);
			}
		}

		for (int i=1; i<=5; i++) wgt[i] = 10.0*o.wwgt[i]/rndiv;
		
		for (int jnode=1; jnode<=5; jnode++) {
			sglbal = sstart+o.nnode[jnode]*(send-sstart);
			T fmod = (*curv)(sglbal);
			FVector<double> fmodv = fmod;
			for (int i=1; i<=3; i++) { 
				intgrlv[i] = intgrlv[i]+fmodv[i]*wgt[jnode];
			}
		}
		double scale = (send-sstart)/2.0;

		return scale*sqrt(intgrlv[1]*intgrlv[1]+intgrlv[2]*intgrlv[2]+intgrlv[3]*intgrlv[3]);
	}


	template<class T>
	 static double IntegralSquared(Curve<T>* c)
	{
		FnCurvObject<T>* curv = new FnCurvObject<T>(c);
		FVector<double> intgrlv(1,3,0.0);
		FVector<double> wgt(1,5);
		FVector<double> snode(1,5);
		Orthonormal o;
		int rndiv = 5;
		int ndiv = 5;
		double ridiv, rjdiv;
		double sstart = c->GetLeftLimit();
		double send = c->GetRightLimit();
		double sglbal;
		o.InitialiseQuadrature();

		for (int jdiv=1; jdiv<=ndiv; jdiv++) {
			rjdiv = (double)(jdiv-1);

			for (int jnode=1; jnode<=5; jnode++) {
				snode[jnode] = (rjdiv/rndiv)+(o.nnode[jnode]/rndiv);
			}
		}

		for (int i=1; i<=5; i++) wgt[i] = 10.0*o.wwgt[i]/rndiv;
		
		for (int jnode=1; jnode<=5; jnode++) {
			sglbal = sstart+o.nnode[jnode]*(send-sstart);
			T fmod = (*curv)(sglbal);
			FVector<double> fmodv = fmod;
			for (int i=1; i<=3; i++) { 
				intgrlv[i] = intgrlv[i]+fmodv[i]*fmodv[i]*wgt[jnode];
			}
		}
		double scale = (send-sstart)/2.0;

		return scale*sqrt(intgrlv[1]*intgrlv[1]+intgrlv[2]*intgrlv[2]+intgrlv[3]*intgrlv[3]);
	}


	template<class T>
	 static double IntegralSquared(Surf<T>* f)
	{
		FnSurfObject<T>* surf = new FnSurfObject<T>(f);
		FVector<double> intgrlv(1,3,0.0);
		FVector<double> wgt(1,5);
		FVector<double> snode(1,5), tnode(1,5);
		Orthonormal o;
		int rndiv = 5;
		int ndiv = 5;
		double ridiv, rjdiv;
		double sstart = f->GetLeftLimitU();
		double send = f->GetRightLimitU();
		double tstart = f->GetLeftLimitV();
		double tend = f->GetRightLimitV();
		double sglbal, tglbal;
		o.InitialiseQuadrature();
	
		std::cout << "o.nnode " << o.nnode;

		for (int idiv=1; idiv<=ndiv; idiv++) {
			ridiv = (double)(idiv-1);
			for (int inode=1; inode<=5; inode++) {
				tnode[inode] = (ridiv/rndiv)+(o.nnode[inode]/rndiv);
			}
		}
	
		for (int jdiv=1; jdiv<=ndiv; jdiv++) {
			rjdiv = (double)(jdiv-1);

			for (int jnode=1; jnode<=5; jnode++) {
				snode[jnode] = (rjdiv/rndiv)+(o.nnode[jnode]/rndiv);
			}
		}
		
		for (int i=1; i<=5; i++) wgt[i] = 10.0*o.wwgt[i]/rndiv;
	
		for (int inode=1; inode<=5; inode++) {
			tglbal = tstart + o.nnode[inode]*(tend-tstart);
			for (int jnode=1; jnode<=5; jnode++) {
				sglbal = sstart+o.nnode[jnode]*(send-sstart);
			
				T fmod = (*surf)(sglbal,tglbal);
			
				FVector<double> fmodv = fmod;
			
				for (int i=1; i<=3; i++) { 
					intgrlv[i] = intgrlv[i]+fmodv[i]*fmodv[i]*wgt[inode]*wgt[jnode];
		
				}
			}
		}
	
	
		double scale = (send-sstart)*(tend-tstart)/4.0;

		return scale*sqrt(intgrlv[1]*intgrlv[1]+intgrlv[2]*intgrlv[2]+intgrlv[3]*intgrlv[3]);
	}


	template<class T>
	 static double IntegralSquared(Vol<T>* f)
	{
		FnVolObject<T>* vol = new FnVolObject<T>(f);
		FVector<double> intgrlv(1,3,0.0);
		FVector<double> wgt(1,5);
		FVector<double> snode(1,5), tnode(1,5), unode(1,5);
		Orthonormal o;
		int rndiv = 5;
		int ndiv = 5;
		double ridiv, rjdiv, rkdiv;
		double sstart = f->GetLeftLimitU();
		double send = f->GetRightLimitU();
		double tstart = f->GetLeftLimitV();
		double tend = f->GetRightLimitV();
		double ustart = f->GetLeftLimitW();
		double uend = f->GetRightLimitW();
		double sglbal, tglbal, uglbal;
		o.InitialiseQuadrature();
	
	
		for (int idiv=1; idiv<=ndiv; idiv++) {
			ridiv = (double)(idiv-1);
			for (int inode=1; inode<=5; inode++) {
				tnode[inode] = (ridiv/rndiv)+(o.nnode[inode]/rndiv);
			}
		}
	
		for (int jdiv=1; jdiv<=ndiv; jdiv++) {
			rjdiv = (double)(jdiv-1);

			for (int jnode=1; jnode<=5; jnode++) {
				snode[jnode] = (rjdiv/rndiv)+(o.nnode[jnode]/rndiv);
			}
		}


		for (int kdiv=1; kdiv<=ndiv; kdiv++) {
			rkdiv = (double)(kdiv-1);

			for (int knode=1; knode<=5; knode++) {
				unode[knode] = (rjdiv/rndiv)+(o.nnode[knode]/rndiv);
			}
		}
		
		for (int i=1; i<=5; i++) wgt[i] = 10.0*o.wwgt[i]/rndiv;
		
		for (int knode=1; knode<=5; knode++) {
			uglbal = ustart + o.nnode[knode]*(uend-ustart);
			for (int inode=1; inode<=5; inode++) {
				tglbal = tstart + o.nnode[inode]*(tend-tstart);
				for (int jnode=1; jnode<=5; jnode++) {
					sglbal = sstart+o.nnode[jnode]*(send-sstart);
					T fmod = (*vol)(sglbal,tglbal,uglbal);
					FVector<double> fmodv = fmod;
					for (int i=1; i<=3; i++) { 
						intgrlv[i] = intgrlv[i]+fmodv[i]*fmodv[i]*wgt[inode]*wgt[jnode]*wgt[knode];
					}
				}
			}
		}
	
		double scale = (send-sstart)*(tend-tstart)*(uend-ustart)/8.0;

		return scale*sqrt(intgrlv[1]*intgrlv[1]+intgrlv[2]*intgrlv[2]+intgrlv[3]*intgrlv[3]);
	}


	 static double Laplace3D1(double x, double y, double z, int num)
	{
		double sum=0.0;
		for (int i=0; i<num; i++)
			for (int j=0; j<num; j++)
				sum = sum + (sin((2.0*i+1)*pi*x)/(2.0*i+1))*(sin((2.0*j+1)*pi*y)/(2.0*j+1))*sinh(pi*z*sqrt((2.0*i+1)*(2.0*i+1)+(2.0*j+1)*(2.0*j+1)))/sinh(pi*sqrt((2.0*i+1)*(2.0*i+1)+(2.0*j+1)*(2.0*j+1)));
		return sum*16.0/pow(pi,2.0);
	}
	

	 static double Laplace3D2(double x, double y, double z, int num)
	{
		double sum=0.0;
		for (int i=0; i<num; i++)
			for (int j=0; j<num; j++)
				sum = sum + (sin((2.0*i+1)*pi*x)/(2.0*i+1))*(sin((2.0*j+1)*pi*y)/(2.0*j+1))*sinh(pi*z*sqrt((2.0*i+1)*(2.0*i+1)+(2.0*j+1)*(2.0*j+1)))/(sqrt((2.0*i+1)*(2.0*i+1)+(2.0*j+1)*(2.0*j+1))*cosh(pi*sqrt((2.0*i+1)*(2.0*i+1)+(2.0*j+1)*(2.0*j+1))));
		return sum*16.0/pow(pi,3.0);
	}


	 static double Laplace3D1HeatFlowZ0(int num)
	{
		double sum=0.0;
		for (int i=0; i<num; i++)
			for (int j=0; j<num; j++)
				sum = sum + (1.0/(2.0*i+1)*(2.0*i+1))*(1.0/(2.0*j+1)*(2.0*j+1))*sqrt((2.0*i+1)*(2.0*i+1)+(2.0*j+1)*(2.0*j+1))/sinh(pi*sqrt((2.0*i+1)*(2.0*i+1)+(2.0*j+1)*(2.0*j+1)));
		return sum*64.0/pow(pi,3.0);
	}
	

	 static double Laplace3D2HeatFlowZ0(int num)
	{
		double sum=0.0;
		for (int i=0; i<num; i++)
			for (int j=0; j<num; j++)
				sum = sum + (1.0/(2.0*i+1)*(2.0*i+1))*(1.0/(2.0*j+1)*(2.0*j+1))/cosh(pi*sqrt((2.0*i+1)*(2.0*i+1)+(2.0*j+1)*(2.0*j+1)));
		return sum*64.0/pow(pi,4.0);
	}


	 static Point3D Cross(const Point3D& p1, const Point3D& p2)
	{
		return Point3D(p1.GetY()*p2.GetZ()-p1.GetZ()*p2.GetY(),-(p1.GetX()*p2.GetZ()-p1.GetZ()*p2.GetX()),p1.GetX()*p2.GetY()-p1.GetY()*p2.GetX());
	}


	 static double Determinant(const Matrix<double>& mat)
	{
		return mat[0][0]*(mat[1][1]*mat[2][2]-mat[1][2]*mat[2][1])-mat[0][1]*(mat[1][0]*mat[2][2]-mat[2][0]*mat[1][2])+
			mat[0][2]*(mat[1][0]*mat[2][1]-mat[2][0]*mat[1][1]);
	}


	 static Point3D Normalise(const Point3D& p)
	{
		double d = sqrt(p.GetX()*p.GetX()+p.GetY()*p.GetY()+p.GetZ()*p.GetZ());

		return Point3D(p.GetX()/d,p.GetY()/d,p.GetZ()/d);
	}


	 static Matrix<double> ComputeIdentityMatrix(int ord) 
	{
		Matrix<double> mat(ord, ord, 0.0);

		for (int i=0; i<ord; i++) mat[i][i]=1.0;
		
		return mat;
	}


	// matrix to convert poly to Bezier

	 static Matrix<double> ComputePolyCurvMatrix(int ord) 
	{
		Matrix<double> mat(ord, ord, 0.0);

		// matrix from Sherar/Goult book
		for (int i=0; i<ord; i++)
			for (int j=0; j<=i; j++)
				mat[i][j] = pow(-1.0,double(i-j))*Combin(ord-1, i)*Combin(i,j);
		return mat;
	}

	// transpose of above

	 static Matrix<double> ComputePolyCurvMatrixTranspose(int ord) 
	{
		Matrix<double> mat(ord, ord, 0.0);

		// matrix from Sherar/Goult book
		for (int i=0; i<ord; i++)
			for (int j=0; j<=i; j++)
				mat[j][i] = pow(-1.0,double(i-j))*Combin(ord-1, i)*Combin(i,j);
		return mat;
	}

	// computes the matrix inverse form of the Bezier curve

	 static Matrix<double> ComputePolyCurvMatrixInverse(int ord) 
	{
		Matrix<double> mat(ord, ord, 0.0);
		
		for (int i=0; i<ord; i++)
			for (int j=0; j<=i; j++)
				mat[i][j] = Combin(i,j)/Combin(ord-1,j);
		return mat;
	}

	// computes the matrix inverse form of the Bezier curve

	 static Matrix<double> ComputePolyCurvMatrixInverseTranspose(int ord) 
	{
		Matrix<double> mat(ord, ord, 0.0);
		
		for (int i=0; i<ord; i++)
			for (int j=0; j<=i; j++)
				mat[j][i] = Combin(i,j)/Combin(ord-1,j);
		return mat;
	}


	template<class T>
	 static Matrix<T> add(const Matrix<T>& t1, const Matrix<T>& t2)
	{
		if (Matrix<T>::CheckIndicesAdd(t1,t2)) {
			Matrix<T> res(t1.GetNumRows(),t1.GetNumCols());
			for (int i=0; i<t1.GetNumRows(); i++)
					for (int j=0; j<t1.GetNumCols(); j++)
						res[i][j]=t1[i][j]+t2[i][j];
			return res;
		} else return t1;
	}

	
	template<class T>
	 static Matrix3D<T> add(const Matrix3D<T>& t1, const Matrix3D<T>& t2)
	{
			Matrix3D<T> res(t1.GetNumRows(),t1.GetNumCols(),t1.GetNumPols());
			for (int k=0; k<t1.GetNumPols(); k++)
				for (int i=0; i<t1.GetNumRows(); i++)
					for (int j=0; j<t1.GetNumCols(); j++)		
						res[k][i][j]=t1[k][i][j]+t2[k][i][j];
			return res;
		//} else return t1;
	}


	template<class T>
	 static Matrix<T> subtract(const Matrix<T>& t1, const Matrix<T>& t2)
	{
		if (Matrix<T>::CheckIndicesAdd(t1,t2)) {
			Matrix<T> res(t1.GetNumRows(),t1.GetNumCols());
			for (int i=0; i<t1.GetNumRows(); i++)
					for (int j=0; j<t1.GetNumCols(); j++)
						res[i][j]=t1[i][j]-t2[i][j];
			return res;
		} else return t1;
	}

	template<class T>
	 static Matrix3D<T> subtract(const Matrix3D<T>& t1, const Matrix3D<T>& t2)
	{
		Matrix3D<T> res(t1.GetNumRows(),t1.GetNumCols(),t1.GetNumPols());
		for (int k=0; k<t1.GetNumPols(); k++) 
			for (int i=0; i<t1.GetNumRows(); i++)
					for (int j=0; j<t1.GetNumCols(); j++)
						res[k][i][j]=t1[k][i][j]-t2[k][i][j];
		return res;
	}

	
	// matrix multiplication functions
	template<class T>
	 static Vector<T> vsubtract(const Vector<T>& t1, const Vector<double>& t2)
	{
		if (t1.GetNum() == t2.GetNum()) {
			Vector<T> sum(t1.GetNum());
			for (int i=0; i<t1.GetNum(); i++)
				sum[i]=t1[i]-t2[i];
			return sum;;
		} else return Vector<T>();
	}


		// matrix multiplication functions
	template<class T>
	 static Vector<T> vadd(const Vector<T>& t1, const Vector<double>& t2)
	{
		if (t1.GetNum() == t2.GetNum()) {
			Vector<T> sum(t1.GetNum());
			for (int i=0; i<t1.GetNum(); i++)
				sum[i]=t1[i]+t2[i];
			return sum;;
		} else return Vector<T>();
	}

			// matrix multiplication functions
	template<class T>
	 static Vector<T> vmult(double d, const Vector<T>& t1)
	{
		Vector<T> sum(t1.GetNum());
		for (int i=0; i<t1.GetNum(); i++)
			sum[i]=d*t1[i];
		return sum;;
	}


	template<class T>
	 static Vector<T> vmultnew(T d, const Vector<double>& t1)
	{
		Vector<T> sum(t1.GetNum());
		for (int i=0; i<t1.GetNum(); i++)
			sum[i]=d*t1[i];
		return sum;;
	}

			// matrix multiplication functions
	template<class T>
	 static Matrix<T> mmult(double d, const Matrix<T>& t1)
	{
		Matrix<T> sum(t1.GetNumRows(),t1.GetNumCols());
		for (int i=0; i<t1.GetNumRows(); i++)
			for (int j=0; j<t1.GetNumCols(); j++)
				sum[i][j]=d*t1[i][j];

		return sum;;
	}

	// matrix multiplication functions
	template<class T>
	 static Matrix3D<T> mmult(double d, const Matrix3D<T>& t1)
	{
		Matrix3D<T> sum(t1.GetNumRows(),t1.GetNumCols(),t1.GetNumPols());
		for (int k=0; k<t1.GetNumPols(); k++) 
			for (int i=0; i<t1.GetNumRows(); i++)
				for (int j=0; j<t1.GetNumCols(); j++)
					sum[k][i][j]=d*t1[k][i][j];

		return sum;;
	}

	// matrix multiplication functions
	template<class T>
	static T mult0(const Vector<T>& t1, const Vector<double>& t2)
	{
		T sum=0.0;
		if (t1.GetNum() == t2.GetNum()) {
			for (int i=0; i<t1.GetNum(); i++)
				sum=sum+t1[i]*t2[i];
			return sum;;
		} else return sum;
	}


	// matrix multiplication functions


	 static Matrix<double> mult00(const Vector<double>& t1, const Vector<double>& t2)
	{
		Matrix<double> mat(t1.GetNum(),t2.GetNum());
		for (int i=0; i<t1.GetNum(); i++)
			for (int j=0; j<t2.GetNum(); j++)
				mat[i][j]=t1[i]*t2[j];
		return mat;;
	}

	
	 static Vector<double> ComputeParameterisation1(int num)
	{
		Vector<double> param(num);

		for (int i=0; i<num; i++) param[i] = (double)i/(double)num;

		return param;
	}

	
	template<class T>
	 static Vector<double> ComputeParameterisation2(const Vector<T>& pts, int num)
	{
		Vector<double> param(num);

		param[0]=0.0;
		double sum=0.0;
		for (int i=0; i<num-1; i++) {
			double temp = sqrt((pts[i+1]-pts[i])*(pts[i+1]-pts[i]));
			sum = sum + temp;
			param[i+1] = sum;
		}

		for (int i=0; i<num; i++) param[i] = param[i]/param[num-1];

		return param;
	}

	
	template<class T>
	 static Vector<double> ComputeParameterisation3(const Vector<T>& pts, int num)
	{
		Vector<double> param(num);

		param[0]=0.0;
		double sum=0.0;
		for (int i=0; i<num-1; i++) {
			double temp = sqrt((pts[i+1]-pts[i])*(pts[i+1]-pts[i]));
			sum = sum + sqrt(temp);
			param[i+1] = sum;
		}

		for (int i=0; i<num; i++) param[i] = param[i]/param[num-1];

		return param;
	}

	template<class T>
	 static Matrix<T> MatrixFromVector(Vector<T>& vec, int numu, int numv)
	{

		Matrix<T> mat(numu,numv);

		int count=0;
		for (int j=0; j<numv; j++)
			for (int i=0; i<numu; i++) {
				mat[i][j] = vec[count];
				count++;
			}

		return mat;
	}

	template<class T>
	 static Matrix3D<T> CreateMatrix3DTensor(const Vector<T>& vecx, const Vector<T>& vecy, const Vector<T>& vecz)
	{

		Matrix3D<T> mat(vecx.GetNum(),vecy.GetNum(),vecz.GetNum());

		for (int k=0; k<vecz.GetNum(); k++) 
			for (int j=0; j<vecy.GetNum(); j++)
				for (int i=0; i<vecx.GetNum(); i++) 
					mat[k][i][j] = vecx[i]*vecy[j]*vecz[k];

		return mat;
	}


	template<class T>
	 static Matrix3D<T> Matrix3DFromVector(Vector<T>& vec, int numu, int numv, int numw)
	{

		Matrix3D<T> mat(numu,numv,numw);

		int count=0;
		for (int k=0; k<numw; k++) 
			for (int j=0; j<numv; j++)
				for (int i=0; i<numu; i++) {
					mat[k][i][j] = vec[count];
					count++;
				}

		return mat;
	}

	template<class T>
	 static Vector<T> CreateKroneckerVector(const Matrix<T>& mat) 
	{
		int num = mat.GetNumRows()*mat.GetNumCols();

		Vector<T> v(num,0.0);

		int count=0;
		for (int j=0; j<mat.GetNumCols(); j++) 
			for (int i=0; i<mat.GetNumRows(); i++) {
				v[count]= mat[i][j];
				count++;
			}

		return v;
	}


	template<class T>
	 static Vector<T> CreateKroneckerVector(const Matrix3D<T>& mat) 
	{
		int num = mat.GetNumRows()*mat.GetNumCols()*mat.GetNumPols();

		Vector<T> v(num,0.0);

		int count=0;
		for (int k=0; k<mat.GetNumPols(); k++) 
			for (int j=0; j<mat.GetNumCols(); j++) 
				for (int i=0; i<mat.GetNumRows(); i++) {
					v[count]= mat[k][i][j];
					count++;
				}

		return v;
	}


	
	 static Matrix<double> kronecker(Matrix<double>& A, Matrix<double>& B)
	{
		Matrix<double> res(A.GetNumRows()*B.GetNumRows(),A.GetNumCols()*B.GetNumCols());

		int count1=0;
		int count2=0;
		int rowstart, colstart;

		for (int i=0; i<A.GetNumRows(); i++) {
			rowstart = i*B.GetNumRows();
			count1=rowstart;
			for (int j=0; j<A.GetNumCols(); j++) {
				colstart = j*B.GetNumCols();
				count2=colstart;
				for (int k=0; k<B.GetNumRows(); k++) {
					for (int l=0; l<B.GetNumCols(); l++)  {
						res[count1][count2] = A[i][j]*B[k][l];
						count2++;
					}
					count2=colstart;
					count1++;
				}
				count1=rowstart;
			}
		}
				

		return res;
	}

	template<class T>
	 static Matrix<T> Kronecker(Matrix<double>& A, Matrix<double>& B, Matrix<double>& C, Matrix<double>& D, Matrix<T>& E, Matrix<T>& X, Matrix<T>& Y)
	{
		
		Matrix<T> m1 = mult2(A,X);
		Matrix<T> m2 = mult1(m1,transpose(B));
	
		Matrix<T> n1 = mult2(C,X);
		Matrix<T> n2 = mult1(n1,transpose(D));
		
		Matrix<T> p1 = add(n2,m2);

		Matrix<T> q1 = subtract(p1,E);

		
		return q1;
	}



	// matrix multiplication functions
	template<class T>
	 static Matrix<T> mult1(const Matrix<T>& t1, const Matrix<double>& t2)
	{
		T sum=0.0;
		if (Matrix<T>::CheckIndicesMult(t1,t2)) {
			Matrix<T> res(t1.GetNumRows(),t2.GetNumCols());
			for (int i=0; i<t1.GetNumRows(); i++)
				for (int k=0; k<t2.GetNumCols(); k++) {
					for (int j=0; j<t1.GetNumCols(); j++)
						sum=sum+t1[i][j]*t2[j][k];
					res[i][k]=sum;
					sum=0.0;
				}
			return res;
		} else return t1;
	}

	template<class T>
	 static Matrix<T> mult2(const Matrix<double>& t1, const Matrix<T>& t2)
	{
		T sum=0.0;
		if (Matrix<T>::CheckIndicesMult(t1,t2)) {
			Matrix<T> res(t1.GetNumRows(),t2.GetNumCols());
			for (int i=0; i<t1.GetNumRows(); i++)
				for (int k=0; k<t2.GetNumCols(); k++) {
					for (int j=0; j<t1.GetNumCols(); j++)
						sum=sum+t1[i][j]*t2[j][k];
					res[i][k]=sum;
					sum=0.0;
				}
			return res;
		} else return t1;
	}


	
	 static Point3D mult31(const Matrix<double>& m, const Point4D& pt)
	{
		Vector<double> v(3);
		double sum=0.0;

		for (int i=0; i<3; i++) {
			sum+=m[i][0]*pt.GetX()+m[i][1]*pt.GetY()+m[i][2]*pt.GetZ()+m[i][3]*pt.GetW();			
			v[i]=sum;
			sum=0.0;
		}
		return Point3D(v[0],v[1],v[2]);
	}

	template<class T>
	 static Vector<T> mult3(const Matrix<double>& t1, const Vector<T>& t2)
	{
		T sum=0.0;
		if (t1.GetNumCols() == t2.GetNum()) {
			Vector<T> res(t1.GetNumRows());
			for (int i=0; i<t1.GetNumRows(); i++) {
				for (int k=0; k<t1.GetNumCols(); k++) 
					sum=sum+t1[i][k]*t2[k];
				res[i]=sum;
				sum=0.0;
			}
			return res;
		} else return t2;
	}

	template<class T>
	 static Vector<T> mult4(const Vector<T>& t1, const Matrix<double>& t2)
	{	
		T sum=0.0;
		if (t1.GetNum() == t2.GetNumRows()) {
			Vector<T> res(t2.GetNumCols());
			for (int i=0; i<t2.GetNumCols(); i++) {
				for (int k=0; k<t2.GetNumRows(); k++) 
					sum=sum+t1[k]*t2[k][i];
				res[i]=sum;
				sum=0.0;
			}
			return res;
		} else return t1;
	}

	template<class T>
	 static Matrix3D<T> mult5(const Matrix3D<T>& t1, const Matrix<double>& t2)
	{
		T sum=0.0;
			Matrix3D<T> res(t1.GetNumRows(),t2.GetNumCols(),t1.GetNumPols());
			for (int l=0; l<t1.GetNumPols(); l++) {
				for (int i=0; i<t1.GetNumRows(); i++)
					for (int k=0; k<t2.GetNumCols(); k++) {
						for (int j=0; j<t1.GetNumCols(); j++)
							sum=sum+t1[l][i][j]*t2[j][k];
						res[l][i][k]=sum;
						sum=0.0;
					}
			}
			return res;
	
	} 	


	template<class T>
	 static Matrix3D<T> mult10(const Matrix3D<T>& t1, const Matrix<double>& t2)
	{
		T sum=0.0;
			Matrix3D<T> res(t1.GetNumRows(),t1.GetNumCols(),t2.GetNumCols());
			for (int i=0; i<t1.GetNumRows(); i++) {
				for (int j=0; j<t1.GetNumCols(); j++)
					for (int k=0; k<t2.GetNumCols(); k++) {
						for (int l=0; l<t1.GetNumPols(); l++)
							sum=sum+t1[l][i][j]*t2[l][k];
						res[k][i][j]=sum;
						sum=0.0;
					}
			}
			return res;
	
	} 	
	template<class T>
	 static Matrix3D<T> mult6(const Matrix<double>& t1, const Matrix3D<T>& t2)
	{
		T sum=0.0;
		if (Matrix3D<T>::CheckIndicesMult(t1,t2)) {
			Matrix3D<T> res(t1.GetNumRows(),t2.GetNumCols(),t2.GetNumPols());
			for (int l=0; l<t2.GetNumPols(); l++) {
				for (int i=0; i<t1.GetNumRows(); i++)
					for (int k=0; k<t2.GetNumCols(); k++) {
						for (int j=0; j<t1.GetNumCols(); j++)
							sum=sum+t1[i][j]*t2[l][j][k];
						res[l][i][k]=sum;
						sum=0.0;
					}
			}
			return res;
		} else return t2;
	} 	


	template<class T>
	 static Matrix<T> mult7(const Vector<double>& t1, const Matrix3D<T>& t2)
	{
		T sum=0.0;
			Matrix<T> res(t2.GetNumPols(),t2.GetNumCols());
			for (int l=0; l<t2.GetNumPols(); l++) {
					for (int k=0; k<t2.GetNumCols(); k++) {
						for (int j=0; j<t1.GetNum(); j++)
							sum=sum+t1[j]*t2[l][j][k];
						res[l][k]=sum;
						sum=0.0;
					}
			}
			return res;
	} 	


	template<class T>
	 static Matrix<T> mult8(const Matrix3D<T>& t1, const Vector<double>& t2)
	{
		T sum=0.0;
			Matrix<T> res(t1.GetNumPols(),t2.GetNum());
			for (int l=0; l<t1.GetNumPols(); l++) {
				for (int i=0; i<t1.GetNumRows(); i++) {
						for (int j=0; j<t1.GetNumCols(); j++)
							sum=sum+t1[l][i][j]*t2[j];
						res[l][i]=sum;
						sum=0.0;
					}
			}
			return res;
	} 	



	template<class T>
	 static Matrix3D<T> mult9(const Matrix<double>& t1, const Matrix3D<T>& t2)
	{
		T sum=0.0;
			Matrix3D<T> res(t2.GetNumRows(),t2.GetNumCols(),t1.GetNumRows());
			for (int l=0; l<t2.GetNumRows(); l++) {
				for (int i=0; i<t2.GetNumCols(); i++)
					for (int k=0; k<t1.GetNumRows(); k++) {
						for (int j=0; j<t1.GetNumCols(); j++)
							sum=sum+t1[k][j]*t2[j][l][i];
						res[k][l][i]=sum;
						sum=0.0;
					}
			}
			return res;
	} 	

		template<class T>
	 static Matrix<T> VM3DI_1(const Vector<double>& t1, const Matrix3D<T>& t2)
	{
		T sum=0.0;
			Matrix<T> res(t2.GetNumCols(),t2.GetNumPols());
			for (int k=0; k<t2.GetNumPols(); k++) {
				for (int j=0; j<t2.GetNumCols(); j++) {
					for (int i=0; i<t2.GetNumRows(); i++)
						sum=sum+t2[k][i][j]*t1[i];
					res[j][k]=sum;
					sum=0.0;
					}
			}
			return res;
	} 	

	template<class T>
	 static Matrix<T> VM3DI_2(const Vector<double>& t1, const Matrix3D<T>& t2)
	{
		T sum=0.0;
			Matrix<T> res(t2.GetNumCols(),t2.GetNumPols());
			for (int j=0; j<t2.GetNumCols(); j++) {
				for (int k=0; k<t2.GetNumPols(); k++) {
					for (int i=0; i<t2.GetNumRows(); i++)
						sum=sum+t2[k][i][j]*t1[i];
					res[j][k]=sum;
					sum=0.0;
					}
			}
			return res;
	} 	


		template<class T>
	 static Matrix<T> VM3DJ_1(const Vector<double>& t1, const Matrix3D<T>& t2)
	{
		T sum=0.0;
			Matrix<T> res(t2.GetNumRows(),t2.GetNumPols());
			for (int k=0; k<t2.GetNumPols(); k++) {
				for (int i=0; i<t2.GetNumRows(); i++) {
					for (int j=0; j<t2.GetNumCols(); j++)
						sum=sum+t2[k][i][j]*t1[j];
					res[i][k]=sum;
					sum=0.0;
					}
			}
			return res;
	} 	

	template<class T>
	 static Matrix<T> VM3DJ_2(const Vector<double>& t1, const Matrix3D<T>& t2)
	{
		T sum=0.0;
			Matrix<T> res(t2.GetNumRows(),t2.GetNumPols());
			for (int i=0; i<t2.GetNumRows(); i++) {
				for (int k=0; k<t2.GetNumPols(); k++) {
					for (int j=0; j<t2.GetNumCols(); j++)
						sum=sum+t2[k][i][j]*t1[j];
					res[i][k]=sum;
					sum=0.0;
					}
			}
			return res;
	} 	



		template<class T>
	 static Matrix<T> VM3DK_1(const Vector<double>& t1, const Matrix3D<T>& t2)
	{
		T sum=0.0;
			Matrix<T> res(t2.GetNumRows(),t2.GetNumCols());
			for (int i=0; i<t2.GetNumRows(); i++) {
				for (int j=0; j<t2.GetNumCols(); j++) {
					for (int k=0; k<t2.GetNumPols(); k++)
						sum=sum+t2[k][i][j]*t1[k];
					res[i][j]=sum;
					sum=0.0;
					}
			}
			return res;
	} 	

	template<class T>
	 static Matrix<T> VM3DK_2(const Vector<double>& t1, const Matrix3D<T>& t2)
	{
		T sum=0.0;
			Matrix<T> res(t2.GetNumRows(),t2.GetNumPols());
			for (int j=0; j<t2.GetNumCols(); j++) {
				for (int i=0; i<t2.GetNumRows(); i++) {
					for (int k=0; k<t2.GetNumPols(); k++)
						sum=sum+t2[k][i][j]*t1[k];
					res[i][j]=sum;
					sum=0.0;
					}
			}
			return res;

	} 	




	template<class T>
	 static Matrix3D<T> MM3DI_1(const Matrix<double>& t1, const Matrix3D<T>& t2)
	{
		T sum=0.0;
			Matrix3D<T> res(t1.GetNumCols(),t2.GetNumCols(),t2.GetNumPols());
			for (int k=0; k<t2.GetNumPols(); k++) {
				for (int l=0; l<t1.GetNumCols(); l++)
					for (int j=0; j<t2.GetNumCols(); j++) {
						for (int i=0; i<t2.GetNumRows(); i++)
							sum=sum+t1[i][l]*t2[k][i][j];
						res[k][l][j]=sum;
						sum=0.0;
					}
			}
			return res;
	} 	

		template<class T>
	 static Matrix3D<T> MM3DI_2(const Matrix<double>& t1, const Matrix3D<T>& t2)
	{
		T sum=0.0;
			Matrix3D<T> res(t1.GetNumCols(),t2.GetNumCols(),t2.GetNumPols());
			for (int j=0; j<t2.GetNumCols(); j++) {
				for (int l=0; l<t1.GetNumCols(); l++)
					for (int k=0; k<t2.GetNumPols(); k++) {
						for (int i=0; i<t2.GetNumRows(); i++)
							sum=sum+t1[i][l]*t2[k][i][j];
						res[k][l][j]=sum;
						sum=0.0;
					}
			}
			return res;

	} 	

	
	template<class T>
	 static Matrix3D<T> MM3DJ_1(const Matrix<double>& t1, const Matrix3D<T>& t2)
	{
		T sum=0.0;

			Matrix3D<T> res(t2.GetNumRows(),t1.GetNumCols(),t2.GetNumPols());
			for (int k=0; k<t2.GetNumPols(); k++) {
				for (int l=0; l<t1.GetNumCols(); l++)
					for (int i=0; i<t2.GetNumRows(); i++) {
						for (int j=0; j<t2.GetNumCols(); j++)
							sum=sum+t1[j][l]*t2[k][i][j];
						res[k][i][l]=sum;
						sum=0.0;
					}
			}
			return res;

	} 	

		template<class T>
	 static Matrix3D<T> MM3DJ_2(const Matrix<double>& t1, const Matrix3D<T>& t2)
	{
		T sum=0.0;

			Matrix3D<T> res(t2.GetNumRows(),t1.GetNumCols(),t2.GetNumPols());
			for (int i=0; i<t2.GetNumRows(); i++) {
				for (int l=0; l<t1.GetNumCols(); l++)
					for (int k=0; k<t2.GetNumPols(); k++) {
						for (int j=0; j<t2.GetNumCols(); j++)
							sum=sum+t1[j][l]*t2[k][i][j];
						res[k][i][l]=sum;
						sum=0.0;
					}
			}
			return res;

	} 	

	
	template<class T>
	 static Matrix3D<T> MM3DK_1(const Matrix<double>& t1, const Matrix3D<T>& t2)
	{
		T sum=0.0;

			Matrix3D<T> res(t2.GetNumRows(),t2.GetNumCols(),t1.GetNumCols());
			for (int j=0; j<t2.GetNumCols(); j++) {
				for (int l=0; l<t1.GetNumCols(); l++)
					for (int i=0; i<t2.GetNumRows(); i++) {
						for (int k=0; k<t2.GetNumPols(); k++)
							sum=sum+t1[k][l]*t2[k][i][j];
						res[l][i][j]=sum;
						sum=0.0;
					}
			}
			return res;
/
	} 	

		template<class T>
	 static Matrix3D<T> MM3DK_2(const Matrix<double>& t1, const Matrix3D<T>& t2)
	{
		T sum=0.0;

			Matrix3D<T> res(t2.GetNumRows(),t2.GetNumCols(),t1.GetNumCols());
			for (int i=0; i<t2.GetNumRows(); i++) {
				for (int l=0; l<t1.GetNumCols(); l++)
					for (int j=0; j<t2.GetNumCols(); j++) {
						for (int k=0; k<t2.GetNumPols(); k++)
							sum=sum+t1[k][l]*t2[k][i][j];
						res[l][i][j]=sum;
						sum=0.0;
					}
			}
			return res;

	} 	


	template<class T>
	 static Matrix<T> swapmat(const Matrix<T>& mat, int ind1, int ind2)
	{
		Matrix<T> m(mat);

		Vector<T> v1 = m.GetRow(ind1);
		Vector<T> v2 = m.GetRow(ind2);

		for (int j=0; j<m.GetNumCols(); j++) m[ind1][j] = v2[j];
		for (int j=0; j<m.GetNumCols(); j++) m[ind2][j] = v1[j];

		return m;
	}

	
	 static Matrix<double> swapmat1(const Matrix<double>& mat, int ind1, int ind2)
	{
		Matrix<double> m(mat);

		Vector<double> v1 = m.GetRow(ind1);
		Vector<double> v2 = m.GetRow(ind2);

		for (int j=0; j<m.GetNumCols(); j++) m[ind1][j] = v2[j];
		for (int j=0; j<m.GetNumCols(); j++) m[ind2][j] = v1[j];

		return m;
	}


	template<class T>
	 static Vector<T> swapvec(const Vector<T>& vec, int ind1, int ind2)
	{
		Vector<T> v(vec);

		T t1 = v[ind1];
		T t2 = v[ind2];

		v[ind1] = t2;
		v[ind2] = t1;

		return v;
	}


	template<class T>
	 static Matrix<T> Transform(const Matrix<T>& cpts, const Matrix<double>& m)
	{
		Matrix<T> ncpts(cpts.GetNumRows(),cpts.GetNumCols());

		for (int i=0; i<cpts.GetNumRows(); i++)
			for (int j=0; j<cpts.GetNumCols(); j++) ncpts[i][j] = mult31(m,Point4D(cpts[i][j]));

		return ncpts;
	}

	
	

	template<class T>
	static Matrix<T> transpose(const Matrix<T>& t1)
	{
		Matrix<T> mat(t1.GetNumCols(),t1.GetNumRows());

		for (int i=0; i<t1.GetNumRows(); i++)
			for (int j=0; j<t1.GetNumCols(); j++)
				mat[j][i]=t1[i][j];
		return mat;
	} 	


	template<class T>
	 static FMatrix<T> transpose(const FMatrix<T>& t1)
	{
		FMatrix<T> mat(t1.GetMinCol(),t1.GetMaxCol(),t1.GetMinRow(),t1.GetMaxRow());

		for (int i=t1.GetMinRow(); i<=t1.GetMaxRow(); i++)
			for (int j=t1.GetMinCol(); j<=t1.GetMaxCol(); j++)
				mat[j][i]=t1[i][j];
		return mat;
	} 	

	
	 static Matrix<double> EliminateMat(const Matrix<double>& mat, const Vector<double>& eqn, int n, int& num)
	{
		// find first non-zero component from eqn vector and then eliminate the
		// equation corresponding to this element found in the matrix equation

		// find index
		int i=-1;
		Matrix<double> newmat(mat);
		int index;
		double Tol=0.000001;
		do {
			i++;
		} while (eqn[i]<Tol && i < n);
		index = i;



		num = index;
		
		for (int j=0; j<n; j++) newmat[index][j] = eqn[j];
		return newmat;
	}

	template<class T>
	 static Vector<double> EliminateVec(const Matrix<double>& mat, const Vector<double>& vec, const Vector<double>& eqn, T rhs, int n)
	{
		int i=-1;
		Vector<double> newvec(vec);
		int index;
		double Tol=0.000001;
		do {
			i++;
		} while (eqn[i]<Tol && i < n);
		index = i;

		newvec[index] = rhs;
		return newvec;
	}


	
	 static Vector<int> GetEliminateIndices(const Matrix<double>& mat1, const Matrix<double>& mat2)
	{
		// 1. take a row of the matrix mat2
		// 2. find the first non-zero element
		// 3. find vector expressing that index element in terms of all others (that index missing)
		// 4. for each row of mat1 excluding index found in 2 replace that index with vector from 3
		// 5. collect together correspoding elements to get new row for matrix


		double tol = 0.0000001;

		Vector<int> v1(mat2.GetNumRows(),0);
		Matrix<double> mat3(mat2);

		bool flag1 = false;
		int swap1, swap2;
		do {
			flag1 = false;
			for (int i=0; i<mat2.GetNumRows(); i++) {
				Vector<double> v = mat2.GetRow(i);
				bool flag = false;
				int j=-1;
				do {
					flag = false;
					do {
						j++;
					} while (j < mat2.GetNumCols()-1 && fabs(v[j]) < tol);
					v1[i]=j;
					if (j == mat2.GetNumCols()-2 && fabs(v[j]) < tol) {
						// swap rows k and i
						mat3 = swapmat1(mat3,swap1,swap2);
						flag1 = true;
						i = mat2.GetNumRows(); // exit loop
					} else {
						for (int k=0; k<=i-1; k++) 
							if (v1[i] == v1[k]) {
								flag = true;
								swap1 = k;
								swap2 = i;
							}
					}
				} while (flag);
			}
		} while (flag1);

	
		return v1;
	}


	
	 static Matrix<double> GetReducedMatrix(Matrix<double>& mat2)
	{
		return mat2;
	}

	



	 static Matrix<double> EliminateMVariables(Matrix<double>& mat1, Matrix<double>& mat5, Vector<int>& v1)
	{
		// 1. take a row of the matrix mat2
		// 2. find the first non-zero element
		// 3. find vector expressing that index element in terms of all others (that index missing)
		// 4. for each row of mat1 excluding index found in 2 replace that index with vector from 3
		// 5. collect together correspoding elements to get new row for matrix

		Matrix<double> mat2(mat5);

		// eliminate matrix elements
		Matrix<double> mat3(mat1);
		for (int k=0; k<mat2.GetNumRows(); k++) {
			for (int i=0; i<mat1.GetNumRows(); i++) {
				double temp = mat3[i][v1[k]];
				for (int j=0; j<mat1.GetNumCols(); j++) {
					mat3[i][j] = mat3[i][j]-temp*mat2[k][j]/mat2[k][v1[k]];
				}
			}
		}
		Matrix<double> mat4 = EliminateRowCol(mat3,mat3.GetNumRows(),v1,mat2.GetNumRows());
		return mat4;
	}


	template<class T>
	 static Vector<T> EliminateVVariables(Matrix<double>& mat1, Matrix<double>& mat4, Vector<T>& rhs1, Vector<int>& v1)
	{
		// 1. take a row of the matrix mat2
		// 2. find the first non-zero element
		// 3. find vector expressing that index element in terms of all others (that index missing)
		// 4. for each row of mat1 excluding index found in 2 replace that index with vector from 3
		// 5. collect together correspoding elements to get new row for matrix

		double tol = 0.0000001;

		Matrix<double> mat2(mat4);
		Matrix<double> mat3(mat1);
		Vector<double> rhs(rhs1);
	
		
		Vector<T> v(rhs);

		
	
		for (int k=0; k<mat2.GetNumRows(); k++) 
			for (int i=0; i<mat1.GetNumRows(); i++) {
				v[i] = v[i]-mat2[k][mat2.GetNumCols()-1]*mat1[i][v1[k]];///mat2[k][v1[k]];
			
			}

		Vector<T> v2 = EliminateRow(v,mat1.GetNumRows(), v1,mat2.GetNumRows());
		
		
		return v2;
	}



	template<class T>
	 static Vector<T> EliminateVVariables(Matrix<double>& mat1, Matrix<double>& mat4, Vector<T>& rhs2, Vector<T>& rhs1, Vector<int>& v1)
	{
		// 1. take a row of the matrix mat2
		// 2. find the first non-zero element
		// 3. find vector expressing that index element in terms of all others (that index missing)
		// 4. for each row of mat1 excluding index found in 2 replace that index with vector from 3
		// 5. collect together correspoding elements to get new row for matrix

		double tol = 0.0000001;

		Matrix<double> mat2(mat4);
		Matrix<double> mat3(mat1);
		Vector<T> rhs(rhs1);
	
		
		
	
		// rhs vector
		Vector<T> v(rhs);

		
	
		for (int k=0; k<mat2.GetNumRows(); k++) 
			for (int i=0; i<mat1.GetNumRows(); i++) {
				v[i] = v[i]-rhs2[k]*mat1[i][v1[k]];
		
			}

		Vector<T> v2 = EliminateRow(v,mat1.GetNumRows(), v1,mat2.GetNumRows());
		
		
		return v2;
	}


	template<class T>
	 static Vector<T> EliminateRow(const Vector<T>& v, int n, Vector<int>& list, int m)
	{
		std::set<int> set1(list.begin(),list.end());
		
		Vector<int> v1(n);
		
		for (int i=0; i<n; i++) v1[i] = i;

		std::set<int> set2(v1.begin(),v1.end());
		std::set<int> set3;

		std::set_difference(set2.begin(),set2.end(),set1.begin(),set1.end(),std::inserter(set3,set3.begin()));
		Vector<int> diff(set3.size());
		std::copy(set3.begin(),set3.end(),diff.begin());

		Vector<T> res(set3.size());

		for (unsigned int i=0; i<set3.size(); i++) res[i]=v[diff[i]];

		return res;

	}

	
	 static Matrix<double> EliminateRowCol(const Matrix<double>& mat, int n, Vector<int>& list, int m)
	{
		std::set<int> set1(list.begin(),list.end());
		
		Vector<int> v(n);
		
		for (int i=0; i<n; i++) v[i] = i;

		std::set<int> set2(v.begin(),v.end());
		std::set<int> set3;

		std::set_difference(set2.begin(),set2.end(),set1.begin(),set1.end(),std::inserter(set3,set3.begin()));
		Vector<int> diff(set3.size());
		std::copy(set3.begin(),set3.end(),diff.begin());

		Matrix<double> res(set3.size(),set3.size());

		for (unsigned int i=0; i<set3.size(); i++) 
			for (unsigned int j=0; j<set3.size(); j++) 
				res[i][j]=mat[diff[i]][diff[j]];

		return res;
	}
	
	template<class T>
	 static Vector<T> Solve(const Matrix<double>& mat, const Vector<T>& rhs)
	{
		T sum;

		Matrix<double> LU(mat);
		Vector<T> x(rhs);
    
		int n = rhs.GetNum();
		// Call factorization program. 

	
		LU=FactLU(LU);

		// Forward substitute. 

		for(int j=0;j<n;j++)
		{
			sum=x[j];
			for(int i=0;i<=j-1;i++)
				sum=sum-LU[j][i]*x[i];
			x[j]=sum/LU[j][j];
		}
      
		/* Back substitute. */

		for(int j=n-1;j>=0;j--)
		{
			sum=x[j];
			for(int i=j+1;i<n;i++)
				sum=sum-LU[j][i]*x[i];
			x[j]=sum;
		}
		return(x);
	}                                                 
     
    template<class T>
	 static Vector<T> Tridiag(const Vector<double>& a, const Vector<T>& b, const Vector<double>& c, const Vector<double>& d, int n)
	{
		Vector<double> d1(d);
		Vector<T> b1(b);

		double z;

		for (int i=1; i<n; i++) {
			z = a[i-1]/d[i-1];
			d1[i] = d1[i] - z * c[i-1];
			b1[i] = b1[i] - z * b1[i-1];
		}

		b1[n-1] = b1[n-1]/d1[n-1];

		for (int i=0; i<n-1; i++) {
			b1[n-i-1] = (b1[n-i-1] - c[n-i-1]*b1[n-i-1])/d1[n-2-i];
		}

		return b1;
	}                

	
	 static Matrix<double> FactLU(const Matrix<double>& mat)
	{
		double piv,sum;
		int n=mat.GetNumRows();
		Matrix<double> LU(mat);
	
    
		/* Column 1 of L is Column 1 of A. */
		/* Find row 1 of U. */
      
		/* Now work through the remainder. */

	    for(int k=0;k<n;k++)
		{
      
			/* Find column k of L. */

			for(int j=k;j<n;j++)
			{
				sum=LU[j][k];            
				for(int i=0; i<=k-1;i++)
					sum-=LU[j][i]*LU[i][k];
				LU[j][k]=sum;
			}
         
			/* Find row k of U except for k=n (diagonal elements are assumed). */

			piv=LU[k][k];
	
			for(int j=k+1;j<n;j++)
			{
				sum=LU[k][j];
				for (int i=0;i<=k-1;i++)
					sum-=LU[k][i]*LU[i][j];
				LU[k][j]=sum/piv;
			}
		}
		return(LU);  
	}

	
	 static Vector<double> CreateKnots(int Num, int Ord, const Vector<double>& Limits)
	{
		// assume all segments of same order	
		int count=0;

		Vector<double> knots((Num+1)*Ord);

		// assign knots using increasing integers
		for (int j=0; j<=Num; j++) 
			for (int i=0; i<Ord; i++) {
				knots[count] = Limits[j];
				count++;
			}
		return knots;
	}

	
	 static Vector<double> CreateKnots(int Ord) 
	{
		// assign knots 0 and 1 of multiplicity ord
		Vector<double> Kts(2*Ord);
	
		for (int i=0; i<Ord; i++) {
			Kts[i]=0.0;
			Kts[i+Ord]=1.0;
		}
		return Kts;
	}	
	
	
	 static Vector<double> CreateKnots(int Ord, int Num) 
	{
		// assume all segments of same order	
		int count = 0;

		Vector<double> Kts(Num+Ord);
		// assign knots using increasing integers
		for (int i=0; i<Ord; i++) {
			Kts[count] = 0.0;
			count++;
		}
	
		for (int i=0; i<Num-Ord+1; i++) {
			Kts[count] = (double)(i+1);
			count++;
		}

		for (int i=1; i<Ord; i++) {
			Kts[count] = (double)(Num-Ord+1);
			count++;
		}
		return Kts;
	}

	template<class T>
	 static T** CStyleMatrix(const Matrix<T>& m) 
	{
		T **mat;
		int nrl = 0;
		int nrh = m.GetNumRows()-1;
		int ncl = 0;
		int nch = m.GetNumCols()-1;

		mat=(T **) malloc((unsigned) (nrh-nrl+1)*sizeof(T*));
		if (!mat) std::cerr << "allocation failure 1 in dmatrix()";
	
		for(int i=0;i<=nrh;i++) {
			mat[i]=(T *) malloc((unsigned) (nch-ncl+1)*sizeof(T));
			if (!mat[i]) std::cerr << "allocation failure 2 in dmatrix()";
		}

		for (int i=0; i<=nrh; i++)
			for (int j=0; j<=nch; j++) mat[i][j] = m[i][j];
		return mat;
	}

	template<class T>
	 static T* CStyleVector(const Vector<T>& v) 
	{
		T *v1;

		v1=(T *)malloc((unsigned) 200*sizeof(T));
		if (!v1) std::cerr << "allocation failure in dvector()";

		for (int i=0; i<v.GetNum(); i++) v1[i] = v[i];
		return v1;
	}

	template<class T>
	 static Vector<T> gauss(Matrix<double>& a, Vector<T>& b, int n )
	{
		int i, j, k;
		T sum;
		double c;
		Vector<int> icol(n);
		Vector<T> x(n);
		Vector<T> b1(b);
		Matrix<double> a1(a);
	

		// Initialize icol 

		for ( k=0; k<n; k++ ) icol[k]=k;

		// Outer loop - eliminate column

		for (int k=0; k<n-1; k++ ) {

		// Find the best pivotal element

			pivot( a1, b1, n, k, icol );
		
			

		// This loop - row to operate on

			for ( i=k+1; i<n; i++ ) {

		// This loop actually does it

				c = a1[i][k] / a1[k][k];
				for ( j=k; j<n; j++ ) 
					a1[i][j] -= a1[k][j] * c;

		// Operate also on b.                                                

				b1[i] = b1[i] - b1[k] * c;
			}
		}
	
		// Back substitute.                                                      

		x[n-1] = b1[n-1] / a1[n-1][n-1];

		// Take rows in reverse order.                                           

		for ( k=n-1; k>=0; k-- ) {

		// Sum over known results.                                              

		    sum = b1[k];
			for ( j=k+1; j<n; j++ )
				sum = sum- x[j] * a1[k][j];

			x[k] = sum / a1[k][k];
		}

		// Finally undo reordering caused by chosing best pivot.            

		reorder( x, b1, n, icol );

		return x;
	}

/****************************************************************************/
/* Function to pivot equations.                                             */
/*                                                                          */
/* Search n x n real array a for largest absolute value a[i][j]             */
/* for i,j  >=  k.  Then switch rows and columns to bring this              */
/* to a[k][k].  Also switch right hand side vector b and column             */
/* memory array icol.                                                       */
/****************************************************************************/

	template<class T>
		static void pivot(Matrix<double>& a, Vector<T>& b, int n, int k, Vector<int>& icol )
	{
		int imax, jmax, i, j;
		double absa, absmx, z;

		if(k >= n) return;


		// Search a for largest a[i][j]                                           

		absmx = 0.0;
		imax = k;
		jmax = k;

		for ( j=k; j<n; j++ )	{
			for ( i=k; i<n; i++ ) {
				if((absa = fabs( a[i][j] )) > absmx) {
					absmx = absa;
					imax = i;
					jmax = j;
				}
			}
		}

		// Switch rows k and imax in a,b                                          

		for ( j=0; j<n; j++ ) {
			swap( a[k][j], a[imax][j] );
		}
		swap( b[k], b[imax] );

		// Switch columns k,jmax in a,icol.                                       

		for ( i=0; i<n; i++ ) {
			swap( a[i][k], a[i][jmax] );
		}
		swap( icol[k], icol[jmax] );

		// Scale row k.                                                           

		z = a[k][k];
		for ( j=k; j<n; j++ )	{
			a[k][j] /= z;
	    }

		b[k] = b[k]/ z;
	}

/****************************************************************************/
/* Undo column rearrangements.                                              */
/*                                                                          */
/* On input x is real solution vector - length n - which has been reordered. */
/* int  array icol on input holds desired order. b is used as working;      */
/* space.  On output x is the solution in the correct order.                */
/****************************************************************************/

	template<class T>
		static void reorder(Vector<T>& x, Vector<T>& b, int n, Vector<int>& icol )
	{
		int j, itrue;

		for ( j=0; j<n; j++ )	{
			itrue = icol[j];
			b[itrue] = x[j];
		}

		for ( j=0; j<n; j++ ) {
			x[j] = b[j];
		}
	}

};




#endif