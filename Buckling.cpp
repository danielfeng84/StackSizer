// This file is part of StackSizer, a kit for stack design of successive laminates.
//
// Copyright (C) 2018 Jianwen Feng <jianwenfeng@hotmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.



#include "Buckling.h"
#include "Composite.h"
#include <math.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;





struct PlateBuckleState PlateBuckling_UC_C_SSSS_Ortho(struct Laminate plate)
{
	struct PlateBuckleState BucklingLoad ;
	double m, n, t, epsi, candidate;
	struct CompsiteE E;
	double temp1, temp2;

	epsi = 1.0e19;
	BucklingLoad = PlateBuckleState();
	t = PlateThick(plate);
	E = PlateStiff(plate, 0);

	for (m=1.;m<=6.;m=m+1.)
	{
		for (n=1.; n<=5; n=n+1.)
		{
			temp1 = (m*PI/plate.a)*(m*PI/plate.a);
			temp2 = (n*plate.a/m/plate.b)*(n*plate.a/m/plate.b);
			candidate = temp1*(plate.stack.abdm.D11 + 2.*temp2*(plate.stack.abdm.D12+2.*plate.stack.abdm.D66) + temp2*temp2*plate.stack.abdm.D22 )/E.Ex/t + E.Ey/temp1/plate.r/plate.r/(E.Ex - temp2*(2.*E.MIUxy*E.Ey-E.Ex*E.Ey/E.Gxy) + E.Ey*temp2*temp2);
			if (candidate<epsi)
			{
				epsi = candidate;
				BucklingLoad.m=(int)m;
				BucklingLoad.n=(int)n;
			}
		}
	}


	BucklingLoad.CrBuckleStrain.InPlane.EpsiX = -epsi;
	BucklingLoad.CrBuckleLoad.InPlane.Nxx = BucklingLoad.CrBuckleStrain.InPlane.EpsiX * E.Ex * t;
	return(BucklingLoad);
}

struct PlateBuckleState PlateBuckling_UC_F_SSSS_Ortho(struct Laminate plate)
{

	struct PlateBuckleState BucklingLoad ;
	struct PlateBuckleState Candidate ;
	double AR, theta,m;


    BucklingLoad = PlateBuckleState();
    Candidate=PlateBuckleState();
	BucklingLoad.CrBuckleLoad.InPlane.Nxx=1.0E20;

	AR = plate.a/plate.b;
	for (m=1.;m<=8.;m=m+1.)
	{
		theta = m*m*m*m*plate.stack.abdm.D11 + 2.0*m*m*AR*AR*(plate.stack.abdm.D12 + 2.0*plate.stack.abdm.D66) + AR*AR*AR*AR*plate.stack.abdm.D22;

		Candidate.CrBuckleLoad.InPlane.Nxx = PI*PI*theta/plate.a/plate.a/m/m;
		Candidate.m = (int)m;
		Candidate.n = 1;

		if(Candidate.CrBuckleLoad.InPlane.Nxx<BucklingLoad.CrBuckleLoad.InPlane.Nxx)


				BucklingLoad = Candidate;
	}


	BucklingLoad.CrBuckleLoad.InPlane.Nxx = - BucklingLoad.CrBuckleLoad.InPlane.Nxx;// change into negative value
	BucklingLoad.CrBuckleStrain = Disp(plate, BucklingLoad.CrBuckleLoad);

	return(BucklingLoad);

}




struct PlateBuckleState PlateBuckling_UC_F_CCCC_Ortho(struct Laminate plate)
{
         struct PlateBuckleState BucklingLoad ;
	     struct PlateBuckleState Candidate ;
         double kx,AR,C1,C2, lamda;
         double m;

         BucklingLoad = PlateBuckleState();
         Candidate = PlateBuckleState();

         AR = plate.a/plate.b;

         lamda = plate.a/plate.b*sqrt(sqrt(plate.stack.abdm.D22/plate.stack.abdm.D11));

         C1=plate.stack.abdm.D11/plate.stack.abdm.D22;
         C2=(plate.stack.abdm.D12+2*plate.stack.abdm.D66)/plate.stack.abdm.D22;
         BucklingLoad.CrBuckleLoad.InPlane.Nxx=1.0E20;

         if (lamda<1.11)
         {
                   kx=4*PI*PI/AR/AR*C1+8/3*C2+4*AR*AR;
                   BucklingLoad.CrBuckleLoad.InPlane.Nxx = kx/plate.b/plate.b*(plate.stack.abdm.D22);
                   BucklingLoad.m=1;

         }
         else
         {
                   for (m=2.;m<=20.;m=m+1.)
                  {
                            kx=PI*PI*(((m-1)*(m-1)*(m-1)*(m-1)+(m+1)*(m+1)*(m+1)*(m+1))/((m-1)*(m-1)+(m+1)*(m+1))/AR/AR*C1+8/3*C2+10.66*AR*AR/((m-1)*(m-1)+(m+1)*(m+1)));
                            Candidate.CrBuckleLoad.InPlane.Nxx = kx/plate.b/plate.b*(plate.stack.abdm.D22);
                            Candidate.m = (int)m;
                            Candidate.n = 1;

                            if(Candidate.CrBuckleLoad.InPlane.Nxx < BucklingLoad.CrBuckleLoad.InPlane.Nxx)
                                     BucklingLoad = Candidate;
                   }
         }




	BucklingLoad.CrBuckleLoad.InPlane.Nxx = - BucklingLoad.CrBuckleLoad.InPlane.Nxx;// change into negative value
	BucklingLoad.CrBuckleStrain = Disp(plate, BucklingLoad.CrBuckleLoad);



         return(BucklingLoad);
}

struct PlateBuckleState PlateBuckling_UC_F_SSSF_Ortho(struct Laminate plate)  // contains bug
{

	struct PlateBuckleState BucklingLoad ;
	struct PlateBuckleState Candidate ;
	double m, alpha, beta, a, b, D11, D12, D22, D66, N, can, target, temp1, temp2, k0, pace, count;

	BucklingLoad = PlateBuckleState();
    Candidate = PlateBuckleState();

	a =  plate.a;
	b =  plate.b;
	D11 = plate.stack.abdm.D11;
	D12 = plate.stack.abdm.D12;
	D22 = plate.stack.abdm.D22;
	D66 = plate.stack.abdm.D66;

	can = 1.e20;

	for (m=1.;m<10.;m=m+1.)
	{
		N = D11*m*m*PI*PI/a/a;

		alpha = sqrt(N*m*m*PI*PI/D22/a/a+(D12*D12-D11*D22+4*D12*D66+4*D66*D66)/D22/D22*m*m*m*m*PI*PI*PI*PI/a/a/a/a)+(D12+2.*D66)/D22*m*m*PI*PI/a/a;
		beta = sqrt(N*m*m*PI*PI/D22/a/a+(D12*D12-D11*D22+4*D12*D66+4*D66*D66)/D22/D22*m*m*m*m*PI*PI*PI*PI/a/a/a/a)-(D12+2.*D66)/D22*m*m*PI*PI/a/a;
		alpha = sqrt(alpha);
		beta = sqrt(beta);

		temp1 =alpha*(a*a*alpha*alpha*D22-(D12+4.*D66)*m*m*PI*PI)*(a*a*beta*beta*D22+D12*m*m*PI*PI)*tan(b*beta);
		temp2 =beta*(a*a*alpha*alpha*D22-D12*m*m*PI*PI)*(a*a*beta*beta*D22+(D12+4.*D66)*m*m*PI*PI)*tanh(b*alpha);

		target = temp1/temp2;
		k0=target/N;

		

		pace = 10;
		count = 0 ;

		for(count=0;count<1000;count++)
		{
			
			N = N + (1.-target)/k0/pace;

			alpha = sqrt(N*m*m*PI*PI/D22/a/a+(D12*D12-D11*D22+4*D12*D66+4*D66*D66)/D22/D22*m*m*m*m*PI*PI*PI*PI/a/a/a/a)+(D12+2.*D66)/D22*m*m*PI*PI/a/a;
			beta  = sqrt(N*m*m*PI*PI/D22/a/a+(D12*D12-D11*D22+4*D12*D66+4*D66*D66)/D22/D22*m*m*m*m*PI*PI*PI*PI/a/a/a/a)-(D12+2.*D66)/D22*m*m*PI*PI/a/a;
			
			
			alpha = sqrt(alpha);
			beta = sqrt(beta);

			temp1 =alpha*(a*a*alpha*alpha*D22-(D12+4.*D66)*m*m*PI*PI)*(a*a*beta*beta*D22+D12*m*m*PI*PI)*tan(b*beta);
			temp2 =beta*(a*a*alpha*alpha*D22-D12*m*m*PI*PI)*(a*a*beta*beta*D22+(D12+4.*D66)*m*m*PI*PI)*tanh(b*alpha);

			target = temp1/temp2;
			
			k0 = target/N;

			if(fabs(target-1.)<1.0e-3)
				break;
			
		}

		if(count>=999)
		{
			printf("fail to converge\n");
			return(BucklingLoad);
		}



		if (N<can)
		{
			BucklingLoad.CrBuckleLoad.InPlane.Nxx=-N;
			BucklingLoad.m=m;
			can = N;
		}

	}


	BucklingLoad.CrBuckleStrain = Disp(plate, BucklingLoad.CrBuckleLoad);
	return(BucklingLoad);
}




struct PlateBuckleState PlateBuckling_UC_F_SBSB_Ortho(struct Laminate plate, double EA, double EI, double GJ)  
{

	struct PlateBuckleState BucklingLoad ;
	struct PlateBuckleState Candidate ;
	double m,  t, alpha, beta, a, b, AR, D11, D12, D22, D66, can, N0, N1, N2, temp1, temp2, sigma, term1, term2, f0, f1, f2;
	int i;
	BucklingLoad = PlateBuckleState();
    Candidate = PlateBuckleState();

	a =  plate.a;
	b =  plate.b;
	AR =a/b;
	D11 = plate.stack.abdm.D11;
	D12 = plate.stack.abdm.D12;
	D22 = plate.stack.abdm.D22;
	D66 = plate.stack.abdm.D66;

	can = -1.e20;

	for (m=1.;m<10.;m=m+1.)
	{

		N2 = -1.e20;
		
		for(i=0;i<1000;i++)
		{
			t = 2./1000.*((double)(i));
			N0 = - PI*PI/a/a/(m*m)*(m*m*m*m*D11 + 2.0*m*m*AR*AR*(D12 + 2.0*D66)*t*t + AR*AR*AR*AR*D22*t*t*t*t);

			alpha = sqrt(sqrt(-N0*m*m*PI*PI/D22/a/a+(D12*D12-D11*D22+4.*D12*D66+4.*D66*D66)/D22/D22*m*m*m*m*PI*PI*PI*PI/a/a/a/a)+(D12+2.*D66)/D22*m*m*PI*PI/a/a);
			beta = sqrt(sqrt(-N0*m*m*PI*PI/D22/a/a+(D12*D12-D11*D22+4.*D12*D66+4.*D66*D66)/D22/D22*m*m*m*m*PI*PI*PI*PI/a/a/a/a)-(D12+2.*D66)/D22*m*m*PI*PI/a/a);

			sigma = N0/(ElastMod(plate))/PlateThick(plate)*EA;
		
			temp1 = (a*a*a*a*a*a)*(beta*beta*beta*beta)*D22*D22 + 2.*(a*a*a*a)*(beta*beta)*D12*D22*(m*m)*(PI*PI)+EI*(-GJ)*(m*PI)*(m*PI)*(m*PI)*(m*PI)*(m*PI)*(m*PI) + (a*a)*(m*PI)*(m*PI)*(m*PI)*(m*PI)*(D12*D12+(-GJ)*sigma);
			temp2 = (a*a*a*a*a*a)*(alpha*alpha*alpha*alpha)*D22*D22 - 2.*(a*a*a*a)*(alpha*alpha)*D12*D22*(m*m)*(PI*PI)+EI*(-GJ)*(m*PI)*(m*PI)*(m*PI)*(m*PI)*(m*PI)*(m*PI) + (a*a)*(m*PI)*(m*PI)*(m*PI)*(m*PI)*(D12*D12+(-GJ)*sigma);


			term1 = alpha*sinh(alpha*b/2.)*cos(b*beta/2.)*temp1 + beta*temp2*cosh(alpha*b/2.)*sin(b*beta/2.) - a*a*D22*m*m*PI*PI*(alpha*alpha+beta*beta)*(EI*m*m*PI*PI + a*a*sigma)*cosh(alpha*b/2)*cos(b*beta/2.) - a*a*a*a*alpha*beta*(alpha*alpha+beta*beta)*D22*(-GJ)*m*m*PI*PI*sinh(alpha*b/2.)*sin(b*beta/2.);
			term2 = beta*temp2*sinh(alpha*b/2.)*cos(b*beta/2.) - alpha*temp1*cosh(alpha*b/2.)*sin(b*beta/2.) + a*a*D22*m*m*PI*PI*(alpha*alpha+beta*beta)*(EI*m*m*PI*PI + a*a*sigma)*sinh(alpha*b/2)*sin(b*beta/2.) - a*a*a*a*alpha*beta*(alpha*alpha+beta*beta)*D22*(-GJ)*m*m*PI*PI*cosh(alpha*b/2.)*cos(b*beta/2.);


			f0 = term1*term2;

			t = 2./1000.*((double)(i+1));
			
			N1 = - PI*PI/a/a/(m*m)*(m*m*m*m*D11 + 2.0*m*m*AR*AR*(D12 + 2.0*D66)*t*t + AR*AR*AR*AR*D22*t*t*t*t);

			alpha = sqrt(sqrt(-N1*m*m*PI*PI/D22/a/a+(D12*D12-D11*D22+4*D12*D66+4*D66*D66)/D22/D22*m*m*m*m*PI*PI*PI*PI/a/a/a/a)+(D12+2.*D66)/D22*m*m*PI*PI/a/a);
			beta = sqrt(sqrt(-N1*m*m*PI*PI/D22/a/a+(D12*D12-D11*D22+4*D12*D66+4*D66*D66)/D22/D22*m*m*m*m*PI*PI*PI*PI/a/a/a/a)-(D12+2.*D66)/D22*m*m*PI*PI/a/a);

			sigma = N1/(ElastMod(plate))/PlateThick(plate)*EA;
		
			temp1 = (a*a*a*a*a*a)*(beta*beta*beta*beta)*D22*D22 + 2.*(a*a*a*a)*(beta*beta)*D12*D22*(m*m)*(PI*PI)+EI*(-GJ)*(m*PI)*(m*PI)*(m*PI)*(m*PI)*(m*PI)*(m*PI) + (a*a)*(m*PI)*(m*PI)*(m*PI)*(m*PI)*(D12*D12+(-GJ)*sigma);
			temp2 = (a*a*a*a*a*a)*(alpha*alpha*alpha*alpha)*D22*D22 - 2.*(a*a*a*a)*(alpha*alpha)*D12*D22*(m*m)*(PI*PI)+EI*(-GJ)*(m*PI)*(m*PI)*(m*PI)*(m*PI)*(m*PI)*(m*PI) + (a*a)*(m*PI)*(m*PI)*(m*PI)*(m*PI)*(D12*D12+(-GJ)*sigma);


			term1 = alpha*sinh(alpha*b/2.)*cos(b*beta/2.)*temp1 + beta*cosh(alpha*b/2.)*sin(b*beta/2.)*temp2 - a*a*D22*m*m*PI*PI*(alpha*alpha+beta*beta)*(EI*m*m*PI*PI + a*a*sigma)*cosh(alpha*b/2)*cos(b*beta/2.) - a*a*a*a*alpha*beta*(alpha*alpha+beta*beta)*D22*(-GJ)*m*m*PI*PI*sinh(alpha*b/2.)*sin(b*beta/2.);
			term2 = beta*temp2*sinh(alpha*b/2.)*cos(b*beta/2.) - alpha*temp1*cosh(alpha*b/2.)*sin(b*beta/2.) + a*a*D22*m*m*PI*PI*(alpha*alpha+beta*beta)*(EI*m*m*PI*PI + a*a*sigma)*sinh(alpha*b/2)*sin(b*beta/2.) - a*a*a*a*alpha*beta*(alpha*alpha+beta*beta)*D22*(-GJ)*m*m*PI*PI*cosh(alpha*b/2.)*cos(b*beta/2.);


			f1 = term1*term2;

			if (f0*f1<0.)
			{
				N2 = (N1+N0)/2.;
				break;
			}


		}


		if (f1*f0>0.)
		{
			printf("error, m=%f\n",m);
			return(BucklingLoad);
		}


		if (N2>can)
		{
			BucklingLoad.CrBuckleLoad.InPlane.Nxx = N2;
			BucklingLoad.m=m;
			can = N2;
		}

	}


	BucklingLoad.CrBuckleStrain = Disp(plate, BucklingLoad.CrBuckleLoad);
	return(BucklingLoad);
}




struct PlateBuckleState PlateBuckling_UC_F_SESF_Ortho(struct Laminate plate, double C0)
{

	struct PlateBuckleState BucklingLoad ;
	struct PlateBuckleState Candidate;

	double m, alpha, beta, a, b, D11, D12, D22, D66, N, can, target, temp1, temp2, temp3, temp4, k0;
	double k1,k2,k3,k4;
	BucklingLoad = PlateBuckleState();
    Candidate = PlateBuckleState();

	a =  plate.a;
	b =  plate.b;
	D11 = plate.stack.abdm.D11;
	D12 = plate.stack.abdm.D12;
	D22 = plate.stack.abdm.D22;
	D66 = plate.stack.abdm.D66;

	can = 1.e20;

	for (m=1.;m<=20.;m=m+1.)
	{
		N = D11*m*m*PI*PI/a/a;


		do
		{
			alpha = sqrt(N*m*m*PI*PI/D22/a/a+(D12*D12-D11*D22+4*D12*D66+4*D66*D66)/D22/D22*m*m*m*m*PI*PI*PI*PI/a/a/a/a)+(D12+2.*D66)/D22*m*m*PI*PI/a/a;
			beta = sqrt(N*m*m*PI*PI/D22/a/a+(D12*D12-D11*D22+4*D12*D66+4*D66*D66)/D22/D22*m*m*m*m*PI*PI*PI*PI/a/a/a/a)-(D12+2.*D66)/D22*m*m*PI*PI/a/a;
			alpha = sqrt(alpha);
			beta = sqrt(beta);

			temp1 = -a*a*alpha*alpha*D22+(D12+4.*D66)*m*m*PI*PI;
	 	    temp2 = a*a*beta*beta*D22 + (D12+4.*D66)*m*m*PI*PI;
			temp3 = a*a*alpha*alpha*D22-D12*m*m*PI*PI;
			temp4 = a*a*beta*beta*D22+D12*m*m*PI*PI;

			k1 = (1.+exp(2.*alpha*b))*temp1-2.*exp(alpha*b)*temp2*cos(b*beta);
			k1=k1*beta;

			k2 = alpha*C0*m*m*PI*PI*temp4*cos(b*beta)+temp3*(alpha*C0*m*m*PI*PI*cosh(alpha*b)+a*a*(alpha*alpha+beta*beta)*D22*sinh(alpha*b));

			k3 = beta*(-1.+exp(2.*alpha*b))*temp3 + 2.*alpha*exp(alpha*b)*temp4*sin(b*beta);
			k4 = beta*C0*m*m*PI*PI*temp2*sin(b*beta)+temp1*(a*a*(alpha*alpha+beta*beta)*D22*cosh(alpha*b)+alpha*C0*m*m*PI*PI*sinh(alpha*b));

			//target = k1*k2/k3/k4;
			target = k3*k4/k1/k2;

			k0=target/N;

			//N = N*target;
			N = N + (1-target)/k0;

		}while(fabs(target-1.)>1.0e-6);

		if (N<can)
		{
			BucklingLoad.CrBuckleLoad.InPlane.Nxx=-N;
			BucklingLoad.m=m;
			can = N;
		}

	}
	BucklingLoad.CrBuckleStrain = Disp(plate, BucklingLoad.CrBuckleLoad);


	return(BucklingLoad);
}

struct PlateBuckleState PlateBuckling_UC_F_SCSC(struct Laminate plate)
{

	struct PlateBuckleState BucklingLoad ;
	struct PlateBuckleState Candidate;

	double m, alpha, beta, a, b, D11, D12, D22, D66, N0, N1, N2, can, f1, f0, f2, k0, temp1, temp2;


	BucklingLoad = PlateBuckleState();
    Candidate = PlateBuckleState();

	a =  plate.a;
	b =  plate.b;
	D11 = plate.stack.abdm.D11;
	D12 = plate.stack.abdm.D12;
	D22 = plate.stack.abdm.D22;
	D66 = plate.stack.abdm.D66;

	can = -1.e20;



	for (m=1.;m<=10.;m=m+1.)
	{
		N0 = -PI*PI/a/a*(D11*m*m + 2*(D12+2*D66)*(a/b)*(a/b) + D22*(a/b)*(a/b)*(a/b)*(a/b)/m/m);
		N1 = 2*N0;

		alpha = sqrt(sqrt(-a*a*N0*D22 +D12*D12*m*m*PI*PI-D11*D22*m*m*PI*PI+4*D12*D66*m*m*PI*PI+4*D66*D66*m*m*PI*PI)*PI*m/a/a/D22 + (D12+2.*D66)/D22*m*m*PI*PI/a/a);
		beta = sqrt(sqrt(-a*a*N0*D22 +D12*D12*m*m*PI*PI-D11*D22*m*m*PI*PI+4*D12*D66*m*m*PI*PI+4*D66*D66*m*m*PI*PI)*PI*m/a/a/D22 - (D12+2.*D66)/D22*m*m*PI*PI/a/a);
		f0 = 2.*(1.-cos(beta*b)*cosh(alpha*b)) + (alpha/beta-beta/alpha)*sin(beta*b)*sinh(alpha*b);

		do
		{
			N1=N1*2.;
			alpha = sqrt(sqrt(-a*a*N1*D22 +D12*D12*m*m*PI*PI-D11*D22*m*m*PI*PI+4*D12*D66*m*m*PI*PI+4*D66*D66*m*m*PI*PI)*PI*m/a/a/D22 + (D12+2.*D66)/D22*m*m*PI*PI/a/a);
			beta = sqrt(sqrt(-a*a*N1*D22 +D12*D12*m*m*PI*PI-D11*D22*m*m*PI*PI+4*D12*D66*m*m*PI*PI+4*D66*D66*m*m*PI*PI)*PI*m/a/a/D22 - (D12+2.*D66)/D22*m*m*PI*PI/a/a);
			f1 = 2.*(1.-cos(beta*b)*cosh(alpha*b)) + (alpha/beta-beta/alpha)*sin(beta*b)*sinh(alpha*b);
		
		}while((f0*f1)>0);

		if ((f0*f1)>0)
		{
			printf("The bisectional method will not converge\n");
			return(BucklingLoad);
		}
		
		do
		{
			N2 = (N0+N1)/2.;
			alpha = sqrt(sqrt(-a*a*N2*D22 +D12*D12*m*m*PI*PI-D11*D22*m*m*PI*PI+4*D12*D66*m*m*PI*PI+4*D66*D66*m*m*PI*PI)*PI*m/a/a/D22 + (D12+2.*D66)/D22*m*m*PI*PI/a/a);
			beta = sqrt(sqrt(-a*a*N2*D22 +D12*D12*m*m*PI*PI-D11*D22*m*m*PI*PI+4*D12*D66*m*m*PI*PI+4*D66*D66*m*m*PI*PI)*PI*m/a/a/D22 - (D12+2.*D66)/D22*m*m*PI*PI/a/a);
			f2 = 2.*(1.-cos(beta*b)*cosh(alpha*b)) + (alpha/beta-beta/alpha)*sin(beta*b)*sinh(alpha*b);

			if ((f2*f1)>0)
			{
				f1=f2;
				N1=N2;
			}
			else
			{
				f0=f2;
				N0=N2;

			}

		}while((fabs(f2)>1.0e-6)&&(fabs(N1-N0)>1.));

		
		if (N2>can)
		{
			BucklingLoad.CrBuckleLoad.InPlane.Nxx = N2;
			BucklingLoad.m=m;
			can = N2;
		}


	}

	BucklingLoad.CrBuckleStrain = Disp(plate, BucklingLoad.CrBuckleLoad);


	return(BucklingLoad);

}


struct PlateBuckleState PlateBuckling_UC_F_SESE_Ortho(struct Laminate plate, double C0)
{
	struct PlateBuckleState BucklingLoad ;
	struct PlateBuckleState Candidate;

	double m, alpha, beta, a, b, D11, D12, D22, D66, N0, N1, N2, can, f1, f0, f2, k0, temp1, temp2, t0, t1, t2;

	BucklingLoad = PlateBuckleState();
    Candidate = PlateBuckleState();

	a =  plate.a;
	b =  plate.b;
	D11 = plate.stack.abdm.D11;
	D12 = plate.stack.abdm.D12;
	D22 = plate.stack.abdm.D22;
	D66 = plate.stack.abdm.D66;

	can = -1.e20;

	

	for (m=1.;m<=10.;m=m+1.)
	{
		N0 = -PI*PI/a/a*(D11*m*m + 2*(D12+2*D66)*(a/b)*(a/b) + D22*(a/b)*(a/b)*(a/b)*(a/b)/m/m);
		alpha = sqrt(sqrt(-a*a*N0*D22 +D12*D12*m*m*PI*PI-D11*D22*m*m*PI*PI+4*D12*D66*m*m*PI*PI+4*D66*D66*m*m*PI*PI)*PI*m/a/a/D22 + (D12+2.*D66)/D22*m*m*PI*PI/a/a);
		beta = sqrt(sqrt(-a*a*N0*D22 +D12*D12*m*m*PI*PI-D11*D22*m*m*PI*PI+4*D12*D66*m*m*PI*PI+4*D66*D66*m*m*PI*PI)*PI*m/a/a/D22 - (D12+2.*D66)/D22*m*m*PI*PI/a/a);
		t0 = a*a*a*a*(alpha*alpha + beta*beta)*(alpha*alpha + beta*beta)*sin(beta*b)*sinh(alpha*b)*D22*D22;
		t1 = 2.*D22*a*a*m*m*PI*PI*(alpha*alpha + beta*beta)*(alpha*cosh(alpha*b)*sin(beta*b) - beta*sinh(alpha*b)*cos(beta*b));
		t2 =(m*m*m*m)*(PI*PI*PI*PI)*alpha*beta*(2.*(1.-cos(beta*b)*cosh(alpha*b)) + (alpha/beta-beta/alpha)*sin(beta*b)*sinh(alpha*b));
		f0 = t0 + C0*t1 + C0*C0*t2 ;
		
		

		N1 = -PI*PI/a/a*(D11*m*m + 2*(D12+2*D66)*(a/b)*(a/b)*4. + D22*(a/b)*(a/b)*(a/b)*(a/b)/m/m*16.);
		alpha = sqrt(sqrt(-a*a*N1*D22 +D12*D12*m*m*PI*PI-D11*D22*m*m*PI*PI+4*D12*D66*m*m*PI*PI+4*D66*D66*m*m*PI*PI)*PI*m/a/a/D22 + (D12+2.*D66)/D22*m*m*PI*PI/a/a);
		beta = sqrt(sqrt(-a*a*N1*D22 +D12*D12*m*m*PI*PI-D11*D22*m*m*PI*PI+4*D12*D66*m*m*PI*PI+4*D66*D66*m*m*PI*PI)*PI*m/a/a/D22 - (D12+2.*D66)/D22*m*m*PI*PI/a/a);
		t0 = a*a*a*a*(alpha*alpha + beta*beta)*(alpha*alpha + beta*beta)*sin(beta*b)*sinh(alpha*b)*D22*D22;
		t1 = 2.*D22*a*a*m*m*PI*PI*(alpha*alpha + beta*beta)*(alpha*cosh(alpha*b)*sin(beta*b) - beta*sinh(alpha*b)*cos(beta*b));
		t2 =(m*m*m*m)*(PI*PI*PI*PI)*alpha*beta*(2.*(1.-cos(beta*b)*cosh(alpha*b)) + (alpha/beta-beta/alpha)*sin(beta*b)*sinh(alpha*b));
		f1 = t0 + C0*t1 + C0*C0*t2 ;

	
		
		do
		{
			
			N2 = (N0+N1)/2.;
			alpha = sqrt(sqrt(-a*a*N2*D22 +D12*D12*m*m*PI*PI-D11*D22*m*m*PI*PI+4*D12*D66*m*m*PI*PI+4*D66*D66*m*m*PI*PI)*PI*m/a/a/D22 + (D12+2.*D66)/D22*m*m*PI*PI/a/a);
			beta = sqrt(sqrt(-a*a*N2*D22 +D12*D12*m*m*PI*PI-D11*D22*m*m*PI*PI+4*D12*D66*m*m*PI*PI+4*D66*D66*m*m*PI*PI)*PI*m/a/a/D22 - (D12+2.*D66)/D22*m*m*PI*PI/a/a);
			t0 = a*a*a*a*(alpha*alpha + beta*beta)*(alpha*alpha + beta*beta)*sin(beta*b)*sinh(alpha*b)*D22*D22;
			t1 = 2.*D22*a*a*m*m*PI*PI*(alpha*alpha + beta*beta)*(alpha*cosh(alpha*b)*sin(beta*b) - beta*sinh(alpha*b)*cos(beta*b));
			t2 =(m*m*m*m)*(PI*PI*PI*PI)*alpha*beta*(2.*(1.-cos(beta*b)*cosh(alpha*b)) + (alpha/beta-beta/alpha)*sin(beta*b)*sinh(alpha*b));
			f2 = t0 + C0*t1 + C0*C0*t2 ;

			if ((f2*f1)>0)
			{
				f1=f2;
				N1=N2;
			}
			else
			{
				f0=f2;
				N0=N2;

			}

			
			
		
		}while(fabs(N1-N0)>0.001);

		if (N2>can)
		{
			BucklingLoad.CrBuckleLoad.InPlane.Nxx = N2;
			BucklingLoad.m=m;
			can = N2;
		}
	}

	BucklingLoad.CrBuckleStrain = Disp(plate, BucklingLoad.CrBuckleLoad);


	return(BucklingLoad);

}







struct PlateBuckleState PlateBuckling_BC_F_SSSS_Ortho(struct Laminate plate, double k)
{

	struct PlateBuckleState BucklingLoad ;
	struct PlateBuckleState Candidate;
	double AR,theta, m, n;

	BucklingLoad=PlateBuckleState();
	Candidate=PlateBuckleState();
	BucklingLoad.CrBuckleLoad.InPlane.Nxx=1.0E20;


	AR=plate.a/plate.b;

	for (m=1.;m<=20.;m=m+1.)
	{
		for (n=1.;n<20.;n=n+1.)
		{
			theta = m*m*m*m*plate.stack.abdm.D11 + 2.0*m*m*n*n*AR*AR*(plate.stack.abdm.D12 + 2.0*plate.stack.abdm.D66) + n*n*n*n*AR*AR*AR*AR*plate.stack.abdm.D22;
			Candidate.CrBuckleLoad.InPlane.Nxx = PI*PI*theta/plate.a/plate.a/(m*m+k*n*n*AR*AR);
			Candidate.m = (int)m;
			Candidate.n = (int)n;

			if(Candidate.CrBuckleLoad.InPlane.Nxx < BucklingLoad.CrBuckleLoad.InPlane.Nxx)
				BucklingLoad = Candidate;
		}
	}

	BucklingLoad.CrBuckleLoad.InPlane.Nxx = - BucklingLoad.CrBuckleLoad.InPlane.Nxx;// change into negative value
	BucklingLoad.CrBuckleLoad.InPlane.Nyy = BucklingLoad.CrBuckleLoad.InPlane.Nxx*k;
	BucklingLoad.CrBuckleStrain = Disp(plate, BucklingLoad.CrBuckleLoad);
	return(BucklingLoad);
}



struct PlateBuckleState PlateBuckling_S_F_SSSS_Ortho(struct Laminate plate)
{
	/*simply supported is assumed*/
	/*double A, B;*/
	double K, beta, theta, BucklingLoad;
	struct PlateBuckleState BLoad;

	BLoad = PlateBuckleState();

    theta = sqrt(plate.stack.abdm.D11*plate.stack.abdm.D22)/(plate.stack.abdm.D12+2.0*plate.stack.abdm.D66);
	beta = sqrt(sqrt(plate.stack.abdm.D11/plate.stack.abdm.D22))*plate.b/plate.a;
    K = (3.32 + 2.17/theta - 0.163/theta/theta) + MIN(beta , 1.0/beta)*MIN(beta , 1.0/beta)*(1.54 + 2.36/theta + 0.1/theta/theta);
    BucklingLoad = K * PI *PI * sqrt(sqrt(plate.stack.abdm.D11*plate.stack.abdm.D22*plate.stack.abdm.D22*plate.stack.abdm.D22))/plate.b/plate.b;
	if(beta>1)
		BucklingLoad = K * PI *PI * sqrt(sqrt(plate.stack.abdm.D22*plate.stack.abdm.D11*plate.stack.abdm.D11*plate.stack.abdm.D11))/plate.a/plate.a;


	BLoad.CrBuckleLoad.InPlane.Nxy = BucklingLoad;
	BLoad.CrBuckleStrain = Disp(plate, BLoad.CrBuckleLoad);

	return(BLoad);
}

struct PlateBuckleState PlateBuckling_S_F_CCCC_Ortho(struct Laminate plate)
{
	/*simply supported is assumed*/
	/*double A, B;*/
	double K, beta, theta, BucklingLoad;
	struct PlateBuckleState BLoad;

    BLoad = PlateBuckleState();
    theta = sqrt(plate.stack.abdm.D11*plate.stack.abdm.D22)/(plate.stack.abdm.D12+2.0*plate.stack.abdm.D66);
	beta = sqrt(sqrt(plate.stack.abdm.D11/plate.stack.abdm.D22))*plate.b/plate.a;
    K = (6.21 + 2.36/theta - 0.116/theta/theta) + MIN(beta , 1.0/beta)*MIN(beta , 1.0/beta)*(3.21+ 2.11/theta + 0.045/theta/theta);
    BucklingLoad = K * PI *PI * sqrt(sqrt(plate.stack.abdm.D11*plate.stack.abdm.D22*plate.stack.abdm.D22*plate.stack.abdm.D22))/plate.b/plate.b;
	if(beta>1)
		BucklingLoad = K * PI *PI * sqrt(sqrt(plate.stack.abdm.D22*plate.stack.abdm.D11*plate.stack.abdm.D11*plate.stack.abdm.D11))/plate.a/plate.a;

	BLoad.CrBuckleLoad.InPlane.Nxy = BucklingLoad;
	BLoad.CrBuckleStrain = Disp(plate, BLoad.CrBuckleLoad);

	return(BLoad);
}


struct PlateBuckleState PlateBuckling_S_F_SCSC_Ortho(struct Laminate plate)
{
	/*simply supported is assumed*/
	/*double A, B;*/
	double K, beta, theta, BucklingLoad;
	struct PlateBuckleState BLoad ;

    BLoad = PlateBuckleState();

    theta = sqrt(plate.stack.abdm.D11*plate.stack.abdm.D22)/(plate.stack.abdm.D12+2.0*plate.stack.abdm.D66);
	beta = sqrt(sqrt(plate.stack.abdm.D11/plate.stack.abdm.D22))*plate.b/plate.a;
    K = (6.21 + 2.36/theta - 0.116/theta/theta) + MIN(beta , 1.0/beta)*MIN(beta , 1.0/beta)*(3.21+ 2.11/theta + 0.045/theta/theta);
    K=8.0;

	BucklingLoad = K * PI *PI * sqrt(sqrt(plate.stack.abdm.D11*plate.stack.abdm.D22*plate.stack.abdm.D22*plate.stack.abdm.D22))/plate.b/plate.b;
	if(beta>1)
		BucklingLoad = K * PI *PI * sqrt(sqrt(plate.stack.abdm.D22*plate.stack.abdm.D11*plate.stack.abdm.D11*plate.stack.abdm.D11))/plate.a/plate.a;

	BLoad.CrBuckleLoad.InPlane.Nxy = BucklingLoad;
	BLoad.CrBuckleStrain = Disp(plate, BLoad.CrBuckleLoad);
	return(BLoad);
}


struct PlateBuckleState PlateBuckling_NS_F_SSSS_Aniso(struct Laminate plate, struct PlateLoad pl)
{
	int mnlim, ind_mn, ind_ij, m, n,  i, j, num_odd, num_even, **ind;
	double a,b,D11,D12,D22,D66,D16,D26, K_d, P_d, K_off_d, P_off_d, lambda;
	struct PlateBuckleState BucklingLoad;
	BucklingLoad = PlateBuckleState();

	a = plate.a;
	b = plate.b;
	D11 = plate.stack.abdm.D11;
	D12 = plate.stack.abdm.D12;
	D22 = plate.stack.abdm.D22;
	D66 = plate.stack.abdm.D66;
	D16 = plate.stack.abdm.D16;
	D26 = plate.stack.abdm.D26;

	mnlim = 25;
	num_odd=0;
	num_even=0;


	ind = (int**) malloc(sizeof(int*)*mnlim*mnlim);




	for (i=0;i<mnlim*mnlim;i++)
	{
		ind[i] = (int*) malloc(sizeof(int)*mnlim*mnlim);
		for (j=0;j<mnlim*mnlim;j++)
		{
				ind[i][j]=0;
		}
	}


	for (m=1;m<=mnlim ;m++)
	{
		for (n=1;n<=mnlim ;n++)
		{
			if((m+n)%2==0)
			{
				num_even =  num_even + 1;
				ind[m-1][n-1] = num_even;
			}
			else
				{
					num_odd =  num_odd + 1;
					ind[m-1][n-1] =  - num_odd;
				}
		}

	}

	MatrixXd  K_eig_odd(num_odd,  num_odd);
	MatrixXd  P_eig_odd(num_odd,  num_odd);

	MatrixXd  K_eig_even(num_even,  num_even);
	MatrixXd  P_eig_even(num_even,  num_even);

	for (i=0;i<num_odd;i++)
	{
		for (j=0;j<num_odd;j++)
		{
			K_eig_odd(i,j)=0.;
			P_eig_odd(i,j)=0.;
		}

	}

	for (i=0;i<num_even;i++)
	{
		for (j=0;j<num_even;j++)
		{
			K_eig_even(i,j)=0.;
			P_eig_even(i,j)=0.;
		}

	}
	for (m=1; m<=mnlim; m++)
	{
		for (n=1; n<=mnlim; n++)
		{
			ind_mn = ind[m-1][n-1];
			K_d = PI*PI*PI*PI*(b*b*b*b*m*m*m*m*D11+2.*a*a*b*b*m*m*n*n*(D12+2.*D66)+a*a*a*a*n*n*n*n*D22)/4./a/a/a/b/b/b;
			P_d = -b*pl.InPlane.Nxx*m*m*PI*PI/4./a - a*pl.InPlane.Nyy*n*n*PI*PI/4./b;

			if (ind_mn>0)
			{
				K_eig_even(ind_mn-1,  ind_mn-1) = K_d;
				P_eig_even(ind_mn-1,  ind_mn-1) = P_d;
			}
			else
			{
				K_eig_odd(-ind_mn-1,  -ind_mn-1) =  K_d;
				P_eig_odd(-ind_mn-1,  -ind_mn-1) =  P_d;
			}



			for (i=(m%2+1); i<=mnlim; i=i+2)
			{
				for (j=(n%2+1); j<=mnlim; j=j+2)
				{
					ind_ij = ind[i-1][j-1];

					K_off_d = -8.*m*n*i*j*PI*PI*(b*b*(m*m+i*i)*D16+a*a*(n*n+j*j)*D26)/a/a/b/b/(m-i)/(m+i)/(n-j)/(n+j);

					P_off_d = 8.*m*n*i*j/(m-i)/(m+i)/(n-j)/(n+j)*pl.InPlane.Nxy;

					if (ind_mn>0)
					{
						K_eig_even(ind_mn-1,  ind_ij-1) = K_off_d;


						P_eig_even(ind_mn-1,  ind_ij-1) = P_off_d;


					}
					else
					{
						K_eig_odd(-ind_mn-1,  -ind_ij-1) =  K_off_d;


						P_eig_odd(-ind_mn-1,  -ind_ij-1) =  P_off_d;

					}
				}
			}



		}
	}




	GeneralizedSelfAdjointEigenSolver<MatrixXd> es1(P_eig_odd, K_eig_odd);
	GeneralizedSelfAdjointEigenSolver<MatrixXd> es2(P_eig_even, K_eig_even);

	lambda = MAX(es1.eigenvalues()[num_odd-1], es2.eigenvalues()[num_even-1]);

	BucklingLoad.CrBuckleLoad.InPlane.Nxx = pl.InPlane.Nxx/lambda;
	BucklingLoad.CrBuckleLoad.InPlane.Nxy = pl.InPlane.Nxy/lambda;
	BucklingLoad.CrBuckleLoad.InPlane.Nyy = pl.InPlane.Nyy/lambda;

	BucklingLoad.CrBuckleStrain = Disp(plate, BucklingLoad.CrBuckleLoad);

	if (lambda<0)
	{
		lambda=1.0e-20;
		BucklingLoad.CrBuckleLoad.InPlane.Nxx = pl.InPlane.Nxx/lambda;
		BucklingLoad.CrBuckleLoad.InPlane.Nxy = pl.InPlane.Nxy/lambda;
		BucklingLoad.CrBuckleLoad.InPlane.Nyy = pl.InPlane.Nyy/lambda;
		printf("The plate will not buckle in that loading direction\n");
	}



	return(BucklingLoad);

}

struct PlateBuckleState ShearBuckling_SSSS_Ortho(struct Laminate plate)
{
	int mnlim, m, n , i,  j, num_odd, num_even, ind_mn, ind_ij, **ind, x, y;
	double AR, **M, *V, BuckForce;
	struct PlateBuckleState BucklingLoad;

	BucklingLoad = PlateBuckleState();
	FILE *fp;

	fp=fopen("M.txt","w");

	BuckForce = 0.;
	mnlim = 10;
	num_odd=0;
	num_even=0;

	AR=plate.a/plate.b;


	V = (double*) malloc(sizeof(double)*mnlim*mnlim);
	M = (double**) malloc(sizeof(double*)*mnlim*mnlim);
	ind = (int**) malloc(sizeof(int*)*mnlim*mnlim);

	for (i=0;i<mnlim*mnlim;i++)
	{
		M[i] = (double*) malloc(sizeof(double)*mnlim*mnlim);
		ind[i] = (int*) malloc(sizeof(int)*mnlim*mnlim);
		V[i] = 0.0;
		for (j=0;j<mnlim*mnlim;j++)
		{
				M[i][j]=0.0;
				ind[i][j]=0;
		}
	}


	for (m=1;m<=mnlim ;m++)
	{
		for (n=1;n<=mnlim ;n++)
		{
			ind_mn = (m-1) + (n-1)*mnlim;
			V[ind_mn] = PI*PI*PI*PI*(m*m*m*m*plate.stack.abdm.D11 + 2.0*m*m*n*n*AR*AR*(plate.stack.abdm.D12 + 2.0*plate.stack.abdm.D66) + n*n*n*n*AR*AR*AR*AR*plate.stack.abdm.D22);


			if((m+n)%2==0)
			{

				num_even =  num_even + 1;
				ind[m-1][n-1] = num_even;
			}
			else
				{

					num_odd =  num_odd + 1;
					ind[m-1][n-1] =  - num_odd;
				}


			for (i=1;i<=mnlim;i++)
			{
				for (j=1;j<=mnlim;j++)
				{
					ind_ij = (i-1) + (j-1)*mnlim;
					if (((m+i)%2!=0)&&((n+j)%2!=0))
						M[ind_mn][ind_ij]= 32.0*m*n*i*j*AR*AR*AR*plate.b*plate.b/(m*m-i*i)/(n*n-j*j);
				}
			}


		}

	}


	for (i=0;i<mnlim*mnlim;i++)
	{
		V[i]=1.0/sqrt(V[i]);
	}
	for (i=0;i<mnlim*mnlim;i++)
	{
		for (j=0;j<mnlim*mnlim;j++)
		{
			M[i][j] = M[i][j] * V[i] * V[j];
		}
	}


	MatrixXd  M_eig_odd(num_odd,  num_odd);
	MatrixXd  M_eig_even(num_even,  num_even);


	for (m=1; m<=mnlim; m++)
	{
		for (n=1; n<=mnlim; n++)
		{
			for (i=1; i<=mnlim; i++)
			{
				for (j=1; j<=mnlim; j++)
				{
					x = ind[m-1][n-1];
					y = ind[i-1][j-1];

					if ((x>0)&&(y>0))
						M_eig_even(x-1, y-1) = M[(m-1) + (n-1)*mnlim][(i-1) + (j-1)*mnlim];
					else if  ((x<0)&&(y<0))
						M_eig_odd(-x-1, -y-1) = M[(m-1) + (n-1)*mnlim][(i-1) + (j-1)*mnlim];
				}
			}

		}
	}
	for (i=0;i<num_odd;i++)
	{
		for (j=0;j<num_odd;j++)
		{
			fprintf(fp,"%e, ",M_eig_odd(i,j));
		}
		fprintf(fp,"\n");
	}

	fclose(fp);

	SelfAdjointEigenSolver<MatrixXd> es_odd, es_even;
	es_odd.compute(M_eig_odd);
	es_even.compute(M_eig_even);

	for (i=0;i<num_odd;i++)
	{
		if(es_odd.eigenvalues()[i]>BuckForce)
			BuckForce=es_odd.eigenvalues()[i];
	}

	for (i=0;i<num_even;i++)
	{
		if(es_even.eigenvalues()[i]>BuckForce)
			BuckForce=es_even.eigenvalues()[i];
	}

	BuckForce = 1.0/BuckForce;

	BucklingLoad.CrBuckleLoad.InPlane.Nxy = BuckForce;
	BucklingLoad.CrBuckleStrain = Disp(plate, BucklingLoad.CrBuckleLoad);


	return(BucklingLoad);
}







struct PlateBuckleState NSBuckling(double Nxcrit, double Nxycrit, double NSAngle)
{
	double c, s, a, b, BuckForce;
	struct PlateBuckleState BLoad={0};


	s = sin(PI*NSAngle/180.);
	c = cos(PI*NSAngle/180.);



	if (fabs(s)<1.0e-6)//uniaxial tension or compression
	{
		
		
		if(c<0)
		{
			BLoad.CrBuckleLoad.InPlane.Nxy = 0.;
			BLoad.CrBuckleLoad.InPlane.Nxx = -fabs(Nxcrit);
			return(BLoad);
		}
		else
		{
			BLoad.CrBuckleLoad.InPlane.Nxy = 0.;
			BLoad.CrBuckleLoad.InPlane.Nxx = 1.0e14;
			return(BLoad);

		}
	}





    a = s*s/Nxycrit/Nxycrit;
	b = c/abs(Nxcrit);

	BuckForce = (b + sqrt(b*b + 4.0*a))/2.0/a; //if b>0, tension, BuckForce is amplified; elsewise BuckForce is diminished
	BLoad.CrBuckleLoad.InPlane.Nxx = BuckForce*c;
	BLoad.CrBuckleLoad.InPlane.Nxy = BuckForce*s;

	return(BLoad);
}










