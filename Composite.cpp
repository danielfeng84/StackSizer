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


#include "Composite.h"
#include <math.h>
#include <stdio.h>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

void ABDCalc(struct Stack *dims)
{
		struct ABDM abdmatrix;
		int i;
		double Q11_global, Q12_global, Q16_global, Q22_global, Q26_global, Q66_global;
		double temp, plythichness, LayerP;
		struct LaminaQ Q;

		abdmatrix.A11=0.0;
		abdmatrix.A12=0.0;
		abdmatrix.A16=0.0;
		abdmatrix.A22=0.0;
		abdmatrix.A26=0.0;
		abdmatrix.A66=0.0;

		abdmatrix.B11=0.0;
		abdmatrix.B12=0.0;
		abdmatrix.B16=0.0;
		abdmatrix.B22=0.0;
		abdmatrix.B26=0.0;
		abdmatrix.B66=0.0;

		abdmatrix.D11=0.0;
		abdmatrix.D12=0.0;
		abdmatrix.D16=0.0;
		abdmatrix.D22=0.0;
		abdmatrix.D26=0.0;
		abdmatrix.D66=0.0;


	for (i=0;i<dims->NumLayers;i++)
	{


		plythichness = dims->ply[i].thick;

		LayerP = LayerHeight(*dims, i);

		Q = GlobalLaminaQ((dims->ply[i]), (double) dims->PlyDirection[i]);


		Q11_global = Q.Q11;
		Q12_global = Q.Q12;
		Q16_global = Q.Q16;
		Q22_global = Q.Q22;
		Q26_global = Q.Q26;
		Q66_global = Q.Q66;



		abdmatrix.A11 = abdmatrix.A11 + Q11_global*plythichness;
		abdmatrix.A12 = abdmatrix.A12 + Q12_global*plythichness;
		abdmatrix.A16 = abdmatrix.A16 + Q16_global*plythichness;
		abdmatrix.A22 = abdmatrix.A22 + Q22_global*plythichness;
		abdmatrix.A26 = abdmatrix.A26 + Q26_global*plythichness;
		abdmatrix.A66 = abdmatrix.A66 + Q66_global*plythichness;

		temp = plythichness*(2.0*LayerP)/2.0;
		abdmatrix.B11 = abdmatrix.B11 + Q11_global*temp;
		abdmatrix.B12 = abdmatrix.B12 + Q12_global*temp;
		abdmatrix.B16 = abdmatrix.B16 + Q16_global*temp;
		abdmatrix.B22 = abdmatrix.B22 + Q22_global*temp;
		abdmatrix.B26 = abdmatrix.B26 + Q26_global*temp;
		abdmatrix.B66 = abdmatrix.B66 + Q66_global*temp;

		temp = plythichness*(plythichness*plythichness/4.0 + 3.0*LayerP*LayerP)/3.0;
		abdmatrix.D11 = abdmatrix.D11 + Q11_global*temp;
		abdmatrix.D12 = abdmatrix.D12 + Q12_global*temp;
		abdmatrix.D16 = abdmatrix.D16 + Q16_global*temp;
		abdmatrix.D22 = abdmatrix.D22 + Q22_global*temp;
		abdmatrix.D26 = abdmatrix.D26 + Q26_global*temp;
		abdmatrix.D66 = abdmatrix.D66 + Q66_global*temp;
	}
	dims->abdm = abdmatrix;

}


struct ABDM ABD(struct Stack stack)
{
		struct ABDM abdmatrix;
		int i;

		double Q11_global, Q12_global, Q16_global, Q22_global, Q26_global, Q66_global;
		double temp, plythichness, LayerP;
		struct LaminaQ Q;

		abdmatrix.A11=0.0;
		abdmatrix.A12=0.0;
		abdmatrix.A16=0.0;
		abdmatrix.A22=0.0;
		abdmatrix.A26=0.0;
		abdmatrix.A66=0.0;

		abdmatrix.B11=0.0;
		abdmatrix.B12=0.0;
		abdmatrix.B16=0.0;
		abdmatrix.B22=0.0;
		abdmatrix.B26=0.0;
		abdmatrix.B66=0.0;

		abdmatrix.D11=0.0;
		abdmatrix.D12=0.0;
		abdmatrix.D16=0.0;
		abdmatrix.D22=0.0;
		abdmatrix.D26=0.0;
		abdmatrix.D66=0.0;


	for (i=0;i<stack.NumLayers;i++)
	{



		plythichness = stack.ply[i].thick;


		LayerP = LayerHeight(stack, i);

		Q = GlobalLaminaQ(stack.ply[i], (double) stack.PlyDirection[i]);


		Q11_global = Q.Q11;
		Q12_global = Q.Q12;
		Q16_global = Q.Q16;
		Q22_global = Q.Q22;
		Q26_global = Q.Q26;
		Q66_global = Q.Q66;



		abdmatrix.A11 = abdmatrix.A11 + Q11_global*plythichness;
		abdmatrix.A12 = abdmatrix.A12 + Q12_global*plythichness;
		abdmatrix.A16 = abdmatrix.A16 + Q16_global*plythichness;
		abdmatrix.A22 = abdmatrix.A22 + Q22_global*plythichness;
		abdmatrix.A26 = abdmatrix.A26 + Q26_global*plythichness;
		abdmatrix.A66 = abdmatrix.A66 + Q66_global*plythichness;

		temp = plythichness*(2.0*LayerP)/2.0;
		abdmatrix.B11 = abdmatrix.B11 + Q11_global*temp;
		abdmatrix.B12 = abdmatrix.B12 + Q12_global*temp;
		abdmatrix.B16 = abdmatrix.B16 + Q16_global*temp;
		abdmatrix.B22 = abdmatrix.B22 + Q22_global*temp;
		abdmatrix.B26 = abdmatrix.B26 + Q26_global*temp;
		abdmatrix.B66 = abdmatrix.B66 + Q66_global*temp;

		temp = plythichness*(plythichness*plythichness/4.0 + 3.0*LayerP*LayerP)/3.0;
		abdmatrix.D11 = abdmatrix.D11 + Q11_global*temp;
		abdmatrix.D12 = abdmatrix.D12 + Q12_global*temp;
		abdmatrix.D16 = abdmatrix.D16 + Q16_global*temp;
		abdmatrix.D22 = abdmatrix.D22 + Q22_global*temp;
		abdmatrix.D26 = abdmatrix.D26 + Q26_global*temp;
		abdmatrix.D66 = abdmatrix.D66 + Q66_global*temp;

	}

	return(abdmatrix);

}

struct LaminaQ LocalLaminaQ(struct CompsiteMat ply)
{
	LaminaQ Q;
	double dominator;
	dominator = (ply.Ex - ply.MIUxy *ply.MIUxy*ply.Ey);
	Q.Q11 = (ply.Ex)*(ply.Ex)/dominator;
	Q.Q12 =  (ply.MIUxy)*(ply.Ex)* (ply.Ey)/dominator;
	Q.Q22 = (ply.Ex)* (ply.Ey)/dominator;
	Q.Q66 = ply.Gxy;
	Q.Q16 = 0.0;
	Q.Q26 = 0.0;
	return(Q);
}

struct CompsiteE Q2E(struct LaminaQ Q)
{
	struct CompsiteE E;
	E.Ex = (Q.Q11*Q.Q22-Q.Q12*Q.Q12)/Q.Q22;
	E.Ey = (Q.Q11*Q.Q22-Q.Q12*Q.Q12)/Q.Q11;
	E.Gxy = Q.Q66;
	E.MIUxy = Q.Q12/Q.Q22;
	return(E);
}
struct TwoDStrain StrainG2L(struct TwoDStrain StrainG, double angle)
{
	struct TwoDStrain StrainL;
	double C, S;
	C = cos(angle*PI/180.);
	S = sin(angle*PI/180.);

	StrainL.EpsiX =   C*C*StrainG.EpsiX +   S*S*StrainG.EpsiY +       C*S*StrainG.EpsiXY;
	StrainL.EpsiY =   S*S*StrainG.EpsiX +   C*C*StrainG.EpsiY -       C*S*StrainG.EpsiXY;
	StrainL.EpsiXY   =-2*C*S*StrainG.EpsiX + 2*C*S*StrainG.EpsiY + (C*C-S*S)*StrainG.EpsiXY;

	return(StrainL);
}

struct TwoDStress Strain2Stress(struct LaminaQ Q, struct TwoDStrain strain)
{
	struct TwoDStress stress;
	stress.SigmaX = strain.EpsiX *Q.Q11 +   strain.EpsiY*Q.Q12;
	stress.SigmaY = strain.EpsiY*Q.Q22 +  strain.EpsiX*Q.Q12;
	stress.SigmaXY= strain.EpsiXY  *Q.Q66;

	return(stress);
}
struct LaminaQ LaminaQL2G(struct LaminaQ QL, double angle)
{
	double C, S;
	struct LaminaQ QG;
	C = cos(angle*PI/180.);
	S = sin(angle*PI/180.);
	QG.Q11 = QL.Q11 * C*C*C*C  +  (QL.Q12 + 2.0*QL.Q66) *  2.0*C*C*S*S +  QL.Q22 * S*S*S*S;
	QG.Q12 = QL.Q12 * ( C*C*C*C + S*S*S*S) + (QL.Q11 + QL.Q22 - 4*QL.Q66 ) * C*C*S*S;
	QG.Q16 = (QL.Q11 -QL.Q12- 2*QL.Q66)*C*C*C*S - (QL.Q22 -QL.Q12- 2*QL.Q66)*C*S*S*S;
	QG.Q22 = QL.Q11 * S*S*S*S + (QL.Q12 + 2*QL.Q66) * 2*C*C*S*S + QL.Q22 * C*C*C*C;
	QG.Q26 = (QL.Q11 -QL.Q12- 2*QL.Q66)*C*S*S*S - (QL.Q22 -QL.Q12- 2*QL.Q66)*C*C*C*S;
	QG.Q66 = (QL.Q11 + QL.Q22 - 2*QL.Q12 - 2*QL.Q66)*C*C*S*S + QL.Q66*(C*C*C*C + S*S*S*S);

	return(QG);

}

struct LaminaQ GlobalLaminaQ(struct CompsiteMat ply, double angle)
{
	double Q11_local, Q12_local, Q22_local, Q66_local, C, S, dominator;
	LaminaQ Q;
	dominator = (ply.Ex - ply.MIUxy *ply.MIUxy*ply.Ey);
	C = cos(angle*PI/180.);
	S = sin(angle*PI/180.);


		Q11_local = (ply.Ex)*(ply.Ex)/dominator;
		Q12_local =  (ply.MIUxy)*(ply.Ex)* (ply.Ey)/dominator;
		Q22_local = (ply.Ex)* (ply.Ey)/dominator;
		Q66_local = ply.Gxy;

		Q.Q11 = Q11_local * C*C*C*C  +  (Q12_local + 2.0*Q66_local) *  2.0*C*C*S*S +  Q22_local * S*S*S*S;
		Q.Q12 = Q12_local * ( C*C*C*C + S*S*S*S) + (Q11_local + Q22_local - 4*Q66_local ) * C*C*S*S;
		Q.Q16 = (Q11_local -Q12_local- 2*Q66_local)*C*C*C*S - (Q22_local -Q12_local- 2*Q66_local)*C*S*S*S;
		Q.Q22 = Q11_local * S*S*S*S + (Q12_local + 2*Q66_local) * 2*C*C*S*S + Q22_local * C*C*C*C;
		Q.Q26 = (Q11_local -Q12_local- 2*Q66_local)*C*S*S*S - (Q22_local -Q12_local- 2*Q66_local)*C*C*C*S;
		Q.Q66 = (Q11_local + Q22_local - 2*Q12_local - 2*Q66_local)*C*C*S*S + Q66_local*(C*C*C*C + S*S*S*S);

		return(Q);
}

struct CompsiteE PlateStiff(struct Laminate plate, double angle)
{
	struct LaminaQ Q, QL;
	struct CompsiteE Eangle;
	double t;
	t = PlateThick(plate);

	Q.Q11 = plate.stack.abdm.A11/t;
	Q.Q12 = plate.stack.abdm.A12/t;
	Q.Q16 = plate.stack.abdm.A16/t;
	Q.Q22 = plate.stack.abdm.A22/t;
	Q.Q26 = plate.stack.abdm.A26/t;
	Q.Q66 = plate.stack.abdm.A66/t;

	QL = LaminaQL2G(Q, -angle);

	Eangle = Q2E(QL);

	return(Eangle);

}


double ElastMod(struct Laminate plate)
{
	double em;

		em = (plate.stack.abdm.A11-plate.stack.abdm.A12*plate.stack.abdm.A12/plate.stack.abdm.A22)/ PlateThick(plate);


	return(em);
}

double StackThick(struct Stack stack)
{
	double t;
	int i;

	t=0.0;

	for (i=0;i<stack.NumLayers; i++)
		t =t + stack.ply[i].thick;
	return(t);
}


double PlateThick(struct Laminate plate)
{
	return(StackThick(plate.stack));
}

double SecArea(struct Laminate plate)
{
	return(StackThick(plate.stack)*plate.b);
}



double LayerHeight(struct Stack stack, int lay_id)
{
	double position;
	int i;
	position = -0.5*StackThick(stack);
	for (i=0;i<lay_id;i++)
	{
		position = position + stack.ply[i].thick;

	}

	position  = position + 0.5*stack.ply[lay_id].thick;

	return(position);
}






struct PlateDeform Disp(struct Laminate plate, struct PlateLoad pl)
{
	struct PlateDeform deform;
	/*send to EIGEN*/
	MatrixXd  abd(6,6);
	VectorXd loadv(6), strain;

	abd(0,0)=plate.stack.abdm.A11;
	abd(0,1)=plate.stack.abdm.A12;
	abd(0,2)=plate.stack.abdm.A16;
	abd(0,3)=plate.stack.abdm.B11;
	abd(0,4)=plate.stack.abdm.B12;
	abd(0,5)=plate.stack.abdm.B16;

	abd(1,0)=plate.stack.abdm.A12;
	abd(1,1)=plate.stack.abdm.A22;
	abd(1,2)=plate.stack.abdm.A26;
	abd(1,3)=plate.stack.abdm.B12;
	abd(1,4)=plate.stack.abdm.B22;
	abd(1,5)=plate.stack.abdm.B26;

	abd(2,0)=plate.stack.abdm.A16;
	abd(2,1)=plate.stack.abdm.A26;
	abd(2,2)=plate.stack.abdm.A66;
	abd(2,3)=plate.stack.abdm.B16;
	abd(2,4)=plate.stack.abdm.B26;
	abd(2,5)=plate.stack.abdm.B66;

	abd(3,0)=plate.stack.abdm.B11;
	abd(3,1)=plate.stack.abdm.B12;
	abd(3,2)=plate.stack.abdm.B16;
	abd(3,3)=plate.stack.abdm.D11;
	abd(3,4)=plate.stack.abdm.D12;
	abd(3,5)=plate.stack.abdm.D16;

	abd(4,0)=plate.stack.abdm.B12;
	abd(4,1)=plate.stack.abdm.B22;
	abd(4,2)=plate.stack.abdm.B26;
	abd(4,3)=plate.stack.abdm.D12;
	abd(4,4)=plate.stack.abdm.D22;
	abd(4,5)=plate.stack.abdm.D26;

	abd(5,0)=plate.stack.abdm.B16;
	abd(5,1)=plate.stack.abdm.B26;
	abd(5,2)=plate.stack.abdm.B66;
	abd(5,3)=plate.stack.abdm.D16;
	abd(5,4)=plate.stack.abdm.D26;
	abd(5,5)=plate.stack.abdm.D66;

	loadv(0) = pl.InPlane.Nxx;
	loadv(1) = pl.InPlane.Nyy;
	loadv(2) = pl.InPlane.Nxy;
	loadv(3) = pl.OutofPlane.Mxx;
	loadv(4) = pl.OutofPlane.Myy;
	loadv(5) = pl.OutofPlane.Mxy;

	//strain = abd.inverse()*loadv;
	strain = abd.llt().solve(loadv);
	
	/*recerved data from EIGEN*/
	deform.InPlane.EpsiX = strain(0);
	deform.InPlane.EpsiY = strain(1);
	deform.InPlane.EpsiXY = strain(2);
	deform.OutofPlane.KappaX = strain(3);
	deform.OutofPlane.KappaY = strain(4);
	deform.OutofPlane.KappaXY = strain(5);

	return(deform);

}










struct PlateLoad Load(struct Laminate plate, struct PlateDeform deform)
{
	struct PlateLoad pl;

	MatrixXd  abd(6,6);
	VectorXd  strain(6), loadv;

	abd(0,0)=plate.stack.abdm.A11;
	abd(0,1)=plate.stack.abdm.A12;
	abd(0,2)=plate.stack.abdm.A16;
	abd(0,3)=plate.stack.abdm.B11;
	abd(0,4)=plate.stack.abdm.B12;
	abd(0,5)=plate.stack.abdm.B16;

	abd(1,0)=plate.stack.abdm.A12;
	abd(1,1)=plate.stack.abdm.A22;
	abd(1,2)=plate.stack.abdm.A26;
	abd(1,3)=plate.stack.abdm.B12;
	abd(1,4)=plate.stack.abdm.B22;
	abd(1,5)=plate.stack.abdm.B26;

	abd(2,0)=plate.stack.abdm.A16;
	abd(2,1)=plate.stack.abdm.A26;
	abd(2,2)=plate.stack.abdm.A66;
	abd(2,3)=plate.stack.abdm.B16;
	abd(2,4)=plate.stack.abdm.B26;
	abd(2,5)=plate.stack.abdm.B66;

	abd(3,0)=plate.stack.abdm.B11;
	abd(3,1)=plate.stack.abdm.B12;
	abd(3,2)=plate.stack.abdm.B16;
	abd(3,3)=plate.stack.abdm.D11;
	abd(3,4)=plate.stack.abdm.D12;
	abd(3,5)=plate.stack.abdm.D16;

	abd(4,0)=plate.stack.abdm.B12;
	abd(4,1)=plate.stack.abdm.B22;
	abd(4,2)=plate.stack.abdm.B26;
	abd(4,3)=plate.stack.abdm.D12;
	abd(4,4)=plate.stack.abdm.D22;
	abd(4,5)=plate.stack.abdm.D26;

	abd(5,0)=plate.stack.abdm.B16;
	abd(5,1)=plate.stack.abdm.B26;
	abd(5,2)=plate.stack.abdm.B66;
	abd(5,3)=plate.stack.abdm.D16;
	abd(5,4)=plate.stack.abdm.D26;
	abd(5,5)=plate.stack.abdm.D66;

	strain(0) = deform.InPlane.EpsiX ;
	strain(1) = deform.InPlane.EpsiY ;
	strain(2) = deform.InPlane.EpsiXY;
	strain(3) = deform.OutofPlane.KappaX;
	strain(4) = deform.OutofPlane.KappaY;
	strain(5) = deform.OutofPlane.KappaXY;

	loadv = abd*strain;

	pl.InPlane.Nxx = loadv(0);
	pl.InPlane.Nyy = loadv(1);
	pl.InPlane.Nxy = loadv(2);
	pl.OutofPlane.Mxx = loadv(3);
	pl.OutofPlane.Myy = loadv(4);
	pl.OutofPlane.Mxy = loadv(5);

	return(pl);
}

struct ABDM ComplianceABD(struct ABDM Stiffness)
{
	MatrixXd  abd(6,6), abd_inv(6,6);
	struct ABDM Compliance;


	abd(0,0)=Stiffness.A11;
	abd(0,1)=Stiffness.A12;
	abd(0,2)=Stiffness.A16;
	abd(0,3)=Stiffness.B11;
	abd(0,4)=Stiffness.B12;
	abd(0,5)=Stiffness.B16;

	abd(1,0)=Stiffness.A12;
	abd(1,1)=Stiffness.A22;
	abd(1,2)=Stiffness.A26;
	abd(1,3)=Stiffness.B12;
	abd(1,4)=Stiffness.B22;
	abd(1,5)=Stiffness.B26;

	abd(2,0)=Stiffness.A16;
	abd(2,1)=Stiffness.A26;
	abd(2,2)=Stiffness.A66;
	abd(2,3)=Stiffness.B16;
	abd(2,4)=Stiffness.B26;
	abd(2,5)=Stiffness.B66;

	abd(3,0)=Stiffness.B11;
	abd(3,1)=Stiffness.B12;
	abd(3,2)=Stiffness.B16;
	abd(3,3)=Stiffness.D11;
	abd(3,4)=Stiffness.D12;
	abd(3,5)=Stiffness.D16;

	abd(4,0)=Stiffness.B12;
	abd(4,1)=Stiffness.B22;
	abd(4,2)=Stiffness.B26;
	abd(4,3)=Stiffness.D12;
	abd(4,4)=Stiffness.D22;
	abd(4,5)=Stiffness.D26;

	abd(5,0)=Stiffness.B16;
	abd(5,1)=Stiffness.B26;
	abd(5,2)=Stiffness.B66;
	abd(5,3)=Stiffness.D16;
	abd(5,4)=Stiffness.D26;
	abd(5,5)=Stiffness.D66;

	abd_inv = abd.inverse();

	Compliance.A11 = abd_inv(0,0);
	Compliance.A12 = abd_inv(0,1);
	Compliance.A16 = abd_inv(0,2);
	Compliance.B11 = abd_inv(0,3);
	Compliance.B12 = abd_inv(0,4);
	Compliance.B16 = abd_inv(0,5);

	Compliance.A12 = abd_inv(1,0);
	Compliance.A22 = abd_inv(1,1);
	Compliance.A26 = abd_inv(1,2);
	Compliance.B12 = abd_inv(1,3);
	Compliance.B22 = abd_inv(1,4);
	Compliance.B26 = abd_inv(1,5);

	Compliance.A16 = abd_inv(2,0);
	Compliance.A26 = abd_inv(2,1);
	Compliance.A66 = abd_inv(2,2);
	Compliance.B16 = abd_inv(2,3);
	Compliance.B26 = abd_inv(2,4);
	Compliance.B66 = abd_inv(2,5);

	Compliance.B11 = abd_inv(3,0);
	Compliance.B12 = abd_inv(3,1);
	Compliance.B16 = abd_inv(3,2);
	Compliance.D11 = abd_inv(3,3);
	Compliance.D12 = abd_inv(3,4);
	Compliance.D16 = abd_inv(3,5);

	Compliance.B12 = abd_inv(4,0);
	Compliance.B22 = abd_inv(4,1);
	Compliance.B26 = abd_inv(4,2);
	Compliance.D12 = abd_inv(4,3);
	Compliance.D22 = abd_inv(4,4);
	Compliance.D26 = abd_inv(4,5);

	Compliance.B16 = abd_inv(5,0);
	Compliance.B26 = abd_inv(5,1);
	Compliance.B66 = abd_inv(5,2);
	Compliance.D16 = abd_inv(5,3);
	Compliance.D26 = abd_inv(5,4);
	Compliance.D66 = abd_inv(5,5);


	return(Compliance);


}



struct TwoDStress PrincipalStress(struct TwoDStress stress)
{
	struct TwoDStress p_stress;
	double sig1, sig2;
	sig1 = 0.5*(stress.SigmaX + stress.SigmaY);
	sig2 = sqrt(0.25*(stress.SigmaX - stress.SigmaY)*(stress.SigmaX - stress.SigmaY) + stress.SigmaXY * stress.SigmaXY);
	p_stress.SigmaX = sig1 + sig2;
	p_stress.SigmaY = sig1 - sig2;
	p_stress.SigmaXY = 0.;
	return(p_stress);
}


struct TwoDStrain PrincipalStrain(struct TwoDStrain strain)
{
	struct TwoDStrain p_strain;
	double epsi1, epsi2;

	epsi1 = 0.5*(strain.EpsiX + strain.EpsiY);
	epsi2 = 0.5*sqrt((strain.EpsiX - strain.EpsiY)*(strain.EpsiX - strain.EpsiY)+ strain.EpsiXY*strain.EpsiXY);

	p_strain.EpsiX = epsi1 + epsi2;
	p_strain.EpsiY = epsi1 - epsi2;
	p_strain.EpsiXY = 0.;
	return(p_strain);

}

struct TwoDStress StressTransfer(struct TwoDStress stress, double angle)
{
	struct TwoDStress stress_local;
	double C, S;
	C = cos(angle*PI/180.);
	S = sin(angle*PI/180.);
	
	stress_local.SigmaX =	C*C*stress.SigmaX  +   S*S*stress.SigmaY    +    2.*C*S*stress.SigmaXY;
	stress_local.SigmaY =	S*S*stress.SigmaX  +   C*C*stress.SigmaY    -    2.*C*S*stress.SigmaXY;
	stress_local.SigmaXY = -C*S*stress.SigmaX  +   C*S*stress.SigmaY    + (C*C-S*S)*stress.SigmaXY;

	return(stress_local);
}

struct TwoDStrain StrainTransfer(struct TwoDStrain strain, double angle)
{
	struct TwoDStrain StrainL;
	double C,S;
	C = cos(angle*PI/180.);
	S = sin(angle*PI/180.);

	StrainL.EpsiX =   C*C*strain.EpsiX +   S*S*strain.EpsiY  +       C*S*strain.EpsiXY;
	StrainL.EpsiY =   S*S*strain.EpsiX +   C*C*strain.EpsiY  -       C*S*strain.EpsiXY;
	StrainL.EpsiXY   =-2*C*S*strain.EpsiX + 2*C*S*strain.EpsiY + (C*C-S*S)*strain.EpsiXY;

	return(StrainL);
}


