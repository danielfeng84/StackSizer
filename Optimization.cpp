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
#include "Optimization.h"
#include "PlyCheck.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <omp.h>
#include <Eigen/Dense>
using namespace Eigen;

int CompareID (const void * a, const void * b)
{
	 ud *p1 = ( ud *)a;

	 ud * p2 = (ud*)b;

	 if(p1->order < p2->order)
		 return -1;
	 else
	 {
		 if(p1->order == p2->order)
			return 0;
		 else
			 return 1;
	 }
}

long long powint (int a,int b)
{

	int i;
	long long c, a_long;
	a_long = (long long)a;

	c=1;
	for (i=1;i<=b;i++)
		c=c*a_long;
	return(c);
}


int *MaxStiffEstim(int num_ply_web)
{
	static int num[3];

	num[0] = MAX(0.5*num_ply_web,3);
	num[1] = MAX(2*(int)(0.2*num_ply_web),2);
	num[2] = MAX(ceil(0.1*num_ply_web),1);

	return(num);
}





inline double RandN(double std)
{
	double U1,U2,X1;
	U1 = (double)rand()/(double)RAND_MAX;
	U2 = (double)rand()/(double)RAND_MAX;
	X1 = std * sqrt(-2*log(U1)) * cos(2*PI*U2);
	return(X1);

}


int Binomial2(double Pr)
{
    double coin;
    coin = (double) rand();
    if (coin<= (double) RAND_MAX*(1.-Pr)*(1.-Pr))
        return(0);
    else
	{
		if  (coin<=(double) RAND_MAX*(1.-Pr*Pr))
			return(1);
		else
			return(2);
	}
}


inline int RandU(int num)
{
	//srand((int)time(NULL));
	return(rand()%num);
}


inline int RandWander(int num)
{
	//srand((int)time(NULL));
	return(rand()%(2*num)-num);
}





void EsimatePercentFirst(struct PlateCopula* plate)
{
	int i, j,  num_comb, list[100][4];

	double Candidate, PrincipStrain;
	struct Laminate  plate_test;
	struct PlateDeform deform;
	struct PlateLoad pl;
	pl = plate->pl;

	plate_test = plate->plate;
	Candidate =-1.;


	if (plate_test.stack.NumLayers==12)//the only alaowable stacking sequence for 12 plies is [45/-45/90/0/45/-45]s
	{
		plate->dirc[0] = 2;
		plate->dirc[1] = 4;
		plate->dirc[2] = 4;
		plate->dirc[3] = 2;


	}

	else
	{
		num_comb = PlyBundleList(list, plate_test.stack.NumLayers);
		
		for (i=0;i<num_comb;i++)
		{
			for (j=0;j<list[i][0];j++)
				plate_test.stack.PlyDirection[j] = 0;
			for (j=list[i][0];j<list[i][0]+list[i][1];j++)
				plate_test.stack.PlyDirection[j] = 45;
			for (j=list[i][0]+list[i][1];j<list[i][0]+list[i][1]+list[i][2];j++)
				plate_test.stack.PlyDirection[j] = -45;
			for (j=list[i][0]+list[i][1]+list[i][2];j<list[i][0]+list[i][1]+list[i][2]+list[i][3];j++)
				plate_test.stack.PlyDirection[j] = 90;


			plate_test.stack.abdm = ABD(plate_test.stack);

			deform = Disp(plate_test,pl);


			PrincipStrain = (deform.InPlane.EpsiX + deform.InPlane.EpsiY - sqrt((deform.InPlane.EpsiX - deform.InPlane.EpsiY)*(deform.InPlane.EpsiX - deform.InPlane.EpsiY)+deform.InPlane.EpsiXY*deform.InPlane.EpsiXY))/2.;


			if (PrincipStrain>Candidate)
			{
				Candidate = PrincipStrain;
				for (j=0;j<4;j++)
					plate->dirc[j]=list[i][j];
			}


		}

	}


	

}






void EsimatePercent(struct PlateCopula *plates, int id)
{
	int layer_num_new, layer_num_old, PercentInc[30];
	int i, j, k, flag, total_num,  num_fail, dist[4];
	struct Laminate  plate_test;
	double frac[4], PrincipStrain, Candidate;
	struct PlateDeform deform;
	struct PlateLoad pl;
	pl = plates[id].pl;

	Candidate = -1.;

	layer_num_new = plates[id].plate.stack.NumLayers;
	layer_num_old = plates[id-1].plate.stack.NumLayers;


	if (layer_num_old == layer_num_new)
	{
		for (j=0;j<layer_num_new;j++)
			plates[id].dirc[j] = plates[id-1].dirc[j];
		return;
	}

	plate_test = plates[id].plate;


	for (j=0;j<plates[id-1].dirc[0];j++)
		plate_test.stack.PlyDirection[j] = 0;
	for (j=plates[id-1].dirc[0];j<plates[id-1].dirc[0]+plates[id-1].dirc[1];j++)
		plate_test.stack.PlyDirection[j] = 45;
	for (j=plates[id-1].dirc[0]+plates[id-1].dirc[1];j<plates[id-1].dirc[0]+plates[id-1].dirc[1]+plates[id-1].dirc[2];j++)
		plate_test.stack.PlyDirection[j] = -45;
	for (j=plates[id-1].dirc[0]+plates[id-1].dirc[1]+plates[id-1].dirc[2];j<plates[id-1].dirc[0]+plates[id-1].dirc[1]+plates[id-1].dirc[2]+plates[id-1].dirc[3];j++)
		plate_test.stack.PlyDirection[j] = 90;



	
	
	total_num = powint (4, (layer_num_new-layer_num_old)/4);

	for (j=0;j<4;j++)
		plates[id].dirc[j] = plates[id-1].dirc[j];



	for(j=0;j<(layer_num_new-layer_num_old)/4;j++)
	{
		num_fail =0;
		
		
		for (flag=0;flag<4;flag++)
		{
			for (i=0;i<4;i++)
			{
				PercentInc[i] = plates[id].dirc[i];
				dist[i]=0; 
			}

			switch (flag)
			{
			case 0:// add 0,0,0,0
				{
					plate_test.stack.PlyDirection[layer_num_old + 2*j] = 0;
					plate_test.stack.PlyDirection[layer_num_new-1-layer_num_old-2*j] = 0;
					
					
					plate_test.stack.PlyDirection[layer_num_old + 2*j + 1] = 0;
					plate_test.stack.PlyDirection[layer_num_new-2-layer_num_old-2*j] = 0;
					
					
					PercentInc[0] = PercentInc[0]+4;

					
					for (i=0;i<4;i++)
					{
						dist[0] = dist[0] + (PercentInc[i]-plate_test.stack.NumLayers/4)*(PercentInc[i]-plate_test.stack.NumLayers/4);
					}
					
					
						
					break;
				}

			case 1://add 90,90,90,90
				{
					plate_test.stack.PlyDirection[layer_num_old + 2*j] = 90;
					plate_test.stack.PlyDirection[layer_num_new-1-layer_num_old - 2*j] = 90;
					
					
					
					plate_test.stack.PlyDirection[layer_num_old + 2*j + 1] = 90;
					plate_test.stack.PlyDirection[layer_num_new-2-layer_num_old - 2*j] = 90;
					
					
					
					
					PercentInc[3] = PercentInc[3] + 4;

					for (i=0;i<4;i++)
					{
						dist[1] = dist[1] + (PercentInc[i]-plate_test.stack.NumLayers/4)*(PercentInc[i]-plate_test.stack.NumLayers/4);
					}
					break;
				}
			
			case 2:
				{
					plate_test.stack.PlyDirection[layer_num_old + 2*j] = 90;
					plate_test.stack.PlyDirection[layer_num_new-1-layer_num_old-2*j] = 90;
					
					
					plate_test.stack.PlyDirection[layer_num_old + 2*j + 1] = 0;
					plate_test.stack.PlyDirection[layer_num_new-2-layer_num_old- 2*j] = 0;
					
					
					
					PercentInc[3] = PercentInc[3] +2;
					
					PercentInc[0] = PercentInc[0] + 2;
					
					
					for (i=0;i<4;i++)
					{
						dist[2] = dist[2] + (PercentInc[i]-plate_test.stack.NumLayers/4)*(PercentInc[i]-plate_test.stack.NumLayers/4);
					}
					
					
					break;
				}

			default:
				{
					plate_test.stack.PlyDirection[layer_num_old + 2*j] = 45;
					plate_test.stack.PlyDirection[layer_num_new-1-layer_num_old- 2*j] = 45;
					
					
					plate_test.stack.PlyDirection[layer_num_old + 2*j + 1] = -45;
					plate_test.stack.PlyDirection[layer_num_new-2-layer_num_old- 2*j] = -45;
					
					
					
					
					PercentInc[1] = PercentInc[1] + 2;
					PercentInc[2] = PercentInc[2] + 2;

					for (i=0;i<4;i++)
					{
						dist[3] = dist[3] + (PercentInc[i]-plate_test.stack.NumLayers/4)*(PercentInc[i]-plate_test.stack.NumLayers/4);
					}

					break;
				}
			}
			
			for (i=0;i<4;i++)
				frac[i] = (double) PercentInc[i]/layer_num_new;
			if ((((frac[0]>0.45)||(frac[0]<0.15))||(((frac[3]>0.45)||(frac[3]<0.15))))||(((frac[1]>0.65)||(frac[1]<0.15))))
			{
				num_fail++;
				continue;
			}
			
			plate_test.stack.abdm = ABD(plate_test.stack);
			deform = Disp(plate_test, pl);
			
			PrincipStrain = (deform.InPlane.EpsiX + deform.InPlane.EpsiY - sqrt(( deform.InPlane.EpsiX - deform.InPlane.EpsiY)*(deform.InPlane.EpsiX - deform.InPlane.EpsiY)+deform.InPlane.EpsiXY*deform.InPlane.EpsiXY))/2.;

			if (PrincipStrain>Candidate)
			{
				Candidate = PrincipStrain;
				k=flag;
			}
	}


	if (num_fail==4)
	{

		k=0;

		for (i=1;i<4;i++)
		{
			if (dist[i]<dist[k])
			{
				k=i;
			}
		}

	}

	switch (k)
		{
		case 0: plates[id].dirc[0] = plates[id].dirc[0]+4; break;
		case 1: plates[id].dirc[3] = plates[id].dirc[3] + 4;break;
		case 2:  plates[id].dirc[3] = plates[id].dirc[3] +2; plates[id].dirc[0] =plates[id].dirc[0] + 2;break;
		default:	plates[id].dirc[1] = plates[id].dirc[1] + 2;plates[id].dirc[2] =plates[id].dirc[2] + 2;break;
		}







}



if (plates[id].plate.stack.NumLayers!= (plates[id].dirc[0] + plates[id].dirc[1] + plates[id].dirc[2] +plates[id].dirc[3]))
{
	printf("There is no solution satisfying the percentage requirement, id:%d,  thick:%d \n",id, plates[id].plate.stack.NumLayers);

	printf("check\n");
}



}





void InitialPlyDeploy(struct PlateCopula plates[], int num_plates)
{
	int i, j, pos, ii, flag, p0, p45, p_45, p90, act;

	int index[10000];


	//ALL deactivated
	for (i=0;i<num_plates;i++)
	{
		for (j=0;j<MAXLAYEROFPLY;j++)
		{
			plates[i].posi0[j][0]=MAXLAYEROFPLY+1;
			plates[i].posi0[j][1]=0;
		}
		for (j=0;j<MAXLAYEROFPLY;j++)
		{
			plates[i].posi45[j][0]=MAXLAYEROFPLY+1;
			plates[i].posi45[j][1]=0;
		}
		for (j=0;j<MAXLAYEROFPLY;j++)
		{
			plates[i].posi_45[j][0]=MAXLAYEROFPLY+1;
			plates[i].posi_45[j][1]=0;
		}
		for (j=0;j<MAXLAYEROFPLY;j++)
		{
			plates[i].posi90[j][0]=MAXLAYEROFPLY+1;
			plates[i].posi90[j][1]=0;
		}
	}



	plates[num_plates-1].posi45[0][0]=0;
	plates[num_plates-1].posi45[0][1]=2;
	plates[num_plates-1].posi45[0][3]=45;


	plates[num_plates-1].posi_45[0][0]=1;
	plates[num_plates-1].posi_45[0][1]=2;
	plates[num_plates-1].posi_45[0][3]=-45;

	plates[num_plates-1].posi0[0][0]=2;
	plates[num_plates-1].posi0[0][1]=2;
	plates[num_plates-1].posi0[0][3]=0;

	plates[num_plates-1].posi90[0][0]=3;
	plates[num_plates-1].posi90[0][1]=2;
	plates[num_plates-1].posi90[0][3]=90;


	plates[num_plates-1].posi45[1][0]=4;
	plates[num_plates-1].posi45[1][1]=2;
	plates[num_plates-1].posi45[1][3]=45;

	plates[num_plates-1].posi_45[1][0]=5;
	plates[num_plates-1].posi_45[1][1]=2;
	plates[num_plates-1].posi_45[1][3]=-45;

	
	ii=0;
	
	p0=1;
	p45=2;
	p_45=2;
	p90=1;




	for (pos=6;pos<plates[num_plates-1].plate.stack.NumLayers/2;pos++)
	{
		act=0;

		do
		{
			flag = ii%4;

			if (flag==0)
			{
				if (p0<plates[num_plates-1].half_nums[0])
				{
					plates[num_plates-1].posi0[p0][0]=pos;
					plates[num_plates-1].posi0[p0][1]=1;
					plates[num_plates-1].posi0[p0][3]=0;
					
					ii++;
					p0++;
					break;
				}
				else
					ii++;
			}

			
			if (flag==1)
			{
				if (p45<plates[num_plates-1].half_nums[1])
				{
					plates[num_plates-1].posi45[p45][0]=pos;
					plates[num_plates-1].posi45[p45][1]=1;
					plates[num_plates-1].posi45[p45][3]=45;
					
					ii++;
					p45++;
					break;
				}
				else
					ii++;
			}


			if (flag==2)
			{
				if (p_45<plates[num_plates-1].half_nums[2])
				{
					plates[num_plates-1].posi_45[p_45][0]=pos;
					plates[num_plates-1].posi_45[p_45][1]=1;
					plates[num_plates-1].posi_45[p_45][3]=-45;
					
					ii++;
					p_45++;
					break;
				}
				else
					ii++;
			}


			if (flag==3)
			{
				if (p90<plates[num_plates-1].half_nums[3])
				{
					plates[num_plates-1].posi90[p90][0]=pos;
					plates[num_plates-1].posi90[p90][1]=1;
					plates[num_plates-1].posi90[p90][3]=90;
					
					ii++;
					p90++;
					break;
				}
				else
					ii++;
			}




		}while(act==0);
		


	
	}


	if (p0!= plates[num_plates-1].half_nums[0])
		printf("Error.....................\n");
	if (p45!= plates[num_plates-1].half_nums[1])
		printf("Error.....................\n");
	if (p_45!= plates[num_plates-1].half_nums[2])
		printf("Error.....................\n");
	if (p90!= plates[num_plates-1].half_nums[3])
		printf("Error.....................\n");


	initialize_and_permute(index, plates[num_plates-1].plate.stack.NumLayers/2-2);

	for (i=num_plates-2;i>=0;i--)
	{

		for (j=0;j<plates[num_plates-1].half_nums[0];j++)
		{
			plates[i].posi0[j][0] = plates[i+1].posi0[j][0];
			plates[i].posi0[j][3] = 0;

			if (j<plates[i].half_nums[0])
				plates[i].posi0[j][1] = 1;
			else
				plates[i].posi0[j][1] = 0;
		}


		for (j=0;j<plates[num_plates-1].half_nums[1];j++)
		{

			plates[i].posi45[j][0] = plates[i+1].posi45[j][0];
			plates[i].posi45[j][3] = 45;

			if (j<plates[i].half_nums[1])
				plates[i].posi45[j][1] = 1;
			else
				plates[i].posi45[j][1] = 0;

		}



		for (j=0;j<plates[num_plates-1].half_nums[2];j++)
		{
			plates[i].posi_45[j][0] = plates[i+1].posi_45[j][0];
			plates[i].posi_45[j][3] = -45;

			if (j<plates[i].half_nums[2])
				plates[i].posi_45[j][1] = 1;
			else
				plates[i].posi_45[j][1] = 0;

		}


		for (j=0;j<plates[num_plates-1].half_nums[3];j++)
		{
			plates[i].posi90[j][0] = plates[i+1].posi90[j][0];
			plates[i].posi90[j][3] = 90;

			if (j<plates[i].half_nums[3])
				plates[i].posi90[j][1] = 1;
			else
				plates[i].posi90[j][1] = 0;
		}

	}


}








ud Find(const struct PlateCopula* plate, int layID)
{
	int i;
	ud ply;


	for (i=0;i<plate->half_nums[0];i++)
	{
		if ((plate->posi0[i][0])==layID)
		{
			ply.angle = 0;
			ply.order = i;
			return(ply);
		}
	}


	for (i=0;i<plate->half_nums[1];i++)
	{
		if ((plate->posi45[i][0])==layID)
		{
			ply.angle = 45;
			ply.order = i;
			return(ply);
		}
	}


	for (i=0;i<plate->half_nums[2];i++)
	{
		if ((plate->posi_45[i][0])==layID)
		{
			ply.angle = -45;
			ply.order = i;
			return(ply);
		}
	}


	for (i=0;i<plate->half_nums[3];i++)
	{
		if ((plate->posi90[i][0])==layID)
		{
			ply.angle = 90;
			ply.order = i;
			return(ply);
		}
	}




	printf("error to find a play");
	exit(0);

}


struct Laminate LamGen(struct  PlateCopula plate)
{
	struct Laminate laminate;
	int i, count, lay, num_layer;
	ud record[MAXLAYEROFPLY];
	count=0;

	laminate = plate.plate;
	num_layer = laminate.stack.NumLayers;

	for (i=0;i<MAXLAYEROFPLY;i++)
		record[i].order = 1e6;
	for (i=0;i<plate.half_nums[0];i++)
	{
		lay = plate.posi0[i][0];
		if (plate.posi0[i][1]!=0)
		{
			record[count].order = lay;
			record[count].angle = plate.posi0[i][3];
			//record[count].angle = 0;
			//if (record[count].angle!=0)
			//	printf("incorrect 0 \n");
			if(((record[count].angle!=45)&&(record[count].angle!=-45))&&((record[count].angle!=0)&&(record[count].angle!=90)))
				printf("ply angle incorrect\n");

			count++;
		}
	}

	for (i=0;i<plate.half_nums[1];i++)
	{
		lay = plate.posi45[i][0];
		if (plate.posi45[i][1]!=0)
		{
			record[count].order = lay;
			record[count].angle = plate.posi45[i][3];
			//record[count].angle = 45;
			//if (record[count].angle!=45)
			//	printf("incorrect 45 \n");

			if(((record[count].angle!=45)&&(record[count].angle!=-45))&&((record[count].angle!=0)&&(record[count].angle!=90)))
				printf("ply angle incorrect\n");
			count++;
		}
	}

	for (i=0;i<plate.half_nums[2];i++)
	{
		lay = plate.posi_45[i][0];
		if (plate.posi_45[i][1]!=0)
		{
			record[count].order = lay;
			record[count].angle = plate.posi_45[i][3];
			//if (record[count].angle!=-45)
			//	printf("incorrect -45 %d\n",record[count].angle);
			//record[count].angle = -45;
			if(((record[count].angle!=45)&&(record[count].angle!=-45))&&((record[count].angle!=0)&&(record[count].angle!=90)))
				printf("ply angle incorrect\n");
			count++;
		}
	}

	for (i=0;i<plate.half_nums[3];i++)
	{
		lay = plate.posi90[i][0];
		if (plate.posi90[i][1]!=0)
		{
			record[count].order = lay;
			record[count].angle = plate.posi90[i][3];
			
			if(((record[count].angle!=45)&&(record[count].angle!=-45))&&((record[count].angle!=0)&&(record[count].angle!=90)))
				printf("ply angle incorrect\n");
			count++;
		}
	}

	qsort(record, num_layer, sizeof (ud), CompareID);


	for (i=0;i<num_layer/2;i++)
	{
		laminate.stack.PlyDirection[i] = record[i].angle;

		laminate.stack.PlyDirection[num_layer-1-i] =laminate.stack.PlyDirection[i];
	}

	if ((plate.half_nums[1]!=plate.half_nums[2])&&(laminate.stack.PlyDirection[num_layer/2]==45))
		laminate.stack.PlyDirection[num_layer/2]=-45;

	return(laminate);



}





void RestorePly(struct PlateCopula plates[], int num_plates, int lay1)
{
	ud ply1;
	int max_pos , i, j;


	max_pos =  plates[num_plates-1].plate.stack.NumLayers/2-1;
	ply1 = Find(&plates[num_plates-1], max_pos);



	for (i=0;i<num_plates;i++)
	{
		for (j=0;j<MAXLAYEROFPLY;j++)
		{

			if (plates[i].posi0[j][0]>=lay1)
				plates[i].posi0[j][0]++;

			if (plates[i].posi90[j][0]>=lay1)
				plates[i].posi90[j][0]++;

			if (plates[i].posi45[j][0]>=lay1)
				plates[i].posi45[j][0]++;

			if (plates[i].posi_45[j][0]>=lay1)
				plates[i].posi_45[j][0]++;

		}

	}

	for (i=0;i<num_plates;i++)
	{
		switch(ply1.angle)
		{
		case 0: plates[i].posi0[ply1.order][0] = lay1; break;

		case 45: plates[i].posi45[ply1.order][0] = lay1; break;

		case -45: plates[i].posi_45[ply1.order][0] = lay1; break;

		case 90: plates[i].posi90[ply1.order][0] = lay1; break;

		}
	}

}





int SinkPly(struct PlateCopula plates[], int num_plates, int lay1)
{

	ud ply1;
	int i, j, max_pos;

	ply1 = Find(&plates[num_plates-1], lay1);
	max_pos =  plates[num_plates-1].plate.stack.NumLayers/2-1;

	for (i=0;i<num_plates;i++)
	{
		for (j=0;j<MAXLAYEROFPLY;j++)
		{
			if (plates[i].posi0[j][0]>lay1)
				plates[i].posi0[j][0]--;

			if (plates[i].posi90[j][0]>lay1)
				plates[i].posi90[j][0]--;

			if (plates[i].posi45[j][0]>lay1)
				plates[i].posi45[j][0]--;

			if (plates[i].posi_45[j][0]>lay1)
				plates[i].posi_45[j][0]--;
		}

	}

	for (i=0;i<num_plates;i++)
	{
		switch(ply1.angle)
		{
		case 0: plates[i].posi0[ply1.order][0] = max_pos; break;

		case 45: plates[i].posi45[ply1.order][0] = max_pos; break;

		case -45: plates[i].posi_45[ply1.order][0] = max_pos; break;

		case 90: plates[i].posi90[ply1.order][0] = max_pos; break;

		}
	}

#pragma omp parallel for
    for (int para_j=0; para_j<num_plates; para_j++)
	{
        plates[para_j].plate = LamGen(plates[para_j]);
	}

	return(1);

}




int Exchange2Plies(struct PlateCopula plates[], int num_plates, int lay1, int lay2)
{
	ud ply1, ply2;
	int i;

	ply1 = Find(&plates[num_plates-1], lay1);
	ply2 = Find(&plates[num_plates-1], lay2);



	switch(ply1.angle)
	{
	case 0:
		{
			for (i=0;i<num_plates;i++)
				plates[i].posi0[ply1.order][0] = lay2;
			break;
		}
	case 45:
		{
			for (i=0;i<num_plates;i++)
				plates[i].posi45[ply1.order][0] = lay2;
			break;
		}
	case -45:
		{
			for (i=0;i<num_plates;i++)
				plates[i].posi_45[ply1.order][0] = lay2;
			break;
		}
	case 90:
		{
			for (i=0;i<num_plates;i++)
				plates[i].posi90[ply1.order][0] = lay2;
			break;
		}
	}

	switch(ply2.angle)
	{
	case 0:
		{
			for (i=0;i<num_plates;i++)
				plates[i].posi0[ply2.order][0] = lay1;
			break;
		}
	case 45:
		{
			for (i=0;i<num_plates;i++)
				plates[i].posi45[ply2.order][0] = lay1;
			break;
		}
	case -45:
		{
			for (i=0;i<num_plates;i++)
				plates[i].posi_45[ply2.order][0] = lay1;
			break;
		}
	case 90:
		{
			for (i=0;i<num_plates;i++)
				plates[i].posi90[ply2.order][0] = lay1;
			break;
		}
	}


#pragma omp parallel for
    for (int para_j=0; para_j<num_plates; para_j++)
	{
        plates[para_j].plate = LamGen(plates[para_j]);
	}

	/*
	if ((abs(ply1.angle)==45)&&(abs(ply2.angle)!=45))
		printf("45 and -45 are leaving\n");

	if ((abs(ply1.angle)!=45)&&(abs(ply2.angle)==45))
		printf("45 and -45 are leaving\n");
*/

	return(1);
}







void Exchange2Plies0(struct PlateCopula plates[], int num_plates, int lay1, int lay2)
{
	ud ply1, ply2;
	int i;

	ply1 = Find(&plates[num_plates-1], lay1);
	ply2 = Find(&plates[num_plates-1], lay2);


	switch(ply1.angle)
	{
	case 0:
		{
			for (i=0;i<num_plates;i++)
				plates[i].posi0[ply1.order][0] = lay2;
			break;
		}
	case 45:
		{
			for (i=0;i<num_plates;i++)
				plates[i].posi45[ply1.order][0] = lay2;
			break;
		}
	case -45:
		{
			for (i=0;i<num_plates;i++)
				plates[i].posi_45[ply1.order][0] = lay2;
			break;
		}
	case 90:
		{
			for (i=0;i<num_plates;i++)
				plates[i].posi90[ply1.order][0] = lay2;
			break;
		}
	}

	switch(ply2.angle)
	{
	case 0:
		{
			for (i=0;i<num_plates;i++)
				plates[i].posi0[ply2.order][0] = lay1;
			break;
		}
	case 45:
		{
			for (i=0;i<num_plates;i++)
				plates[i].posi45[ply2.order][0] = lay1;
			break;
		}
	case -45:
		{
			for (i=0;i<num_plates;i++)
				plates[i].posi_45[ply2.order][0] = lay1;
			break;
		}
	case 90:
		{
			for (i=0;i<num_plates;i++)
				plates[i].posi90[ply2.order][0] = lay1;
			break;
		}
	}


#pragma omp parallel for
    for (int para_j=0; para_j<num_plates; para_j++)
	{
        plates[para_j].plate = LamGen(plates[para_j]);
	}




}


int Reversal(struct PlateCopula plates[], int num_plates, int *clusters, int num_clusters, int cluster_id)
{
	int i, j, opera;

	if (cluster_id>=num_clusters)
		return(0);
	if ((clusters[cluster_id]==45)||(clusters[cluster_id]==-45))
		return(0);

	opera=0;
	for (i=0;i<num_plates;i++)
	{
		for (j=1;j<MAXLAYEROFPLY;j++)
		{
			if ((plates[i].posi0[j][1]!=0)&&(plates[i].posi0[j][2]==cluster_id))
			{
				opera=1;
				plates[i].posi0[j][3]=abs(90-clusters[cluster_id]);
			//	printf("%d\n",plates[i].posi0[j][3] );
			}

			if ((plates[i].posi90[j][1]!=0)&&(plates[i].posi90[j][2]==cluster_id))
			{
				opera=1;

				plates[i].posi90[j][3]=abs(90-clusters[cluster_id]);
			}

		}

	}

	if (opera==1)
	{
		clusters[cluster_id] = abs(90-clusters[cluster_id]);
#pragma omp parallel for
    for (int para_j=0; para_j<num_plates; para_j++)
	{
        plates[para_j].plate = LamGen(plates[para_j]);
	}

	}

	return(opera);
}



struct MarginOfSafety MarginSafety(struct Laminate plate, struct PlateLoad pl, double t_allow, double c_allow, double PB_coef_c, double PB_coef_s, double Smear_Ratio)
{
	struct PlateDeform disp;
	struct TwoDStrain PStrain;
	struct PlateBuckleState bl_c, bl_s;
	double temp, smear;
	struct MarginOfSafety ms;

	disp = Disp(plate, pl);
	PStrain = PrincipalStrain(disp.InPlane);

	ms.subcase = pl.LoadID;
	
	if ((PStrain.EpsiX>0)&&(PStrain.EpsiY>0))
	{
		ms.ms = t_allow/PStrain.EpsiX-1;
		ms.mode = 1;
		return(ms);
	}

	else
	{
		bl_c = PlateBuckling_UC_F_SSSS_Ortho(plate);// two directional compression
		bl_s = PlateBuckling_S_F_SSSS_Ortho(plate);

		//PB_coef_c = 1.5;
		//PB_coef_s = 1.25;

		//smear =(double) (plate.stack.NumLayers*140+9*70)/(plate.stack.NumLayers*140);

		smear = 1.+ Smear_Ratio;
		smear = smear * smear *smear;

		if (pl.InPlane.Nxx<0)
			smear=1.;


		temp = ((pl.InPlane.Nxx/bl_c.CrBuckleLoad.InPlane.Nxx/smear*PB_coef_c) + (pl.InPlane.Nxy/bl_s.CrBuckleLoad.InPlane.Nxy*PB_coef_s)*(pl.InPlane.Nxy/bl_s.CrBuckleLoad.InPlane.Nxy*PB_coef_s));



		if (temp>0)
		{

			ms.ms = 1/temp-1; //buckling is allowed to occur at 1.2LL


			if (PStrain.EpsiX>0)
				ms.ms = MIN(ms.ms, fabs(t_allow/PStrain.EpsiX)-1);
			else
				ms.ms = MIN(ms.ms, fabs(c_allow/PStrain.EpsiX)-1);


			if (PStrain.EpsiY<0)
				ms.ms = MIN(ms.ms, fabs(c_allow/PStrain.EpsiY)-1);
			else
				ms.ms = MIN(ms.ms, fabs(t_allow/PStrain.EpsiY)-1);


			if (ms.ms==(1/temp-1))
				ms.mode = 2;
			else
				ms.mode = 1;


			return(ms);
		}

		else
		{
			//buckling will not occure
			ms.mode =3;
			ms.ms = MIN(fabs(t_allow/PStrain.EpsiX)-1, fabs(c_allow/PStrain.EpsiY)-1);
			return(ms);
		}


	}




}




int IsNeighbor(int num_plates,  int *ConnectVec, int id1, int id2)
{

	int pos;

	pos = id1 + id2*num_plates;

	if ((pos<0)||(pos>=num_plates*num_plates))
	{
		printf("error in computing index \n");
		//exit(0);
	}

	return(ConnectVec[pos]);



}








void PartitionCalc(int *clusters, int *num_clusters,  struct PlateCopula *plates, int num_plates, int max_layer, int *ConnectVec)
{
	int i, j, lay1, flag, isolate, ci, cj, count;
	ud ply1;
	struct Clusters cluster;



	count=0;

	for (lay1=0; lay1<max_layer/2; lay1++)
	{
		ply1 = Find(&plates[num_plates-1], lay1);

		switch(ply1.angle)
		{

		case 0:
			{
				cluster.num_clusters = 0;


				for (i=0;i<num_plates;i++)
				{
					if (plates[i].posi0[ply1.order][1]!=0)
					{
						cluster.cluster_size[cluster.num_clusters] =1;
						cluster.container[cluster.num_clusters] = (int *)malloc(sizeof(int));
						cluster.container[cluster.num_clusters][0] = i;
						cluster.num_clusters++;
					}

				}
				for (int temp=0; temp<5;temp++)
				{
						do
				{

					isolate = cluster.num_clusters;

					for (ci=0;ci<cluster.num_clusters;ci++)
					{
						for (cj=ci+1;cj<cluster.num_clusters;cj++)
						{

							for (i=0;i<cluster.cluster_size[ci]; i++)
							{

								for (j=0;j<cluster.cluster_size[cj]; j++)
								{

									flag = IsNeighbor(num_plates,  ConnectVec, plates[ cluster.container[ci][i] ].PlateID, plates[ cluster.container[cj][j] ].PlateID);
									if (flag==1)
									{
										 MergeTwoClusters(&cluster, ci, cj);
									//	printf("%d, %d, %d\n",cluster.num_clusters, ci, cj);
										break;
									}
								}
								if (flag==1)
									break;
							}
						}
					}
				} while((cluster.num_clusters<isolate)&&(cluster.num_clusters>1));

				}



				/*
				for (i=0;i<cluster.num_clusters;i++)
				{
					for (j=0;j<cluster.cluster_size[i];j++)
						plates[cluster.container[i][j]].posi0[ply1.order][2] = i;

				}
				*/

				//printf("%d clusters for layer %d\n", cluster.num_clusters, lay1);

				for (i=0;i<cluster.num_clusters;i++)
				{
					clusters[count] = 0;
						for (j=0;j<cluster.cluster_size[i];j++)
							plates[cluster.container[i][j]].posi0[ply1.order][2] = count;
					count++;

					//printf("%d, ", cluster.cluster_size[i]);
				}
				//printf("\n");


				break;
			}

		case 45:
			{

				cluster.num_clusters = 0;


				for (i=0;i<num_plates;i++)
				{
					if (plates[i].posi45[ply1.order][1]!=0)
					{
						cluster.cluster_size[cluster.num_clusters] =1;
						cluster.container[cluster.num_clusters] = (int *)malloc(sizeof(int));
						cluster.container[cluster.num_clusters][0] = i;
						cluster.num_clusters++;
					}

				}


				for (int temp=0; temp<5;temp++)
				{
					do
					{
						isolate = cluster.num_clusters;
						for (ci=0;ci<cluster.num_clusters;ci++)
						{
							for (cj=ci+1;cj<cluster.num_clusters;cj++)
							{

							//	printf("%d, %d, %d\n",cluster.num_clusters, ci, cj);
								for (i=0;i<cluster.cluster_size[ci]; i++)
								{
									for (j=0;j<cluster.cluster_size[cj]; j++)
									{
										int id1 = cluster.container[ci][i];
										int id2 = cluster.container[cj][j];
										
										flag = IsNeighbor(num_plates,  ConnectVec, plates[ cluster.container[ci][i] ].PlateID, plates[ cluster.container[cj][j] ].PlateID);
										if (flag==1)
										{
											MergeTwoClusters(&cluster, ci, cj);
											//printf("after merge: %d, %d, %d\n",cluster.num_clusters, ci, cj);
											break;
										}
									}
									if (flag==1)
										break;
								}
							//	printf("%d, %d\n",cluster.num_clusters, cj);
							}
						}
					} while((cluster.num_clusters<isolate)&&(cluster.num_clusters>1));
				}





				//printf("%d clusters for layer %d\n", cluster.num_clusters, lay1);

				for (i=0;i<cluster.num_clusters;i++)
				{
					clusters[count] = 45;
						for (j=0;j<cluster.cluster_size[i];j++)
							plates[cluster.container[i][j]].posi45[ply1.order][2] = count;
					count++;

					//printf("%d, ", cluster.cluster_size[i]);
				}


				//printf("\n");

				break;


			}

		case -45:
			{

				cluster.num_clusters = 0;


				for (i=0;i<num_plates;i++)
				{
					if (plates[i].posi_45[ply1.order][1]!=0)
					{
						cluster.cluster_size[cluster.num_clusters] =1;
						cluster.container[cluster.num_clusters] = (int *)malloc(sizeof(int));
						cluster.container[cluster.num_clusters][0] = i;
						cluster.num_clusters++;
					}

				}





				do
				{

					isolate = cluster.num_clusters;

					for (ci=0;ci<cluster.num_clusters;ci++)
					{
						for (cj=ci+1;cj<cluster.num_clusters;cj++)
						{
							for (i=0;i<cluster.cluster_size[ci]; i++)
							{

								for (j=0;j<cluster.cluster_size[cj]; j++)
								{

									flag = IsNeighbor(num_plates,  ConnectVec, plates[ cluster.container[ci][i] ].PlateID, plates[ cluster.container[cj][j] ].PlateID);
									if (flag==1)
									{
										MergeTwoClusters(&cluster, ci, cj);
								//		printf("%d, %d, %d\n",cluster.num_clusters, ci, cj);
										break;
									}
								}
								if (flag==1)
									break;
							}
						}
					}
				} while((cluster.num_clusters<isolate)&&(cluster.num_clusters>1));






				//printf("%d clusters for layer %d\n", cluster.num_clusters, lay1);
				for (i=0;i<cluster.num_clusters;i++)
				{
					clusters[count] = -45;
						for (j=0;j<cluster.cluster_size[i];j++)
							plates[cluster.container[i][j]].posi45[ply1.order][2] = count;
					count++;

					//printf("%d, ", cluster.cluster_size[i]);
				}
				//printf("\n");
				break;

			}

		case 90:
			{

				cluster.num_clusters = 0;


				for (i=0;i<num_plates;i++)
				{
					if (plates[i].posi90[ply1.order][1]!=0)
					{
						cluster.cluster_size[cluster.num_clusters] =1;
						cluster.container[cluster.num_clusters] = (int *)malloc(sizeof(int));
						cluster.container[cluster.num_clusters][0] = i;
						cluster.num_clusters++;
					}

				}

				for (int temp=0;temp<5;temp++)
				{
					do
				{

					isolate = cluster.num_clusters;

					for (ci=0;ci<cluster.num_clusters;ci++)
					{
						for (cj=ci+1;cj<cluster.num_clusters;cj++)
						{
							if (cluster.num_clusters<=1)
								break;

							for (i=0;i<cluster.cluster_size[ci]; i++)
							{

								for (j=0;j<cluster.cluster_size[cj]; j++)
								{

									flag = IsNeighbor(num_plates,  ConnectVec, plates[ cluster.container[ci][i] ].PlateID, plates[ cluster.container[cj][j] ].PlateID);
									if (flag==1)
									{
										//printf("%d, %d, %d\n",cluster.num_clusters, ci, cj);
										MergeTwoClusters(&cluster, ci, cj);

										break;
									}
								}
								if (flag==1)
									break;
							}
						}
					}
				} while((cluster.num_clusters<isolate)&&(cluster.num_clusters>1));

				}


				//printf("%d clusters for layer %d\n", cluster.num_clusters, lay1);
				for (i=0;i<cluster.num_clusters;i++)
				{
					clusters[count] = 90;
						for (j=0;j<cluster.cluster_size[i];j++)
							plates[cluster.container[i][j]].posi90[ply1.order][2] = count;
					count++;

					//printf("%d, ", cluster.cluster_size[i]);
				}
				//printf("\n");
				break;

			}

		}




	}


	*num_clusters = count;
}


void MergeTwoClusters(struct Clusters *cluster0, int cluster_i, int cluster_j)
{
	int i0, j0, i, count1, count2;
	int *temp;



	i0 = MIN(cluster_i,cluster_j);
	j0 = MAX(cluster_i,cluster_j);



	count1 = cluster0->cluster_size[i0];
	count2 = cluster0->cluster_size[j0];

	temp = (int *)malloc(sizeof(int)*count1);
	for (i=0;i<count1;i++)
		temp[i] = cluster0->container[i0][i];


	cluster0->num_clusters = cluster0->num_clusters-1;



	cluster0->cluster_size[i0] = count1 + count2;
	cluster0->container[i0] = (int *)malloc(sizeof(int)* cluster0->cluster_size[i0]);
	for (i=0;i<count1;i++)
		cluster0->container[i0][i] = temp[i];
	for (i=0;i<count2;i++)
		cluster0->container[i0][i+count1] = cluster0->container[j0][i];


	for (i=j0;i<cluster0->num_clusters;i++)
	{
		cluster0->cluster_size[i] = cluster0->cluster_size[i+1];
		cluster0->container[i] = cluster0->container[i+1];
	}



	cluster0->container[cluster0->num_clusters]=NULL;

	free(temp);


//	printf("success %4d,%4d,%4d\n", cluster0->num_clusters+1, cluster_i, cluster_j );
	/*
	for (i=0;i<cluster0->num_clusters;i++)
	{
		printf("\n%4d,%4d,", i,cluster0->cluster_size[i] );
		for (j=0;j<cluster0->cluster_size[i];j++)
			printf("%d,",cluster0->container[i][j] );
	}
	*/

}
struct ABDM IsoABD(struct LaminaQ Q, double t)
{
	struct ABDM abd;
	abd = ABDM();

	abd.A11 = Q.Q11 * t;
	abd.A12 = Q.Q12 * t;
	abd.A16 = Q.Q16 * t;
	abd.A22 = Q.Q22 * t;
	abd.A26 = Q.Q26 * t;
	abd.A66 = Q.Q66 * t;

	abd.D11 = Q.Q11 * t*t*t/12.;
	abd.D12 = Q.Q12 * t*t*t/12.;
	abd.D16 = Q.Q16 * t*t*t/12.;
	abd.D22 = Q.Q22 * t*t*t/12.;
	abd.D26 = Q.Q26 * t*t*t/12.;
	abd.D66 = Q.Q66 * t*t*t/12.;

	return(abd);
}











void EstimateThicknessIso(struct PlateCopula* plates, struct PlateLoad *pl, int num_subcase, double epsi_al_t, double epsi_al_c)
{
	int i, j, k, flag;
	double thick, min_ms;
	struct LaminaQ Q, q[4];
	struct MarginOfSafety *ms;


	ms = (struct MarginOfSafety*)malloc(sizeof(struct MarginOfSafety)*num_subcase);

	q[0] = GlobalLaminaQ(plates->plate.stack.ply[0], 0.);
	q[1] = GlobalLaminaQ(plates->plate.stack.ply[0], 45.);
	q[2] = GlobalLaminaQ(plates->plate.stack.ply[0], -45.);
	q[3] = GlobalLaminaQ(plates->plate.stack.ply[0], 90.);

	Q =LaminaQ();
	for (i=0;i<4;i++)
	{
		Q.Q11 = Q.Q11 + q[i].Q11/4.;
		Q.Q12 = Q.Q12 + q[i].Q12/4.;
		Q.Q16 = Q.Q16 + q[i].Q16/4.;
		Q.Q22 = Q.Q22 + q[i].Q22/4.;
		Q.Q26 = Q.Q26 + q[i].Q26/4.;
		Q.Q66 = Q.Q66 + q[i].Q66/4.;
	}


	for (k=12; k<MAXLAYEROFPLY; k=k+4)
	{
		thick =(double)k* plates->plate.stack.ply[0].thick; 
			plates->plate.stack.abdm = IsoABD(Q,thick);
				
		flag=1;	
		
		for (j=0; j<num_subcase; j++)
		{
			ms[j] = MarginSafety(plates->plate, pl[j], epsi_al_t, epsi_al_c, 1., 1., 0.);
			if (ms[j].ms <0.)
			{
				flag=0;
				break;
			}
		}

		
		
		if (flag==1)// ALL ms are positive
		{
			min_ms = 10000.;

			for(j=0; j<num_subcase; j++)
			{
				if (ms[j].ms<min_ms)
				{
					plates->pl =  pl[j];
					plates->ms = ms[j];
					min_ms = ms[j].ms;
					i=j;
				}
			}





			plates->pl = pl[i];
			plates->plate.stack.NumLayers = k;
			
			return;
		}
		


	}



	free(ms);
}



void UpdateCriticalCase(struct PlateCopula* plates, struct PlateLoad *pl, int num_subcase, double epsi_al_t, double epsi_al_c, double PB_coef_c, double PB_coef_s, double Smear_Ratio)
{
	int j, temp;
	struct MarginOfSafety ms;
	


	temp = plates->pl.LoadID;

	for(j=0; j<num_subcase; j++)
	{
		if (pl[j].LoadID == plates->pl.LoadID)
			continue;
		ms = MarginSafety(plates->plate, pl[j], epsi_al_t, epsi_al_c, PB_coef_c, PB_coef_s, Smear_Ratio);
		
		if (ms.ms< plates->ms.ms )
		{
			plates->pl =  pl[j];
			plates->ms = ms;
			printf("Critical subcase for plate %d has been updated from %d to %d\n", plates->FEM_ID, temp,plates->pl.LoadID );
		}
	}
	return;
}


int  PlyBundleList(int list[][4], int num_layer)
{
	int i, count;
	int lb, ub;
	int num0,num45;

	count=0;
	for (i=2*((int)(0.2*num_layer))+2; i<=2*((int)(0.4*num_layer)); i=i+2)
	{
		num45 = i;

		if (i>0.7*num_layer)
		{
			lb = ((int)(0.05*num_layer))*2;
			ub = num_layer - num45 - (int)(0.1*num_layer);
		}
		else
		{
			lb = MAX(2*(int)((num_layer - (0.6*num_layer) - num45)/2), ((int)(0.05*num_layer))*2) ;
			ub = MIN(0.6*num_layer, num_layer - (0.1*num_layer) - num45);
		}

		for (num0=lb;num0<=ub;num0=num0+2)
		{
			list[count][0] = num0;
			list[count][1] = num45/2;
			list[count][2] = num45/2;
			list[count][3] = num_layer-num0-num45;
			count++;
		}

	}
	return(count);
}


void initialize_and_permute(int index[], int n)
{
	int i, j, temp;
	for (i=0;i<n;i++)
		index[i] = i;

	for (i=1;i<n-2;i++)
	{
		j = RandU(i);

		temp = index[i];

		index[i] = index[j];

		index[j] = temp;
	}


}




void Adjust45Pair(struct PlateCopula plates[], int num_plates, int cover_thick)
{
	int pos_index, i, j, lay1;
	ud ply, ply2;
	
	pos_index = cover_thick;




		for (i=cover_thick;i<(plates[num_plates-1].plate.stack.NumLayers/2);i++)
		{
		
			ply = Find(&plates[num_plates-1], i);


			if ((ply.angle==0)||(ply.angle==90))
			{

				continue;
			}
		
			
			
			
			else
			{
				for (j=i+1;j<(plates[num_plates-1].plate.stack.NumLayers/2);j++)
				{
					ply2 = Find(&plates[num_plates-1], j);

					if ((ply2.angle==0)||(ply2.angle==90))
						continue;

					else
					{

						if (ply.angle==45)
							lay1 = plates[num_plates-1].posi_45[ply.order][0];
						else
							lay1 = plates[num_plates-1].posi45[ply.order][0];
						
						if (lay1!=j)
							Exchange2Plies0(plates, num_plates, lay1, j);

						//printf("Exchange, %d ,%d\n",  lay1, j);
						
						i=j;

						break;

					}

				}

			}




		
		
		}





}