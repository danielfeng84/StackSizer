// This file is part of StackSizer, a kit for stack design of successive laminates.
//
// Copyright (C) 2018 Jianwen Feng <jianwenfeng@hotmail.com>
//
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// 
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

// ����Ȩ����Ȩ��ʹ�á����ơ��޸ġ��ϲ������淢�С�ɢ��������Ȩ���������������ĸ�����
// ����Ȩ�˿ɸ��ݳ�ʽ����Ҫ�޸���Ȩ����Ϊ�ʵ������ݡ�
// ���������������и����ж����������Ȩ���������������


// EIGEN3 is an open source project subject to the terms of the Mozilla Public License v. 2.0. 

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
// DEALINGS IN THE SOFTWARE.


#include "Buckling.h"
#include "Composite.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "Buckling.h"
#include "Optimization.h"
#include <time.h>
#include "PlyCheck.h"
#include <omp.h>

//// Please redefine the macro MAXLAYEROFPLY a prior  which defines the max. num of plies.
///  ��ȡ.f06�ļ���CQUAD4��Ԫ�����ĺ���ΪIO.h�ļ��ж����   int ExtractStress(char Quad4File[MAXLEN], char f06file[MAXLEN], char outputfile[MAXLEN])
///  Ŀ�꺯���趨Ϊ�� num_safe = num_safe - ((double)quality)/10.;�����Ը�����Ҫ�޸�



/// Compare two thicknesses

int CompareThick (const void * a, const void * b)
{
	struct PlateCopula * p1 = (struct PlateCopula*)a;
	struct PlateCopula * p2 = (struct PlateCopula*)b;

	if(p1->plate.stack.NumLayers < p2->plate.stack.NumLayers)
		return -1;
	else
	{
		if(p1->plate.stack.NumLayers == p2->plate.stack.NumLayers)
			return 0;
		else
			return 1;
	}
}


/// Tensional strain allowable
double epsi_al_t(struct Laminate plate)
{
	/*
	allowable tensional strain for plates is 4500
	*/

	return(4000.e-6);

}


/// Compressive strain allowable
/// If the laminate thickness is thinner than 8mm, then the allowed compressive strain is -3300, otherwise -3800
double epsi_al_c(struct Laminate plate)
{
	/*
	allowable compressive strain for plates thinner than 8 mm is 3300, for plates thicker than 8 mm is 3800
	*/
	double t;
	t = PlateThick(plate);

	if (t<10.)
		return(-3000.e-6);
	else
		return(-3500.e-6);

}



/// Compressive postbuckling allowable
/// If the laminate is thinner than 2mm, then compressive buckling is allowed to occur at 67% ultimate load level;
/// If the laminate is thicker than 5mm, then compressive buckling is not allowed up to ultimate load
double PB_coef_c(struct Laminate plate)
{
	double t;
	t = PlateThick(plate);

	if (t<=2.)
		return(0.67);
	else
		if (t>=5.)
			return(1.);
		else
			return(t/9.+4./9.);
}



/// Shear postbuckling allowable
/// If the laminate is thinner than 2mm, then shear buckling is allowed to occur at 80% ultimate load level;
/// If the laminate is thicker than 5mm, then shear buckling is not allowed up to ultimate load
double PB_coef_s(struct Laminate plate)
{
	double t;
	t = PlateThick(plate);

	if (t<=2)
		return(0.8);
	else
		if (t>=5.)
			return(1.);
		else
			return(t/15.+ 2./3.);
}



/// ��ѯID���ɹ��򷵻ص�Ԫ�ţ����ɹ��򷵻�-1
///
///
/// @param
/// @return ID
int QueryID(int id, const struct PlateCopula *bay, int num_bay)
{
	int i;
	for (i=0;i<num_bay;i++)
	{
		if (bay[i].FEM_ID==id)
			return(i);
	}

	printf("Element ID cannot be found for %d\n", id);
	return(-1);
}








int main()
{

	struct PlateCopula *candidate, *candidate2, *optimal_solution;
	struct CompsiteMat ply;
	struct PlateLoad **LoadM;

	int i, j, k, t, RoundNo, id, num_bay, lay1, lay2, num_try, max_layer,  min_layer, count, quality, best, list[MAXLAYEROFPLY], opera, elem1, elem2, num_good, num_field,  flag_size, pace, num_size, num_lock, switch_trim;
	int *ConnectVec, *clusters, *v_case_id;
	int num_round, num_subcases, num_outer_loop, num_inner_loop, num_outer_loop_inp, num_inner_loop_inp, num_clusters, cover_thick;
	double num_safe, num_safe_trim, num_safe_max, dice, a, b, reserved_ms;
	char buffer[50];
	char OutFileName[50];
	FILE *fp, *fp2, *fp_in, *fp_log;


	min_layer = 12;                 // Minimum gage  ��С���Ϊ12��
	cover_thick = 2;				// The outermost two cover layers �������ι涨Ϊ45 ��-45



	//�������Զ���
	ply.Ex = 154.e3;
	ply.Ey = 8.5e3;
	ply.Gxy = 4.2e3;
	ply.MIUxy = 0.35;
	ply.thick = 0.184;




	fp_log = fopen("log.log","w"); //Log�ļ�
	if (fp_log == NULL)
	{
		printf("ERROR: Fail to open the log file.\n");
		return(0);
	}





	//���������ļ�
	printf("Input the input file name\n");
	scanf("%s", buffer);
	fp_in = fopen(buffer,"r");//fp_in:�����ļ���ָ��
	if (fp_in == NULL)
	{
		printf("ERROR: Fail to open the input file %s\n", buffer);
		fprintf(fp_log, "ERROR: Fail to open the input file %s\n", buffer);
		return(0);
	}




	fscanf(fp_in, "%d%d%d%d%lf", &num_bay, &num_subcases, &num_size, &num_lock, &reserved_ms); //�����һ��
	fscanf(fp_in, "%d%d%d%d", &num_round, &num_inner_loop_inp, &num_outer_loop_inp, &switch_trim);//����ڶ���




	printf("Number of skin bays: %d\n", num_bay);  //num_bay: ��Ƥ������Ŀ
	printf("Number of Subcases: %d\n", num_subcases); //num_subcase��������Ŀ
	printf("Reserved margin of safety: %f\n", reserved_ms);//reserved_ms�� ����ԣ��


	num_round = MAX(num_round, 3); //num_round: ����������������Ҫ��������




	/// Mutural connection matrix ����ConnectVec��ŵ�Ԫ�������ӵ���Ϣ������Ϊnum_bay*num_bay
	ConnectVec = (int *)malloc(sizeof(int)*num_bay*num_bay);
	if (ConnectVec==NULL)
	{
		printf("ERROR: Connection matrix Mallocation Error\n");
		fprintf(fp_log,"ERROR: Connection matrix Mallocation Error\n" );
		return(0);
	}
	for (i=0;i<num_bay*num_bay;i++)
		ConnectVec[i]=0;
	for (i=0;i<num_bay;i++)
		ConnectVec[i+i*num_bay]=1; //��Ԫ�Լ����Լ�����






	//LoadM�غɾ���num_bay��*num_subcases
	//����ָ�����飬
	LoadM = (struct PlateLoad**)malloc(sizeof(struct PlateLoad*)*num_bay);
	if( LoadM == NULL )
	{
		printf( "ERROR: Loading matrix Allocation Error\n" );
		fprintf(fp_log, "ERROR: Loading matrix Allocation Error\n"  );
		return(0);
	}


	//��������һά����
	for(i=0;i<num_bay;i++)
	{
		LoadM[i] = (struct PlateLoad*)malloc(sizeof(struct PlateLoad)*num_subcases);
		if( LoadM[i] == NULL )
		{
			printf( "ERROR: Loading matrix Allocation Error\n" );
			fprintf(fp_log, "ERROR: Loading matrix Allocation Error\n");
			return(0);
		}

	}


	//v_case_id:��Ź������ƣ�����������Ϊnum_subcases
	v_case_id = (int *)malloc(sizeof(int)*num_subcases);
	if (v_case_id==NULL)
	{
		printf("ERROR: Malloc Error\n");
		fprintf(fp_log, "ERROR: Malloc Error\n");
		return(0);
	}


	clusters = (int *)malloc(sizeof(int)*num_bay);
	if (clusters==NULL)
	{
		printf("ERROR: Malloc Error\n");
		fprintf(fp_log, "ERROR: Malloc Error\n");
		return(0);
	}



	candidate = (struct PlateCopula*) malloc(sizeof(struct PlateCopula)*num_bay);
	if (candidate==NULL)
	{
		printf("ERROR: Malloc Error\n");
		fprintf(fp_log, "ERROR: Malloc Error\n" );
		return(0);
	}


	candidate2 = (struct PlateCopula*) malloc(sizeof(struct PlateCopula)*num_bay);
	if (candidate2==NULL)
	{
		printf("ERROR: Malloc Error\n");
		fprintf(fp_log, "ERROR: Malloc Error\n" );
		return(0);
	}


	optimal_solution = (struct PlateCopula*) malloc(sizeof(struct PlateCopula)*num_bay);
	if (candidate==NULL)
	{
		printf("ERROR: Malloc Error\n");
		fprintf(fp_log, "ERROR: Malloc Error\n" );
		return(0);
	}




	//��ʼ��candidate
	for (i=0;i<num_bay;i++)
	{
		candidate[i] = PlateCopula();

	}

	for (j=0;j<num_subcases;j++)
	{
		fscanf(fp_in,"%s%d",buffer,&(v_case_id[j]));//���빤������

		for (i=0;i<num_bay;i++)
		{
			fscanf(fp_in, "%d%lf%lf%lf",  &k, &(LoadM[i][j].InPlane.Nxx),&(LoadM[i][j].InPlane.Nyy),&(LoadM[i][j].InPlane.Nxy) );//�����˳��Ϊ�� ��Ԫ�� FX FY, FXY
			LoadM[i][j].LoadID = v_case_id[j];

			if (j==0)
			{
				candidate[i].FEM_ID = k;//��ŵ�Ԫ��
			}
			else
			{
				if(candidate[i].FEM_ID != k)//.f06�ļ����չ���˳�򱨳����е�Ԫ���������ڲ�ͬ�����µ�Ԫ��˳����Ҫһ�¡��������鵥Ԫ���һ����
				{
					printf("ERROR: Element order is incorrect\n");
					fprintf(fp_log, "ERROR: Element order is incorrect\n" );
					return(0);
				}
			}
		}
	}



	for (j=0;j<num_size;j++)//num_size: �ߴ綨�������Ŀ��
	{
		flag_size = fscanf(fp_in,"%s%lf%lf\n",buffer, &a, &b);
		//printf("%s  %lf  %lf;\n",buffer, a, b);

		num_field = sscanf(buffer, "%d:%d:%d", &elem1, &elem2,  &pace);

		switch (num_field)
		{
		case 1://���1�� ���˵�һ����
			i = QueryID(elem1, candidate, num_bay);
			if (i==-1)//����Ҳ��������ĵ�Ԫ�����������
				continue;
			candidate[i].plate.a = a;
			candidate[i].plate.b = b;

			break;

		case 2://���2��������ʼ��Ԫ�ͽ�����Ԫ������Ϊ1
			for(k=elem1;k<=elem2;k++)
			{
				i = QueryID(k, candidate, num_bay);
				if (i==-1)//����Ҳ��������ĵ�Ԫ�����������
					continue;
				candidate[i].plate.a = a;
				candidate[i].plate.b = b;
			}

			break;

		case 3://���3��������ʼ��Ԫ��������Ԫ�Լ�����pace
			for(k=elem1;k<=elem2;k=k+pace)
			{
				i = QueryID(k, candidate, num_bay);
				if (i==-1)//����Ҳ��������ĵ�Ԫ�����������
					continue;
				candidate[i].plate.a = a;
				candidate[i].plate.b = b;
			}
			break;

		default:
			{
				printf("ERROR: Unrecognized block\n");
				fprintf(fp_log,"ERROR: Unrecognized block\n" );
				return(0);
			}

		}


	}



	for (j=0;j<num_lock;j++)//num_lock: ������ȵ������Ŀ��
	{
		flag_size = fscanf(fp_in,"%s%d", buffer, &t);
		num_field = sscanf(buffer, "%d:%d:%d", &elem1, &elem2,  &pace);
		switch (num_field)
		{
		case 1:
			i = QueryID(elem1, candidate, num_bay);
			if (i==-1)
				continue;
			candidate[i].lock=1;
			candidate[i].plate.stack.NumLayers= 4*int(t/4);//�����Ҫ4�ı�������������Ĳ���4�ı�����������

			break;

		case 2:
			for(k=elem1;k<=elem2;k++)
			{
				i = QueryID(k, candidate, num_bay);
				if (i==-1)
					continue;
				candidate[i].lock=1;
				candidate[i].plate.stack.NumLayers= 4*int(t/4);//�����Ҫ4�ı���
			}

			break;

		case 3:
			for(k=elem1;k<=elem2;k=k+pace)
			{
				i = QueryID(k, candidate, num_bay);
				if (i==-1)
					continue;
				candidate[i].lock=1;
				candidate[i].plate.stack.NumLayers= 4*int(t/4);//�����Ҫ4�ı���
			}
			break;
		default:
			{
				printf("ERROR: Unrecognized block\n");
				fprintf(fp_log, "ERROR: Unrecognized block\n" );
				return(0);
			}
		}
	}





	//���ݼ�飬�鿴�Ƿ����е�Ԫ�������˳ߴ�
	for(i=0;i<num_bay;i++)
	{
		if ((candidate[i].plate.a<1.0e-6)||(candidate[i].plate.b<1.0e-6))
		{
			printf("WARNING: Miss size property of element %d\n",candidate[i].FEM_ID );
			fprintf(fp_log, "WARNING: Miss size property of element %d\n",candidate[i].FEM_ID );
			candidate[i].plate.a = 600.;//���û�и����ߴ磬������Ϊ��635����150��������������Ϣ
			candidate[i].plate.b = 100.;
		}

		if ((candidate[i].lock==1)&&(candidate[i].plate.stack.NumLayers==0))//��������˺�ȣ�����û�и����ֵ����������ֵΪ12��������������Ϣ
		{
			printf("WARNING: Miss thickness property of element %d\n",candidate[i].FEM_ID );
			fprintf(fp_log,  "WARNING: Miss thickness property of element %d\n",candidate[i].FEM_ID);
			candidate[i].plate.stack.NumLayers = min_layer;

		}
	}




	//����������Ϣ��ֱ���ļ�������������Ϣ����Ŀ������Ҫ����
	while(!feof(fp_in))
	{
		fscanf(fp_in, "%d %d",  &elem1, &elem2);
		i = QueryID(elem1, candidate, num_bay);
		j = QueryID(elem2, candidate, num_bay);
		ConnectVec[i+j*num_bay]=1;
		ConnectVec[j+i*num_bay]=1;
	}


	fclose(fp_in);//�ر������ļ�






	///Pickup crital subcases for each bay and estimate thickness �²������Ƥ�������ع����������Ƶ�һ�ֺ��
	for (j=0;j<num_bay;j++)
	{
		candidate[j].PlateID  = j; //PlateIDΪ���ţ���0��num_bay-1
		candidate[j].DropCount =0; //DropCount��¼�������

		for (k=0;k<MAXLAYEROFPLY;k++)
			candidate[j].plate.stack.ply[k] = ply;

		if (candidate[j].lock ==0)
			//Estimate thickness,���ú���EstimateThicknessIso���ٶ�����ͬ�ԣ�����������������ѹ����Ӧ��
			EstimateThicknessIso(&(candidate[j]), LoadM[j], num_subcases, 4500e-6, -3300e-6);

		candidate[j].plate.stack.NumLayers = MAX(candidate[j].plate.stack.NumLayers, min_layer);//������Ҫ12��
		if ((candidate[j].plate.stack.NumLayers%4)!=0)
			candidate[j].plate.stack.NumLayers = candidate[j].plate.stack.NumLayers + 2;//�����Ҫ4�ı���


	}


	//// Start optimization process
	for (RoundNo=1; RoundNo<=num_round; RoundNo++)
	{
		printf("Starting round %d\n", RoundNo);
		fprintf(fp_log, "Starting round %d\n", RoundNo);
		num_good = 0;

		if (RoundNo==(num_round))
		{
			num_outer_loop = num_outer_loop_inp; //���һ�ۼ���ʱ����Ҫ��ϸһЩ��ʹ���û�������ѭ����
			num_inner_loop = num_inner_loop_inp;
		}
		else
		{
			num_outer_loop = 100;//100�Ǹ��ݾ������
			num_inner_loop = 100;//100�Ǹ��ݾ������
		}


		if (RoundNo!=1)
		{
			for (j=0;j<num_bay;j++)
			{
				//optimal_solution�ǵõ��ĺ��ʵ��̲㷽�����˴��൱�ڰ�optimal_solution����һ�ݸ�candidate
				candidate[j] = optimal_solution[j];


				//update subcase here
				//��һ���̲������һ�ֵ����ع�������õ��������п����ڱ�Ĺ�����ԣ�ȸ���
				UpdateCriticalCase(&(candidate[j]), LoadM[j], num_subcases, epsi_al_t(candidate->plate), epsi_al_c(candidate->plate), PB_coef_c(candidate->plate) , PB_coef_s(candidate->plate), 0.);


				if ((candidate[j].ms.ms < reserved_ms)&&(candidate[j].lock==0)) // less than reserved ms�� ���ԣ�ȵ��ڸ����ı���ԣ�ȣ�������4��
				{
					candidate[j].plate.stack.NumLayers = MIN(candidate[j].plate.stack.NumLayers + 4, MAXLAYEROFPLY);

				}


				if ((candidate[j].ms.ms > 1.0)&&(candidate[j].lock==0))//over conservative��ԣ�ȸ���1��Ϊ���ڱ��أ���ȥ4�㣬����Ϊ�˱��ⷴ������-��ԣ�Ȳ���-���Ӻ�-��ԣ�ȹ��ߵ�ѭ��������ʹ��DropCount����
				{
					if(candidate[j].DropCount<2)
					{
						candidate[j].plate.stack.NumLayers = MAX(candidate[j].plate.stack.NumLayers - 4, min_layer);
						candidate[j].DropCount++;
					}

				}
			}
		}



		//���ݺ�ȶ԰嵥Ԫ��������
		qsort (candidate, num_bay, sizeof(struct PlateCopula), CompareThick);
		max_layer = candidate[num_bay-1].plate.stack.NumLayers;//�����


		EsimatePercentFirst(candidate); // ������С��Ӧ���ԭ�򣬹������Ԫ���̲����

		for (j=1;j<num_bay;j++)
		{
			EsimatePercent(candidate, j);//���մӱ������˳����Ƹ�����Ƥ���̲����
		}


		for (j=0;j<num_bay;j++)
		{
			candidate[j].half_nums[0] = candidate[j].dirc[0]/2;
			candidate[j].half_nums[3] = candidate[j].dirc[3]/2;

			if (candidate[j].dirc[1]%2==0)
			{
				candidate[j].half_nums[1] = (candidate[j].dirc[1]/2);
				candidate[j].half_nums[2] = candidate[j].half_nums[1];
			}

			else
			{
				candidate[j].half_nums[1] = (candidate[j].dirc[1]+1)/2;
				candidate[j].half_nums[2] = candidate[j].half_nums[1]-1;
			}
		}

		InitialPlyDeploy(candidate, num_bay);//�����̲㣬���մӺ񵽱���˳�򣬴˴��̲�����������Ҫ�󣬵��ǲ������̲����ԭ��



		//�������
		PartitionCalc(clusters, &num_clusters, candidate, num_bay, max_layer, ConnectVec);




		for (i=0;i<num_bay;i++)
			optimal_solution[i]=candidate[i];


		num_safe_max = 0.;
		num_safe = 0.;

		for (count=0;count<num_outer_loop;count++)
		{
			
			
			
			
			if (num_safe > num_safe_max)//�������Ž�
			{
				num_safe_max = num_safe;
				for (i=0;i<num_bay;i++)
					optimal_solution[i]=candidate2[i];  //update the candidate
				printf("strength satidfied %f,   %d\n",  num_safe_max, count+1);
				fprintf(fp_log, "strength satidfied %f,   %d\n",  num_safe_max, count+1 );
			}

			else//����һ�����ʸ������Ž�
			{
				dice = MAX(((double)(rand()%10000)/10000.),1.0e-10);
				double temp = pow(3.,(num_safe-num_safe_max)/0.005);

				if ((dice >= temp)||(count==0))//rej
				{
					for (i=0;i<num_bay;i++)
						candidate2[i]=optimal_solution[i];
				}
			}



			srand(time(NULL));//�����������������
			best = num_bay*10;
			num_try =0;

			while ((best!=0)&&(num_try<num_inner_loop))
			{
				//
				
				
				num_try++;


				//���ѡ�����㣨�������������cover_thick��
				lay1 = rand()%(max_layer/2-cover_thick)+cover_thick;
				lay2 = rand()%(max_layer/2)+cover_thick;




				if (lay2<max_layer/2)
				{
					opera =  Exchange2Plies(candidate2, num_bay, lay1, lay2);//����lay1�� lay2�����λ��
					quality=0;

					for (i=0;i<num_bay;i++)
					{
						quality = quality + plycheck1(candidate2[i].plate.stack) + plycheck5(candidate2[i].plate.stack)+ plycheck9(candidate2[i].plate.stack);
						//plycheck1: �������Ĳ������̲�
						//plycheck5: ʹ45���̲�ʹ-45���̲㾡���ɶԳ��֣��������45/45/-45/-45
						//plycheck9: ���������̲��ĵ����㲻����ǰ��������
						//���Ը�����Ҫ��Ӹ���׼���̲���Ĵ�����PlyCheck.h�ļ���
						// �̲���ͨ������0������qualityԽСԽ��
						//�̲���������������㣬�����̲������Ϊһ��Ȩ�ؼ��ɵ�ԣ����
					}
					if (quality<=best)
					{
						best=quality;
					}
					else
					{
						if (opera==1)
						{
							Exchange2Plies0(candidate2, num_bay, lay2, lay1); //�µõ����̲����ڲ������̲����ԭ���ٻ���ȥ��

						}
					}

				}

				else
				{

					opera =  SinkPly(candidate2, num_bay, lay1);//��lay1�������洦
					quality=0;

					for (i=0;i<num_bay;i++)
					{
						quality = quality + plycheck1(candidate2[i].plate.stack) + plycheck5(candidate2[i].plate.stack)+ plycheck9(candidate2[i].plate.stack);
					}

					if (quality<=best)
					{
						best=quality;
					}
					else
					{
						if (opera==1)
						{

							RestorePly(candidate2, num_bay, lay1);//�µõ����̲����ڲ������̲����ԭ���ٻ���ȥ��

						}
					}
				}

				if (quality==0)//�����̲㶼�������׼��
					break;


			}
			printf("Obtained a candidate after %d tries, SQI %d\n",num_try, quality);
			fprintf(fp_log, "Obtained a candidate after %d tries, SQI %d\n",num_try, quality );


			Adjust45Pair(candidate2, num_bay, cover_thick);

			//�̲��������Ҫ���ж�ԣ��
			num_safe = 0.;
#pragma omp parallel for
			for (int para_i=0;  para_i<num_bay;  para_i++)
			{
				candidate2[para_i].plate.stack.abdm = ABD(candidate2[para_i].plate.stack);
				candidate2[para_i].ms = MarginSafety( candidate2[para_i].plate, candidate2[para_i].pl, epsi_al_t(candidate2[para_i].plate), epsi_al_c(candidate2[para_i].plate), PB_coef_c(candidate2[para_i].plate) , PB_coef_s(candidate2[para_i].plate), 0.);
				
			}
			
			for (i=0;  i<num_bay;  i++)
			{
				if (candidate2[i].ms.ms > reserved_ms)
					num_safe = num_safe +1.;
			}


			if (num_safe>num_safe_max)
			{

				if ((switch_trim==1)&&(RoundNo==num_round))//΢����
				{
					for (j=2;j<num_clusters;j++)
					{
						opera = Reversal(candidate2, num_bay, clusters, num_clusters, j);//0��Ϊ90��90��Ϊ0������ò�Ϊ45����45���򷵻�0
						if (opera==0)
							continue;

						quality=0;
						for (i=0;i<num_bay;i++)
						{
							quality = quality + plycheck1(candidate2[i].plate.stack) + plycheck5(candidate2[i].plate.stack)+ plycheck9(candidate2[i].plate.stack);
						}

						if (quality> (num_bay/20))
						{
							opera = Reversal(candidate2, num_bay, clusters, num_clusters, j);

							for (i=0;i<num_bay;i++)
							{
								candidate2[i].plate.stack.abdm = ABD(candidate2[i].plate.stack);
								candidate2[i].ms = MarginSafety( candidate2[i].plate, candidate2[i].pl, epsi_al_t(candidate2[i].plate), epsi_al_c(candidate2[i].plate), PB_coef_c(candidate2[i].plate), PB_coef_s(candidate2[i].plate), 0.);
							}

							continue;
						}



						num_safe_trim = 0.;
						for (i=0;i<num_bay;i++)
						{
							candidate2[i].plate.stack.abdm = ABD(candidate2[i].plate.stack);
							candidate2[i].ms = MarginSafety( candidate2[i].plate, candidate2[i].pl, epsi_al_t(candidate2[i].plate), epsi_al_c(candidate2[i].plate), PB_coef_c(candidate2[i].plate), PB_coef_s(candidate2[i].plate), 0.);
							num_safe_trim =num_safe_trim+ candidate2[i].ms.ms;
							if (candidate2[i].ms.ms>reserved_ms)
								num_safe_trim = num_safe_trim+1.;
						}
						num_safe_trim = num_safe_trim - ((double)quality)/10.;


						if (num_safe_trim<num_safe)
						{
							opera = Reversal(candidate2, num_bay, clusters, num_clusters, j);

							for (i=0;i<num_bay;i++)
							{
								candidate2[i].plate.stack.abdm = ABD(candidate2[i].plate.stack);
								candidate2[i].ms = MarginSafety( candidate2[i].plate, candidate2[i].pl, epsi_al_t(candidate2[i].plate), epsi_al_c(candidate2[i].plate), PB_coef_c(candidate2[i].plate), PB_coef_s(candidate2[i].plate), 0.);
							}
						}

						else
						{

							printf("Trim accepted: %d,   %d ----> %d,  %f,  %f\n", j, abs(90-clusters[j]),clusters[j], num_safe, num_safe_trim);
							fprintf(fp_log, "Trim accepted: %d,   %d ----> %d,  %f,  %f\n", j, abs(90-clusters[j]),clusters[j], num_safe, num_safe_trim);
							num_safe = num_safe_trim;
						}
					}
				}
				//΢������





				//������


				sprintf(OutFileName,"Solution_%d_%d.txt", RoundNo, ++num_good);

				fp =fopen(OutFileName,"w");
				if (fp == NULL)
				{
					printf("ERROR: Fail to open the output file.\n");
					fprintf(fp_log, "ERROR: Fail to open the output file.\n");
					return(0);
				}
				/*
				Plate No:         ��ţ����ݶ���˳������
				FEM ELEMENT ID:  FEM�嵥Ԫ��
				RANK:			  ���ݺ�ȵ�����
				Subcase:		  ��Σ�չ���
				Layers:			  �̲���
				M.S:              ��Σ�չ�����Ӧ��ԣ��
				FailMode:         ʧЧģʽ��1��3ΪӦ������ֵ���ƣ�2Ϊ�ȶ��Կ���
				*/



				sprintf(OutFileName,"StackSequence_%d_%d.txt", RoundNo, ++num_good);
				fp2=fopen(OutFileName,"w");
				if (fp2 == NULL)
				{
					printf("ERROR: Fail to open the output file.\n");
					fprintf(fp_log, "ERROR: Fail to open the output file.\n");
					return(0);
				}




				for (id=0;id<num_bay;id++)
				{
					for (i=0;i<num_bay;i++)
					{
						if (candidate2[i].PlateID == id)
							break;
					}

					fprintf(fp2,"%8d,%4d,",candidate2[i].FEM_ID,candidate2[i].plate.stack.NumLayers);
					for (j=0;j<candidate2[i].plate.stack.NumLayers;j++)
						fprintf(fp2,"%4d,",candidate2[i].plate.stack.PlyDirection[j]);
					fprintf(fp2,"\n");


					fprintf(fp,"Plate No.: %4d, FEM ELEMENT ID: %8d, RANK: %8d, Subcase: %8d, Layers: %4d, M.S: %10.2f, FailMode: %d,", id+1, candidate2[i].FEM_ID, i+1, candidate2[i].ms.subcase, candidate2[i].plate.stack.NumLayers, candidate2[i].ms.ms, candidate2[i].ms.mode);
					for (j=0;j<MAXLAYEROFPLY;j++)
						list[j] = DROPEDPLY;

					for (j=0;j<MAXLAYEROFPLY;j++)
					{
						if (candidate2[i].posi0[j][1]!=0)
						{
							list[candidate2[i].posi0[j][0]] = candidate2[i].posi0[j][3];
							list[max_layer-1-candidate2[i].posi0[j][0]] = candidate2[i].posi0[j][3];
						}



						if (candidate2[i].posi45[j][1]!=0)
						{
							list[candidate2[i].posi45[j][0]]=45;
							if (candidate2[i].posi45[j][1] == 1)
							{
								list[max_layer-1-candidate2[i].posi45[j][0]]=45;
							}


							else
							{
								if (candidate2[i].half_nums[1]!=candidate2[i].half_nums[2])
									list[max_layer-candidate2[i].posi45[j][0]]=-45;
								else
									list[max_layer-1-candidate2[i].posi45[j][0]]=45;//
							}
						}



						if (candidate2[i].posi_45[j][1]!=0)
						{
							list[candidate2[i].posi_45[j][0]]=-45;
							list[max_layer-1-candidate2[i].posi_45[j][0]]=-45;
						}




						if (candidate2[i].posi90[j][1]!=0)
						{
							list[candidate2[i].posi90[j][0]] = candidate2[i].posi90[j][3];
							list[max_layer-1-candidate2[i].posi90[j][0]] = candidate2[i].posi90[j][3];
						}
					}

					for (j=0;j<max_layer/2;j++)
					{
						if (list[j]==DROPEDPLY)
							fprintf(fp, "    ,");
						else
							fprintf(fp,"%4d,", list[j]);

					}
					fprintf(fp,"\n");
				}

				fclose(fp);
				fclose(fp2);
			}

		}

	}







	//�ͷ��ڴ�

	free(clusters);
	free(candidate);
	free(candidate2);
	free(optimal_solution) ;
	free(ConnectVec);
	free(v_case_id);

	for(i=0;i<num_bay;i++)
		free (LoadM[i]);
	free(LoadM);


	printf("Finished!");
	fprintf(fp_log, "Finished!");



	fclose(fp_log);
	return(1);



}
