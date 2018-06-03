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

// 被授权人有权利使用、复制、修改、合并、出版发行、散布、再授权及贩售软件及软件的副本。
// 被授权人可根据程式的需要修改授权条款为适当的内容。
// 在软件和软件的所有副本中都必须包含版权声明和许可声明。


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
///  提取.f06文件中CQUAD4单元线力的函数为IO.h文件中定义的   int ExtractStress(char Quad4File[MAXLEN], char f06file[MAXLEN], char outputfile[MAXLEN])
///  目标函数设定为： num_safe = num_safe - ((double)quality)/10.;，可以根据需要修改



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



/// 查询ID，成功则返回单元号，不成功则返回-1
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


	min_layer = 12;                 // Minimum gage  最小厚度为12层
	cover_thick = 2;				// The outermost two cover layers 最外两次规定为45 ，-45



	//材料属性定义
	ply.Ex = 154.e3;
	ply.Ey = 8.5e3;
	ply.Gxy = 4.2e3;
	ply.MIUxy = 0.35;
	ply.thick = 0.184;




	fp_log = fopen("log.log","w"); //Log文件
	if (fp_log == NULL)
	{
		printf("ERROR: Fail to open the log file.\n");
		return(0);
	}





	//读入输入文件
	printf("Input the input file name\n");
	scanf("%s", buffer);
	fp_in = fopen(buffer,"r");//fp_in:输入文件的指针
	if (fp_in == NULL)
	{
		printf("ERROR: Fail to open the input file %s\n", buffer);
		fprintf(fp_log, "ERROR: Fail to open the input file %s\n", buffer);
		return(0);
	}




	fscanf(fp_in, "%d%d%d%d%lf", &num_bay, &num_subcases, &num_size, &num_lock, &reserved_ms); //读入第一行
	fscanf(fp_in, "%d%d%d%d", &num_round, &num_inner_loop_inp, &num_outer_loop_inp, &switch_trim);//读入第二行




	printf("Number of skin bays: %d\n", num_bay);  //num_bay: 蒙皮格子数目
	printf("Number of Subcases: %d\n", num_subcases); //num_subcase：工况数目
	printf("Reserved margin of safety: %f\n", reserved_ms);//reserved_ms： 保留裕度


	num_round = MAX(num_round, 3); //num_round: 迭代轮数，最少需要迭代三轮




	/// Mutural connection matrix 向量ConnectVec存放单元两两连接的信息，长度为num_bay*num_bay
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
		ConnectVec[i+i*num_bay]=1; //单元自己与自己相连






	//LoadM载荷矩阵，num_bay行*num_subcases
	//申请指针数组，
	LoadM = (struct PlateLoad**)malloc(sizeof(struct PlateLoad*)*num_bay);
	if( LoadM == NULL )
	{
		printf( "ERROR: Loading matrix Allocation Error\n" );
		fprintf(fp_log, "ERROR: Loading matrix Allocation Error\n"  );
		return(0);
	}


	//逐行申请一维数组
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


	//v_case_id:存放工况名称，向量，长度为num_subcases
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




	//初始化candidate
	for (i=0;i<num_bay;i++)
	{
		candidate[i] = PlateCopula();

	}

	for (j=0;j<num_subcases;j++)
	{
		fscanf(fp_in,"%s%d",buffer,&(v_case_id[j]));//读入工况名称

		for (i=0;i<num_bay;i++)
		{
			fscanf(fp_in, "%d%lf%lf%lf",  &k, &(LoadM[i][j].InPlane.Nxx),&(LoadM[i][j].InPlane.Nyy),&(LoadM[i][j].InPlane.Nxy) );//读入的顺序为： 单元号 FX FY, FXY
			LoadM[i][j].LoadID = v_case_id[j];

			if (j==0)
			{
				candidate[i].FEM_ID = k;//存放单元号
			}
			else
			{
				if(candidate[i].FEM_ID != k)//.f06文件按照工况顺序报出所有单元的内力，在不同工况下单元的顺序需要一致。这段语句检查单元编号一致性
				{
					printf("ERROR: Element order is incorrect\n");
					fprintf(fp_log, "ERROR: Element order is incorrect\n" );
					return(0);
				}
			}
		}
	}



	for (j=0;j<num_size;j++)//num_size: 尺寸定义语句条目数
	{
		flag_size = fscanf(fp_in,"%s%lf%lf\n",buffer, &a, &b);
		//printf("%s  %lf  %lf;\n",buffer, a, b);

		num_field = sscanf(buffer, "%d:%d:%d", &elem1, &elem2,  &pace);

		switch (num_field)
		{
		case 1://情况1： 给了第一个点
			i = QueryID(elem1, candidate, num_bay);
			if (i==-1)//如果找不到所给的单元，程序会跳过
				continue;
			candidate[i].plate.a = a;
			candidate[i].plate.b = b;

			break;

		case 2://情况2：给了起始单元和结束单元，步长为1
			for(k=elem1;k<=elem2;k++)
			{
				i = QueryID(k, candidate, num_bay);
				if (i==-1)//如果找不到所给的单元，程序会跳过
					continue;
				candidate[i].plate.a = a;
				candidate[i].plate.b = b;
			}

			break;

		case 3://情况3：给了起始单元，结束单元以及步长pace
			for(k=elem1;k<=elem2;k=k+pace)
			{
				i = QueryID(k, candidate, num_bay);
				if (i==-1)//如果找不到所给的单元，程序会跳过
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



	for (j=0;j<num_lock;j++)//num_lock: 给定厚度的语句条目数
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
			candidate[i].plate.stack.NumLayers= 4*int(t/4);//厚度需要4的倍数，如果给定的不是4的倍数，则舍入

			break;

		case 2:
			for(k=elem1;k<=elem2;k++)
			{
				i = QueryID(k, candidate, num_bay);
				if (i==-1)
					continue;
				candidate[i].lock=1;
				candidate[i].plate.stack.NumLayers= 4*int(t/4);//厚度需要4的倍数
			}

			break;

		case 3:
			for(k=elem1;k<=elem2;k=k+pace)
			{
				i = QueryID(k, candidate, num_bay);
				if (i==-1)
					continue;
				candidate[i].lock=1;
				candidate[i].plate.stack.NumLayers= 4*int(t/4);//厚度需要4的倍数
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





	//数据检查，查看是否所有单元都定义了尺寸
	for(i=0;i<num_bay;i++)
	{
		if ((candidate[i].plate.a<1.0e-6)||(candidate[i].plate.b<1.0e-6))
		{
			printf("WARNING: Miss size property of element %d\n",candidate[i].FEM_ID );
			fprintf(fp_log, "WARNING: Miss size property of element %d\n",candidate[i].FEM_ID );
			candidate[i].plate.a = 600.;//如果没有给定尺寸，则设置为长635，宽150，并给出警告信息
			candidate[i].plate.b = 100.;
		}

		if ((candidate[i].lock==1)&&(candidate[i].plate.stack.NumLayers==0))//如果锁定了厚度，但是没有给厚度值，则给定厚度值为12，并给出警告信息
		{
			printf("WARNING: Miss thickness property of element %d\n",candidate[i].FEM_ID );
			fprintf(fp_log,  "WARNING: Miss thickness property of element %d\n",candidate[i].FEM_ID);
			candidate[i].plate.stack.NumLayers = min_layer;

		}
	}




	//读入连接信息，直至文件结束，连接信息的条目数不需要给出
	while(!feof(fp_in))
	{
		fscanf(fp_in, "%d %d",  &elem1, &elem2);
		i = QueryID(elem1, candidate, num_bay);
		j = QueryID(elem2, candidate, num_bay);
		ConnectVec[i+j*num_bay]=1;
		ConnectVec[j+i*num_bay]=1;
	}


	fclose(fp_in);//关闭输入文件






	///Pickup crital subcases for each bay and estimate thickness 猜测各个蒙皮的最严重工况，并估计第一轮厚度
	for (j=0;j<num_bay;j++)
	{
		candidate[j].PlateID  = j; //PlateID为代号，从0至num_bay-1
		candidate[j].DropCount =0; //DropCount记录丢层次数

		for (k=0;k<MAXLAYEROFPLY;k++)
			candidate[j].plate.stack.ply[k] = ply;

		if (candidate[j].lock ==0)
			//Estimate thickness,调用函数EstimateThicknessIso，假定各项同性，控制因素是受拉受压许用应变
			EstimateThicknessIso(&(candidate[j]), LoadM[j], num_subcases, 4500e-6, -3300e-6);

		candidate[j].plate.stack.NumLayers = MAX(candidate[j].plate.stack.NumLayers, min_layer);//至少需要12层
		if ((candidate[j].plate.stack.NumLayers%4)!=0)
			candidate[j].plate.stack.NumLayers = candidate[j].plate.stack.NumLayers + 2;//厚度需要4的倍数


	}


	//// Start optimization process
	for (RoundNo=1; RoundNo<=num_round; RoundNo++)
	{
		printf("Starting round %d\n", RoundNo);
		fprintf(fp_log, "Starting round %d\n", RoundNo);
		num_good = 0;

		if (RoundNo==(num_round))
		{
			num_outer_loop = num_outer_loop_inp; //最后一论计算时候需要精细一些，使用用户给定的循环数
			num_inner_loop = num_inner_loop_inp;
		}
		else
		{
			num_outer_loop = 100;//100是根据经验而来
			num_inner_loop = 100;//100是根据经验而来
		}


		if (RoundNo!=1)
		{
			for (j=0;j<num_bay;j++)
			{
				//optimal_solution是得到的合适的铺层方案，此处相当于把optimal_solution复制一份给candidate
				candidate[j] = optimal_solution[j];


				//update subcase here
				//上一轮铺层根据上一轮的严重工况计算得到，但是有可能在别的工况下裕度更低
				UpdateCriticalCase(&(candidate[j]), LoadM[j], num_subcases, epsi_al_t(candidate->plate), epsi_al_c(candidate->plate), PB_coef_c(candidate->plate) , PB_coef_s(candidate->plate), 0.);


				if ((candidate[j].ms.ms < reserved_ms)&&(candidate[j].lock==0)) // less than reserved ms， 如果裕度低于给定的保留裕度，则增厚4层
				{
					candidate[j].plate.stack.NumLayers = MIN(candidate[j].plate.stack.NumLayers + 4, MAXLAYEROFPLY);

				}


				if ((candidate[j].ms.ms > 1.0)&&(candidate[j].lock==0))//over conservative，裕度高于1认为过于保守，减去4层，但是为了避免反复减薄-》裕度不够-》加厚-》裕度过高的循环，所以使用DropCount控制
				{
					if(candidate[j].DropCount<2)
					{
						candidate[j].plate.stack.NumLayers = MAX(candidate[j].plate.stack.NumLayers - 4, min_layer);
						candidate[j].DropCount++;
					}

				}
			}
		}



		//根据厚度对板单元升序排列
		qsort (candidate, num_bay, sizeof(struct PlateCopula), CompareThick);
		max_layer = candidate[num_bay-1].plate.stack.NumLayers;//最大厚度


		EsimatePercentFirst(candidate); // 根据最小主应变的原则，估计最薄单元的铺层比例

		for (j=1;j<num_bay;j++)
		{
			EsimatePercent(candidate, j);//按照从薄到后的顺序估计各块蒙皮的铺层比例
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

		InitialPlyDeploy(candidate, num_bay);//给出铺层，按照从厚到薄的顺序，此处铺层满足连续性要求，但是不满足铺层设计原则



		//聚类分析
		PartitionCalc(clusters, &num_clusters, candidate, num_bay, max_layer, ConnectVec);




		for (i=0;i<num_bay;i++)
			optimal_solution[i]=candidate[i];


		num_safe_max = 0.;
		num_safe = 0.;

		for (count=0;count<num_outer_loop;count++)
		{
			
			
			
			
			if (num_safe > num_safe_max)//跟新最优解
			{
				num_safe_max = num_safe;
				for (i=0;i<num_bay;i++)
					optimal_solution[i]=candidate2[i];  //update the candidate
				printf("strength satidfied %f,   %d\n",  num_safe_max, count+1);
				fprintf(fp_log, "strength satidfied %f,   %d\n",  num_safe_max, count+1 );
			}

			else//按照一定概率更新最优解
			{
				dice = MAX(((double)(rand()%10000)/10000.),1.0e-10);
				double temp = pow(3.,(num_safe-num_safe_max)/0.005);

				if ((dice >= temp)||(count==0))//rej
				{
					for (i=0;i<num_bay;i++)
						candidate2[i]=optimal_solution[i];
				}
			}



			srand(time(NULL));//随机数发生器的种子
			best = num_bay*10;
			num_try =0;

			while ((best!=0)&&(num_try<num_inner_loop))
			{
				//
				
				
				num_try++;


				//随机选择两层（不包括最外面的cover_thick）
				lay1 = rand()%(max_layer/2-cover_thick)+cover_thick;
				lay2 = rand()%(max_layer/2)+cover_thick;




				if (lay2<max_layer/2)
				{
					opera =  Exchange2Plies(candidate2, num_bay, lay1, lay2);//交换lay1， lay2两层的位置
					quality=0;

					for (i=0;i<num_bay;i++)
					{
						quality = quality + plycheck1(candidate2[i].plate.stack) + plycheck5(candidate2[i].plate.stack)+ plycheck9(candidate2[i].plate.stack);
						//plycheck1: 不超过四层连续铺层
						//plycheck5: 使45度铺层使-45度铺层尽量成对出现，不会出现45/45/-45/-45
						//plycheck9: 两层连续铺层后的第三层不能与前两层正交
						//可以根据需要添加更多准则，铺层检查的代码在PlyCheck.h文件中
						// 铺层检测通过返回0，所以quality越小越好
						//铺层很难做到绝对满足，所有铺层最后作为一个权重集成到裕度中
					}
					if (quality<=best)
					{
						best=quality;
					}
					else
					{
						if (opera==1)
						{
							Exchange2Plies0(candidate2, num_bay, lay2, lay1); //新得到的铺层由于不满足铺层设计原则，再换回去。

						}
					}

				}

				else
				{

					opera =  SinkPly(candidate2, num_bay, lay1);//把lay1放在中面处
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

							RestorePly(candidate2, num_bay, lay1);//新得到的铺层由于不满足铺层设计原则，再换回去。

						}
					}
				}

				if (quality==0)//所有铺层都满足设计准则
					break;


			}
			printf("Obtained a candidate after %d tries, SQI %d\n",num_try, quality);
			fprintf(fp_log, "Obtained a candidate after %d tries, SQI %d\n",num_try, quality );


			Adjust45Pair(candidate2, num_bay, cover_thick);

			//铺层满足设计要求，判断裕度
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

				if ((switch_trim==1)&&(RoundNo==num_round))//微调，
				{
					for (j=2;j<num_clusters;j++)
					{
						opera = Reversal(candidate2, num_bay, clusters, num_clusters, j);//0换为90，90换为0，如果该层为45或者45，则返回0
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
				//微调结束





				//输出结果


				sprintf(OutFileName,"Solution_%d_%d.txt", RoundNo, ++num_good);

				fp =fopen(OutFileName,"w");
				if (fp == NULL)
				{
					printf("ERROR: Fail to open the output file.\n");
					fprintf(fp_log, "ERROR: Fail to open the output file.\n");
					return(0);
				}
				/*
				Plate No:         编号，根据读入顺序排列
				FEM ELEMENT ID:  FEM板单元号
				RANK:			  根据厚度的排序
				Subcase:		  最危险工况
				Layers:			  铺层数
				M.S:              最危险工况对应的裕度
				FailMode:         失效模式，1，3为应变许用值控制，2为稳定性控制
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







	//释放内存

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
