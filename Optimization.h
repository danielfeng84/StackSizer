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

#ifndef _OPTIMIZATIONFUNCTIONS
#define _OPTIMIZATIONFUNCTIONS
#include "Composite.h"


struct MarginOfSafety
{
	double ms;
	int subcase;
	int mode;
};


struct PlateCopula
{
	struct Laminate  plate;     //蒙皮板 
	struct PlateLoad pl;        //最危险工况的内力
	int half_nums[4];			//0,45,-45,90度铺层数目
	int dirc[4];
	int DropCount;
	
	// first digit: position;  second digit: 0,null 1,2 antive; third digit: group
	int posi0[MAXLAYEROFPLY][4];
	int posi45[MAXLAYEROFPLY][4];
	int posi_45[MAXLAYEROFPLY][4];
	int posi90[MAXLAYEROFPLY][4];


	int PlateID;				//
	int FEM_ID;                //蒙皮单元对应的GFEM单元号
	int lock;                   //1代表蒙皮厚度由用户给定，程序将不考虑修改厚度
	
	struct MarginOfSafety ms;   //
};







struct Clusters
{
	int num_clusters;
	int cluster_size[5000];
	int* container[5000];
};

struct ud
{ 
	int angle;
	int order;
};



struct ABDM IsoABD(struct LaminaQ Q, double t);







/** \brief
 *
 * \param
 * \param
 * \return
 *
 */

struct Stack SeachMaxOEF(struct CompsiteMat ply, int num_0, int num_45, int num_90);





int *MaxStiffEstim(int num_ply_web);



void TrimWeb(struct CompositePanel *panel, double TargetStrain);
void TrimFlange(struct CompositePanel *panel, double TargetStrain);


/** \brief generate a uniformly distributed random interger with interval [0, num]
 *
 * \param num an integer
 * \param
 * \return
 *
 */
inline int RandU(int num);

inline int RandWander(int num);

inline double RandN(double std);

int Binomial2(double Pr);

/** \brief
 *
 * \param
 * \param
 * \return
 *
 */
void AttribParam(struct CompositePanel *panel, struct CompsiteMat ply, struct PanelParam param);

void InitialPlyDeploy(struct PlateCopula plates[], int num_plates);
void EsimatePercent(struct PlateCopula* plates, int id);
void EsimatePercentFirst(struct PlateCopula* plate);

void EstimateThicknessIso(struct PlateCopula* plates, struct PlateLoad *pl, int num_subcase, double epsi_al_t, double epsi_al_c);

ud Find(const struct PlateCopula* plates, int layID);

void RestorePly(struct PlateCopula plates[], int num_plates, int lay1);
struct Laminate LamGen(struct  PlateCopula plate);
int SinkPly(struct PlateCopula plates[], int num_plates, int lay1);
int Exchange2Plies(struct PlateCopula plates[], int num_plates, int lay1, int lay2);
void Exchange2Plies0(struct PlateCopula plates[], int num_plates, int lay1, int lay2);
int Reversal(struct PlateCopula plates[], int num_plates, int *clusters, int num_clusters, int cluster_id);
long long powint (int a,int b);
int CompareID (const void * a, const void * b);
struct MarginOfSafety MarginSafety(struct Laminate plate, struct PlateLoad pl, double t_allow, double c_allow, double PB_coef_c, double PB_coef_s, double Smear_Ratio);
int IsNeighbor(int num_plates,  int *ConnectVec, int id1, int id2);
void PartitionCalc(int *clusters, int *num_clusters, struct PlateCopula *plates, int num_plates, int max_layer, int *ConnectVec);
void MergeTwoClusters(struct Clusters *cluster0, int cluster_i, int cluster_j);
void UpdateCriticalCase(struct PlateCopula* plates, struct PlateLoad *pl, int num_subcase, double epsi_al_t, double epsi_al_c, double PB_coef_c, double PB_coef_s, double Smear_Ratio);
int  PlyBundleList(int list[][4], int num_layer);
void initialize_and_permute(int index[], int n);
void Adjust45Pair(struct PlateCopula plates[], int num_plates, int cover_thick);
#endif
