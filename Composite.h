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

/*
this header file is to define the macros and data structures
Last modified: 14-oct-2012
*/
#ifndef _COMPOSITELAMINATE
#define _COMPOSITELAMINATE
#define DROPEDPLY 999
#define MAXLAYEROFPLY 152   /*maximum layers*/
#define PI 3.141592653589793238462643383279502884
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#define	MAX(A, B)	((A) > (B) ? (A) : (B))


/// ABD matrix ABD矩阵
struct ABDM
{
	double	A11;/**< In-plane membrane stiffness, Unit: N/m 拉伸刚度系数，单位： N/m*/
	double	A12;/**< In-plane membrane stiffness, Unit: N/m 拉伸刚度系数，单位： N/m*/
	double	A16;/**< In-plane membrane stiffness, Unit: N/m 拉伸刚度系数，单位： N/m*/
	double	A22;/**< In-plane membrane stiffness, Unit: N/m 拉伸刚度系数，单位： N/m*/
	double	A26;/**< In-plane membrane stiffness, Unit: N/m 拉伸刚度系数，单位： N/m*/ 
	double	A66;/**< In-plane membrane stiffness, Unit: N/m 拉伸刚度系数，单位： N/m*/
	double	B11;/**< In-plane and bending coupled terms, Unit: N 耦合刚度系数，单位：N*/
	double	B12;/**< In-plane and bending coupled terms, Unit: N 耦合刚度系数，单位：N*/
	double	B16;/**< In-plane and bending coupled terms, Unit: N 耦合刚度系数，单位：N*/
	double	B22;/**< In-plane and bending coupled terms, Unit: N 耦合刚度系数，单位：N*/
	double	B26;/**< In-plane and bending coupled terms, Unit: N 耦合刚度系数，单位：N*/
	double	B66;/**< In-plane and bending coupled terms, Unit: N 耦合刚度系数，单位：N*/
	double	D11;/**< Flexural stiffness matrix entry， Unit: N*m 弯曲刚度系数，单位：N*m*/ 
	double	D12;/**< Flexural stiffness matrix entry， Unit: N*m 弯曲刚度系数，单位：N*m*/ 
	double	D16;/**< Flexural stiffness matrix entry， Unit: N*m 弯曲刚度系数，单位：N*m*/
	double	D22;/**< Flexural stiffness matrix entry， Unit: N*m 弯曲刚度系数，单位：N*m*/
	double	D26;/**< Flexural stiffness matrix entry， Unit: N*m 弯曲刚度系数，单位：N*m*/
	double	D66;/**< Flexural stiffness matrix entry， Unit: N*m 弯曲刚度系数，单位：N*m*/
};


/// Q reduced stiffness matrix of a single lamina 单层材料平面应力状态刚度矩阵
struct LaminaQ
{
	double		 Q11; /**< Unit: N/m/m 单位：帕 */
	double		 Q12; /**< Unit: N/m/m 单位：帕 */
	double		 Q16; /**< Unit: N/m/m 单位：帕 */
	double		 Q22; /**< Unit: N/m/m 单位：帕 */
	double		 Q26; /**< Unit: N/m/m 单位：帕 */
	double		 Q66; /**< Unit: N/m/m 单位：帕 */
};


/// Elastic modulus of a unidirectional lamina 正交各项异性材料刚度
struct CompsiteE
{
	double Ex; /**< Longitudinal elastic modulus. Unit: N/m/m. 1主方向弹性模量。单位：N/m/m*/
	double Ey; /**< Transverse elastic modulus. Unit: N/m/m.   2主方向弹性模量。单位：N/m/m*/
	double Gxy;/**< Shear modulus. Unit: N/m/m. 剪切模量，单位：N/m/m */
	double MIUxy; /**< Poison's ratio 泊松比*/
};

/// 2D strain state 2维应变
struct TwoDStrain
{
	double EpsiX;/**< 沿x方向应变*/
	double EpsiY;/**< 沿y方向应变 */
	double EpsiXY;/**< 剪应变*/
};


/// vector of curvature 层合板曲率
struct Curvature
{
	double		KappaX;	 /**< D(w,x,x) 曲率，单位：1/m*/
	double		KappaY;	 /**< D(w,y,y) 曲率，单位：1/m*/
	double		KappaXY; /**< D(w,x,y) 曲率，单位：1/m*/
};



/// 2D stress state 平面内应力状态
struct TwoDStress
{
	double SigmaX;/**< In-plane stress in N/m/m.  第一方向正应力， 单位：N/m/m*/
	double SigmaY;/**< In-plane stress in N/m/m.  第二方向正应力， 单位：N/m/m*/
	double SigmaXY;/**< In-plane stress in N/m/m. 剪切应力， 单位：N/m/m*/
};



/// Membrane Force 层合板膜内力
struct PlateMembraneForce
{
	double		Nxx;  /**< In-plane normal force along the X-direction with the unit of N/m  层合板横截面单位宽度上拉伸内力，单位：N/m*/   
	double		Nyy;  /**< In-plane normal force along the Y-direction with the unit of N/m  层合板横截面单位长度上拉伸内力，单位：N/m*/  
	double		Nxy;  /**< In-plane shear force with the unit of N/m 层合板横截面单位长度上剪切内力，单位：N/m*/ 
};


/// Bending moment vector 层合板力矩
struct PlateBendMoment
{
	double		Mxx;  /**< Distributed bending moment (moment per unit length). Unit: N  层合板横截面单位宽度上拉伸内力矩，单位：N*/
	double		Myy;  /**< Distributed bending moment (moment per unit length). Unit: N  层合板横截面单位宽度上拉伸内力矩，单位：N*/
	double		Mxy;  /**< Distributed torques with the unit of N  层合板横截面单位宽度上拉伸内力矩，单位：N*/
};


struct TransverseShear
{
	double Qx;
	double Qy;
};

/// Prototype of plate load, including three in-plane forces and three bending moments层合板内力，包括薄膜力和力矩
struct PlateLoad
{
	struct PlateMembraneForce InPlane; /**< In-plane force in N/m. 薄膜内力，单位：N/m*/
	struct PlateBendMoment    OutofPlane; /**< Distributed bending moment in N 层合板力矩，单位：N*/
	struct TransverseShear    TShear;
	int LoadID;
	double ms;
};



/// Prototype of plate deformation, including three membrane strains and three curvatures  层合板应变
struct PlateDeform
{
	
	struct TwoDStrain InPlane;/**< mid-plane strain 层合板中面应变*/
	struct Curvature OutofPlane;/**< mid-plane curvature in 1/m 层合板中面曲率，单位：1/m*/
	
};



/// Prototype of beam force 
struct BeamLoad
{
	double N;  /**< axial force, unit: N */
	double My; /**< bending moment, unit: N*m */
	double Mz; /**< bending moment, unit: N*m */
	double T; /**< torque, unit: N*m */
};


/// Prototype of plate deformation, including three membrane strains and three curvatures  层合板应变
struct BeamDeform
{
	
	double EpsiX;
	double KappaY;
	double KappaZ;
	double Rot; /**< angle of twist,unit: 1/m */
	
};



/// Prototype of unidirectional ply material model 单铺层材料
struct CompsiteMat
{
	double		Ex;    /**< Longitudinal elastic modulus 单铺层轴向模量，单位：帕*/
	double		Ey;	   /**< Transverse elastic modulus 单铺层横向模量，单位：帕*/
	double		Gxy;   /**< Shear modulus 单铺层剪切刚度，单位：帕 */
	double		MIUxy;   /**< Poison's ratio 单铺层纵向泊松比*/
	double		thick; /**< Ply thickness 单铺层厚度，单位：m*/
	double		rou;   /**< Density 单铺层密度，单位：kg/m/m/m*/

	double		Sxt;     /**< Tensional failure STRESS of an unidirectional tape 单铺层纵向拉伸强度，单位：帕*/
	double		Sxc;    /**< Compressive failure STRESS of an unidirectional tape 单铺层纵向压缩强度，单位：帕*/
	double		Syt;    /**< Tensional failure STRESS of an unidirectional tape 单铺层横向拉伸强度，单位：帕*/
	double		Syc;   /**< Compressive failure STRESS of an unidirectional tape 单铺层横向压缩强度，单位：帕*/
	double		Ss;	   /**< Shear failure STRESS of an unidirectional tape 铺层剪切强度，单位：帕**/

};

/// stack of a laminate 铺层
struct Stack
{
	int      NumLayers;									/**< number of layers 层数*/
	int      stat;
	int      PlyDirection[MAXLAYEROFPLY];               /**< ply directions, an array 铺层角度*/
	struct   CompsiteMat ply[MAXLAYEROFPLY];            /**< the ply 铺层材料*/
	struct   ABDM abdm;									/**< The ABD matrix ABD矩阵*/	
};


/// Prototype of laminate 层合板
struct Laminate
{
	double   a;     /**< laminate length in m 层合板长度，单位：m*/
	double   b;     /**< laminate width in m 层合板宽度，单位：m*/
	double   r;		/**< radius of curvature in 层合板曲率半径，单位m*/
	int id;         /**< 层合板编号*/
	struct Stack stack; /**< stack 层合板铺层*/
};

/// Prototype of unidirectional ply material model 层合板屈曲状态
struct PlateBuckleState
{
	struct  PlateLoad	  CrBuckleLoad;		  /**< Critical buckling force 临界屈曲载荷*/
	struct  PlateDeform   CrBuckleStrain;	  /**< Critical buckling strain 临界屈曲应变*/
	int			m;						  /**< Number of half waves along the x direction 沿x方向半波数*/
	int			n;                        /**< Number of half waves along the y direction 沿y方向半波数*/
};




/// Lamina failure status
struct LaminaFail
{
	double ms; /**< Margin of safety */
	char   mode; /**< Failure mode */
};





/// Lamina failure status 层合板单铺层失效状态
struct LayerFail
{
	struct PlateLoad CrLoad;  /**< Critical failure load 铺层失效时层合板临界应力状态*/
	struct PlateDeform CrDeform; /**< Critical failure strain 铺层失效时层合板临界变形状态*/
	int    FailedLayerId;  /**< Failed layer No., starting from 0 失效层号（层编号从0开始）*/
	char   FailureMode;   /**< Laminar failure mode: F-fibre failure, M-matrix failure 单铺层失效模式F-纤维失效；M-基体失效*/
};



/// a point (y,z) on beam's cross-section 梁横截面内一点
struct PlanePoint 
{
	double y;  /**< y y坐标*/
	double z; /**< z z坐标*/
};


/// Prototype of beam stiffness 梁截面刚度
struct BeamSecStiff
{
	double		 A;     /**< Cross-sectional area, unit:m*m. 横截面面积，单位：m*m*/
	double		 EA;    /**< Axial stiffness,unit:N. 轴向刚度，单位：N*/
	
	double		 EIyy;  /**< Bending stiffness along the y-axis, unit:N*m*m.  沿第一坐标轴(y轴)弯曲刚度，单位:N*m*m */
	double		 EIzz;  /**< Bending stiffness along the z-axis, unit:N*m*m.  沿第二坐标轴(z轴)弯曲刚度，单位:N*m*m*/
	double       EIxy;  /**< Coupling term between bending about y-axis and bending about z-axis, unit:N*m*m. */
	
	double		 GJ;   /**< Torsional stiffness, unit:N*m*m. */
	double		 WarpStiff;   /**< warp stiffness of the cross section, unit:N*m*m*m*m. */
	double       PolarRadiusGyration; /**< Polar Radius of Gyration of beam cross section, unit:m*m*m*m */
	
	
	
	struct PlanePoint Centroid; /**< centroid, a point(y,z) 弯曲中心*/

	struct PlanePoint ShearCentre; /**< shear centre, a point(y,z) 剪切中心*/
};






/// Prototype of beam stiffness 梁截面刚度
struct BeamSection
{
	
	char	  BeamType;
	struct    Laminate flange_up;    /**< upper flange 长桁上缘条*/
	struct    Laminate flange_low;	 /**< lower flange 长桁下缘条*/
	struct    Laminate web;			 /**< web 长桁腹板*/
	struct    BeamSecStiff stiff;     /**< beam cross-sectional stiffness 梁截面刚度*/

};

/// End FIXITY COEFFICIENTS 
struct EndSupportCoef
{
	double Cy; 		 /**< Fixty coefficient for buckling in the plane perpendicular to y-axis*/  
	double Cz;       /**< Fixty coefficient for buckling in the plane perpendicular to y-axis*/  
	double Ct;       /**< Fixty coefficient for buckling in the z-plane*/
};



struct LateralSupport
{
	struct PlanePoint fix;

	double Ky   /**<Translational Restraints spring stiffness, unit: N/m/m  平动弹簧刚度，单位N/m/m */;
	double Kz;   /**<Translational Restraints spring stiffness, unit: N/m/m  平动弹簧刚度，单位N/m/m */;
	double Kt;  /**<Rotational Restraints spring stiffness, unit: N/rad/m 转动弹簧刚度，单位N/rad/m */;
};




struct Beam
{
	double a;         /**< beam length, unit:m*/
	
	struct EndSupportCoef C;
	
	
	struct LateralSupport LateralSpring;

	struct BeamSection CrossSection;
};



/// Prototype stringer 长桁
struct StringerProfile
{
	char      StrType;               /**< a char, stringer type 长桁形状*/
	struct    Laminate flange_up;    /**< upper flange 长桁上缘条*/
	struct    Laminate flange_low;	 /**< lower flange 长桁下缘条*/
	struct    Laminate web;			 /**< web 长桁腹板*/
	struct    Laminate skin;			/**< skin 蒙皮*/
	struct    BeamSecStiff BeamSecStiff;		/**< stiffness of stringer 长桁截面刚度*/
};




/// frame type define
struct		FrameProfile
{
	char	  FrameType;
	struct    Laminate flange_up;   /**< upper flange */
	struct    Laminate flange_low;	/**< lower flange */
	struct    Laminate web;			/**< web */
	struct    Laminate skin;        /**< the skin */
	struct	  BeamSecStiff FrameStiff;   /**< the frame stiffness */
};




/// prototype of stiffened panel 复合材料加筋壁板
struct CompositePanel
{
	double      b_bay;					/**< width of the panel bay 长桁间距，单位：m*/
	double		a_bay;					/**< length of the panel bay 壁板长度，单位：m*/
	double      R;						/**< radius of curvature 壁板曲率半径，单位：m*/
	int			num_frame;				/**< number of frame, a_p/num_frame is the frame pitch */
	int			num_stringer;				/**< number of bays, b_p/num_bay is the stringer pitch */
	struct		StringerProfile SuperStringer;	/**< the stringer 长桁*/
	struct		FrameProfile frame;     /**< The frame 框*/
	
};



/// prototype of panel load 加筋板载荷
struct PanelLoad
{
	double		Nxx;	     	/**< the uniaxial compressive/tensional load in N/m 轴向载荷，单位：N/m*/
	double		Nxy;        /**< the shear load in N/m 剪切载荷，单位：N/m*/
	double		EpsiX; /**< the uniaxial strain 轴向应变*/
	double		EpsiXY; /**< the uniaxial strain 剪切应变*/

};


/// prototype of stiffener strain in the postbuckling regiem
struct StiffnerStrain
{
	double eff;      /**< effective skin width 蒙皮有效宽度，单位：m*/
	double dtf;      /**< diagonal tension factor 对角张拉系数*/
	double intensity;     /**< */
	double epsi_ave; /**< average compressive strain */
	double epsi_idt_ave; /**< average compressive strain caused by diagonal tension */
	double epsi_idt_max; /**< maximum compressive strain caused by diagonal tension */
};





/// The in-plane two principal strains of a plate 计算平面内主应变
/// @param strain 应变
/// @return principal strains 返回平面内两个主应变
struct TwoDStrain PrincipalStrain(struct TwoDStrain strain);






/// The in-plane two principal stress of a plate 平面内两个主应力
/// @param stress 应力
/// @return principal stress p_stressp_stress.SigmaX: the first principal stress, p_stress.SigmaY: the second principal stress, p_stress.SigmaXY: N/A 返回平面内两个主应力
struct TwoDStress PrincipalStress(struct TwoDStress stress);





/// The Transfer the stress tensor from the old coordinate system to a new coordinate system平面应力状态下应力随坐标系转换关系
///
/// Last modified: 14-Oct-2012
/// @param stress stress tensor under the OLD coordinate system 旧坐标系下应力
/// @param angle the ply angle 转轴角度
/// @return stress tensor under the NEW coordinate system 返回新坐标系下应力
struct TwoDStress StressTransfer(struct TwoDStress stress, double angle);




/// The Q matrix of a lamina 单铺层在铺层材料坐标系下刚度矩阵
///
/// Last modified: 14-Oct-2012
/// @param ply 单铺层
/// @return Q matrix in the ply's material coordinate system 返回单铺层刚度矩阵
struct LaminaQ LocalLaminaQ(struct CompsiteMat ply);







/// The global Q matrix of a lamina with a given ply angle 单铺层在总体坐标下的刚度矩阵
///
/// Last modified: 14-Oct-2012
/// @param ply 单铺层
/// @param angle the ply angle 铺层角度
/// @return struct LaminaQ QG the Q matrix of the ply in the global coordinate system 返回单铺层在总体坐标下的刚度矩阵
struct LaminaQ GlobalLaminaQ(struct CompsiteMat ply, double angle);





/// transfer local Q matrix to global Q matrix 单铺层在材料坐标系下的刚度矩阵获得单铺层在总体坐标系下的刚度矩阵
///
/// Last modified: 14-Oct-2012
/// @param QL the Q matrix of the ply in its material coordinate system 单铺层在材料坐标系下的刚度矩阵
/// @param angle the ply angle 铺层角度
/// @return struct LaminaQ QG the Q matrix of the ply in the global coordinate system 返回单铺层在总体坐标系下的刚度矩阵
struct LaminaQ LaminaQL2G(struct LaminaQ QL, double angle);




/// calculate the stress of a lamina from its strain plane stress status applied 通过单层材料应变状态计算单层材料应力状态
///
/// Last modified: 14-Oct-2012
/// @param Q  单铺层刚度矩阵
/// @param strain 单铺层应变
/// @return stress 单铺层应力
struct TwoDStress Strain2Stress(struct LaminaQ Q, struct TwoDStrain strain);




///Transfer the strain from global coordinate system to local coordinate system 从层合板平面内应变计算单铺层应变
/// @param StrainG  strain tensor under the global coordinate system 层合板应变
/// @param angle Ply angle 单铺层应变
/// @return the ply's strain under its matrical coordinate system
struct TwoDStrain StrainG2L(struct TwoDStrain StrainG, double angle);



///Calculate ABD matrix, passing the address as argument 计算ABD矩阵（地址传递）
/// @param stack stack 铺层
void ABDCalc(struct Stack *stack);


///Calculate ABD matrix, passing the value as argument  计算ABD矩阵（值传递）
/// @param stack stack 铺层
/// @return ABD matrix 返回ABD矩阵
struct ABDM ABD(struct Stack stack);



/// Calculate compliance matrix of a stack from stiffness matrix 
/// @param Stiffness
/// @return compliance matrix
struct ABDM ComplianceABD(struct ABDM Stiffness);




///Calculate equivalent Young's modulus of a composite plate
/// @param plate laminate
/// @return the laminate's Young's modulus
double ElastMod(struct Laminate plate);



///Calculate a laminate's elestic modulus along a specific direction 层合板沿某一方向刚度
/// @param plate laminate
/// @param angle the given rotation angle 
/// @return the laminate's Young's modulus
struct CompsiteE PlateStiff(struct Laminate plate, double angle);





///Calculate thickness of a stack 计算铺层厚度
/// @param stack 铺层
/// @return stackthickness 返回铺层厚度
double StackThick(struct Stack stack);


///Calculate thickness of a laminate 层合板厚度
/// @param plate 层合板
/// @return plate thicknes 返回层合板厚度
double PlateThick(struct Laminate plate);




///Calculate thickness of a laminate 层合板厚度
/// @param plate 层合板
/// @return plate thicknes 返回层合板厚度
double SecArea(struct Laminate plate);



///hightness of the lay_id+1-th layer
/// @param stack
/// @param lay_id layer ID, starting from 0
/// @return layer height
double LayerHeight(struct Stack stack, int lay_id);



///displacement of a composite plate subjected to load
/// @param plate
/// @param pl palte load
/// @return plate displacement vector
struct PlateDeform Disp(struct Laminate plate, struct PlateLoad pl);


///calculate the plate load from a given plate displacement
/// @param plate
/// @param disp palte displacement
/// @return plate load
struct PlateLoad Load(struct Laminate plate, struct PlateDeform disp);




/// Calculate the elastic modulus form the Q matrix of a lamina
/// @param Q
/// @return Modules
struct CompsiteE Q2E(struct LaminaQ Q);


///Transfer the strain from one coordinate system to another coordinate system 从层合板平面内应变计算单铺层应变
/// @param strain  strain tensor under the global coordinate system 层合板应变
/// @param angle Ply angle 单铺层应变
/// @return the ply's strain under its matrical coordinate system
struct TwoDStrain StrainTransfer(struct TwoDStrain strain, double angle);

#endif
