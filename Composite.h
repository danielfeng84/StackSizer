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


/// ABD matrix ABD����
struct ABDM
{
	double	A11;/**< In-plane membrane stiffness, Unit: N/m ����ն�ϵ������λ�� N/m*/
	double	A12;/**< In-plane membrane stiffness, Unit: N/m ����ն�ϵ������λ�� N/m*/
	double	A16;/**< In-plane membrane stiffness, Unit: N/m ����ն�ϵ������λ�� N/m*/
	double	A22;/**< In-plane membrane stiffness, Unit: N/m ����ն�ϵ������λ�� N/m*/
	double	A26;/**< In-plane membrane stiffness, Unit: N/m ����ն�ϵ������λ�� N/m*/ 
	double	A66;/**< In-plane membrane stiffness, Unit: N/m ����ն�ϵ������λ�� N/m*/
	double	B11;/**< In-plane and bending coupled terms, Unit: N ��ϸն�ϵ������λ��N*/
	double	B12;/**< In-plane and bending coupled terms, Unit: N ��ϸն�ϵ������λ��N*/
	double	B16;/**< In-plane and bending coupled terms, Unit: N ��ϸն�ϵ������λ��N*/
	double	B22;/**< In-plane and bending coupled terms, Unit: N ��ϸն�ϵ������λ��N*/
	double	B26;/**< In-plane and bending coupled terms, Unit: N ��ϸն�ϵ������λ��N*/
	double	B66;/**< In-plane and bending coupled terms, Unit: N ��ϸն�ϵ������λ��N*/
	double	D11;/**< Flexural stiffness matrix entry�� Unit: N*m �����ն�ϵ������λ��N*m*/ 
	double	D12;/**< Flexural stiffness matrix entry�� Unit: N*m �����ն�ϵ������λ��N*m*/ 
	double	D16;/**< Flexural stiffness matrix entry�� Unit: N*m �����ն�ϵ������λ��N*m*/
	double	D22;/**< Flexural stiffness matrix entry�� Unit: N*m �����ն�ϵ������λ��N*m*/
	double	D26;/**< Flexural stiffness matrix entry�� Unit: N*m �����ն�ϵ������λ��N*m*/
	double	D66;/**< Flexural stiffness matrix entry�� Unit: N*m �����ն�ϵ������λ��N*m*/
};


/// Q reduced stiffness matrix of a single lamina �������ƽ��Ӧ��״̬�նȾ���
struct LaminaQ
{
	double		 Q11; /**< Unit: N/m/m ��λ���� */
	double		 Q12; /**< Unit: N/m/m ��λ���� */
	double		 Q16; /**< Unit: N/m/m ��λ���� */
	double		 Q22; /**< Unit: N/m/m ��λ���� */
	double		 Q26; /**< Unit: N/m/m ��λ���� */
	double		 Q66; /**< Unit: N/m/m ��λ���� */
};


/// Elastic modulus of a unidirectional lamina �����������Բ��ϸն�
struct CompsiteE
{
	double Ex; /**< Longitudinal elastic modulus. Unit: N/m/m. 1��������ģ������λ��N/m/m*/
	double Ey; /**< Transverse elastic modulus. Unit: N/m/m.   2��������ģ������λ��N/m/m*/
	double Gxy;/**< Shear modulus. Unit: N/m/m. ����ģ������λ��N/m/m */
	double MIUxy; /**< Poison's ratio ���ɱ�*/
};

/// 2D strain state 2άӦ��
struct TwoDStrain
{
	double EpsiX;/**< ��x����Ӧ��*/
	double EpsiY;/**< ��y����Ӧ�� */
	double EpsiXY;/**< ��Ӧ��*/
};


/// vector of curvature ��ϰ�����
struct Curvature
{
	double		KappaX;	 /**< D(w,x,x) ���ʣ���λ��1/m*/
	double		KappaY;	 /**< D(w,y,y) ���ʣ���λ��1/m*/
	double		KappaXY; /**< D(w,x,y) ���ʣ���λ��1/m*/
};



/// 2D stress state ƽ����Ӧ��״̬
struct TwoDStress
{
	double SigmaX;/**< In-plane stress in N/m/m.  ��һ������Ӧ���� ��λ��N/m/m*/
	double SigmaY;/**< In-plane stress in N/m/m.  �ڶ�������Ӧ���� ��λ��N/m/m*/
	double SigmaXY;/**< In-plane stress in N/m/m. ����Ӧ���� ��λ��N/m/m*/
};



/// Membrane Force ��ϰ�Ĥ����
struct PlateMembraneForce
{
	double		Nxx;  /**< In-plane normal force along the X-direction with the unit of N/m  ��ϰ����浥λ�����������������λ��N/m*/   
	double		Nyy;  /**< In-plane normal force along the Y-direction with the unit of N/m  ��ϰ����浥λ������������������λ��N/m*/  
	double		Nxy;  /**< In-plane shear force with the unit of N/m ��ϰ����浥λ�����ϼ�����������λ��N/m*/ 
};


/// Bending moment vector ��ϰ�����
struct PlateBendMoment
{
	double		Mxx;  /**< Distributed bending moment (moment per unit length). Unit: N  ��ϰ����浥λ��������������أ���λ��N*/
	double		Myy;  /**< Distributed bending moment (moment per unit length). Unit: N  ��ϰ����浥λ��������������أ���λ��N*/
	double		Mxy;  /**< Distributed torques with the unit of N  ��ϰ����浥λ��������������أ���λ��N*/
};


struct TransverseShear
{
	double Qx;
	double Qy;
};

/// Prototype of plate load, including three in-plane forces and three bending moments��ϰ�������������Ĥ��������
struct PlateLoad
{
	struct PlateMembraneForce InPlane; /**< In-plane force in N/m. ��Ĥ��������λ��N/m*/
	struct PlateBendMoment    OutofPlane; /**< Distributed bending moment in N ��ϰ����أ���λ��N*/
	struct TransverseShear    TShear;
	int LoadID;
	double ms;
};



/// Prototype of plate deformation, including three membrane strains and three curvatures  ��ϰ�Ӧ��
struct PlateDeform
{
	
	struct TwoDStrain InPlane;/**< mid-plane strain ��ϰ�����Ӧ��*/
	struct Curvature OutofPlane;/**< mid-plane curvature in 1/m ��ϰ��������ʣ���λ��1/m*/
	
};



/// Prototype of beam force 
struct BeamLoad
{
	double N;  /**< axial force, unit: N */
	double My; /**< bending moment, unit: N*m */
	double Mz; /**< bending moment, unit: N*m */
	double T; /**< torque, unit: N*m */
};


/// Prototype of plate deformation, including three membrane strains and three curvatures  ��ϰ�Ӧ��
struct BeamDeform
{
	
	double EpsiX;
	double KappaY;
	double KappaZ;
	double Rot; /**< angle of twist,unit: 1/m */
	
};



/// Prototype of unidirectional ply material model ���̲����
struct CompsiteMat
{
	double		Ex;    /**< Longitudinal elastic modulus ���̲�����ģ������λ����*/
	double		Ey;	   /**< Transverse elastic modulus ���̲����ģ������λ����*/
	double		Gxy;   /**< Shear modulus ���̲���иնȣ���λ���� */
	double		MIUxy;   /**< Poison's ratio ���̲������ɱ�*/
	double		thick; /**< Ply thickness ���̲��ȣ���λ��m*/
	double		rou;   /**< Density ���̲��ܶȣ���λ��kg/m/m/m*/

	double		Sxt;     /**< Tensional failure STRESS of an unidirectional tape ���̲���������ǿ�ȣ���λ����*/
	double		Sxc;    /**< Compressive failure STRESS of an unidirectional tape ���̲�����ѹ��ǿ�ȣ���λ����*/
	double		Syt;    /**< Tensional failure STRESS of an unidirectional tape ���̲��������ǿ�ȣ���λ����*/
	double		Syc;   /**< Compressive failure STRESS of an unidirectional tape ���̲����ѹ��ǿ�ȣ���λ����*/
	double		Ss;	   /**< Shear failure STRESS of an unidirectional tape �̲����ǿ�ȣ���λ����**/

};

/// stack of a laminate �̲�
struct Stack
{
	int      NumLayers;									/**< number of layers ����*/
	int      stat;
	int      PlyDirection[MAXLAYEROFPLY];               /**< ply directions, an array �̲�Ƕ�*/
	struct   CompsiteMat ply[MAXLAYEROFPLY];            /**< the ply �̲����*/
	struct   ABDM abdm;									/**< The ABD matrix ABD����*/	
};


/// Prototype of laminate ��ϰ�
struct Laminate
{
	double   a;     /**< laminate length in m ��ϰ峤�ȣ���λ��m*/
	double   b;     /**< laminate width in m ��ϰ��ȣ���λ��m*/
	double   r;		/**< radius of curvature in ��ϰ����ʰ뾶����λm*/
	int id;         /**< ��ϰ���*/
	struct Stack stack; /**< stack ��ϰ��̲�*/
};

/// Prototype of unidirectional ply material model ��ϰ�����״̬
struct PlateBuckleState
{
	struct  PlateLoad	  CrBuckleLoad;		  /**< Critical buckling force �ٽ������غ�*/
	struct  PlateDeform   CrBuckleStrain;	  /**< Critical buckling strain �ٽ�����Ӧ��*/
	int			m;						  /**< Number of half waves along the x direction ��x����벨��*/
	int			n;                        /**< Number of half waves along the y direction ��y����벨��*/
};




/// Lamina failure status
struct LaminaFail
{
	double ms; /**< Margin of safety */
	char   mode; /**< Failure mode */
};





/// Lamina failure status ��ϰ嵥�̲�ʧЧ״̬
struct LayerFail
{
	struct PlateLoad CrLoad;  /**< Critical failure load �̲�ʧЧʱ��ϰ��ٽ�Ӧ��״̬*/
	struct PlateDeform CrDeform; /**< Critical failure strain �̲�ʧЧʱ��ϰ��ٽ����״̬*/
	int    FailedLayerId;  /**< Failed layer No., starting from 0 ʧЧ��ţ����Ŵ�0��ʼ��*/
	char   FailureMode;   /**< Laminar failure mode: F-fibre failure, M-matrix failure ���̲�ʧЧģʽF-��άʧЧ��M-����ʧЧ*/
};



/// a point (y,z) on beam's cross-section ���������һ��
struct PlanePoint 
{
	double y;  /**< y y����*/
	double z; /**< z z����*/
};


/// Prototype of beam stiffness ������ն�
struct BeamSecStiff
{
	double		 A;     /**< Cross-sectional area, unit:m*m. ������������λ��m*m*/
	double		 EA;    /**< Axial stiffness,unit:N. ����նȣ���λ��N*/
	
	double		 EIyy;  /**< Bending stiffness along the y-axis, unit:N*m*m.  �ص�һ������(y��)�����նȣ���λ:N*m*m */
	double		 EIzz;  /**< Bending stiffness along the z-axis, unit:N*m*m.  �صڶ�������(z��)�����նȣ���λ:N*m*m*/
	double       EIxy;  /**< Coupling term between bending about y-axis and bending about z-axis, unit:N*m*m. */
	
	double		 GJ;   /**< Torsional stiffness, unit:N*m*m. */
	double		 WarpStiff;   /**< warp stiffness of the cross section, unit:N*m*m*m*m. */
	double       PolarRadiusGyration; /**< Polar Radius of Gyration of beam cross section, unit:m*m*m*m */
	
	
	
	struct PlanePoint Centroid; /**< centroid, a point(y,z) ��������*/

	struct PlanePoint ShearCentre; /**< shear centre, a point(y,z) ��������*/
};






/// Prototype of beam stiffness ������ն�
struct BeamSection
{
	
	char	  BeamType;
	struct    Laminate flange_up;    /**< upper flange ������Ե��*/
	struct    Laminate flange_low;	 /**< lower flange ������Ե��*/
	struct    Laminate web;			 /**< web ���츹��*/
	struct    BeamSecStiff stiff;     /**< beam cross-sectional stiffness ������ն�*/

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

	double Ky   /**<Translational Restraints spring stiffness, unit: N/m/m  ƽ�����ɸնȣ���λN/m/m */;
	double Kz;   /**<Translational Restraints spring stiffness, unit: N/m/m  ƽ�����ɸնȣ���λN/m/m */;
	double Kt;  /**<Rotational Restraints spring stiffness, unit: N/rad/m ת�����ɸնȣ���λN/rad/m */;
};




struct Beam
{
	double a;         /**< beam length, unit:m*/
	
	struct EndSupportCoef C;
	
	
	struct LateralSupport LateralSpring;

	struct BeamSection CrossSection;
};



/// Prototype stringer ����
struct StringerProfile
{
	char      StrType;               /**< a char, stringer type ������״*/
	struct    Laminate flange_up;    /**< upper flange ������Ե��*/
	struct    Laminate flange_low;	 /**< lower flange ������Ե��*/
	struct    Laminate web;			 /**< web ���츹��*/
	struct    Laminate skin;			/**< skin ��Ƥ*/
	struct    BeamSecStiff BeamSecStiff;		/**< stiffness of stringer �������ն�*/
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




/// prototype of stiffened panel ���ϲ��ϼӽ�ڰ�
struct CompositePanel
{
	double      b_bay;					/**< width of the panel bay �����࣬��λ��m*/
	double		a_bay;					/**< length of the panel bay �ڰ峤�ȣ���λ��m*/
	double      R;						/**< radius of curvature �ڰ����ʰ뾶����λ��m*/
	int			num_frame;				/**< number of frame, a_p/num_frame is the frame pitch */
	int			num_stringer;				/**< number of bays, b_p/num_bay is the stringer pitch */
	struct		StringerProfile SuperStringer;	/**< the stringer ����*/
	struct		FrameProfile frame;     /**< The frame ��*/
	
};



/// prototype of panel load �ӽ���غ�
struct PanelLoad
{
	double		Nxx;	     	/**< the uniaxial compressive/tensional load in N/m �����غɣ���λ��N/m*/
	double		Nxy;        /**< the shear load in N/m �����غɣ���λ��N/m*/
	double		EpsiX; /**< the uniaxial strain ����Ӧ��*/
	double		EpsiXY; /**< the uniaxial strain ����Ӧ��*/

};


/// prototype of stiffener strain in the postbuckling regiem
struct StiffnerStrain
{
	double eff;      /**< effective skin width ��Ƥ��Ч��ȣ���λ��m*/
	double dtf;      /**< diagonal tension factor �Խ�����ϵ��*/
	double intensity;     /**< */
	double epsi_ave; /**< average compressive strain */
	double epsi_idt_ave; /**< average compressive strain caused by diagonal tension */
	double epsi_idt_max; /**< maximum compressive strain caused by diagonal tension */
};





/// The in-plane two principal strains of a plate ����ƽ������Ӧ��
/// @param strain Ӧ��
/// @return principal strains ����ƽ����������Ӧ��
struct TwoDStrain PrincipalStrain(struct TwoDStrain strain);






/// The in-plane two principal stress of a plate ƽ����������Ӧ��
/// @param stress Ӧ��
/// @return principal stress p_stressp_stress.SigmaX: the first principal stress, p_stress.SigmaY: the second principal stress, p_stress.SigmaXY: N/A ����ƽ����������Ӧ��
struct TwoDStress PrincipalStress(struct TwoDStress stress);





/// The Transfer the stress tensor from the old coordinate system to a new coordinate systemƽ��Ӧ��״̬��Ӧ��������ϵת����ϵ
///
/// Last modified: 14-Oct-2012
/// @param stress stress tensor under the OLD coordinate system ������ϵ��Ӧ��
/// @param angle the ply angle ת��Ƕ�
/// @return stress tensor under the NEW coordinate system ����������ϵ��Ӧ��
struct TwoDStress StressTransfer(struct TwoDStress stress, double angle);




/// The Q matrix of a lamina ���̲����̲��������ϵ�¸նȾ���
///
/// Last modified: 14-Oct-2012
/// @param ply ���̲�
/// @return Q matrix in the ply's material coordinate system ���ص��̲�նȾ���
struct LaminaQ LocalLaminaQ(struct CompsiteMat ply);







/// The global Q matrix of a lamina with a given ply angle ���̲������������µĸնȾ���
///
/// Last modified: 14-Oct-2012
/// @param ply ���̲�
/// @param angle the ply angle �̲�Ƕ�
/// @return struct LaminaQ QG the Q matrix of the ply in the global coordinate system ���ص��̲������������µĸնȾ���
struct LaminaQ GlobalLaminaQ(struct CompsiteMat ply, double angle);





/// transfer local Q matrix to global Q matrix ���̲��ڲ�������ϵ�µĸնȾ����õ��̲�����������ϵ�µĸնȾ���
///
/// Last modified: 14-Oct-2012
/// @param QL the Q matrix of the ply in its material coordinate system ���̲��ڲ�������ϵ�µĸնȾ���
/// @param angle the ply angle �̲�Ƕ�
/// @return struct LaminaQ QG the Q matrix of the ply in the global coordinate system ���ص��̲�����������ϵ�µĸնȾ���
struct LaminaQ LaminaQL2G(struct LaminaQ QL, double angle);




/// calculate the stress of a lamina from its strain plane stress status applied ͨ���������Ӧ��״̬���㵥�����Ӧ��״̬
///
/// Last modified: 14-Oct-2012
/// @param Q  ���̲�նȾ���
/// @param strain ���̲�Ӧ��
/// @return stress ���̲�Ӧ��
struct TwoDStress Strain2Stress(struct LaminaQ Q, struct TwoDStrain strain);




///Transfer the strain from global coordinate system to local coordinate system �Ӳ�ϰ�ƽ����Ӧ����㵥�̲�Ӧ��
/// @param StrainG  strain tensor under the global coordinate system ��ϰ�Ӧ��
/// @param angle Ply angle ���̲�Ӧ��
/// @return the ply's strain under its matrical coordinate system
struct TwoDStrain StrainG2L(struct TwoDStrain StrainG, double angle);



///Calculate ABD matrix, passing the address as argument ����ABD���󣨵�ַ���ݣ�
/// @param stack stack �̲�
void ABDCalc(struct Stack *stack);


///Calculate ABD matrix, passing the value as argument  ����ABD����ֵ���ݣ�
/// @param stack stack �̲�
/// @return ABD matrix ����ABD����
struct ABDM ABD(struct Stack stack);



/// Calculate compliance matrix of a stack from stiffness matrix 
/// @param Stiffness
/// @return compliance matrix
struct ABDM ComplianceABD(struct ABDM Stiffness);




///Calculate equivalent Young's modulus of a composite plate
/// @param plate laminate
/// @return the laminate's Young's modulus
double ElastMod(struct Laminate plate);



///Calculate a laminate's elestic modulus along a specific direction ��ϰ���ĳһ����ն�
/// @param plate laminate
/// @param angle the given rotation angle 
/// @return the laminate's Young's modulus
struct CompsiteE PlateStiff(struct Laminate plate, double angle);





///Calculate thickness of a stack �����̲���
/// @param stack �̲�
/// @return stackthickness �����̲���
double StackThick(struct Stack stack);


///Calculate thickness of a laminate ��ϰ���
/// @param plate ��ϰ�
/// @return plate thicknes ���ز�ϰ���
double PlateThick(struct Laminate plate);




///Calculate thickness of a laminate ��ϰ���
/// @param plate ��ϰ�
/// @return plate thicknes ���ز�ϰ���
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


///Transfer the strain from one coordinate system to another coordinate system �Ӳ�ϰ�ƽ����Ӧ����㵥�̲�Ӧ��
/// @param strain  strain tensor under the global coordinate system ��ϰ�Ӧ��
/// @param angle Ply angle ���̲�Ӧ��
/// @return the ply's strain under its matrical coordinate system
struct TwoDStrain StrainTransfer(struct TwoDStrain strain, double angle);

#endif
