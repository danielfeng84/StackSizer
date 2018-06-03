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

#ifndef _BUCKLINGFUNCTIONS
#define _BUCKLINGFUNCTIONS
#include "Composite.h"












/// Unidirectional Compressive Buckling load of a plate, S stands for simply supported C stands for clamp
/// @param plate: laminate
/// @return BucklingLoad.BuckForce: the critical buckling load in N/m
///			BucklingLoad.m: N/A
///			BucklingLoad.n: N/A
struct PlateBuckleState PlateBuckling_UC_F_SCSC_Ortho(struct Laminate plate);



/// Unidirectional Compressive Buckling load of a plate, S stands for simply supported C stands for clamp
/// @param plate:		laminate
/// @return BucklingLoad.BuckForce: the critical buckling load in N/m
///			BucklingLoad.m: number of half waves along the loading direction
///			BucklingLoad.n: 1
struct PlateBuckleState PlateBuckling_UC_F_SSSS_Ortho(struct Laminate plate);



/// Unidirectional Compressive Buckling load of a plate, S stands for simply supported C stands for clamp
/// @param plate:		laminate
/// @return BucklingLoad.BuckForce: the critical buckling load in N/m
///			BucklingLoad.m: number of half waves along the loading direction
///			BucklingLoad.n: 1
struct PlateBuckleState PlateBuckling_UC_F_CCCC_Ortho(struct Laminate plate);



/// Unidirectional Compressive Buckling load of a Curved Plate with four sides simply supported
/// @param plate laminate
/// @return BucklingLoad.BuckForce: the critical buckling load in N/m
///			BucklingLoad.m: number of half waves along the loading direction
///			BucklingLoad.n: 1
struct PlateBuckleState PlateBuckling_UC_C_SSSS_Ortho(struct Laminate plate);




/// Unidirectional Compressive Buckling load of a plate, S stands for simply supported C stands for clamp
/// @param plate laminate
/// @param k:		End fix coefficient
/// @return BucklingLoad.BuckForce: the critical buckling load in N/m
///			BucklingLoad.m: number of half waves along the loading direction
///			BucklingLoad.n: number of half waves along the unloaded direction
struct PlateBuckleState PlateBuckling_BC_F_SSSS_Ortho(struct Laminate plate, double k);




struct PlateBuckleState PlateBuckling_UC_F_SBSB_Ortho(struct Laminate plate, double EA, double EI, double GJ);

/// Unidirectional Compressive Buckling load of a plate. The two loading edges are simply supported while one unloaded edge is simplily supported and the other edge free
/// @param plate laminate
/// @return BucklingLoad.BuckForce: the critical buckling load in N/m
///			BucklingLoad.m: number of half waves along the loading direction
///			BucklingLoad.n: 1
struct PlateBuckleState PlateBuckling_UC_F_SSSF_Ortho(struct Laminate plate);


/// Unidirectional Compressive Buckling load of a plate. The two loading edges are simply supported while one unloaded edge is simplily supported and the other edge free
/// @param plate laminate
/// @param C0:		End fix coefficient
/// @return BucklingLoad.BuckForce: the critical buckling load in N/m
///			BucklingLoad.m: number of half waves along the loading direction
///			BucklingLoad.n: 1
struct PlateBuckleState PlateBuckling_UC_F_SESF_Ortho(struct Laminate plate, double C0);



/// Unidirectional Compressive Buckling load of a plate. The two loading edges are simply supported while one unloaded edge is simplily supported and the other edge free
/// @param plate laminate
/// @param C0:		End fix coefficient
/// @return BucklingLoad.BuckForce: the critical buckling load in N/m
///			BucklingLoad.m: number of half waves along the loading direction
///			BucklingLoad.n: 1
struct PlateBuckleState PlateBuckling_UC_F_SESE_Ortho(struct Laminate plate, double C0);



/// Shear buckling load of a plate with four edges clamped. This function uses the empirical formula recommend by NASA
///
/// Last modified: 14-Oct-2012
/// @param plate
/// @return struct PlateBuckleState BucklingLoad
///
///			BucklingLoad.m: N/A
///			BucklingLoad.n: N/A
struct PlateBuckleState PlateBuckling_S_F_CCCC_Ortho(struct Laminate plate);




/// Shear buckling load of a plate with four edges simply supported. This function uses the empirical formula recommend by NASA
///
/// Last modified: 14-Oct-2012
/// @param plate
/// @return struct PlateBuckleState BucklingLoad
///
///			BucklingLoad.m: N/A
///			BucklingLoad.n: N/A
struct PlateBuckleState PlateBuckling_S_F_SSSS_Ortho(struct Laminate plate);




/// Shear buckling load. This function uses the empirical formula recommend by NASA
///
/// Last modified: 14-Oct-2012. THIS FUNCTION IS INCORRECT
/// @param plate
/// @return struct PlateBuckleState BucklingLoad
///
///			BucklingLoad.m: N/A
///			BucklingLoad.n: N/A
struct PlateBuckleState PlateBuckling_S_F_SCSC_Ortho(struct Laminate plate);



/// Shear buckling load. This function uses the Galekin method
///
/// Last modified: 14-Oct-2012
/// @param plate
/// @return struct PlateBuckleState BucklingLoad
///
///			BucklingLoad.BuckForce: the critical shear buckling load in N/m
struct PlateBuckleState ShearBuckling_SSSS_Ortho(struct Laminate plate);



/// Critical buckling load of an anisotropic plate
///
/// Last modified: 10-Dec-2012
/// @param plate
/// @param pl plate load prescribing the the load direction
/// @return Output: struct PlateBuckleState BucklingLoad
///
///			BucklingLoad.BuckForce: the critical shear buckling load in N/m
struct PlateBuckleState PlateBuckling_NS_F_SSSS_Aniso(struct Laminate plate, struct PlateLoad pl);






/// Compressive and shear buckling load
///
/// Last modified: 14-Oct-2012
/// @param Nxcrit critical compressive buckling load in N/m
/// @param Nxycrit: critical shear buckling load in N/m
/// @param NSAngle: tan(NSAngle)=Nxy/Nx
struct PlateBuckleState NSBuckling(double Nxcrit, double Nxycrit, double NSAngle);






#endif
