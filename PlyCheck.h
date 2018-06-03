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

#ifndef _PLYCHECK
#define _PLYCHECK
#include "Composite.h"


/*
Reference:  NASA Contract NAS1-19347,  A SUMMARY AND REVIEW OF COMPOSITE LAMINATE DESIGN GUIDELINES
*/

/// Guideline 1:  no more than 4 successive plies
/// @param stack
/// @return return 0 if passes the test; 1 if fails the test 
int plycheck1(struct Stack stack);



/// Guideline 2: Laminates Are Required to Be Balanced.
/// @param stack
/// @return return 0 if passes the test; 2 if fails the test
int plycheck2(struct Stack stack);
 


/// Guideline 3: no successive 90¡ãplies
/// @param stack
/// @return return 0 if passes the test; else if fails the test
int plycheck3(struct Stack stack);


/// Guideline 4:    Laminates Will Be Fiber Dominated, Having at Least 10% of Their Plies in Each of the 0¡ã, ¡À45¡ã, and 90¡ã Directions.
///
/// Guideline 17:   The Maximum Percentage of Plies in Any Direction Will Be 60%.
/// @param stack
/// @return return 0 if passes the test; else if fails the test
int plycheck4(struct Stack stack);


///NOT AVAILABLE YET
int plycheck5(struct Stack stack);



///NOT AVAILABLE YET
/// @param stack
/// @return return 0 if passes the test; else if fails the test
int plycheck6(struct Stack stack);


/// The out most layer is in 45 degree
/// @param stack
/// @return return 0 if passes the test; else if fails the test
int plycheck7(struct Stack stack);




/// The outmost 4 layers cover all the four directions
/// @param stack
/// @return return 0 if passes the test; else if fails the test
int plycheck8(struct Stack stack);



int plycheck9(struct Stack stack);




#endif
