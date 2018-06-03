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

#include "PlyCheck.h"
#include "Composite.h"
#include <math.h>



/// Guideline 1:  no more than 4 successive plies
/// @param stack
int plycheck1(struct Stack stack)
{
    int i, num, quality;

	quality = 0;

	num=1;
	for(i=1; i<stack.NumLayers; i++)
    {
        if (stack.PlyDirection[i]==stack.PlyDirection[i-1])
        {
            num = num+1;
            if (num>4)
			{
               quality = quality + 1;
			   num=1;
			}
        }
        else
            num=1;
    }

    return(quality);
}




int plycheck2(struct Stack stack)
{
    int i;
    for(i=1; i<stack.NumLayers-1; i++)
    {
        if ((stack.PlyDirection[i]!=0)&&(stack.PlyDirection[i]!=90))
        {

            if ((stack.PlyDirection[i+1]!=-stack.PlyDirection[i])&&(stack.PlyDirection[i-1]!=-stack.PlyDirection[i]))
                return(2);

        }
    }

    return(0);
}



/// Guideline 3: no successive 90¡ãplies
/// @param stack
int plycheck3(struct Stack stack)
{
    int i, quality;
	quality =0;

    for(i=1; i<stack.NumLayers-1; i++)
    {
        if ((stack.PlyDirection[i]==stack.PlyDirection[i-1])&&(stack.PlyDirection[i]==90))
        {
			if (stack.PlyDirection[i+1]==90)
                quality++;
        }
    }

    return(quality);
}




int plycheck4(struct Stack stack)
{
    int i, num0, num45,num_45, num90;
	double frac0;
	double BelowBound, Upperbound;

	BelowBound = 0.08;
	Upperbound = 0.65;

    num0 = 0;
    num45 = 0;
    num_45 = 0;
    num90 = 0;

    for(i=0; i<stack.NumLayers; i++)
    {
        if (stack.PlyDirection[i]==0)
            num0++;
        else
            if ((stack.PlyDirection[i]==45))
            num45++;
        else
            if (stack.PlyDirection[i]==90)
            num90++;
        else
            if  ((stack.PlyDirection[i]==-45))
            num_45++;
    }
	frac0=((double)num0)/((double)stack.NumLayers);

	if ((frac0<BelowBound)||(frac0>Upperbound))
    	return(4);
	if ((((double)num45)/((double)stack.NumLayers)<BelowBound)|| (((double)num45)/((double)stack.NumLayers)>Upperbound))
    	return(45);
	if ((((double)num90)/((double)stack.NumLayers)<BelowBound)||(((double)num90)/((double)stack.NumLayers)>Upperbound))
    	return(90);
	if ((((double)num_45)/((double)stack.NumLayers)<BelowBound)||(((double)num_45)/((double)stack.NumLayers)>Upperbound))
    	return(-45);

return(0);
}




int plycheck5(struct Stack stack)
{
    int i;
	int num45,num_45,cum_ub, quality;
	quality =0;


	num45 = 0;
	num_45=0;
	cum_ub=0;

    for(i=0; i<(stack.NumLayers); i++)
    {
		if (stack.PlyDirection[i]==45)
        {
			num45++;
        }
       else if (stack.PlyDirection[i]==-45)
	   {
		   num_45++;
	   }

		cum_ub = cum_ub + num45 - num_45;

	   if ((cum_ub>1)||(cum_ub<-1))
	   {
		   quality ++;
		   cum_ub=0;
	   }

    }

	return(quality);
}



int plycheck6(struct Stack stack)
{
	int i, count;
	int num45,num_45;
	count=0;
	num45 =0;
	num_45=0;
	/*
	for(i=0; i<(stack.NumLayers); i++)
	{
		if (stack.PlyDirection[i]==45)
        {
			num45++;
        }
       if (stack.PlyDirection[i]==-45)
	   {
		   num_45++;
	   }
	}
	if (num45!=num_45)
		return(0);
	*/

    for(i=0; i<(stack.NumLayers); i++)
	{
		if (stack.PlyDirection[i]==45)
        {
			num45++;
        }
       if (stack.PlyDirection[i]==-45)
	   {
		   num_45++;
	   }

	}

	if (num45==num_45)
		return(0);
	else
		return(6);
}

int plycheck7(struct Stack stack)
{
	if ((stack.PlyDirection[0]==45)&&(stack.PlyDirection[1]==-45))
		return(0);
	else
		return(7);
}


int plycheck8(struct Stack stack)
{
	if ((stack.PlyDirection[2]!=0)&&(stack.PlyDirection[3]!=0))
		return(8);
	if ((stack.PlyDirection[2]!=90)&&(stack.PlyDirection[3]!=90))
		return(8);

	return(0);
}



int plycheck9(struct Stack stack)
{
	int i, quality;
	quality =0;

	for(i=2;i<stack.NumLayers;i++)
	{
		if(stack.PlyDirection[i-2]==stack.PlyDirection[i-1])
		{
			if(fabs((float)stack.PlyDirection[i]-stack.PlyDirection[i-1])==90)
				quality = quality + 1;
		}

		if(stack.PlyDirection[i]==stack.PlyDirection[i-1])
		{
			if(fabs((float)stack.PlyDirection[i]-stack.PlyDirection[i-1])==90)
				quality = quality + 1;
		}

	}

	return(quality);

}
