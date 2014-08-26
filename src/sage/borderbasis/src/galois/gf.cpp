/* gf.c
 * Yang Yuwang, Zhao Wei, Wang Lei
 * June, 2012

Fast Galois Field Arithmetic in C/C++
Copright (C) 2012 Yang Yuwang, Zhao Wei, Wang Lei

These files are free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

These files are distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

Yang Yuwang, Zhao Wei, Wang Lei
Department of Computer Science
Nanjing University of Science and Technology

yuwangyang@mail.njust.edu.cn
http://www.sensor608.com/gf.html
 */

#include "gf.h"
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>

//
GFType gfmul(state* s, GFType a, GFType b);
GFType gfdiv(state* s, GFType a, GFType b);
//

GFType prim_poly[13] =
{
/*	0 */	0x00000000,
/*  1 */    0x00000001,
/*  2 */    0x00000007,
/*  3 */    0x0000000b,
/*  4 */    0x00000013,
/*  5 */    0x00000025,
/*  6 */    0x00000043,
/*  7 */    0x00000089,
/*  8 */    0x00000187,
/*  9 */    0x00000211,
/* 10 */    0x00000409,
/* 11 */    0x00000805,
/* 12 */    0x00001053,
 };

state* gf_init(unsigned int m, unsigned int prim)// GF(2^m), primitive polymonial
{
    state* result = (state*)malloc(sizeof(state));

	int i=0,j=0;

	if (m > 12)	// the field size is supported from GF(2^1) to GF(2^12).
		return result;

	result->gFieldSize = 1<<m;

	if (0 == prim)
		prim = prim_poly[m];


	result->table_alpha = (GFType*)malloc(sizeof(GFType)*result->gFieldSize);
	result->table_index = (GFType*)malloc(sizeof(GFType)*result->gFieldSize);
	result->table_mul = (GFType**)malloc(sizeof(GFType*)*result->gFieldSize);
	result->table_div = (GFType**)malloc(sizeof(GFType*)*result->gFieldSize);
	for(i=0; i<result->gFieldSize; i++)
	{
		result->table_mul[i] = (GFType *)malloc(sizeof(GFType) * result->gFieldSize);
		result->table_div[i] = (GFType *)malloc(sizeof(GFType) * result->gFieldSize);
	}


	result->table_alpha[0]=1;
	result->table_index[0]=-1;

	for (i=1; i<result->gFieldSize; i++)
	{
		result->table_alpha[i] = result->table_alpha[i-1]<<1;
		if (result->table_alpha[i]>=result->gFieldSize)
		{
			result->table_alpha[i]^=prim;
		}

		result->table_index[result->table_alpha[i]]=i;
	}

	result->table_index[1]=0;

	// create the tables of mul and div
	for (i=0; i<result->gFieldSize; i++)
		for (j=0; j<result->gFieldSize; j++)
		{
			result->table_mul[i][j]=gfmul(result,i,j);
			result->table_div[i][j]=gfdiv(result,i,j);

		}

    return result;
}
void gf_uninit(state* s){
	int i = 0;

	free(s->table_alpha);
	free(s->table_index);

	for(i=0; i<s->gFieldSize; i++)
	{
		free(s->table_mul[i]);
		free(s->table_div[i]);
	}
	free(s->table_mul);
	free(s->table_div);
    free(s);

}

GFType gfmul(state* s,GFType a, GFType b)
{
	if (0==a || 0==b)
		return 0;

	return s->table_alpha[(s->table_index[a]+s->table_index[b])%(s->gFieldSize-1)];
}

GFType gfdiv(state* s,GFType a, GFType b)
{
	if (0==a || 0==b)
		return 0;

	return s->table_alpha[(s->table_index[a]-s->table_index[b]+(s->gFieldSize-1))%(s->gFieldSize-1)];
}


GFType gf_exp(state* s,GFType a, GFType n)
{
	return s->table_alpha[s->table_index[a]*n%(s->gFieldSize-1)];
}

