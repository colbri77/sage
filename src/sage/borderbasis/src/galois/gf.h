/* gf.h
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

#ifndef _GF_H
#define _GF_H

typedef unsigned int GFType;

struct state {
    int gFieldSize;
    GFType* table_alpha;
    GFType* table_index;
    GFType** table_mul;
    GFType** table_div;
};

state* gf_init(unsigned int m, unsigned int prim);
void gf_uninit(state* s);
void gf_print(state* s);

#define  gf_alpha(s,n)   (s->table_alpha[n])
#define  gf_index(s,n)   (s->table_index[n])

#define  gf_add(s,a,b)	(a^b)
#define  gf_sub(s,a,b)	(a^b)

#define  gf_mul(s,a,b)	(s->table_mul[a][b])
#define  gf_div(s,a,b)	(s->table_div[a][b])

#define  gf_neg(s,a)		(a)
#define  gf_inv(s,a)		(s->table_div[1][a])
GFType gf_exp(state* s,GFType a, GFType n);


#endif
