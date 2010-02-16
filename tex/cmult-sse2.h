/* SSE2 complex multiplication routines
   Copyright (C) 2010 John C. Bowman, University of Alberta

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA. */

#ifndef __cmult_sse2_h__
#define __cmult_sse2_h__ 1

#ifdef __SSE2__      

#include <emmintrin.h>

typedef __m128d V;

#define UNPACKL _mm_unpacklo_pd
#define UNPACKH _mm_unpackhi_pd

union uvec {
  unsigned u[4];
  V v;
};
  
// TODO: move to cc file.
const union uvec sse2_pm = {
  { 0x00000000,0x00000000,0x00000000,0x80000000 }
};

static inline V FLIP(const V& z)
{
  return _mm_shuffle_pd(z,z,1);
}

static inline V CONJ(const V& z)
{
  return _mm_xor_pd(sse2_pm.v,z);
}

// Return the complex product of iz.
static inline V ZMULTI(const V& z)
{
  return FLIP(CONJ(z));
}

// Return the complex product of z and w.
static inline V ZMULT(const V& z, const V& w)
{
  return w*UNPACKL(z,z)+UNPACKH(z,z)*ZMULTI(w);
}

// Return the complex product of z and I*w.
static inline V ZMULTI(const V& z, const V& w)
{
  return ZMULTI(w)*UNPACKL(z,z)-UNPACKH(z,z)*w;
}

static inline V ZMULT(const V& t0, const V& t1, V w)
{
  return t0*w+t1*FLIP(w);
}

static inline V LOAD(const Complex *z)
{
  return *(const V *) z;
}

static inline void STORE(Complex *z, const V& v)
{
  *(V *) z = v;
}

#endif

#endif
