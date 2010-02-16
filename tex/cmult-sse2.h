/*
 * SSE2 complex multiplication routines derived from FFTW simd-sse2.h
 * by John C. Bowman (University of Alberta)
 *
 * Copyright (c) 2003, 2007-8 Matteo Frigo
 * Copyright (c) 2003, 2007-8 Massachusetts Institute of Technology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

typedef __m128d V;

#define VADD _mm_add_pd
#define VSUB _mm_sub_pd
#define VMUL _mm_mul_pd
#define VXOR _mm_xor_pd
#define SHUFPD _mm_shuffle_pd
#define UNPCKL _mm_unpacklo_pd
#define UNPCKH _mm_unpackhi_pd
#define STOREH _mm_storeh_pd
#define STOREL _mm_storel_pd

union uvec {
  unsigned u[4];
  V v;
};
  
const union uvec sse2_pm = {
  { 0x00000000,0x00000000,0x00000000,0x80000000 }
};

static inline V FLIP_RI(const V& x)
{
  return SHUFPD(x,x,1);
}

static inline V VCONJ(const V& x)
{
  return VXOR(sse2_pm.v,x);
}

static inline V VBYI(const V& x)
{
  return FLIP_RI(VCONJ(x));
}

static inline V VZMUL(const V& tx, const V& sr)
{
  return VADD(VMUL(sr,UNPCKL(tx,tx)),VMUL(UNPCKH(tx,tx),VBYI(sr)));
}

static inline V VZMUL(const V& t0, const V& t1, V sr)
{
  return VADD(VMUL(t0,sr),VMUL(t1,FLIP_RI(sr)));
}

static inline V LDA(const Complex *x)
{
  return *(const V *) x;
}

static inline void STA(Complex *x, const V& v)
{
  *(V *) x = v;
}

