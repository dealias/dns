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
     { 0x00000000, 0x00000000, 0x00000000, 0x80000000 }
};

static inline V FLIP_RI(V x)
{
     return SHUFPD(x, x, 1);
}

static inline V VCONJ(V x)
{
     return VXOR(sse2_pm.v, x);
}

static inline V VBYI(V x)
{
     x = VCONJ(x);
     x = FLIP_RI(x);
     return x;
}

static inline V VZMUL(V tx, V sr)
{
     V tr = UNPCKL(tx, tx);
     V ti = UNPCKH(tx, tx);
     tr = VMUL(sr, tr);
     sr = VBYI(sr);
     return VADD(tr, VMUL(ti, sr));
}

static inline V LDA(const Complex *x)
{
  return *(const V *)x;
}

static inline void STA(Complex *x, V v)
{
  *(V *)x = v;
}

