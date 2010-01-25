pair[][] operator *(pair[][] a, pair b)
{
  int n=a.length;
  pair[][] m=new pair[n][];
  for(int i=0; i < n; ++i)
    m[i]=a[i]*b;
  return m;

}

pair[][] operator *(pair b, pair[][] a)
{
  return a*b;
}

pair[][] operator cast(real[][] f) 
{
  pair[][] F=new pair[f.length][];
  for(int i=0; i < f.length; ++i)
    F[i]=(pair[]) f[i];
  return F;
}

real maxerror(pair[][] f, pair[][] g) {
  real e=0;
  for (int i=0; i < f.length; ++i) {
    for (int j=0; j < f[i].length; ++j) {
      if (e < abs(f[i][j]-g[i][j]))
	e=abs(f[i][j]-g[i][j]);
    }
  }
  return e;
}
