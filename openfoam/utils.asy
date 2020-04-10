string run;

real[][][] read(real i, string file) {
  file in=input(run+"/"+string(i)+"/"+file).line();
  while(!eof(in)) {
    string[] s=split(in);
    if(s.length > 0 && s[0] == "internalField") break;
  }
  int n2=in;
  int n=(int) sqrt(n2);
  assert(n*n == n2);

  string s=in;
  assert(s == "(");

  in.line(false);
  triple[][] z=in.dimension(n,n);
  
  real[][] u=sequence(new real[](int i) {
      return
        sequence(new real(int j) {return z[j][i].x;},n);},n);
  real[][] v=sequence(new real[](int i) {
      return
        sequence(new real(int j) {return z[j][i].y;},n);},n);

  u.cyclic=true;
  v.cyclic=true;

  for(int i=0; i < n; ++i)
    u[i].cyclic=v[i].cyclic=true;

  return new real[][][]{u,v};
}

