void multadvection2(double **F, unsigned int m,
                    const unsigned int indexsize,
                    const unsigned int *index,
                    unsigned int r, unsigned int threads)
{
  double* F0=F[0];
  double* F1=F[1];
  
  for(unsigned int j=0; j < m; ++j) {
    double u=F0[j];
    double v=F1[j];
    F0[j]=v*v-u*u;
    F1[j]=u*v;
  }
}

