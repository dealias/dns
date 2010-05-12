size(8cm,0);

int f=3;
int c=5;
void cross(pair z=(0,0), real w=0.1)
{
  draw((z+w*(1,1))--(z-w*(1,1)),blue);
  draw((z+w*(1,-1))--(z-w*(1,-1)),blue);
}

for (int i=-f+1; i < f; ++i ) {
  for (int j=-f+1; j < f; ++j ) {
    if ((i,j) != (0,0)) 
      dot((i,j),red);
  }
}




for (int i=-c+1; i < c; ++i ) {
  for (int j=-c+1; j < c; ++j ) {
    if (i % 2 == 0 && j % 2 == 1)
      cross((i,j));
    if (i % 2 == 1 && j % 2 == 0)
      cross((i,j));
  }
}
