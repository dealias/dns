size(16cm,0);

string[] outnames={"lambda1","lambdar2","lambdar2rot","lambda2"};
real[] lambda={1,sqrt(2),sqrt(2),2};
pair[] R={(1,0),(1,0),exp(-pi*I/4),(1,0)};

int m=4; // really should be called mx, and a power of two.
int Ngrids=2; // the following three arrays only work up to ngrids=3

//bool drawring=true;

pen[] dotpen={black,blue,deepgreen};
filltype[] dotfill={Fill,NoFill,NoFill};
pen[] dotfillpen={black,invisible,invisible};
int n;
path g;
picture pic;
pen p, fillpen;
real L, s;
pair r;
filltype F;

usersetting();

int Invisible;

void drawdots() {
  for(int i=-n+1; i <n; ++i) {
    for(int j= i >= 0? 0 : 1; j <n; ++j) {
      pair a=L*r*(i,j);
      g=scale(s)*unitcircle;
      if((j >= Invisible) || (abs(i) >= Invisible)) {
	//label(pic,"("+(string) i +"," + (string) j+")",a,NE);
	filldraw(pic,shift(a)*g,fillpen,p);
      } else {
	filldraw(pic,shift(a)*g,fillpen,p+dashed);
      }
    }
  }
}


for(int i=0; i < lambda.length; ++i) {
  pic = new picture;
  size(pic,30cm,0);
  
  for(int G=0; G < Ngrids; ++G) {
    n=m;
    if(i==0) n=m*2^G;
    Invisible= G==0 ? 0 : floor(n/2);
    p=dotpen[G];
    fillpen=dotfillpen[G];
    F=dotfill[G];
    r=R[i]^G;
    L=lambda[i]^G;
    s=(G+1)*0.1;
    drawdots();
  }

  if(i==3) {
    // FIXME: this stuff really only applies to radix-4 grids

    int glast=Ngrids-1;
    real radlast=0.5;
    path semi=(1,0){N}..(0,1){W}..{S}(-1,0);
    for(int G=0; G < Ngrids; ++G) {
      int L=2^G;
      real istart=G==0 ? 0 : quotient(m,2);
      real rad = istart*L;
      int gmy=m;
      real kmy=L*gmy;
      real stoprad=G ==glast ? sqrt(2)*kmy-L : kmy;
      
      while(rad < stoprad) {
	path[] paths={scale(rad+L/2)*semi,scale(radlast)*semi};
	pen fpen=dotpen[G]+opacity(0.3);
	pen[] pens={fpen,fpen};
	draw(pic,paths,pens);
	
	draw(pic,scale(rad+L/2)*semi,dotpen[G]+dashed);
	radlast=rad+L/2;
	rad += L;
      }
      /*
      int start=G==0? 0 : quotient(n0,L);
      real stop=G==glast ? ceil((kmax +1/2)/L) : n0-1;

      for(int j=start; j < stop; j += 1) {
	draw(pic,scale(rad)*unitcircle,dotpen[G]+dashed);
	//draw(pic,scale(rad+w)*unitcircle,dotpen[G]+Dotted);
	write(rad);
	lastrad=rad;
      }
      */
    }
  }
  shipout(outnames[i],pic);

  
}

