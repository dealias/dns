size(32cm,0);

// asy -f pdf -u "m=8" -u "Ngrids=3" lambda.asy

string[] outnames={"lambda1","lambdar2","lambdar2rot","lambda2"};
real[] lambda={1,sqrt(2),sqrt(2),2};
pair[] R={(1,0),(1,0),exp(pi*I/4),(1,0)};
int[] radices={1,2,2,4};

int m=4; // really should be called mx, and a power of two.
int Ngrids=2; // the following three arrays only work up to ngrids=3

//bool drawring=true;

pen[] dotpen={black,blue,deepgreen};
filltype[] dotfill={Fill,NoFill,NoFill};
pen[] dotfillpen={black,invisible,invisible};
bool circular=false;
int n,G;
int radix;
path g;
picture pic;
pen p, fillpen;
real L, s;
pair r;
filltype F;
int nold=0;
usersetting();

real size=20cm;
if(m==8)
  size=30cm;
if(m==16)
  size=70cm;


int Invisible;

bool visible(int i, int j, int Invisible) {
  if(circular) {
    real a=max((nold-1.0),0);
    real minr2=a*a/radix;
    int maxr2=(n-1)*(n-1);
    int r2=i*i+j*j;
    if(Invisible ==0)
      return (r2 <= maxr2);
    if(radix==2) {
      return((r2 >= minr2) && (r2 < maxr2));
    }
    if(radix==4) {
      return((r2 >= minr2) && (r2 < maxr2));
    }
  } else {
    if(radix==1 || radix == 4)
      if((abs(j) >= Invisible) || (abs(i) >= Invisible))
	return true;
    if(radix==2)
      if(abs(i) + abs(j)  >= Invisible)
	return true;
  }
  return false;
}

void drawdots() {  
  for(int i=-n+1; i <n; ++i) {
    //for(int i=0; i <n; ++i) {
    for(int j= i >= 0? 0 : 1; j <n; ++j) {
      //for(int j=-n+1; j <n; ++j) {
      pair a=L*r*(i,j);
      g=scale(s)*unitcircle;
      if(visible(i,j,Invisible)) {
	filldraw(pic,shift(a)*g,fillpen,p);
	//label(pic,"("+(string) i +"," + (string) j+")",a,NE);
	//label(pic,"("+(string) i +"," + (string) j+")",a,NE,dotfillpen[G]);
	label(pic,(string) ((radix^G)*(i*i+j*j)),a,NE,p);
      } else {
	//label(pic,"("+(string) i +"," + (string) j+")",a,NE,dotfillpen[G]);
	filldraw(pic,shift(a)*g,fillpen,p+dashed);
	//label(pic,(string) ((radix^G)*(i*i+j*j)),a,NE);
      }
    }
  }
}


for(int i=0; i < lambda.length; ++i) {
  pic = new picture;
  size(pic,size,0);
  radix=radices[i];
  nold=0;
  for(G=0; G < Ngrids; ++G) {
    n=m;
    if(radix != 4) n=m*2^G;
    Invisible= G==0 ? 0 : floor(n/2);
    if(radix==2) {
      Invisible= G==0 ? 0 : floor(n/2);
    }
    p=dotpen[G];
    fillpen=dotfillpen[G];
    F=dotfill[G];
    r=R[i]^G;
    L=lambda[i]^G;
    s=(G+1)*0.1;
    drawdots();
    
    if(radix==2) {// centers picture correctly
      dot(pic,(n*L,0),invisible);
      dot(pic,(-n*L,0),invisible); 
      dot(pic,(0,n*L),invisible);
      dot(pic,(0,-n*L),invisible); 
    }
    nold=n;
  }
  
  if(radix==4 && false) {

    int glast=Ngrids-1;
    real radlast=0.5;

    
    //path semi=(1,0){N}..(0,1){W}..{(0,-1)}(-1,0);
    for(int G=0; G < Ngrids; ++G) {
      int L=2^G;
      real istart=G==0 ? 0 : quotient(m,2);
      real rad = istart*L;
      int gmy=m;
      real kmy=L*gmy;
      real stoprad=G ==glast ? sqrt(2)*kmy-L : kmy;

      real[] a={G*pi/Ngrids,(2G+1)*pi/(Ngrids*2),(2G+2)*pi/(Ngrids*2)};
      pair[] p={(cos(a[0]),sin(a[0])),
		(cos(a[1]),sin(a[1])),
		(cos(a[2]),sin(a[2]))};
      path semi=p[0]{(-p[0].y,p[0].x)}
      ..p[1]{(-p[1].y,p[1].x)}
      ..p[2]{(-p[2].y,p[2].x)};
      path Semi=(1,0){N}..(0,1){W}..{S}(-1,0);
      while(rad < stoprad) {

	// only draw separate shells
	//path[] paths={scale(rad+L/2)*semi,scale(radlast)*semi};

	// draw overlapping shells
	path[] paths={scale(rad+L/2)*semi,scale(max(L/2,rad-L/2))*semi};

	pen fpen=dotpen[G]+opacity(0.3);
	pen[] pens={fpen,fpen};
	draw(pic,paths,pens);
	
	//draw(pic,scale(rad+L/2)*Semi,dotpen[G]+dashed);
	radlast=rad+L/2;
	rad += L;
      }


      int imax=floor(sqrt(2)*m);
      for(int i=1; i <= imax; ++ i) {
	draw(pic,scale(L*(i-0.5))*semi,dotpen[G]+dashed);
      }
      dot(pic,(L*(imax+0.5),0),invisible); // centres picture correctly
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

