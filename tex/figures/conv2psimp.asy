size(11cm,0);


pair center(path g) {return 0.5*(min(g)+max(g));}


void drawboxes(pair z=(0,0), pair p=(1,1))
{
  path g=box((0,0),p);
  draw(shift(z)*g);
}
path g, p;
pair b=(1,1), a=(0,0), z=(0,0);
real wx=0.01;
real dd=1.2*wx;
real dh=1.1*wx;

g=box(a,(wx,wx));
//draw(shift(z)*g);

pen datapen=yellow, convpen=green, labelpen=black;
int here=0;
void iframe(int i) {
  if (here == i) {
    datapen=invisible;
    convpen=invisible;
    labelpen=invisible;
    currentpen=invisible;
  }
  ++here;
}

pen[] b1={datapen,datapen,datapen,datapen,convpen,convpen,convpen};
pen[] b2={invisible,datapen,datapen,datapen,convpen,convpen,invisible};
pen[] b3={datapen,datapen,datapen,datapen,convpen,convpen,invisible};
pen[] b4={invisible,datapen,datapen,datapen,convpen,convpen,invisible};
pen[] u1={invisible,invisible,datapen,convpen,convpen,invisible,invisible};
pen[] u2={invisible,invisible,datapen,convpen,convpen,invisible,invisible};

string[] l1N={"$F$",
	      "${\rm FFT}^{-1}_x\{F\}$",
	      "","","",
	      "${\rm FFT}^{-1}_x\{F*G\}$",
	      "$F*G$"};
string[] l1S={"",
	      "$n_x$ \rm{even}",
	      "","","",
	      "$n_x$ \rm{even}",
	      ""};
string[] l2N={"",
	      "${\rm FFT}^{-1}_x\{F\}$",
	      "","","",
	      "${\rm FFT}^{-1}_x\{F*G\}$",
	      ""};
string[] l2S={"",
	      "$n_x$ \rm{odd}",
	      "","","",
	      "$n_x$ \rm{odd}",
	      ""};
string[] l3N={"$G$",
	      "${\rm FFT}^{-1}_x\{G\}$",
	      "","","",
	      "${\rm FFT}^{-1}_x\{F*G\}$",
	      ""};
string[] l3S={"",
	      "$n_x$ \rm{even}",
	      "","","",
	      "$n_x$ \rm{even}"
	      ,""};
string[] l4N={"",
	      "${\rm FFT}^{-1}_x\{G\}$",
	      "","","",
	      "${\rm FFT}^{-1}_x\{F*G\}$",
	      ""};
string[] l4S={"",
	      "$n_x$ \rm{odd}",
	      "","","",
	      "$n_x$ \rm{odd}",
	      ""};

for(int i=0; i < b1.length; ++i){
  picture pic;
  size(pic,16cm,0);
  
  here=0;
  datapen=yellow;
  convpen=lightgreen;
  labelpen=black;
  currentpen=black;
  real ddh=0.1dh+2dd;
  
  
  
  z=(0,0);
  p=shift(z)*g;
  filldraw(pic,p,b1[i]);
  label(pic,l1N[i],center(p),N,labelpen);
  label(pic,l1S[i],center(p),S,labelpen);
  
  z=(dh,0);
  p=shift(z)*g;
  filldraw(pic,p,b2[i]);
  label(pic,l2N[i],center(p),N,labelpen);
  label(pic,l2S[i],center(p),S,labelpen);
  
  z=(ddh,0);
  p=shift(z)*g;
  filldraw(pic,p,b3[i]);
  label(pic,l3N[i],center(p),N,labelpen);
  label(pic,l3S[i],center(p),S,labelpen);
  
  z=(ddh+dh,0);
  p=shift(z)*g;
  filldraw(pic,p,b4[i]);
  label(pic,l4N[i],center(p),N,labelpen);
  label(pic,l4S[i],center(p),S,labelpen);
  

  real x=0;
  real wu=wx/8;
  real X= i < 4 ? 0 : dh+wx-wu;
  p=shift(x+X,dh)*box(a,(wu,wx));
  filldraw(pic,p,u1[i]);
  p=shift(0,0)*box(a,(wu,wx));
  fill(pic,p,u1[i]);
  //restore the overwritten box edge
  draw(pic,shift((x,0))*g);

  x=ddh;
  p=shift(x+X,dh)*box(a,(wu,wx));
  filldraw(pic,p,u2[i]);
  p=shift(ddh,0)*box(a,(wu,wx));
  fill(pic,p,u1[i]);

  //restore the overwritten box edge
  draw(pic,shift((x,0))*g);
  
  //label(pic,"asdf",center(p),labelpen);
  
  shipout("conv2psimp"+(string) i,pic);
}
