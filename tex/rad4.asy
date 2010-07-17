
string[] outnames={"rad1","rad2cross","rad4row","rad4cross"};
pair[][] al={{(0,0)},
	     {(0,1),(1,0),(0,-1),(-1,0)},
	     {(-1,0),(1,0)},
	     {(1,1),(-1,1),(-1,-1),(1,-1)}};
string[][] L={{"asdf"},
	      {"$\alpha$","$\beta$","$1-\alpha$","$1-\beta$"},
	      {"$\alpha$","$1-\alpha$"},
	      {"$\alpha$","$\beta$","$1-\alpha$","$1-\beta$"}};
path g;
real s;
pen p, fillpen;
picture pic;
pair a;
real loffset=4;

for(int i=0; i < outnames.length; ++i) {
  pic = new picture;
  size(pic,6cm,0);
  
  // resolved grid
  s=0.1;
  if(i==0)
    s=0.1;
  g=scale(s)*unitcircle;
  
  p=black;
  fillpen=black;

  a=(0,0);
  filldraw(pic,shift(a)*g,fillpen,p);
  //  label(pic,"$K$",a,loffset*S);
  
  // decimated grid
  p=blue;
  fillpen=invisible;
  s=0.2;
  g=scale(s)*unitcircle;


  for(int j=0; j < al[i].length; ++j) {
    a=al[i][j];
    //draw(pic,L[i][j],(0,0)..(1-s)*a,S,EndArrow);
    draw(pic,(0,0)..(1-s)*a,S,EndArrow);
    filldraw(pic,shift(a)*g,fillpen,p);
    //    label(pic,"$k_"+(string)(j+1)+"$",a,loffset*S);
    // FIXME: label aligns should not all be S.
  }
  shipout(outnames[i],pic);
}
