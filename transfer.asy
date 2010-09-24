include getparam;
include averages;

scale(Log,Linear);

pen p=linewidth(1);

string Pitext="$\Pi$";
string epstext="$\eta$";

real[] Pi,Epsilon;
real[][] M2,Tk;

while(nextrun()) {
  gettime(n == 0);
  M2=moment2();
  Tk=transfer();
  Pi=-partialsum(Tk[NL]);
  Pi.insert(0,0);
  Epsilon=reverse(partialsum(reverse(M2[LIN]-Tk[LIN])));
  Epsilon.push(0);
  kb[0]=k0;
  
  string runtext=" ("+run+")";
  draw(graph(kb,Pi),p+Pen(2*n+1),Pitext+runtext);
  if(!all(Epsilon == 0)) 
    draw(graph(kb,Epsilon),p+Pen(2*n)+dashed,epstext+runtext);
};

if(n == 1) {
  currentpicture.legend[0].label=Pitext;
  if(currentpicture.legend.length > 1)
    currentpicture.legend[1].label=epstext;
}

xaxis("$k$",BottomTop,LeftTicks);
yaxis("Cumulative enstropy transfer",LeftRight,RightTicks(trailingzero));

yequals(0,Dotted);

attach(legend(),point(E),20E);
