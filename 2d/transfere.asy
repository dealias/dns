include getparam;
include averages;

scale(Log,Linear);

pen p=linewidth(1);

string Pitext="$\Pi_E$";
string etatext="$\epsilon$";

real[] Pi,Eps;
real[][] M2,Tk;

while(nextrun()) {
  //  gettime(n == 0);
  gettime();
  M2=moment2();
  Tk=transfer();

  write("  Enery injection rate=",2*sum(Tk[EPS]));
  write("Energy dissipation rate=",2*sum(Tk[DE]));

  Pi=-2*partialsum(Tk[TE]);
  Pi.insert(0,0);
  Eps=2*reverse(partialsum(reverse(Tk[DE]-Tk[EPS])));
  Eps.push(0);
  kb[0]=k0;
  
  string runtext=" ("+run+")";
  draw(graph(kb,Pi),p+Pen(2*n+1),Pitext+runtext);
  if(!all(Eps == 0)) 
    draw(graph(kb,Eps),p+Pen(2*n)+dashed,etatext+runtext);
};

if(n == 1) {
  currentpicture.legend[0].label=Pitext;
  if(currentpicture.legend.length > 1)
    currentpicture.legend[1].label=etatext;
}

xaxis("$k$",BottomTop,LeftTicks);
yaxis("Cumulative energy transfer",LeftRight,RightTicks(trailingzero));

yequals(0,Dotted);

attach(legend(),point(E),20E);
