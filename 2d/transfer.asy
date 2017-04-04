include getparam;
include averages;

scale(Log,Linear);

pen p=linewidth(1);

string Pitext="$\Pi$";
string etatext="$\eta$";

real[] Pi,Eta;
real[][] M2,Tk;

while(nextrun()) {
  //  gettime(n == 0);
  gettime();
  M2=moment2();
  Tk=transfer();

  write("  Enstrophy injection rate=",2*sum(Tk[ETA]));
  write("Enstrophy dissipation rate=",2*sum(Tk[DZ]));

  Pi=-2*partialsum(Tk[TZ]);
  Pi.insert(0,0);
  Eta=2*reverse(partialsum(reverse(Tk[DZ]-Tk[ETA])));
  Eta.push(0);
  kb[0]=k0;
  
  string runtext=" ("+run+")";
  draw(graph(kb,Pi),p+Pen(2*n+1),Pitext+runtext);
  if(!all(Eta == 0)) 
    draw(graph(kb,Eta),p+Pen(2*n)+dashed,etatext+runtext);
};

if(n == 1) {
  currentpicture.legend[0].label=Pitext;
  if(currentpicture.legend.length > 1)
    currentpicture.legend[1].label=etatext;
}

xaxis("$k$",BottomTop,LeftTicks);
yaxis("Cumulative enstropy transfer",LeftRight,RightTicks(trailingzero));

yequals(0,Dotted);

attach(legend(),point(E),20E);
