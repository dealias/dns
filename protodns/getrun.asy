public string runs;
public string run;
public int lastpos;

int nruns=-1;

void getrun(string default="") 
{
  runs=getstring("run",default);
  run=runs;
  lastpos=0;
} 

string rundir(string dir="")
{
  string s=run+"/";
  if(dir == "") return s;
  return s+dir+"/";
}

bool nextrun()
{
  ++nruns;
  int pos=find(runs,",",lastpos);
  if(lastpos == -1) {run=""; return false;}
  run=substr(runs,lastpos,pos-lastpos);
  lastpos=pos > 0 ? pos+1 : -1;
  return true;
}

getrun();

pen p=linewidth(1);
marker mark=marker(scale(0.8mm)*cross(4),p,above=false);
defaultpen(fontsize(10));
