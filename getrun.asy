public string runs;
public string run;

public int lastpos;

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

getrun();

size(230,200,IgnoreAspect);
pen p=linewidth(1);
marker mark=marker(scale(0.8mm)*cross(4),p,above=false);
defaultpen(fontsize(10));
