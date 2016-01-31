import graph;
include getparam;
include averages;
size(200,150,IgnoreAspect);

string fieldname=getstring("field","dt");

bool dt=fieldname == "dt";

if(dt)
  scale(Linear,Log);


while(nextrun()) {
  real[][] a;
  real[] data;
  file fin=input(run+"/stat",comment="").line();
  
  string[] names=fin.word();
  
  int field=find(names == fieldname)-1;
  if(field < 0) abort("No such field: "+fieldname);
  a=fin.word(false);
  a=transpose(a);
  data=a[field];

  real recount=0;
  real[] stops;
  for (int i=1; i < a[0].length; ++i) {
    if (a[0][i] == 1) {
      recount += a[0][i-1];
      stops.push(a[0][i-1]);
    }
    a[0][i] += recount;
  }

  // if there are quite a lot of data points, thin the data for the graph
  int maxplotlength=1000;
  int skip = ceil(data.length/maxplotlength);
  real[] aS, dataS;
  for(int i=0; i < data.length; ++i) {
    if(i % skip == 0) {
      aS.push(a[0][i]);
      dataS.push(data[i]);
    }
  }

  bool[] allgood=array(dataS.length,true);
  draw(graph(aS,dataS,dt ? dataS > 0 : allgood), p+Pen(n),texify(run));

  // draw vertical lines when the run was restarted
  for (int i=1; i < stops.length; ++i)
    xequals(stops[i],Pen(n)+dashed);
  
  real[] data0=copy(data);
  data0.delete(0);
  write(run+":");
  if(data0.length > 0) {
    write(" max=",max(data));
    write(" min=",min(data));
    write(" min0=",min(data0));
    write(" mean=",sum(data)/data.length);
    write();
  } else {
    write("no progress");
  }
  
}

if(n > 1) attach(legend(),point(E),20E);

xaxis("it",BottomTop,LeftTicks);
yaxis(texify(fieldname),LeftRight,RightTicks(trailingzero));
