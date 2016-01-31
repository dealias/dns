import graph;
include getparam;
include averages;
size(200,150,IgnoreAspect);

while(nextrun()) {
  real[][] a;
  real[] time;
  file fin=input(run+"/stat",comment="").line();
  
  string[] names=fin.word();
    
  int field=find(names == "t")-1;
  a=fin.word(false);
  a=transpose(a);
  time=a[field];
  
  real[] stops;
  real[] it=a[0];
  real[] cpu=a[5];

  for (int i=1; i < it.length; ++i) {
    if (it[i] == 1) {
      stops.push(cpu[i-1]);
    }
  }

  int field2=find(names == "CPU")-1;
  
  // if there are quite a lot of data points, thin the data for the graph
  int maxplotlength=1000;
  int skip = ceil(time.length/maxplotlength);
  real[] aS, timeS;
  for(int i=0; i < time.length; ++i) {
    if(i % skip == 0) {
      aS.push(cpu[i]);
      timeS.push(time[i]);
    }
  }
  
  //  for (int i=0; i < stops.length; ++i)
  //    xequals(stops[i],Pen(n)+dashed);
  //draw(graph(a[field2],time),p+Pen(n),texify(run));
  draw(graph(aS,timeS),p+Pen(n),texify(run));

  if (cpu.length > 1) {
    write("slope=",(time[time.length-1]-time[0])/(cpu[cpu.length-1]-cpu[0]));
  } else
    write("insufficient progress.");

}

if(n > 1) attach(legend(),point(E),20E);

xaxis("CPU",BottomTop,LeftTicks);
yaxis("$t$",LeftRight,RightTicks(trailingzero));
