import graph;
size(500,100);

real f(real t) {
  if(t > -1 && t <1)
    return 1.3;
  return 0;

};


draw(graph(f,-6,6),invisible);

draw((-6,0)..(-1,0),blue);
draw((-1,0)..(-1,1),blue);
draw((-1,1)..(1,1),blue);
draw((1,0)..(1,1),blue);
draw((1,0)..(6,0),blue);

//axes();

//xaxis("$t$",BottomTop,LeftTicks);
//yaxis("$y(t)$",LeftRight,RightTicks);