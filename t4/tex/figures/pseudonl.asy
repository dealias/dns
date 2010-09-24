size(0,2cm);
usepackage("amsmath,mathdef,pde");

draw((-1,0)..(-0.7,0),EndArrow);
draw((0.7,0)..(1,0),EndArrow);

label("$\frac{\partial }{\partial t}\hat{\v{u}}_{\v k}
+ \[(\hat{\v{u}}\cdot \v{k} ) \hat{\v{u}}\]_{\v k}
= \frac{1}{\rho}i \v{k} P - \nu k^2 \hat{\v{u}}_{\v k}$",(0,0));

//dot((-2,0),invisible);
//dot((1,0),invisible);
real x=-0.18;
real y=-0.05;
real dy=-0.2;
//dot((x,y),red);
draw((x,y)..(x,y+dy),Arrows);
label("$(\v{u} \cdot \del)\v{u}$",(x,2y+dy));




