import math;
import graph;

size(15cm,0);
int mx=4;
int my=mx;
int Nx=2*mx-1;
int Ny=2*my-1;
pen p=linewidth(0.75mm);

fill((-mx+0.5,0.5)--(-mx+0.5,my-0.5)--(mx-0.5,my-0.5)--(mx-0.5,-0.5)--(0.5,-0.5)--(0.5,0.5)--cycle,grey);

add(shift(-mx+0.5,-my+0.5)*grid(Nx,Ny));

dot((mx-1,0),red);
dot((-mx+1,0),red);
dot((mx-1,my-1),red);
dot((-mx+1,my-1),red);
label("(0,0)",(0,0),red);

label("$(-m_x+1,0)$",(-mx+0.5,0),W,red);
label("$(m_x-1,0)$",(mx-0.5,0),E,red);
label("$(m_x-1,m_y-1)$",(mx-1,my-0.5),N,red);
label("$(-m_x+1,m_y-1)$",(-mx+1,my-0.5),N,red);

draw("$N_x=2m_x-1$",brace((mx-0.5,-my+0.5),(-mx+0.5,-my+0.5)),S,red);
int o=-2;
draw("$N_y=2m_y-1$",brace((-mx+0.5+o,-my+0.5),(-mx+0.5+o,my-0.5)),W,red);

draw(box((-0.5,-0.5),(0.5,0.5)),p);
draw((-mx+0.5,-0.5)--(-mx+0.5,my-0.5)--(mx-0.5,my-0.5)--(mx-0.5,-0.5)--cycle,
     blue+p);

//axes(Arrow);
