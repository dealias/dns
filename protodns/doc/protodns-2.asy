if(!settings.multipleView) settings.batchView=false;
settings.tex="pdflatex";
defaultfilename="protodns-2";
if(settings.render < 0) settings.render=4;
settings.outformat="";
settings.inlineimage=true;
settings.embed=true;
settings.toolbar=false;
viewportmargin=(2,2);

size(13cm);

import flowchart;

block init=rectangle("Initialize $\omega_k$",(0,3),minwidth=100,palered);
block Sk=rectangle("Calculate  $S_k$",(0,2),minwidth=100,palered);
block wk=rectangle("Update $\omega_k$",(0,1),minwidth=100,palered);
block Ek=rectangle("Calculate $E(k)$",(0,0),minwidth=100,palered);
real x=1.222;
block uv=rectangle("Calculate $u,v $ based on $\omega_k$",(x,2.5),minwidth=150,palered);
block conv=rectangle("Calculate convolutions",(x,1.5),minwidth=150,palered);

draw(init);
draw(Sk);
draw(wk);
draw(Ek);
draw(uv);
draw(conv);

add(new void(picture pic, transform t) {
    blockconnector operator --=blockconnector(pic,t);
    draw(pic,brace(conv.left(t),uv.left(t)));
    init--Down--Arrow--Sk;
    Sk--Down--Arrow--wk;
    wk--Down--Arrow--Ek;
    Ek--Left--Up--Right--Arrow--Sk;
    uv--Down--Arrow--conv--Down;
  });

