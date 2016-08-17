size(0,300);

import flowchart;

block block1=rectangle("Initialize $w_r$",(0,3),palered);
block block2=rectangle("Calculate $S_r$",(0,2),palered);
block block3=roundrectangle("Do something",(-1,1));
block block4=bevel("Don't do something",(1,1));
block block5=circle("End",(0,0));

draw(block1);
draw(block2);
draw(block3);
draw(block4);
draw(block5);

add(new void(picture pic, transform t) {
    blockconnector operator --=blockconnector(pic,t);
    //    draw(pic,block1.right(t)--block2.top(t));
    block1--Down--Arrow--block2;
    block2--Label("Yes",0.5,NW)--Left--Down--Arrow--block3;
    block2--Right--Label("No",0.5,NE)--Down--Arrow--block4;
    block4--Down--Left--Arrow--block5;
    block3--Down--Right--Arrow--block5;
  });
