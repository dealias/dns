size(300,0);

import flowchart;

block block1=rectangle("Initialize $\omega_k$",(0,2.5),minwidth=100,palered);
block block2=rectangle("Calculate  $S_k$",(0,2.0),minwidth=100,palered);
block block3=rectangle("Update $\omega_k$",(0,1.5),minwidth=100,palered);
block block4=rectangle("Calculate $E_k, Z_k$",(0,1.0),minwidth=100,palered);
block block5=rectangle("Calculate $u,v $ based on $\omega_k$",(0.8,2.25),minwidth=150,palered);
block block6=rectangle("Calculate convolutions",(0.8,2),minwidth=150,palered);
block block7=rectangle("No mean flow",(0.8,1.75),minwidth=150,palered);

//block block3=roundrectangle("Do something",(-1,1));
//block block4=bevel("Don't do something",(1,1));
//block block5=circle("End",(0,0));

draw(block1);
draw(block2);
draw(block3);
draw(block4);
draw(block5);
draw(block6);
draw(block7);

draw(brace((0.375,1.70),(0.375,2.30)),E,black);

add(new void(picture pic, transform t) {
    blockconnector operator --=blockconnector(pic,t);
    //draw(pic,block1.right(t)--block2.top(t));
    block1--Down--Arrow--block2;
    block2--Down--Arrow--block3;
    block3--Down--Arrow--block4;
    block4--Left--Up--Right--Arrow--block2;
    block5--Down--Arrow--block6--Down--Arrow--block7;


    //block2--Right--Label("No",0.5,NE)--Down--Arrow--block4;
    //block4--Down--Left--Arrow--block5;
    //block3--Down--Right--Arrow--block5;
    
  });
