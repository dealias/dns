size(12cm);

import flowchart;

int width = 100;
int height = 45;
block init=rectangle("Initialize $\omega_k$.",(-1,2),minheight = height, minwidth=width,palered);
block uk=rectangle(pack("Calculate  $u_k$ and","$v_k$ from $\omega_k$."),(0,2),minheight = height,minwidth=width,palered);
block fft=rectangle(pack("Calculate the","advection term."),(0,1),minheight = height,minwidth=width,palered);
block wk=rectangle("Update $\omega_k$.",(0,0),minheight = height,minwidth=width,palered);
block Ek=rectangle("Calculate $E$, $Z$, and $P$.",(-1,0),minheight = height,minwidth=width,palered);

draw(init);
draw(uk);
draw(fft);
draw(wk);
draw(Ek);

add(new void(picture pic, transform t) {
    blockconnector operator --=blockconnector(pic,t);
    init--Right--Arrow--uk;
    uk--Down--Arrow--fft;
    fft--Down--Arrow--wk;
    wk--Right--Up--Left--Arrow--uk;
    wk--Left--Arrow--Ek;
  });
