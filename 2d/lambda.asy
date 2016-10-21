real kforce=getreal("kforce");
real deltaf=getreal("deltaf");
real kfm=kforce-deltaf/2;
real kfp=kforce+deltaf/2;

int K=ceil(kfp);

real[] Lambda;
int[] Count;

for(int i=-K; i <= K; ++i) {
  for(int j=-K; j <= K; ++j)  {
    real lambda=hypot(i,j);
    int i=search(Lambda,lambda);
    if(i >= 0 && Lambda[i] == lambda)
        Count[i] += 1;
    else {
      Lambda.insert(i+1,lambda);
      Count.insert(i+1,1);
    }
  }
}

write(Lambda,Count);
write();
write("average modes/ring=",sum(Count)/Count.length);
write();

int im=search(Lambda,kfm);
int ip=search(Lambda,kfp);

if(im < 0 || Lambda[im] < kfm) ++im;
if(ip < 0 || Lambda[ip] < kfp) ++ip;

real sum=0;

for(int i=im; i < ip; ++i) {
  sum += Count[i];
  write(Lambda[i],Count[i]);
}

write(sum,none);
write(" modes.");
