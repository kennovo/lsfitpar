#!/usr/bin/awk -f
BEGIN{nam="illcond";n=360;pi=3.141592653589793238;deg=pi/180;pin=pi*2/n;degn=360/n;
  ntarg=1;mult[1]=pin;k[1]=1;d[1]=2*deg;
  nmeas=2;measd[1]=0.5;measd[2]=180;
  for(j=1;j<=nmeas;j++) printf("a%i b%i c%i d%i\n",j,j,j,j) > nam ".meas" j;
  for(i=0;i<n;i++){c=0;for(j=1;j<=ntarg;j++)c+=k[j]*cos(i*mult[j]+d[j]);
    printf("%.15g\n",c) > nam ".qme";print 0 > nam ".mme";
    for(j=1;j<=nmeas;j++) printf("%.15g\n",(i*degn+measd[j])%360) > nam ".meas" j};
  exit}
