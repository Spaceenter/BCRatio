import math, sys;

eka_list = [0.14791,0.192765,0.253927,0.335611,0.439944,0.5637,0.703865,0.875894,1.07782,1.30478,1.57625,1.88956,2.23265,2.63464,3.101,3.62237,4.21492,4.91096,5.70801,6.60399,7.6715,8.91966,10.331,12.0382,14.0708,16.4138,19.2986,22.7862,26.9574,32.1969,38.6047,46.6328,56.4995,68.6876,84.322,104.197,130.658,165.985,214.797,287.611,408.689,669.728,1634.63];
mc = 11.1779;

num = 0;
for eka in eka_list:
  rig = math.sqrt((eka*12+mc)*(eka*12+mc)-mc*mc)/6;
  sys.stdout.write("%4.3f, " % rig);
  num = num+1;
print ;
print "number of bin edges = "+str(num);
print "number of bins= "+str(num-1);
print ;

num = 0;
truenum = 0;
for eka in eka_list:
  rig = math.sqrt((eka*12+mc)*(eka*12+mc)-mc*mc)/6;
  if((num%2)==1 or eka==1634.63): 
    sys.stdout.write("%4.3f, " % rig);
    truenum = truenum+1;
  num = num+1;
print ;
print "number of bin edges = "+str(truenum);
print "number of bins= "+str(truenum-1);
print ;
