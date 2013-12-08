import math;

MC12 = 11.1779;

Beta = 0.5;
vek = (1/math.sqrt(1-Beta*Beta)-1)*MC12/12;
vrig = MC12/6/math.sqrt(1/Beta/Beta-1);
print("TOF Beta = "+str(Beta)+" <--> Ek/A = "+str(vek)+" <--> Rig = "+str(vrig));
Beta = 0.9;
vek = (1/math.sqrt(1-Beta*Beta)-1)*MC12/12;
vrig = MC12/6/math.sqrt(1/Beta/Beta-1);
print("TOF Beta = "+str(Beta)+" <--> Ek/A = "+str(vek)+" <--> Rig = "+str(vrig));
Beta = 0.96;
vek = (1/math.sqrt(1-Beta*Beta)-1)*MC12/12;
vrig = MC12/6/math.sqrt(1/Beta/Beta-1);
print("RICH Beta = "+str(Beta)+" <--> Ek/A = "+str(vek)+" <--> Rig = "+str(vrig));
Beta = 0.99;
vek = (1/math.sqrt(1-Beta*Beta)-1)*MC12/12;
vrig = MC12/6/math.sqrt(1/Beta/Beta-1);
print("RICH Beta = "+str(Beta)+" <--> Ek/A = "+str(vek)+" <--> Rig = "+str(vrig));
