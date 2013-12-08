import os, random, time;

BaseCut = [15, 15, 0.9, 0.5, 0.5, 0.5, 15, 0.15, 0.5, 0.8];
Lower = [5, 5, 0.6, 0.3, 0.3, 0.3, 5, 0.1, 0.3, 0.5];
Upper = [50, 50, 0.96, 0.8, 0.7, 0.7, 50, 0.3, 0.7, 0.9];

for file_index in range(500):

  input = "input"+str(file_index)+".txt";
  output = "output"+str(file_index)+".root";

  file = open(input,"w");
  for iCut in range(10):
    file.write(str(BaseCut[iCut])+"\n");
    #Update BaseCut
    step = (Upper[iCut]-Lower[iCut])*(random.random());
    BaseCut[iCut] = BaseCut[iCut]+step;
    if BaseCut[iCut]>Upper[iCut]: 
      BaseCut[iCut] = BaseCut[iCut]-Upper[iCut]+Lower[iCut];
  file.close();

  job = str(file_index)+".job";
  file = open(job,"w");
  file.write("#BSUB-q 8nh\n");
  file.write("#BSUB-o job.output."+str(file_index)+".txt\n");
  file.write("#BSUB-J "+str(file_index)+".job\n");
  file.write("#BSUB-c 480\n");
  file.write("cd "+os.getcwd()+"\n");
  file.write("export ROOTSYS=/afs/cern.ch/exp/ams/Offline/AMSsoft/linux_slc5_gcc64/root_v5.27ams\n");
  file.write("source $ROOTSYS/bin/thisroot.sh\n");
  file.write("./BCRatio.exe "+input+" "+output+"\n");
  file.close();

  os.system("bsub < "+job);
  time.sleep(1);
  print(job+" is submitted.");
