# AllHadronicSUSY

## Instructions

```
cmsrel CMSSW_7_4_6_patch6

cd CMSSW_7_4_6_patch6/src/

cmsenv

git cms-merge-topic -u cms-met:METCorUnc74X

git clone https://github.com/TreeMaker/TreeMaker.git -b Run2

git clone https://github.com/chrosa/AllHadronicSUSY.git

scram b -j 8
```
