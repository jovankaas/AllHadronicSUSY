# AllHadronicSUSY

## Instructions

```
cmsrel CMSSW_7_4_15
cd CMSSW_7_4_15/src/
cmsenv


git cms-merge-topic -u kpedro88:METfix7415

git clone https://github.com/TreeMaker/TreeMaker.git -b Run2

git clone https://github.com/jovankaas/AllHadronicSUSY.git

scram b -j 8
```
