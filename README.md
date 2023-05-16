# bph-hlt-tools
Ntuples producer + Trigger plots for CMS BPAG

Instructions for running with 2022 data
```
cmsrel CMSSW_12_4_3
cd CMSSW_12_4_3/src && cmsenv
git clone https://github.com/vjmastra/bph-hlt-tools.git myAnalyzers/bph-hlt-tools
scram b -j8
```

Instructions for running with 2023 data
```
cmsrel CMSSW_13_0_3
cd CMSSW_13_0_3/src && cmsenv
git clone -b 2023data https://github.com/horace-cl/bph-hlt-tools.git myAnalyzers/bph-hlt-tools
scram b -j8
```

Run a simple test to produce ntuples
```
cmsRun myAnalyzers/bph-hlt-tools/test/PsikaonRootupler.py
```
