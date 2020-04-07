HiggsAnalysis-CombinedLimit
===========================

### Official documentation

[Manual to run combine](http://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/)

### Standalone compilation in `lxplus`
```
official:
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit

Tao's: 
git clone https://github.com/tahuang1991/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
git checkout TaoLimits_102x
cd HiggsAnalysis/CombinedLimit
source env_standalone.sh 
make -j 8; make # second make fixes compilation error of first
```
