cd /nfs/fanae/user/juanr/CMSSW_7_4_7_patch1/StopLimits2016
source /cms/cmsset_default.sh
eval `scramv1 runtime -sh`
python LimitsSR.py  900 350 10000.0
