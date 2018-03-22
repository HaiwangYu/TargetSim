# TargetSim

## libraries
[Modified sPHENIX g4main and g4detector](https://github.com/HaiwangYu/coresoftware/tree/TargetSim)

```
git clone git@github.com:HaiwangYu/coresoftware.git
git checkout -b TargetSim remotes/origin/TargetSim
```
Then compile the g4main and g4detector package

If on RCF, could also use my libraries by
```
export LD_LIBRARY_PATH=/phenix/hhj/yuhw/GitHub/HaiwangYu/install_TargetSim/lib/:$LD_LIBRARY_PATH
```

## analysis module
'TruthEval' contains the analysis module.
Need to compile it before running 'Fun4All_G4_Target.C'

If on RCF, could also use my libraries by
```
export LD_LIBRARY_PATH=/phenix/hhj/yuhw/GitHub/HaiwangYu/install_TargetSim/lib/:$LD_LIBRARY_PATH
```

## How to run
```
git clone git@github.com:HaiwangYu/TargetSim.git
cd TargetSim
root -l Fun4All_G4_Target.C
```

- 'Fun4All_G4_Target.C' is the main macro to run the simulation
- 'ana.C' is the ploting macro
