#!/opt/sphenix/core/bin/python3

from subprocess import call
import sys
import fileinput
import os.path

jobname = sys.argv[1]
#input_lst_name = sys.argv[2]
njobs = int(sys.argv[2])
nevents = int(sys.argv[3])
do_submition = bool(sys.argv[4])

sphenix_build_tag = 'new'
local_build_tag = 'install'

condor_base_path = '/gpfs/mnt/gpfs04/sphenix/user/yuhw/workspace/TargetSim/TargetSim'
macros_path = condor_base_path
output_path = condor_base_path  + '/output/' + '/{}/'.format(jobname)
condor_work_dir_base = condor_base_path + '/work/' + '/{}/'.format(jobname)

condor_local = condor_base_path + '/temp/'

if do_submition :
	condor_local = '${_CONDOR_SCRATCH_DIR}'

call(["mkdir","-p",output_path])


#ijob = 0
#input_file_names = [line.rstrip('\n') for line in open(input_lst_name)]
#for input_file_name in input_file_names:
for ijob in range(0, njobs) : 
# index is the input file name (base name) without extension
	#index = os.path.splitext(os.path.basename(input_name_abs))[0];
	#output_name_abs = output_path+'/'+index+".root"
	#print output_name_abs
	condor_work_dir = condor_work_dir_base + '/{}/'.format(ijob)
	#print condor_work_dir

	call(["mkdir","-p",condor_work_dir])
	#call(["ln","-sf",C_macro_abs_path,condor_work_dir])

# make executable: run.csh
	runscript = condor_work_dir + '/run.csh'
	RUNSCRIPT = open(runscript,'w')
	RUNSCRIPT.write('#!/bin/tcsh -f\n')
	RUNSCRIPT.write('#automagically generated by submit.py\n')
	RUNSCRIPT.write('unsetenv OFFLINE_MAIN \n')
	RUNSCRIPT.write('unsetenv ONLINE_MAIN\n')
	RUNSCRIPT.write('unsetenv ROOTSYS\n')
	RUNSCRIPT.write('unsetenv LD_LIBRARY_PATH\n')
	RUNSCRIPT.write('source /opt/sphenix/core/bin/sphenix_setup.csh -n {} \n'.format(sphenix_build_tag))
	RUNSCRIPT.write('setenv LD_LIBRARY_PATH $GITHUB_ROOT/HaiwangYu/{}/lib/:$LD_LIBRARY_PATH \n'.format(local_build_tag))
	RUNSCRIPT.write('source /opt/phenix/bin/odbcini_setup.csh \n')
	RUNSCRIPT.write('echo $OFFLINE_MAIN \n')

	RUNSCRIPT.write('cp {}/*.C {} \n'.format(macros_path,condor_local))
	RUNSCRIPT.write('cp {}/*.cfg {} \n'.format(macros_path,condor_local))

	RUNSCRIPT.write('cd {} \n'.format(condor_local))

	cmd = 'time root -l -b -q Fun4All_G4_E1039_R1.C\\({}\\)'.format(nevents)
	RUNSCRIPT.write(cmd + ' | tee {}.out \n'.format(ijob))

	cmd = 'time root -l -b -q ana.C\\(\\)'
	RUNSCRIPT.write(cmd + ' | tee -a {}.out \n'.format(ijob))

	RUNSCRIPT.write('cp DSTReader.root {}/dstr_{:06d}.root \n'.format(output_path,ijob))
	RUNSCRIPT.write('cp eval.root {}/eval_{:06d}.root \n'.format(output_path,ijob))
	RUNSCRIPT.write('cp hist.root {}/hist_{:06d}.root \n'.format(output_path,ijob))
	RUNSCRIPT.write('cp *.out {} \n'.format(output_path))

	RUNSCRIPT.close()

	call(['chmod','0755',runscript])
# make condor job file: condor.job
	condorfile = condor_work_dir + '/condor.job'
	logfile = condor_work_dir + '/log.txt'
	outfile = condor_work_dir + '/out.txt'
	errfile = condor_work_dir + '/err.txt'

	CONDORFILE = open(condorfile,'w')
	CONDORFILE.write('#automagically generated by submit.py\n')
	CONDORFILE.write('Executable = ' + runscript + '\n')
	CONDORFILE.write('Universe   = vanilla\n')
	CONDORFILE.write('GetEnv     = True\n')
	CONDORFILE.write('Requirements = (CPU_Experiment == "phenix")\n')
	CONDORFILE.write('+Experiment  = \"phenix\"\n')
	CONDORFILE.write('+Job_Type  = \"cas\"\n')
	#CONDORFILE.write('+Job_Type  = \"highmem\"\n')
	CONDORFILE.write('Initialdir = ' + condor_work_dir + '\n')
	CONDORFILE.write('Output  = ' + outfile + '\n')
	CONDORFILE.write('Error = ' + errfile + '\n')
	CONDORFILE.write('Log  = ' + logfile + '\n')
	CONDORFILE.write('Notification = Never\n')
	CONDORFILE.write('Queue \n')
	CONDORFILE.close()

	if do_submition :
		call(['/usr/bin/condor_submit',condorfile])
		#print("");
	else :
		call([runscript])

	#ijob = ijob+1;
	
print('SUCCESS!')
