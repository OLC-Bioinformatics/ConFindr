import glob
import os

miseq_runs = glob.glob('/mnt/nas/MiSeq_Backup/??????')
miseq_runs = sorted(miseq_runs)
for run in miseq_runs:
	cmd = 'python3 Detector.py ' + run + ' ' + run.split('/')[-1] + ' -tr'
	os.system(cmd)
	# print(cmd)

miseq_runs = glob.glob('/mnt/nas/External_MiSeq_Backup/???/??????')
miseq_runs = sorted(miseq_runs)
for run in miseq_runs:
	cmd = 'python3 Detector.py ' + run + ' ' + run.split('/')[-2] + run.split('/')[-1] + ' -tr'
	os.system(cmd)
	# print(cmd)
