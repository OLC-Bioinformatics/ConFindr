import glob
import os

miseq_runs = glob.glob('/mnt/nas/MiSeq_Backup/??????')
miseq_runs = sorted(miseq_runs)
for run in miseq_runs:
	cmd = 'python3 New_Detector.py ' + run + ' ' + run.split('/')[-1] + ' database.fasta'
	print(cmd)
	if not os.path.exists(run.split('/')[-1] + '.csv'):
		os.system(cmd)


miseq_runs = glob.glob('/mnt/nas/External_MiSeq_Backup/???/??????')
miseq_runs = sorted(miseq_runs)
for run in miseq_runs:
	cmd = 'python3 New_Detector.py ' + run + ' ' + run.split('/')[-2] + run.split('/')[-1] + ' database.fasta'
	print(cmd)
	os.system(cmd)

