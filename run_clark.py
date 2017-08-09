import os
import csv
import subprocess


def classify_metagenome(dir_db, fastq_files, cpus):
    cwd = os.getcwd()
    cmd = 'which set_targets.sh'
    clark_dir = subprocess.check_output(cmd.split())
    clark_dir = clark_dir.split('/')[:-1]
    clark_dir = '/'.join(clark_dir)
    os.chdir(clark_dir)
    with open(cwd + '/junk.txt', 'w') as outjunk:
        cmd = './set_targets.sh ' + dir_db + ' bacteria'
        subprocess.call(cmd, shell=True, stderr=outjunk, stdout=outjunk)
        cmd = './classify_metagenome.sh -P ' + fastq_files[0] + ' ' + fastq_files[1] + ' -R results -n ' + str(cpus) + ' --light'
        if '.gz' in fastq_files[0]:
            cmd += ' --gzipped'
        subprocess.call(cmd, shell=True, stderr=outjunk, stdout=outjunk)
        cmd = './estimate_abundance.sh -F results.csv -D ' + dir_db + ' > abundance.csv'
        subprocess.call(cmd, shell=True, stderr=outjunk)
    os.rename('abundance.csv', cwd + '/abundance.csv')
    # os.remove(cwd + '/junk.txt')
    os.chdir(cwd)


def read_clark_output(result_file):
    results = ''
    with open(result_file) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            try:
                if float(row['Proportion_Classified(%)']) > 1:
                    results += row['Name'] + '_' + row['Proportion_Classified(%)'] + ':'
            except ValueError:
                pass
    return results

