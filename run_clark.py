import os
import csv
import subprocess


# TODO: Get things converted to subprocess calls so stderr can get redirected into the void of nothingness.
def classify_metagenome(dir_db, fastq_files, cpus):
    cwd = os.getcwd()
    cmd = 'which set_targets.sh'
    clark_dir = subprocess.check_output(cmd.split())
    clark_dir = clark_dir.split('/')[:-1]
    clark_dir = '/'.join(clark_dir)
    os.chdir(clark_dir)
    cmd = './set_targets.sh ' + dir_db + ' bacteria'
    os.system(cmd)
    cmd = './classify_metagenome.sh -P ' + fastq_files[0] + ' ' + fastq_files[1] + ' -R results -n ' + str(cpus) + ' --light'
    os.system(cmd)
    cmd = './estimate_abundance.sh -F results.csv -D ' + dir_db + ' > abundance.csv'
    os.system(cmd)
    os.rename('abundance.csv', cwd + '/abundance.csv')
    os.chdir(cwd)


def read_clark_output(result_file):
    results = ''
    with open(result_file) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            try:
                if float(row['Proportion_Classified(%)']) > 1:
                    print(row)
                    results += row['Name'] + '_' + row['Proportion_Classified(%)'] + ':'
            except ValueError:
                pass
    return results

