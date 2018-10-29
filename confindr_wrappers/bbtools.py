import os
import subprocess
from subprocess import Popen, PIPE


def run_subprocess(command):
    """
    command is the command to run, as a string.
    runs a subprocess, returns stdout and stderr from the subprocess as strings.
    """
    x = Popen(command, shell=True, stdout=PIPE, stderr=PIPE)
    out, err = x.communicate()
    out = out.decode('utf-8')
    err = err.decode('utf-8')
    if x.returncode != 0:
        raise subprocess.CalledProcessError(x.returncode, cmd=command)
    return out, err


def kwargs_to_string(kwargs):
    """
    Given a set of kwargs, turns them into a string which can then be passed to a command.
    :param kwargs: kwargs from a function call.
    :return: outstr: A string, which is '' if no kwargs were given, and the kwargs in string format otherwise.
    """
    outstr = ''
    for arg in kwargs:
        outstr += ' {}={}'.format(arg, kwargs[arg])
    return outstr


def bbmap(reference, forward_in, out_bam, reverse_in='NA', returncmd=False, **kwargs):
    """
    Wrapper for bbmap. Assumes that bbmap executable is in your $PATH.
    :param reference: Reference fasta. Won't be written to disk by default. If you want it to be, add nodisk='t' as an arg.
    :param forward_in: Input reads. Should be in fastq format.
    :param out_bam: Output file. Should end in .sam or .bam
    :param returncmd: If set to true, function will return the cmd string passed to subprocess as a third value.
    :param reverse_in: If your reverse reads are present and normal conventions (_R1 for forward, _R2 for reverse) are
     followed, the reverse reads will be followed automatically. If you want to specify reverse reads, you may do so.
    :param kwargs: Other arguments to give to bbmap in parameter=argument format. See bbmap documentation for full list.
    :return: out and err: stdout string and stderr string from running bbmap.
    """
    options = kwargs_to_string(kwargs)
    if os.path.isfile(forward_in.replace('_R1', '_R2')) and reverse_in == 'NA' and '_R1' in forward_in:
        reverse_in = forward_in.replace('_R1', '_R2')
        cmd = 'bbmap.sh ref={} in={} in2={} out={} nodisk{}'.format(reference, forward_in, reverse_in, out_bam, options)
    elif reverse_in == 'NA':
        cmd = 'bbmap.sh ref={} in={} out={} nodisk{}'.format(reference, forward_in, out_bam, options)
    else:
        cmd = 'bbmap.sh ref={} in={} in2={} out={} nodisk{}'.format(reference, forward_in, reverse_in, out_bam, options)
    out, err = run_subprocess(cmd)
    if returncmd:
        return out, err, cmd
    else:
        return out, err


def bbduk_trim(forward_in, forward_out, reverse_in='NA', reverse_out='NA', returncmd=False, **kwargs):
    """
    Wrapper for using bbduk to quality trim reads. Contains arguments used in OLC Assembly Pipeline, but these can
    be overwritten by using keyword parameters.
    :param forward_in: Forward reads you want to quality trim.
    :param returncmd: If set to true, function will return the cmd string passed to subprocess as a third value.
    :param forward_out: Output forward reads.
    :param reverse_in: Reverse input reads. Don't need to be specified if _R1/_R2 naming convention is used.
    :param reverse_out: Reverse output reads. Don't need to be specified if _R1/_R2 convention is used.
    :param kwargs: Other arguments to give to bbduk in parameter=argument format. See bbduk documentation for full list.
    :return: out and err: stdout string and stderr string from running bbduk.
    """
    options = kwargs_to_string(kwargs)
    cmd = 'which bbduk.sh'
    try:
        subprocess.check_output(cmd.split()).decode('utf-8')
    except subprocess.CalledProcessError:
        print('ERROR: Could not find bbduk. Plase check that the bbtools package is installed and on your $PATH.\n\n')
        raise FileNotFoundError
    if os.path.isfile(forward_in.replace('_R1', '_R2')) and reverse_in == 'NA' and '_R1' in forward_in:
        reverse_in = forward_in.replace('_R1', '_R2')
        if reverse_out == 'NA':
            if '_R1' in forward_out:
                reverse_out = forward_out.replace('_R1', '_R2')
            else:
                raise ValueError('If you do not specify reverse_out, forward_out must contain R1.\n\n')
        cmd = 'bbduk.sh in1={f_in} in2={r_in} out1={f_out} out2={r_out} qtrim=w trimq=20 k=25 minlength=50 ' \
              'forcetrimleft=15 ref=adapters overwrite hdist=1 tpe tbo{optn}'\
            .format(f_in=forward_in,
                    r_in=reverse_in,
                    f_out=forward_out,
                    r_out=reverse_out,
                    optn=options)
    elif reverse_in == 'NA':
        cmd = 'bbduk.sh in={f_in} out={f_out} qtrim=w trimq=20 k=25 minlength=50 forcetrimleft=15' \
              ' ref=adapters overwrite hdist=1 tpe tbo{optn}'\
            .format(f_in=forward_in,
                    f_out=forward_out,
                    optn=options)
    else:
        if reverse_out == 'NA':
            raise ValueError('Reverse output reads must be specified.')
        cmd = 'bbduk.sh in1={f_in} in2={r_in} out1={f_out} out2={r_out} qtrim=w trimq=20 k=25 minlength=50 ' \
              'forcetrimleft=15 ref=adapters overwrite hdist=1 tpe tbo{optn}'\
            .format(f_in=forward_in,
                    r_in=reverse_in,
                    f_out=forward_out,
                    r_out=reverse_out,
                    optn=options)
    out, err = run_subprocess(cmd)
    if returncmd:
        return out, err, cmd
    else:
        return out, err


def tadpole(forward_in, forward_out, reverse_in='NA', returncmd=False, reverse_out='NA', mode='correct', **kwargs):
    """
    Runs tadpole. Default is to run in correction mode, but other modes ('contig', 'extend') can also be specified.
    :param forward_in: Forward input reads.
    :param forward_out: Forward output reads.
    :param returncmd: If set to true, function will return the cmd string passed to subprocess as a third value.
    :param reverse_in: Reverse reads. Only specify if not following _R1/_R2 convention/not in same folder as input.
    :param reverse_out: Reverse output reads. Automatically generated unless specified.
    :param mode: Mode to run tadpole in. Default is 'correct'.
    :param kwargs: Other arguments to give to tadpole in parameter='argument' format. See tadpole documentation for full list.
    :return: out and err: stdout string and stderr string from running tadpole.
    """
    options = kwargs_to_string(kwargs)
    if os.path.isfile(forward_in.replace('_R1', '_R2')) and reverse_in == 'NA' and '_R1' in forward_in:
        reverse_in = forward_in.replace('_R1', '_R2')
        if reverse_out == 'NA':
            if '_R1' in forward_out:
                reverse_out = forward_out.replace('_R1', '_R2')
            else:
                raise ValueError('If you do not specify reverse_out, forward_out must contain _R1.\n\n')
        cmd = 'tadpole.sh in1={} in2={} out1={} out2={} mode={} {}'.format(forward_in, reverse_in,
                                                                           forward_out, reverse_out,
                                                                           mode, options)
    elif reverse_in == 'NA':
        cmd = 'tadpole.sh in={} out={} mode={} {}'.format(forward_in, forward_out, mode, options)
    else:
        if reverse_out == 'NA':
            raise ValueError('Reverse output reads must be specified.')
        cmd = 'tadpole.sh in1={} in2={} out1={} out2={} mode={} {}'.format(forward_in, reverse_in,
                                                                           forward_out, reverse_out,
                                                                           mode, options)
    if not os.path.isfile(forward_out):
        out, err = run_subprocess(cmd)
    else:
        out = str()
        err = str()
    if returncmd:
        return out, err, cmd
    else:
        return out, err


def bbnorm(forward_in, forward_out, returncmd=False, reverse_in='NA', reverse_out='NA', **kwargs):
    """
    Runs bbnorm to normalize read depth. Default target kmer depth is left at bbnorm's default, which is 100.
    :param forward_in: Forward input reads.
    :param forward_out: Forward output reads.
    :param returncmd: If set to true, function will return the cmd string passed to subprocess as a third value.
    :param reverse_in: Reverse reads. Only specify if not following _R1/_R2 convention/not in same folder as input.
    :param reverse_out: Reverse output reads. Automatically generated unless specified.
    :param kwargs: Other arguments to give to bbnorm in parameter='argument' format. See bbnorm documentation for full list.
    :return: out and err: stdout string and stderr string from running bbnorm.
    """
    options = kwargs_to_string(kwargs)
    if os.path.isfile(forward_in.replace('_R1', '_R2')) and reverse_in == 'NA' and '_R1' in forward_in:
        reverse_in = forward_in.replace('_R1', '_R2')
        if reverse_out == 'NA':
            if '_R1' in forward_out:
                reverse_out = forward_out.replace('_R1', '_R2')
            else:
                raise ValueError('If you do not specify reverse_out, forward_out must contain _R1.\n\n')
        cmd = 'bbnorm.sh in1={} in2={} out={} out2={} {}'.format(forward_in, reverse_in,
                                                                  forward_out, reverse_out,
                                                                  options)
    elif reverse_in == 'NA':
        cmd = 'bbnorm.sh in={} out={} {}'.format(forward_in, forward_out, options)
    else:
        if reverse_out == 'NA':
            raise ValueError('Reverse output reads must be specified.')
        cmd = 'bbnorm.sh in1={} in2={} out1={} out2={} {}'.format(forward_in, reverse_in,
                                                                  forward_out, reverse_out,
                                                                  options)
    if not os.path.isfile(forward_out):
        out, err = run_subprocess(cmd)
    else:
        out = str()
        err = str()
    if returncmd:
        return out, err, cmd
    else:
        return out, err


def bbmerge(forward_in, merged_reads, returncmd=False, reverse_in='NA', **kwargs):
    """
    Runs bbmerge.
    :param forward_in: Forward input reads. Reverse reads automatically detected if present in the same folder.
    :param merged_reads: Output file to write merged reads to.
    :param returncmd: If set to true, function will return the cmd string passed to subprocess as a third value.
    :param reverse_in: Reverse input file, if you don't want it autodetected.
    :param kwargs: Other arguments to give to bbmerge in parameter='argument' format. See bbmerge documentation for full list.
    :return: out and err: stdout string and stderr string from running bbmerge.
    """
    options = kwargs_to_string(kwargs)
    if os.path.isfile(forward_in.replace('_R1', '_R2')) and reverse_in == 'NA' and '_R1' in forward_in:
        reverse_in = forward_in.replace('_R1', '_R2')
        cmd = 'bbmerge.sh in={} in2={} out={} {}'.format(forward_in, reverse_in, merged_reads, options)
    elif reverse_in == 'NA':
        cmd = 'bbmerge.sh in={} out={} {}'.format(forward_in, merged_reads, options)
    else:
        cmd = 'bbmerge.sh in={} in2={} out={} {}'.format(forward_in, reverse_in, merged_reads, options)
    if not os.path.isfile(merged_reads):
        out, err = run_subprocess(cmd)
    else:
        out = str()
        err = str()
    if returncmd:
        return out, err, cmd
    else:
        return out, err


def bbduk_bait(reference, forward_in, forward_out, returncmd=False, reverse_in='NA', reverse_out='NA', **kwargs):
    """
    Uses bbduk to bait out reads that have kmers matching to a reference.
    :param reference: Reference you want to pull reads out for. Should be in fasta format.
    :param forward_in: Forward reads you want to quality trim.
    :param returncmd: If set to true, function will return the cmd string passed to subprocess as a third value.
    :param forward_out: Output forward reads.
    :param reverse_in: Reverse input reads. Don't need to be specified if _R1/_R2 naming convention is used.
    :param reverse_out: Reverse output reads. Don't need to be specified if _R1/_R2 convention is used.
    :param kwargs: Other arguments to give to bbduk in parameter=argument format. See bbduk documentation for full list.
    :return: out and err: stdout string and stderr string from running bbduk.
    """
    options = kwargs_to_string(kwargs)
    if os.path.isfile(forward_in.replace('_R1', '_R2')) and reverse_in == 'NA' and '_R1' in forward_in:
        reverse_in = forward_in.replace('_R1', '_R2')
        if reverse_out == 'NA':
            if '_R1' in forward_out:
                reverse_out = forward_out.replace('_R1', '_R2')
            else:
                raise ValueError('If you do not specify reverse_out, forward_out must contain _R1.\n\n')
        cmd = 'bbduk.sh in={} in2={} outm={} outm2={} ref={}{}'.format(forward_in, reverse_in,
                                                                       forward_out, reverse_out,
                                                                       reference, options)
    elif reverse_in == 'NA':
        cmd = 'bbduk.sh in={} outm={} ref={}{}'.format(forward_in, forward_out, reference, options)
    else:
        if reverse_out == 'NA':
            raise ValueError('Reverse output reads must be specified.')
        cmd = 'bbduk.sh in={} in2={} outm={} outm2={} ref={}{}'.format(forward_in, reverse_in,
                                                                       forward_out, reverse_out,
                                                                       reference, options)
    out, err = run_subprocess(cmd)
    if returncmd:
        return out, err, cmd
    else:
        return out, err


def bbduk_filter(reference, forward_in, forward_out, returncmd=False, reverse_in='NA', reverse_out='NA', **kwargs):
    """
    Uses bbduk to filter out reads that have kmers matching to a reference.
    :param reference: Reference you want to pull reads out for. Should be in fasta format.
    :param forward_in: Forward reads you want to quality trim.
    :param returncmd: If set to true, function will return the cmd string passed to subprocess as a third value.
    :param forward_out: Output forward reads.
    :param reverse_in: Reverse input reads. Don't need to be specified if _R1/_R2 naming convention is used.
    :param reverse_out: Reverse output reads. Don't need to be specified if _R1/_R2 convention is used.
    :param kwargs: Other arguments to give to bbduk in parameter=argument format. See bbduk documentation for full list.
    :return: out and err: stdout string and stderr string from running bbduk.
    """
    options = kwargs_to_string(kwargs)
    if os.path.isfile(forward_in.replace('_R1', '_R2')) and reverse_in == 'NA' and '_R1' in forward_in:
        reverse_in = forward_in.replace('_R1', '_R2')
        if reverse_out == 'NA':
            if '_R1' in forward_out:
                reverse_out = forward_out.replace('_R1', '_R2')
            else:
                raise ValueError('If you do not specify reverse_out, forward_out must contain _R1.\n\n')
        cmd = 'bbduk.sh in={} in2={} out={} out2={} ref={}{}'.format(forward_in, reverse_in,
                                                                     forward_out, reverse_out,
                                                                     reference, options)
    elif reverse_in == 'NA':
        cmd = 'bbduk.sh in={} out={} ref={}{}'.format(forward_in, forward_out, reference, options)
    else:
        if reverse_out == 'NA':
            raise ValueError('Reverse output reads must be specified.')
        cmd = 'bbduk.sh in={} in2={} out={} out2={} ref={}{}'.format(forward_in, reverse_in,
                                                                     forward_out, reverse_out,
                                                                     reference, options)
    out, err = run_subprocess(cmd)
    if returncmd:
        return out, err, cmd
    else:
        return out, err


def dedupe(input_file, output_file, returncmd=False, **kwargs):
    """
    Runs dedupe from the bbtools package.
    :param input_file: Input file.
    :param returncmd: If set to true, function will return the cmd string passed to subprocess as a third value.
    :param output_file: Output file.
    :param kwargs: Arguments to give to dedupe in parameter=argument format. See dedupe documentation for full list.
    :return: out and err: stdout string and stderr string from running dedupe.
    """
    options = kwargs_to_string(kwargs)
    cmd = 'dedupe.sh in={} out={}{}'.format(input_file, output_file, options)
    out, err = run_subprocess(cmd)
    if returncmd:
        return out, err, cmd
    else:
        return out, err


def seal(reference, forward_in, output_file, reverse_in='NA', returncmd=False, **kwargs):
    """
    Runs seal from the bbtools package.
    :param reference: Reference file, in fasta format.
    :param returncmd: If set to true, function will return the cmd string passed to subprocess as a third value.
    :param forward_in: Forward reads, fastq format.
    :param output_file: Output file to put rpkm statistics into.
    :param reverse_in: Reverse reads. Not necessary to specify if in same folder and follow _R1/_R2 convention.
    :param kwargs: Arguments to give to seal in parameter=argument format. See seal documentation for full list.
    :return: out and err: stdout string and stderr string from running seal.
    """
    options = kwargs_to_string(kwargs)
    if os.path.isfile(forward_in.replace('_R1', '_R2')) and reverse_in == 'NA' and '_R1' in forward_in:
        reverse_in = forward_in.replace('_R1', '_R2')
        cmd = 'seal.sh ref={} in={} in2={} rpkm={} nodisk{}'.format(reference, forward_in, reverse_in, output_file, options)
    elif reverse_in == 'NA':
        cmd = 'seal.sh ref={} in={} rpkm={} nodisk{}'.format(reference, forward_in, output_file, options)
    else:
        cmd = 'seal.sh ref={} in={} in2={} rpkm={} nodisk{}'.format(reference, forward_in, reverse_in, output_file, options)
    out, err = run_subprocess(cmd)
    if returncmd:
        return out, err, cmd
    else:
        return out, err


def kmercountexact(forward_in, reverse_in='NA', returncmd=False, **kwargs):
    """
    Wrapper for kmer count exact.
    :param forward_in: Forward input reads.
    :param reverse_in: Reverse input reads. Found automatically for certain conventions.
    :param returncmd: If set to true, function will return the cmd string passed to subprocess as a third value.
    :param kwargs: Arguments to give to kmercountexact in parameter='argument' format.
    See kmercountexact documentation for full list.
    :return: out and err: stdout string and stderr string from running kmercountexact.
    """
    options = kwargs_to_string(kwargs)
    if os.path.isfile(forward_in.replace('_R1', '_R2')) and reverse_in == 'NA' and '_R1' in forward_in:
        reverse_in = forward_in.replace('_R1', '_R2')
        cmd = 'kmercountexact.sh in={} in2={} {}'.format(forward_in, reverse_in, options)
    elif reverse_in == 'NA':
        cmd = 'kmercountexact.sh in={} {}'.format(forward_in, options)
    else:
        cmd = 'kmercountexact.sh in={} in2={} {}'.format(forward_in, reverse_in, options)
    out, err = run_subprocess(cmd)
    if returncmd:
        return out, err, cmd
    else:
        return out, err


def genome_size(peaks_file, haploid=True):
    """
    Finds the genome size of an organsim, based on the peaks file created by kmercountexact.sh
    :param peaks_file: Path to peaks file created by kmercountexact.
    :param haploid: Set to True if organism of interest is haploid, False if not. Default True.
    :return: size of genome, as an int. If size could not be found, return will be 0.
    """
    size = 0
    with open(peaks_file) as peaks:
        lines = peaks.readlines()
    for line in lines:
        if haploid:
            if '#haploid_genome_size' in line:
                size = int(line.split()[1])
        else:
            if '#genome_size' in line:
                size = int(line.split()[1])
    return size


def subsample_reads(forward_in, forward_out, num_bases, returncmd=False, reverse_in='NA', reverse_out='NA',
                    **kwargs):
    options = kwargs_to_string(kwargs)
    if os.path.isfile(forward_in.replace('_R1', '_R2')) and reverse_in == 'NA' and '_R1' in forward_in:
        reverse_in = forward_in.replace('_R1', '_R2')
        if reverse_out == 'NA':
            if '_R1' in forward_out:
                reverse_out = forward_out.replace('_R1', '_R2')
            else:
                raise ValueError('If you do not specify reverse_out, forward_out must contain _R1.\n\n')
        cmd = 'reformat.sh in1={} in2={} out1={} out2={} samplebasestarget={} {}'.format(forward_in, reverse_in,
                                                                                         forward_out, reverse_out,
                                                                                         str(num_bases), options)
    elif reverse_in == 'NA':
        cmd = 'reformat.sh in={} out={} samplebasestarget={} {}'.format(forward_in, forward_out,
                                                                       str(num_bases), options)
    else:
        if reverse_out == 'NA':
            raise ValueError('Reverse output reads must be specified.')
        cmd = 'reformat.sh in1={} in2={} out1={} out2={} samplebasestarget={} {}'.format(forward_in, reverse_in,
                                                                                         forward_out, reverse_out,
                                                                                         str(num_bases), options)
    if not os.path.isfile(forward_out):
        out, err = run_subprocess(cmd)
    else:
        out = str()
        err = str()
    if returncmd:
        return out, err, cmd
    else:
        return out, err


def validate_reads(forward_in, returncmd=False, reverse_in='NA'):
    if os.path.isfile(forward_in.replace('_R1', '_R2')) and reverse_in == 'NA' and '_R1' in forward_in:
        reverse_in = forward_in.replace('_R1', '_R2')
        cmd = 'reformat.sh in1={} in2={} vpair'.format(forward_in, reverse_in)
    elif reverse_in == 'NA':
        cmd = 'reformat.sh in={}'.format(forward_in)
    out, err = run_subprocess(cmd)
    if returncmd:
        return out, err, cmd
    else:
        return out, err


def reformat_reads(forward_in, forward_out, returncmd=False, reverse_in='NA', reverse_out='NA'):
    if os.path.isfile(forward_in.replace('_R1', '_R2')) and reverse_in == 'NA' and '_R1' in forward_in:
        reverse_in = forward_in.replace('_R1', '_R2')
        if reverse_out == 'NA':
            if '_R1' in forward_out:
                reverse_out = forward_out.replace('_R1', '_R2')
            else:
                raise ValueError('If you do not specify reverse_out, forward_out must contain _R1.\n\n')
        cmd = 'reformat.sh in1={} in2={} out1={} out2={} tossbrokenreads=t ow=t'\
            .format(forward_in, reverse_in, forward_out, reverse_out)
    elif reverse_in == 'NA':
        cmd = 'reformat.sh in={} out={} tossbrokenreads=t ow=t'.format(forward_in, forward_out)
    else:
        if reverse_out == 'NA':
            raise ValueError('Reverse output reads must be specified.')
        cmd = 'reformat.sh in1={} in2={} out1={} out2={} tossbrokenreads=t ow=t'\
            .format(forward_in, reverse_in, forward_out, reverse_out)
    if not os.path.isfile(forward_out):
        out, err = run_subprocess(cmd)
    else:
        out = str()
        err = str()
    if returncmd:
        return out, err, cmd
    else:
        return out, err


def repair_reads(forward_in, forward_out, returncmd=False, reverse_in='NA', reverse_out='NA'):
    if os.path.isfile(forward_in.replace('_R1', '_R2')) and reverse_in == 'NA' and '_R1' in forward_in:
        reverse_in = forward_in.replace('_R1', '_R2')
        if reverse_out == 'NA':
            if '_R1' in forward_out:
                reverse_out = forward_out.replace('_R1', '_R2')
            else:
                raise ValueError('If you do not specify reverse_out, forward_out must contain _R1.\n\n')
        cmd = 'repair.sh in1={} in2={} out1={} out2={} tossbrokenreads=t repair=t overwrite=t'\
            .format(forward_in, reverse_in, forward_out, reverse_out)
    else:
        if reverse_out == 'NA':
            raise ValueError('Reverse output reads must be specified.')
        cmd = 'repair.sh in1={} in2={} out1={} out2={} tossbrokenreads=t repair=t overwrite=t'\
            .format(forward_in, reverse_in, forward_out, reverse_out)
    if not os.path.isfile(forward_out):
        out, err = run_subprocess(cmd)
    else:
        out = str()
        err = str()
    if returncmd:
        return out, err, cmd
    else:
        return out, err


