#!/share/pkg/python/3.6.4/install/bin/python3.6


# This script is for taking a text file of bash commands and running each line
# as a seperate arrayjob.
# to be run on the SCC (shared computing cluster)
# Written for python3 runs on

# composed by Jacob Pessin <jpessin@bu.edu>
# for help contact Research Computing Services help@scc.bu.edu


"""Create a qsub script, submit it to qsub, if -concat create a second job to concat results."""

import argparse
import re
import subprocess
import os  # , sys

#from __future__ import print_function as print, division as with_ as

########################
# submitter
########################

def qsubmit(qsub_filename):
    """Submit a job with qsub, return what qsub returns."""
    proc = subprocess.run(
        ['qsub', qsub_filename],
        stdout=subprocess.PIPE
    )
    return proc

########################
# qsub templates
########################


SCC_QSUB_PARAMS = """
#!/bin/bash -l
#$ -N {name}
#$ -P {project}
{emailaddy_line}
{emailwhen_line}
#$ -j {join_outputs}
#$ -V
#$ -cwd
{run_time}
#$ -pe omp {ncores}
{mem_per_core}
{array_arg_t_line}
{array_arg_tc_line}

{module_line}

{array_cmd_line}
"""


SCC_CONCAT = """
#$ -N {name}_concat
#$ -P {project}
{emailaddy_line}
{emailwhen_line}
#$ -j {join_outputs}
#$ -V
#$ -cwd
#$ -pe omp 1
#$ -l mem_per_core=1g
{hold_jid_line}

{concat_out_line}
{concat_err_line}

cp $TMPDIR/*concat* .
"""

####################
# fill out templates
####################


def generate_scc_qsub_array(results, arg_t_line, array_cmd_line):
    """Write a scc SGE array script."""

    def module_format(contain):
        """Format module load line."""
        # puting this hear is a functional, yet ugly hack
        modcount = len(contain)
        if modcount > 1:
            contain = " ".join(contain)
            return 'module load {}'.format(contain)
        elif modcount == 1:
            contain = contain[0]
            return 'module load {}'.format(contain)
        else:
            return ''


    fname = results.name + "_array.qsub"
    if not os.path.exists(fname):
        with open(fname, 'wt') as qsub_array:
            qsub_array.write(
                SCC_QSUB_PARAMS.format(
                    name=results.name,
                    project=results.project,
                    emailaddy_line=results.emailaddy_line,
                    emailwhen_line=results.emailwhen_line,
                    join_outputs=results.join_stdoe,
                    run_time=results.run_time,
                    ncores=results.ncores,
                    mem_per_core=results.mem_per_core,
                    array_arg_t_line=arg_t_line,
                    array_arg_tc_line=results.array_arg_tc_line,

                    module_line=module_format(results.module_line),
                    #~ module_line=results.module_line,

                    array_cmd_line=array_cmd_line
                )
            )
    else:
        msg = "{} already exists, change some names and try again".format(fname)
        raise FileExistsError(msg)
    return fname


def generate_scc_qsub_concat(results, concat_out_line, concat_err_line):
    """Write a scc SGE script to concatanate an output from the array script."""
    fname = results.name + "_concat.qsub"
    with open(fname, "wt") as qsub_concat:
        qsub_concat.write(
            SCC_CONCAT.format(
                name=results.name,
                project=results.project,
                emailaddy_line=results.emailaddy_line,
                emailwhen_line=results.emailwhen_line,
                join_outputs=results.join_stdoe,
                mem_per_core=results.mem_per_core,
                hold_jid_line='#$ -hold_jid {}'.format(results.name),

                concat_out_line=concat_out_line,
                concat_err_line=concat_err_line
            )
        )
    return fname


def form_array_submission(results):
    """Parse jobsfile, create <name>_array commands, return Formated arg_t line."""
    def comment_or_empty(string):
        if string.startswith('#'):
            return True
        if string in ("", " ", "\n", " \n"):
            return True
        return False

    with open(results.jobs_file, 'rt') as jobsfile:
        ingst = jobsfile.readlines()
        runlines = [line for line in ingst if not comment_or_empty(line)]
    array_cmds_name = os.path.join(os.getcwd(), results.name + "_array_commands.txt")
    if not os.path.exists(array_cmds_name):
        with open(array_cmds_name, "wt") as cmdtxt:
            cmdtxt.writelines(runlines)
    return array_cmds_name, "#$ -t 1:{}".format(str(len(runlines)))


def form_array_cmd_line(array_cmds_name):
    """Form sed read line in file for TASK_ID."""
    return 'sed -n -e "$SGE_TASK_ID p" {} | bash'.format(array_cmds_name)


def form_hold_jid_line(results):
    """format line for qsub to wait dependently on previous job by name."""
    return '#$ -hold_jid {}'.format(results.name)


def form_concat_out_line(results):
    """Format concat_output line."""
    concat_out_line = 'cat {name}.o* >> $TMPDIR/{name}_concat_stdout'.format(name=results.name)
    return concat_out_line


def form_concat_err_line(results):
    """Format concat_error line."""
    if results.join_stdoe in ['n', 'no']:
        concat_err_line = 'cat {name}.e* >> $TMPDIR/{name}_concat_stderr'.format(name=results.name)
    else:
        concat_err_line = ''
    return concat_err_line


##############################
# argpasre: validate and format
###############################


def valid_job_name(string):
    """Validate the name string for job inputs for the qsub script."""
    stlen = '{' + str(len(string)) + '}'

    def match(string):
        return re.match('[a-zA-Z0-9_-]{}'.format(stlen), string)

    if not string.isprintable() or not match(string):
        msg = "Job names can only contain a-zA-Z0-9 underscore and hyphen"
        raise argparse.ArgumentTypeError(msg)
    return string


def valid_project_name(string):
    """Validate allow character types in a project name."""
    stlen = '{' + str(len(string) - 1) + '}'
    if not re.match('[a-z][a-zA-Z0-9_-]{}'.format(stlen), string):
        msg = "Project string is not in the valid format for a project name"
        raise argparse.ArgumentTypeError(msg)
    return string


def form_emailaddy_line(string):
    """generate the line to tell qsub where to send notification emails."""
    if string:
        return "#$ -M {}".format(string)
    return string


# JP notewhere: validating an email address? Where will it end?
# def form_validate_emailaddy(string):
    # if not re.match([a-zA-Z0-9._-]+@[a-zA-Z0-9.-]+\.[a-zA-Z]*, string):
        # msg = '"{}" does not appear to be a valid email address'.format(string)
        # raise raise argparse.ArgumentTypeError(msg)
    # return form_email_line(string)

def validate_form_emailwhen_line(string):
    """Validate the inputs for and format for qsub -m, when to send notify emails."""
    if string == '':
        return ''
    for char in string:
        if char not in 'beasn':
            msg = """Invalid option with optional -m argument, many only contain
            'b', 'e', 'a', 's', or 'n' """
            raise argparse.ArgumentTypeError(msg)
    return "#$ -m {}".format(string)


def validate_join_stoe(string):
    """Validate options for qsub -j."""
    if string not in {'y', 'yes', 'n', 'no'}:
        msg = """Invalid option with optional -j argument,
        must be one of "y", "yes", "n", "no" """
        raise argparse.ArgumentTypeError(msg)
    return string


def validate_form_runtime(string):
    """Validate runtime string for -l h_rt=<string>, as hh:mm:ss."""
    splt = string.split(':')
    msg = """-h_rt must be in the form hh:mm:ss,
        where hh is hours, mm is minutes, ss is seconds:
            minutes and seconds required"""
    if len(splt) > 3:
        raise argparse.ArgumentTypeError(msg)
    for time in splt:
        if not 1 <= len(time) <= 2:
            raise argparse.ArgumentTypeError(msg)
        elif not re.match('[0-9]{1,2}', time):
            raise argparse.ArgumentTypeError(msg)
    return "#$ -l h_rt={}".format(string)


def form_mem_per_core(string):
    """Format mem_per_core line."""
    #  if string not None, else return ""
    if not string:
        return ''
    stlen = str(len(string) - 1)
    if not string.endswith(('m', 'g')):
        msg = """Missing unit for mem_per_core, last char should be one of
            k, m or g"""
        raise argparse.ArgumentTypeError(msg)
    if not re.match('[0-9]{}'.format('{' + stlen + '}'), string):
        msg = """mem_per_core should be a whole number with a using ex:
            4g
            4096m"""
        raise argparse.ArgumentTypeError(msg)
    return "#$ -l mem_per_core={}".format(string)


def form_tc_line(string):
    """format line for qsub -tc."""
    if int(string) not in range(1, 501):
        msg = "Invalid concerting job limit '-tc', must be an integer between 1-500"
        raise argparse.ArgumentTypeError(msg)
    return "#$ -tc {}".format(string)


# When using argparses nargs, this runs on members, argparse returns list
#~ def format_module_line(contain):
    #~ """Format module load line."""
    #~ # puting this hear is a functional, yet ugly hack
    #~ modcount = len(contain)
    #~ if modcount > 1:
        #~ contain = " ".join(contain)
        #~ return 'module load {}'.format(contain)
    #~ elif modcount == 1:
        #~ contain = contain[0]
        #~ return 'module load {}'.format(contain)
    #~ else:
        #~ return ''


def validate_jobsfile(string):
    """Check jobsfile for valid path/readable, and is text using posix's file."""
    realpath = os.path.realpath(string)
    if not os.path.isfile(realpath):
        if not os.path.exists(realpath):
            msg = "file not found: {}".format(realpath)
        else:
            msg = "{} is not a file".format(realpath)
        raise argparse.ArgumentTypeError(msg)
    else:
        proc = subprocess.run(
            ['file', '-b', realpath],
            stdout=subprocess.PIPE
        )
        if 'text' not in str(proc.stdout):
            msg = "{} is not a text file."
            raise argparse.ArgumentTypeError(msg)
    return realpath


def argparser():
    """Argparser functions."""
    parser = argparse.ArgumentParser(description='''
Create an array job and dependent concatenation script for running on SCC \n
    Non-required option use system defaults if not specified
        ''')
    parser.add_argument('-N',
                        action='store',
                        dest='name',
                        required=True,
                        type=valid_job_name,
                        help='The name of your job.')
    parser.add_argument('-P',
                        action='store',
                        dest='project',
                        required=True,
                        type=valid_project_name,
                        help='Project Name to run the job under')
    parser.add_argument('-M',
                        action='store',
                        dest='emailaddy_line',
                        required=False,
                        default="",
                        type=form_emailaddy_line,
                        help='Prefered email for job notifications')
    parser.add_argument('-m',
                        action='store',
                        dest='emailwhen_line',
                        required=False,
                        default="",
                        type=validate_form_emailwhen_line,
                        help='''Email notification conditions,
                        any combination of b e a s n
                                b begins
                                e ends
                                a aborted
                                s suspended
                                n never''')
    parser.add_argument('-j',
                        action='store',
                        dest='join_stdoe',
                        required=False,
                        default='yes',
                        type=validate_join_stoe,
                        choices=('y', 'yes', 'n', 'no'),
                        help='''Join standard out and standard error,
                            y[es] or n[o], default is yes''')
# # JP: note not using args for -cwd or -V they were requested on by default,
# # no clear way to map the syntax of choosing -cwd or -V  in a qsub script and
# # have on by default
#    parser.add_argument('-V',
#                        action='store_true',
#                        dest='V_opt',
#                        required=False,
#                        help=????????)
#    parser.add_argument('-cwd',
#                        action='store_true',
#                        dest='cwd_opt',
#                        required=False,
#                        help= ??????)
# # JP: note choose -omp as it seems like it might be less confusing
# # than -pe as we are only mapping -pe omp $ncore
    parser.add_argument('-h_rt',
                        action='store',
                        dest='run_time',
                        default='11:59:59',
                        type=validate_form_runtime,
                        help='''The time maximum runtime for each individual job.
                            Format: hh:mm:ss
                            Default 12 hrs''')
    parser.add_argument('-omp',
                        action='store',
                        dest='ncores',
                        required=False,
                        choices=['1', '2', '3', '4', '8', '16', '20', '28'],
                        default=1,
                        help='''Number of processor cores for individual jobs,
                                if the job uses more than this it will be terminated early''')
    parser.add_argument('-mem_per_core',
                        action='store',
                        dest='mem_per_core',
                        required=False,
                        type=form_mem_per_core,
                        default="",
                        help='''Memory provided per core per job,
                            ex with "-omp 4 -mem_per_core 4g"
                            each line in the jobfile will run with 16 Gigabytes of available RAM''')
    parser.add_argument('-tc',
                        action='store',
                        dest='array_arg_tc_line',
                        required=False,
                        # possible issue mixing type with I/O as
                        # string and range with numbers
                        type=form_tc_line,
                        default="200",
                        help='''Maximum number of concurrently running, independent jobs,
                            any number from 1 and 500, default 200.''')
    parser.add_argument('-modload',
                        nargs='*',
                        dest='module_line',
                        required=False,
                        #~ type=format_module_line,
                        default="",
                        help='''Module(s) to load during jobs,
                                Can have multiple module listed use:
                                -modload "python3 R/3.5.2"
                                not
                                -modload python3 -modload R/3.5.2''')
    parser.add_argument('-concat',
                        action='store_true',
                        dest='concat',
                        required=False,
                        help='''If uses, concatenate out and error to final files
                        as jobname_cat.out and jobname_cat.err
                        if "-j y" the combined file will be jobname_concat.out'''
                        )
    parser.add_argument('-jobsfile',
                        action='store',
                        dest='jobs_file',
                        type=validate_jobsfile,
                        required=True,
                        help='Path to the file to generate the jobs')

    results = parser.parse_args()

    return results


def main():
    """Do the thing."""
    results = argparser()
    array_cmds_name, arg_t_line = form_array_submission(results)
    array_cmd_line = form_array_cmd_line(array_cmds_name)
    run_name = generate_scc_qsub_array(results, arg_t_line, array_cmd_line)
    #~ run_qsub_stdoe = qsubmit(run_name)
    # print(run_qsub_stdoe)

    if results.concat:
        concat_out_line = form_concat_out_line(results)
        concat_err_line = form_concat_err_line(results)
        concat_name = generate_scc_qsub_concat(
            results, concat_out_line, concat_err_line)
        #~ concat_qsub_stdoe = qsubmit(concat_name)
        # print(concat_qsub_stdoe)
    # parse results and print jobnames, job id,
    # TODO JP: print for user:
    # to see all your running jobs
    # to see these jobs when running
    # to see all your finished jobs,
    # to see these finished jobs


#~ def main():
    #~ """For argparse testing : just print argparse results namespace"""
    #~ results = argparser()
    #~ print(results)

if __name__ == '__main__':
    main()
