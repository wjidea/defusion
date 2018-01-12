import os, logging, subprocess, shlex

def which(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None


def bioperl_loaded():
    cmd1 = 'perl -MBio::Seq -e 0'
    cmd1L = shlex.split(cmd1)
    p1 = subprocess.Popen(cmd1L, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p1.communicate()
    # print('BioPerl', stdout, stderr)
    if not stderr:
        print("bioPerl not installed")
        return(None)
    else:
        return(True)


def fix_path_slash(path):
    return (os.path.join(os.path.abspath(path), ''))


def log_err(task, seqID, stdout, stderr):
    logging.warning('stderr: {0}'.format(stderr))
    logging.debug('stdout: {0}'.format(stdout))
    logging.info('Finish {0} {1}'.format(task, seqID))


def parse_SeqID(seqID):
    # seqID: Chr10_12968107_12975977
    # prepare the coordinates
    seqID_list = seqID.rstrip().split('_')
    
    end = seqID_list.pop()
    start = seqID_list.pop()
    seqid = '_'.join(seqID_list)
    return (seqid, int(start), int(end))


if __name__ == "__main__":
    print("load utility functions")
    bioperl_loaded()