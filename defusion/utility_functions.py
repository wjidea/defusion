import os, logging, subprocess, shlex, sys


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
        logging.error("bioPerl not installed")
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


def check_coord(prefix, seq_id, start, end):
    
    gene_interval = int(end) - int(start)
    if gene_interval > 100000:
        fo = open(prefix + "/wrong_coordinates.err", 'ab')
        fo.write("{}\t{}\t{}\t{}\n".format(seq_id, start, end, gene_interval))
        logging.error("gene range greater than 100kb detected, please check wrong_coordinates.txt")
        fo.close()
        return(None)
    
    return(True)


def coord_error(error):
    if error:
        cont = raw_input("coords errors detected, do you want to ignore and move on? (y/n)")
    
        if cont == "y" or cont == "Y":
            logging.error("ignore coords error and move on")
        elif cont == "n" or cont == "N":
            logging.error("coord error detected, program terminated and please fix coordnated error acoording to "
                          "prefix/wrong_coordinates.err")
            sys.exit()
        else:
            logging.error(" Please type 'y' or 'n' ")
            coord_error(True)


def input_validate(in_list):
    missing_files = []
    file_valid_list = [os.path.isfile(file) for file in in_list]
    
    if all(file_valid_list):
        return(True)
    else:
        for idx, bl in enumerate(file_valid_list):
            if not bl:
                missing_files.append(in_list[idx])
                logging.error("missing input files {}".format(missing_files))
                sys.exit()
        # return missing_files


if __name__ == "__main__":
    logging.info("load utility functions")
    bioperl_loaded()