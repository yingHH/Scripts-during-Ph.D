"""
Esaily download NCBI data using Aspera

@author: yinghuang
"""
# self pkgs
from pkgkit import runshell
import os


asp_id = '~/.asp_id/asperaweb_id_dsa.openssh'
ascp_path = '/home/yingh/.aspera/connect/bin/ascp'
odir = './'

def check_dir(directory):
    if not os.path.isdir(directory):
        os.mkdir(directory)
    return os.path.abspath(directory)

def easy_asp(data_id, odir=odir, asp_id=asp_id, ascp_path=ascp_path):
    """
    data_id: input a GEO data ID that start with "SRR".
    odir: ouput directory, default "./".
    asp_id: the path to ".openssh" file, default "~/.asp_id/asperaweb_id_dsa.openssh".
    ascp_path: the path to ascp script, default "/home/yingh/.aspera/connect/bin/ascp".
    """

    odir = check_dir(odir)

    id_type = data_id[:3]
    data_pdir = data_id[:6]
    cmd = '{ascp_path} \
        -i {asp_id} \
        -T \
        -l 200m \
        anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/{id_type}/{data_pdir}/{data_id}/{data_id}.sra \
        {odir}'.format(ascp_path=ascp_path, asp_id=asp_id, id_type=id_type, data_pdir=data_pdir, data_id=data_id, odir=odir)

    print('[Running]: ', cmd)
    out = runshell(cmd)
    return out
