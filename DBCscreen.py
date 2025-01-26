#!/usr/bin/env python3
"""
Author: Jiazheng Xie.
Github: https://github.com/xiebio
"""
import getopt
import sys
import os
import gzip
import re
program_dir = os.path.dirname(sys.argv[0])
sys.path.append(f"{program_dir}/scripts")
from util import parse_col_file

def parse_config (config_file):
    parameters_dict = {
        'run_gx':'',#The path of the executable file fcs-gx
        'gxdb':'',#DBCsdb path
    }
    with open (config_file) as fh:
        for line in fh.readlines():
            if line.startswith('#'):continue
            line = line.strip()
            name,path = line.split('=')
            if name in parameters_dict:
                if os.path.isabs(path):
                    parameters_dict[name] = path
                else:
                    parameters_dict[name] = os.path.join(config_dir,path)
    return parameters_dict


def usage():
    '''
    help
    '''
    print("Usage: DBCscreen.py -h help -c config_file -i input_file -o out_dir")
    print("-i --input_file two columns file: fasta taxid")
    print("-o --out_dir outdir prefix[default:./]")


if __name__ == "__main__":
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hc:i:o:",
                ["help", "config_file=", "input_file=","out_dir=",])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    config_file = program_dir+'/DBCscreen.config'
    outdir = './'
    infas_list = ''
    
    for o, a in opts:
        if o in ('-h', '--help'):
            usage()
            sys.exit()
        elif o in ("-c", "--config_file"):
            config_file = os.path.abspath(a)
        elif o in ("-i", "--input_file"):
            infas_list = os.path.abspath(a)
        elif o in ("-o", "--out_dir"):
            outdir = os.path.abspath(a)
    if not infas_list:
        print('-i/--input_file must be set!')
        usage()
        sys.exit(2)

    config_dir = os.path.dirname(config_file)
    params_dict = parse_config(config_file)
    outdir_list = [
        outdir,
        outdir+'/0_gx',
        outdir+'/1_gx_filter',
            ]
    for i in outdir_list:
        if not os.path.exists(i):os.mkdir(i)
    list_info = parse_col_file(infas_list)
    script = ''
    #--------------------------------  0.gx  -----------------------------------
    for infas,taxid in list_info.items():
        id_name = os.path.basename(infas)
        id_name = re.sub(r"\.\w+$", "", id_name)
        if not os.path.exists(f'{outdir}/0_gx/{id_name}.{taxid}.taxonomy.rpt'):
            script += \
'''echo {1} gx begin
{4} --action-report 0 --fasta={2} --tax-id={3} --gx-db={5} --out-dir={0}/0_gx
echo {1} gx done
'''.format(outdir,id_name,infas,taxid,params_dict['run_gx'],params_dict['gxdb'])
    #-------------------------------- 1.gx_filter -----------------------------------
    script += \
'''echo gx_filter begin
python {1}/scripts/gx_filter.py {2} {0}/0_gx {0}/1_gx_filter {3} {4}
echo gx_filter done
'''.format(outdir,program_dir,infas_list,params_dict['run_gx'],params_dict['gxdb'])
    #-------------------------------- Write and run script ---------------------------
    fout = open(outdir+'/run_DBCscreen.sh','w')
    fout.write(script)
    fout.close()
    os.system("sh {}/run_DBCscreen.sh".format(outdir))
