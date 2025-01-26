"""
Author: Jiazheng Xie.
Github: https://github.com/xiebio
"""
import gzip
import re
def parse_col_file(file,info_col=1,sep='\s+'):
    dict_col = {}
    with (
        gzip.open(file, "rt", encoding="utf8") if file.endswith(".gz")
        else open(file, "rt", encoding="utf8")
    ) as fin:
        for line in fin:
            if line and not line.startswith("#"):
                i = line.rstrip()
                i = re.split(sep,i)
                dict_col[ i[0] ] = i[info_col]
    return dict_col
