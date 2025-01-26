#!/usr/bin/env python3
"""
Author: Jiazheng Xie.
Github: https://github.com/xiebio
"""
import sys
from sys import argv
import os
from Bio import SeqIO
import re
import gzip
from collections import defaultdict
from util import parse_col_file

def main():
    infas_list,inrpt_dir,outdir,run_gx,gxdb = argv[1:]
    list_info = parse_col_file(infas_list)
    sys.path.append(f"{os.path.dirname(run_gx)}/../scripts")
    import classify_taxonomy
    gi2div = parse_col_file(f'{gxdb}/db.blast_div.tsv.gz',sep='\t')
    label_virus = defaultdict(bool)
    for infas,intaxid in list_info.items():
        assembly_div = gi2div[intaxid]
        id_name = os.path.basename(infas)
        id_name = re.sub(r"\.\w+$", "", id_name)
        print(f'gx_filter: {id_name}')
        contaminated_seq_id = []
        primary_split_seq_id = []
        split_seq_id_cvg = defaultdict(int)
        inrpt = f'{inrpt_dir}/{id_name}.{intaxid}.taxonomy.rpt'
        fout = open(f"{outdir}/{id_name}.filter.rpt", "wt", encoding="utf8")
        assert os.path.exists( inrpt ),f'Caution: {inrpt} does not exist!'

        for row in classify_taxonomy.Record.get_rows(inrpt):
            if isinstance(row, str) and row.startswith("#"):
                fout.write(row)
                continue
            r = classify_taxonomy.Record.from_row(row)
            contig_acc = r.seq_id.split("~")[0]
            if not r.taxa: continue
            marker_split = "~" in r.seq_id
            if marker_split and (contig_acc in primary_split_seq_id): continue
            t0 = r.taxa[0]            
            if t0.score < 40:
                continue
            if t0.div == assembly_div:
                if marker_split: primary_split_seq_id.append(contig_acc)
                continue
            taxa = [t for t in r.taxa if (t.score >40) and t.score * t.cvg_by_div > t0.score * t0.cvg_by_div * 0.8]
            divs = [t.div for t in taxa]
            marker_primary = False
            for div_gx in divs:
                if ( div_gx == assembly_div):
                    marker_primary = True
                    break
            if marker_primary:
                if marker_split: primary_split_seq_id.append(contig_acc)
                continue

            if (t0.cvg_by_div/r.len > 0.001) and (t0.cvg_by_div > 300 or t0.score * t0.score > 15 * t0.cvg_by_div):
                print(*row, sep="\t", file=fout)
                contaminated_seq_id.append(contig_acc)
                if marker_split:
                    split_seq_id_cvg[contig_acc] += t0.cvg_by_div
                label_virus[contig_acc] = True if t0.div.startswith('virs:') else False
        fout.close()
        
        if len(contaminated_seq_id) > 0:
            seq_records_vir = []
            seq_records_dbc = []
            with (
                gzip.open(infas, "rt", encoding="utf8") if infas.endswith(".gz")
                else open(infas, "rt", encoding="utf8")
            ) as fin:
                for record in SeqIO.parse(fin, "fasta"):
                    record_id = record.id.strip('|').split('|')[-1]
                    if (record_id in contaminated_seq_id) and (record_id not in primary_split_seq_id):
                        if record_id in split_seq_id_cvg and split_seq_id_cvg[contig_acc]/len(record.seq) < 0.001: #Coverage too low
                            continue
                        if label_virus[record_id]:
                            seq_records_vir.append(record)
                        else:seq_records_dbc.append(record)
            SeqIO.write(seq_records_vir, f"{outdir}/{id_name}.vir.fas", "fasta")
            SeqIO.write(seq_records_dbc, f"{outdir}/{id_name}.dbc.fas", "fasta")

if __name__ == '__main__':  
    main()
