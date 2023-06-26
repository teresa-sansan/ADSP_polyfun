### This is from Chirag
import sys
import pyBigWig
import pandas as pd
import numpy as np
import pyfasta
import pyranges as pr
import math
import yaml
import torch
from torch.nn import functional as F
from torch import nn
from pytorch_lightning.core.lightning import LightningModule
from axial_positional_embedding import AxialPositionalEmbedding

#from vectorize_sum import one_hot_C, rev_comp_C, merge_read_counts_C

#from numba import jit

import time
from datetime import date
import datetime

import math

import random
from torch.utils.checkpoint import checkpoint

import os


import statistics

import kipoiseq
from kipoiseq import Interval, Variant

import pyfaidx

#from tqdm import tqdm

from pytorch_lightning.callbacks import ModelSummary

from einops import rearrange, reduce, repeat

from pytorch_lightning.loggers import CSVLogger

from scipy.stats import pearsonr

from math import pi, log

import re
from kipoiseq.extractors import VariantSeqExtractor
from cyvcf2 import VCF
def anti_join(x, y, on):
    """Return rows in x which are not present in y"""
    ans = pd.merge(left=x, right=y, how='left', indicator=True, on=on)
    ans = ans.loc[ans._merge == 'left_only', :].drop(columns='_merge')
    return ans

#24490
idx_num = int(sys.argv[1]) - 1

fasta_file = '/sc/arion/work/lakhac01/microglia_dl/data/hg38.fa'

metadata = '/sc/arion/work/lakhac01/ADSP_reguloML/train_dl_models/unsupervised_pretraining/test_kipoiseq/1000_metadata.tsv'
sample_id_unrelated = '/sc/arion/projects/ad-omics/clakhani/1KG_phased/1000G_2504_high_coverage.sequence.index'
cols = ['ENA_FILE_PATH',
        'MD5SUM',
        'RUN_ID',
        'STUDY_ID',
        'STUDY_NAME',
        'CENTER_NAME',
        'SUBMISSION_ID',
        'SUBMISSION_DATE',
        'SAMPLE_ID',
        'SAMPLE_NAME',
        'POPULATION',
        'EXPERIMENT_ID',
        'INSTRUMENT_PLATFORM',
        'INSTRUMENT_MODEL',            
        'LIBRARY_NAME',
        'RUN_NAME',
        'INSERT_SIZE',
        'LIBRARY_LAYOUT',
        'PAIRED_FASTQ',
        'READ_COUNT',
        'BASE_COUNT',
        'ANALYSIS_GROUP'
        ]

environment = yaml.safe_load(open('../../../../environment.yml'))


genome = pyfasta.Fasta(fasta_file)


metadata = '/sc/arion/work/lakhac01/ADSP_reguloML/train_dl_models/unsupervised_pretraining/test_kipoiseq/1000_metadata.tsv'
metadata_df = pd.read_csv(metadata, sep='\t')

metadata_df = metadata_df[['Sample name','Superpopulation name','Population name','Superpopulation code','Population code']]

metadata_df = metadata_df.rename(columns={"Sample name": "SAMPLE_NAME",
                                          "Superpopulation name": "superpopulation_name",
                                          "Population name": "population_name",
                                          "Superpopulation code": "superpopulation_code",
                                          "Population code": "population_code",
                                          })

populations = metadata_df.population_code.unique().tolist()

test_population = ['YRI']
validation_population = ['CLM']
test_validation_population = ['YRI','CLM']

train_population = list(set(populations).difference(test_validation_population))



black_list_bed = environment['minerva']['blacklist_regions_hg38']



test_chromosomes = ["chr9"]
validation_chromosomes = ["chr5"]
train_chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr6", "chr7", "chr8", "chr10", "chr11", "chr12","chr13",
                     "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22"]



all_chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
                   "chr11", "chr12","chr13",
                   "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22"]


sample_df = pd.read_csv(sample_id_unrelated, sep='\t', comment='#',
                        header=None,
                        names=cols)

sample_df = sample_df[['SAMPLE_NAME']]


metadata_df = pd.read_csv(metadata, sep='\t')

metadata_df = metadata_df[['Sample name','Superpopulation name','Population name','Superpopulation code','Population code']]

metadata_df = metadata_df.rename(columns={"Sample name": "SAMPLE_NAME",
                                          "Superpopulation name": "superpopulation_name",
                                          "Population name": "population_name",
                                          "Superpopulation code": "superpopulation_code",
                                          "Population code": "population_code",
                                          })

populations = metadata_df.population_code.unique().tolist()

test_population = ['YRI']
validation_population = ['CLM']
test_validation_population = ['YRI','CLM']

train_population = list(set(populations).difference(test_validation_population))



def scrub_string(s):
    s = s.strip()
    s = re.sub('[^0-9a-zA-Z]+', '', s)
    s = re.sub('[0-9]+', '', s)
    s = re.sub('[^ACGTN]', 'N', s)
    return s

class SampleSeqExtractor(VariantSeqExtractor):
    def __init__(self, fasta_file, vcf_file):
        """Sequence extractor which can extract an alternate sequence for a
        given interval and the variants corresponding to a given
        sample and phase.
        Args:
          fasta_file: Path to the fasta file containing the reference
            sequence (can be gzipped)
          vcf_file: Path to the VCF file containing phased genotype information
        """
        self.vcf = VCF(vcf_file)
        self._sample_indices = dict(zip(self.vcf.samples,
                                        range(len(self.vcf.samples))))
        super().__init__(fasta_file)
    def extract(self, interval, sample, phase, anchor,
                fixed_len=True, **kwargs):
        """Extracts an alternate sequence for a given interval and the
        variants corresponding to a given sample.
        Args:
          interval: `kipoiseq.dataclasses.Interval`, Region of
            interest from which to query the sequence. 0-based.
          sample: `str`, Sample from the VCF file for which variants should be
            extracted.
          phase: `0` or `1`, Phase for which sequence should be extracted
          anchor: `int`, Absolution position w.r.t. the interval
            start. (0-based).  E.g. for an interval of `chr1:10-20`
            the anchor of 10 denotes the point chr1:10 in the 0-based
            coordinate system.
          fixed_len: `bool`, If True, the return sequence will have the
            same length as the `interval` (e.g. `interval.end -
            interval.start`)
          kwargs: Additional keyword arguments to pass to
            `SampleSeqExtractor.extract`
        Returns:
          A single sequence (`str`) with all the variants applied.
        """
        variants = []
        if sample is not None:
            if sample not in self.vcf.samples:
                raise ValueError(f'Sample {sample} not present in VCF file')
            if phase not in (0, 1):
                raise ValueError('phase argument must be in (0, 1) if sample is not None')
            # Interval is  0-based, cyvcf2 positions are 1-based: need to add 1
            variants = self._get_sample_variants(
                self.vcf(f'{interval.chrom}:'
                    + f'{interval.start + 1}-{interval.end + 1}'
                ),
                sample,
                phase
            )
        return super(SampleSeqExtractor, self).extract(
            interval, variants, anchor, fixed_len, **kwargs)
    def _get_sample_variants(self, variants, sample, phase):
        """Given a list of `cyvcf2.Variant`, returns all those present for a
        given sample and phase and converts them to
        `kipoiseq.dataclasses.Variant`
        Args:
          variants: List of `cyvcf2.Variant`, Variants of interest
          sample: `str`, Sample for which to filter genotypes
          phase: `0` or `1`, Phase for which to filter genotypes
        Returns:
          List of `kipoiseq.dataclasses.Variant`
        """
        sample_index = self._sample_indices[sample]
        return [
            Variant.from_cyvcf(v) for v in variants
            if v.genotypes[sample_index][phase]
        ]


## Change
chromosomes = train_chromosomes
context_length = 2 ** 19
scaling_factor = 1
#shift_range = 75000
#train_list = []
#val_list = []
tile_list = []
context_fourths = context_length//4
shift_range = context_fourths // 2
print(f'Context Length: {context_length}')
for shift in [0, context_fourths * 1 , context_fourths * 2 , context_fourths * 3]:
    print(shift)
    blacklist_bed = pr.read_bed(black_list_bed)
    chromsizes = pr.data.chromsizes()
    tile = pr.gf.tile_genome(chromsizes, context_length)
    tile.Start = tile.Start + shift
    tile.End = tile.End + shift + shift_range
    tile_black_list = tile.overlap(blacklist_bed, how='first')
    #tile_df = tile.dfs
    #tile_df = pd.concat(tile_df.values(), ignore_index=True)
    tile = tile[tile.Chromosome.isin(chromosomes)].dfs
    #val = tile[tile.Chromosome.isin(validation_chromosomes)].dfs
    tile = pd.concat(tile.values(), ignore_index=True)
    #val = pd.concat(val.values(), ignore_index=True)
    tile_black_list_df = pd.concat(tile_black_list.dfs.values(), ignore_index=True)[['Chromosome', 'Start', 'End']]
    tile = anti_join(tile, tile_black_list_df, on=['Chromosome', 'Start', 'End'])
    tile_list.append(tile)
    #val = anti_join(val, tile_black_list_df, on=['Chromosome', 'Start', 'End'])
    #train_list.append(train)
    #val_list.append(val)

tile_all = pd.concat(tile_list, axis=0)
len_tile_all = len(tile_all)

tile_all['length'] = tile_all['End'] - tile_all['Start']

tile_all = tile_all[tile_all['length'] == (context_length + shift_range)]
tile_all.reset_index(inplace=True, drop=True)


#train = pd.concat(train_list, axis=0)
#val = pd.concat(val_list, axis=0)
#train.reset_index(inplace=True, drop=True)
#val.reset_index(inplace=True, drop=True)
#len_train = len(train)
#len_val = len(val)

print("All Length: {}".format(len_tile_all))


###Change
population = train_population
sample_df = sample_df.merge(metadata_df, on='SAMPLE_NAME')
sample_df = sample_df[sample_df['population_code'].isin(population)]
###Change
split = 'train'
sample_df
chr_string = subset['Chromosome']
vals = []
count = 0
len_df = len(sample_df)
vcf = f"/sc/arion/projects/ad-omics/clakhani/1KG_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_{chr_string}.filtered.shapeit2-duohmm-phased.vcf.gz"
df_1KG = SampleSeqExtractor(fasta_file, vcf)

interval = Interval(subset['Chromosome'], subset['Start'], subset['End'])
center = interval.center() - interval.start
print(interval)
print(f'Count: {count} out of {len_df}')
for sample_index, sample_row in sample_df.iterrows():
    start = time.time()
    print(sample_row['SAMPLE_NAME'])
    seq1 = df_1KG.extract(interval, sample_row['SAMPLE_NAME'], 0, center, fixed_len=False)
    seq2 = df_1KG.extract(interval, sample_row['SAMPLE_NAME'], 1, center, fixed_len=False)
    seq1 = seq1.upper()
    seq2 = seq2.upper()
    seq1 = scrub_string(seq1)
    seq2 = scrub_string(seq2)
    dict_seq = {'Chromosome': subset['Chromosome'],
                 'Start': subset['Start'],
                 'End': subset['End'],
                 'SAMPLE_NAME': sample_row['SAMPLE_NAME'],
                 'super_population': sample_row['superpopulation_name'],
                 'population': sample_row['population_name'],
                 'superpopulation_code': sample_row['superpopulation_code'],
                 'population_code': sample_row['population_code'],
                 'seq_0': seq1,
                 'len_seq_0': len(seq1),
                 'pct_N_seq_0': seq1.count('N') / len(seq1),
                 'seq_1': seq2,
                 'len_seq_1': len(seq2),
                 'pct_N_seq_1': seq2.count('N') / len(seq2)
            }
    vals.append(dict_seq)
    end = time.time()
    elapsed = end-start
    print(f"count: {count} out of {len_df} in {elapsed} seconds")
    count += 1

output_df = pd.DataFrame(vals)

output_df = output_df[output_df.pct_N_seq_0 <= .05]
output_df = output_df[output_df.pct_N_seq_1 <= .05]

#output_df = output_df[output_df.len_seq_0 >= context_length]
#output_df = output_df[output_df.len_seq_1 >= context_length]

output_df = output_df[["Chromosome", "Start","End","SAMPLE_NAME","superpopulation_code","population_code","seq_0","seq_1"]]

out_dir = f'/sc/arion/projects/ad-omics/clakhani/1KG_phased/training_data/{split}/'

if not os.path.exists(out_dir):
    os.makedirs(out_dir)

output_file = f"{out_dir}/{split}_{chr_string}_chunk_{idx_num:05d}.tsv.gz"

print(output_df.shape)

print(output_df.head())

if len(output_df) != 0:
    output_df.to_csv(output_file, sep='\t', index=False, compression = 'gzip')
