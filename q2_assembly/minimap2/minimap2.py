# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
import shutil
import subprocess
import tempfile
from typing import Union

import pandas as pd
from q2_types.feature_data import DNAFASTAFormat
from q2_types.per_sample_sequences import (
    SingleLanePerSamplePairedEndFastqDirFmt,
    SingleLanePerSampleSingleEndFastqDirFmt,
)
from q2_types_genomics.per_sample_data._format import BAMDirFmt

from q2_assembly._utils import run_command, run_commands_with_pipe


def _process_sample(sample, fwd, rev, reference, common_args, out):
    """ """
    with tempfile.TemporaryDirectory() as tmp:
        results_fp = os.path.join(tmp, f"{sample}.bam")
        cmd1 = ["minimap2", *common_args, "-a", str(reference.path), fwd]
        if rev:
            cmd1.append(rev)

        cmd2 = ["samtools", "view", "-O", "BAM", "-o", results_fp]

        try:
            run_commands_with_pipe(cmd1, cmd2, verbose=True)
        except subprocess.CalledProcessError as e:
            raise Exception(
                "An error was encountered while running Minimap2, "
                f"(return code {e.returncode}), please inspect "
                "stdout and stderr to learn more."
            )

        shutil.move(
            results_fp,
            os.path.join(str(out), f"{sample}.bam"),
        )


def _align_minimap2(seqs, reference, common_args) -> BAMDirFmt:
    paired = isinstance(seqs, SingleLanePerSamplePairedEndFastqDirFmt)
    manifest = seqs.manifest.view(pd.DataFrame)
    result = BAMDirFmt()

    for samp in list(manifest.index):
        fwd = manifest.loc[samp, "forward"]
        rev = manifest.loc[samp, "reverse"] if paired else None

        _process_sample(samp, fwd, rev, reference, common_args, result)
    return result


def align_minimap2(
    seqs: Union[
        SingleLanePerSamplePairedEndFastqDirFmt, SingleLanePerSampleSingleEndFastqDirFmt
    ],
    reference: DNAFASTAFormat,
    coalign: bool = False,
    preset: str = "sr",
    t: int = 3,
) -> BAMDirFmt:
    common_args = ["-t", str(t), "-x", preset]
    return _align_minimap2(seqs=seqs, reference=reference, common_args=common_args)
