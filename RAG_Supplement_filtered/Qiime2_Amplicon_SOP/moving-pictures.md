---
source_url: https://amplicon-docs.qiime2.org/en/latest/tutorials/moving-pictures.html
download_date: 2025-11-20
---

Downloads

# Moving Pictures tutorial 🎥

## Background(#oR84T9baq5 "Link to this Section")

In this tutorial you’ll use the amplicon distribution of [QIIME 2](https://qiime2.org) to perform an analysis of human microbiome samples from two individuals at four body sites and five timepoints, the first of which immediately followed antibiotic usage.
The samples used here were originally published in Caporaso et al. (2011).
The data used here is a small subset of that generated for the paper, specifically selected so that this tutorial can be run quickly on a personal computer.

The data used in this tutorial was generated using the Earth Microbiome Project protocol.
Specifically, the hypervariable region 4 (V4) of the 16S rRNA gene was amplified using the F515-R806 primers - a broad-coverage primer pair for Bacteria that also picks up some Archaea.
Single-end sequencing was performed on an Illumina HiSeq.

The primers used here are still popular, though more modern ones have been developed since the paper was published.
Typically an Illumina amplicon sequencing run will generate paired end reads (as opposed to the single-end reads used here), and a more modern instrument would be used.
All of that said, this is a great data set for learning QIIME 2; you shouldn’t aim to exactly model the primers and protocol used for data generation though.

## Sample metadata(#Uil4JXXp7a "Link to this Section")

Before starting the analysis, explore the sample metadata to familiarize yourself with the samples used in this study.
The following command will download the sample metadata as tab-separated text and save it in the file `sample-metadata.tsv`.
This `sample-metadata.tsv` file is used throughout the rest of the tutorial.

[Command Line]

[Python API]

[Galaxy]

[R API]

[View Source]

```
wget -O 'sample-metadata.tsv' \
  'https://moving-pictures-tutorial.readthedocs.io/en/latest/data/moving-pictures/sample-metadata.tsv'
```

```
from qiime2 import Metadata
from urllib import request

url = 'https://moving-pictures-tutorial.readthedocs.io/en/latest/data/moving-pictures/sample-metadata.tsv'
fn = 'sample-metadata.tsv'
request.urlretrieve(url, fn)
sample_metadata_md = Metadata.load(fn)
```

Using the `Upload Data` tool:
:   1. On the first tab (**Regular**), press the `Paste/Fetch` data button at the bottom.
       1. Set *"Name"* (first text-field) to: `sample-metadata.tsv`
       2. In the larger text-area, copy-and-paste: <https://moving-pictures-tutorial.readthedocs.io/en/latest/data/moving-pictures/sample-metadata.tsv>
       3. (*"Type"*, *"Genome"*, and *"Settings"* can be ignored)
    2. Press the `Start` button at the bottom.

```
library(reticulate)

Metadata <- import("qiime2")$Metadata
request <- import("urllib")$request

url <- 'https://moving-pictures-tutorial.readthedocs.io/en/latest/data/moving-pictures/sample-metadata.tsv'
fn <- 'sample-metadata.tsv'
request$urlretrieve(url, fn)
sample_metadata_md <- Metadata$load(fn)
```

```
sample_metadata = use.init_metadata_from_url(
   'sample-metadata',
   'https://data.qiime2.org/2025.4/tutorials/moving-pictures/sample_metadata.tsv')
```

* `sample-metadata.tsv` | [download](https://moving-pictures-tutorial.readthedocs.io/en/latest/data/moving-pictures/sample-metadata.tsv)

QIIME 2’s metadata plugin provides a Visualizer called `tabulate` that generates a convenient view of a sample metadata file.
To learn more about `metadata tabulate`, see [here](https://library.qiime2.org/myst/amplicon/_/references/plugins/metadata#q2-action-metadata-tabulate).
Let’s run this, and then we’ll look at the result.
Here’s the first QIIME 2 command that you should run in this tutorial:

[Command Line]

[Python API]

[Galaxy]

[R API]

[View Source]

```
qiime metadata tabulate \
  --m-input-file sample-metadata.tsv \
  --o-visualization sample-metadata-viz.qzv
```

```
import qiime2.plugins.metadata.actions as metadata_actions

sample_metadata_viz_viz, = metadata_actions.tabulate(
    input=sample_metadata_md,
)
```

Using the `qiime2 metadata tabulate` tool:
:   1. For *"input"*:
       * Perform the following steps.
         1. Leave as `Metadata from TSV`
         2. Set *"Metadata Source"* to `sample-metadata.tsv`
    2. Press the `Execute` button.

Once completed, for the new entry in your history, use the `Edit` button to set the name as follows:
:   (Renaming is optional, but it will make any subsequent steps easier to complete.)

    | History Name | "Name" to set (be sure to press [Save]) |
    | --- | --- |
    | `#: qiime2 metadata tabulate [...] : visualization.qzv` | `sample-metadata-viz.qzv` |

```
metadata_actions <- import("qiime2.plugins.metadata.actions")

action_results <- metadata_actions$tabulate(
    input=sample_metadata_md,
)
sample_metadata_viz_viz <- action_results$visualization
```

```
sample_metadata_viz, = use.action(
  use.UsageAction(plugin_id='metadata',
                  action_id='tabulate'),
  use.UsageInputs(input=sample_metadata),
  use.UsageOutputNames(visualization='sample_metadata_viz')
)
```

* `sample-metadata-viz.qzv` | [download](https://moving-pictures-tutorial.readthedocs.io/en/latest/data/moving-pictures/sample-metadata-viz.qzv) | [view](https://view.qiime2.org?src=https%3A%2F%2Fmoving-pictures-tutorial.readthedocs.io%2Fen%2Flatest%2Fdata%2Fmoving-pictures%2Fsample-metadata-viz.qzv)

This will generate a QIIME 2 Visualization.
Visualizations can be viewed by loading them with [QIIME 2 View](https://view.qiime2.org).
Navigate to QIIME 2 View, and drag and drop the visualization that was created to view it.
(You can learn more about viewing Visualizations, including alternatives to QIIME 2 View if you can’t use that for any reason, [here](https://use.qiime2.org/en/latest/how-to-guides/view-visualizations.html).)

## Obtaining and importing data(#mojLSdAFl1 "Link to this Section")

Download the sequence reads that we’ll use in this analysis.
In this tutorial we’ll work with a small subset of the complete sequence data so that the commands will run quickly.

[Command Line]

[Python API]

[Galaxy]

[R API]

[View Source]

```
wget -O 'emp-single-end-sequences.zip' \
  'https://moving-pictures-tutorial.readthedocs.io/en/latest/data/moving-pictures/emp-single-end-sequences.zip'

unzip -d emp-single-end-sequences emp-single-end-sequences.zip
```

```
import zipfile

url = 'https://moving-pictures-tutorial.readthedocs.io/en/latest/data/moving-pictures/emp-single-end-sequences.zip'
fn = 'emp-single-end-sequences.zip'
request.urlretrieve(url, fn)
with zipfile.ZipFile(fn) as zf:
    zf.extractall('emp-single-end-sequences')
```

Using the `Upload Data` tool:
:   * Steps to setup `emp-single-end-sequences`:
      1. On the fourth tab (**Rule-based**):
         1. Set *"Upload data as"* to `Datasets`
         2. Set *"Load tabular data from"* to `Pasted Table`
         3. Paste the following contents into the large text area:

            ```
            emp-single-end-sequences:sequences	https://moving-pictures-tutorial.readthedocs.io/en/latest/data/moving-pictures/emp-single-end-sequences/sequences.fastq.gz
            emp-single-end-sequences:barcodes	https://moving-pictures-tutorial.readthedocs.io/en/latest/data/moving-pictures/emp-single-end-sequences/barcodes.fastq.gz
            ```
         4. Press the `build` button at the bottom.
      2. In the resulting UI, do the following:
         1. Add a rule by pressing the `+ Rules` button and choosing `Add / Modify Column Definitions`.
         2. In the sidebar:
            1. Press `+Add Definition` and select `Name`. (This will choose column "A" by default. You should see a new annotation on the column header.)
            2. Press `+Add Definition` and select `URL`.
            3. Change the dropdown above the button to be `B`. (You should see the table headers list `A (Name)` and `B (URL)`.)
            4. Press the `Apply` button.
      3. Press the `Upload` button at the bottom right.

```
zipfile <- import("zipfile")

url <- 'https://moving-pictures-tutorial.readthedocs.io/en/latest/data/moving-pictures/emp-single-end-sequences.zip'
fn <- 'emp-single-end-sequences.zip'
request$urlretrieve(url, fn)
zf <- zipfile$ZipFile(fn)
zf$extractall('emp-single-end-sequences')
zf$close()
```

```
def emp_factory():
    import os
    import tempfile
    from urllib import request

    from q2_types.multiplexed_sequences import EMPSingleEndDirFmt
    from q2_types.per_sample_sequences import FastqGzFormat

    base_url = 'https://data.qiime2.org/2025.4/tutorials/moving-pictures/'
    bc_url = base_url + 'emp-single-end-sequences/barcodes.fastq.gz'
    seqs_url = base_url + 'emp-single-end-sequences/sequences.fastq.gz'

    fmt = EMPSingleEndDirFmt(mode='w')

    with tempfile.TemporaryDirectory() as tmpdir:
        bc_fp = os.path.join(tmpdir, 'barcodes.fastq.gz')
        bc_fn, _ = request.urlretrieve(bc_url, bc_fp)

        seqs_fp = os.path.join(tmpdir, 'sequences.fastq.gz')
        seqs_fn, _ = request.urlretrieve(seqs_url, seqs_fp)

        fmt.barcodes.write_data(bc_fn, FastqGzFormat)
        fmt.sequences.write_data(seqs_fn, FastqGzFormat)

    fmt.validate()
    return fmt

raw_seqs = use.init_format('emp-single-end-sequences', emp_factory)
```

* `emp-single-end-sequences.zip` | [download](https://moving-pictures-tutorial.readthedocs.io/en/latest/data/moving-pictures/emp-single-end-sequences.zip)

All data that is used as input to QIIME 2 is in form of QIIME 2 artifacts, which contain information about the type of data and the source of the data.
So, the first thing we need to do is import these sequence data files into a QIIME 2 artifact.

The semantic type of this QIIME 2 artifact is `EMPSingleEndSequences`.
`EMPSingleEndSequences` QIIME 2 artifacts contain sequences that are multiplexed, meaning that the sequences have not yet been assigned to samples (where the `barcodes.fastq.gz` contains the barcode read associated with each sequence in `sequences.fastq.gz`.)
To learn about how to import sequence data in other formats, see *[How to import data for use with QIIME 2](https://amplicon-docs.readthedocs.io/en/latest/how-to-guides/how-to-import.html)*.

[Command Line]

[Python API]

[Galaxy]

[R API]

[View Source]

```
qiime tools import \
  --type 'EMPSingleEndSequences' \
  --input-path emp-single-end-sequences \
  --output-path emp-single-end-sequences.qza
```

```
from qiime2 import Artifact

emp_single_end_sequences = Artifact.import_data(
    'EMPSingleEndSequences',
    'emp-single-end-sequences',
)
```

Using the `qiime2 tools import` tool:
:   1. Set *"Type of data to import"* to `EMPSingleEndSequences`
    2. Set *"QIIME 2 file format to import from"* to `EMP Single End Directory Format`
    3. For `import_sequences`, do the following:
       1. Leave *"name"* as `sequences.fastq.gz`
       2. Set *"data"* to `#: emp-single-end-sequences:sequences`
    4. For `import_barcodes`, do the following:
       1. Leave *"name"* as `barcodes.fastq.gz`
       2. Set *"data"* to `#: emp-single-end-sequences:barcodes`
    5. Press the `Execute` button.

Once completed, for the new entry in your history, use the `Edit` button to set the name as follows:
:   (Renaming is optional, but it will make any subsequent steps easier to complete.)

    | History Name | "Name" to set (be sure to press [Save]) |
    | --- | --- |
    | `#: qiime2 tools import [...]` | `emp-single-end-sequences.qza` |

```
Artifact <- import("qiime2")$Artifact

emp_single_end_sequences <- Artifact$import_data(
    'EMPSingleEndSequences',
    'emp-single-end-sequences',
)
```

```
emp_single_end_sequences = use.import_from_format(
    'emp_single_end_sequences',
    'EMPSingleEndSequences',
    raw_seqs)
```

* `emp-single-end-sequences.qza` | [download](https://moving-pictures-tutorial.readthedocs.io/en/latest/data/moving-pictures/emp-single-end-sequences.qza) | [view](https://view.qiime2.org?src=https%3A%2F%2Fmoving-pictures-tutorial.readthedocs.io%2Fen%2Flatest%2Fdata%2Fmoving-pictures%2Femp-single-end-sequences.qza)

## Demultiplexing sequences(#kfkxeS2wMe "Link to this Section")

To demultiplex sequences we need to know which barcode sequence is associated with each sample.
This information is contained in the [sample metadata](https://data.qiime2.org/2025.4/tutorials/moving-pictures/sample_metadata) file.
You can run the following commands to demultiplex the sequences (the `demux emp-single` command refers to the fact that these sequences are barcoded according to the [Earth Microbiome Project](http://earthmicrobiome.org) protocol, and are single-end reads).
The `demux.qza` QIIME 2 artifact will contain the demultiplexed sequences.
The second output (`demux-details.qza`) presents Golay error correction details, and will not be explored in this tutorial (you can visualize these data using `qiime metadata tabulate`).

[Command Line]

[Python API]

[Galaxy]

[R API]

[View Source]

```
qiime demux emp-single \
  --i-seqs emp-single-end-sequences.qza \
  --m-barcodes-file sample-metadata.tsv \
  --m-barcodes-column barcode-sequence \
  --o-per-sample-sequences demux.qza \
  --o-error-correction-details demux-details.qza
```

```
import qiime2.plugins.demux.actions as demux_actions

barcode_sequence_mdc = sample_metadata_md.get_column('barcode-sequence')
demux, demux_details = demux_actions.emp_single(
    seqs=emp_single_end_sequences,
    barcodes=barcode_sequence_mdc,
)
```

Using the `qiime2 demux emp-single` tool:
:   1. Set *"seqs"* to `#: emp-single-end-sequences.qza`
    2. For *"barcodes"*:
       1. Leave as `Metadata from TSV`
       2. Set *"Metadata Source"* to `sample-metadata.tsv`
       3. Set *"Column Name"* to `barcode-sequence`
    3. Press the `Execute` button.

Once completed, for each new entry in your history, use the `Edit` button to set the name as follows:
:   (Renaming is optional, but it will make any subsequent steps easier to complete.)

    | History Name | "Name" to set (be sure to press [Save]) |
    | --- | --- |
    | `#: qiime2 demux emp-single [...] : per_sample_sequences.qza` | `demux.qza` |
    | `#: qiime2 demux emp-single [...] : error_correction_details.qza` | `demux-details.qza` |

```
demux_actions <- import("qiime2.plugins.demux.actions")

barcode_sequence_mdc <- sample_metadata_md$get_column('barcode-sequence')
action_results <- demux_actions$emp_single(
    seqs=emp_single_end_sequences,
    barcodes=barcode_sequence_mdc,
)
demux <- action_results$per_sample_sequences
demux_details <- action_results$error_correction_details
```

```
barcode_sequence = use.get_metadata_column(
    'barcode_sequence', 'barcode-sequence', sample_metadata)

demux, demux_details = use.action(
    use.UsageAction(plugin_id='demux',
                    action_id='emp_single'),
    use.UsageInputs(seqs=emp_single_end_sequences,
                    barcodes=barcode_sequence),
    use.UsageOutputNames(per_sample_sequences='demux',
                         error_correction_details='demux_details'))
```

* `demux.qza` | [download](https://moving-pictures-tutorial.readthedocs.io/en/latest/data/moving-pictures/demux.qza) | [view](https://view.qiime2.org?src=https%3A%2F%2Fmoving-pictures-tutorial.readthedocs.io%2Fen%2Flatest%2Fdata%2Fmoving-pictures%2Fdemux.qza)
* `demux-details.qza` | [download](https://moving-pictures-tutorial.readthedocs.io/en/latest/data/moving-pictures/demux-details.qza) | [view](https://view.qiime2.org?src=https%3A%2F%2Fmoving-pictures-tutorial.readthedocs.io%2Fen%2Flatest%2Fdata%2Fmoving-pictures%2Fdemux-details.qza)

After demultiplexing, it’s useful to generate a summary of the demultiplexing results.
This allows you to determine how many sequences were obtained per sample, and also to get a summary of the distribution of sequence qualities at each position in your sequence data.

[Command Line]

[Python API]

[Galaxy]

[R API]

[View Source]

```
qiime demux summarize \
  --i-data demux.qza \
  --o-visualization demux.qzv
```

```
demux_viz, = demux_actions.summarize(
    data=demux,
)
```

Using the `qiime2 demux summarize` tool:
:   1. Set *"data"* to `#: demux.qza`
    2. Press the `Execute` button.

Once completed, for the new entry in your history, use the `Edit` button to set the name as follows:
:   (Renaming is optional, but it will make any subsequent steps easier to complete.)

    | History Name | "Name" to set (be sure to press [Save]) |
    | --- | --- |
    | `#: qiime2 demux summarize [...] : visualization.qzv` | `demux.qzv` |

```
action_results <- demux_actions$summarize(
    data=demux,
)
demux_viz <- action_results$visualization
```

```
use.action(
    use.UsageAction(plugin_id='demux',
                    action_id='summarize'),
    use.UsageInputs(data=demux),
    use.UsageOutputNames(visualization='demux'))
```

* `demux.qzv` | [download](https://moving-pictures-tutorial.readthedocs.io/en/latest/data/moving-pictures/demux.qzv) | [view](https://view.qiime2.org?src=https%3A%2F%2Fmoving-pictures-tutorial.readthedocs.io%2Fen%2Flatest%2Fdata%2Fmoving-pictures%2Fdemux.qzv)

## Sequence quality control and feature table construction(#yyK0f2Ksls "Link to this Section")

QIIME 2 plugins are available for several quality control methods, including DADA2, Deblur, and basic quality-score-based filtering.
In this tutorial we present this step using DADA2 and Deblur.
These steps are interchangeable, so you can use whichever of these you prefer.
A comparison of these is presented in Nearing *et al.* (2018).
The result of both of these methods will be a `FeatureTable[Frequency]` QIIME 2 artifact, which contains counts (frequencies) of each unique sequence in each sample in the dataset, and a `FeatureData[Sequence]` QIIME 2 artifact, which maps feature identifiers in the `FeatureTable` to the sequences they represent.

### Option 1: DADA2(#E5jTBD9Ne9 "Link to this Section")

DADA2 is a pipeline for detecting and correcting (where possible) Illumina amplicon sequence data.
As implemented in the `q2-dada2` plugin, this quality control process will additionally filter any phiX reads (commonly present in marker gene Illumina sequence data) that are identified in the sequencing data, and will filter chimeric sequences.

The `dada2 denoise-single` method requires two parameters that are used in quality filtering: `--p-trim-left m`, which trims off the first `m` bases of each sequence, and `--p-trunc-len n` which truncates each sequence at position `n`.
This allows the user to remove low quality regions of the sequences.
To determine what values to pass for these two parameters, you should review the *Interactive Quality Plot* tab in the `demux.qzv` file that was generated by `qiime demux summarize` above.

In the `demux.qzv` quality plots, we see that the quality of the initial bases seems to be high, so we won’t trim any bases from the beginning of the sequences.
The quality seems to drop off around position 120, so we’ll truncate our sequences at 120 bases.
This next command may take up to 10 minutes to run, and is the slowest step in this tutorial.

[Command Line]

[View Source]

```
qiime dada2 denoise-single \
  --i-demultiplexed-seqs demux.qza \
  --p-trim-left 0 \
  --p-trunc-len 120 \
  --o-representative-sequences rep-seqs.qza \
  --o-table table.qza \
  --o-denoising-stats stats.qza
```

```
rep_seqs_dada2, table_dada2, stats_dada2 = use.action(
    use.UsageAction(plugin_id='dada2',
                    action_id='denoise_single'),
    use.UsageInputs(demultiplexed_seqs=demux,
                    trim_left=0,
                    trunc_len=120),
    use.UsageOutputNames(representative_sequences='rep_seqs',
                         table='table',
                         denoising_stats='stats'))
```

```
Traceback (most recent call last):
  File "/home/docs/checkouts/readthedocs.org/user_builds/moving-pictures-tutorial/checkouts/latest/.env/lib/python3.10/site-packages/q2doc/transforms/transform_usage.py", line 92, in run
    exec(source, interface['driver'].scope)
  File "<string>", line 1, in <module>
  File "/home/docs/checkouts/readthedocs.org/user_builds/moving-pictures-tutorial/checkouts/latest/.env/lib/python3.10/site-packages/qiime2/plugins.py", line 344, in action
    self._template_action(action, inputs, variables)
  File "/home/docs/checkouts/readthedocs.org/user_builds/moving-pictures-tutorial/checkouts/latest/.env/lib/python3.10/site-packages/qiime2/plugins.py", line 363, in _template_action
    output_vars = self._template_outputs(action, variables)
  File "/home/docs/checkouts/readthedocs.org/user_builds/moving-pictures-tutorial/checkouts/latest/.env/lib/python3.10/site-packages/qiime2/plugins.py", line 390, in _template_outputs
    variable = getattr(variables, output)
AttributeError: 'UsageOutputs' object has no attribute 'base_transition_stats'
```

[Command Line]

[View Source]

```
qiime metadata tabulate \
  --m-input-file stats.qza \
  --o-visualization stats.qzv
```

```
stats_as_md = use.view_as_metadata('stats_dada2_md', stats_dada2)

use.action(
    use.UsageAction(plugin_id='metadata',
                    action_id='tabulate'),
    use.UsageInputs(input=stats_as_md),
    use.UsageOutputNames(visualization='stats'))
```

```
Traceback (most recent call last):
  File "/home/docs/checkouts/readthedocs.org/user_builds/moving-pictures-tutorial/checkouts/latest/.env/lib/python3.10/site-packages/q2doc/transforms/transform_usage.py", line 92, in run
    exec(source, interface['driver'].scope)
  File "<string>", line 1, in <module>
NameError: name 'stats_dada2' is not defined
```

### Option 2: Deblur(#hjmAjuORoh "Link to this Section")

Deblur uses sequence error profiles to associate erroneous sequence reads with the true biological sequence from which they are derived, resulting in high quality sequence variant data.
This is applied in two steps.
First, an initial quality filtering process based on quality scores is applied.
This method is an implementation of the quality filtering approach described by Bokulich et al. (2013).

[Command Line]

[Python API]

[Galaxy]

[R API]

[View Source]

```
qiime quality-filter q-score \
  --i-demux demux.qza \
  --o-filtered-sequences demux-filtered.qza \
  --o-filter-stats demux-filter-stats.qza
```

```
import qiime2.plugins.dada2.actions as dada2_actions
import qiime2.plugins.quality_filter.actions as quality_filter_actions

demux_filtered, demux_filter_stats = quality_filter_actions.q_score(
    demux=demux,
)
```

Using the `qiime2 quality-filter q-score` tool:
:   1. Set *"demux"* to `#: demux.qza`
    2. Press the `Execute` button.

Once completed, for each new entry in your history, use the `Edit` button to set the name as follows:
:   (Renaming is optional, but it will make any subsequent steps easier to complete.)

    | History Name | "Name" to set (be sure to press [Save]) |
    | --- | --- |
    | `#: qiime2 quality-filter q-score [...] : filtered_sequences.qza` | `demux-filtered.qza` |
    | `#: qiime2 quality-filter q-score [...] : filter_stats.qza` | `demux-filter-stats.qza` |

```
quality_filter_actions <- import("qiime2.plugins.quality_filter.actions")

action_results <- quality_filter_actions$q_score(
    demux=demux,
)
demux_filtered <- action_results$filtered_sequences
demux_filter_stats <- action_results$filter_stats
```

```
filtered_seqs, filter_stats = use.action(
    use.UsageAction(plugin_id='quality_filter',
                    action_id='q_score'),
    use.UsageInputs(demux=demux),
    use.UsageOutputNames(filtered_sequences='demux_filtered',
                         filter_stats='demux_filter_stats'))
```

* `demux-filtered.qza` | [download](https://moving-pictures-tutorial.readthedocs.io/en/latest/data/moving-pictures/demux-filtered.qza) | [view](https://view.qiime2.org?src=https%3A%2F%2Fmoving-pictures-tutorial.readthedocs.io%2Fen%2Flatest%2Fdata%2Fmoving-pictures%2Fdemux-filtered.qza)
* `demux-filter-stats.qza` | [download](https://moving-pictures-tutorial.readthedocs.io/en/latest/data/moving-pictures/demux-filter-stats.qza) | [view](https://view.qiime2.org?src=https%3A%2F%2Fmoving-pictures-tutorial.readthedocs.io%2Fen%2Flatest%2Fdata%2Fmoving-pictures%2Fdemux-filter-stats.qza)

Next, the Deblur workflow is applied using the `qiime deblur denoise-16S` method.
This method requires one parameter that is used in quality filtering, `--p-trim-length n` which truncates the sequences at position `n`.
In general, the Deblur developers recommend setting this value to a length where the median quality score begins to drop too low.
On these data, the quality plots (prior to quality filtering) suggest a reasonable choice is in the 115 to 130 sequence position range.
This is a subjective assessment.
One situation where you might deviate from that recommendation is when performing a meta-analysis across multiple sequencing runs.
In this type of meta-analysis, it is critical that the read lengths be the same for all of the sequencing runs being compared to avoid introducing a study-specific bias.
Since we already using a trim length of 120 for `qiime dada2 denoise-single`, and since 120 is reasonable given the quality plots, we’ll pass `--p-trim-length 120`.
This next command may take up to 10 minutes to run.

[Command Line]

[Python API]

[Galaxy]

[R API]

[View Source]

```
qiime deblur denoise-16S \
  --i-demultiplexed-seqs demux-filtered.qza \
  --p-trim-length 120 \
  --p-sample-stats \
  --o-representative-sequences rep-seqs-deblur.qza \
  --o-table table-deblur.qza \
  --o-stats deblur-stats.qza
```

```
import qiime2.plugins.deblur.actions as deblur_actions

table_deblur, rep_seqs_deblur, deblur_stats = deblur_actions.denoise_16S(
    demultiplexed_seqs=demux_filtered,
    trim_length=120,
    sample_stats=True,
)
```

Using the `qiime2 deblur denoise-16S` tool:
:   1. Set *"demultiplexed\_seqs"* to `#: demux-filtered.qza`
    2. Set *"trim\_length"* to `120`
    3. Expand the `additional options` section
       * Set *"sample\_stats"* to `Yes`
    4. Press the `Execute` button.

Once completed, for each new entry in your history, use the `Edit` button to set the name as follows:
:   (Renaming is optional, but it will make any subsequent steps easier to complete.)

    | History Name | "Name" to set (be sure to press [Save]) |
    | --- | --- |
    | `#: qiime2 deblur denoise-16S [...] : table.qza` | `table-deblur.qza` |
    | `#: qiime2 deblur denoise-16S [...] : representative_sequences.qza` | `rep-seqs-deblur.qza` |
    | `#: qiime2 deblur denoise-16S [...] : stats.qza` | `deblur-stats.qza` |

```
deblur_actions <- import("qiime2.plugins.deblur.actions")

action_results <- deblur_actions$denoise_16S(
    demultiplexed_seqs=demux_filtered,
    trim_length=120L,
    sample_stats=TRUE,
)
rep_seqs_deblur <- action_results$representative_sequences
table_deblur <- action_results$table
deblur_stats <- action_results$stats
```

```
rep_seqs_deblur, table_deblur, stats_deblur = use.action(
    use.UsageAction(plugin_id='deblur',
                    action_id='denoise_16S'),
    use.UsageInputs(demultiplexed_seqs=filtered_seqs,
                    trim_length=120,
                    sample_stats=True),
    use.UsageOutputNames(representative_sequences='rep_seqs_deblur',
                         table='table_deblur',
                         stats='deblur_stats'))
```

* `rep-seqs-deblur.qza` | [download](https://moving-pictures-tutorial.readthedocs.io/en/latest/data/moving-pictures/rep-seqs-deblur.qza) | [view](https://view.qiime2.org?src=https%3A%2F%2Fmoving-pictures-tutorial.readthedocs.io%2Fen%2Flatest%2Fdata%2Fmoving-pictures%2Frep-seqs-deblur.qza)
* `table-deblur.qza` | [download](https://moving-pictures-tutorial.readthedocs.io/en/latest/data/moving-pictures/table-deblur.qza) | [view](https://view.qiime2.org?src=https%3A%2F%2Fmoving-pictures-tutorial.readthedocs.io%2Fen%2Flatest%2Fdata%2Fmoving-pictures%2Ftable-deblur.qza)
* `deblur-stats.qza` | [download](https://moving-pictures-tutorial.readthedocs.io/en/latest/data/moving-pictures/deblur-stats.qza) | [view](https://view.qiime2.org?src=https%3A%2F%2Fmoving-pictures-tutorial.readthedocs.io%2Fen%2Flatest%2Fdata%2Fmoving-pictures%2Fdeblur-stats.qza)

The two commands used in this section generate QIIME 2 artifacts containing summary statistics.
To view those summary statistics, you can visualize them using `qiime metadata tabulate` and `qiime deblur visualize-stats`, respectively:

[Command Line]

[Python API]

[Galaxy]

[R API]

[View Source]

```
qiime metadata tabulate \
  --m-input-file demux-filter-stats.qza \
  --o-visualization demux-filter-stats.qzv
qiime deblur visualize-stats \
  --i-deblur-stats deblur-stats.qza \
  --o-visualization deblur-stats.qzv
```

```
filter_stats_md = demux_filter_stats.view(Metadata)
demux_filter_stats_viz, = metadata_actions.tabulate(
    input=filter_stats_md,
)
deblur_stats_viz, = deblur_actions.visualize_stats(
    deblur_stats=deblur_stats,
)
```

Using the `qiime2 metadata tabulate` tool:
:   1. For *"input"*:
       * Perform the following steps.
         1. Change to `Metadata from Artifact`
         2. Set *"Metadata Source"* to `demux-filter-stats.qza`
    2. Press the `Execute` button.

Once completed, for the new entry in your history, use the `Edit` button to set the name as follows:
:   (Renaming is optional, but it will make any subsequent steps easier to complete.)

    | History Name | "Name" to set (be sure to press [Save]) |
    | --- | --- |
    | `#: qiime2 metadata tabulate [...] : visualization.qzv` | `demux-filter-stats.qzv` |

Using the `qiime2 deblur visualize-stats` tool:
:   1. Set *"deblur\_stats"* to `#: deblur-stats.qza`
    2. Press the `Execute` button.

Once completed, for the new entry in your history, use the `Edit` button to set the name as follows:
:   (Renaming is optional, but it will make any subsequent steps easier to complete.)

    | History Name | "Name" to set (be sure to press [Save]) |
    | --- | --- |
    | `#: qiime2 deblur visualize-stats [...] : visualization.qzv` | `deblur-stats.qzv` |

```
filter_stats_md <- demux_filter_stats$view(Metadata)
action_results <- metadata_actions$tabulate(
    input=filter_stats_md,
)
demux_filter_stats_viz <- action_results$visualization
action_results <- deblur_actions$visualize_stats(
    deblur_stats=deblur_stats,
)
deblur_stats_viz <- action_results$visualization
```

```
filter_stats_as_md = use.view_as_metadata('filter_stats', filter_stats)

use.action(
    use.UsageAction(plugin_id='metadata',
                    action_id='tabulate'),
    use.UsageInputs(input=filter_stats_as_md),
    use.UsageOutputNames(visualization='demux_filter_stats'))

use.action(
    use.UsageAction(plugin_id='deblur',
                    action_id='visualize_stats'),
    use.UsageInputs(deblur_stats=stats_deblur),
    use.UsageOutputNames(visualization='deblur_stats'))
```

* `demux-filter-stats.qzv` | [download](https://moving-pictures-tutorial.readthedocs.io/en/latest/data/moving-pictures/demux-filter-stats.qzv) | [view](https://view.qiime2.org?src=https%3A%2F%2Fmoving-pictures-tutorial.readthedocs.io%2Fen%2Flatest%2Fdata%2Fmoving-pictures%2Fdemux-filter-stats.qzv)
* `deblur-stats.qzv` | [download](https://moving-pictures-tutorial.readthedocs.io/en/latest/data/moving-pictures/deblur-stats.qzv) | [view](https://view.qiime2.org?src=https%3A%2F%2Fmoving-pictures-tutorial.readthedocs.io%2Fen%2Flatest%2Fdata%2Fmoving-pictures%2Fdeblur-stats.qzv)

If you’d like to continue the tutorial using this `FeatureTable` (as opposed to the DADA2 feature table generated in *Option 1*), run the following commands.

[Command Line]

[Python API]

```
mv rep-seqs-deblur.qza rep-seqs.qza
mv table-deblur.qza table.qza
```

```
table = table_deblur
rep_seqs = rep_seqs_deblur
```

## FeatureTable and FeatureData summaries(#fKUK4bsiQh "Link to this Section")

After the quality filtering step completes, you’ll want to explore the resulting data.
You can do this using the following two commands, which will create visual summaries of the data.
The `feature-table summarize` command will give you information on how many sequences are associated with each sample and with each feature, histograms of those distributions, and some related summary statistics.
The `feature-table tabulate-seqs` command will provide a mapping of feature IDs to sequences, and provide links to easily BLAST each sequence against the NCBI nt database.
The latter visualization will be very useful later in the tutorial, when you want to learn more about specific features that are important in the data set.

[Command Line]

[View Source]

```
qiime feature-table summarize \
  --i-table table.qza \
  --m-sample-metadata-file sample-metadata.tsv \
  --o-visualization table.qzv
qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv
```

```
use.action(
    use.UsageAction(plugin_id='feature_table',
                    action_id='summarize'),
    use.UsageInputs(table=table_dada2,
                    sample_metadata=sample_metadata),
    use.UsageOutputNames(visualization='table'))

use.action(
    use.UsageAction(plugin_id='feature_table',
                    action_id='tabulate_seqs'),
    use.UsageInputs(data=rep_seqs_dada2),
    use.UsageOutputNames(visualization='rep_seqs'),
)
```

```
Traceback (most recent call last):
  File "/home/docs/checkouts/readthedocs.org/user_builds/moving-pictures-tutorial/checkouts/latest/.env/lib/python3.10/site-packages/q2doc/transforms/transform_usage.py", line 92, in run
    exec(source, interface['driver'].scope)
  File "<string>", line 4, in <module>
NameError: name 'table_dada2' is not defined
```

## Generate a tree for phylogenetic diversity analyses(#lZcW4o56U3 "Link to this Section")

QIIME supports several phylogenetic diversity metrics, including Faith’s Phylogenetic Diversity and weighted and unweighted UniFrac.
In addition to counts of features per sample (i.e., the data in the `FeatureTable[Frequency]` QIIME 2 artifact), these metrics require a rooted phylogenetic tree relating the features to one another.
This information will be stored in a `Phylogeny[Rooted]` QIIME 2 artifact.
To generate a phylogenetic tree we will use `align-to-tree-mafft-fasttree` pipeline from the `q2-phylogeny` plugin.

First, the pipeline uses the `mafft` program to perform a multiple sequence alignment of the sequences in our `FeatureData[Sequence]` to create a `FeatureData[AlignedSequence]` QIIME 2 artifact.
Next, the pipeline masks (or filters) the alignment to remove positions that are highly variable.
These positions are generally considered to add noise to a resulting phylogenetic tree.
Following that, the pipeline applies FastTree to generate a phylogenetic tree from the masked alignment.
The FastTree program creates an unrooted tree, so in the final step in this section midpoint rooting is applied to place the root of the tree at the midpoint of the longest tip-to-tip distance in the unrooted tree.

[Command Line]

[View Source]

```
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza
```

```
_, _, _, rooted_tree = use.action(
    use.UsageAction(plugin_id='phylogeny',
                    action_id='align_to_tree_mafft_fasttree'),
    use.UsageInputs(sequences=rep_seqs_dada2),
    use.UsageOutputNames(alignment='aligned_rep_seqs',
                         masked_alignment='masked_aligned_rep_seqs',
                         tree='unrooted_tree',
                         rooted_tree='rooted_tree'))
```

```
Traceback (most recent call last):
  File "/home/docs/checkouts/readthedocs.org/user_builds/moving-pictures-tutorial/checkouts/latest/.env/lib/python3.10/site-packages/q2doc/transforms/transform_usage.py", line 92, in run
    exec(source, interface['driver'].scope)
  File "<string>", line 4, in <module>
NameError: name 'rep_seqs_dada2' is not defined. Did you mean: 'rep_seqs_deblur'?
```

## Alpha and beta diversity analysis(#ovLK34Yjq3 "Link to this Section")

QIIME 2’s diversity analyses are available through the `q2-diversity` plugin, which supports computing alpha and beta diversity metrics, applying related statistical tests, and generating interactive visualizations.
We’ll first apply the `core-metrics-phylogenetic` method, which rarefies a `FeatureTable[Frequency]` to a user-specified depth, computes several alpha and beta diversity metrics, and generates principle coordinates analysis (PCoA) plots using Emperor for each of the beta diversity metrics.
The metrics computed by default are:

* Alpha diversity
  + Shannon’s diversity index (a quantitative measure of community
    richness)
  + Observed Features (a qualitative measure of community richness)
  + Faith’s Phylogenetic Diversity (a qualitative measure of
    community richness that incorporates phylogenetic relationships
    between the features)
  + Evenness (or Pielou’s Evenness; a measure of community
    evenness)
* Beta diversity
  + Jaccard distance (a qualitative measure of community
    dissimilarity)
  + Bray-Curtis distance (a quantitative measure of community
    dissimilarity)
  + unweighted UniFrac distance (a qualitative measure of community
    dissimilarity that incorporates phylogenetic relationships
    between the features)
  + weighted UniFrac distance (a quantitative measure of community
    dissimilarity that incorporates phylogenetic relationships
    between the features)

An important parameter that needs to be provided to this script is `--p-sampling-depth`, which is the even sampling (i.e. rarefaction) depth.
Because most diversity metrics are sensitive to different sampling depths across different samples, this script will randomly subsample the counts from each sample to the value provided for this parameter.
For example, if you provide `--p-sampling-depth 500`, this step will subsample the counts in each sample without replacement so that each sample in the resulting table has a total count of 500.
If the total count for any sample(s) are smaller than this value, those samples will be dropped from the diversity analysis.
Choosing this value is tricky.
We recommend making your choice by reviewing the information presented in the `table.qzv` file that was created above.
Choose a value that is as high as possible (so you retain more sequences per sample) while excluding as few samples as possible.

[Command Line]

[View Source]

```
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table table.qza \
  --p-sampling-depth 1103 \
  --m-metadata-file sample-metadata.tsv \
  --output-dir diversity-core-metrics-phylogenetic
```

```
core_metrics_results = use.action(
    use.UsageAction(plugin_id='diversity',
                    action_id='core_metrics_phylogenetic'),
    use.UsageInputs(phylogeny=rooted_tree,
                    table=table_dada2,
                    sampling_depth=1103,
                    metadata=sample_metadata),
    use.UsageOutputNames(rarefied_table='rarefied_table',
                         faith_pd_vector='faith_pd_vector',
                         observed_features_vector='observed_features_vector',
                         shannon_vector='shannon_vector',
                         evenness_vector='evenness_vector',
                         unweighted_unifrac_distance_matrix='unweighted_unifrac_distance_matrix',
                         weighted_unifrac_distance_matrix='weighted_unifrac_distance_matrix',
                         jaccard_distance_matrix='jaccard_distance_matrix',
                         bray_curtis_distance_matrix='bray_curtis_distance_matrix',
                         unweighted_unifrac_pcoa_results='unweighted_unifrac_pcoa_results',
                         weighted_unifrac_pcoa_results='weighted_unifrac_pcoa_results',
                         jaccard_pcoa_results='jaccard_pcoa_results',
                         bray_curtis_pcoa_results='bray_curtis_pcoa_results',
                         unweighted_unifrac_emperor='unweighted_unifrac_emperor',
                         weighted_unifrac_emperor='weighted_unifrac_emperor',
                         jaccard_emperor='jaccard_emperor',
                         bray_curtis_emperor='bray_curtis_emperor'))
faith_pd_vec = core_metrics_results.faith_pd_vector
evenness_vec = core_metrics_results.evenness_vector
unweighted_unifrac_dm = core_metrics_results.unweighted_unifrac_distance_matrix
unweighted_unifrac_pcoa = core_metrics_results.unweighted_unifrac_pcoa_results
bray_curtis_pcoa=core_metrics_results.bray_curtis_pcoa_results
```

```
Traceback (most recent call last):
  File "/home/docs/checkouts/readthedocs.org/user_builds/moving-pictures-tutorial/checkouts/latest/.env/lib/python3.10/site-packages/q2doc/transforms/transform_usage.py", line 92, in run
    exec(source, interface['driver'].scope)
  File "<string>", line 4, in <module>
NameError: name 'rooted_tree' is not defined
```

Here we set the `--p-sampling-depth` parameter to 1103.
This value was chosen based on the number of sequences in the `L3S313` sample because it’s close to the number of sequences in the next few samples that have higher sequence counts, and because it is considerably higher (relatively) than the number of sequences in the samples that have fewer sequences.
This will allow us to retain most of our samples.
The three samples that have fewer sequences will be dropped from the `core-metrics-phylogenetic` analyses and anything that uses these results.
It is worth noting that all three of these samples are “right palm” samples.
Losing a disproportionate number of samples from one metadata category is not ideal.
However, we are dropping a small enough number of samples here that this felt like the best compromise between total sequences analyzed and number of samples retained.

After computing diversity metrics, we can begin to explore the microbial composition of the samples in the context of the sample metadata.
This information is present in the [sample metadata](https://data.qiime2.org/2025.4/tutorials/moving-pictures/sample_metadata) file that was downloaded earlier.

We’ll first test for associations between categorical metadata columns and alpha diversity data.
We’ll do that here for the Faith Phylogenetic Diversity (a measure of community richness) and evenness metrics.

[Command Line]

[View Source]

```
qiime diversity alpha-group-significance \
  --i-alpha-diversity diversity-core-metrics-phylogenetic/faith_pd_vector.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization faith-pd-group-significance.qzv
qiime diversity alpha-group-significance \
  --i-alpha-diversity diversity-core-metrics-phylogenetic/evenness_vector.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization evenness-group-significance.qzv
```

```
use.action(
    use.UsageAction(plugin_id='diversity',
                    action_id='alpha_group_significance'),
    use.UsageInputs(alpha_diversity=faith_pd_vec,
                    metadata=sample_metadata),
    use.UsageOutputNames(visualization='faith_pd_group_significance'))

use.action(
    use.UsageAction(plugin_id='diversity',
                    action_id='alpha_group_significance'),
    use.UsageInputs(alpha_diversity=evenness_vec,
                    metadata=sample_metadata),
    use.UsageOutputNames(visualization='evenness_group_significance'))
```

```
Traceback (most recent call last):
  File "/home/docs/checkouts/readthedocs.org/user_builds/moving-pictures-tutorial/checkouts/latest/.env/lib/python3.10/site-packages/q2doc/transforms/transform_usage.py", line 92, in run
    exec(source, interface['driver'].scope)
  File "<string>", line 4, in <module>
NameError: name 'faith_pd_vec' is not defined
```

In this data set, no continuous sample metadata columns (e.g., `days-since-experiment-start`) are correlated with alpha diversity, so we won’t test for those associations here.
If you’re interested in performing those tests (for this data set, or for others), you can use the `qiime diversity alpha-correlation` command.

Next we’ll analyze sample composition in the context of categorical metadata using PERMANOVA (first described in Anderson (2001)) using the `beta-group-significance` command.
The following commands will test whether distances between samples within a group, such as samples from the same body site (e.g., gut), are more similar to each other then they are to samples from the other groups (e.g., tongue, left palm, and right palm).
If you call this command with the `--p-pairwise` parameter, as we’ll do here, it will also perform pairwise tests that will allow you to determine which specific pairs of groups (e.g., tongue and gut) differ from one another, if any.
This command can be slow to run, especially when passing `--p-pairwise`, since it is based on permutation tests.
So, unlike the previous commands, we’ll run `beta-group-significance` on specific columns of metadata that we’re interested in exploring, rather than all metadata columns to which it is applicable.
Here we’ll apply this to our unweighted UniFrac distances, using two sample metadata columns, as follows.

[Command Line]

[View Source]

```
qiime diversity beta-group-significance \
  --i-distance-matrix diversity-core-metrics-phylogenetic/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file sample-metadata.tsv \
  --m-metadata-column body-site \
  --p-pairwise \
  --o-visualization unweighted-unifrac-body-site-group-significance.qzv
qiime diversity beta-group-significance \
  --i-distance-matrix diversity-core-metrics-phylogenetic/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file sample-metadata.tsv \
  --m-metadata-column subject \
  --p-pairwise \
  --o-visualization unweighted-unifrac-subject-group-significance.qzv
```

```
body_site_col = use.get_metadata_column('body_site', 'body-site', sample_metadata)

use.action(
    use.UsageAction(plugin_id='diversity',
                    action_id='beta_group_significance'),
    use.UsageInputs(distance_matrix=unweighted_unifrac_dm,
                    metadata=body_site_col, pairwise=True),
    use.UsageOutputNames(visualization='unweighted_unifrac_body_site_group_significance'))

subject_col = use.get_metadata_column('subject', 'subject', sample_metadata)

use.action(
    use.UsageAction(plugin_id='diversity',
                    action_id='beta_group_significance'),
    use.UsageInputs(distance_matrix=unweighted_unifrac_dm,
                    metadata=subject_col,
                    pairwise=True),
    use.UsageOutputNames(visualization='unweighted_unifrac_subject_group_significance'))
```

```
Traceback (most recent call last):
  File "/home/docs/checkouts/readthedocs.org/user_builds/moving-pictures-tutorial/checkouts/latest/.env/lib/python3.10/site-packages/q2doc/transforms/transform_usage.py", line 92, in run
    exec(source, interface['driver'].scope)
  File "<string>", line 6, in <module>
NameError: name 'unweighted_unifrac_dm' is not defined
```

Again, none of the continuous sample metadata that we have for this data set are correlated with sample composition, so we won’t test for those associations here.
If you’re interested in performing those tests, you can use the `qiime metadata distance-matrix` in combination with `qiime diversity mantel` and `qiime diversity bioenv` commands.

Finally, ordination is a popular approach for exploring microbial community composition in the context of sample metadata.
We can use the Emperor tool to explore principal coordinates (PCoA) plots in the context of sample metadata.
While our `core-metrics-phylogenetic` command did already generate some Emperor plots, we want to pass an optional parameter, `--p-custom-axes`, which is very useful for exploring time series data.
The PCoA results that were used in `core-metrics-phylogeny` are also available, making it easy to generate new visualizations with Emperor.
We will generate Emperor plots for unweighted UniFrac and Bray-Curtis so that the resulting plot will contain axes for principal coordinate 1, principal coordinate 2, and days since the experiment start.
We will use that last axis to explore how these samples changed over time.

[Command Line]

[View Source]

```
qiime emperor plot \
  --i-pcoa diversity-core-metrics-phylogenetic/unweighted_unifrac_pcoa_results.qza \
  --m-metadata-file sample-metadata.tsv \
  --p-custom-axes days-since-experiment-start \
  --o-visualization unweighted-unifrac-emperor-days-since-experiment-start.qzv
qiime emperor plot \
  --i-pcoa diversity-core-metrics-phylogenetic/bray_curtis_pcoa_results.qza \
  --m-metadata-file sample-metadata.tsv \
  --p-custom-axes days-since-experiment-start \
  --o-visualization bray-curtis-emperor-days-since-experiment-start.qzv
```

```
use.action(
    use.UsageAction(plugin_id='emperor',
                    action_id='plot'),
    use.UsageInputs(pcoa=unweighted_unifrac_pcoa,
                    metadata=sample_metadata,
                    custom_axes=['days-since-experiment-start']),
    use.UsageOutputNames(visualization='unweighted-unifrac-emperor-days-since-experiment-start'))

use.action(
    use.UsageAction(plugin_id='emperor',
                    action_id='plot'),
    use.UsageInputs(pcoa=bray_curtis_pcoa,
                    metadata=sample_metadata,
                    custom_axes=['days-since-experiment-start']),
    use.UsageOutputNames(visualization='bray-curtis-emperor-days-since-experiment-start'))
```

```
Traceback (most recent call last):
  File "/home/docs/checkouts/readthedocs.org/user_builds/moving-pictures-tutorial/checkouts/latest/.env/lib/python3.10/site-packages/q2doc/transforms/transform_usage.py", line 92, in run
    exec(source, interface['driver'].scope)
  File "<string>", line 4, in <module>
NameError: name 'unweighted_unifrac_pcoa' is not defined
```

## Alpha rarefaction plotting(#StmGMuGJRY "Link to this Section")

In this section we’ll explore alpha diversity as a function of sampling depth using the `qiime diversity alpha-rarefaction` visualizer.
This visualizer computes one or more alpha diversity metrics at multiple sampling depths, in steps between 1 (optionally controlled with `--p-min-depth`) and the value provided as `--p-max-depth`.
At each sampling depth step, 10 rarefied tables will be generated, and the diversity metrics will be computed for all samples in the tables.
The number of iterations (rarefied tables computed at each sampling depth) can be controlled with `--p-iterations`.
Average diversity values will be plotted for each sample at each even sampling depth, and samples can be grouped based on metadata in the resulting visualization if sample metadata is provided with the `--m-metadata-file` parameter.

[Command Line]

[View Source]

```
qiime diversity alpha-rarefaction \
  --i-table table.qza \
  --i-phylogeny rooted-tree.qza \
  --p-max-depth 4000 \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization alpha-rarefaction.qzv
```

```
use.action(
    use.UsageAction(plugin_id='diversity',
                    action_id='alpha_rarefaction'),
    use.UsageInputs(table=table_dada2,
                    phylogeny=rooted_tree,
                    max_depth=4000,
                    metadata=sample_metadata),
    use.UsageOutputNames(visualization='alpha_rarefaction'))
```

```
Traceback (most recent call last):
  File "/home/docs/checkouts/readthedocs.org/user_builds/moving-pictures-tutorial/checkouts/latest/.env/lib/python3.10/site-packages/q2doc/transforms/transform_usage.py", line 92, in run
    exec(source, interface['driver'].scope)
  File "<string>", line 4, in <module>
NameError: name 'table_dada2' is not defined
```

The visualization will have two plots.
The top plot is an alpha rarefaction plot, and is primarily used to determine if the richness of the samples has been fully observed or sequenced.
If the lines in the plot appear to “level out” (i.e., approach a slope of zero) at some sampling depth along the x-axis, that suggests that collecting additional sequences beyond that sampling depth would not be likely to result in the observation of additional features.
If the lines in a plot don’t level out, this may be because the richness of the samples hasn’t been fully observed yet (because too few sequences were collected), or it could be an indicator that a lot of sequencing error remains in the data (which is being mistaken for novel diversity).

The bottom plot in this visualization is important when grouping samples by metadata.
It illustrates the number of samples that remain in each group when the feature table is rarefied to each sampling depth.
If a given sampling depth `d` is larger than the total frequency of a sample `s` (i.e., the number of sequences that were obtained for sample `s`), it is not possible to compute the diversity metric for sample `s` at sampling depth `d`.
If many of the samples in a group have lower total frequencies than `d`, the average diversity presented for that group at `d` in the top plot will be unreliable because it will have been computed on relatively few samples.
When grouping samples by metadata, it is therefore essential to look at the bottom plot to ensure that the data presented in the top plot is reliable.

## Taxonomic analysis(#mK1xGXss11 "Link to this Section")

In the next sections we’ll begin to explore the taxonomic composition of the samples, and again relate that to sample metadata.
The first step in this process is to assign taxonomy to the sequences in our `FeatureData[Sequence]` QIIME 2 artifact.
We’ll do that using a pre-trained Naive Bayes classifier and the `q2-feature-classifier` plugin.
This classifier was trained on the Greengenes 13\_8 99% OTUs, where the sequences have been trimmed to only include 250 bases from the region of the 16S that was sequenced in this analysis (the V4 region, bound by the 515F/806R primer pair).
We’ll apply this classifier to our sequences, and we can generate a visualization of the resulting mapping from sequence to taxonomy.

[Command Line]

[Python API]

[Galaxy]

[R API]

[View Source]

```
wget -O 'gg-13-8-99-515-806-nb-classifier.qza' \
  'https://moving-pictures-tutorial.readthedocs.io/en/latest/data/moving-pictures/gg-13-8-99-515-806-nb-classifier.qza'
```

```
body_site_mdc = sample_metadata_md.get_column('body-site')
url = 'https://moving-pictures-tutorial.readthedocs.io/en/latest/data/moving-pictures/gg-13-8-99-515-806-nb-classifier.qza'
fn = 'gg-13-8-99-515-806-nb-classifier.qza'
request.urlretrieve(url, fn)
gg_13_8_99_515_806_nb_classifier = Artifact.load(fn)
```

Using the `Upload Data` tool:
:   1. On the first tab (**Regular**), press the `Paste/Fetch` data button at the bottom.
       1. Set *"Name"* (first text-field) to: `gg-13-8-99-515-806-nb-classifier.qza`
       2. In the larger text-area, copy-and-paste: <https://moving-pictures-tutorial.readthedocs.io/en/latest/data/moving-pictures/gg-13-8-99-515-806-nb-classifier.qza>
       3. (*"Type"*, *"Genome"*, and *"Settings"* can be ignored)
    2. Press the `Start` button at the bottom.

```
url <- 'https://moving-pictures-tutorial.readthedocs.io/en/latest/data/moving-pictures/gg-13-8-99-515-806-nb-classifier.qza'
fn <- 'gg-13-8-99-515-806-nb-classifier.qza'
request$urlretrieve(url, fn)
gg_13_8_99_515_806_nb_classifier <- Artifact$load(fn)
```

```
def classifier_factory():
    from urllib import request
    from qiime2 import Artifact
    fp, _ = request.urlretrieve(
        'https://data.qiime2.org/classifiers/sklearn-1.4.2/greengenes/gg-13-8-99-515-806-nb-classifier.qza')

    return Artifact.load(fp)

classifier = use.init_artifact('gg-13-8-99-515-806-nb-classifier', classifier_factory)
```

* `gg-13-8-99-515-806-nb-classifier.qza` | [download](https://moving-pictures-tutorial.readthedocs.io/en/latest/data/moving-pictures/gg-13-8-99-515-806-nb-classifier.qza) | [view](https://view.qiime2.org?src=https%3A%2F%2Fmoving-pictures-tutorial.readthedocs.io%2Fen%2Flatest%2Fdata%2Fmoving-pictures%2Fgg-13-8-99-515-806-nb-classifier.qza)

[Command Line]

[View Source]

```
qiime feature-classifier classify-sklearn \
  --i-classifier gg-13-8-99-515-806-nb-classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza
qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv
```

```
taxonomy, = use.action(
    use.UsageAction(plugin_id='feature_classifier',
                    action_id='classify_sklearn'),
    use.UsageInputs(classifier=classifier,
                    reads=rep_seqs_dada2),
    use.UsageOutputNames(classification='taxonomy'))

taxonomy_as_md = use.view_as_metadata('taxonomy_as_md', taxonomy)

use.action(
    use.UsageAction(plugin_id='metadata',
                    action_id='tabulate'),
    use.UsageInputs(input=taxonomy_as_md),
    use.UsageOutputNames(visualization='taxonomy'))
```

```
Traceback (most recent call last):
  File "/home/docs/checkouts/readthedocs.org/user_builds/moving-pictures-tutorial/checkouts/latest/.env/lib/python3.10/site-packages/q2doc/transforms/transform_usage.py", line 92, in run
    exec(source, interface['driver'].scope)
  File "<string>", line 5, in <module>
NameError: name 'rep_seqs_dada2' is not defined. Did you mean: 'rep_seqs_deblur'?
```

Next, we can view the taxonomic composition of our samples with interactive bar plots.
Generate those plots with the following command and then open the visualization.

[Command Line]

[View Source]

```
qiime taxa barplot \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization taxa-bar-plots.qzv
```

```
use.action(
    use.UsageAction(plugin_id='taxa',
                    action_id='barplot'),
    use.UsageInputs(table=table_dada2,
                    taxonomy=taxonomy,
                    metadata=sample_metadata),
    use.UsageOutputNames(visualization='taxa_bar_plots'))
```

```
Traceback (most recent call last):
  File "/home/docs/checkouts/readthedocs.org/user_builds/moving-pictures-tutorial/checkouts/latest/.env/lib/python3.10/site-packages/q2doc/transforms/transform_usage.py", line 92, in run
    exec(source, interface['driver'].scope)
  File "<string>", line 4, in <module>
NameError: name 'table_dada2' is not defined
```

## Differential abundance testing with ANCOM-BC(#srGSTaILQA "Link to this Section")

ANCOM-BC can be applied to identify features that are differentially abundant (i.e. present in different abundances) across sample groups.
As with any bioinformatics method, you should be aware of the assumptions and limitations of ANCOM-BC before using it.
We recommend reviewing the ANCOM-BC paper before using this method.

ANCOM-BC is a compositionally-aware linear regression model that allows for testing differentially abundant features across groups while also implementing bias correction, and is currently implemented in the `q2-composition` plugin.

Because we expect a lot of features to change in abundance across body sites, in this tutorial we’ll filter our full feature table to only contain gut samples.
We’ll then apply ANCOM-BC to determine which, if any, sequence variants and genera are differentially abundant across the gut samples of our two subjects.

We’ll start by creating a feature table that contains only the gut samples.

[Command Line]

[View Source]

```
qiime feature-table filter-samples \
  --i-table table.qza \
  --m-metadata-file sample-metadata.tsv \
  --p-where '[body-site]="gut"' \
  --o-filtered-table gut-table.qza
```

```
gut_table, = use.action(
    use.UsageAction(plugin_id='feature_table',
                    action_id='filter_samples'),
    use.UsageInputs(table=table_dada2,
                    metadata=sample_metadata,
                    where='[body-site]="gut"'),
    use.UsageOutputNames(filtered_table='gut_table'))
```

```
Traceback (most recent call last):
  File "/home/docs/checkouts/readthedocs.org/user_builds/moving-pictures-tutorial/checkouts/latest/.env/lib/python3.10/site-packages/q2doc/transforms/transform_usage.py", line 92, in run
    exec(source, interface['driver'].scope)
  File "<string>", line 4, in <module>
NameError: name 'table_dada2' is not defined
```

ANCOM-BC operates on a `FeatureTable[Frequency]` QIIME 2 artifact.
We can run ANCOM-BC on the subject column to determine what features differ in abundance across gut samples of the two subjects.

[Command Line]

[View Source]

```
qiime composition ancombc \
  --i-table gut-table.qza \
  --m-metadata-file sample-metadata.tsv \
  --p-formula subject \
  --o-differentials ancombc-subject.qza
qiime composition da-barplot \
  --i-data ancombc-subject.qza \
  --p-significance-threshold 0.001 \
  --o-visualization da-barplot-subject.qzv
```

```
ancombc_subject, = use.action(
    use.UsageAction(plugin_id='composition',
                    action_id='ancombc'),
    use.UsageInputs(table=gut_table,
                    metadata=sample_metadata,
                    formula='subject'),
    use.UsageOutputNames(differentials='ancombc_subject'))

use.action(
    use.UsageAction(plugin_id='composition',
                    action_id='da_barplot'),
    use.UsageInputs(data=ancombc_subject,
                    significance_threshold=0.001),
    use.UsageOutputNames(visualization='da_barplot_subject'))
```

```
Traceback (most recent call last):
  File "/home/docs/checkouts/readthedocs.org/user_builds/moving-pictures-tutorial/checkouts/latest/.env/lib/python3.10/site-packages/q2doc/transforms/transform_usage.py", line 92, in run
    exec(source, interface['driver'].scope)
  File "<string>", line 4, in <module>
NameError: name 'gut_table' is not defined
```

We’re also often interested in performing a differential abundance test at a specific taxonomic level.
To do this, we can collapse the features in our `FeatureTable[Frequency]` at the taxonomic level of interest, and then re-run the above steps.
In this tutorial, we collapse our feature table at the genus level (i.e. level 6 of the Greengenes taxonomy).

[Command Line]

[View Source]

```
qiime taxa collapse \
  --i-table gut-table.qza \
  --i-taxonomy taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table gut-table-l6.qza
qiime composition ancombc \
  --i-table gut-table-l6.qza \
  --m-metadata-file sample-metadata.tsv \
  --p-formula subject \
  --o-differentials l6-ancombc-subject.qza
qiime composition da-barplot \
  --i-data l6-ancombc-subject.qza \
  --p-significance-threshold 0.001 \
  --o-visualization l6-da-barplot-subject.qzv
```

```
l6_gut_table, = use.action(
    use.UsageAction(plugin_id='taxa',
                    action_id='collapse'),
    use.UsageInputs(table=gut_table,
                    taxonomy=taxonomy,
                    level=6),
    use.UsageOutputNames(collapsed_table='gut_table_l6'))

l6_ancombc_subject, = use.action(
    use.UsageAction(plugin_id='composition',
                    action_id='ancombc'),
    use.UsageInputs(table=l6_gut_table,
                    metadata=sample_metadata,
                    formula='subject'),
    use.UsageOutputNames(differentials='l6_ancombc_subject'))

use.action(
    use.UsageAction(plugin_id='composition',
                    action_id='da_barplot'),
    use.UsageInputs(data=l6_ancombc_subject,
                    significance_threshold=0.001),
    use.UsageOutputNames(visualization='l6_da_barplot_subject'))
```

```
Traceback (most recent call last):
  File "/home/docs/checkouts/readthedocs.org/user_builds/moving-pictures-tutorial/checkouts/latest/.env/lib/python3.10/site-packages/q2doc/transforms/transform_usage.py", line 92, in run
    exec(source, interface['driver'].scope)
  File "<string>", line 4, in <module>
NameError: name 'gut_table' is not defined
```

## Next steps(#dw3ls6gVl9 "Link to this Section")

You’ve reached the end of the QIIME 2 *Moving Pictures* tutorial.
Many users who are new to the platform will next adapt the commands presented in this tutorial to their own data, so that might be a good next step for you.
If you need help, head over to the [QIIME 2 Forum](https://forum.qiime2.org).

[Microbiome marker gene analysis with QIIME 2

Microbiome marker gene analysis with QIIME 2](/en/latest/)[Tutorials

Gut-to-soil axis tutorial 💩🌱](/en/latest/tutorials/gut-to-soil.html)
