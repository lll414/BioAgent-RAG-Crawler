---
source_url: https://amplicon-docs.qiime2.org/en/latest/tutorials/gut-to-soil.html
download_date: 2025-11-20
---

Downloads

# Gut-to-soil axis tutorial 💩🌱

## Background(#M8WOP8paii "Link to this Section")

In this tutorial you’ll learn an end-to-end microbiome marker-gene data science workflow, building on data presented in Meilander *et al.* (2024): Upcycling Human Excrement: The Gut Microbiome to Soil Microbiome Axis.
The data used here is a subset (a single sequencing run) of that generated for the paper, specifically selected so that this tutorial can be run quickly on a personal computer.
The full data set for the paper can be found in the study’s Artifact Repository.
In the final step (**not yet written, as of 17 April 2025**), you’ll learn how to adapt the workflow for use in analyzing your own data using Provenance Replay.

The data used in this tutorial was generated using the Earth Microbiome Project protocol.
Specifically, the hypervariable region 4 (V4) of the 16S rRNA gene was amplified using the F515-R806 primers - a broad-coverage primer pair for Bacteria that also amplifies some Archaea.
Paired-end sequencing was performed on an Illumina MiSeq.
Full details are presented in Meilander *et al.* (2024).

## Sample metadata(#FaxCmiQaO1 "Link to this Section")

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
  'https://gut-to-soil-tutorial.readthedocs.io/en/latest/data/gut-to-soil/sample-metadata.tsv'
```

```
from qiime2 import Metadata
from urllib import request

url = 'https://gut-to-soil-tutorial.readthedocs.io/en/latest/data/gut-to-soil/sample-metadata.tsv'
fn = 'sample-metadata.tsv'
request.urlretrieve(url, fn)
sample_metadata_md = Metadata.load(fn)
```

Using the `Upload Data` tool:
:   1. On the first tab (**Regular**), press the `Paste/Fetch` data button at the bottom.
       1. Set *"Name"* (first text-field) to: `sample-metadata.tsv`
       2. In the larger text-area, copy-and-paste: <https://gut-to-soil-tutorial.readthedocs.io/en/latest/data/gut-to-soil/sample-metadata.tsv>
       3. (*"Type"*, *"Genome"*, and *"Settings"* can be ignored)
    2. Press the `Start` button at the bottom.

```
library(reticulate)

Metadata <- import("qiime2")$Metadata
request <- import("urllib")$request

url <- 'https://gut-to-soil-tutorial.readthedocs.io/en/latest/data/gut-to-soil/sample-metadata.tsv'
fn <- 'sample-metadata.tsv'
request$urlretrieve(url, fn)
sample_metadata_md <- Metadata$load(fn)
```

```
sample_metadata = use.init_metadata_from_url(
   'sample-metadata',
   'https://zenodo.org/records/15390940/files/gut-to-soil-tutorial-sample-metadata.tsv?download=1')
```

* `sample-metadata.tsv` | [download](https://gut-to-soil-tutorial.readthedocs.io/en/latest/data/gut-to-soil/sample-metadata.tsv)

QIIME 2’s [metadata plugin](https://amplicon-docs.qiime2.org/en/latest/references/plugins/metadata#q2-plugin-metadata) provides a [Visualizer called `tabulate`](https://amplicon-docs.qiime2.org/en/latest/references/plugins/metadata#q2-action-metadata-tabulate) that generates a convenient view of a sample metadata file.
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
  --o-visualization sample-metadata.qzv
```

```
import qiime2.plugins.metadata.actions as metadata_actions

sample_metadata_viz, = metadata_actions.tabulate(
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
    | `#: qiime2 metadata tabulate [...] : visualization.qzv` | `sample-metadata.qzv` |

```
metadata_actions <- import("qiime2.plugins.metadata.actions")

action_results <- metadata_actions$tabulate(
    input=sample_metadata_md,
)
sample_metadata_viz <- action_results$visualization
```

```
use.action(
  use.UsageAction(plugin_id='metadata',
                  action_id='tabulate'),
  use.UsageInputs(input=sample_metadata),
  use.UsageOutputNames(visualization='sample_metadata')
)
```

* `sample-metadata.qzv` | [download](https://gut-to-soil-tutorial.readthedocs.io/en/latest/data/gut-to-soil/sample-metadata.qzv) | [view](https://view.qiime2.org?src=https%3A%2F%2Fgut-to-soil-tutorial.readthedocs.io%2Fen%2Flatest%2Fdata%2Fgut-to-soil%2Fsample-metadata.qzv)

This will generate a QIIME 2 [Visualization](https://amplicon-docs.readthedocs.io/en/latest/explanations/getting-started.html#getting-started-artifacts-and-visualizations).
Visualizations can be viewed by loading them with [QIIME 2 View](https://view.qiime2.org).
Navigate to QIIME 2 View, and drag and drop the visualization that was created to view it.

## Access already-imported QIIME 2 data(#z6eef3dRr2 "Link to this Section")

This tutorial begins with paired-end read sequencing data that has already been demultiplexed and imported into a QIIME 2 Artifact.
Because sequence data can be delivered to you in many different forms, it’s not possible to cover the varieties here.
Instead we refer you to [*How to import data for use with QIIME 2*](https://amplicon-docs.readthedocs.io/en/latest/how-to-guides/how-to-import.html) to learn how to import your data.
If you want to learn why importing is necessary, refer to [Why importing is necessary](https://amplicon-docs.readthedocs.io/en/latest/explanations/why-importing.html).

[Command Line]

[Python API]

[Galaxy]

[R API]

[View Source]

```
wget -O 'demux.qza' \
  'https://gut-to-soil-tutorial.readthedocs.io/en/latest/data/gut-to-soil/demux.qza'
```

```
from qiime2 import Artifact

url = 'https://gut-to-soil-tutorial.readthedocs.io/en/latest/data/gut-to-soil/demux.qza'
fn = 'demux.qza'
request.urlretrieve(url, fn)
demux = Artifact.load(fn)
```

Using the `Upload Data` tool:
:   1. On the first tab (**Regular**), press the `Paste/Fetch` data button at the bottom.
       1. Set *"Name"* (first text-field) to: `demux.qza`
       2. In the larger text-area, copy-and-paste: <https://gut-to-soil-tutorial.readthedocs.io/en/latest/data/gut-to-soil/demux.qza>
       3. (*"Type"*, *"Genome"*, and *"Settings"* can be ignored)
    2. Press the `Start` button at the bottom.

```
Artifact <- import("qiime2")$Artifact

url <- 'https://gut-to-soil-tutorial.readthedocs.io/en/latest/data/gut-to-soil/demux.qza'
fn <- 'demux.qza'
request$urlretrieve(url, fn)
demux <- Artifact$load(fn)
```

```
demux = use.init_artifact_from_url(
   'demux',
   'https://zenodo.org/records/15390940/files/gut-to-soil-tutorial-nano2-demux-10p.qza?download=1')
```

* `demux.qza` | [download](https://gut-to-soil-tutorial.readthedocs.io/en/latest/data/gut-to-soil/demux.qza) | [view](https://view.qiime2.org?src=https%3A%2F%2Fgut-to-soil-tutorial.readthedocs.io%2Fen%2Flatest%2Fdata%2Fgut-to-soil%2Fdemux.qza)

## Summarize demultiplexed sequences(#aGk1YDQ7Vf "Link to this Section")

When you have demultiplexed sequence data, the next step is typically to generate a visual summary of it.
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
import qiime2.plugins.demux.actions as demux_actions

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
demux_actions <- import("qiime2.plugins.demux.actions")

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

* `demux.qzv` | [download](https://gut-to-soil-tutorial.readthedocs.io/en/latest/data/gut-to-soil/demux.qzv) | [view](https://view.qiime2.org?src=https%3A%2F%2Fgut-to-soil-tutorial.readthedocs.io%2Fen%2Flatest%2Fdata%2Fgut-to-soil%2Fdemux.qzv)

## Upstream data analysis(#ch0dM5aUvM "Link to this Section")

Generally, the term “upstream” is used to refer to data analysis pre-feature-asv\_table, and “downstream” is used to refer to data analysis post-feature-asv\_table.
Let’s jump into our upstream analysis.

### Sequence quality control and feature table construction(#wxHGadXHVj "Link to this Section")

QIIME 2 plugins are available for several quality control methods, including DADA2, Deblur, and basic quality-score-based filtering.
In this tutorial we present this step using [DADA2](https://www.ncbi.nlm.nih.gov/pubmed/27214047).
The result of this method will be a `FeatureTable[Frequency]` QIIME 2 artifact, which contains counts (frequencies) of each unique sequence in each sample in the dataset, and a `FeatureData[Sequence]` QIIME 2 artifact, which maps feature identifiers in the `FeatureTable` to the sequences they represent.

DADA2 is a pipeline for detecting and correcting (where possible) Illumina amplicon sequence data.
As implemented in the [dada2 plugin](https://amplicon-docs.qiime2.org/en/latest/references/plugins/dada2#q2-plugin-dada2), this quality control process will additionally filter any phiX reads (commonly present in marker gene Illumina sequence data) that are identified in the sequencing data, filter chimeric sequences, and merge paired end reads.

The [`denoise-paired` action](https://amplicon-docs.qiime2.org/en/latest/references/plugins/dada2#q2-action-dada2-denoise-paired), which we’ll use here, requires four parameters that are used in quality filtering:

* `trim-left-f a`, which trims off the first `a` bases of each forward read
* `trunc-len-f b` which truncates each forward read at position `b`
* `trim-left-r c`, which trims off the first `c` bases of each forward read
* `trunc-len-r d` which truncates each forward read at position `d`
  This allows the user to remove low quality regions of the sequences.
  To determine what values to pass for these parameters, you should review the *Interactive Quality Plot* tab in the [`demux.qzv`](https://gut-to-soil-tutorial.readthedocs.io/en/latest/gut-to-soil-16s-tutorial#demux-summary-viz) file that was generated above.

Solution to [Exercise 2](https://gut-to-soil-tutorial.readthedocs.io/en/latest/gut-to-soil-16s-tutorial#dada2-trim-trunc)

The quality of the initial bases seems to be high, so I choose not to trim any bases from the beginning of the sequences.
The quality also seems good all the way out to the end, though maybe dropping off after 250 bases.
I’ll therefore truncate at 250.
I’ll keep these values the same for both the forward and reverse reads, though that is not a requirement.

Now run your DADA2 command.
This step may take up to 10 minutes to complete - it’s the longest running step in this tutorial.

[Command Line]

[View Source]

```
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux.qza \
  --p-trim-left-f 0 \
  --p-trunc-len-f 250 \
  --p-trim-left-r 0 \
  --p-trunc-len-r 250 \
  --o-representative-sequences asv-seqs.qza \
  --o-table asv-table.qza \
  --o-denoising-stats stats.qza
```

```
asv_seqs, asv_table, stats = use.action(
    use.UsageAction(plugin_id='dada2',
                    action_id='denoise_paired'),
    use.UsageInputs(demultiplexed_seqs=demux,
                    trim_left_f=0,
                    trunc_len_f=250,
                    trim_left_r=0,
                    trunc_len_r=250),
    use.UsageOutputNames(representative_sequences='asv_seqs',
                         table='asv_table',
                         denoising_stats='stats'))
```

```
Traceback (most recent call last):
  File "/home/docs/checkouts/readthedocs.org/user_builds/gut-to-soil-tutorial/checkouts/latest/.env/lib/python3.10/site-packages/q2doc/transforms/transform_usage.py", line 92, in run
    exec(source, interface['driver'].scope)
  File "<string>", line 1, in <module>
  File "/home/docs/checkouts/readthedocs.org/user_builds/gut-to-soil-tutorial/checkouts/latest/.env/lib/python3.10/site-packages/qiime2/plugins.py", line 344, in action
    self._template_action(action, inputs, variables)
  File "/home/docs/checkouts/readthedocs.org/user_builds/gut-to-soil-tutorial/checkouts/latest/.env/lib/python3.10/site-packages/qiime2/plugins.py", line 363, in _template_action
    output_vars = self._template_outputs(action, variables)
  File "/home/docs/checkouts/readthedocs.org/user_builds/gut-to-soil-tutorial/checkouts/latest/.env/lib/python3.10/site-packages/qiime2/plugins.py", line 390, in _template_outputs
    variable = getattr(variables, output)
AttributeError: 'UsageOutputs' object has no attribute 'base_transition_stats'
```

One of the outputs created by DADA2 is a summary of the denoising run.
That is generated as an [Artifact](https://amplicon-docs.readthedocs.io/en/latest/explanations/getting-started.html#getting-started-artifacts-and-visualizations), so can’t be viewed directly.
However this is one of many QIIME 2 types that can be [viewed as Metadata](https://use.qiime2.org/en/latest/how-to-guides/artifacts-as-metadata.html) - a very powerful concept that we’ll use again later in this tutorial.
Learning to view artifacts as Metadata creates nearly infinite possibilities for how you can explore your microbiome data with QIIME 2.

Here, we’ll again use the [metadata plugins `tabulate` visualizer](https://amplicon-docs.qiime2.org/en/latest/references/plugins/metadata#q2-action-metadata-tabulate), but this time we’ll apply it to the DADA2 statistics.

[Command Line]

[View Source]

```
qiime metadata tabulate \
  --m-input-file stats.qza \
  --o-visualization stats.qzv
```

```
stats_as_md = use.view_as_metadata('stats_as_md', stats)

use.action(
    use.UsageAction(plugin_id='metadata',
                    action_id='tabulate'),
    use.UsageInputs(input=stats_as_md),
    use.UsageOutputNames(visualization='stats'))
```

```
Traceback (most recent call last):
  File "/home/docs/checkouts/readthedocs.org/user_builds/gut-to-soil-tutorial/checkouts/latest/.env/lib/python3.10/site-packages/q2doc/transforms/transform_usage.py", line 92, in run
    exec(source, interface['driver'].scope)
  File "<string>", line 1, in <module>
NameError: name 'stats' is not defined
```

Solution to [Exercise 4](https://gut-to-soil-tutorial.readthedocs.io/en/latest/gut-to-soil-16s-tutorial#merge-metadata)

[Command Line]

[View Source]

```
qiime metadata tabulate \
  --m-input-file sample-metadata.tsv stats.qza \
  --o-visualization sample-metadata-w-dada2-stats.qzv
```

```
sample_metadata_and_dada2_stats_md = use.merge_metadata('sample_metadata_and_dada2_stats_md', sample_metadata, stats_as_md)

use.action(
    use.UsageAction(plugin_id='metadata',
                    action_id='tabulate'),
    use.UsageInputs(input=sample_metadata_and_dada2_stats_md),
    use.UsageOutputNames(visualization='sample_metadata_w_dada2_stats'))
```

```
Traceback (most recent call last):
  File "/home/docs/checkouts/readthedocs.org/user_builds/gut-to-soil-tutorial/checkouts/latest/.env/lib/python3.10/site-packages/q2doc/transforms/transform_usage.py", line 92, in run
    exec(source, interface['driver'].scope)
  File "<string>", line 1, in <module>
NameError: name 'stats_as_md' is not defined
```

### Feature table and feature data summaries(#yGZ0wj8ibk "Link to this Section")

After DADA2 completes, you’ll want to explore the resulting data.
You can do this using the following two commands, which will create visual summaries of the data.
The [`feature-table summarize` action](https://amplicon-docs.qiime2.org/en/latest/references/plugins/feature-table#q2-action-feature-table-summarize) command will give you information on how many sequences are associated with each sample and with each feature, histograms of those distributions, and some related summary statistics.

[Command Line]

[View Source]

```
qiime feature-table summarize-plus \
  --i-table asv-table.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-summary asv-table.qzv \
  --o-sample-frequencies sample-frequencies.qza \
  --o-feature-frequencies asv-frequencies.qza
```

```
_, _, asv_frequencies = use.action(
    use.UsageAction(plugin_id='feature_table',
                    action_id='summarize_plus'),
    use.UsageInputs(table=asv_table,
                    metadata=sample_metadata),
    use.UsageOutputNames(summary='asv_table',
                         sample_frequencies='sample_frequencies',
                         feature_frequencies='asv_frequencies'))
```

```
Traceback (most recent call last):
  File "/home/docs/checkouts/readthedocs.org/user_builds/gut-to-soil-tutorial/checkouts/latest/.env/lib/python3.10/site-packages/q2doc/transforms/transform_usage.py", line 92, in run
    exec(source, interface['driver'].scope)
  File "<string>", line 4, in <module>
NameError: name 'asv_table' is not defined
```

The [`feature-table tabulate-seqs` action](https://amplicon-docs.qiime2.org/en/latest/references/plugins/feature-table#q2-action-feature-table-tabulate-seqs) command will provide a mapping of feature IDs to sequences, and provide links to easily BLAST each sequence against the NCBI nt database.
We can also include the feature frequency information in this visualization by passing it as metadata, similar to how we [merged metadata](https://gut-to-soil-tutorial.readthedocs.io/en/latest/gut-to-soil-16s-tutorial#merge-metadata) in the exercise above.
In this case, however, we’re looking a *feature metadata*, as opposed to *sample metadata*.
As far as QIIME 2 is concerned, there is no difference between these two - in our case, it’ll only be the identifiers that differ.

This visualization will be very useful later in the tutorial, when you want to learn more about specific features that are important in the data set.

[Command Line]

[View Source]

```
qiime feature-table tabulate-seqs \
  --i-data asv-seqs.qza \
  --m-metadata-file asv-frequencies.qza \
  --o-visualization asv-seqs.qzv
```

```
asv_frequencies_as_md = use.view_as_metadata('asv_frequencies_md',
                                                 asv_frequencies)

use.action(
    use.UsageAction(plugin_id='feature_table',
                    action_id='tabulate_seqs'),
    use.UsageInputs(data=asv_seqs,
                    metadata=asv_frequencies_as_md),
    use.UsageOutputNames(visualization='asv_seqs'),
)
```

```
Traceback (most recent call last):
  File "/home/docs/checkouts/readthedocs.org/user_builds/gut-to-soil-tutorial/checkouts/latest/.env/lib/python3.10/site-packages/q2doc/transforms/transform_usage.py", line 92, in run
    exec(source, interface['driver'].scope)
  File "<string>", line 2, in <module>
NameError: name 'asv_frequencies' is not defined
```

### Filtering features from a feature table(#zAWX75ujyJ "Link to this Section")

If you review the tabulated feature sequences, or the feature detail table of the feature table summary, you’ll notice that there are many sequences that are observed in only a single sample.
Let’s filter those out to reduce the number of sequences we’re working with - this will speed up several slower steps that are coming up.

This is a two-step process.
First we filter our feature table, and then we use the new feature table to filter our sequences to only the ones that are contained in the new table.

[Command Line]

[View Source]

```
qiime feature-table filter-features \
  --i-table asv-table.qza \
  --p-min-samples 2 \
  --o-filtered-table asv-table-ms2.qza
```

```
asv_table_ms2, = use.action(
    use.UsageAction(plugin_id='feature_table',
                    action_id='filter_features'),
    use.UsageInputs(table=asv_table,
                    min_samples=2),
    use.UsageOutputNames(filtered_table='asv_table_ms2'),
)
```

```
Traceback (most recent call last):
  File "/home/docs/checkouts/readthedocs.org/user_builds/gut-to-soil-tutorial/checkouts/latest/.env/lib/python3.10/site-packages/q2doc/transforms/transform_usage.py", line 92, in run
    exec(source, interface['driver'].scope)
  File "<string>", line 4, in <module>
NameError: name 'asv_table' is not defined
```

[Command Line]

[View Source]

```
qiime feature-table filter-seqs \
  --i-data asv-seqs.qza \
  --i-table asv-table-ms2.qza \
  --o-filtered-data asv-seqs-ms2.qza
```

```
asv_seqs_ms2, = use.action(
    use.UsageAction(plugin_id='feature_table',
                    action_id='filter_seqs'),
    use.UsageInputs(data=asv_seqs,
                    table=asv_table_ms2),
    use.UsageOutputNames(filtered_data='asv_seqs_ms2'),
)
```

```
Traceback (most recent call last):
  File "/home/docs/checkouts/readthedocs.org/user_builds/gut-to-soil-tutorial/checkouts/latest/.env/lib/python3.10/site-packages/q2doc/transforms/transform_usage.py", line 92, in run
    exec(source, interface['driver'].scope)
  File "<string>", line 4, in <module>
NameError: name 'asv_seqs' is not defined
```

Solution to [Exercise 7](https://gut-to-soil-tutorial.readthedocs.io/en/latest/gut-to-soil-16s-tutorial#summarize-asv-table-ms2)

Here’s the command you would use:

[Command Line]

[View Source]

```
qiime feature-table summarize-plus \
  --i-table asv-table-ms2.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-summary asv-table-ms2.qzv \
  --o-sample-frequencies sample-frequencies-ms2.qza \
  --o-feature-frequencies asv-frequencies-ms2.qza
```

```
_, _, asv_frequencies_ms2 = use.action(
    use.UsageAction(plugin_id='feature_table',
                    action_id='summarize_plus'),
    use.UsageInputs(table=asv_table_ms2,
                    metadata=sample_metadata),
    use.UsageOutputNames(summary='asv_table_ms2',
                         sample_frequencies='sample_frequencies_ms2',
                         feature_frequencies='asv_frequencies_ms2'))
```

```
Traceback (most recent call last):
  File "/home/docs/checkouts/readthedocs.org/user_builds/gut-to-soil-tutorial/checkouts/latest/.env/lib/python3.10/site-packages/q2doc/transforms/transform_usage.py", line 92, in run
    exec(source, interface['driver'].scope)
  File "<string>", line 4, in <module>
NameError: name 'asv_table_ms2' is not defined
```

Be sure to run this as we’re going to use one of the results below.

### Taxonomic annotation(#UIi09DapR3 "Link to this Section")

Before we complete our upstream analysis steps, we’ll generate taxonomic annotations for our sequences using the [feature-classifier plugin](https://amplicon-docs.qiime2.org/en/latest/references/plugins/feature-classifier#q2-plugin-feature-classifier).

First, we’ll download a pre-trained classifier artifact.

[Command Line]

[Python API]

[Galaxy]

[R API]

[View Source]

```
wget -O 'suboptimal-16S-rRNA-classifier.qza' \
  'https://gut-to-soil-tutorial.readthedocs.io/en/latest/data/gut-to-soil/suboptimal-16S-rRNA-classifier.qza'
```

```
import qiime2.plugins.dada2.actions as dada2_actions

url = 'https://gut-to-soil-tutorial.readthedocs.io/en/latest/data/gut-to-soil/suboptimal-16S-rRNA-classifier.qza'
fn = 'suboptimal-16S-rRNA-classifier.qza'
request.urlretrieve(url, fn)
suboptimal_16S_rRNA_classifier = Artifact.load(fn)
```

Using the `Upload Data` tool:
:   1. On the first tab (**Regular**), press the `Paste/Fetch` data button at the bottom.
       1. Set *"Name"* (first text-field) to: `suboptimal-16S-rRNA-classifier.qza`
       2. In the larger text-area, copy-and-paste: <https://gut-to-soil-tutorial.readthedocs.io/en/latest/data/gut-to-soil/suboptimal-16S-rRNA-classifier.qza>
       3. (*"Type"*, *"Genome"*, and *"Settings"* can be ignored)
    2. Press the `Start` button at the bottom.

```
url <- 'https://gut-to-soil-tutorial.readthedocs.io/en/latest/data/gut-to-soil/suboptimal-16S-rRNA-classifier.qza'
fn <- 'suboptimal-16S-rRNA-classifier.qza'
request$urlretrieve(url, fn)
suboptimal_16S_rRNA_classifier <- Artifact$load(fn)
```

```
def classifier_factory():
    from urllib import request
    from qiime2 import Artifact
    fp, _ = request.urlretrieve(
        'https://data.qiime2.org/classifiers/sklearn-1.4.2/greengenes/gg-13-8-99-515-806-nb-classifier.qza')

    return Artifact.load(fp)

classifier = use.init_artifact('suboptimal-16S-rRNA-classifier', classifier_factory)
```

* `suboptimal-16S-rRNA-classifier.qza` | [download](https://gut-to-soil-tutorial.readthedocs.io/en/latest/data/gut-to-soil/suboptimal-16S-rRNA-classifier.qza) | [view](https://view.qiime2.org?src=https%3A%2F%2Fgut-to-soil-tutorial.readthedocs.io%2Fen%2Flatest%2Fdata%2Fgut-to-soil%2Fsuboptimal-16S-rRNA-classifier.qza)

Then, we’ll apply it to our sequences using [`classify-sklearn`](https://amplicon-docs.qiime2.org/en/latest/references/plugins/feature-classifier#q2-action-feature-classifier-classify-sklearn).

[Command Line]

[View Source]

```
qiime feature-classifier classify-sklearn \
  --i-classifier suboptimal-16S-rRNA-classifier.qza \
  --i-reads asv-seqs-ms2.qza \
  --o-classification taxonomy.qza
```

```
taxonomy, = use.action(
    use.UsageAction(plugin_id='feature_classifier',
                    action_id='classify_sklearn'),
    use.UsageInputs(classifier=classifier,
                    reads=asv_seqs_ms2),
    use.UsageOutputNames(classification='taxonomy'))
```

```
Traceback (most recent call last):
  File "/home/docs/checkouts/readthedocs.org/user_builds/gut-to-soil-tutorial/checkouts/latest/.env/lib/python3.10/site-packages/q2doc/transforms/transform_usage.py", line 92, in run
    exec(source, interface['driver'].scope)
  File "<string>", line 5, in <module>
NameError: name 'asv_seqs_ms2' is not defined
```

Then, to get an initial look at our taxonomic classifications, let’s integrate taxonomy in the sequence summary, like the one we generated above.

[Command Line]

[View Source]

```
qiime feature-table tabulate-seqs \
  --i-data asv-seqs-ms2.qza \
  --i-taxonomy Greengenes_13_8:taxonomy.qza \
  --m-metadata-file asv-frequencies-ms2.qza \
  --o-visualization asv-seqs-ms2.qzv
```

```
asv_frequencies_ms2_as_md = use.view_as_metadata('asv_frequencies',
                                                 asv_frequencies_ms2)

taxonomy_collection = use.construct_artifact_collection(
    'taxonomy_collection', {'Greengenes_13_8': taxonomy}
)

use.action(
    use.UsageAction(plugin_id='feature_table',
                    action_id='tabulate_seqs'),
    use.UsageInputs(data=asv_seqs_ms2,
                    taxonomy=taxonomy_collection,
                    metadata=asv_frequencies_ms2_as_md),
    use.UsageOutputNames(visualization='asv_seqs_ms2'),
)
```

```
Traceback (most recent call last):
  File "/home/docs/checkouts/readthedocs.org/user_builds/gut-to-soil-tutorial/checkouts/latest/.env/lib/python3.10/site-packages/q2doc/transforms/transform_usage.py", line 92, in run
    exec(source, interface['driver'].scope)
  File "<string>", line 2, in <module>
NameError: name 'asv_frequencies_ms2' is not defined
```

### Building a tree for phylogenetic diversity calculations(#JCBZPhaSQg "Link to this Section")

QIIME supports several phylogenetic diversity metrics, including Faith’s Phylogenetic Diversity and weighted and unweighted UniFrac.
In addition to counts of features per sample (i.e., the data in the `FeatureTable[Frequency]` QIIME 2 artifact), these metrics require a rooted phylogenetic tree relating the features to one another.
The amplicon distribution offers a few ways to build these trees, including a reference-based approach in the [fragment-insertion plugin](https://amplicon-docs.qiime2.org/en/latest/references/plugins/fragment-insertion#q2-plugin-fragment-insertion) and *de novo* (i.e., reference-free) approaches in the [phylogeny plugin](https://amplicon-docs.qiime2.org/en/latest/references/plugins/phylogeny#q2-plugin-phylogeny).

The reference based approach, by default, is specific to 16S rRNA marker gene analysis.
We could use that here, but the runtime is too long for our documentation.[[1]](#fn-build-requirements-exceed-resources "Link to Footnote")
If you’d like to see this demonstrated, you can refer to the [*Parkinson’s Mouse* tutorial](https://docs.qiime2.org/2024.10/tutorials/pd-mice/).

The *de novo* approach is known to generate low quality trees when very short sequences are used as input, but can be used with any phylogenetically informative marker gene (not just 16S).
If you’d like to see this demonstrated, you can refer to the [*Moving Pictures* tutorial](https://amplicon-docs.readthedocs.io/en/latest/tutorials/moving-pictures.html#generate-a-tree-for-phylogenetic-diversity-analyses).

For those reasons, we’re going to skip building phylogenetic trees and instead use an analog of phylogenetic diversity metrics here.

## Downstream data analysis(#gEfdm8ocZo "Link to this Section")

As mentioned above, we tend to think of “downstream” analysis as beginning with a feature table, taxonomic annotation of our features, and optionally a phylogenetic tree.
Now that we have those (with the exception of the tree, [which we won’t use here](https://gut-to-soil-tutorial.readthedocs.io/en/latest/gut-to-soil-16s-tutorial#phylogenetic-tree-building)), let’s jump in.
This is where it starts to get fun! ⛷️

### Kmer-based diversity analysis(#NRj2muy9oD "Link to this Section")

As mentioned above, we’re going to skip building phylogenetic trees and instead use an analog of phylogenetic diversity metrics here.
This will use two [stand-alone QIIME 2 plugins](https://amplicon-docs.qiime2.org/en/latest/back-matter/glossary#term-stand-alone-plugin), [q2-boots](https://library.qiime2.org/plugins/caporaso-lab/q2-boots) and [q2-kmerizer](https://library.qiime2.org/plugins/bokulich-lab/q2-kmerizer), which are integrated through the [`kmer-diversity`](https://q2-boots.readthedocs.io/en/latest/plugin-reference/plugins/boots#q2-action-boots-kmer-diversity) action in q2-boots.
To learn about kmerization of features and how this relates to phylogenetic diversity metrics, read the q2-kmerizer paper.
[q2-boots](https://library.qiime2.org/plugins/caporaso-lab/q2-boots) provides actions that mirror the interface of the diversity metric calculation actions in the diversity plugin, but generates more robust results because it integrates rarefaction and/or bootstrapping.
You can learn more about this in the q2-boots paper.

[`kmer-diversity`](https://q2-boots.readthedocs.io/en/latest/plugin-reference/plugins/boots#q2-action-boots-kmer-diversity) is a [`qiime2.Pipeline`](https://amplicon-docs.qiime2.org/en/latest/back-matter/glossary): a type of action that links multiple other QIIME 2 actions together for convenience.
As such, it does a lot of work.
Here are the steps that it takes:

1. Resample the input feature table (i.e., the ASV table) to contain exactly `sampling-depth` sequences per sample, either by bootstrapping or rarefaction, `n` times.
   Samples with fewer than `sampling-depth` sequences will be removed from the feature table and not included in the subsequent steps.
   This will result in `n` feature tables.
2. For each feature table resulting from step 1, using the input sequences (i.e., the ASV sequences), kmerize all sequences into kmers of length `kmer-size`. [[2]](#fn-iab-database-searching "Link to Footnote")
   Use this information to create one kmer table per resampled feature table.
   This will result in `n` feature tables, where the features are kmers (instead of ASVs, as in the input feature table).
3. Compute the user-requested alpha- and beta-diversity metrics on each of the kmer tables resulting from step 2. [[3]](#fn-forum-diversity-metrics "Link to Footnote")
   The metrics computed by default are:
   * Alpha diversity
     + Shannon’s diversity index (a quantitative measure of community richness)
     + Observed Features (a qualitative measure of community richness)
     + Evenness (i.e., Pielou’s Evenness; a measure of community evenness)
   * Beta diversity
     + Jaccard distance (a qualitative measure of community dissimilarity)
     + Bray-Curtis distance (a quantitative measure of community dissimilarity)
4. For each diversity metric, average the results computed across the `n` kmer tables.
   These results can be used in subsequent analysis steps (e.g., ordination, statistical modeling, machine learning).
5. Perform PCoA ordination on the averaged beta diversity distance matrices resulting from Step 4. [[4]](#fn-iab-machine-learning "Link to Footnote")
6. Generate an interactive [q2-vizard scatter plot](https://amplicon-docs.qiime2.org/en/latest/references/plugins/vizard#q2-action-vizard-scatterplot-2d) that contains all user-provided sample metadata, all averaged alpha diversity metrics, and the first three ordination axes for all PCoA matrices computed in step 5.

A key parameter that needs to be provided to [`kmer-diversity`](https://q2-boots.readthedocs.io/en/latest/plugin-reference/plugins/boots#q2-action-boots-kmer-diversity) is `sampling-depth`, which is the even sampling (i.e., bootstrapping or rarefaction) depth.
Because most diversity metrics are sensitive to different sampling depths (i.e., sequence counts) across different samples, the tables are randomly subsampled such that the total frequency for each sample is the user-specified sampling depth.
For example, if you set `sampling-depth=500`, this step will subsample the sequence counts in each sample so that each sample in the resulting table has a total frequency of 500.
If the total frequency (i.e., the number of sequences observed) for any sample(s) is smaller than this value, those samples will be dropped from the diversity analysis.
Choosing this value is tricky.
We recommend making your choice by reviewing the information presented in the `asv-table-ms2.qzv` file that [you created above](https://gut-to-soil-tutorial.readthedocs.io/en/latest/gut-to-soil-16s-tutorial#summarize-asv-table-ms2).

I’m going to choose values that is around the first quartile of the sample total frequencies.

[Command Line]

[View Source]

```
qiime boots kmer-diversity \
  --i-table asv-table-ms2.qza \
  --i-sequences asv-seqs-ms2.qza \
  --m-metadata-file sample-metadata.tsv \
  --p-sampling-depth 96 \
  --p-n 10 \
  --p-replacement \
  --p-alpha-average-method median \
  --p-beta-average-method medoid \
  --output-dir boots-kmer-diversity
```

```
use.action(
    use.UsageAction(plugin_id='boots',
                    action_id='kmer_diversity'),
    use.UsageInputs(table=asv_table_ms2,
                    sequences=asv_seqs_ms2,
                    metadata=sample_metadata,
                    sampling_depth=96,
                    n=10,
                    replacement=True,
                    alpha_average_method='median',
                    beta_average_method='medoid'),
    use.UsageOutputNames(
        resampled_tables='bootstrap_tables',
        kmer_tables='kmer_tables',
        alpha_diversities='bootstrap_alpha_diversities',
        distance_matrices='bootstrap_distance_matrices',
        pcoas='bootstrap_pcoas',
        scatter_plot='kmer_diversity_scatter_plot')
)
```

```
Traceback (most recent call last):
  File "/home/docs/checkouts/readthedocs.org/user_builds/gut-to-soil-tutorial/checkouts/latest/.env/lib/python3.10/site-packages/q2doc/transforms/transform_usage.py", line 92, in run
    exec(source, interface['driver'].scope)
  File "<string>", line 4, in <module>
NameError: name 'asv_table_ms2' is not defined
```

After computing diversity metrics, we can begin to explore the microbial composition of the samples in the context of the sample metadata.
You can review the sample metadata using one of the tabulated views of this file that [we created above](https://gut-to-soil-tutorial.readthedocs.io/en/latest/gut-to-soil-16s-tutorial#sample-metadata-tabulate-viz).

### Alpha rarefaction plotting(#nvW5Kgmqfp "Link to this Section")

In this section we’ll explore alpha diversity as a function of sampling depth using the[`alpha-rarefaction` action](https://amplicon-docs.qiime2.org/en/latest/references/plugins/diversity#q2-action-diversity-alpha-rarefaction).
This visualizer computes one or more alpha diversity metrics at multiple sampling depths, in steps between 1 (optionally controlled with `min-depth`) and the value provided as `max-depth`.
At each sampling depth step, 10 rarefied tables will be generated, and the diversity metrics will be computed for all samples in the tables.
The number of iterations (rarefied tables computed at each sampling depth) can be controlled with the `iterations` parameter.
Average diversity values will be plotted for each sample at each even sampling depth, and samples can be grouped based on metadata in the resulting visualization if sample metadata is provided.

The value that you provide for `max-depth` should be determined by reviewing the “Frequency per sample” information presented in the `asv-table-ms2.qzv` file.
In general, choosing a value that is somewhere around the median frequency seems to work well, but you may want to increase that value if the lines in the resulting rarefaction plot don’t appear to be leveling out, or decrease that value if you seem to be losing many of your samples due to low total frequencies closer to the minimum sampling depth than the maximum sampling depth.

[Command Line]

[View Source]

```
qiime diversity alpha-rarefaction \
  --i-table asv-table-ms2.qza \
  --p-max-depth 260 \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization alpha-rarefaction.qzv
```

```
use.action(
    use.UsageAction(plugin_id='diversity',
                    action_id='alpha_rarefaction'),
    use.UsageInputs(table=asv_table_ms2,
                    max_depth=260,
                    metadata=sample_metadata),
    use.UsageOutputNames(visualization='alpha_rarefaction'))
```

```
Traceback (most recent call last):
  File "/home/docs/checkouts/readthedocs.org/user_builds/gut-to-soil-tutorial/checkouts/latest/.env/lib/python3.10/site-packages/q2doc/transforms/transform_usage.py", line 92, in run
    exec(source, interface['driver'].scope)
  File "<string>", line 4, in <module>
NameError: name 'asv_table_ms2' is not defined
```

The visualization will have two plots.
The top plot is an alpha rarefaction plot, and is primarily used to determine if the richness of the samples has been fully observed or sequenced.
If the lines in the plot appear to “level out” (i.e., approach a slope of zero) at some sampling depth along the x-axis, that suggests that collecting additional sequences beyond that sampling depth would not be likely to result in the observation of additional features.
If the lines in the plot don’t level out, this may be because the richness of the samples hasn’t been fully observed yet (because too few sequences were collected), or it could be an indicator that a lot of sequencing error remains in the data (which is being mistaken for novel diversity).

The bottom plot in this visualization is important when grouping samples by metadata.
It illustrates the number of samples that remain in each group when the feature table is rarefied to each sampling depth.
If a given sampling depth `d` is larger than the total frequency of a sample `s` (i.e., the number of sequences that were obtained for sample `s`), it is not possible to compute the diversity metric for sample `s` at sampling depth `d`.
If many of the samples in a group have lower total frequencies than `d`, the average diversity presented for that group at `d` in the top plot will be unreliable because it will have been computed on relatively few samples.
When grouping samples by metadata, it is therefore essential to look at the bottom plot to ensure that the data presented in the top plot is reliable.

### Taxonomic analysis(#BKAsL1UzW7 "Link to this Section")

Next, we can view the taxonomic composition of our samples with interactive bar plots.
Generate those plots with the following command and then open the visualization.

[Command Line]

[View Source]

```
qiime taxa barplot \
  --i-table asv-table-ms2.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization taxa-bar-plots.qzv
```

```
use.action(
    use.UsageAction(plugin_id='taxa',
                    action_id='barplot'),
    use.UsageInputs(table=asv_table_ms2,
                    taxonomy=taxonomy,
                    metadata=sample_metadata),
    use.UsageOutputNames(visualization='taxa_bar_plots'))
```

```
Traceback (most recent call last):
  File "/home/docs/checkouts/readthedocs.org/user_builds/gut-to-soil-tutorial/checkouts/latest/.env/lib/python3.10/site-packages/q2doc/transforms/transform_usage.py", line 92, in run
    exec(source, interface['driver'].scope)
  File "<string>", line 4, in <module>
NameError: name 'asv_table_ms2' is not defined
```

Solution to [Exercise 13](https://gut-to-soil-tutorial.readthedocs.io/en/latest/gut-to-soil-16s-tutorial#why-asv-table)

We use the ASV table here because the feature ids in that table are the same as the ones used in the `FeatureData[Taxonomy]` artifact.
Additionally, kmerization of our data is a tool used for computing diversity metrics - not something we generally intend to use throughout our analyses.

### Differential abundance testing with ANCOM-BC2(#BBVbhnpd1X "Link to this Section")

ANCOM-BC2 is a compositionally-aware linear regression model that allows testing for differentially abundant features across sample groups while also implementing bias correction.
This can be accessed using the [`ancombc2` action](https://amplicon-docs.qiime2.org/en/latest/references/plugins/composition#q2-action-composition-ancombc2) in the [composition plugin](https://amplicon-docs.qiime2.org/en/latest/references/plugins/composition#q2-plugin-composition).

We’ll perform this analysis in a few steps.

#### Filter samples from the feature table(#hBXClDB2eV "Link to this Section")

First, we’ll filter our samples from our feature table such that we only have the three groups that we have the most samples for.

[Command Line]

[View Source]

```
qiime feature-table filter-samples \
  --i-table asv-table-ms2.qza \
  --m-metadata-file sample-metadata.tsv \
  --p-where '[SampleType] IN ("Human Excrement Compost", "Human Excrement", "Food Compost")' \
  --o-filtered-table asv-table-ms2-dominant-sample-types.qza
```

```
asv_table_ms2_dominant_sample_types, = use.action(
    use.UsageAction(plugin_id='feature_table',
                    action_id='filter_samples'),
    use.UsageInputs(table=asv_table_ms2,
                    metadata=sample_metadata,
                    where='[SampleType] IN ("Human Excrement Compost", "Human Excrement", "Food Compost")'),
    use.UsageOutputNames(filtered_table='asv_table_ms2_dominant_sample_types'))
```

```
Traceback (most recent call last):
  File "/home/docs/checkouts/readthedocs.org/user_builds/gut-to-soil-tutorial/checkouts/latest/.env/lib/python3.10/site-packages/q2doc/transforms/transform_usage.py", line 92, in run
    exec(source, interface['driver'].scope)
  File "<string>", line 4, in <module>
NameError: name 'asv_table_ms2' is not defined
```

#### Apply differential abundance testing(#EbMYAdk1px "Link to this Section")

Then, we’ll apply ANCOM-BC2 to see which ASV are differentially abundant across those sample types.
I specify a reference level here as this defines what each group is compared against.
Since the focus of this study is HEC, I choose that as my reference level.
That will let us see what ASVs are over- or under-represented in the other two sample groups (*Human Excrement* and *Food Compost*) relative to HEC, as HEC defines the “global intercept” that will be measured against.

[Command Line]

[View Source]

```
qiime composition ancombc2 \
  --i-table asv-table-ms2-dominant-sample-types.qza \
  --m-metadata-file sample-metadata.tsv \
  --p-fixed-effects-formula SampleType \
  --p-reference-levels 'SampleType::Human Excrement Compost' \
  --o-ancombc2-output ancombc2-results.qza
```

```
ancombc2_results, = use.action(
    use.UsageAction(plugin_id='composition',
                    action_id='ancombc2'),
    use.UsageInputs(table=asv_table_ms2_dominant_sample_types,
                    metadata=sample_metadata,
                    fixed_effects_formula='SampleType',
                    reference_levels=['SampleType::Human Excrement Compost']),
    use.UsageOutputNames(ancombc2_output='ancombc2_results'))
```

```
Traceback (most recent call last):
  File "/home/docs/checkouts/readthedocs.org/user_builds/gut-to-soil-tutorial/checkouts/latest/.env/lib/python3.10/site-packages/q2doc/transforms/transform_usage.py", line 92, in run
    exec(source, interface['driver'].scope)
  File "<string>", line 4, in <module>
NameError: name 'asv_table_ms2_dominant_sample_types' is not defined
```

Finally, we’ll visualize the results.
Taxonomic annotations can optionally be provided here: this helps make ASVs ids more interpretable.

[Command Line]

[View Source]

```
qiime composition ancombc2-visualizer \
  --i-data ancombc2-results.qza \
  --i-taxonomy taxonomy.qza \
  --o-visualization ancombc2-barplot.qzv
```

```
use.action(
    use.UsageAction(plugin_id='composition',
                    action_id='ancombc2_visualizer'),
    use.UsageInputs(data=ancombc2_results,
                    taxonomy=taxonomy),
    use.UsageOutputNames(visualization='ancombc2-barplot'))
```

```
Traceback (most recent call last):
  File "/home/docs/checkouts/readthedocs.org/user_builds/gut-to-soil-tutorial/checkouts/latest/.env/lib/python3.10/site-packages/q2doc/transforms/transform_usage.py", line 92, in run
    exec(source, interface['driver'].scope)
  File "<string>", line 4, in <module>
NameError: name 'ancombc2_results' is not defined
```

Solution to [Exercise 16](https://gut-to-soil-tutorial.readthedocs.io/en/latest/gut-to-soil-16s-tutorial#ancombc2-genera)

To collapse our ASVs into genera (i.e. level 6 of the Greengenes taxonomy), we can use the following command.

[Command Line]

[View Source]

```
qiime taxa collapse \
  --i-table asv-table-ms2-dominant-sample-types.qza \
  --i-taxonomy taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table genus-table-ms2-dominant-sample-types.qza
```

```
genus_table_ms2_dominant_sample_types, = use.action(
    use.UsageAction(plugin_id='taxa',
                    action_id='collapse'),
    use.UsageInputs(table=asv_table_ms2_dominant_sample_types,
                    taxonomy=taxonomy,
                    level=6),
    use.UsageOutputNames(collapsed_table='genus_table_ms2_dominant_sample_types'))
```

```
Traceback (most recent call last):
  File "/home/docs/checkouts/readthedocs.org/user_builds/gut-to-soil-tutorial/checkouts/latest/.env/lib/python3.10/site-packages/q2doc/transforms/transform_usage.py", line 92, in run
    exec(source, interface['driver'].scope)
  File "<string>", line 4, in <module>
NameError: name 'asv_table_ms2_dominant_sample_types' is not defined
```

We can then provide the resulting table as the input to ANCOM-BC2.

[Command Line]

[View Source]

```
qiime composition ancombc2 \
  --i-table genus-table-ms2-dominant-sample-types.qza \
  --m-metadata-file sample-metadata.tsv \
  --p-fixed-effects-formula SampleType \
  --p-reference-levels 'SampleType::Human Excrement Compost' \
  --o-ancombc2-output genus-ancombc2-results.qza
```

```
genus_ancombc2_results, = use.action(
    use.UsageAction(plugin_id='composition',
                    action_id='ancombc2'),
    use.UsageInputs(table=genus_table_ms2_dominant_sample_types,
                    metadata=sample_metadata,
                    fixed_effects_formula='SampleType',
                    reference_levels=['SampleType::Human Excrement Compost']),
    use.UsageOutputNames(ancombc2_output='genus_ancombc2_results'))
```

```
Traceback (most recent call last):
  File "/home/docs/checkouts/readthedocs.org/user_builds/gut-to-soil-tutorial/checkouts/latest/.env/lib/python3.10/site-packages/q2doc/transforms/transform_usage.py", line 92, in run
    exec(source, interface['driver'].scope)
  File "<string>", line 4, in <module>
NameError: name 'genus_table_ms2_dominant_sample_types' is not defined
```

And finally, we can visualize the results.
Notice that in this case we’re not providing the taxonomy, because we’ve already intergrated that information by collapsing at the genus level.

[Command Line]

[View Source]

```
qiime composition ancombc2-visualizer \
  --i-data genus-ancombc2-results.qza \
  --o-visualization genus-ancombc2-barplot.qzv
```

```
use.action(
    use.UsageAction(plugin_id='composition',
                    action_id='ancombc2_visualizer'),
    use.UsageInputs(data=genus_ancombc2_results),
    use.UsageOutputNames(visualization='genus-ancombc2-barplot'))
```

```
Traceback (most recent call last):
  File "/home/docs/checkouts/readthedocs.org/user_builds/gut-to-soil-tutorial/checkouts/latest/.env/lib/python3.10/site-packages/q2doc/transforms/transform_usage.py", line 92, in run
    exec(source, interface['driver'].scope)
  File "<string>", line 4, in <module>
NameError: name 'genus_ancombc2_results' is not defined
```

## That’s it for now, but more is coming soon!(#CvUsCK1jPE "Link to this Section")

In the near future (as of 17 April 2025) we plan to integrate analyses using the [sample-classifier plugin](https://amplicon-docs.qiime2.org/en/latest/references/plugins/sample-classifier#q2-plugin-sample-classifier) and [longitudinal plugin](https://amplicon-docs.qiime2.org/en/latest/references/plugins/longitudinal#q2-plugin-longitudinal).
In the meantime, here are some suggestions to continue your learning:

1. Build a machine learning classifier that classifies samples accordining to the three dominant sample types in the feature table that we used with ANCOM-BC.
   (Hint: see [`classify-samples`](https://amplicon-docs.qiime2.org/en/latest/references/plugins/sample-classifier#q2-action-sample-classifier-classify-samples).)
2. Perform a longitudinal analysis that tracks samples from different buckets over time. Which taxa change most over time? (Hint: see [`feature-volatility`](https://amplicon-docs.qiime2.org/en/latest/references/plugins/longitudinal#q2-action-longitudinal-feature-volatility).)
3. Remember that the full data set (five sequencing runs) are available in the gut-to-soil Artifact Repository.
   Grab one of the larger sequencing runs (we worked with a small sequencing run that was generated as a preliminary test), and adapt the commands in this tutorial to work on a bigger data set.

We’re also in the process of refactoring our statistical methods for assessing alpha and beta diversity across groups, using the new [stats plugin](https://amplicon-docs.qiime2.org/en/latest/references/plugins/stats#q2-plugin-stats).
We’re therefore holding off on integrating statistical analysis until we have that ready.
In the meantime, you can refer to you can refer to the [*Moving Pictures*](https://amplicon-docs.readthedocs.io/en/latest/tutorials/moving-pictures.html) tutorial, as well as the [sample-classifier](https://docs.qiime2.org/2024.10/tutorials/sample-classifier/) and [longitudinal](https://docs.qiime2.org/2024.10/tutorials/longitudinal/) tutorials.

## Replay provenance (work in progress!)(#vkatrBHcFX "Link to this Section")

You might next want to try to adapt the commands presented in this tutorial to your own data, adjusting parameter settings and metadata column headers as is relevant.
QIIME 2’s provenance replay functionality can help with this.
Assuming that you ran all of the steps above in a directory called `gut-to-soil/`, run the following command to generate a template script that you can adapt for your workflow:

```
qiime tools replay-provenance \
  --in-fp gut-to-soil/taxa-bar-plots.qzv \
  --out-fp g2s-replayed.bash
```

If you need help, head over to the [QIIME 2 Forum](https://forum.qiime2.org).

Footnotes(#footnotes "Link to Footnotes")

1. The resource requirements exceed those provided by the [*Read the Docs* (RTD) build system](https://docs.readthedocs.com/platform/stable/builds.html#build-resources), which is used to build the documentation that you’re reading.
   RTD provides systems with 7GB of RAM for 30 minutes maximum to build documentation.
   That’s a very reasonable (and generous) allocation for building documentation, so we choose to work within those contraints rather than creating our own documentation build system like we’ve had in the past (e.g., for `https://docs.qiime2.org`).

   [](#fnref-Ry6x5KyZj8 "Link to Content")[](#fnref-ryOTr7NNFK "Link to Content")
2. kmerization of biological sequences is described in the [*Database Searching* chapter of *An Introduction to Applied Bioinformatics*](https://readiab.org/database-searching.html#kmer-content).

   [](#fnref-KfD19Y7i9G "Link to Content")
3. Learn more about the available metrics in [this QIIME 2 Forum post](https://forum.qiime2.org/t/alpha-and-beta-diversity-explanations-and-commands/2282).

   [](#fnref-EUusfiLe0S "Link to Content")
4. This process is discussed in the [*Machine Learning in Bioinformatics* chapter of *An Introduction to Applied Bioinformatics*](https://readiab.org/machine-learning.html#unsupervised-learning).

   [](#fnref-xkxuydolLx "Link to Content")

[Tutorials

Moving Pictures tutorial 🎥](/en/latest/tutorials/moving-pictures.html)[How To Guides

How to install QIIME 2](/en/latest/how-to-guides/install.html)
