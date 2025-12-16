Source URL: https://biocpy.github.io/tutorial/chapters/workflow.html
Date Scraped: 2025-12-15

---

In this section, we will illustrate a workflow that utilizes either language-agnostic representations for storing genomic data or reading RDS files directly in Python, to facilitate seamless access to datasets and analysis results.

Note

Check out

* the [interop with R](https://biocpy.github.io/tutorial/chapters/interop.html) section for reading RDS files directly in Python or
* the [language agnostic](https://biocpy.github.io/tutorial/chapters/language_agnostic.html) representations for storing genomic data

To begin, we will download the “zilionis lung” dataset from the [scRNAseq](https://bioconductor.org/packages/release/data/experiment/html/scRNAseq.html) package. Subsequently, we will store this dataset as an RDS file.

```
library(scRNAseq)

sce <- ZeiselBrainData()
sub <- sce[,1:2000]
saveRDS(sub, "../assets/data/zilinois-lung-subset.rds")
```

To demonstrate this workflow, we will use the models from the [CellTypist](https://github.com/Teichlab/celltypist) package to annotate cell types for this dataset. CellTypist operates on an `AnnData` representation.

```
from rds2py import read_rds, as_summarized_experiment
import numpy as np

r_object = read_rds("../assets/data/zilinois-lung-subset.rds")
sce = as_summarized_experiment(r_object)
adata, _ = sce.to_anndata()
adata.X = np.log1p(adata.layers["counts"])
adata.var.index = adata.var["genes"].tolist()
print(adata)
```

```
/opt/hostedtoolcache/Python/3.9.18/x64/lib/python3.9/site-packages/genomicranges/SeqInfo.py:348: UserWarning: 'seqnames' is deprecated, use 'get_seqnames' instead
  warn("'seqnames' is deprecated, use 'get_seqnames' instead", UserWarning)
```

```
AnnData object with n_obs × n_vars = 2000 × 41861
    obs: 'Library', 'Barcode', 'Patient', 'Tissue', 'Used', 'Barcoding emulsion', 'Total counts', 'Percent counts from mitochondrial genes', 'Most likely LM22 cell type', 'Major cell type', 'Minor subset', 'used_in_NSCLC_all_cells', 'used_in_NSCLC_and_blood_immune', 'used_in_NSCLC_immune', 'used_in_NSCLC_non_immune', 'used_in_blood'
    var: 'genes'
    layers: 'counts'
```

Before inferring cell types, let’s download the “human lung atlas” model from CellTypist.

```
import celltypist
from celltypist import models

models.download_models()
model_name = "Human_Lung_Atlas.pkl"
model = models.Model.load(model = model_name)
print(model)
```

```
📜 Retrieving model list from server https://celltypist.cog.sanger.ac.uk/models/models.json
📚 Total models in list: 44
📂 Storing models in /home/runner/.celltypist/data/models
💾 Downloading model [1/44]: Immune_All_Low.pkl
💾 Downloading model [2/44]: Immune_All_High.pkl
💾 Downloading model [3/44]: Adult_CynomolgusMacaque_Hippocampus.pkl
💾 Downloading model [4/44]: Adult_Human_PancreaticIslet.pkl
💾 Downloading model [5/44]: Adult_Human_Skin.pkl
💾 Downloading model [6/44]: Adult_Mouse_Gut.pkl
💾 Downloading model [7/44]: Adult_Mouse_OlfactoryBulb.pkl
💾 Downloading model [8/44]: Adult_Pig_Hippocampus.pkl
💾 Downloading model [9/44]: Adult_RhesusMacaque_Hippocampus.pkl
💾 Downloading model [10/44]: Autopsy_COVID19_Lung.pkl
💾 Downloading model [11/44]: COVID19_HumanChallenge_Blood.pkl
💾 Downloading model [12/44]: COVID19_Immune_Landscape.pkl
💾 Downloading model [13/44]: Cells_Fetal_Lung.pkl
💾 Downloading model [14/44]: Cells_Intestinal_Tract.pkl
💾 Downloading model [15/44]: Cells_Lung_Airway.pkl
💾 Downloading model [16/44]: Developing_Human_Brain.pkl
💾 Downloading model [17/44]: Developing_Human_Gonads.pkl
💾 Downloading model [18/44]: Developing_Human_Hippocampus.pkl
💾 Downloading model [19/44]: Developing_Human_Organs.pkl
💾 Downloading model [20/44]: Developing_Human_Thymus.pkl
💾 Downloading model [21/44]: Developing_Mouse_Brain.pkl
💾 Downloading model [22/44]: Developing_Mouse_Hippocampus.pkl
💾 Downloading model [23/44]: Fetal_Human_AdrenalGlands.pkl
💾 Downloading model [24/44]: Fetal_Human_Pancreas.pkl
💾 Downloading model [25/44]: Fetal_Human_Pituitary.pkl
💾 Downloading model [26/44]: Fetal_Human_Retina.pkl
💾 Downloading model [27/44]: Fetal_Human_Skin.pkl
💾 Downloading model [28/44]: Healthy_Adult_Heart.pkl
💾 Downloading model [29/44]: Healthy_COVID19_PBMC.pkl
💾 Downloading model [30/44]: Healthy_Human_Liver.pkl
💾 Downloading model [31/44]: Healthy_Mouse_Liver.pkl
💾 Downloading model [32/44]: Human_AdultAged_Hippocampus.pkl
💾 Downloading model [33/44]: Human_Developmental_Retina.pkl
💾 Downloading model [34/44]: Human_Embryonic_YolkSac.pkl
💾 Downloading model [35/44]: Human_IPF_Lung.pkl
💾 Downloading model [36/44]: Human_Longitudinal_Hippocampus.pkl
💾 Downloading model [37/44]: Human_Lung_Atlas.pkl
💾 Downloading model [38/44]: Human_PF_Lung.pkl
💾 Downloading model [39/44]: Lethal_COVID19_Lung.pkl
💾 Downloading model [40/44]: Mouse_Dentate_Gyrus.pkl
💾 Downloading model [41/44]: Mouse_Isocortex_Hippocampus.pkl
💾 Downloading model [42/44]: Mouse_Postnatal_DentateGyrus.pkl
💾 Downloading model [43/44]: Nuclei_Lung_Airway.pkl
💾 Downloading model [44/44]: Pan_Fetal_Human.pkl
```

```
CellTypist model with 61 cell types and 5017 features
    date: 2023-05-17 19:51:55.661237
    details: integrated Human Lung Cell Atlas (HLCA) combining multiple datasets of the healthy respiratory system
    source: https://doi.org/10.1101/2022.03.10.483747
    version: v2
    cell types: AT0, AT1, ..., pre-TB secretory
    features: TSPAN6, FGR, ..., RP1-34B20.21
```

Now, let’s annotate our dataset.

```
predictions = celltypist.annotate(adata, model = model_name, majority_voting = True)
print(predictions.predicted_labels)
```

```
⚠️ Warning: invalid expression matrix, expect ALL genes and log1p normalized expression to 10000 counts per cell. The prediction result may not be accurate
🔬 Input data has 2000 cells and 41861 genes
🔗 Matching reference genes in the model
🧬 4829 features used for prediction
⚖️ Scaling input data
🖋️ Predicting labels
✅ Prediction done!
👀 Can not detect a neighborhood graph, will construct one before the over-clustering
⛓️ Over-clustering input data with resolution set to 5
🗳️ Majority voting the predictions
✅ Majority voting done!
```

```
            predicted_labels over_clustering       majority_voting
bcIIOD        Neuroendocrine               4          Hillock-like
bcHTNA  Monocyte-derived Mph              28  Monocyte-derived Mph
bcDLAV  EC general capillary               4          Hillock-like
bcHNVA  Monocyte-derived Mph              28  Monocyte-derived Mph
bcALZN  Monocyte-derived Mph              26  Monocyte-derived Mph
...                      ...             ...                   ...
bcCXKF         Basal resting              16                   DC2
bcHCHL  EC general capillary               7           CD4 T cells
bcESWD         Basal resting              31         Basal resting
bcEKPJ           CD4 T cells              22         Basal resting
bcEFKS           CD4 T cells               2         Basal resting

[2000 rows x 3 columns]
```

Note

This CellTypist workflow is based on the tutorial described [here](https://colab.research.google.com/github/Teichlab/celltypist/blob/main/docs/notebook/celltypist_tutorial.ipynb#scrollTo=postal-chicken).

Next, let’s retrieve the `AnnData` object with the predicted labels embedded into the `obs` dataframe.

```
adata = predictions.to_adata()
adata
```

```
AnnData object with n_obs × n_vars = 2000 × 41861
    obs: 'Library', 'Barcode', 'Patient', 'Tissue', 'Used', 'Barcoding emulsion', 'Total counts', 'Percent counts from mitochondrial genes', 'Most likely LM22 cell type', 'Major cell type', 'Minor subset', 'used_in_NSCLC_all_cells', 'used_in_NSCLC_and_blood_immune', 'used_in_NSCLC_immune', 'used_in_NSCLC_non_immune', 'used_in_blood', 'predicted_labels', 'over_clustering', 'majority_voting', 'conf_score'
    var: 'genes'
    uns: 'neighbors', 'leiden'
    obsm: 'X_pca'
    layers: 'counts'
    obsp: 'connectivities', 'distances'
```

We can now reverse the workflow and save this object into an Artifactdb format from Python. However, the object needs to be converted to a `SingleCellExperiment` class first. Read more about our experiment representations [here](https://biocpy.github.io/tutorial/chapters/experiments/single_cell_experiment.html).

```
from singlecellexperiment import SingleCellExperiment

sce = SingleCellExperiment.from_anndata(adata)
print(sce)
```

```
class: SingleCellExperiment
dimensions: (41861, 2000)
assays(2): ['counts', 'X']
row_data columns(1): ['genes']
row_names(41861): ['5S_rRNA', '5_8S_rRNA', '7SK', ..., 'snoZ6', 'snosnR66', 'uc_338']
column_data columns(20): ['Library', 'Barcode', 'Patient', ..., 'over_clustering', 'majority_voting', 'conf_score']
column_names(2000): ['bcIIOD', 'bcHTNA', 'bcDLAV', ..., 'bcESWD', 'bcEKPJ', 'bcEFKS']
main_experiment_name:  
reduced_dims(0): []
alternative_experiments(0): []
row_pairs(0): []
column_pairs(0): []
metadata(2): neighbors leiden
```

We use the [dolomite suite](https://github.com/ArtifactDB/dolomite-base) to save it into a language-agnostic format from Python.

```
import dolomite_base
import dolomite_sce

dolomite_base.save_object(sce, "./zilinois_lung_with_celltypist")
```

Finally, read the object back in R.

```
sce_with_celltypist <- readObject(path=paste(getwd(), "zilinois_lung_with_celltypist", sep="/"))
sce_with_celltypist
```

And that concludes the workflow. Leveraging the generic **read** functions `readObject` (R) and `read_object` (Python), along with the **save** functions `saveObject` (R) and `save_object` (Python), you can seamlessly store most Bioconductor objects in language-agnostic formats.

---

* ArtifactDB GitHub organization - <https://github.com/ArtifactDB>.