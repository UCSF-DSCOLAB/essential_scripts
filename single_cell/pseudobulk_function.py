from ..python_utils.ts_log import ts_log
from scanpy.get import aggregate
import anndata as ad
from pandas.api.types import is_numeric_dtype

def dsco_pseudobulk(
    object, sample_by, cell_by,
    layer = None, use_raw = True,
    modality = "rna",
    min_cells = 10,
    cell_targets = None,
    # features = None,
    meta_targets = None,
    meta_num_summary_method = 'median',
    output_style = 'adata', # 'adata' or 'raw'
    output_metadata_cell_count = 'cells_in_pseudobulk',
    verbose = True):
    """
The primary role of this function is standardizing the pseudobulking process.
Note that it is under construction still and subject to change!.

Pseudobulking is performed per the groupings given by cell_by and any number of sample_by metadata using the scanpy.get.aggregate method to sums counts across cells.
Beforehand, for MuData 'object's, the targeted modality will be extracted.

Raw counts data is normally the desired data for pseudobulking.
If 'layer' is given, the matrix used for pseudobulking will be the 'object'.layers['layer'].
If no 'layer' is given and (default) 'use_raw' = True, the 'object'.raw.X slot element is used for pseudobulking.
Otherwise, 'object'.X is used in pseudobulking.

The number of cells that each pseudbulk represents are added to a metadata named 'cells_in_pseudobulk' by default, but this column name is adjustable with the 'output_metadata_cell_count' argument.

Pseudobulks representing fewer than 'min.cells' cells are removed.

'meta_targets', or all metadata of the target 'object' (or the targeted 'modality' for a muon 'object') are then extracted for each pseudobulk.
For discrete metadata, the column will be ignored if data are not consistent within ALL pseudobulks.
For numeric metadata, values from each cell in the pseudobulks are summarized by the 'meta_num_summary_method' method, ("median" by default).

Finally, the output structure is determined by the output_style argument, either ('adata') as an AnnData or ('raw') as a dict with keys ['counts', 'feature_names', 'metadata'] 

Prior to pseudobulking, particular cell identities or genes can be targeted using the \code{cell.targets} and \code{features} inputs, respectively.

Arguments:
    object                      An AnnData or MuData. For MuData objects, only elements inside the internal AnnData of the target 'modality' will remain accessible to the function.
    samply_by                   A string or list of strings naming columns of 'object'.obs to use for assigning sample identities of cells. A biospecimen, timepoint, or subject name is the common target here, and a batch identifier may be desired as well.
    cell_by                     String naming the metadata column of 'object'.obs to use for assigning annotation or cluster identities of cells.
    layer                       Optionally, a string naming as 'object'.layers key holding the count data.
    use_raw                     Boolean, True by default, denoting whether to use the object's .raw.X matrix instead of its .X matrix unless a 'layer' is given.
    modality                    String, "rna" by default, giving the modality name to use when 'object' is a MuData.
    min_cells                   Number, 10 by default, which sets the minimum number of cells that a pseudobulk should represent in order to be retained.
    cell_targets                Optionally, a string list naming particular 'cell_by' data values to target. Only these cell identities will be retained and pseudobulked.
    meta_targets                Optionally, a string list naming particular columns of 'object'.obs to target for retention. By default, as much a possible will be retained.
    meta_num_summary_method     String like "median" (default) or "mean" naming how to summarize numeric metadata describing how to summarize numeric metadata from all cells of each pseudobulk.
    output_style                'adata' (default) or 'raw'. Determines how data should be returned. Either as an AnnData object, or as a 'raw' dict of the 'counts', 'feature_names', and 'metadata'.
    output_metadata_cell_count  String ('cells_in_pseudobulk' by default) denoting the metadata column name added to hold the number of original cells going into each pseudobulk.
    verbose                     Logical which controls whether timestamped log messages should be shown

Example call:
    dsco_pseudobulk(
        object = adata,
        sample_by = ['biospecimen_id', 'orig.ident'],
        cell_by = 'leiden_res1',
        layer = 'counts'
    )

    """
    if verbose:
        msg_if = ts_log
    else:
        def msg_if(*args, **kwargs):
            return

    # To anndata with target data in 'layer' or .X
    if not isinstance(object, ad.AnnData):
        if modality in object.mod_names:
            object = object[modality]
        else:
            raise Exception("'object' is not an AnnData or not a MuData with 'modality' in its mod_names.")
    if layer is None and use_raw:
        object.X = object.raw.X

    # Ensure all of cell_by and sample_by exist as metadata
    orig_metas = object.obs.columns
    group_metas = list(sample_by) + [cell_by]
    if not all([i in orig_metas for i in group_metas]):
        raise Exception("A 'cell_by' or 'sample_by' element does not exist as a column of 'object.obs'.")

    # Trim to 'cell_targets'
    if not cell_targets is None:
        object = object[object.obs[cell_by] in cell_targets].copy()
        msg_if(
            "Trimming to 'cell_targets' retained ", object.n_obs,
            " cells of identities:\n", object.obs[cell_by].unique.join(", ")
        )

    ### Pseudobulk
    # (Trim to features internally)
    msg_if("Initiating pseudobulking with scanpy.get.aggregate...")
    psobject = aggregate(
        adata = object,
        by = group_metas,
        func = 'sum',
        layer = layer)
    msg_if("Moving counts to .X")
    psobject.X = psobject.layers['sum'].copy()
    del psobject.layers['sum']

    ### Add cell counts and Trim too small pseudobulks
    msg_if("Adding cell counts as metadata '", output_metadata_cell_count, "'.")
    counts = object.obs.value_counts(subset=group_metas)
    pseudo_vals = counts.index.copy()
    counts.index = ['_'.join(map(str, x)) for x in counts.index]
    psobject.obs[output_metadata_cell_count] = counts
    if min_cells > 0:
        too_small = sum(psobject.obs[output_metadata_cell_count] <= min_cells)
        if too_small > 0:
            msg_if(f"Trimming {too_small} pseudobulks built from fewer then {min_cells} cells.")
            psobject = psobject[psobject.obs[output_metadata_cell_count] >= min_cells].copy()

    ### Ensure meta_targets, or as many metadata as possible, are retained.
    meta_ignored = []
    meta_kept = psobject.obs.columns
    if meta_targets is None:
        meta_targets = orig_metas
    meta_targets = [i for i in meta_targets if not i in meta_kept]
    if len(meta_targets) > 0:
        msg_if(f"Grabbing additional metadata '{', '.join(meta_targets)}'.")
        numeric = [i for i in meta_targets if is_numeric_dtype(object.obs[i])]
        discrete = [i for i in meta_targets if not is_numeric_dtype(object.obs[i])]
        if len(discrete)>0:
            for col in discrete:
                if len(object.obs.value_counts(subset=group_metas+[col]).index) > len(pseudo_vals):
                    meta_ignored.append(col)
                    meta_targets.remove(col)
        if len(meta_ignored)>0:
            msg_if(f"\tignoring discrete metadata column(s) '{', '.join(meta_ignored)}' because of inconsistency within some pseudobulks")
        dtypes = object.obs.dtypes.astype(str).to_dict()
        for col in meta_targets:
            method = 'first' if not col in numeric else meta_num_summary_method
            new_dat = object.obs.groupby(group_metas, observed = True).agg(x=(col, method)).loc[pseudo_vals,:]
            new_dat.index = counts.index
            psobject.obs[col] = new_dat
            try:
                psobject.obs[col] = psobject.obs[col].astype(dtypes[col])
            except (KeyError, ValueError):
                pass

    # Output in requested style
    msg_if("Pseudobulk function COMPLETE.")
    if (output_style == 'adata'):
        return psobject
    return {
        'counts': psobject.X,
        'feature_names': psobject.var_names,
        'metadata': psobject.obs
    }
