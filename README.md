ChromBPNet Variant-Effect Prediction Pipeline


<img width="1282" height="1079" alt="synthetic_shap_logo_1" src="https://github.com/user-attachments/assets/42f069a3-eea8-4e5e-b6c9-1338a48b8929" />
Synthetic Example Output Plot from Variant Effect Scorer


A SLURM-ready Bash pipeline that scores genetic variants with ChromBPNet models,
averages their effects across cross-validation folds, annotates them with peaks,
nearest genes, and motif hits, and optionally computes per-variant SHAP scores.
Built to run unattended on an HPC GPU node: every stage is idempotent, input-checked,
and individually toggleable, so a partial run can be resumed without recomputing
finished work.

Stages (each toggleable):
1. Score    - run the ChromBPNet variant scorer once per model fold.
2. Average  - collapse per-fold scores into one mean score per variant, per condition.
3. Annotate - intersect averaged variants with peaks, nearest genes, and a union of
              motif hits.
4. SHAP     - (optional) per-fold counts and profile SHAP scores per variant.

Design notes:
- All inputs are validated before any GPU compute starts.
- Resumable: stages skip existing outputs unless OVERWRITE_OUTPUTS=true.
- Config-only edits: all site-specific values live in one CONFIGURATION block; paths
  expand from {cell_type}/{condition}/{fold} templates.


Acknowledgements:
Please cite ChromBPNet and the Kundaje Lab for original software creation if used in a published work. 
