ChromBPNet Variant-Effect Prediction Pipeline



<img width="800" height="600" alt="synthetic_shap_logo_1" src="https://github.com/user-attachments/assets/0534a396-3978-452d-8c24-622f88e3af2f" />

*Example Output plot after running "VEP_Scorer.sh" and "Motifplot.R"

A SLURM-ready Bash pipeline that scores genetic variants with ChromBPNet models,
averages their effects across cross-validation folds, annotates them with peaks,
nearest genes, and motif hits, optionally computes per-variant SHAP scores, and
visualizes variant effects on TF motifs. Built to run unattended on an HPC GPU
node: every stage is idempotent, input-checked, and individually toggleable, so a
partial run can be resumed without recomputing finished work.

Stages (each toggleable):
1. Score     - run the ChromBPNet variant scorer once per model fold.
2. Average   - collapse per-fold scores into one mean score per variant, per condition.
3. Annotate  - intersect averaged variants with peaks, nearest genes, and a union of
               motif hits.
4. SHAP      - (optional) per-fold counts and profile SHAP scores per variant.
5. Visualize - plot SHAP contribution logos and PWM motif logos, and score variants
               for predicted gain/loss of TF motif binding.

Design notes:
- All inputs are validated before any GPU compute starts.
- Resumable: stages skip existing outputs unless OVERWRITE_OUTPUTS=true.
- Config-only edits: all site-specific values live in one CONFIGURATION block; paths
  expand from {cell_type}/{condition}/{fold} templates.


Acknowledgements:
Please cite ChromBPNet and the Kundaje Lab for original software creation if used in a
published work.
