# Metrics Reference - Loading Assay Quantification

Quick reference for the computed metrics in the loading analysis scripts.
Use this when deciding which metric to report in a paper or presentation.

## Percent Wildtype

- **Formula**: Raw value, computed upstream in Excel. WT = 100 per replicate.
- **Interpretation**: Direct measure of loading activity as a fraction of
  wildtype. A value of 25 means 25% of WT loading.
- **When to use**: Primary metric for figures and results text. Standard for
  loading assays. Report this on the y-axis of bar charts.
- **Example**: "+4sofa loaded at 23.1 ñ 10.9% of wildtype levels."

## Percent Change from WT

- **Formula**: `(condition - WT) / WT x 100`
- **Interpretation**: Directional. Negative values mean reduced loading
  relative to WT. Mathematically equivalent to `Percent Wildtype - 100`.
- **When to use**: When framing results as loss-of-function. Useful in text
  when emphasizing how much loading was lost.
- **Example**: "ORC4R reduced loading by 94.2 ñ 6.6% relative to wildtype."
- **Caveat**: Redundant with Percent Wildtype (just shifted by 100). Choose
  one framing or the other, not both.

## Percent Difference from WT / from ORC4R

- **Formula**: `|A - B| / ((A + B) / 2) x 100`
- **Interpretation**: Symmetric. Neither value is treated as the reference.
  Ranges from 0% (identical) to 200% (one value is zero).
- **When to use**: When comparing two conditions on equal footing, such as
  two mutants against each other where neither is a defined standard.
  Appropriate in discussion sections for pairwise mutant comparisons.
- **Example**: "The percent difference between +4sofa and +1sofa was X%,
  indicating comparable rescue activity."
- **Caveat**: Less common in loading assay literature. Reviewers may not
  immediately recognize the metric. Define it explicitly if used.

## Percent Change from ORC4R

- **Formula**: `(condition - ORC4R) / ORC4R x 100`
- **Interpretation**: Directional. Positive values mean the condition loads
  more than ORC4R alone.
- **When to use**: Rarely recommended for this dataset.
- **Caveat**: ORC4R values are very small (mean ~5.8% of WT) and variable
  across replicates (1.3 to 13.4). Dividing by these small numbers produces
  inflated, unstable values (e.g., WT shows 3936 ñ 3419%). Use fold change
  instead.

## Fold Change from ORC4R

- **Formula**: `condition / ORC4R` (paired within replicate)
- **Interpretation**: How many times more loading the condition achieves
  compared to ORC4R alone. A value of 1 means identical to ORC4R.
- **When to use**: When the narrative is about suppressor rescue relative to
  the mutant baseline. Standard phrasing in suppressor genetics.
- **Example**: "+4sofa restored loading approximately 10-fold relative to
  ORC4R alone."
- **Caveat**: Same small-denominator issue as percent change from ORC4R, but
  the values are more interpretable (single-digit fold changes vs thousands
  of percent). Still expect high SDs. Consider reporting median instead of
  mean, or note the variability explicitly. With n=3, individual replicate
  values may be more informative than summary statistics.

## Choosing a Metric for Publication

| Narrative framing                        | Recommended metric          |
|------------------------------------------|-----------------------------|
| How active is each condition?            | Percent Wildtype            |
| How much loading was lost vs WT?         | Percent Change from WT      |
| How much did the suppressor rescue?      | Fold Change from ORC4R      |
| How different are two mutants?           | Percent Difference          |
| Avoid                                    | Percent Change from ORC4R   |

When in doubt, use Percent Wildtype for figures and fold change from ORC4R
for text describing suppressor rescue. These are the two most conventional
and interpretable metrics for this type of experiment.
