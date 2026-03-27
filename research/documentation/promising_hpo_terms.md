# Promising HPO terms 


## HP:0001022 (Albinism)
is the clear best overall: best F-score ~0.71 with balanced precision/recall (0.36 / 0.80) on 20 positives; 16/20 positives captured, 29 false positives.

## HP:0003387 (Decreased large myelinated fibers)
 is highly precise but very low recall: precision 0.83, recall 0.25; 5/20 positives captured, only 1 false positive. A good starting rule that needs recall boosts.


## HP:0003215 (Dicarboxylic aciduria) and HP:0004756 (Ventricular tachycardia) 
lean heavily on recall at the expense of precision (prec ~0.08–0.10, rec ~0.54–0.77), with large non-HPO hit counts (619 and 355 respectively). Rules are long (~3.8 terms).