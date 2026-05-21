def rounddecimals:
  .overdispersion=(.overdispersion*100.0 + 0.5 | floor/100.0) |
  .avg_per_cell=(.avg_per_cell*100.0 + 0.5 | floor/100.0) |
  .avg_per_gene=(.avg_per_gene*100.0 + 0.5 | floor/100.0) |
  .density=(.density*100.00 + 0.50 | floor/100.00)
;

(rounddecimals)
