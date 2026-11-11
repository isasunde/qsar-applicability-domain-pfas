| step                                |   records_in |   records_out |   removed | note                                 |
|:------------------------------------|-------------:|--------------:|----------:|:-------------------------------------|
| Raw export                          |         4200 |          4200 |         0 | ChEMBL activity export               |
| Type == Small molecule              |         4200 |          4195 |         5 |                                      |
| Standard Type == LogD7.4            |         4195 |          4195 |         0 |                                      |
| Numeric Standard Value              |         4195 |          4195 |         0 |                                      |
| RDKit parse + Standardize           |         4195 |          4195 |         0 | Unparsable or Failed standardization |
| Merge duplicate structures          |         4195 |          4190 |         5 | 5 structures had repeats             |
| Drop conflicting repeats (SD > 0.5) |         4190 |          4190 |         0 |                                      |
| Final curated set                   |         4190 |          4190 |         0 |                                      |