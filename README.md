# repo_p123
Project p123 (SCD growth). Algorithm repository

Analysis Lead:  Christina Howitt
Contributor:    Ian Hambleton

Initial DO files provided by Ian Hambleton, from an initial analysis performed in (ahem...) 2004.

scdgrowth_001.do / scdgrowth_001.smcl
  (a) Load static (demographics) and height/weight information from the Jamaican cohort study.
  (b) Merge with WHO, US, and UK growth standards
  (c) Create individual-level z-scores based on WHO growth studies
  (d) Graphic. FOR MEN. Growth-chart of WHO standard growth, with SS and AA genotypes overlaid
  (e) Graphic. FOR WOMEN. Growth-chart of WHO standard growth, with SS and AA genotypes overlaid

  ToDo
  (a) There are some obvious Ht and Wt errors. We have "trimmed" the dataset for now - in an informal twoway
      See line 132 for this restriction, that we should discuss with Jamaica
