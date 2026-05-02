# References: Runahead Loop/SC and Block-Level Loop/SC

Papers sourced from references sections of the 19 PDFs in `../papers/`.

---

## Runahead / Branch Precomputation

- S. Pruett, A. Zangeneh, A. Fakhrzadehgan, B. Lin, Y. Patt — **"Branch Runahead: An Alternative to Branch Prediction for Impossible to Predict Branches"** — MICRO-54, 2021
- A. Deshmukh, L. Cai, Y. Patt — **"Timely, Efficient, and Accurate Branch Precomputation"** — MICRO, 2024
- A. Deshmukh, Y. Patt — **"Criticality Driven Fetch"** — MICRO-54, 2021
- A. Nathani, J. Feliu, A. Adileh, L. Eeckhout — **"Precise Runahead Execution"** — HPCA, 2020
- A. Deshmukh, L. Cai, Y. Patt — **"Alternate Path Fetch"** — ISCA, 2024
- D. Schall, A. Sandberg, B. Grot — **"The Last-Level Branch Predictor"** — MICRO, 2024
- (Filtered Runahead — cited as [14] in ahead-prediction.pdf, full title/authors not resolved)

---

## Block-Level / Multiple Branches Per Cycle / Ahead Prediction

- A. Seznec, S. Jourdan, P. Sainrat, P. Michaud — **"Multiple-Block Ahead Branch Predictors"** — ASPLOS-VII, 1996 *(cited in BranchNet, EV8, AheadPred)*
- T-Y. Yeh, D. T. Marr, Y. N. Patt — **"Increasing the Instruction Fetch Rate via Multiple Branch Prediction and a Branch Address Cache"** — ICS '93, 1993
- G. Reinman, B. Calder, T. Austin — **"Fetch Directed Instruction Prefetching"** — MICRO-32, 1999
- A. Seznec, S. Felix, V. Krishnan, Y. Sazeides — **"Design Tradeoffs for the EV8 Branch Predictor"** — ISCA-29, 2002
- P. Michaud, A. Seznec, S. Jourdan — **"Exploration of Instruction Fetch Requirement in Out-of-Order Superscalar Processors"** — IJPP, 2001
- G. Giacalone, J. Edmonson — **"Method and Apparatus for Predicting Multiple Conditional Branches"** — US Patent 6,272,624, 2001
- J. L. Aragón, J. González, A. González, J. E. Smith — **"Dual Path Instruction Processing"** — ICS '02, 2002
- J. Pierce, T. Mudge — **"Wrong-Path Instruction Prefetching"** — MICRO-29, 1996
- A. Perais, R. Sheikh — **"Branch Target Buffer Organizations"** — MICRO, 2019

---

## Loop Predictor / Statistical Corrector

- A. Seznec, J. San Miguel, J. Albericio — **"The Innermost Loop Iteration Counter: A New Dimension in Branch History"** — MICRO, 2015 *(IMLI)*
- A. Seznec — **"A 256 Kbits L-TAGE Branch Predictor"** — CBP-2, 2007 *(introduces loop predictor + ahead pipelining)*
- A. Seznec — **"A New Case for the TAGE Branch Predictor"** — MICRO-44, 2011 *(introduces SC component)*
- A. Seznec — **"A 64 Kbits ISL-TAGE Branch Predictor"** — CBP-3, 2011
- A. Seznec, P. Michaud — **"A Case for (Partially)-Tagged Geometric History Length Predictors"** — JILP, 2006 *(original TAGE)*
- T. Heil, Z. Smith, J. E. Smith — **"Improving Branch Predictors by Correlating on Data Values"** — MICRO, 1999
