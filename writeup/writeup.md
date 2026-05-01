# Introduction

## Motivation

Frontier models write working code for real lab tasks, but they cost cents to dollars per call. Asking Opus to write a hundred-line bash script every time a postdoc aligns a few FASTQ files spends the budget that should go to harder problems. We ask a simple question: can a frontier model write the recipe once, and a free, small open-weight model running on the lab's own hardware turn that recipe into a working script every time after? If yes, the per-task cost drops from dollars to power and time. We test the split on one concrete task: calling mitochondrial variants from short-read sequencing data with BWA, samtools, LoFreq, and bcftools. We run it on four kinds of machine a real lab has — a workstation with one consumer GPU, a Jetson edge device, an Apple-silicon laptop, and a two-GPU workstation. The tools live in one pinned conda environment so the model can only use what we already have.
