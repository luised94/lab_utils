## Introduction
Data from next-generation-sequencing experiments is stored in a linux cluster. Code from lab_utils repository/next-generation-sequecing is meant to be used on that data. 

Some of the utils are designed to address the provenance of the data since I work with the MIT Biomicrocenter to generate the libraries, do the sequencing and some additional preprocessing. Furthermore, the samples are named according to their convention, which is YYMMDD<name_of_lab> for example 201001Bel for year, month, day and Bel for the Bell Lab. 

The following is a brief description of the data folders in my data directory.

## 221024Bel
Chromatin immunoprecipitation samples from W303 bar1 delta alpha factor and Nocodazole arrested cells with AID degradation module of ORC4. First troubleshooting experiment to analyze effect of suppressors on ORC and MCM binding. 

## 240324Bel
Second chromatin immunoprecipitation samples after some troubleshooting using CHIP-PCR. Still shows problems with MCM IP using polyclonal antibodies.

## 240630Bel

Chromatin immunoprecipitation after troubleshooting the MCM IPs further. In this case, the UM174 antibody was adjusted to 10ug (5 times more than reported in Miller 2023).

## EatonBel
Reference samples from Eaton 2010. Used for comparison to previous experiments.

## 250207Bel
Description: Chromatin immunoprecipitation after arresting cells in either alpha factor or nocodaloze. IP with HM1108 and UM174. Two repeats.
```{conclusions}
Preliminary: Very interesting results but unfortunately one of the repeats did not work due to alpha factor arrest problems
```

## 250324Bel
Received: 2025-04-29
Description:
  > Chromatin immunoprecipitation after arresting cells in either alpha factor or nocodaloze. IP with HM1108 and UM174. Two repeats.
  > Repeat to have at least two biological replicates
Related: 250207Bel
Jobs: 10105558
