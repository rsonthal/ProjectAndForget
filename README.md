# Code for Project and Forget paper: https://jmlr.org/papers/v23/20-1424.html

Please cite using 

@article{JMLR:v23:20-1424,
  author  = {Rishi Sonthalia and Anna C. Gilbert},
  title   = {Project and Forget: Solving Large-Scale Metric Constrained Problems},
  journal = {Journal of Machine Learning Research},
  year    = {2022},
  volume  = {23},
  number  = {326},
  pages   = {1--54},
  url     = {http://jmlr.org/papers/v23/20-1424.html}
}

Here we have the following files

Notebooks: All code is written in Julia 1.1.0. The include is normally present at the top of each notebook.

1) Metric Nearness.ipynb -- This has the code to run the metric nearness experiment. 

2) Correlation Clustering.ipynb --- This has the code to the run the Correlation Clustering experiment for both the dense and sparse case

3) L2 SVM --- This has the code to run the L2 SVM experiment

4) Plotting.ipynb -- This notebook plots all of the figures that we see in the paper. 

5) ParallelCC.jl, ParallelMetricOpt_1.0.jl, MetricOpt_helper.jl -- These run the code for the method we compare against for the CC problem. 

6) Information Theoretic Metric Learning.ipynb -- This contains the code for the information theoretic metric learning experiments. 

Inputs: 

1) Graphs --- This folder contains a variety of graphs used as input. 

2) Signed Graphs --- This contains the signed graphs

3) CA-GrQc.lgz, CA-HepTh.lgz, CA-HepPh.lgz ---- these contain some of the inputs in formats that we can directly use. (the bigger inputs can be found at the link in the paper)

4) ITML Data -- This folder contains all the data files for ITML experiments. 
