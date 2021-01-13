# TAI_V3H

This is the repository for the revised paper **V3H: View Variation and View Heredity for Incomplete Multi-view Clustering** resubmitted to **IEEE Transactions on Artificial Intelligence (TAI)**  by Xiang Fang, Yuchong Hu, Pan Zhou, and Dapeng Oliver Wu.

## View Variation and View Heredity

Multi-view clustering has wide applications in many real-world applications. In these applications, original image data often contain missing instances. In this repository, we implement a novel approach  View Variation and View Heredity (V3H) for incomplete multi-view clustering. 

We conduct extensive experiments on fifteen real-world datasets, and experimental results demonstrate V3H's superior advantages over other state-of-the-art clustering algorithms.
The codes of the compared methods can be found on the authors'  claimed websites.


## File directory

```bash
.
├── run_V3H_R1.m				                                               # DEMO file of V3H
├── V3H_R1.m				                                           # core function of V3H
├── YaleB.mat				                                             # data mat files
├── splitDigitData.m			                                       # construction of incomplete multi-view data
├── solveF.m				                                               # the initialization of F
├── NormalizeFea.m				                                       # regularization of data
├── ClusteringMeasure.m		                                       # clustering performance
└── constructW.m, EuDist2.m, L2_distance_1.m, and readsparse.m			 # intermediate functions 
```

## Usage

## Recommended operating environment

MATLAB R2020a, Windows 10, 3.30 GHz E3-1225 CPU, and 64 GB main memory.

### Download the V3H repository

0. Install the MATLAB. The scripts have been verified in Matlab 2020a.

1. Download this repository via git
    ```bash
    git clone https://github.com/ZeusDavide/TAI_V3H.git
    ```
    or download the [zip file](https://github.com/ZeusDavide/TAI_V3H/archive/master.zip) manually.
    
2. Get multi-view dataset: 
We provide the YaleB dataset "YaleB.mat" in this repository as an example. For the other datasets in the experiments, please refer to the corresponding links or articles.

3. Add the root folder to the Matlab path before running the scripts.

### Run V3H on incomplete multi-view data

To reproduce the experimental results in Section V-D of the paper, we need to run the scripts `run_V3H.m`.   


### Parameter tuning tips:

- For $\eta$, we set $\eta=10^{-3}$ (i.e., relatively small $\eta$) for $||\bm{M}||_{\eta}$ and $\tau=10^{-2}$ (i.e., relatively small $\tau$) for $||\bm{E}^{(v)}||_{\tau}$.
- In general, increasing iteration number `iter` will promote the clustering performance and consume more time. We recommend its maximum value is 30.



## Contact

[Xiang Fang, HUST](xfang9508@gmail.com)
