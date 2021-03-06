These results broadly validate the utility of mass-conserved flow-laws applied to reach-averaged data across a mass-conserved set of reaches. Errors introduced using simplified representations of physical processes lead to only minor degradation in performance metrics when replicating discharge estimates from distributed numerical models. Fully optimized, only a single case's estimates were beyond the peformance threshold of 35% rRMSE; no cases were outside the nRMSE threshold. 

The simplified flow-laws are variants of the long-established Manning's equation. However, their success is notable since Manning's equation explicitly describes relationships at a cross-section under steady-state conditions. This study's results support the extension of Manning's equation to dynamic flow conditions at the reach-scale (averaging many cross-sections) with only minor errors (<10% on average)

As reflected in the results, the mass-conservation constraint can actually improve flow-law estimates, even though it constitutes a simplification. This is because mass-conservation has two competing effects on flow-law estimates--one that creates errors by degrading the ability for estimates to reproduce differences across locations, and one that reduces errors by smoothing over flow-law errors present at individual reaches. In most of the validation datasets considered here, discharge is highly similar across reaches and so the second, error-reducing effect dominates. An exception to this (reflected in the rRMSE metric) is the Platte River, in which flood waves propagate quickly relative to the spacing of reaches along the river. 

As described in section ####(methods), this study validated three different parameter estimation methods--fully optimized, partially optimized, and physical/statistical parameter estimation. Which of these provides the most truthful representation of best performance? The answer depends on whether one is interested in physical or merely empirical truth. It is not inexcusable to use an empirical truth value, as doing so may be correcting some latent model bias arising from model misspecification (e.g. geometry simplification, flow-resistance free-surface relationships). However this is an assumption not determinable from the data, and may be considered overfitting the data. 

Mcfl parameters in this study were of two varieties--cross-sectional area parameters that have an exact physical interpretation, and flow-resistance parameters that do not, and merely arise out of simplified open-channel hydraulic equations. The parially optimized results differentiate between these parameter varieties, representing the geometric as its physical truth, such that $A_0 + A' = A$ for all times and locations. 

Differences between this study's results and those reported in the McFLI literature (Durand et al., 2016, Hagemann et al., 2017, ####) suggest that parameter estimation is the primary difficulty and source of error in producing inversions. Here we discuss reasons for this added difficulty, noting differences between mcfl and McFLI. 

The fundamental difference between the forward optimization taken in this paper and the inversion approach reported elsewhere, is the inclusion of discharge as a measured quantity indexed in time and space. This transforms the problem from
$$
\text{Given } W_{st}, A'_{st}, S_{st}, \text{ optimize } \theta \\
\text{Use optimal } \theta \text{ to estimate } \hat{Q}_{st}
$$
to 
$$
\text{Given } Q_{st}, W_{st}, A'_{st}, S_{st}, \text{ optimize } \theta \\
\text{Use optimal } \theta \text{ to estimate } \hat{Q}_{st}
$$

The first problem is underconstrained--there is not sufficient information in the given variables to uniquely determine the desired quantities. Their estimation therefore requires external constraints (e.g. via Bayesian priors), and is performed with substantial uncertainty. The second problem, however, is fully constrained and results in a unique optimum. It should not be surprising, then, that forward application of mass-conserved flow laws achieves higher performance than their "inversion", and that parameter error is the main culprit for the disparity. 