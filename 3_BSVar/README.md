# Between-study variability

## Introduction

In this section, we provide an empirical assessment of the amount of between-study variability that can be observed in fMRI meta-analyses.  
We will focus on image-based meta-analyses where we use whole brain images containing a test statistics in each voxel. The test statistics are transformed to standardized effect sizes and between-study heterogeneity is estimated using several methods.

We obtain data using [neurovault](www.neurovault.org) where we focus on the experience of *"pain"*.

Data is obtained by searching for *"pain"* (at 11/01/2017) in the [neurovault](www.neurovault.org) database and manually checking all results.

## Mask
Our goal is to estimate between-study heterogeneity using whole brain statistical parametric maps on the topic *Pain* versus *No pain*.
We will restrict our analysis to brain areas known to be involved in pain processing. To get these areas, we create a mask using [neurosynth](www.neurosynth.org) searching for the term **pain**.
We obtain an automated meta-analysis of 420 studies. From here we use forward inference with a FDR at 0.01. Forward inference is equal to: P(Activation|Term)<sup>1</sup>.
This is the map called *pain_pAgF_z_FDR_0.01_forward_mask.nii*. After estimating the between-study heterogeneity, we will use this mask to obtain the distribution.

## Database

The *raw* data is stored in the *1_Data/1_Raw* folder. However, this folders is not pushed to the Github repository as they are too large. <br>
To resolve, we push the processed data to *1_Data/2_Processed*. These are **R** objects that can be read in to reproduce the figures of the paper.

The database contains the following studies:


| Study        | Sample size           | Type  |  Contrast |
| ------------- |:-------------|:-----:|:-----:|
|Braboszcz 2017 | 17            | T-map | Painful images >Painless images (Normal state) |
| Hebestreit 2017      | 23      |   T-map | Ammonia stimulation (trigeminal pain) > Air in both sessions (medication and placebo) |
| Tamm 2017 | 86      |  T-map | Pain > no pain [no covariates] |
| Karjalainen 2017 | 35      |  T-map | Main effect of vicarious pain |
| Atlas 2010 | 15      |  beta-map<sup>2</sup> | Thermal high vs low stimulation pain |
| Wager 2013 | 15      |  beta-map<sup>2</sup> | Somatic pain vs baseline |
| Kano 2017 | 15      |  beta-map<sup>2</sup> | Visceral pain vs baseline |
| Rubio 2015 | 15      |  beta-map<sup>2</sup> | Visceral pain vs baseline |
| Unpublished | 15      |  beta-map<sup>2</sup> | Mechanical high pressure pain vs baseline |
| Unpublished | 15      |  beta-map<sup>2</sup> | Mechanical medium pressure pain vs baseline |
| Patil 2017 | 46      |  T-map | Outcome of pain induced to others versus baseline |
| Maumet 2016 | Total = 334<sup>3</sup>      |  T-maps | Pain versus baseline |
| Chang 2015 | 28      |  T-maps | High pain versus low pain |


* Note that we only include studies doing whole-brain analyses.
* <sup>2</sup>Data comes from Kragel et al, 2017. These maps contain the contrasts of the parameter estimates for each subject. Hence we first pool these subjects using OLS.
* <sup>3</sup>Study of Maumet et al., 2016 is about the *nidm* data structure. However, it contains **21** studies about pain. Note that these come from the same site! We will need to control for this when estimating between-study variability.

## Structure



<sup>1</sup> http://neurosynth.org/faq/#q15
