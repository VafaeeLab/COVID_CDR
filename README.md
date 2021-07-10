# COVID_CDR
**COVID-19 drug combination and their network-based Mechanism-of-Action (MOA) explorer**

## Usage
- Data availability: [```data directory```](https://github.com/VafaeeLab/COVID_CDR/tree/main/data)
- Code availability: [```script directory```](https://github.com/VafaeeLab/COVID_CDR/tree/main/scripts)

## Introduction
<p>The war against COVID-19 is far from over. Many Countries if immune to the first wave of COVID-19 impact are now facing the second, but stronger wave of these infections and much higher death rates. As per projected data there will be 300 million infections and more than 2 million deaths caused by COVID-19 globally by the end of March 2021. Although considerable scientific attention has been focused on identifying a cure for covid19 there has been a limited progress, few drugs are showing some promise in clinical trials, while others are raising controversy. </p>

<p>Overall, there is yet no licensed treatment or vaccine available to prevent or treat this deadly virus. The available therapeutic options are merely the off-label or repurposed drugs making clinicians and researchers more confused to make best choices as per the symptoms of the patient. Here we represent a web-based platform which is designed to help the personnel's in choosing the possible choice of drugs for treating COVID-19 patients. The platform displays the network-based drug mechanism of action thus simultaneously giving a visual of the cellular pathways involved in the mode of action of the chosen drugs. The platform also allows the freedom to choose more than two drugs and allows the options to understand the similarity score, contraindications of using these drugs</p>

<p>Here, we developed COVID-CDR, a unique platform to observe the network-based Mechanism of action (MOA) of <b>867 drugs </b> that are reported to treat COVID-19 symptoms.</p>

## Demo (click to watch)
[![IMAGE ALT TEXT HERE](https://img.youtube.com/vi/UsWSHu-wNCM/0.jpg)](https://youtu.be/tAd03VqbnXE)

## Acknowledgement
If you find COVID-CDR as useful for your research, please cite our work by including the following citation:
- <b>An Integrative Resource for Network-Based Investigation of COVID-19 Combinatorial Drug Repositioning and Mechanism of Action.</b> <i>ChemRxiv</i> [```Link to the paper```](https://chemrxiv.org/engage/chemrxiv/article-details/60c7523ff96a00638428817e)
- Citation:
```
@article{azad_fatima_vafaee_2020, 
  place={Cambridge}, 
  title={An Integrative Resource for Network-Based Investigation of COVID-19 Combinatorial Drug Repositioning and Mechanism of Action}, 
  DOI={10.26434/chemrxiv.13271096.v1}, 
  journal={ChemRxiv}, 
  publisher={Cambridge Open Engage}, 
  author={Azad, AKM and Fatima, Shadma and Vafaee, Fatemeh}, year={2020}
}
```
## Feature
- Pre-clinical exploration a network-based MOA of a COVID-19 drug combination
- Network-pharamacology based quantification of a drug combination (e.g. disease proximity, network separation, functional relevance, etc.)
- Numerous meta-data related to drugs, their target proteins or their PPI-partners of interest, that are either derived in COVID-CDR or retried from other databases through RESTful APIs.
- Highly queryble and user-friendly GUI
- Explore and download the whole COVID-CDR database including intermediate data (for developers)
- Explore and download the multi-modal similarity database of COVID-19 drugs
- Explore the COVID-CDR variant of other popular databases (e.g. DrugComDB and FDA-approaved combination)

## Dataset
- A set of [```867 drugs```](https://github.com/VafaeeLab/COVID_CDR/blob/main/data/allCOVID_drug_collection_full_v5.csv) that are reported to treat COVID-19 patients by various sources (i.e. ClinicalTrials, DrugBank, and Wikipedia)
- A list of other datasets
