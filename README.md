**A Bias Aware Multimodal Workflow for Seed Zone Maps and Conservation of Vitellaria-Paradoxa** 
Author: Gabriel Salako et al. 2025, PhD	Affiliation: Senckenberg Museum of Natural History, Görlitz, Germany
Correspondence: gabriel.saIako@kwasu.edu.ng
**Abstract**
Integrating heterogeneous ecological information—remote sensing, climate surfaces, plot surveys, and citizen■science contributions—poses major analytical challenges for species distribution modelling (SDM). Many workflows do not properly address sampling bias, spatial autocorrelation, predictor collinearity, or uncertainty communication, resulting in overconfident and poorly generalizable predictions. Here, I present a reproducible muIti■modaI SDM framework combining: (1) ensemble machine■Iearning models, (2) Remote sensing Vegetation index, (3) spatial cross■vaIidation, and (4) muIti■criteria thresholding and seed■zone identification. The workflow is demonstrated using Vitellaria paradoxa, an economically and ecologically important savanna tree species. Although tailored to this case, the pipeline is generalizable to any ecological modelling context requiring robust predictions under data heterogeneity.
**1.	Introduction**
Species distribution modelling has become fundamental to biodiversity conservation, ecological forecasting, and restoration planning. Yet the availability of large muIti■source datasets introduces new modelling complexities. Ecological occurrence data frequently combines structured field surveys (systematic, high quality), crowdsourced data (GBIF, iNaturalist; high volume but biased), environmental sensors, climate models, and remote sensing imagery. This heterogeneity increases the risk of bias misestimation, overfitting, and inflated accuracy metrics when simplistic cross■vaIidation schemes are used. To address these issues, this workflow merges ensemble ML predictive strength, and bias■aware spatial validation, implemented in R with fully reproducible materials.
**2.	Data Sources & Preprocessing**
Occurrence data: field observations (2021—2024) and GBIF records were cleaned (duplicate removal, coordinate checks) and spatially thinned at 1 km to reduce clustering.
Predictors: CHELSA climate variables; Landsat 8/9■derived indices (NDVI, GNDVI, SAVI, NDWI/MNDWI); Iand■surface and Iand■cover layers. All rasters were resampled to 30 m and aligned on a common grid to ensure consistency across data modalities.
**3.	Modelling Framework**
Ensemble ML: Random Forest (robust to collinearity; variable importance) and Boosted Regression Trees (nonlinear responses; flexible interactions) provide spatially rich suitability indices and fast screening of predictor relevance.
**4.	Outputs**
The pipeline generates: ensemble suitability rasters (RF/BRT), partial dependence plots, ecological response curves, and SZPI maps. Figures are exported as PNG/TIFF and can be used directly in decision support.
**5.	Reproducibility & Transparency**
The repository includes modular scripts, an RMarkdown report, and an renv.lock file for strict dependency control. All code and example outputs are openly available.
**6.	Discussion**
By combining ML prediction strength and spatially aware validation, this workflow reduces overconfidence and improves transferability—two common problems in SDM. It is well suited for ecological forecasting, monitoring, and climate adaptation and conservation planning. Future extensions include joint species models, eDNA integration, trait■based modelling, and deep learning for remote sensing trait extraction.

9.	Availability
Code and documentation: https://github.com/gabsaIako/Vitellaria-Paradoxa
