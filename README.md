# MetaSpread: a cancer and metastasis simulation package

MetaSpread is a cancer and metastasis simulator written in Python, based on the work of Franssen et al [[1]](#1) and the MESA framework [[2]](#2). With this package it is possible to simulate how cancer cells reproduce and diffuse locally, and generate metastases by spreading to other parts of the body. The underlying model is a reaction-diffusion type model with deterministic and stochastic aspects, tracking the spatiotemporal evolution of: i) numbers of epithelial and mesenchymal cancer cells (modelled as discrete agents); and ii) the densities of the extracellular matrix and the diffusible matrix metalloproteinase-2  (MMP-2), in both the primary site and all secondary sites. Other processes represented in the simulator are the circulation of cancer cell clusters and single cancer cells in the vasculature, including cell death and disaggregation of clusters, and the extravasation of clusters and single cells to distant metastatic sites.

# Installing the package

MetaSpread is available as an official PyPi package. To install, simply run:

```
pip install metaspread
```

If you want to install manually, you can also download the source code as a zip file, unzip in an appropiate directory, and run

```
python -m metaspread
```

# References

><a id="1">[1]</a>  Franssen, L. C., Lorenzi, T., Burgess, A. E., & Chaplain, M. A. (2019). A mathematical framework for modelling the metastatic spread of cancer. Bulletin of Mathematical Biology, 81, 1965–2010
><a id="2">[2]</a>  Masad, D., & Kazil, J. (2015). MESA: An agent-based modeling framework. 14th PYTHON in Science Conference, 2015, 53–60.
><a id="3">[3]</a>  Hernández-Inostroza, A. and Gjini, E., MetaSpread: A cancer growth and metastatic spread simulation program in Python, 2024 (submitted).


