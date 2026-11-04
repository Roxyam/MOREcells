![](https://capsule-render.vercel.app/api?type=waving&color=0:98A6D4,100:a9d498&customColorList=2&height=320&section=header&text=MOREcells&fontSize=58&animation=fadeIn&fontAlignY=32&desc=Regulatory%20network%20inference%20in%20single-cell%20multiomic%20scenarios&descSize=20&descAlignY=52&descAlign=50&fontColor=ffffff&descColor=ffffff)

<table>
  <tr>
    <td>
      <img src="man/figures/Logo_MOREcells.png" width="250">
    </td>
    <td>
      <h1>MOREcells</h1>
      <p>Regulatory network inference in single-cell multiomic scenarios</p>
    </td>
  </tr>
</table>

`MOREcells` is an R framework designed to analyze gene regulation and multi-omic single-cell data by modeling target omic expression as a function of experimental variables (such as biological conditions, treatments, and cell types) alongside potential regulators of a given target feature, while fully maintaining single-cell resolution and hierarchical data structure.  

To achieve this, MOREcells utilizes a two-step approach that combines Generalized Linear Mixed Models (GLMMs) to account for donor-specific random effects and prevent pseudoreplication, with flexible post-hoc inference options that allow researchers to infer phenotype-specific regulatory networks and perform custom statistical contrasts across different biological groups. 

<div align="center">
  <img src="man/figures/Abstract.png" alt="Graphical Abstract" width="700" />
</div>

## Installation

Currently, the package can be installed directly from GitHub:  

```r
# Install from GitHub (requires devtools) 
devtools::install_github("BiostatOmics/MOREcells")
```

or 

```r
# Install from GitHub (requires remotes)
remotes::install_github("BiostatOmics/MOREcells")
```

---

![](https://capsule-render.vercel.app/api?type=waving&color=0:98A6D4,100:a9d498&customColorList=2&height=200&section=footer&desc=Regulatory%20network%20inference%20in%20single-cell%20multiomic%20scenarios&descSize=20&descAlignY=81&descAlign=62&fontColor=ffffff&descColor=ffffff)
