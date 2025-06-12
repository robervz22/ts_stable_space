# ts_stable_space

A code repository to reproduce the results presented in the paper *Over the Stability Space of a Multivariate Time Series*.

## Table of Contents

-   [Description](#description)
-   [Usage](#usage)
-   [Example](#example)
-   [Contributing](#contributing)
-   [License](#license)

## Description

The repository contains four main directories:

-   `source`: This directory contains the source code for the algorithms implemented based on the article *Over the Stability Space of a Multivariate Time Series*. It includes the following files:

> `vectorial_methods.R`: Implements the three methodologies proposed in the article for estimating the stable space of a multivariate time series.  
>
> `simulations.R`: Contains the main functions required to reproduce the simulation study presented in the article.  
>
> `auxiliar_methods.R`: Includes auxiliary functions for visualisation and error estimation.  

-   `notebooks`: This directory contains Jupyter notebooks that demonstrate the use of the source code and provide step-by-step instructions to reproduce the various tables and figures from the paper.

-   `databases`: Contains all datasets used in the `notebooks` directory, including those for the practical examples.

-   `images`: Contains the exact images used in the article. These can be reproduced using the notebooks, with some light post-processing. For instance, `Figure_1.png` and `Figure_2.png` are the result of running the notebook `stability_space.ipynb` four times and then combining the corresponding images.

## Usage

To run the notebooks in the `notebooks` directory and the scripts in the `source` directory, use `R` version `4.1.2`.

Enjoy it ☕

## Example

Suppose we want to estimate the stability space of the following dataset using the three methodologies proposed in the article:

```R
dt_inflation <- data.table::fread("../databases/variables_inflation.csv")
```

First, we define the input matrices for each method:

```R
X_inflation <- as.matrix(dt_inflation, rownames = "Date")
rownames(X_inflation) <- NULL
X_inflation <- scale(ts(X_inflation))
XX_inflation <- X_inflation[1:(nrow(X_inflation)-1),]
Y_inflation <- X_inflation[2:nrow(X_inflation),]
```

Then, we run the following code lines—each corresponds to one of the methods:

```R
basis_inflation_PLS <- basis_stable(X_inflation, method = "pls")       # PLS method
basis_inflation_PCA <- basis_stable(XX_inflation, method = "pca")      # PCA method
basis_inflation_johansen <- basis_stable(XX_inflation, method = "johansen")  # Johansen procedure
```

Each `basis_inflation_...` object contains two elements: `basis_S` and `basis_N`, which represent the estimators of the stable space and the non-stable space basis, respectively.

## Contributing

1.  Fork the repository.  
2.  Create a new branch (`git checkout -b feature-branch`).  
3.  Commit your changes (`git commit -m 'Add a new feature'`).  
4.  Push to the branch (`git push origin feature-branch`).  
5.  Open a Pull Request.

## License

This project is licensed under the MIT License – see the [LICENSE](LICENSE) file for details.
