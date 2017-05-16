# gretaDFA

Using the `greta` package to fit Dynamic Factor Models.

## Background

Dynamic Factor Analysis (DFA) is a dimension reduction technique specific to time series analysis. The general idea is to model $N$ time series as a linear combination of $M$ "hidden" temporal trends, where $N >> M$. For an $N \times T$ matrix $\mathbf{y}$, the DFA model is

$$
\begin{gathered}
\mathbf{y}_t = \mathbf{Z}\mathbf{x}_t+\mathbf{a}+\mathbf{v}_t \text{ where } \mathbf{v}_t \sim \text{MVN}(0,\mathbf{R}) \\
\mathbf{x}_t = \mathbf{x}_{t-1}+\mathbf{w}_t \text{ where } \mathbf{w}_t \sim \text{MVN}(0,\mathbf{Q})
\end{gathered}   
$$

The matrix $\mathbf{Z}$ maps the latent trends onto the observed data. For convenience, the variance-covariance matrix of the process errors governing the multivariate random walk, $\mathbf{Q}$, is generally assumed to be an $M \times M$ identity matrix, $\mathbf{I}_M$. In addition, to make the model identifiable when $M > 1$, the upper right section of $\mathbf{Z}$ must be set equal to zero. For example, if $M = 3$, then 

$$
\mathbf{Z} = 
\begin{bmatrix}
    z_{11} & 0      & 0 \\
    z_{21} & z_{22} & 0 \\
    z_{31} & z_{32} & z_{33} \\
    z_{41} & z_{42} & z_{43} \\
    z_{51} & z_{52} & z_{53}
\end{bmatrix}   
$$

## Why `greta`?

Solving for the parameters in a DFA model is not trivial. There are several available methods for doing so, including an EM algorithm as implemented in the [`MARSS`](https://cran.r-project.org/web/packages/MARSS/index.html) package, or an MCMC algorithm as implemented in the __Stan__ language via the [`statss`](https://github.com/nwfsc-timeseries/statss) package. However, the fitting process can be extremely slow when working with really large data sets (e.g., 50+ time series).

The [`greta`](https://github.com/goldingn/greta) package was designed by [Nick Golding](https://scholar.google.co.uk/citations?user=peoal7wAAAAJ&hl=en) to do [Hamiltonian Monte Carlo](http://onlinelibrary.wiley.com/doi/10.1111/2041-210X.12681/full) via Google's [TensorFlow](https://www.tensorflow.org/) computational engine, which makes it really fast on big datasets.

## Installation

You can install `greta` from GitHub using `devtools`

```{r}
devtools::install_github('goldingn/greta')
```

You will also need to install TensorFlow (version 1.0.0 or higher) before `greta` will work. Check [here](https://www.tensorflow.org/install/) for instructions on installing TensorFlow.

## Example


