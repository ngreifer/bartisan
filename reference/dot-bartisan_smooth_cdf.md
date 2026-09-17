# Smoothed empirical distribution function, and the sums a bandwidth selector needs

The empirical distribution function convolved with a kernel, evaluated
on a grid, together with two sums over the same kernel weights that the
bandwidth selectors of Bergmann and Zaehle (2026) are built from.

## Usage

``` r
.bartisan_smooth_cdf(x, grid, h, kernel)
```

## Arguments

- x:

  `numeric`; the observations, **sorted ascending**.

- grid:

  `numeric`; where to evaluate, **sorted ascending**.

- h:

  `numeric`; the bandwidth, strictly positive.

- kernel:

  `integer`; 0 for Epanechnikov, 1 for Gaussian.

## Value

A matrix with one row per grid point and three columns: the smoothed
distribution function, the sum of squared kernel weights, and the sum of
squared differences between each kernel weight and the step it smooths.
The last two are the terms that distinguish the two selectors.

## Details

This is in compiled code because the obvious way to write it is
quadratic: one pass over every observation for every grid point is
`n * m` kernel evaluations, which is tens of millions per predictor on a
large fit, and a bandwidth selector repeats the whole thing at every
candidate bandwidth.

Almost all of that work is avoidable, because almost every weight is a 0
or a 1 rather than something in between. An observation further below a
grid point than the kernel reaches has passed it entirely and
contributes a weight of 1; one further above has not been reached and
contributes 0. Only those within reach need evaluating. With the data
and the grid both sorted, the two ends of that window move forward
monotonically as the grid advances, so each is found by a pointer that
never goes back, and the sweep costs `O(n + m)` plus one kernel
evaluation per observation actually inside a window.

The Epanechnikov kernel makes this exact rather than approximate, its
support being `[-1, 1]`, so the window is `h` either side and nothing
outside it is discarded. The Gaussian reaches everywhere and is cut off
at five bandwidths, beyond which the weight being rounded to 0 or 1 is
under `3e-7`.
