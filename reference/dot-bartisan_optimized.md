# Whether the installed shared library was compiled with optimization

Compilers define `__OPTIMIZE__` when they are optimizing, so this is
exact rather than a guess. It exists because an unoptimized build of
this package is between five and twenty times slower, and nothing else
about it looks wrong – which makes it very easy to spend a long time
drawing conclusions from the wrong numbers.

## Usage

``` r
.bartisan_optimized()
```

## Value

`TRUE` if the library was optimized.
