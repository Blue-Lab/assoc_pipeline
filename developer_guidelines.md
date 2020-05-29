# Developer Guidelines

This document is intended to specify guidelines, style and
best practices for contributors. 

## Use RDS for input/output

Formatting and compatability can get lost converting to and from formats
like .txt and .csv. Use RDS with the functions `saveRDS()` and `readRDS()`.

## Use the `argparser` package for command line arguments

It allows for some expanded functionality not offered by `commandArgs()`, like
positional arguments.

## Use standardized script header

Each script should start with a few elements:

1. The standard R shebang: `#! /usr/bin/env Rscript`
1. Load `argparser` and other packages needed to parse the command line
arguments
1. Parse command line arguments with `arg_parser()`, `add_argument()` and
`parse_args()`
1. Load other required packages. These are further down because the scripts
help page will run more quickly if this is after the arguments are parsed.
1. Call `sessionInfo()` and print the arguments.

See a toy example below:

```r
#! /usr/bin/env Rscript
library(argparser)
library(magrittr)

# read arguments
argp <- arg_parser("Here is an example script header") %>%
  add_argument("arg1", help = "This is a positional argument") %>%
  add_argument("--arg2", help = "This is an optional argument with a default",
               default = "") %>%
  add_argument("--arg3", help = "This is the third argument")
argv <- parse_args(argp)

library(GENESIS)

sessionInfo()
print(argv)
```

## Follow [Hadley's style guide](http://adv-r.had.co.nz/Style.html)
