
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mitoXpress

<!-- badges: start -->
<!-- badges: end -->

Vincent de Boer, Dec 12th, 2022

🎅 🔬 🤖 🥼

This is an import and analysis script written for TRF mitoXpress
experiments ran on a SpectraMax. Please have a good look at the data
file and key file.

## Install

For the script to run the following libraries are needed:

``` r
library(tidyverse)
library(broom)
library(here)
library(janitor)
library(lubridate)
library(readxl)
```

Also the renv lockfile in this repository is up to date with all
packages that are used in this project.

This repository can be cloned to your own github page or the files can
be copied to your own Rstudio project.

## Test data

The data and key files are living in the folder `inst\extdata`.

Experimental file:

- evelien_test_expt6.txt

Key file:

- key_file_template.xlsx

## Data file

A data file starts its data column names on row three. The script skips
the first two lines. The file contains a `Time` column and a
`Temperature` column followed by one column for each `well`. The script
removes rows that contain `NaN` and also removes the closing rows of the
documents.

The data file is a tab-delimited text file. Apparently our Spectramax
spits out txt files in UTF16 unicode format, which is somewhat outdated,
since practically all unicode text files now are encoded as UTF8. This
is likely the reason why Excel has difficulty reading the exported files
from the spectramax properly.

The tab-delimited text column names are also renamed. Please notice that
well names need to be three characters, so `A1` becomes `A01`. This
conversion is also incorporated in the script since the plate reader
outputs two and three digit well names. This is annoying when sorting
data and plotting data.

The imported and processed data looks like this in R:

``` r
mito_express_data
#> # A tibble: 610 × 5
#>    time     time_sec  temp well  fluorescence
#>    <Period>    <dbl> <dbl> <chr>        <dbl>
#>  1 0S              0    27 C08         42318.
#>  2 0S              0    27 C09         50590.
#>  3 0S              0    27 D08         27063.
#>  4 0S              0    27 D09         27544.
#>  5 0S              0    27 E08         13598.
#>  6 0S              0    27 E09         12513.
#>  7 0S              0    27 F08          9162.
#>  8 0S              0    27 F09          8408.
#>  9 0S              0    27 G08          7690.
#> 10 0S              0    27 G09          8232.
#> # … with 600 more rows
```

## Key file

The `key file` is needed for the user to manually add meta info to the
analysis. The most important is the choice for time frame for `slope`
calculation. The time frame is given in seconds (`start` and `end`).
Also `sample_type` and `group` information should be added.

The key file is an excel file (.xlsx file).

The key file in R looks like:

``` r
key_file
#> # A tibble: 10 × 6
#>    well  sample_type      group              window  start   end
#>    <chr> <chr>            <chr>              <chr>   <dbl> <dbl>
#>  1 C08   sample           GOx 10U            window1     0   120
#>  2 C09   sample           GOx 10U            window1     0   120
#>  3 D08   sample           GOx 5U             window1     0   120
#>  4 D09   sample           GOx 5U             window1     0   120
#>  5 E08   sample           GOx 1U             window2   300   780
#>  6 E09   sample           GOx 1U             window2   300   780
#>  7 F08   blank no enzyme  GOx 0U             window2   300   780
#>  8 F09   blank no enzyme  GOx 0U             window2   300   780
#>  9 G08   blank no glucose GOx 10U no glucose window2   300   780
#> 10 G09   blank no glucose GOx 10U no glucose window2   300   780
```

## Plots

Two examples of very basic plots are also in the script

This is the output for the example data:

``` r
df <- readRDS(here::here("data", "test_data.rds"))

df %>%
  ggplot(aes(x = time_sec, y = fluorescence, group = well, color = well))+
  geom_point()+
  labs(y = "fluorescence (AU)",
       x = "time (sec)")+
  theme_bw()
```

![](README_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
slopes <- readRDS(here::here("data", "test_slopes.rds"))

slopes %>%
  ggplot(aes(x = well, y = slope))+
  geom_point()+
  labs(y = "slope (dAU/dt)",
       x = "well name")+
  theme_bw()
```

![](README_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

And an example of somewhat more complex ggplot, with the time frame used
for slope calc also plotted.

![](README_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

## Export

The imported data and the calculated slopes can be exported to excel (as
.csv file) using the commands described in the script. Please notice
that the files are written to a folder named `output`.

``` r
df_wider <- df %>%
  pivot_wider(names_from = well, values_from = fluorescence) %>%
  select(-time)

readr::write_csv(df_wider,
                 here::here("output",
                            paste0(str_sub(basename(filepath_mXp), end=-5), ".csv")))

readr::write_csv(slopes,
                 here::here("output",
                            paste0("slopes", str_sub(basename(filepath_mXp), end=-5), ".csv")))
```

## Points of attention

- Please pay attention to the location of the data files. I have the
  data files in the `inst\extdata` folder and files are written to
  `output`.

- The exported files are named according to the input file name.

- One function `get_range_with_key` is needed to calculate the slopes.
  This function needs to be in the global environment, meaning that it
  should be run first before the slopes pipe can run.
