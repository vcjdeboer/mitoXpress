library(tidyverse)
library(broom)
library(here)
library(janitor)
library(lubridate)
library(readxl)

# set the data file -------------------------------------------------------

filepath_mXp <- here::here("inst", "extdata", "evelien_test_expt6.txt")

# key file ----------------------------------------------------------------

# key file needs the following columns
# well, sample_type, group, window, start, end
# see example key template for details

filepath_keyfile <- here::here("inst", "extdata", "key_file_template.xlsx")

key_file <- readxl::read_excel(filepath_keyfile)

# import data file and preprocess -----------------------------------------

#import is specifically for TRF experiment on Spectramax

df <- filepath_mXp %>%
  readr::read_tsv(locale = readr::locale(encoding = "UTF-16LE"),
                  col_names = TRUE,
                  show_col_types = FALSE,
                  skip = 2) %>%
  filter(Time != "NaN",
         Time != "~End") %>%
  filter(!str_detect(Time, "Original")) %>%
  janitor::clean_names(case = "none") %>%
  rename(time = Time,
         temp = Temperature_C) %>%
  pivot_longer(c(-time, -temp),
               names_to = "well",
               values_to = "fluorescence") %>%
  mutate(well = well %>%
           map_chr(~paste0(
             str_extract(.x, "\\D"),
             str_extract_all(.x, "\\d") %>%
               pluck(1) %>%
               str_flatten() %>%
               str_pad(width = 2, pad = "0")))) %>%
  mutate(time = lubridate::hms(time),
         time_sec = lubridate::period_to_seconds(time)) %>%
  select(time, time_sec, temp, well, fluorescence) %>%
  drop_na()

# export to excel (as csv file) -------------------------------------------

df_wider <- df %>%
  pivot_wider(names_from = well, values_from = fluorescence) %>%
  select(-time)

readr::write_csv(df_wider,
                 here::here("inst", "extdata",
                            paste0(str_sub(basename(filepath_mXp), end=-5), ".csv")))

# plot the data for all wells ---------------------------------------------

df %>%
  ggplot(aes(x = time_sec, y = fluorescence, group = well, color = well))+
  geom_point()

# calculate the slopes ----------------------------------------------------

# y = ax+b
# fluorescence = a*time_sec + b
# a = 24.9
# b = intercept = 22857.7

# manually set the time windows for each well
window1 <- c(0:120)
window2 <- c(300:780)

windows <- tibble(
  well = df$well %>% unique(),
  window = c("window1", "window1", "window1", "window1",
             "window2", "window2", "window2", "window2",
             "window2", "window2"))

# range filter functions
get_range_manual <- function(well_df, my_well, windows){

  my_window <- windows %>%
    filter(well == my_well) %>%
    pull(window)

  well_df <- well_df %>%
    filter(time_sec %in% get(my_window))

  return(well_df)

}

get_range_with_key <- function(well_df, my_well, my_key){

  my_range <- my_key %>%
    filter(well == my_well) %>%
    select(start, end)

  well_df <- well_df %>%
    filter(time_sec >= my_range$start,
           time_sec <= my_range$end)

  return(well_df)

}


# get slope with manual windows
df %>%
  group_by(well) %>%
  nest() %>%
  mutate(model = map2(.x = data,
                      .y = well,
                      ~.x %>% get_range(.y, windows) %>%
                       lm(formula = fluorescence ~ time_sec, data = .)),
         summaries = map(model, broom::glance),
         model_coef = map(model, broom::tidy)) %>%
  mutate(slope = map_dbl(model_coef, ~.x %>% pluck("estimate", 2))) %>%
  select(well, slope)

#get slope with key file
df %>%
  group_by(well) %>%
  nest() %>%
  mutate(model = map2(.x = data,
                      .y = well,
                      ~.x %>% get_range_with_key(.y, key_file) %>%
                        lm(formula = fluorescence ~ time_sec, data = .)),
         summaries = map(model, broom::glance),
         model_coef = map(model, broom::tidy)) %>%
  mutate(slope = map_dbl(model_coef, ~.x %>% pluck("estimate", 2))) %>%
  select(well, slope)






