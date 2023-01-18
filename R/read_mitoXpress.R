library(tidyverse)
library(broom)
library(here)
library(janitor)
library(lubridate)
library(readxl)
library(assertthat)

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
#
# you can ignore the warning:
#
# • `` -> `...99`
# Warning message:
#   One or more parsing issues, call `problems()` on your data frame for details,
# e.g.:
#   dat <- vroom(...)
# problems(dat)
#
# other warnings cannot be ignored

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

readr::write_txt(df_wider,
                 here::here("output",
                            paste0(str_sub(basename(filepath_mXp), end=-5), ".txt")))

# plot the data for all wells ---------------------------------------------

df %>%
  ggplot(aes(x = time_sec, y = fluorescence, group = well, color = well))+
  geom_point()+
  labs(y = "fluorescence (AU)",
       x = "time (sec)")+
  theme_bw()

# calculate the slopes ----------------------------------------------------

# y = ax+b
# fluorescence = a*time_sec + b
# a = 24.9
# b = intercept = 22857.7

#function needed for window selection
get_range_with_key <- function(well_df, my_well, my_key){

  my_range <- my_key %>%
    filter(well == my_well) %>%
    select(start, end)

  well_df <- well_df %>%
    filter(time_sec >= my_range$start,
           time_sec <= my_range$end)

  return(well_df)

}

#get slope with key file
slopes <- df %T>% #please note the Tee pipe here "%T>%" works for the assertion here (from https://github.com/hadley/assertthat/issues/41)
  assert_that(identical(well %>% unique(),
                        key_file$well %>%  unique()), env = .) %>%
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


# plot slopes -------------------------------------------------------------

slopes %>%
  ggplot(aes(x = well, y = slope))+
  geom_point()+
  labs(y = "slope (dAU/dt)",
       x = "well name")+
  theme_bw()


# export slopes -----------------------------------------------------------

# filename is set automatically using the initial input filename
readr::write_delim(slopes,
                 here::here("output",
                            paste0("slopes", str_sub(basename(filepath_mXp), end=-5), ".txt")))



# plot slopes in raw data -------------------------------------------------

df %>%
  ggplot( mapping = aes(x = time_sec, y = fluorescence,
                       group = well, color = well))+
  geom_point(alpha = 0.9) +
  geom_smooth(data = . %>%
                filter(well %in% c("E08", "E09", "F08", "F09", "G08", "G09")) %>%
                filter(time_sec >= 300,
                       time_sec <= 780),
              method = lm,
              se = FALSE,
              formula = y ~ x,
              color = "grey50",
              size = 1.2)+
  geom_smooth(data = . %>%
                filter(well %in% c("C08", "C09", "D08", "D09")) %>%
                filter(time_sec >= 0,
                       time_sec <= 120),
              method = lm,
              se = FALSE,
              formula = y ~ x,
              color = "grey50",
              size = 1.2)+
  scale_color_manual(values =c(colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(10)))+
  labs(y = "fluorescence (AU)",
       x = "time (sec)")+
  theme_bw()+
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 16))





