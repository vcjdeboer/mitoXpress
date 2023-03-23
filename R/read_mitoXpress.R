library(tidyverse)
library(broom)
library(here)
library(janitor)
library(lubridate)
library(readxl)
library(assertthat)
library(bayestestR)
library(pspline)
library(colorspace)
library(ggrepel)


library(conflicted)

conflicts_prefer(dplyr::filter())

# set the data file -------------------------------------------------------

filepath_mXp <- here::here("inst", "extdata", "evelien_test_expt6.txt")

filepath_mXp <- here::here("inst", "extdata", "20230306_5dpfoptimalisation.txt")


# key file ----------------------------------------------------------------

# key file needs the following columns
# well, sample_type, group, window, start, end
# see example key template for details

filepath_keyfile <- here::here("inst", "extdata", "key_file_template.xlsx")

filepath_keyfile <- here::here("inst", "extdata", "key_file_exp20230306_5dpfoptimalisation.xlsx")

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

readr::write_delim(df_wider,
                 here::here("output",
                            paste0(str_sub(basename(filepath_mXp), end=-5), ".txt")))

# plot the data for all wells ---------------------------------------------

df %>%
  left_join(key_file, by = c("well")) %>%
  ggplot(aes(x = time_sec, y = fluorescence, group = group, color = group))+
  geom_point()+
  labs(y = "fluorescence (AU)",
       x = "time (sec)")+
  theme_bw()+
  facet_wrap(~group)

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

# issue output folder not on my comp

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
              linewidth = 1.2)+
  geom_smooth(data = . %>%
                filter(well %in% c("C08", "C09", "D08", "D09")) %>%
                filter(time_sec >= 0,
                       time_sec <= 120),
              method = lm,
              se = FALSE,
              formula = y ~ x,
              color = "grey50",
              linewidth = 1.2)+
  scale_color_manual(values =c(colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(10)))+
  labs(y = "fluorescence (AU)",
       x = "time (sec)")+
  theme_bw()+
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 16))


#AUC blubs

background_drift <- function(x,y){
  y <- y/min(y)-1
  auc_func <- bayestestR::area_under_curve(x, y, method = "trapezoid")
  return(auc_func)
}

background_drift2 <- function(x,y, baseline){
  y <- y-baseline
  auc_func <- bayestestR::area_under_curve(x, y, method = "trapezoid")
  return(auc_func)
}

background_drift3 <- function(x,y){
  #y <- y-y[1]
  auc_func <- bayestestR::area_under_curve(x, y, method = "trapezoid")
  return(auc_func)
}

df %>%
  filter(well == "E09") %>%
  ggplot(aes(x = time_sec, y = fluorescence))+
  geom_point()

test_data <- df %>%
  filter(well == "E09")

test_data2 <-  df %>%
  left_join(key_file, by = c("well")) %>%
  filter(well == "E02")

background_drift3(test_data$time_sec, test_data$fluorescence)/max(test_data$time_sec)


y = test_data$fluorescence

test_data %>% ggplot(aes(x = time_sec, y = y))+
  geom_point()

test_data %>% mutate(scaled = scale(fluorescence, center = TRUE, scale = FALSE)) %>%
  ggplot(aes(x = time_sec, y = scaled[,1]))+geom_point()

df %>%
  left_join(key_file, by = c("well")) %>%
  ggplot(aes(x = time_sec, y = fluorescence, group = group, color = well))+
  geom_point()+
  labs(y = "fluorescence (AU)",
       x = "time (sec)")+
  theme_bw()+
  facet_wrap(~group)

scaling_factors <- df %>%
  left_join(key_file, by = c("well")) %>%
  filter(time_sec == 0) %>%
  group_by(group) %>%
  mutate(group_mean = mean(fluorescence)) %>%
  ungroup() %>%
  mutate(relative = fluorescence/group_mean) %>%
  select(well, factor = relative)

df %>%
  left_join(key_file, by = c("well")) %>%
  left_join(scaling_factors, by = c("well")) %>%
  mutate(new_fl = fluorescence/factor) %>%
  ggplot(aes(x = time_sec, y = new_fl, group = group, color = group))+
    geom_point()+
    labs(y = "fluorescence (AU)",
         x = "time (sec)")+
    theme_bw()+
    facet_wrap(~group)

df %>%
  left_join(key_file, by = c("well")) %>%
  filter(group == "5dpf + 1hour inc.") %>%
  ggplot(aes(x = time_sec, y = fluorescence, group = group, color = group))+
  geom_point()+
  geom_smooth(formula = y~ x, method = "loess")+
  labs(y = "fluorescence (AU)",
       x = "time (sec)")+
  theme_bw()+
  facet_wrap(~well, scales = "fixed")

splineFN_new <- function(x, n_order = 0) {
  x <- x %>% select(time_sec, fluorescence)
  colnames(x) <- c("timescale", "param")
  y <- as.vector(predict(sm.spline(x$timescale, x$param),
                         x$timescale,
                         n_order))
  df <- tibble(timescale = x$timescale, parameter = y, fluorescence = x$param)
  return(df)
}

splineFN_new2 <- function(x, n_order = 0, spar = 0.8, method = 2) {
  x <- x %>% select(time_sec, fluorescence)
  colnames(x) <- c("timescale", "param")
  y <- as.vector(predict(smooth.Pspline(x$timescale, x$param,
                                        spar = spar,
                                        method = method),
                         x$timescale,
                         n_order))
  df <- tibble(timescale = x$timescale, parameter = y)
  return(df)
}


df %>%
  left_join(key_file, by = c("well")) %>%
  filter(well == "D02") %>%
  splineFN_new2(n_order = 1) %>%
  ggplot(aes(x = timescale, y = parameter))+
    geom_point(color = "red", alpha = 0.9)

df %>%
  left_join(key_file, by = c("well")) %>%
  filter(well == "G02") %>%
  splineFN_new2() %>%
  ggplot(aes(x = timescale, y = parameter))+
  geom_point(color = "red", alpha = 0.9)+
  geom_point(aes(y = fluorescence), alpha = 0.5)

transformed <- df %>%
  left_join(key_file, by = c("well")) %>%
  nest(.by = c(group, well)) %>%
  mutate(spline_0 = map(.x = data,
                            .f = ~splineFN_new2(.x) %>%
                                  select(zero_order = parameter))) %>%
  mutate(spline_1 = map(.x = data,
                        .f = ~splineFN_new2(.x,
                                            n_order = 1) %>%
                              select(first_order = parameter))) %>%
  mutate(auc_1 = map2_dbl(.x = data,
                     .y = spline_1,
                     .f = ~(background_drift3(.x$time_sec,.y$first_order)))) %>%
  unnest(c(data, spline_0, spline_1)) %>%
  select(group, well, time_sec, fluorescence, zero_order, first_order, auc_1)





transformed %>%
  distinct(well, .keep_all = TRUE) %>%
  ggplot(aes(x = group, y = auc_1, color = well))+
  geom_jitter()

transformed %>%
  distinct(well, .keep_all = TRUE) %>%
  #filter(!group %in% c("empty", "cellfree neg.")) %>%
  ggplot(aes(x = group, y = auc_1, color = well))+
  geom_point()+
  expand_limits(x = 0, y = 0)

transformed %>%
  summarize(mean_1 = mean(first_order), .by = c(well, group)) %>%
  ggplot(aes(x = group, y = mean_1, color = well))+
  geom_point()


transformed %>%
  nest(.by = c(well, group)) %>%
  mutate(max_time = map_dbl(.x = data,
                        .f = ~.x$time_sec[which.max(.x$first_order)])) %>%
  mutate(u = map2_dbl(.x = data,
                                    .y = max_time,
                                   .f = ~.x %>%
                                     filter(time_sec %in% c(0:.y)) %>%
                                     pull(first_order) %>%
                                     mean() #%>%
                                     #'/'(.y)
                                     )) %>%
  ggplot(aes(x = group, y = u, color = well))+
  geom_point() +
  theme_bw()

transformed %>%
  nest(.by = c(well, group)) %>%
  mutate(max_time = map_dbl(.x = data,
                            .f = ~.x$time_sec[which.max(.x$first_order)])) %>%
  mutate(u_per_time = map2_dbl(.x = data,
                      .y = max_time,
                      .f = ~.x %>%
                        filter(time_sec %in% c(0:.y)) %>%
                        pull(first_order) %>%
                        mean() %>%
                      '/'(.y)
  )) %>%
  ggplot(aes(x = group, y = u_per_time, color = well))+
  geom_point() +
  theme_bw()

transformed %>%
  nest(.by = c(well, group)) %>%
  mutate(max_time = map_dbl(.x = data,
                            .f = ~.x$time_sec[which.max(.x$first_order)])) %>%
  mutate(auc = map2_dbl(.x = data,
                               .y = max_time,
                               .f = ~.x %>%
                                 filter(time_sec %in% c(0:.y)) %>%
                                 select(time_sec, first_order)  %>%
                                 with(background_drift3(x = .$time_sec, y = .$first_order))
  )) %>%
  ggplot(aes(x = group, y = auc, color = well))+
  geom_point() +
  theme_bw()

transformed %>%
  nest(.by = c(well, group)) %>%
  mutate(max_time = map_dbl(.x = data,
                            .f = ~.x$time_sec[which.max(.x$first_order)])) %>%
  mutate(auc_time = map2_dbl(.x = data,
                        .y = max_time,
                        .f = ~.x %>%
                          filter(time_sec %in% c(0:.y)) %>%
                          select(time_sec, first_order)  %>%
                          with(background_drift3(x = .$time_sec, y = .$first_order) %>%
                          '/'(.y))
  )) %>%
  ggplot(aes(x = group, y = auc_time, color = well))+
  geom_point() +
  theme_bw()






get_processsed_mXp <- function(data_filename, key_filename ){

  #read keyfile
  filepath_keyfile <- here::here("inst", "extdata", key_filename )
  key_file <- readxl::read_excel(filepath_keyfile)

  #read data
  filepath_mXp <- here::here("inst", "extdata", data_filename)

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

  df <- df %>%
    left_join(key_file, by = c("well")) %>%

  return(df)

}
get_transformed_mXp <- function(df){

  transformed <- df %>%
  #df %>%
    nest(.by = c(group, well)) %>%
    mutate(spline_0 = map(.x = data,
                          .f = ~splineFN_new2(.x) %>%
                            select(zero_order = parameter))) %>%
    mutate(spline_1 = map(.x = data,
                          .f = ~splineFN_new2(.x,
                                              n_order = 1) %>%
                            select(timescale, first_order = parameter))) %>%
    mutate(auc_1 = map2_dbl(.x = data,
                            .y = spline_1,
                            .f = ~(background_drift3(.x$time_sec,
                                                     .y$first_order)))) %>%
    mutate(model =  map(.x = data,
                        ~.x %>%
                          filter(time_sec %in% c(800:2800)) %>%
                          lm(formula = fluorescence ~ time_sec, data = .)),
           summaries = map(model, broom::glance),
           model_coef = map(model, broom::tidy)) %>%
    mutate(slope = map_dbl(model_coef, ~.x %>% pluck("estimate", 2))) %>%

    mutate(max_time = map_dbl(.x = spline_1,
                              .f = ~.x$timescale[which.max(.x$first_order)])) %>%

    mutate(log_model_till_max = map2(.x = data,
                                 .y = max_time,
                                 .f = ~.x %>%
                                   filter(time_sec %in% c(0:.y)) %>%
                                   lm(formula = log(fluorescence) ~ time_sec, data = .)),
                                       summaries_log = map(log_model_till_max, broom::glance),
                                       model_coef_log = map(log_model_till_max, broom::tidy)) %>%
    mutate(log_slope_till_max = map_dbl(model_coef_log, ~.x %>% pluck("estimate", 2) %>% '*'(1000))) %>%


    unnest(c(data, spline_0, spline_1)) %>%
    select(group, well, time_sec, fluorescence, zero_order, first_order, auc_1, slope, log_slope_till_max)

  return(transformed)
}

df <- get_processsed_mXp(data_filename, key_filename)

data_filename <- "20230306_5dpfoptimalisation.txt"
key_filename <- "key_file_exp20230306_5dpfoptimalisation.xlsx"

data_filename <- "20230313_3dpfwildtypeOCR.txt"
key_filename <- "key_file_exp20230313_3dpfwildtypeOCR.xlsx"

data_filename <- "20230314_4dpfwildtypeOCR.txt"
key_filename <- "key_file_exp20230314_4dpfwildtypeOCR_AUC.xlsx"

data_filename <- "20230315_5dpfwildtypeOCRfail.txt"
key_filename <- "key_file_exp20230315_5dpfwildtypeOCRfail_AUC.xlsx"

data_filename <- "20230315_5dpfwildtypeOCR2.txt"
key_filename <- "key_file_exp20230315_5dpfwildtypeOCR2_AUC.xlsx"

transformed <-
  get_processsed_mXp(data_filename, key_filename) %>%
  get_transformed_mXp()

#plot data
transformed %>%
  ggplot(aes(x=time_sec, y = (fluorescence), color = well, group = group))+
  geom_point()+
  geom_line(aes(y = (zero_order), color = well, group = well),
            alpha = 0.5,
            linewidth = 2)+
  ggrepel::geom_label_repel(data = . %>% filter(time_sec == max(time_sec)),
                            aes(label = well),
                            show.legend = FALSE)+
  colorspace::scale_colour_discrete_divergingx(palette = "Geyser", rev = FALSE)+
  labs(y = "fluorescence (AU)",
       x = "time (sec)") +
  facet_wrap(~group, scales = "fixed")+
  theme_bw(base_size = 24)

#plot first order
transformed %>%
  ggplot(aes(x=time_sec, y = first_order, color = well, group = group))+
  geom_point()+
  ggrepel::geom_label_repel(data = . %>% filter(time_sec == max(time_sec)),
                            aes(label = well),
                            show.legend = FALSE)+
  colorspace::scale_colour_discrete_divergingx(palette = "Geyser", rev = FALSE)+
  labs(y = "change in fluorescence (1st der.)",
       x = "time (sec)") +
  facet_wrap(~group, scales = "fixed")+
  theme_bw(base_size = 24)

transformed %>%
  filter(!group %in% c("empty", "cellfree neg.")) %>%
  ggplot(aes(x=time_sec, y = log(fluorescence), color = well, group = group))+
  geom_point()+
  geom_line(aes(y = log(zero_order), color = well, group = well),
            alpha = 0.5,
            linewidth = 2)+
  ggrepel::geom_label_repel(data = . %>% filter(time_sec == max(time_sec)),
                            aes(label = well),
                            show.legend = FALSE)+
  colorspace::scale_colour_discrete_divergingx(palette = "Geyser", rev = FALSE)+
  labs(y = "log fluorescence (AU)",
       x = "time (sec") +
  facet_wrap(~group, scales = "fixed")+
  theme_bw(base_size = 24)


transformed %>%
  filter(!group %in% c("empty")) %>%
  slice(1, .by = well) %>%
  ggplot(aes(x=group, y = log_slope_till_max, color = well, group = group))+
  geom_point(size = 4)+
  ggdist::stat_pointinterval(position = position_dodge(width = -0.4, preserve = "single"))+
  colorspace::scale_colour_discrete_divergingx(palette = "Geyser", rev = FALSE)+
  labs(y = "log slope (AU/1000sec)",
       x = "") +
  theme_bw(base_size = 24)+
  theme(panel.border = element_blank(),
        axis.line.y = element_line(),
        axis.line.x = element_line())







