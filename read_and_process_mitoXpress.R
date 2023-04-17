#Vincent de Boer
#Thu, Mar 23rd 2023


# packages ----------------------------------------------------------------

library(tidyverse)
library(broom)
library(here)
library(janitor)
library(lubridate)
library(readxl)
library(bayestestR)
library(pspline)
library(colorspace)
library(ggrepel)
library(ggdist)

library(conflicted)
conflicts_prefer(dplyr::filter())

# functions ---------------------------------------------------------------

background_drift3 <- function(x,y){
  #y <- y-y[1]
  auc_func <- bayestestR::area_under_curve(x, y, method = "trapezoid")
  return(auc_func)
}

splineFN_new2 <- function(x, n_order = 0, method = 2, degrees = 10) {
  x <- x %>% select(time_sec, fluorescence)
  colnames(x) <- c("timescale", "param")
  y <- as.vector(predict(smooth.Pspline(x$timescale, x$param,
                                        method = method,
                                        df = degrees),
                         x$timescale,
                         n_order))
  df <- tibble(timescale = x$timescale, parameter = y)
  return(df)
}


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

get_transformed_mXp <- function(df, max_cutoff, degrees_of_freedom){

  # df <- df_test
  # max_cutoff = 8*60
  # degrees_of_freedom = 10

  transformed <- df %>%
    #slice(5:n(), .by = well) %>% #remove first 5 rows of each well
    nest(.by = c(group, well)) %>%
    mutate(spline_0 = map(.x = data,
                          .f = ~splineFN_new2(.x,
                                              degrees = degrees_of_freedom) %>%
                            select(zero_order = parameter))) %>%
    mutate(spline_1 = map(.x = data,
                          .f = ~splineFN_new2(.x,
                                              degrees = degrees_of_freedom,
                                              n_order = 1) %>%
                            select(timescale, first_order = parameter))) %>%
    mutate(spline_2 = map(.x = data,
                          .f = ~splineFN_new2(.x,
                                              degrees = degrees_of_freedom,
                                              n_order = 2) %>%
                            select(timescale, second_order = parameter))) %>%
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
    mutate(max_time_2 = map_dbl(.x = spline_2,
                              .f = ~.x$timescale[which.min(.x$second_order)])) %>%
    mutate(time_end = map2_dbl(.x = max_time,
                              .y = max_time_2,
                              .f = ~ifelse(.x>(max_cutoff-1), .x, .y))) %>%
    mutate(log_model_till_max = map2(.x = data,
                                     .y = time_end,
                                     .f = ~.x %>%
                                       filter(time_sec %in% c(0:.y)) %>%
                                       lm(formula = log(fluorescence) ~ time_sec, data = .)),
           summaries_log = map(log_model_till_max, broom::glance),
           model_coef_log = map(log_model_till_max, broom::tidy)) %>%
    mutate(log_slope_till_max = map_dbl(model_coef_log, ~.x %>% pluck("estimate", 2) %>% '*'(1000))) %>%
    #this is done a second time because including time_scale, messess up the unnesting
    mutate(spline_2 = map(.x = data,
                          .f = ~splineFN_new2(.x,
                                              n_order = 2) %>%
                            select(second_order = parameter))) %>%
    unnest(c(data, spline_0, spline_1, spline_2)) %>%
    select(group, well, time_sec, fluorescence,
           zero_order, first_order, second_order,
           auc_1, slope, log_slope_till_max,
           max_time, max_time_2, time_end)

  return(transformed)
}

# input -------------------------------------------------------------------

#assign the input data and key file name to the variable
#choose a pair or add your own
#data files live in the /inst/extdata folder

#data_filename <- "evelien_test_expt6.txt"
#key_filename <- "key_file_template.xlsx"

# data_filename <- "20230306_5dpfoptimalisation.txt"
# key_filename <- "key_file_exp20230306_5dpfoptimalisation.xlsx"

# data_filename <- "20230313_3dpfwildtypeOCR.txt"
# key_filename <- "key_file_exp20230313_3dpfwildtypeOCR.xlsx"

data_filename <- "20230314_4dpfwildtypeOCR.txt"
key_filename <- "key_file_exp20230314_4dpfwildtypeOCR_AUC.xlsx"

# data_filename <- "20230315_5dpfwildtypeOCRfail.txt"
# key_filename <- "key_file_exp20230315_5dpfwildtypeOCRfail_AUC.xlsx"

# data_filename <- "20230315_5dpfwildtypeOCR2.txt"
# key_filename <- "key_file_exp20230315_5dpfwildtypeOCR2_AUC.xlsx"

data_filename <- "20230403_FCCP4dpf.txt"
key_filename <- "key_file_exp20230403_FCCP4dpf.xlsx"



# transform the input data ------------------------------------------------

transformed <-
  get_processsed_mXp(data_filename, key_filename) %>%
  get_transformed_mXp(max_cutoff = 8*60, degrees_of_freedom = 10)


df_test <- get_processsed_mXp(data_filename, key_filename)
#saveRDS(transformed, file = here::here("data", "transformed_4dpf_example.rds"))


# explore the "transformed" tibble ----------------------------------------

transformed %>%
  slice(1, .by = well) %>%
  summarise(count = sum(is.na(log_slope_till_max)))

view(transformed %>% slice(1, .by = well))

# export output -----------------------------------------------------------

#export log slopes
readr::write_delim(transformed %>%
                     group_by(well) %>%
                     slice(1)%>%
                     select(well, group, time_sec, log_slope_till_max),
                   here::here("output",
                              paste0("output_slopes_", str_sub(data_filename, end=-4), ".txt")))

#export total data
readr::write_delim(transformed,
                   here::here("output",
                              paste0("output_total_", str_sub(data_filename, end=-4), ".txt")))

# plot the data -----------------------------------------------------------

#plot raw data with fit smoothed line
transformed %>%
  #filter(well %in% c("C02", "B02", "B03", "C03")) %>%
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

#plot first order data
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

#plot second order data
transformed %>%
  ggplot(aes(x=time_sec, y = second_order, color = well, group = group))+
  geom_point()+
  ggrepel::geom_label_repel(data = . %>% filter(time_sec == max(time_sec)),
                            aes(label = well),
                            show.legend = FALSE)+
  colorspace::scale_colour_discrete_divergingx(palette = "Geyser", rev = FALSE)+
  labs(y = "change in dAU/dt (2nd der.)",
       x = "time (sec)") +
  facet_wrap(~group, scales = "fixed")+
  theme_bw(base_size = 24)

#plot log transformed data
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

#plot the log slopes
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






