
library(tidyverse)
library(NHSRdatasets)
library(NHSRplotthedots)
library(gt)
library(gtExtras)
library(purrr)
library(glue)
library(cli)

cli_inform('Loading data')
data <- read.csv('MHSdata.csv')


cli_alert_success('Loading data - complete')

cli_alert ('Wrangling the datas')

# add metric ids to metrics
# remove duplicates
# reformat dates to dates
dat <- data |>
  mutate(metric = paste0(OrgCode, MetricName),
         metric_id = case_when(MetricName == "Bed occupancy classed as clinically ready for discharge (weekly average patient count, acute)" ~ 1,
                               MetricName == "Patients discharged with length of stay of LOS 21+ days (%)"   ~ 2,
                               MetricName == "Patients discharged with length of stay of LOS 21+ Days (count)"   ~ 3,
                               MetricName == "Bed occupancy classed as clinically ready for discharge (%, acute)"   ~ 4,
                               MetricName == "Virtual ward capacity" ~ 5,
                               MetricName == "Virtual ward occupancy, %" ~ 6,
                               MetricName == "The number of patients with a length of stay of 21 days and over who do not meet the criteria to reside that day" ~ 7,
                               MetricName == "Average number of cases per 4 hour session" ~ 8,
                               MetricName == "PIFU utilisation rate"  ~ 9,
                               MetricName == "Utilisation rate for 'all types of specialist advice' in all specialties" ~ 10,
                               MetricName == "BADS Emergency surgery: Day case and outpatient % of total procedures (inpatient, day case and outpatient) (3mths to period end)" ~ 11,
                               MetricName == "Ratio of follow up : first appointment" ~ 12,
                               MetricName == "Missed outpatient appointments (DNAs) rate" ~ 13,
                               .default = NA
                               ),
         imp = ImprovementDirection,
         ReportingDate = as.Date(ReportingDate, '%d/%m/%Y')) |>
  filter(!Domain %in% c("Daycase Sentinel",
                       "Specialty Overview" ,
                       "Operational and Clinical Productivity"))

# runs a check that all metrics have an id
if(sum(is.na(dat$metric_id)) > 0) 
  {cli_alert("Warning: Metric name out of bounds - check list")}

cli_alert_success('All metrics sucessfully given ID')

cli_alert('Defining SPC function')

spc_icons <- function (df, met, assu_or_var) {

    # filter data to just the metric
  dat <- df |>
    filter (metric == met)
  
  # pull the improvement direction
  imp <- dat$imp[1]
  
  # run the spc dataframe
  spc_dat <- ptd_spc(dat,
                     value_field = Value,
                     date_field = ReportingDate,
                     #target = tg,
                     improvement_direction = 'increase')
  
  # find the latest value, upl, lpl, and latest special point type
  latest_val <- spc_dat$y[spc_dat$x == max(spc_dat$x, na.rm = TRUE)]
  upl <- spc_dat$upl[1]
  lpl <- spc_dat$lpl[1]
  latest_pt <- spc_dat$point_type[spc_dat$x == max(spc_dat$x, na.rm = TRUE)]
  
  # calculate which icon to use
  
  # variation icons
  icon_var <- case_when(latest_pt == 'special_cause_improvement' ~ 'Improved',  
                        latest_pt == 'special_cause_concern'  ~ 'Deteriorated', 
                        .default = 'Maintained') # common case

  # narr_varr = case_when (icon_var == "CCV" ~ "The measure is within common cause variation, with no significant change.", 
  #                        icon_var == "SCH" ~ "There is evidence of special cause variation of a concerning nature.", 
  #                        icon_var == "SCL" ~ "There is evidence of special cause variation of a concerning nature.",
  #                        icon_var == "SIH" ~ "There is evidence of special cause variation of a improving nature.", 
  #                        icon_var == "SIL" ~ "There is evidence of special cause variation of a improving nature.",
  #                        icon_var == "BLANK" ~ " ",
  #                        .default =  "Error - please check")
  
  # return the icon or narrative
  res <- case_when (#assu_or_var == 'assurance' ~ icon_assu, 
                    assu_or_var == 'variation' ~ icon_var,
                    .default = '')  #glue('<span style="font-size:0.7em;">{narr_assur} <br/> {narr_varr}</span>'))
  
  res
}

#spc_icons(dat, 'RQM', 'assurance')
#spc_icons(dat, 'RQM','variation')

# that's the functions set up, now want to run functions across each of the metrics

dat_f <- dat |>
  filter (ReportingDate == max(ReportingDate),
          .by = metric)

# make a list of the metrics 
metrics_list <- unique(dat_f$metric)

cli_alert('Running SPC function across all metrics')

# calculate variation icon for each metric
spc_ic_var  <- map(.x = metrics_list, 
                   .f = ~spc_icons(df = dat, 
                                   met = .x,
                                   assu_or_var = 'variation'))

cli_alert_success('Running SPC function across all metrics - complete')

cli_alert('Mangling data into table format')

# create a dataframe of just the latest result
dat_f <- dat |>
  filter (ReportingDate == max(ReportingDate),
          .by = metric) |>
  mutate(ic_var = spc_ic_var,
         chk_id = metrics_list) |>
  select(OrgName,
         metric_id,
         ic_var,
         StandardizedQuartile)

grouped_data_list_tidy <- dat_f %>% 
  mutate(lst = list(unique(metric_id)),
         .by = c(OrgName, 
                 ic_var, 
                 StandardizedQuartile
                 )) |>
  mutate(tes = paste(unique(sort(unlist(lst))), collapse = ", "),
         .by = c(OrgName, 
                 ic_var, 
                 StandardizedQuartile,
                 OrgName
         )) |>
  select(-lst,
         -metric_id)


dat_fill <- data.frame(ic_var_t = c('Improved', 'Deteriorated', 'Maintained'))

dat_tab <- cross_join(grouped_data_list_tidy, dat_fill) |>
  filter(StandardizedQuartile  != 'Undefined') |>
  mutate(metric_id = if_else(ic_var_t == ic_var, tes, NA),
         StandardizedQuartile = if_else(ic_var_t == ic_var , StandardizedQuartile, NA),
         ic_var_t = factor(ic_var_t, levels = c('Improved',
                                                   'Maintained',
                                                   'Deteriorated'))) |>
  select(-tes,
         -ic_var,
         StandardizedQuartile) |>
  unique() |>
  pivot_wider(names_from = StandardizedQuartile, 
              names_prefix = 'quart', 
              values_from = metric_id) |>
  select(OrgName,
         ic_var_t,
         quart4,
         quart3,
         quart2,
         quart1) |>
  arrange(OrgName,
          ic_var_t)

cli_alert_success('Mangling data into table format - complete')

cli_alert('Making pretty table')

dat_tab |> gt(groupname_col = 'OrgName') |>
  tab_style(
    style = list(
      cell_fill(color = '#ED8B00', alpha = 0.5)
      ),
      locations = cells_body(
        columns = c('quart4', 'quart3'),
        rows = ic_var_t %in% c('Maintained', 'Deteriorated') 
      )
    ) |>
  tab_style(
    style = list(
      cell_fill(color = '#768692', alpha = 0.5)
    ),
    locations = cells_body(
      columns = c('quart4', 'quart3'),
      rows = ic_var_t %in% c('Improved') 
    )
  ) |>
  tab_style(
    style = list(
      cell_fill(color = '#768692', alpha = 0.5)
    ),
    locations = cells_body(
      columns = c('quart2', 'quart1'),
      rows = ic_var_t %in% c('Deteriorated') 
    )
  ) |>
  tab_style(
    style = list(
      cell_fill(color = '#00A9CE', alpha = 0.5)
    ),
    locations = cells_body(
      columns = c('quart2', 'quart1'),
      rows = ic_var_t %in% c('Improved', 'Maintained') 
    )
  ) |>
  tab_style(
    style = list(
      cell_fill(color = '#0072CE', alpha = 0.5)
    ),
    locations = cells_body(
      columns = c('ic_var_t'),
      rows = ic_var_t %in% c('Improved') 
    )
  ) |>
  tab_style(
    style = list(
      cell_fill(color = '#768692', alpha = 0.5)
    ),
    locations = cells_body(
      columns = c('ic_var_t'),
      rows = ic_var_t %in% c('Maintained') 
    )
  ) |>
  tab_style(
    style = list(
      cell_fill(color = '#DA291C', alpha = 0.5)
    ),
    locations = cells_body(
      columns = c('ic_var_t'),
      rows = ic_var_t %in% c('Deteriorated') 
    )
  ) |>
  tab_style(
    style = list(
      cell_fill(color = '#DA291C', alpha = 0.5)
    ),
    locations = cells_column_labels(
      columns = c('quart4')
    )
  ) |>
  tab_style(
    style = list(
      cell_fill(color = '#ED8B00', alpha = 0.5)
    ),
    locations = cells_column_labels(
      columns = c('quart3')
    )
  ) |>
  tab_style(
    style = list(
      cell_fill(color = '#FAE100', alpha = 0.5)
    ),
    locations = cells_column_labels(
      columns = c('quart2')
    )
  ) |>
  tab_style(
    style = list(
      cell_fill(color = '#009639', alpha = 0.5)
    ),
    locations = cells_column_labels(
      columns = c('quart1')
    )
  ) |>
  tab_style(
    style = list(
      cell_text(weight = 'bolder',
                style = "italic")
    ),
    locations = cells_column_labels(
      columns = c('ic_var_t')
    )
  ) |>
  sub_missing(missing_text = '-') |>
  cols_align(align = 'center',
             columns = starts_with('quart')) |>
  tab_spanner(label = 'QUARTILE',
              columns = c(quart4,
                          quart3,
                          quart2,
                          quart1)) |>
  cols_label(ic_var_t = ' ',
             quart4 = '4',
             quart3 = '3',
             quart2 = '2',
             quart1 = '1') |>
  tab_options(data_row.padding = px(30)) |>
  tab_style(
    style = cell_borders(
      sides = c("all"),
      color = "white",
      weight = px(1.5),
      style = "solid"
    ),
    locations = cells_body()
  )
  
cli_alert_success('Making pretty table - complete')

cli_alert_success('Ta da!')
  

# helper function to eyeball results
plot_spc_test <- function(org = 'Dorset County Hospital NHS Foundation Trust', id) {
  p_data <- dat |>
    filter(OrgName == org,
           metric_id == id)
  
  imp <- unique(p_data$imp)
  
  ptd_spc(p_data,
          value_field = Value,
          date_field = ReportingDate,
          improvement_direction = imp)

}

# check function
# plot_spc_test('Gloucestershire Hospitals NHS Foundation Trust', id = 5)


