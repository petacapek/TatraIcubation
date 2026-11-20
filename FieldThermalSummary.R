####################################################################################################
###                        Tatra Mountains soil cores incubation experiment                      ###
######                             background field thermal data                             #######
####################################################################################################
source('packages.R')
# ---- DATA ----
tatraT <- read.xlsx('SourceDataXLSX/Tatry_soil_temperature_winter_summary.xlsx', 8)
tatraT$Time <- convertToDateTime(tatraT$Time, tz = "Europe/Prague")
summary(tatraT)

# ---- winter conditions ----
season_of <- function(d) ifelse(month(d) >= 7, year(d), year(d) - 1)

# 0. wide -> long (all columns except datetime are localities)
long <- tatraT |>
  pivot_longer(-Time, names_to = "locality", values_to = "temp") |>
  mutate(date = date(Time), season = season_of(date))

# 1. daily means + validity, per locality
daily <- long |>
  group_by(locality, date) |>
  summarise(n_obs  = sum(!is.na(temp)),
            mean_t = mean(temp, na.rm = TRUE), .groups = "drop") |>
  group_by(locality) |>
  complete(date = seq(min(date), max(date), by = "day"),
           fill = list(n_obs = 0)) |>
  ungroup() |>
  mutate(valid  = n_obs >= 20,
         cold   = coalesce(valid & mean_t < 5, FALSE),
         season = season_of(date))

# 2. winter bounds per locality and season
winter_bounds <- function(date, cold, min_run = 5) {
  r      <- rle(cold)
  ends   <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1
  ok     <- which(r$values & r$lengths >= min_run)
  if (length(ok) == 0) return(tibble(start = as.Date(NA), end = as.Date(NA)))
  tibble(start = date[starts[min(ok)]],
         end   = date[ends[max(ok)]])    # use starts[max(ok)] for first day of last run
}

bounds <- daily |>
  arrange(locality, date) |>
  group_by(locality, season) |>
  reframe(winter_bounds(date, cold))
bounds

# 3. subset hourly data
winter <- long |>
  inner_join(bounds, by = c("locality", "season")) |>
  filter(date >= start, date <= end)

summary(winter)

# ---- temperature summaries for the winter period ----
winter |> group_by(locality) |> summarise(MeanT = mean(temp, na.rm = T),
                                          MedianT = median(temp, na.rm = T),
                                          MinT = min(temp, na.rm = T),
                                          MaxT = max(temp, na.rm = T))

# ---- number of freezing periods ----
# frozen days within each winter
frozen_days <- daily |>
  inner_join(bounds, by = c("locality", "season")) |>
  filter(date >= start, date <= end) |>
  arrange(locality, date) |>
  group_by(locality, season) |>
  mutate(frozen = coalesce(valid & mean_t < 0, FALSE),
         run_id = consecutive_id(frozen)) |>      # new id whenever frozen switches
  ungroup()

# one row per freezing period
periods <- frozen_days |>
  filter(frozen) |>
  group_by(locality, season, run_id) |>
  summarise(start    = min(date),
            end      = max(date),
            n_days   = n(),
            min_mean = min(mean_t), .groups = "drop")

# number of freezing periods per locality and winter
n_periods <- bounds |>
  left_join(count(periods, locality, season, name = "n_freezing"),
            by = c("locality", "season")) |>
  mutate(n_freezing = if_else(is.na(start), NA_integer_, coalesce(n_freezing, 0L)))
n_periods

# ================= Scree soils
screeT <- read.xlsx('SourceDataXLSX/Tatry_soil_temperature_winter_summary.xlsx', 9)
screeT$Time <- convertToDateTime(screeT$Time, tz = "Europe/Prague")
summary(screeT)

# ---- winter period ----
# 1. daily means + validity
dailyScree <- screeT |>
  mutate(date = date(Time)) |>
  group_by(date) |>
  summarise(n_obs  = sum(!is.na(temp)),
            mean_t = mean(temp, na.rm = TRUE), .groups = "drop") |>
  complete(date = seq(min(date), max(date), by = "day"),
           fill = list(n_obs = 0)) |>              # missing days become invalid
  mutate(valid  = n_obs >= 20,
         cold   = coalesce(valid & mean_t < 5, FALSE),
         season = season_of(date))

# 2. start/end of winter from runs of cold days
boundsScree <- dailyScree |>
  arrange(date) |>
  group_by(season) |>
  reframe(winter_bounds(date, cold))
boundsScree

# 3. subset the hourly data
winterScree <- screeT |>
  mutate(date = date(Time), season = season_of(date)) |>
  inner_join(boundsScree, by = "season") |>
  filter(date >= start, date <= end)

# ---- temperature summaries for the winter period ----
winterScree |> group_by(season) |> summarise(MeanT = mean(temp, na.rm = T),
                                          MedianT = median(temp, na.rm = T),
                                          MinT = min(temp, na.rm = T),
                                          MaxT = max(temp, na.rm = T))

# ---- number of freezing periods ----
# frozen days within each winter
frozen_days_scree <- dailyScree |>
  inner_join(boundsScree, by = c("season")) |>
  filter(date >= start, date <= end) |>
  arrange(date) |>
  group_by(season) |>
  mutate(frozen = coalesce(valid & mean_t < 0, FALSE),
         run_id = consecutive_id(frozen)) |>      # new id whenever frozen switches
  ungroup()

# one row per freezing period
periods_scree <- frozen_days_scree |>
  filter(frozen) |>
  group_by(season, run_id) |>
  summarise(start    = min(date),
            end      = max(date),
            n_days   = n(),
            min_mean = min(mean_t), .groups = "drop")

# number of freezing periods per locality and winter
n_periods_scree <- periods_scree |>
  left_join(count(periods_scree, season, name = "n_freezing"),
            by = c("season")) |>
  mutate(n_freezing = if_else(is.na(start), NA_integer_, coalesce(n_freezing, 0L)))
n_periods_scree
