## code to prepare `DATASET` dataset goes here

sim_params = list(
  # A short generation time
  mean_gt = 4,
  sd_gt = 2,
  R0 = 2,
  I0 = 10,
  # Add a longish and very variable delay to symptoms
  p_symptomatic = 0.3,
  mean_onset = 7,
  sd_onset = 5,
  # and a slightly shorter delay to observation:
  # Only symptomatic cases are observed
  p_detected_given_symptoms = 0.7,
  mean_obs = 5,
  sd_obs = 3,
  # Observation cutoff:
  T_obs = 40
)

# Run simulation ----

sim_ip = ggoutbreak::make_gamma_ip(
  median_of_mean = sim_params$mean_gt,
  median_of_sd = sim_params$sd_gt
)

sim_params$r0 = ggoutbreak::inv_wallinga_lipsitch(sim_params$R0, sim_ip)

truth = ggoutbreak::sim_branching_process(
  fn_Rt = ~ sim_params$R0,
  fn_ip = ~sim_ip,
  fn_imports = \(t) ifelse(t == 1, sim_params$I0, 0),
  max_time = 40,
  seed = 123
)

delayed = truth %>%
  ggoutbreak::sim_delay(
    p_fn = ~ sim_params$p_symptomatic,
    delay_fn = ~ ggoutbreak::rgamma2(
      .x,
      sim_params$mean_onset,
      sim_params$sd_onset
    ),
    input = "time",
    output = "symptom",
    seed = 123
  ) %>%
  ggoutbreak::sim_delay(
    p_fn = \(t, symptom) {
      ifelse(symptom, sim_params$p_detected_given_symptoms, 0)
    },
    delay_fn = ~ ggoutbreak::rgamma2(
      .x,
      sim_params$mean_obs,
      sim_params$sd_obs
    ),
    input = "symptom_time",
    output = "observation",
    seed = 123
  ) %>%
  dplyr::mutate(
    dplyr::across(dplyr::where(ggoutbreak::is.time_period), as.numeric)
  ) %>%
  dplyr::filter(observation_time < sim_params$T)


observed = delayed %>%
  dplyr::left_join(
    delayed %>% dplyr::transmute(infector = id, contact_id = id),
    by = "infector"
  ) %>%
  dplyr::transmute(
    id,
    contact_id = contact_id,
    onset_time = as.integer(floor(symptom_time)),
    obs_time = as.integer(floor(observation_time))
  ) %>%
  glimpse()

sim_outbreak = list(
  parameters = sim_params,
  outbreak_truth = delayed,
  contact_tracing = observed
)

usethis::use_data(sim_outbreak, overwrite = TRUE)
