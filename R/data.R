## documentation block for `sim_outbreak` ----

#' The `sim_outbreak` dataset
#'
#' This is a generated data set of an outbreak based on a branching process
#' model with fixed parameters. There is a both the full infection line list
#' as it happened in the simulation and the subset of this data that might have
#' been observed in a contact tracing exercise where identification is done
#' through symptoms.
#'
#' ## `sim_outbreak` is a named list with 3 items
#'
#' \tabular{lll}{
#'     \bold{Item} \tab \bold{Type} \tab \bold{Description} \cr
#'     `parameters` \tab `list[dbl]` \tab the ground truth of the simulation parameters \cr
#'     `outbreak_truth` \tab `df[outbreak_truth]`* \tab the full infection line list \cr
#'     `contact_tracing` \tab `df[contact_tracing]`* \tab the "observed" contact tracing subset \cr
#'  }
#'
#'
#' ## `df[outbreak_truth]` dataframe with 663 rows and 11 columns
#'
#' The simulation details
#'
#' \tabular{lll}{
#'     \bold{Column} \tab \bold{Type} \tab \bold{Description} \cr
#'     `time` \tab `dbl` \tab The true time of infection \cr
#'     `id` \tab `int` \tab Person unique id \cr
#'     `generation_interval` \tab `dbl` \tab The time since infector's infection \cr
#'     `infector` \tab `int` \tab The unique id of the infector \cr
#'     `generation` \tab `dbl` \tab Which generation is this infection since the simulation start \cr
#'     `symptom` \tab `lgl` \tab Did this person experience symptoms \cr
#'     `symptom_delay` \tab `dbl` \tab How long after infection were their symptoms?\cr
#'     `symptom_time` \tab `dbl` \tab When? (from the simulation start)\cr
#'     `observation` \tab `lgl` \tab Was this person detected (only if symptoms) \cr
#'     `observation_delay` \tab `dbl` \tab How long after symptoms were they observed? \cr
#'     `observation_time` \tab `dbl` \tab When? (from the simulation start) \cr
#'  }
#'
#'
#' ## `df[contact_tracing]` dataframe with 663 rows and 4 columns
#'
#' A minimal set of data that might be collected in a contact tracing
#' exercise.
#'
#' \tabular{lll}{
#'     \bold{Column} \tab \bold{Type} \tab \bold{Description} \cr
#'     `id` \tab `int` \tab Unique person id \cr
#'     `contact_id` \tab `int` \tab Unique id of infectious contact \cr
#'     `onset_time` \tab `int` \tab Time of symptom onset (from the simulation start) \cr
#'     `obs_time` \tab `int` \tab Time of observation (from the simulation start) \cr
#'  }
#'
#' @source <https://ai4ci.github.io/ggoutbreak/articles/simulation-test-models.html#line-list-simulations>
#'
"sim_outbreak"

## end of documentation block for `sim_outbreak`
