# The `sim_outbreak` dataset

This is a generated data set of an outbreak based on a branching process
model with fixed parameters. There is a both the full infection line
list as it happened in the simulation and the subset of this data that
might have been observed in a contact tracing exercise where
identification is done through symptoms.

## Usage

``` r
sim_outbreak
```

## Format

An object of class `list` of length 3.

## Source

<https://ai4ci.github.io/ggoutbreak/articles/simulation-test-models.html#line-list-simulations>

## Details

### `sim_outbreak` is a named list with 3 items

|  |  |  |
|----|----|----|
| **Item** | **Type** | **Description** |
| `parameters` | `list[dbl]` | the ground truth of the simulation parameters |
| `outbreak_truth` | `df[outbreak_truth]`\* | the full infection line list |
| `contact_tracing` | `df[contact_tracing]`\* | the "observed" contact tracing subset |

### `df[outbreak_truth]` dataframe with 663 rows and 11 columns

The simulation details

|  |  |  |
|----|----|----|
| **Column** | **Type** | **Description** |
| `time` | `dbl` | The true time of infection |
| `id` | `int` | Person unique id |
| `generation_interval` | `dbl` | The time since infector's infection |
| `infector` | `int` | The unique id of the infector |
| `generation` | `dbl` | Which generation is this infection since the simulation start |
| `symptom` | `lgl` | Did this person experience symptoms |
| `symptom_delay` | `dbl` | How long after infection were their symptoms? |
| `symptom_time` | `dbl` | When? (from the simulation start) |
| `observation` | `lgl` | Was this person detected (only if symptoms) |
| `observation_delay` | `dbl` | How long after symptoms were they observed? |
| `observation_time` | `dbl` | When? (from the simulation start) |

### `df[contact_tracing]` dataframe with 663 rows and 4 columns

A minimal set of data that might be collected in a contact tracing
exercise.

|              |          |                                                   |
|--------------|----------|---------------------------------------------------|
| **Column**   | **Type** | **Description**                                   |
| `id`         | `int`    | Unique person id                                  |
| `contact_id` | `int`    | Unique id of infectious contact                   |
| `onset_time` | `int`    | Time of symptom onset (from the simulation start) |
| `obs_time`   | `int`    | Time of observation (from the simulation start)   |
