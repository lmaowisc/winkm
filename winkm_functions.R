

# Utilities ---------------------------------------------------------------
# A simple 'lag' function that defaults to 1 for the first entry
lag_surv <- function(x) {
  dplyr::lag(x, default = 1)
}

# A simple interpolation function that uses approx()
# rule=2 prevents NA at extremes
interpol <- function(x, y, xout) {
  approx(x = x, y = y, xout = xout, rule = 2)$y
}


# Reading and Pivoting Digitized KM Data -------------------------------
# Suppose you have digitized data in a wide format with columns such as:
#   "time_placebo", "surv_placebo", "time_camr", "surv_camr"
# You can use the following function to pivot them into a long format
# Time and survival columns are expected to be named in the pattern:
#   "time_group", "surv_group"
# It automatically adds (time=0, surv=1) for each group if zero_time = TRUE.
prepare_km_data <- function(km_data,
                            time_cols,
                            surv_cols,
                            group_labels = NULL,
                            ref = NULL,    # specify the reference (control) group
                            zero_time = TRUE) {
  # 1) Pivot from wide to long using a name pattern:
  #    e.g. time_placebo, surv_placebo -> .value = "time"/"surv", group = "placebo"
  df_long <- km_data |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(c(time_cols, surv_cols)),
      names_to = c(".value", "group"),
      names_pattern = "(time|surv)_(.*)"
    ) |>
    dplyr::arrange(group, time) |>
    dplyr::mutate(
      surv = pmax(pmin(surv, 1), 0)  # Ensure 0 <= surv <= 1
      ) |>
    group_by(group) |>
    mutate(surv = cummin(surv)) # Cumulative min to ensure survival is non-increasing





  # 2) Optionally add a row (time=0, surv=1) to represent 100% survival at t=0
  if (isTRUE(zero_time)) {
    unique_groups <- unique(df_long$group)
    zero_rows <- tibble::tibble(
      group = unique_groups,
      time  = 0,
      surv  = 1
    )
    df_long <- dplyr::bind_rows(df_long, zero_rows) |>
      dplyr::arrange(group, time)
  }

  # 3) If 'ref' is specified and is a valid group, reorder factor levels
  #    so that the reference group comes first.
  if (!is.null(ref) && ref %in% df_long$group) {
    all_groups <- unique(df_long$group)
    # Place 'ref' first, keep any others in their existing order
    new_levels <- c(ref, setdiff(all_groups, ref))
    df_long$group <- factor(df_long$group, levels = new_levels)
  }

  # 4) If user provides group_labels, rename groups according to the *current* levels
  #    (which already puts 'ref' first if ref was given).
  if (!is.null(group_labels) && length(group_labels) == length(levels(df_long$group))) {
    df_long$group <- factor(
      df_long$group,
      levels = levels(df_long$group),   # preserve the newly ordered levels
      labels = group_labels
    )
  }

  return(df_long)
}



# Merging and Interpolating Multiple Endpoints ----------------------------
# Merge OS and PFS data into one table so that each has the same set of time points
merge_endpoints <- function(os_data_long, pfs_data_long, extra_times = numeric()) {
  # Gather all relevant time points
  times_all <- sort(unique(c(os_data_long$time, pfs_data_long$time, extra_times)))

  # Interpolate OS at the combined time grid
  os_wide <- os_data_long |>
    dplyr::group_by(group) |>
    dplyr::reframe(
      os   = interpol(time, surv, times_all),
      time = times_all
    )

  # Interpolate PFS at the combined time grid
  pfs_wide <- pfs_data_long |>
    dplyr::group_by(group) |>
    dplyr::reframe(
      pfs  = interpol(time, surv, times_all),
      time = times_all
    )

  # Join OS and PFS side by side
  df_combined <- dplyr::full_join(os_wide, pfs_wide, by = c("group", "time")) |>
    dplyr::arrange(group, time) |>
    dplyr::select(group, time, os, pfs)

  return(df_combined)
}




# Computing Incremental Quantities (dos, dpfs, etc.) ----------------------
# Once you have a combined dataset with os and pfs
# you can compute the incremental decrements:
compute_increments <- function(df_km) {
  df_km |>
    dplyr::group_by(group) |>
    dplyr::mutate(
      os_p   = lag_surv(os),
      dos    = pmax(os_p - os, 0),
      pfs_p  = lag_surv(pfs),
      dpfs   = pmax(pfs_p - pfs, 0),
      dLpfs  = dplyr::if_else(pfs_p > 0, dpfs / pfs_p, 0)
    ) |>
    dplyr::ungroup()
}


###############################################################################
# FUNCTION: compute_followup()
#
# Purpose:
#   - Take a risk table in "wide" format (multiple group columns) and compute
#     total follow-up time by endpoint and group.
#
# Arguments:
#   risk_table   - Data frame/tibble with columns specifying time, endpoint,
#                  and at-risk numbers for each group.
#   time_col     - Name of the column holding time points (default "time").
#   endpoint_col - Name of the column holding endpoint labels (default "endpoint").
#   group_cols   - Character vector specifying the wide-format group columns
#                  (default c("trt", "ctr")).
#   group_labels - Optional character vector to rename the pivoted groups.
#                  Length must match length(group_cols).
#
# Returns:
#   A tibble/data frame with columns: <endpoint> <group> total_followup
###############################################################################
compute_followup <- function(risk_table,
                             time_col     = "time",
                             endpoint_col = "endpoint",
                             group_cols   = c("trt", "ctr"),
                             group_labels = NULL) {
  # 1) Pivot from wide to long: group_cols -> single "group" column
  df_long <- risk_table |>
    tidyr::pivot_longer(
      cols      = dplyr::all_of(group_cols),
      names_to  = "group",
      values_to = "n_at_risk"
    )

  # 2) Optionally rename (factor) the group levels if group_labels provided
  if (!is.null(group_labels) && length(group_labels) == length(unique(df_long$group))) {
    df_long$group <- factor(
      df_long$group,
      levels = group_cols,
      labels = group_labels
    )
  }

  # 3) For each endpoint & group, sort by time, then compute person-time
  df_followup <- df_long |>
    # Group by the user-specified endpoint_col and the pivoted "group"
    dplyr::group_by(
      dplyr::across(dplyr::all_of(endpoint_col)),
      group
    ) |>
    # Sort within each group by the user-specified time column
    dplyr::arrange(
      dplyr::across(dplyr::all_of(time_col)),
      .by_group = TRUE
    ) |>
    dplyr::mutate(
      # Next time & at-risk in sequence
      t_next = dplyr::lead(.data[[time_col]], default = NA_real_),
      n_next = dplyr::lead(n_at_risk, default = NA_real_),
      # Number lost in [time_i, time_{i+1}]
      lost   = n_at_risk - n_next,
      # Midpoint of interval
      midpoint = (.data[[time_col]] + t_next) / 2,
      # Interval-specific follow-up
      interval_fu = lost * midpoint
    ) |>
    # 4) Summarize total follow-up for each group
    dplyr::summarize(
      total_followup = sum(interval_fu, na.rm = TRUE),
      .groups = "drop"
    )

  # Return one row per (endpoint, group) with total_followup
  df_followup
}


###############################################################################
# FUNCTION: compute_theta
#
# Purpose:
#   - Computes two "theta" parameters: one for the treatment arm (theta_trt)
#     and one for the reference/control arm (theta_ctr).
#   - Uses event counts and total follow-up times for OS and PFS to derive:
#       rD  = ND / LD   (rate of death events)
#       rs  = Ns / Ls   (rate of composite events: death or PD)
#       rE  = NP / Ls   (rate of progression events)
#
#     Then applies the formula:
#       theta = max( log(1 - rE / rs) / log(rD / rs), 1 )
#
# Arguments:
#   event_nums_trt - Named numeric vector for the treatment arm with keys:
#       * "ND" = number of OS (death) events
#       * "Ns" = number of composite events (death or PD)
#       * "NP" = number of progression events
#
#   event_nums_ctr - Named numeric vector for the control arm with the same keys.
#
#   fl_times       - Data frame/tibble containing:
#       * "group"          : group label (treatment or control)
#       * "endpoint"       : "os" or "pfs"
#       * "total_followup" : total person-time for that (group, endpoint)
#
#   ref            - Character string to identify the reference group (control).
#                    If NULL or not found in fl_times, user must ensure
#                    data structure aligns with what is needed.
#
# Returns:
#   A two-element numeric vector with names:
#       * theta_trt
#       * theta_ctr
#   Each value is constrained to be at least 1.
###############################################################################
compute_theta <- function(
    event_nums_trt,
    event_nums_ctr,
    fl_times,
    ref = NULL
) {

  # If ref is NULL or not found in fl_times$group, default to the first unique group.
  if (is.null(ref) || !ref %in% fl_times$group) {
    first_grp <- unique(fl_times$group)[1]
    message(
      "No valid reference group specified; defaulting to the first group: ",
      first_grp
    )
    ref <- first_grp
  }
  #---------------------------------------------------------------------------
  # 1) Extract total follow-up for OS and PFS, splitting by group == ref or not
  #---------------------------------------------------------------------------
  LD_trt <- fl_times$total_followup[
    fl_times$group != ref & fl_times$endpoint == "os"
  ]
  LD_ctr <- fl_times$total_followup[
    fl_times$group == ref & fl_times$endpoint == "os"
  ]
  Ls_trt <- fl_times$total_followup[
    fl_times$group != ref & fl_times$endpoint == "pfs"
  ]
  Ls_ctr <- fl_times$total_followup[
    fl_times$group == ref & fl_times$endpoint == "pfs"
  ]
  #---------------------------------------------------------------------------
  # 2) Compute event rates for the treatment arm:
  #      rD = ND / LD   (death event rate)
  #      rs = Ns / Ls   (composite event rate, death or PD)
  #      rE = NP / Ls   (progression event rate)
  #---------------------------------------------------------------------------
  rD_trt <- event_nums_trt["ND"] / LD_trt
  rs_trt <- event_nums_trt["Ns"] / Ls_trt
  rE_trt <- event_nums_trt["NP"] / Ls_trt
  #---------------------------------------------------------------------------
  # 3) Calculate theta for treatment arm:
  #      theta_trt = max( log(1 - rE_trt / rs_trt) / log(rD_trt / rs_trt), 1 )
  #---------------------------------------------------------------------------
  theta_trt <- max(
    log(1 - rE_trt / rs_trt) / log(rD_trt / rs_trt),
    1
  )
  #---------------------------------------------------------------------------
  # 4) Compute event rates for the control arm, same logic
  #---------------------------------------------------------------------------
  rD_ctr <- event_nums_ctr["ND"] / LD_ctr
  rs_ctr <- event_nums_ctr["Ns"] / Ls_ctr
  rE_ctr <- event_nums_ctr["NP"] / Ls_ctr

  theta_ctr <- max(
    log(1 - rE_ctr / rs_ctr) / log(rD_ctr / rs_ctr),
    1
  )
  #---------------------------------------------------------------------------
  # 5) Return a named vector with both thetas
  #---------------------------------------------------------------------------
  return(c(theta_trt = theta_trt, theta_ctr = theta_ctr))
}




# Function to compute the output factor for PD-based win/loss
compute_pd_fctr <- function(theta, Lpd, os, dos){
  # Relative reduction
  reduction <- ifelse(Lpd * os > 0,
                      (theta - 1) * dos / (Lpd * os), 0)
  # Output factor
  pmax(1 - reduction, 0)

}




###############################################################################
# MERGED FUNCTION: compute_win_loss
#
# Purpose:
#   - Provide a single interface to compute win/loss for OS and PD.
#   - If 'theta' is NULL -> run the original "compute_win_loss" logic (verbatim).
#   - If 'theta' is provided -> run the original "compute_win_loss1" logic (verbatim).
#
# Arguments:
#   df_inc         : A data frame/tibble containing (time, group, os, dLpfs, etc.).
#   event_nums_trt : Named numeric vector for the treatment arm
#                    with c("NP" = #progression_events, "Ns" = #subjects_at_risk).
#   event_nums_ctr : Same as above for the control arm.
#   ref            : Character specifying which 'group' is reference (control).
#                    If NULL or invalid, defaults to the first group found in df_inc.
#   theta          : If NULL, the simpler approach (original compute_win_loss) is used.
#                    If non-NULL, the advanced approach (original compute_win_loss1) applies.
#
# Returns:
#   A data frame containing:
#       time,
#       win_os, loss_os,   # OS-based win/loss
#       win_pd, loss_pd,   # PD-based win/loss
#       win, loss          # total win/loss = OS + PD
###############################################################################
compute_win_loss <- function(df_inc,
                             event_nums_trt,
                             event_nums_ctr,
                             ref = NULL,
                             theta = NULL) {


  ###########################################################################
  # Decide which code block to run based on 'theta':
  #   - If is.null(theta) --> Run "original compute_win_loss"
  #   - Else --> Run "original compute_win_loss1"
  ###########################################################################
  if (is.null(theta)) {
    ###########################################################################
    # BEGIN: ORIGINAL compute_win_loss CODE (verbatim)
    ###########################################################################

    # 1) Determine reference (control) group:
    #    - If ref is NULL or not in df_inc$group, use the first group found in df_inc
    if (is.null(ref) || !ref %in% df_inc$group) {
      message("No valid reference group specified; defaulting to the first group: ",
              unique(df_inc$group)[1])
      ref <- unique(df_inc$group)[1]
    }

    # 2) Extract progression rates from event numbers:
    #    event_nums_xxx = c("NP" = progression_events, "Ns" = total_at_risk)
    trt_pdr <- event_nums_trt["NP"] / event_nums_trt["Ns"]
    ctr_pdr <- event_nums_ctr["NP"] / event_nums_ctr["Ns"]

    # 3) OS part: pivot to (os_trt, dos_trt, os_ctr, dos_ctr)
    wl_os <- df_inc |>
      dplyr::select(group, time, os, dos) |>
      dplyr::mutate(
        # Mark 'ref' group as control (ctr); all others as 'trt'
        group = dplyr::if_else(group == ref, "ctr", "trt")
      ) |>
      tidyr::pivot_wider(
        id_cols     = "time",
        names_from  = "group",
        values_from = c(os, dos)
      ) |>
      dplyr::mutate(
        win_os  = cumsum(os_trt * dos_ctr),
        loss_os = cumsum(os_ctr * dos_trt)
      ) |>
      dplyr::select(time, win_os, loss_os)

    # 4) PD part: pivot to (os_trt, dLpfs_trt, os_ctr, dLpfs_ctr)
    wl_pd <- df_inc |>
      dplyr::select(group, time, os, dLpfs) |>
      dplyr::mutate(
        group = dplyr::if_else(group == ref, "ctr", "trt")
      ) |>
      tidyr::pivot_wider(
        id_cols     = "time",
        names_from  = "group",
        values_from = c(os, dLpfs)
      ) |>
      dplyr::mutate(
        # Weighted progression hazard
        dLpd_trt = dLpfs_trt * trt_pdr,
        dLpd_ctr = dLpfs_ctr * ctr_pdr,

        # Convert hazard increments to survival probabilities for PD
        pd_trt   = cumprod(1 - dLpd_trt),
        pd_ctr   = cumprod(1 - dLpd_ctr),

        pd_p_trt = lag_surv(pd_trt),
        pd_p_ctr = lag_surv(pd_ctr),

        dpd_trt  = pmax(pd_p_trt - pd_trt, 0),
        dpd_ctr  = pmax(pd_p_ctr - pd_ctr, 0),

        # Win / Loss on PD
        win_pd  = os_trt * os_ctr * cumsum(pd_trt * dpd_ctr),
        loss_pd = os_trt * os_ctr * cumsum(pd_ctr * dpd_trt)
      ) |>
      dplyr::select(time, win_pd, loss_pd)

    # 5) Merge OS and PD components, compute total win/loss
    df_wl <- dplyr::left_join(wl_os, wl_pd, by = "time") |>
      dplyr::mutate(
        win  = win_os + win_pd,
        loss = loss_os + loss_pd
      )

    ###########################################################################
    # END: ORIGINAL compute_win_loss CODE
    ###########################################################################

  } else {
    ###########################################################################
    # BEGIN: ORIGINAL compute_win_loss1 CODE (verbatim)
    ###########################################################################
    theta[is.infinite(theta)] <- 10^4


    # 1) Determine reference (control) group:
    #    - If ref is NULL or not in df_inc$group, use the first group found in df_inc
    if (is.null(ref) || !ref %in% df_inc$group) {
      message("No valid reference group specified; defaulting to the first group: ",
              unique(df_inc$group)[1])
      ref <- unique(df_inc$group)[1]
    }

    # 2) Extract progression rates from event numbers:
    #    event_nums_xxx = c("NP" = progression_events, "Ns" = total_at_risk)
    trt_pdr <- event_nums_trt["NP"] / event_nums_trt["Ns"]
    ctr_pdr <- event_nums_ctr["NP"] / event_nums_ctr["Ns"]

    # 3) OS part: pivot to (os_trt, dos_trt, os_ctr, dos_ctr)
    wl_os <- df_inc |>
      dplyr::select(group, time, os, dos) |>
      dplyr::mutate(
        # Mark 'ref' group as control (ctr); all others are 'trt'
        group = dplyr::if_else(group == ref, "ctr", "trt")
      ) |>
      tidyr::pivot_wider(
        id_cols     = "time",
        names_from  = "group",
        values_from = c(os, dos)
      ) |>
      dplyr::mutate(
        win_os  = cumsum(os_trt * dos_ctr),
        loss_os = cumsum(os_ctr * dos_trt)
      ) |>
      dplyr::select(time, win_os, loss_os)

    # 4) PD part: pivot to (os_trt, dLpfs_trt, os_ctr, dLpfs_ctr)
    #    Extract the theta values
    theta_trt <- theta["theta_trt"]
    theta_ctr <- theta["theta_ctr"]

    df_pd <- df_inc |>
      dplyr::select(group, time, os, pfs, dos, dLpfs) |>
      dplyr::mutate(
        group = dplyr::if_else(group == ref, "ctr", "trt")
      ) |>
      tidyr::pivot_wider(
        id_cols     = "time",
        names_from  = "group",
        values_from = c(os, dos, dLpfs, pfs)
      ) |>
      dplyr::mutate(
        dLpd_ctr = dLpfs_ctr * ctr_pdr,
        dLpd_trt = dLpfs_trt * trt_pdr,
        # Compute cumulative hazard
        Lpfs_ctr = cumsum(dLpfs_ctr),
        Lpfs_trt = cumsum(dLpfs_trt),
        pd_fctr_ctr = compute_pd_fctr(theta_ctr, Lpfs_ctr, os_ctr, dos_ctr),
        pd_fctr_trt = compute_pd_fctr(theta_trt, Lpfs_trt, os_trt, dos_trt),
        os2 = os_trt * os_ctr
      ) |>
      dplyr::select(time, os2, dLpd_ctr, dLpd_trt, pd_fctr_ctr, pd_fctr_trt, os_trt, os_ctr,
                    pfs_trt, pfs_ctr)

    # A helper function to compute PD-based win/loss for each time
    compute_pd_wl <- function(tau, os2, os_trt, os_ctr) {
      df_pd |>
        dplyr::select(-c(os2, os_trt, os_ctr)) |>
        dplyr::filter(time <= tau) |>
        dplyr::mutate(
          cum_fctr_ctr = rev(cumprod(rev(pd_fctr_ctr))),
          cum_fctr_trt = rev(cumprod(rev(pd_fctr_trt))),
          dLpd_ctr = dLpd_ctr * cum_fctr_ctr,
          dLpd_trt = dLpd_trt * cum_fctr_trt,
          pd_trt = cumprod(1 - dLpd_trt),
          pd_ctr = cumprod(1 - dLpd_ctr),
          pd_trt = pmin(pd_trt, pfs_trt / os_trt),
          pd_ctr = pmin(pd_ctr, pfs_ctr / os_ctr),
          pd_p_trt = lag_surv(pd_trt),
          pd_p_ctr = lag_surv(pd_ctr),
          dpd_trt  = pmax(pd_p_trt - pd_trt, 0),
          dpd_ctr  = pmax(pd_p_ctr - pd_ctr, 0)
        ) |>
        dplyr::summarize(
          win_pd  = os2 * sum(pd_trt * dpd_ctr),
          loss_pd = os2 * sum(pd_ctr * dpd_trt)
        )
    }

    wl_pd <- tibble::tibble(
      time = df_pd$time,
      os2  = df_pd$os2,
      os_trt = df_pd$os_trt,
      os_ctr = df_pd$os_ctr,
      wl   = purrr::pmap(list(time, os2, os_trt, os_ctr), compute_pd_wl)
    ) |>
      tidyr::unnest(wl)

    # 5) Merge OS and PD components, compute total win/loss
    df_wl <- dplyr::left_join(wl_os, wl_pd, by = "time") |>
      dplyr::mutate(
        win  = win_os + win_pd,
        loss = loss_os + loss_pd
      )

    ###########################################################################
    # END: ORIGINAL compute_win_loss1 CODE
    ###########################################################################
  }

  df_wl <- df_wl |> mutate(
    win_os = cummax(win_os),
    loss_os = cummax(loss_os),
    win = cummax(win),
    loss = cummax(loss)
  )

  # Finally, return df_wl from whichever block was executed
  return(df_wl)
}


###############################################################################
# FUNCTION: run_win_loss_workflow
#
# Purpose:
#   - Provide a single "workflow" for digitized KM analysis:
#       * Read + pivot OS & PFS data (prepare_km_data)
#       * Optionally merge OS & PFS => increments (merge_endpoints + compute_increments)
#       * Optionally read or use a risk table to compute follow-up => compute_theta
#       * Finally, compute win/loss (compute_win_loss) with or without theta
#
# Arguments:
#   os_file, pfs_file    : Paths to CSV for OS/PFS data (optional if os_data, pfs_data given).
#   risk_table_file      : Path to CSV for the risk table (only if compute_theta=TRUE).
#   risk_table_data      : Already-loaded data frame for the risk table (optional).
#
#   os_data, pfs_data    : Already-loaded data frames for OS/PFS (optional).
#
#   event_nums_trt, event_nums_ctr : Named numeric vectors for treatment/control arms,
#                                    e.g. c(ND=..., NP=..., Ns=...).
#
#   ref_group            : Name of the reference group, used in pivoting OS/PFS.
#   ref_group_labelled   : Name of the reference group once labeled
#                          (e.g., "Placebo" if original is "placebo").
#                          Passed to compute_win_loss & compute_theta as 'ref'.
#   group_labels         : (Optional) vector for labeling arms in OS/PFS pivot.
#   merge_os_pfs         : Logical; if TRUE, merges OS & PFS & computes increments.
#   compute_theta        : Logical; if TRUE, reads or uses a risk table to compute follow-up, then theta.
#
#   col_names_os, col_names_pfs : Column names for reading OS/PFS CSV
#                                 (e.g. c("time_placebo","surv_placebo","time_camr","surv_camr")).
#   skip_lines_os, skip_lines_pfs : Integers for # lines to skip in OS/PFS CSV.
#
#   group_cols_at_risk   : (Only if compute_theta=TRUE) wide-format group columns in risk table
#                          (e.g. c("placebo_chemo", "camr_chemo")).
#   group_labels_at_risk : (Optional) labels for pivoting the risk table columns.
#
#   ... : Passed to prepare_km_data (e.g. zero_time=TRUE).
#
# Returns:
#   A list with:
#       os_long, pfs_long : pivoted OS & PFS data
#       df_merged         : merged OS & PFS (if merge_os_pfs=TRUE)
#       df_increments     : increments from compute_increments (if merge_os_pfs=TRUE)
#       theta             : numeric vector c(theta_trt, theta_ctr) or NULL
#       df_wl             : final data frame from compute_win_loss()
###############################################################################
run_win_loss_workflow <- function(
    # File paths for OS/PFS (optional)
  os_file          = NULL,
  pfs_file         = NULL,

  # Risk table: either a file or an already-loaded data frame
  risk_table_file  = NULL,
  risk_table_data  = NULL,

  # Alternatively, already-loaded data frames for OS/PFS
  os_data          = NULL,
  pfs_data         = NULL,

  # Named event numbers
  event_nums_trt,
  event_nums_ctr,

  # Reference group & labeling
  ref_group,
  ref_group_labelled = NULL,  # used as 'ref' in compute_theta & compute_win_loss
  group_labels     = NULL,    # for pivoting OS/PFS
  merge_os_pfs     = TRUE,
  compute_theta    = FALSE,

  # CSV reading specifics
  col_names_os     = NULL,
  col_names_pfs    = NULL,
  skip_lines_os    = 0,
  skip_lines_pfs   = 0,

  # For the risk table step (if compute_theta=TRUE)
  group_cols_at_risk   = NULL, # must specify if compute_theta=TRUE
  group_labels_at_risk = NULL, # optional labels for at-risk pivot

  # Additional arguments for prepare_km_data()
  ...
) {
  #############################################################################
  # 1) Load OS & PFS data from files or use existing data frames
  #############################################################################
  if (!is.null(os_file)) {
    os_data <- readr::read_csv(
      file      = os_file,
      skip      = skip_lines_os,
      col_names = col_names_os
    )
  }
  if (!is.null(pfs_file)) {
    pfs_data <- readr::read_csv(
      file      = pfs_file,
      skip      = skip_lines_pfs,
      col_names = col_names_pfs
    )
  }

  #############################################################################
  # 2) Pivot OS & PFS to long format (prepare_km_data), if they exist
  #############################################################################
  os_long  <- NULL
  pfs_long <- NULL

  if (!is.null(os_data)) {
    time_cols_os <- if (!is.null(col_names_os)) col_names_os[grep("^time_", col_names_os)] else NULL
    surv_cols_os <- if (!is.null(col_names_os)) col_names_os[grep("^surv_", col_names_os)] else NULL

    os_long <- prepare_km_data(
      km_data      = os_data,
      time_cols    = time_cols_os,
      surv_cols    = surv_cols_os,
      ref          = ref_group,
      group_labels = group_labels,
      ...
    )
  }

  if (!is.null(pfs_data)) {
    time_cols_pfs <- if (!is.null(col_names_pfs)) col_names_pfs[grep("^time_", col_names_pfs)] else NULL
    surv_cols_pfs <- if (!is.null(col_names_pfs)) col_names_pfs[grep("^surv_", col_names_pfs)] else NULL

    pfs_long <- prepare_km_data(
      km_data      = pfs_data,
      time_cols    = time_cols_pfs,
      surv_cols    = surv_cols_pfs,
      ref          = ref_group,
      group_labels = group_labels,
      ...
    )
  }

  #############################################################################
  # 3) Merge OS & PFS if requested; compute increments
  #############################################################################
  df_merged     <- NULL
  df_increments <- NULL

  if (merge_os_pfs && !is.null(os_long) && !is.null(pfs_long)) {
    df_merged     <- merge_endpoints(os_long, pfs_long)
    df_increments <- compute_increments(df_merged)
  }

  #############################################################################
  # 4) If compute_theta=TRUE, read/use risk table & compute follow-up => theta
  #############################################################################
  theta_vec <- NULL
  if (isTRUE(compute_theta)) {
    if (is.null(group_cols_at_risk)) {
      stop("group_cols_at_risk must be specified if compute_theta=TRUE,
            to indicate the wide-format group columns in the risk table.")
    }

    # If user did NOT supply risk_table_data but gave risk_table_file, read it
    if (is.null(risk_table_data) && !is.null(risk_table_file)) {
      risk_table_data <- readr::read_csv(risk_table_file)
    }

    # If we still have no risk_table_data, error out
    if (is.null(risk_table_data)) {
      stop("No risk table data available. Please provide either risk_table_file or risk_table_data.")
    }

    # 4a) Compute total follow-up using compute_followup
    fl_times <- compute_followup(
      risk_table       = risk_table_data,
      time_col         = "time",
      endpoint_col     = "endpoint",
      group_cols       = group_cols_at_risk,
      group_labels     = group_labels_at_risk
      # By design, we do NOT pass 'ref' here unless you want to reorder factor levels
      # inside follow-up data. If needed, add ref=ref_group_labelled or ref=ref_group.
    )

    # 4b) Then compute theta
    theta_vec <- compute_theta(
      event_nums_trt = event_nums_trt,
      event_nums_ctr = event_nums_ctr,
      fl_times       = fl_times,
      ref            = ref_group_labelled  # e.g. "Placebo" if pivoting gave that label
    )
  }

  #############################################################################
  # 5) Final call to compute_win_loss() with or without theta
  #############################################################################
  # Decide which data to pass: typically df_increments if merged, else OS or PFS alone
  df_for_win_loss <- df_increments %||% os_long %||% pfs_long
  if (is.null(df_for_win_loss)) {
    stop("No data available for compute_win_loss. Provide or merge OS/PFS.")
  }

  df_wl <- compute_win_loss(
    df_inc         = df_for_win_loss,
    event_nums_trt = event_nums_trt,
    event_nums_ctr = event_nums_ctr,
    ref            = ref_group_labelled,
    theta          = theta_vec
  )

  #############################################################################
  # 6) Return everything as a list
  #############################################################################
  out <- list(
    os_long        = os_long,
    pfs_long       = pfs_long,
    df_merged      = df_merged,
    df_increments  = df_increments,
    theta          = theta_vec,
    df_wl          = df_wl
  )
  return(out)
}


###############################################################################
# FUNCTION: analyze_ipd_data
#
# Purpose:
#   - Perform a comprehensive analysis on a data frame containing:
#       (id, time, status, <treatment_col>).
#   - The user specifies:
#       * Which column holds the treatment-group variable (e.g., "rx").
#       * The control level (control_level) vs. treatment level (treat_level).
#   - Produces KM curves for OS & PFS, merges them, computes increments,
#     event counts, and optionally applies compute_theta() logic to re-run
#     compute_win_loss() with an advanced factor.
#
# Arguments:
#   df               : A data frame with columns at least (id, time, status) plus
#                      one additional column (treatment_col) to specify arm assignment.
#   treatment_col    : Character name of the column in `df` that denotes
#                      the treatment group (e.g., "rx").
#   control_level    : Character specifying which value in `df[[treatment_col]]`
#                      is considered the control arm.
#   treat_level      : Character specifying which value in `df[[treatment_col]]`
#                      is considered the treatment arm.
#   compute_theta    : Logical; if TRUE, we calculate total follow-up times
#                      (simple approach) and compute \(\theta\), then re-run
#                      compute_win_loss() with that factor.
#
# Details of the Analysis Steps:
#   1. Data Preparation: group by 'id', arrange by 'id' and 'time', then ungroup.
#   2. KM Estimation (OS & PFS):
#      - OS excludes rows with status == 1 (optional choice).
#      - PFS uses first row per 'id', measuring if status > 0.
#   3. Convert these survfit() results into tidy data frames (os_long, pfs_long).
#   4. Merge endpoints => df_km, compute increments => df_increments.
#   5. Compute event counts (ND, NP, Ns) for each arm.
#   6. Compute basic Win/Loss with compute_win_loss (no \(\theta\)).
#   7. If compute_theta=TRUE:
#        - Summarize total OS & PFS follow-up times for each arm,
#          call compute_theta() => \(\theta\).
#        - Re-run compute_win_loss() with \(\theta\).
#
# Returns:
#   A named list containing:
#       - df_clean        : the cleaned, arranged data
#       - km_os_fit,
#         km_pfs_fit      : survfit summaries for OS & PFS
#       - figs_km         : (optional) combined plot for OS & PFS (edit as needed)
#       - os_long,
#         pfs_long        : the tidy KM data for OS & PFS
#       - df_km           : merged OS & PFS
#       - df_increments   : the increments from compute_increments(df_km)
#       - event_nums_trt,
#         event_nums_ctr  : named numeric vectors with ND, NP, Ns
#       - df_wl_no_theta  : basic Win/Loss (no \(\theta\))
#       - theta           : \(\theta\) vector if compute_theta=TRUE, else NULL
#       - df_wl_theta     : Win/Loss with \(\theta\) if compute_theta=TRUE, else NULL
###############################################################################
analyze_ipd_data <- function(df,
                             treatment_col    = "rx",
                             control_level    = "Control",
                             treat_level      = "Lev+5FU",
                             compute_theta    = TRUE) {
  #--------------------------------------------------------------------------
  # 1) PREPARE THE DATA
  #--------------------------------------------------------------------------
  # Group by 'id', arrange by (id, time), then ungroup
  df_clean <- df %>%
    dplyr::group_by(id) %>%
    dplyr::arrange(id, time) %>%
    dplyr::ungroup()

  #--------------------------------------------------------------------------
  # 2) KM ESTIMATES FOR OS & PFS
  #--------------------------------------------------------------------------
  # OS excludes rows with status == 1 (optional design)
  # We interpret status>0 as an event, ignoring status==1 rows
  os_formula <- stats::as.formula(paste0(
    "Surv(time, status > 0) ~ ", treatment_col
  ))

  km_os_fit <- summary(
    survival::survfit(
      formula = os_formula,
      data = df_clean %>% dplyr::filter(status != 1)
    )
  )

  # PFS uses only the first row per id
  pfs_formula <- stats::as.formula(paste0(
    "Surv(time, status > 0) ~ ", treatment_col
  ))
  # group_by(id) => slice(1)
  km_pfs_fit <- summary(
    survival::survfit(
      formula = pfs_formula,
      data = df_clean %>%
        dplyr::group_by(id) %>%
        dplyr::slice(1) %>%
        dplyr::ungroup()
    )
  )

  #--------------------------------------------------------------------------
  # 3) CONVERT survfit OUTPUT => TIDY (os_long, pfs_long)
  #--------------------------------------------------------------------------
  os_long <- tibble::tibble(
    group = stringr::str_remove(km_os_fit$strata, paste0(treatment_col, "=")),
    time  = km_os_fit$time,
    surv  = km_os_fit$surv
  )
  pfs_long <- tibble::tibble(
    group = stringr::str_remove(km_pfs_fit$strata, paste0(treatment_col, "=")),
    time  = km_pfs_fit$time,
    surv  = km_pfs_fit$surv
  )

  #--------------------------------------------------------------------------
  # 4) MERGE OS/PFS => df_km; COMPUTE INCREMENTS => df_increments
  #--------------------------------------------------------------------------
  df_km         <- merge_endpoints(os_long, pfs_long)
  df_increments <- compute_increments(df_km)

  #--------------------------------------------------------------------------
  # 5) COMPUTE EVENT NUMBERS (ND=2, NP=1, Ns= at risk for PD)
  #--------------------------------------------------------------------------
  # We'll interpret 'control_level' & 'treat_level' as the relevant groups
  # for ND, NP, Ns. Adjust logic as needed.
  df_NDP <- df_clean %>%
    dplyr::count(
      !!rlang::sym(treatment_col),
      status
    ) %>%
    tidyr::pivot_wider(
      names_from  = !!rlang::sym(treatment_col),
      values_from = "n",
      values_fill = 0
    )

  # Next, for PFS at-risk (# with status>0?), keep first row per id
  df_tfe <- df_clean %>%
    dplyr::group_by(id) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::count(status > 0, !!rlang::sym(treatment_col)) %>%
    dplyr::filter(`status > 0` == TRUE)

  # Extract counts from pivoted data
  event_nums_trt <- c(
    ND = df_NDP[[treat_level]][df_NDP$status == 2],
    NP = df_NDP[[treat_level]][df_NDP$status == 1],
    Ns = df_tfe$n[df_tfe[[treatment_col]] == treat_level]
  )
  event_nums_ctr <- c(
    ND = df_NDP[[control_level]][df_NDP$status == 2],
    NP = df_NDP[[control_level]][df_NDP$status == 1],
    Ns = df_tfe$n[df_tfe[[treatment_col]] == control_level]
  )

  #--------------------------------------------------------------------------
  # 6) BASIC WIN/LOSS (NO THETA)
  #--------------------------------------------------------------------------
  df_wl_no_theta <- compute_win_loss(
    df_inc         = df_increments,
    event_nums_trt = event_nums_trt,
    event_nums_ctr = event_nums_ctr,
    ref            = control_level
  )

  #--------------------------------------------------------------------------
  # 7) OPTIONAL: COMPUTE THETA => df_wl_theta
  #--------------------------------------------------------------------------
  df_wl_theta <- NULL
  thetas      <- NULL

  if (isTRUE(compute_theta)) {
    # Summarize total OS & PFS follow-up times (simplistic approach)
    fl_times <- df_clean %>%
      dplyr::group_by(id, !!rlang::sym(treatment_col)) %>%
      dplyr::summarize(
        os  = max(time),
        pfs = min(time),
        .groups = "drop"
      ) %>%
      dplyr::group_by(!!rlang::sym(treatment_col)) %>%
      dplyr::summarize(
        dplyr::across(c(os, pfs), sum),
        .groups = "drop"
      ) %>%
      tidyr::pivot_longer(
        cols = c(os, pfs),
        names_to = "endpoint",
        values_to = "total_followup"
      ) %>%
      dplyr::rename(
        group = !!rlang::sym(treatment_col)
      )

    # Compute theta
    thetas <- compute_theta(
      event_nums_trt = event_nums_trt,
      event_nums_ctr = event_nums_ctr,
      fl_times       = fl_times,
      ref            = control_level
    )

    # Re-run compute_win_loss with that theta
    df_wl_theta <- compute_win_loss(
      df_inc         = df_increments,
      event_nums_trt = event_nums_trt,
      event_nums_ctr = event_nums_ctr,
      theta          = thetas,
      ref            = control_level
    )
  }

  #--------------------------------------------------------------------------
  # 8) RETURN RESULTS
  #--------------------------------------------------------------------------
  out_list <- list(
    # The cleaned data
    df_clean      = df_clean,

    # Tidy OS/PFS data
    os_long       = os_long,
    pfs_long      = pfs_long,
    df_km         = df_km,
    df_increments = df_increments,

    # Event numbers
    event_nums_trt = event_nums_trt,
    event_nums_ctr = event_nums_ctr,

    # Basic Win/Loss
    df_wl_no_theta = df_wl_no_theta,

    # Theta-based Win/Loss (if compute_theta=TRUE)
    theta       = thetas,       # c(theta_trt, theta_ctr) or NULL
    df_wl_theta = df_wl_theta   # or NULL
  )
  return(out_list)
}


# IPCW to estimate win-loss probabilities
# using individual patient data
ipcw_win_loss <- function(df, ref = NULL, time_grid = NULL) {

  system.time({
  # IPCW to estimate win-loss probabilities ---------------------------------

  df_C <- df |>
    rename(
      trt = rx
    ) |>
    filter(status != 1) |>
    mutate(
      status = ifelse(status == 0, 1, 0)
    )


  # Calculate KM estimates by trt
  objC <- survfit(Surv(time, status) ~ trt, df_C) |>
    summary()

  # Extract the survival probabilities G_C

  G_C <- tibble(
    trt = str_remove(objC$strata, "trt="),
    time = objC$time,
    surv = objC$surv
  )


 if (is.null(time_grid)) {
  # Unique event times
  times <- df |> filter(status != 0) |> pull(time) |> unique() |> sort()
  # times <- seq(0, 6, by = 0.1)
} else {
 times <- time_grid
}

  df_Cw <- df_C |>
    mutate(
      X = time
    ) |>
    group_by(id, trt, X) |>
    reframe(
      time = times,
      iweight = approx(
        G_C$time[G_C$trt == trt], G_C$surv[G_C$trt == trt], pmin(X, time), method = "constant", f = 0, rule = 2
      )$y,
      weight = ((status == 0) | (X >= time)) / iweight
    ) |>
    ungroup()


  # Stablize weight
  # df_Cw <- df_Cw |>
  #   group_by(trt, time) |>
  #   mutate(
  #     # weight = (1 / iweight) / mean(1 / iweight)
  #     weight = weight / mean(weight)
  #   ) |>
  #   ungroup()


  # Compute pairwise win-loss -----------------------------------------------


  df_wide <- df |>
    rename(trt = rx) |>
    pivot_wider(
      names_from = status,
      values_from = time,
      values_fill = Inf
    ) |>
    rename(
      pd = `1`,
      dth = `2`,
      censor = `0`
    ) |>
    select(id, trt, pd, dth, censor)

  # df_wide |>
  #   filter(
  #     dth == pd,
  #     dth < Inf
  #   )


  # Pairwise data

  df_pw <- crossing(
    i = df_wide |> filter(trt != ref),
    j = df_wide |> filter(trt == ref)
  ) |>
    unnest(c(i, j), names_sep = "_") |>
    select(-c(i_trt, j_trt)) |>
    mutate(
      # obj = pmap(list(i_dth, i_pd, i_censor, j_dth, j_pd, j_censor), compare_bivar_time)
      # Follow-up time for dth
      i_x = pmin(i_dth, i_censor),
      j_x = pmin(j_dth, j_censor),
      # Follow-up time for pd
      i_xp = pmin(i_pd, i_x),
      j_xp = pmin(j_pd, j_x),

      # compare on pd
      delta_1 = (j_pd < i_xp) - (i_pd < j_xp),
      t_1 = ifelse(delta_1 == 0, Inf, pmin(i_pd, j_pd)),
      # compare on dth
      delta_2 = (j_dth < i_x) - (i_dth < j_x),
      t_2 = ifelse(delta_2 == 0, Inf, pmin(i_dth, j_dth))
    ) |>
    select(i_id, j_id, delta_1, delta_2,  t_1, t_2)


  # View(df_pw)

  df_pw_wl <- df_pw |>
    pivot_longer(
      cols = c(delta_1, delta_2,   t_1,   t_2),
      names_to = c(".value", "cpn"),
      names_sep = "_"
    ) |>
    group_by(i_id, j_id, t) |>
    slice_max(cpn) |>
    group_by(
      i_id, j_id
    ) |>
    reframe(
      time = times,
      wl = approx(t, delta, time, method = "constant", f = 0, yleft = 0, rule = 2)$y
    ) |>
    # impute na to 0
    mutate(
      wl = replace_na(wl, 0)
    )

  # View(df_pw_wl)

  ## Compute pairwise inverse weight

  df_Cw_trt <- df_Cw |> filter(trt != ref) |> select(id, time, weight)
  df_Cw_ctr <- df_Cw |> filter(trt == ref) |> select(id, time, weight)

  df_pw_wl_wt <- df_pw_wl |>
    left_join(
      df_Cw_trt,
      by = c("i_id" = "id", "time" = "time")
    ) |> left_join(
      df_Cw_ctr,
      by = c("j_id" = "id", "time" = "time")
    ) |>
    mutate(
      weight = weight.x * weight.y
    )

  final_wl <- df_pw_wl_wt |>
    group_by(time) |>
    summarize(
      win = mean((wl == 1) *  weight),
      loss = mean((wl == -1) *  weight)
    ) |>
    mutate(
      win = cummax(win),
      loss = cummax(loss)
    )

})

  final_wl

}



# Model fitting code ------------------------------------------------------

fit_km_ipd <- function(df){

res_theta <- analyze_ipd_data(
  df            = df,
  treatment_col    = "rx",
  control_level    = "0",
  treat_level      = "1",
  compute_theta    = TRUE
)

# res_theta$theta

df_wl <- res_theta$df_wl_theta
df_wl0 <- res_theta$df_wl_no_theta

# IPCW ---------

res_ipcw <- ipcw_win_loss(
  df            = df,
  ref = "0",
  time_grid = seq(0, 3, by = 0.1)
)


final_wl <- res_ipcw


obj <- list(
  df_wl = df_wl,
  df_wl0 = df_wl0,
  final_wl = final_wl
)

obj
}


# Summarize model fit ------------------------------------------------------

summarize_mod_fit <- function(obj, t0){

  # Extract WL estimates
  df_wl <- obj$df_wl # KM w theta
  df_wl0 <- obj$df_wl0 #KM wo theta
  final_wl <- obj$final_wl # IPCW

  # Harmonize with t0

  tibble(
    time = t0,
    # KM w theta
    win_os = approx(x = df_wl$time, y = df_wl$win_os, xout = time, f = 0, rule = 2)$y,
    loss_os = approx(x = df_wl$time, y = df_wl$loss_os, xout = time, f = 0,rule = 2)$y,
    win = approx(x = df_wl$time, y = df_wl$win, xout = time, f = 0,rule = 2)$y,
    loss = approx(x = df_wl$time, y = df_wl$loss, xout = time, f = 0,rule = 2)$y,
    # KM wo theta
    win0 = approx(x = df_wl0$time, y = df_wl0$win, xout = time, f = 0,rule = 2)$y,
    loss0 = approx(x = df_wl0$time, y = df_wl0$loss, xout = time, f = 0,rule = 2)$y,
    # IPCW
    win_ipcw = approx(x = final_wl$time, y = final_wl$win, xout = time, f = 0,rule = 2)$y,
    loss_ipcw = approx(x = final_wl$time, y = final_wl$loss, xout = time, f = 0,rule = 2)$y
  )

}



###############################################################################
# FUNCTION: plot_win_probability
#
# Purpose:
#   - Given a data frame (summary_results) with summary statistics for different
#     scenarios (e.g., "KM Naive," "KM Adjust," "IPCW"), create a faceted ggplot
#     of win probability over time, including ribbons for standard deviations
#     and lines for means.
#
# Arguments:
#   summary_results : A data frame containing columns such as:
#       - time
#       - win_ipcw_mean, win_ipcw_sd
#       - win0_mean,      win0_sd
#       - win_mean,       win_sd
#       - lambdaD,        kappa      (for faceting)
#
# Returns:
#   A ggplot object, which can be printed or further modified.
###############################################################################
plot_win_probability <- function(summary_results) {
  summary_results |>
    # 1) Convert 'lambdaD' and 'kappa' to parseable label expressions,
    #    so facet labels appear as e.g. "lambda[D]==0.2" instead of "0.2".
    dplyr::mutate(
      lambdaD = paste0("lambda[D]==", lambdaD),
      kappa   = paste0("kappa==", kappa) |> as.factor()
    ) |>

    # 2) Initialize ggplot with 'time' on the x-axis
    ggplot(aes(x = time)) +

    # 3) IPCW ribbons & line:
    #    - ribbons use (mean ± sd) for shading
    #    - line uses the mean
    geom_ribbon(aes(
      ymin = win_ipcw_mean - win_ipcw_sd,
      ymax = win_ipcw_mean + win_ipcw_sd,
      fill = "IPCW"
    ), alpha = 0.4) +
    geom_line(aes(y = win_ipcw_mean, color = "IPCW"), linewidth = 1.2) +

    # 4) KM Naive ribbons & line:
    geom_ribbon(aes(
      ymin = win0_mean - win0_sd,
      ymax = win0_mean + win0_sd,
      fill = "KM Naive"
    ), alpha = 0.4) +
    geom_line(aes(y = win0_mean, color = "KM Naive"), linewidth = 1.2) +

    # 5) KM Adjust ribbons & line:
    geom_ribbon(aes(
      ymin = win_mean - win_sd,
      ymax = win_mean + win_sd,
      fill = "KM Adjust"
    ), alpha = 0.4) +
    geom_line(aes(y = win_mean, color = "KM Adjust"), linewidth = 1.2) +

    # 6) Facet by 'kappa' and 'lambdaD'. Use label_parsed so the facet text
    #    (e.g., "kappa==0.1") is parsed as an expression and displayed nicely
    facet_grid(kappa ~ lambdaD, labeller = label_parsed) +

    # 7) Define color and fill scales, linking "KM Naive", "KM Adjust", "IPCW"
    #    to custom colors. The 'limits' ensures the legend ordering.
    scale_color_manual(
      name   = NULL,
      values = c("KM Naive" = "#FF7F00",
                 "KM Adjust"= "#4DAF4A",
                 "IPCW"     = "#00A9FF"),
      limits = c("KM Naive", "KM Adjust", "IPCW")
    ) +
    scale_fill_manual(
      name   = NULL,
      values = c("KM Naive" = "#FF7F00",
                 "KM Adjust"= "#4DAF4A",
                 "IPCW"     = "#00A9FF"),
      limits = c("KM Naive", "KM Adjust", "IPCW")
    ) +

    # 8) Label axes and set minimal theme
    labs(
      x = expression(tau),
      y = expression(w[1*","*0](tau))
    )
    theme_minimal() +

    # 9) Place legend at the top of the plot
    theme(
      legend.position = "top"
    )
}


###############################################################################
# FUNCTION: plot_loss_probability
#
# Purpose:
#   - Given a data frame (summary_results) with summary statistics for different
#     scenarios (e.g., "KM Naive," "KM Adjust," "IPCW"), create a faceted ggplot
#     of loss probability over time, including ribbons for standard deviations
#     and lines for means.
#
# Arguments:
#   summary_results : A data frame containing columns such as:
#       - time
#       - loss_ipcw_mean, loss_ipcw_sd
#       - loss0_mean,      loss0_sd
#       - loss_mean,       loss_sd
#       - lambdaD,        kappa      (for faceting)
#
# Returns:
#   A ggplot object, which can be printed or further modified.
###############################################################################
plot_loss_probability <- function(summary_results) {
  summary_results |>
    # 1) Convert 'lambdaD' and 'kappa' to parseable label expressions,
    #    so facet labels appear as e.g. "lambda[D]==0.2" instead of "0.2".
    dplyr::mutate(
      lambdaD = paste0("lambda[D]==", lambdaD),
      kappa   = paste0("kappa==", kappa) |> as.factor()
    ) |>

    # 2) Initialize ggplot with 'time' on the x-axis
    ggplot(aes(x = time)) +

    # 3) IPCW ribbons & line:
    #    - ribbons use (mean ± sd) for shading
    #    - line uses the mean
    geom_ribbon(aes(
      ymin = loss_ipcw_mean - loss_ipcw_sd,
      ymax = loss_ipcw_mean + loss_ipcw_sd,
      fill = "IPCW"
    ), alpha = 0.4) +
    geom_line(aes(y = loss_ipcw_mean, color = "IPCW"), linewidth = 1.2) +

    # 4) KM Naive ribbons & line:
    geom_ribbon(aes(
      ymin = loss0_mean - loss0_sd,
      ymax = loss0_mean + loss0_sd,
      fill = "KM Naive"
    ), alpha = 0.4) +
    geom_line(aes(y = loss0_mean, color = "KM Naive"), linewidth = 1.2) +

    # 5) KM Adjust ribbons & line:
    geom_ribbon(aes(
      ymin = loss_mean - loss_sd,
      ymax = loss_mean + loss_sd,
      fill = "KM Adjust"
    ), alpha = 0.4) +
    geom_line(aes(y = loss_mean, color = "KM Adjust"), linewidth = 1.2) +

    # 6) Facet by 'kappa' and 'lambdaD'. Use label_parsed so the facet text
    #    (e.g., "kappa==0.1") is parsed as an expression and displayed nicely
    facet_grid(kappa ~ lambdaD, labeller = label_parsed) +

    # 7) Define color and fill scales, linking "KM Naive", "KM Adjust", "IPCW"
    #    to custom colors. The 'limits' ensures the legend ordering.
    scale_color_manual(
      name   = NULL,
      values = c("KM Naive" = "#FF7F00",
                 "KM Adjust"= "#4DAF4A",
                 "IPCW"     = "#00A9FF"),
      limits = c("KM Naive", "KM Adjust", "IPCW")
    ) +
    scale_fill_manual(
      name   = NULL,
      values = c("KM Naive" = "#FF7F00",
                 "KM Adjust"= "#4DAF4A",
                 "IPCW"     = "#00A9FF"),
      limits = c("KM Naive", "KM Adjust", "IPCW")
    ) +

    # 8) Label axes and set minimal theme
    labs(
      x = expression(tau),
      y = expression(w[0*","*1](tau))
    )
  theme_minimal() +

    # 9) Place legend at the top of the plot
    theme(
      legend.position = "top"
    )
}

