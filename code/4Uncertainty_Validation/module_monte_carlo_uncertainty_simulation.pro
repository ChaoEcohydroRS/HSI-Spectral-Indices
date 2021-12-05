pro module_monte_carlo_uncertainty_simulation, x, y, y_e, monte_carlo_size, slope_e = slope_e, intercept_e = intercept_e, confidence_band = confidence_band
  ; Module coded by Dr. J. Wang, 5/2016
  ; jidawang@ksu.edu

  ; Inputs (please note: NULL values are not permitted in any input variable. Please remove null values in the simulation)
  ; x: X series 
  ; y: Y series
  ; y_e: Y uncertainty series
  ; monte_carlo_size: the number/interation of simulation (recommended to be bigger than 1000)
 
  ; Outputs 
  ; slope_e: estimated uncertain for fitting slope
  ; intercept_e: estimated uncertain for intercept
  ; confidence_band: estimated confidence intervals (half interval size)

  Trends_monte_carlo = []
  nn =n_elements(x)
  e_simulated = fltarr(nn, monte_carlo_size)
  for i = 0, nn-1 do begin
    e_simulated[i, *] = transpose(RANDOMN(seed, monte_carlo_size) * y_e[i])
  endfor

  for i = 0, monte_carlo_size-1 do begin
    Y_i = y + e_simulated(*,i)
    coefficient_values =LINFIT(x, Y_i)
    Trends_monte_carlo = [[Trends_monte_carlo], [coefficient_values]]
  endfor
  all_intercepts = Trends_monte_carlo[0,*]
  all_slopes = Trends_monte_carlo[1,*]
  slope_e = stddev(all_slopes,/NaN);*1.96
  intercept_e = stddev(all_intercepts, /NaN);*1.96 ;all 95% interval

  ; return confidence bands
  confidence_band = fltarr(nn,1)
  for i = 0, nn-1 do begin
    confidence_band[i] = stddev(x[i]*all_slopes + all_intercepts);*1.96
  endfor

  return
end