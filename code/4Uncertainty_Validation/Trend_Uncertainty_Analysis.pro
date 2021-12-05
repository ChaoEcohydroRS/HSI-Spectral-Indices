work_folder = 'G:\Default\Uncertainty_analysis'
date_file = work_folder + '\Caspian.txt'

; Convert 8-day series to monthly series, and then to monthly anomaly series
module_conversion_from_8day_to_monthly_anomalies, date_file, $
  lake_volume_anomaly_monthly = lake_volume_anomaly_monthly, lake_volume_e_monthly = lake_volume_e_monthly, $
  valid_array_series_key = valid_array, valid_volume_key = valid_volume, valid_volume_e_key = valid_volume_e, $
  lake_volume_monthly_original = lake_volume_monthly_original  ;valid_x: null-value free

; Linear fitting
g_LSR = LINFIT(valid_array, valid_volume)  
; Calculate fitting residual of monthly anomalies
residual_std = stddev( (valid_array*g_LSR[1] + g_LSR[0])-valid_volume, /NaN)  

; Combine volume estimation errors and fitting residuals to propagate total error series (updated valid_volume_e)
valid_volume_e =  sqrt(valid_volume_e^2 + residual_std^2)   

; Monte Carlo uncertainty simulation: please read the module .pro file for input and output variables
simulation_number = 10000 ;recommended to be greater than 1000
module_monte_carlo_uncertainty_simulation, valid_array, valid_volume, valid_volume_e, simulation_number, $
  slope_e = slope_e, intercept_e = intercept_e,  confidence_band = confidence_band

; Plot result
window,4
plot, lake_volume_anomaly_monthly, title = 'Fitting with uncertainty interval'
; Add error bars
ERRPLOT, lake_volume_anomaly_monthly-lake_volume_e_monthly, lake_volume_anomaly_monthly+lake_volume_e_monthly
oplot, valid_array, (valid_array*g_LSR[1] + g_LSR[0])
; Add confidence bands (interval)
oplot, valid_array, (valid_array*g_LSR[1] + g_LSR[0])+confidence_band ;upper band
oplot, valid_array, (valid_array*g_LSR[1] + g_LSR[0])-confidence_band ;lower band

; Print resul
print, 'Trend: per month'
print, string(g_LSR[1]), ' +/- ', string(slope_e)

end