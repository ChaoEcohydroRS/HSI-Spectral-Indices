pro module_conversion_from_8day_to_monthly_anomalies, date_file, $
  lake_volume_anomaly_monthly = lake_volume_anomaly_monthly, lake_volume_e_monthly = lake_volume_e_monthly, $
  valid_array_series_key = valid_array, valid_volume_key = valid_volume, valid_volume_e_key = valid_volume_e, $
  lake_volume_monthly_original =  lake_volume_monthly
  ; Coded by Dr. J. Wang, 5/2016
  ; jidawang@ksu.edu

  
  OPENR,1,date_file
  line_value = strarr(1, 1) ;this is required
  lake_years = []
  lake_months = []
  lake_days = []
  lake_volume = []
  lake_volume_e = []
  while ~ EOF(1) do begin
    Readf, 1, line_value
    ;split the line_value by space
    line_value_array = STRSPLIT(line_value,  /extract) ;default by tab
    this_date = long(line_value_array[0])
    this_year = Fix(StrMid(StrTrim(this_date,2), 0, 4))
    dayofyear = Fix(StrMid(StrTrim(this_date,2), 4, 3))
    CALDAT, JULDAY(1, dayofyear, this_year), this_month, this_day
    ;Print, this_month, this_day
    lake_years = [lake_years, this_year]
    lake_months = [lake_months, this_month]
    lake_days = [lake_days, this_day]

    lake_volume = [lake_volume, float(line_value_array[7])]
    lake_volume_e = [lake_volume_e, float(line_value_array[8])]
  endwhile
  close,1
  window,0
  plot, lake_volume, title = '8-day series';, colo = 255;, yrange = [-5000, 5000]
  ERRPLOT, lake_volume-lake_volume_e, lake_volume+lake_volume_e

  years_GLDAS = []
  months_GLDAS = []
  for year_i = 2002, 2016 do begin ;although GLDAS original files start from 2000/01, the processed .sav starts from 2002/01
    for month_i = 1, 12 do begin
      years_GLDAS = [years_GLDAS, year_i]
      months_GLDAS = [months_GLDAS, month_i]
    endfor
  endfor
  start_ind_GLDAS = where(years_GLDAS eq 2002 and months_GLDAS eq 4)
  end_ind_GLDAS = where(years_GLDAS eq 2016 and months_GLDAS eq 3)
  years_GLDAS = years_GLDAS[start_ind_GLDAS:end_ind_GLDAS]
  months_GLDAS = months_GLDAS[start_ind_GLDAS:end_ind_GLDAS]
  n = n_elements(years_GLDAS)
  lake_volume_monthly = []
  lake_volume_e_monthly = []
  for i = 0, n-1 do begin  ;now trim down to 2002/4 to 2016/3
    this_month = months_GLDAS[i]
    this_year = years_GLDAS[i]
    this_indices = where(lake_years eq this_year and lake_months eq this_month, count)
    if count GT 0 then begin
      lake_volume_monthly = [lake_volume_monthly, mean(lake_volume[this_indices], /NaN)]
      lake_volume_e_monthly = [lake_volume_e_monthly, sqrt(mean(lake_volume_e[this_indices]^2, /NaN))]
    endif else begin
      lake_volume_monthly = [lake_volume_monthly, !Valus.F_NaN]
      lake_volume_e_monthly = [lake_volume_e_monthly, !Values.F_NaN]
    endelse
  endfor
  window,1
  plot, lake_volume_monthly, title = 'monthly series';, colo = 255;, yrange = [-5000, 5000]
  ERRPLOT, lake_volume_monthly-lake_volume_e_monthly, lake_volume_monthly+lake_volume_e_monthly
  ;now convert to monthly anomalies
  lake_volume_anomaly_monthly = []
  for i = 0, n-1 do begin
    this_month = months_GLDAS[i]
    ind = where(months_GLDAS eq this_month)
    lake_volume_anomaly_monthly = [lake_volume_anomaly_monthly, lake_volume_monthly[i] - mean(lake_volume_monthly[ind],/NaN)]
  endfor
  window,2
  plot, lake_volume_anomaly_monthly, title = 'monthly anomaly series';, colo = 255;, yrange = [-5000, 5000]
  ERRPLOT, lake_volume_anomaly_monthly-lake_volume_e_monthly, lake_volume_anomaly_monthly+lake_volume_e_monthly


  valid_array = []
  valid_volume = []
  valid_volume_e = []
  for i = 0, n-1 do begin
    if finite(lake_volume_anomaly_monthly[i]) EQ 1 then begin
      valid_array = [valid_array, i]
      valid_volume = [valid_volume, lake_volume_anomaly_monthly[i]]
      valid_volume_e = [valid_volume_e, lake_volume_e_monthly[i]]
    endif
  endfor
  return
end
