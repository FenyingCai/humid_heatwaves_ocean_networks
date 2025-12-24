load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
begin

;;------------------------------------------------ read data ------------------------------------------------;;
f1 = addfile("/public/home/fcai/data0/topography/topo_adapt_era5_air2m.nc", "r")
topo0 = (f1->topo({-90:90},:))





;;------------------------------------------- calculate wbt95 intensity --------------------------------------------;;
wbt = new((/42,92,181,360/), float)
wbt95 = new((/42,92,181,360/), float)

f1 = addfile("/public/home/fcai/data-observation/ERA5_daily/humid_heat/WetBulbTemp_365days_daymax_1979_2023.nc", "r")
wbt(:,:,90:180,:) = (f1->wbt(3:44,151:151+91,90:180,:))   ;; JJA, 0-90N
wbt(:,0:58,0:89,:) = (f1->wbt(3:44,0:58,0:89,:))   ;; DJF, 90S-1S
wbt(:,59:89,0:89,:) = (f1->wbt(3:44,334:364,0:89,:))   ;; DJF, 90S-1S

lat=fspan(-90.0, 90.0, 181)
lat@units="degrees_north"
wbt&latitude = lat


f4 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/1_threshold_1x1/ERA5_wbt95_DJF_1982_2023.nc", "r")
wbt95(:,0:89,0:89,:) = (f4->wbt95(:,:,0:89,:))
f5 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/1_threshold_1x1/ERA5_wbt95_JJA_1982_2023.nc", "r")
wbt95(:,:,90:180,:) = (f5->wbt95)


wbt95_clm = dim_sum_n_Wrap(dim_sum_n_Wrap(wbt*wbt95, 1), 0) / dim_sum_n_Wrap(dim_sum_n_Wrap(wbt95, 1), 0)
wbt95_clm = where(topo0.ge.0.0, wbt95_clm, 1.0e20)
wbt95_clm@_FillValue = 1.0e20
copy_VarCoords(topo0, wbt95_clm)

printVarSummary(wbt95_clm)
printMinMax(wbt95_clm, 1)





;;--------------------------- to nc file -----------------------------;;
output="/public/home/fcai/extreme5_degree/25_humid_heat/3_analysis1982_2023/fig1_similarity/"
system("/bin/rm -f "+output+"ERA5_wbt95_clm_1982_2023.nc")
ncdf=addfile(output+"ERA5_wbt95_clm_1982_2023.nc","c")
ncdf->wbt95_clm=wbt95_clm







end







