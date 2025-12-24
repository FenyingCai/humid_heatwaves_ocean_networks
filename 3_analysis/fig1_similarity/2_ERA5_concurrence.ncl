load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
begin

;;------------------------------------------------ read data ------------------------------------------------;;
f1 = addfile("/public/home/fcai/data0/topography/topo_adapt_era5_air2m.nc", "r")
topo0 = (f1->topo({-90:90},:))





;;------------------------------------------------ read data ------------------------------------------------;;
tmax95 = new((/42,92,91+90,360/), float)
f1 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/1_threshold_1x1/ERA5_tmax95_DJF_1982_2023.nc", "r")
tmax95(:,0:89,0:89,:) = (f1->tmax95(:,:,0:89,:))
f2 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/1_threshold_1x1/ERA5_tmax95_JJA_1982_2023.nc", "r")
tmax95(:,:,90:180,:) = (f2->tmax95)

wbt95 = new((/42,92,91+90,360/), float)
f3 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/1_threshold_1x1/ERA5_wbt95_DJF_1982_2023.nc", "r")
wbt95(:,0:89,0:89,:) = (f3->wbt95(:,:,0:89,:))
f4 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/1_threshold_1x1/ERA5_wbt95_JJA_1982_2023.nc", "r")
wbt95(:,:,90:180,:) = (f4->wbt95)

q95 = new((/42,92,91+90,360/), float)
f5 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/1_threshold_1x1/ERA5_qsurface95_DJF_1982_2023.nc", "r")
q95(:,0:89,0:89,:) = (f5->q95(:,:,0:89,:))
f6 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/1_threshold_1x1/ERA5_qsurface95_JJA_1982_2023.nc", "r")
q95(:,:,90:180,:) = (f6->q95)

printMinMax(tmax95, 1)
printMinMax(wbt95, 1)
printMinMax(q95, 1)



per_tmax = dim_sum_n_Wrap(dim_sum_n_Wrap(tmax95*wbt95, 1), 0) / dim_sum_n_Wrap(dim_sum_n_Wrap(wbt95, 1), 0) 
per_q = dim_sum_n_Wrap(dim_sum_n_Wrap(q95*wbt95, 1), 0) / dim_sum_n_Wrap(dim_sum_n_Wrap(wbt95, 1), 0) 
per_tmax = where(topo0.gt.0.0, per_tmax, 1.0e20)
per_tmax@_FillValue = 1.0e20
per_q = where(topo0.gt.0.0, per_q, 1.0e20)
per_q@_FillValue = 1.0e20
copy_VarCoords(topo0, per_tmax)
copy_VarCoords(topo0, per_q)
printVarSummary(per_tmax)

per_diff = per_tmax - per_q
copy_VarCoords(topo0, per_diff)

print("  ")
printMinMax(per_tmax, 1)
printMinMax(per_q, 1)
printMinMax(per_diff, 1)





;;--------------------------- to nc file -----------------------------;;
output="/public/home/fcai/extreme5_degree/25_humid_heat/3_analysis1982_2023/fig1_similarity/"
system("/bin/rm -f "+output+"ERA5_top5_concurrence_1982_2023.nc")
ncdf=addfile(output+"ERA5_top5_concurrence_1982_2023.nc","c")
ncdf->per_tmax=per_tmax
ncdf->per_q=per_q



; f3 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/3_analysis/fig1_data/ERA5_top5_concurrence.nc", "r")
; per_tmax = (f3->per_tmax)
; per_q = (f3->per_q)

; per_diff = per_tmax - per_q
; copy_VarCoords(topo0, per_diff)

; cor_diff = smth9_Wrap(cor_diff, 0.5, 0.25, True)
; ; cor_diff = smth9_Wrap(cor_diff, 0.5, 0.25, True)






end







