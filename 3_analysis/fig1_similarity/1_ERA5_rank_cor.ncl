load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
begin

;;------------------------------------------------ read data ------------------------------------------------;;
f1 = addfile("/public/home/fcai/data0/topography/topo_adapt_era5_air2m.nc", "r")
topo0 = (f1->topo({-90:90},:))





;;-------------------------------------------------- read data --------------------------------------------------;;
tmax_stz = new((/42,92,91+90,360/), float)
f1 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/1_threshold_1x1/ERA5_tmax_stz_DJF_1982_2023.nc", "r")
tmax_stz(:,0:89,0:89,:) = (f1->tmax_stz(:,:,0:89,:))
f2 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/1_threshold_1x1/ERA5_tmax_stz_JJA_1982_2023.nc", "r")
tmax_stz(:,:,90:180,:) = (f2->tmax_stz)

wbt_stz = new((/42,92,91+90,360/), float)
f3 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/1_threshold_1x1/ERA5_wbt_stz_DJF_1982_2023.nc", "r")
wbt_stz(:,0:89,0:89,:) = (f3->wbt_stz(:,:,0:89,:))
f4 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/1_threshold_1x1/ERA5_wbt_stz_JJA_1982_2023.nc", "r")
wbt_stz(:,:,90:180,:) = (f4->wbt_stz)

q_stz = new((/42,92,91+90,360/), float)
f5 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/1_threshold_1x1/ERA5_qsurface_stz_DJF_1982_2023.nc", "r")
q_stz(:,0:89,0:89,:) = (f5->q_stz(:,:,0:89,:))
f6 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/1_threshold_1x1/ERA5_qsurface_stz_JJA_1982_2023.nc", "r")
q_stz(:,:,90:180,:) = (f6->q_stz)



tmax_stz3D = new((/42*92,181,360/), float)
wbt_stz3D = new((/42*92,181,360/), float)
q_stz3D = new((/42*92,181,360/), float)
do iyear = 0,41
print("  iyear  =  "+iyear)
tmax_stz3D(iyear*92:iyear*92+91,:,:) = tmax_stz(iyear,:,:,:)
wbt_stz3D(iyear*92:iyear*92+91,:,:) = wbt_stz(iyear,:,:,:)
q_stz3D(iyear*92:iyear*92+91,:,:) = q_stz(iyear,:,:,:)
end do

cor_wbt_tmax = spcorr_n(wbt_stz3D, tmax_stz3D, 0)
cor_wbt_tmax = where(topo0.gt.0.0, cor_wbt_tmax, 1.0e20)
cor_wbt_tmax@_FillValue = 1.0e20
copy_VarCoords(topo0, cor_wbt_tmax)
printVarSummary(cor_wbt_tmax)

cor_wbt_q = spcorr_n(wbt_stz3D, q_stz3D, 0)
cor_wbt_q = where(topo0.gt.0.0, cor_wbt_q, 1.0e20)
cor_wbt_q@_FillValue = 1.0e20
copy_VarCoords(topo0, cor_wbt_q)


delete([/tmax_stz, wbt_stz, q_stz/])
delete([/tmax_stz3D, wbt_stz3D, q_stz3D/])



;;--------------------------- to nc file -----------------------------;;
output="/public/home/fcai/extreme5_degree/25_humid_heat/3_analysis1982_2023/fig1_similarity/"
system("/bin/rm -f "+output+"ERA5_rank_cor_1982_2023.nc")
ncdf=addfile(output+"ERA5_rank_cor_1982_2023.nc","c")
ncdf->cor_wbt_tmax=cor_wbt_tmax
ncdf->cor_wbt_q=cor_wbt_q


; f3 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/3_analysis/fig1_data/ERA5_rank_cor.nc", "r")
; cor_wbt_tmax = (f3->cor_wbt_tmax)
; cor_wbt_q = (f3->cor_wbt_q)

; cor_diff = cor_wbt_tmax - cor_wbt_q
; copy_VarCoords(topo0, cor_diff)
; cor_diff = smth9_Wrap(cor_diff, 0.5, 0.25, True)
; cor_diff = smth9_Wrap(cor_diff, 0.5, 0.25, True)






end







