load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
begin


;;------------------------------------------------------------------------------------------------------;;
;;------------------------------------------------  NH  ------------------------------------------------;;
f1 = addfile("/public/data2/model/fcai/data0/topography/topo_2x2.nc", "r")
topo0 := (f1->topo({0:90},:))
dim_0 = dimsizes(topo0)

q = new((/42,92,45,180/), float)
u = new((/42,92,45,180/), float)
v = new((/42,92,45,180/), float)
do i = 0,41
print(" year = " + (i+1982))
f1 = addfile("/public/data2/model/fcai/data0/ERA5_q2m_daily2x2/q2m."+flt2string(i+1982)+".nc", "r")
q(i,:,:,:) = (f1->q2m(151:151+91,{0:90},:))   ;; JJA
f2 = addfile("/public/data2/model/fcai/data0/ERA5_u10m_daily2x2/u10m."+flt2string(i+1982)+".nc", "r")
u(i,:,:,:) = (f2->u10m(151:151+91,{0:90},:))
f3 = addfile("/public/data2/model/fcai/data0/ERA5_v10m_daily2x2/v10m."+flt2string(i+1982)+".nc", "r")
v(i,:,:,:) = (f3->v10m(151:151+91,{0:90},:))
end do
qus = q*u
qvs = q*v
copy_VarCoords(q, qus)
copy_VarCoords(q, qvs)
delete([/q,u,v/])

qu_yearly_NH = dim_avg_n_Wrap(qus, 1)
qv_yearly_NH = dim_avg_n_Wrap(qvs, 1)





;;------------------------------------------------------------------------------------------------------;;
;;------------------------------------------------  SH  ------------------------------------------------;;
f1 = addfile("/public/data2/model/fcai/data0/topography/topo_2x2.nc", "r")
topo0 := (f1->topo({-90:0},:))
dim_0 = dimsizes(topo0)

q = new((/42,90,45,180/), float)
u = new((/42,90,45,180/), float)
v = new((/42,90,45,180/), float)
do i = 0,41
print(" year = " + (i+1982))
f1 = addfile("/public/data2/model/fcai/data0/ERA5_q2m_daily2x2/q2m."+flt2string(i+1982)+".nc", "r")
q(i,0:58,:,:) = (f1->q2m(0:58,{-90:0},:))   ;; DJF
f1 = addfile("/public/data2/model/fcai/data0/ERA5_q2m_daily2x2/q2m."+flt2string(i+1982)+".nc", "r")
q(i,59:89,:,:) = (f1->q2m(334:364,{-90:0},:)) 
f2 = addfile("/public/data2/model/fcai/data0/ERA5_u10m_daily2x2/u10m."+flt2string(i+1982)+".nc", "r")
u(i,0:58,:,:) = (f2->u10m(0:58,{-90:0},:))
f2 = addfile("/public/data2/model/fcai/data0/ERA5_u10m_daily2x2/u10m."+flt2string(i+1982)+".nc", "r")
u(i,59:89,:,:) = (f2->u10m(334:364,{-90:0},:))
f3 = addfile("/public/data2/model/fcai/data0/ERA5_v10m_daily2x2/v10m."+flt2string(i+1982)+".nc", "r")
v(i,0:58,:,:) = (f3->v10m(0:58,{-90:0},:))
f3 = addfile("/public/data2/model/fcai/data0/ERA5_v10m_daily2x2/v10m."+flt2string(i+1982)+".nc", "r")
v(i,59:89,:,:) = (f3->v10m(334:364,{-90:0},:))
end do
qus := q*u
qvs := q*v
copy_VarCoords(q, qus)
copy_VarCoords(q, qvs)
delete([/q,u,v/])

qu_yearly_SH = dim_avg_n_Wrap(qus, 1)
qv_yearly_SH = dim_avg_n_Wrap(qvs, 1)



;;------------------------------------------------------------------------------------------------------;;
;;---------------------------------------------  NH + SH  ----------------------------------------------;;
f1 = addfile("/public/data2/model/fcai/data0/topography/topo_2x2.nc", "r")
topo0 := (f1->topo({-90:90},:))
dim_0 := dimsizes(topo0)

qu_yearly = new((/42,90,180/), float)
qv_yearly = new((/42,90,180/), float)
qu_yearly(:,0:44,:) = qu_yearly_SH
qu_yearly(:,45:89,:) = qu_yearly_NH
qv_yearly(:,0:44,:) = qv_yearly_SH
qv_yearly(:,45:89,:) = qv_yearly_NH

qu_trend = regCoef_n(fspan(1,42,42), qu_yearly, 0, 0) * 10.0 * 1000.0
qv_trend = regCoef_n(fspan(1,42,42), qv_yearly, 0, 0) * 10.0 * 1000.0
qu_trend_cor = escorc_n(fspan(1,42,42), qu_yearly, 0, 0)
qv_trend_cor = escorc_n(fspan(1,42,42), qv_yearly, 0, 0)
qu_trend = where(abs(qu_trend_cor).ge.0.25.or.abs(qv_trend_cor).ge.0.25, qu_trend, 1.0e20)
qv_trend = where(abs(qu_trend_cor).ge.0.25.or.abs(qv_trend_cor).ge.0.25, qv_trend, 1.0e20)
qu_trend@_FillValue = 1.0e20
qv_trend@_FillValue = 1.0e20
copy_VarCoords(topo0, qu_trend)
copy_VarCoords(topo0, qv_trend)

printMinMax(qv_trend, 1)




distance_to_land = new((/dim_0(0),dim_0(1)/), float)
distance_to_ocean = new((/dim_0(0),dim_0(1)/), float)
f3 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/1_threshold_2x2/inshore_uv_SH2x2.nc", "r")
distance_to_land(0:44,:) = f3->distance_to_land
distance_to_ocean(0:44,:) = f3->distance_to_ocean
f4 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/1_threshold_2x2/inshore_uv_NH2x2.nc", "r")
distance_to_land(45:89,:) = f4->distance_to_land
distance_to_ocean(45:89,:) = f4->distance_to_ocean

qu_trend = where(distance_to_land.le.1800.0, qu_trend, 1.0e20)
qv_trend = where(distance_to_land.le.1800.0, qv_trend, 1.0e20)
qu_trend = where(distance_to_ocean.le.1800.0, qu_trend, 1.0e20)
qv_trend = where(distance_to_ocean.le.1800.0, qv_trend, 1.0e20)
quv_trend = sqrt(qu_trend*qu_trend  + qv_trend*qv_trend)
qu_trend = where(quv_trend.ge.1.5, qu_trend, 1.0e20)
qv_trend = where(quv_trend.ge.1.5, qv_trend, 1.0e20)
copy_VarCoords(topo0, qu_trend)
copy_VarCoords(topo0, qv_trend)










;;------------------------------------------------ long-term trends in wbt ------------------------------------------------;;
f1 = addfile("/public/data2/model/fcai/data0/topography/topo_adapt_era5_air2m.nc", "r")
topo1x1 = (f1->topo({-90:90},:))
dim_1x1 = dimsizes(topo1x1)

f1 = addfile("/public/data2/model/fcai/data0/topography/topo_2x2.nc", "r")
topo0 := (f1->topo({-90:90},:))
dim_0 := dimsizes(topo0)

wbt95 = new((/42,92,dim_1x1(0),dim_1x1(1)/), float)
f3 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/1_threshold_1x1/ERA5_wbt95_DJF_1982_2023.nc", "r")
wbt95(:,0:89,0:89,:) = (f3->wbt95(:,:,0:89,:))
f4 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/1_threshold_1x1/ERA5_wbt95_JJA_1982_2023.nc", "r")
wbt95(:,:,90:180,:) = (f4->wbt95(:,:,0:90,:))

wbt95_annual = dim_sum_n_Wrap(wbt95, 1)
wbt95_trend = regCoef_n(fspan(1,42,42), wbt95_annual, 0, 0) * 10.0
wbt95_trend_cor = escorc_n(fspan(1,42,42), wbt95_annual, 0, 0)
wbt95_trend = where(topo1x1.ge.0.0, wbt95_trend, 1.0e20)
wbt95_trend_cor = where(topo1x1.ge.0.0, wbt95_trend_cor, 1.0e20)
wbt95_trend@_FillValue = 1.0e20
wbt95_trend_cor@_FillValue = 1.0e20
wbt95_trend = where(wbt95_trend_cor.ge.0.25.or.wbt95_trend_cor.le.-0.25, wbt95_trend, 1.0e20)
copy_VarCoords(topo1x1, wbt95_trend)
copy_VarCoords(topo1x1, wbt95_trend_cor)
print("  wbt95 trend  =  ")
printMinMax(wbt95_trend, 1)


























;;------------------------------------------------------------------------------------------------------------------------------------;;
;;----------------------------------------------------------------- fig (b) ----------------------------------------------------------;;
;;------------------------ read lat, lon -------------------------;;
f00 = addfile("/public/data2/model/fcai/data0/topography/topo_2x2.nc","r")
topo0 := f00->topo({-90:90},:)
lat := topo0&lat
lon := topo0&lon
dim_0 := dimsizes(topo0)
topo_01 = where(topo0.ge.0.0, 1.0, 0.0)
copy_VarCoords(topo0, topo_01)



;;-------------------------------------------------------------------------------------------------------------------------------------;;
wbt_yearly = new((/42,dim_0(0),dim_0(1)/), float)
f1 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/2_network/part1_wbt95_sst95_NH/threshold/ERA5_wbt95_JJA_1982_2023.nc","r")
wbt_yearly(:,45:89,:) = int2flt(dim_sum_n_Wrap(f1->wbt95(:,:,:,:), 1))
f2 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/2_network/part1_wbt95_sst95_NH/threshold/ERA5_wbt95_DJF_1982_2023.nc","r")
wbt_yearly(:,0:44,:) = int2flt(dim_sum_n_Wrap(f2->wbt95(:,:,:,:), 1))
wbt_yearly@_FillValue = 1.0e20
value = 1.0e20
if (any(isnan_ieee(wbt_yearly))) then
      replace_ieeenan (wbt_yearly, value, 0)
      wbt_yearly@_FillValue = value
end if
wbt_yearly@_FillValue = 1.0e20


wbt_yearly!1 = "lat"
wbt_yearly!2 = "lon"
wbt_yearly&lat = lat
wbt_yearly&lon = lon
wbt_yearly_detrend = dtrend_msg_n(fspan(1, 42, 42), wbt_yearly, False, False, 0)
wbt_yearly_detrend@_FillValue = 1.0e20
copy_VarCoords(wbt_yearly, wbt_yearly_detrend)

wbt_trend = regCoef_n(fspan(1,42,42), wbt_yearly, 0, 0) * 10.0
wbt_trend_cor = escorc_n(fspan(1,42,42), wbt_yearly, 0, 0)
wbt_trend = where(topo0.ge.0.0, wbt_trend, 1.0e20)
wbt_trend_cor = where(topo0.ge.0.0, wbt_trend_cor, 1.0e20)
wbt_trend@_FillValue = 1.0e20
wbt_trend_cor@_FillValue = 1.0e20
wbt_trend_original = wbt_trend
wbt_trend = where(wbt_trend_cor.ge.0.25.or.wbt_trend_cor.le.-0.25, wbt_trend, 1.0e20)
copy_VarCoords(topo0, wbt_trend)
copy_VarCoords(topo0, wbt_trend_original)
copy_VarCoords(topo0, wbt_trend_cor)






sst_yearly = new((/42,dim_0(0),dim_0(1)/), float)
f1 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/2_network/part1_wbt95_sst95_NH/threshold/OISST_sst95_JJA_1982_2023.nc","r")
sst_yearly(:,45:89,:) = int2flt(dim_sum_n_Wrap(f1->sst95(:,:,:,:), 1))
f2 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/2_network/part1_wbt95_sst95_NH/threshold/OISST_sst95_DJF_1982_2023.nc","r")
sst_yearly(:,0:44,:) = int2flt(dim_sum_n_Wrap(f2->sst95(:,:,:,:), 1))
sst_yearly@_FillValue = 1.0e20
value = 1.0e20
if (any(isnan_ieee(sst_yearly))) then
      replace_ieeenan (sst_yearly, value, 0)
      sst_yearly@_FillValue = value
end if
sst_yearly@_FillValue = 1.0e20

sst_yearly!1 = "lat"
sst_yearly!2 = "lon"
sst_yearly&lat = lat
sst_yearly&lon = lon
sst_yearly_detrend = dtrend_msg_n(fspan(1, 42, 42), sst_yearly, False, False, 0)
sst_yearly_detrend@_FillValue = 1.0e20
copy_VarCoords(sst_yearly, sst_yearly_detrend)

sst_trend = regCoef_n(fspan(1,42,42), sst_yearly, 0, 0) * 10.0
sst_trend_cor = escorc_n(fspan(1,42,42), sst_yearly, 0, 0)
sst_trend = where(topo0.le.0.0, sst_trend, 1.0e20)
sst_trend_cor = where(topo0.le.0.0, sst_trend_cor, 1.0e20)
sst_trend@_FillValue = 1.0e20
sst_trend_cor@_FillValue = 1.0e20
copy_VarCoords(topo0, sst_trend)
copy_VarCoords(topo0, sst_trend_cor)







Teleconnection_SST = new((/42,dim_0(0),dim_0(1)/), float)
f1 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/3_analysis1982_2023/fig4_grid_based_SST/Teleconnection_SST_extremes_mostconnected_NH.nc","r")
Teleconnection_SST(:,45:89,:) = f1->Teleconnection_SST
f2 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/3_analysis1982_2023/fig4_grid_based_SST/Teleconnection_SST_extremes_mostconnected_SH.nc","r")
Teleconnection_SST(:,0:44,:) = f2->Teleconnection_SST
Teleconnection_SST@_FillValue = 1.0e20
value = 1.0e20
if (any(isnan_ieee(Teleconnection_SST))) then
      replace_ieeenan (Teleconnection_SST, value, 0)
      Teleconnection_SST@_FillValue = value
end if
Teleconnection_SST_detrend = dtrend_msg_n(fspan(1, 42, 42), Teleconnection_SST, False, False, 0)
Teleconnection_SST_detrend@_FillValue = 1.0e20
copy_VarCoords(wbt_yearly, Teleconnection_SST_detrend)





;;------------------------------------------------------------------------------------------;;
;;----------------------- explained variances R2 and long-term trend  ----------------------;;
cor_detrend = escorc_n(wbt_yearly_detrend, Teleconnection_SST_detrend, 0, 0)
; reg_detrend = regCoef_n(wbt_yearly_detrend, Teleconnection_SST_detrend, 0, 0)
reg_detrend = cor_detrend * (dim_stddev_n_Wrap(wbt_yearly_detrend,0)/dim_stddev_n_Wrap(Teleconnection_SST_detrend,0))
R2_detrend = cor_detrend * cor_detrend * 100.0
cor_detrend = where(topo0.ge.0.0, cor_detrend, 1.0e20)
reg_detrend = where(topo0.ge.0.0, reg_detrend, 1.0e20)
R2_detrend = where(topo0.ge.0.0, R2_detrend, 1.0e20)
cor_detrend@_FillValue = 1.0e20
reg_detrend@_FillValue = 1.0e20
R2_detrend@_FillValue = 1.0e20
copy_VarCoords(wbt_yearly(0,:,:), cor_detrend)
copy_VarCoords(wbt_yearly(0,:,:), reg_detrend)
copy_VarCoords(wbt_yearly(0,:,:), R2_detrend)
R2_detrend_lat = dim_median_n(R2_detrend, 1)





trend_Teleconnection_SST = regCoef_n(fspan(1,42,42), Teleconnection_SST, 0, 0) * 10.0
copy_VarCoords(wbt_yearly(0,:,:), trend_Teleconnection_SST)

estimated_trend = reg_detrend*trend_Teleconnection_SST
estimated_trend = where(topo0.ge.0.0, estimated_trend, 1.0e20)
estimated_trend@_FillValue = 1.0e20
copy_VarCoords(wbt_yearly(0,:,:), estimated_trend)


estimated_percentage = estimated_trend / wbt_trend
copy_VarCoords(wbt_yearly(0,:,:), estimated_percentage)












;;--------------------------------------------------------------------------------------------;;
;;----------------------------------  Local  Concurrences  -----------------------------------;;
f3 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/3_analysis1982_2023/fig1_local_concurrence/wbt95_wbt95_LCle2000_42yrs.nc","r")
LC_wbt95_wbt95 := f3->degree

LC_trend = regCoef_n(fspan(1,42,42), LC_wbt95_wbt95, 0, 0) * 10.0 * 10.0
LC_trend_cor = escorc_n(fspan(1,42,42), LC_wbt95_wbt95, 0, 0)
LC_trend = where(topo0.ge.0.0, LC_trend, 1.0e20)
LC_trend_cor = where(topo0.ge.0.0, LC_trend_cor, 1.0e20)
LC_trend@_FillValue = 1.0e20
LC_trend_cor@_FillValue = 1.0e20
LC_trend_original = LC_trend
copy_VarCoords(topo0, LC_trend_original)
LC_trend = where(LC_trend_cor.ge.0.25.or.LC_trend_cor.le.-0.25, LC_trend, 1.0e20)
copy_VarCoords(topo0, LC_trend)

LC_wbt95_wbt95_detrend = dtrend_msg_n(fspan(1, 42, 42), LC_wbt95_wbt95, False, False, 0)
cor_LC_detrended = escorc_n(LC_wbt95_wbt95_detrend, Teleconnection_SST_detrend, 0, 0)
copy_VarCoords(wbt_yearly(0,:,:), cor_LC_detrended)
R2_LC_detrend = cor_LC_detrended * cor_LC_detrended * 100.0
copy_VarCoords(wbt_yearly(0,:,:), R2_LC_detrend)

reg_LC_detrend = cor_LC_detrended * (dim_stddev_n_Wrap(LC_wbt95_wbt95_detrend,0)/dim_stddev_n_Wrap(Teleconnection_SST_detrend,0))

estimated_LC_trend = reg_LC_detrend*trend_Teleconnection_SST * 10.0
estimated_LC_trend = where(topo0.ge.0.0, estimated_LC_trend, 1.0e20)
estimated_LC_trend@_FillValue = 1.0e20
copy_VarCoords(wbt_yearly(0,:,:), estimated_LC_trend)

estimated_LC_percentage = estimated_LC_trend / LC_trend
copy_VarCoords(wbt_yearly(0,:,:), estimated_LC_percentage)








;;************************************************** plots **********************************************;;
res = True
res@tmBorderThicknessF = 1.2
res@mpPerimOn = False
res@tmXTBorderOn = False
res@tmXBBorderOn = False
res@tmYLBorderOn = False
res@tmYRBorderOn = False
res@gsnDraw      =  False
res@gsnFrame     =  False
res@gsnAddCyclic =  True
res@gsnRightString       = ""
res@gsnLeftString        = ""
res@gsnLeftStringFontHeightF   = 0.04
res@gsnRightStringFontHeightF  = 0.04
res@gsnCenterStringFontHeightF  = 0.018
res@tmXTLabelsOn  = False
res@tmYRLabelsOn  = False
vcres = res
res@tmXTOn        = False
res@tmYROn        = False
res@tmXBOn        = False
res@tmYLOn        = False
res@tmXBLabelFontHeightF = 0.013
res@tmYLLabelFontHeightF = 0.013
res@tmXBTickSpacingF = 45.0
res@tmYLTickSpacingF = 20.0




res@tmYLMajorThicknessF = 1.5
res@tmXBMajorThicknessF = 1.5
res@tmYLMinorThicknessF = 1.0
res@tmXBMinorThicknessF = 1.0
res@tmXBMajorLengthF = 0.011
res@tmYLMinorLengthF = 0.011
res@tmXBMinorLengthF = 0.008
res@tmYLMinorLengthF = 0.008
res@tmXBLabelDeltaF = -0.35
res@tmYLLabelDeltaF = -0.35
res@tmXBMajorOutwardLengthF = 0.011
res@tmYLMajorOutwardLengthF = 0.011
res@tmXBMinorOutwardLengthF = 0.008
res@tmYLMinorOutwardLengthF = 0.008





res@mpMinLatF = -60
res@mpMaxLatF = 90
res@mpMinLonF = 0
res@mpMaxLonF = 360
res@mpCenterLonF = 180-165-3

res@tmXBMode = "Explicit"
res@tmXBValues = (/-135,-90,-45,0,45,90,135,180/) + 165+3
res@tmXBLabels = (/"135W","90W","45W","0","45E","90E","135E","180"/)
res@tmXBMinorOn = True
res@tmXBMinorValues = fspan(0, 360, 25)+3

res@mpFillOn                    = False
res@mpOutlineOn                 = True
res@mpGeophysicalLineThicknessF = 1.5
res@mpGeophysicalLineColor      = "gray10";"gray20"
;res@mpGridAndLimbOn = True
;res@mpGridLatSpacingF = 90
;res@mpGridLonSpacingF = 360
;res@gsnMaximize = True
res@cnFillOn             = True
res@cnLinesOn            = False
;res@cnLineThicknessF     = 6.0
;res@cnLineColor          = "red"
res@cnLineLabelsOn       = False
res@lbLabelBarOn         = False
res@lbOrientation        = "Vertical"
res@pmLabelBarWidthF     = 0.045
res@pmLabelBarHeightF    = 0.14
res@pmLabelBarOrthogonalPosF = 0.012
res@lbLabelFontHeightF  = 0.015
;res@cnMissingValFillColor = "white"
;res@gsnYRefLine = 0.0
;res@gsnYRefLineColor = "black"
;res@gsnYRefLineThicknessF = 5.0
res@cnLevelSelectionMode = "ExplicitLevels"              
; res@cnLevels             = (/140.0,155.0,170.0,185.0,200.0,215.0,230.0,245.0/)
; res@cnFillColors         = (/16,15,13,12,8,6,5,3,-1/)
res@cnFillPalette        = "precip2_17lev" 



;vcres@mpFillDrawOrder         = "PostDraw"
vcres@vcRefAnnoOrthogonalPosF = -0.255
;vcres@vcRefAnnoSide           = "TopRight"
vcres@vcGlyphStyle            = "LineArrow"
vcres@vcRefAnnoArrowLineColor   = "black"         ; change ref vector color
vcres@vcRefAnnoArrowUseVecColor = False           ; don't use vec color for ref
vcres@vcMinDistanceF          = 0.0105             ; thin out vectors
vcres@vcLineArrowColor        = "black"           ; change vector color
vcres@vcRefAnnoOn             = True
vcres@vcRefLengthF            = 0.02            ; ref vec length
vcres@vcRefAnnoFontHeightF    = 0.012   ;参考箭头字体大小
vcres@vcRefAnnoString1On      = True
vcres@vcRefAnnoString2On      = False
vcres@vcLineArrowThicknessF   = 1.0            ; make vectors larger
vcres@vcVectorDrawOrder = "PostDraw"
vcres@vcRefAnnoOrthogonalPosF = -0.19 ;;正向下
vcres@vcRefAnnoParallelPosF   = 0.15  ;;正向右
vcres@vcRefAnnoString1        = "5"
vcres@vcRefMagnitudeF         = 5



sres = True
sres@cnLineLabelsOn   = False
sres@cnConstFLabelOn  = False
sres@cnInfoLabelOn    = False 
sres@gsnDraw          = False                   
sres@gsnFrame         = False
sres@gsnLeftString    = ""
sres@gsnRightString   = ""
sres@lbLabelBarOn = False
sres@cnFillOn  = False
sres@cnLinesOn = True
sres2 = sres
sres3 = sres
sres@cnLineColor = "red"
sres@cnLevelSelectionMode = "ExplicitLevels"
sres@cnLevels         = (/2700.0/)
sres@cnFillColors     = (/"red","gray20"/)
sres@cnMissingValFillColor= -1
sres@cnLineThicknessF = 3.0
sres@cnLineDashPattern = 14


sres2@cnLevelSelectionMode = "ExplicitLevels"
sres2@cnLevels         = (/2/)
sres2@cnLineColor = "blue"
sres2@cnLineDashPattern = 14
sres2@cnLineThicknessF = 4
; sres2@cnLineColors     = (/"red","red","red","red"/)
; sres2@cnLineDashPatterns = (/16,16,16,0,0,0,0/)
; sres2@cnLineThicknesses = (/4,4,4,4,4,4,4/)






lires = True
lires@tmBorderThicknessF = 2.2
lires@gsnDraw        = False
lires@gsnFrame       = False
lires@tiXAxisOn      = False
lires@tiYAxisOn      = False
lires@tmXTLabelsOn   = False
lires@tmXTOn         = False

lires@trYMinF        = 10
lires@trYMaxF        = 85

lires@tmXBLabelFontHeightF = 0.0225
lires@tmYLLabelFontHeightF = 0.0225
lires@tmYRLabelFontHeightF = 0.0225
lires@gsnLeftStringFontHeightF   = 0.022
lires@gsnRightStringFontHeightF  = 0.022
lires@tmYLLabelsOn   = False
lires@tmYLOn         = False

lires@tmYROn         = True
lires@tmYRLabelsOn   = True
lires@tmYRMode = "Explicit"
lires@tmYRValues = fspan(0, 80, 5)
lires@tmYRLabels = (/"0","20N","40N","60N","80N"/)
lires@tmYRMinorOn = True
lires@tmYRMinorValues = fspan(0, 80, 9)

lires@tmYLValues = fspan(0, 80, 5)
lires@tmYLLabels = (/"0","20N","40N","60N","80N"/)
lires@tmYLMinorOn = True
lires@tmYLMinorValues = fspan(0, 80, 9);lires@tmYLMode = "Explicit"
;lires@tmYLValues = fspan(1, 71, 8)
;lires@tmYLLabels = (/"1-April","11-April","21-April","1-May","11-May","21-May","31-May","10-June"/)

lires@xyMarkLineModes = (/"Lines","Lines","Lines"/)
;lires@xyMarkers = (/4,1,1/)
;lires@xyMarkerSizes = (/0.015,0.001,0.001/)
;lires@xyMarkerColors= (/"black","black","black"/)
;lires@xyMarkerThicknesses = (/4.0,2.0,1.0/)

lires@lgPerimOn = False
lires@lgLabelFontHeightF = 0.018
lires@vpWidthF  = 0.3
lires@vpHeightF = 0.26
lires@vpXF      = 0.06
lires@vpYF      = 0.75

lires@gsnLeftString  = ""
lires@gsnRightString = ""
lires@xyMarkLineModes = (/"Lines","Lines","Lines"/)
; lires@xyMarkers = (/1,1,1/)
; lires@xyMarkerSizes = (/0.015,0.001,0.001/)
; lires@xyMarkerColors= (/"black","black","black"/)
; lires@xyMarkerThicknesses = (/4.0,2.0,1.0/)
lires@xyLineThicknesses  = (/4.5,8.0,8.0/)
lires@xyLineColors       = (/"red","blue","orange"/)
lires@xyDashPatterns     = (/0,4,2/)

lires@tmYLMajorThicknessF = 2.5
lires@tmYRMajorThicknessF = 2.5
lires@tmXBMajorThicknessF = 2.5
lires@tmYLMinorThicknessF = 1.5
lires@tmYRMinorThicknessF = 1.5
lires@tmXBMinorThicknessF = 1.5
lires@tmYLMajorLengthF = 0.015
lires@tmYRMajorLengthF = 0.015
lires@tmXBMajorLengthF = 0.015
lires@tmYLMinorLengthF = 0.006
lires@tmYRMinorLengthF = 0.006
lires@tmXBMinorLengthF = 0.006
lires@tmXBLabelDeltaF = -0.25
lires@tmYLLabelDeltaF = -0.25
lires@tmXBMajorOutwardLengthF = 0.015
lires@tmYLMajorOutwardLengthF = 0.015
lires@tmYRMajorOutwardLengthF = 0.015
lires@tmXBMinorOutwardLengthF = 0.006
lires@tmYLMinorOutwardLengthF = 0.006
lires@tmYRMinorOutwardLengthF = 0.006






;************************* dotting ***************************
ores                 = True            
ores@gsnDraw         = False               
ores@gsnFrame        = False                 
ores@cnLineLabelsOn  = False               
ores@cnLinesOn       = False                 
ores@cnInfoLabelOn   = False                                                                                                        
ores@cnFillOn        = True                               
ores@lbLabelBarOn    = False                                                                                                            
ores@cnLevelSelectionMode = "ExplicitLevels"                                                                                                                                                                                
; ores@cnLevels        = (/14.007/)   ;;90%
; ores@cnLevels        = (/-0.312,0.312/)   ;;95%
ores@cnMonoFillPattern    = False            ; want multiple patterns                                                               
; ores@cnFillPatterns       = (/-1,11/)     ; the patterns                                                                         
ores@cnMonoFillScale      = False            ; want different densities                                                             
; ores@cnFillScales         = (/0.6,0.6/)    ; change densities                                                                         
ores@cnMonoFillColor      =True                                                                                                     
ores@cnFillDotSizeF       = 0.003    
ores@cnFillDrawOrder ="postDraw"
ores@cnFillColor = "gray30"

ores@cnLevels        = (/0.05/)   ;;90%
; ores@cnLevels        = (/-0.312,0.312/)   ;;95%
ores@cnMonoFillPattern    = False            ; want multiple patterns                                                               
ores@cnFillPatterns       = (/5,-1/)     ; the patterns                                                                         
ores@cnMonoFillScale      = False            ; want different densities                                                             
ores@cnFillScales         = (/0.8,0.8/)    ; change densities   




;**********************plot**************************
plot = new(14, "graphic")
vector = new(14, "graphic")
contour1 = new(14, "graphic")
contour2 = new(14, "graphic")
contour3 = new(14, "graphic")
contour4 = new(14, "graphic")
contour5 = new(14, "graphic")
contour6 = new(14, "graphic")
contour7 = new(14, "graphic")
contour8 = new(14, "graphic")
topooo = new(14, "graphic") 











;;---------------------------------------------------------------------------------------------------;;
pltType = "eps"
pltName = "/public/home/fcai/extreme5_degree/25_humid_heat/3_analysis1982_2023/Figure5_trend_1982_2023"
wks = gsn_open_wks(pltType, pltName)




res@mpProjection = "Robinson"
res@mpLimitMode = "LatLon" ; "Angles" ;
; res@mpBottomAngleF = 60
; res@mpTopAngleF = 90
res@mpMinLonF = 0
res@mpMaxLonF = 360
res@mpCenterLonF = 180-165-3
res@cnInfoLabelOn = False

res@cnFillMode = "CellFill"
res@cnFillPalette = "WhiteBlueGreenYellowRed" ; "WhiteBlue"
; res@cnLevels     := (/5,10,20,30,50,70/)
res@cnLevels     := (/1,2,3,4,5/)
res@cnFillColors := (/-1,35,70,120,170,235/)
res@tmXBLabelsOn  = False
res@mpGeophysicalLineThicknessF = 0.8
; plot(0) = gsn_csm_map(wks, res)
; contour1(0) = gsn_csm_contour(wks, degree, res)
; overlay(plot(0), contour1(0))



delete(res@cnFillMode)
res@cnFillPalette = "MPL_RdBu" ; "WhiteBlue"
res@cnLevels     := (/-2.0,-1.0,0,1.0,2.0/)
res@cnFillColors := (/100,85,73,53,40,30/)
res@tmXBLabelsOn  = True
res@mpGeophysicalLineThicknessF = 0.5



res@cnFillPalette = "WhiteBlueGreenYellowRed" ; "WhiteBlue"
res@cnLevels     := (/10.0,20.0,30.0,40.0,50.0/)
res@cnFillColors := (/15,65,105,140,175,210/)
plot(0) = gsn_csm_map(wks, res)
contour1(0) = gsn_csm_contour(wks, R2_detrend, res)
overlay(plot(0), contour1(0))


res@cnFillPalette = "WhiteBlueGreenYellowRed" ; "WhiteBlue"
res@cnLevels     := (/10.0,20.0,30.0,40.0,50.0/)
res@cnFillColors := (/15,65,105,140,175,210/)
plot(1) = gsn_csm_map(wks, res)
contour1(1) = gsn_csm_contour(wks, R2_LC_detrend, res)
overlay(plot(1), contour1(1))




res@cnFillPalette = "MPL_RdBu" ; "WhiteBlue"
res@cnLevels     := (/-2.0,-1.0,0,1.0,2.0/)
res@cnLevels     := (/-1.5,-0.5,0,0.5,1.5/)
; res@cnFillColors := (/100,85,73,53,40,30/)
res@cnFillColors := (/100,87,73,53,38,30/)
plot(2) = gsn_csm_map(wks, res)
contour1(2) = gsn_csm_contour(wks, wbt95_trend, res)
overlay(plot(2), contour1(2))

colors = (/"#33968E","#33968E"/) ;#C0C0C0"
vcres@vcLevelPalette = colors ;"MPL_BrBG" ; "WhiteBlue"
vcres@vcLevelSelectionMode = "ExplicitLevels"
vcres@vcLevels     := (/-999.0,-0.01,0,0.01,999/)
vcres@lbLabelBarOn = False
vector(2) = gsn_csm_vector_scalar(wks, qu_trend, qv_trend, qu_trend, vcres)
overlay(plot(2), vector(2))


res@cnFillPalette = "MPL_RdBu" ; "WhiteBlue"
res@cnLevels     := (/-4.0,-2,0,2,4/)
res@cnLevels     := (/-3.0,-1.5,0,1.5,3/)
; res@cnFillColors := (/100,85,73,53,40,30/)
plot(3) = gsn_csm_map(wks, res)
contour1(3) = gsn_csm_contour(wks, LC_trend, res)
overlay(plot(3), contour1(3))




res@cnFillPalette = "MPL_RdBu" ; "WhiteBlue"
; res@cnLevels     := (/-1.0,-0.5,0,0.5,1.0/)
; res@cnFillColors := (/100,85,73,53,40,25/)
res@cnLevels     := (/-1.5,-0.5,0,0.5,1.5/)
; res@cnFillColors := (/100,85,73,53,40,30/)

plot(4) = gsn_csm_map(wks, res)
contour1(4) = gsn_csm_contour(wks, estimated_trend, res)
overlay(plot(4), contour1(4))

ores@cnLevels        = (/0.3/)   ;;90%                                                            
ores@cnFillPatterns       = (/-1,17/)     ; the patterns                                                                                                                                   
ores@cnFillScales         = (/0.8,0.8/)    ; change densities 
ores@cnFillColor = "red4"
ores@cnFillColors = (/"gray30","red4"/) 
contour2(4) = gsn_csm_contour(wks, estimated_percentage, ores)
overlay(plot(4), contour2(4))



res@cnFillPalette = "MPL_RdBu" ; "WhiteBlue"
; res@cnLevels     := (/-2.0,-1,0,1,2/)
; res@cnFillColors := (/100,85,73,53,40,25/)
res@cnLevels     := (/-3.0,-1.5,0,1.5,3/)
; res@cnFillColors := (/100,85,73,53,40,30/)

plot(5) = gsn_csm_map(wks, res)
contour1(5) = gsn_csm_contour(wks, estimated_LC_trend, res)
overlay(plot(5), contour1(5))

ores@cnLevels        = (/0.3/)   ;;90%                                                            
ores@cnFillPatterns       = (/-1,17/)     ; the patterns                                                                                                                                   
ores@cnFillScales         = (/0.8,0.8/)    ; change densities 
ores@cnFillColor = "red4"
ores@cnFillColors = (/"gray30","red4"/) 
contour2(5) = gsn_csm_contour(wks, estimated_LC_percentage, ores)
overlay(plot(5), contour2(5))








lnres   =    True
lnres@gsLineThicknessF = 1.0
lnres@gsLineDashPattern=0

SAsia_lat1 = 8
SAsia_lat2 = 28
SAsia_lon1 = 70
SAsia_lon2 = 90
y1=SAsia_lat1
y2=SAsia_lat2
x1=SAsia_lon1
x2=SAsia_lon2
x=(/x1,x2,x2,x1,x1/)
y=(/y2,y2,y1,y1,y2/)
lnres@gsLineColor="green4"
  dum11 = gsn_add_polyline(wks,plot(4),x,y,lnres)
  dum21 = gsn_add_polyline(wks,plot(5),x,y,lnres)

WAsia_lat1 = 15
WAsia_lat2 = 35
WAsia_lon1 = 35
WAsia_lon2 = 55
y1=WAsia_lat1
y2=WAsia_lat2
x1=WAsia_lon1
x2=WAsia_lon2
x=(/x1,x2,x2,x1,x1/)
y=(/y2,y2,y1,y1,y2/)
  dum12 = gsn_add_polyline(wks,plot(4),x,y,lnres)
  dum22 = gsn_add_polyline(wks,plot(5),x,y,lnres)

NEAsia_lat1 = 40
NEAsia_lat2 = 65
NEAsia_lon1 = 120
NEAsia_lon2 = 140
y1=NEAsia_lat1
y2=NEAsia_lat2
x1=NEAsia_lon1
x2=NEAsia_lon2
x=(/x1,x2,x2,x1,x1/)
y=(/y2,y2,y1,y1,y2/)
  dum13 = gsn_add_polyline(wks,plot(4),x,y,lnres)
  dum23 = gsn_add_polyline(wks,plot(5),x,y,lnres)

SWAfrica_lat1 = -20
SWAfrica_lat2 = -5
SWAfrica_lon1 = 10
SWAfrica_lon2 = 40
y1=SWAfrica_lat1
y2=SWAfrica_lat2
x1=SWAfrica_lon1
x2=SWAfrica_lon2
x=(/x1,x2,x2,x1,x1/)
y=(/y2,y2,y1,y1,y2/)
  dum14 = gsn_add_polyline(wks,plot(4),x,y,lnres)
  dum24 = gsn_add_polyline(wks,plot(5),x,y,lnres)

EEur_lat1 = 45
EEur_lat2 = 70
EEur_lon1 = 20
EEur_lon2 = 40
y1=EEur_lat1
y2=EEur_lat2
x1=EEur_lon1
x2=EEur_lon2
x=(/x1,x2,x2,x1,x1/)
y=(/y2,y2,y1,y1,y2/)
  dum15 = gsn_add_polyline(wks,plot(4),x,y,lnres)
  dum25 = gsn_add_polyline(wks,plot(5),x,y,lnres)

Greenland_lat1 = 65
Greenland_lat2 = 85
Greenland_lon1 = -55 + 360
Greenland_lon2 = -20 + 360
y1=Greenland_lat1
y2=Greenland_lat2
x1=Greenland_lon1
x2=Greenland_lon2
x=(/x1,x2,x2,x1,x1/)
y=(/y2,y2,y1,y1,y2/)
  dum16 = gsn_add_polyline(wks,plot(4),x,y,lnres)
  dum26 = gsn_add_polyline(wks,plot(5),x,y,lnres)

NAmer_lat1 = 30
NAmer_lat2 = 40
NAmer_lon1 = -105 + 360
NAmer_lon2 = -85 + 360
y1=NAmer_lat1
y2=NAmer_lat2
x1=NAmer_lon1
x2=NAmer_lon2
x=(/x1,x2,x2,x1,x1/)
y=(/y2,y2,y1,y1,y2/)
  dum17 = gsn_add_polyline(wks,plot(4),x,y,lnres)
  dum27 = gsn_add_polyline(wks,plot(5),x,y,lnres)

SAmer_lat1 = 3
SAmer_lat2 = 13
SAmer_lon1 = -75 + 360
SAmer_lon2 = -55 + 360
y1=SAmer_lat1
y2=SAmer_lat2
x1=SAmer_lon1
x2=SAmer_lon2
x=(/x1,x2,x2,x1,x1/)
y=(/y2,y2,y1,y1,y2/)
  dum18 = gsn_add_polyline(wks,plot(4),x,y,lnres)
  dum28 = gsn_add_polyline(wks,plot(5),x,y,lnres)

WAfrica_lat1 = 5
WAfrica_lat2 = 15
WAfrica_lon1 = -15 + 360
WAfrica_lon2 = 20 + 360
y1=WAfrica_lat1
y2=WAfrica_lat2
x1=WAfrica_lon1
x2=WAfrica_lon2
x=(/x1,x2,x2,x1,x1/)
y=(/y2,y2,y1,y1,y2/)
  dum19 = gsn_add_polyline(wks,plot(4),x,y,lnres)
  dum29 = gsn_add_polyline(wks,plot(5),x,y,lnres)



WAsia_observed = dim_sum(dim_sum(wbt_trend_original({WAsia_lat1:WAsia_lat2},{WAsia_lon1:WAsia_lon2}))) / dim_sum(dim_sum(topo_01({WAsia_lat1:WAsia_lat2},{WAsia_lon1:WAsia_lon2})))
WAsia_estimated = dim_sum(dim_sum(estimated_trend({WAsia_lat1:WAsia_lat2},{WAsia_lon1:WAsia_lon2}))) / dim_sum(dim_sum(topo_01({WAsia_lat1:WAsia_lat2},{WAsia_lon1:WAsia_lon2})))
print(" WAsia (%) =  "+ WAsia_estimated/WAsia_observed )

SAsia_observed = dim_sum(dim_sum(wbt_trend_original({SAsia_lat1:SAsia_lat2},{SAsia_lon1:SAsia_lon2}))) / dim_sum(dim_sum(topo_01({SAsia_lat1:SAsia_lat2},{SAsia_lon1:SAsia_lon2})))
SAsia_estimated = dim_sum(dim_sum(estimated_trend({SAsia_lat1:SAsia_lat2},{SAsia_lon1:SAsia_lon2}))) / dim_sum(dim_sum(topo_01({SAsia_lat1:SAsia_lat2},{SAsia_lon1:SAsia_lon2})))
print(" SAsia (%) =  "+ SAsia_estimated/SAsia_observed +" "+SAsia_estimated+" "+SAsia_observed)

NEAsia_observed = dim_sum(dim_sum(wbt_trend_original({NEAsia_lat1:NEAsia_lat2},{NEAsia_lon1:NEAsia_lon2}))) / dim_sum(dim_sum(topo_01({NEAsia_lat1:NEAsia_lat2},{NEAsia_lon1:NEAsia_lon2})))
NEAsia_estimated = dim_sum(dim_sum(estimated_trend({NEAsia_lat1:NEAsia_lat2},{NEAsia_lon1:NEAsia_lon2}))) / dim_sum(dim_sum(topo_01({NEAsia_lat1:NEAsia_lat2},{NEAsia_lon1:NEAsia_lon2})))
print(" NEAsia (%) =  "+ NEAsia_estimated/NEAsia_observed  +" "+NEAsia_estimated+" "+NEAsia_observed)

SWAfrica_observed = dim_sum(dim_sum(wbt_trend_original({SWAfrica_lat1:SWAfrica_lat2},{SWAfrica_lon1:SWAfrica_lon2}))) / dim_sum(dim_sum(topo_01({SWAfrica_lat1:SWAfrica_lat2},{SWAfrica_lon1:SWAfrica_lon2})))
SWAfrica_estimated = dim_sum(dim_sum(estimated_trend({SWAfrica_lat1:SWAfrica_lat2},{SWAfrica_lon1:SWAfrica_lon2}))) / dim_sum(dim_sum(topo_01({SWAfrica_lat1:SWAfrica_lat2},{SWAfrica_lon1:SWAfrica_lon2})))
print(" SWAfrica (%) =  "+ SWAfrica_estimated/SWAfrica_observed )

EEur_observed = dim_sum(dim_sum(wbt_trend_original({EEur_lat1:EEur_lat2},{EEur_lon1:EEur_lon2}))) / dim_sum(dim_sum(topo_01({EEur_lat1:EEur_lat2},{EEur_lon1:EEur_lon2})))
EEur_estimated = dim_sum(dim_sum(estimated_trend({EEur_lat1:EEur_lat2},{EEur_lon1:EEur_lon2}))) / dim_sum(dim_sum(topo_01({EEur_lat1:EEur_lat2},{EEur_lon1:EEur_lon2})))
print(" EEur (%) =  "+ EEur_estimated/EEur_observed )

Greenland_observed = dim_sum(dim_sum(wbt_trend_original({Greenland_lat1:Greenland_lat2},{Greenland_lon1:Greenland_lon2}))) / dim_sum(dim_sum(topo_01({Greenland_lat1:Greenland_lat2},{Greenland_lon1:Greenland_lon2})))
Greenland_estimated = dim_sum(dim_sum(estimated_trend({Greenland_lat1:Greenland_lat2},{Greenland_lon1:Greenland_lon2}))) / dim_sum(dim_sum(topo_01({Greenland_lat1:Greenland_lat2},{Greenland_lon1:Greenland_lon2})))
print(" Greenland (%) =  "+ Greenland_estimated/Greenland_observed )

NAmer_observed = dim_sum(dim_sum(wbt_trend_original({NAmer_lat1:NAmer_lat2},{NAmer_lon1:NAmer_lon2}))) / dim_sum(dim_sum(topo_01({NAmer_lat1:NAmer_lat2},{NAmer_lon1:NAmer_lon2})))
NAmer_estimated = dim_sum(dim_sum(estimated_trend({NAmer_lat1:NAmer_lat2},{NAmer_lon1:NAmer_lon2}))) / dim_sum(dim_sum(topo_01({NAmer_lat1:NAmer_lat2},{NAmer_lon1:NAmer_lon2})))
print(" NAmer (%) =  "+ NAmer_estimated/NAmer_observed )

SAmer_observed = dim_sum(dim_sum(wbt_trend_original({SAmer_lat1:SAmer_lat2},{SAmer_lon1:SAmer_lon2}))) / dim_sum(dim_sum(topo_01({SAmer_lat1:SAmer_lat2},{SAmer_lon1:SAmer_lon2})))
SAmer_estimated = dim_sum(dim_sum(estimated_trend({SAmer_lat1:SAmer_lat2},{SAmer_lon1:SAmer_lon2}))) / dim_sum(dim_sum(topo_01({SAmer_lat1:SAmer_lat2},{SAmer_lon1:SAmer_lon2})))
print(" SAmer (%) =  "+ SAmer_estimated/SAmer_observed )

; WAfrica_lon1 = -15 + 360
; WAfrica_lon2 = 20 + 360
WAfrica_lon123 = (/1,3,5,7,9,11,13,15,17,19,345,347,349,351,353,355,357,359/)
WAfrica_observed = dim_sum(dim_sum(wbt_trend_original({WAfrica_lat1:WAfrica_lat2},{WAfrica_lon123}))) / dim_sum(dim_sum(topo_01({WAfrica_lat1:WAfrica_lat2},{WAfrica_lon123})))
WAfrica_estimated = dim_sum(dim_sum(estimated_trend({WAfrica_lat1:WAfrica_lat2},{WAfrica_lon123}))) / dim_sum(dim_sum(topo_01({WAfrica_lat1:WAfrica_lat2},{WAfrica_lon123})))
print(" WAfrica (%) =  "+ WAfrica_estimated/WAfrica_observed )

; printVarSummary(wbt_trend_original)
; printVarSummary(estimated_trend)







WAsia_observed = dim_sum(dim_sum(LC_trend_original({WAsia_lat1:WAsia_lat2},{WAsia_lon1:WAsia_lon2}))) / dim_sum(dim_sum(topo_01({WAsia_lat1:WAsia_lat2},{WAsia_lon1:WAsia_lon2})))
WAsia_estimated = dim_sum(dim_sum(estimated_LC_trend({WAsia_lat1:WAsia_lat2},{WAsia_lon1:WAsia_lon2}))) / dim_sum(dim_sum(topo_01({WAsia_lat1:WAsia_lat2},{WAsia_lon1:WAsia_lon2})))
print(" WAsia (% LC) =  "+ WAsia_estimated/WAsia_observed )

SAsia_observed = dim_sum(dim_sum(LC_trend_original({SAsia_lat1:SAsia_lat2},{SAsia_lon1:SAsia_lon2}))) / dim_sum(dim_sum(topo_01({SAsia_lat1:SAsia_lat2},{SAsia_lon1:SAsia_lon2})))
SAsia_estimated = dim_sum(dim_sum(estimated_LC_trend({SAsia_lat1:SAsia_lat2},{SAsia_lon1:SAsia_lon2}))) / dim_sum(dim_sum(topo_01({SAsia_lat1:SAsia_lat2},{SAsia_lon1:SAsia_lon2})))
print(" SAsia (% LC) =  "+ SAsia_estimated/SAsia_observed )

NEAsia_observed = dim_sum(dim_sum(LC_trend_original({NEAsia_lat1:NEAsia_lat2},{NEAsia_lon1:NEAsia_lon2}))) / dim_sum(dim_sum(topo_01({NEAsia_lat1:NEAsia_lat2},{NEAsia_lon1:NEAsia_lon2})))
NEAsia_estimated = dim_sum(dim_sum(estimated_LC_trend({NEAsia_lat1:NEAsia_lat2},{NEAsia_lon1:NEAsia_lon2}))) / dim_sum(dim_sum(topo_01({NEAsia_lat1:NEAsia_lat2},{NEAsia_lon1:NEAsia_lon2})))
print(" NEAsia (% LC) =  "+ NEAsia_estimated/NEAsia_observed )

SWAfrica_observed = dim_sum(dim_sum(LC_trend_original({SWAfrica_lat1:SWAfrica_lat2},{SWAfrica_lon1:SWAfrica_lon2}))) / dim_sum(dim_sum(topo_01({SWAfrica_lat1:SWAfrica_lat2},{SWAfrica_lon1:SWAfrica_lon2})))
SWAfrica_estimated = dim_sum(dim_sum(estimated_LC_trend({SWAfrica_lat1:SWAfrica_lat2},{SWAfrica_lon1:SWAfrica_lon2}))) / dim_sum(dim_sum(topo_01({SWAfrica_lat1:SWAfrica_lat2},{SWAfrica_lon1:SWAfrica_lon2})))
print(" SWAfrica (% LC) =  "+ SWAfrica_estimated/SWAfrica_observed )

EEur_observed = dim_sum(dim_sum(LC_trend_original({EEur_lat1:EEur_lat2},{EEur_lon1:EEur_lon2}))) / dim_sum(dim_sum(topo_01({EEur_lat1:EEur_lat2},{EEur_lon1:EEur_lon2})))
EEur_estimated = dim_sum(dim_sum(estimated_LC_trend({EEur_lat1:EEur_lat2},{EEur_lon1:EEur_lon2}))) / dim_sum(dim_sum(topo_01({EEur_lat1:EEur_lat2},{EEur_lon1:EEur_lon2})))
print(" EEur (% LC) =  "+ EEur_estimated/EEur_observed )

Greenland_observed = dim_sum(dim_sum(LC_trend_original({Greenland_lat1:Greenland_lat2},{Greenland_lon1:Greenland_lon2}))) / dim_sum(dim_sum(topo_01({Greenland_lat1:Greenland_lat2},{Greenland_lon1:Greenland_lon2})))
Greenland_estimated = dim_sum(dim_sum(estimated_LC_trend({Greenland_lat1:Greenland_lat2},{Greenland_lon1:Greenland_lon2}))) / dim_sum(dim_sum(topo_01({Greenland_lat1:Greenland_lat2},{Greenland_lon1:Greenland_lon2})))
print(" Greenland (% LC) =  "+ Greenland_estimated/Greenland_observed )

NAmer_observed = dim_sum(dim_sum(LC_trend_original({NAmer_lat1:NAmer_lat2},{NAmer_lon1:NAmer_lon2}))) / dim_sum(dim_sum(topo_01({NAmer_lat1:NAmer_lat2},{NAmer_lon1:NAmer_lon2})))
NAmer_estimated = dim_sum(dim_sum(estimated_LC_trend({NAmer_lat1:NAmer_lat2},{NAmer_lon1:NAmer_lon2}))) / dim_sum(dim_sum(topo_01({NAmer_lat1:NAmer_lat2},{NAmer_lon1:NAmer_lon2})))
print(" NAmer (% LC) =  "+ NAmer_estimated/NAmer_observed )

SAmer_observed = dim_sum(dim_sum(LC_trend_original({SAmer_lat1:SAmer_lat2},{SAmer_lon1:SAmer_lon2}))) / dim_sum(dim_sum(topo_01({SAmer_lat1:SAmer_lat2},{SAmer_lon1:SAmer_lon2})))
SAmer_estimated = dim_sum(dim_sum(estimated_LC_trend({SAmer_lat1:SAmer_lat2},{SAmer_lon1:SAmer_lon2}))) / dim_sum(dim_sum(topo_01({SAmer_lat1:SAmer_lat2},{SAmer_lon1:SAmer_lon2})))
print(" SAmer (% LC) =  "+ SAmer_estimated/SAmer_observed )

; WAfrica_lon1 = -15 + 360
; WAfrica_lon2 = 20 + 360
WAfrica_lon123 = (/1,3,5,7,9,11,13,15,17,19,345,347,349,351,353,355,357,359/)
WAfrica_observed = dim_sum(dim_sum(LC_trend_original({WAfrica_lat1:WAfrica_lat2},{WAfrica_lon123}))) / dim_sum(dim_sum(topo_01({WAfrica_lat1:WAfrica_lat2},{WAfrica_lon123})))
WAfrica_estimated = dim_sum(dim_sum(estimated_LC_trend({WAfrica_lat1:WAfrica_lat2},{WAfrica_lon123}))) / dim_sum(dim_sum(topo_01({WAfrica_lat1:WAfrica_lat2},{WAfrica_lon123})))
print(" WAfrica (% LC) =  "+ WAfrica_estimated/WAfrica_observed )









txres               = True
txres@txPerimOn = False
txres@txFontHeightF = 0.012
txres@txFontColor = "green4"
dtext0 = gsn_add_text(wks, plot(4), "40%", 54, 8.5, txres)
dtext1 = gsn_add_text(wks, plot(4), "41%", 85, 2, txres)
dtext2 = gsn_add_text(wks, plot(4), "57%", 157, 45, txres)
dtext3 = gsn_add_text(wks, plot(4), "44%", 5, -25, txres)
dtext4 = gsn_add_text(wks, plot(4), "75%", 54, 53, txres)
dtext5 = gsn_add_text(wks, plot(4), "47%", -25, 59, txres)
dtext6 = gsn_add_text(wks, plot(4), "53%", -72, 30, txres)
dtext7 = gsn_add_text(wks, plot(4), "35%", -43, 10, txres)
dtext8 = gsn_add_text(wks, plot(4), "36%", -8, -1, txres)

dtext10 = gsn_add_text(wks, plot(5), "58%", 54, 8.5, txres)
dtext11 = gsn_add_text(wks, plot(5), "66%", 82, 2, txres)
dtext12 = gsn_add_text(wks, plot(5), "57%", 157, 45, txres)
dtext13 = gsn_add_text(wks, plot(5), "56%", 5, -25, txres)
dtext14 = gsn_add_text(wks, plot(5), "97%", 54, 53, txres)
dtext15 = gsn_add_text(wks, plot(5), "64%", -25, 59, txres)
dtext16 = gsn_add_text(wks, plot(5), "56%", -72, 30, txres)
dtext17 = gsn_add_text(wks, plot(5), "52%", -43, 10, txres)
dtext18 = gsn_add_text(wks, plot(5), "62%", -8, -1, txres)




pres = True
pres@gsnFrame = False
pres@gsnPanelLabelBar = True
pres@lbOrientation        = "Vertical"
pres@pmLabelBarWidthF     = 0.03
pres@pmLabelBarHeightF    = 0.14
; pres@pmLabelBarParallelPosF = 0.11
pres@pmLabelBarOrthogonalPosF = 0.005
pres@lbLabelFontHeightF  = 0.009

pres@gsnPanelTop  = 0.91
pres@gsnPanelBottom = 0.61
pres@gsnPanelLeft  = 0.01
pres@gsnPanelRight = 0.53
pres@pmLabelBarOrthogonalPosF = -0.002
gsn_panel(wks,plot(0),(/1,1/),pres)

pres@gsnPanelLeft  = 0.47
pres@gsnPanelRight = 0.99
pres@pmLabelBarOrthogonalPosF = -0.002
gsn_panel(wks,plot(1),(/1,1/),pres)


pres@gsnPanelTop  = 0.65
pres@gsnPanelBottom = 0.35
pres@gsnPanelLeft  = 0.01
pres@gsnPanelRight = 0.53
pres@pmLabelBarOrthogonalPosF = -0.002
gsn_panel(wks,plot(2),(/1,1/),pres)

pres@gsnPanelLeft  = 0.47
pres@gsnPanelRight = 0.99
pres@pmLabelBarOrthogonalPosF = -0.002
gsn_panel(wks,plot(3),(/1,1/),pres)


pres@gsnPanelTop  = 0.39
pres@gsnPanelBottom = 0.09
pres@gsnPanelLeft  = 0.01
pres@gsnPanelRight = 0.53
pres@pmLabelBarOrthogonalPosF = -0.002
gsn_panel(wks,plot(4),(/1,1/),pres)

pres@gsnPanelLeft  = 0.47
pres@gsnPanelRight = 0.99
pres@pmLabelBarOrthogonalPosF = -0.002
gsn_panel(wks,plot(5),(/1,1/),pres)



txres               = True
txres@txFontHeightF = 0.0115
txres@txFont = "helvetica"
txres@txFontColor = "black"
gsn_text_ndc(wks, "Explained    variances    of    humid    heatwaves",  0.23, 0.87, txres)
gsn_text_ndc(wks, "Explained    variances    of    spatial    aggregation",  0.70, 0.87, txres)
gsn_text_ndc(wks, "Long-term    trend    in    humid    heatwaves",  0.23, 0.61, txres)
gsn_text_ndc(wks, "Long-term    trend    in    spatial    aggregation",  0.70, 0.61, txres)
gsn_text_ndc(wks, "Estimated    trend    in    humid    heatwaves",  0.23, 0.35, txres)
gsn_text_ndc(wks, "Estimated    trend    in    spatial    aggregation",  0.70, 0.35, txres)

txres@txFontHeightF = 0.0155
txres@txFont = "helvetica-bold"
gsn_text_ndc(wks, "a",  0.045, 0.895, txres)
gsn_text_ndc(wks, "c",  0.045, 0.63, txres)
gsn_text_ndc(wks, "e",  0.045, 0.37, txres)

gsn_text_ndc(wks, "b",  0.510, 0.895, txres)
gsn_text_ndc(wks, "d",  0.510, 0.63, txres)
gsn_text_ndc(wks, "f",  0.510, 0.37, txres)

frame(wks)
   delete(wks)  
   system("convert -geometry 1600x1600 -density 1200x1200 -trim " + pltName + "."+pltType + " " + pltName + ".png")
   system("rm " + pltName + "." + pltType + " -f")

end












