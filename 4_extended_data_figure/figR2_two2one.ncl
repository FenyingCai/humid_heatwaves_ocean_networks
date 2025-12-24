load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
begin


;;-------------------------------- read data ----------------------------------;;
f1 = addfile("/public/data2/model/fcai/data0/topography/topo_adapt_era5_air2m.nc", "r")
topo0 = (f1->topo({-90:90},:))




tmax_yearly = new((/42,181,360/), float)
do iyear = 0,41
print("iyear = "+iyear)

tmax_NH0 = new((/92,91,360/), float)
f1 = addfile("/public/data2/model/fcai/data-observation/ERA5_daily/t2m_max/t2m."+(iyear+1982)+"-06.daily.nc", "r")
tmax_NH0(0:29,:,:) = short2flt(f1->t2m(:,{0:90},:))
f2 = addfile("/public/data2/model/fcai/data-observation/ERA5_daily/t2m_max/t2m."+(iyear+1982)+"-07.daily.nc", "r")
tmax_NH0(30:60,:,:) = short2flt(f2->t2m(:,{0:90},:))
f3 = addfile("/public/data2/model/fcai/data-observation/ERA5_daily/t2m_max/t2m."+(iyear+1982)+"-08.daily.nc", "r")
tmax_NH0(61:91,:,:) = short2flt(f3->t2m(:,{0:90},:))
tmax_yearly(iyear,90:180,:) = dim_avg_n_Wrap(tmax_NH0, 0)


tmax_SH0 = new((/92,90,360/), float)
f1 = addfile("/public/data2/model/fcai/data-observation/ERA5_daily/t2m_max/t2m."+(iyear+1982)+"-01.daily.nc", "r")
tmax_SH0(0:30,:,:) = short2flt(f1->t2m(:,{-90:-1},:))
f2 = addfile("/public/data2/model/fcai/data-observation/ERA5_daily/t2m_max/t2m."+(iyear+1982)+"-02.daily.nc", "r")
tmax_SH0(31:58,:,:) = short2flt(f2->t2m(0:27,{-90:-1},:))
f3 = addfile("/public/data2/model/fcai/data-observation/ERA5_daily/t2m_max/t2m."+(iyear+1982)+"-12.daily.nc", "r")
tmax_SH0(59:89,:,:) = short2flt(f3->t2m(:,{-90:-1},:))
tmax_yearly(iyear,0:89,:) = dim_avg_n_Wrap(tmax_SH0, 0)
end do
tmax_clim = dim_avg_n_Wrap(tmax_yearly, 0)
tmax_clim = tmax_clim - 273.15
tmax_clim = where(topo0.ge.0.0, tmax_clim, 1.0e20)
tmax_clim@_FillValue = 1.0e20
copy_VarCoords(topo0, tmax_clim)




q_yearly = new((/42,181,360/), float)
do iyear = 0,41
print("iyear = "+iyear)

q_NH0 = new((/92,91,360/), float)
f1 = addfile("/public/data2/model/fcai/data0/ERA5_q2m_daily1x1/q2m."+(iyear+1982)+".nc", "r")
q_NH0(0:91,:,:) = f1->q2m(151:151+91,{0:90},:)
q_yearly(iyear,90:180,:) = dim_avg_n_Wrap(q_NH0, 0)


q_SH0 = new((/92,90,360/), float)
q_SH0(0:58,:,:) = f1->q2m(0:58,{-90:-1},:)
q_SH0(59:89,:,:) = f1->q2m(334:364,{-90:-1},:)
q_yearly(iyear,0:89,:) = dim_avg_n_Wrap(q_SH0, 0)
end do
q_clim = dim_avg_n_Wrap(q_yearly, 0)
q_clim = q_clim * 1000.0
q_clim = where(topo0.ge.0.0, q_clim, 1.0e20)
q_clim@_FillValue = 1.0e20
copy_VarCoords(topo0, q_clim)













;;----------------------------------------------------  read networks (NH)  ----------------------------------------------------------;;
f1 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/2_network/part21_wbt95_tmax95_NH/wbt95_tmax95_NH2x2_1982_2023_sign99.nc","r")
networks1 = (f1->networks0)
networks1 = where(networks1.gt.0.0, networks1, 0.0)
printVarSummary(networks1)
printMinMax(networks1, 1)
dim_0 := dimsizes(networks1)

;;------------- remove ocean ---------------;;
f00 = addfile("/public/data2/model/fcai/data0/topography/topo_2x2.nc","r")
topo := f00->topo({0:90},:)
topo_networks1 = conform(networks1, topo, (/0,1/))
topo_networks2 = conform(networks1, topo, (/2,3/))

networks1 = where(topo_networks1.ge.0.0, networks1, 1.0e20)
networks1 = where(topo_networks2.ge.0.0, networks1, 1.0e20)
network_ref = where(topo_networks1.ge.0.0, 1.0, 1.0e20)
network_ref = where(topo_networks2.ge.0.0, network_ref, 1.0e20)
delete([/topo_networks1,topo_networks2/])
networks1@_FillValue = 1.0e20
printMinMax(networks1, 0)


f0 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/1_threshold_2x2/distance_i_j_2x2_NH.nc","r")
distance = f0->distance(:,:,:,:)
networks1 = where(distance.le.2000.0, networks1, 1.0e20)
delete(distance)

;;------------------------  lat weighted  ------------------------;;
lat_matrix = new((/dim_0(0),dim_0(1), dim_0(0),dim_0(1)/), float)
do i = 0,dim_0(0)-1
  do j = 0,dim_0(0)-1
    lat_matrix(i,:,j,:) = cos((i*2+1+j*2+1) /2.0 /180.0 * 3.14159)
  end do
end do

networks1__ = where(networks1.gt.0.0, 1.0, 0.0)
ref1 = dim_sum_n_Wrap(dim_sum_n_Wrap(networks1__(:,:,:,:)*network_ref, 1), 0) 
ref2 = dim_sum_n_Wrap(dim_sum_n_Wrap(networks1__(:,:,:,:)*network_ref, 3), 2) 
ref1 = where(ref1.gt.0.0, ref1, 1.0e20)
ref2 = where(ref2.gt.0.0, ref2, 1.0e20)
ref1@_FillValue = 1.0e20
ref2@_FillValue = 1.0e20
degree_land = dim_sum_n_Wrap(dim_sum_n_Wrap(networks1(:,:,:,:)*lat_matrix, 3), 2) / ref2
delete(networks1)
degree_NH = where(topo.ge.0.0, degree_land, 0.0)



;;----------------------------------------------------  read networks (SH)  ----------------------------------------------------------;;
f1 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/2_network/part22_wbt95_tmax95_SH/wbt95_tmax95_SH2x2_1982_2023_sign99.nc","r")
networks1 = (f1->networks0)
networks1(0:14,:,:,:) = 0.0
networks1(:,:,0:14,:) = 0.0
networks1 = where(networks1.gt.0.0, networks1, 0.0)
printVarSummary(networks1)

;;------------- remove ocean ---------------;;
f00 = addfile("/public/data2/model/fcai/data0/topography/topo_2x2.nc","r")
topo = f00->topo({-90:0},:)
topo_networks1 = conform(networks1, topo, (/0,1/))
topo_networks2 = conform(networks1, topo, (/2,3/))

networks1 = where(topo_networks1.ge.0.0, networks1, 1.0e20)
networks1 = where(topo_networks2.ge.0.0, networks1, 1.0e20)
network_ref = where(topo_networks1.ge.0.0, 1.0, 1.0e20)
network_ref = where(topo_networks2.ge.0.0, network_ref, 1.0e20)
delete([/topo_networks1,topo_networks2/])
networks1@_FillValue = 1.0e20
printMinMax(networks1, 0)


f0 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/1_threshold_2x2/distance_i_j_2x2_SH.nc","r")
distance = f0->distance(:,:,:,:)
networks1 = where(distance.le.2000.0, networks1, 1.0e20)
delete(distance)

;;------------------------  lat weighted  ------------------------;;
lat_matrix = new((/dim_0(0),dim_0(1), dim_0(0),dim_0(1)/), float)
do i = 0,dim_0(0)-1
  do j = 0,dim_0(0)-1
    lat_matrix(i,:,j,:) = cos((i*2-89+j*2-89) /2.0 /180.0 * 3.14159)
  end do
end do

networks1__ = where(networks1.gt.0.0, 1.0, 0.0)
ref1 = dim_sum_n_Wrap(dim_sum_n_Wrap(networks1__(:,:,:,:)*network_ref, 1), 0) 
ref2 = dim_sum_n_Wrap(dim_sum_n_Wrap(networks1__(:,:,:,:)*network_ref, 3), 2) 
ref1 = where(ref1.gt.0.0, ref1, 1.0e20)
ref2 = where(ref2.gt.0.0, ref2, 1.0e20)
ref1@_FillValue = 1.0e20
ref2@_FillValue = 1.0e20
degree_land = dim_sum_n_Wrap(dim_sum_n_Wrap(networks1(:,:,:,:)*lat_matrix, 3), 2) / ref2

degree_SH = where(topo.ge.0.0, degree_land, 0.0)

;;---------------------  Degree (NH + SH)  ----------------------;;
f00 = addfile("/public/data2/model/fcai/data0/topography/topo_2x2.nc","r")
topo := f00->topo({-90:90},:)
dim_0 := dimsizes(topo)
degree_HWs_Tmax = new((/dim_0(0),dim_0(1)/), float)
degree_HWs_Tmax(0:44,:) = degree_SH * 1.0
degree_HWs_Tmax(45:89,:) = degree_NH * 1.0
copy_VarCoords(topo, degree_HWs_Tmax)
printMinMax(degree_HWs_Tmax,1)











;;---------------------------------------------------  plots  -----------------------------------------------------
res = True
res@tmBorderThicknessF = 0.8
res@gsnDraw      =  False
res@gsnFrame     =  False
res@gsnAddCyclic =  True
res@gsnRightString       = ""
res@gsnLeftString        = ""
res@gsnLeftStringFontHeightF   = 0.04
res@gsnRightStringFontHeightF  = 0.04
res@gsnCenterStringFontHeightF  = 0.04
res@tmXTLabelsOn  = False
res@tmYRLabelsOn  = False
vcres = res
cnres = res
cnres@gsnAddCyclic =  False
res@tmXTOn        = False
res@tmYROn        = False
res@tmXBOn        = True
res@tmYLOn        = True
; res@tmXBLabelFontHeightF = 0.038
; res@tmYLLabelFontHeightF = 0.038
res@tmXBTickSpacingF = 45.0
res@tmYLTickSpacingF = 20.0





res@mpMinLatF = -60
res@mpMaxLatF = 90
res@mpMinLonF = 0
res@mpMaxLonF = 360
res@mpCenterLonF = 180-25
res@tmXBMode = "Explicit"
res@tmXBValues = (/-135,-90,-45,0,45,90,135,180,225,270,315/) + 25
res@tmXBLabels = (/"135W","90W","45W","0","45E","90E","135E","180","135W","90W","45W"/)
res@tmXBMinorOn = True
res@tmXBMinorValues = fspan(-360, 360, 49) + 25
res@tmYLMode = "Explicit"
res@tmYLValues = (/-90,-60,-30,0,30,60,90/)
res@tmYLLabels = (/"90S","60S","30S","0","30N","60N","90N"/)
res@tmYLMinorValues = fspan(-90, 90, 19)
res@tmXBLabelsOn = True
res@tmXBLabelFontHeightF = 0.017
res@tmYLLabelsOn = True
res@tmYLLabelFontHeightF = 0.017

res@tmYLMajorThicknessF = 1.1
res@tmXBMajorThicknessF = 1.1
res@tmYLMinorThicknessF = 0.8
res@tmXBMinorThicknessF = 0.8
res@tmYLMajorLengthF = 0.015
res@tmXBMajorLengthF = 0.015
res@tmYLMinorLengthF = 0.011
res@tmXBMinorLengthF = 0.011
res@tmXBLabelDeltaF = -0.35
res@tmYLLabelDeltaF = -0.35
res@tmXBMajorOutwardLengthF = 0.015
res@tmYLMajorOutwardLengthF = 0.015
res@tmXBMinorOutwardLengthF = 0.011
res@tmYLMinorOutwardLengthF = 0.011





res@mpFillOn                    = False
res@mpOutlineOn                 = True
res@mpGeophysicalLineThicknessF = 0.6
res@mpGeophysicalLineColor      = "black";"gray20"
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
res@cnMissingValFillColor = -1
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
vcres@vcMinDistanceF          = 0.03             ; thin out vectors
vcres@vcLineArrowColor        = "black"           ; change vector color
vcres@vcRefAnnoOn             = True
vcres@vcRefLengthF            = 0.022            ; ref vec length
vcres@vcRefAnnoFontHeightF    = 0.014   ;参考箭头字体大小
vcres@vcRefAnnoString1On      = True
vcres@vcRefAnnoString2On      = False
vcres@vcLineArrowThicknessF   = 3.0            ; make vectors larger
vcres@vcVectorDrawOrder = "PostDraw"
vcres@vcRefAnnoOrthogonalPosF = -0.26 ;;正向下
vcres@vcRefAnnoParallelPosF   = 1.0  ;;正向右
vcres@vcRefAnnoString1        = "1"
vcres@vcRefMagnitudeF         = 1



sres = True
; sres@gsnAddCyclic =  True
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
; sres@cnLineColor = "red"
sres@cnLevelSelectionMode = "ExplicitLevels"
; sres@cnLevels         = (/2700.0/)
; sres@cnFillColors     = (/"red","gray20"/)
sres@cnMissingValFillColor= -1
; sres@cnLineThicknessF = 6.0
; sres@cnLineDashPattern = 14


sres2@cnLevelSelectionMode = "ExplicitLevels"
sres2@cnLevels         = (/2/)
sres2@cnLineColor = "blue"
sres2@cnLineDashPattern = 14
sres2@cnLineThicknessF = 4
; sres2@cnLineColors     = (/"red","red","red","red"/)
; sres2@cnLineDashPatterns = (/16,16,16,0,0,0,0/)
; sres2@cnLineThicknesses = (/4,4,4,4,4,4,4/)




cnres@tmBorderThicknessF = 2.2
cnres@gsnDraw        = False
cnres@gsnFrame       = False
cnres@lgPerimOn = False
cnres@lgLabelFontHeightF = 0.018
cnres@vpWidthF  = 0.80
cnres@vpHeightF = 0.32
cnres@vpXF      = 0.06
cnres@vpYF      = 0.5

cnres@trXMinF = -60.0
cnres@trXMaxF = 60.0
cnres@trYMinF = 0.0
cnres@trYMaxF = 85.0


cnres@tmXTOn        = False
cnres@tmYROn        = False
cnres@tmXBOn        = True
cnres@tmYLOn        = True
cnres@tmXBLabelFontHeightF = 0.0215
cnres@tmYLLabelFontHeightF = 0.0215
cnres@tmYRLabelFontHeightF = 0.0215
cnres@gsnLeftStringFontHeightF   = 0.022
cnres@gsnRightStringFontHeightF  = 0.022
cnres@gsnCenterStringFontHeightF  = 0.02
cnres@tmYRLabelsOn   = False
; cnres@tmYROn         = True


cnres@tmYLMode = "Explicit"
cnres@tmYLValues = fspan(1980.0, 2020.0, 5)
cnres@tmYLLabels = fspan(1980.0, 2020.0, 5)
cnres@tmYLMinorOn = True
cnres@tmYLMinorValues = fspan(1980.0, 2020.0, 9)

cnres@tmXBMode = "Explicit"
cnres@tmXBValues = fspan(-60, 120.0, 7)
cnres@tmXBLabels = fspan(-60, 120.0, 7)
cnres@tmXBMinorOn = True
cnres@tmXBMinorValues = fspan(-75, 120.0, 14)

cnres@cnFillOn       = True
cnres@cnLinesOn      = False
cnres@cnLineLabelsOn = False
cnres@lbLabelBarOn   = False
cnres@cnInfoLabelOn  = False
cnres@cnLevelSelectionMode = "ExplicitLevels"

cnres@tmYLMajorThicknessF = 2.5
cnres@tmXBMajorThicknessF = 2.5
cnres@tmYLMinorThicknessF = 1.5
cnres@tmXBMinorThicknessF = 1.5
cnres@tmYLMajorLengthF = 0.014
cnres@tmXBMajorLengthF = 0.014
cnres@tmYLMinorLengthF = 0.01
cnres@tmXBMinorLengthF = 0.01
cnres@tmXBLabelDeltaF = -0.4
cnres@tmYLLabelDeltaF = -0.4
cnres@tmXBMajorOutwardLengthF = 0.014
cnres@tmYLMajorOutwardLengthF = 0.014
cnres@tmXBMinorOutwardLengthF = 0.01
cnres@tmYLMinorOutwardLengthF = 0.01






lires = True
lires@tmBorderThicknessF = 1.5
lires@gsnDraw        = False
lires@gsnFrame       = False
lires@tiXAxisOn      = False
lires@tiYAxisOn      = False
lires@tmXTLabelsOn   = False
lires@tmXTOn         = False

lires@trYMinF        = 0.0
lires@trYMaxF        = 0.5
lires@trXMinF        = 0.5
lires@trXMaxF        = 12.5

lires@tmXBLabelFontHeightF = 0.042
lires@tmXTLabelFontHeightF = 0.042
lires@tmYLLabelFontHeightF = 0.042
lires@tmYRLabelFontHeightF = 0.042
lires@gsnLeftStringFontHeightF   = 0.04
lires@gsnRightStringFontHeightF  = 0.04
lires@tmYRLabelsOn   = False
lires@tmYROn         = False

lires@tmYLOn         = True
lires@tmYLLabelsOn   = True
lires@tmYLMode = "Explicit"
lires@tmXBMode = "Explicit"
; lires@tmYLValues = fspan(0, 80, 5)
; lires@tmYLLabels = (/"0","20N","40N","60N","80N"/)
; lires@tmYLMinorOn = True
; lires@tmYLMinorValues = fspan(0, 90, 10)

; lires@tmYLValues = fspan(0, 80, 5)
; lires@tmYLLabels = (/"0","20N","40N","60N","80N"/)
; lires@tmYLMinorOn = True
; lires@tmYLMinorValues = fspan(0, 90, 10)
;lires@tmYLMode = "Explicit"
;lires@tmYLValues = fspan(1, 71, 8)
;lires@tmYLLabels = (/"1-April","11-April","21-April","1-May","11-May","21-May","31-May","10-June"/)

lires@xyMarkLineModes = (/"Lines","Lines","Lines"/)
;lires@xyMarkers = (/4,1,1/)
;lires@xyMarkerSizes = (/0.015,0.001,0.001/)
;lires@xyMarkerColors= (/"black","black","black"/)
;lires@xyMarkerThicknesses = (/4.0,2.0,1.0/)

lires@lgPerimOn = False
lires@lgLabelFontHeightF = 0.018
lires@vpWidthF  = 0.8
lires@vpHeightF = 0.7
lires@vpXF      = 0.06
lires@vpYF      = 0.75

lires@gsnLeftString  = ""
lires@gsnRightString = ""
lires@xyMarkLineModes = (/"Lines","Lines","Lines"/)
; lires@xyMarkers = (/1,1,1/)
; lires@xyMarkerSizes = (/0.015,0.001,0.001/)
; lires@xyMarkerColors= (/"black","black","black"/)
; lires@xyMarkerThicknesses = (/4.0,2.0,1.0/)
lires@xyLineThicknesses  = (/3.5,8.0,8.0/)
lires@xyLineColors       = (/"red","blue","orange"/)
lires@xyDashPatterns     = (/0,4,2/)

lires@tmYLMajorThicknessF = 2.5
lires@tmXTMajorThicknessF = 2.5
lires@tmXBMajorThicknessF = 2.5
lires@tmYLMinorThicknessF = 1.5
lires@tmXTMinorThicknessF = 1.5
lires@tmXBMinorThicknessF = 1.5
lires@tmYLMajorLengthF = 0.038
lires@tmXTMajorLengthF = 0.038
lires@tmXBMajorLengthF = 0.038
lires@tmYLMinorLengthF = 0.028
lires@tmXTMinorLengthF = 0.028
lires@tmXBMinorLengthF = 0.028
lires@tmXBLabelDeltaF = -0.25
lires@tmXTLabelDeltaF = -0.25
lires@tmYLLabelDeltaF = -0.25
lires@tmXBMajorOutwardLengthF = 0.038
lires@tmYLMajorOutwardLengthF = 0.038
lires@tmXTMajorOutwardLengthF = 0.038
lires@tmXBMinorOutwardLengthF = 0.028
lires@tmYLMinorOutwardLengthF = 0.028
lires@tmXTMinorOutwardLengthF = 0.028






;*************************打点***************************
ores                 = True   
ores@gsnAddCyclic    =  True         
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
ores@cnFillDrawOrder ="postDraw"
ores@cnFillColor = "white"

ores@cnLevels        = (/0.05/)   ;;90%
; ores@cnLevels        = (/-0.312,0.312/)   ;;95%
ores@cnMonoFillPattern    = False            ; want multiple patterns                                                               
ores@cnFillPatterns       = (/-1,6/)     ; the patterns                                                                         
ores@cnMonoFillScale      = False            ; want different densities                                                             
ores@cnFillScales         = (/0.45,0.45/)    ; change densities   
ores@cnFillDotSizeF       = 0.0025



;**********************plot**************************
plot = new(14, "graphic")
plot2 = new(20, "graphic")
plot3 = new(20, "graphic")
plot4 = new(20, "graphic")
vector = new(14, "graphic")
contour1 = new(14, "graphic")
contour2 = new(24, "graphic")
contour3 = new(14, "graphic")
contour4 = new(14, "graphic")
contour5 = new(14, "graphic")
contour6 = new(14, "graphic")
contour7 = new(14, "graphic")
contour8 = new(14, "graphic")
topooo = new(14, "graphic") 






;;-------------------------------------------------------------------------------------------
pltType = "eps"
pltName = "/public/home/fcai/extreme5_degree/25_humid_heat/7_response/FigureR2_two2one"
wks = gsn_open_wks(pltType, pltName)





res@cnFillPalette = "WhiteBlueGreenYellowRed" ; "WhiteBlue"
res@cnLevels     := (/14,18,22,26,30,34/)
res@cnFillColors := (/-1,42,84,126,168,210,240/)
res@tmXBLabelsOn  = True
res@mpGeophysicalLineThicknessF = 0.3
res@mpGeophysicalLineColor = "black"
plot(0) = gsn_csm_contour_map_ce(wks, tmax_clim, res)




res@cnFillPalette = "WhiteBlueGreenYellowRed" ; "WhiteBlue"
res@cnLevels     := (/3,6,9,12,15,18/)
res@cnFillColors := (/-1,42,84,126,168,210,235/)
res@tmXBLabelsOn  = True
res@tmYLLabelsOn  = False
; res@mpGeophysicalLineThicknessF = 0.8
plot(1) = gsn_csm_contour_map_ce(wks, q_clim, res)




res@mpGeophysicalLineThicknessF = 0.3
res@mpGeophysicalLineColor = "black"
res@cnFillPalette = "WhiteBlueGreenYellowRed" ; "WhiteBlue"
res@cnLevels     := (/17,22,27,32,37,42/)
res@cnFillColors := (/-1,42,84,126,168,210,235/)
res@tmXBLabelsOn  = True
res@tmYLLabelsOn  = True
res@cnLevels     := (/17,22,27,32,37,42/)
plot(2) = gsn_csm_contour_map_ce(wks, degree_HWs_Tmax, res)






  pres = True
  pres@gsnFrame = False
  pres@gsnPanelLabelBar = True
  pres@lbOrientation        = "Vertical"
  pres@pmLabelBarWidthF     = 0.035
  pres@pmLabelBarHeightF    = 0.15
  pres@pmLabelBarOrthogonalPosF = 0.003
  ; pres@pmLabelBarParallelPosF = 0.01
  pres@lbLabelFontHeightF  = 0.008

  pres@gsnPanelTop  = 0.71
  pres@gsnPanelBottom = 0.51
  pres@gsnPanelLeft  = 0.01
  pres@gsnPanelRight = 0.525
  gsn_panel(wks,plot(0),(/1,1/),pres)

  pres@gsnPanelLeft  = 0.495
  pres@gsnPanelRight = 0.99
  gsn_panel(wks,plot(1),(/1,1/),pres)


  pres@gsnPanelTop  = 0.46
  pres@gsnPanelBottom = 0.26
  pres@gsnPanelLeft  = 0.242
  pres@gsnPanelRight = 0.758
  gsn_panel(wks,plot(2),(/1,1/),pres)  



txres               = True
txres@txFontColor = "black"
txres@txFontHeightF = 0.0105
txres@txFont = "helvetica"
gsn_text_ndc(wks, "Near-surface  air  temperature",  0.235, 0.715, txres)
gsn_text_ndc(wks, "Near-surface  specific  humidity",  0.697, 0.715, txres)
gsn_text_ndc(wks, "Local   concurrences   with   heat   events",  0.47, 0.465, txres)

txres@txFontHeightF = 0.0135
txres@txFont = "helvetica-bold"
gsn_text_ndc(wks, "a",  0.03, 0.73, txres)
gsn_text_ndc(wks, "b",  0.51, 0.73, txres)
gsn_text_ndc(wks, "c",  0.26, 0.48, txres)



frame(wks)
   delete(wks)  
   system("convert -geometry 2400x2400 -density 1600x1600 -trim " + pltName + "."+pltType + " " + pltName + ".png")     
  system("rm " + pltName + "." + pltType + " -f")

end







