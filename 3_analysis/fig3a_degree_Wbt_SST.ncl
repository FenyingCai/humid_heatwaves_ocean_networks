load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
begin
; ;;-----------------  read lat, lon  ----------------;;
lat=fspan(1,89.0,45)
lon=fspan(1,359.0,180)
lat@units="degrees_north"
lon@units="degrees_east"
dim_0 = (/45,180, 45,180/)





;;-------------  read networks (NH)  ------------;;
f1 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/2_network/part1_wbt95_sst95_NH/wbt95_sst95_NH2x2_1982_2023_sign99.nc","r")
networks1 = (f1->networks0)
networks1 = where(networks1.gt.0.0, 1.0, 0.0)
printVarSummary(networks1)

;;------------- remove ocean ---------------;;
f00 = addfile("/public/data2/model/fcai/data0/topography/topo_2x2.nc","r")
topo = f00->topo({0:90},:)
topo_networks1 = conform(networks1, topo, (/0,1/))
topo_networks2 = conform(networks1, topo, (/2,3/))

networks1 = where(topo_networks1.ge.0.0, networks1, 1.0e20)
networks1 = where(topo_networks2.le.0.0, networks1, 1.0e20)
network_ref = where(topo_networks1.ge.0.0, 1.0, 1.0e20)
network_ref = where(topo_networks2.le.0.0, network_ref, 1.0e20)
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

ref1 = dim_sum_n_Wrap(dim_sum_n_Wrap(networks1(:,:,:,:)*network_ref, 1), 0) 
ref2 = dim_sum_n_Wrap(dim_sum_n_Wrap(networks1(:,:,:,:)*network_ref, 3), 2) 
ref1 = where(ref1.gt.0.0, ref1, 1.0e20)
ref2 = where(ref2.gt.0.0, ref2, 1.0e20)
ref1@_FillValue = 1.0e20
ref2@_FillValue = 1.0e20
degree_ocean = dim_sum_n_Wrap(dim_sum_n_Wrap(networks1(:,:,:,:)*lat_matrix, 1), 0) / ref1
degree_land = dim_sum_n_Wrap(dim_sum_n_Wrap(networks1(:,:,:,:)*lat_matrix, 3), 2) / ref2
delete(networks1)

degree_ocean!0 = "lat"
degree_ocean!1 = "lon"
degree_ocean&lat = lat
degree_ocean&lon = lon
printVarSummary(degree_ocean)
copy_VarCoords(degree_ocean, degree_land)

degree_NH = where(topo.ge.0.0, degree_land, degree_ocean)
copy_VarCoords(degree_ocean, degree_NH)






;;-------------  read networks (SH)  ------------;;
f1 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/2_network/part2_wbt95_sst95_SH/wbt95_sst95_SH2x2_1982_2023_sign99.nc","r")
networks1 = (f1->networks0)
networks1(0:14,:,:,:) = 0.0
networks1(:,:,0:14,:) = 0.0
networks1 = where(networks1.gt.0.0, 1.0, 0.0)
printVarSummary(networks1)

;;------------- remove ocean ---------------;;
f00 = addfile("/public/data2/model/fcai/data0/topography/topo_2x2.nc","r")
topo = f00->topo({-90:0},:)
topo_networks1 = conform(networks1, topo, (/0,1/))
topo_networks2 = conform(networks1, topo, (/2,3/))

networks1 = where(topo_networks1.ge.0.0, networks1, 1.0e20)
networks1 = where(topo_networks2.le.0.0, networks1, 1.0e20)
network_ref = where(topo_networks1.ge.0.0, 1.0, 1.0e20)
network_ref = where(topo_networks2.le.0.0, network_ref, 1.0e20)
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

ref1 = dim_sum_n_Wrap(dim_sum_n_Wrap(networks1(:,:,:,:)*network_ref, 1), 0) 
ref2 = dim_sum_n_Wrap(dim_sum_n_Wrap(networks1(:,:,:,:)*network_ref, 3), 2) 
ref1 = where(ref1.gt.0.0, ref1, 1.0e20)
ref2 = where(ref2.gt.0.0, ref2, 1.0e20)
ref1@_FillValue = 1.0e20
ref2@_FillValue = 1.0e20
degree_ocean = dim_sum_n_Wrap(dim_sum_n_Wrap(networks1(:,:,:,:)*lat_matrix, 1), 0) / ref1 
degree_land = dim_sum_n_Wrap(dim_sum_n_Wrap(networks1(:,:,:,:)*lat_matrix, 3), 2) / ref2

degree_ocean!0 = "lat"
degree_ocean!1 = "lon"
degree_ocean&lat = lat
degree_ocean&lon = lon
printVarSummary(degree_ocean)
copy_VarCoords(degree_ocean, degree_land)

degree_SH = where(topo.ge.0.0, degree_land, degree_ocean)
copy_VarCoords(degree_ocean, degree_SH)





;;-------------  Degree (NH + SH)  ------------;;
f00 = addfile("/public/data2/model/fcai/data0/topography/topo_2x2.nc","r")
topo := f00->topo({-90:90},:)
dim_0 := dimsizes(topo)
degree = new((/dim_0(0),dim_0(1)/), float)
degree(0:44,:) = degree_SH * 100.0
degree(45:89,:) = degree_NH * 100.0
copy_VarCoords(topo, degree)
printMinMax(degree,1)


















;;------------------------------------------- read composite -------------------------------------------;;
f1 = addfile("/public/data2/model/fcai/data0/topography/topo_adapt_era5_air2m.nc", "r")
topo0 := (f1->topo({-90:90},:))
dim_0 := dimsizes(topo0)
lat := topo0&lat
lon := topo0&lon

f1 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/6_response/composite_analysis/SAsia_composite.nc","r")
q_SAsia_diff=f1->q_diff
u_SAsia_diff=f1->u_diff
v_SAsia_diff=f1->v_diff
qu_SAsia_diff=f1->qu_diff
qv_SAsia_diff=f1->qv_diff
sst_SAsia_diff=f1->sst_diff
wbt_SAsia_diff=f1->wbt_diff

q_SAsia_t=f1->q_t
u_SAsia_t=f1->u_t
v_SAsia_t=f1->v_t
qu_SAsia_t=f1->qu_t
qv_SAsia_t=f1->qv_t
sst_SAsia_t=f1->sst_t
wbt_SAsia_t=f1->wbt_t

qu_SAsia_diff = qu_SAsia_diff*1000.0
qv_SAsia_diff = qv_SAsia_diff*1000.0
qu_SAsia_diff = where(qu_SAsia_t.le.0.05.or.qv_SAsia_t.le.0.05, qu_SAsia_diff, 1.0e20)
qv_SAsia_diff = where(qu_SAsia_t.le.0.05.or.qv_SAsia_t.le.0.05, qv_SAsia_diff, 1.0e20)
qu_SAsia_diff@_FillValue = 1.0e20
qv_SAsia_diff@_FillValue = 1.0e20
copy_VarCoords(u_SAsia_diff, qu_SAsia_diff)
copy_VarCoords(v_SAsia_diff, qv_SAsia_diff)
q_SAsia_diff = q_SAsia_diff*1000.0
copy_VarCoords(u_SAsia_diff, q_SAsia_diff)

; q_SAsia_diff = where(q_SAsia_t.le.0.05, q_SAsia_diff, 1.0e20)
q_SAsia_diff@_FillValue = 1.0e20
copy_VarCoords(q_SAsia_t, q_SAsia_diff)
; sst_SAsia_diff = where(sst_SAsia_t.le.0.05, sst_SAsia_diff, 1.0e20)
; wbt_SAsia_diff = where(wbt_SAsia_t.le.0.05, wbt_SAsia_diff, 1.0e20)
sst_SAsia_diff@_FillValue = 1.0e20
wbt_SAsia_diff@_FillValue = 1.0e20
copy_VarCoords(sst_SAsia_t, sst_SAsia_diff)
copy_VarCoords(sst_SAsia_t, wbt_SAsia_diff)






;;---------------------------------------------------;;
SA_lat1 = 8
SA_lat2 = 28
SA_lon1 = 70
SA_lon2 = 90
key_region3 = topo0
key_region3 = 0.0
key_region3({SA_lat1:SA_lat2},{SA_lon1:SA_lon2}) = 1.0
key_region3 = where(topo0.ge.0.0, key_region3, -1.0)
copy_VarCoords(topo0, key_region3)





;************************************************** plots **********************************************
res = True
res@tmBorderThicknessF = 0.9
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
res@tmXBOn        = True
res@tmYLOn        = True
res@tmXBLabelFontHeightF = 0.015
res@tmYLLabelFontHeightF = 0.015
res@tmXBTickSpacingF = 45.0
res@tmYLTickSpacingF = 20.0




res@tmYLMajorThicknessF = 1.0
res@tmXBMajorThicknessF = 1.0
res@tmYLMinorThicknessF = 0.9
res@tmXBMinorThicknessF = 0.9
res@tmXBMajorLengthF = 0.014
res@tmYLMinorLengthF = 0.014
res@tmXBMinorLengthF = 0.010
res@tmYLMinorLengthF = 0.010
res@tmXBLabelDeltaF = -0.35
res@tmYLLabelDeltaF = -0.35
res@tmXBMajorOutwardLengthF = 0.014
res@tmYLMajorOutwardLengthF = 0.014
res@tmXBMinorOutwardLengthF = 0.010
res@tmYLMinorOutwardLengthF = 0.010





res@mpMinLatF = -60
res@mpMaxLatF = 90
res@mpMinLonF = 0
res@mpMaxLonF = 360
res@mpCenterLonF = 180-165-0

res@tmXBMode = "Explicit"
res@tmXBValues = (/-135,-90,-45,0,45,90,135,180/) + 165+0
res@tmXBLabels = (/"135W","90W","45W","0","45E","90E","135E","180"/)
res@tmYLMode = "Explicit"
res@tmYLValues = (/-90,-60,-30,0,30,60,90/)
res@tmYLLabels = (/"90S","60S","30S","0","30N","60N","90N"/)
res@tmXBMinorOn = True
res@tmXBMinorValues = fspan(0, 360, 25) + 0

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
vcres@vcMinDistanceF          = 0.03             ; thin out vectors
vcres@vcLineArrowColor        = "black"           ; change vector color
vcres@vcRefAnnoOn             = True
vcres@vcRefLengthF            = 0.022            ; ref vec length
vcres@vcRefAnnoFontHeightF    = 0.018   ;参考箭头字体大小
vcres@vcRefAnnoString1On      = True
vcres@vcRefAnnoString2On      = False
vcres@vcLineArrowThicknessF   = 1.2            ; make vectors larger
vcres@vcVectorDrawOrder = "PostDraw"
vcres@vcRefAnnoOrthogonalPosF = -0.19 ;;正向下
vcres@vcRefAnnoParallelPosF   = 1.0  ;;正向右
vcres@vcRefAnnoString1        = "10"
vcres@vcRefMagnitudeF         = 10





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

lires@tmYLMajorThicknessF = 1.5
lires@tmYRMajorThicknessF = 1.5
lires@tmXBMajorThicknessF = 1.5
lires@tmYLMinorThicknessF = 1.0
lires@tmYRMinorThicknessF = 1.0
lires@tmXBMinorThicknessF = 1.0
lires@tmYLMajorLengthF = 0.018
lires@tmYRMajorLengthF = 0.018
lires@tmXBMajorLengthF = 0.018
lires@tmYLMinorLengthF = 0.009
lires@tmYRMinorLengthF = 0.009
lires@tmXBMinorLengthF = 0.009
lires@tmXBLabelDeltaF = -0.25
lires@tmYLLabelDeltaF = -0.25
lires@tmXBMajorOutwardLengthF = 0.018
lires@tmYLMajorOutwardLengthF = 0.018
lires@tmYRMajorOutwardLengthF = 0.018
lires@tmXBMinorOutwardLengthF = 0.009
lires@tmYLMinorOutwardLengthF = 0.009
lires@tmYRMinorOutwardLengthF = 0.009






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

ores1 = ores
ores1@cnFillColor = "black"
ores1@cnLevels        = (/0.05/)   ;;95%
ores1@cnFillPatterns  = (/3,-1/)     ; the patterns   
 
ores2 = ores
ores2@cnFillColor = "white"
ores2@cnLevels        = (/0.05/)   ;;95%, 0.294
ores2@cnFillPatterns  = (/3,-1/)     ; the patterns  





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











;;-----------------------------------------------------------------------------------------;;
pltType = "eps"
pltName = "/public/home/fcai/extreme5_degree/25_humid_heat/3_analysis1982_2023/Figure3a_degree_Wbt_SST"
wks = gsn_open_wks(pltType, pltName)




res@cnFillMode = "CellFill"
res@cnFillPalette = "WhiteBlueGreenYellowRed" ; "WhiteBlue"
res@cnLevels     := (/10,30,50,70,90/)
res@cnFillColors := (/-1,50,100,150,180,230/)
res@tmXBLabelsOn  = True
res@mpGeophysicalLineThicknessF = 0.8

plot(0) = gsn_csm_contour_map_ce(wks, degree, res)
txres               = True
txres@txPerimOn = False
txres@txFontHeightF = 0.015
; txres@txFont = "helvetica"
txres@txFontColor = "black"
dtext1 = gsn_add_text(wks, plot(0), "< 2,000 km", -136, -48, txres)






res@tmXBValues := (/-150,-120,-90,-60,-30,0,30,60,90,120,150,180/)
res@tmXBLabels := (/"150W","120W","90W","60W","30W","0","30E","60E","90E","120E","150E","180"/)
res@tmYLValues := (/0,20,40/)
res@tmYLLabels := (/"0","20N","40N"/)
res@tmXBLabelsOn = True
res@tmXBLabelFontHeightF = 0.022
res@tmYLLabelsOn = True
res@tmYLLabelFontHeightF = 0.022

res@tmYLMajorThicknessF = 1.2
res@tmXBMajorThicknessF = 1.2
res@tmYLMinorThicknessF = 1.0
res@tmXBMinorThicknessF = 1.0
res@tmXBMajorLengthF = 0.018
res@tmYLMinorLengthF = 0.018
res@tmXBMinorLengthF = 0.015
res@tmYLMinorLengthF = 0.015
res@tmXBLabelDeltaF = -0.35
res@tmYLLabelDeltaF = -0.35
res@tmXBMajorOutwardLengthF = 0.018
res@tmYLMajorOutwardLengthF = 0.018
res@tmXBMinorOutwardLengthF = 0.015
res@tmYLMinorOutwardLengthF = 0.015

res@mpMinLonF = 30
res@mpMaxLonF = 135
res@mpMinLatF = -10
res@mpMaxLatF = 50
res@mpCenterLonF = 0

res@cnFillPalette = "MPL_RdBu"
res@cnLevels     := (/-0.5,-0.3,-0.1,0.1,0.3,0.5/)
res@cnFillColors := (/100,88,77,-1,48,37,25/)
res@cnFillPalette = "BlueWhiteOrangeRed"
res@cnLevels     := (/-0.5,-0.3,-0.1,0.1,0.3,0.5/)
res@cnFillColors := (/50,70,95,-1,160,185,205/)
plot(3) = gsn_csm_contour_map_ce(wks, sst_SAsia_diff, res)
contour1(3) = gsn_csm_contour(wks, sst_SAsia_t, ores1)
overlay(plot(3), contour1(3))
keyres = res
keyres@cnFillPalette = "GMT_gray"
keyres@cnLevels     := (/0.5/)
keyres@cnFillColors := (/-1,7/)
keyres@cnInfoLabelOn = False
contour1(3) = gsn_csm_contour(wks, key_region3, keyres)
overlay(plot(3), contour1(3))



res@cnFillPalette = "MPL_RdBu"
res@cnLevels     := (/-1.0,-0.5,-0.1,0.1,0.5,1.0/)
res@cnFillColors := (/100,88,77,-1,48,37,27/)
plot(4) = gsn_csm_contour_map_ce(wks, wbt_SAsia_diff, res)
contour1(4) = gsn_csm_contour(wks, wbt_SAsia_t, ores1)
overlay(plot(4), contour1(4))

res@cnFillPalette = "MPL_BrBG"
res@cnLevels     := (/-1.0,-0.5,-0.1,0.1,0.5,1.0/)
res@cnFillColors := (/28,37,50,-1,75,89,98/)
plot(5) = gsn_csm_contour_map_ce(wks, q_SAsia_diff({850},:,:), res)
contour1(5) = gsn_csm_contour(wks, q_SAsia_t({850},:,:), ores2)
overlay(plot(5), contour1(5))
vector(5) = gsn_csm_vector(wks, qu_SAsia_diff({850},:,:), qv_SAsia_diff({850},:,:), vcres)
overlay(plot(5), vector(5))



txres               = True
txres@txFont = "helvetica-bold"
txres@txPerimOn = False
txres@txFontHeightF = 0.035
txres@txFontColor = "red4"
dtext1 = gsn_add_text(wks, plot(5), "C", 65, 17.5, txres)






pres = True
pres@gsnFrame = False
pres@gsnPanelLabelBar = True
pres@lbOrientation        = "Vertical"
pres@pmLabelBarWidthF     = 0.03
pres@pmLabelBarHeightF    = 0.16
; pres@pmLabelBarParallelPosF = 0.11
pres@pmLabelBarOrthogonalPosF = 0.002
pres@lbLabelFontHeightF  = 0.008

pres@gsnPanelTop  = 0.73
pres@gsnPanelBottom = 0.5
pres@gsnPanelLeft  = 0.01
pres@gsnPanelRight = 0.555
pres@amJust  ="BottomRight"
pres@gsnPanelFigureStringsFontHeightF = 0.0115
gsn_panel(wks,plot(0),(/1,1/),pres)




pres@pmLabelBarHeightF = 0.13
pres@pmLabelBarOrthogonalPosF = 0.004
pres@gsnPanelTop  = 0.435
pres@gsnPanelBottom = 0.235
pres@gsnPanelLeft  = 0.01
pres@gsnPanelRight = 0.34
gsn_panel(wks,plot(3),(/1,1/),pres)
pres@gsnPanelLeft  = 0.34
pres@gsnPanelRight = 0.67
gsn_panel(wks,plot(4),(/1,1/),pres)
pres@gsnPanelLeft  = 0.67
pres@gsnPanelRight = 1.00
pres@pmLabelBarOrthogonalPosF = 0.007
gsn_panel(wks,plot(5),(/1,1/),pres)





txres               = True
txres@txFontHeightF = 0.011
txres@txFont = "helvetica"
txres@txFontColor = "black"
gsn_text_ndc(wks, "Linkage    between    humid    heatwaves    and    adjacent    warm    oceans",  0.50, 0.765, txres)
gsn_text_ndc(wks, "Composite    analysis    of    humid    heatwaves    over    South    Asia",  0.50, 0.47, txres)

txres@txFontHeightF = 0.0095
gsn_text_ndc(wks, "SST",  0.16, 0.425, txres)
gsn_text_ndc(wks, "Wet-bulb  temperature",  0.49, 0.425, txres)
gsn_text_ndc(wks, "q850  &  quv850",  0.82, 0.425, txres)


txres@txFontHeightF = 0.0130
txres@txFont = "helvetica-bold"
gsn_text_ndc(wks, "a",  0.03, 0.735, txres)
gsn_text_ndc(wks, "b",  0.56, 0.735, txres)
gsn_text_ndc(wks, "c",  0.03, 0.44, txres)
gsn_text_ndc(wks, "d",  0.36, 0.44, txres)
gsn_text_ndc(wks, "e",  0.69, 0.44, txres)


frame(wks)
   delete(wks)  
   system("convert -geometry 2400x2400 -density 1600x1600 -trim " + pltName + "."+pltType + " " + pltName + ".png")
   system("rm " + pltName + "." + pltType + " -f")

end








