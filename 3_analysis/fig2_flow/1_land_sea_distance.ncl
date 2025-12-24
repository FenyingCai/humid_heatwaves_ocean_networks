load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
begin



;----------------------------------------------------------------  NH  ------------------------------------------------------------------------;;
f1 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/2_network/part15_wbt95_qs_NH_lag5/wbt95_qs95_lag5_NH2x2_1982_2023_sign99.nc","r")
networks0 = f1->networks0(:,:,:,:)
dim_0 = dimsizes(networks0)


f00 = addfile("/public/home/fcai/data0/topography/topo_2x2.nc","r")
topo = f00->topo({0:90},:)
topo_networks1 = conform(networks0, topo, (/0,1/))
topo_networks2 = conform(networks0, topo, (/2,3/))
delete(networks0)


;;--------------------------------------------  read distance  -------------------------------------------;;
f0 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/1_threshold_2x2/distance_i_j_2x2_NH.nc","r")
distance = f0->distance(:,:,:,:)
distance_to_land0 = where(topo_networks1.ge.0.0, distance, 1.0e20)
distance_to_land0@_FillValue = 1.0e20
distance_to_land = dim_min_n_Wrap(dim_min_n_Wrap(distance_to_land0, 1), 0)
distance_to_ocean0 = where(topo_networks2.lt.0.0, distance, 1.0e20)
distance_to_ocean0@_FillValue = 1.0e20
distance_to_ocean = dim_min_n_Wrap(dim_min_n_Wrap(distance_to_ocean0, 3), 2)

distance_landocean = where(topo.ge.0.0, distance_to_ocean, distance_to_land*(-1.0))
copy_VarCoords(topo, distance_landocean)
printMinMax(distance_landocean, 1)


land1 = where(distance_landocean.ge.0.0.and.distance_landocean.lt.500.0, 1.0, 0.0)
land2 = where(distance_landocean.ge.500.0.and.distance_landocean.lt.2600.0, 1.0, 0.0)
ocean1 = where(distance_landocean.le.0.0.and.distance_landocean.gt.-500.0, 1.0, 0.0)
ocean2 = where(distance_landocean.le.-500.0.and.distance_landocean.gt.-1800.0, 1.0, 0.0)
copy_VarCoords(topo, land1)
copy_VarCoords(topo, land2)
copy_VarCoords(topo, ocean1)
copy_VarCoords(topo, ocean2)


lat_weight2D = new((/dim_0(0),dim_0(1)/), float)
do i = 0,dim_0(0)-1
  lat_weight2D(i,:) = cos((i*2+1) /180.0 * 3.14159)
end do
area_land1 = dim_sum_n_Wrap(dim_sum_n_Wrap(land1*lat_weight2D, 1), 0)
area_land2 = dim_sum_n_Wrap(dim_sum_n_Wrap(land2*lat_weight2D, 1), 0)
area_ocean1 = dim_sum_n_Wrap(dim_sum_n_Wrap(ocean1*lat_weight2D, 1), 0)
area_ocean2 = dim_sum_n_Wrap(dim_sum_n_Wrap(ocean2*lat_weight2D, 1), 0)
print("area = "+area_ocean2+"  "+area_ocean1+"  "+area_land1+"  "+area_land2)


;;-------------------------------------------------------------------------------;;
output="/public/home/fcai/extreme5_degree/25_humid_heat/3_analysis1982_2023/fig2_flow/"
name  ="Distance_landocean_NH"

system("/bin/rm -f "+output+name+".nc")
ncdf=addfile(output+name+".nc","c")
ncdf->distance_landocean=distance_landocean
ncdf->land1=land1
ncdf->land2=land2
ncdf->ocean1=ocean1
ncdf->ocean2=ocean2







;----------------------------------------------------------------  SH  ------------------------------------------------------------------------;;
f1 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/2_network/part16_wbt95_qs_SH_lag5/wbt95_qs95_lag5_SH2x2_1982_2023_sign99.nc","r")
networks0 = f1->networks0(:,:,:,:)
dim_0 = dimsizes(networks0)


f00 = addfile("/public/home/fcai/data0/topography/topo_2x2.nc","r")
topo = f00->topo({-90:0},:)
topo_networks1 = conform(networks0, topo, (/0,1/))
topo_networks2 = conform(networks0, topo, (/2,3/))
delete(networks0)


;;--------------------------------------------  read distance  -------------------------------------------;;
f0 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/1_threshold_2x2/distance_i_j_2x2_SH.nc","r")
distance = f0->distance(:,:,:,:)
distance(0:14,:,:,:) = 1.0e20  ;; exclude 90S-60S
distance(:,:,0:14,:) = 1.0e20
distance@_FillValue = 1.0e20
distance_to_land0 = where(topo_networks1.ge.0.0, distance, 1.0e20)
distance_to_land0@_FillValue = 1.0e20
distance_to_land = dim_min_n_Wrap(dim_min_n_Wrap(distance_to_land0, 1), 0)
distance_to_ocean0 = where(topo_networks2.lt.0.0, distance, 1.0e20)
distance_to_ocean0@_FillValue = 1.0e20
distance_to_ocean = dim_min_n_Wrap(dim_min_n_Wrap(distance_to_ocean0, 3), 2)

distance_landocean = where(topo.ge.0.0, distance_to_ocean, distance_to_land*(-1.0))
copy_VarCoords(topo, distance_landocean)
printMinMax(distance_landocean, 1)


land1 = where(distance_landocean.ge.0.0.and.distance_landocean.lt.500.0, 1.0, 0.0)
land2 = where(distance_landocean.ge.500.0.and.distance_landocean.lt.2600.0, 1.0, 0.0)
ocean1 = where(distance_landocean.le.0.0.and.distance_landocean.gt.-500.0, 1.0, 0.0)
ocean2 = where(distance_landocean.le.-500.0.and.distance_landocean.gt.-1100.0, 1.0, 0.0)
copy_VarCoords(topo, land1)
copy_VarCoords(topo, land2)
copy_VarCoords(topo, ocean1)
copy_VarCoords(topo, ocean2)


lat_weight2D = new((/dim_0(0),dim_0(1)/), float)
do i = 0,dim_0(0)-1
  lat_weight2D(i,:) = cos((i*2-89) /180.0 * 3.14159)
end do
area_land1 = dim_sum_n_Wrap(dim_sum_n_Wrap(land1*lat_weight2D, 1), 0)  ;; 60S-0
area_land2 = dim_sum_n_Wrap(dim_sum_n_Wrap(land2*lat_weight2D, 1), 0)
area_ocean1 = dim_sum_n_Wrap(dim_sum_n_Wrap(ocean1*lat_weight2D, 1), 0)
area_ocean2 = dim_sum_n_Wrap(dim_sum_n_Wrap(ocean2*lat_weight2D, 1), 0)
print("area = "+area_ocean2+"  "+area_ocean1+"  "+area_land1+"  "+area_land2)


;;-------------------------------------------------------------------------------;;
output="/public/home/fcai/extreme5_degree/25_humid_heat/3_analysis1982_2023/fig2_flow/"
name  ="Distance_landocean_SH"

system("/bin/rm -f "+output+name+".nc")
ncdf=addfile(output+name+".nc","c")
ncdf->distance_landocean=distance_landocean
ncdf->land1=land1
ncdf->land2=land2
ncdf->ocean1=ocean1
ncdf->ocean2=ocean2



















; f1 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/3_analysis1982_2023/flow/Distance_landocean_NH.nc", "r")
; distance_landocean = f1->distance_landocean
; land1 = f1->land1
; land2 = f1->land2
; ocean1 = f1->ocean1
; ocean2 = f1->ocean2








; ;;---------------------------------------------------  plots  ---------------------------------------------------
; res = True
; res@tmBorderThicknessF = 1.4
; res@gsnDraw      =  False
; res@gsnFrame     =  False
; res@gsnAddCyclic =  True
; res@gsnRightString       = ""
; res@gsnLeftString        = ""
; res@gsnLeftStringFontHeightF   = 0.04
; res@gsnRightStringFontHeightF  = 0.04
; res@gsnCenterStringFontHeightF  = 0.04
; res@tmXTLabelsOn  = False
; res@tmYRLabelsOn  = False
; vcres = res
; cnres = res
; cnres@gsnAddCyclic =  False
; res@tmXTOn        = False
; res@tmYROn        = False
; res@tmXBOn        = True
; res@tmYLOn        = True
; ; res@tmXBLabelFontHeightF = 0.038
; ; res@tmYLLabelFontHeightF = 0.038
; res@tmXBTickSpacingF = 45.0
; res@tmYLTickSpacingF = 20.0





; res@mpMinLatF = -90
; res@mpMaxLatF = 90
; res@mpMinLonF = -180-47.5
; res@mpMaxLonF = 180-47.5
; res@mpCenterLonF = -47.5
; res@tmXBMode = "Explicit"
; res@tmXBValues = (/-225,-180,-135,-90,-45,0,45,90,135,180,225,270,315,360/)
; res@tmXBLabels = (/"135E","180","135W","90W","45W","0","45E","90E","135E","180","135W","90W","45W","0"/)
; res@tmYLMode = "Explicit"
; res@tmYLValues = (/-90,-60,-30,0,30,60,90/)
; res@tmYLLabels = (/"90S","60S","30S","0","30N","60N","90N"/)
; res@tmXBMinorOn = True
; res@tmXBMinorValues = fspan(-360, 360+15, 49+1)
; res@tmXBLabelsOn = True
; res@tmXBLabelFontHeightF = 0.021
; res@tmYLLabelsOn = True
; res@tmYLLabelFontHeightF = 0.021

; res@tmYLMajorThicknessF = 1.3
; res@tmXBMajorThicknessF = 1.3
; res@tmYLMinorThicknessF = 0.9
; res@tmXBMinorThicknessF = 0.9
; res@tmYLMajorLengthF = 0.014
; res@tmXBMajorLengthF = 0.014
; res@tmYLMinorLengthF = 0.01
; res@tmXBMinorLengthF = 0.01
; res@tmXBLabelDeltaF = -0.35
; res@tmYLLabelDeltaF = -0.35
; res@tmXBMajorOutwardLengthF = 0.014
; res@tmYLMajorOutwardLengthF = 0.014
; res@tmXBMinorOutwardLengthF = 0.01
; res@tmYLMinorOutwardLengthF = 0.01





; res@mpFillOn                    = False
; res@mpOutlineOn                 = True
; res@mpGeophysicalLineThicknessF = 1.5
; res@mpGeophysicalLineColor      = "gray45";"gray20"
; ;res@mpGridAndLimbOn = True
; ;res@mpGridLatSpacingF = 90
; ;res@mpGridLonSpacingF = 360
; ;res@gsnMaximize = True
; res@cnFillOn             = True
; res@cnLinesOn            = False
; ;res@cnLineThicknessF     = 6.0
; ;res@cnLineColor          = "red"
; res@cnLineLabelsOn       = False
; res@lbLabelBarOn         = False
; res@lbOrientation        = "Vertical"
; res@pmLabelBarWidthF     = 0.045
; res@pmLabelBarHeightF    = 0.14
; res@pmLabelBarOrthogonalPosF = 0.012
; res@lbLabelFontHeightF  = 0.015
; ;res@cnMissingValFillColor = "white"
; ;res@gsnYRefLine = 0.0
; ;res@gsnYRefLineColor = "black"
; ;res@gsnYRefLineThicknessF = 5.0
; res@cnMissingValFillColor = -1
; res@cnLevelSelectionMode = "ExplicitLevels"              
; ; res@cnLevels             = (/140.0,155.0,170.0,185.0,200.0,215.0,230.0,245.0/)
; ; res@cnFillColors         = (/16,15,13,12,8,6,5,3,-1/)
; res@cnFillPalette        = "precip2_17lev" 



; ;vcres@mpFillDrawOrder         = "PostDraw"
; vcres@vcRefAnnoOrthogonalPosF = -0.255
; ;vcres@vcRefAnnoSide           = "TopRight"
; vcres@vcGlyphStyle            = "LineArrow"
; vcres@vcRefAnnoArrowLineColor   = "black"         ; change ref vector color
; vcres@vcRefAnnoArrowUseVecColor = False           ; don't use vec color for ref
; vcres@vcMinDistanceF          = 0.03             ; thin out vectors
; vcres@vcLineArrowColor        = "black"           ; change vector color
; vcres@vcRefAnnoOn             = True
; vcres@vcRefLengthF            = 0.022            ; ref vec length
; vcres@vcRefAnnoFontHeightF    = 0.014   ;参考箭头字体大小
; vcres@vcRefAnnoString1On      = True
; vcres@vcRefAnnoString2On      = False
; vcres@vcLineArrowThicknessF   = 3.0            ; make vectors larger
; vcres@vcVectorDrawOrder = "PostDraw"
; vcres@vcRefAnnoOrthogonalPosF = -0.26 ;;正向下
; vcres@vcRefAnnoParallelPosF   = 1.0  ;;正向右
; vcres@vcRefAnnoString1        = "1"
; vcres@vcRefMagnitudeF         = 1



; sres = True
; sres@cnLineLabelsOn   = False
; sres@cnConstFLabelOn  = False
; sres@cnInfoLabelOn    = False 
; sres@gsnDraw          = False                   
; sres@gsnFrame         = False
; sres@gsnLeftString    = ""
; sres@gsnRightString   = ""
; sres@lbLabelBarOn = False
; sres@cnFillOn  = False
; sres@cnLinesOn = True
; sres2 = sres
; sres3 = sres
; sres@cnLineColor = "red"
; sres@cnLevelSelectionMode = "ExplicitLevels"
; ; sres@cnLevels         = (/2700.0/)
; ; sres@cnFillColors     = (/"red","gray20"/)
; sres@cnMissingValFillColor= -1
; sres@cnLineThicknessF = 6.0
; sres@cnLineDashPattern = 14


; sres2@cnLevelSelectionMode = "ExplicitLevels"
; sres2@cnLevels         = (/2/)
; sres2@cnLineColor = "blue"
; sres2@cnLineDashPattern = 14
; sres2@cnLineThicknessF = 4
; ; sres2@cnLineColors     = (/"red","red","red","red"/)
; ; sres2@cnLineDashPatterns = (/16,16,16,0,0,0,0/)
; ; sres2@cnLineThicknesses = (/4,4,4,4,4,4,4/)




; cnres@tmBorderThicknessF = 2.2
; cnres@gsnDraw        = False
; cnres@gsnFrame       = False
; cnres@lgPerimOn = False
; cnres@lgLabelFontHeightF = 0.018
; cnres@vpWidthF  = 0.42
; cnres@vpHeightF = 0.42
; cnres@vpXF      = 0.06
; cnres@vpYF      = 0.5

; cnres@trXMinF = -0.5
; cnres@trXMaxF = 18.5
; cnres@trYMinF = 0.0
; cnres@trYMaxF = 90.0


; cnres@tmXTOn        = False
; cnres@tmYROn        = False
; cnres@tmXBOn        = True
; cnres@tmYLOn        = True
; cnres@tmXBLabelFontHeightF = 0.014
; cnres@tmYLLabelFontHeightF = 0.014
; cnres@tmYRLabelFontHeightF = 0.014
; cnres@gsnLeftStringFontHeightF   = 0.022
; cnres@gsnRightStringFontHeightF  = 0.022
; cnres@gsnCenterStringFontHeightF  = 0.022
; cnres@tmYRLabelsOn   = False
; ; cnres@tmYROn         = True

; cnres@tmYLMode := "Explicit"
; cnres@tmYLValues := fspan(0.0, 80.0, 5)
; cnres@tmYLLabels := (/"0","20N","40N","60N","80N"/)
; cnres@tmYLMinorOn := True
; cnres@tmYLMinorValues := fspan(0.0, 80.0, 9)

; cnres@tmXBMode := "Explicit"
; cnres@tmXBValues := fspan(0, 18.0, 19)
; cnres@tmXBLabels := (/"M","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18"/)
; cnres@tmXBMinorOn := False

; cnres@cnFillOn       = True
; cnres@cnLinesOn      = False
; cnres@cnLineLabelsOn = False
; cnres@lbLabelBarOn   = False
; cnres@cnInfoLabelOn  = False
; cnres@cnLevelSelectionMode = "ExplicitLevels"

; cnres@tmYLMajorThicknessF = 2.5
; cnres@tmXBMajorThicknessF = 2.5
; cnres@tmYLMinorThicknessF = 1.5
; cnres@tmXBMinorThicknessF = 1.5
; cnres@tmYLMajorLengthF = 0.012
; cnres@tmXBMajorLengthF = 0.012
; cnres@tmYLMinorLengthF = 0.009
; cnres@tmXBMinorLengthF = 0.009
; cnres@tmXBLabelDeltaF = -0.4
; cnres@tmYLLabelDeltaF = -0.4
; cnres@tmXBMajorOutwardLengthF = 0.012
; cnres@tmYLMajorOutwardLengthF = 0.012
; cnres@tmXBMinorOutwardLengthF = 0.009
; cnres@tmYLMinorOutwardLengthF = 0.009







; lires = True
; lires@tmBorderThicknessF = 1.1
; lires@gsnDraw        = False
; lires@gsnFrame       = False
; lires@tiXAxisOn      = False
; lires@tiYAxisOn      = False
; lires@tmXTLabelsOn   = False
; lires@tmXTOn         = False

; lires@trYMinF        = 0.0
; lires@trYMaxF        = 0.5
; lires@trXMinF        = 0.5
; lires@trXMaxF        = 12.5

; lires@tmXBLabelFontHeightF = 0.072
; lires@tmYLLabelFontHeightF = 0.072
; lires@tmYRLabelFontHeightF = 0.072
; lires@gsnLeftStringFontHeightF   = 0.025
; lires@gsnRightStringFontHeightF  = 0.025
; lires@tmYRLabelsOn   = False
; lires@tmYROn         = False

; lires@tmYLOn         = True
; lires@tmYLLabelsOn   = True
; lires@tmYLMode = "Explicit"
; lires@tmXBMode = "Explicit"
; ; lires@tmYLValues = fspan(0, 80, 5)
; ; lires@tmYLLabels = (/"0","20N","40N","60N","80N"/)
; ; lires@tmYLMinorOn = True
; ; lires@tmYLMinorValues = fspan(0, 90, 10)

; ; lires@tmYLValues = fspan(0, 80, 5)
; ; lires@tmYLLabels = (/"0","20N","40N","60N","80N"/)
; ; lires@tmYLMinorOn = True
; ; lires@tmYLMinorValues = fspan(0, 90, 10)
; ;lires@tmYLMode = "Explicit"
; ;lires@tmYLValues = fspan(1, 71, 8)
; ;lires@tmYLLabels = (/"1-April","11-April","21-April","1-May","11-May","21-May","31-May","10-June"/)

; lires@xyMarkLineModes = (/"Lines","Lines","Lines"/)
; ;lires@xyMarkers = (/4,1,1/)
; ;lires@xyMarkerSizes = (/0.015,0.001,0.001/)
; ;lires@xyMarkerColors= (/"black","black","black"/)
; ;lires@xyMarkerThicknesses = (/4.0,2.0,1.0/)

; lires@lgPerimOn = False
; lires@lgLabelFontHeightF = 0.018
; lires@vpWidthF  = 0.99
; lires@vpHeightF = 0.75
; lires@vpXF      = 0.06
; lires@vpYF      = 0.75

; lires@gsnLeftString  = ""
; lires@gsnRightString = ""
; lires@xyMarkLineModes = (/"Lines","Lines","Lines"/)
; ; lires@xyMarkers = (/1,1,1/)
; ; lires@xyMarkerSizes = (/0.015,0.001,0.001/)
; ; lires@xyMarkerColors= (/"black","black","black"/)
; ; lires@xyMarkerThicknesses = (/4.0,2.0,1.0/)
; lires@xyLineThicknesses  = (/3.5,8.0,8.0/)
; lires@xyLineColors       = (/"red","blue","orange"/)
; lires@xyDashPatterns     = (/0,4,2/)

; lires@tmYLMajorThicknessF = 1.5
; lires@tmYRMajorThicknessF = 1.5
; lires@tmXBMajorThicknessF = 1.5
; lires@tmYLMinorThicknessF = 1.0
; lires@tmYRMinorThicknessF = 1.0
; lires@tmXBMinorThicknessF = 1.0
; lires@tmYLMajorLengthF = 0.06
; lires@tmYRMajorLengthF = 0.06
; lires@tmXBMajorLengthF = 0.06
; lires@tmYLMinorLengthF = 0.04
; lires@tmYRMinorLengthF = 0.04
; lires@tmXBMinorLengthF = 0.04
; lires@tmXBLabelDeltaF = -0.25
; lires@tmYLLabelDeltaF = -0.25
; lires@tmXBMajorOutwardLengthF = 0.06
; lires@tmYLMajorOutwardLengthF = 0.06
; lires@tmYRMajorOutwardLengthF = 0.06
; lires@tmXBMinorOutwardLengthF = 0.04
; lires@tmYLMinorOutwardLengthF = 0.04
; lires@tmYRMinorOutwardLengthF = 0.04






; ;;*************************打点***************************;;
; ores                 = True   
; ores@gsnAddCyclic    =  True         
; ores@gsnDraw         = False               
; ores@gsnFrame        = False                 
; ores@cnLineLabelsOn  = False               
; ores@cnLinesOn       = False                 
; ores@cnInfoLabelOn   = False                                                                                                        
; ores@cnFillOn        = True                               
; ores@lbLabelBarOn    = False                                                                                                            
; ores@cnLevelSelectionMode = "ExplicitLevels"                                                                                                                                                                                
; ; ores@cnLevels        = (/14.007/)   ;;90%
; ; ores@cnLevels        = (/-0.312,0.312/)   ;;95%
; ores@cnMonoFillPattern    = False            ; want multiple patterns                                                               
; ; ores@cnFillPatterns       = (/-1,11/)     ; the patterns                                                                         
; ores@cnMonoFillScale      = False            ; want different densities                                                             
; ; ores@cnFillScales         = (/0.6,0.6/)    ; change densities                                                                         
; ores@cnMonoFillColor      =True                                                                                                     
; ores@cnFillDotSizeF       = 0.003    
; ores@cnFillDrawOrder ="postDraw"
; ores@cnFillColor = "gray30"

; ores@cnLevels        = (/0.05/)   ;;90%
; ; ores@cnLevels        = (/-0.312,0.312/)   ;;95%
; ores@cnMonoFillPattern    = False            ; want multiple patterns                                                               
; ores@cnFillPatterns       = (/5,-1/)     ; the patterns                                                                         
; ores@cnMonoFillScale      = False            ; want different densities                                                             
; ores@cnFillScales         = (/0.8,0.8/)    ; change densities   




; ;**********************plot**************************
; plot = new(14, "graphic")
; plot2 = new(20, "graphic")
; plot3 = new(20, "graphic")
; plot4 = new(20, "graphic")
; vector = new(14, "graphic")
; contour1 = new(14, "graphic")
; contour2 = new(24, "graphic")
; contour3 = new(14, "graphic")
; contour4 = new(14, "graphic")
; contour5 = new(14, "graphic")
; contour6 = new(14, "graphic")
; contour7 = new(14, "graphic")
; contour8 = new(14, "graphic")
; topooo = new(14, "graphic") 




; ;;-------------------------------------------------------------------------------------------
; pltType = "eps"
; pltName = "/public/home/fcai/extreme5_degree/25_humid_heat/3_analysis1982_2023/Figure3_distance_landocean"
; wks = gsn_open_wks(pltType, pltName)


; colors = (/"#a21d2f","#2b2a76","#73aa43"/)

; res@cnMissingValFillColor = -1
; res@cnFillMode = "CellFill"
; res@gsnLeftString = ""
; res@gsnLeftStringFontHeightF = 0.022
; res@gsnCenterStringFontHeightF = 0.022
; res@tmXBLabelFontHeightF = 0.019
; res@tmYLLabelFontHeightF = 0.019

; ; res@mpMinLonF = -180-75
; ; res@mpMaxLonF = 180-75
; ; res@mpCenterLonF = -75
; res@tmXBLabelsOn = True
; res@cnFillPalette := colors
; res@cnLevels      = (/0.5/)
; res@cnFillColors  = (/-1,0/)
; plot(0) = gsn_csm_contour_map_ce(wks, land1, res)





; ; res@cnFillPalette = "WhiteBlue" 
; res@cnFillColors  = (/-1,1/)
; plot(1) = gsn_csm_contour_map_ce(wks, land2, res)
; res@tmXBLabelsOn = True
; ; res@cnFillPalette = "WhiteYellowOrangeRed" 
; res@cnFillColors  = (/-1,2/)
; plot(2) = gsn_csm_contour_map_ce(wks, ocean1, res)
; plot(3) = gsn_csm_contour_map_ce(wks, ocean2, res)




; ;;--------------------------------------------------------;;
; ; txres = True
; ; txres@txFont = "helvetica"
; ; txres@txFontHeightF = 0.022
; ; txres@txFontColor = "#a21d2f"
; ; eur5 = gsn_add_text(wks,plot(0), "Pacific", -27, -45, txres)
; ; txres@txFontColor = "#2b2a76"
; ; eur6 = gsn_add_text(wks,plot(1), "Atlantic", 53, -45, txres)
; ; txres@txFontColor = "#73aa43"
; ; eur7 = gsn_add_text(wks,plot(2), "Indian", -13, -28, txres)
; ; eur8 = gsn_add_text(wks,plot(2), "Ocean", -13, -45, txres)




;   pres = True
;   pres@gsnFrame = False
;   pres@gsnPanelLabelBar = False

;   pres@gsnPanelTop  = 0.95
;   pres@gsnPanelBottom = 0.05
;   pres@gsnPanelLeft  = 0.05
;   pres@gsnPanelRight = 0.95
;   gsn_panel(wks,plot(0:3),(/2,2/),pres)







; txres               = True
; txres@txFont = "helvetica"
; txres@txFontColor = "black"
; txres@txFontHeightF = 0.0175
; gsn_text_ndc(wks, "Four   types   of   grid   cells",  0.515, 0.775, txres)


; txres@txFontHeightF = 0.0195
; txres@txFont = "helvetica-bold"
; ; gsn_text_ndc(wks, "a",  0.20, 0.95, txres)
; ; gsn_text_ndc(wks, "b",  0.20, 0.645, txres)
; ; gsn_text_ndc(wks, "c",  0.20, 0.34, txres)



; frame(wks)
;    delete(wks)  
;    system("convert -geometry 2400x2400 -density 1800x1800 -trim " + pltName + "."+pltType + " " + pltName + ".png")     
;   system("rm " + pltName + "." + pltType + " -f")



end


