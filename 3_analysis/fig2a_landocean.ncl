load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
begin


;;------------------------------------------------------------------------------------------------------------------;;
f00 = addfile("/public/data2/model/fcai/data0/topography/topo_2x2.nc","r")
topo = f00->topo({-90:90},:)

land1 = new((/90,180/), float)
land2 = new((/90,180/), float)
ocean1 = new((/90,180/), float)
ocean2 = new((/90,180/), float)

f1 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/3_analysis1982_2023/flow/Distance_landocean_NH.nc", "r")
land1(45:89,:) = f1->land1
land2(45:89,:) = f1->land2
ocean1(45:89,:) = f1->ocean1
ocean2(45:89,:) = f1->ocean2

f2 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/3_analysis1982_2023/flow/Distance_landocean_SH.nc", "r")
land1(0:44,:) = f2->land1
land2(0:44,:) = f2->land2
ocean1(0:44,:) = f2->ocean1
ocean2(0:44,:) = f2->ocean2


landocean = new((/90,180/), float)
landocean = where(ocean2.eq.1, -2, landocean)
landocean = where(ocean1.eq.1, -1, landocean)
landocean = where(land1.eq.1, 1, landocean)
landocean = where(land2.eq.1, 2, landocean)
copy_VarCoords(topo, landocean)






;;---------------------------------------------------  plots  ---------------------------------------------------
res = True
res@tmBorderThicknessF = 2.0
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
res@mpMinLonF = -180+10+1.5
res@mpMaxLonF = 180+10
res@mpCenterLonF = 10
res@tmXBMode = "Explicit"
res@tmXBValues = (/-225,-180,-135,-90,-45,0,45,90,135,180,225,270,315,360/)
res@tmXBLabels = (/"135E","180","135W","90W","45W","0","45E","90E","135E","180","135W","90W","45W","0"/)
res@tmYLMode = "Explicit"
res@tmYLValues = (/-90,-60,-30,0,30,60,90/)
res@tmYLLabels = (/"90S","60S","30S","0","30N","60N","90N"/)
res@tmXBMinorOn = True
res@tmXBMinorValues = fspan(-360, 360+15, 49+1)
res@tmXBLabelsOn = True
res@tmXBLabelFontHeightF = 0.020
res@tmYLLabelsOn = True
res@tmYLLabelFontHeightF = 0.020

res@tmYLMajorThicknessF = 1.5
res@tmXBMajorThicknessF = 1.5
res@tmYLMinorThicknessF = 1.0
res@tmXBMinorThicknessF = 1.0
res@tmYLMajorLengthF = 0.014
res@tmXBMajorLengthF = 0.014
res@tmYLMinorLengthF = 0.010
res@tmXBMinorLengthF = 0.010
res@tmXBLabelDeltaF = -0.35
res@tmYLLabelDeltaF = -0.35
res@tmXBMajorOutwardLengthF = 0.014
res@tmYLMajorOutwardLengthF = 0.014
res@tmXBMinorOutwardLengthF = 0.010
res@tmYLMinorOutwardLengthF = 0.010





res@mpFillOn                    = False
res@mpOutlineOn                 = True
res@mpGeophysicalLineThicknessF = 1.5
res@mpGeophysicalLineColor      = "gray45";"gray20"
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
; sres@cnLevels         = (/2700.0/)
; sres@cnFillColors     = (/"red","gray20"/)
sres@cnMissingValFillColor= -1
sres@cnLineThicknessF = 6.0
sres@cnLineDashPattern = 14


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
cnres@vpWidthF  = 0.42
cnres@vpHeightF = 0.42
cnres@vpXF      = 0.06
cnres@vpYF      = 0.5

cnres@trXMinF = -0.5
cnres@trXMaxF = 18.5
cnres@trYMinF = 0.0
cnres@trYMaxF = 90.0


cnres@tmXTOn        = False
cnres@tmYROn        = False
cnres@tmXBOn        = True
cnres@tmYLOn        = True
cnres@tmXBLabelFontHeightF = 0.014
cnres@tmYLLabelFontHeightF = 0.014
cnres@tmYRLabelFontHeightF = 0.014
cnres@gsnLeftStringFontHeightF   = 0.022
cnres@gsnRightStringFontHeightF  = 0.022
cnres@gsnCenterStringFontHeightF  = 0.022
cnres@tmYRLabelsOn   = False
; cnres@tmYROn         = True

cnres@tmYLMode := "Explicit"
cnres@tmYLValues := fspan(0.0, 80.0, 5)
cnres@tmYLLabels := (/"0","20N","40N","60N","80N"/)
cnres@tmYLMinorOn := True
cnres@tmYLMinorValues := fspan(0.0, 80.0, 9)

cnres@tmXBMode := "Explicit"
cnres@tmXBValues := fspan(0, 18.0, 19)
cnres@tmXBLabels := (/"M","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18"/)
cnres@tmXBMinorOn := False

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
cnres@tmYLMajorLengthF = 0.012
cnres@tmXBMajorLengthF = 0.012
cnres@tmYLMinorLengthF = 0.009
cnres@tmXBMinorLengthF = 0.009
cnres@tmXBLabelDeltaF = -0.4
cnres@tmYLLabelDeltaF = -0.4
cnres@tmXBMajorOutwardLengthF = 0.012
cnres@tmYLMajorOutwardLengthF = 0.012
cnres@tmXBMinorOutwardLengthF = 0.009
cnres@tmYLMinorOutwardLengthF = 0.009







lires = True
lires@tmBorderThicknessF = 1.1
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

lires@tmXBLabelFontHeightF = 0.072
lires@tmYLLabelFontHeightF = 0.072
lires@tmYRLabelFontHeightF = 0.072
lires@gsnLeftStringFontHeightF   = 0.025
lires@gsnRightStringFontHeightF  = 0.025
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
lires@vpWidthF  = 0.99
lires@vpHeightF = 0.75
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

lires@tmYLMajorThicknessF = 1.5
lires@tmYRMajorThicknessF = 1.5
lires@tmXBMajorThicknessF = 1.5
lires@tmYLMinorThicknessF = 1.0
lires@tmYRMinorThicknessF = 1.0
lires@tmXBMinorThicknessF = 1.0
lires@tmYLMajorLengthF = 0.06
lires@tmYRMajorLengthF = 0.06
lires@tmXBMajorLengthF = 0.06
lires@tmYLMinorLengthF = 0.04
lires@tmYRMinorLengthF = 0.04
lires@tmXBMinorLengthF = 0.04
lires@tmXBLabelDeltaF = -0.25
lires@tmYLLabelDeltaF = -0.25
lires@tmXBMajorOutwardLengthF = 0.06
lires@tmYLMajorOutwardLengthF = 0.06
lires@tmYRMajorOutwardLengthF = 0.06
lires@tmXBMinorOutwardLengthF = 0.04
lires@tmYLMinorOutwardLengthF = 0.04
lires@tmYRMinorOutwardLengthF = 0.04






;;*************************打点***************************;;
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
pltName = "/public/home/fcai/extreme5_degree/25_humid_heat/3_analysis1982_2023/Figure2a_landocean"
wks = gsn_open_wks(pltType, pltName)


colors = (/"#2b2a76","#73aa43","#f3c80d","#a21d2f"/)  ;; 深色
colors = (/"#6A699F","#9CC37B","#f5d33d","#BD606D"/)  ;; 浅色,0.5左右
; colors = (/"#6A699F","#5d8d35","#f3c80d","#8a1a26"/)  ;; 中色,0.8左右

res@cnMissingValFillColor = -1
res@cnFillMode = "CellFill"
res@gsnLeftString = ""
res@gsnLeftStringFontHeightF = 0.022
res@gsnCenterStringFontHeightF = 0.0235
res@tmXBLabelFontHeightF = 0.0175
res@tmYLLabelFontHeightF = 0.0175

; res@mpMinLonF = -180-75
; res@mpMaxLonF = 180-75
; res@mpCenterLonF = -75
res@tmXBLabelsOn = True
res@cnFillPalette := colors
res@cnLevels      = (/-1.5,0.0,1.5/)
res@cnFillColors  = (/0,1,2,3/)
plot(0) = gsn_csm_contour_map_ce(wks, landocean, res)





txres               = True
txres@txPerimOn = False
txres@txFontHeightF = 0.018
; txres@txFont = "helvetica"
txres@txFontColor = "#BD606D"
dtext1 = gsn_add_text(wks, plot(0), "inland", -130, -6, txres)
txres@txFontColor = "#f5d33d"
dtext2 = gsn_add_text(wks, plot(0), "coastland", -130, -19, txres)
txres@txFontColor = "#9CC37B"
dtext3 = gsn_add_text(wks, plot(0), "near-shore", -130, -32, txres)
txres@txFontColor = "#6A699F"
dtext4 = gsn_add_text(wks, plot(0), "off-shore", -130, -45, txres)







  pres = True
  pres@gsnFrame = False
  pres@gsnPanelLabelBar = False

  pres@gsnPanelTop  = 0.8
  pres@gsnPanelBottom = 0.2
  pres@gsnPanelLeft  = 0.05
  pres@gsnPanelRight = 0.95
  gsn_panel(wks,plot(0),(/1,1/),pres)





  



txres               = True
txres@txFont = "helvetica"
txres@txFontColor = "black"
txres@txFontHeightF = 0.0195
gsn_text_ndc(wks, "Four   types   of   grid   cells",  0.51, 0.72, txres)


txres@txFontHeightF = 0.027
txres@txFont = "helvetica-bold"
gsn_text_ndc(wks, "a",  0.10, 0.75, txres)
; gsn_text_ndc(wks, "b",  0.20, 0.645, txres)
; gsn_text_ndc(wks, "c",  0.20, 0.34, txres)



frame(wks)
   delete(wks)  
   system("convert -geometry 1600x1600 -density 1200x1200 -trim " + pltName + "."+pltType + " " + pltName + ".png")     
  system("rm " + pltName + "." + pltType + " -f")



end


