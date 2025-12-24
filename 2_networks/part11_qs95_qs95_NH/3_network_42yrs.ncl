load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
begin

dim_0 = (/45,180, 45,180/)


network_44yrs = new((/dim_0(0),dim_0(1),dim_0(2),dim_0(3)/), float)
network_44yrs = 0.0


do iyear=0,41
print(" iyear = " + (iyear+1982))
;;------------  read concurrence matrix  ------------;;
f1 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/2_network/part11_qs95_qs95_NH/qs95_qs95_NH2x2_"+int2flt(iyear+1982)+".nc","r")
networks0 = (f1->networks0)

networks0@_FillValue = 1.0e20
value = 1.0e20
if (any(isnan_ieee(networks0))) then
      replace_ieeenan (networks0, value, 0)
      networks0@_FillValue = value
end if


network_44yrs = network_44yrs + networks0
delete(networks0)
printMinMax(network_44yrs, 0)
end do

printVarSummary(network_44yrs)

network_44yrs = where(network_44yrs.ge.17.0, network_44yrs, 0.0)
network_44yrs@_FillValue = 1.0e20



;;-------------------- to nc file ---------------------;;
output="/public/home/fcai/extreme5_degree/25_humid_heat/2_network/part11_qs95_qs95_NH/"
name1="qs95_qs95_NH2x2_1982_2023_sign99"
system("/bin/rm -f "+output+name1+".nc")
ncdf1=addfile(output+name1+".nc","c")
ncdf1->networks0=network_44yrs


end

