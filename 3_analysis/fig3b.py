# %%
import numpy as np
import xarray as xr




##--##--##--##--##--##--      --##--##--##--##--##--##
ds1 = xr.open_dataset("/public/home/fcai/extreme5_degree/25_humid_heat/2_network/part1_wbt95_sst95_NH/distance_PDF_1982_2023_wbt95_sst95_NH.nc")
pdf1_NH = (ds1.pdf.data)
ref1_NH = (ds1.ref.data)
pdf_ref1_NH = (ds1.pdf_ref.data)

ds2 = xr.open_dataset("/public/home/fcai/extreme5_degree/25_humid_heat/2_network/part2_wbt95_sst95_SH/distance_PDF_1982_2023_wbt95_sst95_SH.nc")
pdf1_SH = (ds2.pdf.data)
ref1_SH = (ds2.ref.data)
pdf_ref1_SH = (ds2.pdf_ref.data)

pdf_ref1 = (pdf1_NH + pdf1_SH) / (ref1_NH + ref1_SH)
# print(pdf_ref1)
# print(pdf_ref1_NH)



ds3 = xr.open_dataset("/public/home/fcai/extreme5_degree/25_humid_heat/2_network/part3_wbt95_qs95_NH/distance_PDF_1982_2023_wbt95_qs95_NH.nc")
pdf2_NH = (ds3.pdf.data)
ref2_NH = (ds3.ref.data)
pdf_ref2_NH = (ds3.pdf_ref.data)

ds4 = xr.open_dataset("/public/home/fcai/extreme5_degree/25_humid_heat/2_network/part4_wbt95_qs95_SH/distance_PDF_1982_2023_wbt95_qs95_SH.nc")
pdf2_SH = (ds4.pdf.data)
ref2_SH = (ds4.ref.data)
pdf_ref2_SH = (ds4.pdf_ref.data)

pdf_ref2 = (pdf2_NH + pdf2_SH) / (ref2_NH + ref2_SH)



ds5 = xr.open_dataset("/public/home/fcai/extreme5_degree/25_humid_heat/2_network/part5_sst95_qs95_NH/distance_PDF_1982_2023_sst95_qs95_NH.nc")
pdf3_NH = (ds5.pdf.data)
ref3_NH = (ds5.ref.data)
pdf_ref3_NH = (ds5.pdf_ref.data)

ds6 = xr.open_dataset("/public/home/fcai/extreme5_degree/25_humid_heat/2_network/part6_sst95_qs95_SH/distance_PDF_1982_2023_sst95_qs95_SH.nc")
pdf3_SH = (ds6.pdf.data)
ref3_SH = (ds6.ref.data)
pdf_ref3_SH = (ds6.pdf_ref.data)

pdf_ref3 = (pdf3_NH + pdf3_SH) / (ref3_NH + ref3_SH)






# %%
## Plotting Import
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.use('agg')


##--##--##--  create figure  --##--##--##
width_in_inches = 210 / 25.4  ## A4 page
height_in_inches = 297 / 25.4
fig = plt.figure(dpi=400,figsize=(width_in_inches,height_in_inches))





def create_map(ax):
  # plt.axes( xscale="log" )
  plt.xlim(0.0,80.0)
  plt.ylim(0.0,1.0)
  ax.tick_params(axis='both', length=6.0, width=0.7, labelsize=9.5)
  ax.spines['bottom'].set_linewidth(0.8)
  ax.spines['top'].set_linewidth(0.8)
  ax.spines['left'].set_linewidth(0.8)
  ax.spines['right'].set_linewidth(0.8)
  # ax.set_xticks(np.linspace(0.0,100.0,11)); ax.set_xticklabels(["0","1","2","3","4","5","6","7","8","9","10"], color='black')
  ax.set_xticks(np.linspace(0.0,80.0,9)); ax.set_xticklabels(["0","1","2","3","4","5","6","7","8"], color='black')
  # ax.set_yticks(np.linspace(0.0,1.0,11)); ax.set_yticklabels(["0","","20%","","40%","","60%","","80%","","100%"], color='black')
  ax.spines['left'].set_color('black')
  ax.tick_params(axis='y', right=False, length=6.0, width=0.7, labelsize=8.5, color='black')
  return ax

  


x = np.linspace(0,102,52)
ax3 = plt.axes([0.05,0.4, 0.23,0.15])
create_map(ax3)
ax3.tick_params(axis='both', length=5.0, width=0.8, labelsize=7.5)
ax3.set_xlabel('Distance  (x1000 km)', size=9.0)
ax3.set_ylabel('Proportion', size=9.0)
plt.ylim(0.0,0.8)
ax3.set_yticks(np.linspace(0.0,0.8,9)); ax3.set_yticklabels(["0","","20%","","40%","","60%","","80%"], color='black')
ax3.plot(x[1:41], pdf_ref1_NH[1:41], label='NH', color='#2b2a76', markerfacecolor='#2b2a76', markeredgecolor='#2b2a76', markeredgewidth=0.0, linewidth=0.8,  markersize=0.0, linestyle='-')
ax3.plot(x[1:41], pdf_ref1_SH[1:41], label='SH', color='#009E74', markerfacecolor='#009E74', markeredgecolor='#009E74', markeredgewidth=0.0, linewidth=0.8,  markersize=0.0, linestyle='-')
ax3.plot(x[1:41], pdf_ref1[1:41], label='Globe', color='#B83945', markerfacecolor='#DB9BA1', markeredgecolor='#B83945', marker='o', markeredgewidth=0.8, linewidth=0.7,  markersize=4.5, linestyle='-')
ax3.legend(loc=(0.48, 0.36), fontsize=8.5, frameon=False)
ax3.text(24.0,0.69, 'SST   &   Wbt', fontsize=8.5, color='black')
# ax3.set_title('SST   &   Wbt', fontsize=9.5, color='black')




ax5 = plt.axes([0.31,0.4, 0.23,0.15])
create_map(ax5)
ax5.tick_params(axis='both', length=5.0, width=1.0, labelsize=7.5)
ax5.set_xlabel('Distance  (x1000 km)', size=9.0)
ax5.set_ylabel('', size=9.0)
plt.ylim(0.0,0.8)
ax5.set_yticks(np.linspace(0.0,0.8,9)); ax5.set_yticklabels(["","","","","","","","",""], color='black')
ax5.plot(x[1:41], pdf_ref3_NH[1:41], label='NH', color='#2b2a76', markerfacecolor='#2b2a76', markeredgecolor='#2b2a76', markeredgewidth=0.0, linewidth=0.8,  markersize=0.0, linestyle='-')
ax5.plot(x[1:41], pdf_ref3_SH[1:41], label='SH', color='#009E74', markerfacecolor='#009E74', markeredgecolor='#009E74', markeredgewidth=0.0, linewidth=0.8,  markersize=0.0, linestyle='-')
ax5.plot(x[1:41], pdf_ref3[1:41], label='Globe', color='#B83945', markerfacecolor='#DB9BA1', markeredgecolor='#B83945', marker='o', markeredgewidth=0.8, linewidth=0.7,  markersize=4.5, linestyle='-')
ax5.legend(loc=(0.48, 0.36), fontsize=8.5, frameon=False)
ax5.text(25.0,0.69, 'SST   &   Wet', fontsize=8.5, color='black')
# ax5.set_title('SST   &   Wet', fontsize=9.5, color='black')






# ax3.text(-1.0, 0.81, s='b', fontsize=13.8, color='black', fontweight='bold')
# ax3.text(100.0, 0.81, s='c', fontsize=13.8, color='black', fontweight='bold')
# ax3.text(203.5, 0.81, s='d', fontsize=13.8, color='black', fontweight='bold')


fig.savefig('/public/home/fcai/extreme5_degree/25_humid_heat/3_analysis1982_2023/Figure3b.png', bbox_inches = 'tight')
# fig.savefig('/public/home/fcai/extreme1_AT/NC2_1980_2010/3_era5_analysis/Figure1.pdf')
plt.show()


