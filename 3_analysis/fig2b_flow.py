

import plotly.graph_objects as go
from plotly.subplots import make_subplots
import xarray as xr
import numpy as np



labels = [" off-shore", " near-shore", " coastland", " inland", " off-shore", " near-shore", " coastland", " inland"]

source = [0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3]
target = [4, 5, 6, 7, 4, 5, 6, 7, 4, 5, 6, 7, 4, 5, 6, 7]

source_circle = [0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3]
target_circle = [1, 2, 3, 0, 2, 3, 0, 1, 3, 0, 1, 2]
# link_color1 = ["rgba(128,128,128,0.3)","rgba(128,128,128,0.3)","rgba(162,29,47,0.3)",
#               "rgba(128,128,128,0.3)","rgba(128,128,128,0.3)","rgba(128,128,128,0.3)",
#               "rgba(115,170,67,0.4)","rgba(128,128,128,0.3)","rgba(128,128,128,0.3)"]
# link_color2 = ["rgba(128,128,128,0.3)","rgba(162,29,47,0.3)","rgba(128,128,128,0.3)",
#               "rgba(43,42,118,0.4)","rgba(128,128,128,0.3)","rgba(128,128,128,0.3)",
#               "rgba(128,128,128,0.3)","rgba(128,128,128,0.3)","rgba(128,128,128,0.3)"]
# link_color3 = ["rgba(128,128,128,0.3)","rgba(162,29,47,0.3)","rgba(128,128,128,0.3)",
#               "rgba(43,42,118,0.4)","rgba(128,128,128,0.3)","rgba(43,42,118,0.4)",
#               "rgba(128,128,128,0.3)","rgba(115,170,67,0.4)","rgba(128,128,128,0.3)"]

ds1 = xr.open_dataset("/public/home/fcai/extreme5_degree/25_humid_heat/3_analysis1982_2023/fig2_flow/Flow_landocean_qs_qs_lag0_5_NH.nc")

flow0 = ds1.flow_weight[0,:,:].data
flow_0day_NH = np.zeros(16)
for i in range(4): flow0[i,i] = 0
flow_0day_NH[0:4] = flow0[0,:]
flow_0day_NH[4:8] = flow0[1,:]
flow_0day_NH[8:12] = flow0[2,:]
flow_0day_NH[12:16] = flow0[3,:]

flow5 = ds1.flow_weight[1,:,:].data
flow_5day_NH = np.zeros(16)
for i in range(4): flow5[i,i] = 0
flow_5day_NH[0:4] = (flow5[0,:] - flow5[:,0]) * (-1.0)
flow_5day_NH[4:8] = (flow5[1,:] - flow5[:,1]) * (-1.0)
flow_5day_NH[8:12] = (flow5[2,:] - flow5[:,2]) * (-1.0)
flow_5day_NH[12:16] = (flow5[3,:] - flow5[:,3]) * (-1.0)
flow_5day_NH[np.where(flow_5day_NH<0)] = 0




ds1 = xr.open_dataset("/public/home/fcai/extreme5_degree/25_humid_heat/3_analysis1982_2023/fig2_flow/Flow_landocean_qs_qs_lag0_5_SH.nc")

flow0 = ds1.flow_weight[0,:,:].data
flow_0day_SH = np.zeros(16)
for i in range(4): flow0[i,i] = 0
flow_0day_SH[0:4] = flow0[0,:]
flow_0day_SH[4:8] = flow0[1,:]
flow_0day_SH[8:12] = flow0[2,:]
flow_0day_SH[12:16] = flow0[3,:]

flow5 = ds1.flow_weight[1,:,:].data
flow_5day_SH = np.zeros(16)
for i in range(4): flow5[i,i] = 0
flow_5day_SH[0:4] = (flow5[0,:] - flow5[:,0]) * (-1.0)
flow_5day_SH[4:8] = (flow5[1,:] - flow5[:,1]) * (-1.0)
flow_5day_SH[8:12] = (flow5[2,:] - flow5[:,2]) * (-1.0)
flow_5day_SH[12:16] = (flow5[3,:] - flow5[:,3]) * (-1.0)
flow_5day_SH[np.where(flow_5day_SH<0)] = 0
print(flow5[0,2])
print(flow5[2,0])






##--------------------------------------------- Plot  --------------------------------------------##
fig = make_subplots(rows=1, cols=2, specs=[[{'type':'sankey'},{'type':'sankey'}]])


colors = ["#2b2a76","#73aa43","#f3c80d","#a21d2f","#2b2a76","#73aa43","#f3c80d","#a21d2f"]  ## 深色
colors = ["#6A699F","#9CC37B","#f5d33d","#BD606D","#6A699F","#9CC37B","#f5d33d","#BD606D"]  ## 浅色，0.5左右
# colors = ["#6A699F","#5d8d35","#f3c80d","#7f1a26","#6A699F","#5d8d35","#f3c80d","#7f1a26"]  ## 中色，0.8左右

Sankey2 = go.Sankey(
node=dict(pad=15, thickness=20, color=colors,
          line=dict(color="black", width=0.5),label=labels),
link=dict(source=source, target=target, value=flow_5day_NH, color='rgba(128,128,128,0.4)'))


Sankey4 = go.Sankey(
node=dict(pad=15, thickness=20, color=colors,
          line=dict(color="black", width=0.5),label=labels),
link=dict(source=source, target=target, value=flow_5day_SH, color='rgba(128,128,128,0.4)'))


fig.add_trace(Sankey2, row=1, col=1)
fig.add_trace(Sankey4, row=1, col=2)





fig.update_layout(
    font=dict(color='black', size=13),
    width = 500,
    height = 250,
    
    title=dict(
    text='',#'Flow    among    three    oceans',
    x=0.5,
    y=0.995,
    xanchor='center',
    font=dict(color='black', size=20)),
    margin=dict(l=10,r=10,t=50,b=28)
)


fig.add_annotation(text='b', x=0.007, y=1.32, xref='paper',yref='paper',showarrow=False, font=dict(color='black', size=20.5, family='Arial', weight='bold'))
fig.add_annotation(text='c', x=0.57, y=1.32, xref='paper',yref='paper',showarrow=False, font=dict(color='black', size=20.5, family='Arial', weight='bold'))


fig.add_annotation(text='moist  flow  (NH)', x=0.103, y=1.22, xref='paper',yref='paper',showarrow=False, font=dict(color='black', size=16, family='Arial'))
fig.add_annotation(text='moist  flow  (SH)', x=0.907, y=1.22, xref='paper',yref='paper',showarrow=False, font=dict(color='black', size=16, family='Arial'))
# fig.add_annotation(text='lag = 5 day', x=0.501, y=-0.16, xref='paper',yref='paper',showarrow=False, font=dict(color='black', size=15.2, family='Arial'))
fig.add_annotation(text='day0', x=-0.008, y=-0.16, xref='paper',yref='paper',showarrow=False, font=dict(color='black', size=12.0, family='Arial'))
fig.add_annotation(text='day5', x=0.43, y=-0.16, xref='paper',yref='paper',showarrow=False, font=dict(color='black', size=12.0, family='Arial'))
fig.add_annotation(text='day0', x=0.573, y=-0.16, xref='paper',yref='paper',showarrow=False, font=dict(color='black', size=12.0, family='Arial'))
fig.add_annotation(text='day5', x=1.01, y=-0.16, xref='paper',yref='paper',showarrow=False, font=dict(color='black', size=12.0, family='Arial'))




fig.write_image('/public/home/fcai/extreme5_degree/25_humid_heat/3_analysis1982_2023/Figure2b_flow.png', scale=4)

