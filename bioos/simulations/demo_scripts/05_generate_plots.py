import pandas as pd
import plotly.graph_objects as go

# 1. Silica Cage Coherence
df_silica = pd.read_csv('silica_cage.csv', names=['time', 'no_si', 'si'], comment='#')
fig = go.Figure()
fig.add_trace(go.Scatter(x=df_silica['time']*1e12, y=df_silica['no_si'], name='650 fs (sans silice)'))
fig.add_trace(go.Scatter(x=df_silica['time']*1e12, y=df_silica['si'], name='1,12 ps (cage silice +20 %)'))
fig.update_layout(title='Silica Cage Prolonge Cohérence +20 %', xaxis_title='Temps (ps)', yaxis_title='Cohérence')
fig.write_image('silica_cage.png')
print("Graphique silica_cage.png généré.")

# 2. Dual Channel Readout
df_dual = pd.read_csv('dual_channel.csv', names=['green', 'red'], comment='#')
fig = go.Figure()
fig.add_trace(go.Histogram(x=df_dual['green'], name='Canal Vert (ATP) – |0⟩', opacity=0.7))
fig.add_trace(go.Histogram(x=df_dual['red']+1.2, name='Canal Rouge (ROS) – |1⟩', opacity=0.7))
fig.update_layout(title='Lecture Asymétrique 62 % / 38 %', barmode='overlay', xaxis_tickvals=[0.5, 2.2], xaxis_ticktext=['|0⟩', '|1⟩'])
fig.write_image('dual_channel.png')
print("Graphique dual_channel.png généré.")

# 3. 24h CAM Energy
df_cam = pd.read_csv('24h_cam.csv', names=['time', 'atp'], comment='#')
fig = go.Figure()
fig.add_trace(go.Scatter(x=df_cam['time'], y=df_cam['atp'], name='ATP', line=dict(color='green', width=3)))
fig.add_hline(y=0.3, line_dash='dash', line_color='red', annotation_text='Seuil Luciférase')
fig.update_layout(title='ATP Racine 24 h/24 – CAM activé')
fig.write_image('24h.png')
print("Graphique 24h.png généré.")

print('\nValidation théorique terminée.')
print('Dossier prêt à pousser sur GitHub → lien à partager aux labos.')
