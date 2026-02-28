import streamlit as st
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from App_Functions import (
    solve_betaC, calculate_point2, calculate_point3, calculate_pure_moment,
    generate_side_view, cross_section, Icr_function, Moment_Calculation,
    draw_blocks_plotly
)

# 1. Page Config (Force wide layout)
st.set_page_config(page_title="Masonry Walls Design", layout="wide")

# CSS to force specific styling for the report look
st.markdown("""
<style>
    /* Force tables to look like engineering reports */
    table {
        width: 100%;
        border-collapse: collapse;
    }
    th {
        background-color: #f2f2f2;
        color: #333;
        font-weight: bold;
        border: 1px solid #ddd;
        padding: 8px;
        text-align: center;
    }
    td {
        border: 1px solid #ddd;
        padding: 8px;
        text-align: center;
    }
    /* Status tags */
    .pass-box {
        background-color: #d4edda;
        color: #155724;
        padding: 4px 8px;
        border-radius: 4px;
        font-weight: bold;
        display: inline-block;
    }
    .fail-box {
        background-color: #f8d7da;
        color: #721c24;
        padding: 4px 8px;
        border-radius: 4px;
        font-weight: bold;
        display: inline-block;
    }
</style>
""", unsafe_allow_html=True)

st.title("Masonry Walls Design Tool (Out-of-Plane)")

# --- SIDEBAR INPUTS ---
with st.sidebar:
    st.header("1. Wall Properties")
    H_val = st.number_input("Wall Height (H) [m]", value=8.0, step=0.1)
    
    t_map = {0.140: "140 mm", 0.190: "190 mm", 0.240: "240 mm", 0.290: "290 mm"}
    t_val = st.selectbox("Wall Thickness (t)", options=list(t_map.keys()), format_func=lambda x: t_map[x], index=2)
    
    fblock = st.selectbox("Block Strength (f'm) [MPa]", options=[10, 15, 20, 25, 30], index=3)
    
    S_options = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6]
    S_val = st.selectbox("Bar Spacing (S) [m]", options=S_options, index=5, format_func=lambda x: f"{int(x*1000)} mm")
    
    bar_options = [0.010, 0.015, 0.020, 0.025]
    bar_val = st.selectbox("Bar Diameter (db) [m]", options=bar_options, index=3, format_func=lambda x: f"{int(x*1000)} mm")

    st.header("2. Loads")
    P_DL = st.number_input("Dead Load (PDL) [kN/m]", value=10.0, step=1.0)
    P_LL = st.number_input("Live Load (PLL) [kN/m]", value=0.0, step=1.0)
    P_S = st.number_input("Snow Load (Ps) [kN/m]", value=0.0, step=1.0)
    e_val = st.number_input("Eccentricity (e) [mm]", value=0.0, step=5.0)
    W_val = st.number_input("Wind Load (W) [kPa]", value=1.0, step=0.1)

# --- CALCULATIONS ---

# Units & Constants
mm = 0.001
m = 1
faim = 0.6
fais = 0.85
emu = 0.003
d = t_val / 2

# 1. Cross Section
(t_res, beff_m_1, beff_m_2, As, Aseff_m, bg_m, bug_m_1, bug_m_2, A_gr, A_ug_1, 
 A_ug_2, Ae_1, Ae_2, fm_e_1, fm_e_2, I_gross_gr, I_gross_ug_1, I_gross_eff, 
 I_cr_eff_init, kd, n, E_m, ek, rho_SW, rho_g, rho_ug, fm_g, fm_ug, tf) = cross_section(t_val, S_val, bar_val, fblock)

# 2. Interaction Diagram Points
PMax = solve_betaC(0.6, fm_e_1, Ae_1, fm_g, bg_m, fm_ug, tf, bug_m_1, t_res)
betaC1 = float(PMax[0])

point2_results = calculate_point2(betaC1, faim, fm_g, bg_m, fm_ug, tf, bug_m_1, t_res, d, num_points=15)
point3_results, Mr_y, Pr_y, ey = calculate_point3(faim, fais, emu, fm_g, bg_m, fm_ug, tf, bug_m_2, t_res, d, Aseff_m, num_points=40)
pure_moment = calculate_pure_moment(faim, Aseff_m, d, fm_g, bg_m, fm_ug, tf, bug_m_2, t_res)

# Arrays for Envelope
M_env = [0] + [PMax[-1]] + [pt[4] for pt in point2_results] + [pt[4] for pt in point3_results] + [pure_moment[4]]
P_env = [PMax[-2]] + [PMax[-2]] + [pt[3] for pt in point2_results] + [pt[3] for pt in point3_results] + [pure_moment[3]]

# 3. Factored Loads Definitions
P_SW_mid = rho_SW * H_val / 2
# Load combinations descriptions
load_combos_desc = [
    "1.4 D",
    "1.25 D + 1.5 L",
    "0.9 D + 1.5 L",
    "1.25 D + Wind",
    "0.9 D + Wind",
    "1.25 D + 1.5 S",
    "0.9 D + 1.5 S"
]

P_F_list_calcs = [
    1.4 * (P_DL + P_SW_mid),
    1.25 * (P_DL + P_SW_mid) + 1.5 * P_LL,
    0.9 * (P_DL + P_SW_mid) + 1.5 * P_LL,
    1.25 * (P_DL + P_SW_mid),
    0.9 * (P_DL + P_SW_mid),
    1.25 * (P_DL + P_SW_mid) + 1.5 * P_S,
    0.9 * (P_DL + P_SW_mid) + 1.5 * P_S
]

# 4. Icr Iteration
Icr_results = []
for pf in P_F_list_calcs:
    I_cr_val, kd_val = Icr_function(Aseff_m, d, beff_m_2, bg_m, E_m, pf, tf)
    Icr_results.append(I_cr_val)

# 5. Moment Calculations (Extracting more details this time)
(ev, Mt_F_list, M_F_list, P_F_list, Pcr_list, Mag_list, Cm_list, EI_eff_raw, 
 EI_eff_list, betad, Min_EIeff, Max_EIeff, I_cr_array, betad_list) = Moment_Calculation(
    t_res, e_val, H_val, rho_SW, W_val, P_DL, P_LL, P_S, I_gross_eff, E_m, ek, Icr_results
)

# 6. Interpolation for Mr
sorted_indices = np.argsort(P_env)
P_env_sorted = np.array(P_env)[sorted_indices]
M_env_sorted = np.array(M_env)[sorted_indices]
M_interpolated = np.interp(P_F_list, P_env_sorted, M_env_sorted)

# --- VISUALS ---

c1, c2 = st.columns([1, 1])
with c1:
    st.subheader("Wall Cross-Section")
    fig_wall = draw_blocks_plotly(t_mm_actual=t_val*1000, s=S_val*1000, bar_diameter_mm=bar_val*1000)
    st.plotly_chart(fig_wall, use_container_width=True)
with c2:
    st.subheader("Deflection / Side View")
    fig_side = generate_side_view(H_val, P_DL, P_LL, P_S, e_val, W_val)
    st.plotly_chart(fig_side, use_container_width=True)

st.subheader("Interaction Diagram Checks")
fig = go.Figure()
fig.add_trace(go.Scatter(x=M_env, y=P_env, mode="lines", name="Resistance Envelope", line=dict(color='black', width=3)))
fig.add_trace(go.Scatter(x=Mt_F_list, y=P_F_list, mode="markers", name="Applied Loads (Mt)", marker=dict(color='red', size=10, symbol='circle')))
fig.update_layout(xaxis_title="Moment (kN⋅m)", yaxis_title="Axial Force (kN)", height=500, margin=dict(l=40, r=40, t=20, b=20))
st.plotly_chart(fig, use_container_width=True)

# --- DETAILED REPORT SECTION ---
st.markdown("---")
st.subheader("Detailed Design Report")

data_rows = []
for i in range(len(P_F_list)):
    # Calculate deflection (delta) roughly if not explicitly in output lists
    # Or use Mag_list to show Magnification
    
    mt = Mt_F_list[i]
    mr = M_interpolated[i]
    pf = P_F_list[i]
    
    # Status Logic
    if mt < 0:
         status_html = "<span class='fail-box'>Buckling</span>"
    elif mr >= mt:
         status_html = "<span class='pass-box'>Pass</span>"
    else:
         status_html = "<span class='fail-box'>Fail</span>"

    data_rows.append({
        "Load Case": f"<b>{i+1}</b><br><small>{load_combos_desc[i]}</small>",
        "Pf (Axial)<br>[kN]": f"{pf:.2f}",
        "Mf (Primary)<br>[kNm]": f"{M_F_list[i]:.2f}",
        "Mt (Total)<br>[kNm]": f"<b>{mt:.2f}</b>",
        "Mr (Resist)<br>[kNm]": f"<b>{mr:.2f}</b>",
        "Result": status_html
    })

# Convert to DF for display
df_results = pd.DataFrame(data_rows)

# Render as a HTML Table (Static, full view)
# This mimics the "Report" feel better than the interactive dataframe
st.markdown(
    df_results.to_html(escape=False, index=False), 
    unsafe_allow_html=True
)

# Intermediate Data Expander
with st.expander("See Calculation Details (Stiffness & Deflection)"):
    st.write("Intermediate values used in the Moment Magnification method:")
    detail_rows = []
    for i in range(len(P_F_list)):
        detail_rows.append({
            "Case": f"{i+1}",
            "EI effective [kN m²]": f"{EI_eff_list[i]:.0f}",
            "Magnifier (δns)": f"{Mag_list[i]:.3f}",
            "Cm": f"{Cm_list[i]:.2f}",
            "Beta_d": f"{betad_list[i]:.2f}"
        })
    st.table(pd.DataFrame(detail_rows))
