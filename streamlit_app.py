import streamlit as st
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from App_Functions import (
    solve_betaC, calculate_point2, calculate_point3, calculate_pure_moment,
    generate_side_view, cross_section, Icr_function, Moment_Calculation,
    draw_blocks_plotly
)

# 1. Page Config (Wide mode for better tables)
st.set_page_config(page_title="Masonry Walls Design", layout="wide")

st.markdown("""
<style>
    .reportview-container {
        background: #f0f2f6;
    }
    .main-header {
        font-size: 30px; 
        font-weight: bold; 
        color: #333; 
        text-align: center;
        margin-bottom: 20px;
    }
    .pass-tag {
        background-color: #d4edda;
        color: #155724;
        padding: 4px 8px;
        border-radius: 4px;
        font-weight: bold;
    }
    .fail-tag {
        background-color: #f8d7da;
        color: #721c24;
        padding: 4px 8px;
        border-radius: 4px;
        font-weight: bold;
    }
</style>
""", unsafe_allow_html=True)

st.markdown('<div class="main-header">Masonry Walls Design Tool (Out-of-Plane)</div>', unsafe_allow_html=True)

# --- SIDEBAR INPUTS ---
with st.sidebar:
    st.header("1. Wall Properties")
    H_val = st.number_input("Wall Height (H) [m]", value=8.0, step=0.1)
    
    # Options mapped to match your original logic
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

# --- MAIN LOGIC ---

# Units & Constants
mm = 0.001
m = 1
faim = 0.6
fais = 0.85
emu = 0.003
d = t_val / 2
# Note: In your original code, 'd' was roughly t/2. 
# Adjust logic if 'd' is calculated differently in your real function.

# 1. Run Cross Section
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

# 3. Factored Loads
P_SW_mid = rho_SW * H_val / 2
P_F_list_calcs = [
    1.4 * (P_DL + P_SW_mid),
    1.25 * (P_DL + P_SW_mid) + 1.5 * P_LL,
    0.9 * (P_DL + P_SW_mid) + 1.5 * P_LL,
    1.25 * (P_DL + P_SW_mid),
    0.9 * (P_DL + P_SW_mid),
    1.25 * (P_DL + P_SW_mid) + 1.5 * P_S,
    0.9 * (P_DL + P_SW_mid) + 1.5 * P_S
]

# 4. Icr Iteration (Copied from your logic)
Icr_results = []
for pf in P_F_list_calcs:
    I_cr_val, kd_val = Icr_function(Aseff_m, d, beff_m_2, bg_m, E_m, pf, tf)
    Icr_results.append(I_cr_val)

# 5. Moment Calculations
(ev, Mt_F_list, M_F_list, P_F_list, Pcr_list, Mag_list, Cm_list, EI_eff_raw, 
 EI_eff_list, betad, Min_EIeff, Max_EIeff, I_cr_array, betad_list) = Moment_Calculation(
    t_res, e_val, H_val, rho_SW, W_val, P_DL, P_LL, P_S, I_gross_eff, E_m, ek, Icr_results
)

# 6. Interpolation Logic (CRITICAL for Output Table)
# We sort P_env to make interpolation work correctly
sorted_indices = np.argsort(P_env)
P_env_sorted = np.array(P_env)[sorted_indices]
M_env_sorted = np.array(M_env)[sorted_indices]
M_interpolated = np.interp(P_F_list, P_env_sorted, M_env_sorted)

# --- DISPLAY SECTION ---

# Top Visuals
c1, c2 = st.columns([1, 1])
with c1:
    st.subheader("Wall Cross-Section")
    fig_wall = draw_blocks_plotly(t_mm_actual=t_val*1000, s=S_val*1000, bar_diameter_mm=bar_val*1000)
    st.plotly_chart(fig_wall, use_container_width=True)
with c2:
    st.subheader("Deflection / Side View")
    fig_side = generate_side_view(H_val, P_DL, P_LL, P_S, e_val, W_val)
    st.plotly_chart(fig_side, use_container_width=True)

# Interaction Diagram
st.subheader("Interaction Diagram Checks")
fig = go.Figure()
# Envelope
fig.add_trace(go.Scatter(x=M_env, y=P_env, mode="lines", name="Envelope (Resistance)", line=dict(color='black', width=3)))
# Applied Moments (Total)
fig.add_trace(go.Scatter(x=Mt_F_list, y=P_F_list, mode="markers", name="Mt (Applied Total)", marker=dict(color='red', size=10, symbol='circle')))
# Applied Moments (Primary)
fig.add_trace(go.Scatter(x=M_F_list, y=P_F_list, mode="markers", name="Mf (Primary)", marker=dict(color='blue', size=8, symbol='x')))

fig.update_layout(
    xaxis_title="Moment (kN⋅m)",
    yaxis_title="Axial Force (kN)",
    height=600,
    margin=dict(l=40, r=40, t=40, b=40),
    legend=dict(x=0.7, y=0.9)
)
st.plotly_chart(fig, use_container_width=True)


# --- DETAILED OUTPUT TABLE ---
st.subheader("Design Verification Report")

# Create Data for the Table
data_rows = []
for i in range(len(P_F_list)):
    pf = P_F_list[i]
    mt = Mt_F_list[i] # Total Moment
    mr = M_interpolated[i] # Resistance
    
    # Logic from your Dash App
    if mt < 0: # Buckling failure if Moment is negative/undefined in your logic
         status = "❌ Buckling Fail"
         is_safe = False
    elif mr >= mt:
         status = "✅ Pass"
         is_safe = True
    else:
         status = "❌ Fail"
         is_safe = False

    data_rows.append({
        "Load Case": f"Case {i+1}",
        "Axial Load (Pf) [kN]": f"{pf:.2f}",
        "Total Moment (Mt) [kNm]": f"{mt:.2f}",
        "Resistance (Mr) [kNm]": f"{mr:.2f}",
        "Result": status
    })

# Convert to DataFrame
df_results = pd.DataFrame(data_rows)

# Function to color code the output similar to Dash
def color_status(val):
    if "Pass" in val:
        return 'color: green; font-weight: bold;'
    else:
        return 'color: red; font-weight: bold;'

# Apply styling
styled_df = df_results.style.applymap(color_status, subset=['Result'])

# Display
st.dataframe(styled_df, use_container_width=True, hide_index=True)

# Final Summary Box
all_passed = all("Pass" in row["Result"] for row in data_rows)
if all_passed:
    st.success("### ✅ DESIGN OK: All Load Cases Passed")
else:
    st.error("### ❌ DESIGN FAILED: Check red Load Cases above")
