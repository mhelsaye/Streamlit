import streamlit as st
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from PIL import Image
# Import your original functions
from App_Functions import (
    solve_betaC, calculate_point2, calculate_point3, calculate_pure_moment,
    generate_side_view, cross_section, Icr_function, Moment_Calculation,
    draw_blocks_plotly
)

# 1. Page Configuration
st.set_page_config(page_title="Masonry Walls Design Tool", layout="wide")

# CSS to match your original professional report look
st.markdown("""
<style>
    .report-header {
        color: #0d6efd;
        font-weight: bold;
        font-size: 24px;
        margin-top: 20px;
        margin-bottom: 10px;
    }
    .sub-header {
        font-weight: bold;
        font-size: 18px;
        margin-top: 15px;
    }
    /* Force tables to look like the Dash Bootstrap tables */
    table {
        width: 100%;
        text-align: center;
        border-collapse: collapse;
    }
    th {
        background-color: #f8f9fa;
        font-weight: bold;
        padding: 10px;
        border-bottom: 2px solid #dee2e6;
    }
    td {
        padding: 8px;
        border-bottom: 1px solid #dee2e6;
    }
</style>
""", unsafe_allow_html=True)

st.title("Masonry Walls Design Tool")
st.subheader("(Out-of-Plane)")

# --- SIDEBAR INPUTS (Exact same defaults as your app.py) ---
with st.sidebar:
    st.header("Inputs")
    
    # Wall Properties
    H_val = st.number_input("Wall Height [m]", value=8.0, step=1.0)
    
    # T options: 140, 190, 240, 290
    t_options = [0.140, 0.190, 0.240, 0.290]
    t_labels = ["140 mm", "190 mm", "240 mm", "290 mm"]
    t_val = st.selectbox("Wall Thickness", options=t_options, format_func=lambda x: t_labels[t_options.index(x)], index=2)
    
    # Fblock options
    fblock_options = [10, 15, 20, 25, 30]
    fblock = st.selectbox("Block Strength (f'm) [MPa]", options=fblock_options, index=3)
    
    # Spacing options
    S_options = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6]
    S_val = st.selectbox("Bar Spacing [mm]", options=S_options, index=5, format_func=lambda x: f"{int(x*1000)} mm")
    
    # Bar options
    bar_options = [0.010, 0.015, 0.020, 0.025]
    bar_val = st.selectbox("Rebar Diameter", options=bar_options, index=3, format_func=lambda x: f"{int(x*1000)} mm")

    st.header("Loads")
    P_DL = st.number_input("Dead Load [kN/m]", value=10.0, step=0.1)
    P_LL = st.number_input("Live Load [kN/m]", value=0.0, step=0.1)
    P_S = st.number_input("Snow Load [kN/m]", value=0.0, step=0.1)
    e_val = st.number_input("Eccentricity [mm]", value=0.0, step=1.0)
    W_val = st.number_input("Wind Load [kPa]", value=1.0, step=0.1)

# --- CALCULATIONS ---

# Units & Constants
mm = 0.001
m = 1
faim = 0.6
fais = 0.85
emu = 0.003
d = t_val / 2

# 1. Run Cross Section Logic
(t_res, beff_m_1, beff_m_2, As, Aseff_m, bg_m, bug_m_1, bug_m_2, A_gr, A_ug_1, 
 A_ug_2, Ae_1, Ae_2, fm_e_1, fm_e_2, I_gross_gr, I_gross_ug_1, I_gross_eff, 
 I_cr_eff_init, kd, n, E_m, ek, rho_SW, rho_g, rho_ug, fm_g, fm_ug, tf) = cross_section(t_val, S_val, bar_val, fblock)

# 2. Interaction Diagram Data
PMax = solve_betaC(0.6, fm_e_1, Ae_1, fm_g, bg_m, fm_ug, tf, bug_m_1, t_res)
betaC1 = float(PMax[0])

point2_results = calculate_point2(betaC1, faim, fm_g, bg_m, fm_ug, tf, bug_m_1, t_res, d, num_points=15)
point3_results, Mr_y, Pr_y, ey = calculate_point3(faim, fais, emu, fm_g, bg_m, fm_ug, tf, bug_m_2, t_res, d, Aseff_m, num_points=40)
pure_moment = calculate_pure_moment(faim, Aseff_m, d, fm_g, bg_m, fm_ug, tf, bug_m_2, t_res)

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

# 4. Icr Iteration
Icr_results = []
kd_vals = []
for pf in P_F_list_calcs:
    I_cr_eff, kd_val = Icr_function(Aseff_m, d, beff_m_2, bg_m, E_m, pf, tf)
    Icr_results.append(I_cr_eff)
    kd_vals.append(kd_val)

# 5. Moment Calculations (The heavy lifting)
(ev, Mt_F_list, M_F_list, P_F_list, Pcr_list, Mag_list, Cm_list, EI_eff_raw, 
 EI_eff_list, betad, Min_EIeff, Max_EIeff, I_cr_array, betad_list) = Moment_Calculation(
    t_res, e_val, H_val, rho_SW, W_val, P_DL, P_LL, P_S, I_gross_eff, E_m, ek, Icr_results
)

# 6. Interpolation for Resistance (Mr)
sorted_indices = np.argsort(P_env)
P_env_sorted = np.array(P_env)[sorted_indices]
M_env_sorted = np.array(M_env)[sorted_indices]
M_interpolated = np.interp(P_F_list, P_env_sorted, M_env_sorted)


# --- VISUALIZATION (Top Section) ---
col1, col2 = st.columns([1, 1])
with col1:
    fig_wall = draw_blocks_plotly(t_mm_actual=t_val*1000, s=S_val*1000, bar_diameter_mm=bar_val*1000)
    st.plotly_chart(fig_wall, use_container_width=True)
with col2:
    fig_side = generate_side_view(H_val, P_DL, P_LL, P_S, e_val, W_val)
    st.plotly_chart(fig_side, use_container_width=True)

st.header("Interaction Diagram")
fig_id = go.Figure()
fig_id.add_trace(go.Scatter(x=M_env, y=P_env, mode="lines", name="Resistance Envelope", line=dict(color='black', width=3)))
fig_id.add_trace(go.Scatter(x=Mt_F_list, y=P_F_list, mode="markers", name="Total Moment (Mt)", marker=dict(color='red', size=10, symbol='circle')))
fig_id.add_trace(go.Scatter(x=M_F_list, y=P_F_list, mode="markers", name="Primary Moment (Mf)", marker=dict(color='blue', size=8, symbol='x')))
fig_id.update_layout(xaxis_title="Moment (kN⋅m)", yaxis_title="Axial Force (kN)", height=600)
st.plotly_chart(fig_id, use_container_width=True)


# --- DETAILED REPORT GENERATION (Restoring your exact text) ---

# Variable formatting for report
h_over_t = H_val / t_val * 1000 # Your logic uses 1000 scale
t_check_text = f"**Thickness of the block ($t$):** {t_val*1000:.0f} mm  ✅" if (t_val*1000) >= 140 else f"**Thickness of the block ($t$):** {t_val*1000:.0f} mm  ❌"

st.markdown(f"""
## 🧱 **Masonry Wall Summary Report**
This report summarizes the design of a masonry wall based on the provided inputs. The design is based on the **CSA S304-24** standard.
The wall is designed to resist out-of-plane loads, including external dead loads, self weight, live loads, and wind loads. The design also considers the slenderness effects and the interaction between axial loads and moments. The design procedure follows Masonry Structures Behavior and Design 2nd Edition by Bennett & Drysdale.
The inputs loads are unfactored. The Load combinations are obtained from **NBC 2020** as follows:
- **Load Combination 1**: $1.4\\, P_\\text{{DL}}$ 
- **Load Combination 2**: $1.25\\, P_\\text{{DL}} + 1.5\\, P_\\text{{LL}}$
- **Load Combination 3**: $0.9\\, P_\\text{{DL}} +  1.5\\, P_\\text{{LL}}$
- **Load Combination 4**: $1.25\\, P_\\text{{DL}}  + 1.4\\, W$ 
- **Load Combination 5**: $0.9\\, P_\\text{{DL}}  + 1.4\\, W$
- **Load Combination 6**: $1.25\\, P_\\text{{DL}} + 1.5\\, P_\\text{{S}}$
- **Load Combination 7**: $0.9\\, P_\\text{{DL}} + 1.5\\, P_\\text{{S}}$

### **Wall Properties**
- **Height**: ${H_val:.0f}\\,\\text{{m}}$
- **Thickness**: ${t_val*1000:.0f}\\,\\text{{mm}}$
- **Block Strength**: ${fblock:.0f}\\,\\text{{MPa}}$
- **Bar Spacing**: ${S_val*1000:.0f}\\,\\text{{mm}}$
- **Bar Diameter**: ${bar_val*1000:.0f}\\,\\text{{M}}$

### **Applied Loads**
- **Dead Load ($P_\\text{{DL}}$)**: ${P_DL:.2f}\\,\\text{{kN}}$
- **Live Load ($P_\\text{{LL}}$)**: ${P_LL:.2f}\\,\\text{{kN}}$
- **Snow Load ($P_\\text{{S}}$)**: ${P_S:.2f}\\,\\text{{kN}}$
- **Eccentricity ($e$)**: ${e_val:.1f}\\,\\text{{mm}}$
- **Wind Load ($W$)**: ${W_val:.2f}\\,\\text{{kPa}}$

---
### **Slenderness Check**
- $\\dfrac{{h}}{{t}} = {h_over_t:.2f}$

**According to CSA S304-24:** $h/t > 30 \\Rightarrow \\text{{Slenderness effects must be considered}}$

In accordance with relevant provisions:
- The minimum thickness of the block is 140 mm (C11.7.4.6.2).
- Boundary conditions are assumed to be pinned unless a more detailed analysis is provided (C11.7.4.6.3).
- {t_check_text}
""")

st.markdown(f"""
### **Effective Section Properties**
The width of the effective compression zone ($b$) is determined based on CSA S304-24 (C 11.6.1):
$$ b = \\min(6t, s) = {beff_m_2*1000:.0f}\\,\\text{{mm}} $$

The effective compression zone consists of grouted and ungrouted portions:
- **Grouted Width ($b_{{gr}}$):** {bg_m*1000:.0f} mm
- **Ungrouted Width ($b_{{ug}}$):** {bug_m_2*1000:.0f} mm
- **Effective Area ($A_e$):** {Ae_2*1e6:.0f} mm²
""")

st.markdown(f"""
### **Self-Weight Calculation**
For slender walls ($h/t > 30$), self-weight must be included.
- **Density Used:** {rho_g:.2f} kN/m³ (Grouted) / {rho_ug:.2f} kN/m³ (Ungrouted)
- **Calculated Self-Weight at Mid-height:** {P_SW_mid:.2f} kN/m
""")

# --- TABLE 1: Cracked Inertia ---
st.markdown("### **Calculation of Cracked Moment of Inertia ($I_{cr}$)**")
st.markdown("CSA S304-24 (C11.7.4.4) defines $I_{cr}$ based on the factored axial load.")

icr_data = []
for i in range(len(Icr_results)):
    icr_data.append({
        "Load Case": i+1,
        "kd [mm]": f"{kd_vals[i]*1000:.2f}",
        "Icr [mm⁴]": f"{Icr_results[i]:.2e}",
        "Pf [kN]": f"{P_F_list_calcs[i]:.2f}"
    })
st.table(pd.DataFrame(icr_data))

# --- TABLE 2: Moment Magnification ---
st.markdown("### **Determining the Total Applied Moment Magnification (C11.7.4.3)**")
st.markdown(f"""
To determine the total applied moment, the moment magnification factor ($$\\Psi$$) is calculated using the effective stiffness ($$EI_{{eff}}$$).
$$ M_{{tot}} = \\Psi M_f $$
""")

mag_data = []
for i in range(len(P_F_list)):
    mag_data.append({
        "Case": i+1,
        "Mₜ [kNm]": f"{Mt_F_list[i]:.2f}",
        "M₁ [kNm]": f"{M_F_list[i]:.2f}",
        "Pf [kN]": f"{P_F_list[i]:.2f}",
        "Pcr [kN]": f"{Pcr_list[i]:.0f}",
        "Mag. Factor (Ψ)": f"{Mag_list[i]:.3f}",
        "EIeff [kNm²]": f"{EI_eff_list[i]:.0f}",
    })
st.table(pd.DataFrame(mag_data))


# --- TABLE 3: Final Design Verification ---
st.markdown("### **Final Design Verification ($M_r$ vs $M_{tot}$)**")

final_data = []
for i in range(len(P_F_list)):
    mt = Mt_F_list[i]
    mr = M_interpolated[i]
    pf = P_F_list[i]
    
    # Your exact logic for Check
    if mt < 0:
        status = "❌ Buckling Failure"
        is_safe = False
    elif mr >= mt:
        status = "✅ Pass"
        is_safe = True
    else:
        status = "❌ Fail"
        is_safe = False
        
    final_data.append({
        "Load Case": i+1,
        "Pf [kN]": f"{pf:.2f}",
        "Mt (Total) [kNm]": f"{mt:.2f}",
        "Mr (Resisting) [kNm]": f"{mr:.2f}",
        "Result": status
    })

df_final = pd.DataFrame(final_data)

# Color styling function
def style_results(val):
    color = 'red' if 'Fail' in val or 'Buckling' in val else 'green'
    return f'color: {color}; font-weight: bold'

st.table(df_final.style.applymap(style_results, subset=['Result']))
