import streamlit as st
import tempfile
import os
import matplotlib.pyplot as plt
from defect import calculate_E_formation

# Streamlit page config
st.set_page_config(layout="wide")

# Title
st.title("Defect Formation Energy Viewer")

# Horizontal layout: left - editor, right - plot
col1, col2 = st.columns([1, 1.5])

with col1:
    st.subheader("📄 Input File Editor")
    default_text = """&VBM        # (eV)
7.1666

&band_gap        # with respect to the VBM (eV)
2.0091

&Host_type
Zn S

&Vacancies        # For substitutional impurity, the substituted atom too
Zn S

&Impurities
Cu

&Chemical_potentials            # (eV)
Zn  -5446.8099906
S   -277.5947224
Cu  -4968.9433547

&Host_supercell_energy          # (Ry)
-13468.80835818

&Charge_state_range  #two integer numbers // make sure you have a right number of alignment term
-1 1

&Defective_supercell_energy # charge_state(q) and energy (Ry) 
-1	-13412.38748518
0	-13413.00734819
1	-13413.58004343

&Correction_terms #(eV), for each charge state (e.g. -2 to 2): short-range potential energy from sxdefectalign, and E_corr (MP or LZ)
0.110   0.21088
0.000	0.00000
0.450   0.21088"""

    user_input = st.text_area("Edit input file:", value=default_text, height=500)
    submitted = st.button("Plot Formation Energy Diagram")

with col2:
    st.subheader("📊 Formation Energy Plot")
    if submitted:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_input_path = os.path.join(tmpdir, "input_file.txt")
            with open(tmp_input_path, "w") as f:
                f.write(user_input)

            try:
                # Run the defect calculation
                defect_instance = calculate_E_formation(tmp_input_path)

                # Since calculate_E_formation already plots, we'll just reuse the last plot
                st.pyplot(defect_instance.fig)
            except Exception as e:
                st.error(f"Error: {e}")
