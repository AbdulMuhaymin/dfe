import streamlit as st
import tempfile
import os
from defect import calculate_E_formation

# Streamlit page config
st.set_page_config(layout="wide")

# Title
st.title("Defect Formation Energy Viewer")
isCalculated = False

# Horizontal layout: left - editor, right - plot
col1, col2, col3 = st.columns([1, 1.5, 1])

with col1:
    st.subheader("📄 Input File Editor")
    default_text = """&Metadata
Co-Vac in MgO, relaxed, Mg-rich

&VBM        # (eV)
8.3619

&band_gap      # with respect to the VBM (eV)
4.47

&Host_type
Mg O

&Vacancies
Mg O

&Impurities
Co

&Chemical_potentials            # (eV)
Mg  -1478.3036971           #Mg-rich
O   -438.112148             #O-poor
#Mg -1483.7237              #Mg-poor
#O -432.69214758            #O-rich
Co  -3685.37905

&Host_supercell_energy          # (Ry)
-4507.34148788

&Charge_state_range
-2 2

&Defective_supercell_energy     # (Ry)
-2      -4667.53085681
-1	-4668.53056447
0	-4669.37300222
1	-4670.09251159
2       -4670.74789727

&Correction_terms #(eV)
0.000     0.00000
0.000	  0.00000
0.000     0.00000
0.000	  0.00000
0.000     0.00000"""

    user_input = st.text_area("Edit input file:", value=default_text, height=500)
    submitted = st.button("Plot Formation Energy Diagram")

with col2:
    st.subheader("📊 Defect Formation Energy")
    if submitted:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_input_path = os.path.join(tmpdir, "input_file.txt")
            with open(tmp_input_path, "w") as f:
                f.write(user_input)

            try:
                # Run the defect calculation
                defect_instance = calculate_E_formation(tmp_input_path)
                analysis = defect_instance.get_analysis()
                isCalculated = True

                # Since calculate_E_formation already plots, we'll just reuse the last plot
                st.pyplot(defect_instance.fig)
            except Exception as e:
                st.error(f"Error: {e}")

with col3:
    st.subheader("📈 Analysis")
    if isCalculated:
        st.write(analysis[0])
        st.markdown(":blue-badge[Download Formation Energy Data:]")
        st.download_button(
        data=analysis[1],
        file_name= defect_instance.get_output_fname() + '.txt',
        label= defect_instance.get_output_fname() + '.txt',
        mime='text/plain'
    )
