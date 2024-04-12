import streamlit as st
from streamlit_utils import *

st.set_page_config(
  page_title='Target Selection Dashboard',
  page_icon=n23_icon,
  layout='wide'
  )

add_logo()

if not check_password():
    st.stop()  # Do not continue if check_password is not True.

st.markdown(
    """
    # Target Selection Dashboard

    This dashboard was built to help facilitate the selection of new
    targets at Neuron23. By displaying all facets of data in one place,
    including genetic scores, biological function, and chemical druggability,
    we can  more easily select new targets in a data-driven way.\n
    <br/>

    **👈  Click "Target Selection" on the sidebar** to get started!\n
    <br/>

    **👈  Check out the "About"** page for detailed explanations on what info is available.
    """,
    unsafe_allow_html=True,
)