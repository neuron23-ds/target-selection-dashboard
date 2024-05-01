import streamlit as st
from streamlit_utils import n23_icon, add_logo, check_password

st.set_page_config(
  page_title='Target Selection Dashboard',
  page_icon=n23_icon,
  # layout='wide'
  )

add_logo()

if not check_password():
    st.stop()  # Do not continue if check_password is not True.

st.markdown(
    """
    # Target Selection Dashboard

    Welcome to the Target Selection Dashboard!

    This was built to help facilitate the selection of new targets at Neuron23. 
    In it, you will find data on genetic scores, biological function, chemical druggability, and clinical information.
    This was built by the Data Science team in conjunction with both the Biology team and Chemistry team.\n

    See the "About" page for an explanation on how to use the dashboard and details on where information has been sourced from.

    **👈  Click "Target Selection" on the sidebar to get started!**\n
    """,
    unsafe_allow_html=True,
)