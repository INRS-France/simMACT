# -*- coding: utf-8 -*-
"""GUI management for MACT simulations and analyses."""

# =============================
# DEPENDENCIES
# =============================
from os.path import join, dirname, abspath
import tkinter
import logging
import MSM_Processor


# =============================
# CONSTANTS
# =============================
lgg = logging.getLogger(__name__)
DEF_LOG_CONF = './logConfig/GUI_log.conf'
DEF_LOGGER = "GUI_logger"
lgg = logging.getLogger(__name__)   # if module not run as 'main'

DEF_MDL_DIR = join(dirname(dirname(abspath(__file__))), "Models/")
DEF_MDL_FILE = join(DEF_MDL_DIR, "arm26/arm26.osim")

GUI_COL_SZ = 10  # number of item per column in GUI windows


# ==============================================================================
# GUI functions
# ==============================================================================
# ------------------------------------------------------------------------------
def chooseCoordinateAndMuscleForTVA(prc):
    """GUI dialog for TVA simulation.

    GUI dialog to let the user choose the joint and the muscle considered for
    TVA simulation.

    Input :
    -------
    prc (MSM_Processor) : the musculoskeletal processor

    Return :
    --------
    string : name of the chosen muscle
    string : name of the chosen coordinate
    """
    lgg.debug("Asking muscle and coordinate for TVA...")
    root = tkinter.Tk()
    root.title("Muscle TVA simulation :")

    frm1 = tkinter.LabelFrame(root,
                              text='Select muscle for simulation :',
                              padx=10, pady=10)
    frm1.grid(column=0, row=0, padx=10, pady=10)

    frm2 = tkinter.LabelFrame(root,
                              text='Select coordinate for simulation:',
                              padx=10, pady=10)
    frm2.grid(column=1, row=0, padx=10, pady=10)

    # create radio buttons for muscles and coordinates choice
    cnames = [c.getName() for c in prc.model.getCoordinateSet()]
    mnames = [m.getName() for m in prc.model.getMuscleList()]
    cVar = tkinter.StringVar(root, "")
    mVar = tkinter.StringVar(root, "")

    for i, mn in enumerate(mnames):
        rb = tkinter.Radiobutton(frm1,
                                 text=mn,
                                 variable=mVar,
                                 value=mn,
                                 tristatevalue="x")
        # 'tristate' to prevent buttons to be activated by default
        rb.grid(column=i // GUI_COL_SZ, row=i % GUI_COL_SZ,
                ipadx=10, ipady=10)

    for j, cn in enumerate(cnames):
        rb = tkinter.Radiobutton(frm2,
                                 text=cn,
                                 variable=cVar,
                                 value=cn,
                                 tristatevalue="x")
        # 'tristate' to prevent buttons to be activated by default
        rb.grid(column=j // GUI_COL_SZ, row=j % GUI_COL_SZ,
                ipadx=10, ipady=10)

    bOK = tkinter.Button(root, text="Compute TVA", command=root.destroy)
    bOK.grid(column=0, row=1, padx=10, pady=10)

    # compute the size of the main window
    x_size = ((i+1) // GUI_COL_SZ) + \
             ((j+1) // GUI_COL_SZ)
    # number of columns in frm1 and frm2
    y_size = GUI_COL_SZ             # max number of rows
    margin = 10
    x_scale = 300                    # approximated width  of a column
    y_scale = 50                     # approximated height of a row

    sz = ""
    sz = str(x_size*x_scale + 2*margin) + "x" + \
        str((2+y_size)*y_scale + 2*margin)
    root.geometry(sz)
    root.mainloop()

    return mVar.get(), cVar.get()


# ------------------------------------------------------------------------------
def chooseCoordinateForTVA(prc):
    """GUI dialog to choose the joint considered for TVA simulation.

    Input :
    -------
    prc (MSM_Processor) : the musculoskeletal processor

    Return :
    --------
    string : name of the chosen coordinate
    """
    lgg.debug("Asking coordinate for TVA simulation... ")

    root = tkinter.Tk()
    root.title("Coordinate TVA simulation :")

    frm = tkinter.LabelFrame(
        root, text='Select coordinate for TVA:', padx=10, pady=10)
    frm.pack()

    # create choice radio buttons
    cnames = [c.getName() for c in prc.model.getCoordinateSet()]
    cnVar = tkinter.StringVar(root, "")
    for j, cn in enumerate(cnames):
        rb = tkinter.Radiobutton(frm,
                                 text=cn,
                                 variable=cnVar,
                                 value=cn,
                                 tristatevalue="x")
        # 'tristate' to prevent buttons to be activated by default
        rb.grid(column=j // GUI_COL_SZ, row=j % GUI_COL_SZ,
                ipadx=10, ipady=10)

    # create OK button
    bOK = tkinter.Button(root, text="Compute TVA", command=root.destroy)
    bOK.pack(anchor='c')

    x_size = (j+1) // GUI_COL_SZ
    # number of columns in frm
    y_size = GUI_COL_SZ             # max number of rows
    margin = 10
    x_scale = 300                    # approximated width  of a column
    y_scale = 50                     # approximated height of a row
    sz = ""
    sz = str(x_size*x_scale + 2*margin) + "x" + \
        str((2+y_size)*y_scale + 2*margin)
    root.geometry(sz)
    root.mainloop()

    return cnVar.get()


# =====================================
# MAIN CODE
# =====================================
def main():
    """Test the module."""
    lgg.debug("Loading an Opensim model...")
    prc = MSM_Processor.MSM_Processor(DEF_MDL_FILE)

    # Testing GUI
    chooseCoordinateForTVA(prc)
    lgg.info("Done !")


# =====================================
if __name__ == "__main__":
    # logging stuff
    lgg.setLevel(logging.DEBUG)
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    formatter = logging.Formatter(
        '%(name)s - %(levelname)s - %(filename)s:%(lineno)s\n%(message)s')
    ch.setFormatter(formatter)
    lgg.addHandler(ch)

    main()
