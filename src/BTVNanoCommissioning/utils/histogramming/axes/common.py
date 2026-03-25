import hist

axes = {
    "flav": hist.axis.IntCategory([0, 1, 4, 5, 6], name="flav", label="Genflavour"),
    "syst": hist.axis.StrCategory([], name="syst", growth=True),
    "pt": hist.axis.Regular(60, 0, 300, name="pt", label=" $p_{T}$ [GeV]"),
    "jpt": hist.axis.Regular(300, 0, 3000, name="pt", label=" $p_{T}$ [GeV]"),
    "ljpt": hist.axis.Regular(200, 0, 2000, name="pt", label=" $p_{T}$ [GeV]"),
    "softlpt": hist.axis.Regular(25, 0, 25, name="pt", label=" $p_{T}$ [GeV]"),
    "mass": hist.axis.Regular(50, 0, 300, name="mass", label=" mass [GeV]"),
    "bdt": hist.axis.Regular(50, 0, 1, name="bdt", label=" BDT discriminant"),
    "eta": hist.axis.Regular(25, -2.5, 2.5, name="eta", label=" $\eta$"),
    "phi": hist.axis.Regular(30, -3, 3, name="phi", label="$\phi$"),
    "mt": hist.axis.Regular(30, 0, 300, name="mt", label=" $m_{T}$ [GeV]"),
    "iso": hist.axis.Regular(30, 0, 0.05, name="pfRelIso03_all", label="Rel. Iso"),
    "softliso": hist.axis.Regular(
        20, 0.2, 6.2, name="pfRelIso03_all", label="Rel. Iso"
    ),
    "npv": hist.axis.Integer(0, 100, name="npv", label="N PVs"),
    "dr": hist.axis.Regular(20, 0, 8, name="dr", label="$\Delta$R"),
    "dr_s": hist.axis.Regular(20, 0, 0.5, name="dr", label="$\Delta$R"),
    "dr_SV": hist.axis.Regular(20, 0, 1.0, name="dr", label="$\Delta$R"),
    "dxy": hist.axis.Regular(40, -0.05, 0.05, name="dxy", label="$d_{xy}$ [cm]"),
    "dz": hist.axis.Regular(40, -0.01, 0.01, name="dz", label="$d_{z}$ [cm]"),
    "qcddxy": hist.axis.Regular(40, -0.002, 0.002, name="dxy", label="$d_{xy}$ [cm]"),
    "sip3d": hist.axis.Regular(20, 0, 0.2, name="sip3d", label="SIP 3D"),
    "ptratio": hist.axis.Regular(50, 0, 1, name="ratio", label="ratio"),
    "n": hist.axis.Integer(0, 10, name="n", label="N obj"),
    "osss": hist.axis.IntCategory([1, -1], name="osss", label="OS(+)/SS(-)"),
}
