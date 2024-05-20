import matplotlib.pyplot as plt
import mplhep as hep
from coffea.util import load

out = load("hists_ctag_Wc_sf_test/hists_ctag_Wc_sf_test.coffea")
for i in range(4):
    hep.histplot(
        out[list(out.keys())[0]]["weird_MET_phi"][i, :],
        label=out[list(out.keys())[0]]["weird_MET_phi"].axes[0].value(i),
    )
plt.legend()
plt.savefig("all_MET_phi.pdf")
plt.cla()
for i in range(4):
    hep.histplot(
        out[list(out.keys())[0]]["weird_MET_pt"][i, :],
        label=out[list(out.keys())[0]]["weird_MET_pt"].axes[0].value(i),
    )
plt.legend()
plt.savefig("all_MET_pt.pdf")
plt.cla()

print(out[list(out.keys())[0]]["MET_phi"])
hep.histplot(out[list(out.keys())[0]]["MET_phi"][0, sum, :])
plt.savefig("selcted_phi.pdf")
plt.cla()
hep.histplot(out[list(out.keys())[0]]["MET_pt"][0, sum, :])
plt.savefig("selcted_pt.pdf")
