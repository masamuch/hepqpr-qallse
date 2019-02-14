import pickle

results = {}
modelList = ["D0"]
#modelList = ["D0","Doublet"]

try:
    picklename = ".tmp.seedingLoss.pickle"
    with open(picklename,'rb') as f:
        results = pickle.load(f)
except:
    print ("No picklet files")
    exit()

import matplotlib.pyplot as plt
import numpy as np

variableList = [
    "Density",
    "ptThr",
    "nPhi",
    "Nhits",
    "Ntracks",
    "Ndoublets",
    "TReadingHits",
    "TInitialDoubletBuilding",
    "TModelBuilding",
    "TNeal",
    "diffEnergy",
    "OccFrac",
    "Precision",
    "Recall",
    "Missing",
    "TrackmlScore",
    "QuboSize",
    "precision_initDoublet",
    "recall_initDoublet",
    "missing_initDoublet",
    "TQuboBuilding"
]
variables = {}
for k, v in sorted(results.items(), key=lambda x: x[1]["density"]):
    m = k[0]
    #i = k[3]
    if m not in variables:
        variables[m] = {}
        for v1 in variableList:
            variables[m][v1] = []
    
    variables[m]["ptThr"]                   .append(1./(v["ptThr"]*1.2e+3))
    variables[m]["Density"]                 .append(v["density"]*100)
    variables[m]["nPhi"]                    .append(k[3])
    variables[m]["Nhits"]                   .append(v["meta"]["num_hits"])
    variables[m]["Ntracks"]                 .append(v["meta"]["num_tracks"])
    variables[m]["Ndoublets"]               .append(v["Ndoublets"])
    variables[m]["TReadingHits"]            .append(v["TReadingHits"])
    variables[m]["TInitialDoubletBuilding"] .append(v["TInitialDoubletBuilding"])
    variables[m]["TModelBuilding"]          .append(v["TModelBuilding"])
    variables[m]["TNeal"]                   .append(v["TNeal"])
    variables[m]["diffEnergy"]              .append((v["obsEnergy"]-v['idealEnergy'])/abs(v['idealEnergy']))
    variables[m]["OccFrac"]                 .append(v["bestOcc"]/float(v["OccSum"]))
    variables[m]["Precision"]               .append(v["precision"]*100)
    variables[m]["Recall"]                  .append(v["recall"]*100)
    variables[m]["Missing"]                 .append(v["missing"])
    variables[m]["TrackmlScore"]            .append(v["trackmlScore"])
    variables[m]["QuboSize"]                .append(v["QuboSize"])
    variables[m]["TQuboBuilding"]           .append(v["TQuboBuilding"])
    variables[m]["precision_initDoublet"]   .append(v["precision_initDoublet"])
    variables[m]["recall_initDoublet"]      .append(v["recall_initDoublet"])
    variables[m]["missing_initDoublet"]     .append(v["missing_initDoublet"])


plotVar = {
    ("nPhi"  , "TModelBuilding"),
    ("nPhi"  , "TReadingHits"),
    ("nPhi"  , "TInitialDoubletBuilding"),
    ("nPhi"  , "TNeal"),
    ("nPhi"  , "QuboSize"),
    ("nPhi"  , "Precision"),
    ("nPhi"  , "Recall"),
    ("nPhi"  , "Missing"),
    ("nPhi"  , "diffEnergy"),
    ("nPhi"  , "Nhits"),
    ("nPhi"  , "Ntracks"),
    ("nPhi"  , "Ndoublets"),
    ("nPhi"  , "TrackmlScore"),
    ("nPhi"  , "TQuboBuilding"),
    ("nPhi"  , "precision_initDoublet"),
    ("nPhi"  , "recall_initDoublet"),
    ("nPhi"  , "missing_initDoublet"),
}

legendStr = {
    "Density"       :r"Density [%]",
    "nPhi"          :r"N of #phi slices",
    "ptThr"         :r"pT threshold [GeV]",
    "Nhits"         :r"N of hits",
    "Ntracks"       :r"N of tracks",
    "Ndoublets"     :r"N of doublets",
    "TReadingHits"  :r"Time of reading hits",
    "TInitialDoubletBuilding":r"Time of initial doublet building",
    "TModelBuilding":r"Time of model building",
    "TNeal"         :r"Time of SimulatedAnnealingSampler",
    "diffEnergy"    :r"(Observed energy - ideal energy)/|ideal energy|",
    "OccFrac"       :r"Best sample occurrence / sum of occurences",
    "Precision"     :r"Precision (Purity) [%]",
    "Recall"        :r"Recall (Efficienty) [%]",
    "Missing"       :r"N of missing doublets",
    "QuboSize"      :r"N of QUBO",
    "TQuboBuilding" :r"Time of QUBO building",
    "precision_initDoublet":r"Precision (Purity) in initial doublet creation [%]",
    "recall_initDoublet"   :r"Recall (Efficienccy) in initial doublet creation [%]",
    "missing_initDoublet"  :r"N of missing doublets in initial doublet creation",
    "TrackmlScore"  :r"TrackML Score"
}

cmap0 = plt.get_cmap("tab10")
cmap = {
    "Mp"     :cmap0(0),
    "D0"     :cmap0(1),
    "Doublet":cmap0(2)
}
for v in plotVar:
    for m in modelList:
        x = variables[m][v[0]]
        y = variables[m][v[1]]
        plt.plot(x, y, 'ro', color=cmap[m], label=m)
        if (
            v == ("Density"  , "Ndoublets")
         or v == ("Density"  , "TInitialDoubletBuilding")
            ):
            p = np.polyfit(x, y, 2)
            plt.plot(x, np.poly1d(p)(x), color=cmap[m], label=f'{m}, poly2 fit')
        if (
            v == ("Density"  , "TModelBuilding")
            ):
            p = np.polyfit(x, y, 3)
            plt.plot(x, np.poly1d(p)(x), color=cmap[m], label=f'{m}, poly3 fit')
        if (
            v == ("Density"  , "QuboSize")
            ):
            p = np.polyfit(x, y, 6)
            plt.plot(x, np.poly1d(p)(x), color=cmap[m], label=f'{m}, poly6 fit')
        if (
            v == ("QuboSize" , "TNeal")
            ):
            p = np.polyfit(x, y, 1)
            plt.plot(x, np.poly1d(p)(x), color=cmap[m], label=f'{m}, poly1 fit')
        #if (
        #    v == ("Density"  , "diffEnergy")
        #    ):
        #    p = np.polyfit(x, y, 2)
        #    plt.plot(x, np.poly1d(p)(x), label=f'{m}, poly2 fit')
        #plt.plot(variables[v[0]], np.poly1d(np.polyfit(variables[v[0]], variables[v[1]], 3))(variables[v[0]]), label='poly3')
        #plt.plot(variables[v[0]], np.poly1d(np.polyfit(variables[v[0]], variables[v[1]], 4))(variables[v[0]]), label='poly4')
    plt.xlabel(legendStr[v[0]])
    plt.ylabel(legendStr[v[1]])
    plt.title(f'{legendStr[v[0]]} vs {legendStr[v[1]]}')
    plt.legend()
    plt.savefig("figure_nPhi/fig_"+v[0]+"_"+v[1]+".png")
    plt.close()
    #plt.show()
