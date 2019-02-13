import pickle

results = {}
modelList = ["Mp","D0","Doublet"]
for modelName in modelList:
    picklename = ".tmp.testPlot."+modelName+".pickle"
    
    try:
        with open(picklename,'rb') as f:
            results[modelName] = pickle.load(f)
    except:
        print ("No picklet files of ", modelName)
        exit()

import matplotlib.pyplot as plt
import numpy as np

variableList = [
    "Density",
    "ptThr",
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
    "QuboSize"
]
variables = {}
for m in modelList:
    variables[m] = {}
    for v in variableList:
        variables[m][v] = []

    for k, v in sorted(results[m].items(), key=lambda x: x[1]["density"]):
        #if "5.0e-03" not in k:
        if "020perc" not in k:
            continue
        print(v)
        variables[m]["ptThr"]                   .append(1./(v["ptThr"]*1.2e+3))
        variables[m]["Density"]                 .append(v["density"]*100)
        variables[m]["Nhits"]                   .append(v["meta"]["num_hits"])
        variables[m]["Ntracks"]                 .append(v["meta"]["num_tracks"])
        variables[m]["Ndoublets"]               .append(v["Ndoublets"])
        variables[m]["TReadingHits"]            .append(v["TReadingHits"])
        variables[m]["TInitialDoubletBuilding"] .append(v["TInitialDoubletBuilding"])
        variables[m]["TModelBuilding"]          .append(v["TModelBuilding"])
        variables[m]["TNeal"]                   .append(v["TNeal"])
        variables[m]["diffEnergy"]              .append(v["obsEnergy"]-v['idealEnergy'])
        variables[m]["OccFrac"]                 .append(v["bestOcc"]/float(v["OccSum"]))
        variables[m]["Precision"]               .append(v["precision"]*100)
        variables[m]["Recall"]                  .append(v["recall"]*100)
        variables[m]["Missing"]                 .append(v["missing"])
        variables[m]["TrackmlScore"]            .append(v["trackmlScore"])
        variables[m]["QuboSize"]                .append(v["QuboSize"])


plotVar = {
    ("Ndoublets", "TModelBuilding"),
    ("Density"  , "TModelBuilding"),
    ("Density"  , "TReadingHits"),
    ("Density"  , "TInitialDoubletBuilding"),
    ("Density"  , "TNeal"),
    ("QuboSize" , "TNeal"),
    ("Density"  , "QuboSize"),
    ("Density"  , "Precision"),
    ("Density"  , "Recall"),
    ("Density"  , "Missing"),
    ("Density"  , "diffEnergy"),
    ("Density"  , "Nhits"),
    ("Density"  , "Ntracks"),
    ("Density"  , "Ndoublets"),
    ("Density"  , "TrackmlScore"),
    ("ptThr"    , "TModelBuilding"),
    ("ptThr"    , "TInitialDoubletBuilding"),
    ("ptThr"    , "TNeal"),
    ("ptThr"    , "QuboSize"),
    ("ptThr"    , "Precision"),
    ("ptThr"    , "Recall"),
    ("ptThr"    , "Missing"),
    ("ptThr"    , "diffEnergy"),
    ("ptThr"    , "Ndoublets"),
    ("ptThr"    , "TrackmlScore")
}

legendStr = {
    "Density"       :r"Density [%]",
    "ptThr"         :r"pT threshold [GeV]",
    "Nhits"         :r"N of hits",
    "Ntracks"       :r"N of tracks",
    "Ndoublets"     :r"N of doublets",
    "TReadingHits"  :r"Time of reading hits",
    "TInitialDoubletBuilding":r"Time of initial doublet building",
    "TModelBuilding":r"Time of model building",
    "TNeal"         :r"Time of SimulatedAnnealingSampler",
    "diffEnergy"    :r"Observed energy - ideal energy",
    "OccFrac"       :r"Best sample occurrence / sum of occurences",
    "Precision"     :r"Precision (Purity) [%]",
    "Recall"        :r"Recall (Efficienty) [%]",
    "Missing"       :r"N of missing doublets",
    "QuboSize"      :r"N of QUBO",
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
            plt.plot(x, np.poly1d(p)(x), label=f'{m}, poly2 fit')
        if (
            v == ("Density"  , "TModelBuilding")
            ):
            p = np.polyfit(x, y, 3)
            plt.plot(x, np.poly1d(p)(x), label=f'{m}, poly3 fit')
        if (
            v == ("Density"  , "QuboSize")
            ):
            p = np.polyfit(x, y, 3)
            plt.plot(x, np.poly1d(p)(x), label=f'{m}, poly3 fit')
        if (
            v == ("QuboSize" , "TNeal")
            ):
            p = np.polyfit(x, y, 1)
            plt.plot(x, np.poly1d(p)(x), label=f'{m}, poly1 fit')
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
    plt.savefig("fig_"+v[0]+"_"+v[1]+".png")
    plt.close()
    #plt.show()
