import pickle

results = {}
modelList = ["D0"]

try:
    picklename = ".tmp.checkpT_curv.pickle"
    with open(picklename,'rb') as f:
        results = pickle.load(f)
except:
    print ("No picklet files")
    exit()

import matplotlib.pyplot as plt
import numpy as np

variableList = [
    "Density",
    "ptThr_r",
    "ptThr_w",
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
    
    variables[m]["ptThr_r"]                 .append(0.6/(v["ptThr_r"]*1e+3))
    variables[m]["ptThr_w"]                 .append(v["ptThr_w"])
    variables[m]["Density"]                 .append(v["density"]*100)
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


plotVar = [
    ("ptThr_w"  , "ptThr_r", "TrackmlScore"),
    ("ptThr_w"  , "ptThr_r", "Recall")
]

legendStr = {
    "ptThr_w"       :r"p_{T} threshold for trackML score",
    "ptThr_r"       :r"p_{T} threshold for track reconstruction",
    "Density"       :r"Density [%]",
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
        z = variables[m][v[2]]
        #plt.plot(x, y, 'ro', color=cmap[m], label=m)
        plt.scatter(x, y, c=z)
        #plt.contour(x, y, z)
    plt.xlabel(legendStr[v[0]])
    plt.ylabel(legendStr[v[1]])
    plt.title(f'{legendStr[v[0]]} vs {legendStr[v[1]]}')
    #plt.legend()
    plt.savefig("figure_ptThr/fig_"+v[0]+"_"+v[1]+"_"+v[2]+".png")
    plt.close()
    #plt.show()
